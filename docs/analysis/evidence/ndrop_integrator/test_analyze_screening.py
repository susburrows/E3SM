"""Short synthetic tests for the screening-analysis core."""

from __future__ import annotations

import math
import unittest

from analyze_screening import (
    BreakEvenStatus,
    ExclusionReason,
    SelectionStatus,
    ScreeningInputError,
    bootstrap_median_interval,
    deterministic_groups,
    ingest_rows,
    interpolate_break_even,
    select_measured_targets,
    summarize_timing,
)


def _valid_row(**updates: object) -> dict:
    row = {
        "case_id": "7",
        "nlev": "3",
        "transported_fields": "9",
        "batch_size": "16",
        "mixing_index": "1",
        "mixing_scale": "0.5",
        "dt": "50",
        "method": "imex",
        "forced_fe_substeps": "0",
        "validation_ok": "true",
        "eligible": "true",
        "richardson_ok": "true",
        "screening_error": "0.003",
        "median_seconds_per_column": "2",
        "ci95_lower": "1.5",
        "ci95_upper": "2.5",
        "clipping_count": "0",
        "material_negative_count": "0",
        "repetitions": "3",
        "raw_samples_seconds_per_column": "1;3;2",
    }
    row.update(updates)
    return row


def _constant_timing_row(case_id: int, error: float, seconds: float) -> dict:
    return _valid_row(
        case_id=str(case_id),
        screening_error=str(error),
        median_seconds_per_column=str(seconds),
        ci95_lower=str(seconds),
        ci95_upper=str(seconds),
        raw_samples_seconds_per_column="{0};{0};{0}".format(seconds),
    )


def _timing_row(
    case_id: int,
    error: float,
    median: float,
    lower: float,
    upper: float,
    method: str = "imex",
    fe_nsubmix: int = 0,
) -> dict:
    return _valid_row(
        case_id=str(case_id),
        screening_error=str(error),
        method=method,
        forced_fe_substeps=str(fe_nsubmix),
        median_seconds_per_column=str(median),
        ci95_lower=str(lower),
        ci95_upper=str(upper),
        raw_samples_seconds_per_column="{0};{0};{0}".format(median),
    )


class TimingSummaryTests(unittest.TestCase):
    def test_bootstrap_matches_hand_derived_splitmix_fixture(self) -> None:
        # With seed 0, the four 3-draw index sets modulo 3 are 101, 110,
        # 222, and 211.  Their medians are 2, 2, 3, and 2, so the exact
        # percentile bounds are [2, 3] and the relative half-width is 0.25.
        interval = bootstrap_median_interval([1.0, 2.0, 3.0], 4, 0)

        self.assertEqual(interval.median_seconds_per_column, 2.0)
        self.assertEqual(interval.ci95_lower_seconds_per_column, 2.0)
        self.assertEqual(interval.ci95_upper_seconds_per_column, 3.0)
        self.assertEqual(interval.ci95_relative_half_width, 0.25)

    def test_even_sample_median_uses_supplied_confidence_interval(self) -> None:
        # Catches selecting a middle element or recomputing/replacing the CI.
        summary = summarize_timing([4.0, 1.0, 3.0, 2.0], (1.5, 3.5))

        self.assertEqual(summary.median_seconds_per_column, 2.5)
        self.assertEqual(summary.ci95_lower_seconds_per_column, 1.5)
        self.assertEqual(summary.ci95_upper_seconds_per_column, 3.5)

    def test_timing_shape_and_finiteness_fail_closed(self) -> None:
        # Catches accepting empty, nonfinite, negative, or inverted timing data.
        bad_inputs = (
            (None, (0.0, 1.0)),
            ([], (0.0, 1.0)),
            ([1.0, math.inf], (0.0, 1.0)),
            ([-1.0], (0.0, 1.0)),
            ([1.0], (2.0, 1.0)),
            ([1.0], (0.0,)),
            ([1.0], None),
        )
        for samples, interval in bad_inputs:
            with self.subTest(samples=samples, interval=interval):
                with self.assertRaises(ScreeningInputError):
                    summarize_timing(samples, interval)  # type: ignore[arg-type]


class RowIngestionTests(unittest.TestCase):
    def test_validation_failure_is_excluded_not_treated_as_malformed(self) -> None:
        # Catches attempting to analyze a failed row or requiring absent timings.
        row = _valid_row(
            validation_ok="false",
            eligible="false",
            richardson_ok="false",
            screening_error="",
            median_seconds_per_column="",
            ci95_lower="",
            ci95_upper="",
            repetitions="0",
            raw_samples_seconds_per_column="",
        )

        dataset = ingest_rows([row])

        self.assertEqual(dataset.points, ())
        self.assertEqual(len(dataset.exclusions), 1)
        self.assertEqual(
            dataset.exclusions[0].reason, ExclusionReason.VALIDATION_FAILURE
        )

    def test_each_analysis_gate_has_a_specific_exclusion(self) -> None:
        # Catches admitting unsafe points or hiding their specific gate failure.
        rows = [
            _valid_row(case_id="1", eligible="false", clipping_count="1"),
            _valid_row(case_id="2", material_negative_count="1"),
            _valid_row(case_id="3", eligible="false"),
            _valid_row(
                case_id="4", richardson_ok="false", screening_error=""
            ),
        ]

        dataset = ingest_rows(rows)

        self.assertEqual(dataset.points, ())
        self.assertEqual(
            tuple(item.reason for item in dataset.exclusions),
            (
                ExclusionReason.CLIPPING,
                ExclusionReason.MATERIAL_NEGATIVITY,
                ExclusionReason.INELIGIBLE,
                ExclusionReason.MISSING_RICHARDSON_ERROR,
            ),
        )

    def test_included_row_uses_exact_raw_median_and_supplied_ci(self) -> None:
        # Catches trusting an approximate statistic or discarding bootstrap bounds.
        dataset = ingest_rows([_valid_row()])

        self.assertEqual(dataset.exclusions, ())
        self.assertEqual(len(dataset.points), 1)
        point = dataset.points[0]
        self.assertEqual(point.screening_error, 0.003)
        self.assertEqual(point.timing.median_seconds_per_column, 2.0)
        self.assertEqual(
            (
                point.timing.ci95_lower_seconds_per_column,
                point.timing.ci95_upper_seconds_per_column,
            ),
            (1.5, 2.5),
        )

    def test_reported_count_and_median_must_match_raw_samples(self) -> None:
        # Catches silently accepting truncated samples or a stale reported median.
        bad_rows = (
            _valid_row(repetitions="4"),
            _valid_row(median_seconds_per_column="2.1"),
            _valid_row(
                validation_ok="false",
                eligible="false",
                richardson_ok="false",
                screening_error="",
                repetitions="1",
                raw_samples_seconds_per_column="",
            ),
        )
        for row in bad_rows:
            with self.subTest(row=row):
                with self.assertRaises(ScreeningInputError):
                    ingest_rows([row])

    def test_top_level_shape_and_required_schema_fail_closed(self) -> None:
        # Catches leaking TypeError or accepting a partially decoded CSV schema.
        missing_dt = _valid_row()
        del missing_dt["dt"]
        for rows in (None, [missing_dt]):
            with self.subTest(rows=rows):
                with self.assertRaises(ScreeningInputError):
                    ingest_rows(rows)  # type: ignore[arg-type]

    def test_richardson_error_must_be_a_finite_scalar(self) -> None:
        # Catches admitting NaN or leaking a TypeError for a nonscalar value.
        for value in ("nan", []):
            with self.subTest(value=value):
                with self.assertRaises(ScreeningInputError):
                    ingest_rows([_valid_row(screening_error=value)])


class GroupingTests(unittest.TestCase):
    def test_group_and_point_order_does_not_depend_on_input_order(self) -> None:
        # Catches nondeterministic reports caused by source-row ordering.
        rows = [
            _valid_row(case_id="9", nlev="3"),
            _valid_row(case_id="5", nlev="1"),
            _valid_row(
                case_id="7",
                nlev="3",
                method="forward_euler",
                forced_fe_substeps="5",
            ),
            _valid_row(
                case_id="6",
                nlev="3",
                method="forward_euler",
                forced_fe_substeps="2",
            ),
        ]

        forward = deterministic_groups(ingest_rows(rows).points)
        backward = deterministic_groups(ingest_rows(reversed(rows)).points)

        self.assertEqual(forward, backward)
        self.assertEqual(tuple(group.nlev for group, _ in forward), (1, 3))
        self.assertEqual(
            tuple(point.case_id for point in forward[1][1]), (6, 7, 9)
        )

    def test_group_key_includes_numeric_mixing_strength(self) -> None:
        # Catches merging workloads whose index label is stale or reused.
        points = ingest_rows(
            [
                _valid_row(case_id="1", mixing_scale="0.05"),
                _valid_row(case_id="2", mixing_scale="0.5"),
            ]
        ).points

        groups = deterministic_groups(points)

        self.assertEqual(
            tuple(group.mixing_scale for group, _ in groups), (0.05, 0.5)
        )


class TargetSelectionTests(unittest.TestCase):
    def test_targets_choose_only_fastest_qualifying_measured_point(self) -> None:
        # Catches interpolation/extrapolation on error or choosing a slower point.
        points = ingest_rows(
            [
                _constant_timing_row(1, 0.01, 1.0),
                _constant_timing_row(2, 0.003, 2.0),
                _constant_timing_row(3, 0.001, 4.0),
            ]
        ).points

        selections = select_measured_targets(points, [0.005, 0.0005, 0.02])

        self.assertEqual(
            tuple(item.target_error for item in selections),
            (0.0005, 0.005, 0.02),
        )
        self.assertEqual(selections[0].status, SelectionStatus.UNAVAILABLE)
        self.assertIsNone(selections[0].point)
        self.assertEqual(selections[1].point.case_id, 2)  # type: ignore[union-attr]
        self.assertEqual(selections[2].point.case_id, 1)  # type: ignore[union-attr]

    def test_fe_nsubmix_values_remain_separate_measured_axes(self) -> None:
        # Catches collapsing the FE curve to one fastest nsubmix before break-even.
        points = ingest_rows(
            [
                _timing_row(1, 0.003, 2.0, 2.0, 2.0, "forward_euler", 1),
                _timing_row(2, 0.003, 6.0, 6.0, 6.0, "forward_euler", 3),
            ]
        ).points

        selections = select_measured_targets(points, [0.005])

        self.assertEqual(
            tuple(item.fe_nsubmix for item in selections), (1, 3)
        )
        self.assertEqual(
            tuple(item.point.case_id for item in selections),
            (1, 2),
        )


class BreakEvenTests(unittest.TestCase):
    def test_adjacent_interpolation_carries_candidate_and_fe_timing_cis(self) -> None:
        # Catches non-adjacent fits or dropping either timing uncertainty source.
        dataset = ingest_rows(
            [
                _timing_row(1, 0.003, 2.0, 1.8, 2.2, "forward_euler", 1),
                _timing_row(2, 0.003, 6.0, 5.4, 6.6, "forward_euler", 3),
                _timing_row(3, 0.003, 10.0, 9.0, 11.0, "forward_euler", 5),
                _timing_row(4, 0.003, 8.0, 7.0, 9.0),
            ]
        )
        candidate = next(
            point for point in dataset.points if point.method == "imex"
        )
        fe_points = [
            point for point in dataset.points if point.method == "forward_euler"
        ]

        result = interpolate_break_even(fe_points, candidate)

        self.assertEqual(result.status, BreakEvenStatus.AVAILABLE)
        self.assertEqual(result.first_measured_fe_nsubmix, 5)
        self.assertEqual(result.adjacent_bracket, (3, 5))
        self.assertAlmostEqual(result.median_fe_nsubmix, 4.0)
        self.assertAlmostEqual(result.ci95_lower_fe_nsubmix, 3.0 + 2.0 / 11.0)
        self.assertAlmostEqual(result.ci95_upper_fe_nsubmix, 5.0)

    def test_missing_candidate_or_fe_curve_is_unavailable(self) -> None:
        # Catches presenting an estimate when one side has no measured point.
        candidate = ingest_rows([_constant_timing_row(1, 0.003, 8.0)]).points[0]

        no_fe = interpolate_break_even([], candidate)
        no_candidate = interpolate_break_even([], None)

        self.assertEqual(no_fe.status, BreakEvenStatus.UNAVAILABLE)
        self.assertEqual(no_candidate.status, BreakEvenStatus.UNAVAILABLE)
        self.assertIsNone(no_fe.median_fe_nsubmix)
        self.assertIsNone(no_candidate.median_fe_nsubmix)

    def test_measured_fe_curve_below_candidate_reports_no_crossing(self) -> None:
        # Catches extrapolating a crossing beyond the largest measured FE nsubmix.
        dataset = ingest_rows(
            [
                _timing_row(1, 0.003, 2.0, 1.8, 2.2, "forward_euler", 1),
                _timing_row(2, 0.003, 6.0, 5.4, 6.6, "forward_euler", 3),
                _timing_row(3, 0.003, 8.0, 7.0, 9.0),
            ]
        )
        candidate = next(
            point for point in dataset.points if point.method == "imex"
        )
        fe_points = [
            point for point in dataset.points if point.method == "forward_euler"
        ]

        result = interpolate_break_even(fe_points, candidate)

        self.assertEqual(result.status, BreakEvenStatus.NO_CROSSING)
        self.assertIsNone(result.median_fe_nsubmix)
        self.assertIsNone(result.adjacent_bracket)

    def test_exact_measured_crossing_needs_no_fractional_interpolation(self) -> None:
        dataset = ingest_rows(
            [
                _timing_row(1, 0.003, 2.0, 2.0, 2.0, "forward_euler", 1),
                _timing_row(2, 0.003, 6.0, 6.0, 6.0, "forward_euler", 3),
                _timing_row(3, 0.003, 6.0, 6.0, 6.0),
            ]
        )
        candidate = next(
            point for point in dataset.points if point.method == "imex"
        )
        fe_points = [
            point for point in dataset.points if point.method == "forward_euler"
        ]

        result = interpolate_break_even(fe_points, candidate)

        self.assertEqual(result.status, BreakEvenStatus.AVAILABLE)
        self.assertEqual(result.median_fe_nsubmix, 3.0)
        self.assertEqual(result.adjacent_bracket, (3, 3))

    def test_break_even_rejects_mixed_aci_intervals(self) -> None:
        dataset = ingest_rows(
            [
                _timing_row(1, 0.003, 2.0, 2.0, 2.0, "forward_euler", 1),
                _valid_row(
                    case_id="2",
                    dt="75",
                    median_seconds_per_column="2",
                    ci95_lower="2",
                    ci95_upper="2",
                    raw_samples_seconds_per_column="2;2;2",
                ),
            ]
        )
        candidate = next(
            point for point in dataset.points if point.method == "imex"
        )
        fe_points = [
            point for point in dataset.points if point.method == "forward_euler"
        ]

        with self.assertRaises(ScreeningInputError):
            interpolate_break_even(fe_points, candidate)

    def test_malformed_candidate_fails_closed(self) -> None:
        # Catches leaking AttributeError before the break-even shape guard.
        with self.assertRaises(ScreeningInputError):
            interpolate_break_even([], object())  # type: ignore[arg-type]

        fe_point = ingest_rows(
            [_timing_row(1, 0.003, 2.0, 2.0, 2.0, "forward_euler", 1)]
        ).points[0]
        with self.assertRaises(ScreeningInputError):
            interpolate_break_even([fe_point], fe_point)


if __name__ == "__main__":
    unittest.main()
