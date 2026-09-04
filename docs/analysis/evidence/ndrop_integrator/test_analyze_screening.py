"""Short synthetic tests for the screening-analysis core."""

from __future__ import annotations

import math
import unittest

import analyze_screening as screening
from analyze_screening import (
    BreakEvenStatus,
    ExclusionReason,
    PhaseDiagnostics,
    SelectionStatus,
    ScreeningInputError,
    TimingEstimate,
    bootstrap_median_interval,
    deterministic_groups,
    ingest_paired_artifacts,
    ingest_rows,
    interpolate_break_even,
    select_measured_targets,
    summarize_timing,
)


def _valid_row(**updates: object) -> dict:
    row = {
        "artifact_schema": "ndrop_summary_v2",
        "run_id": "task5-remediated-20260825T120000Z-61ae8249af0e",
        "case_id": "7",
        "nlev": "3",
        "nonzero_aerosol_pairs": "4",
        "nonzero_state_fields": "9",
        "executed_aerosol_pairs": "25",
        "executed_transport_fields": "51",
        "batch_size": "16",
        "mixing_index": "1",
        "mixing_scale": "0.5",
        "dt": "50",
        "method": "imex",
        "linear_phase_mode": "combined",
        "forced_fe_substeps": "0",
        "validation_ok": "true",
        "eligible": "true",
        "richardson_ok": "true",
        "screening_error": "0.003",
        "median_seconds_per_column": "2",
        "ci95_lower": "2",
        "ci95_upper": "2",
        "ci95_relative_half_width": "0",
        "clipping_count": "0",
        "material_negative_count": "0",
        "repetitions": "3",
        "raw_samples_seconds_per_column": "2;2;2",
        "instrumented_samples_seconds_per_column": "2;3;1",
        "instrumented_median_seconds_per_column": "2",
        "instrumented_ci95_lower": "1",
        "instrumented_ci95_upper": "3",
        "instrumentation_overhead": "0",
        "stop_reason": "ci_target",
        "minimum_repetitions": "3",
        "maximum_repetitions": "3",
        "bootstrap_resamples": "4",
        "setup_seconds_per_column": "0.1",
        "factorization_seconds_per_column": "",
        "solve_seconds_per_column": "",
        "combined_factorization_solve_seconds_per_column": "0.6",
        "update_seconds_per_column": "0.2",
        "diagnostics_seconds_per_column": "0.3",
        "fence_seconds_per_column": "0.05",
    }
    row.update(updates)
    if "linear_phase_mode" not in updates:
        if row["method"] == "forward_euler":
            row["linear_phase_mode"] = "none"
        elif row["method"] == "backward_euler":
            row["linear_phase_mode"] = "separate"
        else:
            row["linear_phase_mode"] = "combined"
    linear_updates = {
        "factorization_seconds_per_column",
        "solve_seconds_per_column",
        "combined_factorization_solve_seconds_per_column",
    }
    if not linear_updates.intersection(updates):
        if row["method"] == "forward_euler":
            row["factorization_seconds_per_column"] = ""
            row["solve_seconds_per_column"] = ""
            row["combined_factorization_solve_seconds_per_column"] = ""
        elif row["method"] == "backward_euler":
            row["factorization_seconds_per_column"] = "0.25"
            row["solve_seconds_per_column"] = "0.35"
            row["combined_factorization_solve_seconds_per_column"] = ""
    if row["repetitions"] == "0":
        for name in (
            "instrumented_samples_seconds_per_column",
            "instrumented_median_seconds_per_column",
            "instrumented_ci95_lower",
            "instrumented_ci95_upper",
            "instrumentation_overhead",
            "ci95_relative_half_width",
            "stop_reason",
            "setup_seconds_per_column",
            "factorization_seconds_per_column",
            "solve_seconds_per_column",
            "combined_factorization_solve_seconds_per_column",
            "update_seconds_per_column",
            "diagnostics_seconds_per_column",
            "fence_seconds_per_column",
        ):
            if name not in updates:
                row[name] = ""
        if "stop_reason" not in updates:
            row["stop_reason"] = "invalid"
    return row


def _phase_rows(summary_row: dict) -> list[dict]:
    repetitions = int(summary_row["repetitions"])
    if repetitions == 0:
        return []
    instrumented = str(
        summary_row["instrumented_samples_seconds_per_column"]
    ).split(";")
    setup = ("0.1", "0.2", "0")
    factorization = ("", "", "")
    solve = ("", "", "")
    combined = ("0.6", "0.7", "0.5")
    if summary_row["method"] == "forward_euler":
        combined = ("", "", "")
    elif summary_row["method"] == "backward_euler":
        factorization = ("0.25", "0.3", "0.2")
        solve = ("0.35", "0.4", "0.3")
        combined = ("", "", "")
    return [
        {
            "artifact_schema": "ndrop_phase_samples_v1",
            "run_id": summary_row["run_id"],
            "case_id": summary_row["case_id"],
            "method": summary_row["method"],
            "linear_phase_mode": summary_row["linear_phase_mode"],
            "instrumented_repetition": str(repetition),
            "instrumented_seconds_per_column": instrumented[repetition],
            "setup_seconds_per_column": setup[repetition],
            "factorization_seconds_per_column": factorization[repetition],
            "solve_seconds_per_column": solve[repetition],
            "combined_factorization_solve_seconds_per_column": combined[
                repetition
            ],
            "update_seconds_per_column": ("0.2", "0.3", "0.1")[
                repetition
            ],
            "diagnostics_seconds_per_column": ("0.3", "0.4", "0.2")[
                repetition
            ],
            "fence_seconds_per_column": ("0.05", "0.06", "0.04")[
                repetition
            ],
        }
        for repetition in range(repetitions)
    ]


def _constant_timing_row(case_id: int, error: float, seconds: float) -> dict:
    return _valid_row(
        case_id=str(case_id),
        dt=str(25 + 25 * case_id),
        screening_error=str(error),
        median_seconds_per_column=str(seconds),
        ci95_lower=str(seconds),
        ci95_upper=str(seconds),
        ci95_relative_half_width="0",
        raw_samples_seconds_per_column="{0};{0};{0}".format(seconds),
        instrumentation_overhead=str(2.0 / seconds - 1.0),
        stop_reason="ci_target",
    )


def _timing_row(
    case_id: int,
    error: float,
    median: float,
    lower: float,
    upper: float,
    method: str = "imex",
    fe_nsubmix: int = 0,
    dt: float = 50.0,
) -> dict:
    relative_half_width = (upper - lower) / 2.0 / median
    return _valid_row(
        case_id=str(case_id),
        dt=str(dt),
        screening_error=str(error),
        method=method,
        forced_fe_substeps=str(fe_nsubmix),
        median_seconds_per_column=str(median),
        ci95_lower=str(lower),
        ci95_upper=str(upper),
        ci95_relative_half_width=str(relative_half_width),
        raw_samples_seconds_per_column="{};{};{}".format(lower, median, upper),
        bootstrap_resamples="100",
        instrumentation_overhead=str(2.0 / median - 1.0),
        stop_reason=(
            "ci_target"
            if relative_half_width < 0.05
            else "maximum_repetitions"
        ),
    )


def _timing_policy_row(
    samples: tuple[float, ...],
    *,
    median: float,
    lower: float,
    upper: float,
    relative_half_width: float,
    stop_reason: str,
    minimum_repetitions: int,
    maximum_repetitions: int,
    resamples: int = 4,
) -> dict:
    joined = ";".join(str(value) for value in samples)
    return _valid_row(
        case_id=str(0x4E44524F),
        repetitions=str(len(samples)),
        raw_samples_seconds_per_column=joined,
        instrumented_samples_seconds_per_column=joined,
        median_seconds_per_column=str(median),
        ci95_lower=str(lower),
        ci95_upper=str(upper),
        ci95_relative_half_width=str(relative_half_width),
        instrumented_median_seconds_per_column=str(median),
        instrumented_ci95_lower=str(lower),
        instrumented_ci95_upper=str(upper),
        instrumentation_overhead="0",
        stop_reason=stop_reason,
        minimum_repetitions=str(minimum_repetitions),
        maximum_repetitions=str(maximum_repetitions),
        bootstrap_resamples=str(resamples),
    )


class TimingSummaryTests(unittest.TestCase):
    def test_even_median_exactly_matches_cpp_overflow_safe_formula(self) -> None:
        lower = 4.154765028703989e-220
        upper = 3.588780579606149e-219
        cpp_median = lower + (upper - lower) / 2.0

        summary = summarize_timing((lower, upper), (cpp_median, cpp_median))

        self.assertEqual(summary.median_seconds_per_column, cpp_median)

    def test_bootstrap_matches_hand_derived_splitmix_fixture(self) -> None:
        case_id = 0x4E44524F
        # With seed 0, the four 3-draw index sets modulo 3 are 101, 110,
        # 222, and 211.  Their medians are 2, 2, 3, and 2, so the exact
        # percentile bounds are [2, 3] and the relative half-width is 0.25.
        interval = bootstrap_median_interval(
            [1.0, 2.0, 3.0], 4, 0x4E44524F ^ case_id
        )

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
    def test_reconstructs_primary_bootstrap_and_first_ci_target_prefix(self) -> None:
        # Seed 0 gives prefix-3 draws 101, 110, 222, 211: CI [1, 2].
        # Prefix-4 draws 3030, 3210, 3212, 3313: every median is 1.
        row = _timing_policy_row(
            (1.0, 2.0, 1.0, 1.0),
            median=1.0,
            lower=1.0,
            upper=1.0,
            relative_half_width=0.0,
            stop_reason="ci_target",
            minimum_repetitions=3,
            maximum_repetitions=5,
        )

        point = ingest_rows([row]).points[0]

        self.assertEqual(point.timing, TimingEstimate(1.0, 1.0, 1.0))

    def test_maximum_repetitions_replays_every_eligible_prefix(self) -> None:
        row = _timing_policy_row(
            (1.0, 2.0, 3.0, 2.0, 3.0),
            median=2.0,
            lower=1.0,
            upper=3.0,
            relative_half_width=0.5,
            stop_reason="maximum_repetitions",
            minimum_repetitions=3,
            maximum_repetitions=5,
        )

        point = ingest_rows([row]).points[0]

        self.assertEqual(point.timing, TimingEstimate(2.0, 1.0, 3.0))

    def test_stale_primary_timing_policy_fields_fail_closed(self) -> None:
        base = _timing_policy_row(
            (1.0, 2.0, 3.0, 2.0, 3.0),
            median=2.0,
            lower=1.0,
            upper=3.0,
            relative_half_width=0.5,
            stop_reason="maximum_repetitions",
            minimum_repetitions=3,
            maximum_repetitions=5,
        )
        mutations = (
            {"median_seconds_per_column": "2.1"},
            {
                "ci95_lower": "1.5",
                "ci95_upper": "2.5",
                "ci95_relative_half_width": "0.25",
            },
            {"ci95_relative_half_width": "0.4"},
            {"repetitions": "4"},
            {"bootstrap_resamples": "1"},
            {"minimum_repetitions": "6"},
            {"maximum_repetitions": "6"},
            {"stop_reason": "unexpected"},
        )
        for mutation in mutations:
            row = dict(base)
            row.update(mutation)
            with self.subTest(mutation=mutation):
                with self.assertRaises(ScreeningInputError):
                    ingest_rows([row])

    def test_stopping_after_an_earlier_qualifying_prefix_fails_closed(self) -> None:
        # Prefix 4 has CI [1, 1], so the producer cannot collect sample 5.
        row = _timing_policy_row(
            (1.0, 2.0, 1.0, 1.0, 1.0),
            median=1.0,
            lower=1.0,
            upper=2.0,
            relative_half_width=0.5,
            stop_reason="maximum_repetitions",
            minimum_repetitions=3,
            maximum_repetitions=5,
        )

        with self.assertRaises(ScreeningInputError):
            ingest_rows([row])

    def test_all_timing_policy_fields_fail_closed(self) -> None:
        malformed = (
            _valid_row(ci95_relative_half_width="0.4"),
            _valid_row(instrumented_ci95_lower="2.1"),
            _valid_row(instrumented_ci95_upper="nan"),
            _valid_row(instrumentation_overhead="0.1"),
            _valid_row(stop_reason="unexpected"),
            _valid_row(stop_reason="maximum_repetitions"),
            _valid_row(
                stop_reason="maximum_repetitions",
                ci95_lower="1.92",
                ci95_upper="2.08",
                ci95_relative_half_width="0.04",
            ),
            _valid_row(
                validation_ok="false",
                eligible="false",
                richardson_ok="false",
                screening_error="",
                median_seconds_per_column="",
                ci95_lower="",
                ci95_upper="",
                repetitions="0",
                raw_samples_seconds_per_column="",
                instrumentation_overhead="0",
            ),
        )
        for row in malformed:
            with self.subTest(row=row):
                with self.assertRaises(ScreeningInputError):
                    ingest_rows([row])

    def test_v2_schema_width_and_phase_nullability_fail_closed(self) -> None:
        bad_rows = (
            _valid_row(artifact_schema="ndrop_summary_v1"),
            _valid_row(run_id="bad/run"),
            _valid_row(nonzero_state_fields="8"),
            _valid_row(executed_aerosol_pairs="24"),
            _valid_row(executed_transport_fields="50"),
            _valid_row(factorization_seconds_per_column="0"),
        )
        for row in bad_rows:
            with self.subTest(row=row):
                with self.assertRaises(ScreeningInputError):
                    ingest_rows([row])

    def test_old_width_columns_cannot_substitute_for_v2_occupancy(self) -> None:
        row = _valid_row(active_pairs="4", transported_fields="9")
        del row["nonzero_aerosol_pairs"]
        del row["nonzero_state_fields"]

        with self.assertRaises(ScreeningInputError):
            ingest_rows([row])

    def test_paired_phase_diagnostics_preserve_raw_observations(self) -> None:
        row = _valid_row()
        artifacts = ingest_paired_artifacts([row], _phase_rows(row))

        self.assertEqual(len(artifacts.phase_diagnostics), 1)
        diagnostics = artifacts.phase_diagnostics[0]
        self.assertIsInstance(diagnostics, PhaseDiagnostics)
        self.assertNotIsInstance(diagnostics, TimingEstimate)
        self.assertEqual(diagnostics.instrumented_seconds_per_column, (2, 3, 1))
        self.assertEqual(
            diagnostics.combined_factorization_solve_seconds_per_column,
            (0.6, 0.7, 0.5),
        )
        self.assertIsNone(diagnostics.factorization_seconds_per_column)
        self.assertIsNone(diagnostics.solve_seconds_per_column)

    def test_sidecar_coverage_correspondence_and_medians_fail_closed(self) -> None:
        row = _valid_row()
        phases = _phase_rows(row)
        malformed_pairs = (
            phases[:-1],
            phases + [dict(phases[0])],
            [dict(item) for item in phases],
            [dict(item) for item in phases],
        )
        malformed_pairs[2][1]["instrumented_seconds_per_column"] = "9"
        malformed_pairs[3][1]["combined_factorization_solve_seconds_per_column"] = (
            "0.1"
        )
        for phase_rows in malformed_pairs:
            with self.subTest(phase_rows=phase_rows):
                with self.assertRaises(ScreeningInputError):
                    ingest_paired_artifacts([row], phase_rows)

    def test_invalid_rows_have_no_samples_and_timed_exclusions_keep_them(self) -> None:
        invalid = _valid_row(
            case_id="8",
            validation_ok="false",
            eligible="false",
            richardson_ok="false",
            screening_error="",
            median_seconds_per_column="",
            ci95_lower="",
            ci95_upper="",
            repetitions="0",
            raw_samples_seconds_per_column="",
            instrumented_samples_seconds_per_column="",
            instrumented_median_seconds_per_column="",
            setup_seconds_per_column="",
            combined_factorization_solve_seconds_per_column="",
            update_seconds_per_column="",
            diagnostics_seconds_per_column="",
            fence_seconds_per_column="",
        )
        artifacts = ingest_paired_artifacts([invalid], [])
        self.assertEqual(artifacts.phase_diagnostics, ())
        with self.assertRaises(ScreeningInputError):
            ingest_paired_artifacts([invalid], _phase_rows(_valid_row(case_id="8")))

        excluded = _valid_row(case_id="9", eligible="false")
        artifacts = ingest_paired_artifacts([excluded], _phase_rows(excluded))
        self.assertEqual(len(artifacts.phase_diagnostics), 1)
        self.assertEqual(len(artifacts.dataset.exclusions), 1)

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

    def test_invalid_row_still_requires_valid_repetition_policy(self) -> None:
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
            bootstrap_resamples="not-an-integer",
        )

        with self.assertRaises(ScreeningInputError):
            ingest_rows([row])

    def test_each_analysis_gate_has_a_specific_exclusion(self) -> None:
        # Catches admitting unsafe points or hiding their specific gate failure.
        rows = [
            _valid_row(
                case_id="1", dt="33.33", eligible="false", clipping_count="1"
            ),
            _valid_row(case_id="2", dt="50", material_negative_count="1"),
            _valid_row(case_id="3", dt="75", eligible="false"),
            _valid_row(
                case_id="4",
                dt="100",
                richardson_ok="false",
                screening_error="",
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
        # Catches trusting a stale statistic instead of reconstructing the bootstrap.
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
            (2.0, 2.0),
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
        self.assertEqual(forward[0][0].nonzero_state_fields, 9)
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


class MatchedErrorComparisonTests(unittest.TestCase):
    def _compare(
        self,
        dataset,
        first_axis,
        second_axis,
        target_error=0.005,
    ):
        comparison = getattr(screening, "compare_matched_error_cost", None)
        if comparison is None:
            self.fail("matched-error comparison status path is missing")
        return comparison(
            dataset,
            dataset.points[0].group,
            first_axis,
            second_axis,
            target_error,
        )

    def test_status_contract_is_closed_and_constructor_is_fail_closed(self) -> None:
        status_type = getattr(screening, "MatchedErrorStatus", None)
        if status_type is None:
            self.fail("matched-error statuses are not a closed enum")
        self.assertEqual(
            {status.value for status in status_type},
            {
                "determinate",
                "indeterminate",
                "unavailable",
                "unsupported_cost_basis",
            },
        )
        with self.assertRaises(ScreeningInputError):
            screening.MatchedErrorComparison("unknown", None)
        with self.assertRaises(ScreeningInputError):
            screening.MatchedErrorComparison(
                status_type.UNSUPPORTED_COST_BASIS, ("imex", None)
            )
        with self.assertRaises(ScreeningInputError):
            screening.MatchedErrorComparison(status_type.DETERMINATE, None)

    def test_same_timestep_nonoverlapping_primary_timing_is_determinate(self) -> None:
        dataset = ingest_rows(
            [
                _timing_row(1, 0.005, 2.0, 1.8, 2.2, method="imex"),
                _timing_row(
                    2,
                    0.005,
                    3.0,
                    2.7,
                    3.3,
                    method="backward_euler",
                ),
            ]
        )

        result = self._compare(
            dataset, ("imex", None), ("backward_euler", None)
        )

        self.assertEqual(result.status.value, "determinate")
        self.assertEqual(result.winner_axis, ("imex", None))

    def test_tie_and_overlapping_uncertainty_are_indeterminate(self) -> None:
        fixtures = (
            (
                _timing_row(1, 0.005, 2.0, 2.0, 2.0, method="imex"),
                _timing_row(
                    2,
                    0.005,
                    2.0,
                    2.0,
                    2.0,
                    method="backward_euler",
                ),
            ),
            (
                _timing_row(3, 0.005, 2.0, 1.8, 2.2, method="imex"),
                _timing_row(
                    4,
                    0.005,
                    2.3,
                    2.1,
                    2.5,
                    method="backward_euler",
                ),
            ),
        )
        for rows in fixtures:
            with self.subTest(rows=rows):
                result = self._compare(
                    ingest_rows(rows),
                    ("imex", None),
                    ("backward_euler", None),
                )
                self.assertEqual(result.status.value, "indeterminate")
                self.assertIsNone(result.winner_axis)

    def test_unbracketed_or_missing_common_support_is_unavailable(self) -> None:
        dataset = ingest_rows(
            [
                _timing_row(1, 0.001, 2.0, 2.0, 2.0, method="imex"),
                _timing_row(
                    2,
                    0.003,
                    3.0,
                    3.0,
                    3.0,
                    method="backward_euler",
                ),
            ]
        )

        result = self._compare(
            dataset,
            ("imex", None),
            ("backward_euler", None),
            target_error=0.005,
        )

        self.assertEqual(result.status.value, "unavailable")
        self.assertIsNone(result.winner_axis)

    def test_bracketed_or_direct_cross_timestep_cost_is_unsupported(self) -> None:
        bracketed = ingest_rows(
            [
                _timing_row(
                    1, 0.001, 4.0, 4.0, 4.0, method="imex", dt=25.0
                ),
                _timing_row(
                    2, 0.01, 1.0, 1.0, 1.0, method="imex", dt=100.0
                ),
                _timing_row(
                    3,
                    0.001,
                    6.0,
                    6.0,
                    6.0,
                    method="backward_euler",
                    dt=25.0,
                ),
                _timing_row(
                    4,
                    0.01,
                    2.0,
                    2.0,
                    2.0,
                    method="backward_euler",
                    dt=100.0,
                ),
            ]
        )
        direct = ingest_rows(
            [
                _timing_row(
                    5, 0.005, 2.0, 2.0, 2.0, method="imex", dt=50.0
                ),
                _timing_row(
                    6,
                    0.005,
                    3.0,
                    3.0,
                    3.0,
                    method="backward_euler",
                    dt=100.0,
                ),
            ]
        )

        for dataset in (bracketed, direct):
            with self.subTest(dataset=dataset):
                result = self._compare(
                    dataset,
                    ("imex", None),
                    ("backward_euler", None),
                )
                self.assertEqual(
                    result.status.value, "unsupported_cost_basis"
                )
                self.assertIsNone(result.winner_axis)
                self.assertEqual(
                    set(vars(result)), {"status", "winner_axis"}
                )

    def test_phase_diagnostics_cannot_enter_comparative_timing(self) -> None:
        row = _valid_row()
        diagnostics = ingest_paired_artifacts(
            [row], _phase_rows(row)
        ).phase_diagnostics[0]
        dataset = ingest_rows(
            [
                _timing_row(1, 0.005, 2.0, 2.0, 2.0, method="imex"),
                _timing_row(
                    2,
                    0.005,
                    3.0,
                    3.0,
                    3.0,
                    method="backward_euler",
                ),
            ]
        )
        comparison = getattr(screening, "compare_matched_error_cost", None)
        if comparison is None:
            self.fail("matched-error comparison status path is missing")

        with self.assertRaises(ScreeningInputError):
            comparison(
                diagnostics,
                dataset.points[0].group,
                ("imex", None),
                ("backward_euler", None),
                0.005,
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
                    ci95_relative_half_width="0",
                    raw_samples_seconds_per_column="2;2;2",
                    stop_reason="ci_target",
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
