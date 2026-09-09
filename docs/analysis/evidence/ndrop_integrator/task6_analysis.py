"""Fail-closed Task 6 analysis entry points and retained output receipts."""

from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import re
import sys
from typing import Dict, Optional, Tuple

from analyze_screening import (
    AnalysisDataset,
    GroupKey,
    MatchedErrorStatus,
    ScreeningInputError,
    ScreeningPoint,
    compare_matched_error_cost,
)
from ndrop_artifacts import (
    ArtifactError,
    expected_retained_case,
    load_validated_analysis_input,
    sha256_file,
)


RETAINED_RECEIPT_SCHEMA = "ndrop_task6_analysis_receipt_v1"
DIRECT_ANALYSIS_MODE = "direct_same_timestep"
COMPUTED_NSUBMIX_EXACT_MODE = "computed_nsubmix_exact_same_timestep"
PRODUCTION_NSUBMIX_SOURCE_COMMIT = (
    "9af7f38a8f441b346d07ddd1f36eede1926c956f"
)
_LEVELS = (1, 2, 3, 16, 32, 64, 72, 127, 128)
_MIXING_SCALES = (0.05, 0.5)
_INTERVALS = (33.33, 50.0, 75.0, 100.0, 150.0, 300.0)
_MEASURED_FE_NSUBMIX = (1, 2, 3, 5, 8, 16, 21)
_COMPARATORS = ("imex", "backward_euler")
_MIXING_RATE_THRESHOLD = 1.0e-6
_SOURCE_COMMIT = re.compile(r"^[0-9a-f]{40}$")
_RECEIPT_FIELDS = frozenset(
    {
        "artifact_schema",
        "completion_status",
        "identity_manifest_sha256",
        "final_manifest_sha256",
        "analysis_source_commit",
        "analysis_parameters",
        "analysis_command",
        "analysis_output",
        "analysis_output_sha256",
    }
)


@dataclass(frozen=True)
class RetainedAnalysisOutput:
    """A completed retained result and its authoritative receipt."""

    analysis_path: Path
    receipt_path: Path


def _axis_payload(method: str, fe_nsubmix: Optional[int]) -> dict:
    result = {"method": method}
    if fe_nsubmix is not None:
        result["fe_nsubmix"] = fe_nsubmix
    return result


def _group_payload(group: GroupKey) -> dict:
    return {
        "nlev": group.nlev,
        "nonzero_state_fields": group.nonzero_state_fields,
        "batch_size": group.batch_size,
        "mixing_index": group.mixing_index,
        "mixing_scale": group.mixing_scale,
    }


def _direct_comparison_result(
    first: Optional[ScreeningPoint], second: Optional[ScreeningPoint]
) -> dict:
    """Compare primary timings only when both points share one measured step."""

    if first is None or second is None:
        return {"status": MatchedErrorStatus.UNAVAILABLE.value}
    if not isinstance(first, ScreeningPoint) or not isinstance(second, ScreeningPoint):
        raise ScreeningInputError(
            "direct comparison requires validated ScreeningPoint values"
        )
    if first.group != second.group or first.dt_seconds != second.dt_seconds:
        return {"status": MatchedErrorStatus.UNSUPPORTED_COST_BASIS.value}

    first_timing = first.timing
    second_timing = second.timing
    values = (
        first_timing.median_seconds_per_column,
        first_timing.ci95_lower_seconds_per_column,
        first_timing.ci95_upper_seconds_per_column,
        second_timing.median_seconds_per_column,
        second_timing.ci95_lower_seconds_per_column,
        second_timing.ci95_upper_seconds_per_column,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ScreeningInputError(
            "direct runtime ratios require finite positive primary timing bounds"
        )

    overlaps = not (
        first_timing.ci95_upper_seconds_per_column
        < second_timing.ci95_lower_seconds_per_column
        or second_timing.ci95_upper_seconds_per_column
        < first_timing.ci95_lower_seconds_per_column
    )
    result = {
        "status": (
            MatchedErrorStatus.INDETERMINATE.value
            if overlaps
            else MatchedErrorStatus.DETERMINATE.value
        ),
        "runtime_ratio_first_over_second": (
            first_timing.median_seconds_per_column
            / second_timing.median_seconds_per_column
        ),
        "runtime_ratio_ci95_lower": (
            first_timing.ci95_lower_seconds_per_column
            / second_timing.ci95_upper_seconds_per_column
        ),
        "runtime_ratio_ci95_upper": (
            first_timing.ci95_upper_seconds_per_column
            / second_timing.ci95_lower_seconds_per_column
        ),
    }
    if not overlaps:
        winner = first if (
            first_timing.median_seconds_per_column
            < second_timing.median_seconds_per_column
        ) else second
        result["winner"] = _axis_payload(winner.method, winner.fe_nsubmix)
    return result


def _safe_median(values: Tuple[float, ...]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return ordered[middle - 1] + (
        ordered[middle] - ordered[middle - 1]
    ) / 2.0


def build_direct_same_timestep_payload(dataset: AnalysisDataset) -> dict:
    """Build the measured IMEX/BE/FE runtime-comparison matrix."""

    if not isinstance(dataset, AnalysisDataset):
        raise ScreeningInputError(
            "direct comparison matrix requires a validated AnalysisDataset"
        )
    by_measurement = {}
    group_steps = set()
    fe_counts: Dict[GroupKey, set] = {}
    for point in dataset.points:
        if not isinstance(point, ScreeningPoint):
            raise ScreeningInputError(
                "direct comparison matrix requires ScreeningPoint values"
            )
        axis = (point.method, point.fe_nsubmix)
        key = (point.group, point.dt_seconds, axis)
        if key in by_measurement:
            raise ScreeningInputError("direct comparison measurement is duplicated")
        by_measurement[key] = point
        group_steps.add((point.group, point.dt_seconds))
        if point.method == "forward_euler":
            fe_counts.setdefault(point.group, set()).add(point.fe_nsubmix)

    status_counts = {status.value: 0 for status in MatchedErrorStatus}
    comparisons = []
    summary_values = {}
    summary_statuses = {}
    for group, dt_seconds in sorted(group_steps):
        pairings = [
            (
                "imex_vs_backward_euler",
                ("imex", None),
                ("backward_euler", None),
            )
        ]
        for fe_nsubmix in sorted(fe_counts.get(group, set())):
            pairings.extend(
                (
                    (
                        "forward_euler_nsubmix_{}_vs_imex".format(fe_nsubmix),
                        ("forward_euler", fe_nsubmix),
                        ("imex", None),
                    ),
                    (
                        "forward_euler_nsubmix_{}_vs_backward_euler".format(
                            fe_nsubmix
                        ),
                        ("forward_euler", fe_nsubmix),
                        ("backward_euler", None),
                    ),
                )
            )
        for pairing, first_axis, second_axis in pairings:
            result = _direct_comparison_result(
                by_measurement.get((group, dt_seconds, first_axis)),
                by_measurement.get((group, dt_seconds, second_axis)),
            )
            status_counts[result["status"]] += 1
            summary_statuses.setdefault(
                pairing, {status.value: 0 for status in MatchedErrorStatus}
            )[result["status"]] += 1
            if "runtime_ratio_first_over_second" in result:
                summary_values.setdefault(pairing, []).append(
                    result["runtime_ratio_first_over_second"]
                )
            comparisons.append(
                {
                    "pairing": pairing,
                    "group": _group_payload(group),
                    "dt_seconds": dt_seconds,
                    "first_axis": _axis_payload(*first_axis),
                    "second_axis": _axis_payload(*second_axis),
                    "result": result,
                }
            )

    pairing_summaries = []
    for pairing in sorted(summary_statuses):
        summary = {
            "pairing": pairing,
            "comparison_count": sum(summary_statuses[pairing].values()),
            "status_counts": summary_statuses[pairing],
        }
        ratios = tuple(summary_values.get(pairing, ()))
        if ratios:
            summary.update(
                {
                    "measured_ratio_count": len(ratios),
                    "runtime_ratio_first_over_second_minimum": min(ratios),
                    "runtime_ratio_first_over_second_median": _safe_median(ratios),
                    "runtime_ratio_first_over_second_maximum": max(ratios),
                }
            )
        pairing_summaries.append(summary)

    return {
        "analysis_mode": DIRECT_ANALYSIS_MODE,
        "timing_source": "primary_uninstrumented",
        "input_point_count": len(dataset.points),
        "input_exclusion_count": len(dataset.exclusions),
        "status_counts": status_counts,
        "pairing_summaries": pairing_summaries,
        "comparisons": comparisons,
    }


def _direct_analysis_parameters() -> dict:
    return {
        "analysis_mode": DIRECT_ANALYSIS_MODE,
        "comparison_policy": (
            "imex_vs_backward_euler_and_each_fe_nsubmix_vs_imex_and_be"
        ),
        "cost_basis": "same_timestep_only",
        "timing_source": "primary_uninstrumented",
        "uncertainty": "conservative_ratio_of_reconstructed_ci95_bounds",
    }


def _require_binary64() -> None:
    if not (
        sys.float_info.radix == 2
        and sys.float_info.mant_dig == 53
        and sys.float_info.max_exp == 1024
    ):
        raise ArtifactError("production nsubmix reconstruction requires binary64")


def reconstruct_production_nsubmix(
    nlev: int, mixing_scale: float, dt_seconds: float
) -> int:
    """Reproduce retained make_operator/prepare_operator/production_nsubmix."""

    _require_binary64()
    if type(nlev) is not int or nlev <= 0:
        raise ArtifactError("nlev must be a positive integer")
    if (
        isinstance(mixing_scale, bool)
        or not isinstance(mixing_scale, (int, float))
        or not math.isfinite(mixing_scale)
        or mixing_scale < 0.0
    ):
        raise ArtifactError("mixing scale must be finite and nonnegative")
    if (
        isinstance(dt_seconds, bool)
        or not isinstance(dt_seconds, (int, float))
        or not math.isfinite(dt_seconds)
        or dt_seconds < 0.0
    ):
        raise ArtifactError("timestep must be finite and nonnegative")

    # Retained benchmark make_operator(): construct the synthetic column.
    pressure_weight = [
        0.8 + 0.4 * float((17 * level + 3) % 11) / 10.0
        for level in range(nlev)
    ]
    interface_exchange = [0.0] * (nlev + 1)
    for interface in range(1, nlev):
        offset = (7 * interface + 1) % 5 - 2
        modulation = 1.0 + 0.05 * offset
        interface_exchange[interface] = (
            mixing_scale
            * modulation
            * min(pressure_weight[interface - 1], pressure_weight[interface])
        )

    # Retained prepare_operator()/production_nsubmix(): derive rates in source
    # order, apply the stability threshold, and use positive C++ truncation.
    first_level_controller_rate = 0.25 * mixing_scale
    minimum_timescale = sys.float_info.max
    for level in range(nlev):
        rate_above = interface_exchange[level] / pressure_weight[level]
        rate_below = interface_exchange[level + 1] / pressure_weight[level]
        controller_rate = first_level_controller_rate if level == 0 else 0.0
        inverse_time = rate_below + rate_above + controller_rate
        if inverse_time > _MIXING_RATE_THRESHOLD:
            minimum_timescale = min(minimum_timescale, 1.0 / inverse_time)

    candidate_step = 0.9 * minimum_timescale
    substep_ratio = dt_seconds / candidate_step
    if (
        not math.isfinite(substep_ratio)
        or substep_ratio < 0.0
        or substep_ratio > 2147483646.0
    ):
        return -1
    return int(substep_ratio) + 1


def _mixing_condition(mixing_index: int, mixing_scale: float) -> str:
    expected = _MIXING_SCALES[mixing_index] if mixing_index in (0, 1) else None
    if mixing_scale != expected:
        raise ArtifactError("mixing index and scale are inconsistent")
    return (
        "Lower-intensity mixing (scale 0.05)"
        if mixing_index == 0
        else "Higher-intensity mixing (scale 0.5)"
    )


def _support_class(nsubmix: int) -> str:
    if nsubmix in _MEASURED_FE_NSUBMIX:
        return "exact_measured"
    if _MEASURED_FE_NSUBMIX[0] <= nsubmix <= _MEASURED_FE_NSUBMIX[-1]:
        return "not_measured_in_range"
    return "outside_measured_range"


def _direct_group_key(group: dict) -> tuple:
    fields = {
        "nlev",
        "nonzero_state_fields",
        "batch_size",
        "mixing_index",
        "mixing_scale",
    }
    if not isinstance(group, dict) or set(group) != fields:
        raise ArtifactError("direct comparison group fields are invalid")
    if (
        type(group["nlev"]) is not int
        or group["nlev"] not in _LEVELS
        or type(group["nonzero_state_fields"]) is not int
        or group["nonzero_state_fields"] not in (3, 9, 17, 33, 51)
        or type(group["batch_size"]) is not int
        or group["batch_size"] not in (1, 16, 64, 256, 1024)
        or type(group["mixing_index"]) is not int
        or isinstance(group["mixing_scale"], bool)
        or not isinstance(group["mixing_scale"], (int, float))
    ):
        raise ArtifactError("direct comparison group values are invalid")
    _mixing_condition(group["mixing_index"], group["mixing_scale"])
    return (
        group["nlev"],
        group["nonzero_state_fields"],
        group["batch_size"],
        group["mixing_index"],
        group["mixing_scale"],
    )


def _validated_direct_result(
    result: dict, first_axis: dict, second_axis: dict
) -> dict:
    if not isinstance(result, dict) or "status" not in result:
        raise ArtifactError("direct comparison result is malformed")
    status = result["status"]
    if status == MatchedErrorStatus.UNAVAILABLE.value:
        if set(result) != {"status"}:
            raise ArtifactError("unavailable direct result must be status-only")
        return {"status": status}
    if status == MatchedErrorStatus.UNSUPPORTED_COST_BASIS.value:
        raise ArtifactError("same-timestep direct result has unsupported cost basis")
    ratio_fields = {
        "status",
        "runtime_ratio_first_over_second",
        "runtime_ratio_ci95_lower",
        "runtime_ratio_ci95_upper",
    }
    expected_fields = (
        ratio_fields | {"winner"}
        if status == MatchedErrorStatus.DETERMINATE.value
        else ratio_fields
    )
    if status not in (
        MatchedErrorStatus.DETERMINATE.value,
        MatchedErrorStatus.INDETERMINATE.value,
    ) or set(result) != expected_fields:
        raise ArtifactError("direct comparison status or fields are invalid")
    values = tuple(result[name] for name in sorted(ratio_fields - {"status"}))
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0.0
        for value in values
    ):
        raise ArtifactError("direct comparison ratios must be finite and positive")
    ratio = result["runtime_ratio_first_over_second"]
    lower = result["runtime_ratio_ci95_lower"]
    upper = result["runtime_ratio_ci95_upper"]
    if not lower <= ratio <= upper:
        raise ArtifactError("direct comparison ratio is outside its interval")
    if status == MatchedErrorStatus.INDETERMINATE.value:
        if not lower <= 1.0 <= upper:
            raise ArtifactError("indeterminate direct result must overlap parity")
        return dict(result)
    if upper < 1.0:
        expected_winner = first_axis
    elif lower > 1.0:
        expected_winner = second_axis
    else:
        raise ArtifactError("determinate direct result overlaps parity")
    if result["winner"] != expected_winner:
        raise ArtifactError("determinate direct result has an inconsistent winner")
    normalized = dict(result)
    normalized["winner"] = dict(result["winner"])
    return normalized


def _direct_comparison_index(direct_payload: dict) -> dict:
    if not isinstance(direct_payload, dict):
        raise ArtifactError("direct analysis payload must be an object")
    if direct_payload.get("analysis_mode") != DIRECT_ANALYSIS_MODE:
        raise ArtifactError("parent analysis mode is not direct same-timestep")
    if direct_payload.get("timing_source") != "primary_uninstrumented":
        raise ArtifactError("parent analysis is not primary uninstrumented timing")
    if "phase_diagnostics" in direct_payload:
        raise ArtifactError("phase diagnostics cannot enter comparative analysis")
    records = direct_payload.get("comparisons")
    if not isinstance(records, list):
        raise ArtifactError("parent direct comparisons are missing")

    index = {}
    for record in records:
        if not isinstance(record, dict) or set(record) != {
            "pairing",
            "group",
            "dt_seconds",
            "first_axis",
            "second_axis",
            "result",
        }:
            raise ArtifactError("direct comparison record fields are invalid")
        if record["pairing"] == "imex_vs_backward_euler":
            continue
        first_axis = record["first_axis"]
        second_axis = record["second_axis"]
        if (
            not isinstance(first_axis, dict)
            or set(first_axis) != {"method", "fe_nsubmix"}
            or first_axis["method"] != "forward_euler"
            or type(first_axis["fe_nsubmix"]) is not int
            or first_axis["fe_nsubmix"] not in _MEASURED_FE_NSUBMIX
            or not isinstance(second_axis, dict)
            or set(second_axis) != {"method"}
            or second_axis["method"] not in _COMPARATORS
        ):
            raise ArtifactError("direct comparison axes are invalid")
        nsubmix = first_axis["fe_nsubmix"]
        comparator = second_axis["method"]
        expected_pairing = "forward_euler_nsubmix_{}_vs_{}".format(
            nsubmix, comparator
        )
        if record["pairing"] != expected_pairing:
            raise ArtifactError("direct comparison pairing is inconsistent")
        group_key = _direct_group_key(record["group"])
        dt_seconds = record["dt_seconds"]
        if (
            isinstance(dt_seconds, bool)
            or not isinstance(dt_seconds, (int, float))
            or dt_seconds not in _INTERVALS
        ):
            raise ArtifactError("direct comparison timestep is invalid")
        key = group_key + (dt_seconds, nsubmix, comparator)
        if key in index:
            raise ArtifactError("direct comparison record is duplicated")
        index[key] = _validated_direct_result(
            record["result"], first_axis, second_axis
        )
    return index


def _computed_nsubmix_map() -> list:
    entries = []
    for mixing_index, mixing_scale in enumerate(_MIXING_SCALES):
        for nlev in _LEVELS:
            for dt_seconds in _INTERVALS:
                nsubmix = reconstruct_production_nsubmix(
                    nlev, mixing_scale, dt_seconds
                )
                entries.append(
                    {
                        "mixing_index": mixing_index,
                        "mixing_scale": mixing_scale,
                        "mixing_condition": _mixing_condition(
                            mixing_index, mixing_scale
                        ),
                        "nlev": nlev,
                        "dt_seconds": dt_seconds,
                        "computed_fe_nsubmix": nsubmix,
                        "fe_timing_support": _support_class(nsubmix),
                    }
                )
    return entries


def build_computed_nsubmix_exact_payload(direct_payload: dict) -> dict:
    """Extract exact computed-count FE comparisons from validated direct output."""

    direct_index = _direct_comparison_index(direct_payload)
    computed_map = _computed_nsubmix_map()
    support_counts = {
        "exact_measured": 0,
        "not_measured_in_range": 0,
        "outside_measured_range": 0,
    }
    status_counts = {status.value: 0 for status in MatchedErrorStatus}
    comparisons = []
    strata = {}
    exact_workloads = 0
    valid_direct = 0

    for case_id in range(0, 24300, 9):
        case = expected_retained_case(case_id)
        nsubmix = reconstruct_production_nsubmix(
            case["nlev"], case["mixing_scale"], case["dt"]
        )
        support = _support_class(nsubmix)
        if support == "exact_measured":
            exact_workloads += 1
        group = {
            "nlev": case["nlev"],
            "nonzero_state_fields": case["nonzero_state_fields"],
            "batch_size": case["batch_size"],
            "mixing_index": case["mixing_index"],
            "mixing_scale": case["mixing_scale"],
        }
        group_key = _direct_group_key(group)
        for comparator in _COMPARATORS:
            support_counts[support] += 1
            result = {"status": MatchedErrorStatus.UNAVAILABLE.value}
            if support == "exact_measured":
                result = direct_index.get(
                    group_key + (case["dt"], nsubmix, comparator), result
                )
            if "runtime_ratio_first_over_second" in result:
                valid_direct += 1
            status_counts[result["status"]] += 1
            comparison = {
                "group": group,
                "dt_seconds": case["dt"],
                "computed_fe_nsubmix": nsubmix,
                "fe_timing_support": support,
                "comparator": {"method": comparator},
                "result": result,
            }
            comparisons.append(comparison)
            stratum_key = (
                case["mixing_index"],
                case["mixing_scale"],
                case["nlev"],
                case["dt"],
                comparator,
            )
            stratum = strata.setdefault(
                stratum_key,
                {
                    "mixing_index": case["mixing_index"],
                    "mixing_scale": case["mixing_scale"],
                    "mixing_condition": _mixing_condition(
                        case["mixing_index"], case["mixing_scale"]
                    ),
                    "nlev": case["nlev"],
                    "dt_seconds": case["dt"],
                    "computed_fe_nsubmix": nsubmix,
                    "fe_timing_support": support,
                    "comparator": comparator,
                    "designed_workload_count": 0,
                    "status_counts": {
                        status.value: 0 for status in MatchedErrorStatus
                    },
                    "ratios": [],
                },
            )
            stratum["designed_workload_count"] += 1
            stratum["status_counts"][result["status"]] += 1
            if "runtime_ratio_first_over_second" in result:
                stratum["ratios"].append(
                    result["runtime_ratio_first_over_second"]
                )

    stratum_payloads = []
    for key in sorted(strata):
        stratum = strata[key]
        ratios = tuple(stratum.pop("ratios"))
        stratum["valid_direct_comparison_count"] = len(ratios)
        if ratios:
            stratum.update(
                {
                    "runtime_ratio_first_over_second_minimum": min(ratios),
                    "runtime_ratio_first_over_second_median": _safe_median(ratios),
                    "runtime_ratio_first_over_second_maximum": max(ratios),
                }
            )
        stratum_payloads.append(stratum)

    return {
        "analysis_mode": COMPUTED_NSUBMIX_EXACT_MODE,
        "timing_source": "primary_uninstrumented",
        "cost_basis": "same_timestep_only",
        "fe_timing_basis": (
            "measured_forced_substeps_equal_reconstructed_production_nsubmix"
        ),
        "adaptive_selection_path_timed": False,
        "vertical_axis": "nlev_proxy_not_physical_layer_thickness",
        "mixing_conditions": [
            {
                "mixing_index": index,
                "mixing_scale": scale,
                "label": _mixing_condition(index, scale),
            }
            for index, scale in enumerate(_MIXING_SCALES)
        ],
        "measured_fe_nsubmix": list(_MEASURED_FE_NSUBMIX),
        "production_nsubmix_source_commit": PRODUCTION_NSUBMIX_SOURCE_COMMIT,
        "design_counts": {
            "workload_timestep_cases": 2700,
            "pairwise_comparison_opportunities": len(comparisons),
            "exact_match_workload_opportunities": exact_workloads,
            "exact_match_comparison_opportunities": (
                support_counts["exact_measured"]
            ),
            "valid_direct_comparisons": valid_direct,
        },
        "support_counts": support_counts,
        "status_counts": status_counts,
        "computed_nsubmix_map": computed_map,
        "strata": stratum_payloads,
        "comparisons": comparisons,
    }


def _computed_exact_analysis_parameters(
    parent_analysis_hash: str, parent_receipt_hash: str
) -> dict:
    return {
        "analysis_mode": COMPUTED_NSUBMIX_EXACT_MODE,
        "cost_basis": "same_timestep_only",
        "timing_source": "primary_uninstrumented",
        "match_policy": "exact_integer_nsubmix_only_no_interpolation",
        "fe_timing_basis": (
            "measured_forced_substeps_equal_reconstructed_production_nsubmix"
        ),
        "adaptive_selection_path_timed": False,
        "vertical_axis": "nlev_proxy_not_physical_layer_thickness",
        "measured_fe_nsubmix": list(_MEASURED_FE_NSUBMIX),
        "production_nsubmix_source_commit": PRODUCTION_NSUBMIX_SOURCE_COMMIT,
        "parent_analysis_output_sha256": parent_analysis_hash,
        "parent_analysis_receipt_sha256": parent_receipt_hash,
    }


def _comparison_payload(
    dataset: AnalysisDataset,
    group: GroupKey,
    first_axis: Tuple[str, Optional[int]],
    second_axis: Tuple[str, Optional[int]],
    target_error: float,
) -> dict:
    comparison = compare_matched_error_cost(
        dataset,
        group,
        first_axis,
        second_axis,
        target_error,
    )
    result = {"status": comparison.status.value}
    if comparison.status is MatchedErrorStatus.DETERMINATE:
        if comparison.winner_axis is None:
            raise RuntimeError("determinate comparison is missing its winner")
        method, fe_nsubmix = comparison.winner_axis
        winner = {"method": method}
        if fe_nsubmix is not None:
            winner["fe_nsubmix"] = fe_nsubmix
        result["winner"] = winner
    return {"comparison": result}


def _analysis_parameters(
    group: GroupKey,
    first_axis: Tuple[str, Optional[int]],
    second_axis: Tuple[str, Optional[int]],
    target_error: float,
) -> dict:
    return {
        "group": {
            "nlev": group.nlev,
            "nonzero_state_fields": group.nonzero_state_fields,
            "batch_size": group.batch_size,
            "mixing_index": group.mixing_index,
            "mixing_scale": group.mixing_scale,
        },
        "first_axis": {
            "method": first_axis[0],
            "fe_nsubmix": first_axis[1],
        },
        "second_axis": {
            "method": second_axis[0],
            "fe_nsubmix": second_axis[1],
        },
        "target_error": target_error,
    }


def _checked_analysis_provenance(
    analysis_source_commit: str,
    analysis_command: Tuple[str, ...],
) -> Tuple[str, Tuple[str, ...]]:
    if (
        not isinstance(analysis_source_commit, str)
        or _SOURCE_COMMIT.fullmatch(analysis_source_commit) is None
    ):
        raise ArtifactError(
            "analysis source commit must be a lowercase 40-character SHA"
        )
    if (
        not isinstance(analysis_command, tuple)
        or not analysis_command
        or any(not isinstance(item, str) or not item for item in analysis_command)
    ):
        raise ArtifactError("analysis command must be a nonempty tuple of strings")
    return analysis_source_commit, analysis_command


def _output_paths(output_directory: Path) -> Tuple[Path, Path, Path, Path]:
    root = Path(output_directory)
    return (
        root / "analysis.json.partial",
        root / "analysis.json",
        root / "analysis.done.json.partial",
        root / "analysis.done.json",
    )


def _require_paths_absent(paths: Tuple[Path, ...]) -> None:
    for path in paths:
        if os.path.lexists(str(path)):
            raise ArtifactError("refusing to overwrite {}".format(path))


def _fsync_directory(directory: Path) -> None:
    try:
        descriptor = os.open(str(directory), os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except OSError as error:
        raise ArtifactError("cannot fsync output directory: {}".format(error))


def _write_json_exclusive(path: Path, document: dict) -> None:
    payload = (json.dumps(document, sort_keys=True, indent=2) + "\n").encode()
    try:
        with Path(path).open("xb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(path, 0o444)
        _fsync_directory(Path(path).parent)
    except FileExistsError:
        raise ArtifactError("refusing to overwrite {}".format(path))
    except OSError as error:
        raise ArtifactError("cannot write {}: {}".format(path, error))


def _rename_without_overwrite(source: Path, destination: Path) -> None:
    try:
        os.link(source, destination)
        os.unlink(source)
        _fsync_directory(destination.parent)
    except FileExistsError:
        raise ArtifactError("refusing to overwrite {}".format(destination))
    except OSError as error:
        raise ArtifactError(
            "cannot promote {} to {}: {}".format(source, destination, error)
        )


def run_task6_analysis(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    output_path: Path,
    *,
    group: GroupKey,
    first_axis: Tuple[str, Optional[int]],
    second_axis: Tuple[str, Optional[int]],
    target_error: float,
    require_retained_contract: bool = True,
) -> Path:
    """Validate paired input, compare primary timing, and write one JSON result."""

    if require_retained_contract:
        raise ArtifactError(
            "retained-contract execution requires run_retained_task6_analysis"
        )
    paired = load_validated_analysis_input(
        identity_manifest_path,
        final_manifest_path,
        require_retained_contract=False,
    )
    payload = _comparison_payload(
        paired.dataset,
        group,
        first_axis,
        second_axis,
        target_error,
    )
    output = Path(output_path)
    output.write_bytes(
        (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    )
    return output


def run_retained_task6_analysis(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    output_directory: Path,
    *,
    group: GroupKey,
    first_axis: Tuple[str, Optional[int]],
    second_axis: Tuple[str, Optional[int]],
    target_error: float,
    analysis_source_commit: str,
    analysis_command: Tuple[str, ...],
) -> RetainedAnalysisOutput:
    """Run retained analysis through the strict validator and seal one receipt."""

    source_commit, command = _checked_analysis_provenance(
        analysis_source_commit, analysis_command
    )
    paths = _output_paths(output_directory)
    _require_paths_absent(paths)
    identity_path = Path(identity_manifest_path)
    final_path = Path(final_manifest_path)
    identity_hash = sha256_file(identity_path)
    final_hash = sha256_file(final_path)
    paired = load_validated_analysis_input(
        identity_path,
        final_path,
        require_retained_contract=True,
    )
    if (
        sha256_file(identity_path) != identity_hash
        or sha256_file(final_path) != final_hash
    ):
        raise ArtifactError("validated input manifests changed during analysis loading")

    payload = _comparison_payload(
        paired.dataset,
        group,
        first_axis,
        second_axis,
        target_error,
    )
    if (
        sha256_file(identity_path) != identity_hash
        or sha256_file(final_path) != final_hash
    ):
        raise ArtifactError("validated input manifests changed during analysis")

    analysis_partial, analysis_path, receipt_partial, receipt_path = paths
    analysis_partial.parent.mkdir(parents=True, exist_ok=True)
    _require_paths_absent(paths)
    _write_json_exclusive(analysis_partial, payload)
    _rename_without_overwrite(analysis_partial, analysis_path)
    receipt = {
        "artifact_schema": RETAINED_RECEIPT_SCHEMA,
        "completion_status": "complete",
        "identity_manifest_sha256": identity_hash,
        "final_manifest_sha256": final_hash,
        "analysis_source_commit": source_commit,
        "analysis_parameters": _analysis_parameters(
            group, first_axis, second_axis, target_error
        ),
        "analysis_command": list(command),
        "analysis_output": "analysis.json",
        "analysis_output_sha256": sha256_file(analysis_path),
    }
    _write_json_exclusive(receipt_partial, receipt)
    _rename_without_overwrite(receipt_partial, receipt_path)
    return RetainedAnalysisOutput(analysis_path, receipt_path)


def validate_retained_analysis_output(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    output_directory: Path,
    *,
    group: GroupKey,
    first_axis: Tuple[str, Optional[int]],
    second_axis: Tuple[str, Optional[int]],
    target_error: float,
    analysis_source_commit: str,
    analysis_command: Tuple[str, ...],
) -> RetainedAnalysisOutput:
    """Validate one retained result receipt without reopening benchmark CSVs."""

    source_commit, command = _checked_analysis_provenance(
        analysis_source_commit, analysis_command
    )
    analysis_partial, analysis_path, receipt_partial, receipt_path = _output_paths(
        output_directory
    )
    if analysis_partial.exists() or receipt_partial.exists():
        raise ArtifactError("retained analysis has incomplete partial output")
    if not analysis_path.is_file() or not receipt_path.is_file():
        raise ArtifactError("retained analysis output or receipt is missing")
    try:
        with receipt_path.open() as stream:
            receipt = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ArtifactError("cannot read retained analysis receipt: {}".format(error))
    if not isinstance(receipt, dict) or frozenset(receipt) != _RECEIPT_FIELDS:
        raise ArtifactError("retained analysis receipt fields are invalid")
    expected = {
        "artifact_schema": RETAINED_RECEIPT_SCHEMA,
        "completion_status": "complete",
        "identity_manifest_sha256": sha256_file(Path(identity_manifest_path)),
        "final_manifest_sha256": sha256_file(Path(final_manifest_path)),
        "analysis_source_commit": source_commit,
        "analysis_parameters": _analysis_parameters(
            group, first_axis, second_axis, target_error
        ),
        "analysis_command": list(command),
        "analysis_output": "analysis.json",
        "analysis_output_sha256": sha256_file(analysis_path),
    }
    if receipt != expected:
        raise ArtifactError("retained analysis receipt binding is invalid")
    return RetainedAnalysisOutput(analysis_path, receipt_path)


def run_retained_direct_task6_analysis(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    output_directory: Path,
    *,
    analysis_source_commit: str,
    analysis_command: Tuple[str, ...],
) -> RetainedAnalysisOutput:
    """Run the direct measured matrix through the strict retained boundary."""

    source_commit, command = _checked_analysis_provenance(
        analysis_source_commit, analysis_command
    )
    paths = _output_paths(output_directory)
    _require_paths_absent(paths)
    identity_path = Path(identity_manifest_path)
    final_path = Path(final_manifest_path)
    identity_hash = sha256_file(identity_path)
    final_hash = sha256_file(final_path)
    paired = load_validated_analysis_input(
        identity_path,
        final_path,
        require_retained_contract=True,
    )
    if (
        sha256_file(identity_path) != identity_hash
        or sha256_file(final_path) != final_hash
    ):
        raise ArtifactError("validated input manifests changed during analysis loading")

    payload = build_direct_same_timestep_payload(paired.dataset)
    if (
        sha256_file(identity_path) != identity_hash
        or sha256_file(final_path) != final_hash
    ):
        raise ArtifactError("validated input manifests changed during analysis")

    analysis_partial, analysis_path, receipt_partial, receipt_path = paths
    analysis_partial.parent.mkdir(parents=True, exist_ok=True)
    _require_paths_absent(paths)
    _write_json_exclusive(analysis_partial, payload)
    _rename_without_overwrite(analysis_partial, analysis_path)
    receipt = {
        "artifact_schema": RETAINED_RECEIPT_SCHEMA,
        "completion_status": "complete",
        "identity_manifest_sha256": identity_hash,
        "final_manifest_sha256": final_hash,
        "analysis_source_commit": source_commit,
        "analysis_parameters": _direct_analysis_parameters(),
        "analysis_command": list(command),
        "analysis_output": "analysis.json",
        "analysis_output_sha256": sha256_file(analysis_path),
    }
    _write_json_exclusive(receipt_partial, receipt)
    _rename_without_overwrite(receipt_partial, receipt_path)
    return RetainedAnalysisOutput(analysis_path, receipt_path)


def validate_retained_direct_task6_output(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    output_directory: Path,
    *,
    analysis_source_commit: str,
    analysis_command: Tuple[str, ...],
) -> RetainedAnalysisOutput:
    """Validate one direct-matrix receipt without reopening benchmark CSVs."""

    source_commit, command = _checked_analysis_provenance(
        analysis_source_commit, analysis_command
    )
    analysis_partial, analysis_path, receipt_partial, receipt_path = _output_paths(
        output_directory
    )
    if analysis_partial.exists() or receipt_partial.exists():
        raise ArtifactError("retained analysis has incomplete partial output")
    if not analysis_path.is_file() or not receipt_path.is_file():
        raise ArtifactError("retained analysis output or receipt is missing")
    try:
        with receipt_path.open() as stream:
            receipt = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ArtifactError("cannot read retained analysis receipt: {}".format(error))
    if not isinstance(receipt, dict) or frozenset(receipt) != _RECEIPT_FIELDS:
        raise ArtifactError("retained analysis receipt fields are invalid")
    expected = {
        "artifact_schema": RETAINED_RECEIPT_SCHEMA,
        "completion_status": "complete",
        "identity_manifest_sha256": sha256_file(Path(identity_manifest_path)),
        "final_manifest_sha256": sha256_file(Path(final_manifest_path)),
        "analysis_source_commit": source_commit,
        "analysis_parameters": _direct_analysis_parameters(),
        "analysis_command": list(command),
        "analysis_output": "analysis.json",
        "analysis_output_sha256": sha256_file(analysis_path),
    }
    if receipt != expected:
        raise ArtifactError("retained analysis receipt binding is invalid")
    return RetainedAnalysisOutput(analysis_path, receipt_path)


def _load_validated_direct_payload(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    direct_output_directory: Path,
    *,
    direct_analysis_source_commit: str,
    direct_analysis_command: Tuple[str, ...],
) -> Tuple[dict, str, str]:
    direct = validate_retained_direct_task6_output(
        identity_manifest_path,
        final_manifest_path,
        direct_output_directory,
        analysis_source_commit=direct_analysis_source_commit,
        analysis_command=direct_analysis_command,
    )
    analysis_hash = sha256_file(direct.analysis_path)
    receipt_hash = sha256_file(direct.receipt_path)
    try:
        with direct.analysis_path.open() as stream:
            payload = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ArtifactError("cannot read validated direct analysis: {}".format(error))
    if (
        sha256_file(direct.analysis_path) != analysis_hash
        or sha256_file(direct.receipt_path) != receipt_hash
    ):
        raise ArtifactError("validated direct analysis changed while being read")
    return payload, analysis_hash, receipt_hash


def run_retained_computed_nsubmix_exact_analysis(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    direct_output_directory: Path,
    output_directory: Path,
    *,
    direct_analysis_source_commit: str,
    direct_analysis_command: Tuple[str, ...],
    analysis_source_commit: str,
    analysis_command: Tuple[str, ...],
) -> RetainedAnalysisOutput:
    """Create the exact computed-nsubmix extraction from sealed direct output."""

    _checked_analysis_provenance(
        direct_analysis_source_commit, direct_analysis_command
    )
    source_commit, command = _checked_analysis_provenance(
        analysis_source_commit, analysis_command
    )
    paths = _output_paths(output_directory)
    _require_paths_absent(paths)
    identity_path = Path(identity_manifest_path)
    final_path = Path(final_manifest_path)
    identity_hash = sha256_file(identity_path)
    final_hash = sha256_file(final_path)
    direct_payload, parent_analysis_hash, parent_receipt_hash = (
        _load_validated_direct_payload(
            identity_path,
            final_path,
            direct_output_directory,
            direct_analysis_source_commit=direct_analysis_source_commit,
            direct_analysis_command=direct_analysis_command,
        )
    )
    if (
        sha256_file(identity_path) != identity_hash
        or sha256_file(final_path) != final_hash
    ):
        raise ArtifactError("validated manifests changed during parent loading")

    payload = build_computed_nsubmix_exact_payload(direct_payload)
    if (
        sha256_file(identity_path) != identity_hash
        or sha256_file(final_path) != final_hash
    ):
        raise ArtifactError("validated manifests changed during derived analysis")

    analysis_partial, analysis_path, receipt_partial, receipt_path = paths
    analysis_partial.parent.mkdir(parents=True, exist_ok=True)
    _require_paths_absent(paths)
    _write_json_exclusive(analysis_partial, payload)
    _rename_without_overwrite(analysis_partial, analysis_path)
    receipt = {
        "artifact_schema": RETAINED_RECEIPT_SCHEMA,
        "completion_status": "complete",
        "identity_manifest_sha256": identity_hash,
        "final_manifest_sha256": final_hash,
        "analysis_source_commit": source_commit,
        "analysis_parameters": _computed_exact_analysis_parameters(
            parent_analysis_hash, parent_receipt_hash
        ),
        "analysis_command": list(command),
        "analysis_output": "analysis.json",
        "analysis_output_sha256": sha256_file(analysis_path),
    }
    _write_json_exclusive(receipt_partial, receipt)
    _rename_without_overwrite(receipt_partial, receipt_path)
    return RetainedAnalysisOutput(analysis_path, receipt_path)


def validate_retained_computed_nsubmix_exact_output(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    direct_output_directory: Path,
    output_directory: Path,
    *,
    direct_analysis_source_commit: str,
    direct_analysis_command: Tuple[str, ...],
    analysis_source_commit: str,
    analysis_command: Tuple[str, ...],
) -> RetainedAnalysisOutput:
    """Validate one computed-nsubmix extraction and its parent binding."""

    _checked_analysis_provenance(
        direct_analysis_source_commit, direct_analysis_command
    )
    source_commit, command = _checked_analysis_provenance(
        analysis_source_commit, analysis_command
    )
    direct = validate_retained_direct_task6_output(
        identity_manifest_path,
        final_manifest_path,
        direct_output_directory,
        analysis_source_commit=direct_analysis_source_commit,
        analysis_command=direct_analysis_command,
    )
    parent_analysis_hash = sha256_file(direct.analysis_path)
    parent_receipt_hash = sha256_file(direct.receipt_path)
    analysis_partial, analysis_path, receipt_partial, receipt_path = _output_paths(
        output_directory
    )
    if analysis_partial.exists() or receipt_partial.exists():
        raise ArtifactError("retained analysis has incomplete partial output")
    if not analysis_path.is_file() or not receipt_path.is_file():
        raise ArtifactError("retained analysis output or receipt is missing")
    try:
        with receipt_path.open() as stream:
            receipt = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ArtifactError("cannot read retained analysis receipt: {}".format(error))
    if not isinstance(receipt, dict) or frozenset(receipt) != _RECEIPT_FIELDS:
        raise ArtifactError("retained analysis receipt fields are invalid")
    expected = {
        "artifact_schema": RETAINED_RECEIPT_SCHEMA,
        "completion_status": "complete",
        "identity_manifest_sha256": sha256_file(Path(identity_manifest_path)),
        "final_manifest_sha256": sha256_file(Path(final_manifest_path)),
        "analysis_source_commit": source_commit,
        "analysis_parameters": _computed_exact_analysis_parameters(
            parent_analysis_hash, parent_receipt_hash
        ),
        "analysis_command": list(command),
        "analysis_output": "analysis.json",
        "analysis_output_sha256": sha256_file(analysis_path),
    }
    if receipt != expected:
        raise ArtifactError("retained computed-nsubmix receipt binding is invalid")
    return RetainedAnalysisOutput(analysis_path, receipt_path)
