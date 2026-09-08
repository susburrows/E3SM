"""Fail-closed Task 6 analysis entry points and retained output receipts."""

from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import re
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
    load_validated_analysis_input,
    sha256_file,
)


RETAINED_RECEIPT_SCHEMA = "ndrop_task6_analysis_receipt_v1"
DIRECT_ANALYSIS_MODE = "direct_same_timestep"
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
