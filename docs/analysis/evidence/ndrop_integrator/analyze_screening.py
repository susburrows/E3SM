"""Pure analysis primitives for validated N-drop screening measurements.

This module contains no file-system access.  Callers supply already decoded
records, which keeps unit tests synthetic and lets provenance-aware drivers own
all input and output paths.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
import math
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


class ScreeningInputError(ValueError):
    """Raised when a screening record violates the analysis input contract."""


@dataclass(frozen=True)
class TimingEstimate:
    """Median and supplied 95% CI, all in seconds per advanced column."""

    median_seconds_per_column: float
    ci95_lower_seconds_per_column: float
    ci95_upper_seconds_per_column: float

    @property
    def ci95_relative_half_width(self) -> float:
        """Return CI half-width divided by the median."""

        half_width = (
            self.ci95_upper_seconds_per_column
            - self.ci95_lower_seconds_per_column
        ) / 2.0
        if self.median_seconds_per_column > 0.0:
            return half_width / self.median_seconds_per_column
        return 0.0 if half_width == 0.0 else math.inf


@dataclass(frozen=True, order=True)
class GroupKey:
    """Deterministic workload key for matched error-cost comparisons."""

    nlev: int
    nonzero_state_fields: int
    batch_size: int
    mixing_index: int
    mixing_scale: float


@dataclass(frozen=True)
class ScreeningPoint:
    """One eligible measured point; ``screening_error`` is dimensionless."""

    case_id: int
    group: GroupKey
    method: str
    dt_seconds: float
    fe_nsubmix: Optional[int]
    screening_error: float
    timing: TimingEstimate


class ExclusionReason(str, Enum):
    """Reasons a valid record is excluded from matched-error analysis."""

    VALIDATION_FAILURE = "validation_failure"
    CLIPPING = "clipping"
    MATERIAL_NEGATIVITY = "material_negativity"
    INELIGIBLE = "ineligible"
    MISSING_RICHARDSON_ERROR = "missing_richardson_error"


@dataclass(frozen=True)
class Exclusion:
    """A valid screening record omitted from matched-error comparisons."""

    case_id: int
    group: GroupKey
    method: str
    reason: ExclusionReason


@dataclass(frozen=True)
class AnalysisDataset:
    """Immutable included points and explicitly accounted-for exclusions."""

    points: Tuple[ScreeningPoint, ...]
    exclusions: Tuple[Exclusion, ...]


@dataclass(frozen=True)
class PhaseDiagnostics:
    """Raw diagnostic-only instrumented observations for one benchmark case.

    This type deliberately has no :class:`TimingEstimate` conversion and is
    not accepted by comparative functions.
    """

    run_id: str
    case_id: int
    method: str
    linear_phase_mode: str
    instrumented_seconds_per_column: Tuple[float, ...]
    setup_seconds_per_column: Tuple[float, ...]
    factorization_seconds_per_column: Optional[Tuple[float, ...]]
    solve_seconds_per_column: Optional[Tuple[float, ...]]
    combined_factorization_solve_seconds_per_column: Optional[
        Tuple[float, ...]
    ]
    update_seconds_per_column: Tuple[float, ...]
    diagnostics_seconds_per_column: Tuple[float, ...]
    fence_seconds_per_column: Tuple[float, ...]


@dataclass(frozen=True)
class PairedArtifacts:
    """Validated summary analysis plus separate diagnostic phase evidence."""

    dataset: AnalysisDataset
    phase_diagnostics: Tuple[PhaseDiagnostics, ...]


class SelectionStatus(str, Enum):
    """Availability of a measured runtime-versus-error target."""

    MEASURED = "measured"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True)
class TargetSelection:
    """Fastest measured point at or below a dimensionless error target."""

    group: GroupKey
    method: str
    fe_nsubmix: Optional[int]
    target_error: float
    status: SelectionStatus
    point: Optional[ScreeningPoint]


class BreakEvenStatus(str, Enum):
    """Availability of an FE-submix break-even estimate."""

    AVAILABLE = "available"
    UNAVAILABLE = "unavailable"
    NO_CROSSING = "no_crossing"


@dataclass(frozen=True)
class BreakEvenEstimate:
    """Adjacent measured FE-submix break-even and conservative 95% bounds."""

    status: BreakEvenStatus
    group: Optional[GroupKey]
    candidate_method: Optional[str]
    median_fe_nsubmix: Optional[float]
    ci95_lower_fe_nsubmix: Optional[float]
    ci95_upper_fe_nsubmix: Optional[float]
    first_measured_fe_nsubmix: Optional[int]
    adjacent_bracket: Optional[Tuple[int, int]]


_REQUIRED_ROW_FIELDS = frozenset(
    {
        "artifact_schema",
        "run_id",
        "case_id",
        "nlev",
        "nonzero_aerosol_pairs",
        "nonzero_state_fields",
        "executed_aerosol_pairs",
        "executed_transport_fields",
        "batch_size",
        "mixing_index",
        "mixing_scale",
        "dt",
        "method",
        "linear_phase_mode",
        "forced_fe_substeps",
        "validation_ok",
        "eligible",
        "richardson_ok",
        "screening_error",
        "median_seconds_per_column",
        "ci95_lower",
        "ci95_upper",
        "ci95_relative_half_width",
        "clipping_count",
        "material_negative_count",
        "repetitions",
        "minimum_repetitions",
        "maximum_repetitions",
        "bootstrap_resamples",
        "raw_samples_seconds_per_column",
        "instrumented_samples_seconds_per_column",
        "instrumented_median_seconds_per_column",
        "instrumented_ci95_lower",
        "instrumented_ci95_upper",
        "instrumentation_overhead",
        "stop_reason",
        "setup_seconds_per_column",
        "factorization_seconds_per_column",
        "solve_seconds_per_column",
        "combined_factorization_solve_seconds_per_column",
        "update_seconds_per_column",
        "diagnostics_seconds_per_column",
        "fence_seconds_per_column",
    }
)

_REQUIRED_PHASE_FIELDS = frozenset(
    {
        "artifact_schema",
        "run_id",
        "case_id",
        "method",
        "linear_phase_mode",
        "instrumented_repetition",
        "instrumented_seconds_per_column",
        "setup_seconds_per_column",
        "factorization_seconds_per_column",
        "solve_seconds_per_column",
        "combined_factorization_solve_seconds_per_column",
        "update_seconds_per_column",
        "diagnostics_seconds_per_column",
        "fence_seconds_per_column",
    }
)

_RUN_ID = re.compile(r"^[A-Za-z0-9._-]{1,128}$")
_LINEAR_PHASE_MODE = {
    "forward_euler": "none",
    "imex": "combined",
    "backward_euler": "separate",
}


def _finite_nonnegative(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ScreeningInputError("{} must be a finite nonnegative number".format(name))
    try:
        converted = float(value)
    except (TypeError, ValueError):
        raise ScreeningInputError(
            "{} must be a finite nonnegative number".format(name)
        )
    if not math.isfinite(converted) or converted < 0.0:
        raise ScreeningInputError("{} must be a finite nonnegative number".format(name))
    return converted


def _integer(value: object, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ScreeningInputError("{} must be an integer".format(name))
    try:
        converted = int(value)
    except (TypeError, ValueError):
        raise ScreeningInputError("{} must be an integer".format(name))
    if str(converted) != str(value).strip() or converted < minimum:
        raise ScreeningInputError(
            "{} must be an integer >= {}".format(name, minimum)
        )
    return converted


def _boolean(value: object, name: str) -> bool:
    if isinstance(value, bool):
        return value
    if value == "true":
        return True
    if value == "false":
        return False
    raise ScreeningInputError("{} must be true or false".format(name))


def _blank(value: object) -> bool:
    return value is None or value == ""


def _required_or_blank_timing(
    value: object, name: str, required: bool
) -> Optional[float]:
    if required:
        if _blank(value):
            raise ScreeningInputError("{} must be numeric".format(name))
        return _finite_nonnegative(value, name)
    if not _blank(value):
        raise ScreeningInputError("{} must be blank".format(name))
    return None


def _validated_run_id(value: object) -> str:
    if not isinstance(value, str) or _RUN_ID.fullmatch(value) is None:
        raise ScreeningInputError(
            "run_id must match [A-Za-z0-9._-]{1,128}"
        )
    return value


def _validate_linear_contract(
    row: Mapping, method: str, has_samples: bool
) -> None:
    expected_mode = _LINEAR_PHASE_MODE[method]
    if row["linear_phase_mode"] != expected_mode:
        raise ScreeningInputError(
            "linear_phase_mode does not match the benchmark method"
        )
    _required_or_blank_timing(
        row["factorization_seconds_per_column"],
        "factorization_seconds_per_column",
        has_samples and method == "backward_euler",
    )
    _required_or_blank_timing(
        row["solve_seconds_per_column"],
        "solve_seconds_per_column",
        has_samples and method == "backward_euler",
    )
    _required_or_blank_timing(
        row["combined_factorization_solve_seconds_per_column"],
        "combined_factorization_solve_seconds_per_column",
        has_samples and method == "imex",
    )
    for name in (
        "setup_seconds_per_column",
        "update_seconds_per_column",
        "diagnostics_seconds_per_column",
        "fence_seconds_per_column",
    ):
        _required_or_blank_timing(row[name], name, has_samples)


def _row_identity(row: Mapping) -> Tuple[int, GroupKey, str]:
    case_id = _integer(row["case_id"], "case_id")
    nonzero_pairs = _integer(
        row["nonzero_aerosol_pairs"], "nonzero_aerosol_pairs", 1
    )
    if nonzero_pairs not in {1, 4, 8, 16, 25}:
        raise ScreeningInputError(
            "nonzero_aerosol_pairs is outside the approved occupancy set"
        )
    nonzero_fields = _integer(
        row["nonzero_state_fields"], "nonzero_state_fields", 1
    )
    if nonzero_fields != 1 + 2 * nonzero_pairs:
        raise ScreeningInputError(
            "nonzero_state_fields must equal 1 + 2*nonzero_aerosol_pairs"
        )
    if _integer(row["executed_aerosol_pairs"], "executed_aerosol_pairs") != 25:
        raise ScreeningInputError("executed_aerosol_pairs must equal 25")
    if _integer(row["executed_transport_fields"], "executed_transport_fields") != 51:
        raise ScreeningInputError("executed_transport_fields must equal 51")
    group = GroupKey(
        _integer(row["nlev"], "nlev", 1),
        nonzero_fields,
        _integer(row["batch_size"], "batch_size", 1),
        _integer(row["mixing_index"], "mixing_index"),
        _finite_nonnegative(row["mixing_scale"], "mixing_scale"),
    )
    method = row["method"]
    if method not in {"forward_euler", "imex", "backward_euler"}:
        raise ScreeningInputError("method is not a supported screening method")
    return case_id, group, method


def _raw_timing_samples(value: object) -> Tuple[float, ...]:
    if isinstance(value, str):
        if not value:
            return ()
        raw_values = value.split(";")
        if any(item == "" for item in raw_values):
            raise ScreeningInputError("raw timing samples contain an empty value")
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        raw_values = value
    else:
        raise ScreeningInputError(
            "raw timing samples must be semicolon text or a sequence"
        )
    return tuple(
        _finite_nonnegative(item, "raw timing sample") for item in raw_values
    )


def _cpp_sample_median(samples: Sequence[float]) -> float:
    """Match ``sample_median`` in ``ndrop_benchmark.cpp`` exactly."""

    ordered = sorted(samples)
    middle = len(ordered) // 2
    upper = ordered[middle]
    if len(ordered) % 2:
        return upper
    lower = ordered[middle - 1]
    return lower + (upper - lower) / 2.0


def _finite_number(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ScreeningInputError("{} must be finite".format(name))
    try:
        converted = float(value)
    except (TypeError, ValueError):
        raise ScreeningInputError("{} must be finite".format(name))
    if not math.isfinite(converted):
        raise ScreeningInputError("{} must be finite".format(name))
    return converted


def _reconstruct_primary_timing(
    row: Mapping,
    case_id: int,
    raw_samples: Tuple[float, ...],
    repetitions: int,
    minimum: int,
    maximum: int,
    resamples: int,
) -> TimingEstimate:
    """Replay the producer's primary bootstrap and stopping policy exactly."""

    if repetitions < minimum or repetitions > maximum:
        raise ScreeningInputError(
            "valid repetition count must be within the configured policy"
        )
    stop_reason = row["stop_reason"]
    if stop_reason not in {"ci_target", "maximum_repetitions"}:
        raise ScreeningInputError("valid row has an invalid stop_reason")

    seed = 0x4E44524F ^ case_id
    intervals = tuple(
        bootstrap_median_interval(raw_samples[:count], resamples, seed)
        for count in range(minimum, repetitions + 1)
    )
    first_target = next(
        (
            minimum + offset
            for offset, interval in enumerate(intervals)
            if interval.ci95_relative_half_width < 0.05
        ),
        None,
    )
    if first_target is not None:
        if stop_reason != "ci_target" or first_target != repetitions:
            raise ScreeningInputError(
                "timing samples continue after the first eligible CI target"
            )
    elif stop_reason != "maximum_repetitions" or repetitions != maximum:
        raise ScreeningInputError(
            "maximum_repetitions must reach the configured maximum without "
            "meeting the CI target"
        )

    reconstructed = intervals[-1]
    reported = TimingEstimate(
        _finite_nonnegative(
            row["median_seconds_per_column"], "median_seconds_per_column"
        ),
        _finite_nonnegative(row["ci95_lower"], "ci95_lower"),
        _finite_nonnegative(row["ci95_upper"], "ci95_upper"),
    )
    reported_relative_half_width = _finite_nonnegative(
        row["ci95_relative_half_width"], "ci95_relative_half_width"
    )
    if (
        reported != reconstructed
        or reported_relative_half_width
        != reconstructed.ci95_relative_half_width
    ):
        raise ScreeningInputError(
            "reported primary timing must equal its deterministic bootstrap"
        )
    return reconstructed


def ingest_rows(rows: Iterable[Mapping]) -> AnalysisDataset:
    """Validate decoded benchmark rows and account for excluded records.

    Extra source columns are allowed, but all core schema fields must exist.
    Exclusion is an analysis decision and is returned as data; malformed schema,
    types, shapes, or values raise :class:`ScreeningInputError`.
    """

    if isinstance(rows, (str, bytes, Mapping)):
        raise ScreeningInputError("rows must be an iterable of mappings")
    try:
        row_iterator = iter(rows)
    except TypeError:
        raise ScreeningInputError("rows must be an iterable of mappings")

    points = []
    exclusions = []
    seen_case_ids = set()
    seen_composite_keys = set()
    run_id = None
    for index, row in enumerate(row_iterator):
        if not isinstance(row, Mapping):
            raise ScreeningInputError("row {} must be a mapping".format(index))
        missing = sorted(_REQUIRED_ROW_FIELDS.difference(row))
        if missing:
            raise ScreeningInputError(
                "row {} is missing required fields: {}".format(
                    index, ", ".join(missing)
                )
            )
        if row["artifact_schema"] != "ndrop_summary_v2":
            raise ScreeningInputError(
                "artifact_schema must equal ndrop_summary_v2"
            )
        row_run_id = _validated_run_id(row["run_id"])
        if run_id is None:
            run_id = row_run_id
        elif row_run_id != run_id:
            raise ScreeningInputError("summary rows contain mixed run IDs")
        case_id, group, method = _row_identity(row)
        if case_id in seen_case_ids:
            raise ScreeningInputError("summary case_id values must be unique")
        seen_case_ids.add(case_id)
        composite_key = (
            group,
            _finite_nonnegative(row["dt"], "dt"),
            method,
            _integer(row["forced_fe_substeps"], "forced_fe_substeps"),
        )
        if composite_key in seen_composite_keys:
            raise ScreeningInputError("summary composite keys must be unique")
        seen_composite_keys.add(composite_key)
        validation_ok = _boolean(row["validation_ok"], "validation_ok")
        eligible = _boolean(row["eligible"], "eligible")
        richardson_ok = _boolean(row["richardson_ok"], "richardson_ok")
        clipping_count = _integer(row["clipping_count"], "clipping_count")
        material_negative_count = _integer(
            row["material_negative_count"], "material_negative_count"
        )
        raw_samples = _raw_timing_samples(row["raw_samples_seconds_per_column"])
        instrumented_samples = _raw_timing_samples(
            row["instrumented_samples_seconds_per_column"]
        )
        repetitions = _integer(row["repetitions"], "repetitions")
        minimum_repetitions = _integer(
            row["minimum_repetitions"], "minimum_repetitions", 1
        )
        maximum_repetitions = _integer(
            row["maximum_repetitions"], "maximum_repetitions", 1
        )
        bootstrap_resamples = _integer(
            row["bootstrap_resamples"], "bootstrap_resamples", 1
        )
        if maximum_repetitions < minimum_repetitions:
            raise ScreeningInputError(
                "maximum_repetitions must not be less than minimum_repetitions"
            )
        if repetitions != len(raw_samples) or repetitions != len(
            instrumented_samples
        ):
            raise ScreeningInputError(
                "repetitions must equal both raw timing sample counts"
            )
        if validation_ok != (repetitions > 0):
            raise ScreeningInputError(
                "valid rows require samples and invalid rows require zero repetitions"
            )
        _validate_linear_contract(row, method, repetitions > 0)
        dt_seconds = _finite_nonnegative(row["dt"], "dt")
        if dt_seconds == 0.0:
            raise ScreeningInputError("dt must be positive seconds")
        forced_fe_substeps = _integer(
            row["forced_fe_substeps"], "forced_fe_substeps"
        )
        if method == "forward_euler":
            if forced_fe_substeps < 1:
                raise ScreeningInputError(
                    "forward_euler requires positive forced_fe_substeps"
                )
            fe_nsubmix = forced_fe_substeps
        else:
            if forced_fe_substeps != 0:
                raise ScreeningInputError(
                    "non-forward_euler rows require zero forced_fe_substeps"
                )
            fe_nsubmix = None
        timing = None
        if repetitions > 0:
            timing = _reconstruct_primary_timing(
                row,
                case_id,
                raw_samples,
                repetitions,
                minimum_repetitions,
                maximum_repetitions,
                bootstrap_resamples,
            )
            instrumented_median = _cpp_sample_median(instrumented_samples)
            if _finite_nonnegative(
                row["instrumented_median_seconds_per_column"],
                "instrumented_median_seconds_per_column",
            ) != instrumented_median:
                raise ScreeningInputError(
                    "reported instrumented median must equal its raw sample median"
                )
            instrumented_lower = _finite_nonnegative(
                row["instrumented_ci95_lower"],
                "instrumented_ci95_lower",
            )
            instrumented_upper = _finite_nonnegative(
                row["instrumented_ci95_upper"],
                "instrumented_ci95_upper",
            )
            if not (
                instrumented_lower
                <= instrumented_median
                <= instrumented_upper
            ):
                raise ScreeningInputError(
                    "instrumented CI must contain its raw sample median"
                )
            overhead = _finite_number(
                row["instrumentation_overhead"], "instrumentation_overhead"
            )
            expected_overhead = (
                instrumented_median / timing.median_seconds_per_column - 1.0
                if timing.median_seconds_per_column > 0.0
                else 0.0
            )
            if overhead != expected_overhead:
                raise ScreeningInputError(
                    "instrumentation overhead must match the two medians"
                )
        else:
            for name in (
                "median_seconds_per_column",
                "ci95_lower",
                "ci95_upper",
                "ci95_relative_half_width",
                "instrumented_median_seconds_per_column",
                "instrumented_ci95_lower",
                "instrumented_ci95_upper",
                "instrumentation_overhead",
            ):
                if not _blank(row[name]):
                    raise ScreeningInputError(
                        "{} must be blank for invalid rows".format(name)
                    )
            if row["stop_reason"] != "invalid":
                raise ScreeningInputError(
                    "invalid row must have stop_reason=invalid"
                )
        if not validation_ok:
            exclusions.append(
                Exclusion(
                    case_id, group, method, ExclusionReason.VALIDATION_FAILURE
                )
            )
            continue
        if clipping_count:
            exclusions.append(
                Exclusion(case_id, group, method, ExclusionReason.CLIPPING)
            )
            continue
        if material_negative_count:
            exclusions.append(
                Exclusion(
                    case_id, group, method, ExclusionReason.MATERIAL_NEGATIVITY
                )
            )
            continue
        if not eligible:
            exclusions.append(
                Exclusion(case_id, group, method, ExclusionReason.INELIGIBLE)
            )
            continue
        screening_error_value = row["screening_error"]
        if (
            not richardson_ok
            or screening_error_value is None
            or screening_error_value == ""
        ):
            exclusions.append(
                Exclusion(
                    case_id,
                    group,
                    method,
                    ExclusionReason.MISSING_RICHARDSON_ERROR,
                )
            )
            continue
        screening_error = _finite_nonnegative(
            screening_error_value, "screening_error"
        )
        if timing is None:
            raise ScreeningInputError("eligible rows require primary timing")
        points.append(
            ScreeningPoint(
                case_id,
                group,
                method,
                dt_seconds,
                fe_nsubmix,
                screening_error,
                timing,
            )
        )
    return AnalysisDataset(tuple(points), tuple(exclusions))


def ingest_paired_artifacts(
    summary_rows: Iterable[Mapping], phase_rows: Iterable[Mapping]
) -> PairedArtifacts:
    """Validate a paired summary/phase artifact without mixing evidence types.

    The summary must use ``ndrop_summary_v2`` and every instrumented
    observation must be present exactly once in ``ndrop_phase_samples_v1``.
    Phase diagnostics remain raw tuples and cannot become primary timing
    estimates through this interface.
    """

    if isinstance(summary_rows, (str, bytes, Mapping)):
        raise ScreeningInputError("summary_rows must be an iterable of mappings")
    if isinstance(phase_rows, (str, bytes, Mapping)):
        raise ScreeningInputError("phase_rows must be an iterable of mappings")
    try:
        summaries = tuple(summary_rows)
        phases = tuple(phase_rows)
    except TypeError:
        raise ScreeningInputError("paired artifacts must be iterable")
    if not summaries:
        raise ScreeningInputError("paired artifacts require summary rows")
    dataset = ingest_rows(summaries)
    summary_by_case = {
        _integer(row["case_id"], "case_id"): row for row in summaries
    }
    observed: Dict[int, Dict[int, Mapping]] = {}
    for index, row in enumerate(phases):
        if not isinstance(row, Mapping):
            raise ScreeningInputError(
                "phase row {} must be a mapping".format(index)
            )
        missing = sorted(_REQUIRED_PHASE_FIELDS.difference(row))
        if missing:
            raise ScreeningInputError(
                "phase row {} is missing required fields: {}".format(
                    index, ", ".join(missing)
                )
            )
        if row["artifact_schema"] != "ndrop_phase_samples_v1":
            raise ScreeningInputError(
                "phase artifact_schema must equal ndrop_phase_samples_v1"
            )
        phase_run_id = _validated_run_id(row["run_id"])
        case_id = _integer(row["case_id"], "case_id")
        if case_id not in summary_by_case:
            raise ScreeningInputError("phase row has no matching summary row")
        summary = summary_by_case[case_id]
        if phase_run_id != summary["run_id"]:
            raise ScreeningInputError("summary and phase run IDs differ")
        method = row["method"]
        if method != summary["method"]:
            raise ScreeningInputError("summary and phase methods differ")
        if row["linear_phase_mode"] != summary["linear_phase_mode"]:
            raise ScreeningInputError("summary and phase modes differ")
        repetition = _integer(
            row["instrumented_repetition"], "instrumented_repetition"
        )
        case_observed = observed.setdefault(case_id, {})
        if repetition in case_observed:
            raise ScreeningInputError(
                "phase repetition keys must be unique per case"
            )
        for name in (
            "instrumented_seconds_per_column",
            "setup_seconds_per_column",
            "update_seconds_per_column",
            "diagnostics_seconds_per_column",
            "fence_seconds_per_column",
        ):
            _finite_nonnegative(row[name], name)
        _validate_linear_contract(row, method, True)
        case_observed[repetition] = row

    diagnostics = []
    median_fields = (
        "setup_seconds_per_column",
        "factorization_seconds_per_column",
        "solve_seconds_per_column",
        "combined_factorization_solve_seconds_per_column",
        "update_seconds_per_column",
        "diagnostics_seconds_per_column",
        "fence_seconds_per_column",
    )
    for case_id in sorted(summary_by_case):
        summary = summary_by_case[case_id]
        repetitions = _integer(summary["repetitions"], "repetitions")
        case_observed = observed.get(case_id, {})
        if set(case_observed) != set(range(repetitions)):
            raise ScreeningInputError(
                "phase repetitions must be contiguous and cover the summary"
            )
        if repetitions == 0:
            continue
        ordered = tuple(case_observed[index] for index in range(repetitions))
        instrumented = tuple(
            _finite_nonnegative(
                row["instrumented_seconds_per_column"],
                "instrumented_seconds_per_column",
            )
            for row in ordered
        )
        if instrumented != _raw_timing_samples(
            summary["instrumented_samples_seconds_per_column"]
        ):
            raise ScreeningInputError(
                "sidecar instrumented totals must match summary samples elementwise"
            )
        sample_columns: Dict[str, Tuple[float, ...]] = {}
        for name in median_fields:
            if _blank(summary[name]):
                if any(not _blank(row[name]) for row in ordered):
                    raise ScreeningInputError(
                        "sidecar linear nullability differs from the summary"
                    )
                continue
            values = tuple(
                _finite_nonnegative(row[name], name) for row in ordered
            )
            if _finite_nonnegative(summary[name], name) != _cpp_sample_median(
                values
            ):
                raise ScreeningInputError(
                    "summary phase median cannot be reconstructed exactly"
                )
            sample_columns[name] = values
        diagnostics.append(
            PhaseDiagnostics(
                str(summary["run_id"]),
                case_id,
                str(summary["method"]),
                str(summary["linear_phase_mode"]),
                instrumented,
                sample_columns["setup_seconds_per_column"],
                sample_columns.get("factorization_seconds_per_column"),
                sample_columns.get("solve_seconds_per_column"),
                sample_columns.get(
                    "combined_factorization_solve_seconds_per_column"
                ),
                sample_columns["update_seconds_per_column"],
                sample_columns["diagnostics_seconds_per_column"],
                sample_columns["fence_seconds_per_column"],
            )
        )
    return PairedArtifacts(dataset, tuple(diagnostics))


def _point_order(point: ScreeningPoint) -> tuple:
    return (
        point.method,
        -1 if point.fe_nsubmix is None else point.fe_nsubmix,
        point.dt_seconds,
        point.screening_error,
        point.timing.median_seconds_per_column,
        point.case_id,
    )


def deterministic_groups(
    points: Iterable[ScreeningPoint],
) -> Tuple[Tuple[GroupKey, Tuple[ScreeningPoint, ...]], ...]:
    """Group points by workload in deterministic key and point order.

    The key fixes levels, nonzero-state occupancy, batch size, and both the
    mixing-regime label and dimensionless strength. Executed width remains
    fixed at 51 fields. Time step, method, and FE submix count remain
    comparison axes.
    """

    grouped: Dict[GroupKey, List[ScreeningPoint]] = {}
    for point in points:
        if not isinstance(point, ScreeningPoint):
            raise ScreeningInputError("grouped values must be ScreeningPoint instances")
        grouped.setdefault(point.group, []).append(point)
    return tuple(
        (group, tuple(sorted(grouped[group], key=_point_order)))
        for group in sorted(grouped)
    )


def select_measured_targets(
    points: Iterable[ScreeningPoint], targets: Iterable[float]
) -> Tuple[TargetSelection, ...]:
    """Select fastest qualifying measured points without error extrapolation.

    A point qualifies when its dimensionless Richardson screening error is no
    greater than the requested target.  Target values and output axes are
    sorted so identical inputs produce identical report order.
    """

    checked_targets = tuple(
        sorted(_finite_nonnegative(value, "target error") for value in targets)
    )
    if any(value == 0.0 for value in checked_targets):
        raise ScreeningInputError("target errors must be positive")
    by_axis: Dict[Tuple[GroupKey, str, Optional[int]], List[ScreeningPoint]] = {}
    for point in points:
        if not isinstance(point, ScreeningPoint):
            raise ScreeningInputError("target values must be ScreeningPoint instances")
        by_axis.setdefault(
            (point.group, point.method, point.fe_nsubmix), []
        ).append(point)

    selections = []
    ordered_axes = sorted(
        by_axis.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            -1 if item[0][2] is None else item[0][2],
        ),
    )
    for (group, method, fe_nsubmix), axis_points in ordered_axes:
        for target in checked_targets:
            qualifying = [
                point for point in axis_points if point.screening_error <= target
            ]
            if not qualifying:
                selections.append(
                    TargetSelection(
                        group,
                        method,
                        fe_nsubmix,
                        target,
                        SelectionStatus.UNAVAILABLE,
                        None,
                    )
                )
                continue
            selected = min(
                qualifying,
                key=lambda point: (
                    point.timing.median_seconds_per_column,
                    point.screening_error,
                    point.dt_seconds,
                    -1 if point.fe_nsubmix is None else point.fe_nsubmix,
                    point.case_id,
                ),
            )
            selections.append(
                TargetSelection(
                    group,
                    method,
                    fe_nsubmix,
                    target,
                    SelectionStatus.MEASURED,
                    selected,
                )
            )
    return tuple(selections)


def _first_adjacent_crossing(
    fe_points: Sequence[ScreeningPoint],
    fe_timing_attribute: str,
    candidate_seconds_per_column: float,
) -> Optional[Tuple[float, int, Tuple[int, int]]]:
    previous = None
    for point in fe_points:
        current_time = getattr(point.timing, fe_timing_attribute)
        if current_time >= candidate_seconds_per_column:
            current_nsubmix = point.fe_nsubmix
            if current_nsubmix is None:
                raise ScreeningInputError("FE point is missing its submix count")
            if current_time == candidate_seconds_per_column:
                return (
                    float(current_nsubmix),
                    current_nsubmix,
                    (current_nsubmix, current_nsubmix),
                )
            if previous is None:
                return (
                    float(current_nsubmix),
                    current_nsubmix,
                    (current_nsubmix, current_nsubmix),
                )
            previous_time = getattr(previous.timing, fe_timing_attribute)
            previous_nsubmix = previous.fe_nsubmix
            if previous_nsubmix is None:
                raise ScreeningInputError("FE point is missing its submix count")
            fraction = (
                (candidate_seconds_per_column - previous_time)
                / (current_time - previous_time)
            )
            interpolated = previous_nsubmix + fraction * (
                current_nsubmix - previous_nsubmix
            )
            return (
                interpolated,
                current_nsubmix,
                (previous_nsubmix, current_nsubmix),
            )
        previous = point
    return None


def interpolate_break_even(
    fe_points: Iterable[ScreeningPoint],
    candidate: Optional[ScreeningPoint],
) -> BreakEvenEstimate:
    """Interpolate FE-submix break-even only between adjacent measurements.

    The median is the first FE timing crossing of the candidate timing.  The
    conservative 95% lower bound crosses FE upper timing against candidate
    lower timing; the upper bound crosses FE lower timing against candidate
    upper timing.  Missing inputs and median no-crossing cases are explicit,
    and no value is extrapolated outside measured FE submix counts.
    """

    checked_fe = tuple(fe_points)
    if candidate is not None and not isinstance(candidate, ScreeningPoint):
        raise ScreeningInputError("candidate must be a ScreeningPoint or None")
    if candidate is not None and candidate.method not in {
        "imex",
        "backward_euler",
    }:
        raise ScreeningInputError("candidate must be imex or backward_euler")
    if candidate is None or not checked_fe:
        return BreakEvenEstimate(
            BreakEvenStatus.UNAVAILABLE,
            None if candidate is None else candidate.group,
            None if candidate is None else candidate.method,
            None,
            None,
            None,
            None,
            None,
        )
    for point in checked_fe:
        if not isinstance(point, ScreeningPoint) or point.method != "forward_euler":
            raise ScreeningInputError("break-even inputs must be forward_euler points")
        if point.group != candidate.group:
            raise ScreeningInputError("break-even points must share one workload group")
        if point.dt_seconds != candidate.dt_seconds:
            raise ScreeningInputError(
                "break-even points must share one ACI interval"
            )
        if point.fe_nsubmix is None:
            raise ScreeningInputError("FE point is missing its submix count")
    ordered = tuple(sorted(checked_fe, key=lambda point: point.fe_nsubmix))
    nsubmix_values = tuple(point.fe_nsubmix for point in ordered)
    if len(set(nsubmix_values)) != len(nsubmix_values):
        raise ScreeningInputError("FE submix measurements must be unique")

    median_crossing = _first_adjacent_crossing(
        ordered,
        "median_seconds_per_column",
        candidate.timing.median_seconds_per_column,
    )
    if median_crossing is None:
        return BreakEvenEstimate(
            BreakEvenStatus.NO_CROSSING,
            candidate.group,
            candidate.method,
            None,
            None,
            None,
            None,
            None,
        )
    lower_crossing = _first_adjacent_crossing(
        ordered,
        "ci95_upper_seconds_per_column",
        candidate.timing.ci95_lower_seconds_per_column,
    )
    upper_crossing = _first_adjacent_crossing(
        ordered,
        "ci95_lower_seconds_per_column",
        candidate.timing.ci95_upper_seconds_per_column,
    )
    return BreakEvenEstimate(
        BreakEvenStatus.AVAILABLE,
        candidate.group,
        candidate.method,
        median_crossing[0],
        None if lower_crossing is None else lower_crossing[0],
        None if upper_crossing is None else upper_crossing[0],
        median_crossing[1],
        median_crossing[2],
    )


def summarize_timing(
    samples_seconds_per_column: Sequence[float],
    ci95_seconds_per_column: Sequence[float],
) -> TimingEstimate:
    """Return the exact sample median with a caller-supplied 95% timing CI.

    Parameters are in seconds per advanced column.  The sample must be nonempty
    and both the samples and interval must be finite and nonnegative.  The CI
    must contain the exact sample median.
    """

    if isinstance(samples_seconds_per_column, (str, bytes)):
        raise ScreeningInputError("timing samples must be a nonempty sequence")
    try:
        samples = tuple(
            _finite_nonnegative(value, "timing sample")
            for value in samples_seconds_per_column
        )
    except TypeError:
        raise ScreeningInputError("timing samples must be a nonempty sequence")
    if not samples:
        raise ScreeningInputError("timing samples must be a nonempty sequence")
    if isinstance(ci95_seconds_per_column, (str, bytes)):
        raise ScreeningInputError("timing CI must contain exactly two bounds")
    try:
        ci95 = tuple(ci95_seconds_per_column)
    except TypeError:
        raise ScreeningInputError("timing CI must contain exactly two bounds")
    if len(ci95) != 2:
        raise ScreeningInputError("timing CI must contain exactly two bounds")
    lower = _finite_nonnegative(ci95[0], "timing CI lower")
    upper = _finite_nonnegative(ci95[1], "timing CI upper")
    exact_median = _cpp_sample_median(samples)
    if lower > exact_median or exact_median > upper:
        raise ScreeningInputError("timing CI must contain the sample median")
    return TimingEstimate(exact_median, lower, upper)


def bootstrap_median_interval(
    samples_seconds_per_column: Sequence[float],
    resamples: int,
    seed: int = 0x4E44524F,
) -> TimingEstimate:
    """Reproduce the benchmark's deterministic percentile-bootstrap median CI.

    The generator, sampling with replacement, and percentile indices match
    ``ndrop_benchmark.cpp`` exactly.  This function is intended for retained
    result verification and small synthetic fixtures; callers remain
    responsible for checking the recorded case-specific seed.
    """

    checked_resamples = _integer(resamples, "bootstrap resamples", 1)
    checked_seed = _integer(seed, "bootstrap seed")
    if checked_seed > 0xFFFFFFFFFFFFFFFF:
        raise ScreeningInputError("bootstrap seed must fit in 64 bits")
    if isinstance(samples_seconds_per_column, (str, bytes)):
        raise ScreeningInputError("bootstrap samples must be a nonempty sequence")
    try:
        samples = tuple(
            _finite_nonnegative(value, "bootstrap sample")
            for value in samples_seconds_per_column
        )
    except TypeError:
        raise ScreeningInputError("bootstrap samples must be a nonempty sequence")
    if not samples:
        raise ScreeningInputError("bootstrap samples must be a nonempty sequence")

    mask = 0xFFFFFFFFFFFFFFFF
    random_state = checked_seed
    replicate_medians = []
    for _ in range(checked_resamples):
        replicate = []
        for _ in samples:
            random_state = (random_state + 0x9E3779B97F4A7C15) & mask
            value = random_state
            value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & mask
            value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & mask
            value = (value ^ (value >> 31)) & mask
            replicate.append(samples[value % len(samples)])
        replicate_medians.append(_cpp_sample_median(replicate))

    replicate_medians.sort()
    last_index = len(replicate_medians) - 1
    lower_index = math.floor(0.025 * last_index)
    upper_index = math.ceil(0.975 * last_index)
    return TimingEstimate(
        _cpp_sample_median(samples),
        replicate_medians[lower_index],
        replicate_medians[upper_index],
    )
