"""Closed artifacts for the synthetic controller-matched ndrop study.

The accuracy target for each fixture is the error of Forward Euler at the
substep count returned by ``production_nsubmix``.  IMEX and backward Euler are
searched at every integer count, in order.  This module deliberately validates
the resulting evidence more strictly than the benchmark writes it: incomplete
searches, controller mismatches, stale hashes, or unpaired timing blocks remain
visible terminal outcomes and can never become selected timing evidence.

The study is standalone and synthetic.  A dirty or attached source checkout is
labelled exploratory; artifact renaming cannot promote it to retained evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import subprocess
from typing import Any

import ndrop_matched_artifacts as physical_artifacts


METHODS = ("forward_euler", "imex", "backward_euler")
IMPLICIT_METHODS = ("imex", "backward_euler")
FIXTURES = tuple(range(24))
CONTINUITY_FIXTURES = {
    0: (72, 0.05, 50.0, 6),
    1: (72, 0.5, 150.0, 176),
    2: (128, 0.05, 50.0, 6),
    3: (128, 0.5, 150.0, 176),
}
REQUESTED_NSUBMIX = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
STRATIFIED_INTERVAL = 150.0
# Binary64 value produced by the closed synthetic generator at mixing_scale=1.
STRATIFIED_UNIT_MAX_INVERSE_RATE = 2.1000000000000005
STRATIFIED_FIXTURES = {
    4 + 2 * index + level_index: (levels, requested)
    for index, requested in enumerate(REQUESTED_NSUBMIX)
    for level_index, levels in enumerate((72, 128))
}
MAX_STEPS = 8192
ACCURACY_BUDGET_SECONDS = 7200
BOOTSTRAP_SEED = 0x4E44524F

ACCURACY_HEADER = (
    "artifact_schema", "record_kind", "fixture", "fixture_family",
    "requested_nsubmix", "levels", "mixing_scale", "interval_seconds",
    "mixing_rate_t", "activation_rate_t", "method", "steps",
    "production_nsubmix", "reference_status", "reference_reason",
    "reference_disagreement", "fe_target_error", "status",
    "failure_reason", "error_normal", "error_tight", "error_max",
    "normal_weighted_l1_by_field", "normal_linf_by_field",
    "tight_weighted_l1_by_field", "tight_linf_by_field", "clipping_count",
    "material_negative_count", "normalized_residual",
    "normalized_ledger_error", "state_checksum",
)
TIMING_HEADER = (
    "artifact_schema", "fixture", "method", "steps", "block",
    "seconds_per_column", "state_checksum", "stop_reason",
)

ACCURACY_SCHEMA = "ndrop_controller_matched_accuracy_v1"
SELECTION_SCHEMA = "ndrop_controller_matched_selection_v1"
TIMING_SCHEMA = "ndrop_controller_matched_timing_v1"
ANALYSIS_SCHEMA = "ndrop_controller_matched_analysis_v1"
FINAL_SCHEMA = "ndrop_controller_matched_final_v1"
IDENTITY_SCHEMA = "ndrop_controller_matched_identity_v1"
COMPLETION_SCHEMA = "ndrop_controller_matched_completion_v1"

CONTROLLED_FIXTURES = {
    24: (72, 0.05, 150.0, 18), 25: (72, 0.5, 150.0, 176),
    26: (128, 0.05, 150.0, 18), 27: (128, 0.5, 150.0, 176),
}
CONTROLLED_ACTIVATION_MIXING = 0.05
CONTROLLED_ACCURACY_HEADER = (
    *ACCURACY_HEADER[:7], "activation_mixing_scale", *ACCURACY_HEADER[7:])


PHYSICAL_COLUMNS = physical_artifacts.PHYSICAL_FIXTURES
PHYSICAL_GRIDS = physical_artifacts.PHYSICAL_GRIDS
PHYSICAL_METADATA = ("turbulence", "cloud_base_wtke_m_s", "controller_rate_t",
                     "cloud_base_level", "cloud_base_spacing_m", "limiting_level",
                     "limiting_spacing_m")
PHYSICAL_ACCURACY_HEADER = (*ACCURACY_HEADER, *PHYSICAL_METADATA, "activation_limiter_count")

# Internal physical policy discriminator; True keeps the native v1 contract.
GEOMETRY_POLICY = "geometry"
GEOMETRY_INPUTS = {
    28: (72, 1), 29: (72, 1), 30: (128, 1), 31: (128, 1),
    32: (72, 2), 33: (72, 2), 34: (72, 4), 35: (72, 4),
    36: (128, 2), 37: (128, 2), 38: (128, 4), 39: (128, 4),
}
GEOMETRY_COLUMNS = {
    **PHYSICAL_COLUMNS,
    **{f: (levels, "strong" if f % 2 else "moderate",
           30. if f % 2 else 3., 1. if f % 2 else .3)
       for f, levels in ((32, 79), (33, 79), (34, 93), (35, 93),
                         (36, 141), (37, 141), (38, 167), (39, 167))},
}
GEOMETRY_METADATA = (
    "parent_grid_levels", "geometry_subdivision", "refinement_lower_m",
    "refinement_upper_m", "limiting_height_m", "limiting_above_diffusivity_m2_s",
    "limiting_below_diffusivity_m2_s", "limiting_cloud_base_faces",
)
GEOMETRY_APPROVED_ACCURACY_SHA256 = "7bf06451e74ad1b9930e8fd8df00e85ffbcbcb7772798ef5b87144b501d92cb6"
GEOMETRY_PRIOR_SELECTION_SHA256 = "c553ce935b85a3ec675f9099ca900c3cfddcf4c2a975c3f088a16cbdd36ebdf7"


def physical_cases(physical: bool | str) -> dict:
    return GEOMETRY_COLUMNS if physical == GEOMETRY_POLICY else PHYSICAL_COLUMNS


def fixture_policy(physical: bool | str) -> str:
    return ("dycoms_inspired_conservative_cloud_band_refinement_v1"
            if physical == GEOMETRY_POLICY else "dycoms_inspired_isothermal_repository_grid_v1")


def physical_accuracy_option(physical: bool | str) -> str:
    return ("--physical-geometry-controller-accuracy"
            if physical == GEOMETRY_POLICY else "--physical-controller-accuracy")


def study_schema(kind: str, controlled: bool = False, physical: bool | str = False) -> str:
    if controlled and physical:
        raise ArtifactError("physical and controlled policies are mutually exclusive")
    prefix = ("ndrop_physical_geometry_controller" if physical == GEOMETRY_POLICY else
              "ndrop_physical_controller" if physical else
              "ndrop_controlled_mixing" if controlled else "ndrop_controller_matched")
    return f"{prefix}_{kind}_v1"


def physical_identity(identity: dict[str, Any]) -> bool | str:
    physical = (GEOMETRY_POLICY if identity.get("artifact_schema") ==
                study_schema("identity", physical=GEOMETRY_POLICY) else
                identity.get("artifact_schema") == study_schema("identity", physical=True))
    if physical:
        if (identity.get("fixture_policy") != fixture_policy(physical)
                or identity.get("comparison_policy") != "natural_fe_controller_error_v1"
                or "activation_mixing_scale" in identity):
            raise ArtifactError("physical controller identity mismatch")
    elif "comparison_policy" in identity or "grid_files_sha256" in identity:
        raise ArtifactError("physical policy entered a different identity")
    return physical


def controlled_identity(identity: dict[str, Any]) -> bool:
    controlled = identity.get("artifact_schema") == study_schema("identity", True)
    if controlled:
        if (identity.get("fixture_policy") != "controlled_mixing_v1" or
                type(identity.get("activation_mixing_scale")) is not float or
                identity["activation_mixing_scale"] != CONTROLLED_ACTIVATION_MIXING):
            raise ArtifactError("controlled activation identity mismatch")
    elif not physical_identity(identity) and ("fixture_policy" in identity or "activation_mixing_scale" in identity):
        raise ArtifactError("controlled policy entered a legacy identity")
    return controlled


def accuracy_header(controlled: bool = False, physical: bool | str = False) -> tuple[str, ...]:
    if physical:
        if controlled:
            raise ArtifactError("physical and controlled policies are mutually exclusive")
        return (*PHYSICAL_ACCURACY_HEADER, *GEOMETRY_METADATA) \
            if physical == GEOMETRY_POLICY else PHYSICAL_ACCURACY_HEADER
    return CONTROLLED_ACCURACY_HEADER if controlled else ACCURACY_HEADER

SOURCES = (
    "ndrop_integrator.hpp", "ndrop_integrator.cpp", "ndrop_imex.cpp",
    "ndrop_backward_euler.cpp", "ndrop_reference.cpp",
    "ndrop_benchmark.hpp", "ndrop_benchmark.cpp",
    "benchmark_ndrop_matched_error.cpp",
    "ndrop_controller_matched_artifacts.py", "CMakeLists.txt",
)


def source_files(physical: bool | str = False) -> tuple[str, ...]:
    return (*SOURCES, "ndrop_matched_artifacts.py") if physical else SOURCES

SUCCESS_STATUSES = {"fe_baseline", "eligible_above_target",
                    "eligible_match"}
CONFIRMATION_SUCCESS = "eligible_confirmation"
FAILURE_STATUSES = {"numerical_failure", "fe_clipping",
                    "material_negative", "reference_unresolved",
                    "fe_baseline_ineligible", "fixture_controller_mismatch"}
TERMINAL_STATUSES = {"selected", "not_attained_within_cap",
                     "accuracy_budget_exhausted", "reference_unresolved",
                     "fe_baseline_ineligible", "fixture_controller_mismatch"}


class ArtifactError(ValueError):
    """An artifact violates the closed study contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_new(path: Path, content: bytes) -> str:
    temporary = path.with_suffix(path.suffix + ".partial")
    if path.exists() or temporary.exists():
        raise ArtifactError(f"output already exists: {path}")
    try:
        with temporary.open("xb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        # Preserve a partial output as explicit evidence of an interrupted run.
        raise
    return sha256(path)


def canonical_json(data: Any) -> bytes:
    return (json.dumps(data, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("utf-8")


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def csv_rows(path: Path, expected_header: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, strict=True)
        if tuple(reader.fieldnames or ()) != expected_header:
            raise ArtifactError(f"wrong closed CSV header: {path}")
        rows = list(reader)
    if not rows or any(None in row or None in row.values() for row in rows):
        raise ArtifactError(f"empty or malformed CSV: {path}")
    return rows


def finite(text: str, label: str) -> float:
    try:
        value = float(text)
    except (TypeError, ValueError) as exc:
        raise ArtifactError(f"invalid {label}") from exc
    if not math.isfinite(value):
        raise ArtifactError(f"nonfinite {label}")
    return value


def integer(text: str, label: str) -> int:
    try:
        value = int(text)
    except (TypeError, ValueError) as exc:
        raise ArtifactError(f"invalid {label}") from exc
    if str(value) != text:
        raise ArtifactError(f"noncanonical {label}")
    return value


def stratified_mixing_scale(requested: int) -> float:
    """Return the closed generator value used to land inside an N bin."""
    return 0.9 * (requested - 0.5) / (
        STRATIFIED_INTERVAL * STRATIFIED_UNIT_MAX_INVERSE_RATE)


def _blank(row: dict[str, str], names: tuple[str, ...], label: str) -> None:
    if any(row[name] for name in names):
        raise ArtifactError(f"{label} contains inapplicable values")


ERROR_COLUMNS = (
    "error_normal", "error_tight", "error_max",
    "normal_weighted_l1_by_field", "normal_linf_by_field",
    "tight_weighted_l1_by_field", "tight_linf_by_field",
)
NUMERICAL_COLUMNS = (
    "clipping_count", "material_negative_count", "normalized_residual",
    "normalized_ledger_error", "state_checksum",
)


def _validate_success(row: dict[str, str], levels: int) -> float:
    errors = [finite(row[name], name) for name in
              ("error_normal", "error_tight", "error_max")]
    if any(value < 0 for value in errors) or errors[2] != max(errors[:2]):
        raise ArtifactError("candidate error mismatch")
    maxima = []
    for name in ("normal_weighted_l1_by_field", "normal_linf_by_field",
                 "tight_weighted_l1_by_field", "tight_linf_by_field"):
        values = [finite(item, name) for item in row[name].split(";")]
        if len(values) != 51 or any(value < 0 for value in values):
            raise ArtifactError("field error vector is invalid")
        maxima.append(max(values))
    if errors != [max(maxima[:2]), max(maxima[2:]), max(maxima)]:
        raise ArtifactError("field errors do not reproduce maxima")
    residual = finite(row["normalized_residual"], "residual")
    ledger = finite(row["normalized_ledger_error"], "ledger error")
    ledger_limit = max(1e-11, 5 * levels * math.ulp(1.0))
    if (not 0 <= residual <= 1e-11 or not 0 <= ledger <= ledger_limit or
            integer(row["clipping_count"], "clipping count") != 0 or
            integer(row["material_negative_count"], "negative count") != 0):
        raise ArtifactError("ineligible numerical state admitted")
    return errors[2]


def _validate_measurement_numerics(row: dict[str, str]) -> int:
    if integer(row["clipping_count"], "clipping count") < 0 or \
            integer(row["material_negative_count"], "negative count") < 0:
        raise ArtifactError("negative diagnostic count")
    if finite(row["normalized_residual"], "residual") < 0 or \
            finite(row["normalized_ledger_error"], "ledger error") < 0:
        raise ArtifactError("negative numerical diagnostic")
    checksum = integer(row["state_checksum"], "state checksum")
    if not 0 <= checksum < 2**64:
        raise ArtifactError("state checksum is outside uint64")
    return checksum


def validate_accuracy(
    rows: list[dict[str, str]],
    fixtures: tuple[int, ...] = FIXTURES,
    controlled: bool = False,
    physical: bool | str = False,
) -> dict[tuple[int, str, int], dict[str, str]]:
    """Validate the complete controller-matched search and return measurements.

    Confirmation records are omitted from the returned search index, so they
    cannot affect selection. A separately reconstructed admission record can
    withhold timing when refinement evidence is missing or fails.
    """
    fixture_rows: dict[int, list[tuple[int, dict[str, str]]]] = {}
    metadata: dict[int, tuple[Any, ...]] = {}
    production_counts: dict[int, int] = {}
    if controlled and fixtures != tuple(CONTROLLED_FIXTURES):
        raise ArtifactError("controlled mixing requires exact four-case coverage")
    if physical and fixtures != tuple(physical_cases(physical)):
        raise ArtifactError("physical controller requires exact policy coverage")
    for position, row in enumerate(rows):
        if set(row) != set(accuracy_header(controlled, physical)):
            raise ArtifactError("wrong closed accuracy fields")
        if row["artifact_schema"] != study_schema("accuracy", controlled, physical):
            raise ArtifactError("old or mixed accuracy schema")
        fixture = integer(row["fixture"], "fixture")
        if fixture not in (physical_cases(physical) if physical else CONTROLLED_FIXTURES if controlled else FIXTURES):
            raise ArtifactError("unknown fixture")
        family = row["fixture_family"]
        expected_family = ("physical_geometry" if physical == GEOMETRY_POLICY else
                           "physical_column" if physical else "controlled_mixing" if controlled else
                           "continuity" if fixture in CONTINUITY_FIXTURES
                           else "controller_stratified")
        if family != expected_family:
            raise ArtifactError("fixture family mismatch")
        requested_text = row["requested_nsubmix"]
        requested = integer(requested_text, "requested nsubmix") \
            if requested_text else None
        levels = integer(row["levels"], "levels")
        mixing = finite(row["mixing_scale"], "mixing scale")
        interval = finite(row["interval_seconds"], "interval")
        scales = (finite(row["mixing_rate_t"], "mixing rate times interval"),
                  finite(row["activation_rate_t"],
                         "activation rate times interval"))
        production = integer(row["production_nsubmix"], "production nsubmix")
        if mixing <= 0 or interval <= 0 or any(value < 0 for value in scales) or \
                production <= 0:
            raise ArtifactError("fixture scales or controller count are invalid")
        if physical:
            expected_levels, strength, diffusivity, velocity = physical_cases(physical)[fixture]
            if (requested is not None or (levels, mixing, interval) !=
                    (expected_levels, diffusivity, 150.) or row["turbulence"] != strength or
                    finite(row["cloud_base_wtke_m_s"], "cloud-base velocity") != velocity):
                raise ArtifactError("physical fixture identity mismatch")
            rate = finite(row["controller_rate_t"], "controller stiffness")
            if rate <= 0 or production != math.floor(rate / .9) + 1:
                raise ArtifactError("physical controller count inconsistent with stiffness")
            for name in ("cloud_base_level", "limiting_level"):
                level = integer(row[name], name)
                if not 0 <= level < levels or name == "cloud_base_level" and level == levels - 1:
                    raise ArtifactError("physical level outside grid")
            for name in ("cloud_base_spacing_m", "limiting_spacing_m"):
                if finite(row[name], name) <= 0:
                    raise ArtifactError("nonpositive physical spacing")
            if row["record_kind"] == "search_terminal":
                _blank(row, ("activation_limiter_count",), "search terminal")
            elif integer(row["activation_limiter_count"], "activation limiter count") < 0:
                raise ArtifactError("negative activation limiter count")
            if physical == GEOMETRY_POLICY:
                if ((integer(row["parent_grid_levels"], "parent grid"),
                     integer(row["geometry_subdivision"], "subdivision")) != GEOMETRY_INPUTS[fixture]
                        or finite(row["refinement_lower_m"], "band lower") != 400.
                        or finite(row["refinement_upper_m"], "band upper") != 1000.):
                    raise ArtifactError("geometry policy metadata mismatch")
                if (any(finite(row[name], name) < 0 for name in
                        ("limiting_height_m", "limiting_above_diffusivity_m2_s",
                         "limiting_below_diffusivity_m2_s"))
                        or integer(row["limiting_cloud_base_faces"], "limiting cloud faces") not in (0, 1, 2)):
                    raise ArtifactError("invalid geometry limiting metadata")
        elif controlled:
            activation = finite(row["activation_mixing_scale"], "activation identity")
            if (requested is not None or activation != CONTROLLED_ACTIVATION_MIXING or
                    (levels, mixing, interval, production) != CONTROLLED_FIXTURES[fixture]):
                raise ArtifactError("controlled fixture identity mismatch")
        elif fixture in CONTINUITY_FIXTURES:
            expected_levels, expected_mixing, expected_interval, expected_n = \
                CONTINUITY_FIXTURES[fixture]
            if requested is not None or (levels, mixing, interval, production) != \
                    (expected_levels, expected_mixing, expected_interval, expected_n):
                raise ArtifactError("continuity fixture identity mismatch")
        else:
            expected_levels, expected_requested = STRATIFIED_FIXTURES[fixture]
            expected_mixing = stratified_mixing_scale(expected_requested)
            if ((levels, requested) != (expected_levels, expected_requested) or
                    interval != STRATIFIED_INTERVAL or
                    mixing != expected_mixing):
                raise ArtifactError("stratified fixture identity mismatch")
        current_metadata = (family, requested, levels, mixing, interval,
                            scales[0], scales[1], production,
                            row["reference_status"], row["reference_reason"],
                            row["reference_disagreement"])
        if physical:
            current_metadata += tuple(row[name] for name in PHYSICAL_METADATA)
            if physical == GEOMETRY_POLICY:
                current_metadata += tuple(row[name] for name in GEOMETRY_METADATA)
        if fixture in metadata and metadata[fixture] != current_metadata:
            raise ArtifactError("fixture metadata is inconsistent")
        metadata.setdefault(fixture, current_metadata)
        production_counts.setdefault(fixture, production)
        fixture_rows.setdefault(fixture, []).append((position, row))
    if set(fixture_rows) != set(fixtures):
        raise ArtifactError("accuracy fixture coverage is incomplete")
    if controlled and any(metadata[a][6] != metadata[b][6]
                          for a, b in ((24, 25), (26, 27))):
        raise ArtifactError("controlled activation stiffness differs within a grid pair")

    indexed: dict[tuple[int, str, int], dict[str, str]] = {}
    for fixture in fixtures:
        records = fixture_rows[fixture]
        reference_status = records[0][1]["reference_status"]
        reference_reason = records[0][1]["reference_reason"]
        if reference_status == "certified":
            if reference_reason:
                raise ArtifactError("certified reference has a failure reason")
            disagreement = finite(records[0][1]["reference_disagreement"],
                                  "reference disagreement")
            if not 0 <= disagreement <= 2e-11:
                raise ArtifactError("reference disagreement exceeds gate")
        elif reference_status == "unresolved":
            if not reference_reason:
                raise ArtifactError("unresolved reference lacks a failure reason")
            if records[0][1]["reference_disagreement"]:
                value = finite(records[0][1]["reference_disagreement"],
                               "reference disagreement")
                if value < 0:
                    raise ArtifactError("negative reference disagreement")
        else:
            raise ArtifactError("reference status is invalid")

        measurements: dict[str, list[tuple[int, dict[str, str]]]] = {
            method: [] for method in METHODS}
        confirmations: dict[str, list[tuple[int, dict[str, str]]]] = {
            method: [] for method in METHODS}
        terminals: dict[str, list[tuple[int, dict[str, str]]]] = {
            method: [] for method in IMPLICIT_METHODS}
        for position, row in records:
            kind = row["record_kind"]
            method = row["method"]
            if method not in METHODS:
                raise ArtifactError("unknown method")
            if kind == "measurement":
                measurements[method].append((position, row))
            elif kind == "confirmation":
                confirmations[method].append((position, row))
            elif kind == "search_terminal" and method in IMPLICIT_METHODS:
                terminals[method].append((position, row))
            else:
                raise ArtifactError("record kind is incompatible with method")

        if len(measurements["forward_euler"]) != 1:
            raise ArtifactError("fixture needs exactly one FE measurement")
        _, fe = measurements["forward_euler"][0]
        fe_steps = integer(fe["steps"], "FE steps")
        if fe_steps != production_counts[fixture]:
            raise ArtifactError("FE did not use exact production nsubmix")
        fe_status = fe["status"]
        _validate_measurement_numerics(fe)
        fe_target: float | None = None
        if fe_status == "fe_baseline":
            if reference_status != "certified" or fe["failure_reason"]:
                raise ArtifactError("eligible FE baseline lacks certified reference")
            fe_target = _validate_success(fe, integer(fe["levels"], "levels"))
            if finite(fe["fe_target_error"], "FE target error") != fe_target:
                raise ArtifactError("FE target does not equal FE maximum error")
            if fixture in STRATIFIED_FIXTURES and \
                    production_counts[fixture] != STRATIFIED_FIXTURES[fixture][1]:
                raise ArtifactError("stratified controller mismatch entered selection")
        elif fe_status in FAILURE_STATUSES:
            _blank(fe, ERROR_COLUMNS, "failed FE baseline")
            if not fe["failure_reason"]:
                raise ArtifactError("failed FE baseline lacks its first reason")
            if fe["fe_target_error"]:
                raise ArtifactError("failed FE baseline fabricates a target")
        else:
            raise ArtifactError("invalid FE baseline status")

        if physical:
            for _, record in records:
                if fe_target is None:
                    _blank(record, ("fe_target_error",), "ineligible FE case")
                elif finite(record["fe_target_error"], "FE target error") != fe_target:
                    raise ArtifactError("physical record FE target is stale")

        for method in IMPLICIT_METHODS:
            if len(terminals[method]) != 1:
                raise ArtifactError("implicit search needs exactly one terminal")
            terminal_position, terminal = terminals[method][0]
            terminal_status = terminal["status"]
            if terminal_status not in TERMINAL_STATUSES:
                raise ArtifactError("unknown search terminal status")
            _blank(terminal, ERROR_COLUMNS + NUMERICAL_COLUMNS,
                   "search terminal")
            candidates = measurements[method]
            steps = [integer(row["steps"], "candidate steps")
                     for _, row in candidates]
            if any(step > MAX_STEPS for step in steps):
                raise ArtifactError("implicit candidate exceeds count cap")
            if steps != list(range(1, len(steps) + 1)):
                raise ArtifactError("implicit candidate counts are not consecutive")
            if candidates and terminal_position <= candidates[-1][0]:
                raise ArtifactError("search terminal precedes its candidates")
            for _, candidate in candidates:
                candidate_step = integer(candidate["steps"], "candidate steps")
                _validate_measurement_numerics(candidate)
                status = candidate["status"]
                if status in SUCCESS_STATUSES:
                    error = _validate_success(
                        candidate, integer(candidate["levels"], "levels"))
                    if candidate["failure_reason"]:
                        raise ArtifactError("eligible candidate has a failure reason")
                    if fe_target is None or finite(candidate["fe_target_error"],
                                                   "FE target error") != fe_target:
                        raise ArtifactError("candidate FE target is stale")
                    if status == "eligible_above_target" and not error > fe_target:
                        raise ArtifactError("above-target candidate meets target")
                    if status == "eligible_match" and not error <= fe_target:
                        raise ArtifactError("matched candidate exceeds target")
                    if status == "fe_baseline":
                        raise ArtifactError("implicit candidate labelled FE baseline")
                elif status in FAILURE_STATUSES:
                    _blank(candidate, ERROR_COLUMNS, "failed candidate")
                    if not candidate["failure_reason"]:
                        raise ArtifactError("failed candidate lacks its first reason")
                    expected_target = "" if fe_target is None else str(fe_target)
                    if candidate["fe_target_error"] and \
                            finite(candidate["fe_target_error"], "FE target") != fe_target:
                        raise ArtifactError("failed candidate FE target is stale")
                else:
                    raise ArtifactError("unknown candidate status")
                key = (fixture, method, candidate_step)
                if key in indexed:
                    raise ArtifactError("duplicate measurement key")
                indexed[key] = candidate

            if terminal_status == "selected":
                if not candidates or terminal["failure_reason"]:
                    raise ArtifactError("selected terminal is incomplete")
                terminal_step = integer(terminal["steps"], "terminal steps")
                if terminal_step != steps[-1] or \
                        candidates[-1][1]["status"] != "eligible_match" or \
                        any(row["status"] == "eligible_match"
                            for _, row in candidates[:-1]):
                    raise ArtifactError("selected terminal is not the first match")
                if any(row["status"] not in
                       ({"eligible_above_target"} | FAILURE_STATUSES)
                       for _, row in candidates[:-1]):
                    raise ArtifactError("preselection candidate status is invalid")
            else:
                if not terminal["failure_reason"]:
                    raise ArtifactError("nonselected terminal lacks its first reason")
                if any(row["status"] == "eligible_match" for _, row in candidates):
                    raise ArtifactError("nonselected search contains a match")
                if terminal_status == "not_attained_within_cap":
                    if steps != list(range(1, MAX_STEPS + 1)) or \
                            integer(terminal["steps"], "terminal steps") != MAX_STEPS:
                        raise ArtifactError("cap terminal did not exhaust the cap")
                elif terminal_status == "accuracy_budget_exhausted":
                    if candidates:
                        if integer(terminal["steps"], "terminal steps") != steps[-1]:
                            raise ArtifactError("budget terminal step mismatch")
                    elif terminal["steps"] not in {"", "0"}:
                        raise ArtifactError("empty budget search has a terminal step")
                else:
                    if candidates or terminal["steps"] not in {"", "0"}:
                        raise ArtifactError("pre-search terminal contains candidates")
                    expected_fe_status = {
                        "reference_unresolved": "reference_unresolved",
                        "fe_baseline_ineligible": None,
                        "fixture_controller_mismatch":
                            "fixture_controller_mismatch",
                    }[terminal_status]
                    if terminal_status == "reference_unresolved" and \
                            reference_status != "unresolved":
                        raise ArtifactError("reference terminal has certified reference")
                    if terminal_status == "fe_baseline_ineligible" and \
                            fe_status not in FAILURE_STATUSES:
                        raise ArtifactError("FE-ineligible terminal has eligible FE")
                    if expected_fe_status and fe_status != expected_fe_status:
                        raise ArtifactError("terminal and FE failure disagree")

            confirmation_steps: set[int] = set()
            for position, confirmation in confirmations[method]:
                if terminal_status != "selected" or position <= terminal_position:
                    raise ArtifactError("confirmation can affect search selection")
                step = integer(confirmation["steps"], "confirmation steps")
                if step > MAX_STEPS:
                    raise ArtifactError("confirmation exceeds count cap")
                selected_step = steps[-1]
                if step not in {2 * selected_step, 4 * selected_step} or \
                        step in confirmation_steps:
                    raise ArtifactError("confirmation step is invalid")
                confirmation_steps.add(step)
                _validate_measurement_numerics(confirmation)
                status = confirmation["status"]
                if status in ((SUCCESS_STATUSES - {"fe_baseline"}) |
                              {CONFIRMATION_SUCCESS}):
                    error = _validate_success(
                        confirmation, integer(confirmation["levels"], "levels"))
                    if confirmation["failure_reason"]:
                        raise ArtifactError("eligible confirmation has failure reason")
                    if finite(confirmation["fe_target_error"], "FE target") != fe_target:
                        raise ArtifactError("confirmation FE target is stale")
                    if status == "eligible_match" and error > fe_target:
                        raise ArtifactError("matched confirmation exceeds target")
                    if status == "eligible_above_target" and error <= fe_target:
                        raise ArtifactError("above-target confirmation meets target")
                elif status in FAILURE_STATUSES:
                    _blank(confirmation, ERROR_COLUMNS, "failed confirmation")
                    if not confirmation["failure_reason"]:
                        raise ArtifactError("failed confirmation lacks first reason")
                else:
                    raise ArtifactError("invalid confirmation status")

        fe_key = (fixture, "forward_euler", fe_steps)
        confirmation_steps = set()
        for position, confirmation in confirmations["forward_euler"]:
            _validate_measurement_numerics(confirmation)
            step = integer(confirmation["steps"], "FE confirmation steps")
            if (fe_target is None or position <= measurements["forward_euler"][0][0]
                    or step not in {2 * fe_steps, 4 * fe_steps}
                    or step > MAX_STEPS or step in confirmation_steps
                    or finite(confirmation["fe_target_error"], "FE target") != fe_target):
                raise ArtifactError("FE confirmation is inconsistent")
            confirmation_steps.add(step)
            if confirmation["status"] == CONFIRMATION_SUCCESS:
                _validate_success(confirmation, integer(fe["levels"], "levels"))
                if confirmation["failure_reason"]:
                    raise ArtifactError("eligible FE confirmation has failure reason")
            elif confirmation["status"] in FAILURE_STATUSES:
                _blank(confirmation, ERROR_COLUMNS, "failed FE confirmation")
                if not confirmation["failure_reason"]:
                    raise ArtifactError("failed FE confirmation lacks first reason")
            else:
                raise ArtifactError("invalid FE confirmation status")
        if fe_key in indexed:
            raise ArtifactError("duplicate FE measurement")
        indexed[fe_key] = fe
    return indexed


def git_checkout(path: Path) -> dict[str, Any]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=path, check=True,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=path,
                                check=True, capture_output=True,
                                text=True).stdout)
    return {"head": head, "dirty": dirty}


def build_identity(binary: Path) -> dict[str, Any]:
    cache = binary.parent / "CMakeCache.txt"
    if not cache.is_file():
        raise ArtifactError("the exact timing binary needs its CMake cache")
    entries = {}
    for line in cache.read_text(encoding="utf-8").splitlines():
        if ":" in line and "=" in line and not line.startswith(("#", "//")):
            key_type, value = line.split("=", 1)
            key, _ = key_type.split(":", 1)
            entries[key] = value
    required = ("EKAT_SOURCE", "MAM4XX_SOURCE", "Kokkos_DIR",
                "CMAKE_CXX_COMPILER", "CMAKE_BUILD_TYPE",
                "NDROP_ENABLE_SANITIZERS")
    if any(not entries.get(key) for key in required):
        raise ArtifactError("the exact timing build identity is incomplete")
    dependencies = {name: git_checkout(Path(entries[key]))
                    for name, key in (("ekat", "EKAT_SOURCE"),
                                      ("mam4xx", "MAM4XX_SOURCE"))}
    if dependencies["ekat"]["head"] != \
            "3a1485aecfed7a036e4fa30ce845498fd37975e2":
        raise ArtifactError("the EKAT source is not the locked current target")
    if dependencies["mam4xx"]["head"] != \
            "0d076d8ebb004a18f05f51f97b08edd9def6826f":
        raise ArtifactError("the MAM4xx source is not the locked current target")
    if any(record["dirty"] for record in dependencies.values()):
        raise ArtifactError("a locked dependency checkout is dirty")
    kokkos_config = Path(entries["Kokkos_DIR"]) / "KokkosConfig.cmake"
    if not kokkos_config.is_file():
        raise ArtifactError("the Kokkos build configuration is absent")
    return {"cmake_cache_sha256": sha256(cache),
            "compiler": entries["CMAKE_CXX_COMPILER"],
            "build_type": entries["CMAKE_BUILD_TYPE"],
            "sanitizers": entries["NDROP_ENABLE_SANITIZERS"],
            "kokkos_dir": entries["Kokkos_DIR"],
            "kokkos_config_sha256": sha256(kokkos_config),
            "dependencies": dependencies}


def validate_identity(directory: Path) -> dict[str, Any]:
    identity = load_json(directory / "identity.json")
    if identity.get("artifact_schema") not in {IDENTITY_SCHEMA, study_schema("identity", True),
                                              study_schema("identity", physical=True),
                                              study_schema("identity", physical=GEOMETRY_POLICY)}:
        raise ArtifactError("identity schema mismatch")
    controlled = controlled_identity(identity)
    physical = physical_identity(identity)
    if controlled and set(identity) != {
            "artifact_schema", "binary_sha256", "source_head", "source_root",
            "source_status", "source_attached_head", "source_files_sha256", "build",
            "design_sha256", "python_version", "platform", "fixtures",
            "accuracy_command", "accuracy_protocol", "timing_protocol",
            "fixture_policy", "activation_mixing_scale"}:
        raise ArtifactError("wrong closed controlled identity keys")
    if physical and set(identity) != {
            "artifact_schema", "binary_sha256", "source_head", "source_root",
            "source_status", "source_attached_head", "source_files_sha256", "build",
            "design_sha256", "python_version", "platform", "fixtures",
            "accuracy_command", "accuracy_protocol", "timing_protocol",
            "fixture_policy", "comparison_policy", "grid_files_sha256"}:
        raise ArtifactError("wrong closed physical identity keys")
    binary = directory / "inputs" / "benchmark_ndrop_matched_error"
    if sha256(binary) != identity.get("binary_sha256"):
        raise ArtifactError("timing binary differs from identity")
    if identity.get("python_version") != platform.python_version():
        raise ArtifactError("executing Python version differs from identity")
    root_text = identity.get("source_root")
    if not isinstance(root_text, str):
        raise ArtifactError("identity lacks its source root")
    if physical and not Path(root_text).is_absolute():
        raise ArtifactError("physical source root must be absolute")
    source_dir = Path(root_text) / "docs/analysis/evidence/ndrop_integrator"
    expected_sources = identity.get("source_files_sha256")
    if not isinstance(expected_sources, dict) or \
            set(expected_sources) != set(source_files(physical)):
        raise ArtifactError("identity source file set is incomplete")
    for name, expected in expected_sources.items():
        if sha256(source_dir / name) != expected:
            raise ArtifactError(f"source file differs from identity: {name}")
    design = Path(root_text) / "docs/design/eamxx_mam4xx_ndrop.md"
    if sha256(design) != identity.get("design_sha256"):
        raise ArtifactError("design contract differs from identity")
    if physical:
        if set(identity["grid_files_sha256"]) != set(PHYSICAL_GRIDS):
            raise ArtifactError("physical grid set incomplete")
        for name, digest in identity["grid_files_sha256"].items():
            if (sha256(directory / "inputs/grids" / name) != digest or
                    sha256(Path(root_text) / "components/homme/test/vcoord" / name) != digest):
                raise ArtifactError("physical grid copy/source hash mismatch")
        if (identity["accuracy_command"] != ["inputs/benchmark_ndrop_matched_error",
                physical_accuracy_option(physical), "--physical-grid-dir", "inputs/grids",
                "--max-steps", "8192", "--budget-seconds", "7200"] or
                identity["accuracy_protocol"] != {"integer_search": True,
                    "maximum_steps": MAX_STEPS, "budget_seconds": ACCURACY_BUDGET_SECONDS,
                    "confirmation_multipliers": [2, 4]} or
                identity["timing_protocol"] != {"warmups": 10, "minimum_blocks": 30,
                    "maximum_blocks": 100, "budget_seconds": 7200,
                    "complete_interleaved_blocks": True}):
            raise ArtifactError("physical protocol changed")
        dependencies = identity["build"]["dependencies"]
        for name, commit in (("ekat", "3a1485aecfed7a036e4fa30ce845498fd37975e2"),
                             ("mam4xx", "0d076d8ebb004a18f05f51f97b08edd9def6826f")):
            if dependencies[name] != {"head": commit, "dirty": False}:
                raise ArtifactError("physical locked dependency mismatch")
        if (identity["source_status"] not in {"clean_detached_source", "exploratory_source"}
                or identity["source_status"] == "clean_detached_source" and identity["source_attached_head"]):
            raise ArtifactError("physical source status invalid")
    identity_fixtures(identity)
    return identity


def identity_fixtures(identity: dict[str, Any]) -> tuple[int, ...]:
    physical = physical_identity(identity)
    if physical:
        values = identity.get("fixtures")
        if (not isinstance(values, list) or any(type(v) is not int for v in values)
                or values != list(physical_cases(physical))):
            raise ArtifactError("physical controller requires exact policy coverage")
        return tuple(values)
    if controlled_identity(identity):
        values = identity.get("fixtures")
        if (not isinstance(values, list) or any(type(value) is not int for value in values)
                or values != list(CONTROLLED_FIXTURES)):
            raise ArtifactError("controlled mixing requires exact four-case coverage")
        return tuple(values)
    fixtures = identity.get("fixtures", list(FIXTURES))
    if (not isinstance(fixtures, list) or not fixtures or
            any(type(value) is not int or value not in FIXTURES for value in fixtures) or
            fixtures != sorted(set(fixtures))):
        raise ArtifactError("identity fixture set is invalid")
    return tuple(fixtures)


def completed(directory: Path, stem: str, controlled: bool = False, physical: bool | str = False) -> Path:
    artifact = directory / f"{stem}.csv"
    identity = directory / "identity.json"
    record = load_json(directory / f"{stem}.complete.json")
    if record != {"artifact_schema": study_schema("completion", controlled, physical),
                  "artifact_sha256": sha256(artifact),
                  "identity_sha256": sha256(identity),
                  "path": artifact.name}:
        raise ArtifactError(f"{stem} completion is missing or stale")
    return artifact


def _completion(directory: Path, stem: str, controlled: bool = False, physical: bool | str = False) -> None:
    artifact = directory / f"{stem}.csv"
    record = {"artifact_schema": study_schema("completion", controlled, physical),
              "artifact_sha256": sha256(artifact),
              "identity_sha256": sha256(directory / "identity.json"),
              "path": artifact.name}
    atomic_new(directory / f"{stem}.complete.json", canonical_json(record))


def run_accuracy(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    if directory.exists():
        raise ArtifactError("study directory must be new")
    root = args.source_root.resolve()
    binary = args.binary.resolve()
    if not root.is_dir() or not binary.is_file():
        raise ArtifactError("source or binary is absent")
    source_dir = root / "docs/analysis/evidence/ndrop_integrator"
    physical = GEOMETRY_POLICY if getattr(args, "physical_geometry", False) else \
        getattr(args, "physical_columns", False)
    source_hashes = {name: sha256(source_dir / name) for name in source_files(physical)}
    build = build_identity(binary)
    controlled = args.controlled_mixing
    if physical and (controlled or args.fixtures):
        raise ArtifactError("physical controller requires its exact case policy")
    fixture_identity = {"fixtures": [int(value) for value in args.fixtures.split(",")]
                        if args.fixtures else list(physical_cases(physical) if physical else
                                                  CONTROLLED_FIXTURES if controlled else FIXTURES)}
    if physical:
        fixture_identity.update(artifact_schema=study_schema("identity", physical=physical),
                                fixture_policy=fixture_policy(physical),
                                comparison_policy="natural_fe_controller_error_v1")
    if controlled:
        fixture_identity.update(artifact_schema=study_schema("identity", True),
                                fixture_policy="controlled_mixing_v1",
                                activation_mixing_scale=CONTROLLED_ACTIVATION_MIXING)
    fixtures = identity_fixtures(fixture_identity)
    status = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                            check=True, capture_output=True, text=True).stdout
    attached = subprocess.run(["git", "symbolic-ref", "-q", "HEAD"], cwd=root,
                              capture_output=True).returncode == 0
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, check=True,
                          capture_output=True, text=True).stdout.strip()
    directory.mkdir(parents=True)
    (directory / "inputs").mkdir()
    retained = directory / "inputs" / "benchmark_ndrop_matched_error"
    with binary.open("rb") as source, retained.open("xb") as target:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            target.write(chunk)
        target.flush()
        os.fsync(target.fileno())
    retained.chmod(binary.stat().st_mode & 0o777)
    if sha256(retained) != sha256(binary):
        raise ArtifactError("binary copy changed")
    grid_hashes = {}
    if physical:
        (directory / "inputs/grids").mkdir()
        for name in PHYSICAL_GRIDS:
            source = root / "components/homme/test/vcoord" / name
            target = directory / "inputs/grids" / name
            with target.open("xb") as stream:
                stream.write(source.read_bytes()); stream.flush(); os.fsync(stream.fileno())
            grid_hashes[name] = sha256(source)
            if sha256(target) != grid_hashes[name]:
                raise ArtifactError("physical grid copy changed")
    accuracy_command = ["inputs/benchmark_ndrop_matched_error",
                        "--controller-accuracy", "--max-steps", str(MAX_STEPS),
                        "--budget-seconds", str(ACCURACY_BUDGET_SECONDS),
                        "--fixtures", args.fixtures]
    if physical:
        accuracy_command = ["inputs/benchmark_ndrop_matched_error", physical_accuracy_option(physical),
                            "--physical-grid-dir", "inputs/grids", "--max-steps", str(MAX_STEPS),
                            "--budget-seconds", str(ACCURACY_BUDGET_SECONDS)]
    elif controlled:
        accuracy_command = ["inputs/benchmark_ndrop_matched_error",
                            "--controlled-mixing-accuracy", "--max-steps", str(MAX_STEPS),
                            "--budget-seconds", str(ACCURACY_BUDGET_SECONDS)]
    elif not args.fixtures:
        accuracy_command[-1] = ",".join(map(str, fixtures))
    identity = {
        "artifact_schema": study_schema("identity", controlled, physical),
        "binary_sha256": sha256(retained), "source_head": head,
        "source_root": str(root),
        "source_status": "clean_detached_source" if not status and not attached
                         else "exploratory_source",
        "source_attached_head": attached,
        "source_files_sha256": source_hashes, "build": build,
        "design_sha256": sha256(root / "docs/design/eamxx_mam4xx_ndrop.md"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "fixtures": list(fixtures),
        "accuracy_command": accuracy_command,
        "accuracy_protocol": {"integer_search": True,
                              "maximum_steps": MAX_STEPS,
                              "budget_seconds": ACCURACY_BUDGET_SECONDS,
                              "confirmation_multipliers": [2, 4]},
        "timing_protocol": {"warmups": 10, "minimum_blocks": 30,
                            "maximum_blocks": 100, "budget_seconds": 7200,
                            "complete_interleaved_blocks": True},
    }
    if physical:
        identity.update(fixture_policy=fixture_policy(physical),
                        comparison_policy="natural_fe_controller_error_v1", grid_files_sha256=grid_hashes)
    elif controlled:
        identity.update(fixture_policy="controlled_mixing_v1",
                        activation_mixing_scale=CONTROLLED_ACTIVATION_MIXING)
    atomic_new(directory / "identity.json", canonical_json(identity))
    if physical:
        validate_identity(directory)
    process = subprocess.run([str(retained), *accuracy_command[1:]],
                             capture_output=True, check=True, cwd=directory if physical else None)
    output = directory / "accuracy.csv"
    atomic_new(output, process.stdout)
    validate_accuracy(csv_rows(output, accuracy_header(controlled, physical)), fixtures, controlled, physical)
    _completion(directory, "accuracy", controlled, physical)


def expected_selection(
    indexed: dict[tuple[int, str, int], dict[str, str]],
    identity_sha: str,
    accuracy_sha: str,
    controlled: bool = False,
    physical: bool | str = False,
) -> dict[str, Any]:
    entries = []
    unique: set[tuple[int, str, int]] = set()
    for fixture in sorted({key[0] for key in indexed}):
        fe_rows = [(key, row) for key, row in indexed.items()
                   if key[0] == fixture and key[1] == "forward_euler"]
        if len(fe_rows) != 1:
            raise ArtifactError("selection lacks unique FE baseline")
        fe_key, fe = fe_rows[0]
        target = float(fe["fe_target_error"]) if fe["fe_target_error"] else None
        fe_selected = fe["status"] == "fe_baseline"
        fe_error = float(fe["error_max"]) if fe_selected else None
        entries.append({"fixture": fixture, "method": "forward_euler",
                        "status": "selected" if fe_selected else fe["status"],
                        "steps": fe_key[2] if fe_selected else None,
                        "fe_target_error": target,
                        "achieved_error": fe_error,
                        "achieved_error_to_fe_target_ratio":
                            fe_error / target if target and fe_error is not None else None})
        if fe_selected:
            unique.add(fe_key)
        for method in IMPLICIT_METHODS:
            candidates = [(key, row) for key, row in indexed.items()
                          if key[0] == fixture and key[1] == method]
            matched = [(key, row) for key, row in candidates
                       if row["status"] == "eligible_match"]
            if len(matched) > 1:
                raise ArtifactError("selection contains multiple first matches")
            if matched:
                key, row = matched[0]
                error = float(row["error_max"])
                status, steps = "selected", key[2]
                unique.add(key)
            else:
                # The validator proved there is exactly one terminal.  Recover
                # its fail-closed status from the candidate-free FE condition or
                # the last candidate's failure mode is insufficient, so callers
                # should use expected_selection_from_rows when freezing.
                error, status, steps = None, "not_selected", None
            entries.append({"fixture": fixture, "method": method,
                            "status": status, "steps": steps,
                            "fe_target_error": target,
                            "achieved_error": error,
                            "achieved_error_to_fe_target_ratio":
                                error / target if target and error is not None else None})
    return {"artifact_schema": study_schema("selection", controlled, physical),
            "identity_sha256": identity_sha,
            "accuracy_sha256": accuracy_sha,
            "selected": entries,
            "unique_configs": [list(config) for config in sorted(unique)]}


def expected_selection_from_rows(
    rows: list[dict[str, str]],
    indexed: dict[tuple[int, str, int], dict[str, str]],
    identity_sha: str,
    accuracy_sha: str,
    controlled: bool = False,
    physical: bool | str = False,
) -> dict[str, Any]:
    frozen = expected_selection(indexed, identity_sha, accuracy_sha, controlled, physical)
    terminals = {(int(row["fixture"]), row["method"]): row["status"]
                 for row in rows if row["record_kind"] == "search_terminal"}
    for entry in frozen["selected"]:
        if entry["method"] in IMPLICIT_METHODS and entry["status"] == "not_selected":
            entry["status"] = terminals[(entry["fixture"], entry["method"])]
    return frozen


def freeze(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    identity = validate_identity(directory)
    controlled = controlled_identity(identity)
    physical = physical_identity(identity)
    accuracy_path = completed(directory, "accuracy", controlled, physical)
    rows = csv_rows(accuracy_path, accuracy_header(controlled, physical))
    indexed = validate_accuracy(rows, identity_fixtures(identity), controlled, physical)
    frozen = expected_selection_from_rows(
        rows, indexed, sha256(directory / "identity.json"), sha256(accuracy_path), controlled, physical)
    atomic_new(directory / "selection.json", canonical_json(frozen))
    atomic_new(directory / "timing_admission.json", canonical_json(
        expected_admission(rows, indexed, frozen,
                           sha256(directory / "selection.json"), controlled, physical)))


def selection(directory: Path) -> dict[str, Any]:
    identity = validate_identity(directory)
    controlled = controlled_identity(identity)
    physical = physical_identity(identity)
    frozen = load_json(directory / "selection.json")
    accuracy_path = completed(directory, "accuracy", controlled, physical)
    rows = csv_rows(accuracy_path, accuracy_header(controlled, physical))
    indexed = validate_accuracy(rows, identity_fixtures(identity), controlled, physical)
    expected = expected_selection_from_rows(
        rows, indexed, sha256(directory / "identity.json"), sha256(accuracy_path), controlled, physical)
    if frozen != json.loads(canonical_json(expected)):
        raise ArtifactError("selection is stale or malformed")
    return frozen


def expected_admission(rows, indexed, frozen, selection_sha, controlled=False, physical=False):
    """Refinements cannot select work; they can only withhold timing.

    Strict decrease against BOTH certified references is required in each
    aggregate norm. Orders are reported, not forced into a new slope window.
    Formal first order is established separately by the unchanged unit gates.
    """
    confirmations = {(int(row["fixture"]), row["method"], int(row["steps"])): row
                     for row in rows if row["record_kind"] == "confirmation"}
    entries, configs = [], []
    for config in frozen["unique_configs"]:
        fixture, method, steps = config
        base = indexed[tuple(config)]
        refined = [confirmations.get((fixture, method, factor * steps))
                   for factor in (2, 4)]
        reason = ""
        orders = {}
        if any(row is None for row in refined):
            reason = "missing_refinement_evidence"
        elif any(row["status"] != CONFIRMATION_SUCCESS for row in refined):
            reason = "ineligible_refinement"
        else:
            for norm in ("normal_weighted_l1_by_field", "normal_linf_by_field",
                         "tight_weighted_l1_by_field", "tight_linf_by_field"):
                errors = [max(map(float, row[norm].split(";")))
                          for row in (base, *refined)]
                orders[norm] = {"errors_at_n_2n_4n": errors,
                                "observed_orders": [math.log2(errors[i] / errors[i + 1])
                                                    if errors[i] > 0 and errors[i + 1] > 0
                                                    else None for i in (0, 1)]}
                if not errors[0] > errors[1] > errors[2] > 0:
                    reason = "reference_error_does_not_decrease"
        status = "admitted" if not reason else "excluded"
        entries.append({"fixture": fixture, "method": method, "steps": steps,
                        "status": status, "reason": reason, "convergence": orders})
        if not reason:
            configs.append(config)
    return {"artifact_schema": study_schema("timing_admission", controlled, physical),
            "selection_sha256": selection_sha,
            "accuracy_sha256": frozen["accuracy_sha256"],
            "rule": "strict_decrease_both_references_both_aggregate_norms_n_2n_4n",
            "entries": entries, "unique_configs": configs}


def admission(directory: Path, frozen: dict[str, Any]) -> dict[str, Any]:
    identity = validate_identity(directory)
    controlled = controlled_identity(identity)
    physical = physical_identity(identity)
    rows = csv_rows(completed(directory, "accuracy", controlled, physical), accuracy_header(controlled, physical))
    indexed = validate_accuracy(rows, identity_fixtures(identity), controlled, physical)
    expected = expected_admission(rows, indexed, frozen,
                                  sha256(directory / "selection.json"), controlled, physical)
    actual = load_json(directory / "timing_admission.json")
    if actual != expected:
        raise ArtifactError("timing admission is stale or malformed")
    return actual


def validate_timing(
    rows: list[dict[str, str]],
    frozen: dict[str, Any],
    accuracy: dict[tuple[int, str, int], dict[str, str]],
    controlled: bool = False,
    physical: bool | str = False,
) -> None:
    expected = {tuple(config) for config in frozen["unique_configs"]}
    groups: dict[tuple[int, str, int], list[int]] = {}
    checksums: dict[tuple[int, str, int], set[int]] = {}
    reasons: dict[int, set[str]] = {}
    fixture_blocks: dict[tuple[int, int], set[tuple[int, str, int]]] = {}
    for row in rows:
        if row["artifact_schema"] != study_schema("timing", controlled, physical):
            raise ArtifactError("timing schema mismatch")
        key = (integer(row["fixture"], "fixture"), row["method"],
               integer(row["steps"], "steps"))
        if key not in expected:
            raise ArtifactError("timed configuration was not frozen")
        block = integer(row["block"], "block")
        if block < 0 or finite(row["seconds_per_column"], "time") <= 0:
            raise ArtifactError("timing block or duration is invalid")
        groups.setdefault(key, []).append(block)
        checksum = integer(row["state_checksum"], "checksum")
        if checksum != integer(accuracy[key]["state_checksum"],
                               "accuracy checksum"):
            raise ArtifactError("timed state differs from qualified accuracy state")
        checksums.setdefault(key, set()).add(checksum)
        reasons.setdefault(key[0], set()).add(row["stop_reason"])
        block_key = (key[0], block)
        if key in fixture_blocks.setdefault(block_key, set()):
            raise ArtifactError("duplicate timed block")
        fixture_blocks[block_key].add(key)
    if set(groups) != expected:
        raise ArtifactError("missing selected timed configuration")
    for key, blocks in groups.items():
        if blocks != list(range(len(blocks))) or len(checksums[key]) != 1 or \
                not 30 <= len(blocks) <= 100:
            raise ArtifactError("incomplete or inconsistent timing blocks")
    for fixture in sorted({key[0] for key in expected}):
        fixture_expected = {key for key in expected if key[0] == fixture}
        if not fixture_expected:
            # A reference/controller/FE failure intentionally contributes no
            # timing configuration.  Its unavailable comparisons remain in
            # the frozen selection and analysis rather than being fabricated.
            continue
        counts = {len(blocks) for key, blocks in groups.items()
                  if key[0] == fixture}
        if len(counts) != 1 or \
                len(reasons.get(fixture, set())) != 1 or \
                not reasons[fixture] <= {"median_precision", "maximum_blocks",
                                         "budget_exhausted"}:
            raise ArtifactError("timing methods are not paired by complete blocks")
        count = next(iter(counts))
        if any(fixture_blocks.get((fixture, block), set()) != fixture_expected
               for block in range(count)):
            raise ArtifactError("timing block is not a complete method pair")


def timing_command(directory: Path, admitted: dict, controlled=False, physical=False) -> list[str]:
    encoded = ";".join(f"{f}:{m}:{n}" for f, m, n in admitted["unique_configs"])
    binary = directory / "inputs" / "benchmark_ndrop_matched_error"
    command = [str(binary), "--physical-geometry-controller-timing" if physical == GEOMETRY_POLICY else
               "--physical-controller-timing" if physical else
               "--controlled-mixing-timing" if controlled else "--controller-timing",
               "--selection", encoded, "--warmups", "10", "--min-blocks", "30",
               "--max-blocks", "100", "--budget-seconds", "7200"]
    if physical:
        command.extend(["--physical-grid-dir", str(directory / "inputs/grids")])
    return command


def expected_geometry_timing_approval(directory: Path, frozen: dict, admitted: dict) -> dict:
    """Bind this bounded approval, not a general authorization mechanism."""
    if (sha256(directory / "accuracy.csv") != GEOMETRY_APPROVED_ACCURACY_SHA256
            or len(admitted["unique_configs"]) != 36
            or admitted["unique_configs"] != frozen["unique_configs"]):
        raise ArtifactError("geometry numerical evidence changed; renewed approval required")
    return {"artifact_schema": study_schema("timing_approval", physical=GEOMETRY_POLICY),
            "approval_reference": "user_in_thread_2026-10-12_geometry_timing",
            "authority": "bounded_local_12_case_geometry_timing_and_analysis",
            "identity_sha256": sha256(directory / "identity.json"),
            "accuracy_sha256": sha256(directory / "accuracy.csv"),
            "selection_sha256": sha256(directory / "selection.json"),
            "timing_admission_sha256": sha256(directory / "timing_admission.json"),
            "prior_accuracy_sha256": GEOMETRY_APPROVED_ACCURACY_SHA256,
            "prior_selection_sha256": GEOMETRY_PRIOR_SELECTION_SHA256,
            "command": timing_command(directory, admitted, physical=GEOMETRY_POLICY)}


def geometry_timing_approval(directory: Path, frozen: dict, admitted: dict) -> None:
    try:
        record = load_json(directory / "timing_approval.json")
    except FileNotFoundError as exc:
        raise ArtifactError("geometry timing requires its separate approval record") from exc
    if canonical_json(record) != canonical_json(expected_geometry_timing_approval(directory, frozen, admitted)):
        raise ArtifactError("geometry timing approval is stale or outside scope")


def expected_geometry_execution(directory: Path, admitted: dict, rows: list[dict]) -> dict:
    return {"artifact_schema": study_schema("execution", physical=GEOMETRY_POLICY),
            "identity_sha256": sha256(directory / "identity.json"),
            "selection_sha256": sha256(directory / "selection.json"),
            "timing_approval_sha256": sha256(directory / "timing_approval.json"),
            "timing_sha256": sha256(directory / "timing.csv"),
            "command": timing_command(directory, admitted, physical=GEOMETRY_POLICY),
            "process_exit_code": 0, "timing_rows": len(rows),
            "warmups_by_configuration": [{"fixture": f, "method": m, "steps": n, "warmups": 10}
                                         for f, m, n in admitted["unique_configs"]]}


def validate_geometry_execution(directory: Path, admitted: dict, rows: list[dict]) -> None:
    if canonical_json(load_json(directory / "timing_execution.json")) != canonical_json(
            expected_geometry_execution(directory, admitted, rows)):
        raise ArtifactError("geometry timing execution evidence is stale or incomplete")


def evidence_scope(controlled=False, physical=False) -> str:
    return ("physical_geometry_controller_idealization_unweighted" if physical == GEOMETRY_POLICY else
            "physical_controller_idealization_unweighted" if physical else
            "synthetic_controlled_mixing_unweighted" if controlled else "synthetic_controller_strata_unweighted")


def run_timing(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    identity = validate_identity(directory)
    controlled = controlled_identity(identity)
    physical = physical_identity(identity)
    frozen = selection(directory)
    admitted = admission(directory, frozen)
    if physical == GEOMETRY_POLICY:
        geometry_timing_approval(directory, frozen, admitted)
        if any((directory / name).exists() for name in
               ("timing.csv", "timing.csv.partial", "timing.complete.json", "timing_execution.json")):
            raise ArtifactError("geometry timed pass already exists; renewed bounded approval required")
    accuracy = validate_accuracy(csv_rows(completed(directory, "accuracy", controlled, physical),
                                          accuracy_header(controlled, physical)), identity_fixtures(identity), controlled, physical)
    if identity["source_status"] not in {"clean_detached_source",
                                         "exploratory_source"}:
        raise ArtifactError("source status is invalid")
    if not admitted["unique_configs"]:
        raise ArtifactError("no eligible configuration to time")
    command = timing_command(directory, admitted, controlled, physical)
    output = directory / "timing.csv"
    if physical == GEOMETRY_POLICY:
        # Reserve before launch. A process/write/validation failure leaves an
        # inadmissible partial and consumes this directory's bounded attempt.
        partial = directory / "timing.csv.partial"
        with partial.open("xb") as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.PIPE, check=True)
            stream.flush()
            os.fsync(stream.fileno())
        rows = csv_rows(partial, TIMING_HEADER)
        validate_timing(rows, admitted, accuracy, controlled, physical)
        if output.exists():
            raise ArtifactError("geometry timing destination appeared during execution")
        os.replace(partial, output)
    else:
        process = subprocess.run(command, capture_output=True, check=True)
        atomic_new(output, process.stdout)
        rows = csv_rows(output, TIMING_HEADER)
        validate_timing(rows, admitted, accuracy, controlled, physical)
    _completion(directory, "timing", controlled, physical)
    if physical == GEOMETRY_POLICY:
        atomic_new(directory / "timing_execution.json", canonical_json(
            expected_geometry_execution(directory, admitted, rows)))


def median(values: list[float]) -> float:
    if not values:
        raise ArtifactError("median needs samples")
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return ordered[middle - 1] + (ordered[middle] - ordered[middle - 1]) / 2


def paired_ratio(first: list[float], second: list[float], seed: int) -> dict[str, float]:
    if len(first) != len(second) or len(first) < 30:
        raise ArtifactError("paired ratio needs 30 complete blocks")
    ratio = median(first) / median(second)
    rng = random.Random(seed)
    replicates = []
    for _ in range(10000):
        indexes = [rng.randrange(len(first)) for _ in first]
        replicates.append(median([first[index] for index in indexes]) /
                          median([second[index] for index in indexes]))
    replicates.sort()
    return {"ratio": ratio, "lower": replicates[249],
            "upper": replicates[9749], "blocks": len(first)}


def report(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    identity = validate_identity(directory)
    controlled = controlled_identity(identity)
    physical = physical_identity(identity)
    frozen = selection(directory)
    admitted = admission(directory, frozen)
    if physical == GEOMETRY_POLICY:
        geometry_timing_approval(directory, frozen, admitted)
    accuracy = validate_accuracy(csv_rows(completed(directory, "accuracy", controlled, physical),
                                          accuracy_header(controlled, physical)), identity_fixtures(identity), controlled, physical)
    rows = csv_rows(completed(directory, "timing", controlled, physical), TIMING_HEADER)
    validate_timing(rows, admitted, accuracy, controlled, physical)
    if physical == GEOMETRY_POLICY:
        validate_geometry_execution(directory, admitted, rows)
    by_key: dict[tuple[int, str, int], list[float]] = {}
    stop_reasons: dict[int, str] = {}
    for row in rows:
        key = (int(row["fixture"]), row["method"], int(row["steps"]))
        by_key.setdefault(key, []).append(float(row["seconds_per_column"]))
        stop_reasons[key[0]] = row["stop_reason"]
    configurations = [
        {"fixture": fixture, "method": method, "steps": steps,
         "achieved_error": float(accuracy[(fixture, method, steps)]["error_max"]),
         "fe_target_error": float(
             accuracy[(fixture, method, steps)]["fe_target_error"]),
         "achieved_error_to_fe_target_ratio": float(
             accuracy[(fixture, method, steps)]["error_max"]) / float(
                 accuracy[(fixture, method, steps)]["fe_target_error"]),
         "median_seconds_per_column": median(samples),
         "timing_blocks": len(samples), "stop_reason": stop_reasons[fixture]}
        for (fixture, method, steps), samples in sorted(by_key.items())
    ]
    chosen = {(f, m): n for f, m, n in admitted["unique_configs"]}
    comparisons = []
    for fixture in identity_fixtures(identity):
        for numerator, denominator in (("forward_euler", "imex"),
                                       ("forward_euler", "backward_euler"),
                                       ("imex", "backward_euler")):
            left = chosen.get((fixture, numerator))
            right = chosen.get((fixture, denominator))
            if left is None or right is None:
                comparisons.append({"fixture": fixture, "numerator": numerator,
                                    "denominator": denominator,
                                    "status": "missing_eligible_pair",
                                    "reason": "accuracy_or_refinement_ineligible"})
                continue
            ratio = paired_ratio(by_key[(fixture, numerator, left)],
                                 by_key[(fixture, denominator, right)],
                                 BOOTSTRAP_SEED)
            comparisons.append({"fixture": fixture, "numerator": numerator,
                                "denominator": denominator,
                                "numerator_steps": left,
                                "denominator_steps": right,
                                "status": "inconclusive" if ratio["lower"] <= 1 <= ratio["upper"]
                                else "interval_excludes_one", **ratio})
    analysis = {
        "artifact_schema": study_schema("analysis", controlled, physical),
        "identity_sha256": sha256(directory / "identity.json"),
        "accuracy_sha256": sha256(directory / "accuracy.csv"),
        "selection_sha256": sha256(directory / "selection.json"),
        "timing_admission_sha256": sha256(directory / "timing_admission.json"),
        "timing_sha256": sha256(directory / "timing.csv"),
        "source_status": identity["source_status"],
        "evidence_scope": evidence_scope(controlled, physical),
        "reference": {"relative_tolerances": [1e-12, 1e-13],
                      "self_consistency_limit": 2e-11,
                      "absolute_tolerance": "rtol*max(field_initial_max_abs,1e-30)"},
        "fixtures": [
            {"fixture": f, **{name: accuracy[(f, "forward_euler", next(
                key[2] for key in accuracy if key[:2] == (f, "forward_euler")))][name]
                for name in ("fixture_family", "requested_nsubmix", "levels",
                             "mixing_scale", "interval_seconds", "mixing_rate_t",
                             "activation_rate_t", "production_nsubmix",
                             "reference_status", "reference_disagreement") +
                             (("activation_mixing_scale",) if controlled else ()) +
                             (PHYSICAL_METADATA + ("activation_limiter_count",) if physical else ()) +
                             (GEOMETRY_METADATA if physical == GEOMETRY_POLICY else ())}}
            for f in identity_fixtures(identity)],
        "timing_protocol": identity["timing_protocol"],
        "refinement_admission": admitted["entries"],
        "selection_outcomes": frozen["selected"],
        "bootstrap": {"resamples": 10000, "seed": BOOTSTRAP_SEED,
                      "rng": "Python random.Random integer seed",
                      "median": "lower+(upper-lower)/2",
                      "endpoints": "zero_based_sorted_replicates_249_9749"},
        "configurations": configurations, "comparisons": comparisons,
    }
    if physical == GEOMETRY_POLICY:
        analysis.update(timing_approval_sha256=sha256(directory / "timing_approval.json"),
                        timing_execution_sha256=sha256(directory / "timing_execution.json"))
    atomic_new(directory / "analysis.json", canonical_json(analysis))
    files = ("identity.json", "accuracy.csv", "accuracy.complete.json",
             "selection.json", "timing_admission.json", "timing.csv", "timing.complete.json",
             "analysis.json", "inputs/benchmark_ndrop_matched_error")
    if physical:
        files += tuple("inputs/grids/" + name for name in PHYSICAL_GRIDS)
    if physical == GEOMETRY_POLICY:
        files += ("timing_approval.json", "timing_execution.json")
    final = {"artifact_schema": study_schema("final", controlled, physical),
             "source_status": identity["source_status"],
             "evidence_scope": analysis["evidence_scope"],
             "files_sha256": {name: sha256(directory / name) for name in files},
             "fixtures": list(identity_fixtures(identity)), "completed_timing_rows": len(rows)}
    atomic_new(directory / "final_manifest.json", canonical_json(final))


def validate_final(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    identity = validate_identity(directory)
    controlled = controlled_identity(identity)
    physical = physical_identity(identity)
    frozen = selection(directory)
    admitted = admission(directory, frozen)
    if physical == GEOMETRY_POLICY:
        geometry_timing_approval(directory, frozen, admitted)
    accuracy = validate_accuracy(csv_rows(completed(directory, "accuracy", controlled, physical),
                                          accuracy_header(controlled, physical)), identity_fixtures(identity), controlled, physical)
    rows = csv_rows(completed(directory, "timing", controlled, physical), TIMING_HEADER)
    validate_timing(rows, admitted, accuracy, controlled, physical)
    if physical == GEOMETRY_POLICY:
        validate_geometry_execution(directory, admitted, rows)
    final = load_json(directory / "final_manifest.json")
    if (final.get("artifact_schema") != study_schema("final", controlled, physical) or
            final.get("source_status") != identity["source_status"] or
            final.get("evidence_scope") != evidence_scope(controlled, physical) or
            final.get("completed_timing_rows") != len(rows) or
            final.get("fixtures") != list(identity_fixtures(identity))):
        raise ArtifactError("final manifest coverage or identity mismatch")
    required = {"identity.json", "accuracy.csv", "accuracy.complete.json",
                "selection.json", "timing_admission.json", "timing.csv", "timing.complete.json",
                "analysis.json", "inputs/benchmark_ndrop_matched_error"}
    if physical:
        required.update("inputs/grids/" + name for name in PHYSICAL_GRIDS)
    if physical == GEOMETRY_POLICY:
        required.update(("timing_approval.json", "timing_execution.json"))
    if set(final.get("files_sha256", {})) != required:
        raise ArtifactError("final manifest file set is incomplete")
    for name, expected in final["files_sha256"].items():
        if name.startswith("/") or ".." in Path(name).parts or \
                sha256(directory / name) != expected:
            raise ArtifactError("final manifest file hash mismatch")


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run-accuracy", "freeze", "run-timing", "report",
                 "validate-final"):
        command = commands.add_parser(name)
        command.add_argument("--output-dir", type=Path, required=True)
        if name == "run-accuracy":
            command.add_argument("--binary", type=Path, required=True)
            command.add_argument("--source-root", type=Path, required=True)
            command.add_argument("--fixtures", default="")
            policy = command.add_mutually_exclusive_group()
            policy.add_argument("--controlled-mixing", action="store_true",
                                 help="the exact four-case controlled mixing policy")
            policy.add_argument("--physical-columns", action="store_true",
                                help="unchanged physical columns matched to natural FE accuracy")
            policy.add_argument("--physical-geometry", action="store_true",
                                help="12 paired native/refined physical cases, untimed only")
    arguments = parser.parse_args()
    try:
        {"run-accuracy": run_accuracy, "freeze": freeze,
         "run-timing": run_timing, "report": report,
         "validate-final": validate_final}[arguments.command](arguments)
    except (ArtifactError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
