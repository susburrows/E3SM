"""Closed local artifacts for the standalone, synthetic matched-error study.

This tool never admits Task 5R or production replay artifacts. A study remains
exploratory if its source tree is dirty; it cannot acquire retained status by
renaming a directory or editing a status field.
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


TARGETS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
METHODS = ("forward_euler", "imex", "backward_euler")
FIXTURES = (0, 1, 2, 3)
ACCURACY_HEADER = (
    "artifact_schema", "fixture", "levels", "mixing_scale",
    "interval_seconds", "mixing_rate_t", "activation_rate_t",
    "method", "steps", "reference_status", "reference_reason",
    "reference_disagreement", "status", "failure_reason",
    "error_normal", "error_tight",
    "error_max", "normal_weighted_l1_by_field", "normal_linf_by_field",
    "tight_weighted_l1_by_field", "tight_linf_by_field",
    "richardson_to_next", "observed_order", "clipping_count",
    "material_negative_count", "normalized_residual",
    "normalized_ledger_error", "state_checksum",
)
TIMING_HEADER = (
    "artifact_schema", "fixture", "method", "steps", "block",
    "seconds_per_column", "state_checksum", "stop_reason",
)
SOURCES = (
    "ndrop_integrator.hpp", "ndrop_integrator.cpp", "ndrop_imex.cpp",
    "ndrop_backward_euler.cpp", "ndrop_reference.cpp", "ndrop_benchmark.hpp",
    "ndrop_benchmark.cpp", "benchmark_ndrop_matched_error.cpp",
    "ndrop_matched_artifacts.py", "CMakeLists.txt",
)

# Independent fixture/encoding family: none of these definitions admit old
# synthetic artifacts or production captures.
PHYSICAL_TARGETS = (1e-2, 1e-3)
PHYSICAL_FIXTURES = {28: (72, "moderate", 3., .3),
                     29: (72, "strong", 30., 1.),
                     30: (128, "moderate", 3., .3),
                     31: (128, "strong", 30., 1.)}
PHYSICAL_GRIDS = ("acme-72i.ascii", "acme-72m.ascii",
                  "scream-128i.ascii", "scream-128m.ascii")
PHYSICAL_ERROR_VECTORS = ("normal_weighted_l1_by_field", "normal_linf_by_field",
                          "tight_weighted_l1_by_field", "tight_linf_by_field")
PHYSICAL_ACCURACY_HEADER = (
    "artifact_schema", "record_kind", "fixture", "levels", "turbulence",
    "background_diffusivity_m2_s", "cloud_base_wtke_m_s", "interval_seconds",
    "mixing_rate_t", "activation_rate_t", "controller_rate_t",
    "cloud_base_level", "cloud_base_spacing_m", "limiting_level",
    "limiting_spacing_m", "method", "steps", "production_nsubmix",
    "reference_status", "reference_reason", "reference_disagreement", "status",
    "failure_reason", "error_normal", "error_tight", "error_max",
    *PHYSICAL_ERROR_VECTORS, "clipping_count", "material_negative_count",
    "normalized_residual", "normalized_ledger_error", "state_checksum", "activation_limiter_count")
PHYSICAL_NUMERICAL = ("clipping_count", "material_negative_count",
                       "normalized_residual", "normalized_ledger_error", "state_checksum")


def validate_physical_accuracy(rows):
    """Validate fabricated test fixtures or executed offline evidence alike."""
    indexed, identities, references = {}, {}, {}
    for row in rows:
        if set(row) != set(PHYSICAL_ACCURACY_HEADER) or \
                row["artifact_schema"] != "ndrop_physical_column_accuracy_v1":
            raise ArtifactError("wrong physical-column accuracy schema")
        fixture = integer(row["fixture"], "fixture")
        if fixture not in PHYSICAL_FIXTURES:
            raise ArtifactError("unknown physical-column fixture")
        levels, strength, diffusivity, velocity = PHYSICAL_FIXTURES[fixture]
        if (integer(row["levels"], "levels"), row["turbulence"],
                finite(row["background_diffusivity_m2_s"], "diffusivity"),
                finite(row["cloud_base_wtke_m_s"], "velocity")) != \
                (levels, strength, diffusivity, velocity) or \
                finite(row["interval_seconds"], "interval") != 150:
            raise ArtifactError("physical fixture identity mismatch")
        identity_names = PHYSICAL_ACCURACY_HEADER[3:15] + ("production_nsubmix", "activation_limiter_count")
        identity = tuple(row[name] for name in identity_names)
        if identities.setdefault(fixture, identity) != identity:
            raise ArtifactError("physical fixture characteristics changed between rows")
        if integer(row["activation_limiter_count"], "activation limiter count") < 0:
            raise ArtifactError("negative activation limiter count")
        rate_t = finite(row["controller_rate_t"], "controller stiffness")
        count = integer(row["production_nsubmix"], "controller count")
        if rate_t <= 0 or count != math.floor(rate_t / .9) + 1:
            raise ArtifactError("physical controller count inconsistent with stiffness")
        for name in ("mixing_rate_t", "activation_rate_t"):
            if finite(row[name], name) < 0:
                raise ArtifactError("negative physical stiffness")
        for name in ("cloud_base_level", "limiting_level"):
            if not 0 <= integer(row[name], name) < levels:
                raise ArtifactError("physical level outside grid")
        if integer(row["cloud_base_level"], "cloud base") == levels - 1:
            raise ArtifactError("cloud base has no lower interface")
        for name in ("cloud_base_spacing_m", "limiting_spacing_m"):
            if finite(row[name], name) <= 0:
                raise ArtifactError("nonpositive physical spacing")
        kind, method = row["record_kind"], row["method"]
        steps = integer(row["steps"], "steps")
        if kind not in {"baseline", "measurement", "confirmation", "terminal"} or \
                method not in METHODS or steps < (0 if kind == "terminal" else 1):
            raise ArtifactError("physical record kind/method/count invalid")
        if kind == "baseline" and (method != "forward_euler" or steps != count):
            raise ArtifactError("physical baseline is not controller-selected FE")
        if kind == "measurement" and (steps > 8192 or
                method == "forward_euler" and steps < count):
            raise ArtifactError("physical search outside bounds")
        key = kind, fixture, method, steps
        if key in indexed:
            raise ArtifactError("duplicate physical accuracy key")
        certified = row["reference_status"] == "certified"
        ref = tuple(row[name] for name in
                    ("reference_status", "reference_reason", "reference_disagreement"))
        if references.setdefault(fixture, ref) != ref:
            raise ArtifactError("physical reference evidence changed between rows")
        if certified:
            if row["reference_reason"] or not 0 <= finite(
                    row["reference_disagreement"], "reference self-check") <= 2e-11:
                raise ArtifactError("physical reference certification invalid")
        elif row["reference_status"] != "unresolved" or not row["reference_reason"]:
            raise ArtifactError("physical reference outcome invalid")
        terminal_statuses = {"selected", "not_attained_within_cap",
                             "accuracy_budget_exhausted", "reference_unresolved"}
        failure_statuses = {"numerical_failure", "fe_clipping", "material_negative",
                            "reference_unresolved", "fixture_controller_mismatch",
                            "accuracy_budget_exhausted", "confirmation_count_exceeds_cap"}
        allowed = terminal_statuses if kind == "terminal" else {"eligible"} | failure_statuses
        if row["status"] not in allowed or not certified and \
                row["status"] != "reference_unresolved" or certified and \
                row["status"] == "reference_unresolved":
            raise ArtifactError("physical result/reference applicability mismatch")
        error_names = ("error_normal", "error_tight", "error_max", *PHYSICAL_ERROR_VECTORS)
        if row["status"] == "eligible":
            if row["failure_reason"]:
                raise ArtifactError("eligible physical result has failure reason")
            if not 0 <= integer(row["state_checksum"], "eligible checksum") < 2**64:
                raise ArtifactError("eligible physical checksum outside uint64")
            maxima = []
            for name in PHYSICAL_ERROR_VECTORS:
                values = [finite(value, name) for value in row[name].split(";")]
                if len(values) != 51 or min(values) < 0:
                    raise ArtifactError("invalid 51-field physical error vector")
                maxima.append(max(values))
            if [finite(row[name], name) for name in error_names[:3]] != \
                    [max(maxima[:2]), max(maxima[2:]), max(maxima)]:
                raise ArtifactError("physical error maxima inconsistent")
            if integer(row["clipping_count"], "clipping") != 0 or \
                    integer(row["material_negative_count"], "negative count") != 0 or \
                    not 0 <= finite(row["normalized_residual"], "residual") <= 1e-11 or \
                    not 0 <= finite(row["normalized_ledger_error"], "ledger") <= \
                    max(1e-11, 5 * levels * math.ulp(1.0)):
                raise ArtifactError("ineligible physical numerical state admitted")
        elif any(row[name] for name in error_names):
            raise ArtifactError("failed/terminal physical row has fabricated errors")
        if row["status"] != "eligible" and row["status"] != "selected" and \
                not row["failure_reason"]:
            raise ArtifactError("physical failure lacks reason")
        if kind == "terminal" and any(row[name] for name in PHYSICAL_NUMERICAL):
            raise ArtifactError("physical terminal has candidate diagnostics")
        for name in PHYSICAL_NUMERICAL[:2]:
            if row[name] and integer(row[name], name) < 0:
                raise ArtifactError("negative diagnostic count")
        for name in PHYSICAL_NUMERICAL[2:4]:
            if row[name] and finite(row[name], name) < 0:
                raise ArtifactError("negative diagnostic metric")
        if row["state_checksum"] and not 0 <= integer(row["state_checksum"], "checksum") < 2**64:
            raise ArtifactError("physical checksum outside uint64")
        indexed[key] = row
    for fixture in PHYSICAL_FIXTURES:
        baselines = [r for (kind, f, _, _), r in indexed.items()
                     if kind == "baseline" and f == fixture]
        if len(baselines) != 1:
            raise ArtifactError("missing physical controller baseline")
        count = int(baselines[0]["production_nsubmix"])
        for method in METHODS:
            terminals = [r for (kind, f, m, _), r in indexed.items()
                         if kind == "terminal" and (f, m) == (fixture, method)]
            if len(terminals) != 1:
                raise ArtifactError("missing physical search outcome")
            terminal = terminals[0]
            points = sorted((n, r) for (kind, f, m, n), r in indexed.items()
                            if kind == "measurement" and (f, m) == (fixture, method))
            start = count if method == "forward_euler" else 1
            last = points[-1][0] if points else 0
            if [n for n, _ in points] != list(range(start, last + 1)) or \
                    int(terminal["steps"]) != last:
                raise ArtifactError("physical integer search incomplete")
            tight = [n for n, r in points if r["status"] == "eligible" and
                     float(r["error_max"]) <= .001]
            if (terminal["status"] == "selected") != bool(tight) or tight and tight[0] != last:
                raise ArtifactError("physical search did not stop at first tight match")
            if terminal["status"] == "not_attained_within_cap" and \
                    last != (8192 if start <= 8192 else 0):
                raise ArtifactError("physical cap outcome lacks complete search")
            if method == "forward_euler" and points:
                baseline = baselines[0]
                first = points[0][1]
                for name in ("status", "failure_reason", *error_names, *PHYSICAL_NUMERICAL):
                    if baseline[name] != first[name]:
                        raise ArtifactError("physical FE baseline differs from controller measurement")
            selected_counts = {min(n for n, r in points if r["status"] == "eligible" and
                                   float(r["error_max"]) <= target)
                               for target in PHYSICAL_TARGETS
                               if any(r["status"] == "eligible" and float(r["error_max"]) <= target
                                      for _, r in points)}
            permitted = {n * mult for n in selected_counts for mult in (2, 4)}
            for (kind, f, m, n), row in indexed.items():
                if kind != "confirmation" or (f, m) != (fixture, method):
                    continue
                if n not in permitted:
                    raise ArtifactError("unselected physical confirmation count")
                measured = indexed.get(("measurement", f, m, n))
                if measured and row["status"] not in {"accuracy_budget_exhausted",
                                                       "confirmation_count_exceeds_cap"}:
                    for name in ("status", "failure_reason", *error_names, *PHYSICAL_NUMERICAL):
                        if measured[name] != row[name]:
                            raise ArtifactError("physical confirmation disagrees with measured state")
    return indexed


def physical_selection(indexed, identity_sha, accuracy_sha):
    chosen, baselines, configs = [], [], set()
    for fixture in PHYSICAL_FIXTURES:
        baseline = next(r for (kind, f, _, _), r in indexed.items()
                        if kind == "baseline" and f == fixture)
        count = int(baseline["production_nsubmix"])
        baselines.append({"fixture": fixture, "steps": count,
                          "status": baseline["status"],
                          "achieved_error": float(baseline["error_max"])
                          if baseline["status"] == "eligible" else None})
        for method in METHODS:
            for target in PHYSICAL_TARGETS:
                eligible = [(n, r) for (kind, f, m, n), r in indexed.items()
                            if kind == "measurement" and (f, m) == (fixture, method)
                            and r["status"] == "eligible" and float(r["error_max"]) <= target]
                n, coarse = min(eligible) if eligible else (None, None)
                reason = ""
                if coarse:
                    confirmation = [indexed.get(("confirmation", fixture, method, n * mult))
                                    for mult in (2, 4)]
                    if any(r is None for r in confirmation):
                        reason = "missing_confirmation"
                    elif any(r["status"] != "eligible" for r in confirmation):
                        reason = "ineligible_confirmation"
                    else:
                        evidence = [coarse, *confirmation]
                        for name in PHYSICAL_ERROR_VECTORS:
                            maxima = [max(map(float, r[name].split(";"))) for r in evidence]
                            if not maxima[0] > maxima[1] > maxima[2]:
                                reason = "nondecreasing_confirmation"
                else:
                    terminal = next(r for (kind, f, m, _), r in indexed.items()
                                    if kind == "terminal" and (f, m) == (fixture, method))
                    reason = terminal["status"]
                admitted = coarse is not None and not reason
                chosen.append({"fixture": fixture, "method": method, "target": target,
                               "steps": n, "production_nsubmix": count,
                               "fe_work_role": ("controller_selected" if n == count
                                                else "accuracy_refined")
                               if method == "forward_euler" and n is not None else None,
                               "achieved_error": float(coarse["error_max"]) if coarse else None,
                               "status": "selected" if coarse else "unavailable",
                               "timing_admitted": admitted, "exclusion_reason": reason})
                if admitted:
                    configs.add((fixture, method, n))
    return {"artifact_schema": "ndrop_physical_column_selection_v1",
            "identity_sha256": identity_sha, "accuracy_sha256": accuracy_sha,
            "targets": list(PHYSICAL_TARGETS), "baselines": baselines,
            "selected": chosen, "unique_configs": [list(key) for key in sorted(configs)]}


class ArtifactError(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_new(path: Path, content: bytes) -> str:
    if path.exists() or path.with_suffix(path.suffix + ".partial").exists():
        raise ArtifactError(f"output already exists: {path}")
    temporary = path.with_suffix(path.suffix + ".partial")
    try:
        with temporary.open("xb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        # A failed partial remains visibly inadmissible for later inspection.
        raise
    return sha256(path)


def canonical_json(data: Any) -> bytes:
    return (json.dumps(data, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("utf-8")


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def git_checkout(path: Path) -> dict[str, Any]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=path,
                          check=True, capture_output=True, text=True).stdout.strip()
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
    if dependencies["mam4xx"]["head"] != \
            "0d076d8ebb004a18f05f51f97b08edd9def6826f":
        raise ArtifactError("the MAM4xx source is not the locked current target")
    if any(record["dirty"] for record in dependencies.values()):
        raise ArtifactError("a locked dependency checkout is dirty")
    return {"cmake_cache_sha256": sha256(cache),
            "compiler": entries["CMAKE_CXX_COMPILER"],
            "build_type": entries["CMAKE_BUILD_TYPE"],
            "sanitizers": entries["NDROP_ENABLE_SANITIZERS"],
            "kokkos_dir": entries["Kokkos_DIR"],
            "dependencies": dependencies}


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
    except ValueError as exc:
        raise ArtifactError(f"invalid {label}") from exc
    if not math.isfinite(value):
        raise ArtifactError(f"nonfinite {label}")
    return value


def integer(text: str, label: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise ArtifactError(f"invalid {label}") from exc
    if str(value) != text:
        raise ArtifactError(f"noncanonical {label}")
    return value


def validate_accuracy(rows: list[dict[str, str]]) -> dict[tuple[int, str, int], dict[str, str]]:
    indexed: dict[tuple[int, str, int], dict[str, str]] = {}
    fixture_scales: dict[int, tuple[float, float]] = {}
    expected = set()
    for fixture in FIXTURES:
        for method in METHODS:
            for power in range(14):
                steps = 1 << power
                if method == "forward_euler" and steps < (6 if fixture % 2 == 0 else 176):
                    continue
                expected.add((fixture, method, steps))
    for row in rows:
        if row["artifact_schema"] != "ndrop_matched_accuracy_v1":
            raise ArtifactError("old or mixed accuracy schema")
        fixture = integer(row["fixture"], "fixture")
        steps = integer(row["steps"], "steps")
        method = row["method"]
        key = (fixture, method, steps)
        if key not in expected or key in indexed:
            raise ArtifactError("unknown or duplicate accuracy key")
        levels, mixing, interval = ((72, .05, 50), (72, .5, 150),
                                    (128, .05, 50), (128, .5, 150))[fixture]
        if (integer(row["levels"], "levels") != levels or
                finite(row["mixing_scale"], "mixing") != mixing or
                finite(row["interval_seconds"], "interval") != interval):
            raise ArtifactError("fixture identity mismatch")
        scales = (finite(row["mixing_rate_t"], "mixing rate times interval"),
                  finite(row["activation_rate_t"],
                         "activation rate times interval"))
        if any(value < 0 for value in scales) or \
                scales != fixture_scales.setdefault(fixture, scales):
            raise ArtifactError("fixture operator scales are inconsistent")
        if row["reference_status"] not in {"certified", "unresolved"}:
            raise ArtifactError("reference status is invalid")
        if row["reference_status"] == "certified":
            if row["reference_reason"]:
                raise ArtifactError("certified reference has a failure reason")
            disagreement = finite(row["reference_disagreement"],
                                  "reference disagreement")
            if not 0 <= disagreement <= 2e-11:
                raise ArtifactError("reference disagreement exceeds gate")
        elif row["status"] != "reference_unresolved":
            raise ArtifactError("unresolved reference entered candidate analysis")
        elif not row["reference_reason"]:
            raise ArtifactError("unresolved reference lacks a failure reason")
        if row["status"] not in {"eligible", "reference_unresolved",
                                 "numerical_failure", "fe_clipping",
                                 "material_negative"}:
            raise ArtifactError("unknown candidate status")
        if row["status"] == "eligible":
            if row["reference_status"] != "certified":
                raise ArtifactError("eligible without reference")
            if row["failure_reason"]:
                raise ArtifactError("eligible candidate has a failure reason")
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
            if not 0 <= residual <= 1e-11 or not 0 <= ledger <= ledger_limit or \
                    integer(row["clipping_count"], "clipping count") != 0 or \
                    integer(row["material_negative_count"], "negative count") != 0:
                raise ArtifactError("ineligible numerical state admitted")
        elif any(row[name] for name in ("error_normal", "error_tight", "error_max",
                                        "normal_weighted_l1_by_field",
                                        "normal_linf_by_field",
                                        "tight_weighted_l1_by_field",
                                        "tight_linf_by_field")):
            raise ArtifactError("failed candidate has fabricated error")
        elif not row["failure_reason"]:
            raise ArtifactError("failed candidate lacks its first reason")
        if row["richardson_to_next"] and \
                finite(row["richardson_to_next"], "Richardson difference") < 0:
            raise ArtifactError("negative Richardson difference")
        if row["observed_order"]:
            finite(row["observed_order"], "observed order")
        checksum = integer(row["state_checksum"], "state checksum")
        if not 0 <= checksum < 2**64:
            raise ArtifactError("state checksum is outside uint64")
        indexed[key] = row
    if set(indexed) != expected:
        raise ArtifactError("accuracy matrix is incomplete")
    for (fixture, method, steps), row in indexed.items():
        next_row = indexed.get((fixture, method, steps * 2))
        pair_available = row["status"] == "eligible" and \
            next_row is not None and next_row["status"] == "eligible"
        if bool(row["richardson_to_next"]) != pair_available:
            raise ArtifactError("Richardson pair applicability mismatch")
        next_difference = next_row["richardson_to_next"] if next_row else ""
        order_available = bool(row["richardson_to_next"] and next_difference and
                               float(row["richardson_to_next"]) > 0 and
                               float(next_difference) > 0)
        if bool(row["observed_order"]) != order_available:
            raise ArtifactError("observed order applicability mismatch")
    return indexed


def validate_physical_identity(directory):
    identity = load_json(directory / "identity.json")
    required = {"artifact_schema", "fixture_policy", "fixtures", "targets",
                "source_head", "source_root", "source_status", "source_attached_head",
                "source_files_sha256", "binary_sha256", "build", "design_sha256",
                "python_version", "platform", "grid_files_sha256", "accuracy_command",
                "accuracy_protocol", "timing_protocol"}
    if set(identity) != required or identity["artifact_schema"] != \
            "ndrop_physical_column_identity_v1" or identity["fixture_policy"] != \
            "dycoms_inspired_isothermal_repository_grid_v1" or \
            identity["fixtures"] != list(PHYSICAL_FIXTURES) or \
            identity["targets"] != list(PHYSICAL_TARGETS):
        raise ArtifactError("physical identity schema/fixture policy mismatch")
    root = Path(identity["source_root"])
    if not root.is_absolute() or set(identity["source_files_sha256"]) != set(SOURCES):
        raise ArtifactError("physical source identity incomplete")
    for name, digest in identity["source_files_sha256"].items():
        if sha256(root / "docs/analysis/evidence/ndrop_integrator" / name) != digest:
            raise ArtifactError("physical source changed since preregistered accuracy")
    if sha256(root / "docs/design/eamxx_mam4xx_ndrop.md") != identity["design_sha256"]:
        raise ArtifactError("physical design changed since accuracy")
    if set(identity["grid_files_sha256"]) != set(PHYSICAL_GRIDS):
        raise ArtifactError("physical grid set incomplete")
    for name, digest in identity["grid_files_sha256"].items():
        if sha256(directory / "inputs/grids" / name) != digest or \
                sha256(root / "components/homme/test/vcoord" / name) != digest:
            raise ArtifactError("physical grid copy/source hash mismatch")
    if sha256(directory / "inputs/benchmark_ndrop_matched_error") != identity["binary_sha256"]:
        raise ArtifactError("physical executable hash mismatch")
    if identity["python_version"] != platform.python_version() or \
            identity["accuracy_command"] != ["inputs/benchmark_ndrop_matched_error",
                "--physical-column-accuracy", "--physical-grid-dir", "inputs/grids",
                "--max-steps", "8192", "--budget-seconds", "7200"] or \
            identity["accuracy_protocol"] != {"integer_search": True,
                "maximum_steps": 8192, "budget_seconds": 7200,
                "confirmation_multipliers": [2, 4]} or \
            identity["timing_protocol"] != {"warmups": 10, "minimum_blocks": 30,
                "maximum_blocks": 100, "budget_seconds": 7200,
                "complete_interleaved_blocks": True}:
        raise ArtifactError("physical protocol/version changed")
    build = identity["build"]
    if build["dependencies"]["ekat"]["head"] != "3a1485aecfed7a036e4fa30ce845498fd37975e2" or \
            build["dependencies"]["mam4xx"]["head"] != "0d076d8ebb004a18f05f51f97b08edd9def6826f" or \
            any(item["dirty"] for item in build["dependencies"].values()):
        raise ArtifactError("physical locked dependencies mismatch")
    if identity["source_status"] not in {"exploratory_source", "clean_detached_source"} or \
            identity["source_status"] == "clean_detached_source" and identity["source_attached_head"]:
        raise ArtifactError("physical source status invalid")
    return identity


def run_physical_accuracy(args):
    directory, root, binary = args.output_dir.resolve(), args.source_root.resolve(), args.binary.resolve()
    if directory.exists():
        raise ArtifactError("physical study directory must be new")
    build = build_identity(binary)
    source_dir = root / "docs/analysis/evidence/ndrop_integrator"
    source_hashes = {name: sha256(source_dir / name) for name in SOURCES}
    status = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                            check=True, capture_output=True, text=True).stdout
    attached = subprocess.run(["git", "symbolic-ref", "-q", "HEAD"], cwd=root,
                              capture_output=True).returncode == 0
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, check=True,
                          capture_output=True, text=True).stdout.strip()
    (directory / "inputs/grids").mkdir(parents=True)
    copied = directory / "inputs/benchmark_ndrop_matched_error"
    atomic_new(copied, binary.read_bytes())
    copied.chmod(binary.stat().st_mode & 0o777)
    grids = {}
    for name in PHYSICAL_GRIDS:
        source = root / "components/homme/test/vcoord" / name
        grids[name] = atomic_new(directory / "inputs/grids" / name, source.read_bytes())
        if grids[name] != sha256(source):
            raise ArtifactError("physical grid copy changed")
    command = ["inputs/benchmark_ndrop_matched_error", "--physical-column-accuracy",
               "--physical-grid-dir", "inputs/grids", "--max-steps", "8192",
               "--budget-seconds", "7200"]
    identity = {"artifact_schema": "ndrop_physical_column_identity_v1",
                "fixture_policy": "dycoms_inspired_isothermal_repository_grid_v1",
                "fixtures": list(PHYSICAL_FIXTURES), "targets": list(PHYSICAL_TARGETS),
                "source_head": head, "source_root": str(root),
                "source_status": "clean_detached_source" if not status and not attached
                                 else "exploratory_source", "source_attached_head": attached,
                "source_files_sha256": source_hashes, "binary_sha256": sha256(copied),
                "build": build, "grid_files_sha256": grids,
                "design_sha256": sha256(root / "docs/design/eamxx_mam4xx_ndrop.md"),
                "python_version": platform.python_version(), "platform": platform.platform(),
                "accuracy_command": command, "accuracy_protocol": {
                    "integer_search": True, "maximum_steps": 8192,
                    "budget_seconds": 7200, "confirmation_multipliers": [2, 4]},
                "timing_protocol": {"warmups": 10, "minimum_blocks": 30,
                    "maximum_blocks": 100, "budget_seconds": 7200,
                    "complete_interleaved_blocks": True}}
    atomic_new(directory / "identity.json", canonical_json(identity))
    validate_physical_identity(directory)
    process = subprocess.run([str(copied), *command[1:]], cwd=directory,
                             capture_output=True, check=True)
    output = directory / "accuracy.csv"
    atomic_new(output, process.stdout)
    validate_physical_accuracy(csv_rows(output, PHYSICAL_ACCURACY_HEADER))
    completion = {"artifact_schema": "ndrop_physical_column_completion_v1",
                  "artifact_sha256": sha256(output), "path": output.name,
                  "identity_sha256": sha256(directory / "identity.json")}
    atomic_new(directory / "accuracy.complete.json", canonical_json(completion))


def accuracy_for_study(directory, physical):
    path = completed(directory, "accuracy")
    if physical:
        indexed = validate_physical_accuracy(csv_rows(path, PHYSICAL_ACCURACY_HEADER))
        return {(f, m, n): r for (kind, f, m, n), r in indexed.items() if kind == "measurement"}
    return validate_accuracy(csv_rows(path, ACCURACY_HEADER))


def validate_accuracy_study(args):
    identity = validate_identity(args.output_dir)
    accuracy_for_study(args.output_dir, identity["artifact_schema"] == "ndrop_physical_column_identity_v1")
    selection(args.output_dir)


def validate_identity(directory: Path) -> dict[str, Any]:
    identity_path = directory / "identity.json"
    identity = load_json(identity_path)
    if identity.get("artifact_schema") == "ndrop_physical_column_identity_v1":
        return validate_physical_identity(directory)
    if identity.get("artifact_schema") != "ndrop_matched_identity_v1":
        raise ArtifactError("identity schema mismatch")
    if sha256(directory / "inputs" / "benchmark_ndrop_matched_error") != identity["binary_sha256"]:
        raise ArtifactError("timing binary differs from identity")
    return identity


def completed(directory: Path, stem: str) -> Path:
    identity_path = directory / "identity.json"
    record = load_json(directory / f"{stem}.complete.json")
    physical = load_json(identity_path).get("artifact_schema") == "ndrop_physical_column_identity_v1"
    schema = "ndrop_physical_column_completion_v1" if physical else "ndrop_matched_completion_v1"
    if record != {"artifact_schema": schema,
                  "artifact_sha256": sha256(directory / f"{stem}.csv"),
                  "identity_sha256": sha256(identity_path),
                  "path": f"{stem}.csv"}:
        raise ArtifactError(f"{stem} completion is missing or stale")
    return directory / f"{stem}.csv"


def run_accuracy(args: argparse.Namespace) -> None:
    if getattr(args, "physical", False):
        run_physical_accuracy(args)
        return
    directory = args.output_dir.resolve()
    if directory.exists():
        raise ArtifactError("study directory must be new")
    root = args.source_root.resolve()
    binary = args.binary.resolve()
    if not binary.is_file() or not root.is_dir():
        raise ArtifactError("source or binary is absent")
    source_dir = root / "docs/analysis/evidence/ndrop_integrator"
    source_hashes = {name: sha256(source_dir / name) for name in SOURCES}
    build = build_identity(binary)
    status = subprocess.run(["git", "status", "--porcelain"], cwd=root,
                            check=True, capture_output=True, text=True).stdout
    attached = subprocess.run(["git", "symbolic-ref", "-q", "HEAD"],
                              cwd=root, capture_output=True).returncode == 0
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                          check=True, capture_output=True, text=True).stdout.strip()
    directory.mkdir(parents=True)
    (directory / "inputs").mkdir()
    retained_binary = directory / "inputs" / "benchmark_ndrop_matched_error"
    with binary.open("rb") as source, retained_binary.open("xb") as target:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            target.write(chunk)
        target.flush()
        os.fsync(target.fileno())
    retained_binary.chmod(binary.stat().st_mode & 0o777)
    if sha256(retained_binary) != sha256(binary):
        raise ArtifactError("binary copy changed")
    identity = {
        "artifact_schema": "ndrop_matched_identity_v1",
        "binary_sha256": sha256(retained_binary),
        "source_head": head,
        "source_status": "clean_detached_source" if not status and not attached
                         else "exploratory_source",
        "source_attached_head": attached,
        "source_files_sha256": source_hashes,
        "build": build,
        "design_sha256": sha256(root / "docs/design/eamxx_mam4xx_ndrop.md"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "targets": TARGETS,
        "accuracy_command": ["inputs/benchmark_ndrop_matched_error",
                             "--accuracy"],
        "timing_protocol": {"warmups": 10, "minimum_blocks": 30,
                            "maximum_blocks": 100,
                            "budget_seconds": 7200,
                            "complete_interleaved_blocks": True},
    }
    atomic_new(directory / "identity.json", canonical_json(identity))
    process = subprocess.run([str(retained_binary), "--accuracy"],
                             capture_output=True, check=True)
    output = directory / "accuracy.csv"
    atomic_new(output, process.stdout)
    validate_accuracy(csv_rows(output, ACCURACY_HEADER))
    record = {"artifact_schema": "ndrop_matched_completion_v1",
              "artifact_sha256": sha256(output),
              "identity_sha256": sha256(directory / "identity.json"),
              "path": output.name}
    atomic_new(directory / "accuracy.complete.json", canonical_json(record))


def expected_selection(indexed: dict[tuple[int, str, int], dict[str, str]],
                       identity_sha: str, accuracy_sha: str) -> dict[str, Any]:
    selected = []
    unique = set()
    for fixture in FIXTURES:
        for method in METHODS:
            for target in TARGETS:
                qualifying = [steps for (f, m, steps), row in indexed.items()
                              if f == fixture and m == method and
                              row["status"] == "eligible" and
                              float(row["error_max"]) <= target]
                steps = min(qualifying) if qualifying else None
                achieved = (float(indexed[(fixture, method, steps)]["error_max"])
                            if steps is not None else None)
                selected.append({"fixture": fixture, "method": method,
                                 "target": target, "steps": steps,
                                 "achieved_error": achieved,
                                 "target_to_achieved_error_ratio":
                                     target / achieved if achieved else None,
                                 "status": "selected" if steps is not None
                                 else "not_attained_within_cap"})
                if steps is not None:
                    unique.add((fixture, method, steps))
    return {"artifact_schema": "ndrop_matched_selection_v1",
              "identity_sha256": identity_sha,
              "accuracy_sha256": accuracy_sha,
              "targets": TARGETS, "selected": selected,
              "unique_configs": [list(config) for config in sorted(unique)]}


def freeze(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    identity = validate_identity(directory)
    accuracy_path = completed(directory, "accuracy")
    physical = identity["artifact_schema"] == "ndrop_physical_column_identity_v1"
    indexed = (validate_physical_accuracy(csv_rows(accuracy_path, PHYSICAL_ACCURACY_HEADER))
               if physical else validate_accuracy(csv_rows(accuracy_path, ACCURACY_HEADER)))
    frozen = (physical_selection if physical else expected_selection)(
        indexed, sha256(directory / "identity.json"), sha256(accuracy_path))
    atomic_new(directory / "selection.json", canonical_json(frozen))


def selection(directory: Path) -> dict[str, Any]:
    frozen = load_json(directory / "selection.json")
    accuracy_path = completed(directory, "accuracy")
    physical = load_json(directory / "identity.json").get("artifact_schema") == "ndrop_physical_column_identity_v1"
    indexed = (validate_physical_accuracy(csv_rows(accuracy_path, PHYSICAL_ACCURACY_HEADER))
               if physical else validate_accuracy(csv_rows(accuracy_path, ACCURACY_HEADER)))
    expected = (physical_selection if physical else expected_selection)(
        indexed, sha256(directory / "identity.json"), sha256(accuracy_path))
    if frozen != json.loads(canonical_json(expected)):
        raise ArtifactError("selection is stale or malformed")
    return frozen


def validate_timing(rows: list[dict[str, str]], frozen: dict[str, Any],
                    accuracy: dict[tuple[int, str, int], dict[str, str]],
                    physical: bool = False) -> None:
    expected = {tuple(config) for config in frozen["unique_configs"]}
    groups: dict[tuple[int, str, int], list[int]] = {}
    checksums: dict[tuple[int, str, int], set[int]] = {}
    reasons: dict[int, set[str]] = {}
    for row in rows:
        if row["artifact_schema"] != ("ndrop_physical_column_timing_v1"
                                      if physical else "ndrop_matched_timing_v1"):
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
    if set(groups) != expected:
        raise ArtifactError("missing selected timed configuration")
    for key, blocks in groups.items():
        if blocks != list(range(len(blocks))) or len(checksums[key]) != 1 or \
                not 30 <= len(blocks) <= 100:
            raise ArtifactError("incomplete or inconsistent timing blocks")
    for fixture in (PHYSICAL_FIXTURES if physical else FIXTURES):
        group_counts = {len(blocks) for key, blocks in groups.items()
                        if key[0] == fixture}
        if len(group_counts) > 1 or len(reasons.get(fixture, set())) > 1 or \
                not reasons.get(fixture, {"maximum_blocks"}) <= {
                    "median_precision", "maximum_blocks", "budget_exhausted"}:
            raise ArtifactError("timing methods are not paired by complete blocks")


def run_timing(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    identity = validate_identity(directory)
    physical = identity["artifact_schema"] == "ndrop_physical_column_identity_v1"
    frozen = selection(directory)
    accuracy = accuracy_for_study(directory, physical)
    # A dirty or attached source produces exploratory evidence, bound to the
    # exact copied binary and file hashes; the analysis records that status.
    if identity["source_status"] not in {"clean_detached_source",
                                     "exploratory_source"}:
        raise ArtifactError("source status is invalid")
    if not frozen["unique_configs"]:
        raise ArtifactError("no eligible configuration to time")
    encoded = ";".join(f"{f}:{m}:{n}" for f, m, n in frozen["unique_configs"])
    command = [str(directory / "inputs" / "benchmark_ndrop_matched_error"),
         "--physical-column-timing" if physical else "--timing", "--selection", encoded, "--warmups", "10",
         "--min-blocks", "30", "--max-blocks", "100",
         "--budget-seconds", "7200"]
    if physical:
        command.extend(["--physical-grid-dir", str(directory / "inputs/grids")])
    process = subprocess.run(command, capture_output=True, check=True)
    output = directory / "timing.csv"
    atomic_new(output, process.stdout)
    validate_timing(csv_rows(output, TIMING_HEADER), frozen, accuracy, physical)
    record = {"artifact_schema": "ndrop_physical_column_completion_v1" if physical
              else "ndrop_matched_completion_v1",
              "artifact_sha256": sha256(output),
              "identity_sha256": sha256(directory / "identity.json"),
              "path": output.name}
    atomic_new(directory / "timing.complete.json", canonical_json(record))


def median(values: list[float]) -> float:
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
    physical = identity["artifact_schema"] == "ndrop_physical_column_identity_v1"
    frozen = selection(directory)
    accuracy = accuracy_for_study(directory, physical)
    rows = csv_rows(completed(directory, "timing"), TIMING_HEADER)
    validate_timing(rows, frozen, accuracy, physical)
    by_key: dict[tuple[int, str, int], list[float]] = {}
    stop_reasons: dict[int, str] = {}
    for row in rows:
        key = (int(row["fixture"]), row["method"], int(row["steps"]))
        by_key.setdefault(key, []).append(float(row["seconds_per_column"]))
        stop_reasons[key[0]] = row["stop_reason"]
    configurations = [
        {"fixture": fixture, "method": method, "steps": steps,
         "achieved_error": float(accuracy[(fixture, method, steps)]["error_max"]),
         "median_seconds_per_column": median(samples),
         "timing_blocks": len(samples), "stop_reason": stop_reasons[fixture]}
        for (fixture, method, steps), samples in sorted(by_key.items())
    ]
    results = []
    for fixture in (PHYSICAL_FIXTURES if physical else FIXTURES):
        for target in (PHYSICAL_TARGETS if physical else TARGETS):
            chosen = {entry["method"]: entry["steps"] if not physical or entry["timing_admitted"]
                      else None for entry in frozen["selected"]
                      if entry["fixture"] == fixture and entry["target"] == target}
            for numerator, denominator in (("forward_euler", "imex"),
                                           ("forward_euler", "backward_euler"),
                                           ("imex", "backward_euler")):
                left = chosen[numerator]
                right = chosen[denominator]
                if left is None or right is None:
                    results.append({"fixture": fixture, "target": target,
                                    "numerator": numerator, "denominator": denominator,
                                    "status": "missing_eligible_pair"})
                    continue
                ratio = paired_ratio(by_key[(fixture, numerator, left)],
                                     by_key[(fixture, denominator, right)],
                                     0x4e44524f)
                results.append({"fixture": fixture, "target": target,
                                "numerator": numerator, "denominator": denominator,
                                "numerator_steps": left,
                                "denominator_steps": right,
                                "numerator_achieved_error": float(
                                    accuracy[(fixture, numerator, left)]["error_max"]),
                                "denominator_achieved_error": float(
                                    accuracy[(fixture, denominator, right)]["error_max"]),
                                "status": "inconclusive" if ratio["lower"] <= 1 <= ratio["upper"]
                                else "interval_excludes_one", **ratio})
    analysis = {"artifact_schema": "ndrop_physical_column_analysis_v1" if physical
                else "ndrop_matched_analysis_v1",
                "identity_sha256": sha256(directory / "identity.json"),
                "accuracy_sha256": sha256(directory / "accuracy.csv"),
                "selection_sha256": sha256(directory / "selection.json"),
                "timing_sha256": sha256(directory / "timing.csv"),
                "source_status": identity["source_status"],
                "bootstrap": {"resamples": 10000, "seed": 0x4e44524f,
                              "rng": "Python random.Random integer seed",
                              "median": "lower+(upper-lower)/2",
                              "endpoints": "zero_based_sorted_replicates_249_9749"},
                "configurations": configurations,
                "comparisons": results}
    if physical:
        analysis["controller_baselines"] = frozen["baselines"]
        analysis["scope"] = "common-ceiling standalone idealization; no production fidelity or climatology"
        for item in configurations:
            row = accuracy[(item["fixture"], item["method"], item["steps"])]
            item["production_nsubmix"] = int(row["production_nsubmix"])
            item["fe_work_role"] = ("controller_selected" if item["steps"] ==
                int(row["production_nsubmix"]) else "accuracy_refined") \
                if item["method"] == "forward_euler" else None
    atomic_new(directory / "analysis.json", canonical_json(analysis))
    files = ("identity.json", "accuracy.csv", "accuracy.complete.json",
             "selection.json", "timing.csv", "timing.complete.json",
             "analysis.json", "inputs/benchmark_ndrop_matched_error")
    if physical:
        files += tuple("inputs/grids/" + name for name in PHYSICAL_GRIDS)
    final = {"artifact_schema": "ndrop_physical_column_final_v1" if physical else "ndrop_matched_final_v1",
             "source_status": identity["source_status"],
             "files_sha256": {name: sha256(directory / name) for name in files},
             "fixtures": list(PHYSICAL_FIXTURES if physical else FIXTURES),
             "targets": list(PHYSICAL_TARGETS if physical else TARGETS),
             "completed_timing_rows": len(rows)}
    atomic_new(directory / "final_manifest.json", canonical_json(final))


def validate_final(args: argparse.Namespace) -> None:
    directory = args.output_dir.resolve()
    identity = validate_identity(directory)
    physical = identity["artifact_schema"] == "ndrop_physical_column_identity_v1"
    frozen = selection(directory)
    accuracy = accuracy_for_study(directory, physical)
    rows = csv_rows(completed(directory, "timing"), TIMING_HEADER)
    validate_timing(rows, frozen, accuracy, physical)
    final = load_json(directory / "final_manifest.json")
    if final.get("artifact_schema") != ("ndrop_physical_column_final_v1" if physical else "ndrop_matched_final_v1") or \
            final.get("source_status") != identity["source_status"] or \
            final.get("completed_timing_rows") != len(rows) or \
            final.get("fixtures") != list(PHYSICAL_FIXTURES if physical else FIXTURES) or \
            final.get("targets") != list(PHYSICAL_TARGETS if physical else TARGETS):
        raise ArtifactError("final manifest coverage or identity mismatch")
    for name, expected in final["files_sha256"].items():
        if name.startswith("/") or ".." in Path(name).parts or \
                sha256(directory / name) != expected:
            raise ArtifactError("final manifest file hash mismatch")
    required = {"identity.json", "accuracy.csv", "accuracy.complete.json",
                "selection.json", "timing.csv", "timing.complete.json",
                "analysis.json", "inputs/benchmark_ndrop_matched_error"}
    if physical:
        required.update("inputs/grids/" + name for name in PHYSICAL_GRIDS)
    if set(final["files_sha256"]) != required:
        raise ArtifactError("final manifest file set is incomplete")


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run-accuracy", "freeze", "validate-accuracy", "run-timing", "report",
                 "validate-final"):
        command = commands.add_parser(name)
        command.add_argument("--output-dir", type=Path, required=True)
        if name == "run-accuracy":
            command.add_argument("--physical", action="store_true",
                                 help="separate physically scaled common-ceiling study")
            command.add_argument("--binary", type=Path, required=True)
            command.add_argument("--source-root", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        {"run-accuracy": run_accuracy, "freeze": freeze,
         "validate-accuracy": validate_accuracy_study,
         "run-timing": run_timing, "report": report,
         "validate-final": validate_final}[arguments.command](arguments)
    except (ArtifactError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
