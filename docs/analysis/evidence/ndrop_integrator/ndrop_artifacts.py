#!/usr/bin/env python3
"""Fail-closed Task 5R artifact transactions and provenance validation.

This module keeps instrumented phase evidence paired with, but separate from,
the primary uninstrumented timing summary. Retained runs require a clean
detached source, immutable identity metadata, aligned shard transactions, and
an authoritative final manifest written only after every hash and row-level
contract has been revalidated.
"""

import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from analyze_screening import (
    ScreeningInputError,
    ingest_paired_artifacts,
    ingest_rows,
)


SUMMARY_COLUMNS = (
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
    "forced_fe_substeps",
    "linear_phase_mode",
    "validation_ok",
    "eligible",
    "exclusion_reason",
    "message",
    "repetitions",
    "richardson_ok",
    "richardson_message",
    "screening_error",
    "richardson_maximum_weighted_l1",
    "richardson_maximum_linf",
    "richardson_weighted_l1_by_field",
    "richardson_linf_by_field",
    "refined_checksum",
    "stop_reason",
    "median_seconds_per_column",
    "ci95_lower",
    "ci95_upper",
    "ci95_relative_half_width",
    "instrumented_median_seconds_per_column",
    "instrumented_ci95_lower",
    "instrumented_ci95_upper",
    "instrumentation_overhead",
    "setup_seconds_per_column",
    "factorization_seconds_per_column",
    "solve_seconds_per_column",
    "combined_factorization_solve_seconds_per_column",
    "update_seconds_per_column",
    "diagnostics_seconds_per_column",
    "fence_seconds_per_column",
    "clipping_count",
    "material_negative_count",
    "block_factorizations",
    "block_rhs_solves",
    "scalar_factorizations",
    "scalar_rhs_solves",
    "executed_field_level_updates",
    "workspace_bytes",
    "checksum",
    "compiler",
    "platform",
    "kokkos_execution_space",
    "hardware_threads",
    "warmups",
    "minimum_repetitions",
    "maximum_repetitions",
    "bootstrap_resamples",
    "raw_samples_seconds_per_column",
    "instrumented_samples_seconds_per_column",
)

PHASE_COLUMNS = (
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
)

IDENTITY_SCHEMA = "ndrop_identity_manifest_v1"
COMPLETION_SCHEMA = "ndrop_shard_completion_v1"
FINAL_SCHEMA = "ndrop_final_manifest_v1"
SUMMARY_SCHEMA = "ndrop_summary_v2"
PHASE_SCHEMA = "ndrop_phase_samples_v1"
MINIMUM_FREE_BYTES = 5 * 1024 * 1024 * 1024

_RUN_ID = re.compile(r"^[A-Za-z0-9._-]{1,128}$")
_RETAINED_RUN_ID = re.compile(
    r"^task5-remediated-[0-9]{8}T[0-9]{6}Z-([0-9a-f]{12})$"
)
_SOURCE_COMMIT = re.compile(r"^[0-9a-f]{40}$")

_IDENTITY_FIELDS = frozenset(
    {
        "artifact_schema",
        "created_utc",
        "validation_mode",
        "retained_contract",
        "run_id",
        "source_root",
        "source_commit",
        "source_clean",
        "source_detached",
        "benchmark_binary_path",
        "benchmark_binary_sha256",
        "dependency_provenance_path",
        "dependency_provenance_sha256",
        "compiler",
        "compiler_flags",
        "machine",
        "platform",
        "execution_space",
        "validation_commands",
        "summary_schema",
        "phase_samples_schema",
        "scratch_root",
        "output_root",
        "identity_manifest_path",
        "merged_summary_partial_path",
        "merged_phase_samples_partial_path",
        "merged_summary_path",
        "merged_phase_samples_path",
        "final_manifest_path",
        "total_case_count",
        "shard_case_count",
        "shard_count",
        "worker_count",
        "warmups",
        "minimum_repetitions",
        "maximum_repetitions",
        "bootstrap_resamples",
        "benchmark_commands",
        "shards",
    }
)

_COMPLETION_FIELDS = frozenset(
    {
        "artifact_schema",
        "created_utc",
        "run_id",
        "identity_manifest_path",
        "identity_manifest_sha256",
        "case_start",
        "case_count",
        "case_end_exclusive",
        "summary_path",
        "summary_sha256",
        "summary_rows",
        "phase_samples_path",
        "phase_samples_sha256",
        "phase_rows",
    }
)

_FINAL_FIELDS = frozenset(
    {
        "artifact_schema",
        "created_utc",
        "run_id",
        "identity_manifest_path",
        "identity_manifest_sha256",
        "completion_records",
        "summary_path",
        "summary_sha256",
        "summary_rows",
        "phase_samples_path",
        "phase_samples_sha256",
        "phase_rows",
    }
)


class ArtifactError(ValueError):
    """Raised when provenance or paired artifact validation fails."""


@dataclass(frozen=True)
class PairValidation:
    run_id: str
    summary_rows: int
    phase_rows: int
    case_ids: Tuple[int, ...]


@dataclass(frozen=True)
class CompletedShard:
    completion_path: str
    summary_path: str
    phase_samples_path: str
    summary_rows: int
    phase_rows: int


@dataclass(frozen=True)
class FinalRun:
    final_manifest_path: str
    summary_path: str
    phase_samples_path: str
    summary_rows: int
    phase_rows: int


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _absolute(path: Path) -> str:
    return str(Path(path).expanduser().resolve())


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest for one regular file."""

    candidate = Path(path)
    if not candidate.is_file():
        raise ArtifactError("required file does not exist: {}".format(candidate))
    digest = hashlib.sha256()
    try:
        with candidate.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise ArtifactError("cannot hash {}: {}".format(candidate, error))
    return digest.hexdigest()


def _require_exact_fields(
    document: Mapping, expected: frozenset, artifact_name: str
) -> None:
    observed = frozenset(document)
    missing = sorted(expected.difference(observed))
    extra = sorted(observed.difference(expected))
    if missing or extra:
        raise ArtifactError(
            "{} fields differ: missing={} extra={}".format(
                artifact_name, ",".join(missing), ",".join(extra)
            )
        )


def _require_absolute_path(value: object, name: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ArtifactError("{} must be a nonempty absolute path".format(name))
    path = Path(value)
    if not path.is_absolute() or str(path.resolve()) != value:
        raise ArtifactError("{} must be a normalized absolute path".format(name))
    return path


def _integer(value: object, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ArtifactError("{} must be an integer".format(name))
    try:
        result = int(value)
    except (TypeError, ValueError):
        raise ArtifactError("{} must be an integer".format(name))
    if str(result) != str(value).strip() or result < minimum:
        raise ArtifactError(
            "{} must be an integer >= {}".format(name, minimum)
        )
    return result


def _validate_header(reader: csv.DictReader, columns: Tuple[str, ...], name: str) -> None:
    if tuple(reader.fieldnames or ()) != columns:
        raise ArtifactError("{} header does not exactly match its schema".format(name))


def _validate_decoded_row(row: Mapping, name: str) -> None:
    if None in row or any(value is None for value in row.values()):
        raise ArtifactError("{} row width does not match its header".format(name))


def _validate_summary_metadata(row: Mapping, identity: Mapping) -> None:
    exact = {
        "run_id": identity["run_id"],
        "compiler": identity["compiler"],
        "platform": identity["platform"],
        "kokkos_execution_space": identity["execution_space"],
        "warmups": str(identity["warmups"]),
        "minimum_repetitions": str(identity["minimum_repetitions"]),
        "maximum_repetitions": str(identity["maximum_repetitions"]),
        "bootstrap_resamples": str(identity["bootstrap_resamples"]),
    }
    for name, expected in exact.items():
        if row[name] != expected:
            raise ArtifactError(
                "summary {} does not match the identity manifest".format(name)
            )
    if identity["retained_contract"] and row["validation_ok"] == "true":
        repetitions = _integer(row["repetitions"], "repetitions")
        if not (
            identity["minimum_repetitions"]
            <= repetitions
            <= identity["maximum_repetitions"]
        ):
            raise ArtifactError(
                "valid retained rows must meet the configured repetition range"
            )


def expected_retained_case(case_id: int) -> Dict[str, object]:
    """Return the exact approved matrix identity for one retained case ID."""

    if case_id < 0 or case_id >= 24300:
        raise ArtifactError("retained case_id is outside [0, 24300)")
    levels = (1, 2, 3, 16, 32, 64, 72, 127, 128)
    occupancies = (1, 4, 8, 16, 25)
    batches = (1, 16, 64, 256, 1024)
    mixing_scales = (0.05, 0.5)
    intervals = (33.33, 50.0, 75.0, 100.0, 150.0, 300.0)
    fe_substeps = (1, 2, 3, 5, 8, 16, 21)

    group, method_offset = divmod(case_id, 9)
    group, interval_index = divmod(group, len(intervals))
    group, mixing_index = divmod(group, len(mixing_scales))
    group, batch_index = divmod(group, len(batches))
    group, occupancy_index = divmod(group, len(occupancies))
    level_index = group
    if level_index >= len(levels):
        raise ArtifactError("retained case_id cannot be decoded")
    if method_offset < len(fe_substeps):
        method = "forward_euler"
        forced_fe_substeps = fe_substeps[method_offset]
        linear_phase_mode = "none"
    elif method_offset == len(fe_substeps):
        method = "imex"
        forced_fe_substeps = 0
        linear_phase_mode = "combined"
    else:
        method = "backward_euler"
        forced_fe_substeps = 0
        linear_phase_mode = "separate"
    occupancy = occupancies[occupancy_index]
    return {
        "nlev": levels[level_index],
        "nonzero_aerosol_pairs": occupancy,
        "nonzero_state_fields": 1 + 2 * occupancy,
        "executed_aerosol_pairs": 25,
        "executed_transport_fields": 51,
        "batch_size": batches[batch_index],
        "mixing_index": mixing_index,
        "mixing_scale": mixing_scales[mixing_index],
        "dt": intervals[interval_index],
        "method": method,
        "forced_fe_substeps": forced_fe_substeps,
        "linear_phase_mode": linear_phase_mode,
    }


def _validate_retained_case(row: Mapping) -> None:
    case_id = _integer(row["case_id"], "case_id")
    expected = expected_retained_case(case_id)
    integer_fields = (
        "nlev",
        "nonzero_aerosol_pairs",
        "nonzero_state_fields",
        "executed_aerosol_pairs",
        "executed_transport_fields",
        "batch_size",
        "mixing_index",
        "forced_fe_substeps",
    )
    for name in integer_fields:
        if _integer(row[name], name) != expected[name]:
            raise ArtifactError(
                "case {} has an unexpected {}".format(case_id, name)
            )
    for name in ("mixing_scale", "dt"):
        try:
            observed = float(row[name])
        except (TypeError, ValueError):
            raise ArtifactError("{} must be numeric".format(name))
        if observed != expected[name]:
            raise ArtifactError(
                "case {} has an unexpected {}".format(case_id, name)
            )
    for name in ("method", "linear_phase_mode"):
        if row[name] != expected[name]:
            raise ArtifactError(
                "case {} has an unexpected {}".format(case_id, name)
            )


def validate_csv_pair(
    summary_path: Path,
    phase_samples_path: Path,
    expected_run_id: Optional[str] = None,
    expected_case_ids: Optional[Iterable[int]] = None,
    identity: Optional[Mapping] = None,
) -> PairValidation:
    """Validate exact schemas and stream a paired summary/phase artifact.

    Only one case's phase samples are retained in memory at a time. This keeps
    full retained-run validation bounded even at 100 repetitions per case.
    """

    summary_file = Path(summary_path)
    phase_file = Path(phase_samples_path)
    try:
        with summary_file.open(newline="") as stream:
            reader = csv.DictReader(stream, strict=True)
            _validate_header(reader, SUMMARY_COLUMNS, "summary")
            summaries = []
            for row in reader:
                _validate_decoded_row(row, "summary")
                summaries.append(row)
    except (OSError, csv.Error) as error:
        raise ArtifactError("cannot decode summary CSV: {}".format(error))
    if not summaries:
        raise ArtifactError("summary CSV must contain at least one row")

    try:
        ingest_rows(summaries)
    except ScreeningInputError as error:
        raise ArtifactError("summary validation failed: {}".format(error))
    case_ids = tuple(_integer(row["case_id"], "case_id") for row in summaries)
    if case_ids != tuple(sorted(case_ids)):
        raise ArtifactError("summary case IDs must be strictly increasing")
    if expected_case_ids is not None and case_ids != tuple(expected_case_ids):
        raise ArtifactError("summary case IDs do not match the expected range")
    run_id = str(summaries[0]["run_id"])
    if expected_run_id is not None and run_id != expected_run_id:
        raise ArtifactError("summary run ID does not match the expected run")
    if identity is not None:
        for row in summaries:
            _validate_summary_metadata(row, identity)
            if identity["retained_contract"]:
                _validate_retained_case(row)

    phase_rows = 0
    try:
        with phase_file.open(newline="") as stream:
            reader = csv.DictReader(stream, strict=True)
            _validate_header(reader, PHASE_COLUMNS, "phase sidecar")
            iterator = iter(reader)
            for summary in summaries:
                repetitions = _integer(summary["repetitions"], "repetitions")
                current = []
                for repetition in range(repetitions):
                    try:
                        row = next(iterator)
                    except StopIteration:
                        raise ArtifactError(
                            "phase sidecar ends before summary coverage is complete"
                        )
                    _validate_decoded_row(row, "phase sidecar")
                    if row["case_id"] != summary["case_id"]:
                        raise ArtifactError(
                            "phase sidecar is not ordered by summary case ID"
                        )
                    if _integer(
                        row["instrumented_repetition"],
                        "instrumented_repetition",
                    ) != repetition:
                        raise ArtifactError(
                            "phase repetitions must be zero-based and contiguous"
                        )
                    current.append(row)
                    phase_rows += 1
                try:
                    ingest_paired_artifacts([summary], current)
                except ScreeningInputError as error:
                    raise ArtifactError(
                        "paired phase validation failed: {}".format(error)
                    )
            try:
                extra = next(iterator)
            except StopIteration:
                extra = None
            if extra is not None:
                raise ArtifactError("phase sidecar contains unmatched extra rows")
    except ArtifactError:
        raise
    except (OSError, csv.Error) as error:
        raise ArtifactError("cannot decode phase sidecar CSV: {}".format(error))

    return PairValidation(run_id, len(summaries), phase_rows, case_ids)


def _build_shards(
    run_id: str,
    output_root: Path,
    benchmark_binary: Path,
    total_case_count: int,
    shard_case_count: int,
    warmups: int,
    minimum_repetitions: int,
    maximum_repetitions: int,
    bootstrap_resamples: int,
) -> Tuple[List[Dict[str, object]], List[List[str]]]:
    shard_root = output_root / "shards"
    shards = []
    commands = []
    for command_index, case_start in enumerate(
        range(0, total_case_count, shard_case_count)
    ):
        case_count = min(shard_case_count, total_case_count - case_start)
        case_end = case_start + case_count
        tag = "{:05d}-{:05d}".format(case_start, case_end - 1)
        summary = shard_root / (tag + ".summary.csv")
        phase = shard_root / (tag + ".phase_samples.csv")
        completion = shard_root / (tag + ".done.json")
        command = [
            str(benchmark_binary),
            "--output",
            str(summary) + ".partial",
            "--run-id",
            run_id,
            "--phase-samples-output",
            str(phase) + ".partial",
            "--case-start",
            str(case_start),
            "--case-count",
            str(case_count),
            "--warmups",
            str(warmups),
            "--min-repetitions",
            str(minimum_repetitions),
            "--max-repetitions",
            str(maximum_repetitions),
            "--bootstrap-resamples",
            str(bootstrap_resamples),
            "--progress-every",
            "0",
        ]
        shards.append(
            {
                "case_start": case_start,
                "case_count": case_count,
                "case_end_exclusive": case_end,
                "summary_partial_path": str(summary) + ".partial",
                "phase_samples_partial_path": str(phase) + ".partial",
                "summary_path": str(summary),
                "phase_samples_path": str(phase),
                "completion_path": str(completion),
                "benchmark_command_index": command_index,
            }
        )
        commands.append(command)
    return shards, commands


def build_identity_document(
    *,
    run_id: str,
    source_root: Path,
    source_commit: str,
    benchmark_binary: Path,
    dependency_provenance: Path,
    compiler: str,
    compiler_flags: Sequence[str],
    machine: str,
    platform: str,
    execution_space: str,
    validation_commands: Sequence[str],
    output_root: Path,
    total_case_count: int = 24300,
    shard_case_count: int = 180,
    warmups: int = 10,
    minimum_repetitions: int = 30,
    maximum_repetitions: int = 100,
    bootstrap_resamples: int = 10000,
    retained_contract: bool = True,
    scratch_root: Optional[Path] = None,
) -> Dict[str, object]:
    """Build an identity document whose output paths and commands are exact."""

    binary = Path(_absolute(benchmark_binary))
    provenance = Path(_absolute(dependency_provenance))
    root = Path(_absolute(output_root))
    scratch = Path(_absolute(scratch_root or binary.parent))
    if total_case_count <= 0 or shard_case_count <= 0:
        raise ArtifactError("case counts must be positive")
    if total_case_count % shard_case_count:
        raise ArtifactError("shards must partition the case count exactly")
    shards, commands = _build_shards(
        run_id,
        root,
        binary,
        total_case_count,
        shard_case_count,
        warmups,
        minimum_repetitions,
        maximum_repetitions,
        bootstrap_resamples,
    )
    document = {
        "artifact_schema": IDENTITY_SCHEMA,
        "created_utc": _utc_now(),
        "validation_mode": "STRICT",
        "retained_contract": retained_contract,
        "run_id": run_id,
        "source_root": _absolute(source_root),
        "source_commit": source_commit,
        "source_clean": True,
        "source_detached": True,
        "benchmark_binary_path": str(binary),
        "benchmark_binary_sha256": sha256_file(binary),
        "dependency_provenance_path": str(provenance),
        "dependency_provenance_sha256": sha256_file(provenance),
        "compiler": compiler,
        "compiler_flags": list(compiler_flags),
        "machine": machine,
        "platform": platform,
        "execution_space": execution_space,
        "validation_commands": list(validation_commands),
        "summary_schema": SUMMARY_SCHEMA,
        "phase_samples_schema": PHASE_SCHEMA,
        "scratch_root": str(scratch),
        "output_root": str(root),
        "identity_manifest_path": str(root / "identity.json"),
        "merged_summary_partial_path": str(
            root / "ndrop_summary_v2.merged.csv.partial"
        ),
        "merged_phase_samples_partial_path": str(
            root / "ndrop_phase_samples_v1.merged.csv.partial"
        ),
        "merged_summary_path": str(root / "ndrop_summary_v2.merged.csv"),
        "merged_phase_samples_path": str(
            root / "ndrop_phase_samples_v1.merged.csv"
        ),
        "final_manifest_path": str(root / "final_manifest.json"),
        "total_case_count": total_case_count,
        "shard_case_count": shard_case_count,
        "shard_count": len(shards),
        "worker_count": 1,
        "warmups": warmups,
        "minimum_repetitions": minimum_repetitions,
        "maximum_repetitions": maximum_repetitions,
        "bootstrap_resamples": bootstrap_resamples,
        "benchmark_commands": commands,
        "shards": shards,
    }
    _validate_identity_document(document, require_retained=retained_contract)
    return document


def _validate_string_list(value: object, name: str, allow_empty: bool = False) -> None:
    if not isinstance(value, list) or (not value and not allow_empty):
        raise ArtifactError("{} must be a nonempty list".format(name))
    if any(not isinstance(item, str) or not item for item in value):
        raise ArtifactError("{} must contain nonempty strings".format(name))


def _validate_identity_document(
    document: Mapping, require_retained: bool = False
) -> None:
    if not isinstance(document, Mapping):
        raise ArtifactError("identity manifest must be a JSON object")
    _require_exact_fields(document, _IDENTITY_FIELDS, "identity manifest")
    if document["artifact_schema"] != IDENTITY_SCHEMA:
        raise ArtifactError("identity manifest schema is unsupported")
    if document["validation_mode"] != "STRICT":
        raise ArtifactError("identity validation_mode must be STRICT")
    if not isinstance(document["retained_contract"], bool):
        raise ArtifactError("retained_contract must be boolean")
    if require_retained and document["retained_contract"] is not True:
        raise ArtifactError("a retained-contract identity manifest is required")
    run_id = document["run_id"]
    if not isinstance(run_id, str) or _RUN_ID.fullmatch(run_id) is None:
        raise ArtifactError("identity run_id is invalid")
    source_commit = document["source_commit"]
    if (
        not isinstance(source_commit, str)
        or _SOURCE_COMMIT.fullmatch(source_commit) is None
    ):
        raise ArtifactError("source_commit must be a 40-character lowercase SHA")
    if document["source_clean"] is not True or document["source_detached"] is not True:
        raise ArtifactError("retained source identity must be clean and detached")
    for name in (
        "source_root",
        "benchmark_binary_path",
        "dependency_provenance_path",
        "scratch_root",
        "output_root",
        "identity_manifest_path",
        "merged_summary_partial_path",
        "merged_phase_samples_partial_path",
        "merged_summary_path",
        "merged_phase_samples_path",
        "final_manifest_path",
    ):
        _require_absolute_path(document[name], name)
    for name in (
        "compiler",
        "machine",
        "platform",
        "execution_space",
        "created_utc",
    ):
        if not isinstance(document[name], str) or not document[name]:
            raise ArtifactError("{} must be a nonempty string".format(name))
    _validate_string_list(document["compiler_flags"], "compiler_flags", True)
    _validate_string_list(document["validation_commands"], "validation_commands")
    if document["summary_schema"] != SUMMARY_SCHEMA:
        raise ArtifactError("identity summary schema is unsupported")
    if document["phase_samples_schema"] != PHASE_SCHEMA:
        raise ArtifactError("identity phase schema is unsupported")

    total = _integer(document["total_case_count"], "total_case_count", 1)
    shard_size = _integer(document["shard_case_count"], "shard_case_count", 1)
    shard_count = _integer(document["shard_count"], "shard_count", 1)
    if total % shard_size or shard_count != total // shard_size:
        raise ArtifactError("identity shard geometry does not exactly partition cases")
    for name in (
        "worker_count",
        "warmups",
        "minimum_repetitions",
        "maximum_repetitions",
        "bootstrap_resamples",
    ):
        _integer(document[name], name, 1)
    if document["minimum_repetitions"] > document["maximum_repetitions"]:
        raise ArtifactError("minimum repetitions cannot exceed maximum repetitions")
    if not isinstance(document["shards"], list) or len(document["shards"]) != shard_count:
        raise ArtifactError("identity shard list has the wrong length")
    if (
        not isinstance(document["benchmark_commands"], list)
        or len(document["benchmark_commands"]) != shard_count
    ):
        raise ArtifactError("identity benchmark command list has the wrong length")

    output_root = Path(document["output_root"])
    if document["identity_manifest_path"] != str(output_root / "identity.json"):
        raise ArtifactError("identity manifest path is not canonical")
    canonical_merged_paths = {
        "merged_summary_partial_path": str(
            output_root / "ndrop_summary_v2.merged.csv.partial"
        ),
        "merged_phase_samples_partial_path": str(
            output_root / "ndrop_phase_samples_v1.merged.csv.partial"
        ),
        "merged_summary_path": str(output_root / "ndrop_summary_v2.merged.csv"),
        "merged_phase_samples_path": str(
            output_root / "ndrop_phase_samples_v1.merged.csv"
        ),
        "final_manifest_path": str(output_root / "final_manifest.json"),
    }
    for name, expected in canonical_merged_paths.items():
        if document[name] != expected:
            raise ArtifactError("{} is not canonical".format(name))
    rebuilt_shards, rebuilt_commands = _build_shards(
        run_id,
        output_root,
        Path(document["benchmark_binary_path"]),
        total,
        shard_size,
        document["warmups"],
        document["minimum_repetitions"],
        document["maximum_repetitions"],
        document["bootstrap_resamples"],
    )
    if document["shards"] != rebuilt_shards:
        raise ArtifactError("identity shard paths or indices are not canonical")
    if document["benchmark_commands"] != rebuilt_commands:
        raise ArtifactError("identity benchmark commands are not canonical")

    if document["retained_contract"]:
        match = _RETAINED_RUN_ID.fullmatch(run_id)
        if match is None or match.group(1) != source_commit[:12]:
            raise ArtifactError("retained run ID must contain the source SHA prefix")
        exact = {
            "total_case_count": 24300,
            "shard_case_count": 180,
            "shard_count": 135,
            "worker_count": 1,
            "warmups": 10,
            "minimum_repetitions": 30,
            "maximum_repetitions": 100,
            "bootstrap_resamples": 10000,
        }
        for name, expected in exact.items():
            if document[name] != expected:
                raise ArtifactError(
                    "retained identity {} must equal {}".format(name, expected)
                )
        if shard_size % 9:
            raise ArtifactError("retained shards must align to nine-row groups")
        if not document["compiler_flags"]:
            raise ArtifactError("retained identity requires explicit compiler flags")


def _write_immutable_json(path: Path, document: Mapping) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(document, sort_keys=True, indent=2) + "\n").encode()
    try:
        descriptor = os.open(
            str(destination), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444
        )
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(destination, 0o444)
        _fsync_directory(destination.parent)
    except FileExistsError:
        raise ArtifactError("refusing to overwrite {}".format(destination))
    except OSError as error:
        raise ArtifactError("cannot write {}: {}".format(destination, error))


def _fsync_directory(directory: Path) -> None:
    try:
        descriptor = os.open(str(directory), os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except OSError as error:
        raise ArtifactError("cannot fsync directory {}: {}".format(directory, error))


def seal_identity_manifest(path: Path, document: Mapping) -> None:
    """Validate and exclusively create the pre-run identity manifest."""

    _validate_identity_document(
        document, require_retained=bool(document.get("retained_contract"))
    )
    destination = Path(_absolute(path))
    if str(destination) != document["identity_manifest_path"]:
        raise ArtifactError("identity output path differs from the bound path")
    if sha256_file(Path(document["benchmark_binary_path"])) != document[
        "benchmark_binary_sha256"
    ]:
        raise ArtifactError("benchmark binary hash changed before sealing")
    if sha256_file(Path(document["dependency_provenance_path"])) != document[
        "dependency_provenance_sha256"
    ]:
        raise ArtifactError("dependency provenance hash changed before sealing")
    if document["retained_contract"]:
        commit, clean, detached = _git_identity(Path(document["source_root"]))
        if commit != document["source_commit"] or not clean or not detached:
            raise ArtifactError("source is not the bound clean detached commit")
        if not os.access(document["benchmark_binary_path"], os.X_OK):
            raise ArtifactError("retained benchmark binary is not executable")
    output_root = Path(document["output_root"])
    if output_root.exists() and any(output_root.iterdir()):
        raise ArtifactError("refusing to seal over a nonempty run directory")
    shard_root = Path(document["output_root"]) / "shards"
    shard_root.mkdir(parents=True, exist_ok=True)
    _write_immutable_json(destination, document)


def _read_json(path: Path, artifact_name: str) -> Mapping:
    try:
        with Path(path).open() as stream:
            document = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ArtifactError("cannot read {}: {}".format(artifact_name, error))
    if not isinstance(document, Mapping):
        raise ArtifactError("{} must contain a JSON object".format(artifact_name))
    return document


def _git_identity(source_root: Path) -> Tuple[str, bool, bool]:
    def run(arguments: Sequence[str], allow_failure: bool = False) -> subprocess.CompletedProcess:
        result = subprocess.run(
            list(arguments),
            cwd=str(source_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode and not allow_failure:
            raise ArtifactError(
                "git identity command failed: {}".format(result.stderr.strip())
            )
        return result

    commit = run(("git", "rev-parse", "HEAD")).stdout.strip()
    status = run(
        ("git", "status", "--porcelain", "--untracked-files=all")
    ).stdout
    symbolic = run(("git", "symbolic-ref", "-q", "HEAD"), True)
    detached = symbolic.returncode == 1
    return commit, status == "", detached


def _load_identity(
    path: Path,
    *,
    require_retained: bool = False,
    verify_source: bool = True,
) -> Mapping:
    identity_path = Path(_absolute(path))
    document = _read_json(identity_path, "identity manifest")
    _validate_identity_document(document, require_retained=require_retained)
    if str(identity_path) != document["identity_manifest_path"]:
        raise ArtifactError("identity manifest was moved from its bound path")
    if sha256_file(Path(document["benchmark_binary_path"])) != document[
        "benchmark_binary_sha256"
    ]:
        raise ArtifactError("benchmark binary hash differs from identity")
    if sha256_file(Path(document["dependency_provenance_path"])) != document[
        "dependency_provenance_sha256"
    ]:
        raise ArtifactError("dependency provenance hash differs from identity")
    if verify_source and document["retained_contract"]:
        commit, clean, detached = _git_identity(Path(document["source_root"]))
        if commit != document["source_commit"] or not clean or not detached:
            raise ArtifactError("source is not the bound clean detached commit")
    return document


def _shard_for_start(identity: Mapping, case_start: int) -> Mapping:
    matches = [
        shard
        for shard in identity["shards"]
        if shard["case_start"] == case_start
    ]
    if len(matches) != 1:
        raise ArtifactError("case_start does not identify exactly one shard")
    return matches[0]


def _rename_without_overwrite(source: Path, destination: Path) -> None:
    if destination.exists():
        raise ArtifactError("refusing to overwrite {}".format(destination))
    try:
        os.rename(source, destination)
        _fsync_directory(destination.parent)
    except OSError as error:
        raise ArtifactError(
            "cannot rename {} to {}: {}".format(source, destination, error)
        )


def complete_shard(identity_path: Path, case_start: int) -> CompletedShard:
    """Validate and seal a paired partial shard, writing its marker last."""

    identity = _load_identity(identity_path)
    shard = _shard_for_start(identity, case_start)
    summary_partial = Path(shard["summary_partial_path"])
    phase_partial = Path(shard["phase_samples_partial_path"])
    summary_final = Path(shard["summary_path"])
    phase_final = Path(shard["phase_samples_path"])
    completion = Path(shard["completion_path"])
    completion_partial = Path(str(completion) + ".partial")
    for path in (summary_final, phase_final, completion, completion_partial):
        if path.exists():
            raise ArtifactError("refusing to overwrite existing shard artifact")
    expected_ids = range(
        shard["case_start"], shard["case_end_exclusive"]
    )
    validated = validate_csv_pair(
        summary_partial,
        phase_partial,
        expected_run_id=identity["run_id"],
        expected_case_ids=expected_ids,
        identity=identity,
    )
    summary_hash = sha256_file(summary_partial)
    phase_hash = sha256_file(phase_partial)
    identity_hash = sha256_file(Path(identity["identity_manifest_path"]))
    _rename_without_overwrite(summary_partial, summary_final)
    _rename_without_overwrite(phase_partial, phase_final)
    marker = {
        "artifact_schema": COMPLETION_SCHEMA,
        "created_utc": _utc_now(),
        "run_id": identity["run_id"],
        "identity_manifest_path": identity["identity_manifest_path"],
        "identity_manifest_sha256": identity_hash,
        "case_start": shard["case_start"],
        "case_count": shard["case_count"],
        "case_end_exclusive": shard["case_end_exclusive"],
        "summary_path": str(summary_final),
        "summary_sha256": summary_hash,
        "summary_rows": validated.summary_rows,
        "phase_samples_path": str(phase_final),
        "phase_samples_sha256": phase_hash,
        "phase_rows": validated.phase_rows,
    }
    _write_immutable_json(completion_partial, marker)
    _rename_without_overwrite(completion_partial, completion)
    return CompletedShard(
        str(completion),
        str(summary_final),
        str(phase_final),
        validated.summary_rows,
        validated.phase_rows,
    )


def validate_completed_shard(
    identity_path: Path, case_start: int
) -> CompletedShard:
    """Revalidate a completion marker, both hashes, and the paired rows."""

    identity = _load_identity(identity_path, verify_source=False)
    shard = _shard_for_start(identity, case_start)
    completion_path = Path(shard["completion_path"])
    marker = _read_json(completion_path, "shard completion marker")
    _require_exact_fields(marker, _COMPLETION_FIELDS, "shard completion marker")
    exact = {
        "artifact_schema": COMPLETION_SCHEMA,
        "run_id": identity["run_id"],
        "identity_manifest_path": identity["identity_manifest_path"],
        "identity_manifest_sha256": sha256_file(
            Path(identity["identity_manifest_path"])
        ),
        "case_start": shard["case_start"],
        "case_count": shard["case_count"],
        "case_end_exclusive": shard["case_end_exclusive"],
        "summary_path": shard["summary_path"],
        "phase_samples_path": shard["phase_samples_path"],
    }
    for name, expected in exact.items():
        if marker[name] != expected:
            raise ArtifactError(
                "completion marker {} differs from identity".format(name)
            )
    if not isinstance(marker["created_utc"], str) or not marker["created_utc"]:
        raise ArtifactError("completion marker timestamp is invalid")
    for name in (
        "case_start",
        "case_count",
        "case_end_exclusive",
        "summary_rows",
        "phase_rows",
    ):
        _integer(marker[name], name)
    if sha256_file(Path(marker["summary_path"])) != marker["summary_sha256"]:
        raise ArtifactError("completed summary hash is stale")
    if (
        sha256_file(Path(marker["phase_samples_path"]))
        != marker["phase_samples_sha256"]
    ):
        raise ArtifactError("completed phase sidecar hash is stale")
    validated = validate_csv_pair(
        Path(marker["summary_path"]),
        Path(marker["phase_samples_path"]),
        expected_run_id=identity["run_id"],
        expected_case_ids=range(
            shard["case_start"], shard["case_end_exclusive"]
        ),
        identity=identity,
    )
    if validated.summary_rows != marker["summary_rows"]:
        raise ArtifactError("completion summary count is stale")
    if validated.phase_rows != marker["phase_rows"]:
        raise ArtifactError("completion phase count is stale")
    return CompletedShard(
        str(completion_path),
        marker["summary_path"],
        marker["phase_samples_path"],
        validated.summary_rows,
        validated.phase_rows,
    )


def _merge_csv(inputs: Sequence[Path], output: Path, columns: Tuple[str, ...]) -> None:
    try:
        with output.open("x", newline="") as output_stream:
            writer = csv.DictWriter(output_stream, fieldnames=columns)
            writer.writeheader()
            for source in inputs:
                with source.open(newline="") as input_stream:
                    reader = csv.DictReader(input_stream, strict=True)
                    _validate_header(reader, columns, str(source))
                    for row in reader:
                        _validate_decoded_row(row, str(source))
                        writer.writerow(row)
            output_stream.flush()
            os.fsync(output_stream.fileno())
        _fsync_directory(output.parent)
    except FileExistsError:
        raise ArtifactError("refusing to overwrite {}".format(output))
    except (OSError, csv.Error) as error:
        raise ArtifactError("cannot merge CSV artifacts: {}".format(error))


def merge_completed_shards(identity_path: Path) -> FinalRun:
    """Merge every validated shard and create the final manifest last."""

    identity = _load_identity(identity_path)
    completed = [
        validate_completed_shard(identity_path, shard["case_start"])
        for shard in identity["shards"]
    ]
    summary_partial = Path(identity["merged_summary_partial_path"])
    phase_partial = Path(identity["merged_phase_samples_partial_path"])
    summary_final = Path(identity["merged_summary_path"])
    phase_final = Path(identity["merged_phase_samples_path"])
    final_manifest = Path(identity["final_manifest_path"])
    final_partial = Path(str(final_manifest) + ".partial")
    for path in (
        summary_partial,
        phase_partial,
        summary_final,
        phase_final,
        final_manifest,
        final_partial,
    ):
        if path.exists():
            raise ArtifactError("refusing to overwrite a merged artifact")
    _merge_csv(
        [Path(item.summary_path) for item in completed],
        summary_partial,
        SUMMARY_COLUMNS,
    )
    _merge_csv(
        [Path(item.phase_samples_path) for item in completed],
        phase_partial,
        PHASE_COLUMNS,
    )
    validated = validate_csv_pair(
        summary_partial,
        phase_partial,
        expected_run_id=identity["run_id"],
        expected_case_ids=range(identity["total_case_count"]),
        identity=identity,
    )
    summary_hash = sha256_file(summary_partial)
    phase_hash = sha256_file(phase_partial)
    _rename_without_overwrite(summary_partial, summary_final)
    _rename_without_overwrite(phase_partial, phase_final)
    completion_records = [
        {
            "path": item.completion_path,
            "sha256": sha256_file(Path(item.completion_path)),
        }
        for item in completed
    ]
    document = {
        "artifact_schema": FINAL_SCHEMA,
        "created_utc": _utc_now(),
        "run_id": identity["run_id"],
        "identity_manifest_path": identity["identity_manifest_path"],
        "identity_manifest_sha256": sha256_file(
            Path(identity["identity_manifest_path"])
        ),
        "completion_records": completion_records,
        "summary_path": str(summary_final),
        "summary_sha256": summary_hash,
        "summary_rows": validated.summary_rows,
        "phase_samples_path": str(phase_final),
        "phase_samples_sha256": phase_hash,
        "phase_rows": validated.phase_rows,
    }
    _write_immutable_json(final_partial, document)
    _rename_without_overwrite(final_partial, final_manifest)
    return FinalRun(
        str(final_manifest),
        str(summary_final),
        str(phase_final),
        validated.summary_rows,
        validated.phase_rows,
    )


def validate_final_run(
    identity_path: Path,
    final_manifest_path: Path,
    *,
    require_retained_contract: bool = True,
) -> FinalRun:
    """Revalidate every retained-run binding from identity through merge."""

    identity = _load_identity(
        identity_path,
        require_retained=require_retained_contract,
        verify_source=False,
    )
    final_path = Path(_absolute(final_manifest_path))
    if str(final_path) != identity["final_manifest_path"]:
        raise ArtifactError("final manifest was moved from its bound path")
    document = _read_json(final_path, "final manifest")
    _require_exact_fields(document, _FINAL_FIELDS, "final manifest")
    exact = {
        "artifact_schema": FINAL_SCHEMA,
        "run_id": identity["run_id"],
        "identity_manifest_path": identity["identity_manifest_path"],
        "identity_manifest_sha256": sha256_file(
            Path(identity["identity_manifest_path"])
        ),
        "summary_path": identity["merged_summary_path"],
        "phase_samples_path": identity["merged_phase_samples_path"],
    }
    for name, expected in exact.items():
        if document[name] != expected:
            raise ArtifactError("final manifest {} is stale".format(name))
    if not isinstance(document["created_utc"], str) or not document["created_utc"]:
        raise ArtifactError("final manifest timestamp is invalid")
    _integer(document["summary_rows"], "summary_rows")
    _integer(document["phase_rows"], "phase_rows")
    records = document["completion_records"]
    if not isinstance(records, list) or len(records) != identity["shard_count"]:
        raise ArtifactError("final manifest completion coverage is incomplete")
    for index, shard in enumerate(identity["shards"]):
        record = records[index]
        expected_path = shard["completion_path"]
        if (
            not isinstance(record, Mapping)
            or frozenset(record) != frozenset({"path", "sha256"})
            or record["path"] != expected_path
            or record["sha256"] != sha256_file(Path(expected_path))
        ):
            raise ArtifactError("final completion record hash or order is stale")
        validate_completed_shard(identity_path, shard["case_start"])
    if sha256_file(Path(document["summary_path"])) != document["summary_sha256"]:
        raise ArtifactError("merged summary hash is stale")
    if (
        sha256_file(Path(document["phase_samples_path"]))
        != document["phase_samples_sha256"]
    ):
        raise ArtifactError("merged phase sidecar hash is stale")
    validated = validate_csv_pair(
        Path(document["summary_path"]),
        Path(document["phase_samples_path"]),
        expected_run_id=identity["run_id"],
        expected_case_ids=range(identity["total_case_count"]),
        identity=identity,
    )
    if validated.summary_rows != document["summary_rows"]:
        raise ArtifactError("final summary row count is stale")
    if validated.phase_rows != document["phase_rows"]:
        raise ArtifactError("final phase row count is stale")
    return FinalRun(
        str(final_path),
        document["summary_path"],
        document["phase_samples_path"],
        validated.summary_rows,
        validated.phase_rows,
    )


def _preflight_free_space(identity: Mapping) -> None:
    for name in ("scratch_root", "output_root"):
        path = Path(identity[name])
        try:
            free = shutil.disk_usage(path).free
        except OSError as error:
            raise ArtifactError("cannot check free space for {}: {}".format(path, error))
        if free < MINIMUM_FREE_BYTES:
            raise ArtifactError("{} has less than 5 GiB free".format(path))


@contextmanager
def _exclusive_run_lock(identity: Mapping):
    lock_path = Path(identity["output_root"]) / "benchmark_worker.lock"
    try:
        descriptor = os.open(
            str(lock_path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444
        )
    except FileExistsError:
        raise ArtifactError("another shard worker is active or left a stale lock")
    except OSError as error:
        raise ArtifactError("cannot create shard-worker lock: {}".format(error))
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write("pid={}\n".format(os.getpid()))
            stream.flush()
            os.fsync(stream.fileno())
        _fsync_directory(lock_path.parent)
        yield
    finally:
        try:
            lock_path.unlink()
            _fsync_directory(lock_path.parent)
        except FileNotFoundError:
            pass


def run_shard(identity_path: Path, case_start: int) -> CompletedShard:
    """Reverify provenance, execute one exact command, and seal its pair."""

    identity = _load_identity(identity_path, require_retained=True)
    _preflight_free_space(identity)
    shard = _shard_for_start(identity, case_start)
    paths = (
        shard["summary_partial_path"],
        shard["phase_samples_partial_path"],
        shard["summary_path"],
        shard["phase_samples_path"],
        shard["completion_path"],
    )
    if any(Path(path).exists() for path in paths):
        raise ArtifactError("refusing to overwrite an existing shard artifact")
    command = identity["benchmark_commands"][shard["benchmark_command_index"]]
    with _exclusive_run_lock(identity):
        result = subprocess.run(command, check=False)
        if result.returncode:
            raise ArtifactError(
                "benchmark shard exited with status {}".format(result.returncode)
            )
        return complete_shard(identity_path, case_start)


def _create_identity_from_arguments(arguments: argparse.Namespace) -> None:
    source_root = Path(_absolute(arguments.source_root))
    commit, clean, detached = _git_identity(source_root)
    if not clean or not detached:
        raise ArtifactError("source must be clean and detached before sealing")
    document = build_identity_document(
        run_id=arguments.run_id,
        source_root=source_root,
        source_commit=commit,
        benchmark_binary=Path(arguments.binary),
        dependency_provenance=Path(arguments.dependency_provenance),
        compiler=arguments.compiler,
        compiler_flags=arguments.compiler_flag,
        machine=arguments.machine,
        platform=arguments.platform,
        execution_space=arguments.execution_space,
        validation_commands=arguments.validation_command,
        output_root=Path(arguments.output_root),
        scratch_root=Path(arguments.scratch_root),
        retained_contract=True,
    )
    seal_identity_manifest(Path(document["identity_manifest_path"]), document)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create-identity")
    create.add_argument("--run-id", required=True)
    create.add_argument("--source-root", required=True)
    create.add_argument("--binary", required=True)
    create.add_argument("--dependency-provenance", required=True)
    create.add_argument("--compiler", required=True)
    create.add_argument("--compiler-flag", action="append", default=[])
    create.add_argument("--machine", required=True)
    create.add_argument("--platform", required=True)
    create.add_argument("--execution-space", required=True)
    create.add_argument("--validation-command", action="append", required=True)
    create.add_argument("--scratch-root", required=True)
    create.add_argument("--output-root", required=True)

    for name in ("run-shard", "complete-shard"):
        command = subparsers.add_parser(name)
        command.add_argument("--identity", required=True)
        command.add_argument("--case-start", required=True, type=int)
    merge = subparsers.add_parser("merge")
    merge.add_argument("--identity", required=True)
    validate = subparsers.add_parser("validate-final")
    validate.add_argument("--identity", required=True)
    validate.add_argument("--final-manifest", required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        if arguments.command == "create-identity":
            _create_identity_from_arguments(arguments)
        elif arguments.command == "run-shard":
            run_shard(Path(arguments.identity), arguments.case_start)
        elif arguments.command == "complete-shard":
            complete_shard(Path(arguments.identity), arguments.case_start)
        elif arguments.command == "merge":
            merge_completed_shards(Path(arguments.identity))
        elif arguments.command == "validate-final":
            validate_final_run(
                Path(arguments.identity), Path(arguments.final_manifest)
            )
        else:
            raise ArtifactError("unsupported command")
    except ArtifactError as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
