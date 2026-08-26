"""Fail-closed tests for Task 5R artifact transactions and provenance."""

import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from ndrop_artifacts import (
    ArtifactError,
    PHASE_COLUMNS,
    SUMMARY_COLUMNS,
    build_identity_document,
    complete_shard,
    expected_retained_case,
    merge_completed_shards,
    prepare_retained_inputs,
    seal_identity_manifest,
    validate_completed_shard,
    validate_csv_pair,
    validate_final_run,
)
from test_analyze_screening import _phase_rows, _valid_row


def _summary_row(case_id: int, run_id: str) -> dict:
    row = dict.fromkeys(SUMMARY_COLUMNS, "")
    row.update(
        _valid_row(
            case_id=str(case_id),
            run_id=run_id,
            dt=str(25 + case_id),
        )
    )
    row.update(
        {
            "ci95_relative_half_width": "0.25",
            "instrumented_ci95_lower": "1",
            "instrumented_ci95_upper": "3",
            "instrumentation_overhead": "0",
            "block_factorizations": "0",
            "block_rhs_solves": "0",
            "scalar_factorizations": "2",
            "scalar_rhs_solves": "51",
            "executed_field_level_updates": "153",
            "workspace_bytes": "4096",
            "checksum": "0123456789abcdef",
            "compiler": "clang test",
            "platform": "macos-arm64",
            "kokkos_execution_space": "Serial",
            "hardware_threads": "8",
            "warmups": "10",
            "minimum_repetitions": "30",
            "maximum_repetitions": "100",
            "bootstrap_resamples": "10000",
        }
    )
    return row


def _write_pair(
    summary_path: Path, phase_path: Path, case_ids: range, run_id: str
) -> None:
    summaries = [_summary_row(case_id, run_id) for case_id in case_ids]
    with summary_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        writer.writerows(summaries)
    with phase_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=PHASE_COLUMNS)
        writer.writeheader()
        for summary in summaries:
            writer.writerows(_phase_rows(summary))


class PairValidationTests(unittest.TestCase):
    def test_retained_case_decoder_has_exact_method_and_matrix_corners(self) -> None:
        first = expected_retained_case(0)
        self.assertEqual(
            (
                first["nlev"],
                first["nonzero_aerosol_pairs"],
                first["batch_size"],
                first["mixing_index"],
                first["dt"],
                first["method"],
                first["forced_fe_substeps"],
            ),
            (1, 1, 1, 0, 33.33, "forward_euler", 1),
        )
        self.assertEqual(expected_retained_case(6)["forced_fe_substeps"], 21)
        self.assertEqual(expected_retained_case(7)["method"], "imex")
        self.assertEqual(expected_retained_case(8)["method"], "backward_euler")
        last = expected_retained_case(24299)
        self.assertEqual(
            (
                last["nlev"],
                last["nonzero_aerosol_pairs"],
                last["batch_size"],
                last["mixing_index"],
                last["dt"],
                last["method"],
            ),
            (128, 25, 1024, 1, 300.0, "backward_euler"),
        )
        with self.assertRaises(ArtifactError):
            expected_retained_case(24300)

    def test_retained_identity_has_exact_shard_and_command_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / "benchmark"
            provenance = root / "dependencies.json"
            binary.write_bytes(b"binary")
            binary.chmod(0o755)
            provenance.write_text("{}\n")
            commit = "a" * 40
            retained_binary, retained_provenance = prepare_retained_inputs(
                root / "run", binary, provenance
            )
            document = build_identity_document(
                run_id=(
                    "task5-remediated-20260825T120000Z-" + commit[:12]
                ),
                source_root=root / "source",
                source_commit=commit,
                benchmark_binary=retained_binary,
                dependency_provenance=retained_provenance,
                compiler="clang test",
                compiler_flags=("-O3", "-DNDEBUG"),
                machine="test-machine",
                platform="macos-arm64",
                execution_space="Serial",
                validation_commands=("ctest --output-on-failure",),
                output_root=root / "run",
            )
            self.assertEqual(document["total_case_count"], 24300)
            self.assertEqual(document["shard_case_count"], 180)
            self.assertEqual(document["shard_count"], 135)
            self.assertEqual(len(document["benchmark_commands"]), 135)
            self.assertEqual(document["shards"][-1]["case_start"], 24120)
            self.assertEqual(
                document["shards"][-1]["case_end_exclusive"], 24300
            )
            self.assertEqual(
                Path(document["benchmark_binary_path"]).parent,
                root / "run" / "inputs",
            )
            self.assertEqual(
                Path(document["dependency_provenance_path"]).parent,
                root / "run" / "inputs",
            )
            with self.assertRaises(ArtifactError):
                prepare_retained_inputs(root / "run", binary, provenance)

            with self.assertRaises(ArtifactError):
                build_identity_document(
                    run_id="task5-remediated-bad",
                    source_root=root / "source",
                    source_commit=commit,
                    benchmark_binary=binary,
                    dependency_provenance=provenance,
                    compiler="clang test",
                    compiler_flags=("-O3",),
                    machine="test-machine",
                    platform="macos-arm64",
                    execution_space="Serial",
                    validation_commands=("ctest",),
                    output_root=root / "bad-run",
                )

    def test_exact_headers_pairing_and_expected_range_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary = root / "summary.csv"
            phase = root / "phase.csv"
            _write_pair(summary, phase, range(3), "task5-test")

            result = validate_csv_pair(
                summary,
                phase,
                expected_run_id="task5-test",
                expected_case_ids=range(3),
            )
            self.assertEqual(result.summary_rows, 3)
            self.assertEqual(result.phase_rows, 9)
            self.assertEqual(result.case_ids, (0, 1, 2))

            header = summary.read_text().splitlines()[0]
            summary.write_text(header + ",unexpected\n" + "\n".join(
                line + "," for line in summary.read_text().splitlines()[1:]
            ) + "\n")
            with self.assertRaises(ArtifactError):
                validate_csv_pair(summary, phase)

    def test_mixed_run_and_wrong_case_range_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary = root / "summary.csv"
            phase = root / "phase.csv"
            _write_pair(summary, phase, range(3), "task5-test")

            with self.assertRaises(ArtifactError):
                validate_csv_pair(
                    summary,
                    phase,
                    expected_run_id="another-run",
                    expected_case_ids=range(3),
                )
            with self.assertRaises(ArtifactError):
                validate_csv_pair(
                    summary,
                    phase,
                    expected_run_id="task5-test",
                    expected_case_ids=range(1, 4),
                )


class TransactionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.binary = self.root / "benchmark"
        self.provenance = self.root / "dependencies.json"
        self.binary.write_bytes(b"benchmark-binary")
        self.provenance.write_text('{"locked":true}\n')
        self.run_root = self.root / "run"
        self.run_root.mkdir()
        self.identity_path = self.run_root / "identity.json"
        document = build_identity_document(
            run_id="task5-test",
            source_root=self.root / "source",
            source_commit="1" * 40,
            benchmark_binary=self.binary,
            dependency_provenance=self.provenance,
            compiler="clang test",
            compiler_flags=("-O3", "-DNDEBUG"),
            machine="test-machine",
            platform="macos-arm64",
            execution_space="Serial",
            validation_commands=("ctest --output-on-failure",),
            output_root=self.run_root,
            total_case_count=6,
            shard_case_count=3,
            retained_contract=False,
        )
        seal_identity_manifest(self.identity_path, document)
        self.identity = json.loads(self.identity_path.read_text())

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write_shard(self, case_start: int) -> None:
        shard = next(
            item
            for item in self.identity["shards"]
            if item["case_start"] == case_start
        )
        summary = Path(shard["summary_partial_path"])
        phase = Path(shard["phase_samples_partial_path"])
        summary.parent.mkdir(parents=True, exist_ok=True)
        _write_pair(
            summary,
            phase,
            range(case_start, case_start + shard["case_count"]),
            self.identity["run_id"],
        )

    def test_identity_is_exclusive_and_binds_binary_and_provenance(self) -> None:
        with self.assertRaises(ArtifactError):
            seal_identity_manifest(self.identity_path, self.identity)
        self.assertEqual(
            self.identity["benchmark_binary_sha256"],
            hashlib.sha256(self.binary.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            self.identity["dependency_provenance_sha256"],
            hashlib.sha256(self.provenance.read_bytes()).hexdigest(),
        )
        self.assertEqual(len(self.identity["shards"]), 2)
        self.assertEqual(len(self.identity["benchmark_commands"]), 2)

    def test_completion_marker_is_last_and_hash_tampering_is_rejected(self) -> None:
        self._write_shard(0)
        completion = complete_shard(self.identity_path, 0)
        marker_path = Path(completion.completion_path)
        self.assertTrue(marker_path.is_file())
        validated = validate_completed_shard(self.identity_path, 0)
        self.assertEqual(validated.summary_rows, 3)
        self.assertEqual(validated.phase_rows, 9)

        Path(completion.summary_path).write_text("tampered\n")
        with self.assertRaises(ArtifactError):
            validate_completed_shard(self.identity_path, 0)

    def test_malformed_pair_never_gets_a_completion_marker(self) -> None:
        self._write_shard(0)
        shard = self.identity["shards"][0]
        Path(shard["phase_samples_partial_path"]).write_text("bad-header\n")

        with self.assertRaises(ArtifactError):
            complete_shard(self.identity_path, 0)
        self.assertFalse(Path(shard["completion_path"]).exists())
        self.assertTrue(Path(shard["summary_partial_path"]).exists())

    def test_merge_is_deterministic_and_final_manifest_revalidates(self) -> None:
        for case_start in (0, 3):
            self._write_shard(case_start)
            complete_shard(self.identity_path, case_start)

        final = merge_completed_shards(self.identity_path)
        validated = validate_final_run(
            self.identity_path,
            final.final_manifest_path,
            require_retained_contract=False,
        )
        self.assertEqual(validated.summary_rows, 6)
        self.assertEqual(validated.phase_rows, 18)
        with Path(final.summary_path).open(newline="") as stream:
            case_ids = [int(row["case_id"]) for row in csv.DictReader(stream)]
        self.assertEqual(case_ids, list(range(6)))

        with self.assertRaises(ArtifactError):
            merge_completed_shards(self.identity_path)

    def test_missing_marker_means_shard_does_not_exist_for_merge(self) -> None:
        for case_start in (0, 3):
            self._write_shard(case_start)
        complete_shard(self.identity_path, 0)

        with self.assertRaises(ArtifactError):
            merge_completed_shards(self.identity_path)


if __name__ == "__main__":
    unittest.main()
