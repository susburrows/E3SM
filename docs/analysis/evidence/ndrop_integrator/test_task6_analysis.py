"""Synthetic end-to-end tests for the fail-closed Task 6 analysis path."""

import csv
import importlib
import json
from pathlib import Path
import tempfile
import unittest

from analyze_screening import GroupKey
from ndrop_artifacts import (
    PHASE_COLUMNS,
    SUMMARY_COLUMNS,
    ArtifactError,
    build_identity_document,
    complete_shard,
    load_validated_analysis_input,
    merge_completed_shards,
    seal_identity_manifest,
)
from test_analyze_screening import _phase_rows
from test_ndrop_artifacts import _CsvOpenMonitor, _summary_row


class Task6AnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.binary = self.root / "benchmark"
        self.provenance = self.root / "dependencies.json"
        self.binary.write_bytes(b"synthetic-benchmark")
        self.provenance.write_text('{"locked":true}\n')
        self.run_root = self.root / "run"
        self.run_root.mkdir()
        self.identity_path = self.run_root / "identity.json"
        identity = build_identity_document(
            run_id="task6-synthetic",
            source_root=self.root / "source",
            source_commit="2" * 40,
            benchmark_binary=self.binary,
            dependency_provenance=self.provenance,
            compiler="clang test",
            compiler_flags=("-O3",),
            machine="synthetic-machine",
            platform="macos-arm64",
            execution_space="Serial",
            validation_commands=("synthetic validation",),
            output_root=self.run_root,
            total_case_count=2,
            shard_case_count=2,
            minimum_repetitions=3,
            maximum_repetitions=3,
            bootstrap_resamples=4,
            retained_contract=False,
        )
        seal_identity_manifest(self.identity_path, identity)
        self.identity = json.loads(self.identity_path.read_text())
        shard = self.identity["shards"][0]
        summary_path = Path(shard["summary_partial_path"])
        phase_path = Path(shard["phase_samples_partial_path"])
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        first = _summary_row(0, self.identity["run_id"])
        first["dt"] = "50"
        second = _summary_row(1, self.identity["run_id"])
        second.update(
            {
                "dt": "100",
                "method": "backward_euler",
                "linear_phase_mode": "separate",
                "factorization_seconds_per_column": "0.25",
                "solve_seconds_per_column": "0.35",
                "combined_factorization_solve_seconds_per_column": "",
            }
        )
        summaries = (first, second)
        with summary_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=SUMMARY_COLUMNS)
            writer.writeheader()
            writer.writerows(summaries)
        with phase_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=PHASE_COLUMNS)
            writer.writeheader()
            for summary in summaries:
                writer.writerows(_phase_rows(summary))
        complete_shard(self.identity_path, 0)
        self.final_run = merge_completed_shards(self.identity_path)
        self.group = GroupKey(3, 9, 16, 1, 0.5)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _run(self, final_manifest: Path, output: Path) -> Path:
        try:
            module = importlib.import_module("task6_analysis")
        except ModuleNotFoundError:
            self.fail("production Task 6 analysis callable is missing")
        return module.run_task6_analysis(
            self.identity_path,
            final_manifest,
            output,
            group=self.group,
            first_axis=("imex", None),
            second_axis=("backward_euler", None),
            target_error=0.003,
            require_retained_contract=False,
        )

    def test_validated_fixture_writes_deterministic_unsupported_result(self) -> None:
        final_manifest = Path(self.final_run.final_manifest_path)
        with _CsvOpenMonitor() as validation_opens:
            paired = load_validated_analysis_input(
                self.identity_path,
                final_manifest,
                require_retained_contract=False,
            )
        self.assertEqual(len(paired.dataset.points), 2)
        self.assertEqual(len(paired.phase_diagnostics), 2)

        first_output = self.root / "first-analysis.json"
        second_output = self.root / "second-analysis.json"
        with _CsvOpenMonitor() as analysis_opens:
            returned = self._run(final_manifest, first_output)
        self._run(final_manifest, second_output)

        expected = (
            '{\n  "comparison": {\n'
            '    "status": "unsupported_cost_basis"\n'
            "  }\n}\n"
        )
        self.assertEqual(returned, first_output)
        self.assertEqual(first_output.read_text(), expected)
        self.assertEqual(second_output.read_bytes(), first_output.read_bytes())
        self.assertEqual(analysis_opens, validation_opens)
        comparison = json.loads(first_output.read_text())["comparison"]
        self.assertEqual(comparison, {"status": "unsupported_cost_basis"})

    def test_invalid_manifest_schema_fails_before_csv_open_or_output(self) -> None:
        final_manifest = Path(self.final_run.final_manifest_path)
        document = json.loads(final_manifest.read_text())
        document["artifact_schema"] = "unsupported"
        final_manifest.chmod(0o644)
        final_manifest.write_text(json.dumps(document))
        output = self.root / "analysis.json"

        with _CsvOpenMonitor() as csv_opens:
            with self.assertRaises(ArtifactError):
                self._run(final_manifest, output)

        self.assertEqual(csv_opens, [])
        self.assertFalse(output.exists())

    def test_partial_manifest_path_fails_before_csv_open_or_output(self) -> None:
        partial_manifest = Path(self.final_run.final_manifest_path + ".partial")
        output = self.root / "analysis.json"

        with _CsvOpenMonitor() as csv_opens:
            with self.assertRaises(ArtifactError):
                self._run(partial_manifest, output)

        self.assertEqual(csv_opens, [])
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
