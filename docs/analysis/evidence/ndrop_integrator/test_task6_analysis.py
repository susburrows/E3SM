"""Synthetic end-to-end tests for the fail-closed Task 6 analysis path."""

import csv
import hashlib
import importlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

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
        self.paired = load_validated_analysis_input(
            self.identity_path,
            Path(self.final_run.final_manifest_path),
            require_retained_contract=False,
        )
        self.analysis_source_commit = "3" * 40
        self.analysis_command = (
            "python3",
            "task6_analysis.py",
            "--synthetic-retained-test",
        )

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

    def _run_retained(self, output_directory: Path):
        module = importlib.import_module("task6_analysis")
        with mock.patch.object(
            module,
            "load_validated_analysis_input",
            return_value=self.paired,
        ) as loader, mock.patch.object(
            module,
            "compare_matched_error_cost",
            wraps=module.compare_matched_error_cost,
        ) as comparison:
            result = module.run_retained_task6_analysis(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                output_directory,
                group=self.group,
                first_axis=("imex", None),
                second_axis=("backward_euler", None),
                target_error=0.003,
                analysis_source_commit=self.analysis_source_commit,
                analysis_command=self.analysis_command,
            )
        return module, result, loader, comparison

    def _validate_retained(self, module, output_directory: Path):
        return module.validate_retained_analysis_output(
            self.identity_path,
            Path(self.final_run.final_manifest_path),
            output_directory,
            group=self.group,
            first_axis=("imex", None),
            second_axis=("backward_euler", None),
            target_error=0.003,
            analysis_source_commit=self.analysis_source_commit,
            analysis_command=self.analysis_command,
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

    def test_retained_entry_point_forces_contract_and_uses_only_dataset(self) -> None:
        output_directory = self.root / "retained-output"
        with _CsvOpenMonitor() as csv_opens:
            module, result, loader, comparison = self._run_retained(
                output_directory
            )

        loader.assert_called_once_with(
            self.identity_path,
            Path(self.final_run.final_manifest_path),
            require_retained_contract=True,
        )
        self.assertEqual(csv_opens, [])
        self.assertIs(comparison.call_args.args[0], self.paired.dataset)
        self.assertNotIn(
            self.paired.phase_diagnostics, comparison.call_args.args
        )
        self.assertEqual(result.analysis_path, output_directory / "analysis.json")
        self.assertEqual(
            result.receipt_path, output_directory / "analysis.done.json"
        )

        with self.assertRaises(TypeError):
            module.run_retained_task6_analysis(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                self.root / "disabled-contract",
                group=self.group,
                first_axis=("imex", None),
                second_axis=("backward_euler", None),
                target_error=0.003,
                analysis_source_commit=self.analysis_source_commit,
                analysis_command=self.analysis_command,
                require_retained_contract=False,
            )

    def test_synthetic_entry_point_rejects_retained_mode(self) -> None:
        module = importlib.import_module("task6_analysis")
        output = self.root / "legacy-retained.json"
        with mock.patch.object(
            module, "load_validated_analysis_input"
        ) as loader, self.assertRaises(ArtifactError):
            module.run_task6_analysis(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                output,
                group=self.group,
                first_axis=("imex", None),
                second_axis=("backward_euler", None),
                target_error=0.003,
                require_retained_contract=True,
            )

        loader.assert_not_called()
        self.assertFalse(output.exists())

    def test_output_promotion_is_atomic_and_never_replaces(self) -> None:
        module = importlib.import_module("task6_analysis")
        source = self.root / "promotion.partial"
        destination = self.root / "promotion.json"
        source.write_text("candidate")
        real_link = module.os.link

        def collide(source_path, destination_path):
            Path(destination_path).write_text("competitor")
            return real_link(source_path, destination_path)

        with mock.patch.object(module.os, "link", side_effect=collide):
            with self.assertRaises(ArtifactError):
                module._rename_without_overwrite(source, destination)

        self.assertEqual(destination.read_text(), "competitor")
        self.assertEqual(source.read_text(), "candidate")

        dangling = self.root / "dangling.json"
        dangling.symlink_to(self.root / "missing-target")
        with self.assertRaises(ArtifactError):
            module._rename_without_overwrite(source, dangling)
        self.assertTrue(dangling.is_symlink())
        self.assertEqual(source.read_text(), "candidate")

    def test_retained_invalid_manifest_fails_before_csv_or_receipt(self) -> None:
        module = importlib.import_module("task6_analysis")
        final_manifest = Path(self.final_run.final_manifest_path)
        document = json.loads(final_manifest.read_text())
        document["artifact_schema"] = "unsupported"
        final_manifest.chmod(0o644)
        final_manifest.write_text(json.dumps(document))
        output_directory = self.root / "invalid-retained-output"

        def synthetic_strict_loader(
            identity, final, *, require_retained_contract
        ):
            self.assertTrue(require_retained_contract)
            return load_validated_analysis_input(
                identity, final, require_retained_contract=False
            )

        with mock.patch.object(
            module,
            "load_validated_analysis_input",
            side_effect=synthetic_strict_loader,
        ), _CsvOpenMonitor() as csv_opens:
            with self.assertRaises(ArtifactError):
                module.run_retained_task6_analysis(
                    self.identity_path,
                    final_manifest,
                    output_directory,
                    group=self.group,
                    first_axis=("imex", None),
                    second_axis=("backward_euler", None),
                    target_error=0.003,
                    analysis_source_commit=self.analysis_source_commit,
                    analysis_command=self.analysis_command,
                )

        self.assertEqual(csv_opens, [])
        self.assertFalse((output_directory / "analysis.json").exists())
        self.assertFalse((output_directory / "analysis.done.json").exists())

    def test_retained_partial_manifest_fails_before_csv_or_receipt(self) -> None:
        module = importlib.import_module("task6_analysis")
        partial_manifest = Path(self.final_run.final_manifest_path + ".partial")
        output_directory = self.root / "partial-retained-output"

        def synthetic_strict_loader(
            identity, final, *, require_retained_contract
        ):
            self.assertTrue(require_retained_contract)
            return load_validated_analysis_input(
                identity, final, require_retained_contract=False
            )

        with mock.patch.object(
            module,
            "load_validated_analysis_input",
            side_effect=synthetic_strict_loader,
        ), _CsvOpenMonitor() as csv_opens:
            with self.assertRaises(ArtifactError):
                module.run_retained_task6_analysis(
                    self.identity_path,
                    partial_manifest,
                    output_directory,
                    group=self.group,
                    first_axis=("imex", None),
                    second_axis=("backward_euler", None),
                    target_error=0.003,
                    analysis_source_commit=self.analysis_source_commit,
                    analysis_command=self.analysis_command,
                )

        self.assertEqual(csv_opens, [])
        self.assertFalse((output_directory / "analysis.json").exists())
        self.assertFalse((output_directory / "analysis.done.json").exists())

    def test_analysis_failure_produces_no_final_output_or_receipt(self) -> None:
        module = importlib.import_module("task6_analysis")
        output_directory = self.root / "analysis-failure"
        with mock.patch.object(
            module,
            "load_validated_analysis_input",
            return_value=self.paired,
        ), mock.patch.object(
            module,
            "compare_matched_error_cost",
            side_effect=RuntimeError("synthetic analysis failure"),
        ):
            with self.assertRaisesRegex(RuntimeError, "synthetic analysis failure"):
                module.run_retained_task6_analysis(
                    self.identity_path,
                    Path(self.final_run.final_manifest_path),
                    output_directory,
                    group=self.group,
                    first_axis=("imex", None),
                    second_axis=("backward_euler", None),
                    target_error=0.003,
                    analysis_source_commit=self.analysis_source_commit,
                    analysis_command=self.analysis_command,
                )

        self.assertFalse((output_directory / "analysis.json").exists())
        self.assertFalse((output_directory / "analysis.done.json").exists())

    def test_retained_output_and_receipt_are_deterministic_and_bound(self) -> None:
        first_directory = self.root / "retained-first"
        second_directory = self.root / "retained-second"
        module, first, _, _ = self._run_retained(first_directory)
        _, second, _, _ = self._run_retained(second_directory)

        self.assertEqual(
            first.analysis_path.read_bytes(), second.analysis_path.read_bytes()
        )
        self.assertEqual(
            first.receipt_path.read_bytes(), second.receipt_path.read_bytes()
        )
        self.assertFalse((first_directory / "analysis.json.partial").exists())
        self.assertFalse((first_directory / "analysis.done.json.partial").exists())
        receipt = json.loads(first.receipt_path.read_text())
        self.assertEqual(receipt["completion_status"], "complete")
        self.assertEqual(
            receipt["analysis_output_sha256"],
            hashlib.sha256(first.analysis_path.read_bytes()).hexdigest(),
        )
        validated = self._validate_retained(module, first_directory)
        self.assertEqual(validated, first)

        original_output = first.analysis_path.read_bytes()
        original_receipt = first.receipt_path.read_bytes()
        with self.assertRaises(ArtifactError):
            self._run_retained(first_directory)
        self.assertEqual(first.analysis_path.read_bytes(), original_output)
        self.assertEqual(first.receipt_path.read_bytes(), original_receipt)

    def test_partial_or_tampered_retained_output_fails_closed(self) -> None:
        module = importlib.import_module("task6_analysis")
        partial_directory = self.root / "interrupted"
        partial_directory.mkdir()
        (partial_directory / "analysis.json.partial").write_text("partial")
        with self.assertRaises(ArtifactError):
            self._validate_retained(module, partial_directory)

        mutations = {
            "output": lambda output, receipt: output.write_bytes(
                output.read_bytes() + b"tampered\n"
            ),
            "binding": lambda output, receipt: self._mutate_receipt(
                receipt, "final_manifest_sha256", "0" * 64
            ),
            "receipt": lambda output, receipt: self._mutate_receipt(
                receipt, "completion_status", "incomplete"
            ),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                directory = self.root / ("tamper-" + name)
                _, completed, _, _ = self._run_retained(directory)
                completed.analysis_path.chmod(0o644)
                completed.receipt_path.chmod(0o644)
                mutate(completed.analysis_path, completed.receipt_path)
                with self.assertRaises(ArtifactError):
                    self._validate_retained(module, directory)

    @staticmethod
    def _mutate_receipt(path: Path, name: str, value: str) -> None:
        document = json.loads(path.read_text())
        document[name] = value
        path.write_text(json.dumps(document, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
