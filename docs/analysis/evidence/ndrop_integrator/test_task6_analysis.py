"""Synthetic end-to-end tests for the fail-closed Task 6 analysis path."""

import csv
import hashlib
import importlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from analyze_screening import (
    AnalysisDataset,
    GroupKey,
    ScreeningPoint,
    TimingEstimate,
)
from ndrop_artifacts import (
    PHASE_COLUMNS,
    SUMMARY_COLUMNS,
    ArtifactError,
    build_identity_document,
    complete_shard,
    expected_retained_case,
    load_validated_analysis_input,
    merge_completed_shards,
    seal_identity_manifest,
    sha256_file,
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

    def test_direct_matrix_reports_only_measured_same_timestep_ratios(self) -> None:
        module = importlib.import_module("task6_analysis")

        def point(case_id, method, dt, median, lower, upper, fe_nsubmix=None):
            return ScreeningPoint(
                case_id,
                self.group,
                method,
                dt,
                fe_nsubmix,
                0.003,
                TimingEstimate(median, lower, upper),
            )

        dataset = AnalysisDataset(
            (
                point(1, "imex", 50.0, 2.0, 1.8, 2.2),
                point(2, "backward_euler", 50.0, 4.0, 3.8, 4.2),
                point(3, "forward_euler", 50.0, 2.1, 1.7, 2.4, 1),
                point(4, "imex", 100.0, 2.0, 1.8, 2.2),
                point(5, "forward_euler", 100.0, 1.0, 0.9, 1.1, 2),
            ),
            (),
        )

        payload = module.build_direct_same_timestep_payload(dataset)

        self.assertEqual(
            payload["status_counts"],
            {
                "determinate": 3,
                "indeterminate": 1,
                "unavailable": 6,
                "unsupported_cost_basis": 0,
            },
        )
        imex_be_50 = next(
            item
            for item in payload["comparisons"]
            if item["pairing"] == "imex_vs_backward_euler"
            and item["dt_seconds"] == 50.0
        )
        self.assertEqual(imex_be_50["result"]["status"], "determinate")
        self.assertEqual(
            imex_be_50["result"]["runtime_ratio_first_over_second"], 0.5
        )
        self.assertAlmostEqual(
            imex_be_50["result"]["runtime_ratio_ci95_lower"], 1.8 / 4.2
        )
        self.assertAlmostEqual(
            imex_be_50["result"]["runtime_ratio_ci95_upper"], 2.2 / 3.8
        )
        unavailable = [
            item
            for item in payload["comparisons"]
            if item["result"]["status"] == "unavailable"
        ]
        self.assertTrue(unavailable)
        self.assertTrue(
            all(set(item["result"]) == {"status"} for item in unavailable)
        )
        self.assertEqual(payload["input_exclusion_count"], 0)

    def test_cross_timestep_direct_cost_is_status_only_unsupported(self) -> None:
        module = importlib.import_module("task6_analysis")
        first = ScreeningPoint(
            1,
            self.group,
            "imex",
            50.0,
            None,
            0.003,
            TimingEstimate(2.0, 1.8, 2.2),
        )
        second = ScreeningPoint(
            2,
            self.group,
            "backward_euler",
            100.0,
            None,
            0.003,
            TimingEstimate(4.0, 3.8, 4.2),
        )

        result = module._direct_comparison_result(first, second)

        self.assertEqual(result, {"status": "unsupported_cost_basis"})

    def test_retained_direct_matrix_reuses_receipt_and_validated_boundary(self) -> None:
        module = importlib.import_module("task6_analysis")
        output_directory = self.root / "direct-retained-output"
        with mock.patch.object(
            module,
            "load_validated_analysis_input",
            return_value=self.paired,
        ) as loader, _CsvOpenMonitor() as csv_opens:
            completed = module.run_retained_direct_task6_analysis(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                output_directory,
                analysis_source_commit=self.analysis_source_commit,
                analysis_command=self.analysis_command,
            )

        loader.assert_called_once_with(
            self.identity_path,
            Path(self.final_run.final_manifest_path),
            require_retained_contract=True,
        )
        self.assertEqual(csv_opens, [])
        validated = module.validate_retained_direct_task6_output(
            self.identity_path,
            Path(self.final_run.final_manifest_path),
            output_directory,
            analysis_source_commit=self.analysis_source_commit,
            analysis_command=self.analysis_command,
        )
        self.assertEqual(validated, completed)
        payload = json.loads(completed.analysis_path.read_text())
        self.assertEqual(payload["analysis_mode"], "direct_same_timestep")
        self.assertEqual(payload["status_counts"]["unsupported_cost_basis"], 0)

    @staticmethod
    def _golden_production_nsubmix(nlev, mixing_scale, dt_seconds):
        intervals = (33.33, 50.0, 75.0, 100.0, 150.0, 300.0)
        lower = {
            1: (1, 1, 2, 2, 3, 5),
            2: (3, 4, 6, 8, 11, 22),
            3: (3, 5, 7, 9, 14, 27),
        }
        higher = {
            1: (5, 7, 11, 14, 21, 42),
            2: (25, 37, 55, 73, 109, 217),
            3: (30, 44, 66, 88, 132, 263),
        }
        table = lower if mixing_scale == 0.05 else higher
        fallback = (
            (4, 6, 9, 12, 18, 36)
            if mixing_scale == 0.05
            else (39, 59, 88, 117, 176, 351)
        )
        return table.get(nlev, fallback)[intervals.index(dt_seconds)]

    def _synthetic_direct_payload(self):
        comparisons = []
        measured = {1, 2, 3, 5, 8, 16, 21}
        for case_id in range(0, 24300, 9):
            case = expected_retained_case(case_id)
            nsubmix = self._golden_production_nsubmix(
                case["nlev"], case["mixing_scale"], case["dt"]
            )
            if nsubmix not in measured:
                continue
            group = {
                "nlev": case["nlev"],
                "nonzero_state_fields": case["nonzero_state_fields"],
                "batch_size": case["batch_size"],
                "mixing_index": case["mixing_index"],
                "mixing_scale": case["mixing_scale"],
            }
            first_axis = {
                "method": "forward_euler",
                "fe_nsubmix": nsubmix,
            }
            for comparator in ("imex", "backward_euler"):
                result = (
                    {
                        "status": "determinate",
                        "runtime_ratio_first_over_second": 0.5,
                        "runtime_ratio_ci95_lower": 0.4,
                        "runtime_ratio_ci95_upper": 0.6,
                        "winner": first_axis,
                    }
                    if comparator == "imex"
                    else {
                        "status": "indeterminate",
                        "runtime_ratio_first_over_second": 1.0,
                        "runtime_ratio_ci95_lower": 0.8,
                        "runtime_ratio_ci95_upper": 1.2,
                    }
                )
                comparisons.append(
                    {
                        "pairing": "forward_euler_nsubmix_{}_vs_{}".format(
                            nsubmix, comparator
                        ),
                        "group": group,
                        "dt_seconds": case["dt"],
                        "first_axis": first_axis,
                        "second_axis": {"method": comparator},
                        "result": result,
                    }
                )
        return {
            "analysis_mode": "direct_same_timestep",
            "timing_source": "primary_uninstrumented",
            "comparisons": comparisons,
        }

    def _seal_synthetic_direct_output(self, module, payload, directory):
        directory.mkdir()
        analysis_path = directory / "analysis.json"
        receipt_path = directory / "analysis.done.json"
        analysis_path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")
        receipt = {
            "artifact_schema": module.RETAINED_RECEIPT_SCHEMA,
            "completion_status": "complete",
            "identity_manifest_sha256": sha256_file(self.identity_path),
            "final_manifest_sha256": sha256_file(
                Path(self.final_run.final_manifest_path)
            ),
            "analysis_source_commit": "4" * 40,
            "analysis_parameters": module._direct_analysis_parameters(),
            "analysis_command": ["synthetic-direct-analysis"],
            "analysis_output": "analysis.json",
            "analysis_output_sha256": sha256_file(analysis_path),
        }
        receipt_path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
        return analysis_path, receipt_path

    def test_production_nsubmix_matches_complete_retained_source_map(self) -> None:
        module = importlib.import_module("task6_analysis")
        levels = (1, 2, 3, 16, 32, 64, 72, 127, 128)
        intervals = (33.33, 50.0, 75.0, 100.0, 150.0, 300.0)
        for mixing_scale in (0.05, 0.5):
            for nlev in levels:
                for dt_seconds in intervals:
                    with self.subTest(
                        mixing_scale=mixing_scale,
                        nlev=nlev,
                        dt_seconds=dt_seconds,
                    ):
                        self.assertEqual(
                            module.reconstruct_production_nsubmix(
                                nlev, mixing_scale, dt_seconds
                            ),
                            self._golden_production_nsubmix(
                                nlev, mixing_scale, dt_seconds
                            ),
                        )

    def test_first_level_controller_drives_single_level_submixing(self) -> None:
        module = importlib.import_module("task6_analysis")
        self.assertEqual(module.reconstruct_production_nsubmix(1, 0.05, 75.0), 2)
        self.assertEqual(module.reconstruct_production_nsubmix(1, 0.5, 33.33), 5)

    def test_computed_nsubmix_payload_has_complete_exact_only_denominator(self) -> None:
        module = importlib.import_module("task6_analysis")
        payload = module.build_computed_nsubmix_exact_payload(
            self._synthetic_direct_payload()
        )

        self.assertEqual(len(payload["computed_nsubmix_map"]), 108)
        self.assertEqual(len(payload["comparisons"]), 5400)
        self.assertEqual(len(payload["strata"]), 216)
        self.assertTrue(
            all(item["designed_workload_count"] == 25 for item in payload["strata"])
        )
        self.assertEqual(
            payload["design_counts"],
            {
                "workload_timestep_cases": 2700,
                "pairwise_comparison_opportunities": 5400,
                "exact_match_workload_opportunities": 300,
                "exact_match_comparison_opportunities": 600,
                "valid_direct_comparisons": 600,
            },
        )
        self.assertEqual(
            payload["support_counts"],
            {
                "exact_measured": 600,
                "not_measured_in_range": 1950,
                "outside_measured_range": 2850,
            },
        )
        self.assertEqual(payload["status_counts"]["determinate"], 300)
        self.assertEqual(payload["status_counts"]["indeterminate"], 300)
        self.assertEqual(payload["status_counts"]["unavailable"], 4800)
        unavailable = [
            item
            for item in payload["comparisons"]
            if item["result"]["status"] == "unavailable"
        ]
        self.assertTrue(unavailable)
        self.assertTrue(
            all(set(item["result"]) == {"status"} for item in unavailable)
        )
        self.assertFalse(payload["adaptive_selection_path_timed"])
        self.assertEqual(payload["timing_source"], "primary_uninstrumented")
        self.assertEqual(
            [item["label"] for item in payload["mixing_conditions"]],
            [
                "Lower-intensity mixing (scale 0.05)",
                "Higher-intensity mixing (scale 0.5)",
            ],
        )

        missing = self._synthetic_direct_payload()
        missing["comparisons"].pop(0)
        missing_payload = module.build_computed_nsubmix_exact_payload(missing)
        self.assertEqual(
            missing_payload["design_counts"]["valid_direct_comparisons"], 599
        )
        self.assertEqual(missing_payload["status_counts"]["unavailable"], 4801)

    def test_computed_nsubmix_payload_rejects_nonprimary_parent_data(self) -> None:
        module = importlib.import_module("task6_analysis")
        mutations = {
            "wrong_mode": ("analysis_mode", "matched_error"),
            "instrumented": ("timing_source", "instrumented"),
            "phase_diagnostics": ("phase_diagnostics", []),
        }
        for name, (field, value) in mutations.items():
            with self.subTest(name=name):
                payload = self._synthetic_direct_payload()
                payload[field] = value
                with self.assertRaises(ArtifactError):
                    module.build_computed_nsubmix_exact_payload(payload)

    def test_computed_nsubmix_payload_rejects_bad_selected_records(self) -> None:
        module = importlib.import_module("task6_analysis")
        mutations = {
            "duplicate": lambda records: records.append(dict(records[0])),
            "reversed": lambda records: records[0].update(
                {
                    "first_axis": records[0]["second_axis"],
                    "second_axis": records[0]["first_axis"],
                }
            ),
            "wrong_winner": lambda records: records[0]["result"].update(
                {"winner": {"method": "imex"}}
            ),
            "unsupported": lambda records: records[0].update(
                {"result": {"status": "unsupported_cost_basis"}}
            ),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                payload = self._synthetic_direct_payload()
                mutate(payload["comparisons"])
                with self.assertRaises((ArtifactError, ValueError)):
                    module.build_computed_nsubmix_exact_payload(payload)

    def test_retained_computed_nsubmix_output_is_parent_bound(self) -> None:
        module = importlib.import_module("task6_analysis")
        direct_directory = self.root / "direct-parent"
        parent_analysis, parent_receipt = self._seal_synthetic_direct_output(
            module, self._synthetic_direct_payload(), direct_directory
        )
        output_directory = self.root / "computed-output"
        command = ("task6_analysis.run_retained_computed_nsubmix_exact_analysis",)
        with _CsvOpenMonitor() as csv_opens:
            completed = module.run_retained_computed_nsubmix_exact_analysis(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                direct_directory,
                output_directory,
                direct_analysis_source_commit="4" * 40,
                direct_analysis_command=("synthetic-direct-analysis",),
                analysis_source_commit="5" * 40,
                analysis_command=command,
            )
        self.assertEqual(csv_opens, [])
        validated = module.validate_retained_computed_nsubmix_exact_output(
            self.identity_path,
            Path(self.final_run.final_manifest_path),
            direct_directory,
            output_directory,
            direct_analysis_source_commit="4" * 40,
            direct_analysis_command=("synthetic-direct-analysis",),
            analysis_source_commit="5" * 40,
            analysis_command=command,
        )
        self.assertEqual(validated, completed)
        receipt = json.loads(completed.receipt_path.read_text())
        self.assertEqual(
            receipt["analysis_parameters"]["parent_analysis_output_sha256"],
            sha256_file(parent_analysis),
        )
        self.assertEqual(
            receipt["analysis_parameters"]["parent_analysis_receipt_sha256"],
            sha256_file(parent_receipt),
        )

        derived_bytes = completed.analysis_path.read_bytes()
        completed.analysis_path.chmod(0o644)
        completed.analysis_path.write_bytes(derived_bytes + b"tampered\n")
        with self.assertRaises(ArtifactError):
            module.validate_retained_computed_nsubmix_exact_output(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                direct_directory,
                output_directory,
                direct_analysis_source_commit="4" * 40,
                direct_analysis_command=("synthetic-direct-analysis",),
                analysis_source_commit="5" * 40,
                analysis_command=command,
            )
        completed.analysis_path.write_bytes(derived_bytes)

        parent_receipt_bytes = parent_receipt.read_bytes()
        parent_receipt.chmod(0o644)
        self._mutate_receipt(parent_receipt, "completion_status", "incomplete")
        with self.assertRaises(ArtifactError):
            module.validate_retained_computed_nsubmix_exact_output(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                direct_directory,
                output_directory,
                direct_analysis_source_commit="4" * 40,
                direct_analysis_command=("synthetic-direct-analysis",),
                analysis_source_commit="5" * 40,
                analysis_command=command,
            )
        parent_receipt.write_bytes(parent_receipt_bytes)

        parent_analysis.chmod(0o644)
        parent_analysis.write_bytes(parent_analysis.read_bytes() + b"tampered\n")
        with self.assertRaises(ArtifactError):
            module.validate_retained_computed_nsubmix_exact_output(
                self.identity_path,
                Path(self.final_run.final_manifest_path),
                direct_directory,
                output_directory,
                direct_analysis_source_commit="4" * 40,
                direct_analysis_command=("synthetic-direct-analysis",),
                analysis_source_commit="5" * 40,
                analysis_command=command,
            )

    @staticmethod
    def _mutate_receipt(path: Path, name: str, value: str) -> None:
        document = json.loads(path.read_text())
        document[name] = value
        path.write_text(json.dumps(document, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
