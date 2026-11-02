"""Synthetic failure and pairing fixtures for the standalone study artifacts."""

import copy
import json
import os
from pathlib import Path
import subprocess
import shutil
import tempfile
import unittest

import ndrop_matched_artifacts as artifacts


def accuracy_rows():
    rows = []
    for fixture in artifacts.FIXTURES:
        levels, mixing, interval = ((72, .05, 50), (72, .5, 150),
                                    (128, .05, 50), (128, .5, 150))[fixture]
        for method in artifacts.METHODS:
            for power in range(14):
                steps = 1 << power
                if method == "forward_euler" and steps < (6 if fixture % 2 == 0 else 176):
                    continue
                row = dict.fromkeys(artifacts.ACCURACY_HEADER, "")
                row.update(artifact_schema="ndrop_matched_accuracy_v1",
                           fixture=str(fixture), levels=str(levels),
                           mixing_scale=str(mixing), interval_seconds=str(interval),
                           mixing_rate_t="1", activation_rate_t="0.5",
                           method=method, steps=str(steps),
                           reference_status="certified", reference_reason="",
                           reference_disagreement="0", status="numerical_failure",
                           failure_reason="synthetic failure",
                           clipping_count="0", material_negative_count="0",
                           normalized_residual="0", normalized_ledger_error="0",
                           state_checksum="1")
                rows.append(row)
    return rows


class MatchedArtifactsTest(unittest.TestCase):
    def test_complete_matrix_and_accuracy_selection(self):
        rows = accuracy_rows()
        first = next(row for row in rows if row["fixture"] == "0" and
                     row["method"] == "forward_euler" and row["steps"] == "8")
        first.update(status="eligible", failure_reason="",
                     error_normal="0.001",
                     error_tight="0.001", error_max="0.001",
                     normal_weighted_l1_by_field=";".join(["0.001"] * 51),
                     normal_linf_by_field=";".join(["0.001"] * 51),
                     tight_weighted_l1_by_field=";".join(["0.001"] * 51),
                     tight_linf_by_field=";".join(["0.001"] * 51))
        indexed = artifacts.validate_accuracy(rows)
        frozen = artifacts.expected_selection(indexed, "a" * 64, "b" * 64)
        self.assertIn([0, "forward_euler", 8], frozen["unique_configs"])
        self.assertEqual(len(frozen["selected"]), 60)
        self.assertEqual(sum(item["status"] == "selected"
                             for item in frozen["selected"]), 3)
        selected = next(item for item in frozen["selected"]
                        if item["status"] == "selected")
        self.assertEqual(selected["achieved_error"], 0.001)
        self.assertEqual(selected["target_to_achieved_error_ratio"], 10.)

        corrupt = copy.deepcopy(rows)
        damaged = next(row for row in corrupt if row["fixture"] == "0" and
                       row["method"] == "forward_euler" and row["steps"] == "8")
        damaged["normal_linf_by_field"] = "1"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(corrupt)
        missing = rows[:-1]
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(missing)
        duplicate = rows + [rows[-1]]
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(duplicate)
        mixed = copy.deepcopy(rows)
        mixed[0]["artifact_schema"] = "ndrop_summary_v2"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(mixed)
        bad_gate = copy.deepcopy(rows)
        bad_gate_row = next(row for row in bad_gate if row["fixture"] == "0" and
                            row["method"] == "forward_euler" and
                            row["steps"] == "8")
        bad_gate_row["normalized_residual"] = "2e-11"
        with self.assertRaisesRegex(artifacts.ArtifactError,
                                    "ineligible numerical state"):
            artifacts.validate_accuracy(bad_gate)
        missing_reason = copy.deepcopy(rows)
        failed_row = next(row for row in missing_reason
                          if row["status"] == "numerical_failure")
        failed_row["failure_reason"] = ""
        with self.assertRaisesRegex(artifacts.ArtifactError,
                                    "lacks its first reason"):
            artifacts.validate_accuracy(missing_reason)
        false_reference_reason = copy.deepcopy(rows)
        false_reference_reason[0]["reference_reason"] = "fabricated"
        with self.assertRaisesRegex(artifacts.ArtifactError,
                                    "certified reference"):
            artifacts.validate_accuracy(false_reference_reason)

    def test_paired_blocks_and_bootstrap(self):
        frozen = {"unique_configs": [[0, "forward_euler", 8],
                                     [0, "imex", 16]]}
        accuracy = {(0, "forward_euler", 8): {"state_checksum": "7"},
                    (0, "imex", 16): {"state_checksum": "7"}}
        rows = []
        for method, steps, duration in (("forward_euler", 8, 2.),
                                        ("imex", 16, 1.)):
            for block in range(30):
                rows.append({"artifact_schema": "ndrop_matched_timing_v1",
                             "fixture": "0", "method": method,
                             "steps": str(steps), "block": str(block),
                             "seconds_per_column": str(duration),
                             "state_checksum": "7",
                             "stop_reason": "median_precision"})
        artifacts.validate_timing(rows, frozen, accuracy)
        result = artifacts.paired_ratio([2.] * 30, [1.] * 30, 0x4e44524f)
        self.assertEqual((result["ratio"], result["lower"], result["upper"]),
                         (2., 2., 2.))
        broken = copy.deepcopy(rows)
        broken[-1]["block"] = "29"
        broken[-2]["block"] = "29"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_timing(broken, frozen, accuracy)
        broken = copy.deepcopy(rows)
        broken[-1]["state_checksum"] = "8"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_timing(broken, frozen, accuracy)
        stale_accuracy = copy.deepcopy(accuracy)
        stale_accuracy[(0, "imex", 16)]["state_checksum"] = "8"
        with self.assertRaisesRegex(artifacts.ArtifactError,
                                    "timed state differs"):
            artifacts.validate_timing(rows, frozen, stale_accuracy)

    def test_benchmark_median_contract(self):
        self.assertEqual(artifacts.median([1., 2., 3., 4.]), 2.5)
        self.assertEqual(artifacts.median([3., 1., 2.]), 2.)

    def test_completion_requires_current_identity_and_hash(self):
        with tempfile.TemporaryDirectory(prefix="ndrop-matched-artifact-") as folder:
            directory = Path(folder)
            (directory / "identity.json").write_text("{}\n", encoding="utf-8")
            (directory / "accuracy.csv").write_text("payload\n", encoding="utf-8")
            with self.assertRaises(FileNotFoundError):
                artifacts.completed(directory, "accuracy")
            record = {"artifact_schema": "ndrop_matched_completion_v1",
                      "artifact_sha256": artifacts.sha256(directory / "accuracy.csv"),
                      "identity_sha256": artifacts.sha256(directory / "identity.json"),
                      "path": "accuracy.csv"}
            (directory / "accuracy.complete.json").write_text(
                json.dumps(record), encoding="utf-8")
            self.assertEqual(artifacts.completed(directory, "accuracy"),
                             directory / "accuracy.csv")
            (directory / "accuracy.csv").write_text("corrupted\n", encoding="utf-8")
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.completed(directory, "accuracy")


def physical_rows():
    """Fabricated accuracy fixtures; never production or executed evidence."""
    rows = []
    for fixture, levels, strength, diffusivity, velocity in (
            (28, 72, "moderate", 3, .3), (29, 72, "strong", 30, 1.),
            (30, 128, "moderate", 3, .3), (31, 128, "strong", 30, 1.)):
        def row(kind, method, steps, error, status="eligible"):
            result = dict.fromkeys(artifacts.PHYSICAL_ACCURACY_HEADER, "")
            result.update(
                artifact_schema="ndrop_physical_column_accuracy_v1",
                record_kind=kind, fixture=str(fixture), levels=str(levels),
                turbulence=strength, background_diffusivity_m2_s=str(diffusivity),
                cloud_base_wtke_m_s=str(velocity), interval_seconds="150",
                mixing_rate_t="1.0", activation_rate_t="0.2",
                controller_rate_t="1.1", cloud_base_level="60",
                cloud_base_spacing_m="50", limiting_level="60",
                limiting_spacing_m="50", method=method, steps=str(steps),
                activation_limiter_count="0",
                production_nsubmix="2", reference_status="certified",
                reference_disagreement="0", status=status)
            if kind != "terminal":
                result.update(error_normal=str(error), error_tight=str(error),
                              error_max=str(error), clipping_count="0",
                              material_negative_count="0", normalized_residual="0",
                              normalized_ledger_error="0", state_checksum="7")
                for name in artifacts.PHYSICAL_ERROR_VECTORS:
                    result[name] = ";".join([str(error)] * 51)
            return result
        rows.append(row("baseline", "forward_euler", 2, .04))
        for method in artifacts.METHODS:
            start = 2 if method == "forward_euler" else 1
            for steps, error in zip(range(start, start + 3), (.04, .008, .0008)):
                rows.append(row("measurement", method, steps, error))
            rows.append(row("terminal", method, start + 2, None, "selected"))
            for steps in sorted({2 * (start + 1), 4 * (start + 1),
                                 2 * (start + 2), 4 * (start + 2)}):
                rows.append(row("confirmation", method, steps, .0036 / steps))
    return rows


class PhysicalColumnTest(unittest.TestCase):
    def test_limiter_outcome_is_present_nonnegative_and_fixture_bound(self):
        self.assertIn("activation_limiter_count", artifacts.PHYSICAL_ACCURACY_HEADER)
        rows = physical_rows()
        self.assertTrue(all(r["activation_limiter_count"] == "0" for r in rows))
        rows[0]["activation_limiter_count"] = "-1"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_physical_accuracy(rows)
        rows[0]["activation_limiter_count"] = "1"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_physical_accuracy(rows)

    def test_common_ceiling_does_not_substitute_controller_baseline_error(self):
        indexed = artifacts.validate_physical_accuracy(physical_rows())
        frozen = artifacts.physical_selection(indexed, "a" * 64, "b" * 64)
        self.assertEqual(frozen["targets"], [0.01, 0.001])
        self.assertEqual(len(frozen["selected"]), 24)
        self.assertEqual(len(frozen["baselines"]), 4)
        chosen = [r for r in frozen["selected"] if r["fixture"] == 28]
        self.assertEqual([(r["method"], r["target"], r["steps"]) for r in chosen],
                         [("forward_euler", .01, 3), ("forward_euler", .001, 4),
                          ("imex", .01, 2), ("imex", .001, 3),
                          ("backward_euler", .01, 2), ("backward_euler", .001, 3)])
        self.assertEqual(frozen["baselines"][0]["steps"], 2)
        self.assertEqual(frozen["baselines"][0]["achieved_error"], .04)
        self.assertTrue(all(r["timing_admitted"] for r in frozen["selected"]))

    def test_incomplete_search_bad_reference_and_fabricated_metrics_fail_closed(self):
        original = physical_rows()
        mutations = (
            lambda rows: rows.pop(2),
            lambda rows: rows.append(dict(rows[0])),
            lambda rows: rows[0].update(reference_disagreement="1e-9"),
            lambda rows: rows[0].update(production_nsubmix="256"),
            lambda rows: rows[0].update(normal_linf_by_field="0;0"),
            lambda rows: rows[0].update(artifact_schema="ndrop_summary_v2"),
            lambda rows: rows[0].update(state_checksum=""),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                rows = copy.deepcopy(original)
                mutation(rows)
                with self.assertRaises(artifacts.ArtifactError):
                    artifacts.validate_physical_accuracy(rows)

    def test_missing_or_nonconvergent_confirmation_withholds_timing_only(self):
        rows = physical_rows()
        rows = [r for r in rows if not (r["fixture"] == "28" and
                r["method"] == "imex" and r["record_kind"] == "confirmation")]
        frozen = artifacts.physical_selection(
            artifacts.validate_physical_accuracy(rows), "a" * 64, "b" * 64)
        imex = [r for r in frozen["selected"] if r["fixture"] == 28 and
                r["method"] == "imex"]
        self.assertTrue(all(r["status"] == "selected" for r in imex))
        self.assertTrue(all(not r["timing_admitted"] for r in imex))
        self.assertTrue(all(r["exclusion_reason"] == "missing_confirmation"
                            for r in imex))
        rows = physical_rows()
        for row in rows:
            if row["record_kind"] == "confirmation":
                row.update(error_normal="0.1", error_tight="0.1", error_max="0.1")
                for name in artifacts.PHYSICAL_ERROR_VECTORS:
                    row[name] = ";".join(["0.1"] * 51)
        frozen = artifacts.physical_selection(
            artifacts.validate_physical_accuracy(rows), "a" * 64, "b" * 64)
        self.assertTrue(all(r["exclusion_reason"] == "nondecreasing_confirmation"
                            for r in frozen["selected"]))

    def test_reference_outcome_and_confirmation_identity_cannot_change_between_rows(self):
        rows = physical_rows()
        confirmation = next(r for r in rows if r["record_kind"] == "confirmation")
        confirmation["reference_disagreement"] = "1e-14"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_physical_accuracy(rows)
        rows = physical_rows()
        extra = dict(next(r for r in rows if r["record_kind"] == "confirmation" and
                          r["method"] == "imex"))
        extra["steps"] = "2"
        rows.append(extra)
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_physical_accuracy(rows)

    def test_failed_candidate_metrics_are_blank_and_counts_are_retained(self):
        rows = physical_rows()
        failed = next(r for r in rows if r["record_kind"] == "measurement" and
                      r["method"] == "imex" and r["steps"] == "1")
        failed.update(status="material_negative", failure_reason="synthetic negative fixture",
                      material_negative_count="1")
        for name in ("error_normal", "error_tight", "error_max", *artifacts.PHYSICAL_ERROR_VECTORS):
            failed[name] = ""
        artifacts.validate_physical_accuracy(rows)
        failed["normalized_residual"] = ""
        failed["normalized_ledger_error"] = ""
        artifacts.validate_physical_accuracy(rows)
        failed["normalized_residual"] = "inf"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_physical_accuracy(rows)
        failed["normalized_residual"] = ""
        failed["error_max"] = "0"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_physical_accuracy(rows)

    def test_eligible_rows_cannot_consistently_omit_state_checksums(self):
        rows = physical_rows()
        for row in rows:
            row["state_checksum"] = ""
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_physical_accuracy(rows)

    def test_physical_completion_is_separate_and_hash_bound(self):
        with tempfile.TemporaryDirectory(prefix="ndrop-physical-artifact-") as folder:
            directory = Path(folder)
            identity = {"artifact_schema": "ndrop_physical_column_identity_v1"}
            (directory / "identity.json").write_bytes(artifacts.canonical_json(identity))
            (directory / "accuracy.csv").write_bytes(b"synthetic fixture,not evidence\n")
            with self.assertRaises(FileNotFoundError):
                artifacts.completed(directory, "accuracy")
            record = {"artifact_schema": "ndrop_physical_column_completion_v1",
                      "artifact_sha256": artifacts.sha256(directory / "accuracy.csv"),
                      "identity_sha256": artifacts.sha256(directory / "identity.json"),
                      "path": "accuracy.csv"}
            (directory / "accuracy.complete.json").write_bytes(artifacts.canonical_json(record))
            self.assertEqual(artifacts.completed(directory, "accuracy"), directory / "accuracy.csv")
            record["artifact_schema"] = "ndrop_matched_completion_v1"
            (directory / "accuracy.complete.json").write_bytes(artifacts.canonical_json(record))
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.completed(directory, "accuracy")

    def test_physical_timing_requires_matching_states_and_complete_blocks(self):
        accuracy = {(28, "forward_euler", 3): {"state_checksum": "7"},
                    (28, "imex", 2): {"state_checksum": "7"}}
        frozen = {"unique_configs": [[28, "forward_euler", 3], [28, "imex", 2]]}
        rows = [{"artifact_schema": "ndrop_physical_column_timing_v1",
                 "fixture": "28", "method": method, "steps": str(steps),
                 "block": str(block), "seconds_per_column": str(duration),
                 "state_checksum": "7", "stop_reason": "median_precision"}
                for method, steps, duration in (("forward_euler", 3, 2.), ("imex", 2, 1.))
                for block in range(30)]
        artifacts.validate_timing(rows, frozen, accuracy, physical=True)
        broken = copy.deepcopy(rows)
        broken[-1]["state_checksum"] = "9"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_timing(broken, frozen, accuracy, physical=True)
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_timing(rows[:-1], frozen, accuracy, physical=True)

    @unittest.skipUnless(os.environ.get("NDROP_PHYSICAL_TEST_BINARY"),
                         "requires explicitly supplied untimed fixture-check binary")
    def test_repository_grid_fixture_checks(self):
        root = Path(__file__).resolve().parents[4]
        result = subprocess.run(
            [os.environ["NDROP_PHYSICAL_TEST_BINARY"],
             "--check-physical-column-fixtures", "--physical-grid-dir",
             str(root / "components/homme/test/vcoord")],
            capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("failed=0", result.stdout)

    @unittest.skipUnless(os.environ.get("NDROP_PHYSICAL_TEST_BINARY"),
                         "requires local build for copied-input integrity fixtures")
    def test_physical_identity_rejects_copied_grid_binary_and_source_corruption(self):
        """Local integrity fixtures, not a fabricated execution/approval record."""
        root = Path(__file__).resolve().parents[4]
        binary = Path(os.environ["NDROP_PHYSICAL_TEST_BINARY"])
        with tempfile.TemporaryDirectory(prefix="ndrop-physical-identity-") as folder:
            directory = Path(folder)
            (directory / "inputs/grids").mkdir(parents=True)
            shutil.copyfile(binary, directory / "inputs/benchmark_ndrop_matched_error")
            grids = {}
            for name in artifacts.PHYSICAL_GRIDS:
                path = root / "components/homme/test/vcoord" / name
                shutil.copyfile(path, directory / "inputs/grids" / name)
                grids[name] = artifacts.sha256(path)
            identity = dict(
                artifact_schema="ndrop_physical_column_identity_v1",
                fixture_policy="dycoms_inspired_isothermal_repository_grid_v1",
                fixtures=list(artifacts.PHYSICAL_FIXTURES), targets=list(artifacts.PHYSICAL_TARGETS),
                source_head="synthetic identity fixture, not execution evidence",
                source_root=str(root), source_status="exploratory_source", source_attached_head=True,
                source_files_sha256={name: artifacts.sha256(
                    root / "docs/analysis/evidence/ndrop_integrator" / name)
                    for name in artifacts.SOURCES},
                binary_sha256=artifacts.sha256(binary), build=artifacts.build_identity(binary),
                design_sha256=artifacts.sha256(root / "docs/design/eamxx_mam4xx_ndrop.md"),
                python_version=artifacts.platform.python_version(), platform=artifacts.platform.platform(),
                grid_files_sha256=grids,
                accuracy_command=["inputs/benchmark_ndrop_matched_error", "--physical-column-accuracy",
                                  "--physical-grid-dir", "inputs/grids", "--max-steps", "8192",
                                  "--budget-seconds", "7200"],
                accuracy_protocol=dict(integer_search=True, maximum_steps=8192,
                                      budget_seconds=7200, confirmation_multipliers=[2, 4]),
                timing_protocol=dict(warmups=10, minimum_blocks=30, maximum_blocks=100,
                                     budget_seconds=7200, complete_interleaved_blocks=True))
            manifest = directory / "identity.json"
            manifest.write_bytes(artifacts.canonical_json(identity))
            artifacts.validate_physical_identity(directory)
            for path in (directory / "inputs/grids/acme-72i.ascii",
                         directory / "inputs/benchmark_ndrop_matched_error"):
                original = path.read_bytes()
                path.write_bytes(original + b"corrupt")
                with self.assertRaises(artifacts.ArtifactError):
                    artifacts.validate_physical_identity(directory)
                path.write_bytes(original)
            identity["source_files_sha256"]["ndrop_integrator.cpp"] = "0" * 64
            manifest.write_bytes(artifacts.canonical_json(identity))
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_physical_identity(directory)

    @unittest.skipUnless(os.environ.get("NDROP_PHYSICAL_TEST_BINARY"),
                         "requires explicitly supplied untimed fixture-check binary")
    def test_missing_extent_nonfinite_and_ordered_grid_corruption(self):
        root = Path(__file__).resolve().parents[4] / "components/homme/test/vcoord"
        for damage in ("extent", "nonfinite", "order", "missing", "extra"):
            with self.subTest(damage=damage), tempfile.TemporaryDirectory(
                    prefix="ndrop-grid-corruption-") as folder:
                directory = Path(folder)
                for name in artifacts.PHYSICAL_GRIDS:
                    shutil.copyfile(root / name, directory / name)
                path = directory / "acme-72i.ascii"
                if damage == "missing":
                    path.unlink()
                else:
                    lines = path.read_text().splitlines()
                    if damage == "extent": lines[0] = "74"
                    elif damage == "nonfinite": lines[1] = "nan"
                    elif damage == "order": lines[2] = lines[1]
                    else: lines.append("1")
                    path.write_text("\n".join(lines) + "\n")
                result = subprocess.run([os.environ["NDROP_PHYSICAL_TEST_BINARY"],
                    "--check-physical-column-fixtures", "--physical-grid-dir", str(directory)],
                    capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("ERROR:", result.stderr)


if __name__ == "__main__":
    unittest.main()
