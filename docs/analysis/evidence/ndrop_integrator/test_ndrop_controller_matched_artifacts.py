"""Synthetic closed-contract tests for controller-matched artifacts."""

import copy
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import ndrop_controller_matched_artifacts as artifacts


def _base_row(fixture, method, kind, steps="", controlled=False):
    row = dict.fromkeys(artifacts.accuracy_header(controlled), "")
    if controlled:
        levels, mixing, interval, production = artifacts.CONTROLLED_FIXTURES[fixture]
        family, requested = "controlled_mixing", ""
        row["activation_mixing_scale"] = "0.05"
    elif fixture in artifacts.CONTINUITY_FIXTURES:
        levels, mixing, interval, production = \
            artifacts.CONTINUITY_FIXTURES[fixture]
        family, requested = "continuity", ""
    else:
        levels, production = artifacts.STRATIFIED_FIXTURES[fixture]
        family, requested = "controller_stratified", str(production)
        mixing = artifacts.stratified_mixing_scale(production)
        interval = artifacts.STRATIFIED_INTERVAL
    row.update(
        artifact_schema=artifacts.study_schema("accuracy", controlled),
        record_kind=kind, fixture=str(fixture), fixture_family=family,
        requested_nsubmix=requested, levels=str(levels),
        mixing_scale=str(mixing), interval_seconds=str(interval),
        mixing_rate_t=str(0.1 * production), activation_rate_t="0.5",
        method=method, steps=str(steps), production_nsubmix=str(production),
        reference_status="certified", reference_reason="",
        reference_disagreement="0", clipping_count="0",
        material_negative_count="0", normalized_residual="0",
        normalized_ledger_error="0", state_checksum=str(1000 + fixture),
    )
    return row


def _success(row, status, error, target=0.01):
    vector = ";".join([str(error)] * 51)
    row.update(status=status, failure_reason="", fe_target_error=str(target),
               error_normal=str(error), error_tight=str(error),
               error_max=str(error), normal_weighted_l1_by_field=vector,
               normal_linf_by_field=vector,
               tight_weighted_l1_by_field=vector,
               tight_linf_by_field=vector)
    return row


def _failed(row, status="numerical_failure", reason="synthetic failure",
            target=0.01):
    row.update(status=status, failure_reason=reason,
               fe_target_error="" if target is None else str(target))
    return row


def _terminal(fixture, method, status="selected", steps="2",
              reason="", target=0.01, controlled=False):
    row = _base_row(fixture, method, "search_terminal", steps, controlled)
    for name in artifacts.NUMERICAL_COLUMNS:
        row[name] = ""
    row.update(status=status, failure_reason=reason,
               fe_target_error="" if target is None else str(target))
    return row


def accuracy_rows(confirm=True, controlled=False):
    rows = []
    for fixture in (artifacts.CONTROLLED_FIXTURES if controlled else artifacts.FIXTURES):
        production = (artifacts.CONTROLLED_FIXTURES[fixture][3] if controlled else
                      artifacts.CONTINUITY_FIXTURES[fixture][3]
                      if fixture in artifacts.CONTINUITY_FIXTURES else
                      artifacts.STRATIFIED_FIXTURES[fixture][1])
        rows.append(_success(
            _base_row(fixture, "forward_euler", "measurement", production, controlled),
            "fe_baseline", 0.01))
        if confirm:
            for factor in (2, 4):
                rows.append(_success(
                    _base_row(fixture, "forward_euler", "confirmation", factor * production, controlled),
                    artifacts.CONFIRMATION_SUCCESS, 0.01 / factor))
        for method in artifacts.IMPLICIT_METHODS:
            rows.append(_success(_base_row(fixture, method, "measurement", 1, controlled),
                                 "eligible_above_target", 0.02))
            rows.append(_success(_base_row(fixture, method, "measurement", 2, controlled),
                                 "eligible_match", 0.005))
            rows.append(_terminal(fixture, method, controlled=controlled))
            if confirm:
                for factor in (2, 4):
                    rows.append(_success(
                        _base_row(fixture, method, "confirmation", 2 * factor, controlled),
                        artifacts.CONFIRMATION_SUCCESS, 0.005 / factor))
    return rows


def _replace_search(rows, fixture, method, candidates, terminal):
    rows[:] = [row for row in rows if not (
        row["fixture"] == str(fixture) and row["method"] == method and
        row["record_kind"] in {"measurement", "confirmation", "search_terminal"})]
    rows.extend(candidates)
    rows.append(terminal)


def _make_fixture_presearch_failure(rows, fixture, terminal_status,
                                    fe_status, reason):
    for row in rows:
        if row["fixture"] != str(fixture):
            continue
        if row["method"] in artifacts.IMPLICIT_METHODS:
            continue
        if row["record_kind"] == "measurement":
            row["status"] = fe_status
            row["failure_reason"] = reason
            row["fe_target_error"] = ""
            for name in artifacts.ERROR_COLUMNS:
                row[name] = ""
            if terminal_status == "reference_unresolved":
                row["reference_status"] = "unresolved"
                row["reference_reason"] = reason
                row["reference_disagreement"] = ""
    rows[:] = [row for row in rows if not (
        row["fixture"] == str(fixture) and
        (row["method"] in artifacts.IMPLICIT_METHODS or
         row["record_kind"] == "confirmation"))]
    for method in artifacts.IMPLICIT_METHODS:
        terminal = _terminal(fixture, method, terminal_status, "0", reason, None)
        if terminal_status == "reference_unresolved":
            terminal["reference_status"] = "unresolved"
            terminal["reference_reason"] = reason
            terminal["reference_disagreement"] = ""
        rows.append(terminal)


class ControllerMatchedArtifactsTest(unittest.TestCase):
    def test_controlled_closed_policy_and_activation_corruption(self):
        identity = {"artifact_schema": artifacts.study_schema("identity", True),
                    "fixture_policy": "controlled_mixing_v1",
                    "activation_mixing_scale": 0.05, "fixtures": [24, 25, 26, 27]}
        fixtures = artifacts.identity_fixtures(identity)
        rows = accuracy_rows(controlled=True)
        indexed = artifacts.validate_accuracy(rows, fixtures, True)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", True)
        admitted = artifacts.expected_admission(rows, indexed, frozen, "c", True)
        self.assertEqual(len(frozen["selected"]), 12)
        self.assertEqual(admitted["unique_configs"], frozen["unique_configs"])
        self.assertEqual(indexed[(24, "forward_euler", 18)]["steps"], "18")
        self.assertEqual(indexed[(25, "forward_euler", 176)]["steps"], "176")
        for values in ([24], [24, 25, 26], [24, 25, 26, 27, 0], [25, 24, 26, 27]):
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.identity_fixtures({**identity, "fixtures": values})
        for update in ({"activation_mixing_scale": 0.5}, {"fixture_policy": "other"},
                       {"artifact_schema": artifacts.IDENTITY_SCHEMA}):
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.identity_fixtures({**identity, **update})
        for field, value in (("activation_mixing_scale", "0.5"),
                             ("mixing_scale", "0.1"), ("interval_seconds", "50"),
                             ("production_nsubmix", "256"), ("requested_nsubmix", "18"),
                             ("artifact_schema", artifacts.ACCURACY_SCHEMA),
                             ("fixture_family", "continuity"), ("fixture", "0")):
            bad = copy.deepcopy(rows)
            bad[0][field] = value
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_accuracy(bad, fixtures, True)
        missing = copy.deepcopy(rows)
        del missing[0]["activation_mixing_scale"]
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(missing, fixtures, True)
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(rows)
        bad = copy.deepcopy(rows)
        for row in bad:
            if row["fixture"] == "25":
                row["activation_rate_t"] = "5"
        with self.assertRaisesRegex(artifacts.ArtifactError, "stiffness differs"):
            artifacts.validate_accuracy(bad, fixtures, True)
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy([row for row in rows if row["fixture"] != "27"],
                                        fixtures, True)

    def test_controlled_refinement_failure_withholds_not_reselects(self):
        fixtures = tuple(artifacts.CONTROLLED_FIXTURES)
        rows = accuracy_rows(controlled=True)
        indexed = artifacts.validate_accuracy(rows, fixtures, True)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", True)
        for mode in ("missing", "failed", "nondecreasing"):
            bad = copy.deepcopy(rows)
            row = next(row for row in bad if row["fixture"] == "24" and
                       row["method"] == "imex" and row["record_kind"] == "confirmation")
            if mode == "missing":
                bad.remove(row)
            elif mode == "failed":
                for name in artifacts.ERROR_COLUMNS:
                    row[name] = ""
                _failed(row)
            else:
                _success(row, artifacts.CONFIRMATION_SUCCESS, 0.02)
            new_index = artifacts.validate_accuracy(bad, fixtures, True)
            self.assertEqual(frozen, artifacts.expected_selection_from_rows(
                bad, new_index, "a", "b", True))
            admission = artifacts.expected_admission(bad, new_index, frozen, "c", True)
            self.assertNotIn([24, "imex", 2], admission["unique_configs"])

    def test_controlled_artifact_chain_with_fabricated_samples_only(self):
        """Serialization fixture: no executable, measurement, or production authority."""
        source_root = Path(__file__).resolve().parents[4]
        source_dir = source_root / "docs/analysis/evidence/ndrop_integrator"
        with tempfile.TemporaryDirectory(prefix="ndrop-controlled-fixture-") as folder:
            directory = Path(folder)
            (directory / "inputs").mkdir()
            binary = directory / "inputs/benchmark_ndrop_matched_error"
            binary.write_bytes(b"NONEXECUTABLE synthetic serialization fixture")
            identity = {"artifact_schema": artifacts.study_schema("identity", True),
                        "fixture_policy": "controlled_mixing_v1", "activation_mixing_scale": 0.05,
                        "fixtures": list(artifacts.CONTROLLED_FIXTURES),
                        "binary_sha256": artifacts.sha256(binary),
                        "source_head": "synthetic_fixture_not_a_source_commit",
                        "source_root": str(source_root), "source_status": "exploratory_source",
                        "source_attached_head": True,
                        "source_files_sha256": {name: artifacts.sha256(source_dir / name)
                                                for name in artifacts.SOURCES},
                        "build": {"fixture_only": True},
                        "design_sha256": artifacts.sha256(source_root / "docs/design/eamxx_mam4xx_ndrop.md"),
                        "python_version": artifacts.platform.python_version(),
                        "platform": "synthetic_serialization_fixture",
                        "accuracy_command": [], "accuracy_protocol": {}, "timing_protocol": {}}
            (directory / "identity.json").write_bytes(artifacts.canonical_json(identity))
            rows = accuracy_rows(controlled=True)
            with (directory / "accuracy.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=artifacts.CONTROLLED_ACCURACY_HEADER)
                writer.writeheader()
                writer.writerows(rows)
            artifacts._completion(directory, "accuracy", True)
            args = SimpleNamespace(output_dir=directory)
            artifacts.freeze(args)
            frozen = artifacts.selection(directory)
            admitted = artifacts.admission(directory, frozen)
            indexed = artifacts.validate_accuracy(rows, tuple(artifacts.CONTROLLED_FIXTURES), True)
            fake_samples = [{"artifact_schema": artifacts.study_schema("timing", True),
                             "fixture": str(f), "method": m, "steps": str(n), "block": str(block),
                             "seconds_per_column": "1", "state_checksum": indexed[(f, m, n)]["state_checksum"],
                             "stop_reason": "median_precision"}
                            for block in range(30) for f, m, n in admitted["unique_configs"]]
            artifacts.validate_timing(fake_samples, admitted, indexed, True)
            corrupt = copy.deepcopy(fake_samples)
            corrupt[0]["state_checksum"] = "0"
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_timing(corrupt, admitted, indexed, True)
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_timing(fake_samples[:-1], admitted, indexed, True)
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_timing(fake_samples, admitted, indexed)
            with (directory / "timing.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=artifacts.TIMING_HEADER)
                writer.writeheader()
                writer.writerows(fake_samples)
            artifacts._completion(directory, "timing", True)
            artifacts.report(args)
            artifacts.validate_final(args)
            report = artifacts.load_json(directory / "analysis.json")
            self.assertEqual(report["artifact_schema"], artifacts.study_schema("analysis", True))
            self.assertTrue(all(float(row["activation_mixing_scale"]) == 0.05
                                for row in report["fixtures"]))
            self.assertEqual(len(report["comparisons"]), 12)
            self.assertTrue(all(row["numerator_steps"] > 0 for row in report["comparisons"]))
            stale = artifacts.load_json(directory / "timing.complete.json")
            stale["artifact_schema"] = artifacts.COMPLETION_SCHEMA
            (directory / "timing.complete.json").write_bytes(artifacts.canonical_json(stale))
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_final(args)
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.atomic_new(directory / "identity.json", b"overwrite")

    def test_subset_is_frozen_and_coverage_is_exact(self):
        fixtures = (4, 5, 10, 11, 18, 19, 20, 21)
        rows = [row for row in accuracy_rows() if int(row["fixture"]) in fixtures]
        indexed = artifacts.validate_accuracy(rows, fixtures)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b")
        self.assertEqual({entry["fixture"] for entry in frozen["selected"]}, set(fixtures))
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(rows)
        for values in ([], [4, 4], [5, 4], [24], [True]):
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.identity_fixtures({"fixtures": values})

    def test_refinement_admission_does_not_change_error_selection(self):
        rows = accuracy_rows()
        indexed = artifacts.validate_accuracy(rows)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b")
        admitted = artifacts.expected_admission(rows, indexed, frozen, "c")
        self.assertEqual(admitted["unique_configs"], frozen["unique_configs"])
        for mode in ("missing", "failure", "nondecreasing"):
            bad = copy.deepcopy(rows)
            row = next(row for row in bad if row["fixture"] == "0" and
                       row["method"] == "forward_euler" and
                       row["record_kind"] == "confirmation")
            if mode == "missing":
                bad.remove(row)
            elif mode == "failure":
                for name in artifacts.ERROR_COLUMNS:
                    row[name] = ""
                _failed(row)
            else:
                _success(row, artifacts.CONFIRMATION_SUCCESS, 0.02)
            new_index = artifacts.validate_accuracy(bad)
            self.assertEqual(frozen, artifacts.expected_selection_from_rows(
                bad, new_index, "a", "b"))
            admitted = artifacts.expected_admission(bad, new_index, frozen, "c")
            entry = next(entry for entry in admitted["entries"] if
                         entry["fixture"] == 0 and entry["method"] == "forward_euler")
            self.assertEqual(entry["status"], "excluded")
            self.assertNotIn([0, "forward_euler", 6], admitted["unique_configs"])

    def test_complete_selection_consecutive_search_and_non_dyadic_fe(self):
        rows = accuracy_rows()
        indexed = artifacts.validate_accuracy(rows)
        self.assertIn((0, "forward_euler", 6), indexed)
        self.assertIn((1, "forward_euler", 176), indexed)
        frozen = artifacts.expected_selection_from_rows(
            rows, indexed, "a" * 64, "b" * 64)
        self.assertEqual(len(frozen["selected"]), 72)
        self.assertEqual(len(frozen["unique_configs"]), 72)
        imex = next(entry for entry in frozen["selected"]
                    if entry["fixture"] == 0 and entry["method"] == "imex")
        self.assertEqual(imex["steps"], 2)
        self.assertEqual(imex["achieved_error_to_fe_target_ratio"], 0.5)
        self.assertEqual(indexed[(0, "imex", 1)]["status"],
                         "eligible_above_target")

        nonconsecutive = copy.deepcopy(rows)
        candidate = next(row for row in nonconsecutive
                         if row["fixture"] == "0" and row["method"] == "imex"
                         and row["record_kind"] == "measurement"
                         and row["steps"] == "2")
        candidate["steps"] = "3"
        with self.assertRaisesRegex(artifacts.ArtifactError, "not consecutive"):
            artifacts.validate_accuracy(nonconsecutive)

        premature = copy.deepcopy(rows)
        first = next(row for row in premature
                     if row["fixture"] == "0" and row["method"] == "imex"
                     and row["record_kind"] == "measurement"
                     and row["steps"] == "1")
        first["status"] = "eligible_match"
        first["error_normal"] = first["error_tight"] = first["error_max"] = "0.005"
        for name in ("normal_weighted_l1_by_field", "normal_linf_by_field",
                     "tight_weighted_l1_by_field", "tight_linf_by_field"):
            first[name] = ";".join(["0.005"] * 51)
        with self.assertRaisesRegex(artifacts.ArtifactError, "first match"):
            artifacts.validate_accuracy(premature)

    def test_selected_counts_below_equal_and_above_fe_are_admitted(self):
        rows = accuracy_rows(confirm=False)
        fixture = 0  # Exact FE baseline count is six.
        for method, selected in (("imex", 6), ("backward_euler", 7)):
            candidates = []
            for step in range(1, selected + 1):
                status = ("eligible_match" if step == selected
                          else "eligible_above_target")
                error = 0.005 if step == selected else 0.02
                candidates.append(_success(
                    _base_row(fixture, method, "measurement", step),
                    status, error))
            _replace_search(
                rows, fixture, method, candidates,
                _terminal(fixture, method, "selected", str(selected)))
        indexed = artifacts.validate_accuracy(rows)
        frozen = artifacts.expected_selection_from_rows(
            rows, indexed, "a" * 64, "b" * 64)
        selected = {(entry["method"], entry["steps"])
                    for entry in frozen["selected"]
                    if entry["fixture"] == fixture}
        self.assertIn(("forward_euler", 6), selected)
        self.assertIn(("imex", 6), selected)
        self.assertIn(("backward_euler", 7), selected)

    def test_corrupt_mixed_missing_reason_and_controller_mismatch(self):
        rows = accuracy_rows()
        corrupt = copy.deepcopy(rows)
        success = next(row for row in corrupt if row["record_kind"] == "measurement")
        success["normal_linf_by_field"] = "1"
        with self.assertRaisesRegex(artifacts.ArtifactError, "field error vector"):
            artifacts.validate_accuracy(corrupt)

        mixed = copy.deepcopy(rows)
        mixed[0]["artifact_schema"] = "ndrop_matched_accuracy_v1"
        with self.assertRaisesRegex(artifacts.ArtifactError, "mixed"):
            artifacts.validate_accuracy(mixed)

        wrong_stratified_identity = copy.deepcopy(rows)
        for row in wrong_stratified_identity:
            if row["fixture"] == "4":
                row["interval_seconds"] = "50.0"
        with self.assertRaisesRegex(artifacts.ArtifactError,
                                    "stratified fixture identity"):
            artifacts.validate_accuracy(wrong_stratified_identity)

        missing = copy.deepcopy(rows)
        failed = _failed(_base_row(0, "imex", "measurement", 1), reason="")
        _replace_search(missing, 0, "imex", [failed],
                        _terminal(0, "imex", "accuracy_budget_exhausted", "1",
                                  "budget"))
        with self.assertRaisesRegex(artifacts.ArtifactError, "lacks its first reason"):
            artifacts.validate_accuracy(missing)

        mismatch = accuracy_rows()
        fixture = 4
        fe = next(row for row in mismatch if row["fixture"] == str(fixture)
                  and row["method"] == "forward_euler")
        fe["production_nsubmix"] = "2"
        fe["steps"] = "2"
        _make_fixture_presearch_failure(
            mismatch, fixture, "fixture_controller_mismatch",
            "fixture_controller_mismatch", "controller returned 2, requested 1")
        # The helper rebuilds terminals from requested metadata; set the actual
        # production result consistently across every record in the fixture.
        for row in mismatch:
            if row["fixture"] == str(fixture):
                row["production_nsubmix"] = "2"
        indexed = artifacts.validate_accuracy(mismatch)
        frozen = artifacts.expected_selection_from_rows(
            mismatch, indexed, "a" * 64, "b" * 64)
        self.assertEqual(next(entry["status"] for entry in frozen["selected"]
                              if entry["fixture"] == fixture and
                              entry["method"] == "imex"),
                         "fixture_controller_mismatch")

    def test_budget_reference_and_fe_terminal_paths(self):
        budget = accuracy_rows()
        candidate = _failed(_base_row(0, "imex", "measurement", 1))
        _replace_search(
            budget, 0, "imex", [candidate],
            _terminal(0, "imex", "accuracy_budget_exhausted", "1",
                      "accuracy budget exhausted"))
        artifacts.validate_accuracy(budget)

        unresolved = accuracy_rows()
        _make_fixture_presearch_failure(
            unresolved, 0, "reference_unresolved", "reference_unresolved",
            "reference disagreement")
        artifacts.validate_accuracy(unresolved)

        ineligible = accuracy_rows()
        _make_fixture_presearch_failure(
            ineligible, 0, "fe_baseline_ineligible", "numerical_failure",
            "FE residual gate")
        artifacts.validate_accuracy(ineligible)

    def test_cap_terminal_requires_every_integer_count(self):
        rows = accuracy_rows(confirm=False)
        candidates = [
            _failed(_base_row(0, "imex", "measurement", step))
            for step in range(1, artifacts.MAX_STEPS + 1)
        ]
        _replace_search(
            rows, 0, "imex", candidates,
            _terminal(0, "imex", "not_attained_within_cap",
                      str(artifacts.MAX_STEPS), "cap exhausted"))
        artifacts.validate_accuracy(rows)

        over_cap = accuracy_rows(confirm=False)
        candidates = [_success(
            _base_row(0, "imex", "measurement", artifacts.MAX_STEPS + 1),
            "eligible_match", 0.005)]
        _replace_search(
            over_cap, 0, "imex", candidates,
            _terminal(0, "imex", "selected",
                      str(artifacts.MAX_STEPS + 1)))
        with self.assertRaisesRegex(artifacts.ArtifactError, "count cap"):
            artifacts.validate_accuracy(over_cap)
        broken = rows[:-2] + [rows[-1]]
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(broken)

    def test_timing_checksum_pairing_and_bootstrap(self):
        accuracy_rows_value = accuracy_rows(confirm=False)
        indexed = artifacts.validate_accuracy(accuracy_rows_value)
        frozen = artifacts.expected_selection_from_rows(
            accuracy_rows_value, indexed, "a" * 64, "b" * 64)
        rows = []
        for block in range(30):
            for fixture, method, steps in frozen["unique_configs"]:
                rows.append({"artifact_schema": artifacts.TIMING_SCHEMA,
                             "fixture": str(fixture), "method": method,
                             "steps": str(steps), "block": str(block),
                             "seconds_per_column": "2" if method == "forward_euler"
                                                   else "1",
                             "state_checksum": indexed[(fixture, method, steps)][
                                 "state_checksum"],
                             "stop_reason": "median_precision"})
        artifacts.validate_timing(rows, frozen, indexed)
        ratio = artifacts.paired_ratio([2.] * 30, [1.] * 30,
                                        artifacts.BOOTSTRAP_SEED)
        self.assertEqual((ratio["ratio"], ratio["lower"], ratio["upper"]),
                         (2., 2., 2.))

        stale = copy.deepcopy(rows)
        stale[-1]["state_checksum"] = "999999"
        with self.assertRaisesRegex(artifacts.ArtifactError, "timed state differs"):
            artifacts.validate_timing(stale, frozen, indexed)
        missing_pair = rows[:-1]
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_timing(missing_pair, frozen, indexed)

    def test_completion_requires_current_identity_and_hash(self):
        with tempfile.TemporaryDirectory(prefix="ndrop-controller-artifact-") as folder:
            directory = Path(folder)
            (directory / "identity.json").write_text("{}\n", encoding="utf-8")
            (directory / "accuracy.csv").write_text("payload\n", encoding="utf-8")
            with self.assertRaises(FileNotFoundError):
                artifacts.completed(directory, "accuracy")
            record = {"artifact_schema": artifacts.COMPLETION_SCHEMA,
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

    def test_identity_binds_analyzer_python_and_source_files(self):
        source_root = Path(__file__).resolve().parents[4]
        source_dir = source_root / "docs/analysis/evidence/ndrop_integrator"
        with tempfile.TemporaryDirectory(prefix="ndrop-controller-identity-") as folder:
            directory = Path(folder)
            (directory / "inputs").mkdir()
            binary = directory / "inputs/benchmark_ndrop_matched_error"
            binary.write_bytes(b"synthetic binary")
            identity = {
                "artifact_schema": artifacts.IDENTITY_SCHEMA,
                "binary_sha256": artifacts.sha256(binary),
                "python_version": artifacts.platform.python_version(),
                "source_root": str(source_root),
                "source_files_sha256": {
                    name: artifacts.sha256(source_dir / name)
                    for name in artifacts.SOURCES
                },
                "design_sha256": artifacts.sha256(
                    source_root / "docs/design/eamxx_mam4xx_ndrop.md"),
            }
            (directory / "identity.json").write_bytes(
                artifacts.canonical_json(identity))
            artifacts.validate_identity(directory)

            stale = copy.deepcopy(identity)
            stale["source_files_sha256"]["ndrop_controller_matched_artifacts.py"] = \
                "0" * 64
            (directory / "identity.json").write_bytes(
                artifacts.canonical_json(stale))
            with self.assertRaisesRegex(artifacts.ArtifactError,
                                        "source file differs"):
                artifacts.validate_identity(directory)


def physical_controller_rows():
    """Hand-defined serialization fixtures, not simulated production states."""
    rows = []
    for fixture, levels, strength, mixing, velocity, count, stiffness in (
            (28, 72, "moderate", 3., .3, 1, .45),
            (29, 72, "strong", 30., 1., 5, 3.6),
            (30, 128, "moderate", 3., .3, 3, 1.8),
            (31, 128, "strong", 30., 1., 10, 8.1)):
        def row(method, kind, steps):
            result = dict.fromkeys(artifacts.accuracy_header(physical=True), "")
            result.update(artifact_schema="ndrop_physical_controller_accuracy_v1",
                          fixture=str(fixture), fixture_family="physical_column",
                          record_kind=kind, levels=str(levels), mixing_scale=str(mixing),
                          interval_seconds="150", mixing_rate_t=str(stiffness),
                          activation_rate_t="0.05", production_nsubmix=str(count),
                          method=method, steps=str(steps), reference_status="certified",
                          reference_disagreement="0", clipping_count="0",
                          material_negative_count="0", normalized_residual="0",
                          normalized_ledger_error="0", state_checksum=str(1000 + fixture),
                          turbulence=strength, cloud_base_wtke_m_s=str(velocity),
                          controller_rate_t=str(stiffness), cloud_base_level="20",
                          cloud_base_spacing_m="50", limiting_level="21",
                          limiting_spacing_m="40", activation_limiter_count="0")
            return result
        rows.append(_success(row("forward_euler", "measurement", count),
                             "fe_baseline", .03, .03))
        for factor in (2, 4):
            rows.append(_success(row("forward_euler", "confirmation", count * factor),
                                 "eligible_confirmation", .03 / factor, .03))
        for method in artifacts.IMPLICIT_METHODS:
            rows.append(_success(row(method, "measurement", 1),
                                 "eligible_above_target", .04, .03))
            rows.append(_success(row(method, "measurement", 2),
                                 "eligible_match", .025, .03))
            terminal = row(method, "search_terminal", 2)
            for name in artifacts.NUMERICAL_COLUMNS + ("activation_limiter_count",):
                terminal[name] = ""
            terminal.update(status="selected", fe_target_error="0.03")
            rows.append(terminal)
            for factor in (2, 4):
                rows.append(_success(row(method, "confirmation", 2 * factor),
                                     "eligible_confirmation", .025 / factor, .03))
    return rows


def geometry_controller_rows():
    """Fabricated serialization outcomes only; not executed accuracy evidence."""
    templates = physical_controller_rows()
    rows = []
    for fixture, (levels, strength, mixing, velocity) in artifacts.GEOMETRY_COLUMNS.items():
        parent, subdivision = artifacts.GEOMETRY_INPUTS[fixture]
        template_id = (28 if parent == 72 else 30) + fixture % 2
        count = {28: 1, 29: 5, 30: 3, 31: 10}.get(fixture, 7 if fixture % 2 else 2)
        old_count = {28: 1, 29: 5, 30: 3, 31: 10}[template_id]
        for old in templates:
            if int(old["fixture"]) != template_id:
                continue
            row = copy.deepcopy(old)
            row.update(artifact_schema=artifacts.study_schema("accuracy", physical=artifacts.GEOMETRY_POLICY),
                       fixture=str(fixture), fixture_family="physical_geometry", levels=str(levels),
                       turbulence=strength, mixing_scale=str(mixing), cloud_base_wtke_m_s=str(velocity),
                       production_nsubmix=str(count), controller_rate_t=str(.9 * (count - .5)),
                       parent_grid_levels=str(parent), geometry_subdivision=str(subdivision),
                       refinement_lower_m="400", refinement_upper_m="1000", limiting_height_m="500",
                       limiting_above_diffusivity_m2_s=str(mixing),
                       limiting_below_diffusivity_m2_s=str(mixing), limiting_cloud_base_faces="0")
            if row["method"] == "forward_euler":
                row["steps"] = str(count * (int(row["steps"]) // old_count))
            rows.append(row)
    return rows


class GeometryControllerMatchedTest(unittest.TestCase):
    policy = artifacts.GEOMETRY_POLICY
    fixtures = tuple(artifacts.GEOMETRY_COLUMNS)

    def test_exact_policy_and_selection(self):
        rows = geometry_controller_rows()
        indexed = artifacts.validate_accuracy(rows, self.fixtures, physical=self.policy)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=self.policy)
        self.assertEqual(len(frozen["selected"]), 36)
        self.assertEqual(len(artifacts.expected_admission(
            rows, indexed, frozen, "c", physical=self.policy)["unique_configs"]), 36)
        self.assertEqual({int(r["levels"]) for r in rows}, {72, 79, 93, 128, 141, 167})
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(rows, self.fixtures, physical=True)
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(rows, self.fixtures[:-1], physical=self.policy)
        with self.assertRaisesRegex(artifacts.ArtifactError, "missing selected"):
            artifacts.validate_timing([], frozen, indexed, physical=self.policy)

    def test_paired_geometry_timing_and_ratio_serialization(self):
        """Fabricated time samples only; no execution or approval evidence."""
        rows = geometry_controller_rows()
        indexed = artifacts.validate_accuracy(rows, self.fixtures, physical=self.policy)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=self.policy)
        samples = [{"artifact_schema": artifacts.study_schema("timing", physical=self.policy),
                    "fixture": str(f), "method": m, "steps": str(n), "block": str(b),
                    "seconds_per_column": str({"forward_euler": 2, "imex": 1, "backward_euler": 3}[m]),
                    "state_checksum": indexed[(f, m, n)]["state_checksum"], "stop_reason": "median_precision"}
                   for b in range(30) for f, m, n in frozen["unique_configs"]]
        artifacts.validate_timing(samples, frozen, indexed, physical=self.policy)
        ratio = artifacts.paired_ratio([2.] * 30, [1.] * 30, artifacts.BOOTSTRAP_SEED)
        self.assertEqual((ratio["ratio"], ratio["lower"], ratio["upper"]), (2., 2., 2.))
        for bad in (samples[:-1], samples + [samples[0]],
                    [dict(samples[0], state_checksum="0"), *samples[1:]],
                    [dict(samples[0], artifact_schema=artifacts.TIMING_SCHEMA), *samples[1:]]):
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_timing(bad, frozen, indexed, physical=self.policy)
        command = artifacts.timing_command(Path("/fixture"), frozen, physical=self.policy)
        self.assertEqual(command[1], "--physical-geometry-controller-timing")
        self.assertEqual(command[command.index("--warmups") + 1], "10")
        self.assertEqual(command[command.index("--max-blocks") + 1], "100")
        self.assertEqual(artifacts.evidence_scope(physical=self.policy),
                         "physical_geometry_controller_idealization_unweighted")

    def test_failed_process_and_promotion_preserve_bounded_partial(self):
        """Mocked process/authority plumbing; no approval record or real run."""
        rows = geometry_controller_rows()
        indexed = artifacts.validate_accuracy(rows, self.fixtures, physical=self.policy)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=self.policy)
        admitted = artifacts.expected_admission(rows, indexed, frozen, "c", physical=self.policy)
        identity = {"artifact_schema": artifacts.study_schema("identity", physical=self.policy),
                    "fixture_policy": artifacts.fixture_policy(self.policy),
                    "comparison_policy": "natural_fe_controller_error_v1",
                    "fixtures": list(self.fixtures), "source_status": "exploratory_source"}
        for failure in ("process", "flush", "rename"):
            with tempfile.TemporaryDirectory(prefix="ndrop-geometry-transaction-") as folder:
                directory = Path(folder)
                payload = directory / "accuracy.csv"
                with payload.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, artifacts.accuracy_header(physical=self.policy))
                    writer.writeheader(); writer.writerows(rows)
                def fake_process(*args, **kwargs):
                    if failure == "process":
                        raise subprocess.CalledProcessError(1, args[0])
                    samples = [{"artifact_schema": artifacts.study_schema("timing", physical=self.policy),
                                "fixture": str(f), "method": m, "steps": str(n), "block": str(b),
                                "seconds_per_column": "1", "state_checksum": indexed[(f, m, n)]["state_checksum"],
                                "stop_reason": "median_precision"}
                               for b in range(30) for f, m, n in admitted["unique_configs"]]
                    import io
                    text = io.StringIO()
                    writer = csv.DictWriter(text, artifacts.TIMING_HEADER)
                    writer.writeheader(); writer.writerows(samples)
                    kwargs["stdout"].write(text.getvalue().encode())
                    return SimpleNamespace(returncode=0)
                with mock.patch.object(artifacts, "validate_identity", return_value=identity), \
                     mock.patch.object(artifacts, "selection", return_value=frozen), \
                     mock.patch.object(artifacts, "admission", return_value=admitted), \
                     mock.patch.object(artifacts, "geometry_timing_approval"), \
                     mock.patch.object(artifacts, "completed", return_value=payload), \
                     mock.patch.object(artifacts.subprocess, "run", side_effect=fake_process) as process, \
                     mock.patch.object(artifacts.os, "fsync", side_effect=OSError("synthetic flush failure")
                                       if failure == "flush" else None), \
                     mock.patch.object(artifacts.os, "replace", side_effect=OSError("synthetic rename failure")
                                       if failure == "rename" else None):
                    with self.assertRaises((OSError, subprocess.CalledProcessError)):
                        artifacts.run_timing(SimpleNamespace(output_dir=directory))
                    self.assertTrue((directory / "timing.csv.partial").exists())
                    self.assertFalse((directory / "timing.complete.json").exists())
                    self.assertFalse((directory / "timing_execution.json").exists())
                    with self.assertRaisesRegex(artifacts.ArtifactError, "already exists"):
                        artifacts.run_timing(SimpleNamespace(output_dir=directory))
                    self.assertEqual(process.call_count, 1)

    @unittest.skipUnless(os.environ.get("NDROP_GEOMETRY_QUALIFIED_DIR"), "exact requalified untimed artifact required")
    def test_real_qualified_approval_corruption_without_execution(self):
        directory = Path(os.environ["NDROP_GEOMETRY_QUALIFIED_DIR"])
        frozen = artifacts.selection(directory)
        admitted = artifacts.admission(directory, frozen)
        expected = artifacts.expected_geometry_timing_approval(directory, frozen, admitted)
        # Validate in-memory records only. Never write or issue an approval.
        with mock.patch.object(artifacts, "load_json", return_value=expected):
            artifacts.geometry_timing_approval(directory, frozen, admitted)
        for field in expected:
            bad = copy.deepcopy(expected)
            bad[field] = [] if field == "command" else "corrupt"
            with self.subTest(field=field), mock.patch.object(artifacts, "load_json", return_value=bad), \
                 self.assertRaises(artifacts.ArtifactError):
                artifacts.geometry_timing_approval(directory, frozen, admitted)

    def test_execution_counts_types_and_hashes(self):
        """Hash/count grammar fixture only, not a claimed executed process."""
        with tempfile.TemporaryDirectory(prefix="ndrop-geometry-counts-") as folder:
            directory = Path(folder)
            for name in ("identity.json", "selection.json", "timing_approval.json", "timing.csv"):
                (directory / name).write_bytes(b"FABRICATED SERIALIZATION FIXTURE")
            admitted = {"unique_configs": [[28, "forward_euler", 1], [28, "imex", 1]]}
            rows = [dict(fabricated=True)] * 60
            expected = artifacts.expected_geometry_execution(directory, admitted, rows)
            self.assertEqual(expected["timing_rows"], 60)
            self.assertEqual([v["warmups"] for v in expected["warmups_by_configuration"]], [10, 10])
            with mock.patch.object(artifacts, "load_json", return_value=expected):
                artifacts.validate_geometry_execution(directory, admitted, rows)
            for field, value in (("process_exit_code", False), ("process_exit_code", 0.0),
                                 ("timing_rows", 59), ("warmups_by_configuration", []),
                                 ("identity_sha256", "0" * 64), ("command", [])):
                bad = {**expected, field: value}
                with self.subTest(field=field, value=value), \
                     mock.patch.object(artifacts, "load_json", return_value=bad), \
                     self.assertRaises(artifacts.ArtifactError):
                    artifacts.validate_geometry_execution(directory, admitted, rows)

    def test_closed_geometry_and_corruption(self):
        rows = geometry_controller_rows()
        for name, value in (("parent_grid_levels", "128"), ("geometry_subdivision", "4"),
                            ("levels", "140"), ("refinement_lower_m", "399"),
                            ("refinement_upper_m", "1001"), ("limiting_height_m", "nan"),
                            ("limiting_above_diffusivity_m2_s", "-1"),
                            ("limiting_cloud_base_faces", "3"), ("steps", "256"),
                            ("fe_target_error", "0.01"), ("fixture_family", "physical_column"),
                            ("artifact_schema", "ndrop_physical_controller_accuracy_v1")):
            bad = copy.deepcopy(rows)
            next(r for r in bad if r["fixture"] == "32")[name] = value
            with self.subTest(name=name), self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_accuracy(bad, self.fixtures, physical=self.policy)
        missing = copy.deepcopy(rows)
        del missing[0]["geometry_subdivision"]
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(missing, self.fixtures, physical=self.policy)
        stale = copy.deepcopy(rows)
        next(r for r in stale if r["fixture"] == "32" and r["record_kind"] == "confirmation")["limiting_height_m"] = "501"
        with self.assertRaisesRegex(artifacts.ArtifactError, "metadata is inconsistent"):
            artifacts.validate_accuracy(stale, self.fixtures, physical=self.policy)

    def test_lower_failures_and_missing_confirmations_are_retained(self):
        rows = geometry_controller_rows()
        low = next(r for r in rows if r["fixture"] == "39" and r["method"] == "imex"
                   and r["record_kind"] == "measurement" and r["steps"] == "1")
        for name in artifacts.ERROR_COLUMNS:
            low[name] = ""
        _failed(low, reason="fabricated material failure", target=.03)
        indexed = artifacts.validate_accuracy(rows, self.fixtures, physical=self.policy)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=self.policy)
        self.assertEqual(indexed[(39, "imex", 1)]["failure_reason"], "fabricated material failure")
        rows = [r for r in rows if not (r["fixture"] == "39" and r["method"] == "backward_euler"
                                       and r["record_kind"] == "confirmation")]
        indexed = artifacts.validate_accuracy(rows, self.fixtures, physical=self.policy)
        self.assertEqual(frozen, artifacts.expected_selection_from_rows(
            rows, indexed, "a", "b", physical=self.policy))
        self.assertNotIn([39, "backward_euler", 2], artifacts.expected_admission(
            rows, indexed, frozen, "c", physical=self.policy)["unique_configs"])
        with self.assertRaisesRegex(artifacts.ArtifactError, "not consecutive"):
            artifacts.validate_accuracy([r for r in rows if r is not low], self.fixtures, physical=self.policy)

    def test_untimed_identity_completion_and_no_execution_bypass(self):
        root = Path(__file__).resolve().parents[4]
        source = root / "docs/analysis/evidence/ndrop_integrator"
        with tempfile.TemporaryDirectory(prefix="ndrop-geometry-artifact-") as folder:
            directory = Path(folder)
            (directory / "inputs/grids").mkdir(parents=True)
            binary = directory / "inputs/benchmark_ndrop_matched_error"
            binary.write_bytes(b"NONEXECUTABLE fabricated geometry serialization fixture")
            hashes = {}
            for name in artifacts.PHYSICAL_GRIDS:
                (directory / "inputs/grids" / name).write_bytes(
                    (root / "components/homme/test/vcoord" / name).read_bytes())
                hashes[name] = artifacts.sha256(directory / "inputs/grids" / name)
            identity = {
                "artifact_schema": artifacts.study_schema("identity", physical=self.policy),
                "fixture_policy": artifacts.fixture_policy(self.policy),
                "comparison_policy": "natural_fe_controller_error_v1", "fixtures": list(self.fixtures),
                "binary_sha256": artifacts.sha256(binary), "grid_files_sha256": hashes,
                "source_head": "fabricated_not_a_commit", "source_root": str(root),
                "source_status": "exploratory_source", "source_attached_head": True,
                "source_files_sha256": {n: artifacts.sha256(source / n) for n in artifacts.source_files(self.policy)},
                "build": {"fixture_only": True, "dependencies": {
                    "ekat": {"head": "3a1485aecfed7a036e4fa30ce845498fd37975e2", "dirty": False},
                    "mam4xx": {"head": "0d076d8ebb004a18f05f51f97b08edd9def6826f", "dirty": False}}},
                "design_sha256": artifacts.sha256(root / "docs/design/eamxx_mam4xx_ndrop.md"),
                "python_version": artifacts.platform.python_version(), "platform": "fabricated_fixture",
                "accuracy_command": ["inputs/benchmark_ndrop_matched_error", "--physical-geometry-controller-accuracy",
                    "--physical-grid-dir", "inputs/grids", "--max-steps", "8192", "--budget-seconds", "7200"],
                "accuracy_protocol": {"integer_search": True, "maximum_steps": 8192,
                    "budget_seconds": 7200, "confirmation_multipliers": [2, 4]},
                "timing_protocol": {"warmups": 10, "minimum_blocks": 30, "maximum_blocks": 100,
                    "budget_seconds": 7200, "complete_interleaved_blocks": True},
            }
            path = directory / "identity.json"
            path.write_bytes(artifacts.canonical_json(identity))
            artifacts.validate_identity(directory)
            for update in ({"fixtures": [28, 29, 30, 31]}, {"fixture_policy": artifacts.fixture_policy(True)},
                           {"accuracy_command": []}, {"extra": 1}, {"grid_files_sha256": {}}):
                path.write_bytes(artifacts.canonical_json({**identity, **update}))
                with self.assertRaises(artifacts.ArtifactError):
                    artifacts.validate_identity(directory)
            path.write_bytes(artifacts.canonical_json(identity))
            rows = geometry_controller_rows()
            with (directory / "accuracy.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, artifacts.accuracy_header(physical=self.policy))
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaises(FileNotFoundError):
                artifacts.completed(directory, "accuracy", physical=self.policy)
            artifacts._completion(directory, "accuracy", physical=self.policy)
            args = SimpleNamespace(output_dir=directory)
            artifacts.freeze(args)
            frozen = artifacts.selection(directory)
            self.assertEqual(len(artifacts.admission(directory, frozen)["unique_configs"]), 36)
            for operation in (artifacts.run_timing, artifacts.report, artifacts.validate_final):
                with self.assertRaises(artifacts.ArtifactError):
                    operation(args)
            self.assertFalse((directory / "timing.csv").exists())
            bad_selection = artifacts.load_json(directory / "selection.json")
            bad_selection["unique_configs"][0][2] = 256
            (directory / "selection.json").write_bytes(artifacts.canonical_json(bad_selection))
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.selection(directory)
            (directory / "accuracy.csv").write_bytes(b"corrupt")
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.completed(directory, "accuracy", physical=self.policy)
            (directory / "inputs/grids/scream-128i.ascii").write_bytes(b"corrupt")
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_identity(directory)

    @unittest.skipUnless(os.environ.get("NDROP_GEOMETRY_TEST_BINARY"), "explicit untimed local binary required")
    def test_actual_geometry_binary_uses_emergent_counts(self):
        root = Path(__file__).resolve().parents[4]
        binary = os.environ["NDROP_GEOMETRY_TEST_BINARY"]
        process = subprocess.run([binary, "--physical-geometry-controller-accuracy", "--physical-grid-dir",
                                  str(root / "components/homme/test/vcoord"), "--max-steps", "256",
                                  "--budget-seconds", "300"], capture_output=True, text=True, timeout=330)
        self.assertEqual(process.returncode, 0, process.stderr)
        rows = list(csv.DictReader(process.stdout.splitlines()))
        indexed = artifacts.validate_accuracy(rows, self.fixtures, physical=self.policy)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=self.policy)
        for entry in frozen["selected"]:
            if entry["status"] == "selected":
                self.assertLessEqual(entry["achieved_error"], entry["fe_target_error"])
            if entry["method"] == "forward_euler":
                self.assertEqual(entry["steps"], int(indexed[(entry["fixture"], "forward_euler", entry["steps"])]["production_nsubmix"]))
        self.assertEqual(len(frozen["selected"]), 36)
        for invalid in (["--physical-geometry-controller-timing"],
                        ["--physical-geometry-controller-accuracy", "--fixtures", "32"]):
            self.assertNotEqual(subprocess.run([binary, *invalid], capture_output=True).returncode, 0)


class PhysicalControllerMatchedTest(unittest.TestCase):
    def test_failure_and_partial_search_preserve_exclusions(self):
        rows = physical_controller_rows()
        failed = next(r for r in rows if r["fixture"] == "28" and
                      r["method"] == "imex" and r["record_kind"] == "measurement")
        for name in artifacts.ERROR_COLUMNS:
            failed[name] = ""
        _failed(failed, reason="synthetic partial failure", target=.03)
        indexed = artifacts.validate_accuracy(rows, (28, 29, 30, 31), physical=True)
        self.assertEqual(indexed[(28, "imex", 1)]["failure_reason"], "synthetic partial failure")
        for status in ("reference_unresolved", "numerical_failure", "fe_clipping"):
            rows = physical_controller_rows()
            fe = rows[0]
            terminals = [dict(r) for r in rows if r["fixture"] == "28" and
                         r["record_kind"] == "search_terminal"]
            rows = [r for r in rows if r["fixture"] != "28"]
            for name in artifacts.ERROR_COLUMNS:
                fe[name] = ""
            _failed(fe, status, "synthetic first failure", None)
            if status == "fe_clipping": fe["clipping_count"] = "1"
            for terminal in terminals:
                terminal.update(status="reference_unresolved" if status == "reference_unresolved"
                                else "fe_baseline_ineligible", steps="0", fe_target_error="",
                                failure_reason="synthetic first failure")
            if status == "reference_unresolved":
                for r in [fe, *terminals]:
                    r.update(reference_status="unresolved", reference_reason="synthetic reference failure",
                             reference_disagreement="")
            rows.extend([fe, *terminals])
            indexed = artifacts.validate_accuracy(rows, (28, 29, 30, 31), physical=True)
            frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=True)
            self.assertFalse(any(config[0] == 28 for config in frozen["unique_configs"]))

    def test_physical_artifact_chain_and_grid_corruption(self):
        """NONEXECUTABLE synthetic artifact fixture; no performance measurement."""
        root = Path(__file__).resolve().parents[4]
        source = root / "docs/analysis/evidence/ndrop_integrator"
        with tempfile.TemporaryDirectory(prefix="ndrop-pctrl-artifact-") as folder:
            directory = Path(folder); (directory / "inputs/grids").mkdir(parents=True)
            binary = directory / "inputs/benchmark_ndrop_matched_error"
            binary.write_bytes(b"NONEXECUTABLE synthetic fixture")
            hashes = {}
            for name in artifacts.PHYSICAL_GRIDS:
                data = (root / "components/homme/test/vcoord" / name).read_bytes()
                (directory / "inputs/grids" / name).write_bytes(data)
                hashes[name] = artifacts.sha256(directory / "inputs/grids" / name)
            identity = {"artifact_schema": "ndrop_physical_controller_identity_v1",
                        "fixture_policy": "dycoms_inspired_isothermal_repository_grid_v1",
                        "comparison_policy": "natural_fe_controller_error_v1", "fixtures": [28, 29, 30, 31],
                        "binary_sha256": artifacts.sha256(binary), "grid_files_sha256": hashes,
                        "source_head": "synthetic_not_a_commit", "source_root": str(root),
                        "source_status": "exploratory_source", "source_attached_head": True,
                        "source_files_sha256": {n: artifacts.sha256(source / n) for n in artifacts.source_files(True)},
                        "build": {"fixture_only": True, "dependencies": {
                            "ekat": {"head": "3a1485aecfed7a036e4fa30ce845498fd37975e2", "dirty": False},
                            "mam4xx": {"head": "0d076d8ebb004a18f05f51f97b08edd9def6826f", "dirty": False}}},
                        "design_sha256": artifacts.sha256(root / "docs/design/eamxx_mam4xx_ndrop.md"),
                        "python_version": artifacts.platform.python_version(), "platform": "synthetic_fixture",
                        "accuracy_command": ["inputs/benchmark_ndrop_matched_error", "--physical-controller-accuracy",
                            "--physical-grid-dir", "inputs/grids", "--max-steps", "8192", "--budget-seconds", "7200"],
                        "accuracy_protocol": {"integer_search": True, "maximum_steps": 8192,
                            "budget_seconds": 7200, "confirmation_multipliers": [2, 4]},
                        "timing_protocol": {"warmups": 10, "minimum_blocks": 30, "maximum_blocks": 100,
                            "budget_seconds": 7200, "complete_interleaved_blocks": True}}
            def save_identity(value):
                (directory / "identity.json").write_bytes(artifacts.canonical_json(value))
            save_identity(identity)
            artifacts.validate_identity(directory)
            for update in ({"comparison_policy": "fixed_1_percent"}, {"fixtures": [28, 29]},
                           {"accuracy_command": []}, {"extra": 1}, {"grid_files_sha256": {}},
                           {"source_files_sha256": {}}, {"binary_sha256": "0" * 64}):
                save_identity({**identity, **update})
                with self.assertRaises(artifacts.ArtifactError): artifacts.validate_identity(directory)
            save_identity(identity)
            rows = physical_controller_rows()
            with (directory / "accuracy.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, artifacts.PHYSICAL_ACCURACY_HEADER)
                writer.writeheader(); writer.writerows(rows)
            artifacts._completion(directory, "accuracy", physical=True)
            args = SimpleNamespace(output_dir=directory)
            artifacts.freeze(args)
            frozen = artifacts.selection(directory)
            self.assertEqual(len(artifacts.admission(directory, frozen)["unique_configs"]), 12)
            samples = [{"artifact_schema": "ndrop_physical_controller_timing_v1", "fixture": str(f),
                        "method": m, "steps": str(n), "block": str(b), "seconds_per_column": "1",
                        "state_checksum": str(1000 + f), "stop_reason": "median_precision"}
                       for b in range(30) for f, m, n in frozen["unique_configs"]]
            with (directory / "timing.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, artifacts.TIMING_HEADER)
                writer.writeheader(); writer.writerows(samples)
            artifacts._completion(directory, "timing", physical=True)
            artifacts.report(args); artifacts.validate_final(args)
            report = artifacts.load_json(directory / "analysis.json")
            self.assertEqual(report["evidence_scope"], "physical_controller_idealization_unweighted")
            self.assertTrue(all("turbulence" in f for f in report["fixtures"]))
            (directory / "inputs/grids/acme-72i.ascii").write_bytes(b"corrupt")
            with self.assertRaises(artifacts.ArtifactError): artifacts.selection(directory)

    def test_natural_fe_target_and_independent_schema(self):
        rows = physical_controller_rows()
        indexed = artifacts.validate_accuracy(rows, (28, 29, 30, 31), physical=True)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=True)
        self.assertEqual(frozen["artifact_schema"], "ndrop_physical_controller_selection_v1")
        self.assertEqual([e["steps"] for e in frozen["selected"]
                          if e["method"] == "forward_euler"], [1, 5, 3, 10])
        self.assertTrue(all(e["fe_target_error"] == .03 for e in frozen["selected"]))
        self.assertEqual([e["steps"] for e in frozen["selected"]
                          if e["method"] == "imex"], [2, 2, 2, 2])
        admitted = artifacts.expected_admission(rows, indexed, frozen, "c", physical=True)
        self.assertEqual(admitted["unique_configs"], frozen["unique_configs"])
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(rows)

    def test_corruption_never_refines_fe_or_replaces_its_target(self):
        rows = physical_controller_rows()
        for field, value in (("steps", "256"), ("fe_target_error", "0.01"),
                             ("artifact_schema", "ndrop_physical_column_accuracy_v1"),
                             ("cloud_base_wtke_m_s", "1"), ("controller_rate_t", "10"),
                             ("limiting_level", "72"), ("activation_limiter_count", "-1")):
            bad = copy.deepcopy(rows); bad[0][field] = value
            with self.subTest(field=field), self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_accuracy(bad, (28, 29, 30, 31), physical=True)
        missing = [r for r in rows if not (r["fixture"] == "28" and
                   r["method"] == "imex" and r["record_kind"] == "measurement" and r["steps"] == "1")]
        with self.assertRaisesRegex(artifacts.ArtifactError, "not consecutive"):
            artifacts.validate_accuracy(missing, (28, 29, 30, 31), physical=True)
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(rows[:-1], (28,), physical=True)
        stale_terminal = copy.deepcopy(rows)
        terminal = next(r for r in stale_terminal if r["record_kind"] == "search_terminal")
        terminal["fe_target_error"] = "0.01"
        with self.assertRaises(artifacts.ArtifactError):
            artifacts.validate_accuracy(stale_terminal, (28, 29, 30, 31), physical=True)

    def test_missing_confirmation_withholds_timing_not_selection(self):
        rows = physical_controller_rows()
        indexed = artifacts.validate_accuracy(rows, (28, 29, 30, 31), physical=True)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=True)
        rows = [r for r in rows if not (r["fixture"] == "28" and r["method"] == "imex"
                                       and r["record_kind"] == "confirmation")]
        indexed = artifacts.validate_accuracy(rows, (28, 29, 30, 31), physical=True)
        self.assertEqual(frozen, artifacts.expected_selection_from_rows(
            rows, indexed, "a", "b", physical=True))
        admitted = artifacts.expected_admission(rows, indexed, frozen, "c", physical=True)
        self.assertNotIn([28, "imex", 2], admitted["unique_configs"])

    def test_fabricated_timing_requires_paired_qualified_checksums(self):
        rows = physical_controller_rows()
        indexed = artifacts.validate_accuracy(rows, (28, 29, 30, 31), physical=True)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=True)
        samples = [{"artifact_schema": "ndrop_physical_controller_timing_v1",
                    "fixture": str(f), "method": m, "steps": str(n), "block": str(b),
                    "seconds_per_column": "1", "state_checksum": indexed[(f, m, n)]["state_checksum"],
                    "stop_reason": "median_precision"}
                   for b in range(30) for f, m, n in frozen["unique_configs"]]
        artifacts.validate_timing(samples, frozen, indexed, physical=True)
        for bad in (samples[:-1], [dict(samples[0], state_checksum="0"), *samples[1:]],
                    [dict(samples[0], artifact_schema=artifacts.TIMING_SCHEMA), *samples[1:]]):
            with self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_timing(bad, frozen, indexed, physical=True)

    @unittest.skipUnless(os.environ.get("NDROP_PHYSICAL_TEST_BINARY"), "explicit untimed local binary required")
    def test_actual_binary_uses_natural_controller_and_first_matches(self):
        root = Path(__file__).resolve().parents[4]
        command = [os.environ["NDROP_PHYSICAL_TEST_BINARY"], "--physical-controller-accuracy",
                   "--physical-grid-dir", str(root / "components/homme/test/vcoord"),
                   "--max-steps", "64", "--budget-seconds", "120"]
        process = subprocess.run(command, capture_output=True, text=True, timeout=150)
        self.assertEqual(process.returncode, 0, process.stderr)
        rows = list(csv.DictReader(process.stdout.splitlines()))
        indexed = artifacts.validate_accuracy(rows, (28, 29, 30, 31), physical=True)
        frozen = artifacts.expected_selection_from_rows(rows, indexed, "a", "b", physical=True)
        self.assertEqual([e["steps"] for e in frozen["selected"]
                          if e["method"] == "forward_euler"], [1, 5, 3, 10])
        self.assertTrue(all(e["achieved_error"] <= e["fe_target_error"] for e in frozen["selected"]))
        self.assertEqual(len(frozen["selected"]), 12)
        admitted = artifacts.expected_admission(rows, indexed, frozen, "c", physical=True)
        self.assertEqual(len(admitted["unique_configs"]), 12)
        for option in ([], ["--fixtures", "28"]):
            invalid = subprocess.run([command[0], command[1], *option], capture_output=True, text=True)
            self.assertNotEqual(invalid.returncode, 0)


if __name__ == "__main__":
    unittest.main()
