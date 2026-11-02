"""Synthetic OI-5 schema fixtures only; no production admission is attempted."""

import copy
import hashlib
import json
import math
import os
import subprocess
import unittest

import ndrop_replay_artifacts as artifacts


ZERO = artifacts.real_hex(0.0)
ONE = artifacts.real_hex(1.0)
HASH = "a" * 64


def residual(required=True, complete=True):
    value = ZERO if required and complete else None
    index = 0 if required and complete else None
    return {
        "applicability": "required" if required else "not_applicable",
        "aggregation": "maximum_over_substeps_then_fields",
        "numerator_by_field": [value] * 51,
        "normalizer_by_field": [value] * 51,
        "normalized_by_field": [value] * 51,
        "worst_substep_by_field": [index] * 51,
        "maximum_normalized": value,
    }


def measurement(method="forward_euler", complete=True):
    fe = method == "forward_euler"
    value = ZERO if complete else None
    count = 0 if complete else None
    return {
        "method": method,
        "completion": {
            "status": "complete" if complete else "not_started",
            "expected_substeps": 1 if complete else 0,
            "completed_substeps": 1 if complete else 0,
            "failure_stage": None if complete else "diagnostics",
            "failure_reason": None if complete else "diagnostic_incomplete",
        },
        "total_ledger_delta": [value] * 25,
        "clipping_ledger_delta": [value] * 25,
        "operator_ledger_delta": [value] * 25,
        "raw_minimum_by_field": [value] * 51,
        "roundoff_negative_count_by_field": [count] * 51,
        "material_negative_count_by_field": [count] * 51,
        "fe_clipping_count_by_field": [count if fe else None] * 51,
        "fe_weighted_deficit_by_field": [value if fe else None] * 51,
        "matrix_residual": residual(method != "forward_euler", complete),
        "physical_residual": residual(True, complete),
    }


def context():
    fields = []
    for i in range(51):
        constituent = None if i == 0 else (i - 1) % 25
        phase = "droplet" if i == 0 else (
            "interstitial" if i <= 25 else "cloudborne")
        kind = "droplet_number" if i == 0 else (
            "aerosol_number" if constituent in (0, 8, 13, 21)
            else "aerosol_mass")
        fields.append({
            "name": "qcld" if i == 0 else f"synthetic-{constituent}",
            "phase": phase,
            "constituent_index": constituent,
            "quantity_kind": kind,
            "native_unit": "kg_per_kg" if kind == "aerosol_mass"
            else "number_per_kg",
        })
    result = {
        "run_id": "synthetic.oi5-unit",
        "replay_case_sha256": HASH,
        "source_identity_sha256": HASH,
        "record_role": "census",
        "case_status": "ELIGIBLE",
        "field_order": fields,
        "pair_order": fields[1:26],
    }
    for name in ("preregistration", "preregistration_approval",
                 "replay_run_resource", "replay_execution_approval"):
        result[name + "_path"] = "synthetic/" + name + ".json"
        result[name + "_sha256"] = HASH
    return result


class CanonicalRealTests(unittest.TestCase):
    def test_bit_preservation_and_rejection(self):
        values = [0.0, -0.0, 1.0, 0.5, math.ulp(0.0), -math.ulp(0.0)]
        for value in values:
            encoded = artifacts.real_hex(value)
            self.assertEqual(artifacts.real_hex(artifacts.parse_real_hex(encoded)),
                             encoded)
        self.assertNotEqual(artifacts.real_hex(0.0), artifacts.real_hex(-0.0))
        for bad in ("0x1p+0", "0x1.0p+0", "0x2p-1", "0x0.0000000000000p-1022",
                    "NaN", "0x1.0000000000000p+000"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                artifacts.parse_real_hex(bad)


class DiagnosticsTests(unittest.TestCase):
    def test_actual_cpp_reference_component_and_fail_closed_shapes(self):
        executable = os.environ.get("NDROP_OI5_CPP_UNIT")
        if not executable:
            self.skipTest("set NDROP_OI5_CPP_UNIT for C++/Python integration")
        output = subprocess.check_output([executable, "--emit-reference"],
                                         timeout=30)
        value = json.loads(output)
        artifacts.validate_reference_self_check(value)
        self.assertEqual(value["outcome"], "completed")
        self.assertTrue(value["passed"])
        failed = copy.deepcopy(value)
        failed.update(attempt_status="attempted", outcome="failed",
                      primary_completed=True, tight_completed=False,
                      failure_stage="tight_advance",
                      failure_reason="step_underflow",
                      maximum_normalized_disagreement=None, passed=None)
        artifacts.validate_reference_self_check(failed)
        failed["maximum_normalized_disagreement"] = ZERO
        with self.assertRaises(ValueError):
            artifacts.validate_reference_self_check(failed)

    def test_actual_cpp_fixture_measurements_round_trip(self):
        executable = os.environ.get("NDROP_OI5_CPP_UNIT")
        if not executable:
            self.skipTest("set NDROP_OI5_CPP_UNIT for C++/Python integration")
        for method in sorted(artifacts.METHODS):
            with self.subTest(method=method):
                output = subprocess.check_output(
                    [executable, "--emit-measurement", method], timeout=30)
                measured = json.loads(output)
                artifacts.validate_measurement(measured)
                parsed = artifacts.parse_artifact(
                    artifacts.serialize_artifact(measured, context()))
                self.assertEqual(parsed["method"], method)
                self.assertEqual(parsed["completion"]["status"], "complete")
                for key in ("raw_minimum_by_field",
                            "roundoff_negative_count_by_field",
                            "material_negative_count_by_field"):
                    self.assertEqual(len(parsed[key]), 51)
                    self.assertTrue(all(value is not None for value in
                                        parsed[key]))
                for i in range(25):
                    total = artifacts.parse_real_hex(
                        parsed["total_ledger_delta"][i])
                    clipping = artifacts.parse_real_hex(
                        parsed["clipping_ledger_delta"][i])
                    operator = artifacts.parse_real_hex(
                        parsed["operator_ledger_delta"][i])
                    self.assertTrue(math.isclose(total, clipping + operator,
                                                 abs_tol=1e-12))
                if method == "forward_euler":
                    self.assertTrue(all(value is not None for value in
                                        parsed["fe_clipping_count_by_field"]))
                    self.assertTrue(all(value is None for value in
                                        parsed["matrix_residual"]
                                        ["normalized_by_field"]))
                else:
                    self.assertTrue(all(value is None for value in
                                        parsed["fe_clipping_count_by_field"]))
                    self.assertTrue(all(value is not None for value in
                                        parsed["matrix_residual"]
                                        ["normalized_by_field"]))

    def test_all_methods_complete_round_trip(self):
        for method in sorted(artifacts.METHODS):
            with self.subTest(method=method):
                blob = artifacts.serialize_artifact(measurement(method), context())
                parsed = artifacts.parse_artifact(blob)
                self.assertEqual(parsed["method"], method)
                self.assertEqual(len(parsed["field_order"]), 51)
                self.assertEqual(len(parsed["pair_order"]), 25)
                self.assertEqual(parsed["diagnostics_sha256"],
                                 hashlib.sha256(artifacts._canonical({
                                     key: value for key, value in parsed.items()
                                     if key != "diagnostics_sha256"
                                 })).hexdigest())
                self.assertIsNone(parsed["matrix_residual"]["maximum_normalized"]
                                  if method == "forward_euler" else None)

    def test_partial_uncomputed_and_failed(self):
        m = measurement(complete=False)
        artifacts.validate_measurement(m)
        m["completion"].update(status="partial", expected_substeps=3,
                               completed_substeps=1, failure_stage="raw_state",
                               failure_reason="raw_nonfinite")
        m["raw_minimum_by_field"][0] = ONE
        m["roundoff_negative_count_by_field"][0] = 0
        m["material_negative_count_by_field"][0] = 0
        artifacts.validate_measurement(m)
        failed_context = context()
        failed_context["case_status"] = "NUMERICAL_FAILURE"
        blob = artifacts.serialize_artifact(m, failed_context)
        self.assertIsNone(artifacts.parse_artifact(blob)
                          ["physical_residual"]["maximum_normalized"])
        m["completion"]["failure_reason"] = "step_underflow"
        with self.assertRaises(ValueError):
            artifacts.validate_measurement(m)

    def test_implicit_inapplicability_and_zero_ledger(self):
        m = measurement("imex")
        m["fe_clipping_count_by_field"][0] = 0
        with self.assertRaises(ValueError):
            artifacts.validate_measurement(m)
        m = measurement("backward_euler")
        m["clipping_ledger_delta"][3] = ONE
        with self.assertRaises(ValueError):
            artifacts.validate_measurement(m)

    def test_case_status_requires_matching_diagnostic_completion(self):
        partial = measurement(complete=False)
        for case_status in ("ELIGIBLE", "REFERENCE_UNRESOLVED"):
            with self.subTest(case_status=case_status), self.assertRaises(ValueError):
                selected_context = context()
                selected_context["case_status"] = case_status
                artifacts.serialize_artifact(partial, selected_context)
        failed_context = context()
        failed_context["case_status"] = "NUMERICAL_FAILURE"
        artifacts.parse_artifact(artifacts.serialize_artifact(partial,
                                                               failed_context))
        with self.assertRaises(ValueError):
            artifacts.serialize_artifact(measurement(), failed_context)
        sentinel_context = context()
        sentinel_context["record_role"] = "STRESS_SENTINEL"
        sentinel_context["case_status"] = "REFERENCE_UNRESOLVED"
        with self.assertRaises(ValueError):
            artifacts.serialize_artifact(measurement(), sentinel_context)

    def test_torn_residual_and_wrong_maximum(self):
        m = measurement()
        m["physical_residual"]["normalizer_by_field"][2] = None
        with self.assertRaises(ValueError):
            artifacts.validate_measurement(m)
        m = measurement()
        m["physical_residual"]["maximum_normalized"] = ONE
        with self.assertRaises(ValueError):
            artifacts.validate_measurement(m)

    def test_corruption_and_closed_keys(self):
        blob = artifacts.serialize_artifact(measurement(), context())
        for bad in (blob[:-1], blob + b"\n", blob.replace(b'"artifact_schema"',
                                                        b'"wrong_schema"', 1),
                    blob.replace(b'"diagnostics_sha256":"',
                                 b'"diagnostics_sha256":"b', 1),
                    b'\xef\xbb\xbf' + blob):
            with self.subTest(bad=bad[:30]), self.assertRaises(ValueError):
                artifacts.parse_artifact(bad)
        value = json.loads(blob)
        value["extra"] = 1
        with self.assertRaises(ValueError):
            artifacts.validate_artifact(value)
        duplicate = blob.replace(b'"artifact_schema":',
                                 b'"artifact_schema":"duplicate",'
                                 b'"artifact_schema":', 1)
        with self.assertRaises(ValueError):
            artifacts.parse_artifact(duplicate)

    def test_synthetic_context_is_not_admission(self):
        parsed = artifacts.parse_artifact(
            artifacts.serialize_artifact(measurement(), context()))
        self.assertEqual(parsed["run_id"], "synthetic.oi5-unit")
        self.assertFalse(hasattr(artifacts, "admit_census_replay_input"))
        self.assertFalse(hasattr(artifacts, "validate_replay_execution_manifest"))


if __name__ == "__main__":
    unittest.main()
