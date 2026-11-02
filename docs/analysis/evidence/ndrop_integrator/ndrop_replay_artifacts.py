"""Closed OI-5 diagnostics serialization; schema validation is not admission.

No function here authorizes a capture, replay, timing row, or result consumer.
Production admission requires the separate approved manifest/completion gates.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import struct
from typing import Any


SCHEMA = "ndrop_replay_diagnostics_v1"
RUN_ID = re.compile(r"[A-Za-z0-9._-]{1,128}\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
HEX = re.compile(r"-?(?:0x1\.[0-9a-f]{13}p[+-](?:0|[1-9][0-9]*)|"
                 r"0x0\.[0-9a-f]{13}p-1022|0x0\.0000000000000p\+0)\Z")
METHODS = {"forward_euler", "imex", "backward_euler"}
STAGES = {
    "input": {"input_nonfinite"},
    "assembly": {"assembly_nonfinite"},
    "factorization": {"factorization_unsafe_pivot", "factorization_nonfinite"},
    "solve": {"solve_nonfinite"},
    "raw_state": {"raw_material_negative", "raw_nonfinite"},
    "clipping": {"clipping_nonfinite"},
    "projection": {"projection_nonfinite"},
    "matrix_residual": {"matrix_residual_nonfinite", "matrix_residual_exceeded"},
    "physical_residual": {"physical_residual_nonfinite", "physical_residual_exceeded"},
    "ledger": {"ledger_nonfinite", "ledger_mismatch"},
    "boundary": {"boundary_mismatch"},
    "branch": {"branch_mismatch"},
    "diagnostics": {"diagnostic_incomplete"},
}
MEASUREMENT_KEYS = {
    "method", "completion", "total_ledger_delta", "clipping_ledger_delta",
    "operator_ledger_delta", "raw_minimum_by_field",
    "roundoff_negative_count_by_field", "material_negative_count_by_field",
    "fe_clipping_count_by_field", "fe_weighted_deficit_by_field",
    "matrix_residual", "physical_residual",
}
CONTEXT_KEYS = {
    "run_id", "preregistration_path", "preregistration_sha256",
    "preregistration_approval_path", "preregistration_approval_sha256",
    "replay_run_resource_path", "replay_run_resource_sha256",
    "replay_execution_approval_path", "replay_execution_approval_sha256",
    "replay_case_sha256", "source_identity_sha256", "record_role",
    "case_status", "field_order", "pair_order",
}
ARTIFACT_KEYS = MEASUREMENT_KEYS | CONTEXT_KEYS | {
    "artifact_schema", "diagnostics_sha256"
}
RESIDUAL_KEYS = {
    "applicability", "aggregation", "numerator_by_field",
    "normalizer_by_field", "normalized_by_field", "maximum_normalized",
    "worst_substep_by_field",
}
FIELD_KEYS = {"name", "phase", "constituent_index", "quantity_kind", "native_unit"}
REFERENCE_KEYS = {
    "applicability", "attempt_status", "outcome", "failure_stage",
    "failure_reason", "primary_rtol", "tight_rtol", "primary_completed",
    "tight_completed", "maximum_normalized_disagreement", "threshold",
    "passed",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _closed(value: Any, keys: set[str], name: str) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == keys, f"{name}: closed keys")
    return value


def _uint(value: Any, name: str, bits: int = 64) -> None:
    _require(type(value) is int and 0 <= value < 1 << bits, f"{name}: uint{bits}")


def real_hex(value: float) -> str:
    _require(type(value) is float and math.isfinite(value), "finite binary64 required")
    bits = struct.unpack(">Q", struct.pack(">d", value))[0]
    sign = "-" if bits >> 63 else ""
    exponent = (bits >> 52) & 0x7FF
    fraction = bits & ((1 << 52) - 1)
    if exponent == 0:
        power = -1022 if fraction else 0
        return f"{sign}0x0.{fraction:013x}p{power:+d}"
    return f"{sign}0x1.{fraction:013x}p{exponent - 1023:+d}"


def parse_real_hex(value: Any) -> float:
    _require(type(value) is str and HEX.fullmatch(value) is not None,
             "noncanonical real_hex")
    try:
        number = float.fromhex(value)
    except ValueError as exc:
        raise ValueError("invalid real_hex") from exc
    _require(math.isfinite(number) and real_hex(number) == value,
             "real_hex does not round-trip")
    return number


def validate_reference_self_check(value: Any) -> dict[str, Any]:
    """Validate the approved closed result component, without admitting replay."""
    value = _closed(value, REFERENCE_KEYS, "reference_self_check")
    for key, expected in (("primary_rtol", 1e-12), ("tight_rtol", 1e-13),
                          ("threshold", 2e-11)):
        _require(value[key] == real_hex(expected), f"{key}: frozen value")
    applicability = value["applicability"]
    attempt = value["attempt_status"]
    outcome = value["outcome"]
    stage = value["failure_stage"]
    reason = value["failure_reason"]
    primary = value["primary_completed"]
    tight = value["tight_completed"]
    metric = value["maximum_normalized_disagreement"]
    passed = value["passed"]
    _require(applicability in {"required", "not_applicable"} and
             type(primary) is bool and type(tight) is bool,
             "reference applicability/completion")
    if attempt == "not_attempted" and outcome == "not_attempted":
        _require(not primary and not tight and stage is None and
                 metric is None and passed is None and
                 reason in ({"candidate_failed_first"} if applicability ==
                            "required" else {"sentinel_excluded",
                                               "record_not_replayable"}),
                 "not-attempted reference")
    elif attempt == "attempted" and outcome == "failed":
        _require(applicability == "required" and metric is None and
                 passed is None and type(stage) is str and
                 type(reason) is str and (not tight or primary),
                 "failed reference")
        if stage == "pair_comparison":
            _require(primary and tight and reason == "comparison_nonfinite",
                     "pair-comparison failure")
        else:
            _require(stage in {"primary_setup", "primary_advance",
                               "tight_setup", "tight_advance"},
                     "reference failure stage")
            _require((not primary if stage.startswith("primary_") else
                      primary and not tight), "reference completion order")
            _require(reason in ({"incomplete_reference_pair",
                                 "nonfinite_stage"} if stage.endswith("setup")
                                else {"incomplete_reference_pair",
                                      "nonfinite_stage", "step_underflow",
                                      "step_limit_exhausted"}),
                     "reference stage/reason")
    elif attempt == "attempted" and outcome == "completed":
        _require(applicability == "required" and primary and tight and
                 stage is None and reason is None and type(passed) is bool,
                 "completed reference")
        observed = parse_real_hex(metric)
        _require(observed >= 0 and passed == (observed <= 2e-11),
                 "reference self-consistency")
    else:
        raise ValueError("invalid reference attempt/outcome")
    return value


def _canonical(value: Any) -> bytes:
    def inspect(item: Any) -> None:
        if type(item) is dict:
            _require(all(type(key) is str for key in item), "nonstring JSON key")
            for key, child in item.items():
                _require(not any(0xD800 <= ord(c) <= 0xDFFF for c in key),
                         "surrogate JSON key")
                inspect(child)
        elif type(item) is list:
            for child in item:
                inspect(child)
        elif type(item) is str:
            _require(not any(0xD800 <= ord(c) <= 0xDFFF for c in item),
                     "surrogate string")
        else:
            _require(item is None or type(item) in (int, bool),
                     "JSON floating point or unsupported type")
    inspect(value)
    return (json.dumps(value, ensure_ascii=False, sort_keys=True,
                       separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def _vector(value: Any, count: int, kind: str, name: str) -> list[Any]:
    _require(type(value) is list and len(value) == count, f"{name}: length {count}")
    for item in value:
        if item is None:
            continue
        if kind == "real":
            parse_real_hex(item)
        else:
            _uint(item, name)
    return value


def _residual(value: Any, required: bool, complete: bool,
              expected: int, name: str) -> None:
    value = _closed(value, RESIDUAL_KEYS, name)
    _require(value["applicability"] == ("required" if required else "not_applicable"),
             f"{name}: applicability")
    _require(value["aggregation"] == "maximum_over_substeps_then_fields",
             f"{name}: aggregation")
    rows = [
        _vector(value[key], 51, "count" if key == "worst_substep_by_field"
                else "real", f"{name}.{key}")
        for key in ("numerator_by_field", "normalizer_by_field",
                    "normalized_by_field", "worst_substep_by_field")
    ]
    observed: list[float] = []
    for i, tuple_values in enumerate(zip(*rows)):
        _require(all(x is None for x in tuple_values) or
                 all(x is not None for x in tuple_values),
                 f"{name}: torn residual tuple {i}")
        if tuple_values[0] is not None:
            _require(required, f"{name}: inapplicable tuple")
            numerator, normalizer, normalized = map(parse_real_hex,
                                                      tuple_values[:3])
            _require(numerator >= 0 and normalizer >= 0 and normalized >= 0,
                     f"{name}: negative residual")
            _require(real_hex(numerator / normalizer if normalizer > 0
                              else numerator) == tuple_values[2],
                     f"{name}: normalization mismatch")
            _require(tuple_values[3] < expected, f"{name}: substep range")
            observed.append(normalized)
        elif complete and required:
            raise ValueError(f"{name}: missing complete tuple {i}")
    if not observed:
        _require(value["maximum_normalized"] is None,
                 f"{name}: fabricated maximum")
    else:
        _require(value["maximum_normalized"] == real_hex(max(observed)),
                 f"{name}: maximum mismatch")


def validate_measurement(value: Any) -> dict[str, Any]:
    value = _closed(value, MEASUREMENT_KEYS, "diagnostic measurement")
    method = value["method"]
    _require(method in METHODS, "method")
    c = _closed(value["completion"], {
        "status", "expected_substeps", "completed_substeps", "failure_stage",
        "failure_reason"}, "completion")
    _require(c["status"] in {"not_started", "partial", "complete"},
             "completion status")
    _uint(c["expected_substeps"], "expected_substeps", 32)
    _uint(c["completed_substeps"], "completed_substeps", 32)
    _require(c["completed_substeps"] <= c["expected_substeps"],
             "completed exceeds expected")
    complete = c["status"] == "complete"
    if complete:
        _require(c["failure_stage"] is None and c["failure_reason"] is None and
                 c["completed_substeps"] == c["expected_substeps"] and
                 c["expected_substeps"] > 0, "inconsistent completion")
    else:
        _require(c["failure_stage"] in STAGES and
                 c["failure_reason"] in STAGES[c["failure_stage"]],
                 "failure stage/reason")
    if method != "forward_euler":
        _require(c["expected_substeps"] in (0, 1), "implicit substeps")
    for key in ("total_ledger_delta", "clipping_ledger_delta",
                "operator_ledger_delta"):
        _vector(value[key], 25, "real", key)
        if complete:
            _require(all(x is not None for x in value[key]), f"{key}: incomplete")
    for key in ("raw_minimum_by_field", "roundoff_negative_count_by_field",
                "material_negative_count_by_field", "fe_clipping_count_by_field",
                "fe_weighted_deficit_by_field"):
        kind = "count" if "count" in key else "real"
        _vector(value[key], 51, kind, key)
        if complete and (method == "forward_euler" or not key.startswith("fe_")):
            _require(all(x is not None for x in value[key]), f"{key}: incomplete")
    for i in range(51):
        _require((value["roundoff_negative_count_by_field"][i] is None) ==
                 (value["material_negative_count_by_field"][i] is None),
                 "negative-count pair nullability")
    if method != "forward_euler":
        for key in ("fe_clipping_count_by_field", "fe_weighted_deficit_by_field"):
            _require(all(x is None for x in value[key]), f"{key}: inapplicable")
        for entry in value["clipping_ledger_delta"]:
            if entry is not None:
                _require(entry == real_hex(0.0), "implicit clipping ledger")
    _residual(value["matrix_residual"], method != "forward_euler", complete,
              c["expected_substeps"], "matrix_residual")
    _residual(value["physical_residual"], True, complete,
              c["expected_substeps"], "physical_residual")
    return value


def _field_order(value: Any, pairs: Any) -> None:
    _require(type(value) is list and len(value) == 51 and
             type(pairs) is list and len(pairs) == 25, "field/pair ordering")
    for i, item in enumerate(value):
        _closed(item, FIELD_KEYS, f"field {i}")
        _require(type(item["name"]) is str and bool(item["name"]), "field name")
        _require(item["phase"] == ("droplet" if i == 0 else
                                   "interstitial" if i <= 25 else "cloudborne"),
                 "field phase")
        _require(item["quantity_kind"] in {"droplet_number", "aerosol_number",
                                             "aerosol_mass"}, "quantity kind")
        _require(item["native_unit"] == ("kg_per_kg" if item["quantity_kind"] ==
                                               "aerosol_mass" else "number_per_kg"),
                 "field unit")
        if i == 0:
            _require(item["constituent_index"] is None and
                     item["quantity_kind"] == "droplet_number", "droplet order")
        else:
            _uint(item["constituent_index"], "constituent index")
            _require(item["constituent_index"] == (i - 1) % 25,
                     "constituent order")
    for i, pair in enumerate(pairs):
        _closed(pair, FIELD_KEYS, f"pair {i}")
        _require(pair == value[1 + i] and
                 {**pair, "phase": "cloudborne"} == value[26 + i],
                 "pair/field mapping")


def validate_artifact(value: Any) -> dict[str, Any]:
    value = _closed(value, ARTIFACT_KEYS, "diagnostics artifact")
    _require(value["artifact_schema"] == SCHEMA, "artifact schema")
    _require(type(value["run_id"]) is str and RUN_ID.fullmatch(value["run_id"]),
             "run_id")
    for key in ("replay_case_sha256", "source_identity_sha256"):
        _require(type(value[key]) is str and SHA256.fullmatch(value[key]), key)
    for prefix in ("preregistration", "preregistration_approval",
                   "replay_run_resource", "replay_execution_approval"):
        path, digest = value[prefix + "_path"], value[prefix + "_sha256"]
        _require(type(path) is str and path and not path.startswith("/") and
                 all(component not in ("", ".", "..") for component in path.split("/")),
                 prefix + " path")
        _require(type(digest) is str and SHA256.fullmatch(digest), prefix + " hash")
    _require(value["record_role"] in {"census", "STRESS_SENTINEL"},
             "record role")
    _require(value["case_status"] in {"ELIGIBLE", "NUMERICAL_FAILURE",
                                       "REFERENCE_UNRESOLVED"}, "case status")
    if value["record_role"] == "STRESS_SENTINEL":
        _require(value["case_status"] != "REFERENCE_UNRESOLVED",
                 "sentinel reference is inapplicable")
    _field_order(value["field_order"], value["pair_order"])
    measured = validate_measurement({key: value[key]
                                     for key in MEASUREMENT_KEYS})
    completion = measured["completion"]["status"]
    if value["case_status"] in {"ELIGIBLE", "REFERENCE_UNRESOLVED"}:
        _require(completion == "complete", "passing candidate needs complete diagnostics")
    else:
        _require(completion != "complete", "candidate failure cannot be complete")
    _require(type(value["diagnostics_sha256"]) is str and
             SHA256.fullmatch(value["diagnostics_sha256"]), "diagnostics hash")
    content = {key: val for key, val in value.items()
               if key != "diagnostics_sha256"}
    _require(hashlib.sha256(_canonical(content)).hexdigest() ==
             value["diagnostics_sha256"], "diagnostics self-hash")
    return value


def serialize_artifact(measurement: dict[str, Any],
                       context: dict[str, Any]) -> bytes:
    """Build a schema-valid object; this is not a completed replay transaction."""
    validate_measurement(measurement)
    _closed(context, CONTEXT_KEYS, "diagnostics context")
    value = {**measurement, **context, "artifact_schema": SCHEMA}
    value["diagnostics_sha256"] = hashlib.sha256(_canonical(value)).hexdigest()
    validate_artifact(value)
    return _canonical(value)


def parse_artifact(blob: bytes) -> dict[str, Any]:
    _require(type(blob) is bytes and blob.endswith(b"\n") and
             not blob.startswith(b"\xef\xbb\xbf"), "canonical artifact framing")
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, val in items:
            _require(key not in result, "duplicate JSON key")
            result[key] = val
        return result
    try:
        value = json.loads(blob.decode("utf-8"), object_pairs_hook=pairs,
                           parse_float=lambda _: (_ for _ in ()).throw(
                               ValueError("JSON floating point")),
                           parse_constant=lambda _: (_ for _ in ()).throw(
                               ValueError("nonfinite JSON constant")))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("invalid UTF-8 JSON") from exc
    _require(_canonical(value) == blob, "noncanonical JSON bytes")
    return validate_artifact(value)
