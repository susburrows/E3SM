#include "ndrop_replay_schema.hpp"
#include "ndrop_production_replay.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ndrop_test::production_replay_v1 {
namespace {

std::string quote(const std::string &value) {
  std::string result = "\"";
  for (const unsigned char c : value) {
    switch (c) {
    case '\"': result += "\\\""; break;
    case '\\': result += "\\\\"; break;
    case '\b': result += "\\b"; break;
    case '\t': result += "\\t"; break;
    case '\n': result += "\\n"; break;
    case '\f': result += "\\f"; break;
    case '\r': result += "\\r"; break;
    default:
      if (c < 0x20) {
        const char hex[] = "0123456789abcdef";
        result += "\\u00";
        result += hex[c >> 4];
        result += hex[c & 15];
      } else {
        result += static_cast<char>(c);
      }
    }
  }
  return result + "\"";
}

std::string object(const std::map<std::string, std::string> &fields) {
  std::string result = "{";
  for (const auto &[key, value] : fields) {
    if (result.size() > 1) result += ',';
    result += quote(key) + ":" + value;
  }
  return result + "}";
}

std::string array(const std::vector<std::string> &values) {
  std::string result = "[";
  for (const auto &value : values) {
    if (result.size() > 1) result += ',';
    result += value;
  }
  return result + "]";
}

template <std::size_t N, class Getter>
std::string masked_array(const std::array<bool, N> &valid, Getter get) {
  std::vector<std::string> values;
  values.reserve(N);
  for (std::size_t i = 0; i < N; ++i)
    values.push_back(valid[i] ? get(i) : "null");
  return array(values);
}

std::string real_or_null(Real value) { return quote(canonical_real_hex(value)); }

std::string status_name(ReplayDiagnosticStatus status) {
  switch (status) {
  case ReplayDiagnosticStatus::not_started: return "not_started";
  case ReplayDiagnosticStatus::partial: return "partial";
  case ReplayDiagnosticStatus::complete: return "complete";
  }
  throw std::invalid_argument("unknown replay diagnostic status");
}

std::string stage_name(ReplayDiagnosticStage stage) {
  static const char *const names[] = {
      "input", "assembly", "factorization", "solve", "raw_state",
      "clipping", "projection", "matrix_residual", "physical_residual",
      "ledger", "boundary", "branch", "diagnostics"};
  const auto index = static_cast<std::size_t>(stage);
  if (index >= std::size(names)) throw std::invalid_argument("unknown stage");
  return names[index];
}

std::string reason_name(ReplayDiagnosticReason reason) {
  static const char *const names[] = {
      "input_nonfinite", "assembly_nonfinite", "factorization_unsafe_pivot",
      "factorization_nonfinite", "solve_nonfinite", "raw_material_negative",
      "raw_nonfinite", "clipping_nonfinite", "projection_nonfinite",
      "matrix_residual_nonfinite", "matrix_residual_exceeded",
      "physical_residual_nonfinite", "physical_residual_exceeded",
      "ledger_nonfinite", "ledger_mismatch", "boundary_mismatch",
      "branch_mismatch", "diagnostic_incomplete"};
  const auto index = static_cast<std::size_t>(reason);
  if (index >= std::size(names)) throw std::invalid_argument("unknown reason");
  return names[index];
}

std::string residual_object(
    const std::array<std::optional<ReplayResidualMaximum>, 51> &values,
    bool required) {
  std::vector<std::string> numerator, normalizer, normalized, substeps;
  Real maximum = 0;
  bool any = false;
  for (const auto &value : values) {
    if (value) {
      any = true;
      if (value->numerator < 0 || value->normalizer < 0 ||
          value->normalized < 0)
        throw std::invalid_argument("negative residual evidence");
      maximum = std::max(maximum, value->normalized);
      numerator.push_back(real_or_null(value->numerator));
      normalizer.push_back(real_or_null(value->normalizer));
      normalized.push_back(real_or_null(value->normalized));
      substeps.push_back(std::to_string(value->substep));
    } else {
      numerator.push_back("null");
      normalizer.push_back("null");
      normalized.push_back("null");
      substeps.push_back("null");
    }
  }
  if (!required && any) throw std::invalid_argument("inapplicable residual");
  return object({{"aggregation", quote("maximum_over_substeps_then_fields")},
                 {"applicability", quote(required ? "required" : "not_applicable")},
                 {"maximum_normalized", any ? real_or_null(maximum) : "null"},
                 {"normalized_by_field", array(normalized)},
                 {"normalizer_by_field", array(normalizer)},
                 {"numerator_by_field", array(numerator)},
                 {"worst_substep_by_field", array(substeps)}});
}

} // namespace

std::string canonical_real_hex(Real value) {
  static_assert(sizeof(Real) == sizeof(std::uint64_t), "binary64 required");
  if (!std::isfinite(value)) throw std::invalid_argument("nonfinite real_hex");
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  const bool negative = (bits >> 63) != 0;
  const auto exponent = static_cast<int>((bits >> 52) & 0x7ff);
  const auto fraction = bits & ((std::uint64_t(1) << 52) - 1);
  std::ostringstream out;
  out.imbue(std::locale::classic());
  if (negative) out << '-';
  out << "0x" << (exponent == 0 ? '0' : '1') << '.'
      << std::hex << std::nouppercase << std::setfill('0') << std::setw(13)
      << fraction << std::dec << "p";
  const int power = exponent == 0 ? (fraction == 0 ? 0 : -1022)
                                  : exponent - 1023;
  if (power >= 0) out << '+';
  out << power;
  return out.str();
}

std::string serialize_diagnostic_measurement(const AdvanceResult &result,
                                             const std::string &method) {
  if (method != "forward_euler" && method != "imex" &&
      method != "backward_euler")
    throw std::invalid_argument("unknown replay method");
  if (!result.diagnostics.replay_detail)
    throw std::invalid_argument("replay detail was not collected");
  const auto &d = result.diagnostics;
  const auto &x = *d.replay_detail;
  if (x.status == ReplayDiagnosticStatus::complete && !result.ok)
    throw std::invalid_argument("inconsistent diagnostic completion");
  if (x.completed_substeps > x.expected_substeps)
    throw std::invalid_argument("invalid completed substep count");
  const auto minimum = [&](std::size_t i) {
    if (i == 0) return real_or_null(d.droplet_raw_minimum);
    const int j = i <= 25 ? static_cast<int>(i - 1)
                          : static_cast<int>(i - 26);
    return real_or_null(i <= 25
                            ? d.field_clipping[j].interstitial_raw_minimum
                            : d.field_clipping[j].cloudborne_raw_minimum);
  };
  const auto count = [&](std::size_t i, bool material) {
    if (i == 0)
      return std::to_string(material ? d.droplet_material_negative_count
                                     : d.droplet_roundoff_negative_count);
    const int j = i <= 25 ? static_cast<int>(i - 1)
                          : static_cast<int>(i - 26);
    const auto &f = d.field_clipping[j];
    return std::to_string(
        i <= 25 ? (material ? f.interstitial_material_negative_count
                            : f.interstitial_roundoff_negative_count)
                : (material ? f.cloudborne_material_negative_count
                            : f.cloudborne_roundoff_negative_count));
  };
  const auto fe_count = [&](std::size_t i) {
    if (i == 0) return std::to_string(d.droplet_clipping_count);
    const int j = i <= 25 ? static_cast<int>(i - 1)
                          : static_cast<int>(i - 26);
    return std::to_string(i <= 25 ? d.field_clipping[j].interstitial_count
                                   : d.field_clipping[j].cloudborne_count);
  };
  const auto fe_deficit = [&](std::size_t i) {
    if (i == 0) return real_or_null(d.droplet_weighted_deficit);
    const int j = i <= 25 ? static_cast<int>(i - 1)
                          : static_cast<int>(i - 26);
    return real_or_null(i <= 25
                            ? d.field_clipping[j].interstitial_weighted_deficit
                            : d.field_clipping[j].cloudborne_weighted_deficit);
  };
  const auto ledger = [&](const std::array<bool, 25> &mask,
                          const std::array<Real, 25> &values) {
    return masked_array(mask, [&](std::size_t i) {
      return real_or_null(values[i]);
    });
  };
  const std::string completion = object({
      {"status", quote(status_name(x.status))},
      {"expected_substeps", std::to_string(x.expected_substeps)},
      {"completed_substeps", std::to_string(x.completed_substeps)},
      {"failure_stage", x.failure_stage ? quote(stage_name(*x.failure_stage)) : "null"},
      {"failure_reason", x.failure_reason ? quote(reason_name(*x.failure_reason)) : "null"}});
  return object({
      {"method", quote(method)}, {"completion", completion},
      {"raw_minimum_by_field", masked_array(x.raw_minimum_valid, minimum)},
      {"roundoff_negative_count_by_field", masked_array(x.negative_counts_valid,
                                                   [&](std::size_t i) { return count(i, false); })},
      {"material_negative_count_by_field", masked_array(x.negative_counts_valid,
                                                   [&](std::size_t i) { return count(i, true); })},
      {"fe_clipping_count_by_field", masked_array(x.fe_clipping_count_valid, fe_count)},
      {"fe_weighted_deficit_by_field", masked_array(x.fe_weighted_deficit_valid, fe_deficit)},
      {"total_ledger_delta", ledger(x.total_ledger_valid, d.total_ledger_delta)},
      {"clipping_ledger_delta", ledger(x.clipping_ledger_valid, d.clipping_ledger_delta)},
      {"operator_ledger_delta", ledger(x.operator_ledger_valid, d.operator_ledger_delta)},
      {"matrix_residual", residual_object(x.matrix_residual, method != "forward_euler")},
      {"physical_residual", residual_object(x.physical_residual, true)}});
}

std::string serialize_reference_self_check(const ReferenceSelfCheck &check) {
  const bool not_applicable = check.applicability == "not_applicable";
  const bool required = check.applicability == "required";
  const bool not_attempted = check.attempt_status == "not_attempted" &&
                             check.outcome == "not_attempted";
  const bool attempted_failure = check.attempt_status == "attempted" &&
                                 check.outcome == "failed";
  const bool completed = check.attempt_status == "attempted" &&
                         check.outcome == "completed";
  const bool no_metric = !check.maximum_normalized_disagreement &&
                         !check.passed;
  if ((!not_applicable && !required) ||
      check.primary_rtol != Real(1e-12) ||
      check.tight_rtol != Real(1e-13) ||
      check.threshold != Real(2e-11))
    throw std::invalid_argument("invalid reference self-check constants");
  if (not_attempted) {
    if (!no_metric || check.failure_stage || check.primary_completed ||
        check.tight_completed || !check.failure_reason ||
        (not_applicable &&
         *check.failure_reason != "sentinel_excluded" &&
         *check.failure_reason != "record_not_replayable") ||
        (required && *check.failure_reason != "candidate_failed_first"))
      throw std::invalid_argument("invalid not-attempted reference");
  } else if (attempted_failure) {
    if (!required || !no_metric || !check.failure_stage ||
        !check.failure_reason ||
        (check.tight_completed && !check.primary_completed))
      throw std::invalid_argument("invalid failed reference");
    const auto &stage = *check.failure_stage;
    const auto &reason = *check.failure_reason;
    if (stage != "primary_setup" && stage != "primary_advance" &&
        stage != "tight_setup" && stage != "tight_advance" &&
        stage != "pair_comparison")
      throw std::invalid_argument("unknown reference failure stage");
    if (stage == "pair_comparison") {
      if (!check.primary_completed || !check.tight_completed ||
          reason != "comparison_nonfinite")
        throw std::invalid_argument("invalid pair comparison failure");
    } else {
      if ((stage.rfind("primary_", 0) == 0 && check.primary_completed) ||
          (stage.rfind("tight_", 0) == 0 &&
           (!check.primary_completed || check.tight_completed)) ||
          (stage.find("_setup") != std::string::npos
               ? reason != "nonfinite_stage" &&
                     reason != "incomplete_reference_pair"
               : reason != "nonfinite_stage" &&
                     reason != "step_underflow" &&
                     reason != "step_limit_exhausted" &&
                     reason != "incomplete_reference_pair"))
        throw std::invalid_argument("invalid reference failure state");
    }
  } else if (completed) {
    if (!required || !check.primary_completed || !check.tight_completed ||
        check.failure_stage || check.failure_reason ||
        !check.maximum_normalized_disagreement || !check.passed ||
        !std::isfinite(*check.maximum_normalized_disagreement) ||
        *check.maximum_normalized_disagreement < 0 ||
        *check.passed !=
            (*check.maximum_normalized_disagreement <= check.threshold))
      throw std::invalid_argument("invalid completed reference");
  } else {
    throw std::invalid_argument("invalid reference outcome");
  }
  return object({
      {"applicability", quote(check.applicability)},
      {"attempt_status", quote(check.attempt_status)},
      {"outcome", quote(check.outcome)},
      {"failure_stage", check.failure_stage ? quote(*check.failure_stage) : "null"},
      {"failure_reason", check.failure_reason ? quote(*check.failure_reason) : "null"},
      {"primary_rtol", real_or_null(check.primary_rtol)},
      {"tight_rtol", real_or_null(check.tight_rtol)},
      {"primary_completed", check.primary_completed ? "true" : "false"},
      {"tight_completed", check.tight_completed ? "true" : "false"},
      {"maximum_normalized_disagreement",
       check.maximum_normalized_disagreement
           ? real_or_null(*check.maximum_normalized_disagreement) : "null"},
      {"threshold", real_or_null(check.threshold)},
      {"passed", check.passed ? (*check.passed ? "true" : "false") : "null"}});
}

} // namespace ndrop_test::production_replay_v1
