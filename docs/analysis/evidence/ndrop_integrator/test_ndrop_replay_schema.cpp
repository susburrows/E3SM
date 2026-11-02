#include "ndrop_replay_schema.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

void test_real_hex() {
  using ndrop_test::production_replay_v1::canonical_real_hex;
  check(canonical_real_hex(1.0) == "0x1.0000000000000p+0", "one");
  check(canonical_real_hex(0.5) == "0x1.0000000000000p-1", "half");
  check(canonical_real_hex(0.0) == "0x0.0000000000000p+0", "positive zero");
  check(canonical_real_hex(-0.0) == "-0x0.0000000000000p+0",
        "negative zero");
  check(canonical_real_hex(std::numeric_limits<double>::denorm_min()) ==
            "0x0.0000000000001p-1022", "least subnormal");
  bool rejected = false;
  try {
    static_cast<void>(canonical_real_hex(
        std::numeric_limits<double>::quiet_NaN()));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected, "nonfinite value rejected");
}

void test_measurement() {
  using namespace ndrop_test;
  using production_replay_v1::serialize_diagnostic_measurement;
  AdvanceResult complete;
  complete.ok = true;
  auto &detail = complete.diagnostics.replay_detail.emplace();
  detail.status = ReplayDiagnosticStatus::complete;
  detail.expected_substeps = 1;
  detail.completed_substeps = 1;
  for (int i = 0; i < 51; ++i) {
    detail.raw_minimum_valid[i] = true;
    detail.negative_counts_valid[i] = true;
    detail.fe_clipping_count_valid[i] = true;
    detail.fe_weighted_deficit_valid[i] = true;
    detail.physical_residual[i] = ReplayResidualMaximum{0.0, 0.0, 0.0, 0};
  }
  complete.diagnostics.droplet_raw_minimum = 0;
  for (int j = 0; j < num_constituents; ++j) {
    detail.total_ledger_valid[j] = true;
    detail.clipping_ledger_valid[j] = true;
    detail.operator_ledger_valid[j] = true;
    complete.diagnostics.field_clipping[j].interstitial_raw_minimum = 0;
    complete.diagnostics.field_clipping[j].cloudborne_raw_minimum = 0;
  }
  const std::string serialized =
      serialize_diagnostic_measurement(complete, "forward_euler");
  check(serialized.find("\"not_applicable\"") != std::string::npos,
        "FE matrix applicability");
  check(serialized.find("\"maximum_normalized\":\"0x0.0000000000000p+0\"")
            != std::string::npos, "completed zero is not null");
  detail.status = ReplayDiagnosticStatus::partial;
  detail.failure_stage = ReplayDiagnosticStage::raw_state;
  detail.failure_reason = ReplayDiagnosticReason::raw_nonfinite;
  detail.raw_minimum_valid[0] = false;
  complete.ok = false;
  const std::string partial =
      serialize_diagnostic_measurement(complete, "forward_euler");
  check(partial.find("\"failure_reason\":\"raw_nonfinite\"") !=
            std::string::npos, "OI-5 failure reason");
  check(partial.find("\"raw_minimum_by_field\":[null,") !=
            std::string::npos, "uncomputed entry is null");
  bool rejected = false;
  try {
    complete.diagnostics.replay_detail.reset();
    static_cast<void>(serialize_diagnostic_measurement(complete,
                                                        "forward_euler"));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected, "uncollected diagnostic rejected");
}

void test_distinct_field_and_pair_order() {
  using namespace ndrop_test;
  using production_replay_v1::canonical_real_hex;
  using production_replay_v1::serialize_diagnostic_measurement;
  AdvanceResult result;
  result.ok = true;
  auto &detail = result.diagnostics.replay_detail.emplace();
  detail.status = ReplayDiagnosticStatus::complete;
  detail.expected_substeps = 1;
  detail.completed_substeps = 1;
  std::vector<std::string> minima, roundoff, residuals, totals;
  result.diagnostics.droplet_raw_minimum = 1000;
  result.diagnostics.droplet_roundoff_negative_count = 999;
  for (int i = 0; i < 51; ++i) {
    detail.raw_minimum_valid[i] = true;
    detail.negative_counts_valid[i] = true;
    detail.fe_clipping_count_valid[i] = true;
    detail.fe_weighted_deficit_valid[i] = true;
    detail.physical_residual[i] =
        ReplayResidualMaximum{Real(i + 1), 1, Real(i + 1), 0};
    minima.push_back('"' + canonical_real_hex(
        i == 0 ? 1000.0 : i <= 25 ? Real(1100 + i - 1)
                                  : Real(1200 + i - 26)) + '"');
    roundoff.push_back(std::to_string(
        i == 0 ? 999 : i <= 25 ? i : i + 100));
    residuals.push_back('"' + canonical_real_hex(Real(i + 1)) + '"');
  }
  for (int j = 0; j < 25; ++j) {
    auto &field = result.diagnostics.field_clipping[j];
    field.interstitial_raw_minimum = 1100 + j;
    field.cloudborne_raw_minimum = 1200 + j;
    field.interstitial_roundoff_negative_count = j + 1;
    field.cloudborne_roundoff_negative_count = 126 + j;
    result.diagnostics.total_ledger_delta[j] = 100 + j;
    result.diagnostics.clipping_ledger_delta[j] = 10 + j;
    result.diagnostics.operator_ledger_delta[j] = 90;
    detail.total_ledger_valid[j] = true;
    detail.clipping_ledger_valid[j] = true;
    detail.operator_ledger_valid[j] = true;
    totals.push_back('"' + canonical_real_hex(Real(100 + j)) + '"');
  }
  const auto joined = [](const std::vector<std::string> &values) {
    std::string out = "[";
    for (const auto &value : values) {
      if (out.size() > 1) out += ',';
      out += value;
    }
    return out + "]";
  };
  const std::string blob =
      serialize_diagnostic_measurement(result, "forward_euler");
  check(blob.find("\"raw_minimum_by_field\":" + joined(minima)) !=
            std::string::npos, "51 distinct field minima map in order");
  check(blob.find("\"roundoff_negative_count_by_field\":" +
                  joined(roundoff)) != std::string::npos,
        "51 distinct field counts map in order");
  check(blob.find("\"normalized_by_field\":" + joined(residuals)) !=
            std::string::npos, "51 distinct residual tuples map in order");
  check(blob.find("\"total_ledger_delta\":" + joined(totals)) !=
            std::string::npos, "25 distinct pair ledgers map in order");
}

} // namespace

int main() {
  try {
    test_real_hex();
    test_measurement();
    test_distinct_field_and_pair_order();
    std::cout << "OI-5 schema unit tests passed\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
