#ifdef NDROP_LEGACY_BASELINE
#include NDROP_LEGACY_BASELINE_HEADER
#else
#include "ndrop_integrator.hpp"
#include "ndrop_production_replay.hpp"
#include "ndrop_replay_schema.hpp"
#endif

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
using namespace ndrop_test;

#ifndef NDROP_LEGACY_BASELINE
void check(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}
#endif

std::uint64_t bits(Real value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

struct Fixture {
  State state{3};
  OperatorData op{3};
  Real dt = 0.2;
  Fixture() {
    for (int k = 0; k < 3; ++k) {
      op.pressure_weight[k] = 1 + 0.2 * k;
      op.cloud_fraction[k] = k == 2 ? 0 : 0.7;
      state.droplet_number[k] = 2 + 0.01 * k;
      for (int j = 0; j < num_constituents; ++j) {
        state.interstitial[j][k] = 1 + 0.01 * j + 0.001 * k;
        state.cloudborne[j][k] = 0.2 + 0.005 * j + 0.001 * k;
      }
      for (int mode = 0; mode < num_modes; ++mode) {
        op.number_activation_rate[mode][k] = k == 2 ? 0 : 0.003;
        op.mass_activation_rate[mode][k] = k == 2 ? 0 : 0.002;
      }
    }
    op.interface_exchange[1] = 0.07;
    op.interface_exchange[2] = 0.09;
    op.first_level_controller_rate = 0.02;
  }
};

AdvanceResult advance(const std::string &method, const Fixture &f,
                      const AdvanceControl &control) {
  if (method == "forward_euler") {
    ForwardEulerOptions options;
    options.forced_substeps = 3;
    return advance_forward_euler(f.state, f.op, f.dt, options, control);
  }
  if (method == "imex") return advance_imex(f.state, f.op, f.dt, control);
  return advance_backward_euler(f.state, f.op, f.dt, control);
}

void emit_real(Real value) { std::cout << bits(value) << ','; }

void emit_snapshot(const AdvanceResult &r, const std::string &method) {
  std::cout << method << ':' << r.ok << ':' << r.message << ':';
  const auto &d = r.diagnostics;
  std::cout << d.nsubmix << ',' << d.final_buffer << ','
            << d.activation_limiter_count << ',' << d.clipping_count << ',';
  for (Real value : r.state.droplet_number) emit_real(value);
  for (int j = 0; j < num_constituents; ++j)
    for (Real value : r.state.interstitial[j]) emit_real(value);
  for (int j = 0; j < num_constituents; ++j)
    for (Real value : r.state.cloudborne[j]) emit_real(value);
  for (Real value : {d.raw_minimum, d.droplet_raw_minimum,
                     d.droplet_weighted_deficit, d.max_inverse_mixing_time,
                     d.first_level_controller_rate, d.normalized_residual,
                     d.normalized_matrix_residual,
                     d.normalized_physical_residual,
                     d.max_normalized_ledger_error, d.solution_minimum})
    emit_real(value);
  std::cout << d.droplet_clipping_count << ',' << d.negative_count << ','
            << d.roundoff_negative_count << ',' << d.material_negative_count
            << ',' << d.block_matrix_family_count << ','
            << d.block_factorization_count << ',' << d.block_rhs_solve_count
            << ',' << d.scalar_factorization_count << ','
            << d.scalar_rhs_solve_count << ',' << d.positivity_failure << ',';
  for (int j = 0; j < num_constituents; ++j) {
    const auto &x = d.field_clipping[j];
    if (method != "imex") {
      emit_real(x.interstitial_raw_minimum);
      emit_real(x.cloudborne_raw_minimum);
    }
    emit_real(x.interstitial_weighted_deficit);
    emit_real(x.cloudborne_weighted_deficit);
    std::cout << x.interstitial_count << ',' << x.cloudborne_count << ','
              << x.interstitial_roundoff_negative_count << ','
              << x.interstitial_material_negative_count << ','
              << x.cloudborne_roundoff_negative_count << ','
              << x.cloudborne_material_negative_count << ',';
    emit_real(d.initial_ledger[j]);
    emit_real(d.final_ledger[j]);
    emit_real(d.clipping_ledger_delta[j]);
    emit_real(d.total_ledger_delta[j]);
    emit_real(d.operator_ledger_delta[j]);
  }
  std::cout << '\n';
}

#ifndef NDROP_LEGACY_BASELINE
bool same_state(const State &a, const State &b) {
  if (a.nlev() != b.nlev()) return false;
  const auto same = [](const std::vector<Real> &x,
                       const std::vector<Real> &y) {
    if (x.size() != y.size()) return false;
    for (std::size_t i = 0; i < x.size(); ++i)
      if (bits(x[i]) != bits(y[i])) return false;
    return true;
  };
  if (!same(a.droplet_number, b.droplet_number)) return false;
  for (int j = 0; j < num_constituents; ++j)
    if (!same(a.interstitial[j], b.interstitial[j]) ||
        !same(a.cloudborne[j], b.cloudborne[j])) return false;
  return true;
}

void test_diagnostics() {
  const Fixture fixture;
  for (const std::string method : {"forward_euler", "imex", "backward_euler"}) {
    AdvanceControl with_detail;
    with_detail.collect_replay_diagnostics = true;
    const auto enabled = advance(method, fixture, with_detail);
    const auto existing = advance(method, fixture, {});
    check(enabled.ok && existing.ok, method + " successful fixture");
    check(same_state(enabled.state, existing.state), method + " bitwise state");
    check(enabled.message == existing.message &&
              enabled.diagnostics.nsubmix == existing.diagnostics.nsubmix &&
              bits(enabled.diagnostics.normalized_residual) ==
                  bits(existing.diagnostics.normalized_residual) &&
              bits(enabled.diagnostics.max_normalized_ledger_error) ==
                  bits(existing.diagnostics.max_normalized_ledger_error),
          method + " unchanged existing decisions");
    check(enabled.diagnostics.replay_detail.has_value(), method + " detail");
    const auto &x = *enabled.diagnostics.replay_detail;
    check(x.status == ReplayDiagnosticStatus::complete &&
              x.completed_substeps == x.expected_substeps &&
              x.expected_substeps == (method == "forward_euler" ? 3U : 1U),
          method + " complete substep evidence");
    Real physical_max = 0, matrix_max = 0;
    for (int i = 0; i < 51; ++i) {
      check(x.raw_minimum_valid[i] && x.negative_counts_valid[i] &&
                x.physical_residual[i].has_value(),
            method + " complete field " + std::to_string(i));
      physical_max = std::max(physical_max, x.physical_residual[i]->normalized);
      if (method == "forward_euler") {
        check(!x.matrix_residual[i] && x.fe_clipping_count_valid[i] &&
                  x.fe_weighted_deficit_valid[i], "FE phase applicability");
        check(x.physical_residual[i]->substep < 3, "FE worst substep");
      } else {
        check(x.matrix_residual[i].has_value() &&
                  !x.fe_clipping_count_valid[i] &&
                  !x.fe_weighted_deficit_valid[i],
              method + " implicit phase applicability");
        matrix_max = std::max(matrix_max, x.matrix_residual[i]->normalized);
        check(x.matrix_residual[i]->substep == 0, "implicit substep index");
      }
    }
    check(bits(physical_max) ==
              bits(enabled.diagnostics.normalized_physical_residual),
          method + " physical max reproduces original");
    if (method != "forward_euler")
      check(bits(matrix_max) ==
                bits(enabled.diagnostics.normalized_matrix_residual),
            method + " matrix max reproduces original");
    for (int j = 0; j < num_constituents; ++j)
      check(x.total_ledger_valid[j] && x.operator_ledger_valid[j] &&
                x.clipping_ledger_valid[j], method + " ledger masks");
    if (method == "imex")
      check(std::isfinite(enabled.diagnostics.field_clipping[0]
                              .interstitial_raw_minimum) &&
                std::isinf(existing.diagnostics.field_clipping[0]
                               .interstitial_raw_minimum),
            "IMEX opt-in raw minimum before projection");
    const auto measurement =
        production_replay_v1::serialize_diagnostic_measurement(enabled, method);
    check(measurement.find("\"completion\"") != std::string::npos &&
              measurement.find("\"physical_residual\"") != std::string::npos,
          method + " serialized measurement");
  }
  AdvanceControl lean;
  lean.collect_detailed_diagnostics = false;
  lean.collect_replay_diagnostics = true;
  const auto disabled = advance("forward_euler", fixture, lean);
  check(disabled.ok && !disabled.diagnostics.replay_detail,
        "lean primary path leaves detail disengaged");
  check(!advance("forward_euler", fixture, {}).diagnostics.replay_detail,
        "default detail disengaged");

  AdvanceControl with_detail;
  with_detail.collect_replay_diagnostics = true;
  Fixture stationary;
  stationary.op.interface_exchange[1] = 0;
  stationary.op.interface_exchange[2] = 0;
  stationary.op.first_level_controller_rate = 0;
  for (int mode = 0; mode < num_modes; ++mode)
    for (int k = 0; k < 3; ++k) {
      stationary.op.number_activation_rate[mode][k] = 0;
      stationary.op.mass_activation_rate[mode][k] = 0;
    }
  const auto tied = advance("forward_euler", stationary, with_detail);
  check(tied.ok, "stationary FE tie fixture");
  for (const auto &value : tied.diagnostics.replay_detail->physical_residual)
    check(value && value->normalized == 0 && value->substep == 0,
          "FE exact ties keep earliest substep");

  Fixture bad;
  bad.op.pressure_weight[0] = -1;
  const auto invalid = advance("forward_euler", bad, with_detail);
  check(!invalid.ok && invalid.diagnostics.replay_detail &&
            invalid.diagnostics.replay_detail->status ==
                ReplayDiagnosticStatus::not_started &&
            invalid.diagnostics.replay_detail->failure_reason ==
                ReplayDiagnosticReason::diagnostic_incomplete &&
            !invalid.diagnostics.replay_detail->physical_residual[0],
        "input failure has no fabricated diagnostic");

  State overflowing_state(2);
  OperatorData overflowing_op(2);
  overflowing_op.pressure_weight[0] = std::numeric_limits<Real>::min();
  overflowing_op.pressure_weight[1] = std::numeric_limits<Real>::min();
  overflowing_op.cloud_fraction[0] = 0;
  overflowing_state.interstitial[0][0] = std::numeric_limits<Real>::max();
  overflowing_state.cloudborne[0][0] = std::numeric_limits<Real>::max();
  const auto projected = advance_imex(overflowing_state, overflowing_op, 1,
                                      with_detail);
  check(!projected.ok && projected.message ==
            "zero-cloud aerosol resuspension overflowed" &&
            projected.diagnostics.replay_detail->status ==
                ReplayDiagnosticStatus::partial &&
            projected.diagnostics.replay_detail->completed_substeps == 1 &&
            projected.diagnostics.replay_detail->failure_stage ==
                ReplayDiagnosticStage::projection &&
            !projected.diagnostics.replay_detail->total_ledger_valid[0],
        "post-solve projection failure retains partial evidence: " +
            projected.message);
  State negative_state(2);
  OperatorData negative_op(2);
  negative_state.interstitial[0] = {0.001, 100};
  negative_op.interface_exchange[1] = 0.5;
  negative_op.number_activation_rate[0][0] = 0.5;
  const auto negative = advance_imex(negative_state, negative_op, 100,
                                     with_detail);
  check(!negative.ok && negative.message ==
            "material negative value in IMEX solution" &&
            negative.diagnostics.replay_detail &&
            negative.diagnostics.replay_detail->status ==
                ReplayDiagnosticStatus::partial &&
            negative.diagnostics.replay_detail->completed_substeps == 1 &&
            negative.diagnostics.replay_detail->failure_reason ==
                ReplayDiagnosticReason::raw_material_negative &&
            negative.diagnostics.replay_detail->negative_counts_valid[1] &&
            !negative.diagnostics.replay_detail->total_ledger_valid[0],
        "material-negative failure retains finite computed fields only");
}

void test_reference_classifier() {
  using namespace production_replay_v1;
  Fixture f;
  const auto actual = run_reference_self_check(true, true, false, f.state,
                                                f.op, f.dt);
  check(actual.outcome == "completed" && actual.passed == true &&
            actual.primary_completed && actual.tight_completed,
        "actual two-tolerance reference succeeds");
  check(serialize_reference_self_check(actual).find(
            "\"outcome\":\"completed\"") != std::string::npos,
        "actual reference serializes");
  const auto candidate_failure = classify_reference_self_check(
      false, true, false, f.state, f.op, f.dt, {}, {});
  check(candidate_failure.outcome == "not_attempted" &&
            candidate_failure.failure_reason == "candidate_failed_first" &&
            !candidate_failure.maximum_normalized_disagreement,
        "candidate-first reference exclusion");
  check(serialize_reference_self_check(candidate_failure).find(
            "\"candidate_failed_first\"") != std::string::npos,
        "not-attempted reference serializes");
  const auto sentinel = classify_reference_self_check(
      true, true, true, f.state, f.op, f.dt, {}, {});
  check(sentinel.applicability == "not_applicable" &&
            sentinel.failure_reason == "sentinel_excluded", "sentinel exclusion");
  ReferenceResult complete;
  complete.ok = true;
  complete.message = "ok";
  complete.state = f.state;
  complete.accepted_steps = 1;
  complete.function_evaluations = 7;
  const auto missing = classify_reference_self_check(
      true, true, false, f.state, f.op, f.dt, complete, {});
  check(missing.outcome == "failed" && missing.primary_completed &&
            !missing.tight_completed && missing.failure_stage == "tight_setup" &&
            missing.failure_reason == "incomplete_reference_pair" &&
            !missing.passed, "incomplete pair retained without metric");
  const auto fabricated_complete = classify_reference_self_check(
      true, true, false, f.state, f.op, f.dt, complete, complete);
  check(fabricated_complete.outcome == "completed" &&
            fabricated_complete.maximum_normalized_disagreement == 0,
        "fabricated identical pair classified for test only");
  for (const auto &[message, stage, reason] : {
           std::tuple{"reference exceeded maximum attempted steps",
                      "primary_advance", "step_limit_exhausted"},
           std::tuple{"reference step size underflowed",
                      "primary_advance", "step_underflow"},
          std::tuple{"reference stage 1 is nonfinite",
                      "primary_advance", "nonfinite_stage"},
          std::tuple{"reference embedded solution is nonfinite",
                      "primary_advance", "nonfinite_stage"},
          std::tuple{"reference error norm is nonfinite",
                      "primary_advance", "nonfinite_stage"},
          std::tuple{"reference final state is nonfinite",
                      "primary_advance", "nonfinite_stage"},
          std::tuple{"reference zero-cloud postprocess overflowed",
                      "primary_advance", "nonfinite_stage"},
          std::tuple{"reference options must be finite and positive",
                      "primary_setup", "incomplete_reference_pair"},
           std::tuple{"reference operator contains a nonfinite coefficient",
                      "primary_setup", "nonfinite_stage"},
           std::tuple{"unrecognized fabricated failure", "primary_advance",
                      "incomplete_reference_pair"}}) {
    ReferenceResult failed;
    failed.message = message;
    const auto classified = classify_reference_self_check(
        true, true, false, f.state, f.op, f.dt, failed, {});
    check(classified.failure_stage == stage &&
              classified.failure_reason == reason && !classified.passed,
          std::string("fabricated classifier case: ") + message);
    check(serialize_reference_self_check(classified).find(reason) !=
              std::string::npos,
          "fabricated source-bound failure serializes");
  }
  for (int n = 1; n <= 7; ++n) {
    ReferenceResult failed;
    failed.message = "reference stage " + std::to_string(n) +
                     " is nonfinite";
    const auto primary = classify_reference_self_check(
        true, true, false, f.state, f.op, f.dt, failed, {});
    const auto tight = classify_reference_self_check(
        true, true, false, f.state, f.op, f.dt, complete, failed);
    check(primary.failure_stage == "primary_advance" &&
              tight.failure_stage == "tight_advance" &&
              primary.failure_reason == "nonfinite_stage" &&
              tight.failure_reason == "nonfinite_stage" &&
              tight.primary_completed && !tight.tight_completed,
          "source-bound stage N classifier");
  }
  for (const auto &[message, reason] : {
           std::pair{"reference exceeded maximum attempted steps",
                     "step_limit_exhausted"},
           std::pair{"reference step size underflowed", "step_underflow"},
           std::pair{"reference operator contains a nonfinite coefficient",
                     "nonfinite_stage"},
           std::pair{"unrecognized fabricated failure",
                     "incomplete_reference_pair"}}) {
    ReferenceResult failed;
    failed.message = message;
    const auto tight = classify_reference_self_check(
        true, true, false, f.state, f.op, f.dt, complete, failed);
    check(tight.failure_reason == reason && tight.primary_completed &&
              !tight.tight_completed, "tight-call failure classification");
  }
  ReferenceResult discrepant = complete;
  discrepant.state.droplet_number[0] *= 2;
  const auto above = classify_reference_self_check(
      true, true, false, f.state, f.op, f.dt, complete, discrepant);
  check(above.outcome == "completed" && above.passed == false &&
            above.maximum_normalized_disagreement > above.threshold,
        "finite above-gate reference pair");
  ReferenceResult nonfinite_pair = complete;
  nonfinite_pair.state.droplet_number[0] =
      std::numeric_limits<Real>::max();
  nonfinite_pair.state.droplet_number[1] =
      std::numeric_limits<Real>::max();
  const auto comparison_failure = classify_reference_self_check(
      true, true, false, f.state, f.op, f.dt, complete, nonfinite_pair);
  check(comparison_failure.failure_stage == "pair_comparison" &&
            comparison_failure.failure_reason == "comparison_nonfinite" &&
            !comparison_failure.maximum_normalized_disagreement,
        "finite states with overflowing pair norm fail closed");
  Fixture invalid_input;
  invalid_input.op.pressure_weight[0] = -1;
  std::string validation_message;
  check(!validate(invalid_input.state, invalid_input.op, invalid_input.dt,
                  validation_message), "input fixture is invalid");
  ReferenceResult validation_failure;
  validation_failure.message =
      "reference input validation failed: " + validation_message;
  const auto input_classified = classify_reference_self_check(
      true, true, false, invalid_input.state, invalid_input.op,
      invalid_input.dt, validation_failure, {});
  check(input_classified.failure_stage == "primary_setup" &&
            input_classified.failure_reason == "incomplete_reference_pair",
        "exact input-validation suffix");
  ReferenceResult inconsistent = complete;
  inconsistent.message = "reference stage 1 is nonfinite";
  const auto rejected = classify_reference_self_check(
      true, true, false, f.state, f.op, f.dt, inconsistent, {});
  check(rejected.failure_reason == "incomplete_reference_pair" &&
            !rejected.primary_completed, "inconsistent success fails closed");
  inconsistent = complete;
  inconsistent.accepted_steps = 0;
  check(classify_reference_self_check(true, true, false, f.state, f.op, f.dt,
                                      inconsistent, {})
            .failure_reason == "incomplete_reference_pair",
        "inconsistent counter fails closed");
  inconsistent = complete;
  inconsistent.function_evaluations = 1;
  check(classify_reference_self_check(true, true, false, f.state, f.op, f.dt,
                                      inconsistent, {})
            .failure_reason == "incomplete_reference_pair",
        "impossible evaluation count fails closed");
  inconsistent = complete;
  inconsistent.maximum_error_ratio = -1;
  check(classify_reference_self_check(true, true, false, f.state, f.op, f.dt,
                                      inconsistent, {})
            .failure_reason == "incomplete_reference_pair",
        "negative accepted-step error ratio fails closed");
}
#endif
} // namespace

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  try {
    const Fixture f;
    if (argc == 2 &&
        (std::string(argv[1]) == "--baseline-digest" ||
         std::string(argv[1]) == "--baseline-lean-digest")) {
      AdvanceControl control;
      if (std::string(argv[1]) == "--baseline-lean-digest")
        control.collect_detailed_diagnostics = false;
      for (const std::string method : {"forward_euler", "imex", "backward_euler"})
        emit_snapshot(advance(method, f, control), method);
    } else {
#ifndef NDROP_LEGACY_BASELINE
      if (argc == 3 && std::string(argv[1]) == "--emit-measurement") {
        const std::string method = argv[2];
        check(method == "forward_euler" || method == "imex" ||
                  method == "backward_euler", "unknown fixture method");
        AdvanceControl control;
        control.collect_replay_diagnostics = true;
        const auto result = advance(method, f, control);
        check(result.ok, "actual fixture advance failed");
        std::cout << production_replay_v1::serialize_diagnostic_measurement(
                         result, method) << '\n';
      } else if (argc == 2 &&
                 std::string(argv[1]) == "--emit-reference") {
        const auto check = production_replay_v1::run_reference_self_check(
            true, true, false, f.state, f.op, f.dt);
        std::cout << production_replay_v1::serialize_reference_self_check(
                         check) << '\n';
      } else {
      test_diagnostics();
      test_reference_classifier();
      std::cout << "OI-5 offline unit tests passed\n";
      }
#endif
    }
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    Kokkos::finalize();
    return 1;
  }
  Kokkos::finalize();
  return 0;
}
