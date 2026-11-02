#include "ndrop_production_replay.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace ndrop_test::production_replay_v1 {
namespace {

bool finite_complete_state(const State &state, std::size_t nlev) {
  if (state.nlev() != nlev || state.droplet_number.size() != nlev) return false;
  const auto valid = [nlev](const std::vector<Real> &field) {
    return field.size() == nlev &&
           std::all_of(field.begin(), field.end(),
                       [](Real x) { return std::isfinite(x); });
  };
  if (!valid(state.droplet_number)) return false;
  for (int j = 0; j < num_constituents; ++j) {
    if (!valid(state.interstitial[j]) || !valid(state.cloudborne[j]))
      return false;
  }
  return true;
}

bool completed(const ReferenceResult &result, std::size_t nlev) {
  return result.ok && result.message == "ok" &&
         finite_complete_state(result.state, nlev) &&
         result.accepted_steps > 0 && result.rejected_steps >= 0 &&
         std::int64_t(result.function_evaluations) ==
             7 * (std::int64_t(result.accepted_steps) +
                  std::int64_t(result.rejected_steps)) &&
         std::isfinite(result.maximum_error_ratio) &&
         result.maximum_error_ratio >= 0 && result.maximum_error_ratio <= 1;
}

std::pair<std::string, std::string>
failure(const ReferenceResult &result, const State &initial,
        const OperatorData &data, Real dt, const std::string &which) {
  const std::string &message = result.message;
  if (result.ok || result.accepted_steps < 0 ||
      result.rejected_steps < 0 || result.function_evaluations < 0 ||
      !std::isfinite(result.maximum_error_ratio))
    return {which + "_advance", "incomplete_reference_pair"};
  std::string validation_message;
  if (!validate(initial, data, dt, validation_message) &&
      message == "reference input validation failed: " + validation_message)
    return {which + "_setup", "incomplete_reference_pair"};
  if (message == "reference options must be finite and positive")
    return {which + "_setup", "incomplete_reference_pair"};
  if (message == "reference operator contains a nonfinite coefficient")
    return {which + "_setup", "nonfinite_stage"};
  if (message == "reference exceeded maximum attempted steps")
    return {which + "_advance", "step_limit_exhausted"};
  if (message == "reference step size underflowed")
    return {which + "_advance", "step_underflow"};
  for (int n = 1; n <= 7; ++n) {
    if (message == "reference stage " + std::to_string(n) + " is nonfinite")
      return {which + "_advance", "nonfinite_stage"};
  }
  if (message == "reference embedded solution is nonfinite" ||
      message == "reference error norm is nonfinite" ||
      message == "reference final state is nonfinite" ||
      message == "reference zero-cloud postprocess overflowed")
    return {which + "_advance", "nonfinite_stage"};
  return {which + "_advance", "incomplete_reference_pair"};
}

bool field_error_max(const std::vector<Real> &actual,
                     const std::vector<Real> &reference,
                     const std::vector<Real> &initial,
                     const std::vector<Real> &weight, Real &maximum) {
  Real error_l1 = 0, ref_l1 = 0, initial_l1 = 0;
  Real error_linf = 0, ref_linf = 0, initial_linf = 0;
  for (std::size_t k = 0; k < weight.size(); ++k) {
    if (!std::isfinite(weight[k]) || weight[k] <= 0 ||
        !std::isfinite(initial[k])) return false;
    const Real difference = std::abs(actual[k] - reference[k]);
    error_l1 += weight[k] * difference;
    ref_l1 += weight[k] * std::abs(reference[k]);
    initial_l1 += weight[k] * std::abs(initial[k]);
    error_linf = std::max(error_linf, difference);
    ref_linf = std::max(ref_linf, std::abs(reference[k]));
    initial_linf = std::max(initial_linf, std::abs(initial[k]));
    if (!std::isfinite(error_l1) || !std::isfinite(ref_l1) ||
        !std::isfinite(initial_l1) || !std::isfinite(error_linf) ||
        !std::isfinite(ref_linf) || !std::isfinite(initial_linf)) return false;
  }
  const Real l1_scale = std::max(ref_l1, initial_l1);
  const Real linf_scale = std::max(ref_linf, initial_linf);
  if (l1_scale > 0) error_l1 /= l1_scale;
  if (linf_scale > 0) error_linf /= linf_scale;
  if (!std::isfinite(error_l1) || !std::isfinite(error_linf)) return false;
  maximum = std::max({maximum, error_l1, error_linf});
  return true;
}

std::optional<Real> maximum_disagreement(const State &primary,
                                         const State &tight,
                                         const State &initial,
                                         const OperatorData &data) {
  Real maximum = 0;
  const auto &w = data.pressure_weight;
  if (w.size() != initial.nlev() ||
      !field_error_max(primary.droplet_number, tight.droplet_number,
                       initial.droplet_number, w, maximum)) return std::nullopt;
  for (int j = 0; j < num_constituents; ++j) {
    if (!field_error_max(primary.interstitial[j], tight.interstitial[j],
                         initial.interstitial[j], w, maximum) ||
        !field_error_max(primary.cloudborne[j], tight.cloudborne[j],
                         initial.cloudborne[j], w, maximum)) return std::nullopt;
  }
  return maximum;
}

} // namespace

ReferenceSelfCheck classify_reference_self_check(
    bool candidate_ok, bool replayable, bool stress_sentinel,
    const State &initial, const OperatorData &data, Real dt,
    const std::optional<ReferenceResult> &primary,
    const std::optional<ReferenceResult> &tight) {
  ReferenceSelfCheck check;
  if (stress_sentinel || !replayable) {
    check.applicability = "not_applicable";
    check.failure_reason = stress_sentinel ? "sentinel_excluded"
                                           : "record_not_replayable";
    check.message = *check.failure_reason;
    return check;
  }
  if (!candidate_ok) {
    check.failure_reason = "candidate_failed_first";
    check.message = *check.failure_reason;
    return check;
  }
  check.attempt_status = "attempted";
  check.outcome = "failed";
  if (!primary) {
    check.failure_stage = "primary_setup";
    check.failure_reason = "incomplete_reference_pair";
    check.message = "missing primary reference result";
    return check;
  }
  if (!completed(*primary, initial.nlev())) {
    const auto [stage, reason] = failure(*primary, initial, data, dt, "primary");
    check.failure_stage = stage;
    check.failure_reason = reason;
    check.message = primary->message;
    return check;
  }
  check.primary_completed = true;
  if (!tight) {
    check.failure_stage = "tight_setup";
    check.failure_reason = "incomplete_reference_pair";
    check.message = "missing tight reference result";
    return check;
  }
  if (!completed(*tight, initial.nlev())) {
    const auto [stage, reason] = failure(*tight, initial, data, dt, "tight");
    check.failure_stage = stage;
    check.failure_reason = reason;
    check.message = tight->message;
    return check;
  }
  check.tight_completed = true;
  const auto difference = maximum_disagreement(primary->state, tight->state,
                                                initial, data);
  if (!difference) {
    check.failure_stage = "pair_comparison";
    check.failure_reason = "comparison_nonfinite";
    check.message = "reference pair comparison is nonfinite";
    return check;
  }
  check.outcome = "completed";
  check.maximum_normalized_disagreement = *difference;
  check.passed = *difference <= check.threshold;
  check.message = *check.passed ? "ok" : "reference pair exceeds tolerance";
  return check;
}

ReferenceSelfCheck run_reference_self_check(bool candidate_ok,
                                            bool replayable,
                                            bool stress_sentinel,
                                            const State &initial,
                                            const OperatorData &data, Real dt) {
  if (!candidate_ok || !replayable || stress_sentinel)
    return classify_reference_self_check(candidate_ok, replayable,
                                         stress_sentinel, initial, data, dt,
                                         std::nullopt, std::nullopt);
  ReferenceResult primary;
  try {
    primary = advance_reference_dopri54(initial, data, dt);
  } catch (...) {
    primary.message = "reference primary call threw an exception";
    return classify_reference_self_check(true, true, false, initial, data, dt,
                                         primary, std::nullopt);
  }
  if (!completed(primary, initial.nlev()))
    return classify_reference_self_check(true, true, false, initial, data, dt,
                                         primary, std::nullopt);
  ReferenceOptions tighter;
  tighter.relative_tolerance = 1e-13;
  ReferenceResult tight;
  try {
    tight = advance_reference_dopri54(initial, data, dt, tighter);
  } catch (...) {
    tight.message = "reference tight call threw an exception";
  }
  return classify_reference_self_check(true, true, false, initial, data, dt,
                                       primary, tight);
}

} // namespace ndrop_test::production_replay_v1
