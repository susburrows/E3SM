#include "ndrop_integrator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ndrop_test::AdvanceResult;
using ndrop_test::ForwardEulerOptions;
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::State;
using ndrop_test::num_constituents;
using ndrop_test::num_modes;

// This oracle intentionally duplicates the expected MAM constituent layout.
// It must not call the harness mapping helpers, so a mapping permutation in
// the implementation causes the analytic tests to fail.
constexpr std::array<int, num_constituents> expected_mode = {
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3};
constexpr std::array<bool, num_constituents> expected_number = {
    true,  false, false, false, false, false, false, false, true,
    false, false, false, false, true,  false, false, false, false,
    false, false, false, true,  false, false, false};

int passed = 0;
int failed = 0;

void check(const bool condition, const std::string &message) {
  if (!condition) {
    ++failed;
    std::cerr << "FAIL: " << message << '\n';
  } else {
    ++passed;
  }
}

bool near(const Real lhs, const Real rhs, const Real scale = 1) {
  const Real tolerance =
      200 * std::numeric_limits<Real>::epsilon() *
      std::max({std::abs(lhs), std::abs(rhs), scale});
  return std::abs(lhs - rhs) <= tolerance;
}

State make_state(const std::size_t nlev, const bool uniform = false) {
  State state(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    state.droplet_number[k] = uniform ? 4 : 4 + 0.02 * k;
    for (int j = 0; j < num_constituents; ++j) {
      state.interstitial[j][k] =
          uniform ? 1 + 0.01 * j : 1 + 0.01 * j + 0.003 * k;
      state.cloudborne[j][k] =
          uniform ? 0.2 + 0.002 * j : 0.2 + 0.002 * j + 0.001 * k;
    }
  }
  return state;
}

OperatorData make_operator(const std::size_t nlev,
                           const Real interface_exchange,
                           const Real activation_rate,
                           const Real cloud_base = 1) {
  OperatorData op(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    op.pressure_weight[k] = 0.8 + 0.01 * k;
    op.cloud_fraction[k] =
        std::clamp(cloud_base - 0.003 * static_cast<Real>(k), Real(0.1),
                   Real(1));
  }
  for (std::size_t interface = 1; interface < nlev; ++interface) {
    // Keep cell rates close to interface_exchange even with nonuniform weight.
    op.interface_exchange[interface] =
        interface_exchange *
        std::min(op.pressure_weight[interface - 1],
                 op.pressure_weight[interface]);
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    std::fill(op.number_activation_rate[mode].begin(),
              op.number_activation_rate[mode].end(), activation_rate);
    std::fill(op.mass_activation_rate[mode].begin(),
              op.mass_activation_rate[mode].end(),
              Real(0.5) * activation_rate);
    // A zero bottom rate is the only limiter-inactive input under the
    // production no-flux bottom coefficient.
    op.number_activation_rate[mode].back() = 0;
    op.mass_activation_rate[mode].back() = 0;
  }
  return op;
}

void check_state_equal(const State &actual, const State &expected,
                       const std::string &label) {
  check(actual.nlev() == expected.nlev(), label + " nlev");
  for (std::size_t k = 0; k < expected.nlev(); ++k) {
    check(near(actual.droplet_number[k], expected.droplet_number[k]),
          label + " droplet level " + std::to_string(k));
    for (int j = 0; j < num_constituents; ++j) {
      check(near(actual.interstitial[j][k], expected.interstitial[j][k]),
            label + " interstitial");
      check(near(actual.cloudborne[j][k], expected.cloudborne[j][k]),
            label + " cloudborne");
    }
  }
}

void check_ledgers(const AdvanceResult &result, const std::string &label) {
  for (int j = 0; j < num_constituents; ++j) {
    check(near(result.diagnostics.final_ledger[j],
               result.diagnostics.initial_ledger[j]),
          label + " ledger " + std::to_string(j));
  }
}

void test_invalid_inputs() {
  State empty;
  OperatorData empty_op;
  auto result = ndrop_test::advance_forward_euler(empty, empty_op, 1);
  check(!result.ok && result.message == "nlev must be positive",
        "empty state fails loudly");

  State state = make_state(2);
  OperatorData op = make_operator(2, 0.1, 0.01);
  result = ndrop_test::advance_forward_euler(state, op, 0);
  check(!result.ok && result.message == "dt must be finite and positive",
        "nonpositive dt fails loudly");

  op.interface_exchange.front() = 0.1;
  result = ndrop_test::advance_forward_euler(state, op, 1);
  check(!result.ok &&
            result.message ==
                "top and bottom interface exchange coefficients must be zero",
        "nonzero boundary flux fails loudly");

  op = make_operator(2, 0.1, 0.01);
  op.cloud_fraction[0] = 1.1;
  result = ndrop_test::advance_forward_euler(state, op, 1);
  check(!result.ok &&
            result.message == "cloud fractions must be finite and in [0,1]",
        "invalid cloud fraction fails loudly");

  op = make_operator(2, 0.1, 0.01);
  state.interstitial[3][0] = -1;
  result = ndrop_test::advance_forward_euler(state, op, 1);
  check(!result.ok &&
            result.message == "aerosol state must be finite and nonnegative",
        "negative input state fails loudly");

  state = make_state(2);
  op = make_operator(2, std::numeric_limits<Real>::max(), 0);
  op.pressure_weight[0] = std::numeric_limits<Real>::min();
  result = ndrop_test::advance_forward_euler(state, op, 1);
  check(!result.ok &&
            result.message == "derived mixing rates must be finite",
        "overflowed derived rate fails loudly");

  op = make_operator(2, 1e100, 0);
  result = ndrop_test::advance_forward_euler(state, op, 1e100);
  check(!result.ok &&
            result.message ==
                "calculated substep count is outside the supported range",
        "overflowed substep count fails loudly");
}

void test_identity_and_level_counts() {
  for (const std::size_t nlev : {1u, 2u, 3u, 127u, 128u}) {
    State state = make_state(nlev);
    OperatorData op(nlev);
    const auto result = ndrop_test::advance_forward_euler(state, op, 100);
    check(result.ok, "identity advance succeeds");
    check(result.diagnostics.nsubmix == 1, "identity uses one substep");
    check(result.diagnostics.clipping_count == 0,
          "identity does not clip");
    check_state_equal(result.state, state,
                      "identity nlev=" + std::to_string(nlev));
    check_ledgers(result, "identity nlev=" + std::to_string(nlev));
  }
}

void test_uniform_full_overlap() {
  State state = make_state(128, true);
  OperatorData op = make_operator(128, 0.7, 0, 1);
  std::fill(op.cloud_fraction.begin(), op.cloud_fraction.end(), 1);
  const auto result = ndrop_test::advance_forward_euler(state, op, 10);
  check(result.ok, "uniform full-overlap advance succeeds");
  check(result.diagnostics.clipping_count == 0,
        "uniform full-overlap does not clip");
  check_state_equal(result.state, state, "uniform full-overlap");
  check_ledgers(result, "uniform full-overlap");
}

void test_substeps_limiter_and_parity() {
  struct Case {
    Real exchange;
    int expected_substeps;
  };
  for (const Case test : {Case{0.01, 1}, Case{0.5, 2}, Case{1.0, 3},
                          Case{2.0, 5}}) {
    State state = make_state(128);
    OperatorData op = make_operator(128, test.exchange, 0.001, 1);
    std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1);
    for (std::size_t interface = 1; interface < 128; ++interface) {
      op.interface_exchange[interface] = test.exchange;
    }
    const auto result = ndrop_test::advance_forward_euler(state, op, 1);
    check(result.ok, "substep case succeeds");
    check(result.diagnostics.nsubmix == test.expected_substeps,
          "expected nsubmix=" + std::to_string(test.expected_substeps));
    check(result.diagnostics.final_buffer ==
              (test.expected_substeps % 2 == 0 ? 0 : 1),
          "buffer parity for nsubmix=" +
              std::to_string(test.expected_substeps));
    check(result.diagnostics.clipping_count == 0,
          "adaptive FE does not clip in substep matrix");
    check_ledgers(result,
                  "substeps=" + std::to_string(test.expected_substeps));
  }

  State state = make_state(3);
  OperatorData inactive = make_operator(3, 0.5, 0.1);
  auto result = ndrop_test::advance_forward_euler(state, inactive, 1);
  check(result.ok && result.diagnostics.activation_limiter_count == 0,
        "limiter-inactive case performs no limiting");

  OperatorData active = make_operator(3, 0.5, 10);
  result = ndrop_test::advance_forward_euler(state, active, 1);
  check(result.ok && result.diagnostics.activation_limiter_count > 8,
        "large activation rates engage limiter");
  check_ledgers(result, "limiter-active");
}

void test_fractional_overlap_ledger() {
  State state = make_state(127);
  OperatorData op = make_operator(127, 0.8, 0.02, 0.85);
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    op.cloud_fraction[k] =
        k % 4 == 0 ? 0 : (k % 4 == 1 ? 0.2 : (k % 4 == 2 ? 0.65 : 1));
  }
  const auto result = ndrop_test::advance_forward_euler(state, op, 3);
  check(result.ok, "fractional-overlap advance succeeds");
  check(result.diagnostics.clipping_count == 0,
        "fractional-overlap adaptive FE does not clip");
  check_ledgers(result, "fractional-overlap");
}

Real next_unit_random(unsigned long long &state) {
  state = state * 6364136223846793005ULL + 1442695040888963407ULL;
  return static_cast<Real>(state >> 11) /
         static_cast<Real>(1ULL << 53);
}

void test_adaptive_subcycling_avoids_clipping() {
  unsigned long long seed = 0x4645414441505456ULL;
  for (const std::size_t nlev : {1u, 2u, 3u, 127u, 128u}) {
    for (int sample = 0; sample < 4; ++sample) {
      State state(nlev);
      OperatorData op(nlev);
      for (std::size_t k = 0; k < nlev; ++k) {
        op.pressure_weight[k] = 0.2 + 2 * next_unit_random(seed);
        const Real selector = next_unit_random(seed);
        op.cloud_fraction[k] =
            selector < 0.15 ? 0 : (selector > 0.85 ? 1 : selector);
        state.droplet_number[k] = 0.01 + 20 * next_unit_random(seed);
        for (int j = 0; j < num_constituents; ++j) {
          state.interstitial[j][k] =
              0.001 + 5 * next_unit_random(seed);
          state.cloudborne[j][k] =
              0.001 + 2 * next_unit_random(seed);
        }
      }
      for (std::size_t interface = 1; interface < nlev; ++interface) {
        const Real rate = std::pow(10, -3 + 4 * next_unit_random(seed));
        op.interface_exchange[interface] =
            rate * std::min(op.pressure_weight[interface - 1],
                            op.pressure_weight[interface]);
      }
      for (int mode = 0; mode < num_modes; ++mode) {
        for (std::size_t k = 0; k + 1 < nlev; ++k) {
          const Real rate_below =
              op.interface_exchange[k + 1] / op.pressure_weight[k];
          op.number_activation_rate[mode][k] =
              2 * next_unit_random(seed) * rate_below;
          op.mass_activation_rate[mode][k] =
              2 * next_unit_random(seed) * rate_below;
        }
      }
      const Real dt = 0.1 + 3 * next_unit_random(seed);
      const auto result = ndrop_test::advance_forward_euler(state, op, dt);
      check(result.ok, "adaptive randomized FE column succeeds");
      check(result.diagnostics.clipping_count == 0,
            "production-adaptive FE avoids clipping");
      check(result.diagnostics.raw_minimum >= 0,
            "production-adaptive FE raw state remains nonnegative");
    }
  }
}

void test_forced_step_clipping_accounting() {
  OperatorData op = make_operator(3, 2, 0, 1);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1);
  std::fill(op.cloud_fraction.begin(), op.cloud_fraction.end(), 1);
  op.interface_exchange[1] = 2;
  op.interface_exchange[2] = 2;

  ForwardEulerOptions raw_options;
  raw_options.forced_substeps = 1;
  raw_options.apply_clipping = false;
  ForwardEulerOptions clipped_options = raw_options;
  clipped_options.apply_clipping = true;

  for (const bool clip_interstitial : {true, false}) {
    State state = make_state(3, true);
    state.droplet_number =
        clip_interstitial ? std::vector<Real>{10, 0.1, 0.1}
                          : std::vector<Real>{1, 1, 1};
    for (int j = 0; j < num_constituents; ++j) {
      if (clip_interstitial) {
        state.interstitial[j] = {10, 0.1, 0.1};
        state.cloudborne[j] = {1, 1, 1};
      } else {
        state.interstitial[j] = {1, 1, 1};
        state.cloudborne[j] = {10, 0.1, 0.1};
      }
    }

    const auto raw =
        ndrop_test::advance_forward_euler(state, op, 1, raw_options);
    check(raw.ok && raw.diagnostics.raw_minimum < 0,
          "forced unstable raw step records negative value");
    check(raw.diagnostics.clipping_count == 0,
          "unclipped forced step does not hide negative value");

    const auto clipped =
        ndrop_test::advance_forward_euler(state, op, 1, clipped_options);
    const int expected_total = num_constituents +
                               (clip_interstitial ? 1 : 0);
    check(clipped.ok && clipped.diagnostics.clipping_count == expected_total,
          "forced unstable step has exact clipping count");
    check(clipped.diagnostics.droplet_clipping_count ==
              (clip_interstitial ? 1 : 0),
          "droplet clipping count is exact");
    check(near(clipped.diagnostics.droplet_weighted_deficit,
               clip_interstitial ? Real(9.8) : Real(0)),
          "droplet clipping deficit is exact");
    check(near(clipped.diagnostics.droplet_raw_minimum,
               clip_interstitial ? Real(-9.8) : Real(1)),
          "droplet raw minimum is exact");

    for (int j = 0; j < num_constituents; ++j) {
      const auto &diag = clipped.diagnostics.field_clipping[j];
      check(diag.interstitial_count == (clip_interstitial ? 1 : 0),
            "per-field interstitial clipping count is exact");
      check(diag.cloudborne_count == (clip_interstitial ? 0 : 1),
            "per-field cloudborne clipping count is exact");
      check(near(diag.interstitial_weighted_deficit,
                 clip_interstitial ? Real(9.8) : Real(0)),
            "per-field interstitial clipping deficit is exact");
      check(near(diag.cloudborne_weighted_deficit,
                 clip_interstitial ? Real(0) : Real(9.8)),
            "per-field cloudborne clipping deficit is exact");
      check(near(diag.interstitial_raw_minimum,
                 clip_interstitial ? Real(-9.8) : Real(1)),
            "per-field interstitial raw minimum is exact");
      check(near(diag.cloudborne_raw_minimum,
                 clip_interstitial ? Real(1) : Real(-9.8)),
            "per-field cloudborne raw minimum is exact");
      check(near(clipped.diagnostics.clipping_ledger_delta[j], 9.8),
            "per-field clipping ledger increment is exact");
      const Real total_change = clipped.diagnostics.final_ledger[j] -
                                clipped.diagnostics.initial_ledger[j];
      check(near(total_change,
                 clipped.diagnostics.clipping_ledger_delta[j]),
            "per-field ledger change equals clipping contribution");
    }
  }

  const State uniform = make_state(3, true);
  const auto unclipped =
      ndrop_test::advance_forward_euler(uniform, op, 1, clipped_options);
  check(unclipped.ok && unclipped.diagnostics.clipping_count == 0,
        "stable control case has no clipping");
  check(unclipped.diagnostics.droplet_clipping_count == 0 &&
            unclipped.diagnostics.droplet_weighted_deficit == 0,
        "stable control has zero droplet clipping diagnostics");
  for (int j = 0; j < num_constituents; ++j) {
    const auto &diag = unclipped.diagnostics.field_clipping[j];
    check(diag.interstitial_count == 0 && diag.cloudborne_count == 0 &&
              diag.interstitial_weighted_deficit == 0 &&
              diag.cloudborne_weighted_deficit == 0 &&
              unclipped.diagnostics.clipping_ledger_delta[j] == 0,
          "stable control has exact zero per-field clipping diagnostics");
  }
}

void test_zero_cloud_postprocess() {
  for (const std::size_t nlev : {1u, 3u, 8u}) {
    State state = make_state(nlev);
    OperatorData op(nlev);
    for (std::size_t k = 0; k < nlev; ++k) {
      op.cloud_fraction[k] =
          nlev == 1 || k == 1 || (nlev == 8 && k >= 4 && k <= 6) ? 0 : 1;
    }
    const auto initial_ledger =
        ndrop_test::column_ledgers(state, op.pressure_weight);
    const auto result = ndrop_test::advance_forward_euler(state, op, 1);
    check(result.ok, "zero-cloud FE case succeeds");
    for (std::size_t k = 0; k < nlev; ++k) {
      if (op.cloud_fraction[k] == 0) {
        check(result.state.droplet_number[k] == 0,
              "zero-cloud FE resets droplet number");
        for (int j = 0; j < num_constituents; ++j) {
          check(result.state.cloudborne[j][k] == 0,
                "zero-cloud FE resets cloudborne aerosol");
          check(near(result.state.interstitial[j][k],
                     state.interstitial[j][k] + state.cloudborne[j][k]),
                "zero-cloud FE resuspends cloudborne aerosol");
        }
      }
    }
    for (int j = 0; j < num_constituents; ++j) {
      check(near(result.diagnostics.final_ledger[j], initial_ledger[j]),
            "zero-cloud FE preserves aerosol ledger");
    }
  }
}

void test_finite_failures() {
  State state = make_state(2, true);
  OperatorData op(2);
  state.droplet_number[0] = std::numeric_limits<Real>::quiet_NaN();
  auto result = ndrop_test::advance_forward_euler(state, op, 1);
  check(!result.ok &&
            result.message == "droplet number must be finite and nonnegative",
        "FE nonfinite input fails loudly: " + result.message);

  state = make_state(2, true);
  state.interstitial[0][0] = std::numeric_limits<Real>::max();
  state.cloudborne[0][0] = std::numeric_limits<Real>::max();
  op.pressure_weight[0] = std::numeric_limits<Real>::min();
  op.pressure_weight[1] = std::numeric_limits<Real>::min();
  op.cloud_fraction[0] = 0;
  result = ndrop_test::advance_forward_euler(state, op, 1);
  check(!result.ok &&
            result.message == "zero-cloud aerosol resuspension overflowed",
        "zero-cloud FE resuspension overflow fails loudly: " +
            result.message);

  state = make_state(2, true);
  op = OperatorData(2);
  op.pressure_weight[0] = std::numeric_limits<Real>::max();
  state.interstitial[0][0] = std::numeric_limits<Real>::max();
  result = ndrop_test::advance_forward_euler(state, op, 1);
  check(!result.ok &&
            result.message ==
                "aerosol column ledger exceeds the representable range",
        "FE initial ledger overflow fails loudly: " + result.message);

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e-200);
  op.interface_exchange[1] = 1e-100;
  op.mass_activation_rate[0][0] = 1e100;
  state.interstitial[1][1] = 1e300;
  ForwardEulerOptions forced;
  forced.forced_substeps = 1;
  result = ndrop_test::advance_forward_euler(state, op, 1, forced);
  check(!result.ok &&
            result.message ==
                "aerosol activation source contains a nonfinite value",
        "FE activation-source overflow fails loudly: " + result.message);

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e-200);
  op.interface_exchange[1] = 1e-100;
  op.number_activation_rate[0][0] = 1e100;
  state.interstitial[0][1] = 1e300;
  result = ndrop_test::advance_forward_euler(state, op, 1, forced);
  check(!result.ok &&
            result.message ==
                "droplet activation source contains a nonfinite value",
        "FE droplet-source overflow fails loudly: " + result.message);

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e-300);
  op.interface_exchange[1] = 1e-200;
  state.interstitial[0] = {5e307, 0};
  result = ndrop_test::advance_forward_euler(state, op, 4, forced);
  check(!result.ok &&
            result.message ==
                "Forward Euler update produced a nonfinite value",
        "FE raw-update overflow fails loudly: " + result.message);

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e100);
  op.interface_exchange[1] = 1e200;
  state.interstitial[0] = {1e200, 0};
  result = ndrop_test::advance_forward_euler(state, op, 1, forced);
  check(!result.ok &&
            result.message ==
                "Forward Euler clipping diagnostic overflowed",
        "FE clipping-deficit overflow fails loudly: " + result.message);

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e100);
  op.interface_exchange[1] = 1e200;
  state.interstitial[0] = {1e108, 0};
  state.cloudborne[0] = {1e108, 0};
  result = ndrop_test::advance_forward_euler(state, op, 1, forced);
  check(!result.ok &&
            result.message ==
                "Forward Euler clipping ledger diagnostic overflowed",
        "FE clipping-ledger overflow fails loudly: " + result.message);
}

struct ErrorPair {
  Real weighted_l1 = 0;
  Real linf = 0;
};

ErrorPair normalized_error(const std::vector<Real> &actual,
                           const std::vector<Real> &exact,
                           const std::vector<Real> &initial,
                           const std::vector<Real> &weights) {
  ErrorPair error;
  Real exact_l1 = 0;
  Real initial_l1 = 0;
  Real exact_linf = 0;
  Real initial_linf = 0;
  for (std::size_t k = 0; k < actual.size(); ++k) {
    error.weighted_l1 += weights[k] * std::abs(actual[k] - exact[k]);
    exact_l1 += weights[k] * std::abs(exact[k]);
    initial_l1 += weights[k] * std::abs(initial[k]);
    error.linf = std::max(error.linf, std::abs(actual[k] - exact[k]));
    exact_linf = std::max(exact_linf, std::abs(exact[k]));
    initial_linf = std::max(initial_linf, std::abs(initial[k]));
  }
  error.weighted_l1 /= std::max(exact_l1, initial_l1);
  error.linf /= std::max(exact_linf, initial_linf);
  return error;
}

std::vector<ErrorPair> all_field_errors(
    const State &actual, const State &exact, const State &initial,
    const std::vector<Real> &weights) {
  std::vector<ErrorPair> errors;
  errors.push_back(normalized_error(actual.droplet_number,
                                    exact.droplet_number,
                                    initial.droplet_number, weights));
  for (int j = 0; j < num_constituents; ++j) {
    errors.push_back(normalized_error(actual.interstitial[j],
                                      exact.interstitial[j],
                                      initial.interstitial[j], weights));
    errors.push_back(normalized_error(actual.cloudborne[j],
                                      exact.cloudborne[j],
                                      initial.cloudborne[j], weights));
  }
  return errors;
}

std::pair<Real, Real> fitted_order(const std::vector<Real> &steps,
                                  const std::vector<Real> &errors) {
  const std::size_t first = steps.size() - 3;
  Real mean_x = 0;
  Real mean_y = 0;
  for (std::size_t i = first; i < steps.size(); ++i) {
    mean_x += std::log(steps[i]);
    mean_y += std::log(errors[i]);
  }
  mean_x /= 3;
  mean_y /= 3;
  Real covariance = 0;
  Real variance_x = 0;
  Real variance_y = 0;
  for (std::size_t i = first; i < steps.size(); ++i) {
    const Real dx = std::log(steps[i]) - mean_x;
    const Real dy = std::log(errors[i]) - mean_y;
    covariance += dx * dy;
    variance_x += dx * dx;
    variance_y += dy * dy;
  }
  return {covariance / variance_x,
          covariance * covariance / (variance_x * variance_y)};
}

State integrate_fe(State state, const OperatorData &op, const Real total_time,
                   const int steps) {
  ForwardEulerOptions options;
  options.forced_substeps = 1;
  options.apply_clipping = false;
  const Real dt = total_time / steps;
  for (int step = 0; step < steps; ++step) {
    const auto result =
        ndrop_test::advance_forward_euler(state, op, dt, options);
    if (!result.ok) {
      throw std::runtime_error("FE integration failed: " + result.message);
    }
    state = result.state;
  }
  return state;
}

void check_convergence(const std::array<std::vector<Real>, 51> &l1_errors,
                       const std::array<std::vector<Real>, 51> &linf_errors,
                       const std::vector<Real> &steps,
                       const std::string &label) {
  Real min_order = std::numeric_limits<Real>::infinity();
  Real max_order = -std::numeric_limits<Real>::infinity();
  Real min_r2 = 1;
  for (std::size_t field = 0; field < l1_errors.size(); ++field) {
    const auto [l1_order, l1_r2] =
        fitted_order(steps, l1_errors[field]);
    const auto [linf_order, linf_r2] =
        fitted_order(steps, linf_errors[field]);
    check(l1_order >= 0.85 && l1_order <= 1.15 && l1_r2 >= 0.98,
          label + " weighted-L1 first-order convergence");
    check(linf_order >= 0.85 && linf_order <= 1.15 && linf_r2 >= 0.98,
          label + " Linf first-order convergence");
    min_order = std::min({min_order, l1_order, linf_order});
    max_order = std::max({max_order, l1_order, linf_order});
    min_r2 = std::min({min_r2, l1_r2, linf_r2});
  }
  std::cout << "CONVERGENCE case=" << label << " order_range=" << min_order
            << ':' << max_order << " min_r2=" << min_r2 << '\n';
}

void test_diffusion_convergence() {
  constexpr std::size_t nlev = 32;
  constexpr Real exchange = 0.2;
  constexpr Real total_time = 2;
  constexpr Real pi = 3.141592653589793238462643383279502884;
  State initial(nlev);
  State exact(nlev);
  OperatorData op(nlev);
  for (std::size_t interface = 1; interface < nlev; ++interface) {
    op.interface_exchange[interface] = exchange;
  }
  const Real eigenvalue =
      -4 * exchange * std::pow(std::sin(pi / (2 * nlev)), 2);
  const Real decay = std::exp(eigenvalue * total_time);
  for (std::size_t k = 0; k < nlev; ++k) {
    const Real mode = std::cos(pi * (static_cast<Real>(k) + 0.5) / nlev);
    initial.droplet_number[k] = 10 + mode;
    exact.droplet_number[k] = 10 + decay * mode;
    for (int j = 0; j < num_constituents; ++j) {
      initial.interstitial[j][k] = 2 + 0.01 * j + 0.1 * mode;
      exact.interstitial[j][k] = 2 + 0.01 * j + 0.1 * decay * mode;
      initial.cloudborne[j][k] = 0.5 + 0.002 * j + 0.05 * mode;
      exact.cloudborne[j][k] = 0.5 + 0.002 * j + 0.05 * decay * mode;
    }
  }
  const std::array<int, 5> counts = {8, 16, 32, 64, 128};
  std::vector<Real> steps;
  std::array<std::vector<Real>, 51> l1_errors;
  std::array<std::vector<Real>, 51> linf_errors;
  for (const int count : counts) {
    const auto numerical = integrate_fe(initial, op, total_time, count);
    const auto errors =
        all_field_errors(numerical, exact, initial, op.pressure_weight);
    steps.push_back(total_time / count);
    for (std::size_t field = 0; field < errors.size(); ++field) {
      l1_errors[field].push_back(errors[field].weighted_l1);
      linf_errors[field].push_back(errors[field].linf);
    }
  }
  check_convergence(l1_errors, linf_errors, steps, "fe_diffusion");
}

Real exponential_integral(const Real rate, const Real time) {
  return std::abs(rate) < 1e-14 ? time : std::expm1(rate * time) / rate;
}

std::array<Real, 2>
forced_cloud_pair(const std::array<Real, 2> &initial_cloud, const Real d,
                  const Real a, const Real time, const Real coefficient_a,
                  const Real coefficient_b, const Real rate_a,
                  const Real rate_b) {
  const Real integral =
      coefficient_a * exponential_integral(rate_a, time) +
      coefficient_b * exponential_integral(rate_b, time);
  const Real sum = initial_cloud[0] + initial_cloud[1] + a * integral;
  const Real weighted_integral =
      coefficient_a * exponential_integral(rate_a + 2 * d, time) +
      coefficient_b * exponential_integral(rate_b + 2 * d, time);
  const Real difference =
      std::exp(-2 * d * time) *
      (initial_cloud[0] - initial_cloud[1] + a * weighted_integral);
  return {(sum + difference) / 2, (sum - difference) / 2};
}

struct PairSolution {
  std::array<Real, 2> interstitial;
  std::array<Real, 2> cloudborne;
  Real coefficient_a;
  Real coefficient_b;
  Real rate_a;
  Real rate_b;
};

PairSolution exact_activation_pair(
    const std::array<Real, 2> &initial_interstitial,
    const std::array<Real, 2> &initial_cloudborne, const Real d,
    const Real activation, const Real total_time) {
  const Real s = std::sqrt(d * (d - activation));
  const Real rate_a = -d + s;
  const Real rate_b = -d - s;
  const Real coefficient_a =
      Real(0.5) *
      (initial_interstitial[1] + d * initial_interstitial[0] / s);
  const Real coefficient_b =
      Real(0.5) *
      (initial_interstitial[1] - d * initial_interstitial[0] / s);
  PairSolution solution{
      {std::exp(-d * total_time) *
           (std::cosh(s * total_time) * initial_interstitial[0] +
            (d - activation) * std::sinh(s * total_time) *
                initial_interstitial[1] / s),
       coefficient_a * std::exp(rate_a * total_time) +
           coefficient_b * std::exp(rate_b * total_time)},
      {},
      coefficient_a,
      coefficient_b,
      rate_a,
      rate_b};
  solution.cloudborne =
      forced_cloud_pair(initial_cloudborne, d, activation, total_time,
                        coefficient_a, coefficient_b, rate_a, rate_b);
  return solution;
}

void test_all_field_activation_convergence() {
  constexpr Real d = 0.5;
  constexpr std::array<Real, num_modes> number_activation = {0.11, 0.17,
                                                              0.23, 0.29};
  constexpr std::array<Real, num_modes> mass_activation = {0.07, 0.13, 0.19,
                                                            0.25};
  constexpr Real total_time = 1;
  State initial = make_state(2);
  State exact = initial;
  OperatorData op(2);
  op.interface_exchange[1] = d;
  for (int mode = 0; mode < num_modes; ++mode) {
    op.number_activation_rate[mode][0] = number_activation[mode];
    op.mass_activation_rate[mode][0] = mass_activation[mode];
  }

  std::array<std::array<Real, 2>, num_modes> number_forcing{};
  for (int j = 0; j < num_constituents; ++j) {
    const int mode = expected_mode[j];
    const Real activation = expected_number[j] ? number_activation[mode]
                                                : mass_activation[mode];
    const auto solution = exact_activation_pair(
        {initial.interstitial[j][0], initial.interstitial[j][1]},
        {initial.cloudborne[j][0], initial.cloudborne[j][1]}, d, activation,
        total_time);
    exact.interstitial[j] = {solution.interstitial[0],
                             solution.interstitial[1]};
    exact.cloudborne[j] = {solution.cloudborne[0],
                           solution.cloudborne[1]};
    if (expected_number[j]) {
      number_forcing[mode] =
          forced_cloud_pair({0, 0}, d, number_activation[mode], total_time,
                            solution.coefficient_a, solution.coefficient_b,
                            solution.rate_a, solution.rate_b);
    }
  }
  exact.droplet_number =
      {Real(0.5) *
           (initial.droplet_number[0] + initial.droplet_number[1] +
            std::exp(-2 * d * total_time) *
                (initial.droplet_number[0] - initial.droplet_number[1])),
       Real(0.5) *
           (initial.droplet_number[0] + initial.droplet_number[1] -
            std::exp(-2 * d * total_time) *
                (initial.droplet_number[0] - initial.droplet_number[1]))};
  for (const auto &forcing : number_forcing) {
    exact.droplet_number[0] += forcing[0];
    exact.droplet_number[1] += forcing[1];
  }

  const std::array<int, 5> counts = {8, 16, 32, 64, 128};
  std::vector<Real> steps;
  std::array<std::vector<Real>, 51> l1_errors;
  std::array<std::vector<Real>, 51> linf_errors;
  for (const int count : counts) {
    const auto numerical = integrate_fe(initial, op, total_time, count);
    const auto errors =
        all_field_errors(numerical, exact, initial, op.pressure_weight);
    steps.push_back(total_time / count);
    for (std::size_t field = 0; field < errors.size(); ++field) {
      l1_errors[field].push_back(errors[field].weighted_l1);
      linf_errors[field].push_back(errors[field].linf);
    }
  }
  check_convergence(l1_errors, linf_errors, steps, "fe_activation_all_fields");
}

} // namespace

int main() {
  try {
    test_invalid_inputs();
    test_identity_and_level_counts();
    test_uniform_full_overlap();
    test_substeps_limiter_and_parity();
    test_fractional_overlap_ledger();
    test_adaptive_subcycling_avoids_clipping();
    test_forced_step_clipping_accounting();
    test_zero_cloud_postprocess();
    test_finite_failures();
    test_diffusion_convergence();
    test_all_field_activation_convergence();
  } catch (const std::exception &error) {
    ++failed;
    std::cerr << "UNCAUGHT: " << error.what() << '\n';
  }

  std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
