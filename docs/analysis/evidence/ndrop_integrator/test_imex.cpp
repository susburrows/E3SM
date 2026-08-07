#include "ndrop_integrator.hpp"

#include <Kokkos_Core.hpp>

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
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::State;
using ndrop_test::num_constituents;
using ndrop_test::num_modes;

// Independent oracle mapping; do not replace with constituent_mode or
// constituent_is_number from the implementation under test.
constexpr std::array<int, num_constituents> expected_mode = {
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3};
constexpr std::array<bool, num_constituents> expected_number = {
    true,  false, false, false, false, false, false, false, true,
    false, false, false, false, true,  false, false, false, false,
    false, false, false, true,  false, false, false};
constexpr std::array<int, num_modes> expected_number_index = {0, 8, 13, 21};

int passed = 0;
int failed = 0;

void check(const bool condition, const std::string &message) {
  if (condition) {
    ++passed;
  } else {
    ++failed;
    std::cerr << "FAIL: " << message << '\n';
  }
}

bool near(const Real lhs, const Real rhs, const Real factor = 500) {
  const Real tolerance =
      factor * std::numeric_limits<Real>::epsilon() *
      std::max({std::abs(lhs), std::abs(rhs), Real(1)});
  return std::abs(lhs - rhs) <= tolerance;
}

std::array<Real, 2>
solve_symmetric_diffusion(const std::array<Real, 2> &rhs, const Real rate,
                          const Real dt) {
  const Real diagonal = 1 + dt * rate;
  const Real off_diagonal = -dt * rate;
  const Real determinant =
      diagonal * diagonal - off_diagonal * off_diagonal;
  return {(diagonal * rhs[0] - off_diagonal * rhs[1]) / determinant,
          (-off_diagonal * rhs[0] + diagonal * rhs[1]) / determinant};
}

State make_state(const std::size_t nlev, const bool uniform = false) {
  State state(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    const Real level = static_cast<Real>(k);
    state.droplet_number[k] = uniform ? 4 : 4 + 0.02 * level;
    for (int j = 0; j < num_constituents; ++j) {
      state.interstitial[j][k] =
          uniform ? 1 + 0.01 * j : 1 + 0.01 * j + 0.003 * level;
      state.cloudborne[j][k] =
          uniform ? 0.2 + 0.002 * j : 0.2 + 0.002 * j + 0.001 * level;
    }
  }
  return state;
}

OperatorData make_operator(const std::size_t nlev, const Real exchange,
                           const Real activation, const bool fractional) {
  OperatorData op(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    op.pressure_weight[k] = 0.8 + 0.01 * static_cast<Real>(k);
    op.cloud_fraction[k] =
        fractional
            ? std::array<Real, 4>{0, 0.2, 0.65, 1}[k % 4]
            : Real(1);
  }
  for (std::size_t interface = 1; interface < nlev; ++interface) {
    op.interface_exchange[interface] =
        exchange *
        std::min(op.pressure_weight[interface - 1],
                 op.pressure_weight[interface]);
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    std::fill(op.number_activation_rate[mode].begin(),
              op.number_activation_rate[mode].end(), activation);
    std::fill(op.mass_activation_rate[mode].begin(),
              op.mass_activation_rate[mode].end(), Real(0.5) * activation);
    op.number_activation_rate[mode].back() = 0;
    op.mass_activation_rate[mode].back() = 0;
  }
  return op;
}

void check_ledgers(const AdvanceResult &result, const std::string &label) {
  const Real tolerance =
      std::max(Real(1e-11),
               5 * static_cast<Real>(result.state.nlev()) *
                   std::numeric_limits<Real>::epsilon());
  for (int j = 0; j < num_constituents; ++j) {
    const Real scale =
        std::max(std::abs(result.diagnostics.initial_ledger[j]), Real(1));
    check(std::abs(result.diagnostics.final_ledger[j] -
                   result.diagnostics.initial_ledger[j]) <=
              tolerance * scale,
          label + " ledger " + std::to_string(j));
  }
}

void test_identity_level_counts() {
  for (const std::size_t nlev : {1u, 2u, 3u, 127u, 128u}) {
    const State state = make_state(nlev);
    const OperatorData op(nlev);
    const auto result = ndrop_test::advance_imex(state, op, 100);
    check(result.ok, "identity IMEX succeeds");
    check(result.diagnostics.normalized_residual <= 1e-11,
          "identity residual");
    check(!result.diagnostics.positivity_failure, "identity positivity");
    for (std::size_t k = 0; k < nlev; ++k) {
      check(near(result.state.droplet_number[k], state.droplet_number[k]),
            "identity droplet");
      for (int j = 0; j < num_constituents; ++j) {
        check(near(result.state.interstitial[j][k],
                   state.interstitial[j][k]),
              "identity interstitial");
        check(near(result.state.cloudborne[j][k], state.cloudborne[j][k]),
              "identity cloudborne");
      }
    }
    check_ledgers(result, "identity");
  }
}

void test_uniform_and_fractional_ledgers() {
  State state = make_state(128, true);
  OperatorData op = make_operator(128, 0.7, 0, false);
  auto result = ndrop_test::advance_imex(state, op, 10);
  check(result.ok, "uniform full-overlap IMEX succeeds");
  for (std::size_t k = 0; k < state.nlev(); ++k) {
    check(near(result.state.droplet_number[k], state.droplet_number[k], 2000),
          "uniform droplet");
    for (int j = 0; j < num_constituents; ++j) {
      check(near(result.state.interstitial[j][k], state.interstitial[j][k],
                 2000),
            "uniform interstitial");
      check(near(result.state.cloudborne[j][k], state.cloudborne[j][k], 2000),
            "uniform cloudborne");
    }
  }
  check_ledgers(result, "uniform");

  state = make_state(127);
  op = make_operator(127, 0.8, 0.002, true);
  result = ndrop_test::advance_imex(state, op, 0.25);
  check(result.ok, "fractional-overlap IMEX succeeds");
  check(result.diagnostics.normalized_residual <= 1e-11,
        "fractional residual");
  check_ledgers(result, "fractional");
}

Real next_unit_random(unsigned long long &state) {
  state = state * 6364136223846793005ULL + 1442695040888963407ULL;
  return static_cast<Real>(state >> 11) /
         static_cast<Real>(1ULL << 53);
}

void test_randomized_compatible_columns() {
  unsigned long long seed = 0x4e44524f50554c4cULL;
  for (const std::size_t nlev : {1u, 2u, 3u, 127u, 128u}) {
    for (int sample = 0; sample < 4; ++sample) {
      State state(nlev);
      OperatorData op(nlev);
      for (std::size_t k = 0; k < nlev; ++k) {
        op.pressure_weight[k] = 0.2 + 2 * next_unit_random(seed);
        const Real cloud_selector = next_unit_random(seed);
        op.cloud_fraction[k] =
            cloud_selector < 0.15
                ? 0
                : (cloud_selector > 0.85 ? 1 : cloud_selector);
        state.droplet_number[k] = 1 + 10 * next_unit_random(seed);
        for (int j = 0; j < num_constituents; ++j) {
          state.interstitial[j][k] =
              0.1 + 2 * next_unit_random(seed);
          state.cloudborne[j][k] =
              0.05 + next_unit_random(seed);
        }
      }
      for (std::size_t interface = 1; interface < nlev; ++interface) {
        const Real rate = std::pow(10, -2 + 3 * next_unit_random(seed));
        op.interface_exchange[interface] =
            rate * std::min(op.pressure_weight[interface - 1],
                            op.pressure_weight[interface]);
      }
      for (int mode = 0; mode < num_modes; ++mode) {
        for (std::size_t k = 0; k + 1 < nlev; ++k) {
          const Real rate_below =
              op.interface_exchange[k + 1] / op.pressure_weight[k];
          op.number_activation_rate[mode][k] =
              0.001 * next_unit_random(seed) * rate_below;
          op.mass_activation_rate[mode][k] =
              0.0005 * next_unit_random(seed) * rate_below;
        }
      }
      const auto result = ndrop_test::advance_imex(state, op, 0.5);
      check(result.ok, "random compatible IMEX column succeeds");
      check(result.diagnostics.normalized_residual <= 1e-11,
            "random compatible IMEX residual");
      check(result.diagnostics.normalized_matrix_residual <= 1e-11,
            "random compatible IMEX matrix residual");
      check(result.diagnostics.normalized_physical_residual <= 1e-11,
            "random compatible IMEX independently assembled residual");
      check(!result.diagnostics.positivity_failure,
            "random compatible IMEX positivity");
      check_ledgers(result, "random compatible");
    }
  }
}

struct ErrorPair {
  Real weighted_l1 = 0;
  Real linf = 0;
};

ErrorPair normalized_field_error(const std::vector<Real> &actual,
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

std::vector<ErrorPair>
normalized_field_errors(const State &actual, const State &exact,
                        const State &initial,
                        const std::vector<Real> &weights) {
  std::vector<ErrorPair> errors;
  errors.reserve(1 + 2 * num_constituents);
  errors.push_back(normalized_field_error(
      actual.droplet_number, exact.droplet_number, initial.droplet_number,
      weights));
  for (int j = 0; j < num_constituents; ++j) {
    errors.push_back(normalized_field_error(
        actual.interstitial[j], exact.interstitial[j],
        initial.interstitial[j], weights));
    errors.push_back(normalized_field_error(
        actual.cloudborne[j], exact.cloudborne[j], initial.cloudborne[j],
        weights));
  }
  return errors;
}

std::pair<Real, Real> fit_order(const std::vector<Real> &steps,
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
  const Real slope = covariance / variance_x;
  const Real r2 = covariance * covariance / (variance_x * variance_y);
  return {slope, r2};
}

State integrate_imex(State state, const OperatorData &op, const Real total_time,
                     const int steps) {
  const Real dt = total_time / steps;
  for (int step = 0; step < steps; ++step) {
    auto result = ndrop_test::advance_imex(state, op, dt);
    if (!result.ok) {
      throw std::runtime_error("IMEX integration failed: " + result.message);
    }
    state = std::move(result.state);
  }
  return state;
}

State integrate_fe(State state, const OperatorData &op, const Real total_time,
                   const int steps) {
  ndrop_test::ForwardEulerOptions options;
  options.forced_substeps = 1;
  options.apply_clipping = false;
  const Real dt = total_time / steps;
  for (int step = 0; step < steps; ++step) {
    auto result =
        ndrop_test::advance_forward_euler(state, op, dt, options);
    if (!result.ok) {
      throw std::runtime_error("FE integration failed: " + result.message);
    }
    state = std::move(result.state);
  }
  return state;
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
      exact.cloudborne[j][k] =
          0.5 + 0.002 * j + 0.05 * decay * mode;
    }
  }

  const std::array<int, 5> counts = {8, 16, 32, 64, 128};
  std::vector<Real> steps;
  constexpr std::size_t field_count = 1 + 2 * num_constituents;
  std::array<std::vector<Real>, field_count> l1_errors;
  std::array<std::vector<Real>, field_count> linf_errors;
  for (const int count : counts) {
    const State numerical = integrate_imex(initial, op, total_time, count);
    const auto errors =
        normalized_field_errors(numerical, exact, initial, op.pressure_weight);
    steps.push_back(total_time / count);
    for (std::size_t field = 0; field < field_count; ++field) {
      l1_errors[field].push_back(errors[field].weighted_l1);
      linf_errors[field].push_back(errors[field].linf);
    }
  }
  Real min_l1_order = std::numeric_limits<Real>::infinity();
  Real max_l1_order = -std::numeric_limits<Real>::infinity();
  Real min_linf_order = std::numeric_limits<Real>::infinity();
  Real max_linf_order = -std::numeric_limits<Real>::infinity();
  Real min_r2 = 1;
  for (std::size_t field = 0; field < field_count; ++field) {
    const auto [l1_order, l1_r2] = fit_order(steps, l1_errors[field]);
    const auto [linf_order, linf_r2] = fit_order(steps, linf_errors[field]);
    check(l1_order >= 0.85 && l1_order <= 1.15 && l1_r2 >= 0.98,
          "IMEX diffusion field weighted-L1 first-order convergence");
    check(linf_order >= 0.85 && linf_order <= 1.15 && linf_r2 >= 0.98,
          "IMEX diffusion field Linf first-order convergence");
    min_l1_order = std::min(min_l1_order, l1_order);
    max_l1_order = std::max(max_l1_order, l1_order);
    min_linf_order = std::min(min_linf_order, linf_order);
    max_linf_order = std::max(max_linf_order, linf_order);
    min_r2 = std::min({min_r2, l1_r2, linf_r2});
  }
  std::cout << "CONVERGENCE case=diffusion l1_order_range="
            << min_l1_order << ':' << max_l1_order
            << " linf_order_range=" << min_linf_order << ':'
            << max_linf_order << " min_r2=" << min_r2 << '\n';
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

void test_activation_convergence() {
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
  constexpr std::size_t field_count = 1 + 2 * num_constituents;
  std::array<std::vector<Real>, field_count> l1_errors;
  std::array<std::vector<Real>, field_count> linf_errors;
  for (const int count : counts) {
    const State numerical = integrate_imex(initial, op, total_time, count);
    const auto errors =
        normalized_field_errors(numerical, exact, initial, op.pressure_weight);
    steps.push_back(total_time / count);
    for (std::size_t field = 0; field < field_count; ++field) {
      l1_errors[field].push_back(errors[field].weighted_l1);
      linf_errors[field].push_back(errors[field].linf);
    }
  }
  Real min_l1_order = std::numeric_limits<Real>::infinity();
  Real max_l1_order = -std::numeric_limits<Real>::infinity();
  Real min_linf_order = std::numeric_limits<Real>::infinity();
  Real max_linf_order = -std::numeric_limits<Real>::infinity();
  Real min_r2 = 1;
  for (std::size_t field = 0; field < field_count; ++field) {
    const auto [l1_order, l1_r2] = fit_order(steps, l1_errors[field]);
    const auto [linf_order, linf_r2] =
        fit_order(steps, linf_errors[field]);
    check(l1_order >= 0.85 && l1_order <= 1.15 && l1_r2 >= 0.98,
          "IMEX activation field weighted-L1 first-order convergence");
    check(linf_order >= 0.85 && linf_order <= 1.15 && linf_r2 >= 0.98,
          "IMEX activation field Linf first-order convergence");
    min_l1_order = std::min(min_l1_order, l1_order);
    max_l1_order = std::max(max_l1_order, l1_order);
    min_linf_order = std::min(min_linf_order, linf_order);
    max_linf_order = std::max(max_linf_order, linf_order);
    min_r2 = std::min({min_r2, l1_r2, linf_r2});
  }
  std::cout << "CONVERGENCE case=activation l1_order_range="
            << min_l1_order << ':' << max_l1_order
            << " linf_order_range=" << min_linf_order << ':'
            << max_linf_order << " min_r2=" << min_r2 << '\n';
}

void test_activation_limiter_oracle() {
  constexpr Real exchange = 0.25;
  constexpr Real dt = 0.01;
  State initial = make_state(2);
  OperatorData op(2);
  op.interface_exchange[1] = exchange;
  for (int mode = 0; mode < num_modes; ++mode) {
    op.number_activation_rate[mode][0] = 10 + mode;
    op.mass_activation_rate[mode][0] = 20 + mode;
  }

  const auto result = ndrop_test::advance_imex(initial, op, dt);
  check(result.ok, "limiter-active IMEX succeeds: " + result.message);
  check(result.diagnostics.activation_limiter_count == 2 * num_modes,
        "limiter-active IMEX reports all capped rates");

  Real droplet_source = 0;
  for (const int j : expected_number_index) {
    droplet_source += exchange * initial.interstitial[j][1];
  }
  const auto expected_droplet = solve_symmetric_diffusion(
      {initial.droplet_number[0] + dt * droplet_source,
       initial.droplet_number[1]},
      exchange, dt);
  check(near(result.state.droplet_number[0], expected_droplet[0]) &&
            near(result.state.droplet_number[1], expected_droplet[1]),
        "limiter-active droplet solve uses capped modal sources");

  for (int j = 0; j < num_constituents; ++j) {
    const Real source = exchange * initial.interstitial[j][1];
    const auto expected_cloud = solve_symmetric_diffusion(
        {initial.cloudborne[j][0] + dt * source,
         initial.cloudborne[j][1]},
        exchange, dt);
    const auto expected_interstitial = solve_symmetric_diffusion(
        {initial.interstitial[j][0] - dt * source,
         initial.interstitial[j][1]},
        exchange, dt);
    check(near(result.state.cloudborne[j][0], expected_cloud[0]) &&
              near(result.state.cloudborne[j][1], expected_cloud[1]),
          "limiter-active cloudborne field matches direct solve");
    check(near(result.state.interstitial[j][0], expected_interstitial[0]) &&
              near(result.state.interstitial[j][1], expected_interstitial[1]),
          "limiter-active interstitial field matches direct solve");
  }
  check(result.diagnostics.normalized_residual <= 1e-11,
        "limiter-active IMEX residual");
  check_ledgers(result, "limiter-active IMEX");
}

void test_positivity_failure() {
  State state = make_state(2, true);
  state.interstitial[0] = {0.001, 100};
  state.cloudborne[0] = {0, 0};
  OperatorData op(2);
  op.interface_exchange[1] = 0.5;
  op.number_activation_rate[0][0] = 0.5;
  const auto result = ndrop_test::advance_imex(state, op, 100);
  check(!result.ok, "material-negative IMEX case fails");
  check(result.diagnostics.positivity_failure,
        "material-negative IMEX case is classified");
  check(result.diagnostics.negative_count > 0,
        "material-negative IMEX count is reported");
  check(result.message == "material negative value in IMEX solution",
        "material-negative IMEX message");

  state = make_state(2, true);
  for (int j = 0; j < num_constituents; ++j) {
    state.interstitial[j][0] *= 1e-20;
    state.interstitial[j][1] *= 1e-20;
    state.cloudborne[j][0] *= 1e-20;
    state.cloudborne[j][1] *= 1e-20;
  }
  state.droplet_number[0] *= 1e-20;
  state.droplet_number[1] *= 1e-20;
  state.interstitial[0] = {1e-23, 1e-18};
  state.cloudborne[0] = {0, 0};
  const auto tiny_result = ndrop_test::advance_imex(state, op, 100);
  check(!tiny_result.ok && tiny_result.diagnostics.positivity_failure,
        "scale-relative positivity check rejects tiny material negative");
}

void test_zero_cloud_postprocess() {
  for (const std::size_t nlev : {1u, 3u, 8u}) {
    State state = make_state(nlev);
    OperatorData op(nlev);
    for (std::size_t k = 0; k < nlev; ++k) {
      op.cloud_fraction[k] =
          nlev == 1 || k == 1 || (nlev == 8 && k >= 4 && k <= 6) ? 0 : 1;
    }
    const auto result = ndrop_test::advance_imex(state, op, 1);
    check(result.ok, "zero-cloud IMEX case succeeds");
    for (std::size_t k = 0; k < nlev; ++k) {
      if (op.cloud_fraction[k] == 0) {
        check(result.state.droplet_number[k] == 0,
              "zero-cloud IMEX resets droplet number");
        for (int j = 0; j < num_constituents; ++j) {
          check(result.state.cloudborne[j][k] == 0,
                "zero-cloud IMEX resets cloudborne aerosol");
          check(near(result.state.interstitial[j][k],
                     state.interstitial[j][k] + state.cloudborne[j][k]),
                "zero-cloud IMEX resuspends cloudborne aerosol");
        }
      }
    }
    check_ledgers(result, "zero-cloud IMEX");
  }
}

void test_finite_failures() {
  State state = make_state(2, true);
  OperatorData op(2);
  state.droplet_number[0] = std::numeric_limits<Real>::quiet_NaN();
  auto result = ndrop_test::advance_imex(state, op, 1);
  check(!result.ok &&
            result.message == "droplet number must be finite and nonnegative",
        "IMEX nonfinite input fails loudly");

  state = make_state(2, true);
  op.pressure_weight[0] = std::numeric_limits<Real>::max();
  state.interstitial[0][0] = std::numeric_limits<Real>::max();
  result = ndrop_test::advance_imex(state, op, 1);
  check(!result.ok &&
            result.message ==
                "aerosol column ledger exceeds the representable range",
        "IMEX initial-ledger overflow fails loudly");

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e-200);
  op.interface_exchange[1] = 1e-100;
  op.mass_activation_rate[0][0] = 1e100;
  state.interstitial[1][1] = 1e300;
  result = ndrop_test::advance_imex(state, op, 1);
  check(!result.ok &&
            result.message ==
                "aerosol activation source contains a nonfinite value",
        "IMEX activation-source overflow fails loudly");

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e-200);
  op.interface_exchange[1] = 1e-100;
  op.number_activation_rate[0][0] = 1e100;
  state.interstitial[0][1] = 1e300;
  result = ndrop_test::advance_imex(state, op, 1);
  check(!result.ok &&
            result.message ==
                "droplet activation source contains a nonfinite value",
        "IMEX droplet-source overflow fails loudly");

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e-100);
  op.interface_exchange[1] = 1;
  op.mass_activation_rate[0][0] = 1e100;
  state.interstitial[1][1] = 1e100;
  result = ndrop_test::advance_imex(state, op, 1e200);
  check(!result.ok &&
            result.message ==
                "IMEX cloud right-hand side contains a nonfinite value",
        "IMEX right-hand-side overflow fails loudly");

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(),
            std::numeric_limits<Real>::min());
  op.cloud_fraction[0] = 0;
  state.interstitial[0][0] = std::numeric_limits<Real>::max();
  state.cloudborne[0][0] = std::numeric_limits<Real>::max();
  result = ndrop_test::advance_imex(state, op, 1);
  check(!result.ok &&
            result.message == "zero-cloud aerosol resuspension overflowed",
        "IMEX zero-cloud resuspension overflow fails loudly");
}

void test_nonfinite_checker_probes() {
  for (const auto probe :
       {ndrop_test::ImexFiniteProbe::solution,
        ndrop_test::ImexFiniteProbe::matrix_residual,
        ndrop_test::ImexFiniteProbe::physical_residual,
        ndrop_test::ImexFiniteProbe::ledger_normalization}) {
    check(ndrop_test::imex_nonfinite_probe_is_rejected(probe),
          "IMEX finite-stage checker rejects injected nonfinite value");
  }
}

void test_coupled_cross_scheme_convergence() {
  constexpr std::size_t nlev = 12;
  constexpr Real total_time = 1;
  State initial = make_state(nlev);
  OperatorData op(nlev);
  constexpr std::array<Real, 5> clouds = {0.2, 0.35, 0.55, 0.8, 1};
  for (std::size_t k = 0; k < nlev; ++k) {
    op.pressure_weight[k] = 0.7 + 0.03 * k;
    op.cloud_fraction[k] = clouds[k % clouds.size()];
  }
  for (std::size_t interface = 1; interface < nlev; ++interface) {
    const Real rate = 0.08 + 0.015 * (interface % 4);
    op.interface_exchange[interface] =
        rate * std::min(op.pressure_weight[interface - 1],
                        op.pressure_weight[interface]);
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    for (std::size_t k = 0; k + 1 < nlev; ++k) {
      const Real rate_below =
          op.interface_exchange[k + 1] / op.pressure_weight[k];
      op.number_activation_rate[mode][k] =
          (0.04 + 0.005 * mode) * rate_below;
      op.mass_activation_rate[mode][k] =
          (0.02 + 0.003 * mode) * rate_below;
    }
  }

  constexpr std::array<int, 5> counts = {16, 32, 64, 128, 256};
  std::array<State, counts.size()> fe;
  std::array<State, counts.size()> imex;
  for (std::size_t i = 0; i < counts.size(); ++i) {
    fe[i] = integrate_fe(initial, op, total_time, counts[i]);
    imex[i] = integrate_imex(initial, op, total_time, counts[i]);
  }

  const auto fe_coarse =
      normalized_field_errors(fe[2], fe[3], initial, op.pressure_weight);
  const auto fe_fine =
      normalized_field_errors(fe[3], fe[4], initial, op.pressure_weight);
  const auto imex_coarse =
      normalized_field_errors(imex[2], imex[3], initial, op.pressure_weight);
  const auto imex_fine =
      normalized_field_errors(imex[3], imex[4], initial, op.pressure_weight);
  const auto cross_coarse =
      normalized_field_errors(fe[3], imex[3], initial, op.pressure_weight);
  const auto cross_fine =
      normalized_field_errors(fe[4], imex[4], initial, op.pressure_weight);

  Real min_order = std::numeric_limits<Real>::infinity();
  Real max_order = -std::numeric_limits<Real>::infinity();
  Real max_cross = 0;
  for (std::size_t field = 0; field < fe_fine.size(); ++field) {
    for (const auto [coarse, fine] :
         {std::pair<Real, Real>{fe_coarse[field].weighted_l1,
                               fe_fine[field].weighted_l1},
          std::pair<Real, Real>{fe_coarse[field].linf,
                               fe_fine[field].linf},
          std::pair<Real, Real>{imex_coarse[field].weighted_l1,
                               imex_fine[field].weighted_l1},
          std::pair<Real, Real>{imex_coarse[field].linf,
                               imex_fine[field].linf}}) {
      const Real order = std::log2(coarse / fine);
      check(order >= 0.85 && order <= 1.15,
            "coupled all-field first-order self-convergence");
      min_order = std::min(min_order, order);
      max_order = std::max(max_order, order);
    }
    check(cross_fine[field].weighted_l1 <
              cross_coarse[field].weighted_l1 &&
              cross_fine[field].linf < cross_coarse[field].linf,
          "coupled FE/IMEX difference decreases under refinement");
    check(cross_fine[field].weighted_l1 <=
                  2 * std::max(fe_fine[field].weighted_l1,
                               imex_fine[field].weighted_l1) &&
              cross_fine[field].linf <=
                  2 * std::max(fe_fine[field].linf,
                               imex_fine[field].linf),
          "coupled FE/IMEX finest difference satisfies Richardson bound");
    max_cross =
        std::max({max_cross, cross_fine[field].weighted_l1,
                  cross_fine[field].linf});
  }
  std::cout << "CONVERGENCE case=coupled_all_fields order_range="
            << min_order << ':' << max_order
            << " max_finest_cross_error=" << max_cross << '\n';
}

void test_nonfinite_assembly_failure() {
  const State state = make_state(2, true);
  OperatorData op(2);
  op.interface_exchange[1] = 1e200;
  const auto result = ndrop_test::advance_imex(state, op, 1e200);
  check(!result.ok, "nonfinite matrix assembly fails");
  check(result.message == "IMEX matrix contains a nonfinite coefficient",
        "nonfinite matrix failure message");
}

} // namespace

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  try {
    test_identity_level_counts();
    test_uniform_and_fractional_ledgers();
    test_randomized_compatible_columns();
    test_diffusion_convergence();
    test_activation_convergence();
    test_activation_limiter_oracle();
    test_positivity_failure();
    test_zero_cloud_postprocess();
    test_finite_failures();
    test_nonfinite_checker_probes();
    test_coupled_cross_scheme_convergence();
    test_nonfinite_assembly_failure();
  } catch (const std::exception &error) {
    ++failed;
    std::cerr << "UNCAUGHT: " << error.what() << '\n';
  }
  Kokkos::finalize();
  std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
