#include "ndrop_integrator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using ndrop_test::AdvanceResult;
using ndrop_test::OperatorData;
using ndrop_test::PreparedOperator;
using ndrop_test::Real;
using ndrop_test::State;
using ndrop_test::num_constituents;
using ndrop_test::num_modes;

// Independent oracle mappings. The dense oracle must not call the mapping or
// activation-source helpers used by the implementation under test.
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

bool near(const Real lhs, const Real rhs, const Real factor = 5000) {
  const Real tolerance =
      factor * std::numeric_limits<Real>::epsilon() *
      std::max({std::abs(lhs), std::abs(rhs), Real(1)});
  return std::abs(lhs - rhs) <= tolerance;
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
                           const Real activation, const Real cloud_top = 1) {
  OperatorData op(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    op.pressure_weight[k] = 0.8 + 0.01 * static_cast<Real>(k);
    op.cloud_fraction[k] = k == 0 ? cloud_top : Real(1);
  }
  for (std::size_t interface = 1; interface < nlev; ++interface) {
    op.interface_exchange[interface] =
        exchange * std::min(op.pressure_weight[interface - 1],
                            op.pressure_weight[interface]);
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    for (std::size_t k = 0; k + 1 < nlev; ++k) {
      const Real rate_below =
          op.interface_exchange[k + 1] / op.pressure_weight[k];
      op.number_activation_rate[mode][k] = activation * rate_below;
      op.mass_activation_rate[mode][k] =
          Real(0.5) * activation * rate_below;
    }
  }
  return op;
}

bool state_bitwise_equal(const State &lhs, const State &rhs) {
  if (lhs.nlev() != rhs.nlev() ||
      std::memcmp(lhs.droplet_number.data(), rhs.droplet_number.data(),
                  lhs.nlev() * sizeof(Real)) != 0) {
    return false;
  }
  for (int j = 0; j < num_constituents; ++j) {
    if (std::memcmp(lhs.interstitial[j].data(), rhs.interstitial[j].data(),
                    lhs.nlev() * sizeof(Real)) != 0 ||
        std::memcmp(lhs.cloudborne[j].data(), rhs.cloudborne[j].data(),
                    lhs.nlev() * sizeof(Real)) != 0) {
      return false;
    }
  }
  return true;
}

void check_state_near(const State &actual, const State &expected,
                      const std::string &label) {
  check(actual.nlev() == expected.nlev(), label + " level count");
  for (std::size_t k = 0; k < actual.nlev(); ++k) {
    check(near(actual.droplet_number[k], expected.droplet_number[k]),
          label + " droplet " + std::to_string(k));
    for (int j = 0; j < num_constituents; ++j) {
      check(near(actual.interstitial[j][k], expected.interstitial[j][k]),
            label + " interstitial " + std::to_string(j) + ":" +
                std::to_string(k));
      check(near(actual.cloudborne[j][k], expected.cloudborne[j][k]),
            label + " cloudborne " + std::to_string(j) + ":" +
                std::to_string(k));
    }
  }
}

void check_ledgers(const AdvanceResult &result, const std::string &label) {
  const Real tolerance =
      std::max(Real(1e-11),
               5 * static_cast<Real>(result.state.nlev()) *
                   std::numeric_limits<Real>::epsilon());
  for (int j = 0; j < num_constituents; ++j) {
    const Real direct = result.diagnostics.final_ledger[j] -
                        result.diagnostics.initial_ledger[j];
    const Real scale =
        std::max(std::abs(result.diagnostics.initial_ledger[j]), Real(1));
    check(std::abs(direct) <= tolerance * scale,
          label + " ledger " + std::to_string(j));
    check(near(result.diagnostics.total_ledger_delta[j], direct),
          label + " total delta " + std::to_string(j));
    check(result.diagnostics.clipping_ledger_delta[j] == 0,
          label + " no clipping " + std::to_string(j));
    check(near(result.diagnostics.operator_ledger_delta[j], direct),
          label + " operator delta " + std::to_string(j));
  }
}

std::vector<Real> dense_partial_pivot_solve(std::vector<Real> matrix,
                                            std::vector<Real> rhs) {
  const std::size_t n = rhs.size();
  for (std::size_t column = 0; column < n; ++column) {
    std::size_t pivot = column;
    for (std::size_t row = column + 1; row < n; ++row) {
      if (std::abs(matrix[row * n + column]) >
          std::abs(matrix[pivot * n + column])) {
        pivot = row;
      }
    }
    if (matrix[pivot * n + column] == 0) {
      throw std::runtime_error("singular independent dense oracle");
    }
    if (pivot != column) {
      for (std::size_t entry = column; entry < n; ++entry) {
        std::swap(matrix[column * n + entry], matrix[pivot * n + entry]);
      }
      std::swap(rhs[column], rhs[pivot]);
    }
    for (std::size_t row = column + 1; row < n; ++row) {
      const Real multiplier =
          matrix[row * n + column] / matrix[column * n + column];
      matrix[row * n + column] = 0;
      for (std::size_t entry = column + 1; entry < n; ++entry) {
        matrix[row * n + entry] -=
            multiplier * matrix[column * n + entry];
      }
      rhs[row] -= multiplier * rhs[column];
    }
  }
  std::vector<Real> solution(n, 0);
  for (std::size_t reverse = 0; reverse < n; ++reverse) {
    const std::size_t row = n - 1 - reverse;
    Real value = rhs[row];
    for (std::size_t column = row + 1; column < n; ++column) {
      value -= matrix[row * n + column] * solution[column];
    }
    solution[row] = value / matrix[row * n + row];
  }
  return solution;
}

// Dense oracle assembled directly from the signed physical equations. It does
// not call any block-band, matrix-family, source, or solver helper from BE.
State dense_backward_euler_oracle(const State &initial,
                                  const OperatorData &data, const Real dt) {
  const PreparedOperator op = ndrop_test::prepare_operator(data);
  const std::size_t nlev = initial.nlev();
  State result(nlev);
  for (int j = 0; j < num_constituents; ++j) {
    const std::size_t n = 2 * nlev;
    std::vector<Real> matrix(n * n, 0);
    std::vector<Real> rhs(n, 0);
    const int mode = expected_mode[j];
    const auto &rates = expected_number[j] ? op.number_activation_rate[mode]
                                           : op.mass_activation_rate[mode];
    for (std::size_t k = 0; k < nlev; ++k) {
      const std::size_t cloud_row = 2 * k;
      const std::size_t interstitial_row = cloud_row + 1;
      const Real diagonal =
          1 + dt * (op.rate_above[k] + op.rate_below[k]);
      matrix[cloud_row * n + cloud_row] = diagonal;
      matrix[interstitial_row * n + interstitial_row] = diagonal;
      rhs[cloud_row] = initial.cloudborne[j][k];
      rhs[interstitial_row] = initial.interstitial[j][k];
      if (k > 0) {
        matrix[cloud_row * n + 2 * (k - 1)] =
            -dt * op.rate_above[k] * op.overlap_above[k];
        matrix[interstitial_row * n + 2 * (k - 1)] =
            -dt * op.rate_above[k] * (1 - op.overlap_above[k]);
        matrix[interstitial_row * n + 2 * (k - 1) + 1] =
            -dt * op.rate_above[k];
      }
      if (k + 1 < nlev) {
        matrix[cloud_row * n + 2 * (k + 1)] =
            -dt * op.rate_below[k] * op.overlap_below[k];
        matrix[cloud_row * n + 2 * (k + 1) + 1] = -dt * rates[k];
        matrix[interstitial_row * n + 2 * (k + 1)] =
            -dt * op.rate_below[k] * (1 - op.overlap_below[k]);
        matrix[interstitial_row * n + 2 * (k + 1) + 1] =
            -dt * op.rate_below[k] + dt * rates[k];
      } else {
        const Real bottom_rate = op.number_activation_rate[mode][k];
        matrix[cloud_row * n + cloud_row] -= dt * bottom_rate;
        matrix[cloud_row * n + interstitial_row] -= dt * bottom_rate;
        matrix[interstitial_row * n + cloud_row] += dt * bottom_rate;
        matrix[interstitial_row * n + interstitial_row] += dt * bottom_rate;
      }
    }
    const auto solution = dense_partial_pivot_solve(matrix, rhs);
    for (std::size_t k = 0; k < nlev; ++k) {
      result.cloudborne[j][k] = solution[2 * k];
      result.interstitial[j][k] = solution[2 * k + 1];
    }
  }

  std::vector<Real> matrix(nlev * nlev, 0);
  std::vector<Real> rhs(nlev, 0);
  for (std::size_t k = 0; k < nlev; ++k) {
    matrix[k * nlev + k] =
        1 + dt * (op.rate_above[k] + op.rate_below[k]);
    if (k > 0) {
      matrix[k * nlev + k - 1] =
          -dt * op.rate_above[k] * op.overlap_above[k];
    }
    if (k + 1 < nlev) {
      matrix[k * nlev + k + 1] =
          -dt * op.rate_below[k] * op.overlap_below[k];
    }
    Real source = 0;
    for (int mode = 0; mode < num_modes; ++mode) {
      const int j = expected_number_index[mode];
      const std::size_t kp1 = std::min(k + 1, nlev - 1);
      source += op.number_activation_rate[mode][k] *
                result.interstitial[j][kp1];
      if (k + 1 == nlev) {
        source += op.number_activation_rate[mode][k] *
                  (result.interstitial[j][k] + result.cloudborne[j][k]);
      }
    }
    rhs[k] = initial.droplet_number[k] + dt * source;
  }
  result.droplet_number = dense_partial_pivot_solve(matrix, rhs);
  return result;
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
  return {covariance / variance_x,
          covariance * covariance / (variance_x * variance_y)};
}

State integrate_backward_euler(State state, const OperatorData &op,
                               const Real total_time, const int steps) {
  const Real dt = total_time / steps;
  for (int step = 0; step < steps; ++step) {
    auto result = ndrop_test::advance_backward_euler(state, op, dt);
    if (!result.ok) {
      throw std::runtime_error("BE integration failed: " + result.message);
    }
    state = std::move(result.state);
  }
  return state;
}

void check_first_order(const State &initial, const State &exact,
                       const OperatorData &op, const Real total_time,
                       const std::string &label) {
  constexpr std::array<int, 5> counts = {8, 16, 32, 64, 128};
  std::vector<Real> steps;
  constexpr std::size_t field_count = 1 + 2 * num_constituents;
  std::array<std::vector<Real>, field_count> l1_errors;
  std::array<std::vector<Real>, field_count> linf_errors;
  for (const int count : counts) {
    const State numerical =
        integrate_backward_euler(initial, op, total_time, count);
    const auto errors =
        normalized_field_errors(numerical, exact, initial, op.pressure_weight);
    steps.push_back(total_time / count);
    for (std::size_t field = 0; field < field_count; ++field) {
      l1_errors[field].push_back(errors[field].weighted_l1);
      linf_errors[field].push_back(errors[field].linf);
    }
  }
  Real min_order = std::numeric_limits<Real>::infinity();
  Real max_order = -std::numeric_limits<Real>::infinity();
  Real min_r2 = 1;
  for (std::size_t field = 0; field < field_count; ++field) {
    const auto [l1_order, l1_r2] = fit_order(steps, l1_errors[field]);
    const auto [linf_order, linf_r2] = fit_order(steps, linf_errors[field]);
    check(l1_order >= 0.85 && l1_order <= 1.15 && l1_r2 >= 0.98,
          label + " weighted-L1 first order");
    check(linf_order >= 0.85 && linf_order <= 1.15 && linf_r2 >= 0.98,
          label + " Linf first order");
    min_order = std::min({min_order, l1_order, linf_order});
    max_order = std::max({max_order, l1_order, linf_order});
    min_r2 = std::min({min_r2, l1_r2, linf_r2});
  }
  std::cout << "CONVERGENCE case=" << label << " order_range=" << min_order
            << ':' << max_order << " min_r2=" << min_r2 << '\n';
}

void test_identity_and_boundaries() {
  for (const std::size_t nlev : {1u, 2u, 3u, 127u, 128u}) {
    const State state = make_state(nlev);
    OperatorData op(nlev);
    for (int mode = 0; mode < num_modes; ++mode) {
      op.number_activation_rate[mode].back() = 100;
      op.mass_activation_rate[mode].back() = 200;
    }
    const auto result = ndrop_test::advance_backward_euler(state, op, 100);
    check(result.ok, "identity BE succeeds: " + result.message);
    check(result.diagnostics.nsubmix == 1,
          "BE reports one implicit advance");
    check(result.diagnostics.activation_limiter_count == 2 * num_modes,
          "no-flux bottom activation is capped to zero");
    check(result.diagnostics.normalized_residual <= 1e-11,
          "identity BE residual");
    check_state_near(result.state, state, "identity BE");
    check_ledgers(result, "identity BE");
  }
}

void test_dense_oracle_grid() {
  constexpr std::array<std::size_t, 5> levels = {1, 2, 3, 8, 3};
  constexpr std::array<Real, 5> overlaps = {1, 0.75, 0.5, 0.25, 1};
  constexpr std::array<Real, 5> exchanges = {0.05, 5, 0.05, 5, 0.5};
  constexpr std::array<Real, 5> strengths = {0, 0.1, 0.5, 0.9, 2};
  constexpr std::array<int, 5> active_counts = {1, 4, 8, 16, 25};
  for (std::size_t case_index = 0; case_index < levels.size(); ++case_index) {
    State state = make_state(levels[case_index]);
    for (int j = active_counts[case_index]; j < num_constituents; ++j) {
      std::fill(state.interstitial[j].begin(), state.interstitial[j].end(), 0);
      std::fill(state.cloudborne[j].begin(), state.cloudborne[j].end(), 0);
    }
    OperatorData op = make_operator(levels[case_index], exchanges[case_index],
                                    strengths[case_index],
                                    overlaps[case_index]);
    for (int mode = 0; mode < num_modes; ++mode) {
      for (std::size_t k = 0; k + 1 < levels[case_index]; ++k) {
        op.number_activation_rate[mode][k] *= 1 + Real(0.03) * mode;
        op.mass_activation_rate[mode][k] *= 1 + Real(0.05) * mode;
      }
    }
    const auto result = ndrop_test::advance_backward_euler(state, op, 0.2);
    check(result.ok, "dense-grid BE succeeds: " + result.message);
    const State oracle = dense_backward_euler_oracle(state, op, 0.2);
    check_state_near(result.state, oracle,
                     "independent dense case " + std::to_string(case_index));
    check(result.diagnostics.normalized_matrix_residual <= 1e-11,
          "dense-grid matrix residual");
    check(result.diagnostics.normalized_physical_residual <= 1e-11,
          "dense-grid physical residual");
    if (strengths[case_index] > 1) {
      check(result.diagnostics.activation_limiter_count > 0,
            "limiter-active dense-grid case reports capped rates");
    }
    check_ledgers(result, "dense-grid");
  }
}

void test_matrix_family_reuse() {
  State state = make_state(8);
  OperatorData op = make_operator(8, 1, 0);
  for (int mode = 0; mode < num_modes; ++mode) {
    for (std::size_t k = 0; k + 1 < state.nlev(); ++k) {
      op.number_activation_rate[mode][k] = 0.02 * (mode + 1);
      op.mass_activation_rate[mode][k] = 0.011 * (mode + 1);
    }
  }
  auto result = ndrop_test::advance_backward_euler(state, op, 0.25);
  check(result.ok, "distinct-family BE succeeds");
  check(result.diagnostics.block_matrix_family_count == 8 &&
            result.diagnostics.block_factorization_count == 8 &&
            result.diagnostics.block_rhs_solve_count == num_constituents,
        "BE factors eight distinct mode/type families for 25 right sides");
  check(result.diagnostics.scalar_factorization_count == 1 &&
            result.diagnostics.scalar_rhs_solve_count == 1,
        "BE factors and solves one droplet scalar family");

  op = make_operator(8, 1, 0);
  for (int mode = 0; mode < num_modes; ++mode) {
    for (std::size_t k = 0; k + 1 < state.nlev(); ++k) {
      op.number_activation_rate[mode][k] = 0.02;
      op.mass_activation_rate[mode][k] = 0.02;
    }
  }
  result = ndrop_test::advance_backward_euler(state, op, 0.25);
  check(result.ok, "shared-family BE succeeds");
  check(result.diagnostics.block_matrix_family_count == 1 &&
            result.diagnostics.block_factorization_count == 1 &&
            result.diagnostics.block_rhs_solve_count == num_constituents,
        "identical block matrix is factored once for all aerosol fields");
}

void test_controller_is_not_physical() {
  const State state = make_state(8);
  OperatorData base = make_operator(8, 0.4, 0.2, 0.5);
  OperatorData controlled = base;
  controlled.first_level_controller_rate = 1e6;
  const auto plain = ndrop_test::advance_backward_euler(state, base, 0.25);
  const auto with_controller =
      ndrop_test::advance_backward_euler(state, controlled, 0.25);
  check(plain.ok && with_controller.ok,
        "BE controller comparison succeeds");
  check(state_bitwise_equal(plain.state, with_controller.state),
        "first-level FE controller is absent from BE state equations");
  check(with_controller.diagnostics.first_level_controller_rate == 1e6,
        "BE reports the controller diagnostically");
}

void test_uniform_and_zero_cloud() {
  const State uniform = make_state(128, true);
  const OperatorData uniform_op = make_operator(128, 10, 0, 1);
  const auto uniform_result =
      ndrop_test::advance_backward_euler(uniform, uniform_op, 20);
  check(uniform_result.ok, "strong-mixing uniform BE succeeds");
  check_state_near(uniform_result.state, uniform, "uniform BE");
  check_ledgers(uniform_result, "uniform BE");

  for (const std::size_t nlev : {1u, 3u, 8u}) {
    State state = make_state(nlev);
    OperatorData op(nlev);
    for (std::size_t k = 0; k < nlev; ++k) {
      op.cloud_fraction[k] =
          nlev == 1 || k == 1 || (nlev == 8 && k >= 4 && k <= 6) ? 0 : 1;
    }
    const auto result = ndrop_test::advance_backward_euler(state, op, 1);
    check(result.ok, "zero-cloud BE succeeds");
    for (std::size_t k = 0; k < nlev; ++k) {
      if (op.cloud_fraction[k] == 0) {
        check(result.state.droplet_number[k] == 0,
              "zero-cloud BE resets droplets");
        for (int j = 0; j < num_constituents; ++j) {
          check(result.state.cloudborne[j][k] == 0,
                "zero-cloud BE resets cloudborne aerosol");
          check(near(result.state.interstitial[j][k],
                     state.interstitial[j][k] + state.cloudborne[j][k]),
                "zero-cloud BE resuspends cloudborne aerosol");
        }
      }
    }
    check_ledgers(result, "zero-cloud BE");
  }
}

void test_zero_overlap_with_mixing() {
  State state = make_state(3);
  OperatorData op(3);
  op.cloud_fraction = {1, 0, 1};
  op.interface_exchange[1] = 0.2;
  op.interface_exchange[2] = 0.2;
  const PreparedOperator prepared = ndrop_test::prepare_operator(op);
  check(prepared.overlap_above[1] == 0 &&
            prepared.overlap_below[1] == 0,
        "zero-cloud middle level creates exact zero overlap rows");
  const auto result = ndrop_test::advance_backward_euler(state, op, 0.25);
  check(result.ok, "zero-overlap mixing BE succeeds: " + result.message);
  check(result.state.droplet_number[1] == 0,
        "zero-overlap mixing projection resets droplets");
  for (int j = 0; j < num_constituents; ++j) {
    check(result.state.cloudborne[j][1] == 0,
          "zero-overlap mixing projection resets cloudborne aerosol");
  }
  check(result.diagnostics.normalized_physical_residual <= 1e-11,
        "zero-overlap mixing physical residual");
  check_ledgers(result, "zero-overlap mixing BE");
}

void test_serial_determinism() {
  const State state = make_state(128);
  const OperatorData op = make_operator(128, 0.8, 0.5, 0.25);
  const std::array<AdvanceResult, 3> results = {
      ndrop_test::advance_backward_euler(state, op, 0.25),
      ndrop_test::advance_backward_euler(state, op, 0.25),
      ndrop_test::advance_backward_euler(state, op, 0.25)};
  check(results[0].ok && results[1].ok && results[2].ok,
        "three serial BE advances succeed");
  check(state_bitwise_equal(results[0].state, results[1].state) &&
            state_bitwise_equal(results[1].state, results[2].state),
        "three serial BE output states are bitwise identical");
  check(results[0].diagnostics.normalized_residual ==
                results[1].diagnostics.normalized_residual &&
            results[1].diagnostics.normalized_residual ==
                results[2].diagnostics.normalized_residual,
        "three serial BE residual diagnostics are identical");
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
  check_first_order(initial, exact, op, total_time, "be_diffusion");
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
      {}, coefficient_a, coefficient_b, rate_a, rate_b};
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
    exact.cloudborne[j] = {solution.cloudborne[0], solution.cloudborne[1]};
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
  check_first_order(initial, exact, op, total_time, "be_activation");
}

void test_scaling_zero_ledgers_and_positivity() {
  for (const Real factor : {Real(1e-12), Real(1), Real(1e12)}) {
    State state = make_state(3);
    for (Real &value : state.droplet_number) {
      value *= factor;
    }
    for (int j = 0; j < num_constituents; ++j) {
      for (Real &value : state.interstitial[j]) {
        value *= factor;
      }
      for (Real &value : state.cloudborne[j]) {
        value *= factor;
      }
    }
    const OperatorData op = make_operator(3, 2, 0.9, 0.5);
    const auto result = ndrop_test::advance_backward_euler(state, op, 0.25);
    check(result.ok, "scaled BE succeeds: " + result.message);
    check(!result.diagnostics.positivity_failure &&
              result.diagnostics.material_negative_count == 0,
          "scaled BE remains materially nonnegative");
    check(std::isfinite(result.diagnostics.droplet_raw_minimum) &&
              std::isfinite(result.diagnostics.raw_minimum),
          "scaled BE reports finite raw solution minima");
    check(result.diagnostics.normalized_residual <= 1e-11,
          "scaled BE residual");
    check_ledgers(result, "scaled BE");
  }

  State zero(3);
  zero.droplet_number = {1, 2, 3};
  const auto result = ndrop_test::advance_backward_euler(
      zero, make_operator(3, 0.2, 0), 0.25);
  check(result.ok, "zero-ledger BE succeeds");
  for (int j = 0; j < num_constituents; ++j) {
    check(result.diagnostics.initial_ledger[j] == 0 &&
              result.diagnostics.total_ledger_delta[j] == 0 &&
              result.diagnostics.operator_ledger_delta[j] == 0,
          "zero-ledger BE uses absolute diagnostics");
  }
}

void test_fail_closed_paths() {
  State state = make_state(2, true);
  OperatorData op(2);
  state.droplet_number[0] = std::numeric_limits<Real>::quiet_NaN();
  auto result = ndrop_test::advance_backward_euler(state, op, 1);
  check(!result.ok &&
            result.message == "droplet number must be finite and nonnegative",
        "BE rejects nonfinite input");

  state = make_state(2, true);
  op = OperatorData(2);
  op.pressure_weight[0] = std::numeric_limits<Real>::max();
  state.interstitial[0][0] = std::numeric_limits<Real>::max();
  result = ndrop_test::advance_backward_euler(state, op, 1);
  check(!result.ok &&
            result.message ==
                "aerosol column ledger exceeds the representable range",
        "BE rejects initial-ledger overflow");

  state = make_state(2, true);
  op = OperatorData(2);
  op.interface_exchange[1] = 1e200;
  result = ndrop_test::advance_backward_euler(state, op, 1e200);
  check(!result.ok &&
            result.message ==
                "Backward Euler block matrix contains a nonfinite coefficient",
        "BE rejects matrix overflow");

  state = make_state(2, true);
  op = OperatorData(2);
  std::fill(op.pressure_weight.begin(), op.pressure_weight.end(), 1e-200);
  op.interface_exchange[1] = 1e-100;
  op.number_activation_rate[0][0] = 1e100;
  state.interstitial[0][1] = 1e300;
  result = ndrop_test::advance_backward_euler(state, op, 1e-200);
  check(!result.ok &&
            result.message ==
                "Backward Euler droplet activation source contains a nonfinite value",
        "BE rejects solved-state droplet-source overflow");

  state = make_state(2, true);
  op = OperatorData(2);
  op.interface_exchange[1] = 1;
  for (int mode = 0; mode < num_modes; ++mode) {
    op.number_activation_rate[mode][0] = 1;
    state.interstitial[expected_number_index[mode]][1] = 7e307;
  }
  result = ndrop_test::advance_backward_euler(state, op, 2);
  check(!result.ok &&
            result.message ==
                "Backward Euler droplet right-hand side contains a nonfinite value",
        "BE rejects droplet right-hand-side overflow: " + result.message);

  for (const auto probe :
       {ndrop_test::BackwardEulerProbe::solution,
        ndrop_test::BackwardEulerProbe::matrix_residual,
        ndrop_test::BackwardEulerProbe::physical_residual,
        ndrop_test::BackwardEulerProbe::ledger_normalization,
        ndrop_test::BackwardEulerProbe::ledger_decomposition,
        ndrop_test::BackwardEulerProbe::block_pivot,
        ndrop_test::BackwardEulerProbe::scalar_pivot,
        ndrop_test::BackwardEulerProbe::material_negative,
        ndrop_test::BackwardEulerProbe::roundoff_negative,
        ndrop_test::BackwardEulerProbe::bottom_branch}) {
    check(ndrop_test::backward_euler_probe_is_rejected(probe),
          "BE deterministic probe is rejected");
  }
}

} // namespace

int main() {
  try {
    test_identity_and_boundaries();
    test_dense_oracle_grid();
    test_matrix_family_reuse();
    test_controller_is_not_physical();
    test_uniform_and_zero_cloud();
    test_zero_overlap_with_mixing();
    test_serial_determinism();
    test_diffusion_convergence();
    test_activation_convergence();
    test_scaling_zero_ledgers_and_positivity();
    test_fail_closed_paths();
  } catch (const std::exception &error) {
    ++failed;
    std::cerr << "UNCAUGHT: " << error.what() << '\n';
  }
  std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
