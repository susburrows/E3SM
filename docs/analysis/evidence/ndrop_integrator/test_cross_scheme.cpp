#include "ndrop_integrator.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using ndrop_test::AdvanceResult;
using ndrop_test::AdvanceDiagnostics;
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::ReferenceOptions;
using ndrop_test::State;
using ndrop_test::num_constituents;
using ndrop_test::num_modes;

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
  if (condition) {
    ++passed;
  } else {
    ++failed;
    std::cerr << "FAIL: " << message << '\n';
  }
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

bool real_bitwise_equal(const Real lhs, const Real rhs) {
  return std::memcmp(&lhs, &rhs, sizeof(Real)) == 0;
}

bool field_diagnostics_bitwise_equal(
    const ndrop_test::FieldClipping &lhs,
    const ndrop_test::FieldClipping &rhs) {
  return real_bitwise_equal(lhs.interstitial_raw_minimum,
                            rhs.interstitial_raw_minimum) &&
         real_bitwise_equal(lhs.cloudborne_raw_minimum,
                            rhs.cloudborne_raw_minimum) &&
         real_bitwise_equal(lhs.interstitial_weighted_deficit,
                            rhs.interstitial_weighted_deficit) &&
         real_bitwise_equal(lhs.cloudborne_weighted_deficit,
                            rhs.cloudborne_weighted_deficit) &&
         lhs.interstitial_count == rhs.interstitial_count &&
         lhs.cloudborne_count == rhs.cloudborne_count &&
         lhs.interstitial_roundoff_negative_count ==
             rhs.interstitial_roundoff_negative_count &&
         lhs.interstitial_material_negative_count ==
             rhs.interstitial_material_negative_count &&
         lhs.cloudborne_roundoff_negative_count ==
             rhs.cloudborne_roundoff_negative_count &&
         lhs.cloudborne_material_negative_count ==
             rhs.cloudborne_material_negative_count;
}

bool diagnostics_bitwise_equal(const AdvanceDiagnostics &lhs,
                               const AdvanceDiagnostics &rhs) {
  const bool scalar_equal =
      lhs.nsubmix == rhs.nsubmix && lhs.final_buffer == rhs.final_buffer &&
      lhs.activation_limiter_count == rhs.activation_limiter_count &&
      lhs.clipping_count == rhs.clipping_count &&
      real_bitwise_equal(lhs.raw_minimum, rhs.raw_minimum) &&
      real_bitwise_equal(lhs.droplet_raw_minimum,
                         rhs.droplet_raw_minimum) &&
      real_bitwise_equal(lhs.droplet_weighted_deficit,
                         rhs.droplet_weighted_deficit) &&
      lhs.droplet_clipping_count == rhs.droplet_clipping_count &&
      real_bitwise_equal(lhs.max_inverse_mixing_time,
                         rhs.max_inverse_mixing_time) &&
      real_bitwise_equal(lhs.first_level_controller_rate,
                         rhs.first_level_controller_rate) &&
      real_bitwise_equal(lhs.normalized_residual,
                         rhs.normalized_residual) &&
      real_bitwise_equal(lhs.normalized_matrix_residual,
                         rhs.normalized_matrix_residual) &&
      real_bitwise_equal(lhs.normalized_physical_residual,
                         rhs.normalized_physical_residual) &&
      real_bitwise_equal(lhs.max_normalized_ledger_error,
                         rhs.max_normalized_ledger_error) &&
      real_bitwise_equal(lhs.solution_minimum, rhs.solution_minimum) &&
      lhs.negative_count == rhs.negative_count &&
      lhs.roundoff_negative_count == rhs.roundoff_negative_count &&
      lhs.material_negative_count == rhs.material_negative_count &&
      lhs.droplet_roundoff_negative_count ==
          rhs.droplet_roundoff_negative_count &&
      lhs.droplet_material_negative_count ==
          rhs.droplet_material_negative_count &&
      lhs.block_matrix_family_count == rhs.block_matrix_family_count &&
      lhs.block_factorization_count == rhs.block_factorization_count &&
      lhs.block_rhs_solve_count == rhs.block_rhs_solve_count &&
      lhs.scalar_factorization_count == rhs.scalar_factorization_count &&
      lhs.scalar_rhs_solve_count == rhs.scalar_rhs_solve_count &&
      lhs.positivity_failure == rhs.positivity_failure;
  if (!scalar_equal) {
    return false;
  }
  for (int field = 0; field < num_constituents; ++field) {
    if (!field_diagnostics_bitwise_equal(lhs.field_clipping[field],
                                         rhs.field_clipping[field]) ||
        !real_bitwise_equal(lhs.initial_ledger[field],
                            rhs.initial_ledger[field]) ||
        !real_bitwise_equal(lhs.final_ledger[field],
                            rhs.final_ledger[field]) ||
        !real_bitwise_equal(lhs.clipping_ledger_delta[field],
                            rhs.clipping_ledger_delta[field]) ||
        !real_bitwise_equal(lhs.total_ledger_delta[field],
                            rhs.total_ledger_delta[field]) ||
        !real_bitwise_equal(lhs.operator_ledger_delta[field],
                            rhs.operator_ledger_delta[field])) {
      return false;
    }
  }
  return true;
}

struct FieldError {
  Real weighted_l1 = 0;
  Real linf = 0;
};

FieldError field_error(const std::vector<Real> &actual,
                       const std::vector<Real> &reference,
                       const std::vector<Real> &initial,
                       const std::vector<Real> &weights) {
  if (actual.size() != reference.size() || actual.size() != initial.size() ||
      actual.size() != weights.size()) {
    throw std::logic_error("cross-scheme field error size mismatch");
  }
  FieldError error;
  Real reference_l1 = 0;
  Real initial_l1 = 0;
  Real reference_linf = 0;
  Real initial_linf = 0;
  for (std::size_t k = 0; k < actual.size(); ++k) {
    if (!std::isfinite(actual[k]) || !std::isfinite(reference[k]) ||
        !std::isfinite(initial[k]) || !std::isfinite(weights[k]) ||
        weights[k] < 0) {
      throw std::runtime_error("cross-scheme field error is nonfinite");
    }
    const Real absolute_error = std::abs(actual[k] - reference[k]);
    error.weighted_l1 += weights[k] * absolute_error;
    reference_l1 += weights[k] * std::abs(reference[k]);
    initial_l1 += weights[k] * std::abs(initial[k]);
    error.linf = std::max(error.linf, absolute_error);
    reference_linf = std::max(reference_linf, std::abs(reference[k]));
    initial_linf = std::max(initial_linf, std::abs(initial[k]));
    if (!std::isfinite(error.weighted_l1) ||
        !std::isfinite(reference_l1) || !std::isfinite(initial_l1) ||
        !std::isfinite(error.linf) || !std::isfinite(reference_linf) ||
        !std::isfinite(initial_linf)) {
      throw std::runtime_error("cross-scheme field error is nonfinite");
    }
  }
  const Real l1_scale = std::max(reference_l1, initial_l1);
  const Real linf_scale = std::max(reference_linf, initial_linf);
  if (l1_scale > 0) {
    error.weighted_l1 /= l1_scale;
  }
  if (linf_scale > 0) {
    error.linf /= linf_scale;
  }
  if (!std::isfinite(error.weighted_l1) || !std::isfinite(error.linf)) {
    throw std::runtime_error("cross-scheme field error is nonfinite");
  }
  return error;
}

std::array<FieldError, 51>
all_field_errors(const State &actual, const State &reference,
                 const State &initial, const std::vector<Real> &weights) {
  std::array<FieldError, 51> errors{};
  errors[0] = field_error(actual.droplet_number, reference.droplet_number,
                          initial.droplet_number, weights);
  for (int j = 0; j < num_constituents; ++j) {
    errors[1 + 2 * j] =
        field_error(actual.interstitial[j], reference.interstitial[j],
                    initial.interstitial[j], weights);
    errors[2 + 2 * j] =
        field_error(actual.cloudborne[j], reference.cloudborne[j],
                    initial.cloudborne[j], weights);
  }
  return errors;
}

std::pair<Real, Real>
maximum_errors(const std::array<FieldError, 51> &errors) {
  Real weighted_l1 = 0;
  Real linf = 0;
  for (const auto &error : errors) {
    if (!std::isfinite(error.weighted_l1) || !std::isfinite(error.linf)) {
      throw std::runtime_error("cross-scheme field error is nonfinite");
    }
    weighted_l1 = std::max(weighted_l1, error.weighted_l1);
    linf = std::max(linf, error.linf);
  }
  return {weighted_l1, linf};
}

Real next_unit_random(unsigned long long &state) {
  state = state * 6364136223846793005ULL + 1442695040888963407ULL;
  return static_cast<Real>(state >> 11) /
         static_cast<Real>(1ULL << 53);
}

struct CommonCase {
  int id = 0;
  std::size_t nlev = 0;
  Real overlap = 0;
  Real activation = 0;
  int active_pairs = 0;
  int mixing_index = 0;
  Real mixing_scale = 0;
  Real dt = 0.5;
  State initial;
  OperatorData op;
};

CommonCase make_common_case(const int id, const std::size_t nlev,
                            const Real overlap, const Real activation,
                            const int active_pairs, const int mixing_index,
                            const Real mixing_scale) {
  CommonCase test_case;
  test_case.id = id;
  test_case.nlev = nlev;
  test_case.overlap = overlap;
  test_case.activation = activation;
  test_case.active_pairs = active_pairs;
  test_case.mixing_index = mixing_index;
  test_case.mixing_scale = mixing_scale;
  test_case.initial = State(nlev);
  test_case.op = OperatorData(nlev);
  unsigned long long seed =
      0x4e44524fULL ^ (static_cast<unsigned long long>(id) << 24);

  for (std::size_t k = 0; k < nlev; ++k) {
    test_case.op.pressure_weight[k] =
        0.8 + 0.4 * static_cast<Real>((17 * k + 3) % 11) / 10;
    test_case.op.cloud_fraction[k] =
        k % 2 == 0 ? Real(0.8) : Real(0.8) * overlap;
    test_case.initial.droplet_number[k] =
        4 + 0.01 * static_cast<Real>(k) + 0.002 * next_unit_random(seed);
    for (int j = 0; j < active_pairs; ++j) {
      test_case.initial.interstitial[j][k] =
          1 + 0.01 * j + 0.002 * static_cast<Real>(k) +
          0.01 * next_unit_random(seed);
      test_case.initial.cloudborne[j][k] =
          0.2 + 0.002 * j + 0.001 * static_cast<Real>(k) +
          0.005 * next_unit_random(seed);
    }
  }
  for (std::size_t interface = 1; interface < nlev; ++interface) {
    const int offset = static_cast<int>((7 * interface + 1) % 5) - 2;
    const Real modulation = 1 + 0.05 * static_cast<Real>(offset);
    test_case.op.interface_exchange[interface] =
        mixing_scale * modulation *
        std::min(test_case.op.pressure_weight[interface - 1],
                 test_case.op.pressure_weight[interface]);
  }
  constexpr std::array<Real, num_modes> number_multiplier = {0.7, 0.8, 0.9,
                                                              1.0};
  constexpr std::array<Real, num_modes> mass_multiplier = {0.55, 0.65, 0.75,
                                                            0.85};
  for (int mode = 0; mode < num_modes; ++mode) {
    for (std::size_t k = 0; k + 1 < nlev; ++k) {
      const Real rate_below = test_case.op.interface_exchange[k + 1] /
                              test_case.op.pressure_weight[k];
      test_case.op.number_activation_rate[mode][k] =
          activation * number_multiplier[mode] * rate_below;
      test_case.op.mass_activation_rate[mode][k] =
          activation * mass_multiplier[mode] * rate_below;
    }
  }
  test_case.op.first_level_controller_rate = 0.25 * mixing_scale;
  return test_case;
}

enum class Method { forward_euler, imex, backward_euler };

const char *method_name(const Method method) {
  switch (method) {
  case Method::forward_euler:
    return "forward_euler";
  case Method::imex:
    return "imex";
  case Method::backward_euler:
    return "backward_euler";
  }
  return "unknown";
}

AdvanceResult advance(const Method method, const State &state,
                      const OperatorData &op, const Real dt,
                      const bool fixed_fe_step = false) {
  switch (method) {
  case Method::forward_euler: {
    ndrop_test::ForwardEulerOptions options;
    options.forced_substeps = fixed_fe_step ? 1 : 0;
    return ndrop_test::advance_forward_euler(state, op, dt, options);
  }
  case Method::imex:
    return ndrop_test::advance_imex(state, op, dt);
  case Method::backward_euler:
    return ndrop_test::advance_backward_euler(state, op, dt);
  }
  throw std::logic_error("unrecognized method");
}

struct MatrixRecord {
  int case_id = 0;
  std::size_t nlev = 0;
  Real overlap = 0;
  Real activation = 0;
  int active_pairs = 0;
  int mixing_index = 0;
  Real mixing_scale = 0;
  Real dt = 0;
  Method method = Method::forward_euler;
  bool passed = false;
  bool deterministic = false;
  Real reference_consistency = 0;
  Real weighted_l1_error = 0;
  Real linf_error = 0;
  Real residual = 0;
  Real ledger_error = 0;
  Real raw_minimum = 0;
  int material_negative_count = 0;
  int clipping_count = 0;
  Real clipping_deficit = 0;
  int nsubmix = 0;
  std::array<FieldError, 51> field_errors{};
  AdvanceDiagnostics diagnostics;
  std::string message;
};

struct MethodSummary {
  int cases = 0;
  int failures = 0;
  int nondeterministic = 0;
  int material_negative_count = 0;
  int clipping_cases = 0;
  Real maximum_reference_consistency = 0;
  Real maximum_weighted_l1_error = 0;
  Real maximum_linf_error = 0;
  Real maximum_residual = 0;
  Real maximum_ledger_error = 0;
  Real minimum_raw_value = std::numeric_limits<Real>::infinity();
  Real maximum_clipping_deficit = 0;
};

void update_summary(MethodSummary &summary, const MatrixRecord &record) {
  ++summary.cases;
  if (!record.passed) {
    ++summary.failures;
  }
  if (!record.deterministic) {
    ++summary.nondeterministic;
  }
  summary.material_negative_count += record.material_negative_count;
  if (record.clipping_count > 0) {
    ++summary.clipping_cases;
  }
  summary.maximum_reference_consistency =
      std::max(summary.maximum_reference_consistency,
               record.reference_consistency);
  summary.maximum_weighted_l1_error =
      std::max(summary.maximum_weighted_l1_error, record.weighted_l1_error);
  summary.maximum_linf_error =
      std::max(summary.maximum_linf_error, record.linf_error);
  summary.maximum_residual =
      std::max(summary.maximum_residual, record.residual);
  summary.maximum_ledger_error =
      std::max(summary.maximum_ledger_error, record.ledger_error);
  summary.minimum_raw_value =
      std::min(summary.minimum_raw_value, record.raw_minimum);
  summary.maximum_clipping_deficit =
      std::max(summary.maximum_clipping_deficit, record.clipping_deficit);
}

struct MatrixResults {
  std::vector<MatrixRecord> records;
  std::array<MethodSummary, 3> summaries{};
  int reference_failures = 0;
  Real maximum_reference_consistency = 0;
};

std::size_t method_index(const Method method) {
  return static_cast<std::size_t>(method);
}

MatrixResults run_common_matrix(const bool smoke) {
  const std::vector<std::size_t> levels =
      smoke ? std::vector<std::size_t>{3} :
              std::vector<std::size_t>{1, 2, 3, 127, 128};
  const std::vector<Real> overlaps =
      smoke ? std::vector<Real>{0.5} :
              std::vector<Real>{0, 0.25, 0.5, 0.75, 1};
  const std::vector<Real> activations =
      smoke ? std::vector<Real>{0.9} :
              std::vector<Real>{0, 0.1, 0.5, 0.9, 2};
  const std::vector<int> pair_counts =
      smoke ? std::vector<int>{4} : std::vector<int>{1, 4, 8, 16, 25};
  const std::vector<Real> mixing_scales =
      smoke ? std::vector<Real>{0.5} : std::vector<Real>{0, 0.05, 0.5};
  constexpr std::array<Method, 3> methods = {
      Method::forward_euler, Method::imex, Method::backward_euler};

  MatrixResults results;
  int case_id = 0;
  for (const auto nlev : levels) {
    for (const Real overlap : overlaps) {
      for (const Real activation : activations) {
        for (const int active_pairs : pair_counts) {
          for (std::size_t mixing = 0; mixing < mixing_scales.size();
               ++mixing) {
            const CommonCase test_case = make_common_case(
                case_id, nlev, overlap, activation, active_pairs,
                static_cast<int>(mixing), mixing_scales[mixing]);
            ReferenceOptions normal_options;
            ReferenceOptions tight_options;
            tight_options.relative_tolerance = 1e-13;
            const auto reference = ndrop_test::advance_reference_dopri54(
                test_case.initial, test_case.op, test_case.dt, normal_options);
            const auto tight_reference = ndrop_test::advance_reference_dopri54(
                test_case.initial, test_case.op, test_case.dt, tight_options);
            check(reference.ok,
                  "common-matrix reference succeeds case " +
                      std::to_string(case_id) + ": " + reference.message);
            check(tight_reference.ok,
                  "common-matrix tight reference succeeds case " +
                      std::to_string(case_id) + ": " +
                      tight_reference.message);
            Real reference_consistency = std::numeric_limits<Real>::infinity();
            if (reference.ok && tight_reference.ok) {
              const auto reference_errors = all_field_errors(
                  reference.state, tight_reference.state, test_case.initial,
                  test_case.op.pressure_weight);
              const auto maxima = maximum_errors(reference_errors);
              reference_consistency = std::max(maxima.first, maxima.second);
              results.maximum_reference_consistency =
                  std::max(results.maximum_reference_consistency,
                           reference_consistency);
              check(reference_consistency <= 2e-11,
                    "reference two-tolerance self-check case " +
                        std::to_string(case_id));
            } else {
              ++results.reference_failures;
            }

            for (const Method method : methods) {
              MatrixRecord record;
              record.case_id = case_id;
              record.nlev = nlev;
              record.overlap = overlap;
              record.activation = activation;
              record.active_pairs = active_pairs;
              record.mixing_index = static_cast<int>(mixing);
              record.mixing_scale = mixing_scales[mixing];
              record.dt = test_case.dt;
              record.method = method;
              record.reference_consistency = reference_consistency;
              const auto first = advance(method, test_case.initial,
                                         test_case.op, test_case.dt);
              const auto second = advance(method, test_case.initial,
                                          test_case.op, test_case.dt);
              const auto third = advance(method, test_case.initial,
                                         test_case.op, test_case.dt);
              record.message = first.message;
              record.deterministic =
                  first.ok == second.ok && first.ok == third.ok &&
                  first.message == second.message &&
                  first.message == third.message &&
                  state_bitwise_equal(first.state, second.state) &&
                  state_bitwise_equal(first.state, third.state) &&
                  diagnostics_bitwise_equal(first.diagnostics,
                                            second.diagnostics) &&
                  diagnostics_bitwise_equal(first.diagnostics,
                                            third.diagnostics);
              record.diagnostics = first.diagnostics;
              record.residual = first.diagnostics.normalized_residual;
              record.ledger_error =
                  first.diagnostics.max_normalized_ledger_error;
              record.raw_minimum =
                  method == Method::forward_euler
                      ? first.diagnostics.raw_minimum
                      : first.diagnostics.solution_minimum;
              record.material_negative_count =
                  first.diagnostics.material_negative_count;
              record.clipping_count = first.diagnostics.clipping_count;
              record.nsubmix = first.diagnostics.nsubmix;
              record.clipping_deficit =
                  first.diagnostics.droplet_weighted_deficit;
              for (const auto &field : first.diagnostics.field_clipping) {
                record.clipping_deficit +=
                    field.interstitial_weighted_deficit +
                    field.cloudborne_weighted_deficit;
              }
              if (first.ok && reference.ok) {
                record.field_errors = all_field_errors(
                    first.state, reference.state, test_case.initial,
                    test_case.op.pressure_weight);
                const auto maxima = maximum_errors(record.field_errors);
                record.weighted_l1_error = maxima.first;
                record.linf_error = maxima.second;
              }
              const Real ledger_tolerance =
                  std::max(Real(1e-11),
                           5 * static_cast<Real>(nlev) *
                               std::numeric_limits<Real>::epsilon());
              record.passed =
                  reference.ok && tight_reference.ok &&
                  reference_consistency <= 2e-11 && first.ok &&
                  record.deterministic && record.residual <= 1e-11 &&
                  record.ledger_error <= ledger_tolerance &&
                  first.diagnostics.material_negative_count == 0;
              check(record.passed,
                    std::string(method_name(method)) +
                        " common-matrix gate case " +
                        std::to_string(case_id) + ": " + first.message);
              results.records.push_back(record);
              update_summary(results.summaries[method_index(method)], record);
            }
            ++case_id;
          }
        }
      }
    }
  }
  const int expected_cases = smoke ? 1 : 1875;
  check(case_id == expected_cases, "complete common-matrix case count");
  for (const Method method : methods) {
    check(results.summaries[method_index(method)].cases == expected_cases,
          std::string(method_name(method)) + " complete matrix coverage");
  }
  return results;
}

State make_smooth_initial(const std::size_t nlev) {
  State state(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    state.droplet_number[k] = 4 + 0.03 * k + 0.1 * std::sin(0.4 * k);
    for (int j = 0; j < num_constituents; ++j) {
      state.interstitial[j][k] =
          1 + 0.01 * j + 0.005 * k + 0.02 * std::cos(0.3 * k + 0.1 * j);
      state.cloudborne[j][k] =
          0.3 + 0.003 * j + 0.002 * k +
          0.01 * std::sin(0.2 * k + 0.05 * j);
    }
  }
  return state;
}

OperatorData make_smooth_operator(const std::size_t nlev) {
  OperatorData op(nlev);
  constexpr std::array<Real, 5> clouds = {0.25, 0.4, 0.6, 0.8, 1};
  for (std::size_t k = 0; k < nlev; ++k) {
    op.pressure_weight[k] = 0.7 + 0.03 * k;
    op.cloud_fraction[k] = clouds[k % clouds.size()];
  }
  for (std::size_t interface = 1; interface < nlev; ++interface) {
    const Real rate = 0.08 + 0.015 * static_cast<Real>(interface % 4);
    op.interface_exchange[interface] =
        rate * std::min(op.pressure_weight[interface - 1],
                        op.pressure_weight[interface]);
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    for (std::size_t k = 0; k + 1 < nlev; ++k) {
      const Real rate_below =
          op.interface_exchange[k + 1] / op.pressure_weight[k];
      op.number_activation_rate[mode][k] =
          (0.15 + 0.03 * mode) * rate_below;
      op.mass_activation_rate[mode][k] =
          (0.08 + 0.02 * mode) * rate_below;
    }
  }
  return op;
}

State integrate_method(const Method method, State state,
                       const OperatorData &op, const Real total_time,
                       const int steps) {
  const Real dt = total_time / steps;
  for (int step = 0; step < steps; ++step) {
    const auto result = advance(method, state, op, dt, true);
    if (!result.ok) {
      throw std::runtime_error(std::string(method_name(method)) +
                               " convergence advance failed: " +
                               result.message);
    }
    state = result.state;
  }
  return state;
}

struct Fit {
  Real slope = 0;
  Real r2 = 0;
};

Fit fit_three_finest(const std::array<Real, 5> &steps,
                     const std::array<Real, 5> &errors) {
  Real mean_x = 0;
  Real mean_y = 0;
  for (std::size_t i = 2; i < 5; ++i) {
    if (!(errors[i] > 0) || !std::isfinite(errors[i])) {
      throw std::runtime_error("convergence error reached zero or nonfinite");
    }
    mean_x += std::log(steps[i]);
    mean_y += std::log(errors[i]);
  }
  mean_x /= 3;
  mean_y /= 3;
  Real covariance = 0;
  Real variance_x = 0;
  Real variance_y = 0;
  for (std::size_t i = 2; i < 5; ++i) {
    const Real dx = std::log(steps[i]) - mean_x;
    const Real dy = std::log(errors[i]) - mean_y;
    covariance += dx * dy;
    variance_x += dx * dx;
    variance_y += dy * dy;
  }
  Fit fit;
  fit.slope = covariance / variance_x;
  fit.r2 = covariance * covariance / (variance_x * variance_y);
  return fit;
}

struct ConvergenceSummary {
  Method method = Method::forward_euler;
  Real minimum_slope = std::numeric_limits<Real>::infinity();
  Real maximum_slope = -std::numeric_limits<Real>::infinity();
  Real minimum_r2 = std::numeric_limits<Real>::infinity();
  Real finest_weighted_l1 = 0;
  Real finest_linf = 0;
  std::array<Fit, 51> weighted_l1_fit{};
  std::array<Fit, 51> linf_fit{};
  bool passed = true;
};

std::array<ConvergenceSummary, 3> test_common_convergence() {
  constexpr std::size_t nlev = 12;
  constexpr Real total_time = 1;
  constexpr std::array<int, 5> counts = {16, 32, 64, 128, 256};
  constexpr std::array<Method, 3> methods = {
      Method::forward_euler, Method::imex, Method::backward_euler};
  const State initial = make_smooth_initial(nlev);
  const OperatorData op = make_smooth_operator(nlev);
  ReferenceOptions tight;
  tight.relative_tolerance = 1e-13;
  const auto reference =
      ndrop_test::advance_reference_dopri54(initial, op, total_time);
  const auto tight_reference =
      ndrop_test::advance_reference_dopri54(initial, op, total_time, tight);
  check(reference.ok, "coupled reference succeeds: " + reference.message);
  check(tight_reference.ok,
        "coupled tight reference succeeds: " + tight_reference.message);
  if (!reference.ok || !tight_reference.ok) {
    throw std::runtime_error("coupled reference unavailable");
  }
  const auto consistency = maximum_errors(all_field_errors(
      reference.state, tight_reference.state, initial, op.pressure_weight));
  check(std::max(consistency.first, consistency.second) <= 2e-11,
        "coupled reference two-tolerance consistency");

  std::array<Real, 5> step_sizes{};
  std::array<std::array<State, 5>, 3> solutions;
  std::array<ConvergenceSummary, 3> summaries;
  for (const Method method : methods) {
    auto &summary = summaries[method_index(method)];
    summary.method = method;
    std::array<std::array<Real, 5>, 51> l1{};
    std::array<std::array<Real, 5>, 51> linf{};
    for (std::size_t i = 0; i < counts.size(); ++i) {
      step_sizes[i] = total_time / counts[i];
      solutions[method_index(method)][i] =
          integrate_method(method, initial, op, total_time, counts[i]);
      const auto errors = all_field_errors(
          solutions[method_index(method)][i], reference.state, initial,
          op.pressure_weight);
      for (std::size_t field = 0; field < errors.size(); ++field) {
        l1[field][i] = errors[field].weighted_l1;
        linf[field][i] = errors[field].linf;
      }
    }
    for (std::size_t field = 0; field < 51; ++field) {
      summary.weighted_l1_fit[field] =
          fit_three_finest(step_sizes, l1[field]);
      summary.linf_fit[field] =
          fit_three_finest(step_sizes, linf[field]);
      for (const Fit fit : {summary.weighted_l1_fit[field],
                            summary.linf_fit[field]}) {
        const bool fit_pass = fit.slope >= 0.85 && fit.slope <= 1.15 &&
                              fit.r2 >= 0.98;
        check(fit_pass,
              std::string(method_name(method)) +
                  " coupled first-order field " + std::to_string(field));
        summary.passed = summary.passed && fit_pass;
        summary.minimum_slope = std::min(summary.minimum_slope, fit.slope);
        summary.maximum_slope = std::max(summary.maximum_slope, fit.slope);
        summary.minimum_r2 = std::min(summary.minimum_r2, fit.r2);
      }
      summary.finest_weighted_l1 =
          std::max(summary.finest_weighted_l1, l1[field].back());
      summary.finest_linf =
          std::max(summary.finest_linf, linf[field].back());
    }
  }

  for (std::size_t lhs = 0; lhs < methods.size(); ++lhs) {
    for (std::size_t rhs = lhs + 1; rhs < methods.size(); ++rhs) {
      const auto cross_coarse = all_field_errors(
          solutions[lhs][3], solutions[rhs][3], initial,
          op.pressure_weight);
      const auto cross_fine = all_field_errors(
          solutions[lhs][4], solutions[rhs][4], initial,
          op.pressure_weight);
      const auto richardson_lhs = all_field_errors(
          solutions[lhs][3], solutions[lhs][4], initial,
          op.pressure_weight);
      const auto richardson_rhs = all_field_errors(
          solutions[rhs][3], solutions[rhs][4], initial,
          op.pressure_weight);
      for (std::size_t field = 0; field < 51; ++field) {
        const bool decreases =
            cross_fine[field].weighted_l1 <
                cross_coarse[field].weighted_l1 &&
            cross_fine[field].linf < cross_coarse[field].linf;
        const bool bounded =
            cross_fine[field].weighted_l1 <=
                2 * std::max(richardson_lhs[field].weighted_l1,
                             richardson_rhs[field].weighted_l1) &&
            cross_fine[field].linf <=
                2 * std::max(richardson_lhs[field].linf,
                             richardson_rhs[field].linf);
        check(decreases,
              std::string(method_name(methods[lhs])) + "/" +
                  method_name(methods[rhs]) +
                  " finest difference decreases field " +
                  std::to_string(field));
        check(bounded,
              std::string(method_name(methods[lhs])) + "/" +
                  method_name(methods[rhs]) +
                  " Richardson common-solution bound field " +
                  std::to_string(field));
        summaries[lhs].passed = summaries[lhs].passed && decreases && bounded;
        summaries[rhs].passed = summaries[rhs].passed && decreases && bounded;
      }
    }
  }
  return summaries;
}

State make_analytic_state(const std::size_t nlev) {
  State state(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    state.droplet_number[k] = 3 + 0.4 * k;
    for (int j = 0; j < num_constituents; ++j) {
      state.interstitial[j][k] = 1 + 0.02 * j + 0.1 * k;
      state.cloudborne[j][k] = 0.2 + 0.01 * j + 0.05 * k;
    }
  }
  return state;
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
  Real coefficient_a = 0;
  Real coefficient_b = 0;
  Real rate_a = 0;
  Real rate_b = 0;
};

PairSolution exact_activation_pair(
    const std::array<Real, 2> &initial_interstitial,
    const std::array<Real, 2> &initial_cloudborne, const Real d,
    const Real activation, const Real total_time) {
  const Real s = std::sqrt(d * (d - activation));
  const Real rate_a = -d + s;
  const Real rate_b = -d - s;
  const Real coefficient_a =
      0.5 * (initial_interstitial[1] + d * initial_interstitial[0] / s);
  const Real coefficient_b =
      0.5 * (initial_interstitial[1] - d * initial_interstitial[0] / s);
  PairSolution solution;
  solution.interstitial =
      {std::exp(-d * total_time) *
           (std::cosh(s * total_time) * initial_interstitial[0] +
            (d - activation) * std::sinh(s * total_time) *
                initial_interstitial[1] / s),
       coefficient_a * std::exp(rate_a * total_time) +
           coefficient_b * std::exp(rate_b * total_time)};
  solution.cloudborne =
      forced_cloud_pair(initial_cloudborne, d, activation, total_time,
                        coefficient_a, coefficient_b, rate_a, rate_b);
  solution.coefficient_a = coefficient_a;
  solution.coefficient_b = coefficient_b;
  solution.rate_a = rate_a;
  solution.rate_b = rate_b;
  return solution;
}

void test_reference_against_analytic_solutions() {
  constexpr Real d = 0.5;
  constexpr Real total_time = 1;
  State initial = make_analytic_state(2);
  OperatorData op(2);
  op.interface_exchange[1] = d;
  State exact = initial;
  const Real decay = std::exp(-2 * d * total_time);
  auto diffuse = [decay](const std::vector<Real> &field) {
    return std::vector<Real>{
        0.5 * (field[0] + field[1] + decay * (field[0] - field[1])),
        0.5 * (field[0] + field[1] - decay * (field[0] - field[1]))};
  };
  exact.droplet_number = diffuse(initial.droplet_number);
  for (int j = 0; j < num_constituents; ++j) {
    exact.interstitial[j] = diffuse(initial.interstitial[j]);
    exact.cloudborne[j] = diffuse(initial.cloudborne[j]);
  }
  auto reference =
      ndrop_test::advance_reference_dopri54(initial, op, total_time);
  check(reference.ok, "analytic diffusion reference succeeds");
  auto error = maximum_errors(all_field_errors(
      reference.state, exact, initial, op.pressure_weight));
  check(std::max(error.first, error.second) <= 5e-11,
        "reference agrees with analytic no-flux diffusion");

  constexpr std::array<Real, num_modes> number_activation = {0.11, 0.17,
                                                              0.23, 0.29};
  constexpr std::array<Real, num_modes> mass_activation = {0.07, 0.13, 0.19,
                                                            0.25};
  initial = make_analytic_state(2);
  exact = initial;
  op = OperatorData(2);
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
  exact.droplet_number = diffuse(initial.droplet_number);
  for (const auto &forcing : number_forcing) {
    exact.droplet_number[0] += forcing[0];
    exact.droplet_number[1] += forcing[1];
  }
  reference = ndrop_test::advance_reference_dopri54(initial, op, total_time);
  check(reference.ok, "analytic activation reference succeeds");
  error = maximum_errors(all_field_errors(reference.state, exact, initial,
                                          op.pressure_weight));
  check(std::max(error.first, error.second) <= 5e-11,
        "reference agrees with analytic activation transfer");
}

void test_reference_failure_and_determinism_contract() {
  const State initial = make_analytic_state(2);
  OperatorData op(2);
  op.interface_exchange[1] = 0.5;
  const auto first = ndrop_test::advance_reference_dopri54(initial, op, 1);
  const auto second = ndrop_test::advance_reference_dopri54(initial, op, 1);
  check(first.ok && second.ok && state_bitwise_equal(first.state, second.state) &&
            first.accepted_steps == second.accepted_steps &&
            first.rejected_steps == second.rejected_steps &&
            first.function_evaluations == second.function_evaluations &&
            real_bitwise_equal(first.maximum_error_ratio,
                               second.maximum_error_ratio),
        "reference is bitwise deterministic");

  ReferenceOptions options;
  options.relative_tolerance = 0;
  auto result =
      ndrop_test::advance_reference_dopri54(initial, op, 1, options);
  check(!result.ok &&
            result.message == "reference options must be finite and positive",
        "reference rejects a nonpositive tolerance");

  options = ReferenceOptions{};
  options.maximum_attempted_steps = 1;
  result = ndrop_test::advance_reference_dopri54(initial, op, 1, options);
  check(!result.ok &&
            result.message == "reference exceeded maximum attempted steps",
        "reference enforces the attempted-step limit");

  op = OperatorData(2);
  op.pressure_weight[0] = 1e-200;
  op.interface_exchange[1] = 1e200;
  result = ndrop_test::advance_reference_dopri54(initial, op, 1);
  check(!result.ok &&
            result.message ==
                "reference operator contains a nonfinite coefficient",
        "reference rejects derived-coefficient overflow");

  State extreme = make_analytic_state(2);
  extreme.droplet_number = {1e200, 2e200};
  op = OperatorData(2);
  op.interface_exchange[1] = 1e200;
  result = ndrop_test::advance_reference_dopri54(extreme, op, 1e100);
  check(!result.ok && result.message == "reference stage 1 is nonfinite",
        "reference rejects physical-RHS overflow");

  extreme = State(1);
  extreme.droplet_number[0] = 1;
  extreme.interstitial[0][0] = std::numeric_limits<Real>::max();
  extreme.cloudborne[0][0] = std::numeric_limits<Real>::max();
  op = OperatorData(1);
  op.cloud_fraction[0] = 0;
  result = ndrop_test::advance_reference_dopri54(extreme, op, 1);
  check(!result.ok &&
            result.message == "reference zero-cloud postprocess overflowed",
        "reference rejects projection overflow");
}

void test_error_norm_failure_contract() {
  State actual = make_analytic_state(2);
  const State reference = actual;
  const State initial = actual;
  actual.droplet_number[0] = std::numeric_limits<Real>::quiet_NaN();
  bool rejected = false;
  try {
    static_cast<void>(all_field_errors(actual, reference, initial, {1, 1}));
  } catch (const std::runtime_error &error) {
    rejected = std::string(error.what()) ==
               "cross-scheme field error is nonfinite";
  }
  check(rejected, "cross-scheme field errors fail closed on nonfinite input");
}

std::string csv_escape(const std::string &value) {
  if (value.find_first_of(",\"\n") == std::string::npos) {
    return value;
  }
  std::string escaped = "\"";
  for (const char character : value) {
    if (character == '\"') {
      escaped += "\"\"";
    } else {
      escaped += character;
    }
  }
  escaped += '\"';
  return escaped;
}

std::string real_text(const Real value) {
  std::ostringstream text;
  text << std::setprecision(17) << value;
  return text.str();
}

std::string finite_real_text(const Real value) {
  return std::isfinite(value) ? real_text(value) : std::string{};
}

template <typename Getter>
std::string joined_values(const std::size_t count, Getter getter) {
  std::ostringstream values;
  for (std::size_t index = 0; index < count; ++index) {
    if (index > 0) {
      values << ';';
    }
    values << getter(index);
  }
  return values.str();
}

void write_csv(const std::string &path, const MatrixResults &matrix,
               const std::array<ConvergenceSummary, 3> &convergence) {
  if (path.empty()) {
    return;
  }
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error("cannot open CSV output: " + path);
  }
  const std::vector<std::string> headers = {
      "record_type",
      "case_id",
      "nlev",
      "overlap",
      "activation",
      "active_pairs",
      "transported_fields",
      "mixing_index",
      "mixing_scale",
      "dt",
      "method",
      "passed",
      "deterministic",
      "reference_consistency",
      "maximum_weighted_l1_error",
      "maximum_linf_error",
      "residual",
      "ledger_error",
      "raw_minimum",
      "material_negative_count",
      "clipping_count",
      "clipping_deficit",
      "nsubmix",
      "weighted_l1_by_field",
      "linf_by_field",
      "initial_ledger_by_pair",
      "final_ledger_by_pair",
      "total_ledger_delta_by_pair",
      "clipping_ledger_delta_by_pair",
      "operator_ledger_delta_by_pair",
      "interstitial_raw_minimum_by_pair",
      "cloudborne_raw_minimum_by_pair",
      "interstitial_clipping_count_by_pair",
      "cloudborne_clipping_count_by_pair",
      "interstitial_weighted_deficit_by_pair",
      "cloudborne_weighted_deficit_by_pair",
      "interstitial_roundoff_negative_count_by_pair",
      "interstitial_material_negative_count_by_pair",
      "cloudborne_roundoff_negative_count_by_pair",
      "cloudborne_material_negative_count_by_pair",
      "droplet_raw_minimum",
      "droplet_clipping_count",
      "droplet_weighted_deficit",
      "droplet_roundoff_negative_count",
      "droplet_material_negative_count",
      "weighted_l1_slope_by_field",
      "weighted_l1_r2_by_field",
      "linf_slope_by_field",
      "linf_r2_by_field",
      "message"};
  auto write_row = [&output](const std::vector<std::string> &row) {
    for (std::size_t column = 0; column < row.size(); ++column) {
      if (column > 0) {
        output << ',';
      }
      output << csv_escape(row[column]);
    }
    output << '\n';
  };
  auto set_column = [&headers](std::vector<std::string> &row,
                               const std::string &name,
                               const std::string &value) {
    const auto found = std::find(headers.begin(), headers.end(), name);
    if (found == headers.end()) {
      throw std::logic_error("unknown CSV column: " + name);
    }
    row[static_cast<std::size_t>(found - headers.begin())] = value;
  };
  write_row(headers);
  for (const auto &record : matrix.records) {
    std::vector<std::string> row(headers.size());
    set_column(row, "record_type", "case");
    set_column(row, "case_id", std::to_string(record.case_id));
    set_column(row, "nlev", std::to_string(record.nlev));
    set_column(row, "overlap", real_text(record.overlap));
    set_column(row, "activation", real_text(record.activation));
    set_column(row, "active_pairs", std::to_string(record.active_pairs));
    set_column(row, "transported_fields",
               std::to_string(1 + 2 * record.active_pairs));
    set_column(row, "mixing_index", std::to_string(record.mixing_index));
    set_column(row, "mixing_scale", real_text(record.mixing_scale));
    set_column(row, "dt", real_text(record.dt));
    set_column(row, "method", method_name(record.method));
    set_column(row, "passed", record.passed ? "true" : "false");
    set_column(row, "deterministic",
               record.deterministic ? "true" : "false");
    set_column(row, "reference_consistency",
               real_text(record.reference_consistency));
    set_column(row, "maximum_weighted_l1_error",
               real_text(record.weighted_l1_error));
    set_column(row, "maximum_linf_error", real_text(record.linf_error));
    set_column(row, "residual", real_text(record.residual));
    set_column(row, "ledger_error", real_text(record.ledger_error));
    set_column(row, "raw_minimum", real_text(record.raw_minimum));
    set_column(row, "material_negative_count",
               std::to_string(record.material_negative_count));
    set_column(row, "clipping_count", std::to_string(record.clipping_count));
    set_column(row, "clipping_deficit", real_text(record.clipping_deficit));
    set_column(row, "nsubmix", std::to_string(record.nsubmix));
    set_column(row, "weighted_l1_by_field",
               joined_values(51, [&record](const std::size_t field) {
                 return real_text(record.field_errors[field].weighted_l1);
               }));
    set_column(row, "linf_by_field",
               joined_values(51, [&record](const std::size_t field) {
                 return real_text(record.field_errors[field].linf);
               }));
    const auto real_diagnostic_array =
        [](const std::array<Real, num_constituents> &values) {
          return joined_values(num_constituents,
                               [&values](const std::size_t field) {
                                 return real_text(values[field]);
                               });
        };
    set_column(row, "initial_ledger_by_pair",
               real_diagnostic_array(record.diagnostics.initial_ledger));
    set_column(row, "final_ledger_by_pair",
               real_diagnostic_array(record.diagnostics.final_ledger));
    set_column(row, "total_ledger_delta_by_pair",
               real_diagnostic_array(
                   record.diagnostics.total_ledger_delta));
    set_column(row, "clipping_ledger_delta_by_pair",
               real_diagnostic_array(
                   record.diagnostics.clipping_ledger_delta));
    set_column(row, "operator_ledger_delta_by_pair",
               real_diagnostic_array(
                   record.diagnostics.operator_ledger_delta));
    set_column(
        row, "interstitial_raw_minimum_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return finite_real_text(record.diagnostics.field_clipping[field]
                                      .interstitial_raw_minimum);
        }));
    set_column(
        row, "cloudborne_raw_minimum_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return finite_real_text(record.diagnostics.field_clipping[field]
                                      .cloudborne_raw_minimum);
        }));
    set_column(
        row, "interstitial_clipping_count_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return record.diagnostics.field_clipping[field].interstitial_count;
        }));
    set_column(
        row, "cloudborne_clipping_count_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return record.diagnostics.field_clipping[field].cloudborne_count;
        }));
    set_column(
        row, "interstitial_weighted_deficit_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return real_text(record.diagnostics.field_clipping[field]
                               .interstitial_weighted_deficit);
        }));
    set_column(
        row, "cloudborne_weighted_deficit_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return real_text(record.diagnostics.field_clipping[field]
                               .cloudborne_weighted_deficit);
        }));
    set_column(
        row, "interstitial_roundoff_negative_count_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return record.diagnostics.field_clipping[field]
              .interstitial_roundoff_negative_count;
        }));
    set_column(
        row, "interstitial_material_negative_count_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return record.diagnostics.field_clipping[field]
              .interstitial_material_negative_count;
        }));
    set_column(
        row, "cloudborne_roundoff_negative_count_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return record.diagnostics.field_clipping[field]
              .cloudborne_roundoff_negative_count;
        }));
    set_column(
        row, "cloudborne_material_negative_count_by_pair",
        joined_values(num_constituents, [&record](const std::size_t field) {
          return record.diagnostics.field_clipping[field]
              .cloudborne_material_negative_count;
        }));
    set_column(row, "droplet_raw_minimum",
               finite_real_text(record.diagnostics.droplet_raw_minimum));
    set_column(row, "droplet_clipping_count",
               std::to_string(record.diagnostics.droplet_clipping_count));
    set_column(row, "droplet_weighted_deficit",
               real_text(record.diagnostics.droplet_weighted_deficit));
    set_column(
        row, "droplet_roundoff_negative_count",
        std::to_string(record.diagnostics.droplet_roundoff_negative_count));
    set_column(
        row, "droplet_material_negative_count",
        std::to_string(record.diagnostics.droplet_material_negative_count));
    set_column(row, "message", record.message);
    write_row(row);
  }
  for (const auto &summary : convergence) {
    std::vector<std::string> row(headers.size());
    set_column(row, "record_type", "convergence");
    set_column(row, "case_id", "-1");
    set_column(row, "nlev", "12");
    set_column(row, "overlap", "positive_fractional");
    set_column(row, "activation", "all_modes");
    set_column(row, "active_pairs", "25");
    set_column(row, "transported_fields", "51");
    set_column(row, "mixing_index", "1");
    set_column(row, "method", method_name(summary.method));
    set_column(row, "passed", summary.passed ? "true" : "false");
    set_column(row, "deterministic", "true");
    set_column(row, "maximum_weighted_l1_error",
               real_text(summary.finest_weighted_l1));
    set_column(row, "maximum_linf_error",
               real_text(summary.finest_linf));
    set_column(
        row, "weighted_l1_slope_by_field",
        joined_values(51, [&summary](const std::size_t field) {
          return real_text(summary.weighted_l1_fit[field].slope);
        }));
    set_column(
        row, "weighted_l1_r2_by_field",
        joined_values(51, [&summary](const std::size_t field) {
          return real_text(summary.weighted_l1_fit[field].r2);
        }));
    set_column(row, "linf_slope_by_field",
               joined_values(51, [&summary](const std::size_t field) {
                 return real_text(summary.linf_fit[field].slope);
               }));
    set_column(row, "linf_r2_by_field",
               joined_values(51, [&summary](const std::size_t field) {
                 return real_text(summary.linf_fit[field].r2);
               }));
    set_column(row, "message",
               "minimum_slope=" + real_text(summary.minimum_slope) +
                   ";maximum_slope=" + real_text(summary.maximum_slope) +
                   ";minimum_r2=" + real_text(summary.minimum_r2));
    write_row(row);
  }
}

void write_table(const std::string &path, const MatrixResults &matrix,
                 const std::array<ConvergenceSummary, 3> &convergence) {
  std::ostringstream table;
  table << std::setprecision(8);
  table << "# Task 4 common cross-method gate\n\n";
  table << "Maximum reference two-tolerance change: `"
        << matrix.maximum_reference_consistency << "`\n\n";
  table << "| Method | Matrix cases | Failures | Max weighted L1 | Max Linf | "
           "Max residual | Max ledger | Material negatives | Clipping cases | "
           "Slope range | Min R2 | Gate |\n";
  table << "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n";
  for (const auto &fit : convergence) {
    const auto &summary = matrix.summaries[method_index(fit.method)];
    const bool gate = summary.failures == 0 && fit.passed &&
                      matrix.reference_failures == 0;
    table << "| " << method_name(fit.method) << " | " << summary.cases
          << " | " << summary.failures << " | "
          << summary.maximum_weighted_l1_error << " | "
          << summary.maximum_linf_error << " | " << summary.maximum_residual
          << " | " << summary.maximum_ledger_error << " | "
          << summary.material_negative_count << " | "
          << summary.clipping_cases << " | " << fit.minimum_slope << ':'
          << fit.maximum_slope << " | " << fit.minimum_r2 << " | "
          << (gate ? "PASS" : "FAIL") << " |\n";
  }
  if (path.empty()) {
    std::cout << table.str();
    return;
  }
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error("cannot open Markdown output: " + path);
  }
  output << table.str();
}

struct Arguments {
  bool smoke = false;
  std::string csv_path;
  std::string table_path;
};

Arguments parse_arguments(const int argc, char **argv) {
  Arguments arguments;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--smoke") {
      arguments.smoke = true;
    } else if (argument == "--csv" && i + 1 < argc) {
      arguments.csv_path = argv[++i];
    } else if (argument == "--table" && i + 1 < argc) {
      arguments.table_path = argv[++i];
    } else {
      throw std::invalid_argument("unknown or incomplete argument: " +
                                  argument);
    }
  }
  return arguments;
}

} // namespace

int main(int argc, char **argv) {
  const Arguments arguments = parse_arguments(argc, argv);
  Kokkos::initialize(argc, argv);
  try {
    test_reference_against_analytic_solutions();
    test_reference_failure_and_determinism_contract();
    test_error_norm_failure_contract();
    const MatrixResults matrix = run_common_matrix(arguments.smoke);
    const auto convergence = test_common_convergence();
    write_csv(arguments.csv_path, matrix, convergence);
    write_table(arguments.table_path, matrix, convergence);
  } catch (const std::exception &error) {
    ++failed;
    std::cerr << "UNCAUGHT: " << error.what() << '\n';
  }
  Kokkos::finalize();
  std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
