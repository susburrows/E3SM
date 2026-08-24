#ifndef NDROP_INTEGRATOR_HPP
#define NDROP_INTEGRATOR_HPP

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace ndrop_test {

using Real = double;

constexpr int num_modes = 4;
constexpr int num_constituents = 25;

struct State {
  std::vector<Real> droplet_number;
  std::array<std::vector<Real>, num_constituents> interstitial;
  std::array<std::vector<Real>, num_constituents> cloudborne;

  explicit State(std::size_t nlev = 0);
  std::size_t nlev() const;
};

struct OperatorData {
  // pressure_weight[k] = pdel[k] / g [kg m^-2].
  std::vector<Real> pressure_weight;
  // interface_exchange has nlev+1 entries. Endpoints must be exactly zero.
  // Interior entry k is shared by levels k-1 and k [kg m^-2 s^-1].
  std::vector<Real> interface_exchange;
  std::vector<Real> cloud_fraction;
  // Controller-only duplicate of the first active mixing rate [s^-1]. This
  // participates in production substep selection but never in a tendency or
  // physical operator matrix.
  Real first_level_controller_rate = 0;
  std::array<std::vector<Real>, num_modes> number_activation_rate;
  std::array<std::vector<Real>, num_modes> mass_activation_rate;

  explicit OperatorData(std::size_t nlev = 0);
  std::size_t nlev() const;
};

struct PreparedOperator {
  std::vector<Real> pressure_weight;
  std::vector<Real> rate_below;
  std::vector<Real> rate_above;
  std::vector<Real> overlap_below;
  std::vector<Real> overlap_above;
  std::array<std::vector<Real>, num_modes> number_activation_rate;
  std::array<std::vector<Real>, num_modes> mass_activation_rate;
  Real first_level_controller_rate = 0;
  int activation_limiter_count = 0;
  Real max_inverse_mixing_time = 0;

  std::size_t nlev() const;
};

struct FieldClipping {
  Real interstitial_raw_minimum = std::numeric_limits<Real>::infinity();
  Real cloudborne_raw_minimum = std::numeric_limits<Real>::infinity();
  Real interstitial_weighted_deficit = 0;
  Real cloudborne_weighted_deficit = 0;
  int interstitial_count = 0;
  int cloudborne_count = 0;
  int interstitial_roundoff_negative_count = 0;
  int interstitial_material_negative_count = 0;
  int cloudborne_roundoff_negative_count = 0;
  int cloudborne_material_negative_count = 0;
};

struct AdvanceDiagnostics {
  int nsubmix = 0;
  int final_buffer = 0;
  int activation_limiter_count = 0;
  int clipping_count = 0;
  Real raw_minimum = std::numeric_limits<Real>::infinity();
  Real droplet_raw_minimum = std::numeric_limits<Real>::infinity();
  Real droplet_weighted_deficit = 0;
  int droplet_clipping_count = 0;
  Real max_inverse_mixing_time = 0;
  Real first_level_controller_rate = 0;
  Real normalized_residual = 0;
  Real normalized_matrix_residual = 0;
  Real normalized_physical_residual = 0;
  Real max_normalized_ledger_error = 0;
  Real solution_minimum = std::numeric_limits<Real>::infinity();
  int negative_count = 0;
  int roundoff_negative_count = 0;
  int material_negative_count = 0;
  int droplet_roundoff_negative_count = 0;
  int droplet_material_negative_count = 0;
  int block_matrix_family_count = 0;
  int block_factorization_count = 0;
  int block_rhs_solve_count = 0;
  int scalar_factorization_count = 0;
  int scalar_rhs_solve_count = 0;
  bool positivity_failure = false;
  std::array<FieldClipping, num_constituents> field_clipping{};
  std::array<Real, num_constituents> initial_ledger{};
  std::array<Real, num_constituents> final_ledger{};
  std::array<Real, num_constituents> clipping_ledger_delta{};
  std::array<Real, num_constituents> total_ledger_delta{};
  std::array<Real, num_constituents> operator_ledger_delta{};
};

struct AdvanceResult {
  bool ok = false;
  std::string message;
  State state;
  AdvanceDiagnostics diagnostics;
};

struct ForwardEulerOptions {
  // A positive value forces this many equal substeps. Zero uses production
  // adaptive subcycling.
  int forced_substeps = 0;
  bool apply_clipping = true;
};

enum class ImexFiniteProbe {
  solution,
  matrix_residual,
  physical_residual,
  ledger_normalization,
  ledger_decomposition
};

enum class BackwardEulerProbe {
  solution,
  matrix_residual,
  physical_residual,
  ledger_normalization,
  ledger_decomposition,
  block_pivot,
  scalar_pivot,
  material_negative,
  roundoff_negative,
  bottom_branch
};

int constituent_mode(int constituent);
bool constituent_is_number(int constituent);

bool validate(const State &state, const OperatorData &data, Real dt,
              std::string &message);
PreparedOperator prepare_operator(const OperatorData &data);
int production_nsubmix(const PreparedOperator &op, Real dt);

std::vector<Real> aerosol_activation_source(const State &state,
                                            const PreparedOperator &op,
                                            int constituent);
std::vector<Real> droplet_activation_source(const State &state,
                                            const PreparedOperator &op);
std::vector<Real> apply_cloud_operator(const std::vector<Real> &field,
                                       const PreparedOperator &op);
std::vector<Real>
apply_interstitial_operator(const std::vector<Real> &interstitial,
                            const std::vector<Real> &cloudborne,
                            const PreparedOperator &op);

std::array<Real, num_constituents>
column_ledgers(const State &state, const std::vector<Real> &pressure_weight);
bool checked_column_ledgers(
    const State &state, const std::vector<Real> &pressure_weight,
    std::array<Real, num_constituents> &ledgers, std::string &message);
bool finite_prepared_operator(const PreparedOperator &op);
bool finite_values(const std::vector<Real> &values);
bool apply_zero_cloud_postprocess(State &state,
                                  const std::vector<Real> &cloud_fraction,
                                  std::string &message);

AdvanceResult advance_forward_euler(const State &initial,
                                    const OperatorData &data, Real dt,
                                    const ForwardEulerOptions &options = {});

// Implemented in ndrop_imex.cpp and linked only by the EKAT-backed test target.
AdvanceResult advance_imex(const State &initial, const OperatorData &data,
                           Real dt);

// Test-only access to finite-value checkers whose failure branches cannot be
// reached from a valid M-matrix and finite right-hand side.
bool imex_nonfinite_probe_is_rejected(ImexFiniteProbe probe);

// Implemented in ndrop_backward_euler.cpp by the standalone host-only gate.
AdvanceResult advance_backward_euler(const State &initial,
                                     const OperatorData &data, Real dt);

// Test-only access to deterministic failure branches that cannot be reached
// from a valid conservative M-matrix with finite nonnegative input.
bool backward_euler_probe_is_rejected(BackwardEulerProbe probe);

} // namespace ndrop_test

#endif
