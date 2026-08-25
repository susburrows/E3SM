#include "ndrop_integrator.hpp"
#include "ndrop_phase_timer.hpp"

#include <ekat_tridiag.hpp>

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace ndrop_test {
namespace {

using View1D = Kokkos::View<Real *, Kokkos::LayoutRight>;
using View2D = Kokkos::View<Real **, Kokkos::LayoutRight>;

struct Matrix {
  std::vector<Real> lower;
  std::vector<Real> diagonal;
  std::vector<Real> upper;
};

Matrix cloud_matrix(const PreparedOperator &op, const Real dt) {
  Matrix matrix;
  matrix.lower.assign(op.nlev(), 0);
  matrix.diagonal.assign(op.nlev(), 1);
  matrix.upper.assign(op.nlev(), 0);
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    matrix.diagonal[k] += dt * (op.rate_above[k] + op.rate_below[k]);
    if (k > 0) {
      matrix.lower[k] =
          -dt * op.rate_above[k] * op.overlap_above[k];
    }
    if (k + 1 < op.nlev()) {
      matrix.upper[k] =
          -dt * op.rate_below[k] * op.overlap_below[k];
    }
  }
  return matrix;
}

Matrix interstitial_matrix(const PreparedOperator &op, const Real dt) {
  Matrix matrix;
  matrix.lower.assign(op.nlev(), 0);
  matrix.diagonal.assign(op.nlev(), 1);
  matrix.upper.assign(op.nlev(), 0);
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    matrix.diagonal[k] += dt * (op.rate_above[k] + op.rate_below[k]);
    if (k > 0) {
      matrix.lower[k] = -dt * op.rate_above[k];
    }
    if (k + 1 < op.nlev()) {
      matrix.upper[k] = -dt * op.rate_below[k];
    }
  }
  return matrix;
}

bool finite_matrix(const Matrix &matrix) {
  const auto finite = [](const std::vector<Real> &values) {
    return std::all_of(values.begin(), values.end(),
                       [](const Real value) { return std::isfinite(value); });
  };
  return finite(matrix.lower) && finite(matrix.diagonal) &&
         finite(matrix.upper);
}

bool finite_view(const View2D &view) {
  for (int k = 0; k < view.extent_int(0); ++k) {
    for (int j = 0; j < view.extent_int(1); ++j) {
      if (!std::isfinite(view(k, j))) {
        return false;
      }
    }
  }
  return true;
}

void solve(const Matrix &matrix, View2D rhs) {
  const auto nlev = matrix.diagonal.size();
  View1D lower("imex_lower", nlev);
  View1D diagonal("imex_diagonal", nlev);
  View1D upper("imex_upper", nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    lower(k) = matrix.lower[k];
    diagonal(k) = matrix.diagonal[k];
    upper(k) = matrix.upper[k];
  }
  ekat::tridiag::thomas(lower, diagonal, upper, rhs);
}

bool normalized_matrix_residual(const Matrix &matrix, const View2D &solution,
                                const View2D &rhs, Real &normalized_max) {
  Real overall_max = 0;
  const int nlev = solution.extent_int(0);
  const int nrhs = solution.extent_int(1);
  for (int j = 0; j < nrhs; ++j) {
    Real residual_max = 0;
    Real scale_max = 0;
    for (int k = 0; k < nlev; ++k) {
      Real value_scale =
          std::max(std::abs(solution(k, j)), std::abs(rhs(k, j)));
      Real coefficient_scale =
          std::max(Real(1), std::abs(matrix.diagonal[k]));
      if (k > 0) {
        value_scale =
            std::max(value_scale, std::abs(solution(k - 1, j)));
        coefficient_scale =
            std::max(coefficient_scale, std::abs(matrix.lower[k]));
      }
      if (k + 1 < nlev) {
        value_scale =
            std::max(value_scale, std::abs(solution(k + 1, j)));
        coefficient_scale =
            std::max(coefficient_scale, std::abs(matrix.upper[k]));
      }
      if (value_scale == 0) {
        continue;
      }
      Real applied =
          (matrix.diagonal[k] / coefficient_scale) *
          (solution(k, j) / value_scale);
      Real scale = std::abs(applied);
      if (k > 0) {
        const Real term =
            (matrix.lower[k] / coefficient_scale) *
            (solution(k - 1, j) / value_scale);
        applied += term;
        scale += std::abs(term);
      }
      if (k + 1 < nlev) {
        const Real term =
            (matrix.upper[k] / coefficient_scale) *
            (solution(k + 1, j) / value_scale);
        applied += term;
        scale += std::abs(term);
      }
      const Real scaled_rhs =
          (rhs(k, j) / value_scale) / coefficient_scale;
      const Real residual = std::abs(applied - scaled_rhs);
      scale += std::abs(scaled_rhs);
      if (!std::isfinite(applied) || !std::isfinite(residual) ||
          !std::isfinite(scale)) {
        return false;
      }
      residual_max = std::max(residual_max, residual);
      scale_max = std::max(scale_max, scale);
    }
    const Real normalized =
        scale_max > 0 ? residual_max / scale_max : residual_max;
    if (!std::isfinite(normalized)) {
      return false;
    }
    overall_max = std::max(overall_max, normalized);
  }
  normalized_max = overall_max;
  return std::isfinite(normalized_max);
}

bool update_physical_residual(const std::vector<Real> &old,
                              const std::vector<Real> &solution,
                              const std::vector<Real> &source,
                              const PreparedOperator &op, const Real dt,
                              const bool interstitial,
                              const std::vector<Real> *cloudborne,
                              Real &overall_max) {
  Real residual_max = 0;
  Real scale_max = 0;
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    const std::size_t km1 = k == 0 ? 0 : k - 1;
    const std::size_t kp1 = std::min(k + 1, op.nlev() - 1);
    const Real dt_rate_above = dt * op.rate_above[k];
    const Real dt_rate_below = dt * op.rate_below[k];
    if (!std::isfinite(dt_rate_above) || !std::isfinite(dt_rate_below)) {
      return false;
    }
    Real value_scale =
        std::max({std::abs(old[k]), std::abs(solution[k]),
                  std::abs(solution[km1]), std::abs(solution[kp1])});
    if (interstitial) {
      value_scale =
          std::max({value_scale, std::abs((*cloudborne)[km1]),
                    std::abs((*cloudborne)[kp1])});
    }
    if (value_scale == 0) {
      continue;
    }
    const Real coefficient_scale =
        std::max({Real(1), std::abs(dt_rate_above),
                  std::abs(dt_rate_below)});
    const Real scaled_current = solution[k] / value_scale;
    Real residual =
        scaled_current / coefficient_scale -
        (old[k] / value_scale) / coefficient_scale;
    Real scale = std::abs(scaled_current / coefficient_scale) +
                 std::abs((old[k] / value_scale) / coefficient_scale);
    const Real source_term =
        (dt / coefficient_scale) * (source[k] / value_scale) *
        (interstitial ? -1 : 1);
    residual -= source_term;
    scale += std::abs(source_term);

    Real below_bracket =
        solution[kp1] / value_scale - scaled_current;
    Real above_bracket =
        solution[km1] / value_scale - scaled_current;
    if (interstitial) {
      below_bracket += ((*cloudborne)[kp1] / value_scale) *
                       (1 - op.overlap_below[k]);
      above_bracket += ((*cloudborne)[km1] / value_scale) *
                       (1 - op.overlap_above[k]);
    } else {
      below_bracket =
          op.overlap_below[k] * (solution[kp1] / value_scale) -
          scaled_current;
      above_bracket =
          op.overlap_above[k] * (solution[km1] / value_scale) -
          scaled_current;
    }
    const Real below_term =
        (dt_rate_below / coefficient_scale) * below_bracket;
    const Real above_term =
        (dt_rate_above / coefficient_scale) * above_bracket;
    residual -= below_term + above_term;
    scale += std::abs(below_term) + std::abs(above_term);
    residual = std::abs(residual);
    if (!std::isfinite(residual) || !std::isfinite(scale)) {
      return false;
    }
    residual_max = std::max(residual_max, residual);
    scale_max = std::max(scale_max, scale);
  }
  const Real normalized =
      scale_max > 0 ? residual_max / scale_max : residual_max;
  if (!std::isfinite(normalized)) {
    return false;
  }
  overall_max = std::max(overall_max, normalized);
  return std::isfinite(overall_max);
}

bool normalized_ledger_delta(const Real initial, const Real delta,
                             Real &normalized) {
  const Real scale = std::abs(initial);
  normalized = scale > 0 ? std::abs(delta) / scale : std::abs(delta);
  return std::isfinite(normalized);
}

bool checked_difference(const Real lhs, const Real rhs, Real &difference) {
  const Real maximum = std::numeric_limits<Real>::max();
  if ((rhs > 0 && lhs < -maximum + rhs) ||
      (rhs < 0 && lhs > maximum + rhs)) {
    return false;
  }
  difference = lhs - rhs;
  return std::isfinite(difference);
}

void update_positivity(const std::vector<Real> &initial,
                       const std::vector<Real> &solution,
                       AdvanceDiagnostics &diagnostics,
                       int &roundoff_count, int &material_count) {
  Real scale = 0;
  for (const Real value : initial) {
    scale = std::max(scale, std::abs(value));
  }
  for (const Real value : solution) {
    scale = std::max(scale, std::abs(value));
    diagnostics.solution_minimum =
        std::min(diagnostics.solution_minimum, value);
  }
  const Real threshold =
      10 * std::numeric_limits<Real>::epsilon() * scale;
  for (const Real value : solution) {
    if (value < 0) {
      ++diagnostics.negative_count;
      if (-value > threshold) {
        ++diagnostics.material_negative_count;
        ++material_count;
        diagnostics.positivity_failure = true;
      } else {
        ++diagnostics.roundoff_negative_count;
        ++roundoff_count;
      }
    }
  }
}

} // namespace

AdvanceResult advance_imex(const State &initial, const OperatorData &data,
                           const Real dt, const AdvanceControl &control) {
  PhaseStopwatch stopwatch(control.phase_timings);
  stopwatch.phase(&AdvancePhaseTimings::setup);
  AdvanceResult result;
  result.state = State(initial.nlev());
  if (!validate(initial, data, dt, result.message)) {
    return result;
  }
  const PreparedOperator op = prepare_operator(data);
  if (!finite_prepared_operator(op)) {
    result.message = "derived mixing rates must be finite";
    return result;
  }

  result.diagnostics.nsubmix = 1;
  result.diagnostics.final_buffer = 0;
  result.diagnostics.activation_limiter_count =
      op.activation_limiter_count;
  result.diagnostics.max_inverse_mixing_time =
      op.max_inverse_mixing_time;
  result.diagnostics.first_level_controller_rate =
      op.first_level_controller_rate;
  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    if (!checked_column_ledgers(initial, op.pressure_weight,
                                result.diagnostics.initial_ledger,
                                result.message)) {
      return result;
    }
    stopwatch.phase(&AdvancePhaseTimings::setup);
  }

  const Matrix cloud = cloud_matrix(op, dt);
  const Matrix interstitial = interstitial_matrix(op, dt);
  if (!finite_matrix(cloud) || !finite_matrix(interstitial)) {
    result.message = "IMEX matrix contains a nonfinite coefficient";
    return result;
  }
  const auto droplet_source = droplet_activation_source(initial, op);
  if (!finite_values(droplet_source)) {
    result.message = "droplet activation source contains a nonfinite value";
    return result;
  }
  std::array<std::vector<Real>, num_constituents> sources;
  for (int j = 0; j < num_constituents; ++j) {
    sources[j] = aerosol_activation_source(initial, op, j);
    if (!finite_values(sources[j])) {
      result.message = "aerosol activation source contains a nonfinite value";
      return result;
    }
  }

  constexpr int cloud_rhs_count = num_constituents + 1;
  View2D cloud_rhs("imex_cloud_rhs", op.nlev(), cloud_rhs_count);
  View2D cloud_rhs_original;
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    cloud_rhs(k, 0) = initial.droplet_number[k] + dt * droplet_source[k];
    for (int j = 0; j < num_constituents; ++j) {
      cloud_rhs(k, j + 1) =
          initial.cloudborne[j][k] + dt * sources[j][k];
    }
  }
  if (!finite_view(cloud_rhs)) {
    result.message = "IMEX cloud right-hand side contains a nonfinite value";
    return result;
  }
  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    cloud_rhs_original = View2D("imex_cloud_rhs_original", op.nlev(),
                                cloud_rhs_count);
    Kokkos::deep_copy(cloud_rhs_original, cloud_rhs);
    stopwatch.phase(&AdvancePhaseTimings::setup);
  }
  stopwatch.phase(&AdvancePhaseTimings::combined_factorization_solve);
  solve(cloud, cloud_rhs);
  stopwatch.phase(&AdvancePhaseTimings::update);
  if (!finite_view(cloud_rhs)) {
    result.message = "IMEX cloud solve produced a nonfinite value";
    return result;
  }
  Real matrix_residual = 0;
  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    if (!normalized_matrix_residual(cloud, cloud_rhs, cloud_rhs_original,
                                    matrix_residual)) {
      result.message = "IMEX cloud residual is nonfinite";
      return result;
    }
    result.diagnostics.normalized_matrix_residual = matrix_residual;
    stopwatch.phase(&AdvancePhaseTimings::update);
  }

  for (std::size_t k = 0; k < op.nlev(); ++k) {
    result.state.droplet_number[k] = cloud_rhs(k, 0);
    for (int j = 0; j < num_constituents; ++j) {
      result.state.cloudborne[j][k] = cloud_rhs(k, j + 1);
    }
  }

  View2D interstitial_rhs("imex_interstitial_rhs", op.nlev(),
                          num_constituents);
  View2D interstitial_rhs_original;
  for (int j = 0; j < num_constituents; ++j) {
    const auto return_and_mixing = apply_interstitial_operator(
        std::vector<Real>(op.nlev(), 0), result.state.cloudborne[j], op);
    for (std::size_t k = 0; k < op.nlev(); ++k) {
      interstitial_rhs(k, j) =
          initial.interstitial[j][k] - dt * sources[j][k] +
          dt * return_and_mixing[k];
    }
  }
  if (!finite_view(interstitial_rhs)) {
    result.message =
        "IMEX interstitial right-hand side contains a nonfinite value";
    return result;
  }
  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    interstitial_rhs_original =
        View2D("imex_interstitial_rhs_original", op.nlev(),
               num_constituents);
    Kokkos::deep_copy(interstitial_rhs_original, interstitial_rhs);
    stopwatch.phase(&AdvancePhaseTimings::setup);
  }
  stopwatch.phase(&AdvancePhaseTimings::combined_factorization_solve);
  solve(interstitial, interstitial_rhs);
  stopwatch.phase(&AdvancePhaseTimings::update);
  if (!finite_view(interstitial_rhs)) {
    result.message = "IMEX interstitial solve produced a nonfinite value";
    return result;
  }
  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    matrix_residual = 0;
    if (!normalized_matrix_residual(interstitial, interstitial_rhs,
                                    interstitial_rhs_original,
                                    matrix_residual)) {
      result.message = "IMEX interstitial residual is nonfinite";
      return result;
    }
    result.diagnostics.normalized_matrix_residual =
        std::max(result.diagnostics.normalized_matrix_residual,
                 matrix_residual);
    stopwatch.phase(&AdvancePhaseTimings::update);
  }
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    for (int j = 0; j < num_constituents; ++j) {
      result.state.interstitial[j][k] = interstitial_rhs(k, j);
    }
  }

  update_positivity(initial.droplet_number, result.state.droplet_number,
                    result.diagnostics,
                    result.diagnostics.droplet_roundoff_negative_count,
                    result.diagnostics.droplet_material_negative_count);
  for (int j = 0; j < num_constituents; ++j) {
    update_positivity(initial.interstitial[j],
                      result.state.interstitial[j], result.diagnostics,
                      result.diagnostics.field_clipping[j]
                          .interstitial_roundoff_negative_count,
                      result.diagnostics.field_clipping[j]
                          .interstitial_material_negative_count);
    update_positivity(initial.cloudborne[j], result.state.cloudborne[j],
                      result.diagnostics,
                      result.diagnostics.field_clipping[j]
                          .cloudborne_roundoff_negative_count,
                      result.diagnostics.field_clipping[j]
                          .cloudborne_material_negative_count);
  }

  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    Real physical_residual = 0;
    if (!update_physical_residual(
            initial.droplet_number, result.state.droplet_number,
            droplet_source, op, dt, false, nullptr, physical_residual)) {
      result.message = "IMEX physical residual is nonfinite";
      return result;
    }
    for (int j = 0; j < num_constituents; ++j) {
      if (!update_physical_residual(
              initial.cloudborne[j], result.state.cloudborne[j], sources[j],
              op, dt, false, nullptr, physical_residual) ||
          !update_physical_residual(
              initial.interstitial[j], result.state.interstitial[j],
              sources[j], op, dt, true, &result.state.cloudborne[j],
              physical_residual)) {
        result.message = "IMEX physical residual is nonfinite";
        return result;
      }
    }
    result.diagnostics.normalized_physical_residual = physical_residual;
    result.diagnostics.normalized_residual =
        std::max(result.diagnostics.normalized_matrix_residual,
                 result.diagnostics.normalized_physical_residual);
    if (!std::isfinite(result.diagnostics.normalized_residual)) {
      result.message = "IMEX normalized residual is nonfinite";
      return result;
    }
    stopwatch.phase(&AdvancePhaseTimings::update);
  }
  if (result.diagnostics.positivity_failure) {
    result.message = "material negative value in IMEX solution";
    return result;
  }
  if (!apply_zero_cloud_postprocess(result.state, data.cloud_fraction,
                                    result.message)) {
    return result;
  }
  if (!control.collect_detailed_diagnostics) {
    result.ok = true;
    result.message = "ok";
    return result;
  }

  stopwatch.phase(&AdvancePhaseTimings::diagnostics);
  if (!checked_column_ledgers(result.state, op.pressure_weight,
                              result.diagnostics.final_ledger,
                              result.message)) {
    return result;
  }
  for (int j = 0; j < num_constituents; ++j) {
    Real total_delta = 0;
    Real operator_delta = 0;
    if (!checked_difference(result.diagnostics.final_ledger[j],
                            result.diagnostics.initial_ledger[j],
                            total_delta) ||
        !checked_difference(total_delta,
                            result.diagnostics.clipping_ledger_delta[j],
                            operator_delta)) {
      result.message = "IMEX aerosol ledger decomposition is nonfinite";
      return result;
    }
    result.diagnostics.total_ledger_delta[j] = total_delta;
    result.diagnostics.operator_ledger_delta[j] = operator_delta;
    Real normalized = 0;
    if (!normalized_ledger_delta(result.diagnostics.initial_ledger[j],
                                 operator_delta, normalized)) {
      result.message = "IMEX aerosol ledger normalization is nonfinite";
      return result;
    }
    result.diagnostics.max_normalized_ledger_error =
        std::max(result.diagnostics.max_normalized_ledger_error,
                 normalized);
  }
  const Real ledger_tolerance =
      std::max(Real(1e-11),
               5 * static_cast<Real>(op.nlev()) *
                   std::numeric_limits<Real>::epsilon());
  result.ok = !result.diagnostics.positivity_failure &&
              result.diagnostics.normalized_residual <= 1e-11 &&
              result.diagnostics.max_normalized_ledger_error <=
                  ledger_tolerance;
  if (result.diagnostics.normalized_residual > 1e-11) {
    result.message = "IMEX normalized residual exceeds tolerance";
  } else if (result.diagnostics.max_normalized_ledger_error >
             ledger_tolerance) {
    result.message = "IMEX aerosol ledger error exceeds tolerance";
  } else {
    result.message = "ok";
  }
  return result;
}

bool imex_nonfinite_probe_is_rejected(const ImexFiniteProbe probe) {
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  switch (probe) {
  case ImexFiniteProbe::solution: {
    View2D solution("imex_nonfinite_solution_probe", 1, 1);
    Kokkos::deep_copy(solution, nan);
    return !finite_view(solution);
  }
  case ImexFiniteProbe::matrix_residual: {
    const Matrix matrix{{0}, {1}, {0}};
    View2D solution("imex_nonfinite_matrix_solution_probe", 1, 1);
    View2D rhs("imex_nonfinite_matrix_rhs_probe", 1, 1);
    Kokkos::deep_copy(solution, nan);
    Kokkos::deep_copy(rhs, Real(1));
    Real normalized = 0;
    return !normalized_matrix_residual(matrix, solution, rhs, normalized);
  }
  case ImexFiniteProbe::physical_residual: {
    PreparedOperator op;
    op.pressure_weight = {1};
    op.rate_below = {0};
    op.rate_above = {0};
    op.overlap_below = {1};
    op.overlap_above = {1};
    Real normalized = 0;
    return !update_physical_residual({1}, {1}, {nan}, op, 1, false,
                                     nullptr, normalized);
  }
  case ImexFiniteProbe::ledger_normalization: {
    Real normalized = 0;
    return !normalized_ledger_delta(std::numeric_limits<Real>::denorm_min(),
                                    std::numeric_limits<Real>::max(),
                                    normalized);
  }
  case ImexFiniteProbe::ledger_decomposition: {
    Real difference = 0;
    return !checked_difference(std::numeric_limits<Real>::max(),
                               -std::numeric_limits<Real>::max(),
                               difference);
  }
  }
  return false;
}

} // namespace ndrop_test
