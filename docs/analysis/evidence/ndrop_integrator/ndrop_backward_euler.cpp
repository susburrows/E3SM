#include "ndrop_integrator.hpp"
#include "ndrop_phase_timer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace ndrop_test {
namespace {

struct Vec2 {
  Real x = 0;
  Real y = 0;
};

struct Mat2 {
  Real a00 = 0;
  Real a01 = 0;
  Real a10 = 0;
  Real a11 = 0;
};

struct BlockSystem {
  std::vector<Mat2> lower;
  std::vector<Mat2> diagonal;
  std::vector<Mat2> upper;
};

struct BlockFactorization {
  std::vector<Mat2> inverse_diagonal;
  std::vector<Mat2> reduced_upper;
};

struct ScalarSystem {
  std::vector<Real> lower;
  std::vector<Real> diagonal;
  std::vector<Real> upper;
};

struct ScalarFactorization {
  std::vector<Real> inverse_diagonal;
  std::vector<Real> reduced_upper;
};

struct MatrixFamily {
  std::vector<Real> activation;
  std::vector<int> constituents;
  BlockSystem system;
  BlockFactorization factorization;
};

bool finite(const Real value) { return std::isfinite(value); }

bool finite(const Vec2 &value) { return finite(value.x) && finite(value.y); }

bool finite(const Mat2 &matrix) {
  return finite(matrix.a00) && finite(matrix.a01) && finite(matrix.a10) &&
         finite(matrix.a11);
}

Vec2 subtract(const Vec2 &lhs, const Vec2 &rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y};
}

Mat2 subtract(const Mat2 &lhs, const Mat2 &rhs) {
  return {lhs.a00 - rhs.a00, lhs.a01 - rhs.a01,
          lhs.a10 - rhs.a10, lhs.a11 - rhs.a11};
}

Vec2 multiply(const Mat2 &matrix, const Vec2 &vector) {
  return {matrix.a00 * vector.x + matrix.a01 * vector.y,
          matrix.a10 * vector.x + matrix.a11 * vector.y};
}

Mat2 multiply(const Mat2 &lhs, const Mat2 &rhs) {
  return {lhs.a00 * rhs.a00 + lhs.a01 * rhs.a10,
          lhs.a00 * rhs.a01 + lhs.a01 * rhs.a11,
          lhs.a10 * rhs.a00 + lhs.a11 * rhs.a10,
          lhs.a10 * rhs.a01 + lhs.a11 * rhs.a11};
}

Real block_system_norm(const BlockSystem &system) {
  Real norm = 0;
  for (std::size_t k = 0; k < system.diagonal.size(); ++k) {
    const Real row0 =
        std::abs(system.lower[k].a00) + std::abs(system.lower[k].a01) +
        std::abs(system.diagonal[k].a00) +
        std::abs(system.diagonal[k].a01) +
        std::abs(system.upper[k].a00) + std::abs(system.upper[k].a01);
    const Real row1 =
        std::abs(system.lower[k].a10) + std::abs(system.lower[k].a11) +
        std::abs(system.diagonal[k].a10) +
        std::abs(system.diagonal[k].a11) +
        std::abs(system.upper[k].a10) + std::abs(system.upper[k].a11);
    if (!finite(row0) || !finite(row1)) {
      return std::numeric_limits<Real>::infinity();
    }
    norm = std::max({norm, row0, row1});
  }
  return norm;
}

bool invert_block(const Mat2 &matrix, const Real threshold, Mat2 &inverse,
                  ReplayDiagnosticReason *reason = nullptr) {
  Real pivot00 = matrix.a00;
  Real pivot01 = matrix.a01;
  Real lower10 = matrix.a10;
  Real lower11 = matrix.a11;
  if (std::abs(lower10) > std::abs(pivot00)) {
    std::swap(pivot00, lower10);
    std::swap(pivot01, lower11);
  }
  if (!finite(pivot00) || std::abs(pivot00) <= threshold) {
    if (reason != nullptr)
      *reason = finite(pivot00)
                    ? ReplayDiagnosticReason::factorization_unsafe_pivot
                    : ReplayDiagnosticReason::factorization_nonfinite;
    return false;
  }
  const Real multiplier = lower10 / pivot00;
  const Real pivot11 = lower11 - multiplier * pivot01;
  if (!finite(multiplier) || !finite(pivot11) ||
      std::abs(pivot11) <= threshold) {
    if (reason != nullptr)
      *reason = finite(multiplier) && finite(pivot11)
                    ? ReplayDiagnosticReason::factorization_unsafe_pivot
                    : ReplayDiagnosticReason::factorization_nonfinite;
    return false;
  }
  const Real determinant =
      matrix.a00 * matrix.a11 - matrix.a01 * matrix.a10;
  if (!finite(determinant) || determinant == 0) {
    if (reason != nullptr)
      *reason = finite(determinant)
                    ? ReplayDiagnosticReason::factorization_unsafe_pivot
                    : ReplayDiagnosticReason::factorization_nonfinite;
    return false;
  }
  inverse = {matrix.a11 / determinant, -matrix.a01 / determinant,
             -matrix.a10 / determinant, matrix.a00 / determinant};
  const bool ok = finite(inverse);
  if (!ok && reason != nullptr)
    *reason = ReplayDiagnosticReason::factorization_nonfinite;
  return ok;
}

bool finite_block_system(const BlockSystem &system) {
  if (system.lower.size() != system.diagonal.size() ||
      system.upper.size() != system.diagonal.size()) {
    return false;
  }
  for (std::size_t k = 0; k < system.diagonal.size(); ++k) {
    if (!finite(system.lower[k]) || !finite(system.diagonal[k]) ||
        !finite(system.upper[k])) {
      return false;
    }
  }
  return true;
}

bool factor_block_system(const BlockSystem &system,
                         BlockFactorization &factorization,
                         ReplayDiagnosticReason *reason = nullptr) {
  const std::size_t nlev = system.diagonal.size();
  const Real norm = block_system_norm(system);
  if (!finite(norm)) {
    if (reason != nullptr)
      *reason = ReplayDiagnosticReason::factorization_nonfinite;
    return false;
  }
  const Real threshold =
      100 * std::numeric_limits<Real>::epsilon() * norm;
  factorization.inverse_diagonal.assign(nlev, {});
  factorization.reduced_upper.assign(nlev, {});
  for (std::size_t k = 0; k < nlev; ++k) {
    Mat2 reduced = system.diagonal[k];
    if (k > 0) {
      reduced = subtract(
          reduced, multiply(system.lower[k],
                            factorization.reduced_upper[k - 1]));
    }
    if (!finite(reduced)) {
      if (reason != nullptr)
        *reason = ReplayDiagnosticReason::factorization_nonfinite;
      return false;
    }
    if (!invert_block(reduced, threshold,
                      factorization.inverse_diagonal[k], reason)) {
      return false;
    }
    factorization.reduced_upper[k] =
        multiply(factorization.inverse_diagonal[k], system.upper[k]);
    if (!finite(factorization.reduced_upper[k])) {
      if (reason != nullptr)
        *reason = ReplayDiagnosticReason::factorization_nonfinite;
      return false;
    }
  }
  return true;
}

bool solve_block_system(const BlockSystem &system,
                        const BlockFactorization &factorization,
                        const std::vector<Vec2> &rhs,
                        std::vector<Vec2> &solution) {
  const std::size_t nlev = rhs.size();
  std::vector<Vec2> reduced_rhs(nlev);
  for (std::size_t k = 0; k < nlev; ++k) {
    Vec2 value = rhs[k];
    if (k > 0) {
      value = subtract(value, multiply(system.lower[k],
                                       reduced_rhs[k - 1]));
    }
    reduced_rhs[k] = multiply(factorization.inverse_diagonal[k], value);
    if (!finite(reduced_rhs[k])) {
      return false;
    }
  }
  solution.assign(nlev, {});
  for (std::size_t reverse = 0; reverse < nlev; ++reverse) {
    const std::size_t k = nlev - 1 - reverse;
    solution[k] = reduced_rhs[k];
    if (k + 1 < nlev) {
      solution[k] = subtract(
          solution[k],
          multiply(factorization.reduced_upper[k], solution[k + 1]));
    }
    if (!finite(solution[k])) {
      return false;
    }
  }
  return true;
}

BlockSystem assemble_block_system(const PreparedOperator &op,
                                  const std::vector<Real> &activation,
                                  const Real dt) {
  BlockSystem system;
  system.lower.assign(op.nlev(), {});
  system.diagonal.assign(op.nlev(), {});
  system.upper.assign(op.nlev(), {});
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    const Real up = dt * op.rate_above[k];
    const Real down = dt * op.rate_below[k];
    const Real source = dt * activation[k];
    const Real diagonal = 1 + up + down;
    if (k > 0) {
      system.lower[k] =
          {-up * op.overlap_above[k], 0,
           -up * (1 - op.overlap_above[k]), -up};
    }
    system.diagonal[k] = {diagonal, 0, 0, diagonal};
    if (k + 1 < op.nlev()) {
      system.upper[k] =
          {-down * op.overlap_below[k], -source,
           -down * (1 - op.overlap_below[k]), -down + source};
    } else {
      system.diagonal[k] = {diagonal - source, -source, source,
                            diagonal + source};
    }
  }
  return system;
}

std::vector<Real> activation_coefficients(const PreparedOperator &op,
                                          const int constituent) {
  const int mode = constituent_mode(constituent);
  std::vector<Real> coefficients =
      constituent_is_number(constituent)
          ? op.number_activation_rate[mode]
          : op.mass_activation_rate[mode];
  coefficients.back() = op.number_activation_rate[mode].back();
  return coefficients;
}

bool same_coefficients(const std::vector<Real> &lhs,
                       const std::vector<Real> &rhs) {
  return lhs == rhs;
}

bool normalized_block_residual(const BlockSystem &system,
                               const std::vector<Vec2> &solution,
                               const std::vector<Vec2> &rhs,
                               Real &overall_max,
                               ReplayDiagnosticDetail *detail = nullptr,
                               const int constituent = 0,
                               const std::uint32_t substep = 0) {
  for (std::size_t k = 0; k < solution.size(); ++k) {
    Vec2 applied = multiply(system.diagonal[k], solution[k]);
    Real scale0 =
        std::abs(system.diagonal[k].a00 * solution[k].x) +
        std::abs(system.diagonal[k].a01 * solution[k].y) +
        std::abs(rhs[k].x);
    Real scale1 =
        std::abs(system.diagonal[k].a10 * solution[k].x) +
        std::abs(system.diagonal[k].a11 * solution[k].y) +
        std::abs(rhs[k].y);
    if (k > 0) {
      const Vec2 term = multiply(system.lower[k], solution[k - 1]);
      applied.x += term.x;
      applied.y += term.y;
      scale0 += std::abs(system.lower[k].a00 * solution[k - 1].x) +
                std::abs(system.lower[k].a01 * solution[k - 1].y);
      scale1 += std::abs(system.lower[k].a10 * solution[k - 1].x) +
                std::abs(system.lower[k].a11 * solution[k - 1].y);
    }
    if (k + 1 < solution.size()) {
      const Vec2 term = multiply(system.upper[k], solution[k + 1]);
      applied.x += term.x;
      applied.y += term.y;
      scale0 += std::abs(system.upper[k].a00 * solution[k + 1].x) +
                std::abs(system.upper[k].a01 * solution[k + 1].y);
      scale1 += std::abs(system.upper[k].a10 * solution[k + 1].x) +
                std::abs(system.upper[k].a11 * solution[k + 1].y);
    }
    const Real residual0 = std::abs(applied.x - rhs[k].x);
    const Real residual1 = std::abs(applied.y - rhs[k].y);
    if (!finite(applied) || !finite(scale0) || !finite(scale1) ||
        !finite(residual0) || !finite(residual1)) {
      if (detail != nullptr) {
        detail->matrix_residual[26 + constituent].reset();
        detail->matrix_residual[1 + constituent].reset();
      }
      return false;
    }
    const Real ratio0 = scale0 > 0 ? residual0 / scale0 : residual0;
    const Real ratio1 = scale1 > 0 ? residual1 / scale1 : residual1;
    overall_max = std::max({overall_max, ratio0, ratio1});
    if (detail != nullptr) {
      auto &cloud = detail->matrix_residual[26 + constituent];
      auto &interstitial = detail->matrix_residual[1 + constituent];
      if (!cloud || ratio0 > cloud->normalized) {
        cloud = ReplayResidualMaximum{residual0, scale0, ratio0, substep};
      }
      if (!interstitial || ratio1 > interstitial->normalized) {
        interstitial = ReplayResidualMaximum{residual1, scale1, ratio1,
                                             substep};
      }
    }
  }
  return finite(overall_max);
}

ScalarSystem assemble_cloud_system(const PreparedOperator &op,
                                   const Real dt) {
  ScalarSystem system;
  system.lower.assign(op.nlev(), 0);
  system.diagonal.assign(op.nlev(), 1);
  system.upper.assign(op.nlev(), 0);
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    system.diagonal[k] +=
        dt * (op.rate_above[k] + op.rate_below[k]);
    if (k > 0) {
      system.lower[k] =
          -dt * op.rate_above[k] * op.overlap_above[k];
    }
    if (k + 1 < op.nlev()) {
      system.upper[k] =
          -dt * op.rate_below[k] * op.overlap_below[k];
    }
  }
  return system;
}

bool finite_scalar_system(const ScalarSystem &system) {
  return finite_values(system.lower) && finite_values(system.diagonal) &&
         finite_values(system.upper);
}

Real scalar_system_norm(const ScalarSystem &system) {
  Real norm = 0;
  for (std::size_t k = 0; k < system.diagonal.size(); ++k) {
    const Real row = std::abs(system.lower[k]) +
                     std::abs(system.diagonal[k]) +
                     std::abs(system.upper[k]);
    if (!finite(row)) {
      return std::numeric_limits<Real>::infinity();
    }
    norm = std::max(norm, row);
  }
  return norm;
}

bool factor_scalar_system(const ScalarSystem &system,
                          ScalarFactorization &factorization,
                          ReplayDiagnosticReason *reason = nullptr) {
  const Real norm = scalar_system_norm(system);
  if (!finite(norm)) {
    if (reason != nullptr)
      *reason = ReplayDiagnosticReason::factorization_nonfinite;
    return false;
  }
  const Real threshold =
      100 * std::numeric_limits<Real>::epsilon() * norm;
  const std::size_t nlev = system.diagonal.size();
  factorization.inverse_diagonal.assign(nlev, 0);
  factorization.reduced_upper.assign(nlev, 0);
  for (std::size_t k = 0; k < nlev; ++k) {
    Real pivot = system.diagonal[k];
    if (k > 0) {
      pivot -= system.lower[k] * factorization.reduced_upper[k - 1];
    }
    if (!finite(pivot) || std::abs(pivot) <= threshold) {
      if (reason != nullptr)
        *reason = finite(pivot)
                      ? ReplayDiagnosticReason::factorization_unsafe_pivot
                      : ReplayDiagnosticReason::factorization_nonfinite;
      return false;
    }
    factorization.inverse_diagonal[k] = 1 / pivot;
    factorization.reduced_upper[k] =
        factorization.inverse_diagonal[k] * system.upper[k];
    if (!finite(factorization.inverse_diagonal[k]) ||
        !finite(factorization.reduced_upper[k])) {
      if (reason != nullptr)
        *reason = ReplayDiagnosticReason::factorization_nonfinite;
      return false;
    }
  }
  return true;
}

bool solve_scalar_system(const ScalarSystem &system,
                         const ScalarFactorization &factorization,
                         const std::vector<Real> &rhs,
                         std::vector<Real> &solution) {
  const std::size_t nlev = rhs.size();
  std::vector<Real> reduced_rhs(nlev, 0);
  for (std::size_t k = 0; k < nlev; ++k) {
    Real value = rhs[k];
    if (k > 0) {
      value -= system.lower[k] * reduced_rhs[k - 1];
    }
    reduced_rhs[k] = factorization.inverse_diagonal[k] * value;
    if (!finite(reduced_rhs[k])) {
      return false;
    }
  }
  solution.assign(nlev, 0);
  for (std::size_t reverse = 0; reverse < nlev; ++reverse) {
    const std::size_t k = nlev - 1 - reverse;
    solution[k] = reduced_rhs[k];
    if (k + 1 < nlev) {
      solution[k] -= factorization.reduced_upper[k] * solution[k + 1];
    }
    if (!finite(solution[k])) {
      return false;
    }
  }
  return true;
}

bool normalized_scalar_residual(const ScalarSystem &system,
                                const std::vector<Real> &solution,
                                const std::vector<Real> &rhs,
                                Real &overall_max,
                                ReplayDiagnosticDetail *detail = nullptr,
                                const std::uint32_t substep = 0) {
  for (std::size_t k = 0; k < solution.size(); ++k) {
    Real applied = system.diagonal[k] * solution[k];
    Real scale = std::abs(applied) + std::abs(rhs[k]);
    if (k > 0) {
      const Real term = system.lower[k] * solution[k - 1];
      applied += term;
      scale += std::abs(term);
    }
    if (k + 1 < solution.size()) {
      const Real term = system.upper[k] * solution[k + 1];
      applied += term;
      scale += std::abs(term);
    }
    const Real residual = std::abs(applied - rhs[k]);
    if (!finite(applied) || !finite(scale) || !finite(residual)) {
      if (detail != nullptr) detail->matrix_residual[0].reset();
      return false;
    }
    const Real ratio = scale > 0 ? residual / scale : residual;
    overall_max = std::max(overall_max, ratio);
    if (detail != nullptr) {
      auto &maximum = detail->matrix_residual[0];
      if (!maximum || ratio > maximum->normalized) {
        maximum = ReplayResidualMaximum{residual, scale, ratio, substep};
      }
    }
  }
  return finite(overall_max);
}

void update_positivity(const std::vector<Real> &initial,
                       const std::vector<Real> &solution,
                       AdvanceDiagnostics &diagnostics,
                       int &roundoff_count, int &material_count,
                       Real &field_minimum) {
  Real scale = 0;
  for (const Real value : initial) {
    scale = std::max(scale, std::abs(value));
  }
  for (const Real value : solution) {
    scale = std::max(scale, std::abs(value));
    field_minimum = std::min(field_minimum, value);
    diagnostics.raw_minimum = std::min(diagnostics.raw_minimum, value);
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

bool bottom_branch_status(const std::vector<Real> &initial_interstitial,
                          const std::vector<Real> &initial_cloudborne,
                          const std::vector<Real> &solution_interstitial,
                          const std::vector<Real> &solution_cloudborne,
                          bool &changed) {
  const std::size_t bottom = initial_interstitial.size() - 1;
  const Real initial_sum =
      initial_interstitial[bottom] + initial_cloudborne[bottom];
  const Real solution_sum =
      solution_interstitial[bottom] + solution_cloudborne[bottom];
  if (!finite(initial_sum) || !finite(solution_sum)) {
    return false;
  }
  const Real scale =
      std::max({std::abs(initial_sum), std::abs(solution_sum), Real(1)});
  changed = solution_sum <
            -10 * std::numeric_limits<Real>::epsilon() * scale;
  return finite(scale);
}

bool update_physical_residual(const std::vector<Real> &old,
                              const std::vector<Real> &solution,
                              const std::vector<Real> &source,
                              const PreparedOperator &op, const Real dt,
                              const bool interstitial,
                              const std::vector<Real> *cloudborne,
                              Real &overall_max,
                              ReplayDiagnosticDetail *detail = nullptr,
                              const int field = 0,
                              const std::uint32_t substep = 0) {
  Real residual_max = 0;
  Real scale_max = 0;
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    const std::size_t km1 = k == 0 ? 0 : k - 1;
    const std::size_t kp1 = std::min(k + 1, op.nlev() - 1);
    const Real dt_rate_above = dt * op.rate_above[k];
    const Real dt_rate_below = dt * op.rate_below[k];
    if (!finite(dt_rate_above) || !finite(dt_rate_below)) {
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
      if (source[k] != 0) {
        return false;
      }
      continue;
    }
    const Real coefficient_scale =
        std::max({Real(1), std::abs(dt_rate_above),
                  std::abs(dt_rate_below), std::abs(dt * source[k] /
                                                    value_scale)});
    if (!finite(coefficient_scale)) {
      return false;
    }
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
    if (!finite(residual) || !finite(scale)) {
      return false;
    }
    residual_max = std::max(residual_max, residual);
    scale_max = std::max(scale_max, scale);
  }
  const Real normalized =
      scale_max > 0 ? residual_max / scale_max : residual_max;
  if (!finite(normalized)) {
    return false;
  }
  overall_max = std::max(overall_max, normalized);
  if (detail != nullptr) {
    auto &maximum = detail->physical_residual[field];
    if (!maximum || normalized > maximum->normalized) {
      maximum = ReplayResidualMaximum{residual_max, scale_max, normalized,
                                      substep};
    }
  }
  return finite(overall_max);
}

bool checked_difference(const Real lhs, const Real rhs, Real &difference) {
  const Real maximum = std::numeric_limits<Real>::max();
  if ((rhs > 0 && lhs < -maximum + rhs) ||
      (rhs < 0 && lhs > maximum + rhs)) {
    return false;
  }
  difference = lhs - rhs;
  return finite(difference);
}

bool normalized_ledger_delta(const Real initial, const Real delta,
                             Real &normalized) {
  const Real scale = std::abs(initial);
  normalized = scale > 0 ? std::abs(delta) / scale : std::abs(delta);
  return finite(normalized);
}

} // namespace

AdvanceResult advance_backward_euler(const State &initial,
                                     const OperatorData &data,
                                     const Real dt,
                                     const AdvanceControl &control) {
  return advance_backward_euler_interval(initial, data, dt, 1, control);
}

AdvanceResult advance_backward_euler_interval(
    const State &initial, const OperatorData &data, const Real total_dt,
    const int substeps, const AdvanceControl &control) {
  PhaseStopwatch stopwatch(control.phase_timings);
  stopwatch.phase(&AdvancePhaseTimings::setup);
  AdvanceResult result;
  result.state = State(initial.nlev());
  ReplayDiagnosticDetail *detail = nullptr;
  if (control.collect_detailed_diagnostics &&
      control.collect_replay_diagnostics) {
    detail = &result.diagnostics.replay_detail.emplace();
  }
  const auto fail = [&](ReplayDiagnosticStage stage,
                        ReplayDiagnosticReason reason) {
    if (detail != nullptr) {
      detail->failure_stage = stage;
      detail->failure_reason = reason;
    }
  };
  if (substeps <= 0 || !validate(initial, data, total_dt, result.message)) {
    if (substeps <= 0) result.message = "substeps must be positive";
    fail(ReplayDiagnosticStage::diagnostics,
         ReplayDiagnosticReason::diagnostic_incomplete);
    return result;
  }
  const Real dt = total_dt / substeps;
  if (!(dt > 0) || !finite(dt)) {
    result.message = "internal timestep must be finite and positive";
    fail(ReplayDiagnosticStage::diagnostics,
         ReplayDiagnosticReason::diagnostic_incomplete);
    return result;
  }
  const PreparedOperator op = prepare_operator(data);
  if (!finite_prepared_operator(op)) {
    result.message = "derived mixing rates must be finite";
    fail(ReplayDiagnosticStage::assembly,
         ReplayDiagnosticReason::assembly_nonfinite);
    return result;
  }
  result.diagnostics.nsubmix = substeps;
  if (detail != nullptr)
    detail->expected_substeps = static_cast<std::uint32_t>(substeps);
  result.diagnostics.final_buffer = 0;
  result.diagnostics.activation_limiter_count = op.activation_limiter_count;
  result.diagnostics.max_inverse_mixing_time = op.max_inverse_mixing_time;
  result.diagnostics.first_level_controller_rate =
      op.first_level_controller_rate;
  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    if (!checked_column_ledgers(initial, op.pressure_weight,
                                result.diagnostics.initial_ledger,
                                result.message)) {
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_nonfinite);
      return result;
    }
    stopwatch.phase(&AdvancePhaseTimings::setup);
  }

  std::vector<MatrixFamily> families;
  for (int j = 0; j < num_constituents; ++j) {
    const auto coefficients = activation_coefficients(op, j);
    auto family = std::find_if(
        families.begin(), families.end(), [&](const MatrixFamily &candidate) {
          return same_coefficients(candidate.activation, coefficients);
        });
    if (family == families.end()) {
      MatrixFamily created;
      created.activation = coefficients;
      created.system = assemble_block_system(op, coefficients, dt);
      families.push_back(std::move(created));
      family = families.end() - 1;
    }
    family->constituents.push_back(j);
  }
  result.diagnostics.block_matrix_family_count =
      static_cast<int>(families.size());
  std::optional<ScalarSystem> droplet_system_cache;
  State previous(substeps > 1 ? initial.nlev() : 0);
  for (int step = 0; step < substeps; ++step) {
  const State &old = step == 0 ? initial : previous;
  Real matrix_residual = 0;
  bool branch_changed = false;
  for (auto &family : families) {
    if (!finite_block_system(family.system)) {
      result.message =
          "Backward Euler block matrix contains a nonfinite coefficient";
      fail(ReplayDiagnosticStage::assembly,
           ReplayDiagnosticReason::assembly_nonfinite);
      return result;
    }
    stopwatch.phase(&AdvancePhaseTimings::factorization);
    ReplayDiagnosticReason factor_reason =
        ReplayDiagnosticReason::factorization_unsafe_pivot;
    if (!factor_block_system(family.system, family.factorization,
                             detail != nullptr ? &factor_reason : nullptr)) {
      result.message = "Backward Euler block factorization pivot failure";
      fail(ReplayDiagnosticStage::factorization, factor_reason);
      return result;
    }
    ++result.diagnostics.block_factorization_count;
    for (const int j : family.constituents) {
      stopwatch.phase(&AdvancePhaseTimings::setup);
      std::vector<Vec2> rhs(op.nlev());
      for (std::size_t k = 0; k < op.nlev(); ++k) {
        rhs[k] = {old.cloudborne[j][k], old.interstitial[j][k]};
      }
      std::vector<Vec2> solution;
      stopwatch.phase(&AdvancePhaseTimings::solve);
      if (detail != nullptr) detail->status = ReplayDiagnosticStatus::partial;
      if (!solve_block_system(family.system, family.factorization, rhs,
                              solution)) {
        result.message =
            "Backward Euler block solve produced a nonfinite value";
        fail(ReplayDiagnosticStage::solve,
             ReplayDiagnosticReason::solve_nonfinite);
        return result;
      }
      ++result.diagnostics.block_rhs_solve_count;
      if (control.collect_detailed_diagnostics) {
        stopwatch.phase(&AdvancePhaseTimings::diagnostics);
        if (!normalized_block_residual(family.system, solution, rhs,
                                       matrix_residual, detail, j,
                                       static_cast<std::uint32_t>(step))) {
          result.message = "Backward Euler block residual is nonfinite";
          fail(ReplayDiagnosticStage::matrix_residual,
               ReplayDiagnosticReason::matrix_residual_nonfinite);
          return result;
        }
      }
      stopwatch.phase(&AdvancePhaseTimings::update);
      for (std::size_t k = 0; k < op.nlev(); ++k) {
        result.state.cloudborne[j][k] = solution[k].x;
        result.state.interstitial[j][k] = solution[k].y;
      }
      update_positivity(old.cloudborne[j], result.state.cloudborne[j],
                        result.diagnostics,
                        result.diagnostics.field_clipping[j]
                            .cloudborne_roundoff_negative_count,
                        result.diagnostics.field_clipping[j]
                            .cloudborne_material_negative_count,
                        result.diagnostics.field_clipping[j]
                            .cloudborne_raw_minimum);
      update_positivity(old.interstitial[j], result.state.interstitial[j],
                        result.diagnostics,
                        result.diagnostics.field_clipping[j]
                            .interstitial_roundoff_negative_count,
                        result.diagnostics.field_clipping[j]
                            .interstitial_material_negative_count,
                        result.diagnostics.field_clipping[j]
                            .interstitial_raw_minimum);
      if (detail != nullptr) {
        for (const int f : {1 + j, 26 + j}) {
          detail->raw_minimum_valid[f] = true;
          detail->negative_counts_valid[f] = true;
        }
      }
      bool this_branch_changed = false;
      if (!bottom_branch_status(old.interstitial[j],
                                old.cloudborne[j],
                                result.state.interstitial[j],
                                result.state.cloudborne[j],
                                this_branch_changed)) {
        result.message =
            "Backward Euler bottom branch diagnostic is nonfinite";
        fail(ReplayDiagnosticStage::branch,
             ReplayDiagnosticReason::branch_mismatch);
        return result;
      }
      branch_changed = branch_changed || this_branch_changed;
    }
  }

  stopwatch.phase(&AdvancePhaseTimings::setup);
  const auto droplet_source = droplet_activation_source(result.state, op);
  if (!finite_values(droplet_source)) {
    result.message =
        "Backward Euler droplet activation source contains a nonfinite value";
    fail(ReplayDiagnosticStage::assembly,
         ReplayDiagnosticReason::assembly_nonfinite);
    return result;
  }
  if (!droplet_system_cache) {
    droplet_system_cache.emplace(assemble_cloud_system(op, dt));
  }
  const ScalarSystem &droplet_system = *droplet_system_cache;
  if (!finite_scalar_system(droplet_system)) {
    result.message =
        "Backward Euler droplet matrix contains a nonfinite coefficient";
    fail(ReplayDiagnosticStage::assembly,
         ReplayDiagnosticReason::assembly_nonfinite);
    return result;
  }
  ScalarFactorization droplet_factorization;
  stopwatch.phase(&AdvancePhaseTimings::factorization);
  ReplayDiagnosticReason scalar_factor_reason =
      ReplayDiagnosticReason::factorization_unsafe_pivot;
  if (!factor_scalar_system(
          droplet_system, droplet_factorization,
          detail != nullptr ? &scalar_factor_reason : nullptr)) {
    result.message = "Backward Euler scalar factorization pivot failure";
    fail(ReplayDiagnosticStage::factorization, scalar_factor_reason);
    return result;
  }
  ++result.diagnostics.scalar_factorization_count;
  stopwatch.phase(&AdvancePhaseTimings::setup);
  std::vector<Real> droplet_rhs(op.nlev(), 0);
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    droplet_rhs[k] = old.droplet_number[k] + dt * droplet_source[k];
    if (!finite(droplet_rhs[k])) {
      result.message =
          "Backward Euler droplet right-hand side contains a nonfinite value";
      fail(ReplayDiagnosticStage::assembly,
           ReplayDiagnosticReason::assembly_nonfinite);
      return result;
    }
  }
  stopwatch.phase(&AdvancePhaseTimings::solve);
  if (!solve_scalar_system(droplet_system, droplet_factorization,
                           droplet_rhs, result.state.droplet_number)) {
    result.message =
        "Backward Euler droplet solve produced a nonfinite value";
    fail(ReplayDiagnosticStage::solve,
         ReplayDiagnosticReason::solve_nonfinite);
    return result;
  }
  ++result.diagnostics.scalar_rhs_solve_count;
  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    if (!normalized_scalar_residual(droplet_system,
                                    result.state.droplet_number, droplet_rhs,
                                    matrix_residual, detail,
                                    static_cast<std::uint32_t>(step))) {
      result.message = "Backward Euler droplet residual is nonfinite";
      fail(ReplayDiagnosticStage::matrix_residual,
           ReplayDiagnosticReason::matrix_residual_nonfinite);
      return result;
    }
    result.diagnostics.normalized_matrix_residual =
        std::max(result.diagnostics.normalized_matrix_residual,
                 matrix_residual);
  }
  stopwatch.phase(&AdvancePhaseTimings::update);
  update_positivity(old.droplet_number, result.state.droplet_number,
                    result.diagnostics,
                    result.diagnostics.droplet_roundoff_negative_count,
                    result.diagnostics.droplet_material_negative_count,
                    result.diagnostics.droplet_raw_minimum);
  if (detail != nullptr) {
    detail->raw_minimum_valid[0] = true;
    detail->negative_counts_valid[0] = true;
    ++detail->completed_substeps;
  }

  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    std::array<std::vector<Real>, num_constituents> aerosol_sources;
    for (int j = 0; j < num_constituents; ++j) {
      aerosol_sources[j] = aerosol_activation_source(result.state, op, j);
      if (!finite_values(aerosol_sources[j])) {
        result.message = "Backward Euler aerosol activation source contains a "
                         "nonfinite value";
        fail(ReplayDiagnosticStage::assembly,
             ReplayDiagnosticReason::assembly_nonfinite);
        return result;
      }
    }
    Real physical_residual = 0;
    if (!update_physical_residual(
            old.droplet_number, result.state.droplet_number,
            droplet_source, op, dt, false, nullptr, physical_residual,
            detail, 0, static_cast<std::uint32_t>(step))) {
      result.message = "Backward Euler physical residual is nonfinite";
      fail(ReplayDiagnosticStage::physical_residual,
           ReplayDiagnosticReason::physical_residual_nonfinite);
      return result;
    }
    for (int j = 0; j < num_constituents; ++j) {
      if (!update_physical_residual(
              old.cloudborne[j], result.state.cloudborne[j],
              aerosol_sources[j], op, dt, false, nullptr,
              physical_residual, detail, 26 + j,
              static_cast<std::uint32_t>(step)) ||
          !update_physical_residual(
              old.interstitial[j], result.state.interstitial[j],
              aerosol_sources[j], op, dt, true,
              &result.state.cloudborne[j], physical_residual,
              detail, 1 + j, static_cast<std::uint32_t>(step))) {
        result.message = "Backward Euler physical residual is nonfinite";
        fail(ReplayDiagnosticStage::physical_residual,
             ReplayDiagnosticReason::physical_residual_nonfinite);
        return result;
      }
    }
    result.diagnostics.normalized_physical_residual =
        std::max(result.diagnostics.normalized_physical_residual,
                 physical_residual);
    result.diagnostics.normalized_residual =
        std::max(result.diagnostics.normalized_matrix_residual,
                 result.diagnostics.normalized_physical_residual);
    if (!finite(result.diagnostics.normalized_residual)) {
      result.message = "Backward Euler normalized residual is nonfinite";
      fail(ReplayDiagnosticStage::physical_residual,
           ReplayDiagnosticReason::physical_residual_nonfinite);
      return result;
    }
    stopwatch.phase(&AdvancePhaseTimings::update);
  }
  if (branch_changed) {
    result.message =
        "Backward Euler bottom activation left the approved linear branch";
    fail(ReplayDiagnosticStage::branch,
         ReplayDiagnosticReason::branch_mismatch);
    return result;
  }
  if (result.diagnostics.positivity_failure) {
    result.message = "material negative value in Backward Euler solution";
    fail(ReplayDiagnosticStage::raw_state,
         ReplayDiagnosticReason::raw_material_negative);
    return result;
  }
  if (step + 1 < substeps) std::swap(previous, result.state);
  }
  if (!apply_zero_cloud_postprocess(result.state, data.cloud_fraction,
                                    result.message)) {
    fail(ReplayDiagnosticStage::projection,
         ReplayDiagnosticReason::projection_nonfinite);
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
    fail(ReplayDiagnosticStage::ledger,
         ReplayDiagnosticReason::ledger_nonfinite);
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
      result.message =
          "Backward Euler aerosol ledger decomposition is nonfinite";
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_nonfinite);
      return result;
    }
    result.diagnostics.total_ledger_delta[j] = total_delta;
    result.diagnostics.operator_ledger_delta[j] = operator_delta;
    if (detail != nullptr) {
      detail->clipping_ledger_valid[j] = true;
      detail->total_ledger_valid[j] = true;
      detail->operator_ledger_valid[j] = true;
    }
    Real normalized = 0;
    if (!normalized_ledger_delta(result.diagnostics.initial_ledger[j],
                                 operator_delta, normalized)) {
      result.message =
          "Backward Euler aerosol ledger normalization is nonfinite";
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_nonfinite);
      return result;
    }
    result.diagnostics.max_normalized_ledger_error =
        std::max(result.diagnostics.max_normalized_ledger_error, normalized);
  }
  const Real ledger_tolerance =
      std::max(Real(1e-11),
               5 * static_cast<Real>(op.nlev()) *
                   std::numeric_limits<Real>::epsilon());
  result.ok = result.diagnostics.normalized_residual <= 1e-11 &&
              result.diagnostics.max_normalized_ledger_error <=
                  ledger_tolerance;
  if (result.diagnostics.normalized_residual > 1e-11) {
    result.message = "Backward Euler normalized residual exceeds tolerance";
  } else if (result.diagnostics.max_normalized_ledger_error >
             ledger_tolerance) {
    result.message = "Backward Euler aerosol ledger error exceeds tolerance";
  } else {
    result.message = "ok";
  }
  if (detail != nullptr) {
    if (result.ok) {
      detail->status = ReplayDiagnosticStatus::complete;
    } else if (result.diagnostics.normalized_residual > 1e-11) {
      fail(result.diagnostics.normalized_matrix_residual > 1e-11
               ? ReplayDiagnosticStage::matrix_residual
               : ReplayDiagnosticStage::physical_residual,
           result.diagnostics.normalized_matrix_residual > 1e-11
               ? ReplayDiagnosticReason::matrix_residual_exceeded
               : ReplayDiagnosticReason::physical_residual_exceeded);
    } else {
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_mismatch);
    }
  }
  return result;
}

bool backward_euler_probe_is_rejected(const BackwardEulerProbe probe) {
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  switch (probe) {
  case BackwardEulerProbe::solution: {
    std::vector<Vec2> rhs = {{1, 1}};
    BlockSystem system{{{}}, {{1, 0, 0, 1}}, {{}}};
    BlockFactorization factorization;
    if (!factor_block_system(system, factorization)) {
      return false;
    }
    factorization.inverse_diagonal[0].a00 = nan;
    std::vector<Vec2> solution;
    return !solve_block_system(system, factorization, rhs, solution);
  }
  case BackwardEulerProbe::matrix_residual: {
    BlockSystem system{{{}}, {{1, 0, 0, 1}}, {{}}};
    Real normalized = 0;
    return !normalized_block_residual(system, {{nan, 1}}, {{1, 1}},
                                      normalized);
  }
  case BackwardEulerProbe::physical_residual: {
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
  case BackwardEulerProbe::ledger_normalization: {
    Real normalized = 0;
    return !normalized_ledger_delta(std::numeric_limits<Real>::denorm_min(),
                                    std::numeric_limits<Real>::max(),
                                    normalized);
  }
  case BackwardEulerProbe::ledger_decomposition: {
    Real difference = 0;
    return !checked_difference(std::numeric_limits<Real>::max(),
                               -std::numeric_limits<Real>::max(),
                               difference);
  }
  case BackwardEulerProbe::block_pivot: {
    const Real unsafe = 50 * std::numeric_limits<Real>::epsilon();
    BlockSystem system{{{}}, {{unsafe, 0, 0, 1}}, {{}}};
    BlockFactorization factorization;
    return !factor_block_system(system, factorization);
  }
  case BackwardEulerProbe::scalar_pivot: {
    const Real unsafe = 50 * std::numeric_limits<Real>::epsilon();
    ScalarSystem system{{0}, {unsafe}, {1}};
    ScalarFactorization factorization;
    return !factor_scalar_system(system, factorization);
  }
  case BackwardEulerProbe::material_negative: {
    AdvanceDiagnostics diagnostics;
    int roundoff = 0;
    int material = 0;
    Real minimum = std::numeric_limits<Real>::infinity();
    update_positivity({1}, {-1}, diagnostics, roundoff, material, minimum);
    return diagnostics.positivity_failure && material == 1 && roundoff == 0;
  }
  case BackwardEulerProbe::roundoff_negative: {
    AdvanceDiagnostics diagnostics;
    int roundoff = 0;
    int material = 0;
    Real minimum = std::numeric_limits<Real>::infinity();
    update_positivity({1}, {-std::numeric_limits<Real>::epsilon()},
                      diagnostics, roundoff, material, minimum);
    return !diagnostics.positivity_failure && roundoff == 1 && material == 0;
  }
  case BackwardEulerProbe::bottom_branch: {
    bool changed = false;
    return bottom_branch_status({1}, {0}, {-1}, {0}, changed) && changed;
  }
  }
  return false;
}

std::array<Real, 4>
backward_euler_bottom_block_for_probe(const Real rate_above,
                                      const Real activation,
                                      const Real dt) {
  PreparedOperator op;
  op.pressure_weight = {1, 1};
  op.rate_above = {0, rate_above};
  op.rate_below = {rate_above, 0};
  op.overlap_above = {1, 1};
  op.overlap_below = {1, 1};
  const BlockSystem system =
      assemble_block_system(op, {0, activation}, dt);
  const Mat2 &bottom = system.diagonal.back();
  return {bottom.a00, bottom.a01, bottom.a10, bottom.a11};
}

} // namespace ndrop_test
