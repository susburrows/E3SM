#include "ndrop_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace ndrop_test {
namespace {

constexpr Real cloud_threshold = 1e-10;
constexpr Real mixing_rate_threshold = 1e-6;

constexpr std::array<int, num_constituents> mode_for_constituent = {
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3};

constexpr std::array<bool, num_constituents> number_constituent = {
    true,  false, false, false, false, false, false, false, true,
    false, false, false, false, true,  false, false, false, false,
    false, false, false, true,  false, false, false};

bool finite_nonnegative(const std::vector<Real> &values) {
  return std::all_of(values.begin(), values.end(), [](const Real value) {
    return std::isfinite(value) && value >= 0;
  });
}

bool checked_add(const Real lhs, const Real rhs, Real &sum) {
  if ((rhs > 0 && lhs > std::numeric_limits<Real>::max() - rhs) ||
      (rhs < 0 && lhs < -std::numeric_limits<Real>::max() - rhs)) {
    return false;
  }
  sum = lhs + rhs;
  return std::isfinite(sum);
}

bool checked_product(const Real lhs, const Real rhs, Real &product) {
  if (lhs != 0 && rhs != 0 &&
      std::abs(lhs) > std::numeric_limits<Real>::max() / std::abs(rhs)) {
    return false;
  }
  product = lhs * rhs;
  return std::isfinite(product);
}

bool record_raw_value(const Real raw, Real &field_minimum,
                      AdvanceDiagnostics &diagnostics,
                      std::string &message) {
  if (!std::isfinite(raw)) {
    message = "Forward Euler update produced a nonfinite value";
    return false;
  }
  field_minimum = std::min(field_minimum, raw);
  diagnostics.raw_minimum = std::min(diagnostics.raw_minimum, raw);
  return true;
}

bool clipped_value(const Real raw, const Real weight, Real &weighted_deficit,
                   int &count, Real &field_minimum,
                   AdvanceDiagnostics &diagnostics, Real &value,
                   std::string &message) {
  if (!record_raw_value(raw, field_minimum, diagnostics, message)) {
    return false;
  }
  if (raw >= 0) {
    value = raw;
    return true;
  }
  Real increment = 0;
  Real total = 0;
  if (!checked_product(weight, -raw, increment) ||
      !checked_add(weighted_deficit, increment, total)) {
    message = "Forward Euler clipping diagnostic overflowed";
    return false;
  }
  weighted_deficit = total;
  ++count;
  ++diagnostics.clipping_count;
  value = 0;
  return true;
}

Real cloud_raw(const std::vector<Real> &old, const std::vector<Real> &source,
               const PreparedOperator &op, const std::size_t k,
               const Real dt) {
  const std::size_t km1 = k == 0 ? 0 : k - 1;
  const std::size_t kp1 = std::min(k + 1, old.size() - 1);
  return old[k] +
         dt * (source[k] +
               op.rate_below[k] *
                   (op.overlap_below[k] * old[kp1] - old[k]) +
               op.rate_above[k] *
                   (op.overlap_above[k] * old[km1] - old[k]));
}

Real interstitial_raw(const std::vector<Real> &old,
                      const std::vector<Real> &cloud_old,
                      const std::vector<Real> &source,
                      const PreparedOperator &op, const std::size_t k,
                      const Real dt) {
  const std::size_t km1 = k == 0 ? 0 : k - 1;
  const std::size_t kp1 = std::min(k + 1, old.size() - 1);
  return old[k] +
         dt * (-source[k] +
               op.rate_below[k] *
                   (old[kp1] - old[k] +
                    cloud_old[kp1] * (1 - op.overlap_below[k])) +
               op.rate_above[k] *
                   (old[km1] - old[k] +
                    cloud_old[km1] * (1 - op.overlap_above[k])));
}

} // namespace

State::State(const std::size_t nlev) : droplet_number(nlev, 0) {
  for (auto &field : interstitial) {
    field.assign(nlev, 0);
  }
  for (auto &field : cloudborne) {
    field.assign(nlev, 0);
  }
}

std::size_t State::nlev() const { return droplet_number.size(); }

OperatorData::OperatorData(const std::size_t nlev)
    : pressure_weight(nlev, 1), interface_exchange(nlev + 1, 0),
      cloud_fraction(nlev, 1) {
  for (auto &field : number_activation_rate) {
    field.assign(nlev, 0);
  }
  for (auto &field : mass_activation_rate) {
    field.assign(nlev, 0);
  }
}

std::size_t OperatorData::nlev() const { return pressure_weight.size(); }

std::size_t PreparedOperator::nlev() const { return pressure_weight.size(); }

bool finite_values(const std::vector<Real> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](const Real value) { return std::isfinite(value); });
}

bool finite_prepared_operator(const PreparedOperator &op) {
  if (!finite_values(op.pressure_weight) || !finite_values(op.rate_below) ||
      !finite_values(op.rate_above) || !finite_values(op.overlap_below) ||
      !finite_values(op.overlap_above) ||
      !std::isfinite(op.max_inverse_mixing_time)) {
    return false;
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    if (!finite_values(op.number_activation_rate[mode]) ||
        !finite_values(op.mass_activation_rate[mode])) {
      return false;
    }
  }
  return true;
}

int constituent_mode(const int constituent) {
  if (constituent < 0 || constituent >= num_constituents) {
    throw std::out_of_range("constituent index outside [0,25)");
  }
  return mode_for_constituent[constituent];
}

bool constituent_is_number(const int constituent) {
  if (constituent < 0 || constituent >= num_constituents) {
    throw std::out_of_range("constituent index outside [0,25)");
  }
  return number_constituent[constituent];
}

std::vector<Real> aerosol_activation_source(const State &state,
                                            const PreparedOperator &op,
                                            const int constituent) {
  const auto nlev = state.nlev();
  const int mode = constituent_mode(constituent);
  const auto &rates = constituent_is_number(constituent)
                          ? op.number_activation_rate[mode]
                          : op.mass_activation_rate[mode];
  std::vector<Real> source(nlev, 0);
  for (std::size_t k = 0; k + 1 < nlev; ++k) {
    source[k] = rates[k] * state.interstitial[constituent][k + 1];
  }
  if (nlev > 0) {
    const std::size_t bottom = nlev - 1;
    // Preserve the current C++ bottom formula, including nact for mass.
    const Real bottom_rate = op.number_activation_rate[mode][bottom];
    const Real raw =
        bottom_rate * (state.interstitial[constituent][bottom] +
                       state.cloudborne[constituent][bottom]);
    source[bottom] = std::isfinite(raw) ? std::max(Real(0), raw) : raw;
  }
  return source;
}

std::vector<Real> droplet_activation_source(const State &state,
                                            const PreparedOperator &op) {
  const auto nlev = state.nlev();
  std::vector<Real> source(nlev, 0);
  constexpr std::array<int, num_modes> number_indices = {0, 8, 13, 21};
  for (int mode = 0; mode < num_modes; ++mode) {
    const int number_index = number_indices[mode];
    for (std::size_t k = 0; k < nlev; ++k) {
      const std::size_t kp1 = std::min(k + 1, nlev - 1);
      source[k] += op.number_activation_rate[mode][k] *
                   state.interstitial[number_index][kp1];
    }
    if (nlev > 0) {
      const std::size_t bottom = nlev - 1;
      const Real extra =
          op.number_activation_rate[mode][bottom] *
          (state.interstitial[number_index][bottom] +
           state.cloudborne[number_index][bottom]);
      source[bottom] +=
          std::isfinite(extra) ? std::max(Real(0), extra) : extra;
    }
  }
  return source;
}

std::vector<Real> apply_cloud_operator(const std::vector<Real> &field,
                                       const PreparedOperator &op) {
  if (field.size() != op.nlev()) {
    throw std::invalid_argument("cloud operator field size mismatch");
  }
  std::vector<Real> tendency(op.nlev(), 0);
  const std::vector<Real> zero_source(op.nlev(), 0);
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    tendency[k] = cloud_raw(field, zero_source, op, k, 1) - field[k];
  }
  return tendency;
}

std::vector<Real>
apply_interstitial_operator(const std::vector<Real> &interstitial,
                            const std::vector<Real> &cloudborne,
                            const PreparedOperator &op) {
  if (interstitial.size() != op.nlev() ||
      cloudborne.size() != op.nlev()) {
    throw std::invalid_argument("interstitial operator field size mismatch");
  }
  std::vector<Real> tendency(op.nlev(), 0);
  const std::vector<Real> zero_source(op.nlev(), 0);
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    tendency[k] =
        interstitial_raw(interstitial, cloudborne, zero_source, op, k, 1) -
        interstitial[k];
  }
  return tendency;
}

bool validate(const State &state, const OperatorData &data, const Real dt,
              std::string &message) {
  const auto nlev = state.nlev();
  if (nlev == 0) {
    message = "nlev must be positive";
    return false;
  }
  if (!std::isfinite(dt) || dt <= 0) {
    message = "dt must be finite and positive";
    return false;
  }
  if (data.pressure_weight.size() != nlev ||
      data.interface_exchange.size() != nlev + 1 ||
      data.cloud_fraction.size() != nlev) {
    message = "operator dimensions do not match state nlev";
    return false;
  }
  if (!finite_nonnegative(state.droplet_number)) {
    message = "droplet number must be finite and nonnegative";
    return false;
  }
  for (int j = 0; j < num_constituents; ++j) {
    if (state.interstitial[j].size() != nlev ||
        state.cloudborne[j].size() != nlev) {
      message = "aerosol field dimensions do not match state nlev";
      return false;
    }
    if (!finite_nonnegative(state.interstitial[j]) ||
        !finite_nonnegative(state.cloudborne[j])) {
      message = "aerosol state must be finite and nonnegative";
      return false;
    }
  }
  for (const Real weight : data.pressure_weight) {
    if (!std::isfinite(weight) || weight <= 0) {
      message = "pressure weights must be finite and positive";
      return false;
    }
  }
  if (!finite_nonnegative(data.interface_exchange)) {
    message = "interface exchange coefficients must be finite and nonnegative";
    return false;
  }
  if (data.interface_exchange.front() != 0 ||
      data.interface_exchange.back() != 0) {
    message = "top and bottom interface exchange coefficients must be zero";
    return false;
  }
  for (const Real cloud : data.cloud_fraction) {
    if (!std::isfinite(cloud) || cloud < 0 || cloud > 1) {
      message = "cloud fractions must be finite and in [0,1]";
      return false;
    }
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    if (data.number_activation_rate[mode].size() != nlev ||
        data.mass_activation_rate[mode].size() != nlev) {
      message = "activation-rate dimensions do not match state nlev";
      return false;
    }
    if (!finite_nonnegative(data.number_activation_rate[mode]) ||
        !finite_nonnegative(data.mass_activation_rate[mode])) {
      message = "activation rates must be finite and nonnegative";
      return false;
    }
  }
  message.clear();
  return true;
}

PreparedOperator prepare_operator(const OperatorData &data) {
  const auto nlev = data.nlev();
  PreparedOperator op;
  op.pressure_weight = data.pressure_weight;
  op.rate_below.assign(nlev, 0);
  op.rate_above.assign(nlev, 0);
  op.overlap_below.assign(nlev, 1);
  op.overlap_above.assign(nlev, 1);
  op.number_activation_rate = data.number_activation_rate;
  op.mass_activation_rate = data.mass_activation_rate;

  for (std::size_t k = 0; k < nlev; ++k) {
    op.rate_above[k] = data.interface_exchange[k] / data.pressure_weight[k];
    op.rate_below[k] =
        data.interface_exchange[k + 1] / data.pressure_weight[k];
    const std::size_t km1 = k == 0 ? 0 : k - 1;
    const std::size_t kp1 = std::min(k + 1, nlev - 1);
    op.overlap_below[k] =
        data.cloud_fraction[kp1] > cloud_threshold
            ? std::min(data.cloud_fraction[k] / data.cloud_fraction[kp1],
                       Real(1))
            : Real(1);
    op.overlap_above[k] =
        data.cloud_fraction[km1] > cloud_threshold
            ? std::min(data.cloud_fraction[k] / data.cloud_fraction[km1],
                       Real(1))
            : Real(1);
    op.max_inverse_mixing_time =
        std::max(op.max_inverse_mixing_time,
                 op.rate_below[k] + op.rate_above[k]);

    for (int mode = 0; mode < num_modes; ++mode) {
      if (op.number_activation_rate[mode][k] > op.rate_below[k]) {
        op.number_activation_rate[mode][k] = op.rate_below[k];
        ++op.activation_limiter_count;
      }
      if (op.mass_activation_rate[mode][k] > op.rate_below[k]) {
        op.mass_activation_rate[mode][k] = op.rate_below[k];
        ++op.activation_limiter_count;
      }
    }
  }
  return op;
}

int production_nsubmix(const PreparedOperator &op, const Real dt) {
  Real minimum_timescale = std::numeric_limits<Real>::max();
  for (std::size_t k = 0; k < op.nlev(); ++k) {
    const Real inverse_time = op.rate_below[k] + op.rate_above[k];
    if (inverse_time > mixing_rate_threshold) {
      minimum_timescale =
          std::min(minimum_timescale, Real(1) / inverse_time);
    }
  }
  const Real candidate_step = Real(0.9) * minimum_timescale;
  const Real substep_ratio = dt / candidate_step;
  if (!std::isfinite(substep_ratio) || substep_ratio < 0 ||
      substep_ratio >
          static_cast<Real>(std::numeric_limits<int>::max() - 1)) {
    return -1;
  }
  return static_cast<int>(substep_ratio) + 1;
}

std::array<Real, num_constituents>
column_ledgers(const State &state,
               const std::vector<Real> &pressure_weight) {
  std::array<Real, num_constituents> ledgers{};
  std::string message;
  if (!checked_column_ledgers(state, pressure_weight, ledgers, message)) {
    throw std::overflow_error(message);
  }
  return ledgers;
}

bool checked_column_ledgers(
    const State &state, const std::vector<Real> &pressure_weight,
    std::array<Real, num_constituents> &ledgers, std::string &message) {
  if (state.nlev() != pressure_weight.size()) {
    message = "ledger weight size does not match state";
    return false;
  }
  ledgers.fill(0);
  for (int j = 0; j < num_constituents; ++j) {
    Real ledger = 0;
    for (std::size_t k = 0; k < state.nlev(); ++k) {
      Real interstitial_term = 0;
      Real cloudborne_term = 0;
      Real updated = 0;
      if (!checked_product(pressure_weight[k],
                           state.interstitial[j][k],
                           interstitial_term) ||
          !checked_product(pressure_weight[k], state.cloudborne[j][k],
                           cloudborne_term) ||
          !checked_add(ledger, interstitial_term, updated) ||
          !checked_add(updated, cloudborne_term, ledger)) {
        message = "aerosol column ledger exceeds the representable range";
        return false;
      }
    }
    ledgers[j] = ledger;
  }
  message.clear();
  return true;
}

bool apply_zero_cloud_postprocess(
    State &state, const std::vector<Real> &cloud_fraction,
    std::string &message) {
  if (state.nlev() != cloud_fraction.size()) {
    message = "zero-cloud postprocess dimensions do not match state";
    return false;
  }
  for (std::size_t k = 0; k < state.nlev(); ++k) {
    if (cloud_fraction[k] != 0) {
      continue;
    }
    state.droplet_number[k] = 0;
    for (int j = 0; j < num_constituents; ++j) {
      Real combined = 0;
      if (!checked_add(state.interstitial[j][k],
                       state.cloudborne[j][k], combined)) {
        message = "zero-cloud aerosol resuspension overflowed";
        return false;
      }
      state.interstitial[j][k] = combined;
      state.cloudborne[j][k] = 0;
    }
  }
  message.clear();
  return true;
}

AdvanceResult advance_forward_euler(const State &initial,
                                    const OperatorData &data, const Real dt,
                                    const ForwardEulerOptions &options) {
  AdvanceResult result;
  result.state = State(initial.nlev());
  if (options.forced_substeps < 0) {
    result.message = "forced_substeps must be nonnegative";
    return result;
  }
  if (!validate(initial, data, dt, result.message)) {
    return result;
  }

  const PreparedOperator op = prepare_operator(data);
  if (!finite_prepared_operator(op)) {
    result.message = "derived mixing rates must be finite";
    return result;
  }
  const int nsubmix = options.forced_substeps > 0
                          ? options.forced_substeps
                          : production_nsubmix(op, dt);
  if (nsubmix <= 0) {
    result.message = "calculated substep count is outside the supported range";
    return result;
  }

  result.diagnostics.nsubmix = nsubmix;
  result.diagnostics.activation_limiter_count =
      op.activation_limiter_count;
  result.diagnostics.max_inverse_mixing_time =
      op.max_inverse_mixing_time;
  if (!checked_column_ledgers(initial, op.pressure_weight,
                              result.diagnostics.initial_ledger,
                              result.message)) {
    return result;
  }

  std::array<State, 2> buffers = {initial, State(initial.nlev())};
  int old_index = 0;
  int new_index = 1;
  const Real substep = dt / nsubmix;

  for (int step = 0; step < nsubmix; ++step) {
    if (step > 0) {
      std::swap(old_index, new_index);
    }
    const State &old = buffers[old_index];
    State &next = buffers[new_index];
    const auto nc_source = droplet_activation_source(old, op);
    if (!finite_values(nc_source)) {
      result.message = "droplet activation source contains a nonfinite value";
      return result;
    }
    std::array<std::vector<Real>, num_constituents> sources;
    for (int j = 0; j < num_constituents; ++j) {
      sources[j] = aerosol_activation_source(old, op, j);
      if (!finite_values(sources[j])) {
        result.message = "aerosol activation source contains a nonfinite value";
        return result;
      }
    }

    for (std::size_t k = 0; k < initial.nlev(); ++k) {
      const Real raw =
          cloud_raw(old.droplet_number, nc_source, op, k, substep);
      if (options.apply_clipping) {
        if (!clipped_value(
                raw, op.pressure_weight[k],
                result.diagnostics.droplet_weighted_deficit,
                result.diagnostics.droplet_clipping_count,
                result.diagnostics.droplet_raw_minimum, result.diagnostics,
                next.droplet_number[k], result.message)) {
          return result;
        }
      } else {
        if (!record_raw_value(raw, result.diagnostics.droplet_raw_minimum,
                              result.diagnostics, result.message)) {
          return result;
        }
        next.droplet_number[k] = raw;
      }
    }

    for (int j = 0; j < num_constituents; ++j) {
      auto &field_diag = result.diagnostics.field_clipping[j];
      for (std::size_t k = 0; k < initial.nlev(); ++k) {
        const Real raw_cloud =
            cloud_raw(old.cloudborne[j], sources[j], op, k, substep);
        const Real raw_interstitial =
            interstitial_raw(old.interstitial[j], old.cloudborne[j],
                             sources[j], op, k, substep);
        if (options.apply_clipping) {
          if (!clipped_value(raw_cloud, op.pressure_weight[k],
                             field_diag.cloudborne_weighted_deficit,
                             field_diag.cloudborne_count,
                             field_diag.cloudborne_raw_minimum,
                             result.diagnostics, next.cloudborne[j][k],
                             result.message) ||
              !clipped_value(raw_interstitial, op.pressure_weight[k],
                             field_diag.interstitial_weighted_deficit,
                             field_diag.interstitial_count,
                             field_diag.interstitial_raw_minimum,
                             result.diagnostics, next.interstitial[j][k],
                             result.message)) {
            return result;
          }
        } else {
          if (!record_raw_value(raw_cloud,
                                field_diag.cloudborne_raw_minimum,
                                result.diagnostics, result.message) ||
              !record_raw_value(raw_interstitial,
                                field_diag.interstitial_raw_minimum,
                                result.diagnostics, result.message)) {
            return result;
          }
          next.cloudborne[j][k] = raw_cloud;
          next.interstitial[j][k] = raw_interstitial;
        }
      }
    }
  }

  for (int j = 0; j < num_constituents; ++j) {
    const auto &field_diag = result.diagnostics.field_clipping[j];
    Real delta = 0;
    if (!checked_add(field_diag.interstitial_weighted_deficit,
                     field_diag.cloudborne_weighted_deficit, delta)) {
      result.message = "Forward Euler clipping ledger diagnostic overflowed";
      return result;
    }
    result.diagnostics.clipping_ledger_delta[j] = delta;
  }
  result.diagnostics.final_buffer = new_index;
  result.state = buffers[new_index];
  if (!apply_zero_cloud_postprocess(result.state, data.cloud_fraction,
                                    result.message)) {
    return result;
  }
  if (!checked_column_ledgers(result.state, op.pressure_weight,
                              result.diagnostics.final_ledger,
                              result.message)) {
    return result;
  }
  result.diagnostics.solution_minimum = result.diagnostics.raw_minimum;
  result.ok = true;
  result.message = "ok";
  return result;
}

} // namespace ndrop_test
