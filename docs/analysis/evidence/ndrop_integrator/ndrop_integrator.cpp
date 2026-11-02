#include "ndrop_integrator.hpp"
#include "ndrop_phase_timer.hpp"

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

void classify_negative_values(const std::vector<Real> &previous,
                              const std::vector<Real> &raw,
                              AdvanceDiagnostics &diagnostics,
                              int &roundoff_count, int &material_count) {
  Real scale = 0;
  for (const Real value : previous) {
    scale = std::max(scale, std::abs(value));
  }
  for (const Real value : raw) {
    scale = std::max(scale, std::abs(value));
  }
  const Real threshold =
      10 * std::numeric_limits<Real>::epsilon() * scale;
  for (const Real value : raw) {
    if (value >= 0) {
      continue;
    }
    ++diagnostics.negative_count;
    if (-value > threshold) {
      ++diagnostics.material_negative_count;
      ++material_count;
    } else {
      ++diagnostics.roundoff_negative_count;
      ++roundoff_count;
    }
  }
}

bool update_explicit_physical_residual(
    const std::vector<Real> &old, const std::vector<Real> &raw,
    const std::vector<Real> &source, const PreparedOperator &op,
    const Real dt, const bool interstitial,
    const std::vector<Real> *cloudborne, Real &overall_max,
    ReplayDiagnosticDetail *detail = nullptr, const int field = 0,
    const std::uint32_t substep = 0) {
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
        std::max({std::abs(old[k]), std::abs(raw[k]),
                  std::abs(old[km1]), std::abs(old[kp1])});
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
                  std::abs(dt_rate_below)});
    const Real scaled_current = raw[k] / value_scale;
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

    Real below_bracket = old[kp1] / value_scale - old[k] / value_scale;
    Real above_bracket = old[km1] / value_scale - old[k] / value_scale;
    if (interstitial) {
      below_bracket += ((*cloudborne)[kp1] / value_scale) *
                       (1 - op.overlap_below[k]);
      above_bracket += ((*cloudborne)[km1] / value_scale) *
                       (1 - op.overlap_above[k]);
    } else {
      below_bracket =
          op.overlap_below[k] * (old[kp1] / value_scale) -
          old[k] / value_scale;
      above_bracket =
          op.overlap_above[k] * (old[km1] / value_scale) -
          old[k] / value_scale;
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
  if (detail != nullptr) {
    auto &maximum = detail->physical_residual[field];
    if (!maximum || normalized > maximum->normalized) {
      maximum = ReplayResidualMaximum{residual_max, scale_max, normalized,
                                      substep};
    }
  }
  return std::isfinite(overall_max);
}

bool normalized_ledger_delta(const Real initial, const Real delta,
                             Real &normalized) {
  const Real scale = std::abs(initial);
  normalized = scale > 0 ? std::abs(delta) / scale : std::abs(delta);
  return std::isfinite(normalized);
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
      !std::isfinite(op.max_inverse_mixing_time) ||
      !std::isfinite(op.first_level_controller_rate)) {
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
  if (!std::isfinite(data.first_level_controller_rate) ||
      data.first_level_controller_rate < 0) {
    message = "first-level controller rate must be finite and nonnegative";
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
  op.first_level_controller_rate = data.first_level_controller_rate;

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
    const Real controller_rate =
        k == 0 ? op.first_level_controller_rate : Real(0);
    op.max_inverse_mixing_time =
        std::max(op.max_inverse_mixing_time,
                 op.rate_below[k] + op.rate_above[k] + controller_rate);

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
    const Real controller_rate =
        k == 0 ? op.first_level_controller_rate : Real(0);
    const Real inverse_time =
        op.rate_below[k] + op.rate_above[k] + controller_rate;
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
                                    const ForwardEulerOptions &options,
                                    const AdvanceControl &control) {
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
  if (options.forced_substeps < 0) {
    result.message = "forced_substeps must be nonnegative";
    fail(ReplayDiagnosticStage::diagnostics,
         ReplayDiagnosticReason::diagnostic_incomplete);
    return result;
  }
  if (!validate(initial, data, dt, result.message)) {
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
  const int nsubmix = options.forced_substeps > 0
                          ? options.forced_substeps
                          : production_nsubmix(op, dt);
  if (nsubmix <= 0) {
    result.message = "calculated substep count is outside the supported range";
    fail(ReplayDiagnosticStage::diagnostics,
         ReplayDiagnosticReason::diagnostic_incomplete);
    return result;
  }

  result.diagnostics.nsubmix = nsubmix;
  if (detail != nullptr) {
    detail->expected_substeps = static_cast<std::uint32_t>(nsubmix);
  }
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
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_nonfinite);
      return result;
    }
    stopwatch.phase(&AdvancePhaseTimings::setup);
  }

  std::array<State, 2> buffers = {initial, State(initial.nlev())};
  int old_index = 0;
  int new_index = 1;
  const Real substep = dt / nsubmix;
  stopwatch.phase(&AdvancePhaseTimings::update);

  for (int step = 0; step < nsubmix; ++step) {
    if (detail != nullptr) {
      detail->status = ReplayDiagnosticStatus::partial;
    }
    if (step > 0) {
      std::swap(old_index, new_index);
    }
    const State &old = buffers[old_index];
    State &next = buffers[new_index];
    const auto nc_source = droplet_activation_source(old, op);
    if (!finite_values(nc_source)) {
      result.message = "droplet activation source contains a nonfinite value";
      fail(ReplayDiagnosticStage::assembly,
           ReplayDiagnosticReason::assembly_nonfinite);
      return result;
    }
    std::array<std::vector<Real>, num_constituents> sources;
    for (int j = 0; j < num_constituents; ++j) {
      sources[j] = aerosol_activation_source(old, op, j);
      if (!finite_values(sources[j])) {
        result.message = "aerosol activation source contains a nonfinite value";
        fail(ReplayDiagnosticStage::assembly,
             ReplayDiagnosticReason::assembly_nonfinite);
        return result;
      }
    }

    std::vector<Real> raw_droplet(initial.nlev(), 0);
    if (detail != nullptr) {
      detail->raw_minimum_valid[0] = false;
      detail->negative_counts_valid[0] = false;
      detail->fe_clipping_count_valid[0] = false;
      detail->fe_weighted_deficit_valid[0] = false;
    }
    for (std::size_t k = 0; k < initial.nlev(); ++k) {
      raw_droplet[k] =
          cloud_raw(old.droplet_number, nc_source, op, k, substep);
    }
    if (!finite_values(raw_droplet)) {
      result.message = "Forward Euler update produced a nonfinite value";
      fail(ReplayDiagnosticStage::raw_state,
           ReplayDiagnosticReason::raw_nonfinite);
      return result;
    }
    if (control.collect_detailed_diagnostics) {
      stopwatch.phase(&AdvancePhaseTimings::diagnostics);
      if (!update_explicit_physical_residual(
              old.droplet_number, raw_droplet, nc_source, op, substep,
              false, nullptr,
              result.diagnostics.normalized_physical_residual, detail, 0,
              static_cast<std::uint32_t>(step))) {
        result.message = "Forward Euler physical residual is nonfinite";
        if (detail != nullptr) detail->physical_residual[0].reset();
        fail(ReplayDiagnosticStage::physical_residual,
             ReplayDiagnosticReason::physical_residual_nonfinite);
        return result;
      }
      stopwatch.phase(&AdvancePhaseTimings::update);
    }
    classify_negative_values(
        old.droplet_number, raw_droplet, result.diagnostics,
        result.diagnostics.droplet_roundoff_negative_count,
        result.diagnostics.droplet_material_negative_count);
    for (std::size_t k = 0; k < initial.nlev(); ++k) {
      const Real raw = raw_droplet[k];
      if (options.apply_clipping) {
        if (!clipped_value(
                raw, op.pressure_weight[k],
                result.diagnostics.droplet_weighted_deficit,
                result.diagnostics.droplet_clipping_count,
                result.diagnostics.droplet_raw_minimum, result.diagnostics,
                next.droplet_number[k], result.message)) {
          fail(ReplayDiagnosticStage::clipping,
               ReplayDiagnosticReason::clipping_nonfinite);
          return result;
        }
      } else {
        if (!record_raw_value(raw, result.diagnostics.droplet_raw_minimum,
                              result.diagnostics, result.message)) {
          fail(ReplayDiagnosticStage::raw_state,
               ReplayDiagnosticReason::raw_nonfinite);
          return result;
        }
        next.droplet_number[k] = raw;
      }
    }
    if (detail != nullptr) {
      detail->raw_minimum_valid[0] = true;
      detail->negative_counts_valid[0] = true;
      detail->fe_clipping_count_valid[0] = true;
      detail->fe_weighted_deficit_valid[0] = true;
    }

    for (int j = 0; j < num_constituents; ++j) {
      auto &field_diag = result.diagnostics.field_clipping[j];
      if (detail != nullptr) {
        for (const int f : {1 + j, 26 + j}) {
          detail->raw_minimum_valid[f] = false;
          detail->negative_counts_valid[f] = false;
          detail->fe_clipping_count_valid[f] = false;
          detail->fe_weighted_deficit_valid[f] = false;
        }
      }
      std::vector<Real> raw_cloud(initial.nlev(), 0);
      std::vector<Real> raw_interstitial(initial.nlev(), 0);
      for (std::size_t k = 0; k < initial.nlev(); ++k) {
        raw_cloud[k] =
            cloud_raw(old.cloudborne[j], sources[j], op, k, substep);
        raw_interstitial[k] =
            interstitial_raw(old.interstitial[j], old.cloudborne[j],
                             sources[j], op, k, substep);
      }
      if (!finite_values(raw_cloud) || !finite_values(raw_interstitial)) {
        result.message = "Forward Euler update produced a nonfinite value";
        fail(ReplayDiagnosticStage::raw_state,
             ReplayDiagnosticReason::raw_nonfinite);
        return result;
      }
      if (control.collect_detailed_diagnostics) {
        stopwatch.phase(&AdvancePhaseTimings::diagnostics);
        const bool cloud_residual_ok = update_explicit_physical_residual(
                old.cloudborne[j], raw_cloud, sources[j], op, substep,
                false, nullptr,
                result.diagnostics.normalized_physical_residual, detail,
                26 + j, static_cast<std::uint32_t>(step));
        const bool interstitial_residual_ok = cloud_residual_ok &&
            update_explicit_physical_residual(
                old.interstitial[j], raw_interstitial, sources[j], op,
                substep, true, &old.cloudborne[j],
                result.diagnostics.normalized_physical_residual, detail,
                1 + j, static_cast<std::uint32_t>(step));
        if (!cloud_residual_ok || !interstitial_residual_ok) {
          if (detail != nullptr) {
            detail->physical_residual[cloud_residual_ok ? 1 + j : 26 + j]
                .reset();
          }
          result.message = "Forward Euler physical residual is nonfinite";
          fail(ReplayDiagnosticStage::physical_residual,
               ReplayDiagnosticReason::physical_residual_nonfinite);
          return result;
        }
        stopwatch.phase(&AdvancePhaseTimings::update);
      }
      classify_negative_values(
          old.cloudborne[j], raw_cloud, result.diagnostics,
          field_diag.cloudborne_roundoff_negative_count,
          field_diag.cloudborne_material_negative_count);
      classify_negative_values(
          old.interstitial[j], raw_interstitial, result.diagnostics,
          field_diag.interstitial_roundoff_negative_count,
          field_diag.interstitial_material_negative_count);
      for (std::size_t k = 0; k < initial.nlev(); ++k) {
        if (options.apply_clipping) {
          if (!clipped_value(raw_cloud[k], op.pressure_weight[k],
                             field_diag.cloudborne_weighted_deficit,
                             field_diag.cloudborne_count,
                             field_diag.cloudborne_raw_minimum,
                             result.diagnostics, next.cloudborne[j][k],
                             result.message) ||
              !clipped_value(raw_interstitial[k], op.pressure_weight[k],
                             field_diag.interstitial_weighted_deficit,
                             field_diag.interstitial_count,
                             field_diag.interstitial_raw_minimum,
                             result.diagnostics, next.interstitial[j][k],
                             result.message)) {
            fail(ReplayDiagnosticStage::clipping,
                 ReplayDiagnosticReason::clipping_nonfinite);
            return result;
          }
        } else {
          if (!record_raw_value(raw_cloud[k],
                                field_diag.cloudborne_raw_minimum,
                                result.diagnostics, result.message) ||
              !record_raw_value(raw_interstitial[k],
                                field_diag.interstitial_raw_minimum,
                                result.diagnostics, result.message)) {
            fail(ReplayDiagnosticStage::raw_state,
                 ReplayDiagnosticReason::raw_nonfinite);
            return result;
          }
          next.cloudborne[j][k] = raw_cloud[k];
          next.interstitial[j][k] = raw_interstitial[k];
        }
      }
      if (detail != nullptr) {
        for (const int f : {1 + j, 26 + j}) {
          detail->raw_minimum_valid[f] = true;
          detail->negative_counts_valid[f] = true;
          detail->fe_clipping_count_valid[f] = true;
          detail->fe_weighted_deficit_valid[f] = true;
        }
      }
    }
    if (detail != nullptr) {
      ++detail->completed_substeps;
    }
  }

  if (control.collect_detailed_diagnostics) {
    stopwatch.phase(&AdvancePhaseTimings::diagnostics);
    for (int j = 0; j < num_constituents; ++j) {
      const auto &field_diag = result.diagnostics.field_clipping[j];
      Real delta = 0;
      if (!checked_add(field_diag.interstitial_weighted_deficit,
                       field_diag.cloudborne_weighted_deficit, delta)) {
        result.message =
            "Forward Euler clipping ledger diagnostic overflowed";
        fail(ReplayDiagnosticStage::ledger,
             ReplayDiagnosticReason::ledger_nonfinite);
        return result;
      }
      result.diagnostics.clipping_ledger_delta[j] = delta;
      if (detail != nullptr) detail->clipping_ledger_valid[j] = true;
    }
    stopwatch.phase(&AdvancePhaseTimings::update);
  }
  result.diagnostics.final_buffer = new_index;
  result.state = buffers[new_index];
  if (!apply_zero_cloud_postprocess(result.state, data.cloud_fraction,
                                    result.message)) {
    fail(ReplayDiagnosticStage::projection,
         ReplayDiagnosticReason::projection_nonfinite);
    return result;
  }
  if (!control.collect_detailed_diagnostics) {
    result.diagnostics.solution_minimum = result.diagnostics.raw_minimum;
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
    if (!checked_add(result.diagnostics.final_ledger[j],
                     -result.diagnostics.initial_ledger[j], total_delta) ||
        !checked_add(total_delta,
                     -result.diagnostics.clipping_ledger_delta[j],
                     operator_delta)) {
      result.message = "Forward Euler ledger diagnostic overflowed";
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_nonfinite);
      return result;
    }
    result.diagnostics.total_ledger_delta[j] = total_delta;
    result.diagnostics.operator_ledger_delta[j] = operator_delta;
    if (detail != nullptr) {
      detail->total_ledger_valid[j] = true;
      detail->operator_ledger_valid[j] = true;
    }
    Real normalized = 0;
    if (!normalized_ledger_delta(result.diagnostics.initial_ledger[j],
                                 operator_delta, normalized)) {
      result.message =
          "Forward Euler aerosol ledger normalization is nonfinite";
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_nonfinite);
      return result;
    }
    result.diagnostics.max_normalized_ledger_error =
        std::max(result.diagnostics.max_normalized_ledger_error,
                 normalized);
  }
  result.diagnostics.normalized_residual =
      result.diagnostics.normalized_physical_residual;
  if (!std::isfinite(result.diagnostics.normalized_residual)) {
    result.message = "Forward Euler normalized residual is nonfinite";
    fail(ReplayDiagnosticStage::physical_residual,
         ReplayDiagnosticReason::physical_residual_nonfinite);
    return result;
  }
  result.diagnostics.solution_minimum = result.diagnostics.raw_minimum;
  const Real ledger_tolerance =
      std::max(Real(1e-11),
               5 * static_cast<Real>(op.nlev()) *
                   std::numeric_limits<Real>::epsilon());
  result.ok = result.diagnostics.normalized_residual <= 1e-11 &&
              result.diagnostics.max_normalized_ledger_error <=
                  ledger_tolerance;
  if (result.diagnostics.normalized_residual > 1e-11) {
    result.message = "Forward Euler normalized residual exceeds tolerance";
  } else if (result.diagnostics.max_normalized_ledger_error >
             ledger_tolerance) {
    result.message = "Forward Euler aerosol ledger error exceeds tolerance";
  } else {
    result.message = "ok";
  }
  if (detail != nullptr) {
    if (result.ok) {
      detail->status = ReplayDiagnosticStatus::complete;
    } else if (result.diagnostics.normalized_residual > 1e-11) {
      fail(ReplayDiagnosticStage::physical_residual,
           ReplayDiagnosticReason::physical_residual_exceeded);
    } else {
      fail(ReplayDiagnosticStage::ledger,
           ReplayDiagnosticReason::ledger_mismatch);
    }
  }
  return result;
}

} // namespace ndrop_test
