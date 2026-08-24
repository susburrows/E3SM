#include "ndrop_integrator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace ndrop_test {
namespace {

constexpr Real cloud_threshold = 1e-10;

// This mapping is intentionally local to the reference implementation. The
// oracle does not call the candidate integrators' mapping or source helpers.
constexpr std::array<int, num_constituents> reference_mode = {
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3};
constexpr std::array<bool, num_constituents> reference_number = {
    true,  false, false, false, false, false, false, false, true,
    false, false, false, false, true,  false, false, false, false,
    false, false, false, true,  false, false, false};
constexpr std::array<int, num_modes> reference_number_index = {0, 8, 13, 21};

struct ReferenceOperator {
  std::size_t nlev = 0;
  std::vector<Real> rate_below;
  std::vector<Real> rate_above;
  std::vector<Real> overlap_below;
  std::vector<Real> overlap_above;
  std::array<std::vector<Real>, num_modes> number_activation_rate;
  std::array<std::vector<Real>, num_modes> mass_activation_rate;
};

std::size_t droplet_index(const std::size_t k) { return k; }

std::size_t interstitial_index(const ReferenceOperator &op, const int field,
                               const std::size_t k) {
  return op.nlev + static_cast<std::size_t>(field) * op.nlev + k;
}

std::size_t cloudborne_index(const ReferenceOperator &op, const int field,
                             const std::size_t k) {
  return op.nlev + static_cast<std::size_t>(num_constituents) * op.nlev +
         static_cast<std::size_t>(field) * op.nlev + k;
}

std::size_t state_size(const std::size_t nlev) {
  return static_cast<std::size_t>(1 + 2 * num_constituents) * nlev;
}

bool finite_vector(const std::vector<Real> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](const Real value) { return std::isfinite(value); });
}

bool prepare_reference_operator(const OperatorData &data,
                                ReferenceOperator &op) {
  op.nlev = data.nlev();
  op.rate_below.assign(op.nlev, 0);
  op.rate_above.assign(op.nlev, 0);
  op.overlap_below.assign(op.nlev, 1);
  op.overlap_above.assign(op.nlev, 1);
  op.number_activation_rate = data.number_activation_rate;
  op.mass_activation_rate = data.mass_activation_rate;
  for (std::size_t k = 0; k < op.nlev; ++k) {
    op.rate_above[k] = data.interface_exchange[k] / data.pressure_weight[k];
    op.rate_below[k] =
        data.interface_exchange[k + 1] / data.pressure_weight[k];
    const std::size_t km1 = k == 0 ? 0 : k - 1;
    const std::size_t kp1 = std::min(k + 1, op.nlev - 1);
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
    for (int mode = 0; mode < num_modes; ++mode) {
      op.number_activation_rate[mode][k] =
          std::min(op.number_activation_rate[mode][k], op.rate_below[k]);
      op.mass_activation_rate[mode][k] =
          std::min(op.mass_activation_rate[mode][k], op.rate_below[k]);
    }
  }
  if (!finite_vector(op.rate_below) || !finite_vector(op.rate_above) ||
      !finite_vector(op.overlap_below) || !finite_vector(op.overlap_above)) {
    return false;
  }
  for (int mode = 0; mode < num_modes; ++mode) {
    if (!finite_vector(op.number_activation_rate[mode]) ||
        !finite_vector(op.mass_activation_rate[mode])) {
      return false;
    }
  }
  return true;
}

std::vector<Real> pack_state(const State &state) {
  const std::size_t nlev = state.nlev();
  std::vector<Real> packed(state_size(nlev), 0);
  for (std::size_t k = 0; k < nlev; ++k) {
    packed[droplet_index(k)] = state.droplet_number[k];
  }
  ReferenceOperator indexing;
  indexing.nlev = nlev;
  for (int field = 0; field < num_constituents; ++field) {
    for (std::size_t k = 0; k < nlev; ++k) {
      packed[interstitial_index(indexing, field, k)] =
          state.interstitial[field][k];
      packed[cloudborne_index(indexing, field, k)] =
          state.cloudborne[field][k];
    }
  }
  return packed;
}

State unpack_state(const std::vector<Real> &packed, const std::size_t nlev) {
  State state(nlev);
  ReferenceOperator indexing;
  indexing.nlev = nlev;
  for (std::size_t k = 0; k < nlev; ++k) {
    state.droplet_number[k] = packed[droplet_index(k)];
  }
  for (int field = 0; field < num_constituents; ++field) {
    for (std::size_t k = 0; k < nlev; ++k) {
      state.interstitial[field][k] =
          packed[interstitial_index(indexing, field, k)];
      state.cloudborne[field][k] =
          packed[cloudborne_index(indexing, field, k)];
    }
  }
  return state;
}

bool evaluate_reference_rhs(const ReferenceOperator &op,
                            const std::vector<Real> &state,
                            std::vector<Real> &rhs) {
  rhs.assign(state.size(), 0);
  for (std::size_t k = 0; k < op.nlev; ++k) {
    const std::size_t km1 = k == 0 ? 0 : k - 1;
    const std::size_t kp1 = std::min(k + 1, op.nlev - 1);
    const Real droplet = state[droplet_index(k)];
    const Real droplet_below = state[droplet_index(kp1)];
    const Real droplet_above = state[droplet_index(km1)];
    Real droplet_source = 0;
    for (int mode = 0; mode < num_modes; ++mode) {
      const int number_field = reference_number_index[mode];
      const Real rate = op.number_activation_rate[mode][k];
      droplet_source +=
          rate * state[interstitial_index(op, number_field, kp1)];
      if (k + 1 == op.nlev) {
        const Real bottom =
            rate * (state[interstitial_index(op, number_field, k)] +
                    state[cloudborne_index(op, number_field, k)]);
        droplet_source += std::max(Real(0), bottom);
      }
    }
    rhs[droplet_index(k)] =
        droplet_source +
        op.rate_below[k] *
            (op.overlap_below[k] * droplet_below - droplet) +
        op.rate_above[k] *
            (op.overlap_above[k] * droplet_above - droplet);

    for (int field = 0; field < num_constituents; ++field) {
      const int mode = reference_mode[field];
      const auto &rates = reference_number[field]
                              ? op.number_activation_rate[mode]
                              : op.mass_activation_rate[mode];
      const Real interstitial = state[interstitial_index(op, field, k)];
      const Real cloudborne = state[cloudborne_index(op, field, k)];
      const Real interstitial_below =
          state[interstitial_index(op, field, kp1)];
      const Real interstitial_above =
          state[interstitial_index(op, field, km1)];
      const Real cloudborne_below =
          state[cloudborne_index(op, field, kp1)];
      const Real cloudborne_above =
          state[cloudborne_index(op, field, km1)];
      Real activation_source = rates[k] * interstitial_below;
      if (k + 1 == op.nlev) {
        const Real bottom_rate = op.number_activation_rate[mode][k];
        activation_source =
            std::max(Real(0), bottom_rate * (interstitial + cloudborne));
      }
      rhs[cloudborne_index(op, field, k)] =
          activation_source +
          op.rate_below[k] *
              (op.overlap_below[k] * cloudborne_below - cloudborne) +
          op.rate_above[k] *
              (op.overlap_above[k] * cloudborne_above - cloudborne);
      rhs[interstitial_index(op, field, k)] =
          -activation_source +
          op.rate_below[k] *
              (interstitial_below - interstitial +
               cloudborne_below * (1 - op.overlap_below[k])) +
          op.rate_above[k] *
              (interstitial_above - interstitial +
               cloudborne_above * (1 - op.overlap_above[k]));
    }
  }
  return finite_vector(rhs);
}

void linear_combination(const std::vector<Real> &base, const Real step,
                        const std::array<const std::vector<Real> *, 7> &stages,
                        const std::array<Real, 7> &coefficients,
                        std::vector<Real> &output) {
  output.resize(base.size());
  for (std::size_t i = 0; i < base.size(); ++i) {
    Real increment = 0;
    for (std::size_t stage = 0; stage < stages.size(); ++stage) {
      if (stages[stage] != nullptr && coefficients[stage] != 0) {
        increment += coefficients[stage] * (*stages[stage])[i];
      }
    }
    output[i] = base[i] + step * increment;
  }
}

std::vector<Real> field_absolute_tolerances(const State &initial,
                                            const Real rtol) {
  const std::size_t nlev = initial.nlev();
  ReferenceOperator indexing;
  indexing.nlev = nlev;
  std::vector<Real> tolerances(state_size(nlev), 0);
  Real scale = 0;
  for (const Real value : initial.droplet_number) {
    scale = std::max(scale, std::abs(value));
  }
  scale = std::max(scale, Real(1e-30));
  for (std::size_t k = 0; k < nlev; ++k) {
    tolerances[droplet_index(k)] = rtol * scale;
  }
  for (int field = 0; field < num_constituents; ++field) {
    Real interstitial_scale = 0;
    Real cloudborne_scale = 0;
    for (std::size_t k = 0; k < nlev; ++k) {
      interstitial_scale =
          std::max(interstitial_scale,
                   std::abs(initial.interstitial[field][k]));
      cloudborne_scale =
          std::max(cloudborne_scale,
                   std::abs(initial.cloudborne[field][k]));
    }
    interstitial_scale = std::max(interstitial_scale, Real(1e-30));
    cloudborne_scale = std::max(cloudborne_scale, Real(1e-30));
    for (std::size_t k = 0; k < nlev; ++k) {
      tolerances[interstitial_index(indexing, field, k)] =
          rtol * interstitial_scale;
      tolerances[cloudborne_index(indexing, field, k)] =
          rtol * cloudborne_scale;
    }
  }
  return tolerances;
}

bool apply_reference_zero_cloud_postprocess(
    State &state, const std::vector<Real> &cloud_fraction) {
  for (std::size_t k = 0; k < state.nlev(); ++k) {
    if (cloud_fraction[k] != 0) {
      continue;
    }
    state.droplet_number[k] = 0;
    for (int field = 0; field < num_constituents; ++field) {
      const Real combined =
          state.interstitial[field][k] + state.cloudborne[field][k];
      if (!std::isfinite(combined)) {
        return false;
      }
      state.interstitial[field][k] = combined;
      state.cloudborne[field][k] = 0;
    }
  }
  return true;
}

} // namespace

ReferenceResult advance_reference_dopri54(const State &initial,
                                          const OperatorData &data,
                                          const Real dt,
                                          const ReferenceOptions &options) {
  ReferenceResult result;
  result.state = initial;
  std::string validation_message;
  if (!validate(initial, data, dt, validation_message)) {
    result.message = "reference input validation failed: " +
                     validation_message;
    return result;
  }
  if (!std::isfinite(options.relative_tolerance) ||
      options.relative_tolerance <= 0 ||
      options.maximum_attempted_steps <= 0) {
    result.message = "reference options must be finite and positive";
    return result;
  }
  ReferenceOperator op;
  if (!prepare_reference_operator(data, op)) {
    result.message = "reference operator contains a nonfinite coefficient";
    return result;
  }

  std::vector<Real> state = pack_state(initial);
  const std::vector<Real> atol =
      field_absolute_tolerances(initial, options.relative_tolerance);
  std::array<std::vector<Real>, 7> k;
  std::vector<Real> stage_state;
  std::vector<Real> fifth_order;
  std::vector<Real> fourth_order;
  const std::array<const std::vector<Real> *, 7> stages = {
      &k[0], &k[1], &k[2], &k[3], &k[4], &k[5], &k[6]};

  constexpr std::array<Real, 7> a2 = {Real(1) / 5, 0, 0, 0, 0, 0, 0};
  constexpr std::array<Real, 7> a3 = {Real(3) / 40, Real(9) / 40, 0, 0,
                                      0, 0, 0};
  constexpr std::array<Real, 7> a4 = {Real(44) / 45, Real(-56) / 15,
                                      Real(32) / 9, 0, 0, 0, 0};
  constexpr std::array<Real, 7> a5 = {
      Real(19372) / 6561, Real(-25360) / 2187, Real(64448) / 6561,
      Real(-212) / 729, 0, 0, 0};
  constexpr std::array<Real, 7> a6 = {
      Real(9017) / 3168, Real(-355) / 33, Real(46732) / 5247,
      Real(49) / 176, Real(-5103) / 18656, 0, 0};
  constexpr std::array<Real, 7> fifth = {
      Real(35) / 384, 0, Real(500) / 1113, Real(125) / 192,
      Real(-2187) / 6784, Real(11) / 84, 0};
  constexpr std::array<Real, 7> fourth = {
      Real(5179) / 57600, 0, Real(7571) / 16695, Real(393) / 640,
      Real(-92097) / 339200, Real(187) / 2100, Real(1) / 40};

  Real time = 0;
  Real step = dt / 8;
  int attempted_steps = 0;
  while (time < dt) {
    if (++attempted_steps > options.maximum_attempted_steps) {
      result.message = "reference exceeded maximum attempted steps";
      return result;
    }
    const Real remaining = dt - time;
    step = std::min(step, remaining);
    if (!(step > 0) || time + step == time) {
      result.message = "reference step size underflowed";
      return result;
    }
    if (!evaluate_reference_rhs(op, state, k[0])) {
      result.message = "reference stage 1 is nonfinite";
      return result;
    }
    linear_combination(state, step, stages, a2, stage_state);
    if (!finite_vector(stage_state) ||
        !evaluate_reference_rhs(op, stage_state, k[1])) {
      result.message = "reference stage 2 is nonfinite";
      return result;
    }
    linear_combination(state, step, stages, a3, stage_state);
    if (!finite_vector(stage_state) ||
        !evaluate_reference_rhs(op, stage_state, k[2])) {
      result.message = "reference stage 3 is nonfinite";
      return result;
    }
    linear_combination(state, step, stages, a4, stage_state);
    if (!finite_vector(stage_state) ||
        !evaluate_reference_rhs(op, stage_state, k[3])) {
      result.message = "reference stage 4 is nonfinite";
      return result;
    }
    linear_combination(state, step, stages, a5, stage_state);
    if (!finite_vector(stage_state) ||
        !evaluate_reference_rhs(op, stage_state, k[4])) {
      result.message = "reference stage 5 is nonfinite";
      return result;
    }
    linear_combination(state, step, stages, a6, stage_state);
    if (!finite_vector(stage_state) ||
        !evaluate_reference_rhs(op, stage_state, k[5])) {
      result.message = "reference stage 6 is nonfinite";
      return result;
    }
    linear_combination(state, step, stages, fifth, fifth_order);
    if (!finite_vector(fifth_order) ||
        !evaluate_reference_rhs(op, fifth_order, k[6])) {
      result.message = "reference stage 7 is nonfinite";
      return result;
    }
    linear_combination(state, step, stages, fourth, fourth_order);
    if (!finite_vector(fourth_order)) {
      result.message = "reference embedded solution is nonfinite";
      return result;
    }
    Real error_ratio = 0;
    for (std::size_t i = 0; i < state.size(); ++i) {
      const Real scale =
          atol[i] + options.relative_tolerance *
                        std::max(std::abs(state[i]),
                                 std::abs(fifth_order[i]));
      const Real ratio = std::abs(fifth_order[i] - fourth_order[i]) / scale;
      if (!std::isfinite(ratio)) {
        result.message = "reference error norm is nonfinite";
        return result;
      }
      error_ratio = std::max(error_ratio, ratio);
    }
    result.function_evaluations += 7;
    const Real factor =
        error_ratio == 0
            ? Real(5)
            : std::clamp(Real(0.9) * std::pow(error_ratio, Real(-0.2)),
                         Real(0.2), Real(5));
    if (error_ratio <= 1) {
      state.swap(fifth_order);
      ++result.accepted_steps;
      result.maximum_error_ratio =
          std::max(result.maximum_error_ratio, error_ratio);
      time = step == remaining ? dt : time + step;
      step *= factor;
    } else {
      ++result.rejected_steps;
      step *= std::min(Real(1), factor);
    }
  }

  if (!finite_vector(state)) {
    result.message = "reference final state is nonfinite";
    return result;
  }
  result.state = unpack_state(state, op.nlev);
  if (!apply_reference_zero_cloud_postprocess(result.state,
                                               data.cloud_fraction)) {
    result.message = "reference zero-cloud postprocess overflowed";
    return result;
  }
  result.ok = true;
  result.message = "ok";
  return result;
}

} // namespace ndrop_test
