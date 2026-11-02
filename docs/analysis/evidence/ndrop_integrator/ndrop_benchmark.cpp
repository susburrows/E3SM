#include "ndrop_benchmark.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace ndrop_test::benchmark {
namespace {

std::uint64_t next_random(std::uint64_t &state) {
  state += UINT64_C(0x9e3779b97f4a7c15);
  std::uint64_t value = state;
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

void hash_real(std::uint64_t &hash, const Real value) {
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (bits >> (8 * byte)) & UINT64_C(0xff);
    hash *= UINT64_C(1099511628211);
  }
}

bool state_bitwise_equal(const State &lhs, const State &rhs) {
  if (lhs.nlev() != rhs.nlev() ||
      std::memcmp(lhs.droplet_number.data(), rhs.droplet_number.data(),
                  lhs.nlev() * sizeof(Real)) != 0) {
    return false;
  }
  for (int field = 0; field < num_constituents; ++field) {
    if (std::memcmp(lhs.interstitial[field].data(),
                    rhs.interstitial[field].data(),
                    lhs.nlev() * sizeof(Real)) != 0 ||
        std::memcmp(lhs.cloudborne[field].data(),
                    rhs.cloudborne[field].data(),
                    lhs.nlev() * sizeof(Real)) != 0) {
      return false;
    }
  }
  return true;
}

struct FieldErrors {
  Real weighted_l1 = 0;
  Real linf = 0;
};

FieldErrors normalized_field_errors(const std::vector<Real> &coarse,
                                    const std::vector<Real> &refined,
                                    const std::vector<Real> &initial,
                                    const std::vector<Real> &weights) {
  if (coarse.size() != refined.size() || coarse.size() != initial.size() ||
      coarse.size() != weights.size()) {
    throw std::invalid_argument("Richardson field dimensions do not match");
  }
  FieldErrors errors;
  Real refined_l1 = 0;
  Real initial_l1 = 0;
  Real refined_linf = 0;
  Real initial_linf = 0;
  for (std::size_t level = 0; level < coarse.size(); ++level) {
    if (!std::isfinite(coarse[level]) || !std::isfinite(refined[level]) ||
        !std::isfinite(initial[level]) || !std::isfinite(weights[level]) ||
        weights[level] < 0) {
      throw std::invalid_argument(
          "Richardson field data must be finite with nonnegative weights");
    }
    const Real difference = std::abs(coarse[level] - refined[level]);
    errors.weighted_l1 += weights[level] * difference;
    refined_l1 += weights[level] * std::abs(refined[level]);
    initial_l1 += weights[level] * std::abs(initial[level]);
    errors.linf = std::max(errors.linf, difference);
    refined_linf = std::max(refined_linf, std::abs(refined[level]));
    initial_linf = std::max(initial_linf, std::abs(initial[level]));
    if (!std::isfinite(errors.weighted_l1) ||
        !std::isfinite(refined_l1) || !std::isfinite(initial_l1) ||
        !std::isfinite(errors.linf) || !std::isfinite(refined_linf) ||
        !std::isfinite(initial_linf)) {
      throw std::invalid_argument("Richardson error arithmetic overflowed");
    }
  }
  const Real l1_scale = std::max(refined_l1, initial_l1);
  const Real linf_scale = std::max(refined_linf, initial_linf);
  if (l1_scale > 0) {
    errors.weighted_l1 /= l1_scale;
  }
  if (linf_scale > 0) {
    errors.linf /= linf_scale;
  }
  if (!std::isfinite(errors.weighted_l1) || !std::isfinite(errors.linf)) {
    throw std::invalid_argument("Richardson normalized error is nonfinite");
  }
  return errors;
}

} // namespace

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

const char *linear_phase_mode(const Method method) {
  switch (method) {
  case Method::forward_euler:
    return "none";
  case Method::imex:
    return "combined";
  case Method::backward_euler:
    return "separate";
  }
  return "unknown";
}

std::vector<Case> required_cases(const bool smoke) {
  if (smoke) {
    return {{0, 16, 8, 4, 0, 0.05, 0.5, Method::forward_euler, 3},
            {1, 16, 8, 4, 0, 0.05, 0.5, Method::imex, 0},
            {2, 16, 8, 4, 0, 0.05, 0.5, Method::backward_euler, 0}};
  }

  constexpr std::array<std::size_t, 9> levels = {1,  2,  3,  16, 32,
                                                  64, 72, 127, 128};
  constexpr std::array<int, 5> nonzero_aerosol_pairs = {1, 4, 8, 16, 25};
  constexpr std::array<int, 5> batches = {1, 16, 64, 256, 1024};
  constexpr std::array<Real, 2> mixing = {0.05, 0.5};
  constexpr std::array<Real, 6> intervals = {33.33, 50, 75, 100, 150, 300};
  constexpr std::array<int, 7> fe_substeps = {1, 2, 3, 5, 8, 16, 21};

  std::vector<Case> cases;
  cases.reserve(levels.size() * nonzero_aerosol_pairs.size() * batches.size() *
                mixing.size() * intervals.size() *
                (fe_substeps.size() + 2));
  int id = 0;
  for (const std::size_t nlev : levels) {
    for (const int pairs : nonzero_aerosol_pairs) {
      for (const int batch : batches) {
        for (std::size_t mixing_index = 0; mixing_index < mixing.size();
             ++mixing_index) {
          for (const Real dt : intervals) {
            for (const int substeps : fe_substeps) {
              cases.push_back({id++, nlev, pairs, batch,
                               static_cast<int>(mixing_index),
                               mixing[mixing_index], dt,
                               Method::forward_euler, substeps});
            }
            cases.push_back({id++, nlev, pairs, batch,
                             static_cast<int>(mixing_index),
                             mixing[mixing_index], dt, Method::imex, 0});
            cases.push_back({id++, nlev, pairs, batch,
                             static_cast<int>(mixing_index),
                             mixing[mixing_index], dt,
                             Method::backward_euler, 0});
          }
        }
      }
    }
  }
  return cases;
}

Real sample_median(std::vector<Real> samples) {
  if (samples.empty()) {
    throw std::invalid_argument("median requires a nonempty sample");
  }
  for (const Real sample : samples) {
    if (!std::isfinite(sample) || sample < 0) {
      throw std::invalid_argument("median samples must be finite and nonnegative");
    }
  }
  const std::size_t middle = samples.size() / 2;
  std::nth_element(samples.begin(), samples.begin() + middle, samples.end());
  const Real upper = samples[middle];
  if (samples.size() % 2 == 1) {
    return upper;
  }
  const Real lower = *std::max_element(samples.begin(),
                                       samples.begin() + middle);
  return lower + (upper - lower) / 2;
}

BootstrapInterval bootstrap_median_interval(
    const std::vector<Real> &samples, const int resamples,
    const std::uint64_t seed) {
  if (samples.empty() || resamples <= 0) {
    throw std::invalid_argument(
        "bootstrap requires samples and a positive resample count");
  }
  for (const Real sample : samples) {
    if (!std::isfinite(sample) || sample < 0) {
      throw std::invalid_argument(
          "bootstrap samples must be finite and nonnegative");
    }
  }
  std::uint64_t random_state = seed;
  std::vector<Real> replicate(samples.size());
  std::vector<Real> medians;
  medians.reserve(static_cast<std::size_t>(resamples));
  for (int draw = 0; draw < resamples; ++draw) {
    for (Real &value : replicate) {
      value = samples[next_random(random_state) % samples.size()];
    }
    medians.push_back(sample_median(replicate));
  }
  std::sort(medians.begin(), medians.end());
  const std::size_t lower_index = static_cast<std::size_t>(
      std::floor(0.025 * static_cast<Real>(medians.size() - 1)));
  const std::size_t upper_index = static_cast<std::size_t>(
      std::ceil(0.975 * static_cast<Real>(medians.size() - 1)));
  BootstrapInterval interval;
  interval.median = sample_median(samples);
  interval.lower = medians[lower_index];
  interval.upper = medians[upper_index];
  const Real half_width = (interval.upper - interval.lower) / 2;
  if (interval.median > 0) {
    interval.relative_half_width = half_width / interval.median;
  } else if (half_width == 0) {
    interval.relative_half_width = 0;
  } else {
    interval.relative_half_width = std::numeric_limits<Real>::infinity();
  }
  return interval;
}

std::uint64_t state_checksum(const State &state) {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash_real(hash, static_cast<Real>(state.nlev()));
  for (const Real value : state.droplet_number) {
    hash_real(hash, value);
  }
  for (int field = 0; field < num_constituents; ++field) {
    for (const Real value : state.interstitial[field]) {
      hash_real(hash, value);
    }
    for (const Real value : state.cloudborne[field]) {
      hash_real(hash, value);
    }
  }
  return hash;
}

std::size_t estimated_workspace_bytes(const Case &benchmark_case) {
  const std::size_t levels = benchmark_case.nlev;
  std::size_t scalar_count = 0;
  switch (benchmark_case.method) {
  case Method::forward_euler:
    scalar_count = 3 * 51 * levels;
    break;
  case Method::imex:
    scalar_count = (51 + 51 + 6) * levels;
    break;
  case Method::backward_euler:
    scalar_count = (51 + 8 * (12 + 8) + 8) * levels;
    break;
  }
  return scalar_count * sizeof(Real) *
         static_cast<std::size_t>(benchmark_case.batch_size);
}

bool timing_path_matches_validated(const AdvanceResult &validated,
                                   const AdvanceResult &lean) {
  return (!validated.ok || lean.ok) &&
         validated.diagnostics.clipping_count ==
             lean.diagnostics.clipping_count &&
         validated.diagnostics.material_negative_count ==
             lean.diagnostics.material_negative_count &&
         state_bitwise_equal(validated.state, lean.state);
}

RichardsonErrors richardson_errors(const State &coarse, const State &refined,
                                   const State &initial,
                                   const std::vector<Real> &weights) {
  if (coarse.nlev() != refined.nlev() || coarse.nlev() != initial.nlev() ||
      coarse.nlev() != weights.size()) {
    throw std::invalid_argument("Richardson state dimensions do not match");
  }
  RichardsonErrors result;
  const auto assign = [&result](const std::size_t field,
                                const FieldErrors &errors) {
    result.weighted_l1_by_field[field] = errors.weighted_l1;
    result.linf_by_field[field] = errors.linf;
    result.maximum_weighted_l1 =
        std::max(result.maximum_weighted_l1, errors.weighted_l1);
    result.maximum_linf = std::max(result.maximum_linf, errors.linf);
  };
  assign(0, normalized_field_errors(coarse.droplet_number,
                                    refined.droplet_number,
                                    initial.droplet_number, weights));
  for (int field = 0; field < num_constituents; ++field) {
    assign(1 + 2 * field,
           normalized_field_errors(coarse.interstitial[field],
                                   refined.interstitial[field],
                                   initial.interstitial[field], weights));
    assign(2 + 2 * field,
           normalized_field_errors(coarse.cloudborne[field],
                                   refined.cloudborne[field],
                                   initial.cloudborne[field], weights));
  }
  result.screening_error =
      std::max(result.maximum_weighted_l1, result.maximum_linf);
  return result;
}

ReferenceErrors reference_errors(const State &candidate,
                                 const State &reference,
                                 const State &initial,
                                 const std::vector<Real> &weights) {
  // The fieldwise algebra is shared; here the second state is the independent
  // DOPRI54 reference, so its norms are the approved reference denominators.
  return richardson_errors(candidate, reference, initial, weights);
}

} // namespace ndrop_test::benchmark
