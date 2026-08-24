#ifndef NDROP_BENCHMARK_HPP
#define NDROP_BENCHMARK_HPP

#include "ndrop_integrator.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ndrop_test::benchmark {

enum class Method { forward_euler, imex, backward_euler };

struct Case {
  int id = 0;
  std::size_t nlev = 0;
  int active_pairs = 0;
  int batch_size = 0;
  int mixing_index = 0;
  Real mixing_scale = 0;
  Real dt = 0;
  Method method = Method::forward_euler;
  int forced_fe_substeps = 0;
};

struct BootstrapInterval {
  Real median = 0;
  Real lower = 0;
  Real upper = 0;
  Real relative_half_width = 0;
};

struct RichardsonErrors {
  std::array<Real, 51> weighted_l1_by_field{};
  std::array<Real, 51> linf_by_field{};
  Real maximum_weighted_l1 = 0;
  Real maximum_linf = 0;
  Real screening_error = 0;
};

const char *method_name(Method method);
std::vector<Case> required_cases(bool smoke = false);
Real sample_median(std::vector<Real> samples);
BootstrapInterval bootstrap_median_interval(
    const std::vector<Real> &samples, int resamples,
    std::uint64_t seed = UINT64_C(0x4e44524f));
std::uint64_t state_checksum(const State &state);
std::size_t estimated_workspace_bytes(const Case &benchmark_case);
bool timing_path_matches_validated(const AdvanceResult &validated,
                                   const AdvanceResult &lean);
RichardsonErrors richardson_errors(const State &coarse, const State &refined,
                                   const State &initial,
                                   const std::vector<Real> &weights);

} // namespace ndrop_test::benchmark

#endif
