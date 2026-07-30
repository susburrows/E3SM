#include "ndrop_integrator.hpp"

#include <mam4xx/ndrop.hpp>

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using mam4::ColumnView;
using mam4::DeviceType;
using mam4::ndrop::View1D;
using mam4::ndrop::View2D;
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::State;

constexpr int pver = mam4::ndrop::pver;
constexpr int top_lev = 1;
constexpr int active_levels = pver - top_lev;
constexpr int nmodes = mam4::AeroConfig::num_modes();
constexpr int ncnst = mam4::ndrop::ncnst_tot;
constexpr int nspec_max = mam4::ndrop::nspec_max;

static_assert(pver == 128, "comparison must be configured for L128");
static_assert(nmodes == ndrop_test::num_modes, "mode count mismatch");
static_assert(ncnst == ndrop_test::num_constituents,
              "constituent count mismatch");

struct ColumnStorage {
  View1D owned;
  ColumnView view;

  explicit ColumnStorage(const std::string &name)
      : owned(name, pver), view(owned.data(), pver) {}
};

struct Comparison {
  int expected_nsubmix = 0;
  int exact_final_buffer = -1;
  int canonical_final_buffer = -1;
  Real max_absolute_difference = 0;
  Real max_relative_difference = 0;
  int mismatches = 0;
};

enum class CloudProfile {
  full,
  fractional,
  isolated_zero,
  contiguous_zero,
  all_zero
};

void build_mode_map(int nspec_amode[nmodes],
                    int mam_idx[nmodes][nspec_max]) {
  const int species_per_mode[nmodes] = {7, 4, 7, 3};
  int one_based_constituent = 1;
  for (int mode = 0; mode < nmodes; ++mode) {
    nspec_amode[mode] = species_per_mode[mode];
    for (int slot = 0; slot < nspec_max; ++slot) {
      mam_idx[mode][slot] = 1;
    }
    for (int slot = 0; slot <= species_per_mode[mode]; ++slot) {
      mam_idx[mode][slot] = one_based_constituent++;
    }
  }
  if (one_based_constituent != ncnst + 1) {
    throw std::runtime_error("mode map does not cover all constituents");
  }
}

Real cloud_fraction(const int active_k, const CloudProfile profile) {
  if (profile == CloudProfile::full) {
    return 1;
  }
  if (profile == CloudProfile::isolated_zero) {
    return active_k == 31 ? 0 : Real(0.7);
  }
  if (profile == CloudProfile::contiguous_zero) {
    return active_k >= 40 && active_k <= 55 ? 0 : Real(0.7);
  }
  if (profile == CloudProfile::all_zero) {
    return 0;
  }
  constexpr std::array<Real, 4> cloud = {0.15, 0.4, 0.75, 1};
  return cloud[static_cast<std::size_t>(active_k) % cloud.size()];
}

State initial_state() {
  State state(active_levels);
  for (int h = 0; h < active_levels; ++h) {
    state.droplet_number[h] = 10 + 0.01 * h;
    for (int j = 0; j < ncnst; ++j) {
      state.interstitial[j][h] =
          1 + 0.01 * j + 0.001 * static_cast<Real>(h + 1);
      state.cloudborne[j][h] =
          0.2 + 0.002 * j + 0.0002 * static_cast<Real>(h + 1);
    }
  }
  return state;
}

OperatorData canonical_operator(const Real exchange,
                                const Real activation_rate,
                                const CloudProfile profile) {
  OperatorData op(active_levels);
  for (int h = 0; h < active_levels; ++h) {
    op.pressure_weight[h] = 1;
    op.cloud_fraction[h] = cloud_fraction(h, profile);
  }
  for (int interface = 1; interface < active_levels; ++interface) {
    op.interface_exchange[interface] = exchange;
  }
  for (int mode = 0; mode < nmodes; ++mode) {
    std::fill(op.number_activation_rate[mode].begin(),
              op.number_activation_rate[mode].end(), activation_rate);
    std::fill(op.mass_activation_rate[mode].begin(),
              op.mass_activation_rate[mode].end(),
              Real(0.5) * activation_rate);
    if (activation_rate <= exchange) {
      op.number_activation_rate[mode].back() = 0;
      op.mass_activation_rate[mode].back() = 0;
    }
  }
  return op;
}

void compare_value(const Real exact, const Real canonical,
                   Comparison &comparison) {
  const Real difference = std::abs(exact - canonical);
  const Real scale = std::max({std::abs(exact), std::abs(canonical), Real(1)});
  comparison.max_absolute_difference =
      std::max(comparison.max_absolute_difference, difference);
  comparison.max_relative_difference =
      std::max(comparison.max_relative_difference, difference / scale);
  if (difference > 100 * std::numeric_limits<Real>::epsilon() * scale) {
    ++comparison.mismatches;
  }
}

Comparison run_comparison(const Real exchange, const Real activation_rate,
                          const CloudProfile profile) {
  const State initial = initial_state();
  const OperatorData host_op =
      canonical_operator(exchange, activation_rate, profile);
  const auto canonical =
      ndrop_test::advance_forward_euler(initial, host_op, 1);
  if (!canonical.ok) {
    throw std::runtime_error("canonical Forward Euler failed: " +
                             canonical.message);
  }

  ColumnStorage csbot_store("csbot");
  ColumnStorage cldn_store("cldn");
  ColumnStorage zn_store("zn");
  ColumnStorage zs_store("zs");
  ColumnStorage eddy_store("eddy_diff");
  ColumnStorage qcld_store("qcld");
  ColumnStorage overlapp_store("overlapp");
  ColumnStorage overlapm_store("overlapm");
  ColumnStorage dkp_store("eddy_diff_kp");
  ColumnStorage dkm_store("eddy_diff_km");
  ColumnStorage qncld_store("qncld");
  ColumnStorage srcn_store("srcn");
  ColumnStorage source_store("source");

  auto &csbot = csbot_store.view;
  auto &cldn = cldn_store.view;
  auto &zn = zn_store.view;
  auto &zs = zs_store.view;
  auto &eddy_diff = eddy_store.view;
  auto &qcld = qcld_store.view;
  auto &overlapp = overlapp_store.view;
  auto &overlapm = overlapm_store.view;
  auto &eddy_diff_kp = dkp_store.view;
  auto &eddy_diff_km = dkm_store.view;
  auto &qncld = qncld_store.view;
  auto &srcn = srcn_store.view;
  auto &source = source_store.view;

  Kokkos::deep_copy(csbot, Real(1));
  Kokkos::deep_copy(cldn, Real(1));
  Kokkos::deep_copy(zn, Real(1));
  Kokkos::deep_copy(zs, Real(1));
  Kokkos::deep_copy(eddy_diff, Real(0));
  Kokkos::deep_copy(qcld, Real(0));
  Kokkos::deep_copy(overlapp, Real(0));
  Kokkos::deep_copy(overlapm, Real(0));
  Kokkos::deep_copy(eddy_diff_kp, Real(0));
  Kokkos::deep_copy(eddy_diff_km, Real(0));
  Kokkos::deep_copy(qncld, Real(0));
  Kokkos::deep_copy(srcn, Real(0));
  Kokkos::deep_copy(source, Real(0));

  for (int h = 0; h < active_levels; ++h) {
    const int k = h + top_lev;
    cldn(k) = cloud_fraction(h, profile);
    qcld(k) = initial.droplet_number[h];
    if (h + 1 < active_levels) {
      eddy_diff(k) = exchange;
    }
  }

  View2D nact("nact", pver, nmodes);
  View2D mact("mact", pver, nmodes);
  Kokkos::deep_copy(nact, Real(0));
  Kokkos::deep_copy(mact, Real(0));
  for (int h = 0; h < active_levels; ++h) {
    const int k = h + top_lev;
    for (int mode = 0; mode < nmodes; ++mode) {
      const bool limiter_inactive_bottom =
          h + 1 == active_levels && activation_rate <= exchange;
      nact(k, mode) = limiter_inactive_bottom ? 0 : activation_rate;
      mact(k, mode) =
          limiter_inactive_bottom ? 0 : Real(0.5) * activation_rate;
    }
  }

  View1D raercol[pver][2];
  View1D raercol_cw[pver][2];
  for (int k = 0; k < pver; ++k) {
    for (int buffer = 0; buffer < 2; ++buffer) {
      raercol[k][buffer] =
          View1D("raercol_" + std::to_string(k) + "_" +
                     std::to_string(buffer),
                 ncnst);
      raercol_cw[k][buffer] =
          View1D("raercol_cw_" + std::to_string(k) + "_" +
                     std::to_string(buffer),
                 ncnst);
      Kokkos::deep_copy(raercol[k][buffer], Real(0));
      Kokkos::deep_copy(raercol_cw[k][buffer], Real(0));
    }
  }
  for (int h = 0; h < active_levels; ++h) {
    const int k = h + top_lev;
    for (int j = 0; j < ncnst; ++j) {
      raercol[k][0](j) = initial.interstitial[j][h];
      raercol_cw[k][0](j) = initial.cloudborne[j][h];
      // Poison the destination to detect an incomplete write.
      raercol[k][1](j) = 7 * initial.interstitial[j][h];
      raercol_cw[k][1](j) = 7 * initial.cloudborne[j][h];
    }
  }

  int nspec_amode[nmodes];
  int mam_idx[nmodes][nspec_max];
  build_mode_map(nspec_amode, mam_idx);

  typename DeviceType::template view_1d<int> indices("indices", 2);
  const auto policy = haero::ThreadTeamPolicy(1u, Kokkos::AUTO);
  Kokkos::parallel_for(
      "compare_mam4xx_forward_euler", policy,
      KOKKOS_LAMBDA(const mam4::ThreadTeam &team) {
        int nsav = 0;
        int nnew = 1;
        mam4::ndrop::update_from_explmix(
            team, Real(1), csbot, cldn, zn, zs, eddy_diff, nact, mact, qcld,
            raercol, raercol_cw, nsav, nnew, nspec_amode, mam_idx, true,
            top_lev, overlapp, overlapm, eddy_diff_kp, eddy_diff_km, qncld,
            srcn, source);
        indices(0) = nsav;
        indices(1) = nnew;
      });
  Kokkos::fence();

  Comparison comparison;
  comparison.expected_nsubmix = canonical.diagnostics.nsubmix;
  comparison.exact_final_buffer = indices(1);
  comparison.canonical_final_buffer =
      canonical.diagnostics.final_buffer;
  if (indices(1) != canonical.diagnostics.final_buffer) {
    ++comparison.mismatches;
  }
  for (int h = 0; h < active_levels; ++h) {
    const int k = h + top_lev;
    compare_value(qcld(k), canonical.state.droplet_number[h], comparison);
    for (int j = 0; j < ncnst; ++j) {
      compare_value(raercol[k][indices(1)](j),
                    canonical.state.interstitial[j][h], comparison);
      compare_value(raercol_cw[k][indices(1)](j),
                    canonical.state.cloudborne[j][h], comparison);
    }
  }
  return comparison;
}

} // namespace

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  int failures = 0;
  try {
    struct TestCase {
      const char *name;
      Real exchange;
      Real activation_rate;
      CloudProfile profile;
      int expected_nsubmix;
    };
    constexpr std::array<TestCase, 9> cases = {
        TestCase{"one_full_inactive", 0.01, 0.001, CloudProfile::full, 1},
        TestCase{"two_full_inactive", 0.5, 0.001, CloudProfile::full, 2},
        TestCase{"three_full_inactive", 1, 0.001, CloudProfile::full, 3},
        TestCase{"five_full_inactive", 2, 0.001, CloudProfile::full, 5},
        TestCase{"three_fractional_inactive", 1, 0.001,
                 CloudProfile::fractional, 3},
        TestCase{"two_fractional_limiter", 0.5, 10,
                 CloudProfile::fractional, 2},
        TestCase{"isolated_zero", 0.5, 0.001,
                 CloudProfile::isolated_zero, 2},
        TestCase{"contiguous_zero", 0.5, 0.001,
                 CloudProfile::contiguous_zero, 2},
        TestCase{"all_zero", 0.5, 0.001, CloudProfile::all_zero, 2}};

    for (const auto &test : cases) {
      const Comparison result =
          run_comparison(test.exchange, test.activation_rate,
                         test.profile);
      const bool pass =
          result.expected_nsubmix == test.expected_nsubmix &&
          result.mismatches == 0;
      failures += pass ? 0 : 1;
      std::cout << std::setprecision(17) << "CASE name=" << test.name
                << " status=" << (pass ? "PASS" : "FAIL")
                << " nsubmix=" << result.expected_nsubmix
                << " exact_final_buffer=" << result.exact_final_buffer
                << " canonical_final_buffer="
                << result.canonical_final_buffer
                << " max_abs_difference="
                << result.max_absolute_difference
                << " max_rel_difference="
                << result.max_relative_difference
                << " mismatches=" << result.mismatches << '\n';
    }
  } catch (const std::exception &error) {
    ++failures;
    std::cerr << "ERROR " << error.what() << '\n';
  }
  Kokkos::finalize();
  std::cout << "RESULT failed_cases=" << failures << '\n';
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
