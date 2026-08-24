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
using mam4::ndrop::View3D;
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::State;

constexpr int pver = mam4::ndrop::pver;
constexpr int nmodes = mam4::AeroConfig::num_modes();
constexpr int ncnst = mam4::ndrop::ncnst_tot;

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
  int canonical_nsubmix = 0;
  int exact_nsubmix = 0;
  int physical_only_nsubmix = 0;
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

State initial_state(const int top_lev) {
  const int active_levels = pver - top_lev;
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
                                const CloudProfile profile,
                                const bool distinct_modes,
                                const int top_lev,
                                const bool first_interface_only) {
  const int active_levels = pver - top_lev;
  OperatorData op(active_levels);
  for (int h = 0; h < active_levels; ++h) {
    op.pressure_weight[h] = 1;
    op.cloud_fraction[h] = cloud_fraction(h, profile);
  }
  for (int interface = 1; interface < active_levels; ++interface) {
    op.interface_exchange[interface] =
        first_interface_only && interface > 1 ? Real(0) : exchange;
  }
  // Current MAM4xx duplicates eddy_diff(top_lev) in its first-level
  // timescale controller. The duplicate self-cancels in every physical
  // tendency and therefore remains separate from interface_exchange.
  op.first_level_controller_rate = exchange;
  for (int mode = 0; mode < nmodes; ++mode) {
    const Real number_rate =
        activation_rate * (distinct_modes ? Real(mode + 1) : Real(1));
    const Real mass_rate =
        activation_rate * (distinct_modes ? Real(mode + 5) : Real(0.5));
    std::fill(op.number_activation_rate[mode].begin(),
              op.number_activation_rate[mode].end(), number_rate);
    std::fill(op.mass_activation_rate[mode].begin(),
              op.mass_activation_rate[mode].end(), mass_rate);
    if (number_rate <= exchange && mass_rate <= exchange) {
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
                          const CloudProfile profile,
                          const int top_lev,
                          const bool distinct_modes = false,
                          const bool first_interface_only = false) {
  const int active_levels = pver - top_lev;
  const State initial = initial_state(top_lev);
  const OperatorData host_op =
      canonical_operator(exchange, activation_rate, profile, distinct_modes,
                         top_lev, first_interface_only);
  const auto canonical =
      ndrop_test::advance_forward_euler(initial, host_op, 1);
  if (!canonical.ok) {
    throw std::runtime_error("canonical Forward Euler failed: " +
                             canonical.message);
  }
  OperatorData physical_only = host_op;
  physical_only.first_level_controller_rate = 0;
  const int physical_only_nsubmix = ndrop_test::production_nsubmix(
      ndrop_test::prepare_operator(physical_only), 1);

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

  for (int h = 0; h < active_levels; ++h) {
    const int k = h + top_lev;
    cldn(k) = cloud_fraction(h, profile);
    qcld(k) = initial.droplet_number[h];
    if (h + 1 < active_levels &&
        (!first_interface_only || h == 0)) {
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
      const Real number_rate =
          activation_rate * (distinct_modes ? Real(mode + 1) : Real(1));
      const Real mass_rate =
          activation_rate * (distinct_modes ? Real(mode + 5) : Real(0.5));
      const bool limiter_inactive_bottom = h + 1 == active_levels &&
                                           number_rate <= exchange &&
                                           mass_rate <= exchange;
      nact(k, mode) = limiter_inactive_bottom ? 0 : number_rate;
      mact(k, mode) = limiter_inactive_bottom ? 0 : mass_rate;
    }
  }

  View3D raercol("raercol", pver, 2, ncnst);
  View3D raercol_cw("raercol_cw", pver, 2, ncnst);
  Kokkos::deep_copy(raercol, Real(0));
  Kokkos::deep_copy(raercol_cw, Real(0));
  for (int h = 0; h < active_levels; ++h) {
    const int k = h + top_lev;
    for (int j = 0; j < ncnst; ++j) {
      raercol(k, 0, j) = initial.interstitial[j][h];
      raercol_cw(k, 0, j) = initial.cloudborne[j][h];
      // Poison the destination to detect an incomplete write.
      raercol(k, 1, j) = 7 * initial.interstitial[j][h];
      raercol_cw(k, 1, j) = 7 * initial.cloudborne[j][h];
    }
  }

  typename DeviceType::template view_1d<int> indices("indices", 3);
  const auto policy = mam4::ThreadTeamPolicy(1u, Kokkos::AUTO);
  Kokkos::parallel_for(
      "compare_mam4xx_forward_euler", policy,
      KOKKOS_LAMBDA(const mam4::ThreadTeam &team) {
        int nsav = 0;
        int nnew = 1;
        int nsubmix = 0;
        mam4::ndrop::update_from_explmix(
            team, Real(1), csbot, cldn, zn, zs, eddy_diff, nact, mact, qcld,
            raercol, raercol_cw, nsav, nnew, nsubmix, top_lev, overlapp,
            overlapm, eddy_diff_kp, eddy_diff_km, qncld);
        indices(0) = nsav;
        indices(1) = nnew;
        indices(2) = nsubmix;
      });
  Kokkos::fence();

  Comparison comparison;
  comparison.canonical_nsubmix = canonical.diagnostics.nsubmix;
  comparison.exact_nsubmix = indices(2);
  comparison.physical_only_nsubmix = physical_only_nsubmix;
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
      compare_value(raercol(k, indices(1), j),
                    canonical.state.interstitial[j][h], comparison);
      compare_value(raercol_cw(k, indices(1), j),
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
      int top_lev;
      Real exchange;
      Real activation_rate;
      CloudProfile profile;
      int expected_nsubmix;
      int expected_physical_only_nsubmix;
      bool distinct_modes;
      bool first_interface_only;
    };
    constexpr std::array<TestCase, 15> cases = {
        TestCase{"top0_one_full_inactive", 0, 0.01, 0.001,
                 CloudProfile::full, 1, 1, false, false},
        TestCase{"top0_two_fractional", 0, 0.5, 0.001,
                 CloudProfile::fractional, 2, 2, false, false},
        TestCase{"top1_one_full_inactive", 1, 0.01, 0.001,
                 CloudProfile::full, 1, 1, false, false},
        TestCase{"top1_two_full_inactive", 1, 0.5, 0.001,
                 CloudProfile::full, 2, 2, false, false},
        TestCase{"top1_three_full_inactive", 1, 1, 0.001,
                 CloudProfile::full, 3, 3, false, false},
        TestCase{"top1_five_full_inactive", 1, 2, 0.001,
                 CloudProfile::full, 5, 5, false, false},
        TestCase{"top1_three_fractional_inactive", 1, 1, 0.001,
                 CloudProfile::fractional, 3, 3, false, false},
        TestCase{"top1_two_fractional_distinct_modes", 1, 0.5, 0.001,
                 CloudProfile::fractional, 2, 2, true, false},
        TestCase{"top1_two_fractional_limiter", 1, 0.5, 10,
                 CloudProfile::fractional, 2, 2, false, false},
        TestCase{"top1_isolated_zero", 1, 0.5, 0.001,
                 CloudProfile::isolated_zero, 2, 2, false, false},
        TestCase{"top1_contiguous_zero", 1, 0.5, 0.001,
                 CloudProfile::contiguous_zero, 2, 2, false, false},
        TestCase{"top1_all_zero", 1, 0.5, 0.001,
                 CloudProfile::all_zero, 2, 2, false, false},
        TestCase{"top6_one_full_inactive", 6, 0.01, 0.001,
                 CloudProfile::full, 1, 1, false, false},
        TestCase{"top6_two_fractional", 6, 0.5, 0.001,
                 CloudProfile::fractional, 2, 2, false, false},
        TestCase{"top6_controller_limiting", 6, 0.6, 0,
                 CloudProfile::full, 2, 1, false, true}};

    for (const auto &test : cases) {
      const Comparison result =
          run_comparison(test.exchange, test.activation_rate,
                         test.profile, test.top_lev, test.distinct_modes,
                         test.first_interface_only);
      const bool pass =
          result.canonical_nsubmix == test.expected_nsubmix &&
          result.exact_nsubmix == result.canonical_nsubmix &&
          result.physical_only_nsubmix ==
              test.expected_physical_only_nsubmix &&
          result.mismatches == 0;
      failures += pass ? 0 : 1;
      std::cout << std::setprecision(17) << "CASE name=" << test.name
                << " status=" << (pass ? "PASS" : "FAIL")
                << " top_lev=" << test.top_lev
                << " canonical_nsubmix=" << result.canonical_nsubmix
                << " exact_nsubmix=" << result.exact_nsubmix
                << " physical_only_nsubmix="
                << result.physical_only_nsubmix
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
