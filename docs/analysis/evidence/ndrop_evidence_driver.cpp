#include <mam4xx/ndrop.hpp>

#include <Kokkos_Core.hpp>

#include <algorithm>
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
using mam4::Real;
using mam4::ndrop::View1D;
using mam4::ndrop::View2D;

constexpr int pver = mam4::ndrop::pver;
constexpr int nmodes = mam4::AeroConfig::num_modes();
constexpr int ncnst = mam4::ndrop::ncnst_tot;
constexpr int nspec_max = mam4::ndrop::nspec_max;

struct ColumnStorage {
  View1D owned;
  ColumnView view;

  explicit ColumnStorage(const std::string &name)
      : owned(name, pver), view(owned.data(), pver) {}
};

struct RunResult {
  int calculated_nsubmix = -1;
  int nsav = -1;
  int nnew = -1;
  bool buffer_parity_matches = false;
  Real max_abs_ledger_delta = 0;
  Real max_rel_ledger_delta = 0;
  Real min_post_value = std::numeric_limits<Real>::max();
  Real interstitial_change_l1 = 0;
  Real cloudborne_change_l1 = 0;
  Real droplet_change_l1 = 0;
  Real droplet_source_l1 = 0;
  Real post_fingerprint = 0;
  int ledger_violations = 0;
  bool all_post_finite = true;
  bool all_post_nonnegative = true;
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
    throw std::runtime_error("mode map does not cover 25 constituents");
  }
}

int calculate_nsubmix(const Real dtmicro, const int top_lev,
                      const ColumnView &zn, const ColumnView &zs,
                      const ColumnView &eddy_diff,
                      const ColumnView &csbot) {
  // Do not reproduce the suspected top_lev=0 access in host reporting. Let the
  // exact MAM4xx function provide the bounds evidence for that case.
  if (top_lev <= 0 || top_lev >= pver) {
    return -1;
  }
  Real dtmin = std::numeric_limits<Real>::max();
  for (int k = top_lev; k < pver; ++k) {
    const int km1 = std::max(k - 1, top_lev);
    const Real dkp = zn(k) * eddy_diff(k) * csbot(k) * zs(k);
    const Real dkm =
        zn(k) * eddy_diff(k - 1) * csbot(k - 1) * zs(km1);
    const Real tinv = dkp + dkm;
    if (tinv > 1e-6) {
      dtmin = std::min(dtmin, Real(1) / tinv);
    }
  }
  const Real candidate = Real(0.9) * dtmin;
  return static_cast<int>(dtmicro / candidate) + 1;
}

RunResult run_case(const int top_lev, const bool enabled,
                   const Real interior_diffusivity,
                   const Real destination_scale,
                   const Real upper_boundary_diffusivity) {
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
  Kokkos::deep_copy(eddy_diff, interior_diffusivity);
  Kokkos::deep_copy(qcld, Real(10));
  Kokkos::deep_copy(overlapp, Real(0));
  Kokkos::deep_copy(overlapm, Real(0));
  Kokkos::deep_copy(eddy_diff_kp, Real(0));
  Kokkos::deep_copy(eddy_diff_km, Real(0));
  Kokkos::deep_copy(qncld, Real(0));
  Kokkos::deep_copy(srcn, Real(0));
  Kokkos::deep_copy(source, Real(0));

  eddy_diff(pver - 1) = 0;
  if (top_lev > 0 && top_lev < pver) {
    eddy_diff(top_lev - 1) = upper_boundary_diffusivity;
  }

  View2D nact("nact", pver, nmodes);
  View2D mact("mact", pver, nmodes);
  Real initial_qcld[pver];
  const Real nact_rate = Real(0.2) * interior_diffusivity;
  const Real mact_rate = Real(0.1) * interior_diffusivity;
  for (int k = 0; k < pver; ++k) {
    qcld(k) += Real(0.01) * k;
    initial_qcld[k] = qcld(k);
    for (int mode = 0; mode < nmodes; ++mode) {
      nact(k, mode) = nact_rate;
      mact(k, mode) = mact_rate;
    }
  }

  View1D raercol[pver][2];
  View1D raercol_cw[pver][2];
  Real initial_i[pver][ncnst];
  Real initial_c[pver][ncnst];
  for (int k = 0; k < pver; ++k) {
    for (int time = 0; time < 2; ++time) {
      raercol[k][time] =
          View1D("raercol_" + std::to_string(k) + "_" + std::to_string(time),
                 ncnst);
      raercol_cw[k][time] = View1D(
          "raercol_cw_" + std::to_string(k) + "_" + std::to_string(time),
          ncnst);
    }
    for (int constituent = 0; constituent < ncnst; ++constituent) {
      const Real interstitial = Real(1) + Real(0.01) * constituent +
                                Real(0.001) * (k + 1);
      const Real cloudborne = Real(0.2) + Real(0.002) * constituent +
                              Real(0.0002) * (k + 1);
      initial_i[k][constituent] = interstitial;
      initial_c[k][constituent] = cloudborne;
      raercol[k][0](constituent) = interstitial;
      raercol_cw[k][0](constituent) = cloudborne;
      // Production dropmixnuc only loads buffer zero. Varying buffer one
      // detects reads of a destination that should have been overwritten.
      raercol[k][1](constituent) = destination_scale * interstitial;
      raercol_cw[k][1](constituent) = cloudborne;
    }
  }

  int nspec_amode[nmodes];
  int mam_idx[nmodes][nspec_max];
  build_mode_map(nspec_amode, mam_idx);

  const Real dtmicro = 1;
  RunResult result;
  result.calculated_nsubmix = calculate_nsubmix(
      dtmicro, top_lev, zn, zs, eddy_diff, csbot);

  typename DeviceType::template view_1d<int> indices("indices", 2);
  auto policy = haero::ThreadTeamPolicy(1u, Kokkos::AUTO);
  std::cout << "LAUNCH_UPDATE top_lev=" << top_lev << std::endl;
  Kokkos::parallel_for(
      "ndrop_evidence", policy,
      KOKKOS_LAMBDA(const mam4::ThreadTeam &team) {
        int nsav = 0;
        int nnew = 1;
        mam4::ndrop::update_from_explmix(
            team, dtmicro, csbot, cldn, zn, zs, eddy_diff, nact, mact, qcld,
            raercol, raercol_cw, nsav, nnew, nspec_amode, mam_idx, enabled,
            top_lev, overlapp, overlapm, eddy_diff_kp, eddy_diff_km, qncld,
            srcn, source);
        indices(0) = nsav;
        indices(1) = nnew;
      });
  Kokkos::fence();

  result.nsav = indices(0);
  result.nnew = indices(1);
  if (result.calculated_nsubmix > 0) {
    const int expected_nnew = result.calculated_nsubmix % 2 == 0 ? 0 : 1;
    const int expected_nsav = 1 - expected_nnew;
    result.buffer_parity_matches =
        result.nnew == expected_nnew && result.nsav == expected_nsav;
  }
  const Real tolerance = 100 * std::numeric_limits<Real>::epsilon();
  for (int k = std::max(top_lev, 0); k < pver; ++k) {
    result.all_post_finite =
        result.all_post_finite && std::isfinite(qcld(k));
    result.all_post_nonnegative =
        result.all_post_nonnegative && qcld(k) >= Real(0);
    result.min_post_value = std::min(result.min_post_value, qcld(k));
    result.droplet_change_l1 += std::abs(qcld(k) - initial_qcld[k]);
    result.droplet_source_l1 += std::abs(srcn(k));
  }
  for (int constituent = 0; constituent < ncnst; ++constituent) {
    Real pre = 0;
    Real post = 0;
    for (int k = std::max(top_lev, 0); k < pver; ++k) {
      const Real weight = Real(1) / zn(k);
      pre += weight *
             (initial_i[k][constituent] + initial_c[k][constituent]);
      const Real post_i = raercol[k][result.nnew](constituent);
      const Real post_c = raercol_cw[k][result.nnew](constituent);
      result.all_post_finite = result.all_post_finite &&
                               std::isfinite(post_i) && std::isfinite(post_c);
      result.all_post_nonnegative = result.all_post_nonnegative &&
                                    post_i >= Real(0) && post_c >= Real(0);
      post += weight * (post_i + post_c);
      result.min_post_value =
          std::min(result.min_post_value, std::min(post_i, post_c));
      result.interstitial_change_l1 +=
          std::abs(post_i - initial_i[k][constituent]);
      result.cloudborne_change_l1 +=
          std::abs(post_c - initial_c[k][constituent]);
      result.post_fingerprint +=
          (Real(1) + constituent + Real(0.001) * k) *
          (post_i + Real(3) * post_c);
    }
    const Real abs_delta = std::abs(post - pre);
    const Real rel_delta = abs_delta / std::max(std::abs(pre), Real(1));
    result.max_abs_ledger_delta =
        std::max(result.max_abs_ledger_delta, abs_delta);
    result.max_rel_ledger_delta =
        std::max(result.max_rel_ledger_delta, rel_delta);
    if (abs_delta > tolerance * std::max(std::abs(pre), Real(1))) {
      ++result.ledger_violations;
    }
  }
  return result;
}

void print_result(const std::string &scenario, const int top_lev,
                  const bool enabled, const Real diffusivity,
                  const Real destination_scale, const RunResult &result) {
  std::cout << std::setprecision(17)
            << "RESULT scenario=" << scenario << " pver=" << pver
            << " top_lev=" << top_lev << " enabled=" << enabled
            << " diffusivity=" << diffusivity
            << " destination_scale=" << destination_scale
            << " calculated_nsubmix=" << result.calculated_nsubmix
            << " nsav=" << result.nsav << " nnew=" << result.nnew
            << " buffer_parity_matches=" << result.buffer_parity_matches
            << " max_abs_ledger_delta=" << result.max_abs_ledger_delta
            << " max_rel_ledger_delta=" << result.max_rel_ledger_delta
            << " ledger_violations=" << result.ledger_violations
            << " min_post=" << result.min_post_value
            << " interstitial_change_l1=" << result.interstitial_change_l1
            << " cloudborne_change_l1=" << result.cloudborne_change_l1
            << " droplet_change_l1=" << result.droplet_change_l1
            << " droplet_source_l1=" << result.droplet_source_l1
            << " all_post_finite=" << result.all_post_finite
            << " all_post_nonnegative=" << result.all_post_nonnegative
            << " post_fingerprint=" << result.post_fingerprint << '\n';
}

} // namespace

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  int exit_code = 0;
  try {
    if (argc < 3) {
      throw std::runtime_error(
          "usage: ndrop_evidence boundary TOP [UPPER_DIFF] | "
          "ledger ENABLED DIFFUSIVITY DESTINATION_SCALE");
    }
    const std::string scenario = argv[1];
    if (scenario == "boundary") {
      const int top_lev = std::stoi(argv[2]);
      const Real upper_diff = argc > 3 ? std::stod(argv[3]) : Real(0);
      std::cout << "BEFORE_CALL scenario=boundary pver=" << pver
                << " top_lev=" << top_lev
                << " upper_boundary_diffusivity=" << upper_diff << std::endl;
      const auto result = run_case(top_lev, true, Real(0.01), Real(1),
                                   upper_diff);
      print_result("boundary", top_lev, true, Real(0.01), Real(1), result);
    } else if (scenario == "ledger") {
      if (argc != 5) {
        throw std::runtime_error(
            "usage: ndrop_evidence ledger ENABLED DIFFUSIVITY "
            "DESTINATION_SCALE");
      }
      const bool enabled = std::stoi(argv[2]) != 0;
      const Real diffusivity = std::stod(argv[3]);
      const Real destination_scale = std::stod(argv[4]);
      const int top_lev = 1;
      const auto result = run_case(top_lev, enabled, diffusivity,
                                   destination_scale, Real(0));
      print_result("ledger", top_lev, enabled, diffusivity,
                   destination_scale, result);
    } else {
      throw std::runtime_error("unknown scenario: " + scenario);
    }
  } catch (const std::exception &error) {
    std::cerr << "ERROR " << error.what() << '\n';
    exit_code = 2;
  }
  Kokkos::finalize();
  return exit_code;
}
