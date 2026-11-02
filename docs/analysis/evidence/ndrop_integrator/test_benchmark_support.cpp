#include "ndrop_benchmark.hpp"
#include "ndrop_integrator.hpp"

#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using ndrop_test::AdvanceControl;
using ndrop_test::AdvancePhaseTimings;
using ndrop_test::AdvanceResult;
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::State;
using ndrop_test::num_constituents;
using ndrop_test::num_modes;

int passed = 0;
int failed = 0;

void check(const bool condition, const std::string &message) {
  if (condition) {
    ++passed;
  } else {
    ++failed;
    std::cerr << "FAIL: " << message << '\n';
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

State make_state() {
  State state(3);
  for (std::size_t k = 0; k < state.nlev(); ++k) {
    state.droplet_number[k] = 3 + 0.2 * k;
    for (int field = 0; field < num_constituents; ++field) {
      state.interstitial[field][k] = 1 + 0.01 * field + 0.03 * k;
      state.cloudborne[field][k] = 0.2 + 0.002 * field + 0.01 * k;
    }
  }
  return state;
}

OperatorData make_operator() {
  OperatorData op(3);
  op.cloud_fraction = {0.8, 0.5, 0.7};
  op.interface_exchange = {0, 0.04, 0.06, 0};
  for (int mode = 0; mode < num_modes; ++mode) {
    op.number_activation_rate[mode][0] = 0.005 * (mode + 1);
    op.number_activation_rate[mode][1] = 0.004 * (mode + 1);
    op.mass_activation_rate[mode][0] = 0.003 * (mode + 1);
    op.mass_activation_rate[mode][1] = 0.002 * (mode + 1);
  }
  return op;
}

bool finite_nonnegative(const AdvancePhaseTimings &timings) {
  const Real values[] = {timings.setup, timings.factorization, timings.solve,
                         timings.combined_factorization_solve, timings.update,
                         timings.diagnostics};
  for (const Real value : values) {
    if (!std::isfinite(value) || value < 0) {
      return false;
    }
  }
  return true;
}

template <typename Advance>
void test_method(const std::string &name,
                 const ndrop_test::benchmark::Method method,
                 Advance advance) {
  const State initial = make_state();
  const OperatorData op = make_operator();
  const AdvanceResult detailed = advance(initial, op, 0.5, AdvanceControl{});

  AdvanceControl lean_control;
  lean_control.collect_detailed_diagnostics = false;
  const AdvanceResult lean = advance(initial, op, 0.5, lean_control);
  check(detailed.ok && lean.ok, name + " detailed and lean advances succeed");
  check(state_bitwise_equal(detailed.state, lean.state),
        name + " diagnostics-off state is bitwise identical");

  AdvancePhaseTimings timings;
  AdvanceControl timed_control;
  timed_control.phase_timings = &timings;
  const AdvanceResult timed = advance(initial, op, 0.5, timed_control);
  check(timed.ok && state_bitwise_equal(detailed.state, timed.state),
        name + " phase instrumentation is state-neutral");
  check(finite_nonnegative(timings), name + " phase timings are finite");
  check(timings.setup + timings.factorization + timings.solve +
                timings.combined_factorization_solve + timings.update +
                timings.diagnostics >
            0,
        name + " phase instrumentation records elapsed time");
  check(std::string(ndrop_test::benchmark::linear_phase_mode(method)) ==
            (method == ndrop_test::benchmark::Method::forward_euler
                 ? "none"
                 : method == ndrop_test::benchmark::Method::imex ? "combined"
                                                                  : "separate"),
        name + " has the exact linear phase mode");
  if (method == ndrop_test::benchmark::Method::forward_euler) {
    check(timings.factorization == 0 && timings.solve == 0 &&
              timings.combined_factorization_solve == 0,
          name + " has no linear solver timing");
  } else if (method == ndrop_test::benchmark::Method::imex) {
    check(timings.factorization == 0 && timings.solve == 0 &&
              timings.combined_factorization_solve > 0,
          name + " records only combined factorization-plus-solve timing");
  } else {
    check(timings.factorization > 0 && timings.solve > 0 &&
              timings.combined_factorization_solve == 0,
          name + " retains separate factorization and solve timing");
  }
}

void test_required_matrix() {
  const auto cases = ndrop_test::benchmark::required_cases();
  check(cases.size() == 24300, "required benchmark has 24,300 method rows");
  std::size_t fe = 0;
  std::size_t imex = 0;
  std::size_t be = 0;
  bool ids_are_contiguous = true;
  bool composite_keys_are_unique = true;
  bool occupancies_are_approved = true;
  std::set<std::tuple<std::size_t, int, int, int, Real, Real, int, int>>
      composite_keys;
  std::map<std::size_t, std::size_t> level_counts;
  std::map<int, std::size_t> occupancy_counts;
  std::map<int, std::size_t> batch_counts;
  std::map<int, std::size_t> mixing_counts;
  std::map<Real, std::size_t> interval_counts;
  std::map<int, std::size_t> fe_substep_counts;
  for (std::size_t index = 0; index < cases.size(); ++index) {
    const auto &item = cases[index];
    ids_are_contiguous =
        ids_are_contiguous && item.id == static_cast<int>(index);
    composite_keys_are_unique =
        composite_keys
            .insert(std::make_tuple(
                item.nlev, item.nonzero_aerosol_pairs, item.batch_size,
                item.mixing_index, item.mixing_scale, item.dt,
                static_cast<int>(item.method), item.forced_fe_substeps))
            .second &&
        composite_keys_are_unique;
    ++level_counts[item.nlev];
    ++occupancy_counts[item.nonzero_aerosol_pairs];
    occupancies_are_approved =
        occupancies_are_approved &&
        (item.nonzero_aerosol_pairs == 1 ||
         item.nonzero_aerosol_pairs == 4 ||
         item.nonzero_aerosol_pairs == 8 ||
         item.nonzero_aerosol_pairs == 16 ||
         item.nonzero_aerosol_pairs == 25);
    ++batch_counts[item.batch_size];
    ++mixing_counts[item.mixing_index];
    ++interval_counts[item.dt];
    switch (item.method) {
    case ndrop_test::benchmark::Method::forward_euler:
      ++fe;
      ++fe_substep_counts[item.forced_fe_substeps];
      break;
    case ndrop_test::benchmark::Method::imex:
      ++imex;
      break;
    case ndrop_test::benchmark::Method::backward_euler:
      ++be;
      break;
    }
  }
  check(ids_are_contiguous, "required benchmark case ids are contiguous");
  check(composite_keys_are_unique && composite_keys.size() == cases.size(),
        "required benchmark composite keys are unique");
  check(fe == 18900 && imex == 2700 && be == 2700,
        "required benchmark method coverage is exact");
  const auto &first = cases.front();
  const auto &last = cases.back();
  check(first.nlev == 1 && first.nonzero_aerosol_pairs == 1 &&
            first.batch_size == 1 &&
            first.mixing_index == 0 && first.mixing_scale == 0.05 &&
            first.dt == 33.33 && first.forced_fe_substeps == 1,
        "required benchmark starts at the exact lower matrix corner");
  check(last.nlev == 128 && last.nonzero_aerosol_pairs == 25 &&
            last.batch_size == 1024 && last.mixing_index == 1 &&
            last.mixing_scale == 0.5 && last.dt == 300 &&
            last.method == ndrop_test::benchmark::Method::backward_euler,
        "required benchmark ends at the exact upper matrix corner");

  constexpr std::array<int, 7> expected_substeps = {1, 2, 3, 5, 8, 16, 21};
  bool group_layout_is_exact = true;
  for (std::size_t group = 0; group < cases.size(); group += 9) {
    const auto &prototype = cases[group];
    for (std::size_t offset = 0; offset < 9; ++offset) {
      const auto &item = cases[group + offset];
      group_layout_is_exact =
          group_layout_is_exact && item.nlev == prototype.nlev &&
          item.nonzero_aerosol_pairs == prototype.nonzero_aerosol_pairs &&
          item.batch_size == prototype.batch_size &&
          item.mixing_index == prototype.mixing_index &&
          item.mixing_scale == prototype.mixing_scale &&
          item.dt == prototype.dt;
    }
    for (std::size_t offset = 0; offset < expected_substeps.size(); ++offset) {
      group_layout_is_exact =
          group_layout_is_exact &&
          cases[group + offset].method ==
              ndrop_test::benchmark::Method::forward_euler &&
          cases[group + offset].forced_fe_substeps ==
              expected_substeps[offset];
    }
    group_layout_is_exact =
        group_layout_is_exact &&
        cases[group + 7].method == ndrop_test::benchmark::Method::imex &&
        cases[group + 8].method ==
            ndrop_test::benchmark::Method::backward_euler;
  }
  check(group_layout_is_exact,
        "every benchmark group has the exact interleaved method layout");

  const std::map<std::size_t, std::size_t> expected_level_counts = {
      {1, 2700},   {2, 2700},   {3, 2700},   {16, 2700}, {32, 2700},
      {64, 2700},  {72, 2700},  {127, 2700}, {128, 2700}};
  const std::map<int, std::size_t> expected_occupancy_counts = {
      {1, 4860}, {4, 4860}, {8, 4860}, {16, 4860}, {25, 4860}};
  const std::map<int, std::size_t> expected_batch_counts = {
      {1, 4860}, {16, 4860}, {64, 4860}, {256, 4860}, {1024, 4860}};
  const std::map<int, std::size_t> expected_mixing_counts = {
      {0, 12150}, {1, 12150}};
  const std::map<Real, std::size_t> expected_interval_counts = {
      {33.33, 4050}, {50, 4050},  {75, 4050},
      {100, 4050},   {150, 4050}, {300, 4050}};
  const std::map<int, std::size_t> expected_fe_substep_counts = {
      {1, 2700}, {2, 2700}, {3, 2700}, {5, 2700},
      {8, 2700}, {16, 2700}, {21, 2700}};
  check(level_counts == expected_level_counts &&
            occupancy_counts == expected_occupancy_counts &&
            batch_counts == expected_batch_counts &&
            mixing_counts == expected_mixing_counts &&
            interval_counts == expected_interval_counts &&
            fe_substep_counts == expected_fe_substep_counts,
        "required benchmark has exact marginal matrix coverage");
  check(ndrop_test::benchmark::executed_aerosol_pairs == 25 &&
            ndrop_test::benchmark::executed_transport_fields == 51,
        "benchmark executed width is fixed at 25 pairs and 51 fields");
  check(occupancies_are_approved,
        "benchmark occupancy uses only the approved values");

  const auto smoke = ndrop_test::benchmark::required_cases(true);
  check(smoke.size() == 3, "smoke benchmark has one row per method");
  check(smoke[0].forced_fe_substeps == 3,
        "smoke benchmark exercises three FE substeps");
}

void test_statistics_and_accounting() {
  check(ndrop_test::benchmark::sample_median({1, 4, 2, 3}) == 2.5,
        "even sample median is exact");
  check(ndrop_test::benchmark::sample_median({3, 1, 2}) == 2,
        "odd sample median is exact");
  const auto interval = ndrop_test::benchmark::bootstrap_median_interval(
      std::vector<Real>(30, 2), 10000);
  check(interval.median == 2 && interval.lower == 2 && interval.upper == 2 &&
            interval.relative_half_width == 0,
        "constant bootstrap fixture has an exact zero-width interval");
  const auto repeated = ndrop_test::benchmark::bootstrap_median_interval(
      std::vector<Real>(30, 2), 10000);
  check(interval.median == repeated.median && interval.lower == repeated.lower &&
            interval.upper == repeated.upper,
        "bootstrap interval is deterministic");

  bool rejected = false;
  try {
    static_cast<void>(
        ndrop_test::benchmark::bootstrap_median_interval({}, 10000));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected, "bootstrap rejects an empty sample");

  const State state = make_state();
  const auto checksum = ndrop_test::benchmark::state_checksum(state);
  State changed = state;
  changed.droplet_number[0] += 1;
  check(checksum == ndrop_test::benchmark::state_checksum(state) &&
            checksum != ndrop_test::benchmark::state_checksum(changed),
        "state checksum is deterministic and state-sensitive");

  const auto cases = ndrop_test::benchmark::required_cases(true);
  check(ndrop_test::benchmark::estimated_workspace_bytes(cases[0]) > 0 &&
            ndrop_test::benchmark::estimated_workspace_bytes(cases[1]) > 0 &&
            ndrop_test::benchmark::estimated_workspace_bytes(cases[2]) > 0,
        "workspace estimates are positive for every method");

  State richardson_initial(2);
  richardson_initial.droplet_number = {2, 4};
  State richardson_refined = richardson_initial;
  State richardson_coarse = richardson_initial;
  richardson_coarse.droplet_number = {3, 2};
  const auto richardson = ndrop_test::benchmark::richardson_errors(
      richardson_coarse, richardson_refined, richardson_initial, {1, 3});
  check(richardson.weighted_l1_by_field[0] == 0.5 &&
            richardson.linf_by_field[0] == 0.5 &&
            richardson.maximum_weighted_l1 == 0.5 &&
            richardson.maximum_linf == 0.5 &&
            richardson.screening_error == 0.5,
        "Richardson errors match the hand-derived weighted and maximum norms");
  check(richardson.weighted_l1_by_field[1] == 0 &&
            richardson.linf_by_field[50] == 0,
        "Richardson field ordering retains zero unchanged aerosol fields");
  const auto reference = ndrop_test::benchmark::reference_errors(
      richardson_coarse, richardson_refined, richardson_initial, {1, 3});
  check(reference.screening_error == 0.5 &&
            reference.weighted_l1_by_field[0] == 0.5,
        "candidate-reference errors use the independent reference denominator");

  bool bad_richardson_rejected = false;
  try {
    static_cast<void>(ndrop_test::benchmark::richardson_errors(
        richardson_coarse, richardson_refined, richardson_initial, {1}));
  } catch (const std::invalid_argument &) {
    bad_richardson_rejected = true;
  }
  check(bad_richardson_rejected,
        "Richardson errors reject dimensionally inconsistent weights");

  AdvanceResult validated;
  AdvanceResult lean;
  validated.state = make_state();
  lean.state = validated.state;
  validated.ok = false;
  lean.ok = true;
  validated.diagnostics.clipping_count = 7;
  lean.diagnostics.clipping_count = 7;
  validated.diagnostics.material_negative_count = 7;
  lean.diagnostics.material_negative_count = 7;
  check(ndrop_test::benchmark::timing_path_matches_validated(validated,
                                                              lean),
        "a diagnostics-only failure may match an identical lean state");
  validated.ok = true;
  lean.ok = false;
  check(!ndrop_test::benchmark::timing_path_matches_validated(validated,
                                                               lean),
        "a lean failure cannot replace a passing validated path");
  validated.ok = false;
  lean.ok = true;
  lean.state.droplet_number[0] += 1;
  check(!ndrop_test::benchmark::timing_path_matches_validated(validated,
                                                               lean),
        "a state mismatch never passes timing-path equivalence");
}

void test_implicit_intervals() {
  const State initial = make_state();
  const OperatorData op = make_operator();
  constexpr Real total_dt = 0.5;
  for (const auto &method : {ndrop_test::benchmark::Method::imex,
                             ndrop_test::benchmark::Method::backward_euler}) {
    const std::string name =
        ndrop_test::benchmark::method_name(method);
    const auto single = [&](const State &state, const Real dt,
                            const AdvanceControl &control) {
      return method == ndrop_test::benchmark::Method::imex
                 ? ndrop_test::advance_imex(state, op, dt, control)
                 : ndrop_test::advance_backward_euler(state, op, dt, control);
    };
    const auto interval = [&](const int substeps,
                              const AdvanceControl &control) {
      return method == ndrop_test::benchmark::Method::imex
                 ? ndrop_test::advance_imex_interval(initial, op, total_dt,
                                                     substeps, control)
                 : ndrop_test::advance_backward_euler_interval(
                       initial, op, total_dt, substeps, control);
    };
    const auto existing = single(initial, total_dt, AdvanceControl{});
    const auto one = interval(1, AdvanceControl{});
    check(existing.ok && one.ok &&
              state_bitwise_equal(existing.state, one.state) &&
              existing.message == one.message &&
              existing.diagnostics.block_rhs_solve_count ==
                  one.diagnostics.block_rhs_solve_count,
          name + " interval N=1 preserves the existing solve");

    for (const int substeps : {2, 3}) {
      AdvanceControl detailed_control;
      detailed_control.collect_replay_diagnostics = true;
      const auto completed = interval(substeps, detailed_control);
      AdvanceControl lean;
      lean.collect_detailed_diagnostics = false;
      const auto timed_path = interval(substeps, lean);
      check(completed.ok && timed_path.ok &&
                state_bitwise_equal(completed.state, timed_path.state) &&
                completed.diagnostics.nsubmix == substeps,
            name + " interval diagnostics do not change state");
      const auto &detail = completed.diagnostics.replay_detail;
      check(detail && detail->expected_substeps ==
                          static_cast<std::uint32_t>(substeps) &&
                detail->completed_substeps ==
                    static_cast<std::uint32_t>(substeps) &&
                detail->physical_residual[0] &&
                detail->physical_residual[0]->substep <
                    static_cast<std::uint32_t>(substeps),
            name + " interval retains complete substep diagnostics");

      State sequential = initial;
      bool sequential_ok = true;
      for (int step = 0; step < substeps; ++step) {
        const auto next = single(sequential, total_dt / substeps,
                                 AdvanceControl{});
        sequential_ok = sequential_ok && next.ok;
        sequential = next.state;
      }
      // This fixture has positive cloud at every level; its projection is a
      // no-op, so independent one-step composition is a numerical oracle.
      check(sequential_ok && state_bitwise_equal(completed.state, sequential),
            name + " interval matches independently composed steps");
    }

    const auto invalid = interval(0, AdvanceControl{});
    check(!invalid.ok && invalid.message == "substeps must be positive",
          name + " rejects a zero internal count");

    OperatorData masked = op;
    masked.cloud_fraction[1] = 0;
    const auto with_final_projection =
        method == ndrop_test::benchmark::Method::imex
            ? ndrop_test::advance_imex_interval(initial, masked, total_dt, 2)
            : ndrop_test::advance_backward_euler_interval(initial, masked,
                                                           total_dt, 2);
    const auto first_projected =
        method == ndrop_test::benchmark::Method::imex
            ? ndrop_test::advance_imex(initial, masked, total_dt / 2)
            : ndrop_test::advance_backward_euler(initial, masked,
                                                 total_dt / 2);
    check(with_final_projection.ok && first_projected.ok &&
              first_projected.state.cloudborne[0][1] == 0,
          name + " zero-cloud fixture reaches the post-integration reset");
    if (with_final_projection.ok && first_projected.ok) {
      const auto twice_projected =
          method == ndrop_test::benchmark::Method::imex
              ? ndrop_test::advance_imex(first_projected.state, masked,
                                         total_dt / 2)
              : ndrop_test::advance_backward_euler(
                    first_projected.state, masked, total_dt / 2);
      check(twice_projected.ok &&
                !state_bitwise_equal(with_final_projection.state,
                                     twice_projected.state),
            name + " interval does not project at an internal step");
    }
  }
}

void test_controller_matched_forward_euler() {
  const State initial = make_state();
  OperatorData op(initial.nlev());
  op.cloud_fraction = {0.8, 0.5, 0.7};
  op.first_level_controller_rate = 4.5;
  constexpr Real interval = 0.5;
  const auto prepared = ndrop_test::prepare_operator(op);
  const int controller_count =
      ndrop_test::production_nsubmix(prepared, interval);
  check(controller_count == 3,
        "production controller supplies an exact non-dyadic FE count");

  const auto adaptive = ndrop_test::advance_forward_euler(initial, op, interval);
  ndrop_test::ForwardEulerOptions forced_options;
  forced_options.forced_substeps = controller_count;
  const auto forced = ndrop_test::advance_forward_euler(
      initial, op, interval, forced_options);
  check(adaptive.ok && forced.ok &&
            adaptive.diagnostics.nsubmix == controller_count &&
            forced.diagnostics.nsubmix == controller_count &&
            state_bitwise_equal(adaptive.state, forced.state),
        "adaptive FE executes the exact controller count without rounding");
}

} // namespace

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  test_method(
      "Forward Euler", ndrop_test::benchmark::Method::forward_euler,
      [](const State &state, const OperatorData &op, const Real dt,
         const AdvanceControl &control) {
        ndrop_test::ForwardEulerOptions options;
        options.forced_substeps = 3;
        return ndrop_test::advance_forward_euler(state, op, dt, options,
                                                 control);
      });
  test_method("IMEX", ndrop_test::benchmark::Method::imex,
              ndrop_test::advance_imex);
  test_method("Backward Euler",
              ndrop_test::benchmark::Method::backward_euler,
              ndrop_test::advance_backward_euler);
  test_required_matrix();
  test_statistics_and_accounting();
  test_implicit_intervals();
  test_controller_matched_forward_euler();
  Kokkos::finalize();
  std::cout << "RESULT passed=" << passed << " failed=" << failed << '\n';
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
