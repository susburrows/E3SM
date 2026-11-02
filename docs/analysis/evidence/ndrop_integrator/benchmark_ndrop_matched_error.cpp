#include "ndrop_benchmark.hpp"
#include "ndrop_integrator.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using ndrop_test::AdvanceControl;
using ndrop_test::AdvanceResult;
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::State;
using ndrop_test::benchmark::Method;

struct Fixture {
  int id;
  std::size_t levels;
  Real mixing;
  Real interval;
};

struct ControllerFixture {
  Fixture fixture;
  const char *family;
  int requested_nsubmix;
  Real activation_mixing = 0; // Zero preserves the legacy generator verbatim.
};

constexpr std::array<Fixture, 4> fixtures = {{{0, 72, 0.05, 50},
                                               {1, 72, 0.5, 150},
                                               {2, 128, 0.05, 50},
                                               {3, 128, 0.5, 150}}};
constexpr std::array<Method, 3> methods = {
    Method::forward_euler, Method::imex, Method::backward_euler};

OperatorData operator_data(const Fixture &fixture);

std::vector<ControllerFixture> controller_fixtures() {
  std::vector<ControllerFixture> result;
  result.reserve(24);
  for (const auto &fixture : fixtures) {
    result.push_back({fixture, "continuity", 0});
  }
  constexpr std::array<int, 10> requested = {
      1, 2, 4, 8, 16, 32, 64, 128, 256, 512};
  int id = 4;
  for (const int target : requested) {
    for (const std::size_t levels : {std::size_t(72), std::size_t(128)}) {
      constexpr Real interval = 150;
      const Fixture unit{id, levels, 1, interval};
      const auto prepared = ndrop_test::prepare_operator(operator_data(unit));
      Real inverse_rate = 0;
      for (std::size_t level = 0; level < levels; ++level) {
        inverse_rate = std::max(
            inverse_rate,
            prepared.rate_below[level] + prepared.rate_above[level] +
                (level == 0 ? prepared.first_level_controller_rate : Real(0)));
      }
      if (!(inverse_rate > 0) || !std::isfinite(inverse_rate)) {
        throw std::runtime_error("unit synthetic operator has invalid rate");
      }
      const Real mixing =
          Real(0.9) * (static_cast<Real>(target) - Real(0.5)) /
          (interval * inverse_rate);
      ControllerFixture generated{{id++, levels, mixing, interval},
                                  "controller_stratified", target};
      const auto actual = ndrop_test::production_nsubmix(
          ndrop_test::prepare_operator(operator_data(generated.fixture)),
          interval);
      if (actual != target) {
        throw std::runtime_error(
            "controller-stratified fixture did not reach requested count");
      }
      result.push_back(generated);
    }
  }
  return result;
}

int parse_positive(const std::string &name, const std::string &value) {
  std::size_t used = 0;
  const long number = std::stol(value, &used);
  if (used != value.size() || number <= 0 ||
      number > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<int>(number);
}

std::string field_array(const std::array<Real, 51> &values) {
  std::ostringstream out;
  out << std::setprecision(17);
  for (std::size_t field = 0; field < values.size(); ++field) {
    if (field) out << ';';
    out << values[field];
  }
  return out.str();
}

std::string csv_text(const std::string &value) {
  if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
  std::string quoted = "\"";
  for (const char character : value) {
    if (character == '\"') quoted += '\"';
    quoted += character;
  }
  return quoted + '\"';
}

State initial_state(const Fixture &fixture) {
  State state(fixture.levels);
  for (std::size_t level = 0; level < fixture.levels; ++level) {
    state.droplet_number[level] = 4 + 0.01 * level;
    for (int field = 0; field < ndrop_test::num_constituents; ++field) {
      state.interstitial[field][level] =
          1 + 0.01 * field + 0.002 * level;
      state.cloudborne[field][level] =
          0.2 + 0.002 * field + 0.001 * level;
    }
  }
  return state;
}

OperatorData operator_data(const Fixture &fixture) {
  OperatorData op(fixture.levels);
  for (std::size_t level = 0; level < fixture.levels; ++level) {
    op.pressure_weight[level] =
        0.8 + 0.4 * static_cast<Real>((17 * level + 3) % 11) / 10;
    op.cloud_fraction[level] = level % 2 == 0 ? 0.8 : 0.55;
  }
  for (std::size_t face = 1; face < fixture.levels; ++face) {
    const int offset = static_cast<int>((7 * face + 1) % 5) - 2;
    op.interface_exchange[face] =
        fixture.mixing * (1 + 0.05 * offset) *
        std::min(op.pressure_weight[face - 1], op.pressure_weight[face]);
  }
  constexpr std::array<Real, ndrop_test::num_modes> number = {0.35, 0.4,
                                                                0.45, 0.5};
  constexpr std::array<Real, ndrop_test::num_modes> mass = {0.25, 0.3,
                                                              0.35, 0.4};
  for (int mode = 0; mode < ndrop_test::num_modes; ++mode) {
    for (std::size_t level = 0; level + 1 < fixture.levels; ++level) {
      const Real below = op.interface_exchange[level + 1] /
                         op.pressure_weight[level];
      op.number_activation_rate[mode][level] = number[mode] * below;
      op.mass_activation_rate[mode][level] = mass[mode] * below;
    }
  }
  op.first_level_controller_rate = 0.25 * fixture.mixing;
  return op;
}

std::vector<ControllerFixture> study_fixtures(const bool controlled) {
  if (!controlled) return controller_fixtures();
  return {{{24, 72, 0.05, 150}, "controlled_mixing", 0, 0.05},
          {{25, 72, 0.5, 150}, "controlled_mixing", 0, 0.05},
          {{26, 128, 0.05, 150}, "controlled_mixing", 0, 0.05},
          {{27, 128, 0.5, 150}, "controlled_mixing", 0, 0.05}};
}

OperatorData study_operator(const ControllerFixture &fixture) {
  auto op = operator_data(fixture.fixture);
  if (fixture.activation_mixing == 0) return op;
  auto baseline = fixture.fixture;
  baseline.mixing = fixture.activation_mixing;
  const auto activation = operator_data(baseline);
  op.number_activation_rate = activation.number_activation_rate;
  op.mass_activation_rate = activation.mass_activation_rate;
  const auto prepared = ndrop_test::prepare_operator(op);
  if (prepared.activation_limiter_count != 0 ||
      prepared.number_activation_rate != activation.number_activation_rate ||
      prepared.mass_activation_rate != activation.mass_activation_rate) {
    throw std::runtime_error("controlled mixing changed fixed activation inputs");
  }
  return op;
}

// Physical-column fixture equations/units are frozen in the design document's
// physical-column extension. This factory never changes an integrator operator.
struct PhysicalColumn {
  Fixture fixture;
  const char *turbulence;
  Real velocity;
  State state;
  OperatorData op;
  std::vector<Real> height, thickness, diffusivity;
  int cloud_base = -1;
  int limiting_level = -1;
  Real controller_rate_t = 0;
  int parent_grid_levels = 0;
  int subdivision = 1;
  bool geometry_study = false;
  std::vector<Real> interface_pressure{}, midpoint_pressure{};
  std::vector<std::size_t> parent_cell{};
};

std::vector<Real> read_hybrid_pressure(const std::string &path,
                                       const std::size_t count) {
  std::ifstream input(path);
  if (!input) throw std::invalid_argument("cannot read physical grid: " + path);
  const auto next_line = [&]() {
    std::string line;
    while (std::getline(input, line)) {
      line = line.substr(0, line.find('!'));
      if (line.find_first_not_of(" \t\r") != std::string::npos) return line;
    }
    throw std::invalid_argument("incomplete physical grid: " + path);
  };
  const auto one_value = [&](const std::string &line) {
    std::istringstream values(line);
    Real value;
    std::string extra;
    if (!(values >> value) || values >> extra || !std::isfinite(value))
      throw std::invalid_argument("invalid physical grid value: " + path);
    return value;
  };
  std::vector<Real> a(count), pressure(count);
  for (int block = 0; block < 2; ++block) {
    if (one_value(next_line()) != static_cast<Real>(count))
      throw std::invalid_argument("physical grid extent mismatch: " + path);
    for (std::size_t index = 0; index < count; ++index) {
      const Real value = one_value(next_line());
      if (value < 0 || value > 1)
        throw std::invalid_argument("physical grid coefficient outside [0,1]");
      if (block == 0) a[index] = value;
      else pressure[index] = a[index] * 100000 + value * 101780;
    }
  }
  std::string extra;
  while (std::getline(input, extra)) {
    extra = extra.substr(0, extra.find('!'));
    if (extra.find_first_not_of(" \t\r") != std::string::npos)
      throw std::invalid_argument("extra physical grid content: " + path);
  }
  for (std::size_t index = 0; index < count; ++index) {
    if (!(pressure[index] > 0) || !std::isfinite(pressure[index]) ||
        (index && !(pressure[index] > pressure[index - 1])))
      throw std::invalid_argument("physical grid pressures not positive increasing");
  }
  return pressure;
}

// Same face assembly as the approved native fixtures. Refinement changes only
// fixture inputs; all preparation, controller and advance entry points stay put.
void finish_physical_column(PhysicalColumn &column,
                            const std::vector<Real> &pint,
                            const std::vector<Real> &pmid) {
  constexpr Real rd = 287, temperature = 289, gravity = 9.80665;
  constexpr Real scale_height = rd * temperature / gravity;
  constexpr std::array<Real, 4> fn = {.35, .4, .45, .5};
  constexpr std::array<Real, 4> fm = {.25, .3, .35, .4};
  const Real pi = std::acos(Real(-1));
  for (std::size_t face = 1; face < column.fixture.levels; ++face) {
    const std::size_t k = face - 1;
    const Real spacing = column.height[k] - column.height[k + 1];
    if (!(spacing > 0) || !std::isfinite(spacing))
      throw std::invalid_argument("physical midpoint heights not decreasing");
    const Real rho_face = pint[face] / (rd * temperature);
    const Real delta_cloud = column.op.cloud_fraction[k] -
                             column.op.cloud_fraction[k + 1];
    Real diffusivity = scale_height * std::log(101780 / pint[face]) < 840
                           ? column.fixture.mixing : Real(.01);
    if (column.op.cloud_fraction[k] > .01 && delta_cloud > .01) {
      diffusivity = column.velocity * spacing;
      column.cloud_base = static_cast<int>(k);
      const Real crdz = rho_face / (pmid[k] / (rd * temperature)) /
                        column.thickness[k];
      for (int mode = 0; mode < ndrop_test::num_modes; ++mode) {
        column.op.number_activation_rate[mode][k] =
            fn[mode] * column.velocity / std::sqrt(2 * pi) * crdz * delta_cloud;
        column.op.mass_activation_rate[mode][k] =
            fm[mode] * column.velocity / std::sqrt(2 * pi) * crdz * delta_cloud;
      }
    }
    column.diffusivity[face] = diffusivity;
    column.op.interface_exchange[face] = rho_face * diffusivity / spacing;
  }
  column.op.first_level_controller_rate =
      column.op.interface_exchange[1] / column.op.pressure_weight[0];
  std::string message;
  if (column.cloud_base < 0 ||
      !ndrop_test::validate(column.state, column.op, 150, message))
    throw std::invalid_argument("invalid physical fixture: " + message);
  const auto prepared = ndrop_test::prepare_operator(column.op);
  for (std::size_t k = 0; k < column.fixture.levels; ++k) {
    const Real rate_t = column.fixture.interval *
        (prepared.rate_above[k] + prepared.rate_below[k] +
         (k == 0 ? prepared.first_level_controller_rate : Real(0)));
    if (rate_t > column.controller_rate_t) {
      column.controller_rate_t = rate_t;
      column.limiting_level = static_cast<int>(k);
    }
  }
}

std::vector<PhysicalColumn> physical_columns(const std::string &grid_dir) {
  if (grid_dir.empty())
    throw std::invalid_argument("--physical-grid-dir is required");
  constexpr Real rd = 287, temperature = 289, gravity = 9.80665;
  constexpr Real scale_height = rd * temperature / gravity;
  constexpr std::array<Real, 4> number = {100e6, 200e6, 1e6, 10e6};
  constexpr std::array<Real, 4> diameter = {.12e-6, .04e-6, 1e-6, .10e-6};
  constexpr std::array<Real, 4> sigma = {1.8, 1.6, 1.8, 1.6};
  constexpr std::array<int, 4> species = {7, 4, 7, 3};
  const Real pi = std::acos(Real(-1));
  std::vector<PhysicalColumn> columns;
  for (const std::size_t levels : {std::size_t(72), std::size_t(128)}) {
    const std::string prefix = grid_dir +
        (levels == 72 ? "/acme-72" : "/scream-128");
    const auto pint = read_hybrid_pressure(prefix + "i.ascii", levels + 1);
    const auto pmid = read_hybrid_pressure(prefix + "m.ascii", levels);
    if (pint.back() != 101780)
      throw std::invalid_argument("physical grid surface pressure mismatch");
    for (const bool strong : {false, true}) {
      const int id = (levels == 72 ? 28 : 30) + static_cast<int>(strong);
      PhysicalColumn column{{id, levels, strong ? Real(30) : Real(3), 150},
                            strong ? "strong" : "moderate",
                            strong ? Real(1) : Real(.3), State(levels),
                            OperatorData(levels), {}, {}, {}, -1, -1, 0};
      column.height.resize(levels);
      column.thickness.resize(levels);
      column.diffusivity.assign(levels + 1, 0);
      column.parent_grid_levels = static_cast<int>(levels);
      column.interface_pressure = pint;
      column.midpoint_pressure = pmid;
      for (std::size_t k = 0; k < levels; ++k) {
        if (!(pint[k] < pmid[k] && pmid[k] < pint[k + 1]))
          throw std::invalid_argument("physical midpoint outside interfaces");
        const Real upper = scale_height * std::log(101780 / pint[k]);
        const Real lower = scale_height * std::log(101780 / pint[k + 1]);
        const Real dz = scale_height * std::log(pint[k + 1] / pint[k]);
        const Real rho = pmid[k] / (rd * temperature);
        column.height[k] = scale_height * std::log(101780 / pmid[k]);
        column.thickness[k] = dz;
        column.op.pressure_weight[k] = (pint[k + 1] - pint[k]) / gravity;
        column.op.cloud_fraction[k] = .8 *
            (std::max(Real(0), std::min(upper, Real(840)) -
                               std::max(lower, Real(600))) / (upper - lower));
        const Real cloud = column.op.cloud_fraction[k];
        column.state.droplet_number[k] = 80e6 * cloud / rho;
        for (int j = 0; j < ndrop_test::num_constituents; ++j) {
          const int mode = ndrop_test::constituent_mode(j);
          const Real log_sigma = std::log(sigma[mode]);
          const Real mass = number[mode] * (pi / 6) *
              std::pow(diameter[mode], 3) *
              std::exp(4.5 * log_sigma * log_sigma) * 1500;
          const Real total = (ndrop_test::constituent_is_number(j)
                                  ? number[mode] : mass / species[mode]) / rho;
          column.state.cloudborne[j][k] = total * .1 * cloud;
          column.state.interstitial[j][k] = total - column.state.cloudborne[j][k];
        }
      }
      finish_physical_column(column, pint, pmid);
      columns.push_back(std::move(column));
    }
  }
  return columns;
}

std::vector<PhysicalColumn> physical_geometry_columns(const std::string &grid_dir) {
  auto columns = physical_columns(grid_dir);
  const auto native = columns;
  for (auto &column : columns) column.geometry_study = true;
  constexpr Real scale_height = 287 * 289 / 9.80665;
  int id = 32;
  for (std::size_t pair = 0; pair < native.size(); pair += 2) {
    for (const int subdivision : {2, 4}) {
      for (std::size_t strength = 0; strength < 2; ++strength) {
        const auto &parent = native[pair + strength];
        std::vector<Real> pint{parent.interface_pressure.front()}, pmid;
        std::vector<std::size_t> parent_cell;
        for (std::size_t k = 0; k < parent.fixture.levels; ++k) {
          const Real upper = scale_height * std::log(101780 / parent.interface_pressure[k]);
          const Real lower = scale_height * std::log(101780 / parent.interface_pressure[k + 1]);
          const int children = upper > 400 && lower < 1000 ? subdivision : 1;
          for (int j = 1; j <= children; ++j) {
            const Real bottom = j == children ? parent.interface_pressure[k + 1]
                : parent.interface_pressure[k] * std::exp(
                    j * std::log(parent.interface_pressure[k + 1] /
                                 parent.interface_pressure[k]) / children);
            pmid.push_back(children == 1 ? parent.midpoint_pressure[k]
                                        : std::sqrt(pint.back() * bottom));
            pint.push_back(bottom);
            parent_cell.push_back(k);
          }
        }
        const auto levels = pmid.size();
        PhysicalColumn column{{id++, levels, parent.fixture.mixing, 150},
                              parent.turbulence, parent.velocity, State(levels),
                              OperatorData(levels), {}, {}, {}, -1, -1, 0};
        column.parent_grid_levels = parent.parent_grid_levels;
        column.subdivision = subdivision;
        column.geometry_study = true;
        column.interface_pressure = pint;
        column.midpoint_pressure = pmid;
        column.parent_cell = parent_cell;
        column.height.resize(levels);
        column.thickness.resize(levels);
        column.diffusivity.assign(levels + 1, 0);
        for (std::size_t k = 0; k < levels; ++k) {
          if (!(pint[k] < pmid[k] && pmid[k] < pint[k + 1]))
            throw std::invalid_argument("refined midpoint outside interfaces");
          const Real upper = scale_height * std::log(101780 / pint[k]);
          const Real lower = scale_height * std::log(101780 / pint[k + 1]);
          column.height[k] = scale_height * std::log(101780 / pmid[k]);
          column.thickness[k] = scale_height * std::log(pint[k + 1] / pint[k]);
          column.op.pressure_weight[k] = (pint[k + 1] - pint[k]) / 9.80665;
          column.op.cloud_fraction[k] = .8 *
              (std::max(Real(0), std::min(upper, Real(840)) -
                                 std::max(lower, Real(600))) / (upper - lower));
        }
        // Design: conservative per-parent remap, not new concentrations.
        for (std::size_t first = 0; first < levels;) {
          std::size_t end = first + 1;
          while (end < levels && parent_cell[end] == parent_cell[first]) ++end;
          const auto k = parent_cell[first];
          Real mass = 0, cloudy_mass = 0;
          for (auto child = first; child < end; ++child) {
            mass += column.op.pressure_weight[child];
            cloudy_mass += column.op.pressure_weight[child] * column.op.cloud_fraction[child];
          }
          if (!(mass > 0) || !std::isfinite(mass) || !std::isfinite(cloudy_mass))
            throw std::invalid_argument("invalid refined mass");
          if (cloudy_mass == 0) {
            if (parent.state.droplet_number[k] != 0)
              throw std::invalid_argument("cloudless remap lost parent droplets");
            for (const auto &field : parent.state.cloudborne)
              if (field[k] != 0) throw std::invalid_argument("cloudless remap lost parent aerosol");
          }
          for (auto child = first; child < end; ++child) {
            const Real factor = end - first == 1 ? Real(1)
                : cloudy_mass == 0 ? Real(0) : column.op.cloud_fraction[child] / (cloudy_mass / mass);
            column.state.droplet_number[child] = parent.state.droplet_number[k] * factor;
            for (int j = 0; j < ndrop_test::num_constituents; ++j) {
              column.state.interstitial[j][child] = parent.state.interstitial[j][k];
              column.state.cloudborne[j][child] = parent.state.cloudborne[j][k] * factor;
            }
          }
          first = end;
        }
        finish_physical_column(column, pint, pmid);
        columns.push_back(std::move(column));
      }
    }
  }
  return columns;
}

AdvanceResult advance(const Method method, const State &state,
                      const OperatorData &op, const Real interval,
                      const int steps, const AdvanceControl &control) {
  switch (method) {
  case Method::forward_euler: {
    ndrop_test::ForwardEulerOptions options;
    options.forced_substeps = steps;
    return ndrop_test::advance_forward_euler(state, op, interval, options,
                                             control);
  }
  case Method::imex:
    return ndrop_test::advance_imex_interval(state, op, interval, steps,
                                             control);
  case Method::backward_euler:
    return ndrop_test::advance_backward_euler_interval(
        state, op, interval, steps, control);
  }
  throw std::logic_error("unknown method");
}

AdvanceResult advance_controller_fe(const State &state, const OperatorData &op,
                                    const Real interval,
                                    const AdvanceControl &control) {
  return ndrop_test::advance_forward_euler(
      state, op, interval, ndrop_test::ForwardEulerOptions{}, control);
}

void accuracy() {
  std::cout << std::setprecision(17)
            << "artifact_schema,fixture,levels,mixing_scale,interval_seconds,"
               "mixing_rate_t,activation_rate_t,"
               "method,steps,reference_status,reference_reason,"
               "reference_disagreement,status,failure_reason,"
               "error_normal,error_tight,error_max,"
               "normal_weighted_l1_by_field,normal_linf_by_field,"
               "tight_weighted_l1_by_field,tight_linf_by_field,"
               "richardson_to_next,observed_order,"
               "clipping_count,material_negative_count,normalized_residual,"
               "normalized_ledger_error,state_checksum\n";
  for (const Fixture &fixture : fixtures) {
    const State state = initial_state(fixture);
    const OperatorData op = operator_data(fixture);
    ndrop_test::ReferenceOptions tight_options;
    tight_options.relative_tolerance = 1e-13;
    const auto normal = ndrop_test::advance_reference_dopri54(
        state, op, fixture.interval);
    const auto tight = ndrop_test::advance_reference_dopri54(
        state, op, fixture.interval, tight_options);
    Real disagreement = std::numeric_limits<Real>::infinity();
    if (normal.ok && tight.ok) {
      disagreement = ndrop_test::benchmark::reference_errors(
                         normal.state, tight.state, state, op.pressure_weight)
                         .screening_error;
    }
    const bool reference_ok = normal.ok && tight.ok &&
                              std::isfinite(disagreement) &&
                              disagreement <= 2e-11;
    const std::string reference_reason =
        reference_ok ? "" : !normal.ok ? "normal reference: " + normal.message
                           : !tight.ok ? "tight reference: " + tight.message
                           : "reference disagreement exceeds tolerance";
    const auto prepared = ndrop_test::prepare_operator(op);
    Real mixing_rate_t = 0;
    Real activation_rate_t = 0;
    for (std::size_t level = 0; level < fixture.levels; ++level) {
      mixing_rate_t = std::max(
          mixing_rate_t,
          fixture.interval * (prepared.rate_below[level] +
                              prepared.rate_above[level]));
      for (int mode = 0; mode < ndrop_test::num_modes; ++mode) {
        activation_rate_t = std::max(
            {activation_rate_t,
             fixture.interval * prepared.number_activation_rate[mode][level],
             fixture.interval * prepared.mass_activation_rate[mode][level]});
      }
    }
    const int fe_controller =
        ndrop_test::production_nsubmix(prepared, fixture.interval);
    for (const Method method : methods) {
      std::vector<int> step_counts;
      std::vector<AdvanceResult> candidates;
      for (int steps = 1; steps <= 8192; steps *= 2) {
        if (method == Method::forward_euler && steps < fe_controller) continue;
        step_counts.push_back(steps);
        candidates.push_back(advance(method, state, op, fixture.interval,
                                     steps, AdvanceControl{}));
      }
      const auto eligible_candidate = [&](const std::size_t index) {
        const auto &candidate = candidates[index];
        return reference_ok && candidate.ok &&
               candidate.diagnostics.clipping_count == 0 &&
               candidate.diagnostics.material_negative_count == 0;
      };
      std::vector<std::optional<Real>> richardson(candidates.size());
      for (std::size_t index = 0; index + 1 < candidates.size(); ++index) {
        if (eligible_candidate(index) && eligible_candidate(index + 1)) {
          richardson[index] = ndrop_test::benchmark::richardson_errors(
              candidates[index].state, candidates[index + 1].state,
              state, op.pressure_weight).screening_error;
        }
      }
      for (std::size_t index = 0; index < candidates.size(); ++index) {
        const int steps = step_counts[index];
        const auto &candidate = candidates[index];
        const bool eligible =
            eligible_candidate(index);
        ndrop_test::benchmark::ReferenceErrors error_normal;
        ndrop_test::benchmark::ReferenceErrors error_tight;
        if (eligible) {
          error_normal = ndrop_test::benchmark::reference_errors(
              candidate.state, normal.state, state, op.pressure_weight);
          error_tight = ndrop_test::benchmark::reference_errors(
              candidate.state, tight.state, state, op.pressure_weight);
        }
        const char *status = !reference_ok
                                 ? "reference_unresolved"
                                 : !candidate.ok ? "numerical_failure"
                                 : candidate.diagnostics.clipping_count
                                     ? "fe_clipping"
                                 : candidate.diagnostics.material_negative_count
                                     ? "material_negative"
                                     : "eligible";
        const std::string failure_reason =
            eligible ? "" : !reference_ok ? reference_reason
                                           : candidate.message;
        const Real maximum =
            std::max(error_normal.screening_error, error_tight.screening_error);
        std::cout << "ndrop_matched_accuracy_v1," << fixture.id << ','
                  << fixture.levels << ',' << fixture.mixing << ','
                  << fixture.interval << ',' << mixing_rate_t << ','
                  << activation_rate_t << ','
                  << ndrop_test::benchmark::method_name(method) << ',' << steps
                  << ',' << (reference_ok ? "certified" : "unresolved") << ','
                  << csv_text(reference_reason) << ',';
        if (std::isfinite(disagreement)) std::cout << disagreement;
        std::cout << ',' << status << ',' << csv_text(failure_reason) << ',';
        if (eligible) {
          std::cout << error_normal.screening_error << ','
                    << error_tight.screening_error << ',' << maximum << ','
                    << field_array(error_normal.weighted_l1_by_field) << ','
                    << field_array(error_normal.linf_by_field) << ','
                    << field_array(error_tight.weighted_l1_by_field) << ','
                    << field_array(error_tight.linf_by_field);
        } else {
          std::cout << ",,,,,,";
        }
        std::cout << ',';
        if (richardson[index]) std::cout << *richardson[index];
        std::cout << ',';
        if (index + 1 < richardson.size() && richardson[index] &&
            richardson[index + 1] && *richardson[index] > 0 &&
            *richardson[index + 1] > 0) {
          std::cout << std::log2(*richardson[index] / *richardson[index + 1]);
        }
        std::cout << ',' << candidate.diagnostics.clipping_count << ','
                  << candidate.diagnostics.material_negative_count << ','
                  << candidate.diagnostics.normalized_residual << ','
                  << candidate.diagnostics.max_normalized_ledger_error << ','
                  << ndrop_test::benchmark::state_checksum(candidate.state)
                  << '\n';
      }
    }
  }
}

struct CertifiedReferences {
  ndrop_test::ReferenceResult normal;
  ndrop_test::ReferenceResult tight;
  Real disagreement = std::numeric_limits<Real>::infinity();
  bool ok = false;
  std::string reason;
};

CertifiedReferences references(const State &state, const OperatorData &op,
                               const Real interval) {
  CertifiedReferences result;
  result.normal = ndrop_test::advance_reference_dopri54(state, op, interval);
  ndrop_test::ReferenceOptions tight_options;
  tight_options.relative_tolerance = 1e-13;
  result.tight = ndrop_test::advance_reference_dopri54(
      state, op, interval, tight_options);
  if (result.normal.ok && result.tight.ok) {
    result.disagreement = ndrop_test::benchmark::reference_errors(
                              result.normal.state, result.tight.state, state,
                              op.pressure_weight)
                              .screening_error;
  }
  result.ok = result.normal.ok && result.tight.ok &&
              std::isfinite(result.disagreement) &&
              result.disagreement <= 2e-11;
  result.reason =
      result.ok ? ""
                : !result.normal.ok
                      ? "normal reference: " + result.normal.message
                      : !result.tight.ok
                            ? "tight reference: " + result.tight.message
                            : "reference disagreement exceeds tolerance";
  return result;
}

struct CandidateEvaluation {
  AdvanceResult result;
  ndrop_test::benchmark::ReferenceErrors normal_error;
  ndrop_test::benchmark::ReferenceErrors tight_error;
  Real maximum_error = std::numeric_limits<Real>::quiet_NaN();
  bool eligible = false;
  std::string status;
  std::string reason;
};

CandidateEvaluation evaluate_candidate(
    const Method method, const State &state, const OperatorData &op,
    const Real interval, const int steps, const CertifiedReferences &reference) {
  CandidateEvaluation evaluated;
  evaluated.result =
      advance(method, state, op, interval, steps, AdvanceControl{});
  evaluated.eligible = reference.ok && evaluated.result.ok &&
                       evaluated.result.diagnostics.clipping_count == 0 &&
                       evaluated.result.diagnostics.material_negative_count == 0;
  if (!reference.ok) {
    evaluated.status = "reference_unresolved";
    evaluated.reason = reference.reason;
  } else if (!evaluated.result.ok) {
    evaluated.status = "numerical_failure";
    evaluated.reason = evaluated.result.message;
  } else if (evaluated.result.diagnostics.clipping_count != 0) {
    evaluated.status = "fe_clipping";
    evaluated.reason = "Forward Euler clipping count is nonzero";
  } else if (evaluated.result.diagnostics.material_negative_count != 0) {
    evaluated.status = "material_negative";
    evaluated.reason = "material-negative count is nonzero";
  } else {
    evaluated.normal_error = ndrop_test::benchmark::reference_errors(
        evaluated.result.state, reference.normal.state, state,
        op.pressure_weight);
    evaluated.tight_error = ndrop_test::benchmark::reference_errors(
        evaluated.result.state, reference.tight.state, state,
        op.pressure_weight);
    evaluated.maximum_error =
        std::max(evaluated.normal_error.screening_error,
                 evaluated.tight_error.screening_error);
    evaluated.status = "eligible";
  }
  return evaluated;
}

CandidateEvaluation evaluate_controller_fe(
    const State &state, const OperatorData &op, const Real interval,
    const CertifiedReferences &reference, const int expected_nsubmix) {
  CandidateEvaluation evaluated;
  evaluated.result =
      advance_controller_fe(state, op, interval, AdvanceControl{});
  evaluated.eligible = reference.ok && evaluated.result.ok &&
                       evaluated.result.diagnostics.nsubmix == expected_nsubmix &&
                       evaluated.result.diagnostics.clipping_count == 0 &&
                       evaluated.result.diagnostics.material_negative_count == 0;
  if (!reference.ok) {
    evaluated.status = "reference_unresolved";
    evaluated.reason = reference.reason;
  } else if (!evaluated.result.ok) {
    evaluated.status = "numerical_failure";
    evaluated.reason = evaluated.result.message;
  } else if (evaluated.result.diagnostics.nsubmix != expected_nsubmix) {
    evaluated.status = "fixture_controller_mismatch";
    evaluated.reason = "adaptive Forward Euler count differs from controller";
  } else if (evaluated.result.diagnostics.clipping_count != 0) {
    evaluated.status = "fe_clipping";
    evaluated.reason = "Forward Euler clipping count is nonzero";
  } else if (evaluated.result.diagnostics.material_negative_count != 0) {
    evaluated.status = "material_negative";
    evaluated.reason = "material-negative count is nonzero";
  } else {
    evaluated.normal_error = ndrop_test::benchmark::reference_errors(
        evaluated.result.state, reference.normal.state, state,
        op.pressure_weight);
    evaluated.tight_error = ndrop_test::benchmark::reference_errors(
        evaluated.result.state, reference.tight.state, state,
        op.pressure_weight);
    evaluated.maximum_error =
        std::max(evaluated.normal_error.screening_error,
                 evaluated.tight_error.screening_error);
    evaluated.status = "eligible";
  }
  return evaluated;
}

void physical_row(const PhysicalColumn &column, const char *kind,
                  const Method method, const int steps,
                  const CertifiedReferences &reference,
                  const std::string &status, const std::string &reason,
                  const CandidateEvaluation *candidate) {
  const auto prepared = ndrop_test::prepare_operator(column.op);
  Real mixing_rate_t = 0, activation_rate_t = 0;
  for (std::size_t k = 0; k < column.fixture.levels; ++k) {
    mixing_rate_t = std::max(mixing_rate_t, 150 *
                            (prepared.rate_above[k] + prepared.rate_below[k]));
    for (int mode = 0; mode < ndrop_test::num_modes; ++mode)
      activation_rate_t = std::max({activation_rate_t,
          150 * prepared.number_activation_rate[mode][k],
          150 * prepared.mass_activation_rate[mode][k]});
  }
  const int cb = column.cloud_base, limit = column.limiting_level;
  std::cout << "ndrop_physical_column_accuracy_v1," << kind << ','
            << column.fixture.id << ',' << column.fixture.levels << ','
            << column.turbulence << ',' << column.fixture.mixing << ','
            << column.velocity << ",150," << mixing_rate_t << ','
            << activation_rate_t << ',' << column.controller_rate_t << ','
            << cb << ',' << column.height[cb] - column.height[cb + 1] << ','
            << limit << ',' << column.thickness[limit] << ','
            << ndrop_test::benchmark::method_name(method) << ',' << steps << ','
            << ndrop_test::production_nsubmix(prepared, 150) << ','
            << (reference.ok ? "certified" : "unresolved") << ','
            << csv_text(reference.reason) << ',';
  if (std::isfinite(reference.disagreement)) std::cout << reference.disagreement;
  std::cout << ',' << status << ',' << csv_text(reason) << ',';
  if (candidate && candidate->eligible) {
    std::cout << candidate->normal_error.screening_error << ','
              << candidate->tight_error.screening_error << ','
              << candidate->maximum_error << ','
              << field_array(candidate->normal_error.weighted_l1_by_field) << ','
              << field_array(candidate->normal_error.linf_by_field) << ','
              << field_array(candidate->tight_error.weighted_l1_by_field) << ','
              << field_array(candidate->tight_error.linf_by_field);
  } else std::cout << ",,,,,,";
  std::cout << ',';
  if (candidate) {
    std::cout << candidate->result.diagnostics.clipping_count << ','
              << candidate->result.diagnostics.material_negative_count << ',';
    const auto &diagnostics = candidate->result.diagnostics;
    if (std::isfinite(diagnostics.normalized_residual))
      std::cout << diagnostics.normalized_residual;
    std::cout << ',';
    if (std::isfinite(diagnostics.max_normalized_ledger_error))
      std::cout << diagnostics.max_normalized_ledger_error;
    std::cout << ',' << ndrop_test::benchmark::state_checksum(candidate->result.state);
  } else std::cout << ",,,,";
  std::cout << ',' << prepared.activation_limiter_count << '\n';
}

void physical_accuracy(const std::vector<PhysicalColumn> &columns,
                       const int maximum_steps, const int budget_seconds) {
  const auto start = Clock::now();
  std::cout << std::setprecision(17)
      << "artifact_schema,record_kind,fixture,levels,turbulence,"
         "background_diffusivity_m2_s,cloud_base_wtke_m_s,interval_seconds,"
         "mixing_rate_t,activation_rate_t,controller_rate_t,cloud_base_level,"
         "cloud_base_spacing_m,limiting_level,limiting_spacing_m,method,steps,"
         "production_nsubmix,reference_status,reference_reason,"
         "reference_disagreement,status,failure_reason,error_normal,error_tight,"
         "error_max,normal_weighted_l1_by_field,normal_linf_by_field,"
         "tight_weighted_l1_by_field,tight_linf_by_field,clipping_count,"
         "material_negative_count,normalized_residual,normalized_ledger_error,"
         "state_checksum,activation_limiter_count\n";
  const auto exhausted = [&]() {
    return std::chrono::duration<Real>(Clock::now() - start).count() >= budget_seconds;
  };
  for (const auto &column : columns) {
    const int controller = ndrop_test::production_nsubmix(
        ndrop_test::prepare_operator(column.op), 150);
    const auto reference = references(column.state, column.op, 150);
    const auto baseline = evaluate_controller_fe(
        column.state, column.op, 150, reference, controller);
    physical_row(column, "baseline", Method::forward_euler, controller,
                 reference, baseline.status, baseline.reason, &baseline);
    for (const Method method : methods) {
      std::map<int, CandidateEvaluation> candidates;
      const auto evaluate = [&](const int steps) -> const CandidateEvaluation & {
        auto it = candidates.find(steps);
        if (it == candidates.end()) {
          auto candidate = method == Method::forward_euler && steps == controller
              ? baseline : evaluate_candidate(method, column.state, column.op,
                                               150, steps, reference);
          it = candidates.emplace(steps, std::move(candidate)).first;
        }
        return it->second;
      };
      int last = 0;
      std::string terminal = "not_attained_within_cap";
      std::string reason = "no eligible tight-target match within count cap";
      std::map<Real, int> selected;
      if (!reference.ok) {
        terminal = "reference_unresolved";
        reason = reference.reason;
      } else {
        for (int steps = method == Method::forward_euler ? controller : 1;
             steps <= maximum_steps; ++steps) {
          if (exhausted()) {
            terminal = "accuracy_budget_exhausted";
            reason = "untimed accuracy budget exhausted";
            break;
          }
          const auto &candidate = evaluate(steps);
          last = steps;
          physical_row(column, "measurement", method, steps, reference,
                       candidate.status, candidate.reason, &candidate);
          for (const Real target : {Real(.01), Real(.001)})
            if (!selected.count(target) && candidate.eligible &&
                candidate.maximum_error <= target) selected.emplace(target, steps);
          if (selected.count(.001)) {
            terminal = "selected";
            reason.clear();
            break;
          }
        }
      }
      physical_row(column, "terminal", method, last, reference,
                   terminal, reason, nullptr);
      std::set<int> confirmations;
      for (const auto &selection : selected)
        for (const int multiplier : {2, 4})
          confirmations.insert(selection.second * multiplier);
      for (const int steps : confirmations) {
        if (steps > maximum_steps || exhausted()) {
          physical_row(column, "confirmation", method, steps, reference,
                       steps > maximum_steps ? "confirmation_count_exceeds_cap"
                                             : "accuracy_budget_exhausted",
                       "confirmation unavailable within count/budget bounds", nullptr);
        } else {
          const auto &candidate = evaluate(steps);
          physical_row(column, "confirmation", method, steps, reference,
                       candidate.status, candidate.reason, &candidate);
        }
      }
    }
  }
}

void check_physical_fixtures(const std::vector<PhysicalColumn> &columns) {
  int checks = 0;
  const auto require = [&](const bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
  };
  for (std::size_t index = 0; index < columns.size(); ++index) {
    const auto &column = columns[index];
    const auto prepared = ndrop_test::prepare_operator(column.op);
    const int count = ndrop_test::production_nsubmix(prepared, 150);
    require(count > 0, "physical controller count invalid");
    require(column.op.interface_exchange.front() == 0 &&
            column.op.interface_exchange.back() == 0, "physical no-flux mismatch");
    require(column.cloud_base >= 0 && column.limiting_level >= 0,
            "physical cloud/limiting levels missing");
    if (index % 2) {
      const auto &moderate = columns[index - 1];
      const auto same_bits = [](const std::vector<Real> &a, const std::vector<Real> &b) {
        return a.size() == b.size() &&
            std::memcmp(a.data(), b.data(), a.size() * sizeof(Real)) == 0;
      };
      require(column.height == moderate.height &&
              column.op.pressure_weight == moderate.op.pressure_weight &&
              column.op.cloud_fraction == moderate.op.cloud_fraction &&
              column.state.droplet_number == moderate.state.droplet_number &&
              column.state.interstitial == moderate.state.interstitial &&
              column.state.cloudborne == moderate.state.cloudborne,
              "paired physical column inputs differ beyond turbulence");
      require(same_bits(column.state.droplet_number, moderate.state.droplet_number),
              "paired physical droplet bits differ");
      for (int j = 0; j < ndrop_test::num_constituents; ++j)
        require(same_bits(column.state.interstitial[j], moderate.state.interstitial[j]) &&
                same_bits(column.state.cloudborne[j], moderate.state.cloudborne[j]),
                "paired physical aerosol bits differ");
      require(count >= ndrop_test::production_nsubmix(
                  ndrop_test::prepare_operator(moderate.op), 150),
              "strong physical mixing has fewer controller steps");
    }
    for (std::size_t k = 0; k < column.fixture.levels; ++k) {
      require(column.thickness[k] > 0 && column.op.pressure_weight[k] > 0,
              "invalid physical geometry");
      if (k) {
        const Real mass_rate = column.op.pressure_weight[k] * prepared.rate_above[k];
        require(std::abs(mass_rate - column.op.interface_exchange[k]) <=
                    4 * std::numeric_limits<Real>::epsilon() *
                    column.op.interface_exchange[k], "dimensional exchange mismatch");
      }
      require(column.op.cloud_fraction[k] >= 0 && column.op.cloud_fraction[k] <= .8,
              "physical cloud fraction invalid");
      std::array<Real, 4> mode_number{}, mode_mass{};
      for (int j = 0; j < ndrop_test::num_constituents; ++j) {
        const Real total = column.state.interstitial[j][k] + column.state.cloudborne[j][k];
        (ndrop_test::constituent_is_number(j) ? mode_number : mode_mass)
            [ndrop_test::constituent_mode(j)] += total;
      }
      constexpr std::array<Real, 4> diameter = {.12e-6, .04e-6, 1e-6, .10e-6};
      constexpr std::array<Real, 4> sigma = {1.8, 1.6, 1.8, 1.6};
      constexpr std::array<Real, 4> number = {100e6, 200e6, 1e6, 10e6};
      const Real rho = 101780 * std::exp(-column.height[k] * 9.80665 / (287 * 289)) /
                       (287 * 289);
      for (int mode = 0; mode < 4; ++mode) {
        const Real mass_per_particle = std::acos(Real(-1)) / 6 * 1500 *
            std::exp(3 * std::log(diameter[mode]) +
                     4.5 * std::log(sigma[mode]) * std::log(sigma[mode]));
        if (column.subdivision == 1)
          require(std::abs(mode_number[mode] * rho / number[mode] - 1) <=
                      16 * std::numeric_limits<Real>::epsilon(),
                  "physical modal number/density mismatch");
        require(std::abs(mode_mass[mode] / mode_number[mode] / mass_per_particle - 1) <=
                    64 * std::numeric_limits<Real>::epsilon(),
                "physical modal mass/number lognormal moment mismatch");
      }
      if (column.op.cloud_fraction[k] == 0) {
        require(column.state.droplet_number[k] == 0, "clear-cell droplet not zero");
        for (int j = 0; j < ndrop_test::num_constituents; ++j)
          require(column.state.cloudborne[j][k] == 0, "clear-cell cloudborne not zero");
      }
    }
    for (const Method method : methods) {
      AdvanceControl disabled;
      disabled.collect_detailed_diagnostics = false;
      const int steps = method == Method::forward_euler ? count : 4 * count;
      const auto enabled = method == Method::forward_euler
          ? advance_controller_fe(column.state, column.op, 150, AdvanceControl{})
          : advance(method, column.state, column.op, 150, steps, AdvanceControl{});
      const auto uninstrumented = method == Method::forward_euler
          ? advance_controller_fe(column.state, column.op, 150, disabled)
          : advance(method, column.state, column.op, 150, steps, disabled);
      require(enabled.ok && uninstrumented.ok, "physical diagnostic advance failed");
      require(enabled.message == uninstrumented.message &&
              enabled.diagnostics.nsubmix == uninstrumented.diagnostics.nsubmix &&
              enabled.diagnostics.final_buffer == uninstrumented.diagnostics.final_buffer &&
              enabled.diagnostics.activation_limiter_count ==
                  uninstrumented.diagnostics.activation_limiter_count &&
              enabled.diagnostics.clipping_count == uninstrumented.diagnostics.clipping_count,
              "physical diagnostic decisions/counters changed");
      const auto same = [](const std::vector<Real> &a, const std::vector<Real> &b) {
        return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(Real)) == 0;
      };
      require(same(enabled.state.droplet_number, uninstrumented.state.droplet_number),
              "physical diagnostic droplet state changed");
      for (int j = 0; j < ndrop_test::num_constituents; ++j) {
        require(same(enabled.state.interstitial[j], uninstrumented.state.interstitial[j]) &&
                same(enabled.state.cloudborne[j], uninstrumented.state.cloudborne[j]),
                "physical diagnostic aerosol state changed");
      }
      if (method == Method::forward_euler)
        require(enabled.diagnostics.nsubmix == count, "physical FE bypassed controller");
      for (std::size_t k = 0; k < column.fixture.levels; ++k) {
        if (column.op.cloud_fraction[k] != 0) continue;
        require(enabled.state.droplet_number[k] == 0,
                "physical final clear-cell droplet reset missing");
        for (int j = 0; j < ndrop_test::num_constituents; ++j)
          require(enabled.state.cloudborne[j][k] == 0,
                  "physical final clear-cell resuspension missing");
      }
    }
    std::cout << "PHYSICAL fixture=" << column.fixture.id << " levels="
              << column.fixture.levels << " turbulence=" << column.turbulence
              << " controller_steps=" << count << " cloud_base_level="
              << column.cloud_base << " limiting_level=" << column.limiting_level
              << " max_rate_times_interval=" << column.controller_rate_t << '\n';
  }
  std::cout << "PHYSICAL_FIXTURES checks=" << checks << " failed=0\n";
}

void check_geometry_fixtures(const std::vector<PhysicalColumn> &columns,
                             const std::string &grid_dir) {
  check_physical_fixtures(columns);
  const auto native = physical_columns(grid_dir);
  int checks = 0;
  const auto require = [&](const bool condition, const char *message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
  };
  const auto close = [&](const Real a, const Real b) {
    return b == 0 ? a == 0 : std::abs(a / b - 1) <=
        128 * std::numeric_limits<Real>::epsilon();
  };
  const auto field = [](const State &state, const int f, const std::size_t k) {
    return f == 0 ? state.droplet_number[k] : f % 2
        ? state.interstitial[(f - 1) / 2][k] : state.cloudborne[(f - 2) / 2][k];
  };
  require(columns.size() == 12, "geometry fixture coverage mismatch");
  for (std::size_t index = 0; index < columns.size(); ++index) {
    const auto &column = columns[index];
    require(column.fixture.id == 28 + static_cast<int>(index), "geometry ID mismatch");
    const auto &parent = native[(column.parent_grid_levels == 72 ? 0 : 2) + index % 2];
    require(column.interface_pressure.front() == parent.interface_pressure.front() &&
            column.interface_pressure.back() == parent.interface_pressure.back(),
            "geometry full-column endpoints changed");
    if (column.subdivision == 1) {
      require(column.interface_pressure == parent.interface_pressure &&
              column.midpoint_pressure == parent.midpoint_pressure &&
              column.op.interface_exchange == parent.op.interface_exchange &&
              column.op.number_activation_rate == parent.op.number_activation_rate &&
              column.op.mass_activation_rate == parent.op.mass_activation_rate &&
              column.state.droplet_number == parent.state.droplet_number &&
              column.state.cloudborne == parent.state.cloudborne &&
              column.state.interstitial == parent.state.interstitial,
              "native fixture changed in geometry policy");
      continue;
    }
    const std::size_t expected = column.parent_grid_levels == 72
        ? column.subdivision == 2 ? 79 : 93
        : column.subdivision == 2 ? 141 : 167;
    require(column.fixture.levels == expected, "geometry actual level count mismatch");
    for (std::size_t first = 0; first < column.fixture.levels;) {
      std::size_t end = first + 1;
      while (end < column.fixture.levels && column.parent_cell[end] == column.parent_cell[first]) ++end;
      const auto k = column.parent_cell[first];
      const Real upper = 287 * 289 / 9.80665 * std::log(101780 / parent.interface_pressure[k]);
      const Real lower = 287 * 289 / 9.80665 * std::log(101780 / parent.interface_pressure[k + 1]);
      const auto children = upper > 400 && lower < 1000 ? column.subdivision : 1;
      require(end - first == static_cast<std::size_t>(children), "geometry band membership changed");
      require(column.interface_pressure[first] == parent.interface_pressure[k] &&
              column.interface_pressure[end] == parent.interface_pressure[k + 1],
              "parent interfaces changed");
      Real mass = 0;
      std::array<Real, 51> integrals{};
      for (auto child = first; child < end; ++child) {
        const Real w = column.op.pressure_weight[child];
        mass += w;
        // log(p_bottom/p_top) has cancellation on thin cells: bound the
        // absolute height error from rounded pressure endpoints, not a
        // fixed relative tolerance on the small difference.
        require(std::abs(column.thickness[child] * children - parent.thickness[k]) <=
                    8 * children * std::numeric_limits<Real>::epsilon() * (287 * 289 / 9.80665),
                "nonuniform log-pressure subdivision");
        if (children == 1)
          require(column.midpoint_pressure[child] == parent.midpoint_pressure[k],
                  "unsplit midpoint changed");
        for (int f = 0; f < 51; ++f) {
          const Real value = field(column.state, f, child);
          require(value >= 0 && std::isfinite(value), "invalid remapped state");
          if (children == 1)
            require(value == field(parent.state, f, k), "unsplit field changed");
          integrals[f] += w * value;
        }
      }
      require(close(mass, parent.op.pressure_weight[k]), "remap dry mass not conserved");
      for (int f = 0; f < 51; ++f)
        require(close(integrals[f], parent.op.pressure_weight[k] * field(parent.state, f, k)),
                "remap field integral not conserved");
      first = end;
    }
    Real child_mass = 0, native_mass = 0;
    std::array<Real, 51> child_integrals{}, native_integrals{};
    for (std::size_t k = 0; k < column.fixture.levels; ++k) {
      child_mass += column.op.pressure_weight[k];
      for (int f = 0; f < 51; ++f)
        child_integrals[f] += column.op.pressure_weight[k] * field(column.state, f, k);
    }
    for (std::size_t k = 0; k < parent.fixture.levels; ++k) {
      native_mass += parent.op.pressure_weight[k];
      for (int f = 0; f < 51; ++f)
        native_integrals[f] += parent.op.pressure_weight[k] * field(parent.state, f, k);
    }
    require(close(child_mass, native_mass), "column dry mass not conserved");
    for (int f = 0; f < 51; ++f)
      require(close(child_integrals[f], native_integrals[f]), "column field integral not conserved");
  }
  std::cout << "GEOMETRY_FIXTURES checks=" << checks << " failed=0\n";
}

std::vector<ControllerFixture> controller_study_fixtures(
    const bool controlled, const std::vector<PhysicalColumn> *physical) {
  if (!physical) return study_fixtures(controlled);
  std::vector<ControllerFixture> result;
  for (const auto &column : *physical)
    result.push_back({column.fixture, column.geometry_study ? "physical_geometry" : "physical_column", 0});
  return result;
}

const PhysicalColumn *physical_fixture(
    const int id, const std::vector<PhysicalColumn> *physical) {
  if (!physical) return nullptr;
  for (const auto &column : *physical)
    if (column.fixture.id == id) return &column;
  throw std::logic_error("physical fixture identity missing");
}

void controller_row(const char *record_kind,
                    const ControllerFixture &controller_fixture,
                    const Method method, const int steps,
                    const int production_count,
                    const CertifiedReferences &reference,
                    const Real mixing_rate_t, const Real activation_rate_t,
                    const std::optional<Real> fe_target,
                    const std::string &status, const std::string &reason,
                    const CandidateEvaluation *candidate,
                    const PhysicalColumn *column = nullptr) {
  const auto &fixture = controller_fixture.fixture;
  std::cout << (column && column->geometry_study ? "ndrop_physical_geometry_controller_accuracy_v1,"
                      : column ? "ndrop_physical_controller_accuracy_v1,"
                      : controller_fixture.activation_mixing == 0
                    ? "ndrop_controller_matched_accuracy_v1,"
                    : "ndrop_controlled_mixing_accuracy_v1,") << record_kind << ','
            << fixture.id << ',' << controller_fixture.family << ',';
  if (controller_fixture.requested_nsubmix > 0) {
    std::cout << controller_fixture.requested_nsubmix;
  }
  std::cout << ',' << fixture.levels << ',' << fixture.mixing << ',';
  if (controller_fixture.activation_mixing != 0)
    std::cout << controller_fixture.activation_mixing << ',';
  std::cout << fixture.interval << ',' << mixing_rate_t << ','
            << activation_rate_t << ','
            << ndrop_test::benchmark::method_name(method) << ',' << steps << ','
            << production_count << ','
            << (reference.ok ? "certified" : "unresolved") << ','
            << csv_text(reference.reason) << ',';
  if (std::isfinite(reference.disagreement)) {
    std::cout << reference.disagreement;
  }
  std::cout << ',';
  if (fe_target) std::cout << *fe_target;
  std::cout << ',' << status << ',' << csv_text(reason) << ',';
  if (candidate != nullptr && candidate->eligible) {
    std::cout << candidate->normal_error.screening_error << ','
              << candidate->tight_error.screening_error << ','
              << candidate->maximum_error << ','
              << field_array(candidate->normal_error.weighted_l1_by_field)
              << ',' << field_array(candidate->normal_error.linf_by_field)
              << ','
              << field_array(candidate->tight_error.weighted_l1_by_field)
              << ',' << field_array(candidate->tight_error.linf_by_field);
  } else {
    std::cout << ",,,,,,";
  }
  std::cout << ',';
  if (candidate != nullptr) {
    std::cout << candidate->result.diagnostics.clipping_count << ','
              << candidate->result.diagnostics.material_negative_count << ','
              << candidate->result.diagnostics.normalized_residual << ','
              << candidate->result.diagnostics.max_normalized_ledger_error
              << ','
              << ndrop_test::benchmark::state_checksum(candidate->result.state);
  } else {
    std::cout << ",,,,";
  }
  if (column) {
    std::cout << ',' << column->turbulence << ',' << column->velocity << ','
              << column->controller_rate_t << ',' << column->cloud_base << ','
              << column->height[column->cloud_base] -
                     column->height[column->cloud_base + 1] << ','
              << column->limiting_level << ','
              << column->thickness[column->limiting_level] << ',';
    if (candidate) std::cout << candidate->result.diagnostics.activation_limiter_count;
    if (column->geometry_study) {
      const auto k = static_cast<std::size_t>(column->limiting_level);
      const auto cloud_face = [&](const std::size_t face) {
        return face > 0 && face < column->fixture.levels &&
            column->op.cloud_fraction[face - 1] > .01 &&
            column->op.cloud_fraction[face - 1] - column->op.cloud_fraction[face] > .01;
      };
      std::cout << ',' << column->parent_grid_levels << ',' << column->subdivision
                << ",400,1000," << column->height[k] << ','
                << column->diffusivity[k] << ',' << column->diffusivity[k + 1]
                << ',' << static_cast<int>(cloud_face(k)) + static_cast<int>(cloud_face(k + 1));
    }
  }
  std::cout << '\n';
}

void controller_accuracy(const int maximum_steps, const int budget_seconds,
                         const std::string &fixture_ids,
                         const bool controlled = false,
                         const std::vector<PhysicalColumn> *physical = nullptr) {
  const auto budget_start = Clock::now();
  std::cout << std::setprecision(17)
            << "artifact_schema,record_kind,fixture,fixture_family,"
               "requested_nsubmix,levels,mixing_scale,";
  if (controlled) std::cout << "activation_mixing_scale,";
  std::cout << "interval_seconds,mixing_rate_t,activation_rate_t,method,steps,"
               "production_nsubmix,reference_status,reference_reason,"
               "reference_disagreement,fe_target_error,status,failure_reason,"
               "error_normal,error_tight,error_max,"
               "normal_weighted_l1_by_field,normal_linf_by_field,"
               "tight_weighted_l1_by_field,tight_linf_by_field,"
               "clipping_count,material_negative_count,normalized_residual,"
               "normalized_ledger_error,state_checksum";
  if (physical) std::cout << ",turbulence,cloud_base_wtke_m_s,controller_rate_t,"
                               "cloud_base_level,cloud_base_spacing_m,limiting_level,"
                               "limiting_spacing_m,activation_limiter_count";
  if (physical && physical->front().geometry_study)
    std::cout << ",parent_grid_levels,geometry_subdivision,refinement_lower_m,refinement_upper_m,"
                 "limiting_height_m,limiting_above_diffusivity_m2_s,"
                 "limiting_below_diffusivity_m2_s,limiting_cloud_base_faces";
  std::cout << '\n';
  const auto budget_exhausted = [&]() {
    return std::chrono::duration_cast<std::chrono::seconds>(Clock::now() -
                                                            budget_start)
               .count() >= budget_seconds;
  };
  for (const auto &controller_fixture : controller_study_fixtures(controlled, physical)) {
    const Fixture &fixture = controller_fixture.fixture;
    if (!fixture_ids.empty() &&
        ("," + fixture_ids + ",").find("," + std::to_string(fixture.id) + ",") ==
            std::string::npos) continue;
    const auto *column = physical_fixture(fixture.id, physical);
    const State state = column ? column->state : initial_state(fixture);
    const OperatorData op = column ? column->op : study_operator(controller_fixture);
    const auto prepared = ndrop_test::prepare_operator(op);
    const int production_count =
        ndrop_test::production_nsubmix(prepared, fixture.interval);
    Real mixing_rate_t = 0;
    Real activation_rate_t = 0;
    for (std::size_t level = 0; level < fixture.levels; ++level) {
      mixing_rate_t = std::max(
          mixing_rate_t,
          fixture.interval *
              (prepared.rate_below[level] + prepared.rate_above[level]));
      for (int mode = 0; mode < ndrop_test::num_modes; ++mode) {
        activation_rate_t = std::max(
            {activation_rate_t,
             fixture.interval * prepared.number_activation_rate[mode][level],
             fixture.interval * prepared.mass_activation_rate[mode][level]});
      }
    }
    const auto reference = references(state, op, fixture.interval);
    if (controller_fixture.requested_nsubmix > 0 &&
        production_count != controller_fixture.requested_nsubmix) {
      throw std::runtime_error("fixture controller count changed after creation");
    }
    if (production_count <= 0 || production_count > maximum_steps) {
      for (const Method method : {Method::imex, Method::backward_euler}) {
        controller_row("search_terminal", controller_fixture, method, 0,
                       production_count, reference, mixing_rate_t,
                       activation_rate_t, std::nullopt,
                       "fixture_controller_mismatch",
                       "production controller count is outside the study cap",
                       nullptr, column);
      }
      continue;
    }
    const auto fe = evaluate_controller_fe(state, op, fixture.interval,
                                           reference, production_count);
    const std::optional<Real> fe_target =
        fe.eligible ? std::optional<Real>(fe.maximum_error) : std::nullopt;
    controller_row("measurement", controller_fixture, Method::forward_euler,
                   production_count, production_count, reference,
                   mixing_rate_t, activation_rate_t, fe_target,
                   fe.eligible ? "fe_baseline" : fe.status, fe.reason, &fe, column);
    for (const Method method : {Method::imex, Method::backward_euler}) {
      if (!reference.ok || !fe.eligible) {
        controller_row("search_terminal", controller_fixture, method, 0,
                       production_count, reference, mixing_rate_t,
                       activation_rate_t, fe_target,
                       reference.ok ? "fe_baseline_ineligible"
                                    : "reference_unresolved",
                       reference.ok ? fe.reason : reference.reason, nullptr, column);
        continue;
      }
      int selected = 0;
      int last = 0;
      std::string terminal = "not_attained_within_cap";
      std::string terminal_reason =
          "no eligible matched-error count at or below cap";
      for (int steps = 1; steps <= maximum_steps; ++steps) {
        if (budget_exhausted()) {
          terminal = "accuracy_budget_exhausted";
          terminal_reason = "accuracy-stage wall-clock budget exhausted";
          break;
        }
        last = steps;
        auto candidate = evaluate_candidate(method, state, op,
                                            fixture.interval, steps, reference);
        const bool match = candidate.eligible &&
                           candidate.maximum_error <= *fe_target;
        const std::string status =
            match ? "eligible_match"
                  : candidate.eligible ? "eligible_above_target"
                                       : candidate.status;
        controller_row("measurement", controller_fixture, method, steps,
                       production_count, reference, mixing_rate_t,
                       activation_rate_t, fe_target, status, candidate.reason,
                       &candidate, column);
        if (match) {
          selected = steps;
          terminal = "selected";
          terminal_reason.clear();
          break;
        }
      }
      controller_row("search_terminal", controller_fixture, method,
                     selected > 0 ? selected : last, production_count,
                     reference, mixing_rate_t, activation_rate_t, fe_target,
                     terminal, terminal_reason, nullptr, column);
      if (selected > 0) {
        for (const int multiplier : {2, 4}) {
          if (selected > maximum_steps / multiplier || budget_exhausted()) {
            break;
          }
          const int steps = selected * multiplier;
          auto candidate = evaluate_candidate(method, state, op,
                                              fixture.interval, steps,
                                              reference);
          controller_row("confirmation", controller_fixture, method, steps,
                         production_count, reference, mixing_rate_t,
                         activation_rate_t, fe_target,
                         candidate.eligible ? "eligible_confirmation"
                                            : candidate.status,
                         candidate.reason, &candidate, column);
        }
      }
    }
    // Refinement evidence is diagnostic and never changes the controller
    // baseline or integer-search selection. It can withhold timing admission.
    if (fe.eligible) {
      for (const int multiplier : {2, 4}) {
        if (production_count > maximum_steps / multiplier || budget_exhausted())
          break;
        const int steps = production_count * multiplier;
        auto candidate = evaluate_candidate(Method::forward_euler, state, op,
                                             fixture.interval, steps, reference);
        controller_row("confirmation", controller_fixture,
                       Method::forward_euler, steps, production_count, reference,
                       mixing_rate_t, activation_rate_t, fe_target,
                       candidate.eligible ? "eligible_confirmation"
                                          : candidate.status,
                       candidate.reason, &candidate, column);
      }
    }
  }
}

struct Selected {
  int fixture;
  Method method;
  int steps;
  std::vector<Real> samples;
  std::vector<std::uint64_t> checksums;
};

Method parse_method(const std::string &name) {
  for (const auto method : methods) {
    if (name == ndrop_test::benchmark::method_name(method)) return method;
  }
  throw std::invalid_argument("selection has an unsupported method");
}

std::vector<Selected> parse_selection(const std::string &encoded,
                                      const bool physical = false) {
  std::vector<Selected> selections;
  std::istringstream rows(encoded);
  std::string row;
  while (std::getline(rows, row, ';')) {
    if (row.empty()) throw std::invalid_argument("empty selection row");
    std::istringstream fields(row);
    std::string fixture_text, method_text, steps_text, extra;
    if (!std::getline(fields, fixture_text, ':') ||
        !std::getline(fields, method_text, ':') ||
        !std::getline(fields, steps_text, ':') ||
        std::getline(fields, extra, ':')) {
      throw std::invalid_argument("selection must be fixture:method:steps");
    }
    std::size_t fixture_used = 0;
    const int id = std::stoi(fixture_text, &fixture_used);
    if (fixture_used != fixture_text.size() || std::to_string(id) != fixture_text ||
        id < 0) {
      throw std::invalid_argument("selection fixture is invalid");
    }
    const int steps = parse_positive("steps", steps_text);
    if ((physical ? (id < 28 || id > 31)
                  : (id >= static_cast<int>(fixtures.size()) ||
                     (steps & (steps - 1)) != 0)) || steps > 8192) {
      throw std::invalid_argument("selection is outside the frozen matrix");
    }
    const Method method = parse_method(method_text);
    for (const auto &existing : selections) {
      if (existing.fixture == id && existing.method == method &&
          existing.steps == steps) {
        throw std::invalid_argument("duplicate selection row");
      }
    }
    selections.push_back({id, method, steps, {}, {}});
  }
  if (selections.empty()) throw std::invalid_argument("selection is empty");
  std::sort(selections.begin(), selections.end(), [](const Selected &a,
                                                      const Selected &b) {
    if (a.fixture != b.fixture) return a.fixture < b.fixture;
    if (a.method != b.method) return a.method < b.method;
    return a.steps < b.steps;
  });
  return selections;
}

std::vector<Selected> parse_controller_selection(const std::string &encoded,
                                                 const int maximum_steps,
                                                 const bool controlled = false,
                                                 const bool physical = false,
                                                 const bool geometry = false) {
  if (geometry && (!physical || controlled))
    throw std::invalid_argument("geometry selection requires its physical policy");
  std::vector<Selected> selections;
  std::istringstream rows(encoded);
  std::string row;
  while (std::getline(rows, row, ';')) {
    if (row.empty()) throw std::invalid_argument("empty selection row");
    std::istringstream fields(row);
    std::string fixture_text, method_text, steps_text, extra;
    if (!std::getline(fields, fixture_text, ':') ||
        !std::getline(fields, method_text, ':') ||
        !std::getline(fields, steps_text, ':') ||
        std::getline(fields, extra, ':')) {
      throw std::invalid_argument("selection must be fixture:method:steps");
    }
    std::size_t fixture_used = 0;
    const int id = std::stoi(fixture_text, &fixture_used);
    if (fixture_used != fixture_text.size() || std::to_string(id) != fixture_text ||
        (physical ? (id < 28 || id > (geometry ? 39 : 31))
                  : controlled ? (id < 24 || id > 27) : (id < 0 || id >= 24))) {
      throw std::invalid_argument("controller selection fixture is invalid");
    }
    const int steps = parse_positive("steps", steps_text);
    if (steps > maximum_steps) {
      throw std::invalid_argument("controller selection exceeds the count cap");
    }
    const Method method = parse_method(method_text);
    for (const auto &existing : selections) {
      if (existing.fixture == id && existing.method == method) {
        throw std::invalid_argument("duplicate controller fixture/method");
      }
    }
    selections.push_back({id, method, steps, {}, {}});
  }
  if (selections.empty()) throw std::invalid_argument("selection is empty");
  std::sort(selections.begin(), selections.end(), [](const Selected &a,
                                                      const Selected &b) {
    if (a.fixture != b.fixture) return a.fixture < b.fixture;
    return a.method < b.method;
  });
  return selections;
}

// Parser-only coverage: this test must never call the timing routine.
void check_geometry_selection() {
  const auto selected = parse_controller_selection(
      "39:backward_euler:90;28:forward_euler:1;38:imex:10", 8192, false, true, true);
  if (selected.size() != 3 || selected.front().fixture != 28 ||
      selected.back().fixture != 39 || selected.back().steps != 90)
    throw std::runtime_error("geometry timing selection parser changed");
  int rejected = 0;
  for (const std::string invalid : {"27:imex:1", "40:imex:1", "32:imex:0",
                                   "32:imex:8193", "32:imex:1;32:imex:2", ""}) {
    try {
      (void)parse_controller_selection(invalid, 8192, false, true, true);
    } catch (const std::invalid_argument &) {
      ++rejected;
    }
  }
  if (rejected != 6) throw std::runtime_error("geometry invalid selection admitted");
  try {
    (void)parse_controller_selection("32:imex:1", 8192, false, true);
    throw std::runtime_error("geometry fixture entered native timing policy");
  } catch (const std::invalid_argument &) {
    ++rejected;
  }
  std::cout << "GEOMETRY_SELECTION checks=" << 1 + rejected << " failed=0\n";
}

void timing(std::vector<Selected> selections, const int warmups,
            const int min_blocks, const int max_blocks,
            const int budget_seconds,
            const std::vector<PhysicalColumn> *physical = nullptr) {
  const auto budget_start = Clock::now();
  std::cout << std::setprecision(17)
            << "artifact_schema,fixture,method,steps,block,seconds_per_column,"
               "state_checksum,stop_reason\n";
  std::vector<Fixture> timed_fixtures(fixtures.begin(), fixtures.end());
  if (physical) {
    timed_fixtures.clear();
    for (const auto &column : *physical) timed_fixtures.push_back(column.fixture);
  }
  for (const Fixture &fixture : timed_fixtures) {
    std::vector<Selected *> group;
    for (auto &item : selections) {
      if (item.fixture == fixture.id) group.push_back(&item);
    }
    if (group.empty()) continue;
    const PhysicalColumn *column = nullptr;
    if (physical) {
      for (const auto &item : *physical)
        if (item.fixture.id == fixture.id) column = &item;
    }
    const State state = column ? column->state : initial_state(fixture);
    const OperatorData op = column ? column->op : operator_data(fixture);
    const int controller = column ? ndrop_test::production_nsubmix(
        ndrop_test::prepare_operator(op), fixture.interval) : 0;
    for (const auto *selected : group)
      if (column && selected->method == Method::forward_euler && selected->steps < controller)
        throw std::invalid_argument("physical FE selection below controller count");
    const auto measure = [&](Selected &selected, const bool timed) {
      AdvanceControl control;
      control.collect_detailed_diagnostics = false;
      const auto start = Clock::now();
      const auto result = column && selected.method == Method::forward_euler &&
                                       selected.steps == controller
          ? advance_controller_fe(state, op, fixture.interval, control)
          : advance(selected.method, state, op, fixture.interval, selected.steps, control);
      Kokkos::fence();
      const auto stop = Clock::now();
      if (!result.ok || result.diagnostics.clipping_count != 0 ||
          result.diagnostics.material_negative_count != 0) {
        throw std::runtime_error("timed advance failed numerical eligibility");
      }
      if (timed) {
        selected.samples.push_back(
            std::chrono::duration<Real>(stop - start).count());
        selected.checksums.push_back(
            ndrop_test::benchmark::state_checksum(result.state));
      }
    };
    for (int repeat = 0; repeat < warmups; ++repeat) {
      for (std::size_t offset = 0; offset < group.size(); ++offset) {
        measure(*group[(offset + static_cast<std::size_t>(repeat)) %
                       group.size()], false);
      }
    }
    std::string reason = "maximum_blocks";
    for (int block = 0; block < max_blocks; ++block) {
      for (std::size_t offset = 0; offset < group.size(); ++offset) {
        measure(*group[(offset + static_cast<std::size_t>(block)) %
                       group.size()], true);
      }
      if (block + 1 >= min_blocks) {
        bool precision_met = true;
        for (const Selected *selected : group) {
          const auto interval =
              ndrop_test::benchmark::bootstrap_median_interval(
                  selected->samples, 10000);
          precision_met = precision_met &&
                          interval.relative_half_width < 0.05;
        }
        if (precision_met) {
          reason = "median_precision";
          break;
        }
      }
      if (std::chrono::duration_cast<std::chrono::seconds>(
              Clock::now() - budget_start)
              .count() >= budget_seconds) {
        reason = "budget_exhausted";
        break;
      }
    }
    for (const Selected *selected : group) {
      for (std::size_t block = 0; block < selected->samples.size(); ++block) {
        std::cout << (column ? "ndrop_physical_column_timing_v1,"
                            : "ndrop_matched_timing_v1,") << fixture.id << ','
                  << ndrop_test::benchmark::method_name(selected->method)
                  << ',' << selected->steps << ',' << block << ','
                  << selected->samples[block] << ','
                  << selected->checksums[block] << ',' << reason << '\n';
      }
    }
  }
}

void controller_timing(std::vector<Selected> selections, const int warmups,
                       const int min_blocks, const int max_blocks,
                       const int budget_seconds, const bool controlled = false,
                       const std::vector<PhysicalColumn> *physical = nullptr) {
  const auto budget_start = Clock::now();
  std::cout << std::setprecision(17)
            << "artifact_schema,fixture,method,steps,block,seconds_per_column,"
               "state_checksum,stop_reason\n";
  for (const auto &controller_fixture : controller_study_fixtures(controlled, physical)) {
    const Fixture &fixture = controller_fixture.fixture;
    std::vector<Selected *> group;
    for (auto &item : selections) {
      if (item.fixture == fixture.id) group.push_back(&item);
    }
    if (group.empty()) continue;
    const auto *column = physical_fixture(fixture.id, physical);
    const State state = column ? column->state : initial_state(fixture);
    const OperatorData op = column ? column->op : study_operator(controller_fixture);
    const int controller_count = ndrop_test::production_nsubmix(
        ndrop_test::prepare_operator(op), fixture.interval);
    const auto measure = [&](Selected &selected, const bool timed) {
      AdvanceControl control;
      control.collect_detailed_diagnostics = false;
      const auto start = Clock::now();
      const auto result = selected.method == Method::forward_euler
                              ? advance_controller_fe(state, op,
                                                      fixture.interval, control)
                              : advance(selected.method, state, op,
                                        fixture.interval, selected.steps,
                                        control);
      Kokkos::fence();
      const auto stop = Clock::now();
      if (!result.ok || result.diagnostics.clipping_count != 0 ||
          result.diagnostics.material_negative_count != 0 ||
          (selected.method == Method::forward_euler &&
           (selected.steps != controller_count ||
            result.diagnostics.nsubmix != controller_count))) {
        throw std::runtime_error(
            "controller-matched timed advance failed eligibility");
      }
      if (timed) {
        selected.samples.push_back(
            std::chrono::duration<Real>(stop - start).count());
        selected.checksums.push_back(
            ndrop_test::benchmark::state_checksum(result.state));
      }
    };
    for (int repeat = 0; repeat < warmups; ++repeat) {
      for (std::size_t offset = 0; offset < group.size(); ++offset) {
        measure(*group[(offset + static_cast<std::size_t>(repeat)) %
                       group.size()], false);
      }
    }
    std::string reason = "maximum_blocks";
    for (int block = 0; block < max_blocks; ++block) {
      if (std::chrono::duration_cast<std::chrono::seconds>(
              Clock::now() - budget_start)
              .count() >= budget_seconds) {
        reason = "budget_exhausted";
        break;
      }
      for (std::size_t offset = 0; offset < group.size(); ++offset) {
        measure(*group[(offset + static_cast<std::size_t>(block)) %
                       group.size()], true);
      }
      if (block + 1 >= min_blocks) {
        bool precision_met = true;
        for (const Selected *selected : group) {
          const auto interval =
              ndrop_test::benchmark::bootstrap_median_interval(
                  selected->samples, 10000);
          precision_met = precision_met &&
                          interval.relative_half_width < 0.05;
        }
        if (precision_met) {
          reason = "median_precision";
          break;
        }
      }
    }
    for (const Selected *selected : group) {
      for (std::size_t block = 0; block < selected->samples.size(); ++block) {
        std::cout << (column && column->geometry_study ? "ndrop_physical_geometry_controller_timing_v1,"
                               : physical ? "ndrop_physical_controller_timing_v1,"
                               : controlled ? "ndrop_controlled_mixing_timing_v1,"
                                : "ndrop_controller_matched_timing_v1,")
                  << fixture.id << ','
                  << ndrop_test::benchmark::method_name(selected->method) << ','
                  << selected->steps << ',' << block << ','
                  << selected->samples[block] << ','
                  << selected->checksums[block] << ',' << reason << '\n';
      }
    }
  }
}

// Untimed fixture/entry-point checks, not a benchmark or production replay.
void check_controlled_fixtures() {
  int checks = 0;
  const auto require = [&](const bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
    ++checks;
  };
  const auto bits_equal = [](const std::vector<Real> &a,
                             const std::vector<Real> &b) {
    return a.size() == b.size() &&
           std::memcmp(a.data(), b.data(), a.size() * sizeof(Real)) == 0;
  };
  const auto controlled = study_fixtures(true);
  for (std::size_t index = 0; index < controlled.size(); index += 2) {
    const auto &moderate = controlled[index];
    const auto &strong = controlled[index + 1];
    const auto a = study_operator(moderate);
    const auto b = study_operator(strong);
    const auto pa = ndrop_test::prepare_operator(a);
    const auto pb = ndrop_test::prepare_operator(b);
    const State state = initial_state(moderate.fixture);
    const State other = initial_state(strong.fixture);
    require(bits_equal(state.droplet_number, other.droplet_number), "initial droplet mismatch");
    require(bits_equal(a.pressure_weight, b.pressure_weight), "weight mismatch");
    require(bits_equal(a.cloud_fraction, b.cloud_fraction), "cloud mismatch");
    require(pa.activation_limiter_count == 0 && pb.activation_limiter_count == 0,
            "unexpected limiting");
    require(a.first_level_controller_rate == 0.25 * moderate.fixture.mixing &&
                b.first_level_controller_rate == 0.25 * strong.fixture.mixing,
            "controller input mismatch");
    for (int field = 0; field < ndrop_test::num_constituents; ++field) {
      require(bits_equal(state.interstitial[field], other.interstitial[field]) &&
                  bits_equal(state.cloudborne[field], other.cloudborne[field]),
              "initial aerosol mismatch");
    }
    for (int mode = 0; mode < ndrop_test::num_modes; ++mode) {
      require(bits_equal(a.number_activation_rate[mode], b.number_activation_rate[mode]) &&
                  bits_equal(a.mass_activation_rate[mode], b.mass_activation_rate[mode]) &&
                  bits_equal(pa.number_activation_rate[mode], pb.number_activation_rate[mode]) &&
                  bits_equal(pa.mass_activation_rate[mode], pb.mass_activation_rate[mode]),
              "activation mismatch");
      require(a.number_activation_rate[mode].back() == 0 &&
                  a.mass_activation_rate[mode].back() == 0, "bottom activation mismatch");
    }
    for (const auto *fixture : {&moderate, &strong}) {
      const auto op = study_operator(*fixture);
      const auto legacy = operator_data(fixture->fixture);
      require(bits_equal(op.interface_exchange, legacy.interface_exchange), "exchange mismatch");
      const int count = ndrop_test::production_nsubmix(
          ndrop_test::prepare_operator(op), fixture->fixture.interval);
      for (const Method method : methods) {
        AdvanceControl off;
        off.collect_detailed_diagnostics = false;
        const auto run = [&](const AdvanceControl &control) {
          return method == Method::forward_euler
                     ? advance_controller_fe(state, op, fixture->fixture.interval, control)
                     : advance(method, state, op, fixture->fixture.interval, count, control);
        };
        const auto detailed = run(AdvanceControl{});
        const auto primary = run(off);
        require(detailed.ok && primary.ok && detailed.message == primary.message &&
                    detailed.diagnostics.nsubmix == primary.diagnostics.nsubmix &&
                    detailed.diagnostics.clipping_count == primary.diagnostics.clipping_count &&
                    detailed.diagnostics.material_negative_count == primary.diagnostics.material_negative_count,
                "diagnostics changed a decision");
        require(method != Method::forward_euler || detailed.diagnostics.nsubmix == count,
                "FE bypassed the controller");
        require(bits_equal(detailed.state.droplet_number, primary.state.droplet_number),
                "diagnostics changed droplet state");
        for (int field = 0; field < ndrop_test::num_constituents; ++field) {
          require(bits_equal(detailed.state.interstitial[field], primary.state.interstitial[field]) &&
                      bits_equal(detailed.state.cloudborne[field], primary.state.cloudborne[field]),
                  "diagnostics changed aerosol state");
        }
      }
    }
  }
  std::cout << "CONTROLLED_FIXTURES checks=" << checks << " failed=0\n";
}

} // namespace

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  try {
    bool do_accuracy = false;
    bool do_timing = false;
    bool do_controller_accuracy = false;
    bool do_controller_timing = false;
    bool controlled = false;
    bool do_fixture_check = false;
    bool do_physical_accuracy = false, do_physical_timing = false;
    bool do_physical_check = false, physical_controller = false;
    bool geometry = false;
    std::string physical_grid_dir;
    std::string selection;
    std::string fixture_ids;
    int warmups = 10, min_blocks = 30, max_blocks = 100;
    int budget_seconds = 7200;
    int maximum_steps = 8192;
    for (int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      const auto next = [&]() -> std::string {
        if (index + 1 >= argc) throw std::invalid_argument("missing option value");
        return argv[++index];
      };
      if (option == "--accuracy") do_accuracy = true;
      else if (option == "--physical-column-accuracy") do_physical_accuracy = true;
      else if (option == "--physical-column-timing") do_physical_timing = true;
      else if (option == "--check-physical-column-fixtures") do_physical_check = true;
      else if (option == "--check-physical-geometry-fixtures") {
        do_physical_check = true;
        geometry = true;
      } else if (option == "--physical-geometry-controller-accuracy") {
        do_controller_accuracy = true;
        physical_controller = true;
        geometry = true;
      } else if (option == "--physical-geometry-controller-timing") {
        do_controller_timing = true;
        physical_controller = true;
        geometry = true;
      }
      else if (option == "--physical-grid-dir") physical_grid_dir = next();
      else if (option == "--physical-controller-accuracy") {
        do_controller_accuracy = true;
        physical_controller = true;
      } else if (option == "--physical-controller-timing") {
        do_controller_timing = true;
        physical_controller = true;
      }
      else if (option == "--timing") do_timing = true;
      else if (option == "--controller-accuracy")
        do_controller_accuracy = true;
      else if (option == "--controller-timing")
        do_controller_timing = true;
      else if (option == "--controlled-mixing-accuracy") {
        do_controller_accuracy = true;
        controlled = true;
      } else if (option == "--controlled-mixing-timing") {
        do_controller_timing = true;
        controlled = true;
      } else if (option == "--check-controlled-mixing-fixtures")
        do_fixture_check = true;
      else if (option == "--selection") selection = next();
      else if (option == "--fixtures") fixture_ids = next();
      else if (option == "--warmups") warmups = parse_positive(option, next());
      else if (option == "--min-blocks")
        min_blocks = parse_positive(option, next());
      else if (option == "--max-blocks")
        max_blocks = parse_positive(option, next());
      else if (option == "--budget-seconds")
        budget_seconds = parse_positive(option, next());
      else if (option == "--max-steps")
        maximum_steps = parse_positive(option, next());
      else throw std::invalid_argument("unknown option: " + option);
    }
    const int mode_count = static_cast<int>(do_accuracy) +
                           static_cast<int>(do_timing) +
                           static_cast<int>(do_controller_accuracy) +
                           static_cast<int>(do_controller_timing) +
                           static_cast<int>(do_fixture_check) +
                           static_cast<int>(do_physical_accuracy) +
                           static_cast<int>(do_physical_timing) +
                           static_cast<int>(do_physical_check);
    if (mode_count != 1 || min_blocks > max_blocks) {
      throw std::invalid_argument("choose one mode with valid block limits");
    }
    if (controlled && physical_controller)
      throw std::invalid_argument("physical and controlled policies are mutually exclusive");
    if (!fixture_ids.empty()) {
      if (!do_controller_accuracy)
        throw std::invalid_argument("--fixtures requires --controller-accuracy");
      std::istringstream ids(fixture_ids);
      std::string token;
      int previous = -1;
      while (std::getline(ids, token, ',')) {
        std::size_t consumed = 0;
        const int id = std::stoi(token, &consumed);
        if (consumed != token.size() || id <= previous || id >= 24 ||
            token != std::to_string(id))
          throw std::invalid_argument("fixture IDs must be sorted unique 0..23");
        previous = id;
      }
      if (fixture_ids.back() == ',')
        throw std::invalid_argument("empty fixture ID");
    }
    if ((controlled || physical_controller) && !fixture_ids.empty())
      throw std::invalid_argument("physical/controlled mixing requires exact four-case coverage");
    if (!physical_grid_dir.empty() && !do_physical_accuracy &&
        !do_physical_timing && !do_physical_check && !physical_controller)
      throw std::invalid_argument("physical grids require a physical-column mode");
    if (physical_controller) {
      const auto columns = geometry ? physical_geometry_columns(physical_grid_dir)
                                    : physical_columns(physical_grid_dir);
      if (do_controller_accuracy)
        controller_accuracy(maximum_steps, budget_seconds, "", false, &columns);
      else
        controller_timing(parse_controller_selection(selection, maximum_steps, false, true, geometry),
                          warmups, min_blocks, max_blocks, budget_seconds, false, &columns);
    }
    else if (do_physical_accuracy || do_physical_timing || do_physical_check) {
      const auto columns = geometry ? physical_geometry_columns(physical_grid_dir)
                                    : physical_columns(physical_grid_dir);
      if (do_physical_check && geometry) {
        check_geometry_fixtures(columns, physical_grid_dir);
        check_geometry_selection();
      }
      else if (do_physical_check) check_physical_fixtures(columns);
      else if (do_physical_accuracy)
        physical_accuracy(columns, maximum_steps, budget_seconds);
      else timing(parse_selection(selection, true), warmups, min_blocks,
                  max_blocks, budget_seconds, &columns);
    }
    else if (do_fixture_check) check_controlled_fixtures();
    else if (do_accuracy) accuracy();
    else if (do_timing)
      timing(parse_selection(selection), warmups, min_blocks, max_blocks,
             budget_seconds);
    else if (do_controller_accuracy)
      controller_accuracy(maximum_steps, budget_seconds, fixture_ids, controlled);
    else
      controller_timing(
          parse_controller_selection(selection, maximum_steps, controlled), warmups,
          min_blocks, max_blocks, budget_seconds, controlled);
    Kokkos::finalize();
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    Kokkos::finalize();
    std::cerr << "ERROR: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
