#include "ndrop_benchmark.hpp"
#include "ndrop_integrator.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using ndrop_test::AdvanceControl;
using ndrop_test::AdvanceDiagnostics;
using ndrop_test::AdvancePhaseTimings;
using ndrop_test::AdvanceResult;
using ndrop_test::OperatorData;
using ndrop_test::Real;
using ndrop_test::State;
using ndrop_test::benchmark::BootstrapInterval;
using ndrop_test::benchmark::Case;
using ndrop_test::benchmark::Method;
using ndrop_test::num_constituents;
using ndrop_test::num_modes;

struct Arguments {
  bool smoke = false;
  std::string output_path;
  int case_start = 0;
  int case_count = std::numeric_limits<int>::max();
  int warmups = 10;
  int minimum_repetitions = 30;
  int maximum_repetitions = 100;
  int bootstrap_resamples = 10000;
  int progress_every = 10;
};

struct BatchMeasurement {
  bool ok = true;
  std::string message = "ok";
  Real seconds_per_column = 0;
  Real fence_seconds_per_column = 0;
  AdvancePhaseTimings phases_per_column;
  std::uint64_t checksum = 0;
  AdvanceDiagnostics diagnostics;
};

struct RowResult {
  Case benchmark_case;
  bool validation_ok = false;
  bool eligible = false;
  std::string message;
  std::string exclusion_reason;
  AdvanceDiagnostics validation_diagnostics;
  std::vector<Real> samples;
  std::vector<Real> instrumented_samples;
  std::vector<Real> fence_samples;
  std::array<std::vector<Real>, 5> phase_samples;
  BootstrapInterval primary_interval;
  BootstrapInterval instrumented_interval;
  std::array<Real, 5> phase_medians{};
  Real fence_median = 0;
  Real instrumentation_overhead = 0;
  std::uint64_t checksum = 0;
  std::string stop_reason;
};

std::string csv_escape(const std::string &value) {
  if (value.find_first_of(",\"\n") == std::string::npos) {
    return value;
  }
  std::string escaped = "\"";
  for (const char character : value) {
    escaped += character == '\"' ? "\"\"" : std::string(1, character);
  }
  escaped += '\"';
  return escaped;
}

std::string real_text(const Real value) {
  if (!std::isfinite(value)) {
    return {};
  }
  std::ostringstream text;
  text << std::setprecision(17) << value;
  return text.str();
}

std::string hex_text(const std::uint64_t value) {
  std::ostringstream text;
  text << std::hex << std::setw(16) << std::setfill('0') << value;
  return text.str();
}

std::string joined_samples(const std::vector<Real> &samples) {
  std::ostringstream text;
  text << std::setprecision(17);
  for (std::size_t index = 0; index < samples.size(); ++index) {
    if (index > 0) {
      text << ';';
    }
    text << samples[index];
  }
  return text.str();
}

int parse_nonnegative(const std::string &name, const std::string &value) {
  std::size_t consumed = 0;
  const long parsed = std::stol(value, &consumed);
  if (consumed != value.size() || parsed < 0 ||
      parsed > std::numeric_limits<int>::max()) {
    throw std::invalid_argument(name + " must be a nonnegative integer");
  }
  return static_cast<int>(parsed);
}

Arguments parse_arguments(const int argc, char **argv) {
  Arguments arguments;
  bool repetitions_overridden = false;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    auto next_value = [&]() -> std::string {
      if (index + 1 >= argc) {
        throw std::invalid_argument("missing value for " + argument);
      }
      return argv[++index];
    };
    if (argument == "--smoke") {
      arguments.smoke = true;
    } else if (argument == "--output") {
      arguments.output_path = next_value();
    } else if (argument == "--case-start") {
      arguments.case_start = parse_nonnegative(argument, next_value());
    } else if (argument == "--case-count") {
      arguments.case_count = parse_nonnegative(argument, next_value());
    } else if (argument == "--warmups") {
      arguments.warmups = parse_nonnegative(argument, next_value());
      repetitions_overridden = true;
    } else if (argument == "--min-repetitions") {
      arguments.minimum_repetitions =
          parse_nonnegative(argument, next_value());
      repetitions_overridden = true;
    } else if (argument == "--max-repetitions") {
      arguments.maximum_repetitions =
          parse_nonnegative(argument, next_value());
      repetitions_overridden = true;
    } else if (argument == "--bootstrap-resamples") {
      arguments.bootstrap_resamples =
          parse_nonnegative(argument, next_value());
      repetitions_overridden = true;
    } else if (argument == "--progress-every") {
      arguments.progress_every = parse_nonnegative(argument, next_value());
    } else {
      throw std::invalid_argument("unknown argument: " + argument);
    }
  }
  if (arguments.output_path.empty()) {
    throw std::invalid_argument("--output is required");
  }
  if (arguments.smoke && !repetitions_overridden) {
    arguments.warmups = 1;
    arguments.minimum_repetitions = 3;
    arguments.maximum_repetitions = 3;
    arguments.bootstrap_resamples = 100;
  }
  if (arguments.minimum_repetitions <= 0 ||
      arguments.maximum_repetitions < arguments.minimum_repetitions ||
      arguments.bootstrap_resamples <= 0) {
    throw std::invalid_argument(
        "repetition and bootstrap limits are inconsistent");
  }
  return arguments;
}

State make_initial_state(const Case &benchmark_case) {
  State state(benchmark_case.nlev);
  for (std::size_t level = 0; level < benchmark_case.nlev; ++level) {
    state.droplet_number[level] = 4 + 0.01 * level;
    for (int field = 0; field < benchmark_case.active_pairs; ++field) {
      state.interstitial[field][level] =
          1 + 0.01 * field + 0.002 * level;
      state.cloudborne[field][level] =
          0.2 + 0.002 * field + 0.001 * level;
    }
  }
  return state;
}

OperatorData make_operator(const Case &benchmark_case) {
  OperatorData op(benchmark_case.nlev);
  for (std::size_t level = 0; level < benchmark_case.nlev; ++level) {
    op.pressure_weight[level] =
        0.8 + 0.4 * static_cast<Real>((17 * level + 3) % 11) / 10;
    op.cloud_fraction[level] = level % 2 == 0 ? 0.8 : 0.55;
  }
  for (std::size_t interface = 1; interface < benchmark_case.nlev;
       ++interface) {
    const int offset = static_cast<int>((7 * interface + 1) % 5) - 2;
    const Real modulation = 1 + 0.05 * offset;
    op.interface_exchange[interface] =
        benchmark_case.mixing_scale * modulation *
        std::min(op.pressure_weight[interface - 1],
                 op.pressure_weight[interface]);
  }
  constexpr std::array<Real, num_modes> number_multiplier = {0.35, 0.4,
                                                               0.45, 0.5};
  constexpr std::array<Real, num_modes> mass_multiplier = {0.25, 0.3, 0.35,
                                                             0.4};
  for (int mode = 0; mode < num_modes; ++mode) {
    for (std::size_t level = 0; level + 1 < benchmark_case.nlev; ++level) {
      const Real below = op.interface_exchange[level + 1] /
                         op.pressure_weight[level];
      op.number_activation_rate[mode][level] =
          number_multiplier[mode] * below;
      op.mass_activation_rate[mode][level] = mass_multiplier[mode] * below;
    }
  }
  op.first_level_controller_rate = 0.25 * benchmark_case.mixing_scale;
  return op;
}

AdvanceResult advance_case(const Case &benchmark_case, const State &state,
                           const OperatorData &op,
                           const AdvanceControl &control) {
  switch (benchmark_case.method) {
  case Method::forward_euler: {
    ndrop_test::ForwardEulerOptions options;
    options.forced_substeps = benchmark_case.forced_fe_substeps;
    return ndrop_test::advance_forward_euler(
        state, op, benchmark_case.dt, options, control);
  }
  case Method::imex:
    return ndrop_test::advance_imex(state, op, benchmark_case.dt, control);
  case Method::backward_euler:
    return ndrop_test::advance_backward_euler(state, op, benchmark_case.dt,
                                              control);
  }
  throw std::logic_error("unknown benchmark method");
}

void add_phase(AdvancePhaseTimings &sum,
               const AdvancePhaseTimings &value) {
  sum.setup += value.setup;
  sum.factorization += value.factorization;
  sum.solve += value.solve;
  sum.update += value.update;
  sum.diagnostics += value.diagnostics;
}

void combine_checksum(std::uint64_t &accumulator, const std::uint64_t value) {
  accumulator ^=
      value + UINT64_C(0x9e3779b97f4a7c15) + (accumulator << 6) +
      (accumulator >> 2);
}

BatchMeasurement run_batch(const Case &benchmark_case,
                           const std::vector<State> &inputs,
                           const OperatorData &op, const bool detailed,
                           const bool phases) {
  BatchMeasurement measurement;
  std::vector<State> outputs(static_cast<std::size_t>(benchmark_case.batch_size));
  AdvancePhaseTimings phase_sum;
  AdvanceControl control;
  control.collect_detailed_diagnostics = detailed;
  const auto started = Clock::now();
  for (int column = 0; column < benchmark_case.batch_size; ++column) {
    AdvancePhaseTimings column_phases;
    control.phase_timings = phases ? &column_phases : nullptr;
    AdvanceResult result = advance_case(
        benchmark_case, inputs[static_cast<std::size_t>(column)], op, control);
    if (!result.ok && measurement.ok) {
      measurement.ok = false;
      measurement.message = result.message;
    }
    measurement.diagnostics = result.diagnostics;
    outputs[static_cast<std::size_t>(column)] = std::move(result.state);
    if (phases) {
      add_phase(phase_sum, column_phases);
    }
  }
  const auto before_fence = Clock::now();
  Kokkos::fence();
  const auto stopped = Clock::now();
  const Real batch = static_cast<Real>(benchmark_case.batch_size);
  measurement.seconds_per_column =
      std::chrono::duration<Real>(stopped - started).count() / batch;
  measurement.fence_seconds_per_column =
      std::chrono::duration<Real>(stopped - before_fence).count() / batch;
  if (phases) {
    measurement.phases_per_column.setup = phase_sum.setup / batch;
    measurement.phases_per_column.factorization =
        phase_sum.factorization / batch;
    measurement.phases_per_column.solve = phase_sum.solve / batch;
    measurement.phases_per_column.update = phase_sum.update / batch;
    measurement.phases_per_column.diagnostics =
        phase_sum.diagnostics / batch;
  }
  for (const State &output : outputs) {
    combine_checksum(measurement.checksum,
                     ndrop_test::benchmark::state_checksum(output));
  }
  return measurement;
}

bool same_group(const Case &lhs, const Case &rhs) {
  return lhs.nlev == rhs.nlev && lhs.active_pairs == rhs.active_pairs &&
         lhs.batch_size == rhs.batch_size &&
         lhs.mixing_index == rhs.mixing_index && lhs.dt == rhs.dt;
}

std::vector<std::vector<Case>> selected_groups(const Arguments &arguments) {
  const auto all = ndrop_test::benchmark::required_cases(arguments.smoke);
  const long long end =
      static_cast<long long>(arguments.case_start) + arguments.case_count;
  std::vector<Case> selected;
  for (const Case &benchmark_case : all) {
    if (benchmark_case.id >= arguments.case_start &&
        static_cast<long long>(benchmark_case.id) < end) {
      selected.push_back(benchmark_case);
    }
  }
  if (selected.empty()) {
    throw std::invalid_argument("selected benchmark case range is empty");
  }
  std::vector<std::vector<Case>> groups;
  for (const Case &benchmark_case : selected) {
    if (groups.empty() || !same_group(groups.back().front(), benchmark_case)) {
      groups.push_back({});
    }
    groups.back().push_back(benchmark_case);
  }
  return groups;
}

void validate_row(RowResult &row) {
  const State state = make_initial_state(row.benchmark_case);
  const OperatorData op = make_operator(row.benchmark_case);
  const AdvanceResult validation =
      advance_case(row.benchmark_case, state, op, AdvanceControl{});
  AdvanceControl lean_control;
  lean_control.collect_detailed_diagnostics = false;
  const AdvanceResult lean =
      advance_case(row.benchmark_case, state, op, lean_control);
  if (!ndrop_test::benchmark::timing_path_matches_validated(validation,
                                                             lean)) {
    throw std::runtime_error(
        "diagnostics-off advance disagrees with validated path for case " +
        std::to_string(row.benchmark_case.id) + ": detailed_ok=" +
        (validation.ok ? "true" : "false") + ", detailed_message=" +
        validation.message + ", lean_ok=" + (lean.ok ? "true" : "false") +
        ", lean_message=" + lean.message + ", detailed_clipping=" +
        std::to_string(validation.diagnostics.clipping_count) +
        ", lean_clipping=" +
        std::to_string(lean.diagnostics.clipping_count) +
        ", detailed_material_negative=" +
        std::to_string(validation.diagnostics.material_negative_count) +
        ", lean_material_negative=" +
        std::to_string(lean.diagnostics.material_negative_count) +
        ", detailed_checksum=" + hex_text(ndrop_test::benchmark::state_checksum(
                                     validation.state)) +
        ", lean_checksum=" +
        hex_text(ndrop_test::benchmark::state_checksum(lean.state)));
  }
  row.validation_ok = validation.ok;
  row.validation_diagnostics = validation.diagnostics;
  row.message = validation.message;
  if (!validation.ok) {
    row.eligible = false;
    row.exclusion_reason = validation.message;
  } else if (row.benchmark_case.method == Method::forward_euler &&
             validation.diagnostics.clipping_count > 0) {
    row.eligible = false;
    row.exclusion_reason = "forward_euler_clipping";
  } else {
    row.eligible = true;
  }
}

void measure_group(std::vector<RowResult> &rows, const Arguments &arguments) {
  const Case &prototype = rows.front().benchmark_case;
  const OperatorData op = make_operator(prototype);
  const State initial = make_initial_state(prototype);
  const std::vector<State> inputs(
      static_cast<std::size_t>(prototype.batch_size), initial);

  for (RowResult &row : rows) {
    validate_row(row);
  }
  for (int warmup = 0; warmup < arguments.warmups; ++warmup) {
    for (std::size_t offset = 0; offset < rows.size(); ++offset) {
      RowResult &row = rows[(offset + static_cast<std::size_t>(warmup)) %
                            rows.size()];
      if (row.validation_ok) {
        const auto measured = run_batch(row.benchmark_case, inputs, op, false,
                                        false);
        if (!measured.ok) {
          throw std::runtime_error("lean warmup changed status: " +
                                   measured.message);
        }
        combine_checksum(row.checksum, measured.checksum);
      }
    }
  }

  std::vector<bool> stopped(rows.size(), false);
  for (std::size_t index = 0; index < rows.size(); ++index) {
    if (!rows[index].validation_ok) {
      stopped[index] = true;
      rows[index].stop_reason = "invalid";
    }
  }
  for (int repetition = 0; repetition < arguments.maximum_repetitions;
       ++repetition) {
    bool any_active = false;
    for (std::size_t offset = 0; offset < rows.size(); ++offset) {
      const std::size_t index =
          (offset + static_cast<std::size_t>(repetition)) % rows.size();
      if (stopped[index]) {
        continue;
      }
      any_active = true;
      const auto measured = run_batch(rows[index].benchmark_case, inputs, op,
                                      false, false);
      if (!measured.ok) {
        throw std::runtime_error("lean measurement changed status: " +
                                 measured.message);
      }
      rows[index].samples.push_back(measured.seconds_per_column);
      combine_checksum(rows[index].checksum, measured.checksum);
    }
    if (!any_active) {
      break;
    }
    for (std::size_t index = 0; index < rows.size(); ++index) {
      if (stopped[index] ||
          static_cast<int>(rows[index].samples.size()) <
              arguments.minimum_repetitions) {
        continue;
      }
      const auto interval = ndrop_test::benchmark::bootstrap_median_interval(
          rows[index].samples, arguments.bootstrap_resamples,
          UINT64_C(0x4e44524f) ^
              static_cast<std::uint64_t>(rows[index].benchmark_case.id));
      if (interval.relative_half_width < 0.05) {
        stopped[index] = true;
        rows[index].primary_interval = interval;
        rows[index].stop_reason = "ci_target";
      } else if (static_cast<int>(rows[index].samples.size()) >=
                 arguments.maximum_repetitions) {
        stopped[index] = true;
        rows[index].primary_interval = interval;
        rows[index].stop_reason = "maximum_repetitions";
      }
    }
  }

  std::size_t maximum_instrumented_repetitions = 0;
  for (RowResult &row : rows) {
    if (row.validation_ok && row.stop_reason.empty()) {
      row.stop_reason = "maximum_repetitions";
      row.primary_interval =
          ndrop_test::benchmark::bootstrap_median_interval(
              row.samples, arguments.bootstrap_resamples,
              UINT64_C(0x4e44524f) ^
                  static_cast<std::uint64_t>(row.benchmark_case.id));
    }
    maximum_instrumented_repetitions =
        std::max(maximum_instrumented_repetitions, row.samples.size());
  }
  for (std::size_t repetition = 0;
       repetition < maximum_instrumented_repetitions; ++repetition) {
    for (std::size_t offset = 0; offset < rows.size(); ++offset) {
      const std::size_t index = (offset + repetition) % rows.size();
      RowResult &row = rows[index];
      if (!row.validation_ok || repetition >= row.samples.size()) {
        continue;
      }
      const auto measured = run_batch(row.benchmark_case, inputs, op, true,
                                      true);
      if (!measured.ok) {
        throw std::runtime_error("instrumented measurement changed status: " +
                                 measured.message);
      }
      row.instrumented_samples.push_back(measured.seconds_per_column);
      row.fence_samples.push_back(measured.fence_seconds_per_column);
      row.phase_samples[0].push_back(measured.phases_per_column.setup);
      row.phase_samples[1].push_back(
          measured.phases_per_column.factorization);
      row.phase_samples[2].push_back(measured.phases_per_column.solve);
      row.phase_samples[3].push_back(measured.phases_per_column.update);
      row.phase_samples[4].push_back(measured.phases_per_column.diagnostics);
      combine_checksum(row.checksum, measured.checksum);
    }
  }
  for (RowResult &row : rows) {
    if (!row.validation_ok) {
      continue;
    }
    row.instrumented_interval =
        ndrop_test::benchmark::bootstrap_median_interval(
            row.instrumented_samples, arguments.bootstrap_resamples,
            UINT64_C(0x4e44524f) ^ UINT64_C(0x494e535452) ^
                static_cast<std::uint64_t>(row.benchmark_case.id));
    for (std::size_t phase = 0; phase < row.phase_samples.size(); ++phase) {
      row.phase_medians[phase] =
          ndrop_test::benchmark::sample_median(row.phase_samples[phase]);
    }
    row.fence_median =
        ndrop_test::benchmark::sample_median(row.fence_samples);
    row.instrumentation_overhead =
        row.primary_interval.median > 0
            ? row.instrumented_interval.median /
                  row.primary_interval.median -
                  1
            : 0;
  }
}

std::string compiler_name() {
#if defined(__clang__)
  return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
  return std::string("gcc ") + __VERSION__;
#else
  return "unknown";
#endif
}

std::string platform_name() {
#if defined(__APPLE__) && defined(__aarch64__)
  return "macos-arm64";
#elif defined(__APPLE__)
  return "macos";
#else
  return "unknown";
#endif
}

void write_header(std::ofstream &output) {
  output << "case_id,nlev,active_pairs,transported_fields,batch_size,"
            "mixing_index,mixing_scale,dt,method,forced_fe_substeps,"
            "validation_ok,eligible,exclusion_reason,message,repetitions,"
            "stop_reason,median_seconds_per_column,ci95_lower,ci95_upper,"
            "ci95_relative_half_width,instrumented_median_seconds_per_column,"
            "instrumented_ci95_lower,instrumented_ci95_upper,"
            "instrumentation_overhead,setup_seconds_per_column,"
            "factorization_seconds_per_column,solve_seconds_per_column,"
            "update_seconds_per_column,diagnostics_seconds_per_column,"
            "fence_seconds_per_column,clipping_count,material_negative_count,"
            "block_factorizations,block_rhs_solves,scalar_factorizations,"
            "scalar_rhs_solves,executed_field_level_updates,workspace_bytes,"
            "checksum,compiler,platform,kokkos_execution_space,hardware_threads,"
            "warmups,minimum_repetitions,maximum_repetitions,"
            "bootstrap_resamples,raw_samples_seconds_per_column,"
            "instrumented_samples_seconds_per_column\n";
}

void write_row(std::ofstream &output, const RowResult &row,
               const Arguments &arguments) {
  const Case &item = row.benchmark_case;
  const std::uint64_t updates =
      UINT64_C(51) * item.nlev *
      static_cast<std::uint64_t>(item.method == Method::forward_euler
                                     ? item.forced_fe_substeps
                                     : 1) *
      static_cast<std::uint64_t>(item.batch_size);
  auto field = [&output](const std::string &value) {
    output << ',' << csv_escape(value);
  };
  auto timing = [&row](const Real value) {
    return row.samples.empty() ? std::string{} : real_text(value);
  };
  output << item.id << ',' << item.nlev << ',' << item.active_pairs << ','
         << 1 + 2 * item.active_pairs << ',' << item.batch_size << ','
         << item.mixing_index << ',' << real_text(item.mixing_scale) << ','
         << real_text(item.dt);
  field(ndrop_test::benchmark::method_name(item.method));
  output << ',' << item.forced_fe_substeps << ','
         << (row.validation_ok ? "true" : "false") << ','
         << (row.eligible ? "true" : "false");
  field(row.exclusion_reason);
  field(row.message);
  output << ',' << row.samples.size();
  field(row.stop_reason);
  output << ',' << timing(row.primary_interval.median) << ','
         << timing(row.primary_interval.lower) << ','
         << timing(row.primary_interval.upper) << ','
         << timing(row.primary_interval.relative_half_width) << ','
         << timing(row.instrumented_interval.median) << ','
         << timing(row.instrumented_interval.lower) << ','
         << timing(row.instrumented_interval.upper) << ','
         << timing(row.instrumentation_overhead) << ','
         << timing(row.phase_medians[0]) << ','
         << timing(row.phase_medians[1]) << ','
         << timing(row.phase_medians[2]) << ','
         << timing(row.phase_medians[3]) << ','
         << timing(row.phase_medians[4]) << ','
         << timing(row.fence_median) << ','
         << row.validation_diagnostics.clipping_count << ','
         << row.validation_diagnostics.material_negative_count << ','
         << row.validation_diagnostics.block_factorization_count << ','
         << row.validation_diagnostics.block_rhs_solve_count << ','
         << row.validation_diagnostics.scalar_factorization_count << ','
         << row.validation_diagnostics.scalar_rhs_solve_count << ','
         << updates << ','
         << ndrop_test::benchmark::estimated_workspace_bytes(item);
  field(hex_text(row.checksum));
  field(compiler_name());
  field(platform_name());
  field(Kokkos::DefaultExecutionSpace::name());
  output << ',' << std::thread::hardware_concurrency() << ','
         << arguments.warmups << ',' << arguments.minimum_repetitions << ','
         << arguments.maximum_repetitions << ','
         << arguments.bootstrap_resamples;
  field(joined_samples(row.samples));
  field(joined_samples(row.instrumented_samples));
  output << '\n';
}

} // namespace

int main(int argc, char **argv) {
  try {
    const Arguments arguments = parse_arguments(argc, argv);
    Kokkos::initialize(argc, argv);
    try {
      const auto groups = selected_groups(arguments);
      std::ofstream output(arguments.output_path);
      if (!output) {
        throw std::runtime_error("cannot open benchmark output: " +
                                 arguments.output_path);
      }
      write_header(output);
      std::size_t rows_written = 0;
      std::size_t invalid_rows = 0;
      for (std::size_t group_index = 0; group_index < groups.size();
           ++group_index) {
        std::vector<RowResult> rows;
        rows.reserve(groups[group_index].size());
        for (const Case &benchmark_case : groups[group_index]) {
          RowResult row;
          row.benchmark_case = benchmark_case;
          rows.push_back(std::move(row));
        }
        measure_group(rows, arguments);
        for (const RowResult &row : rows) {
          write_row(output, row, arguments);
          ++rows_written;
          if (!row.validation_ok) {
            ++invalid_rows;
          }
        }
        output.flush();
        if (!output) {
          throw std::runtime_error("benchmark output write failed");
        }
        if (arguments.progress_every > 0 &&
            ((group_index + 1) %
                 static_cast<std::size_t>(arguments.progress_every) ==
             0 ||
             group_index + 1 == groups.size())) {
          std::cerr << "PROGRESS groups=" << group_index + 1 << '/'
                    << groups.size() << " rows=" << rows_written << '\n';
        }
      }
      std::cout << "RESULT rows=" << rows_written
                << " invalid=" << invalid_rows << " failed=0\n";
    } catch (...) {
      Kokkos::finalize();
      throw;
    }
    Kokkos::finalize();
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
