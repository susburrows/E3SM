cmake_minimum_required(VERSION 3.18)

project(ndrop_integrator_gate LANGUAGES CXX)

include(CTest)
find_package(Kokkos CONFIG REQUIRED)

get_filename_component(E3SM_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../.."
                       ABSOLUTE)
set(EKAT_SOURCE "${E3SM_SOURCE_ROOT}/externals/ekat" CACHE PATH
    "Path to the locked EKAT source checkout")
set(MAM4XX_SOURCE "${E3SM_SOURCE_ROOT}/externals/mam4xx" CACHE PATH
    "Path to the locked MAM4xx source checkout")

option(NDROP_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan" OFF)
option(NDROP_BUILD_EXACT_ADAPTER
       "Build the exact-current MAM4xx Forward Euler comparison" ON)

function(configure_ndrop_target target)
  target_compile_features(${target} PRIVATE cxx_std_17)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  target_link_libraries(${target} PRIVATE Kokkos::kokkos)
  if(NDROP_ENABLE_SANITIZERS)
    target_compile_options(
      ${target} PRIVATE -O0 -g -fsanitize=address,undefined
                        -fno-omit-frame-pointer)
    target_link_options(
      ${target} PRIVATE -fsanitize=address,undefined
                        -fno-omit-frame-pointer)
  endif()
endfunction()

add_executable(test_forward_euler ndrop_integrator.cpp test_forward_euler.cpp)
configure_ndrop_target(test_forward_euler)
add_test(NAME ndrop_forward_euler COMMAND test_forward_euler)

add_executable(test_imex ndrop_integrator.cpp ndrop_imex.cpp test_imex.cpp)
target_include_directories(
  test_imex SYSTEM PRIVATE
  "${EKAT_SOURCE}/src/algorithm"
  "${EKAT_SOURCE}/src/core"
  "${EKAT_SOURCE}/src/kokkos"
  "${EKAT_SOURCE}/src/pack")
configure_ndrop_target(test_imex)
add_test(NAME ndrop_imex COMMAND test_imex)

add_executable(test_backward_euler ndrop_integrator.cpp
               ndrop_backward_euler.cpp test_backward_euler.cpp)
configure_ndrop_target(test_backward_euler)
add_test(NAME ndrop_backward_euler COMMAND test_backward_euler)

add_executable(test_cross_scheme ndrop_integrator.cpp ndrop_imex.cpp
               ndrop_backward_euler.cpp ndrop_reference.cpp
               test_cross_scheme.cpp)
target_include_directories(
  test_cross_scheme SYSTEM PRIVATE
  "${EKAT_SOURCE}/src/algorithm"
  "${EKAT_SOURCE}/src/core"
  "${EKAT_SOURCE}/src/kokkos"
  "${EKAT_SOURCE}/src/pack")
configure_ndrop_target(test_cross_scheme)
add_test(NAME ndrop_cross_scheme COMMAND test_cross_scheme)

add_executable(test_benchmark_support ndrop_integrator.cpp ndrop_imex.cpp
               ndrop_backward_euler.cpp ndrop_benchmark.cpp
               test_benchmark_support.cpp)
target_include_directories(
  test_benchmark_support SYSTEM PRIVATE
  "${EKAT_SOURCE}/src/algorithm"
  "${EKAT_SOURCE}/src/core"
  "${EKAT_SOURCE}/src/kokkos"
  "${EKAT_SOURCE}/src/pack")
configure_ndrop_target(test_benchmark_support)
add_test(NAME ndrop_benchmark_support COMMAND test_benchmark_support)

# OI-5 offline unit tests use synthetic/repository inputs and never invoke a
# production capture or the timed benchmark executable.
add_executable(test_ndrop_replay_schema ndrop_integrator.cpp ndrop_replay_schema.cpp
               test_ndrop_replay_schema.cpp)
configure_ndrop_target(test_ndrop_replay_schema)
add_test(NAME ndrop_oi5_schema COMMAND test_ndrop_replay_schema)

add_executable(test_ndrop_production_replay ndrop_integrator.cpp ndrop_imex.cpp
               ndrop_backward_euler.cpp ndrop_reference.cpp
               ndrop_production_replay.cpp ndrop_replay_schema.cpp
               test_ndrop_production_replay.cpp)
target_include_directories(
  test_ndrop_production_replay SYSTEM PRIVATE
  "${EKAT_SOURCE}/src/algorithm"
  "${EKAT_SOURCE}/src/core"
  "${EKAT_SOURCE}/src/kokkos"
  "${EKAT_SOURCE}/src/pack")
configure_ndrop_target(test_ndrop_production_replay)
add_test(NAME ndrop_oi5_production_replay_unit
         COMMAND test_ndrop_production_replay)

set(NDROP_OI5_BASELINE_DIR "" CACHE PATH
    "Exact pre-change offline source copy for OI-5 comparison")
if(NDROP_OI5_BASELINE_DIR)
  foreach(_file ndrop_integrator.hpp ndrop_integrator.cpp ndrop_imex.cpp
                ndrop_backward_euler.cpp ndrop_reference.cpp)
    if(NOT EXISTS "${NDROP_OI5_BASELINE_DIR}/${_file}")
      message(FATAL_ERROR "OI-5 baseline is missing ${_file}")
    endif()
  endforeach()
  add_executable(test_ndrop_oi5_legacy_baseline
                 "${NDROP_OI5_BASELINE_DIR}/ndrop_integrator.cpp"
                 "${NDROP_OI5_BASELINE_DIR}/ndrop_imex.cpp"
                 "${NDROP_OI5_BASELINE_DIR}/ndrop_backward_euler.cpp"
                 "${NDROP_OI5_BASELINE_DIR}/ndrop_reference.cpp"
                 test_ndrop_production_replay.cpp)
  target_compile_definitions(test_ndrop_oi5_legacy_baseline PRIVATE
      NDROP_LEGACY_BASELINE
      NDROP_LEGACY_BASELINE_HEADER="${NDROP_OI5_BASELINE_DIR}/ndrop_integrator.hpp")
  target_include_directories(
    test_ndrop_oi5_legacy_baseline SYSTEM PRIVATE
    "${EKAT_SOURCE}/src/algorithm"
    "${EKAT_SOURCE}/src/core"
    "${EKAT_SOURCE}/src/kokkos"
    "${EKAT_SOURCE}/src/pack")
  configure_ndrop_target(test_ndrop_oi5_legacy_baseline)
endif()

add_executable(benchmark_ndrop_integrators ndrop_integrator.cpp
               ndrop_imex.cpp ndrop_backward_euler.cpp ndrop_benchmark.cpp
               benchmark_ndrop_integrators.cpp)
target_include_directories(
  benchmark_ndrop_integrators SYSTEM PRIVATE
  "${EKAT_SOURCE}/src/algorithm"
  "${EKAT_SOURCE}/src/core"
  "${EKAT_SOURCE}/src/kokkos"
  "${EKAT_SOURCE}/src/pack")
configure_ndrop_target(benchmark_ndrop_integrators)

# Separate accuracy-first study. It is never registered as a CTest test: CTest
# may run numerical unit tests, but must not accidentally start a timed study.
add_executable(benchmark_ndrop_matched_error ndrop_integrator.cpp
               ndrop_imex.cpp ndrop_backward_euler.cpp ndrop_reference.cpp
               ndrop_benchmark.cpp benchmark_ndrop_matched_error.cpp)
target_include_directories(
  benchmark_ndrop_matched_error SYSTEM PRIVATE
  "${EKAT_SOURCE}/src/algorithm"
  "${EKAT_SOURCE}/src/core"
  "${EKAT_SOURCE}/src/kokkos"
  "${EKAT_SOURCE}/src/pack")
configure_ndrop_target(benchmark_ndrop_matched_error)
add_test(NAME ndrop_controlled_mixing_fixtures
         COMMAND benchmark_ndrop_matched_error --check-controlled-mixing-fixtures)
add_test(NAME ndrop_physical_column_fixtures
         COMMAND benchmark_ndrop_matched_error --check-physical-column-fixtures
                 --physical-grid-dir "${E3SM_SOURCE_ROOT}/components/homme/test/vcoord")
add_test(NAME ndrop_physical_geometry_fixtures
         COMMAND benchmark_ndrop_matched_error --check-physical-geometry-fixtures
                 --physical-grid-dir "${E3SM_SOURCE_ROOT}/components/homme/test/vcoord")
add_test(
  NAME ndrop_benchmark_smoke
  COMMAND benchmark_ndrop_integrators --smoke --output
          "${CMAKE_CURRENT_BINARY_DIR}/benchmark_smoke.csv"
          --run-id task5-remediated-smoke
          --phase-samples-output
          "${CMAKE_CURRENT_BINARY_DIR}/benchmark_smoke_phases.csv")
add_test(
  NAME ndrop_benchmark_exclusion
  COMMAND benchmark_ndrop_integrators --output
          "${CMAKE_CURRENT_BINARY_DIR}/benchmark_exclusion.csv"
          --run-id task5-remediated-exclusion
          --phase-samples-output
          "${CMAKE_CURRENT_BINARY_DIR}/benchmark_exclusion_phases.csv"
          --case-start 24291 --case-count 9 --warmups 0
          --min-repetitions 1 --max-repetitions 1
          --bootstrap-resamples 10 --progress-every 0)
function(add_ndrop_expected_failure name expected)
  string(JOIN "|" encoded_args ${ARGN})
  add_test(
    NAME ${name}
    COMMAND "${CMAKE_COMMAND}"
            "-DPROGRAM=$<TARGET_FILE:benchmark_ndrop_integrators>"
            "-DEXPECTED=${expected}"
            "-DENCODED_ARGS=${encoded_args}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/expect_failure.cmake")
endfunction()

add_ndrop_expected_failure(
  ndrop_benchmark_requires_run_id "--run-id is required"
  --smoke --output "${CMAKE_CURRENT_BINARY_DIR}/missing_run_id.csv"
  --phase-samples-output
  "${CMAKE_CURRENT_BINARY_DIR}/missing_run_id_phases.csv")
add_ndrop_expected_failure(
  ndrop_benchmark_requires_phase_output "--phase-samples-output is required"
  --smoke --output "${CMAKE_CURRENT_BINARY_DIR}/missing_phase_output.csv"
  --run-id task5-remediated-missing-phase)
add_ndrop_expected_failure(
  ndrop_benchmark_rejects_invalid_run_id
  "--run-id must match [A-Za-z0-9._-]{1,128}"
  --smoke --output "${CMAKE_CURRENT_BINARY_DIR}/invalid_run_id.csv"
  --run-id invalid/run --phase-samples-output
  "${CMAKE_CURRENT_BINARY_DIR}/invalid_run_id_phases.csv")
add_ndrop_expected_failure(
  ndrop_benchmark_phase_output_failure "cannot open phase-samples output"
  --smoke --output "${CMAKE_CURRENT_BINARY_DIR}/phase_output_failure.csv"
  --run-id task5-remediated-phase-failure
  --phase-samples-output "${CMAKE_CURRENT_BINARY_DIR}")
add_ndrop_expected_failure(
  ndrop_benchmark_rejects_aliased_outputs
  "--output and --phase-samples-output must be different paths"
  --smoke --output "${CMAKE_CURRENT_BINARY_DIR}/aliased_output.csv"
  --run-id task5-remediated-alias --phase-samples-output
  "${CMAKE_CURRENT_BINARY_DIR}/./aliased_output.csv")

if(NDROP_BUILD_EXACT_ADAPTER)
  if(NOT EXISTS "${MAM4XX_SOURCE}/src/mam4xx/mam4_config.hpp.in")
    message(FATAL_ERROR
            "Current MAM4xx source not found at ${MAM4XX_SOURCE}; prepare "
            "the locked stack or set -DMAM4XX_SOURCE=/path/to/mam4xx")
  endif()

  set(MAM4XX_REAL_TYPE double)
  set(MAM4XX_DOUBLE_PRECISION 1)
  set(MAM4XX_ENABLE_GPU OFF)
  set(MAM4XX_NUM_VERTICAL_LEVELS 128)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated/mam4xx")
  configure_file(
    "${MAM4XX_SOURCE}/src/mam4xx/mam4_config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/mam4xx/mam4_config.hpp"
    @ONLY)

  add_executable(compare_mam4xx_forward_euler
                 ndrop_integrator.cpp compare_mam4xx_forward_euler.cpp)
  target_include_directories(compare_mam4xx_forward_euler
                             PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
  target_include_directories(
    compare_mam4xx_forward_euler SYSTEM PRIVATE
    "${MAM4XX_SOURCE}/src"
    "${EKAT_SOURCE}/src/algorithm"
    "${EKAT_SOURCE}/src/core"
    "${EKAT_SOURCE}/src/kokkos"
    "${EKAT_SOURCE}/src/pack")
  configure_ndrop_target(compare_mam4xx_forward_euler)
  add_test(NAME ndrop_mam4xx_exact_forward_euler
           COMMAND compare_mam4xx_forward_euler)
endif()
