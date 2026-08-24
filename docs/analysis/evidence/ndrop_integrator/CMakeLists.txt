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
add_test(
  NAME ndrop_benchmark_smoke
  COMMAND benchmark_ndrop_integrators --smoke --output
          "${CMAKE_CURRENT_BINARY_DIR}/benchmark_smoke.csv")
add_test(
  NAME ndrop_benchmark_exclusion
  COMMAND benchmark_ndrop_integrators --output
          "${CMAKE_CURRENT_BINARY_DIR}/benchmark_exclusion.csv"
          --case-start 24291 --case-count 9 --warmups 0
          --min-repetitions 1 --max-repetitions 1
          --bootstrap-resamples 10 --progress-every 0)

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
