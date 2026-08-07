cmake_minimum_required(VERSION 3.18)

project(ndrop_integrator_gate LANGUAGES CXX)

include(CTest)
find_package(Kokkos CONFIG REQUIRED)

get_filename_component(E3SM_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../.."
                       ABSOLUTE)
set(EKAT_SOURCE "${E3SM_SOURCE_ROOT}/externals/ekat")
set(MAM4XX_SOURCE "${E3SM_SOURCE_ROOT}/externals/mam4xx")
set(HAERO_SOURCE "${E3SM_SOURCE_ROOT}/externals/haero" CACHE PATH
    "Path to the Haero source checkout")

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

if(NDROP_BUILD_EXACT_ADAPTER)
  if(NOT EXISTS "${HAERO_SOURCE}/haero/haero_config.hpp.in")
    message(FATAL_ERROR
            "Haero source not found at ${HAERO_SOURCE}; initialize the "
            "submodule or set -DHAERO_SOURCE=/path/to/haero")
  endif()

  set(HAERO_REAL_TYPE double)
  set(HAERO_DOUBLE_PRECISION 1)
  set(HAERO_ENABLE_GPU OFF)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated/haero")
  configure_file(
    "${HAERO_SOURCE}/haero/haero_config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/haero/haero_config.hpp")

  set(NUM_VERTICAL_LEVELS 128)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated/mam4xx")
  configure_file(
    "${MAM4XX_SOURCE}/src/mam4xx/aero_config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/mam4xx/aero_config.hpp")

  add_executable(compare_mam4xx_forward_euler
                 ndrop_integrator.cpp compare_mam4xx_forward_euler.cpp)
  target_compile_definitions(compare_mam4xx_forward_euler
                             PRIVATE NUM_VERTICAL_LEVELS=128)
  target_include_directories(compare_mam4xx_forward_euler
                             PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
  target_include_directories(
    compare_mam4xx_forward_euler SYSTEM PRIVATE
    "${MAM4XX_SOURCE}/src"
    "${HAERO_SOURCE}"
    "${EKAT_SOURCE}/src/algorithm"
    "${EKAT_SOURCE}/src/core"
    "${EKAT_SOURCE}/src/kokkos"
    "${EKAT_SOURCE}/src/pack")
  configure_ndrop_target(compare_mam4xx_forward_euler)
  add_test(NAME ndrop_mam4xx_exact_forward_euler
           COMMAND compare_mam4xx_forward_euler)
endif()
