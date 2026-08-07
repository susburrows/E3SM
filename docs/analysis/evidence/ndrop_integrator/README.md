# Isolated `ndrop` Integrator Harness

This directory contains test-only code for the approved standalone numerical
experiment in `docs/design/eamxx_mam4xx_ndrop.md`. It is not linked into E3SM
or MAM4xx production targets.

The first gate implements the canonical shared-interface operator and the
current Forward Euler subcycling, activation-rate limiting, alternating work
buffers, and componentwise positivity clipping. The test suite covers all 25
aerosol constituents and droplet number.

Build and run outside the OneDrive source tree:

```bash
SRC=/path/to/E3SM/docs/analysis/evidence/ndrop_integrator
BUILD=/private/tmp/ndrop-int
mkdir -p "$BUILD"

clang++ -std=c++17 -O0 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  "$SRC/ndrop_integrator.cpp" "$SRC/test_forward_euler.cpp" \
  -o "$BUILD/test_forward_euler"

"$BUILD/test_forward_euler"
```

Successful harness tests do not validate the known production `top_lev=0` or
`enable_aero_vertical_mix=false` defect paths. Those paths are explicitly
excluded by the approved design.

`compare_mam4xx_forward_euler.cpp` is the separate exact-current-function
adapter. It compares the canonical harness with the unmodified
`mam4::ndrop::update_from_explmix` at `NUM_VERTICAL_LEVELS=128` and a safe
`top_lev=1`. It covers production substep counts 1, 2, 3, and 5; full and
fractional overlap; limiter-active and limiter-inactive activation; and
isolated, contiguous, and all-zero cloud profiles. The destination work array
is filled with unrelated values before every exact call to detect incomplete
writes.

`ndrop_imex.cpp` implements the approved IMEX split with EKAT's serial Thomas
solver. `test_imex.cpp` checks identity and uniform-state behavior, fractional-
overlap ledgers, residuals, fail-loud positivity behavior, and first-order
convergence against discrete no-flux diffusion and closed-form two-level
activation-plus-mixing solutions.

## Solver structure and scope

EKAT is the E3SM Kokkos Application Toolkit. The harness uses its scalar
tridiagonal Thomas solver; it does not implement a new general-purpose linear
solver. With the harness coefficients frozen, the tested operator is linear.
For each aerosol pair, the complete IMEX system is block lower triangular: the
cloud-borne state is solved first with one scalar tridiagonal matrix, its
cloud-to-interstitial return is added to the next right-hand side, and the
interstitial state is solved with a second scalar tridiagonal matrix. Droplet
number shares the cloud-borne matrix family.

The `ndrop`-specific IMEX source evaluation, matrix/right-hand-side assembly,
solve order, zero-cloud transition, positivity checks, and ledger diagnostics
are implemented in the standalone harness. SUNDIALS and PETSc are not used.
The harness also does not implement EKAT's GPU cyclic-reduction path, Backward
Euler's coupled block-tridiagonal solve, a hybrid FE/IMEX policy, or any
higher-order implicit method. Those are separate, unvalidated stages and are
not implied by the CPU correctness results below.

## Remediation validation

An independent review on 2026-07-29 rejected the original Forward Euler and
IMEX gates. The shared harness omitted the production zero-cloud
resuspension/reset, could accept nonfinite diagnostics after finite-input
overflow, and did not contain all of the approved all-field convergence and
clipping checks. The 2026-07-30 remediation adds:

- exact production zero-cloud handling after each complete integration call;
- checked derived coefficients, activation sources, state updates, matrix and
  right-hand-side values, solutions, residuals, ledgers, and normalizations;
- scaled matrix residuals and a separately assembled physical-equation
  residual that does not reuse the tridiagonal matrix;
- droplet and per-constituent interstitial/cloud-borne raw minima, clipping
  counts, weighted clipping changes, and ledger effects;
- all-51-field analytic diffusion and activation convergence;
- a nonuniform, fractional-overlap, all-field coupled FE/IMEX
  self-convergence and common-limit test; and
- randomized weak-to-strong production-adaptive FE tests confirming that the
  approved substep and activation-rate limits keep raw updates nonnegative.

The last property follows from the production substep bound: the retained
coefficient of the old value is positive, mixing contributions from neighboring
levels are nonnegative, and the limited activation removal cannot exceed its
matching positive interstitial contribution from below. The forced unstable
step remains only to test clipping diagnostics.

### Executed results

The final sanitizer and optimized Forward Euler runs each passed 20,899 checks:

```text
CONVERGENCE case=fe_diffusion order_range=1.00004:1.00004 min_r2=1
CONVERGENCE case=fe_activation_all_fields order_range=1.00555:1.00702 min_r2=0.999998
RESULT passed=20899 failed=0
```

The final sanitizer and optimized IMEX runs each passed 21,537 checks:

```text
CONVERGENCE case=diffusion l1_order_range=0.999957:0.999957 linf_order_range=0.999957:0.999958 min_r2=1
CONVERGENCE case=activation l1_order_range=0.992285:0.992673 linf_order_range=0.993273:0.993517 min_r2=0.999998
CONVERGENCE case=coupled_all_fields order_range=0.99701:1.00319 max_finest_cross_error=9.0418e-05
RESULT passed=21537 failed=0
```

The exact-current adapter passed nine sanitizer and optimized cases with zero
mismatches and zero maximum absolute or relative difference. This includes
isolated, contiguous, and all-zero cloud profiles.

The commands were run from the E3SM source root with all generated executables
under `/private/tmp/ndrop-remediate`. The debug Kokkos/EKAT installation was
`/private/tmp/ni2-install`; the release installation was
`/private/tmp/ni2-rel-install`. Representative build forms were:

```bash
clang++ -std=c++17 -O0 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer \
  docs/analysis/evidence/ndrop_integrator/ndrop_integrator.cpp \
  docs/analysis/evidence/ndrop_integrator/test_forward_euler.cpp \
  -o /private/tmp/ndrop-remediate/test_forward_euler_san

clang++ -std=c++17 -O0 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer \
  -isystem /private/tmp/ni2-install/include/ekat \
  -isystem /private/tmp/ni2-install/include/kokkos \
  docs/analysis/evidence/ndrop_integrator/ndrop_integrator.cpp \
  docs/analysis/evidence/ndrop_integrator/ndrop_imex.cpp \
  docs/analysis/evidence/ndrop_integrator/test_imex.cpp \
  /private/tmp/ni2-install/lib/libekat_kokkosutils.a \
  /private/tmp/ni2-install/lib/libkokkoscontainers.a \
  /private/tmp/ni2-install/lib/libkokkoscore.a -ldl -lpthread \
  -o /private/tmp/ndrop-remediate/test_imex_san
```

The optimized forms replace the sanitizer/debug flags with `-O3 -DNDEBUG` and
use `/private/tmp/ni2-rel-install`. The exact-current adapter additionally uses
`-DNUM_VERTICAL_LEVELS=128`, the checked-out MAM4xx headers, and the generated
MAM4xx/Haero configuration headers.

The debug executables used ASan/UBSan and Kokkos bounds checking. The optimized
executables used `-O3 -DNDEBUG` and the release Kokkos/EKAT libraries. Both
configurations use double precision, Apple clang 21.0.0 on arm64 macOS, Kokkos
4.5.1 Serial, and EKAT's serial Thomas solver. Compilation with
`-Wall -Wextra -Wpedantic -Werror` was clean for harness sources; EKAT/Kokkos
headers were treated as system headers for this check.

### Provenance

- E3SM inspected revision:
  `e12e7d61dba924717e8f35fe8288247c1766ee30`
- MAM4xx inspected revision:
  `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`
- EKAT inspected revision:
  `65cda3ff61abe22add622b56cee6234e29fc51c4`
- Debug Kokkos library SHA-256:
  `a040f738b29bbfb6471486fef4a7096df0e8697667f5be965b83adb51711055e`
- Debug EKAT Kokkos-utils library SHA-256:
  `a0338864d8afecd0cb6538a8e7dace8cb0d4f0a8de6cf15f9d5f6f9a4c181f24`
- Release Kokkos library SHA-256:
  `a5246e431826d861593256c42101828c724a1c8d7a0035965bd5eca814b9e402`
- Release EKAT Kokkos-utils library SHA-256:
  `d3a2f9b15dfe08542f5eb3bd92cef3b5ebf13509a69d8fa9b46880c84f3adbcc`

Retained output records:

| Record | SHA-256 |
|---|---|
| `fe_validation_20260730.txt` | `9e2259c612bae542b0d6c2341807b52a974f2267e83fd5a42ed7cc7ba28de8d6` |
| `imex_validation_20260730.txt` | `9b164825cb3882905963003afa49f0ce66ab9305967f832a9569b2811c20be22` |
| `mam4xx_comparison_20260730.txt` | `04c55fd1dbcbbe7a945152b0a94fd1907240500ad60fbaa0255ab19778b1a870` |

The inspected EKAT revision's macOS floating-point-exception fallback accesses
x86-only `fenv_t` members and does not compile on Apple Silicon. The generated
scratch dependency copy used no-op FPE mask-control functions for the local
build. These test executables do not enable or query FPE traps. No repository
EKAT, MAM4xx, or production physics source was modified by this workaround.

This is local CPU correctness evidence for the test-only operator. It is not
GPU performance evidence, whole-model validation, or authorization for a
production IMEX implementation. The known `top_lev=0` and
`enable_aero_vertical_mix=false` defect paths remain excluded.
