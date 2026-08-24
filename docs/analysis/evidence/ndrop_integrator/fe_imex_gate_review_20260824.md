# Forward Euler and IMEX Task 2 gate review

**Execution date:** 2026-08-24

**Code/test commits:** `e67097b6bc29e919329a729d48776f866a32038e`,
`07c786d82d7571fd01527b26cc3c21de09e09acd`

**Final code/test tree:** `7bef91ee45026912c41f4ae605aee0758fa49f39`

**Scope:** standalone test-only FE and IMEX requalification; no production
EAMxx or MAM4xx physics or gitlink changes

## Decision summary

| Gate | Decision | Executed evidence |
|---|---|---|
| Forward Euler Task 2 gate | **PASS** | 23,195 checks in both debug/sanitizer and optimized builds; analytic slopes and `R^2` satisfy the approved first-order gate; signed ledgers, independent per-substep residuals, clipping accounting, scaling, zero-ledger, finite-failure, negativity-classification, and three-run bitwise-state tests pass. |
| IMEX Task 2 gate | **PASS** | 25,138 checks in both builds; analytic and coupled convergence, matrix and independently assembled physical residuals, signed ledgers, scaling, zero-ledger, finite-failure, material-negativity rejection, and three-run bitwise-state tests pass. |
| Locked current-MAM4xx adapter | **PASS** | 15 cases pass in both builds with exact `nsubmix`, buffer, droplet, and all aerosol field agreement; maximum absolute and relative differences are zero. |
| Archived MAM4xx `1b02f272...` at `top_lev=0` | **EXCLUDED expected failure** | Historical evidence only; it is not built or registered in the current gate. |

These are Task 2 method-local decisions. They do not replace the complete
Task 4 common accuracy matrix, authorize production integration, or constitute
GPU or whole-model evidence.

## Locked inputs and environment

The dependency preparation script reverified four clean detached checkouts
before the retained gate logs were written.

| Input | Revision |
|---|---|
| E3SM baseline | `e2c578e40f6602a108eff577091de86a453624ab` |
| MAM4xx | `0d076d8ebb004a18f05f51f97b08edd9def6826f` |
| EKAT | `3a1485aecfed7a036e4fa30ce845498fd37975e2` |
| Kokkos | `23745d8494203711a6090d0a201c3ae3bf57853b` |

Both builds used Apple clang 21.0.0, CMake 4.4.2, Ninja, double precision,
and Kokkos 4.5.1 Serial on an arm64 Apple M2 Max running macOS 26.6.1. The
debug harness enabled ASan and UBSan with a Kokkos debug build. Kokkos's
separate debug-bounds option was not enabled. The optimized harness and Kokkos
dependency used Release
`-O3 -DNDEBUG`. All harness targets compiled with
`-Wall -Wextra -Wpedantic -Werror`.

## As-built equation and interface audit

- `first_level_controller_rate` has units s^-1, is finite and nonnegative,
  and is added only to the first active level's FE inverse-timescale reduction.
  A forced-one-substep mutation produces bitwise-equal physical FE states with
  the controller enabled and disabled. IMEX matrices and both physical
  operator helpers do not reference it.
- For each aerosol pair, the implementation computes
  `total = final - initial`, `operator = total - clipping`, and checks every
  subtraction for a finite representable result. IMEX clipping is exactly
  zero. FE clipping is the accumulated pressure-weighted raw-to-clipped change
  at each substep.
- The FE physical residual is assembled from the old clipped state, source,
  shared physical interface rates, and overlap terms before the raw state is
  clipped. The IMEX physical residual is assembled independently of the EKAT
  tridiagonal bands. Both normalizations use field-state values and
  dimensionless `dt*rate` factors; a raw source rate is never mixed directly
  into a state-value scale.
- Negative values are classified against
  `10*epsilon*max(initial_field_scale,raw_field_scale)`. Roundoff and material
  counts are separate globally and for droplets and each aerosol reservoir.
  FE retains production clipping; IMEX returns a material-negativity failure
  without clipping.
- The physical top and bottom interface rates remain zero. The first-level
  synthetic rate is absent from tendencies and ledgers, so it is not counted
  twice. Activation source signs, limiter behavior, the current bottom source,
  and zero-cloud reset/resuspension are unchanged.

The current adapter was ported from archived arrays and mapping arguments to
the locked MAM4xx rank-3 aerosol views and current returned-`nsubmix` API. It
uses `top_lev=0,1,6`. The deliberately first-interface-only controller case
reports `canonical_nsubmix=2`, `exact_nsubmix=2`, and
`physical_only_nsubmix=1`, while every physical output remains exact.

## TDD and review record

The pre-change baseline passed its two registered FE/IMEX tests. New tests
were then compiled before implementation. The red build exited 1 with 18 FE
and 11 IMEX compile diagnostics for the intentionally absent controller,
signed-ledger, roundoff/material-classification, and finite-decomposition
interfaces. This established the required red boundary before those
interfaces were added.

The first executable run after implementation passed IMEX but exposed 102 FE
assertion failures. Review found that the new test compared a roundoff-size
signed ledger remainder to a unit scale instead of its column-ledger scale;
the production normalized ledger gate was already using the approved scale.
The same run showed the new residual checker masking the older, more specific
raw-update overflow error. The test scale was corrected and raw finiteness is
now checked before residual assembly. No acceptance tolerance was relaxed.

A subsequent arithmetic review replaced direct IMEX ledger subtraction with
explicit overflow bounds before subtraction. A second brief pass checked the
controller call graph, signs, units, no-flux rows, zero-cloud ordering,
clipping/operator decomposition, and current-adapter API against the approved
design. `git diff --cached --check` was clean, both warnings-as-errors builds
were clean. That pass also found that identical summary output alone was
insufficient determinism evidence, so a focused follow-up added three serial
advances and full-state byte comparisons to both FE and IMEX suites. The final
suites below passed. No production source or gitlink changed.

## Executed commands

The essential retained gate commands were:

```text
cmake -S docs/analysis/evidence/ndrop_integrator \
  -B /Users/burr114/e3sm-scratch/ndrop-screening/build/task2-current-debug \
  -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DKokkos_DIR=/Users/burr114/e3sm-scratch/ndrop-screening/kokkos-install-debug/lib/cmake/Kokkos \
  -DEKAT_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/ekat \
  -DMAM4XX_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/mam4xx \
  -DNDROP_ENABLE_SANITIZERS=ON -DNDROP_BUILD_EXACT_ADAPTER=ON
cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/task2-current-debug --parallel 4
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/task2-current-debug -V \
  --output-log /Users/burr114/e3sm-results/ndrop-screening/task2-fe-imex/ctest-debug-sanitizer.log

cmake -S docs/analysis/evidence/ndrop_integrator \
  -B /Users/burr114/e3sm-scratch/ndrop-screening/build/task2-current-release \
  -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_DIR=/Users/burr114/e3sm-scratch/ndrop-screening/kokkos-install-release/lib/cmake/Kokkos \
  -DEKAT_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/ekat \
  -DMAM4XX_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/mam4xx \
  -DNDROP_ENABLE_SANITIZERS=OFF -DNDROP_BUILD_EXACT_ADAPTER=ON
cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/task2-current-release --parallel 4
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/task2-current-release -V \
  --output-log /Users/burr114/e3sm-results/ndrop-screening/task2-fe-imex/ctest-release.log
```

The lock script was also run with `--verify-only` for each build type from a
clean detached worktree at `07c786d82d` and wrote the two provenance JSON files
listed below. The FE and IMEX tests perform three serial advances and compare
every output-state value byte for byte. Three external optimized invocations
also produced identical FE, IMEX, and adapter standard output.

## Numerical results

Debug/sanitizer and optimized builds produced the same scientific summaries:

```text
Forward Euler: RESULT passed=23195 failed=0
fe_diffusion order_range=1.00004:1.00004 min_r2=1
fe_activation_all_fields order_range=1.00514:1.00704 min_r2=0.999998

IMEX: RESULT passed=25138 failed=0
diffusion l1_order_range=0.999957:0.999957
diffusion linf_order_range=0.999957:0.999958 min_r2=1
activation l1_order_range=0.991946:0.992783
activation linf_order_range=0.993205:0.993723 min_r2=0.999998
coupled_all_fields order_range=0.99701:1.00319
coupled_all_fields max_finest_cross_error=9.0418e-05

Current MAM4xx adapter: RESULT failed_cases=0
cases=15 max_abs_difference=0 max_rel_difference=0 mismatches=0
```

Every result entering these decisions also satisfied residual `<=1e-11`,
ledger error `<=max(1e-11,5*nlev*epsilon)`, and the applicable material-
negativity rule. Scaling factors `1e-12`, `1`, and `1e12`, exact zero aerosol
ledgers, level counts `1,2,3,127,128`, fractional/zero cloud masks, limiter
active/inactive cases, and nonfinite diagnostic probes are asserted inside the
reported check totals.

## Retained records

| Record | SHA-256 |
|---|---|
| `/Users/burr114/e3sm-results/ndrop-screening/task2-fe-imex/ctest-debug-sanitizer.log` | `56f6aa2338cfd2c14cf5fdaca1b83ed08991565c8ee2603b81a769c92d5f73d5` |
| `/Users/burr114/e3sm-results/ndrop-screening/task2-fe-imex/ctest-release.log` | `48d58588ea14e0de2ce06f4b4a9d7052c89c72c504192b49cceba490d1ba69f9` |
| `/Users/burr114/e3sm-results/ndrop-screening/provenance/task2-fe-imex-debug.json` | `af9e5fce4889a4dd2d9cc5eed24c12c29085d45dc0732079ed431720fe8d20fd` |
| `/Users/burr114/e3sm-results/ndrop-screening/provenance/task2-fe-imex-release.json` | `06fd93fd4f1d144260f19bb497679028998feacf14a520396e13e425fe7e093b` |

## Limits and next gate

- Evidence is local Kokkos Serial CPU evidence, not GPU performance or
  portability evidence.
- Task 2 reuses the existing analytic and coupled suites; the complete common
  parameter grid and independent Dormand-Prince oracle belong to Task 4.
- IMEX still factors its two scalar matrix families per call; broader
  matrix-family reuse and performance effects are evaluated only after the
  common gate.
- Backward Euler is not implemented. Task 3 must begin with the approved
  block-stencil check and its own red tests.
- Whole-model EAMxx coupling, production distributions, and supported-machine
  testing remain outside this standalone task.

With these limits explicit, Task 2 is complete and the approved Task 3
standalone Backward Euler implementation may begin.
