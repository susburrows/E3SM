# Backward Euler Task 3 gate review

**Execution date:** 2026-08-24

**Code/test commits:** `2f43c8c1e4ebd6b97152414ab02dca760bc883e6`,
`c6f7f18a5d58e118cf4146a790fe69895a1c4adb`

**Final code/test tree:** `f0448c4f7affe62abb9716160d0c246e6ed7d00f`

**Scope:** standalone host-only Backward Euler implementation and method-local
gate; no production EAMxx or MAM4xx source, interface, equation, or gitlink
change

## Decision summary

| Gate | Decision | Executed evidence |
|---|---|---|
| Backward Euler Task 3 method-local gate | **PASS** | 23,165 checks in both debug/ASan/UBSan and optimized builds; analytic diffusion and activation slopes and `R^2` pass; block Thomas agrees with an independently assembled dense partial-pivoting oracle; matrix and physical residuals, signed ledgers, positivity, scaling, projection, failure behavior, family reuse, and three-run bitwise-state checks pass. |
| Complete four-target standalone regression | **PASS** | FE 23,195 checks, IMEX 25,138 checks, BE 23,165 checks, and all 15 locked current-MAM4xx adapter cases pass in both builds. |
| Production or GPU readiness | **NOT CLAIMED** | This is local Kokkos Serial CPU evidence for a test-only host solver. Task 4 common accuracy, Task 5 timing, supported-machine integration, and any production/GPU design remain separate gates. |

The method-local decision permits the approved common cross-method validation
to begin. It does not make BE a production recommendation or establish that
one full-interval BE step is accurate enough for any EAMxx configuration.

## Overall assessment

The approved standalone design is ready for the Task 4 cross-method phase.
The governing equations, signed block rows, source ordering, boundary
conditions, pivot rule, positivity/projection order, diagnostics, tolerances,
and reproducibility contract are implementable and agree with the as-built
gate. No critical or major design finding remains for Task 3.

## Design-review findings

### Minor: stale current-versus-archived status statements

- **Location:** Scope and Status, canonical test operator, staged gates, and
  provisional method-status table.
- **Problem:** Several statements still said BE had not begun or restricted
  exact comparison to `top_lev=1`, despite the locked repaired-current adapter
  already passing `top_lev=0,1,6` and the new BE gate passing.
- **Why it matters:** The stale text could conflate the archived
  `eddy_diff(-1)` defect with current repaired MAM4xx or obscure the actual
  gate sequence.
- **Resolution:** Corrected in the synchronized design. The archived expected
  failure remains separate; the current method-local BE gate is recorded as
  passed while Task 4, production, and GPU work remain pending.

### Minor: admissible columns cannot exercise a nonzero bottom activation row

- **Location:** signed no-flux operator and implicit rows.
- **Problem:** The canonical activation limiter caps every rate by
  `rate_below`; the explicit no-flux bottom makes that rate zero. Therefore an
  ordinary valid harness call cannot numerically exercise the retained
  current-C++ nonzero bottom source branch.
- **Why it matters:** Dense-oracle agreement on admissible columns alone would
  not detect a sign error in
  `[[lambda-ha,-ha],[ha,lambda+ha]]` or the bottom droplet contribution.
- **Resolution:** A narrow test-only inspection now extracts the actual
  assembled bottom block before canonical limiting and checks all four signed
  entries independently. A separate hand calculation verifies the positive
  bottom droplet source `sum_m a_m*(2*I_m+C_m)`. Branch departure and material
  negativity have deterministic failure probes. This does not change the
  physical advance.

## Equation-to-algorithm map

| Approved equation or rule | As-built location and evidence |
|---|---|
| `x_k=(C_k,I_k)` signed `L`, `B`, and `U` blocks | `assemble_block_system`; independent dense oracle directly assembles the `2*n` scalar matrix from the physical equations. |
| `S_k=a_k I_(k+1)` with opposite cloud/interstitial signs | Upper-block entries `-h*a` and `-h*d_down+h*a`; analytic activation convergence and dense agreement exercise all 25 aerosol fields with independent mode/type maps. |
| Retained bottom `S=a(I+C)` row | Bottom diagonal block inspection verifies `lambda-ha,-ha,ha,lambda+ha`; outside-branch probe fails deterministically. |
| Explicit no-flux top and bottom | Shared interface endpoints are exactly zero; identity levels `1,2,3,127,128`, analytic no-flux modes, and zero-overlap mixing pass. The FE-only first-level controller has a bitwise no-effect BE test. |
| Factor each distinct activation family once | Exact-coefficient grouping; a distinct-rate case reports eight factorizations for 25 aerosol right-hand sides, while an identical-rate case reports one. |
| Host block Thomas with unsafe-pivot threshold | Reduced diagonal blocks are inverted with deterministic partial row selection; scalar pivots at or below `100*epsilon*||A||_infinity` fail. Near-threshold block and scalar probes pass. |
| Droplet solve after modal-number aerosol blocks | One scalar cloud matrix is factored once; its right-hand side uses the solved modal-number state. The dense oracle assembles its matrix and source separately. |
| Raw positivity before zero-cloud projection | Per-field minima and roundoff/material counts are populated before isolated, contiguous, all-zero, and zero-overlap reset/resuspension cases. BE never clips. |
| Signed total/clipping/operator ledgers | All 25 pairs report `clipping=0`, `operator=total`, and the approved relative or zero-ledger absolute normalization. |
| Independent like-unit residual | A physical-equation residual is assembled independently of block bands; aerosol mass, aerosol number, and droplet number scales are never summed. |

## Locked inputs and environment

The dependency verifier confirmed clean detached checkouts before both final
provenance records were written.

| Input | Revision |
|---|---|
| E3SM baseline | `e2c578e40f6602a108eff577091de86a453624ab` |
| MAM4xx | `0d076d8ebb004a18f05f51f97b08edd9def6826f` |
| EKAT | `3a1485aecfed7a036e4fa30ce845498fd37975e2` |
| Kokkos | `23745d8494203711a6090d0a201c3ae3bf57853b` |

Both builds used Apple clang 21.0.0, CMake 4.4.2, Ninja, double precision,
and Kokkos 4.5.1 Serial on an arm64 Apple M2 Max running macOS 26.6.1. The
debug harness enabled ASan and UBSan with a Kokkos Debug build. Kokkos's
separate debug-bounds option was not enabled. The optimized harness and Kokkos
dependency used Release `-O3 -DNDEBUG`. Every harness target compiled with
`-Wall -Wextra -Wpedantic -Werror`.

## TDD and review record

The pre-change current-stack baseline passed its three registered FE, IMEX,
and adapter targets. The complete BE test target and public test-only
interfaces were then compiled before implementation. Compilation succeeded,
but linking exited 1 with only the intentionally undefined
`advance_backward_euler` and `backward_euler_probe_is_rejected` symbols. This
established the red boundary.

The first implementation run passed 23,023 checks. A separate arithmetic and
failure-path review then added unsafe-pivot rather than zero-pivot probes,
finite raw minima, explicit bottom-branch arithmetic checks, roundoff/material
classification probes, zero-overlap mixing, solved-source overflow, and
right-hand-side overflow. The deliberately constructed RHS-overflow case first
failed at the earlier block-residual stage because its test magnitudes also
overflowed that normalization. The input was analytically reduced until every
block residual term remained finite while the sum of four finite modal droplet
sources overflowed only the intended RHS. No tolerance or acceptance rule was
changed.

The formal design review found the nonzero-bottom coverage gap described
above. The follow-up commit added independent signed-entry and bottom-source
checks. A second brief pass rechecked the controller call graph, `L/B/U` signs,
mode/type activation families, no-flux rows, block and scalar pivot scaling,
dense-oracle independence, physical residual units, raw positivity before
projection, and fail-closed ledger decomposition. `git diff --cached --check`
and both warnings-as-errors builds were clean.

## Executed commands

The final clean detached source worktree was
`/Users/burr114/e3sm-scratch/ndrop-screening/worktrees/task3-final` at
`c6f7f18a5d58e118cf4146a790fe69895a1c4adb`. The essential retained commands
were:

```text
cmake -S /Users/burr114/e3sm-scratch/ndrop-screening/worktrees/task3-final/docs/analysis/evidence/ndrop_integrator \
  -B /Users/burr114/e3sm-scratch/ndrop-screening/build/task3-current-debug \
  -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DKokkos_DIR=/Users/burr114/e3sm-scratch/ndrop-screening/kokkos-install-debug/lib/cmake/Kokkos \
  -DEKAT_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/ekat \
  -DMAM4XX_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/mam4xx \
  -DNDROP_ENABLE_SANITIZERS=ON -DNDROP_BUILD_EXACT_ADAPTER=ON
cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/task3-current-debug --parallel 4
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/task3-current-debug -V \
  --output-log /Users/burr114/e3sm-results/ndrop-screening/task3-backward-euler/ctest-debug-sanitizer.log

cmake -S /Users/burr114/e3sm-scratch/ndrop-screening/worktrees/task3-final/docs/analysis/evidence/ndrop_integrator \
  -B /Users/burr114/e3sm-scratch/ndrop-screening/build/task3-current-release \
  -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_DIR=/Users/burr114/e3sm-scratch/ndrop-screening/kokkos-install-release/lib/cmake/Kokkos \
  -DEKAT_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/ekat \
  -DMAM4XX_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/mam4xx \
  -DNDROP_ENABLE_SANITIZERS=OFF -DNDROP_BUILD_EXACT_ADAPTER=ON
cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/task3-current-release --parallel 4
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/task3-current-release -V \
  --output-log /Users/burr114/e3sm-results/ndrop-screening/task3-backward-euler/ctest-release.log
```

The lock script was run with `--verify-only` for both build types from that
clean detached worktree. The BE suite itself performs three serial advances
and compares the complete output state byte for byte. Three additional
optimized process invocations produced identical standard output with SHA-256
`87b385f9089e8fec52b1b6cd367f5ff81680798107cf2f6e7facd29a774dc098`.

## Numerical results

Debug/sanitizer and optimized builds produced the same BE summary:

```text
CONVERGENCE case=be_diffusion order_range=0.999955:0.999957 min_r2=1
CONVERGENCE case=be_activation order_range=0.993111:0.994964 min_r2=0.999998
RESULT passed=23165 failed=0
```

The suite covers identity at `nlev=1,2,3,127,128`; nontrivial and uniform
L128-like columns; overlap `0,0.25,0.5,0.75,1`; typical and strong mixing;
activation strengths `0,0.1,0.5,0.9,2` times the limiter; active aerosol-pair
counts `1,4,8,16,25`; all 25 mapped aerosol pairs; droplets; limiter-active
rates; scaling factors `1e-12,1,1e12`; and exact zero aerosol ledgers. The
method-local category values are represented without claiming the complete
Task 4 Cartesian product.

Every successful decision case satisfied matrix and independently assembled
physical residuals at most `1e-11`, ledger error at most
`max(1e-11,5*nlev*epsilon)`, and no material negative. Deterministic failures
cover nonfinite input and initial ledger, matrix assembly, solved-state source,
droplet RHS, solution checker, matrix and physical residuals, ledger
normalization and decomposition, unsafe block and scalar pivots, material
negativity, and the excluded bottom-source branch change.

## Retained records

| Record | SHA-256 |
|---|---|
| `/Users/burr114/e3sm-results/ndrop-screening/task3-backward-euler/ctest-debug-sanitizer.log` | `63e8f629bedacc9d00b4271cd4d382c3b96c361617780172893e8c0117be40f6` |
| `/Users/burr114/e3sm-results/ndrop-screening/task3-backward-euler/ctest-release.log` | `b2b546c13c59fd9885d9b551dc9f9d1d1ac0b96dfad0ed9837ddea960073524d` |
| `/Users/burr114/e3sm-results/ndrop-screening/provenance/task3-backward-euler-debug.json` | `f318fba580f0f41fd45d0b4cdf2b0b140784c1e2c111d91eed96e61d4f54d378` |
| `/Users/burr114/e3sm-results/ndrop-screening/provenance/task3-backward-euler-release.json` | `e82b0c2270a5496c6a269dad6368277e4f1e2c8bcb4a0b5cc7e251bd021a9ce5` |

## Open questions and residual risks

- The complete common Cartesian matrix and independent adaptive
  Dormand-Prince physical-RHS oracle belong to Task 4. Task 3 cannot establish
  matched-error ranking or full-interval suitability.
- The nonzero bottom source is a retained current-C++ algebraic branch, but the
  canonical no-flux limiter makes it zero for admissible harness columns. Its
  signed block and droplet source are unit-checked; a production-physics owner
  must resolve the C++ versus Fortran bottom formulation before adoption.
- Evidence is local host Kokkos Serial CPU evidence. It is not a device block
  solver, supported-machine EAMxx integration, GPU portability, or performance
  result.
- No fallback clipping exists. A Task 4 admissible case with a material raw
  negative must fail and exclude BE for that case rather than silently repair
  it.
- Whole-model coefficient evolution, outer ACI splitting, and coupling to
  other physics remain outside this frozen-operator gate.

## Recommended next revision

Proceed to the approved Task 4 common cross-method accuracy and convergence
gate from one revision containing all three passing methods. Preserve the
method-local PASS decisions, but allow the common matrix to fail or exclude a
method independently if reference consistency, positivity, ledger, residual,
or common-solution criteria are not met.
