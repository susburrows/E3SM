# EAMxx-MAM4xx `ndrop` Cross-Method Screening Execution Plan

**Approved for execution:** 2026-08-18
**Branch:** `codex/ndrop-cross-method-screening`
**Starting commit:** `e8cf1e2d8183b60a2179926ada8bb225e4c393b8`
**Authoritative design:** `docs/design/eamxx_mam4xx_ndrop.md`
**Scope:** standalone, test-only Forward Euler (FE), IMEX Euler, and Backward
Euler (BE) validation and local CPU screening. Production EAMxx and MAM4xx
physics are out of scope.

## Goal

Qualify each integrator independently for the approved `ndrop` activation and
vertical-mixing operator, compare only methods that pass their own numerical
gate, measure their optimized CPU cost on this Mac, estimate FE/implicit and
FE/IMEX-hybrid break-even conditions, and document the remaining production
evidence needed to revise the estimates.

## Global Constraints

- Do not initialize GSD, create a new project, or change production EAMxx or
  MAM4xx physics. All new methods and adapters remain standalone test code.
- Preserve unrelated working-tree changes. Work only on the isolated branch
  named above. Put dependency checkouts, builds, and other rebuildable output
  under `/Users/burr114/e3sm-scratch/ndrop-screening`; put retained results and
  provenance under `/Users/burr114/e3sm-results/ndrop-screening`.
- Before numerical-code changes, amend the authoritative design with the exact
  operator, boundaries, limiter/controller semantics, solver equations,
  diagnostics, tests, tolerances, and revision history recorded below.
- The physical operator has no-flux top and bottom boundaries for all schemes.
  At the first active level `t`, current repaired MAM4xx may use
  `D_top_controller(t) = zn(t) * eddy_diff(t) * csbot(t) * zs(t)` only in the
  FE stability reduction. It is not a physical flux and must never enter the
  IMEX or BE matrices.
- The current MAM4xx adapter must directly compare returned `nsubmix` for
  repaired `top_lev=0` and generalized `t=0,1,6`, including a deliberately
  top-controller-limiting profile. Legacy MAM4xx `1b02f272...` at `top_lev=0`
  remains expected-failure historical evidence, never a passing current gate.
- Do not change production gitlinks. A checked-in lock manifest and preparation
  script must construct detached clean dependency checkouts under
  `/Users/burr114/e3sm-scratch/ndrop-screening`, verify exact SHAs and clean
  state, and record compiler/build metadata under
  `/Users/burr114/e3sm-results/ndrop-screening`. Fail closed on mismatch.
- Pin the current comparison stack to E3SM
  `e2c578e40f6602a108eff577091de86a453624ab`, MAM4xx
  `0d076d8ebb004a18f05f51f97b08edd9def6826f`, EKAT
  `3a1485aecfed7a036e4fa30ce845498fd37975e2`, and Kokkos
  `23745d8494203711a6090d0a201c3ae3bf57853b`. Preserve the archived evidence
  stack and its E3SM, MAM4xx, EKAT, Kokkos, and Haero SHAs as historical.
- FE must reproduce current production subcycling and per-substep clipping.
  IMEX and BE must not silently clip. Report roundoff negatives separately from
  material failures. Optional fallback clipping is a separate experiment with
  its frequency, magnitude, and ledger effect isolated from operator error.
- For aerosol pair `j`, report total, clipping, and operator ledger changes:
  `delta_L_total = L_post - L_pre`, `delta_L_clip` accumulated per substep, and
  `delta_L_operator = delta_L_total - delta_L_clip`. Use an absolute tolerance
  for a zero ledger. Independently compute the FE per-substep physical residual
  `r_s = q_raw(s+1) - q_clipped(s) - dtmix*F(q_clipped(s))` and IMEX/BE residual
  `r = q_new - q_old - dt*(F_explicit + F_implicit)`, normalized only by
  like-unit field terms. Exercise scaling factors `1e-12`, `1`, and `1e12`.
- BE uses a host-only 2-by-2 block Thomas solve and no clipping. Its independent
  dense partial-pivoting oracle must assemble the physical equations directly,
  not reuse block bands or assembly helpers. Non-bottom coupling is
  `S_k = a_k I_(k+1)`. Number and mass modes are distinct when activation
  differs. Independently check the droplet source. Judge raw positivity before
  the approved zero-cloud projection, then test isolated, contiguous, and
  all-zero reset/resuspension exactly. If the bottom `max(0,tmpa)` branch would
  change, fail and report it as outside the approved linear solve. Solver
  breakdown is deterministic at `100*epsilon*||A||_infinity`.
- Factor once per distinct matrix family per column. IMEX has separate scalar
  cloud/droplet and interstitial families. BE has distinct number/mass block
  families per mode when activation differs, plus droplet scalar. Reuse only
  identical matrices/RHS families; never across columns. Repeated matched-error
  steps may reuse within the same column, family, timestep, coefficients, and
  branch.
- Required accuracy levels are `1,2,3,127,128`; overlap values
  `0,0.25,0.5,0.75,1`; activation strengths `0,0.1,0.5,0.9` of the limiter plus
  limiter-active `2.0`; active aerosol-pair counts `1,4,8,16,25`, corresponding
  to `3,9,17,33,51` transported fields including droplet and both reservoirs;
  zero, typical, and strong mixing; deterministic seed `0x4E44524F`.
- Use analytic exact solutions where available. Otherwise use an independently
  assembled adaptive Dormand-Prince 5(4) physical-RHS oracle at `rtol=1e-12`
  with field-scaled absolute tolerances, rerun ten times tighter, and require
  maximum fieldwise change no greater than `2e-11`.
- A first-order gate requires observed slope in `[0.85,1.15]`, `R^2 >= 0.98`,
  residual `<= 1e-11` in double precision, ledger error
  `<= max(1e-11,5*nlev*epsilon)`, no material negative below
  `-10*epsilon*field_scale`, and three serial runs with bitwise-identical
  outputs. Any required admissible IMEX or BE case with a material negative
  fails that method.
- Performance cases use levels `1,2,3,16,32,64,72,127,128`; pair counts
  `1,4,8,16,25`; batches `1,16,64,256,1024`; typical and strong mixing; ACI
  intervals `33.33,50,75,100,150,300` seconds; and FE `nsubmix`
  `1,2,3,5,8,16,21`.
- Time the full standalone advance: assembly, limiter/controller logic,
  matrix/RHS work, copies, factorizations, solves, positivity checks, FE
  clipping, state update, and fences. Also report separable subtimers and
  matched instrumentation-off/on overhead. Use 10 warmups and 30--100 repeated
  samples until the fixed-seed 10,000-resample bootstrap median 95% confidence
  interval has less than 5% half-width. Interleave methods and consume output
  checksums so optimization cannot remove the work.
- Runtime-versus-error targets are `1e-2,3e-3,1e-3,3e-4,1e-4`. Do not
  extrapolate. Define break-even as the smallest measured FE `nsubmix` whose
  median runtime is at least the candidate runtime; interpolate only between
  adjacent measured values and carry bootstrap uncertainty. Use only points
  without clipping or material negativity. Treat hybrid results as threshold
  and packing-overhead estimates, not production implementation. Make no GPU
  speedup claims from Mac CPU data.
- Each method must receive its own PASS/FAIL gate before entering performance
  comparisons. A failed method remains documented but is excluded. The final
  reviewer must use a clean checkout, rerun debug and optimized suites, and
  issue separate FE, IMEX, BE, and cross-scheme PASS/FAIL decisions.
- Every test must name the behavioral break it catches, derive expected values
  independently, and be observed failing before its implementation is added.
  Commit only logical reviewed chunks; do not push or publish.

## Task 1: Freeze the design, provenance, and reproducibility contract

Update `docs/design/eamxx_mam4xx_ndrop.md` before changing numerical code.
Document the controller-only first-active-level FE rate; signed matrix stencils
and boundary rows; BE solver/oracle and pivot rules; diagnostic ledgers and
residuals; matrix-family reuse; all validation/performance grids and acceptance
criteria; test-only scope; exact revisions; approval on 2026-08-18; and a
revision-history row. Add `dependency_lock.json` and
`prepare_locked_dependencies.sh` in this directory. The script must verify or
construct clean detached dependency checkouts under
`/Users/burr114/e3sm-scratch/ndrop-screening`, validate all locked commits and
gitlinks without changing production gitlinks, and emit a machine-readable
provenance record below `/Users/burr114/e3sm-results/ndrop-screening`
containing repository SHAs, dirty state, compiler identity/version, CMake
generator/build type/flags, CPU/OS identity, and the exact invocation. Add
behavior tests for clean success and closed failure on any SHA/dirty-state
mismatch. Update `README.md` with the locked workflow. Commit as
`docs(ndrop): approve current integrator screening contract`.

## Task 2: Requalify Forward Euler and IMEX

Use strict TDD to repair the diagnostic gap in `ndrop_integrator.hpp`,
`ndrop_integrator.cpp`, `ndrop_imex.cpp`, and their tests. Populate and validate
the total, clipping, and operator ledger changes and the independent physical
residuals; distinguish roundoff from material negativity; apply the scaling
mutation cases; and make nonfinite diagnostics fail closed. Extend the current
MAM4xx adapter for repaired `top_lev=0` and `t=0,1,6`, including a top-limiting
controller profile and direct `nsubmix` equality. Keep the archived legacy
expected failure clearly separated. Run focused tests, then the full debug and
optimized FE/IMEX gates using the locked current stack. Commit code/tests as
`test(ndrop): requalify forward Euler and IMEX gates`, then record executed
commands, counts, hashes, and method-specific PASS/FAIL results in reviewed
human-readable evidence and commit as
`docs(ndrop): record forward Euler and IMEX gate review`.

## Task 3: Implement and gate standalone Backward Euler

Before code, confirm the design contains the complete signed 2-by-2 block
stencil, scalar droplet system, boundary rows, source signs, projection order,
and failure behavior. With strict TDD, add `ndrop_backward_euler.cpp` and
`test_backward_euler.cpp`; extend only the test harness interfaces needed.
Implement host-only block Thomas factorization/solve, distinct matrix-family
reuse, deterministic pivot failure, raw positivity reporting, and approved
post-solve reset/resuspension. Add a dense partial-pivoting oracle assembled
directly from the physical equations. Cover identity/zero, uniform state,
analytic no-flux modes, activation transfer, typical/strong mixing, all overlap
values, small/L128-like columns, strengths/field counts, boundaries, ledger
conservation, positivity, nonfinite inputs, determinism, reset patterns,
scaling, dense-oracle agreement, and refinement convergence. Run debug,
sanitizer where supported, and optimized gates. Commit as
`test(ndrop): add backward Euler validation gate`, then add reviewed evidence
and commit as `docs(ndrop): record backward Euler gate review`.

## Task 4: Add the independent common accuracy and convergence gate

Add `ndrop_reference.cpp` and `test_cross_scheme.cpp`. Implement analytic
references and the independently assembled adaptive Dormand-Prince 5(4)
physical-RHS oracle with the two-resolution self-check. Run the complete common
matrix and produce machine-readable results plus human-readable tables for
accuracy, observed order and `R^2`, residuals, positivity classification,
clipping frequency/magnitude, total/operator/clipping ledger changes,
boundaries, determinism, and reference self-consistency. Enforce method-local
PASS/FAIL gates and exclusion. Commit as
`test(ndrop): add cross-scheme accuracy gate`.

## Task 5: Add and run the optimized Mac performance screen

Add `benchmark_ndrop_integrators.cpp` and required CMake targets. Verify the
benchmark implementation against the validated path. Time the full advance and
separable setup/factorization/solve/update/fence phases for the required matrix,
with warmups, interleaved order, adaptive repetitions, robust statistics,
checksums, and instrumentation-off/on matched runs. Record compiler and machine
metadata. Emit raw machine-readable samples; do not commit large generated
build trees. Include only schemes that passed Task 4. Commit as
`perf(ndrop): add validated local screening benchmark`.

## Task 6: Analyze error-cost and break-even behavior

Add `analyze_screening.py` with tests using synthetic fixtures whose medians,
confidence intervals, exclusions, adjacent interpolation, and break-even
answers are hand-derived. Analyze the executed accuracy and timing outputs for
the required tolerances, FE `nsubmix`, levels, transported fields, batch sizes,
mixing strengths, and ACI intervals. Quantify uncertainty. Evaluate an FE/IMEX
threshold hybrid, including packing and workload-variation sensitivity. Keep
executed CPU measurements separate from GPU-structure estimates.

## Task 7: Write the screening report and recommendation

Create `docs/analysis/eamxx_ndrop_integrator_screening.md`. Include separate
FE/IMEX/BE PASS/FAIL gates; accuracy/convergence tables; positivity, clipping,
and conservation tables; typical/strong performance tables; runtime-versus-
error results; break-even `nsubmix` by level/field count; batching and qualitative
GPU implications; exact provenance; uncertainty, limitations, failed/excluded
cases; and the production distributions/counters needed to revise estimates.
Choose exactly one recommendation: current FE, IMEX, BE, FE/IMEX hybrid, or
insufficient evidence. Commit analysis and report as
`docs(ndrop): report integrator screening results`.

## Task 8: Independent final clean-checkout validation and gate review

Give a strict reviewer the design, this plan, locked provenance, complete branch
diff, raw evidence, and exact commands. The reviewer must create or use a clean
detached checkout, reconstruct locked dependencies, rerun the complete debug
and optimized suites, verify hashes/results, and issue separate FE, IMEX, BE,
and cross-scheme PASS/FAIL decisions. It must also judge whether the report
separates measurements from estimates and whether exclusions and uncertainty
are correct. Address all critical/important findings through the review loop.
Do not accept the final comparison without this review.

## Required Outputs

1. Per-method validation report with independent PASS/FAIL gates.
2. Accuracy and convergence tables.
3. Positivity, clipping, and conservation tables.
4. Typical and strong-mixing performance tables.
5. Runtime-versus-error comparisons.
6. Break-even `nsubmix` estimates by level and transported-variable count.
7. One explicit recommendation from the approved list.
8. A concise production-measurement list for revising offline estimates.
