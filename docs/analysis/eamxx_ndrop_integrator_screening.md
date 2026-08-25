# EAMxx-MAM4xx `ndrop` Integrator Screening

**Draft status:** Task 7 structure prepared on 2026-08-24 before completion of
the retained Task 5 performance run and Task 6 analysis.

**Decision status:** **PENDING.** This draft does not select Forward Euler
(FE), IMEX Euler, Backward Euler (BE), an FE/IMEX hybrid, or insufficient
evidence.

**Evidence boundary:** Tasks 2--4 provide executed numerical evidence. The
Task 5 benchmark implementation is committed, but no partial Task 5 shard,
progress log, or timing value is used here. Every performance value,
runtime-versus-error result, break-even result, hybrid result, and final
recommendation is explicitly pending.

A static review of the committed Task 5 source found no Critical issue, but
identified three Important limitations that remain open: the nominal
transported-field axis does not change executed width, IMEX exposes a combined
factorization-plus-solve phase rather than separate phase timings, and raw
phase/fence samples are not serialized. These findings do not change what the
complete uninstrumented timing region contains, but no pre-remediation artifact
is eligible for Task 6 or the final recommendation. They also block field-count
scaling and phase-level conclusions.

The pending Task 5R design amendment selects exact remedies for those
findings: occupancy-only labels with fixed executed-width metadata, a supported
combined IMEX factorization-plus-solve phase, and a repetition-level phase
sidecar. It also adds paired-output completion records and fail-closed analysis
boundaries. These are approved-plan requirements, not implemented or executed
results. Source changes and a remediated retained run require explicit approval
of the design amendment.

## Executive status

The standalone numerical gates admit FE, IMEX, and BE to the local CPU
performance screen. Each method passed all 1,875 cases in the common Task 4
matrix, converged at first order to the same independently integrated frozen
operator, met residual and aerosol-ledger tolerances, and produced neither a
material raw negative nor an FE clipping event in that matrix.

This result establishes numerical qualification for the screen, not
production suitability. The largest measured coarse one-call normalized
L-infinity error in Task 4 was `0.1129882`. Consequently, equal-interval
timing alone cannot select a method: the final comparison must use measured
runtime at supported error levels, retain positivity and clipping exclusions,
and carry timing uncertainty. Those performance-dependent steps are pending.

| Evidence stage | Status | What it establishes |
|---|---|---|
| FE and IMEX method-local gate (Task 2) | **PASS, executed** | Correctness of the standalone FE and IMEX paths under their method-local contracts |
| BE method-local gate (Task 3) | **PASS, executed** | Correctness of the host-only block-Thomas BE path under its method-local contract |
| Common accuracy gate (Task 4) | **PASS, executed** | Independent-reference agreement, first-order behavior, common solution, conservation, positivity, and deterministic results on the required matrix |
| Optimized Mac screen (Task 5) | **PENDING retained execution evidence** | Local Kokkos Serial CPU timing and phase/accounting results |
| Error-cost and break-even analysis (Task 6) | **PENDING** | Supported error points, break-even thresholds, batching sensitivity, and hybrid estimates |
| Final recommendation (Task 7) | **PENDING** | Exactly one choice from the approved recommendation set |

## Scope and decision question

This report concerns the standalone, test-only integration of the frozen
`ndrop` vertical-mixing and activation operator described in the
[approved design](../design/eamxx_mam4xx_ndrop.md) and
[execution plan](evidence/ndrop_integrator/2026-08-18-cross-method-screening-plan.md).
It compares:

- current production-equivalent, subcycled FE with per-substep clipping;
- IMEX Euler with explicit activation and implicit vertical mixing and
  cloud-to-interstitial return; and
- fully implicit BE with a host-only 2-by-2 block-tridiagonal aerosol solve.

The decision question is whether the completed standalone evidence supports
retaining FE, advancing IMEX or BE toward a separately approved production
prototype, advancing a separately designed FE/IMEX hybrid, or concluding that
the evidence is insufficient.

The following are outside this report's authorization and evidence:

- changes to production EAMxx or MAM4xx equations, interfaces, gitlinks, or
  runtime behavior;
- supported-machine EAMxx integration, coupled climate response, or
  whole-model speedup;
- GPU implementation or GPU performance;
- recomputation of activation, turbulence, cloud overlap, or geometry inside
  an implicit solve;
- higher-order implicit methods, nonlinear solves, and optional fallback
  clipping; and
- resolution of the current C++ versus Fortran bottom-source difference.

## Evidence labels used in this report

The report uses four labels to prevent estimates from being presented as
measurements:

- **Measured:** produced by an executed retained gate or completed benchmark,
  with an identified source revision, command, artifact, and hash.
- **Derived:** calculated only from measured retained inputs by the committed
  analysis procedure. A derived value must identify its inputs and uncertainty
  rule.
- **Estimated:** a model-based quantity, including Richardson error estimates,
  workspace formulas, operation counts, break-even interpolation, or hybrid
  packing sensitivity. It is not direct timing or production evidence.
- **Pending:** not available from complete retained evidence at the time of
  this draft. A placeholder must not be read as zero, not applicable, or a
  negative result.

Unless a table says otherwise, all accuracy values below are **measured Task 4
standalone values**. All performance tables remain **pending**.

## Frozen operator and candidate methods

### Common state and operator

At full field count, the state contains droplet number plus 25 interstitial
and 25 cloud-borne aerosol fields, or 51 transported fields. The operator
freezes pressure weights, one shared nonnegative exchange coefficient per
physical interface, cloud-overlap factors, mode/type activation rates, and
geometry over one ACI call. Shared-interface construction enforces weighted
flux cancellation, with explicit zero flux at the top and bottom.

For aerosol pair `j`, the frozen continuous equations are

\[
\frac{dC_j}{dt}=L_C C_j+S_j(I_j,C_j),
\]

\[
\frac{dI_j}{dt}=L_I I_j+R C_j-S_j(I_j,C_j),
\]

and droplet number follows

\[
\frac{dN_c}{dt}=L_C N_c+S_N(I,C).
\]

The integrator is the only intended algorithmic difference. Every method
uses the same operator data and the same post-integration zero-cloud state
transition: cloud-borne aerosol is returned to the interstitial reservoir and
cloud-borne aerosol and droplet number are reset at exactly zero cloud
fraction. Raw IMEX and BE positivity is judged before that transition.

The repaired current-MAM4xx first-active-level controller rate may increase
the FE substep count. It is not a physical flux and is excluded from every
tendency, ledger, IMEX matrix, BE matrix, and implicit residual. The physical
top and bottom interface rates remain zero for all methods.

### Forward Euler

FE reproduces current MAM4xx adaptive subcycling, activation-rate limiting,
alternating buffers, and componentwise `max(0,qnew)` clipping after every
substep. It retains the raw state, clipping counts and deficits, and the signed
separation of total, clipping, and operator aerosol-ledger change. The fixed
FE `nsubmix` values in the performance screen are controlled comparison
points, not changes to the outer MAC/ACI process-group subcycle count.

### IMEX Euler

IMEX evaluates activation at the old state while treating vertical mixing and
cloud-to-interstitial return implicitly. It uses one scalar cloud/droplet
matrix family and one scalar interstitial matrix family, factors each distinct
family once per column, and applies the factors to multiple right-hand sides.
The standalone CPU path uses the EKAT serial Thomas solver. IMEX does not clip;
a material raw negative is a case failure.

### Backward Euler

BE evaluates mixing, return, and activation at the new state. Each aerosol
pair forms a 2-by-2 block-tridiagonal system. The standalone host path uses a
block-Thomas solve, reusing exact matrix families only within a column; the
droplet scalar solve follows the aerosol-number blocks. Unsafe pivots,
material raw negatives, or departure from the retained bottom-source linear
branch fail deterministically. No production or device block solver exists in
this study.

### Conservation and positivity contract

For every aerosol pair, all methods report

\[
\Delta L_{operator}=\Delta L_{total}-\Delta L_{clip}.
\]

IMEX and BE require `Delta L_clip=0`; FE records the per-substep clipping
contribution separately. Residual normalization combines only like-unit terms.
A negative is material below `-10*epsilon*field_scale`; roundoff and material
counts remain separate. Nonfinite state, coefficient, matrix, right-hand side,
solution, ledger, residual, or normalization fails closed.

## Executed numerical qualification

### Gate chronology

| Method or gate | Executed decision | Revision and retained review | Principal result |
|---|---|---|---|
| FE | **PASS** | Code/test `e67097b6bc`, determinism follow-up `07c786d82d`; [Task 2 review](evidence/ndrop_integrator/fe_imex_gate_review_20260824.md) | 23,195 checks in each debug/sanitizer and optimized run |
| IMEX | **PASS** | Code/test `e67097b6bc`, determinism follow-up `07c786d82d`; [Task 2 review](evidence/ndrop_integrator/fe_imex_gate_review_20260824.md) | 25,138 checks in each debug/sanitizer and optimized run |
| BE | **PASS** | Code/test `2f43c8c1e4`, review follow-up `c6f7f18a5d`; [Task 3 review](evidence/ndrop_integrator/backward_euler_gate_review_20260824.md) | 23,165 checks in each debug/sanitizer and optimized run |
| Common FE/IMEX/BE matrix | **PASS** | Code/test `7807aedfde`; [Task 4 review](evidence/ndrop_integrator/cross_scheme_gate_review_20260824.md) | 1,875 passing cases per method; 11,880 checks in each cross-scheme executable |

The exact-current adapter also passed 15 repaired-current MAM4xx cases with
zero state mismatch and exact returned-`nsubmix` equality. This adapter result
does not validate the archived defective MAM4xx revision.

### Accuracy methodology

The common gate is the Cartesian product of:

- active levels `1,2,3,127,128`;
- cloud overlap `0,0.25,0.5,0.75,1`;
- activation strengths `0,0.1,0.5,0.9,2.0` times the limiter;
- active aerosol-pair counts `1,4,8,16,25`, producing nominal nonzero-state
  labels of `3,9,17,33,51` fields; and
- zero, typical, and strong mixing.

This gives 1,875 columns and 5,625 method records. The deterministic seed is
`0x4E44524F`.

The active-pair dimension varies which aerosol reservoirs are initialized
nonzero, not the width of the executed integrator. FE, IMEX, and BE still
traverse all 25 aerosol pairs and all 51 fields. It is therefore a
nonzero-state-occupancy axis and cannot support cost scaling by transported
field count.

Analytic no-flux diffusion and activation cases check all 51 fields. General
coupled cases use an independently assembled adaptive Dormand--Prince 5(4)
physical-RHS reference at `rtol=1e-12`, repeated at `1e-13`; the allowed
maximum normalized reference change is `2e-11`. Candidate methods do not
supply the reference RHS, matrix, or time integration.

Normalized weighted-L1 and L-infinity errors are calculated separately for
every field. First-order acceptance requires slopes in `[0.85,1.15]` and
`R^2 >= 0.98`. Double-precision residuals must be at most `1e-11`, and aerosol
ledger error must be at most `max(1e-11,5*nlev*epsilon)`. Every required
successful case must have no material raw negative. Each method also performs
three full serial advances with bitwise-identical state and diagnostics.

### Accuracy and convergence results

**Measured Task 4 results:**

| Method | Cases | Failures | Max weighted L1 | Max Linf | Max residual | Max ledger | Material negatives | FE clipping cases | Slope range | Min `R^2` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| FE | 1,875 | 0 | `0.091525245` | `0.1129882` | `1.4389395e-16` | `2.3628458e-15` | 0 | 0 | `1.0011674:1.0016497` | `0.9999999` |
| IMEX | 1,875 | 0 | `0.06688989` | `0.070265567` | `3.2474441e-16` | `2.3628644e-15` | 0 | N/A | `0.99825832:0.99868307` | `0.99999989` |
| BE | 1,875 | 0 | `0.057285903` | `0.072283807` | `4.0419594e-16` | `2.4656757e-15` | 0 | N/A | `0.99835794:0.99883807` | `0.9999999` |

The measured reference change between the two oracle tolerances was
`2.0331196e-13`. All three pairwise method differences decreased at the finest
refinement and satisfied the approved Richardson common-solution bound in both
norms for every field. Debug/sanitizer and optimized CSV and Markdown outputs
were bitwise identical.

These maxima summarize the complete matrix and must not be interpreted as
same-case rankings among methods. In particular, the coarse errors are not
production acceptance thresholds. The common PASS means that all methods are
eligible for performance screening; it does not show that one step over a
production ACI interval is sufficiently accurate.

### Positivity, clipping, and conservation results

| Property | FE | IMEX | BE | Evidence status |
|---|---:|---:|---:|---|
| Material-negative records in Task 4 | 0 | 0 | 0 | **Measured** over 1,875 cases per method |
| Clipping-active Task 4 cases | 0 | Prohibited | Prohibited | **Measured**; FE clipping could not repair a common-gate point |
| Maximum normalized ledger error | `2.3628458e-15` | `2.3628644e-15` | `2.4656757e-15` | **Measured** |
| Maximum normalized residual | `1.4389395e-16` | `3.2474441e-16` | `4.0419594e-16` | **Measured** |
| Production/coupled clipping frequency | **Pending production sampling** | N/A | N/A | Not established by the synthetic gate |

## Optimized Mac performance methodology

The pre-remediation Task 5 implementation is in
`benchmark_ndrop_integrators.cpp`, `ndrop_benchmark.cpp`,
`ndrop_benchmark.hpp`, and the phase-timer support. Its benchmark-support tests
encode checks for the required matrix, statistics, state neutrality, case
ordering, accounting, and timing-path equivalence. The complete retained run,
artifact hashes, and performance review remain pending.

### Pending Task 5R evidence contract

**Status: PENDING EXPLICIT APPROVAL AND IMPLEMENTATION.** No value from the
pre-remediation run bound to
`dba6d00598d86a239a2bf2743d95325821c9378c` may be merged into or substituted
for the remediated evidence.

The remediated summary will declare `artifact_schema=ndrop_summary_v2` and a
validated `run_id`. The varying fields will be `nonzero_aerosol_pairs` and
`nonzero_state_fields`, with
`nonzero_state_fields == 1 + 2*nonzero_aerosol_pairs`. Every row will also
declare the fixed executed width, `executed_aerosol_pairs=25` and
`executed_transport_fields=51`. The report may analyze sensitivity to nonzero
occupancy but may not infer transported-width scaling.

The linear phase contract will be method-specific:

| Method | `linear_phase_mode` | Factorization | Solve | Combined factorization-plus-solve |
|---|---|---|---|---|
| FE | `none` | unavailable | unavailable | unavailable |
| IMEX | `combined` | unavailable | unavailable | measured diagnostically |
| BE | `separate` | measured diagnostically | measured diagnostically | unavailable |

Unavailable values will be blank, never numeric zero. The IMEX combined phase
will include matrix-view allocation/coefficient copy and both supported EKAT
serial Thomas factorization-and-solve calls. No unsupported EKAT helper or new
solver is authorized, and BE's phase boundaries remain unchanged.

Every instrumented repetition for every `validation_ok=true` summary row will
be serialized in an `ndrop_phase_samples_v1` sidecar, including rows later
excluded from comparative analysis. Invalid rows will have zero repetitions
and no sidecar rows. Sidecar instrumented totals must match the summary raw
instrumented list element by element, and all phase/fence medians must be
reconstructible with the benchmark's exact median algorithm.

The retained screen will use an immutable pre-run identity manifest and 135
aligned 180-case shards. Each summary/sidecar pair is accepted only when both
files validate, both hashes match, and its completion record was written last.
The deterministic merge and final immutable manifest must bind both merged
artifacts. Only complete primary uninstrumented timings may enter speedup,
matched-error, break-even, hybrid, or recommendation calculations; phase data
remain diagnostic only.

### Required benchmark matrix

The implementation emits 24,300 method rows:

- 18,900 FE rows and 2,700 rows each for IMEX and BE;
- active levels `1,2,3,16,32,64,72,127,128`;
- nonzero aerosol-pair inputs `1,4,8,16,25`, which vary nonzero-state occupancy
  while all methods still execute 25 pairs and 51 fields;
- batch sizes `1,16,64,256,1024` columns;
- synthetic typical and strong mixing scales `0.05` and `0.5`;
- ACI intervals `33.33,50,75,100,150,300` seconds; and
- FE substep counts `1,2,3,5,8,16,21`.

Each equal-input group contains seven FE rows followed by IMEX and BE. The
measurement order is cyclically rotated across warmups and repetitions so a
method does not always occupy the same timing position.

### Inputs and timed region

Each group uses identical deterministic initial states, pressure weights,
cloud fractions, interface exchange, and activation coefficients. The input
state and frozen `OperatorData` are constructed outside the timed batch. Each
timed method call includes its method-local setup, limiter/controller logic,
matrix and right-hand-side work, required copies, factorizations, solves,
positivity handling, state update, and final Kokkos fence. Output checksums are
consumed after timing so the compiler cannot discard the result.

The primary timing path disables detailed diagnostics but preserves the
validated state and required status/positivity behavior. Before timing, the
code executes both the full diagnostic path and the lean path and requires
bitwise-identical state, identical clipping count, and identical material-
negative count. A separate instrumented pass records setup, method-specific
linear work, update, diagnostics, and fence phase observations and reports its
overhead relative to the primary pass. Under Task 5R, IMEX coefficient copying,
factorization, and solve will be recorded only in the combined
factorization-plus-solve phase; its separate factorization and solve values are
unavailable rather than zero. BE retains separate factorization and solve
phases, while FE's linear-phase values are unavailable.

The batch loop is a local Kokkos Serial CPU loop over homogeneous columns. It
does not reuse a matrix or factorization across columns. Batch-size effects are
therefore host screening evidence, not a measurement of GPU column
concurrency.

### Error estimate retained with each timing row

For an eligible coarse row, the benchmark constructs one refinement:

- FE doubles the forced substep count over the same interval; and
- IMEX and BE halve the timestep and apply two consecutive steps over the same
  interval.

For each of 51 fields it records normalized weighted-L1 and L-infinity
coarse-versus-refined differences. Because the qualified methods are first
order, division by `2^1-1` leaves the same difference as the first-order
Richardson estimate. The row's screening error is the maximum of both norms
over all fields.

This is an **estimated local truncation error indicator**, not a direct
Dormand--Prince error measurement for that performance row. It may be used in
Task 6 only when coarse and refined advances succeed and neither has FE
clipping or an implicit material-negative failure.

### Timing statistics and uncertainty

The committed default procedure uses:

- 10 untimed warmups;
- at least 30 and at most 100 primary samples per row;
- the sample median as the central timing;
- a fixed-seed, 10,000-resample percentile bootstrap 95% interval for the
  median; and
- early stopping only when the interval's relative half-width is less than
  5%.

Rows that reach 100 repetitions without the target precision remain reportable
but must be flagged `maximum_repetitions`; their wider interval must be
carried into comparisons. The instrumented pass uses the same number of
samples obtained for the corresponding primary row. Instrumented and primary
timings are separate distributions, so their ratio describes case-matched
overhead; it is not a paired-sample confidence interval.

Phase medians are descriptive and need not sum exactly to the median complete
advance because medians are nonlinear and instrumentation itself has cost.
The pre-remediation CSV does not serialize raw phase or fence samples, so it
cannot support phase-level conclusions. Task 5R requires raw repetition-level
serialization and independent median reconstruction, but those values remain
diagnostic. Only complete uninstrumented timing supports a local CPU
comparison.

### Direct measurements versus accounting estimates

The final tables must preserve these classifications:

| Quantity | Classification |
|---|---|
| Complete and instrumented time per column-call, repetitions, and bootstrap interval | **Measured**, after the retained run completes |
| Phase and fence samples/medians | **Diagnostic measured observations** only after the Task 5R sidecar validates; never comparative timing |
| State status, clipping count, material-negative count, factorization and solve counters, and checksums | **Executed diagnostic/accounting evidence** |
| `executed_field_level_updates` | **Deterministic count estimate** from fields, levels, batch, and method steps |
| `workspace_bytes` | **Implementation formula estimate**, not a measured peak-memory value |
| Richardson `screening_error` | **Estimated error** from coarse/refined state differences |
| Break-even and hybrid results | **Derived/estimated**, never direct timing |

## Exclusion rules

Exclusions are applied before comparative or recommendation tables. Excluded
rows may be retained as diagnostic evidence but cannot support a speedup,
break-even, hybrid, or recommendation claim.

1. A method that fails its Task 2, Task 3, or Task 4 gate is excluded in full.
   All three methods currently pass; a later regression failure would reopen
   the decision.
2. A benchmark row with invalid input, nonfinite arithmetic, solver failure,
   unsafe pivot, outside-linear-branch status, material raw negativity, or
   state disagreement between diagnostic and lean paths is excluded.
3. An FE row with clipping in either the coarse or refined advance is excluded
   from matched-error and break-even analysis. Its raw/clipped difference and
   clipping ledger remain diagnostic evidence.
4. An IMEX or BE row with a material raw negative is excluded; clipping cannot
   be introduced to repair it.
5. A row without a successful refinement and finite Richardson estimate is
   excluded from runtime-versus-error analysis even if its timing is retained.
6. No error target may use extrapolation. Interpolation is allowed only
   between adjacent measured eligible points and must be labeled derived.
7. The archived MAM4xx `1b02f272...` `top_lev=0` bounds failure is historical
   expected-failure evidence and never a passing current case.
8. A BE state that would change the retained bottom `max(0,tmpa)` branch is
   outside the approved linear solve and is excluded.
9. Partial Task 5 shards, interrupted runs, smoke runs, old/unversioned schema,
   summary/sidecar pairs lacking a valid completion record, and rows lacking
   final retained provenance or hashes cannot enter final aggregate tables.
10. CPU timing cannot support a GPU speedup, occupancy, bandwidth,
    synchronization, or production whole-model claim.
11. The nominal transported-field label cannot support field-count cost or
    break-even conclusions because every method executes all 51 fields.

The final report must include an exclusion manifest with counts by reason and
enough case keys to reproduce every exclusion.

## Pending performance results

> **PENDING - DO NOT INTERPRET PLACEHOLDERS AS RESULTS.** Complete these
> sections only from finalized Task 5 retained artifacts and the tested Task 6
> analysis. Do not transcribe partial shard values.

### Task 5 completion checks

| Check | Required | Executed result |
|---|---:|---|
| Data rows in combined required matrix | 24,300 | **PENDING** |
| FE / IMEX / BE rows | 18,900 / 2,700 / 2,700 | **PENDING** |
| Missing or duplicate case IDs | 0 / 0 | **PENDING** |
| Invalid rows | Report count and reasons | **PENDING** |
| Eligible timing rows | Report by method | **PENDING** |
| Rows meeting `<5%` CI half-width | Report count and fraction | **PENDING** |
| Rows stopped at 100 repetitions | Report count and fraction | **PENDING** |
| Nonfinite or malformed output cells | 0 except defined empty cells | **PENDING** |
| Summary schema and run ID | `ndrop_summary_v2`; one validated run ID | **PENDING** |
| Phase sidecar schema and coverage | `ndrop_phase_samples_v1`; exact repetition coverage | **PENDING** |
| Shard completion records | 135 valid records for aligned 180-case shards | **PENDING** |
| Merged summary and sidecar hashes and sizes | Exact SHA-256 and bytes for both | **PENDING** |
| Clean optimized regression | PASS required | **PENDING** |

### Typical-mixing performance

| Method | Levels | Nonzero aerosol pairs | Batch | ACI interval | FE `nsubmix` | Median time/column-call | 95% CI | Repetitions | Eligibility |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **PENDING** |  |  |  |  |  |  |  |  |  |

### Strong-mixing performance

| Method | Levels | Nonzero aerosol pairs | Batch | ACI interval | FE `nsubmix` | Median time/column-call | 95% CI | Repetitions | Eligibility |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **PENDING** |  |  |  |  |  |  |  |  |  |

### Phase, workspace, and work accounting

| Method | Setup | Factorization | Solve | Combined factorization/solve | Update | Diagnostics | Fence | Instrumentation overhead | Workspace estimate | Work-count estimate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **PENDING** |  |  |  |  |  |  |  |  |  |  |

The final narrative must use complete uninstrumented time for comparisons and
must label workspace and work counts as estimates.

### Runtime versus error

Requested screening-error targets are `1e-2`, `3e-3`, `1e-3`, `3e-4`, and
`1e-4`.

| Error target | Method | Supported measured point(s) | Error estimate | Runtime and 95% CI | Interpolation | Status/exclusion |
|---:|---|---|---:|---:|---|---|
| `1e-2` | **PENDING** |  |  |  |  |  |
| `3e-3` | **PENDING** |  |  |  |  |  |
| `1e-3` | **PENDING** |  |  |  |  |  |
| `3e-4` | **PENDING** |  |  |  |  |  |
| `1e-4` | **PENDING** |  |  |  |  |  |

If a target is not bracketed by eligible measured points, report
`unsupported by measured grid`; do not extrapolate.

### Break-even FE substep count

Break-even is the smallest measured FE `nsubmix` whose median runtime is at
least the candidate runtime for the same levels, nonzero-state occupancy,
batch, mixing profile, and ACI interval. Any interpolation must be between
adjacent measured FE points. Bootstrap uncertainty must be propagated to a
threshold interval or set of supported adjacent values. This occupancy label
must not be interpreted as a variable executed-field width.

| Candidate | Mixing | Levels | Nonzero aerosol pairs | Batch | ACI interval | Lower / point / upper break-even `nsubmix` | Eligible points used | Status |
|---|---|---:|---:|---:|---:|---|---|---|
| IMEX | **PENDING** |  |  |  |  |  |  |  |
| BE | **PENDING** |  |  |  |  |  |  |  |

No break-even or hybrid conclusion is made in this draft.

### Batching and qualitative GPU implications

The host screen can measure how complete per-column-call cost changes across
the prescribed homogeneous batch sizes. It can also report right-hand-side
reuse, factorization counts, workspace estimates, and sequential vertical or
substep depth. These are structural inputs to a GPU design decision.

The following implications are qualitative only:

- FE exposes field and column concurrency but retains sequential vertical
  passes and sequential substeps; variable `nsubmix` can create a workload
  tail.
- IMEX has two scalar matrix families and many right-hand sides. EKAT has
  device-oriented scalar tridiagonal machinery, but no production IMEX
  integration or device benchmark is part of this study.
- BE requires a coupled block solver not supplied by the current scalar EKAT
  interface; device layout, solver choice, staging, and synchronization are
  unresolved.
- The current host batch loop is serial. It cannot measure GPU occupancy,
  compaction cost, warp divergence, memory bandwidth, or host/device staging.

Any FE/IMEX threshold or packing analysis belongs to Task 6 and must be labeled
an estimate. It cannot be reported as an implemented hybrid or a measured GPU
result.

## Uncertainty treatment for final analysis

The final analysis must preserve distinct uncertainty sources:

- **Timing sampling:** report the median, 95% bootstrap interval, relative
  half-width, repetition count, and stop reason for every summarized row.
- **Error estimation:** retain both per-field Richardson norms, their maximum,
  the coarse/refined checksums, and the fact that a single first-order
  refinement is an estimate rather than an independent reference solution.
- **Grid support:** report whether an error target or break-even threshold is
  directly measured, bracketed by adjacent points, or unsupported. No
  extrapolation is allowed.
- **Exclusion uncertainty:** show how many rows are removed for clipping,
  positivity, invalid status, failed refinement, incomplete evidence, or
  precision failure. Do not silently drop them.
- **Hybrid sensitivity:** vary the FE/IMEX threshold, packing/compaction
  overhead, and workload distribution only as estimates. Present a range, not
  a production speedup.
- **Representativeness:** separate variation within the synthetic benchmark
  from uncertainty caused by unknown production distributions. Bootstrap
  timing intervals do not cover the latter.

Accuracy-gate outputs are deterministic numerical results and do not receive
statistical confidence intervals. Their limitations arise from case coverage,
the frozen operator, and floating-point/backend scope.

## Final recommendation framework

> **FINAL RECOMMENDATION: PENDING.** No method is selected in this draft.

After Tasks 5 and 6 complete, choose exactly one of:

1. current FE;
2. IMEX Euler;
3. Backward Euler;
4. an FE/IMEX hybrid; or
5. insufficient evidence.

The choice must be supported by eligible matched-error points, timing
uncertainty, break-even behavior across levels and tested nonzero-state
occupancies, batching sensitivity, implementation complexity, unresolved
physics, and the production measurements still required. The current screen
cannot establish executed-field-count scaling. Equal-interval timing, a coarse
one-call accuracy maximum, an older provisional cost model, or CPU-to-GPU
analogy is insufficient by itself.

The completed recommendation section must include:

- one sentence naming exactly one approved choice;
- the measured and derived evidence that controls the choice;
- the regimes for which the choice is and is not supported;
- exclusions and uncertainty that could reverse it;
- the next separately approved implementation or measurement gate; and
- an explicit statement that no production or GPU change is authorized by
  this standalone report.

## Production evidence needed to revise offline estimates

The following production measurements are needed before translating the local
screen into a production implementation or whole-model claim:

1. Per-call and per-physics-step distributions of actual `nsubmix`, including
   counts, sums, maxima, quantiles, and the fraction of columns at each
   prospective threshold.
2. `dt_phys`, configured MAC/ACI subcycles, derived ACI interval, active-level
   count and first active level, and actual transported-field configuration.
3. Sampled pressure weights, cloud fractions/overlap, interface exchange,
   activation rates, and separate physical versus first-level controller
   limits, sufficient to replay representative and tail columns.
4. FE raw minima, per-field clipping counts and weighted deficits, and signed
   total/operator/clipping aerosol ledgers within the call rather than only at
   history-output times.
5. MAM4xx ACI and `dropmixnuc` time and call counts separated from TMS, SHOC,
   cloud fraction, and P3, plus the fraction of total model time available for
   improvement.
6. Column workload heterogeneity and threshold packing/compaction cost on a
   representative GPU, including occupancy, divergence, bandwidth,
   synchronization, workspace, and host/device staging.
7. Candidate raw-positivity and solver-failure rates on replayed production
   columns; clipping cannot be assumed for IMEX or BE.
8. Short coupled scientific comparisons of aerosol partitioning, droplet
   number, cloud fraction, precipitation, and relevant ACI/P3 responses,
   followed by restart and supported-machine regression evidence.
9. A physics-owner decision on the retained current-C++ versus Fortran bottom
   activation formulation before any production adoption.

## Provenance

### Frozen dependencies

The current stack is the only stack authorized for this screen:

| Repository | Locked revision |
|---|---|
| E3SM baseline | `e2c578e40f6602a108eff577091de86a453624ab` |
| MAM4xx | `0d076d8ebb004a18f05f51f97b08edd9def6826f` |
| EKAT | `3a1485aecfed7a036e4fa30ce845498fd37975e2` |
| Kokkos | `23745d8494203711a6090d0a201c3ae3bf57853b` |

The archived 2026-08-06 stack is historical only and is not substituted for
the current gate or performance run.

### Committed evidence inputs

| Artifact | Revision |
|---|---|
| Approved design synchronized through Task 4 | `7c9e187c5305e082f69c8159b5c2a9aeca2658bb` |
| Approved execution plan | `e449075196df85c0ce11f1c25f997390409c17d7` |
| FE/IMEX gate review | `166ed9bae15edde0d61c01fc82727635744a542a` |
| BE gate review | `f2a9cd800b4543d27690173b6f23bac78d07e947` |
| Cross-scheme gate code/test | `7807aedfdecafe5e30b783e98e81ac4b723f56e3` |
| Cross-scheme retained review and design sync | `7c9e187c5305e082f69c8159b5c2a9aeca2658bb` |
| Benchmark implementation | `5882c458ee` |
| Matched-error estimate fields and tests | `dba6d00598d86a239a2bf2743d95325821c9378c` |

### Executed Task 4 environment

Task 4 used double precision, Apple clang 21.0.0, CMake 4.4.2, Ninja,
Kokkos 4.5.1 Serial, and an arm64 Apple M2 Max running macOS 26.6.1. The debug
configuration enabled ASan and UBSan with a Kokkos Debug build; the optimized
configuration used Release `-O3 -DNDEBUG`. Harness targets compiled with
`-Wall -Wextra -Wpedantic -Werror`. This environment description applies to
the retained accuracy gate, not automatically to the pending Task 5 run.

### Retained Task 4 artifacts

| Record | SHA-256 |
|---|---|
| `cross-scheme-debug.csv` and `cross-scheme-release.csv` | `b28f376cb377e25d5ea04156063a04b48553015ed995813d9b8c8db15cd98f44` |
| `cross-scheme-debug.md` and `cross-scheme-release.md` | `14b7e811ac4dc2b0bd56db24efcc0d259f50c4519424b3fd5195021140f91163` |
| `cross-scheme-debug.log` and `cross-scheme-release.log` | `c6e7e1b57f98ed2fca4d5a684e37c6a58ae58bf23aa9971565b75775e5006aaf` |
| `ctest-debug-aux.log` | `73a046b438768da95a367c43b94c3ac5f805adceb8831e243cff41a4ecd7361c` |
| `ctest-release.log` | `c1c1be482d54d1cec6992153d7dfaff184fe7f6fbfd484b1e2e5693cfe13a7d4` |
| `task4-cross-scheme-debug.json` | `7906f643ec7ee76d8616a599233e589e85eeddeccb771ad8e0794eae8575f8c3` |
| `task4-cross-scheme-release.json` | `d0cd64373c0e58eecd296b80ff1a06de90d2f3a59de17b6e5a6424d1dae65bc1` |

### Task 5 and Task 6 provenance to fill

Do not finalize performance or recommendation sections until all fields below
are present or explicitly marked unavailable with a reason.

| Provenance field | Final value |
|---|---|
| Clean detached source commit used for optimized execution | **PENDING** |
| Source tree and dependency dirty-state checks | **PENDING** |
| Locked-dependency provenance JSON path and SHA-256 | **PENDING** |
| Compiler identity/version and complete flags | **PENDING** |
| CMake version, generator, build type, and exact configure/build commands | **PENDING** |
| Kokkos version, execution space, and precision | **PENDING** |
| Machine model, CPU/architecture, hardware threads, OS version | **PENDING** |
| Exact benchmark invocation(s), shard ranges, start/end time, and exit status | **PENDING** |
| Benchmark executable SHA-256 | **PENDING** |
| Identity-manifest path and SHA-256 | **PENDING** |
| Raw summary/sidecar shard inventory with sizes, row counts, ranges, and SHA-256 | **PENDING** |
| All 135 completion-record paths, counts, and SHA-256 values | **PENDING** |
| Merged summary and sidecar paths, schema/header hashes, counts, sizes, and SHA-256 values | **PENDING** |
| Final immutable manifest path and SHA-256 | **PENDING** |
| Smoke/support/regression logs, commands, counts, and SHA-256 | **PENDING** |
| Task 6 analysis source commit and tests | **PENDING** |
| Exact analysis invocation, input hashes, and output hashes | **PENDING** |
| Exclusion manifest path, counts, and SHA-256 | **PENDING** |
| Independent final clean-checkout review and decision | **PENDING** |

## Limitations and evidence boundaries

- The operator, coefficients, cloud mask, and geometry are frozen over one
  ACI call. Evolving coefficients, outer process splitting, and coupling to
  other physics can dominate whole-model error.
- The deterministic matrix is synthetic. It spans required dimensions and
  regimes but does not establish their production frequency or correlation.
- All retained accuracy evidence is local double-precision Kokkos Serial CPU
  evidence. Portability and numerical behavior on supported CPU and GPU
  backends remain untested here.
- The BE solver is host-only and test-only. The current GPU scalar solver does
  not imply availability or performance of a 2-by-2 block solver.
- The reference independently assembles the continuous RHS and integration but
  shares the validated input structure and common zero-cloud transition.
- The current C++ bottom source differs from the inspected Fortran
  formulation. The screen retains the current C++ behavior and cannot resolve
  the intended physics.
- Task 4 found no clipping or material negatives in its synthetic matrix. This
  does not prove their absence in production distributions.
- The performance screen measures the isolated integrator on one Mac. It does
  not include the rest of MAM4xx ACI, other MAC-group processes, MPI, GPU
  kernels, I/O, or coupled feedbacks.
- Bootstrap intervals quantify timing-sample variation for one case and
  environment. They do not quantify production representativeness,
  implementation risk, or scientific uncertainty.
- Workspace and work-count columns are formulas. They are not peak-memory,
  cache, bandwidth, energy, or hardware-counter measurements.
- The pre-remediation `active_pairs` input changes nonzero-state occupancy,
  but every method executes all 25 aerosol pairs and 51 fields. Task 5R will
  correct the labels and preserve that fixed width; no transported-field-count
  performance conclusion is supported.
- The pre-remediation artifact records IMEX's combined factorization and solve
  work in the solve phase. Task 5R will expose a combined diagnostic phase and
  leave separate values unavailable, but cannot decompose the operation.
- Raw phase and fence timing samples are absent from pre-remediation output.
  Task 5R requires a validated sidecar, yet even validated phase data remain
  diagnostic and unavailable to comparative performance functions.
- Earlier provisional cost ratios and hybrid ideas are hypotheses only and do
  not enter the final decision without the completed retained analysis.

## Revision history

| Date | Author/agent | Summary |
|---|---|---|
| 2026-08-25 | Codex | Recorded the pending Task 5R evidence-remediation contract: fixed executed width with occupancy-only labels, combined IMEX timing, raw phase sidecar, transactional manifests, and primary-timing-only comparative analysis. No performance conclusion was added. |
| 2026-08-24 | Codex | Integrated the static Task 5 review limits: reclassified active-pair cases as nonzero-state occupancy, prohibited field-count scaling, and marked IMEX phase separation and raw phase/fence retention as unresolved without adding performance conclusions. |
| 2026-08-24 | Codex | Drafted the Task 7 report structure from approved design, committed Task 2--5 sources, and retained Task 2--4 reviews; recorded executed Task 2--4 numerical evidence and left Task 5 timing, Task 6 analysis, break-even/hybrid results, and the final recommendation explicitly pending. |
