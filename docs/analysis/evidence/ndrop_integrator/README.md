# Isolated `ndrop` Integrator Harness

This directory contains test-only code for the approved standalone numerical
experiment in `docs/design/eamxx_mam4xx_ndrop.md`. It is not linked into E3SM
or MAM4xx production targets.

## Physically scaled common-accuracy column extension

### Natural-FE-accuracy comparison on the same physical columns

The separate `ndrop_physical_controller_*_v1` workflow applies the existing
controller-matched policy to unchanged fixtures 28--31. FE uses only its
natural adaptive controller count (not accuracy-refined work). Its achieved
maximum 51-field, both-norm, both-reference error defines each case's target;
IMEX and BE select the first eligible positive integer count meeting it.
This answers whether another method can match FE's baseline accuracy at lower
cost, not whether that baseline accuracy is scientifically sufficient.

Prepare and freeze untimed evidence in a new local scratch directory:

```bash
python3 -B ndrop_controller_matched_artifacts.py run-accuracy --physical-columns \
  --source-root /path/to/E3SM --binary /path/to/benchmark_ndrop_matched_error \
  --output-dir /Users/burr114/e3sm-scratch/TASK/accuracy
python3 -B ndrop_controller_matched_artifacts.py freeze \
  --output-dir /Users/burr114/e3sm-scratch/TASK/accuracy
```

The copied executable and four repository grids, design/source/build identity,
physical metadata, complete lower-count searches, numerical gates, selection,
and independent `N,2N,4N` timing-admission evidence are hash-bound. Failures
remain visible; missing confirmations exclude timing without reselecting work.
The earlier fixed-ceiling physical and controller-stratified synthetic studies
are separate and unchanged. `mixing_scale` in this physical-controller CSV is
background diffusivity [m2 s-1], not a dimensionless multiplier.

Stop after reviewing the sealed selection. `run-timing` remains subject to
separate explicit approval and uses only diagnostics-disabled complete
advances plus required fences, ten warmups, 30--100 paired blocks and the
unchanged stopping/bootstrap protocol. Synthetic timing fixtures in unit
tests are not measured evidence. These columns remain idealized standalone
shared-exchange operators; they neither resolve directional-rate production
fidelity nor provide climatological weights or a production recommendation.

### Paired high-work geometry sensitivity (untimed stage)

The separately approved `--physical-geometry` policy includes native controls
28--31 and eight conservative cloud-band subdivisions. Whole native cells
intersecting 400--1000 m are split equally in log pressure; endpoints and the
rest of the full column remain unchanged. The actual level counts are
L72-derived 79/93 and L128-derived **141/167**, correcting the proposal's
140/164 labels without changing its band predicate. These are prescribed-
turbulence resolution sensitivities, not either model's current grid.

Each mesh has identical moderate/strong initial state and geometry, with the
existing K=3/30 m2/s, cloud-base wtke=0.3/1.0 m/s, 150 s interval and activation
construction. Dry mass and all 51 parent-cell state integrals are conserved
under subdivision. Cloud/droplet amounts are distributed by child cloudy mass
with clear-child zeros; interstitial mixing ratios are copied. Refined state
concentrations are conservative representations, not fresh point samples.
All integrators, controller, accuracy and convergence gates are unchanged.

```bash
python3 -B ndrop_controller_matched_artifacts.py run-accuracy --physical-geometry \
  --source-root /path/to/E3SM --binary /path/to/benchmark_ndrop_matched_error \
  --output-dir /Users/burr114/e3sm-scratch/TASK/geometry-accuracy
python3 -B ndrop_controller_matched_artifacts.py freeze \
  --output-dir /Users/burr114/e3sm-scratch/TASK/geometry-accuracy
```

The separate closed `ndrop_physical_geometry_controller_*_v1` identity and CSV
bind all 12 cases, parent grid and actual extents, subdivision and band,
limiting-cell height/spacing and adjacent diffusivities/cloud-base predicates,
natural FE counts, achieved errors, all lower implicit failures, references,
copied input hashes and completion records. Numerical timing admission alone
does not grant execution authority. The user separately approved one bounded
local timing stage on 2026-10-12. Geometry `run-timing` and reporting require a
hash-bound ordinary `timing_approval.json`, unchanged approved accuracy CSV,
and freshly sealed selections from the exact timing binary. Missing or stale
approval and execution records fail closed; a partial/previous timing attempt
blocks another run in that directory. The scope and exact closed records are
defined in the design document's bounded timing-stage subsection.
Keep the native-only and prior synthetic artifacts separate. Neither tail
counts nor an implicit advantage are guaranteed; do not tune inputs after
selection. This does not resolve the deferred directional-rate representation.

After the approved timing record is sealed, the bounded stage uses the same
driver's `run-timing`, `report`, then `validate-final`. The geometry executable
mode is `--physical-geometry-controller-timing`; ten warmups and 30--100 complete
paired blocks use the existing 5% precision rule and 7200 s bound. Only the
complete diagnostics-disabled advance plus required fence is measured. The
distinct geometry analysis preserves actual/parent grid sizes, refinement
and limiter metadata, all achieved errors, raw paired timings, fixed 10,000-
resample ratio intervals, command/count evidence and final hashes. Its scope
is unweighted local idealization, never production frequency or model speedup.

### Existing fixed 1%/0.1% ceiling comparison

`ndrop_matched_artifacts.py run-accuracy --physical` prepares four separate
DYCOMS-inspired idealizations (fixtures 28--31): moderate/strong turbulence
at repository L72/L128 nonuniform coordinate families and a common 150 s
interval. It copies and hashes the exact `acme-72{i,m}.ascii` and
`scream-128{i,m}.ascii` inputs before measurement. The design document's
physical-column extension freezes the dimensional geometry/exchange recipe,
600--840 m cloud layer, physically scaled modal mass/number inputs, prescribed
activation fractions, and turbulence sensitivity values. Geometry is an
explicit isothermal approximation; this is neither a complete DYCOMS sounding
nor evidence that those turbulence profiles are observed climatologies.

The natural adaptive FE baseline is retained separately. Common error ceilings
are 1% and 0.1%, using the existing maximum of all 51 fieldwise normalized
weighted-L1/L-infinity errors against both certified references. Every integer
count is searched independently; FE starts at its controller-selected count.
Larger FE counts are **accuracy-refined**, not production controller work.
Selections are accuracy-only; `N,2N,4N` certification is required before timing.
Missing/failed/non-decreasing confirmation withholds timing but remains visible
in selection. All lower failures, both-reference errors, physical stiffness,
limiting levels/spacings and exclusion reasons are retained.

Example untimed workflow (new task-specific scratch paths required):

```bash
python3 -B ndrop_matched_artifacts.py run-accuracy --physical \
  --source-root /path/to/E3SM --binary /path/to/build/benchmark_ndrop_matched_error \
  --output-dir /Users/burr114/e3sm-scratch/<task>/accuracy
python3 -B ndrop_matched_artifacts.py freeze \
  --output-dir /Users/burr114/e3sm-scratch/<task>/accuracy
python3 -B ndrop_matched_artifacts.py validate-accuracy \
  --output-dir /Users/burr114/e3sm-scratch/<task>/accuracy
```

The untimed CTest `ndrop_physical_column_fixtures` checks repository-grid and
paired-input fidelity, no-flux/clear-cell behavior, controller use, and bitwise
diagnostics-on/off numerical states. Supply `NDROP_PHYSICAL_TEST_BINARY` to the
Python unit suite to execute the additional grid corruption fixtures. Select
numerical/support tests explicitly; do not run benchmark smoke tests through a
blanket CTest invocation.

After **separate timed-run approval**, the existing `run-timing`, `report` and
`validate-final` commands recognize only hash-bound
`ndrop_physical_column_*_v1` evidence. Timing includes the complete
uninstrumented advance, uses adaptive FE for its baseline count and the existing
forced-count test entry point only for refined FE, and preserves the existing
complete paired-block/bootstrap protocol. No timing, speedup, production replay
fidelity, climatological frequency or method recommendation follows from the
untimed preparation. Old synthetic artifacts and production captures are not
admitted or reinterpreted. The directional-rate representation issue remains
deferred; no production source, numerical implementation, or dependency changes
are part of this extension.

## Controlled moderate/strong synthetic mixing extension

The separately approved `controlled_mixing_v1` fixture policy adds four cases,
without changing fixtures 0--23 or their artifacts:

| ID | Levels | Mixing scale | Activation-input scale | Interval |
|---:|---:|---:|---:|---:|
| 24 | 72 | 0.05 (moderate) | 0.05 | 150 s |
| 25 | 72 | 0.5 (strong) | 0.05 | 150 s |
| 26 | 128 | 0.05 (moderate) | 0.05 | 150 s |
| 27 | 128 | 0.5 (strong) | 0.05 | 150 s |

Within each grid pair, the initial 51-field state, geometry, cloud fractions,
and all modal number/mass activation inputs are identical. Activation arrays
are copied from the existing moderate-scale generator; mixing exchanges and
the controller-only first-level input use the case's mixing scale. Ordinary
operator preparation and limiting still run. Unexpected limiting or changed
prepared activation invalidates the fixture. These labels are synthetic,
not calibrated atmospheric mixing regimes or production frequencies.

FE executes its actual production-controller count, not a requested bin.
Each case's certified FE error defines its own target. Consecutive IMEX/BE
integer selection, the 8192 cap, two-reference check, numerical gates, and
`N,2N,4N` timing admission are unchanged. Error matching is within each case;
the targets can differ across mixing regimes. Side-by-side costs do not show
a common-error cross-regime speedup or fieldwise dominance.

Use a new unsynced scratch directory. The following commands perform only
untimed accuracy/refinement and selection preparation, authorized separately
from any timing pass (substitute the approved exact build/output paths):

```bash
python3 -B ndrop_controller_matched_artifacts.py run-accuracy \
  --controlled-mixing --source-root /path/to/E3SM \
  --binary /path/to/approved-build/benchmark_ndrop_matched_error \
  --output-dir /Users/burr114/e3sm-scratch/TASK/accuracy
python3 -B ndrop_controller_matched_artifacts.py freeze \
  --output-dir /Users/burr114/e3sm-scratch/TASK/accuracy
```

This policy requires exactly `[24,25,26,27]`; partial or mixed populations
are rejected. All eight artifact kinds use `ndrop_controlled_mixing_*_v1`,
with a bound `fixture_policy` and `activation_mixing_scale=0.05`. Accuracy
records add that column immediately after `mixing_scale`; analysis fixture
entries retain it. The existing source/binary/dependency, completion/hash,
selection, paired-block, timing-state-checksum, and exclusion checks apply.

`ndrop_controlled_mixing_fixtures` is an untimed CTest entry invoking
`--check-controlled-mixing-fixtures`; it checks fixed inputs, normal limiting,
actual adaptive FE counts, and bitwise identical enabled/disabled diagnostic
states for all methods. Select this and the numerical/support unit tests
explicitly; do not run the benchmark smoke/exclusion tests inadvertently.
Python corruption fixtures use fabricated values, not measured performance.

Timing software supports this policy, but a timing run and its performance
analysis require separate approval after the selection/admission artifacts
validate. Only the complete uninstrumented advance may supply primary time,
with ten warmups and 30--100 complete interleaved blocks. Report all three
within-case paired time ratios, errors, counts, refinement evidence, and
exclusions. No production replay, retained Task 5 access, remote/HPC work, or
production/model-wide recommendation is authorized by this extension.

## Locked dependency workflow

`dependency_lock.json` records two immutable stacks. `current` is the only
stack authorized for the 2026-08-18 cross-method screen; `archived` preserves
the 2026-08-06 FE/IMEX evidence and must not be substituted for the current
gate. Neither workflow changes the production gitlinks in this E3SM worktree.

Prepare clean detached checkouts below the local scratch root and write the
machine-readable provenance record to the retained-results root:

```bash
./prepare_locked_dependencies.sh \
  --stack current \
  --root /Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current \
  --results-root /Users/burr114/e3sm-results/ndrop-screening \
  --record /Users/burr114/e3sm-results/ndrop-screening/provenance/dependencies-current.json \
  --build-type Release \
  --generator Ninja \
  --cmake-flags='-O3 -DNDEBUG'
```

The script obtains the locked E3SM base from the local repository, obtains the
locked MAM4xx, EKAT, and Kokkos revisions in separate directories, verifies
the E3SM and EKAT gitlinks, and rejects a SHA mismatch, attached HEAD, or dirty
checkout without resetting it. Its JSON record contains locked and actual
SHAs, dirty state, the source-worktree revision and dirty state, compiler and
CMake identity, build settings, CPU/OS identity, and the exact invocation.
Rerun with `--verify-only` before each retained gate or timing run.

The preparation script itself has a bounded behavior test that uses temporary
local repositories only:

```bash
./test_prepare_locked_dependencies.sh
```

It proves clean success and fail-closed behavior for dirty and SHA-mismatched
checkouts. The test owns and removes only its `mktemp` directory below
`/Users/burr114/e3sm-scratch/` and its matching temporary provenance directory
below `/Users/burr114/e3sm-results/`.

The first gate implements the canonical shared-interface operator and the
current Forward Euler subcycling, activation-rate limiting, alternating work
buffers, and componentwise positivity clipping. The test suite covers all 25
aerosol constituents and droplet number.

Build and run outside the OneDrive source tree. `Kokkos_DIR` must identify a
Serial Kokkos build matching the selected locked EKAT revision. The current
adapter uses the locked current MAM4xx rank-3-view API and does not depend on
Haero. Archived Haero-based outputs remain historical artifacts only.

```bash
cmake -S /path/to/E3SM/docs/analysis/evidence/ndrop_integrator \
  -B /Users/burr114/e3sm-scratch/ndrop-screening/build/debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DKokkos_DIR=/Users/burr114/e3sm-scratch/ndrop-screening/kokkos-install-debug/lib/cmake/Kokkos \
  -DEKAT_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/ekat \
  -DMAM4XX_SOURCE=/Users/burr114/e3sm-scratch/ndrop-screening/dependencies/current/mam4xx \
  -DNDROP_ENABLE_SANITIZERS=ON \
  -DNDROP_BUILD_EXACT_ADAPTER=ON
cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/debug --parallel 4
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/debug \
  --output-on-failure
```

Use a release Kokkos build, `-DCMAKE_BUILD_TYPE=Release`, and
`-DNDROP_ENABLE_SANITIZERS=OFF` for the optimized gate. The checked-in
`CMakeLists.txt` applies warnings-as-errors to all harness sources and registers
the Forward Euler, IMEX, Backward Euler, common cross-scheme, and exact-current
executables with CTest.

## Separate synthetic matched-error interval study

`benchmark_ndrop_matched_error` implements the 2026-10-07 standalone interval
study. It shares the existing fixed synthetic generator, uses the independently
assembled two-tolerance DOPRI54 reference to check the state actually timed,
and lets FE, IMEX, and BE take different internal steps over one frozen ACI
interval. Its artifact schemas are distinct from Task 5R and from production
replay. Source with a dirty or attached HEAD is labeled exploratory; the tool
does not promote that evidence to a retained result.

The four fixtures are L72/L128 crossed with `(mixing_scale,T)=(0.05,50 s)`
and `(0.5,150 s)`, with all 25 pairs nonzero. `run-accuracy` computes the
dyadic ladder through 8192 without a timing pass. `freeze` verifies the exact
matrix and fieldwise errors, chooses the least qualifying tested count for each
method and tolerance, and writes a separate immutable selection. `run-timing`
uses only those frozen configurations, with 10 warmups, 30--100 complete
interleaved blocks, and a 7200 s study-wide cap. `report` validates completion
records and computes paired-block bootstrap median-time ratios. A missing pair
is reported as unavailable.

Use a new local, unsynced scratch directory for each run. The commands below
are illustrative; build and binary paths must name the exact optimized binary
that passed the numerical suite, and the selected dependency checkouts must
match the locked current stack:

```bash
cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/release \
  --parallel 4 --target benchmark_ndrop_matched_error test_cross_scheme \
  test_benchmark_support
python3 ndrop_matched_artifacts.py run-accuracy \
  --source-root /path/to/E3SM \
  --binary /Users/burr114/e3sm-scratch/ndrop-screening/build/release/benchmark_ndrop_matched_error \
  --output-dir /Users/burr114/e3sm-scratch/ndrop-screening/matched-error/RUN
python3 ndrop_matched_artifacts.py freeze \
  --output-dir /Users/burr114/e3sm-scratch/ndrop-screening/matched-error/RUN
python3 ndrop_matched_artifacts.py run-timing \
  --output-dir /Users/burr114/e3sm-scratch/ndrop-screening/matched-error/RUN
python3 ndrop_matched_artifacts.py report \
  --output-dir /Users/burr114/e3sm-scratch/ndrop-screening/matched-error/RUN
```

Only the uninstrumented complete advance supplies primary timing. The
independent reference and detailed validation work are outside the timed
region. The direct comparison is conditional on these frozen synthetic
operators; physical-column population, GPU, and whole-model conclusions need
their separate evidence.

### Production-controller-matched synthetic extension

The separately versioned `ndrop_controller_matched_*_v1` study uses the exact
substep count returned by the existing Forward Euler production controller as
each fixture's baseline. Forward Euler runs through its adaptive entry point;
IMEX and Backward Euler are restarted from the same state at consecutive
integer counts and stop at the first numerically eligible result whose maximum
two-reference field error is no greater than the Forward Euler error. The
accuracy search streams results rather than retaining candidate states, and
confirmation records at `2N` and `4N` cannot alter that selection. For the
bounded pilot, all three methods retain these refinements. The separately
hashed `timing_admission.json` withholds timing unless both aggregate norms
strictly decrease against both references at `N,2N,4N`, with finite, eligible
states throughout. Observed orders are reported without imposing another
order window; formal first order remains the unchanged numerical-unit gate.
Missing or failed refinements exclude a configuration, never select a new
count, relax an error target, or redistribute a missing pair.

The matrix contains the four earlier L72/L128 continuity fixtures and twenty
unweighted controller-stratified synthetic fixtures targeting exact returned
counts `1,2,4,8,16,32,64,128,256,512` at both level counts. It is a workload
diagnostic, not a climatology or production-state sample. Accuracy and frozen
selection are prepared before any timing authorization:

`run-accuracy --fixtures 4,5,10,11,18,19,20,21` freezes a bounded subset at
L72/L128 and exact controller counts 1, 8, 128, 256. Omit `--fixtures` for
the complete 24-fixture matrix. These bins are weak/moderate/stiff synthetic
workloads, not estimates of typical production counts. The 256 bin is
deliberately constructed, not a search saturation or controller cap.

```bash
python3 ndrop_controller_matched_artifacts.py run-accuracy \
  --source-root /path/to/E3SM \
  --binary /path/to/optimized/benchmark_ndrop_matched_error \
  --output-dir /Users/burr114/e3sm-scratch/ndrop-controller-matched/RUN
python3 ndrop_controller_matched_artifacts.py freeze \
  --output-dir /Users/burr114/e3sm-scratch/ndrop-controller-matched/RUN
```

Stop after inspecting the sealed accuracy and selection artifacts. A separate
explicit approval is required before `run-timing`; neither the accuracy pass
nor its selection authorizes a timed benchmark, performance analysis, or a
production claim. When separately authorized, timing uses only the complete,
uninstrumented advances already selected and the existing interleaved-block
bootstrap protocol.

Each analysis retains fixture/stiffness inputs, both reference error vectors
in the hash-bound accuracy file, the self-check and tolerances, controller and
implicit counts, refinement orders/admission exclusions, raw paired blocks,
median times and all three paired time ratios with intervals. A ratio FE/IMEX
greater than one means FE took more elapsed time at the declared maximum-field
error target; matching this maximum does not imply fieldwise dominance.
Ten warmups and 30--100 complete blocks are bound before timing. The stopping
test runs between complete blocks, outside the complete-advance timer, and
requires every timed method's median interval half-width below 5%; capped
results retain their actual block counts and stop reason. These are local,
exploratory CPU observations from synthetic frozen operators, not retained
Task 5 evidence, production-replay validation, or a method recommendation.

## Task 5R paired evidence workflow

The approved Task 5R implementation emits `ndrop_summary_v2` plus the
repetition-level `ndrop_phase_samples_v1` sidecar. The varying case dimension
is nonzero-state occupancy; every method still executes 25 aerosol pairs and
51 transport fields. IMEX records coefficient-view copy/allocation and both
public EKAT Thomas operations as one combined diagnostic phase. Only the
uninstrumented complete-advance timing may enter comparative analysis.

`ndrop_artifacts.py` is the fail-closed retained-run entry point. Run it only
from a clean detached source checkout. `create-identity` copies the exact
optimized executable and locked-dependency provenance record into a new
results directory below `/Users/burr114/e3sm-results/`, binds every one of the
135 commands to those retained copies, and seals `identity.json`. A typical
invocation has this form; replace the uppercase example tokens, and require the
run ID's final 12 lowercase hexadecimal characters to equal the clean detached
source SHA prefix:

```bash
python3 ndrop_artifacts.py create-identity \
  --run-id task5-remediated-20260826T000000Z-0123456789ab \
  --source-root /Users/burr114/e3sm-scratch/ndrop-screening/worktrees/CLEAN_SOURCE \
  --binary /Users/burr114/e3sm-scratch/ndrop-screening/build/RELEASE_BUILD/benchmark_ndrop_integrators \
  --dependency-provenance /Users/burr114/e3sm-results/ndrop-screening/provenance/dependencies-current.json \
  --compiler 'AppleClang <version>' --compiler-flag=-O3 \
  --compiler-flag=-DNDEBUG --machine '<machine>' \
  --platform '<platform>' --execution-space Serial \
  --validation-command '<exact completed validation command>' \
  --scratch-root /Users/burr114/e3sm-scratch/ndrop-screening \
  --output-root /Users/burr114/e3sm-results/ndrop-screening/RUN_DIRECTORY
```

Run each exact shard start (`0,180,...,24120`) with `run-shard`, then merge and
revalidate the authoritative final manifest:

```bash
python3 ndrop_artifacts.py run-shard \
  --identity /Users/burr114/e3sm-results/ndrop-screening/RUN_DIRECTORY/identity.json \
  --case-start ALIGNED_CASE_START
python3 ndrop_artifacts.py merge \
  --identity /Users/burr114/e3sm-results/ndrop-screening/RUN_DIRECTORY/identity.json
python3 ndrop_artifacts.py validate-final \
  --identity /Users/burr114/e3sm-results/ndrop-screening/RUN_DIRECTORY/identity.json \
  --final-manifest /Users/burr114/e3sm-results/ndrop-screening/RUN_DIRECTORY/final_manifest.json
```

An interrupted or malformed pair has no valid completion record and therefore
does not exist for analysis. No command overwrites retained input, shard,
completion, merge, or final-manifest artifacts.

The 2026-08-25 pre-retained implementation review passed 35 Python tests,
13/13 debug ASan/UBSan tests, 13/13 optimized tests, and the five-test optimized
numerical gate. The reviewed implementation ends at `d55ef90e46`. These are
validation results, not retained timing evidence; the clean-detached run and
all performance conclusions remain pending.

The current adapter validates the repaired current-MAM4xx `top_lev=0` path.
The archived `1b02f272...` `top_lev=0` expected failure and its historical
`enable_aero_vertical_mix=false` defect are explicitly excluded from the
current passing gate.

`compare_mam4xx_forward_euler.cpp` is the separate exact-current-function
adapter. It compares the canonical harness with the unmodified
`mam4::ndrop::update_from_explmix` at `MAM4XX_NUM_VERTICAL_LEVELS=128` and
`top_lev=0,1,6`. It directly compares the returned `nsubmix`, covers production
substep counts 1, 2, 3, and 5, and includes a first-interface-only profile in
which the controller changes `nsubmix` from one to two without entering the
physical tendency. The other cases cover full and fractional overlap,
limiter-active and limiter-inactive activation, and isolated, contiguous, and
all-zero cloud profiles. One case assigns distinct number and mass activation
rates to every mode so a mode-mapping permutation cannot be hidden by equal
coefficients. The destination work array is filled with unrelated values
before every exact call to detect incomplete writes.

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
The harness implements a host-only 2-by-2 block-Thomas Backward Euler baseline
with an independent dense oracle. It does not implement a GPU block solver,
EKAT's GPU cyclic-reduction path, a hybrid FE/IMEX policy, or any higher-order
implicit method. Those are separate, unvalidated stages and are not implied by
the CPU correctness results below.

## Gate validation status

### Executed common cross-scheme gate: 2026-08-24

Task 4 implemented the independent adaptive Dormand--Prince 5(4) physical-RHS
oracle and common FE/IMEX/BE gate at commit `7807aedfde`. The complete exact
Cartesian product contains 1,875 cases and 5,625 method records across level
counts `1,2,3,127,128`, overlap `0,0.25,0.5,0.75,1`, activation strengths
`0,0.1,0.5,0.9,2`, active aerosol-pair counts `1,4,8,16,25`, and zero,
typical, and strong mixing.

The debug/ASan/UBSan and optimized Kokkos Serial executables each passed
11,880 checks and produced bitwise-identical CSV and Markdown results. FE,
IMEX, and BE each passed all 1,875 matrix cases, all-field first-order slope
and `R^2` gates, all three pairwise Richardson common-solution bounds, residual
and ledger tolerances, raw positivity, and three-run bitwise determinism. No
matrix point produced a material negative or FE clipping. The maximum reference
change between `rtol=1e-12` and `1e-13` is `2.0331196e-13`.

Run the gate directly to retain its 50-column machine-readable matrix and
human-readable summary:

```bash
/path/to/build/test_cross_scheme \
  --csv /Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme.csv \
  --table /Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme.md
```

Per-case fields that are not defined for a method are empty rather than
nonfinite sentinels. The final CSV has one header, 5,625 method-case rows, and
three convergence rows; semicolon-delimited arrays retain all 51 field errors
and all 25 aerosol-pair diagnostics. The detailed design review, TDD record,
results, hashes, limitations, and separate method decisions are in
`cross_scheme_gate_review_20260824.md`. Retained artifacts are below
`/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme` and
`/Users/burr114/e3sm-results/ndrop-screening/provenance`.

The largest coarse one-call Linf error is `0.1129882`. The Task 4 PASS therefore
establishes reference consistency, first-order convergence, and a common
solution, not that one full-interval candidate step is already accurate enough
for production. Matched-error cost remains a Task 6 decision.

The current gate artifacts include:

- exact production zero-cloud handling after each complete integration call;
- checked derived coefficients, activation sources, state updates, matrix and
  right-hand-side values, solutions, residuals, ledgers, and normalizations;
- scaled matrix residuals and a separately assembled physical-equation
  residual that does not reuse the tridiagonal matrix;
- droplet and per-constituent interstitial/cloud-borne raw minima, clipping
  counts, weighted clipping changes, and ledger effects;
- all-51-field analytic diffusion and activation convergence;
- a nonuniform, fractional-overlap, all-field coupled FE/IMEX
  self-convergence and common-limit test;
- randomized weak-to-strong production-adaptive FE tests confirming that the
  approved substep and activation-rate limits keep raw updates nonnegative;
- independent, deliberately distinct modal activation coefficients in both
  analytic suites and the exact-current adapter;
- exact droplet and all-constituent interstitial/cloud-borne clipping counts,
  raw minima, weighted deficits, and clipping-ledger changes;
- an IMEX limiter-active direct two-level solve for every field; and
- explicit failure tests for nonfinite input, FE raw updates and clipping
  diagnostics, and IMEX solution, matrix-residual, physical-residual, and
  ledger-normalization checkers.

The last property follows from the production substep bound: the retained
coefficient of the old value is positive, mixing contributions from neighboring
levels are nonnegative, and the limited activation removal cannot exceed its
matching positive interstitial contribution from below. The forced unstable
step remains only to test clipping diagnostics.

An independent extra-high review of candidate commit `2c4427e` rebuilt and
reran both configurations, inspected the equations and durable tests, and used
reviewer-only mode, matrix, clipping, and hand-solve checks. It assigned strict
`PASS` decisions to both the Forward Euler and IMEX gates. The complete review
is retained in `independent_review_20260806.md`.

### Executed Backward Euler gate: 2026-08-24

Task 3 implemented the host-only Backward Euler gate at code/test commits
`2f43c8c1e4` and `c6f7f18a5d` using the current lock manifest. The final
debug build with ASan/UBSan and optimized release build each passed all four
registered CTest executables. Backward Euler passed 23,165 checks in each
configuration; the unchanged FE, IMEX, and current-MAM4xx adapter gates also
passed.

The BE suite checks the signed 2-by-2 block stencil against a separately
assembled dense partial-pivoting oracle, verifies exact matrix-family reuse,
tests unsafe block and scalar pivots, covers small and L128-like columns from
zero through full overlap and typical through strong mixing, and exercises
activation strengths, active-field counts, ledger conservation, raw
positivity, reset/resuspension, scaling, fail-closed arithmetic, and three-run
bitwise state determinism. Analytic diffusion and activation cases both meet
the approved first-order slope and `R^2` gates. A test-only inspection verifies
the retained nonzero bottom activation block and positive bottom droplet source
because the canonical no-flux limiter reduces the admissible bottom rate to
zero.

The detailed equation audit, TDD record, results, hashes, limitations, and
method-local `PASS` decision are in
`backward_euler_gate_review_20260824.md`. Retained verbose logs and provenance
are below
`/Users/burr114/e3sm-results/ndrop-screening/task3-backward-euler` and
`/Users/burr114/e3sm-results/ndrop-screening/provenance`.

### Executed results: 2026-08-24 current locked stack

Task 2 requalified FE and IMEX at code/test commits `e67097b6bc` and
`07c786d82d` using the current lock manifest. Both the debug build with
ASan/UBSan and the optimized release build passed all three registered CTest
executables. FE passed 23,195 checks, IMEX passed 25,138 checks, and the
current-MAM4xx adapter passed all 15
cases with zero field mismatches and exact `nsubmix` equality. The
controller-limiting case returned `canonical_nsubmix=2`, `exact_nsubmix=2`, and
`physical_only_nsubmix=1`. Each FE and IMEX executable also performs three
serial advances and verifies bitwise-identical full output states.

The focused gate added signed total/clipping/operator ledger decompositions,
per-substep FE and independent IMEX physical residuals, explicit
roundoff-versus-material negative classifications, finite diagnostic guards,
zero-ledger absolute checks, and state scaling by `1e-12`, `1`, and `1e12`.
The detailed commands, TDD record, hashes, limitations, and separate method
decisions are in `fe_imex_gate_review_20260824.md`. Retained verbose logs and
locked-stack provenance are below
`/Users/burr114/e3sm-results/ndrop-screening/task2-fe-imex` and
`/Users/burr114/e3sm-results/ndrop-screening/provenance`.

### Executed results: 2026-08-06

The sanitizer/debug and optimized Forward Euler runs each passed 21,284 checks:

```text
CONVERGENCE case=fe_diffusion order_range=1.00004:1.00004 min_r2=1
CONVERGENCE case=fe_activation_all_fields order_range=1.00514:1.00704 min_r2=0.999998
RESULT passed=21284 failed=0
```

The sanitizer/debug and optimized IMEX runs each passed 21,622 checks:

```text
CONVERGENCE case=diffusion l1_order_range=0.999957:0.999957 linf_order_range=0.999957:0.999958 min_r2=1
CONVERGENCE case=activation l1_order_range=0.991946:0.992783 linf_order_range=0.993205:0.993723 min_r2=0.999998
CONVERGENCE case=coupled_all_fields order_range=0.99701:1.00319 max_finest_cross_error=9.0418e-05
RESULT passed=21622 failed=0
```

The exact-current adapter passed ten sanitizer/debug and optimized cases with
zero mismatches and zero maximum absolute or relative difference. This includes
the distinct-mode case and isolated, contiguous, and all-zero cloud profiles.

Both CTest configurations passed all three registered executables. All build
products were under `/private/tmp/ndrop-gate-remediate`; no dependency or
production-physics source was modified.

The historical 2026-08-06 debug executables used ASan/UBSan and Kokkos bounds
checking. The historical optimized executables used `-O3 -DNDEBUG` and the
release Kokkos/EKAT libraries. Both configurations used double precision,
Apple clang 21.0.0 on arm64 macOS, Kokkos 4.5.1 Serial, and EKAT's serial
Thomas solver. Compilation with
`-Wall -Wextra -Wpedantic -Werror` was clean for harness sources; EKAT/Kokkos
headers were treated as system headers for this check.

### Provenance

- E3SM remediation base revision:
  `3f8d0812602730b7bcca2b02e5a4ecae42fd62d8`
- MAM4xx inspected revision:
  `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`
- EKAT inspected revision:
  `65cda3ff61abe22add622b56cee6234e29fc51c4`
- Kokkos inspected revision: `6a19eb9454da490dbf785443d8f021028baeb161`
- Haero inspected revision: `8f9b61aed4252461d76eb7c206927e4a005930f5`

Retained output records:

| Record | SHA-256 |
|---|---|
| `fe_validation_20260806.txt` | `b51dbadd2575affd80e8d81d069bf6d690b030d2989b9bd27b05fa30e2b26490` |
| `imex_validation_20260806.txt` | `868a6ce12850251c0274a7cb85abf00be726a542f44ecdb58232d778f152d695` |
| `mam4xx_comparison_20260806.txt` | `95aa39bbeacfde590ec6e3ec2a3ad1520cd04d4b5eeb300c8c52529d2d1ac83b` |
| `independent_review_20260806.md` | `dda1d5f9e0661898327469c1ab94a1e8cef60d6c6b15dbc4f36e8efb3e6d0b26` |

This is local CPU correctness evidence for the test-only operator. It is not
GPU performance evidence, whole-model validation, or authorization for a
production IMEX implementation. The known `top_lev=0` and
`enable_aero_vertical_mix=false` defect paths remain excluded.
