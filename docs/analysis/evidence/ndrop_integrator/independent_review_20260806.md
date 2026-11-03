# Independent Review of the Corrected `ndrop` Forward Euler and IMEX Gates

## 1. Overview and exact provenance/configuration

### Review scope and independence

I reviewed the standalone EAMxx-MAM4xx `ndrop` Forward Euler (FE) and IMEX
gates directly from the specified clean checkout. I did not use prior gate
summaries as evidence for the decisions. I inspected the approved design,
standalone implementation, tests, exact-current adapter, CMake/CTest driver,
retained output records, and the relevant MAM4xx and EKAT/Kokkos source. I then
performed fresh builds and runs and constructed reviewer-only mutation and
hand-calculation tests under `/private/tmp`. The source checkout was not
modified.

The gate scope is deliberately limited to a frozen-coefficient, standalone
column operator. It does not validate the full activation parameterization,
outer MAC/ACI process splitting, GPU performance, production integration, or
the excluded `top_lev=0` and `enable_aero_vertical_mix=false` defect paths.

### Source provenance

| Item | Verified revision/state |
|---|---|
| Review checkout | `/private/tmp/ndrop-rvw-3f8d081` |
| Candidate | `2c4427eb4153788d27714a38229dfd7a209f568c` |
| Candidate tree | `b55db4841aa3762a3568067a8599084fd5094e72` |
| Base before remediation | `3f8d0812602730b7bcca2b02e5a4ecae42fd62d8` |
| MAM4xx | `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db` |
| EKAT | `65cda3ff61abe22add622b56cee6234e29fc51c4` |
| EKAT Kokkos gitlink | `6a19eb9454da490dbf785443d8f021028baeb161` |
| Haero gitlink/source | `8f9b61aed4252461d76eb7c206927e4a005930f5`, clean |
| Checkout status | Clean before and after review |

The reused Kokkos builds were independently checked. Both were built from the
Kokkos revision recorded by the inspected EKAT gitlink and report Kokkos
4.5.1 with the Serial backend. Their generated include paths refer to the
clean cloud-synced source checkout at the same Kokkos revision; all new build
and test output from this review remained under `/private/tmp`.

### Build configurations

| Configuration | Compiler/platform | Relevant flags and checks |
|---|---|---|
| Debug/sanitized | AppleClang 21.0.0, arm64 macOS | `-O0 -g`, ASan, UBSan, frame pointers, Kokkos debug bounds and DualView checks |
| Optimized | AppleClang 21.0.0, arm64 macOS | `-O3 -DNDEBUG`, release Kokkos |
| Both | C++17, Kokkos Serial | `-Wall -Wextra -Wpedantic -Werror` for harness sources |

The debug executable was verified to link the ASan runtime. CMake 4.4.0 was
used. The checked-in CMake driver configured and built without source edits.
`git diff --check 3f8d081..2c4427e` was clean.

## 2. Forward Euler gate

### Methods assessment (Forward Euler)

The FE implementation is appropriate as the canonical reference for this
study. It derives cell rates from shared physical-interface exchange
coefficients and pressure weights, imposes zero exchange at both boundaries,
forms cloud-overlap factors with the current MAM4xx rules, caps both modal
number and mass activation rates by the downward mixing rate, and uses the
current production `nsubmix` rule. Each substep recomputes state-dependent
activation sources, reads only the saved buffer, writes the alternate buffer,
and applies componentwise production clipping.

The signs and indexing agree with the inspected MAM4xx equations. Cloud-borne
aerosol receives `+S`; interstitial aerosol receives `-S`; droplet number
receives the sum of modal-number activation sources. The non-cloudy portion of
cloud-borne material crossing an interface is returned to interstitial aerosol.
The shared-interface construction satisfies
`w[k]*rate_below[k] = w[k+1]*rate_above[k+1]`, so the weighted paired aerosol
ledger cancels interface-by-interface before clipping. The current C++ bottom
source formula, including use of modal-number activation for aerosol mass, is
retained rather than silently corrected.

The exact-current adapter uses `top_lev=1`, initializes the explicit upper
boundary coefficient to zero, enables aerosol vertical mixing, poisons the
destination work buffer, and compares all active droplet and aerosol values.
Its synthetic `mam_idx` is equivalent to the production mapping in
`mam4::ndrop::get_e3sm_parameters`: mode ranges 0-7, 8-12, 13-20, and 21-24,
with modal-number indices 0, 8, 13, and 21.

### Results and realism assessment (Forward Euler)

Fresh debug/sanitized and optimized FE runs each reported:

```text
CONVERGENCE case=fe_diffusion order_range=1.00004:1.00004 min_r2=1
CONVERGENCE case=fe_activation_all_fields order_range=1.00514:1.00704 min_r2=0.999998
RESULT passed=21284 failed=0
```

These slopes are realistic for Forward Euler and are in the approved
first-order interval. They cover separate weighted-L1 and L-infinity fits for
droplet number and all 50 aerosol fields, not a pooled dimensional norm.

The exact-current adapter passed ten cases in both configurations with zero
mismatches and zero maximum absolute or relative difference. It covers
production-adaptive `nsubmix` values 1, 2, 3, and 5; odd/even buffer parity;
full and fractional overlap; limiter-active/inactive states; isolated,
contiguous, and all-zero cloud profiles; and a distinct-mode activation case.

The randomized production-adaptive tests span level counts 1, 2, 3, 127, and
128 and weak-to-strong exchange rates. No clipping was observed. This is
consistent with the discrete positivity argument: the selected substep leaves
a positive old-state coefficient, mixing inflows are nonnegative, and the
capped activation removal does not exceed its associated interstitial inflow.
The forced unstable step is therefore appropriately retained as a diagnostic
test rather than represented as a production-adaptive physical case.

### Validation evidence and harness/tests assessment (Forward Euler)

The checked-in suite provides durable evidence for:

- input dimensions, finite values, nonnegativity, boundaries, and timestep;
- zero-operator identity and full-overlap uniform-state preservation;
- compatible fractional-overlap ledger conservation;
- production substep selection, limiter activation, and buffer parity;
- exact zero-cloud droplet reset and aerosol resuspension;
- first-order analytic no-flux diffusion and activation for all 51 fields;
- exact per-field clipping counts, raw minima, weighted deficits, clipping-only
  ledger increments, and total ledger effects for droplet number and all 25
  interstitial/cloud-borne pairs, with separate zero-clipping controls; and
- deterministic failure at derived-rate, substep-count, activation-source,
  raw-update, clipping-deficit, clipping-ledger, initial-ledger, and zero-cloud
  overflow stages.

Reviewer-only checks strengthened this assessment:

1. Deliberately assigning one accumulation-mode mass field to the Aitken mode
   caused the FE analytic suite to fail five checks, reducing the fitted order
   for affected fields to approximately zero. The exact-current distinct-mode
   case independently reported 254 mismatches.
2. Replacing pressure-weighted clipping deficit with an unweighted value did
   not pass silently; the committed FE suite failed its clipping finite-stage
   assertions.
3. A non-unit-weight hand check independently reproduced one droplet clipping
   event and one interstitial event for every aerosol constituent, raw minimum
   `-9.8`, weighted deficit `19.6`, and aerosol ledger increase `19.6`.

The CMake/CTest harness is concise, reproducible, warning-clean, and isolates
all build products. Tests use deterministic states and a fixed random seed.

### Identified issues or risks (Forward Euler)

No scoped FE implementation defect or gate-blocking validation-artifact
deficiency was found. Residual limitations are:

- exact-current agreement is intentionally restricted to the safe
  `top_lev=1`, enabled-mixing path;
- only the Kokkos Serial CPU backend and double precision were executed;
- synthetic frozen coefficients are not production-column statistics; and
- the retained current-C++ bottom formula is an explicit scientific TODO, not
  evidence that the formula is physically intended.

### Recommendation: PASS (Forward Euler)

**Forward Euler gate: PASS.** The implementation reproduces the scoped current
MAM4xx state transition exactly under safe boundaries, the convergence and
ledger results match theory, clipping and finite-failure behavior are exposed,
and mutation tests demonstrate that the critical mapping and diagnostic
oracles are not merely restating the implementation.

## 3. IMEX gate

### Methods assessment (IMEX)

The IMEX split matches the approved design. Activation is evaluated explicitly
from the old state. Cloud-overlap-weighted mixing is implicit for cloud-borne
aerosol and droplet number. Full interstitial mixing and cloud-to-interstitial
return are implicit, with the return evaluated from the newly solved
cloud-borne state. This yields a block-lower-triangular aerosol system and the
implemented solve order, cloud-borne first and interstitial second, is correct.

The cloud and interstitial matrices have positive diagonals and nonpositive
off-diagonals. Under the validated nonnegative rate and overlap contracts they
are strictly row-diagonally dominant M-matrices because of the identity term.
EKAT's scalar Thomas solver is therefore appropriate for this CPU correctness
gate. No nonlinear solve is needed because rates, overlap, and geometry are
frozen.

The matrix assembly signs and no-flux endpoints are consistent with the
semi-discrete operator. The independently coded physical-equation residual
uses `rate_above`, `rate_below`, overlap, source, old state, and solved state
rather than the assembled tridiagonal coefficients. Aerosol conservation
follows by adding the cloud-borne and interstitial solved equations: activation
cancels and the weighted cloud loss is exactly the interstitial return.

IMEX does not clip. Its positivity check is performed on the raw solved state
before zero-cloud resuspension. Every negative is counted, while only values
larger than the approved field-relative roundoff tolerance cause failure.

### Results and realism assessment (IMEX)

Fresh debug/sanitized and optimized IMEX runs each reported:

```text
CONVERGENCE case=diffusion l1_order_range=0.999957:0.999957 linf_order_range=0.999957:0.999958 min_r2=1
CONVERGENCE case=activation l1_order_range=0.991946:0.992783 linf_order_range=0.993205:0.993723 min_r2=0.999998
CONVERGENCE case=coupled_all_fields order_range=0.99701:1.00319 max_finest_cross_error=9.0418e-05
RESULT passed=21622 failed=0
```

These values are realistic for first-order IMEX Euler. The diffusion and
activation cases compare all 51 fields with analytic solutions. The coupled
case separately establishes first-order self-convergence for FE and IMEX,
decreasing cross-scheme differences, and the approved finest-step Richardson
bound without treating either candidate as truth.

The limiter-active two-level oracle caps all eight modal activation rates,
checks the exact limiter count, and compares droplet and all 25 aerosol pairs
against direct 2-by-2 solves. Material-negative and tiny-scale material-
negative cases both fail without clipping, as required.

### Validation evidence and harness/tests assessment (IMEX)

The checked-in suite covers identity, uniform-state preservation, fractional
overlap, randomized compatible columns, zero-cloud transitions, matrix and
physical residuals, ledgers, positivity failures, analytic convergence,
coupled cross-scheme convergence, limiter activation, and finite failures.

Finite-stage coverage includes end-to-end rejection of nonfinite input,
initial ledgers, derived operator/source values, matrix assembly, right-hand
sides, and zero-cloud resuspension. Solution, matrix-residual,
physical-residual, and ledger-normalization checkers are also exercised with
focused injected nonfinite values. Checker-level probes are acceptable here:
for a finite right-hand side and the validated finite M-matrix, a nonfinite
Thomas solution is not physically constructible without external corruption,
and the approved design explicitly permits a focused checker test in that
case.

Reviewer-only checks provided independent evidence:

1. The same constituent-mode mutation that failed FE also failed six IMEX
   checks, including analytic convergence and finite-stage expectations.
2. A hand-assembled two-level test with pressure weights 2 and 5, cloud
   fractions 0.25 and 0.75, distinct number/mass activation rates in every
   mode, and fractional cloud return matched the IMEX result for droplet number
   and all 25 aerosol pairs and preserved every weighted aerosol ledger.
3. Removing the cloud overlap factor from one IMEX matrix coefficient caused
   448 failures. The matrix solution remained self-consistent, but the
   independently assembled physical residual, ledgers, fractional/randomized
   cases, and coupled integration rejected it. This demonstrates meaningful
   separation between the matrix and physical residual oracles.

The shared CMake/CTest harness is reproducible and maintainable for this local
screening scope. EKAT and Kokkos headers are external system includes, while
all harness translation units retain strict warnings-as-errors.

### Identified issues or risks (IMEX)

No scoped IMEX implementation defect or gate-blocking validation-artifact
deficiency was found. Residual limitations are:

- this is an isolated frozen-linear-operator result, not validation of
  recomputed nonlinear activation physics;
- CPU Serial results do not establish GPU correctness or performance;
- positivity is demonstrated for tested admissible cases but is not guaranteed
  for all explicit-activation timesteps; fail-loud behavior remains essential;
- zero-cloud projection is tested as a discrete state transition, not as part
  of smooth convergence; and
- production integration, data layout, batching, and synchronization costs are
  outside this gate.

### Recommendation: PASS (IMEX)

**IMEX gate: PASS.** The split, matrix assembly, solve order, physical residual,
ledger, positivity policy, and first-order results are consistent with the
approved equations. The committed suite plus independent direct and mutation
tests provide durable, targeted evidence adequate for this standalone gate.

## 4. Cross-gate findings, limitations, and essential next steps

### Cross-gate findings

- FE and IMEX converge at first order independently and toward a common limit
  in the nonuniform, fractional-overlap, all-field coupled case.
- The source-to-constituent mapping is no longer a shared-oracle blind spot:
  distinct analytic coefficients and exact-current comparison reject a
  deliberate permutation.
- Conservation checks are pressure-weighted and field-specific. No dimensional
  pooling of aerosol mass, aerosol number, and droplet number is used.
- Zero-cloud resuspension/reset is consistently applied after each complete
  advance and preserves aerosol ledgers.
- The nonfinite handling is fail-loud at the relevant reachable stages; focused
  probes cover internal checker branches that valid M-matrix inputs cannot
  naturally reach.
- The tests do not count either known production defect path as successful
  validation.

### Limitations

These PASS decisions authorize progression within the approved standalone
plan only. They do not authorize production physics changes, establish GPU
performance, resolve the bottom-source scientific question, or validate the
known production defect paths. The review also did not assess Backward Euler,
which has not yet been implemented at this snapshot.

### Essential next steps

1. Record the independent PASS decisions in the gate-status documentation.
2. Proceed to the separately gated Backward Euler implementation and validate
   its block system against an independent dense partial-pivoting solve before
   any cross-scheme or performance work.
3. Keep production integration and GPU work blocked until all three standalone
   methods pass from one frozen revision and the separate production design
   gate is approved.

No additional FE or IMEX remediation is required before the Backward Euler
gate. Adding the reviewer-only non-unit clipping and direct fractional-overlap
cases to the committed suite would strengthen long-term regression protection,
but their absence is not a gate blocker because the checked-in mutation-
sensitive tests, source review, and executed independent calculations jointly
establish the required behavior.

## 5. Executed commands/tests and concise decision table

### Primary build and test commands

The following command forms were executed with source
`/private/tmp/ndrop-rvw-3f8d081/docs/analysis/evidence/ndrop_integrator`:

```text
cmake -S <source> -B /private/tmp/ndrop-independent-review-2c4427e/build-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DKokkos_DIR=/private/tmp/ndrop-gate-rereview-run/kokkos-debug \
  -DNDROP_ENABLE_SANITIZERS=ON -DHAERO_SOURCE=<clean-haero-source>
cmake --build /private/tmp/ndrop-independent-review-2c4427e/build-debug --parallel 4
ctest --test-dir /private/tmp/ndrop-independent-review-2c4427e/build-debug --output-on-failure

cmake -S <source> -B /private/tmp/ndrop-independent-review-2c4427e/build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_DIR=/private/tmp/ndrop-gate-rereview-run/kokkos-release \
  -DNDROP_ENABLE_SANITIZERS=OFF -DHAERO_SOURCE=<clean-haero-source>
cmake --build /private/tmp/ndrop-independent-review-2c4427e/build-release --parallel 4
ctest --test-dir /private/tmp/ndrop-independent-review-2c4427e/build-release --output-on-failure
```

Each executable was also run directly in each configuration to capture its
full convergence and case output. Both CTest configurations passed 3/3 tests:

| Executable | Debug/sanitized | Optimized |
|---|---:|---:|
| `test_forward_euler` | 21,284 passed, 0 failed | 21,284 passed, 0 failed |
| `test_imex` | 21,622 passed, 0 failed | 21,622 passed, 0 failed |
| `compare_mam4xx_forward_euler` | 10/10 cases, zero mismatches | 10/10 cases, zero mismatches |

The checked-in retained-record SHA-256 values were independently verified:

| Record | SHA-256 |
|---|---|
| `fe_validation_20260806.txt` | `b51dbadd2575affd80e8d81d069bf6d690b030d2989b9bd27b05fa30e2b26490` |
| `imex_validation_20260806.txt` | `868a6ce12850251c0274a7cb85abf00be726a542f44ecdb58232d778f152d695` |
| `mam4xx_comparison_20260806.txt` | `95aa39bbeacfde590ec6e3ec2a3ad1520cd04d4b5eeb300c8c52529d2d1ac83b` |

### Reviewer-only tests

| Test | Expected/observed result |
|---|---|
| Constituent mode mutation | FE failed 5 checks; IMEX failed 6; exact adapter found 254 mismatches |
| Unweighted clipping-deficit mutation | FE suite failed 2 finite-stage checks |
| Non-unit clipping hand oracle | Passed exact minima, counts, weighted deficits, and ledger effects |
| Fractional-overlap IMEX 2×2 oracle | Passed droplet and all 25 pair solutions and ledgers |
| IMEX cloud-matrix overlap mutation | Rejected with 448 failures, including independent physical residual and ledgers |

### Strict decisions

| Gate | Decision | Basis |
|---|---|---|
| Forward Euler | **PASS** | Exact-current agreement, correct first-order behavior, complete fieldwise ledger/clipping evidence, fail-loud finite handling, and mutation-sensitive independent mappings |
| IMEX | **PASS** | Correct approved split and solve order, analytic and common-limit convergence, independent physical residual and direct solve evidence, ledger/positivity enforcement, and finite-stage coverage |
