# Cross-scheme Task 4 gate review

**Execution date:** 2026-08-24

**Code/test commit:** `7807aedfdecafe5e30b783e98e81ac4b723f56e3`

**Final code/test tree:** `645d2c5193f76cf6637e60ab37ba4bd8c720fef2`

**Scope:** standalone common FE, IMEX Euler, and Backward Euler accuracy and
convergence gate; no production EAMxx or MAM4xx source, interface, equation,
gitlink, or runtime change

## Decision summary

| Gate | Decision | Executed evidence |
|---|---|---|
| Forward Euler common gate | **PASS** | All 1,875 matrix cases pass; no material negative or clipping; first-order all-field slopes `1.0011674:1.0016497`, minimum `R^2=0.9999999`; maximum normalized residual `1.4389395e-16` and ledger error `2.3628458e-15`. |
| IMEX Euler common gate | **PASS** | All 1,875 matrix cases pass; no material negative; first-order all-field slopes `0.99825832:0.99868307`, minimum `R^2=0.99999989`; maximum normalized residual `3.2474441e-16` and ledger error `2.3628644e-15`. |
| Backward Euler common gate | **PASS** | All 1,875 matrix cases pass; no material negative; first-order all-field slopes `0.99835794:0.99883807`, minimum `R^2=0.9999999`; maximum normalized residual `4.0419594e-16` and ledger error `2.4656757e-15`. |
| Independent reference gate | **PASS** | Analytic diffusion and activation checks pass; repeating the adaptive Dormand--Prince 5(4) oracle at `rtol=1e-12` and `1e-13` changes the maximum normalized field error by `2.0331196e-13`. |
| Common-solution and determinism gate | **PASS** | All three pairwise method differences decrease at the finest refinement and satisfy the approved Richardson bound in weighted L1 and Linf for all 51 fields; every matrix method performs three bitwise-identical serial advances. |
| Clean debug/optimized regression | **PASS** | Each cross-scheme executable reports 11,880 checks and zero failures. The debug/ASan/UBSan run and optimized run produce bitwise-identical CSV and Markdown outputs. The other four registered executables pass in debug, and all five pass optimized CTest. |
| Production, GPU, or full-interval accuracy readiness | **NOT CLAIMED** | Evidence is a frozen-operator Kokkos Serial CPU screen. Maximum one-call coarse Linf errors reach `0.1129882`; Task 5 timing and Task 6 matched-error analysis remain required before a method recommendation. |

The PASS decisions authorize only the approved local CPU performance screen.
They do not select a production method or establish that one IMEX or BE step
over a production ACI interval meets an accuracy target.

## Overall assessment

The approved design is ready for Task 5 implementation. The governing
operator, method splits, independent reference, parameter matrix, error norms,
positivity and ledger rules, deterministic output contract, and quantitative
acceptance thresholds are implementable and agree with the executed gate. No
critical or unresolved major design finding remains for Task 4.

## Findings

### Major, resolved: nonfinite reporting errors could be hidden by a maximum

- **Location:** cross-scheme per-field error normalization and aggregation.
- **Problem:** the first implementation used ordinary floating-point
  accumulation and `std::max`; a nonfinite difference could propagate to a
  field error and then be ignored by a later maximum reduction.
- **Why it matters:** this violated the approved fail-closed rule and could
  make a corrupted comparison appear to have zero maximum error.
- **Resolution:** a red regression injected a NaN and produced
  `RESULT passed=635 failed=1`. The norm now rejects nonfinite inputs,
  intermediate sums, differences, scales, and normalized outputs, and the
  aggregation independently rejects nonfinite field errors. The smoke rerun
  passed 636 checks.

### Major, resolved: Forward Euler was exempted from the common material-negative gate

- **Location:** common-matrix method-local acceptance expression.
- **Problem:** the first acceptance expression required zero material raw
  negatives only for IMEX and BE, even though the Task 4 amendment says all
  methods must have no material negative.
- **Why it matters:** production FE clipping could otherwise turn an
  inadmissible raw matrix point into a passing common accuracy point.
- **Resolution:** the common condition now requires zero material negatives
  for every method. FE clipping remains fully reported but cannot repair a
  Task 4 gate point. The complete matrix has zero FE material negatives and
  zero clipping cases.

### Minor, resolved: implicit per-field raw-minimum columns had nonfinite sentinels

- **Location:** CSV diagnostic serialization.
- **Problem:** method-local clipping structures use positive infinity as an
  unused sentinel for IMEX and BE, which initially serialized as `inf`.
- **Why it matters:** a retained machine-readable result should distinguish an
  unused field from a measured nonfinite value.
- **Resolution:** implicit global `raw_minimum` reports the checked solved-state
  minimum, while unused FE-specific per-field raw-minimum cells serialize as
  empty strings. Both final CSV files have 50 columns on every row and contain
  no NaN or infinity token.

## Equation-to-algorithm map

| Approved equation or rule | As-built location and evidence |
|---|---|
| Independent physical right-hand side | `ndrop_reference.cpp` reconstructs pressure-normalized interface rates, overlap, activation limiting, mode/type maps, bottom sources, and all 51 tendencies directly from `OperatorData`; it does not call candidate preparation, matrix, or tendency/source helpers. |
| Adaptive Dormand--Prince 5(4) | Standard seven-stage embedded coefficients, infinity-norm error control, field-scaled absolute tolerances, bounded growth/shrinkage, step-underflow and attempted-step limits, and finite checks at every stage. |
| Exact no-flux diffusion and activation oracles | Separate two-level analytic solutions check all 51 reference fields and four mode-specific droplet sources to `5e-11`. |
| Complete deterministic matrix | Exact Cartesian product of five level counts, five overlaps, five activation strengths, five active-pair counts, and three mixing strengths: 1,875 columns and 5,625 method records. |
| Per-field normalized errors | Weighted L1 and Linf remain separate for droplet number and every interstitial/cloud-borne field; nonfinite arithmetic fails closed. |
| Residual, conservation, and positivity | Every matrix record gates the method-local physical residual, all 25 signed ledger decompositions, and raw material-negative count. FE clipping diagnostics are retained but not accepted as a repair. |
| First-order and common solution | Five timestep refinements, three-finest-point fits for both norms and all 51 fields, plus all three pairwise decreasing-difference and Richardson bounds. |
| Three-run determinism | Full state and all diagnostic scalars and arrays compare bitwise across three serial advances per matrix method. |
| Common zero-cloud transition | The reference applies the canonical operation once after continuous integration; matrix zero-cloud cases compare the same one-call state transition, while convergence uses strictly positive clouds. |

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
and Kokkos Serial on an arm64 Apple M2 Max running macOS 26.6.1. The debug
harness enabled ASan and UBSan with a Kokkos Debug build; the separate Kokkos
debug-bounds option was not enabled. The optimized harness and Kokkos
dependency used Release `-O3 -DNDEBUG`. Every harness target compiled with
`-Wall -Wextra -Wpedantic -Werror`.

## TDD and review record

The full cross-scheme test and public test-only reference interface were first
compiled without an implementation. Linking exited nonzero only for the
intentionally missing `advance_reference_dopri54` symbol. The initial oracle
and matrix then passed a focused smoke run and the complete provisional matrix.

The formal design and source review rechecked oracle independence, no-flux
rows, current-C++ source conventions, all-field norm units, two-tolerance
reference consistency, positivity before projection, ledger decompositions,
determinism, and output schema. It exposed the two resolved major findings
above. A second brief pass rechecked every changed acceptance expression and
finite reduction, and both warnings-as-errors builds remained clean.

The final source was committed before verification. A clean detached worktree
at `/Users/burr114/e3sm-scratch/ndrop-screening/worktrees/task4-final` was
created at `7807aedfdecafe5e30b783e98e81ac4b723f56e3`; both provenance records
identify that clean source and the exact locked dependencies.

## Numerical results

Both final configurations produced identical summaries:

| Method | Cases | Failures | Max weighted L1 | Max Linf | Max residual | Max ledger | Material negatives | Clipping cases | Slope range | Min `R^2` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| FE | 1,875 | 0 | `0.091525245` | `0.1129882` | `1.4389395e-16` | `2.3628458e-15` | 0 | 0 | `1.0011674:1.0016497` | `0.9999999` |
| IMEX | 1,875 | 0 | `0.06688989` | `0.070265567` | `3.2474441e-16` | `2.3628644e-15` | 0 | 0 | `0.99825832:0.99868307` | `0.99999989` |
| BE | 1,875 | 0 | `0.057285903` | `0.072283807` | `4.0419594e-16` | `2.4656757e-15` | 0 | 0 | `0.99835794:0.99883807` | `0.9999999` |

The reference two-tolerance maximum is `2.0331196e-13`. Each CSV has 5,629
rows: one header, 5,625 method-case rows, and three convergence rows. Every
row has 50 columns; each method has exactly 1,875 passing case rows; no row
contains a NaN or infinity token. Debug and optimized CSVs are bitwise
identical, as are their Markdown tables and one-line result logs.

The coarse one-call errors are measurements, not acceptance targets. Their
magnitudes make matched-error screening mandatory: the Task 4 PASS means the
methods converge correctly to the same reference, not that these coarse
points are sufficiently accurate for production.

## Retained records

| Record | SHA-256 |
|---|---|
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme-debug.csv` | `b28f376cb377e25d5ea04156063a04b48553015ed995813d9b8c8db15cd98f44` |
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme-debug.md` | `14b7e811ac4dc2b0bd56db24efcc0d259f50c4519424b3fd5195021140f91163` |
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme-debug.log` | `c6e7e1b57f98ed2fca4d5a684e37c6a58ae58bf23aa9971565b75775e5006aaf` |
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme-release.csv` | `b28f376cb377e25d5ea04156063a04b48553015ed995813d9b8c8db15cd98f44` |
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme-release.md` | `14b7e811ac4dc2b0bd56db24efcc0d259f50c4519424b3fd5195021140f91163` |
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/cross-scheme-release.log` | `c6e7e1b57f98ed2fca4d5a684e37c6a58ae58bf23aa9971565b75775e5006aaf` |
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/ctest-debug-aux.log` | `73a046b438768da95a367c43b94c3ac5f805adceb8831e243cff41a4ecd7361c` |
| `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme/ctest-release.log` | `c1c1be482d54d1cec6992153d7dfaff184fe7f6fbfd484b1e2e5693cfe13a7d4` |
| `/Users/burr114/e3sm-results/ndrop-screening/provenance/task4-cross-scheme-debug.json` | `7906f643ec7ee76d8616a599233e589e85eeddeccb771ad8e0794eae8575f8c3` |
| `/Users/burr114/e3sm-results/ndrop-screening/provenance/task4-cross-scheme-release.json` | `d0cd64373c0e58eecd296b80ff1a06de90d2f3a59de17b6e5a6424d1dae65bc1` |

## Open questions

No open question blocks Task 5. The performance screen and Task 6 analysis
must determine which measured timestep reaches each requested error target;
extrapolation remains prohibited.

## Recommended revisions

Proceed to the approved optimized Mac performance screen using only FE, IMEX,
and BE, because all three passed the common gate. Preserve per-case positivity
and clipping exclusions, time the complete uninstrumented advance, and keep
equal-interval runtime separate from matched-error conclusions.

## Residual risks

- The deterministic matrix is synthetic and freezes coefficients over an ACI
  call. It does not sample whole-model coefficient or `nsubmix` distributions.
- Evidence is local Kokkos Serial CPU evidence, not supported-machine EAMxx,
  GPU, or whole-model performance evidence.
- The reference shares the validated input structure and common zero-cloud
  state-transition definition with the candidates, while its continuous RHS
  and time integration are independently assembled.
- The retained current-C++ bottom source remains physically unresolved versus
  the Fortran formulation and must be reviewed before production adoption.
- Outer process splitting, evolving activation coefficients, and coupled cloud
  response remain outside the frozen-operator gate.
