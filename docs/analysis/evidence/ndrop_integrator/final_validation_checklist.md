# `ndrop` Cross-Method Screening Final Validation Checklist

## Purpose and review scope

This checklist is the strict, clean-checkout Task 8 gate for the standalone
Forward Euler (FE), IMEX Euler, and Backward Euler (BE) screening study. It is
controlled by `docs/design/eamxx_mam4xx_ndrop.md` and
`docs/analysis/evidence/ndrop_integrator/2026-08-18-cross-method-screening-plan.md`.
It does not authorize production EAMxx or MAM4xx changes, GPU-performance
claims, or whole-model speedup claims.

The static post-commit review covered base
`7c9e187c5305e082f69c8159b5c2a9aeca2658bb` through Task 5 head
`dba6d00598d86a239a2bf2743d95325821c9378c`. The reviewer inspected committed
source and controlling documents only. No build, executable, test suite, pilot,
partial timing output, active result, or benchmark result was run or read.

External research: not used — the review is fully determined by the committed
repository design, execution plan, source, and tests.

## Static post-commit review findings

### Critical

None found by the code-and-document-only review.

### Important

1. **The transported-field performance axis is nominal, not an executed-width
   axis.** `make_initial_state` makes only `active_pairs` pairs nonzero
   (`benchmark_ndrop_integrators.cpp:207-218`), and the CSV labels that input as
   `1 + 2*active_pairs` transported fields
   (`benchmark_ndrop_integrators.cpp:648-650`). However, FE still constructs
   sources and advances all 25 pairs (`ndrop_integrator.cpp:654-707`), IMEX
   allocates and solves all 51 fields (`ndrop_imex.cpp:337-396`), and BE builds
   families and solves all 25 pairs (`ndrop_backward_euler.cpp:624-670`). The
   reported executed-update count is correspondingly fixed at 51 fields
   (`benchmark_ndrop_integrators.cpp:636-641`). Therefore rows labeled 3, 9,
   17, or 33 transported fields do not measure those execution widths and
   cannot support Task 6 break-even conclusions by transported-variable count.
   Resolve by implementing and validating a genuine variable-width path, or
   reclassify this axis as nonzero-state occupancy and prohibit field-count
   cost conclusions.

2. **IMEX factorization is not separately timed as required.** The IMEX helper
   copies the tridiagonal coefficients and calls the combined EKAT Thomas
   operation inside `solve` (`ndrop_imex.cpp:81-92`). Both IMEX systems put the
   entire call in the `solve` phase (`ndrop_imex.cpp:367-369` and
   `ndrop_imex.cpp:418-420`), so `factorization_seconds_per_column` cannot
   represent IMEX factorization even though Task 5 requires separable
   setup/factorization/solve/update/fence phases. Resolve the phase boundary or
   explicitly document that the backend exposes only a combined
   factorization-plus-solve measurement and do not compare its zero
   factorization field with BE.

3. **Raw phase and fence samples are retained in memory but discarded from the
   machine-readable output.** `RowResult` holds the fence and five phase sample
   vectors (`benchmark_ndrop_integrators.cpp:76-77`), but the CSV writes only
   their medians (`benchmark_ndrop_integrators.cpp:686-691`) and serializes raw
   samples only for the full uninstrumented and instrumented regions
   (`benchmark_ndrop_integrators.cpp:708-709`). A final reviewer cannot
   independently reproduce the phase medians, inspect phase outliers, or audit
   the fence distribution. Serialize those raw samples (preferably as a
   normalized sample table keyed by run/case/repetition/phase) before accepting
   phase-level results.

### Minor

1. **Workspace is an unqualified estimate in the CSV.** The value comes from a
   fixed scalar-count model (`ndrop_benchmark.cpp:244-259`), while the schema
   calls it `workspace_bytes` (`benchmark_ndrop_integrators.cpp:626-627`) and
   the test checks only that it is positive
   (`test_benchmark_support.cpp:226-230`). Rename it
   `estimated_workspace_bytes`, document inclusions/exclusions, and keep it in
   the estimate—not measurement—portion of the final report unless it is
   replaced by an independently verified memory measurement.

2. **The committed matrix test does not independently prove composite-key
   uniqueness or every dimension's exact marginal coverage.** It verifies the
   total/method counts, endpoint rows, contiguous IDs, and nine-row group
   layout (`test_benchmark_support.cpp:119-185`), but a duplicated interior
   dimension value paired with a missing value could pass. The final evidence
   audit below must perform the stronger uniqueness and exact-set checks.

## Pending Task 5R remediation contract

**Status: PENDING EXPLICIT APPROVAL AND IMPLEMENTATION.** The design amendment
records the selected disposition of all three Important findings, but the
findings remain open until the corresponding source, tests, reviews, and
retained rerun pass this checklist. The pre-remediation run bound to
`dba6d00598d86a239a2bf2743d95325821c9378c` is immutable historical evidence
and must never be merged into the remediated dataset.

- [ ] Require summary rows to declare `artifact_schema=ndrop_summary_v2`, a
  valid `run_id`, `nonzero_aerosol_pairs`, `nonzero_state_fields`,
  `executed_aerosol_pairs=25`, and `executed_transport_fields=51`. Remove the
  ambiguous `active_pairs` and `transported_fields` columns, and prohibit all
  executed-width scaling conclusions.
- [ ] Require `linear_phase_mode=none` for FE, `combined` for IMEX, and
  `separate` for BE. For IMEX, time coefficient-copy/allocation plus both
  supported EKAT serial Thomas operations only as combined
  factorization-plus-solve; unavailable separate values are blank, not zero.
- [ ] Require a normalized `ndrop_phase_samples_v1` sidecar with one row for
  every instrumented repetition of every `validation_ok=true` summary row,
  including timed-but-ineligible rows. Invalid rows have `repetitions=0` and
  no sidecar rows.
- [ ] Require a pre-run identity manifest, paired partial summary/sidecar
  outputs, validation and hashes for both files, and a per-shard completion
  marker written last. A shard is nonexistent for analysis unless its marker
  and both hashes validate.
- [ ] Use 180-case shards aligned to complete nine-row method groups, yielding
  exactly 135 shard pairs and completion records. Merge deterministically by
  case ID and instrumented repetition, then bind both merged artifacts in a
  final immutable manifest.
- [ ] Permit speedup, matched-error, break-even, hybrid, and recommendation
  calculations to consume only complete primary uninstrumented timings.
  Instrumented phase/fence observations remain diagnostic-only evidence.

## Gate rule

- [ ] Do not issue a final PASS while any Critical or Important finding above
  is open. Record each resolution as a commit, design/document amendment when
  behavior or interpretation changes, focused regression test, and clean
  debug/optimized rerun.
- [ ] Treat every unchecked item as either a blocking failure or an explicitly
  named limitation. A limitation cannot waive a numerical, provenance,
  completeness, exclusion, uncertainty, or no-extrapolation requirement.
- [ ] Preserve raw and derived evidence. Never repair a retained result in
  place; create a new run with a new immutable manifest and hashes.

## 1. Clean source and locked-dependency provenance

- [ ] Create a short, unsynced, detached clean checkout for final review. Record
  its absolute path. Do not build or run in the OneDrive source tree.
- [ ] Record `git rev-parse HEAD`; require the reviewed candidate commit (or a
  documented descendant containing only reviewed remediations).
- [ ] Require `git status --porcelain=v1 --untracked-files=all` to be empty
  before configure, before each retained run, and after each retained run.
- [ ] Record and hash the complete diff from
  `7c9e187c5305e082f69c8159b5c2a9aeca2658bb`; confirm production EAMxx/MAM4xx
  physics and production gitlinks are unchanged.
- [ ] Hash the authoritative design, execution plan, dependency lock,
  preparation script, CMake file, all harness sources/tests, analysis program,
  analysis tests, report, and this checklist.
- [ ] Run the locked-dependency preparation workflow with `--stack current` in
  verify-only mode immediately before every retained gate or timing run.
- [ ] Require clean detached dependency checkouts and exact SHAs: E3SM baseline
  `e2c578e40f6602a108eff577091de86a453624ab`, MAM4xx
  `0d076d8ebb004a18f05f51f97b08edd9def6826f`, EKAT
  `3a1485aecfed7a036e4fa30ce845498fd37975e2`, and Kokkos
  `23745d8494203711a6090d0a201c3ae3bf57853b`.
- [ ] Verify every locked gitlink, `detached=true`, `dirty=false`, and
  `actual_sha == locked_sha`. Reject the archived stack as current evidence.
- [ ] Archive and hash the provenance JSON. Require it to record source head and
  cleanliness, dependency paths/hashes, compiler command/version, CMake
  generator/build type/flags, CPU/OS/architecture, and exact preparation
  invocation.
- [ ] Create a separate immutable run manifest binding the provenance hash,
  source/binary hashes, exact configure/build/test/benchmark commands,
  environment affecting execution, start/end timestamps with timezone, output
  paths, run kind (`debug`, `optimized-full`, `pilot`, or `smoke`), and hashes
  of every retained output.

## 2. Independent debug/ASan/UBSan and optimized suites

- [ ] Configure fresh build directories under the approved scratch root. Do
  not reuse Task 5 build products.
- [ ] Debug gate: double precision, Debug, AddressSanitizer and UBSan enabled,
  Kokkos bounds checking where supported, warnings-as-errors, exact adapter ON.
- [ ] Optimized gate: fresh release Kokkos build, Release with the recorded
  optimization and `NDEBUG` flags, sanitizers OFF, exact adapter ON.
- [ ] Record CMake cache and configure/build logs; verify the resolved EKAT,
  MAM4xx, Kokkos, compiler, and flags match the run manifest.
- [ ] In both configurations, run the complete registered suite, not selected
  tests: FE, IMEX, BE, cross-scheme, benchmark support, benchmark smoke,
  benchmark exclusion, and the exact-current adapter.
- [ ] Require zero sanitizer findings, zero failed assertions/tests, zero
  unexpected skips, and deterministic nonzero test counts consistent with the
  reviewed executable set. Preserve full logs and hashes.
- [ ] Rerun the complete numerical acceptance suite with the exact optimized
  executable/binary set used for timing. Any optimization after this rerun
  invalidates the gate and requires another complete rerun.

## 3. Exact-current MAM4xx adapter

- [ ] Confirm `NDROP_BUILD_EXACT_ADAPTER=ON` and that
  `ndrop_mam4xx_exact_forward_euler` is registered and executed in both suites.
- [ ] Confirm the adapter uses the locked current rank-3 MAM4xx API, not the
  archived Haero stack or historical expected-failure path.
- [ ] Require direct equality of returned `nsubmix`, final alternating buffer,
  droplets, and every aerosol reservoir for `top_lev=0,1,6`.
- [ ] Require the controller-limiting profile to return two substeps in current
  MAM4xx and the canonical FE path while the same physical no-flux operator
  without the controller returns one.
- [ ] Require production substeps 1, 2, 3, and 5; full/fractional overlap;
  limiter active/inactive; isolated, contiguous, and all-zero cloud patterns;
  distinct number/mass mode mappings; and poisoned destination buffers.
- [ ] Keep legacy `1b02f272...` `top_lev=0` evidence explicitly historical and
  excluded from the current PASS count.

## 4. Task 5 benchmark implementation gate

- [ ] Verify diagnostics-on validation and diagnostics-off timing paths are
  bitwise state-identical with equal clipping/material-negative counts for each
  method. A passing detailed path must never map to a failed lean path.
- [ ] Verify the uninstrumented full region includes coefficient/source and
  matrix/RHS assembly, controller/limiter logic, required copies,
  factorizations, all solves, positivity handling, FE clipping, state update,
  layout/staging costs, and the final Kokkos fence.
- [ ] Verify identical initial states, active levels, coefficients, frozen ACI
  intervals, and production-feasible layouts across methods. Confirm no matrix
  or factorization is reused across columns.
- [ ] Resolve Important finding 1 by relabeling the varying input as
  `nonzero_aerosol_pairs`, deriving `nonzero_state_fields` as
  `1 + 2*nonzero_aerosol_pairs`, and proving that every method still executes
  exactly 25 aerosol pairs and 51 transport fields. Prohibit width-scaling
  conclusions.
- [ ] Resolve Important finding 2 with the exact phase modes and nullability
  contract above. Prove IMEX's combined phase encloses matrix-view
  allocation/copy and both public EKAT Thomas calls without changing solver
  arithmetic; retain BE's existing separate boundaries.
- [ ] Verify primary time is seconds per column-call and is the only timing used
  for CPU speedup. Instrumented and phase times may diagnose overhead only.
- [ ] Verify operation counters against hand-derived counts for representative
  FE, IMEX, and BE cases, including matrix families, factorizations, RHS solves,
  substeps, levels, fields, and batches.
- [ ] Verify workspace accounting against documented live allocations and peak
  concurrency. Label modeled values as estimates and do not present them as
  measured bytes.

## 5. Complete 24,300-row screen and schema

- [ ] Accept only a completed full run or an exact, manifest-verified partition
  of full-run chunks. Require 24,300 data rows plus one header.
- [ ] Require unique contiguous `case_id` values `0..24299` and a unique
  composite key `(nlev, nonzero_aerosol_pairs, batch_size, mixing_index,
  mixing_scale, dt, method, forced_fe_substeps)`.
- [ ] Require exact method counts: 18,900 FE, 2,700 IMEX, and 2,700 BE.
- [ ] Require exactly 2,700 physical groups. Every group must contain FE
  `nsubmix={1,2,3,5,8,16,21}` exactly once, IMEX exactly once, and BE exactly
  once.
- [ ] Require exact dimension sets with no missing, duplicate, or extra values:
  `nlev={1,2,3,16,32,64,72,127,128}`;
  `nonzero_aerosol_pairs={1,4,8,16,25}`;
  `batch_size={1,16,64,256,1024}`; mixing `(index,scale)={(0,0.05),(1,0.5)}`;
  `dt={33.33,50,75,100,150,300}`; FE substeps as above and zero for IMEX/BE.
- [ ] Verify marginal counts independently: 2,700 rows per level; 4,860 per
  pair count; 4,860 per batch; 12,150 per mixing regime; 4,050 per ACI interval;
  and 2,700 per FE substep.
- [ ] Require `nonzero_state_fields == 1 + 2*nonzero_aerosol_pairs`,
  `executed_aerosol_pairs == 25`, and `executed_transport_fields == 51` on
  every row.
- [ ] Require `artifact_schema=ndrop_summary_v2` and a `run_id` matching
  `[A-Za-z0-9._-]{1,128}` on every summary row. Reject `ndrop_summary_v1`,
  unversioned data, mixed run IDs, and missing or extra rows.
- [ ] Enforce linear-phase nullability exactly: FE has mode `none` and all
  three linear columns blank; IMEX has mode `combined`, only the combined
  column numeric; BE has mode `separate`, factorization and solve numeric and
  combined blank.
- [ ] Freeze both exact schema versions. Reject missing/extra columns, duplicate
  column names, malformed CSV quoting, inconsistent row widths, noncanonical
  booleans/statuses, nonfinite numeric text, or semicolon-list parse failures.
- [ ] Require all row metadata to agree with the one run manifest and optimized
  binary: compiler, platform, Kokkos execution space, hardware threads,
  warmups, repetition limits, bootstrap count, and source/binary identity.
- [ ] Recompute and record byte counts, row counts, SHA-256 values, and
  canonical content hashes for the merged summary and phase sidecar after
  parse-and-reserialize validation.

## 6. Warmup, repetition, bootstrap, and interleaving policy

- [ ] For the retained full screen require exactly 10 untimed warmups per row,
  minimum 30 and maximum 100 measured repetitions, 10,000 bootstrap resamples,
  and the approved fixed seed with the documented case-ID derivation.
- [ ] For each valid timed row require `30 <= repetitions <= 100` and exactly
  that many primary raw samples. Require finite nonnegative samples.
- [ ] Independently recompute the sample median and deterministic percentile
  bootstrap 95% interval. Match the recorded median, lower/upper limits, and
  relative half-width within a declared roundoff tolerance.
- [ ] Accept `stop_reason=ci_target` only when repetitions are at least 30 and
  relative half-width is strictly below 0.05. Accept
  `stop_reason=maximum_repetitions` only at 100 samples. Invalid rows must have
  zero samples and `stop_reason=invalid`.
- [ ] Require the instrumented sample count to equal the primary sample count,
  recompute its median/CI, and recompute instrumentation overhead as
  `instrumented_median/primary_median - 1`.
- [ ] Verify temporal method interleaving from the exact executable and run log:
  cyclically rotated order within each nine-row group for warmups and measured
  repetitions. Reject a serial method-at-a-time replacement run.
- [ ] Verify output checksums are consumed outside the timed region and retained
  per row. Confirm diagnostics-on/off state neutrality and deterministic status;
  investigate any checksum drift before using timings.
- [ ] Quantify bootstrap uncertainty for derived speedups, matched-error
  comparisons, and break-even values; do not treat row medians as exact.

## 7. Raw-sample, status, phase, and metadata integrity

- [ ] Preserve the original CSV read-only and verify its hash before every
  analysis. Analysis must fail closed if the hash differs from the manifest.
- [ ] Require raw primary and instrumented sample-list lengths to match
  `repetitions`; no blank, nonfinite, negative, silently truncated, reordered,
  rounded, or deduplicated values are allowed.
- [ ] Resolve Important finding 3 with `ndrop_phase_samples_v1`. Require unique,
  contiguous `(run_id, case_id, instrumented_repetition)` keys and raw
  instrumented total, setup, factorization, solve, combined
  factorization-plus-solve, update/positivity, diagnostics, and fence seconds
  per column for every required repetition. Require the exact timing fields
  `instrumented_seconds_per_column`, `setup_seconds_per_column`,
  `factorization_seconds_per_column`, `solve_seconds_per_column`,
  `combined_factorization_solve_seconds_per_column`,
  `update_seconds_per_column`, `diagnostics_seconds_per_column`, and
  `fence_seconds_per_column`.
- [ ] Require sidecar coverage for every `validation_ok=true` summary row,
  including timed-but-ineligible rows, and none for invalid rows. Require each
  sidecar instrumented total to equal the corresponding summary raw
  instrumented sample element by element.
- [ ] Independently reconstruct every summary phase/fence median using the
  exact C++ median algorithm and require exact method-specific nullability.
- [ ] Require full time, phase times, and fence time to be finite and
  nonnegative. Investigate impossible ordering or phase sums; account for timer
  overhead and any intentionally unclassified interval.
- [ ] Validate `validation_ok`, `eligible`, `exclusion_reason`, `message`,
  `richardson_ok`, `richardson_message`, `stop_reason`, clipping count,
  material-negative count, factorization/solve counts, update count, workspace,
  checksums, compiler/platform/Kokkos, and all timing-policy fields.
- [ ] Do not infer a successful measurement from process exit alone. Reconcile
  the executable's invalid-row total with parsed row statuses and inspect every
  distinct failure/exclusion message.
- [ ] Bind both outputs to source/dependencies/build with the identity manifest,
  per-shard completion records, and final manifest; compiler and platform
  strings embedded in a CSV are not sufficient provenance.

## 8. Richardson errors and exclusions

- [ ] Require exactly 51 weighted-L1 and 51 Linf values whenever
  `richardson_ok=true`, ordered as droplet number, then interstitial/cloudborne
  for constituent pairs 0 through 24.
- [ ] Independently recompute each field's normalized weighted-L1 and Linf
  difference between the coarse result and the two-half-step/refined result,
  with denominators based on the refined and initial field separately. Confirm
  the first-order denominator `2^1-1` and fail on dimension/nonfinite errors.
- [ ] Confirm `richardson_maximum_weighted_l1` and
  `richardson_maximum_linf` are exact maxima of their arrays and
  `screening_error` is the larger maximum. Verify the refined-state checksum.
- [ ] A matched-error point must satisfy all of: method passed Task 4;
  `validation_ok=true`; coarse status successful; `eligible=true`;
  `richardson_ok=true`; FE coarse and refined clipping counts zero; implicit
  coarse and refined material-negative counts zero; all error/timing values
  finite; and no unresolved failure message.
- [ ] Report excluded FE raw/clipped difference, clipping amount/count, and
  ledger impact. Report excluded IMEX/BE raw-negative magnitude/count and
  failure status. Never silently drop excluded rows.
- [ ] Do not let the CSV's `eligible=true` alone override
  `richardson_ok=false`; analysis must encode and test the full predicate above.

## 9. No pilot, smoke, partial, stale, or mixed evidence

- [ ] Store smoke, exclusion-test, pilot, debug, optimized-full, and analysis
  outputs in separate immutable directories with distinct manifests.
- [ ] Reject smoke defaults (1 warmup, 3 repetitions, 100 resamples), exclusion
  test defaults/overrides, or any other reduced policy from the full dataset.
- [ ] Do not append a new run to an old CSV or overwrite an existing retained
  file. Reject duplicate headers and embedded partial-run boundaries.
- [ ] Require exactly 135 aligned 180-case shards. For each shard, require
  identical source/dependency/binary hashes and timing-policy/environment
  metadata; exact disjoint case ranges; a validated summary/sidecar pair; and
  a completion record created only after both final files and hashes exist.
- [ ] Require complete `0..24299` coverage, no overlap, one deterministic
  merge, and hashes for every original summary, sidecar, completion record,
  merged artifact, and manifest.
- [ ] Reject chunks from different builds, machines, execution spaces,
  compiler flags, thread settings, dates after code changes, pilot parameters,
  or source cleanliness states.
- [ ] Reject a truncated file even if its last written group is internally
  complete. Completion is established only by the full schema, count,
  uniqueness, coverage, and manifest gates above.

## 10. Task 6 analysis and synthetic-test gate

- [ ] Review `analyze_screening.py` and its tests without trusting Task 5's
  derived columns. Require fail-closed schema, provenance, hash, completeness,
  status, raw-sample, and uniqueness validation before analysis.
- [ ] Require hand-derived synthetic fixtures for odd/even medians,
  deterministic bootstrap CI bounds, CI propagation, clipping exclusions,
  implicit-negativity exclusions, failed/missing Richardson rows, adjacent
  interpolation, exact measured crossing, no crossing, and break-even answers.
- [ ] Require tests that fail on duplicate/missing case keys, mixed run
  manifests, wrong raw-sample lengths, nonfinite values, pilot parameters,
  partial matrices, and stale hashes.
- [ ] Analyze exact error targets `{1e-2,3e-3,1e-3,3e-4,1e-4}` only where
  bracketed by eligible measured points. Never extrapolate in error, FE
  `nsubmix`, level, nonzero-state occupancy, batch, mixing regime, ACI
  interval, or hardware.
- [ ] Define break-even as the smallest measured FE `nsubmix` whose median
  runtime is at least the candidate median. Interpolate only between adjacent
  measured FE values from `{1,2,3,5,8,16,21}` and only within a valid bracket.
  Return `no measured break-even` outside the bracket—never a capped or
  extrapolated number.
- [ ] Carry raw-sample/bootstrap uncertainty through runtime ratios,
  runtime-versus-error selection, interpolation, and break-even. Report an
  interval or indeterminate classification when uncertainty changes the answer.
- [ ] Never analyze cost by transported-variable count from this screen. The
  nonzero-state-occupancy axis is not a variable-width execution, even after
  its labels and fixed executed-width metadata are corrected.
- [ ] Keep phase diagnostics in a representation that cannot be converted to
  the primary timing estimate. Comparative functions must reject phase data,
  old schemas, incomplete artifact pairs, and missing/stale manifests.
- [ ] Treat the FE/IMEX hybrid as a threshold and packing/workload-variation
  estimate only. Test packing overhead and workload-sensitivity arithmetic with
  hand-derived fixtures; do not call it an implemented or measured hybrid.

## 11. Report: measurements versus estimates

- [ ] Give every table/figure/claim an evidence class: executed Mac CPU
  measurement, directly derived statistic with uncertainty, source-based
  operation/workspace estimate, hybrid packing estimate, or qualitative GPU
  structure consideration.
- [ ] Use only complete uninstrumented advance time for measured CPU speedup.
  Clearly separate instrumented/subtimer diagnostics from primary timing.
- [ ] Label workspace values as estimated until independently measured and
  audited. State the allocation model and exclusions.
- [ ] Keep operation counts separate from elapsed-time measurements and state
  whether each count is instrumented, hand-derived, or source-model based.
- [ ] Make no GPU speedup, occupancy, warp, bandwidth, synchronization, or
  whole-model performance claim from local Mac timing.
- [ ] Report exact provenance and output hashes; accuracy/convergence,
  positivity, clipping, ledgers, typical/strong timing, runtime-versus-error,
  break-even, batching, uncertainty, limitations, failed/excluded cases, and
  the production measurements needed to revise offline estimates.
- [ ] Choose exactly one recommendation: current FE, IMEX, BE, FE/IMEX hybrid,
  or insufficient evidence. Tie it only to passed gates and eligible evidence.

## 12. Separate final decisions and sign-off

- [ ] **FE decision — PASS/FAIL:** state exact adapter result, debug and
  optimized suite results, current-production semantics, clipping/exclusion
  behavior, determinism, source/dependency hashes, and unresolved limitations.
- [ ] **IMEX decision — PASS/FAIL:** state numerical gate, debug and optimized
  suite results, positivity/exclusions, matrix-family behavior, combined or
  separated factorization/solve timing limitation, hashes, and unresolved
  limitations.
- [ ] **BE decision — PASS/FAIL:** state numerical gate, debug and optimized
  suite results, block/droplet factorization and solve behavior, raw positivity,
  bottom-branch exclusions, hashes, and unresolved limitations.
- [ ] **Cross-scheme decision — PASS/FAIL:** state Task 4 common-solution and
  Richardson results, admissible method set, complete Task 5 matrix result,
  uncertainty and exclusion correctness, no-extrapolation compliance, and
  whether comparisons are scientifically supportable.
- [ ] A method FAIL must remain documented and be excluded from performance and
  recommendation comparisons. Cross-scheme PASS is impossible if a failed
  method remains in a claimed comparison.
- [ ] Record reviewer, date/timezone, reviewed commit, base commit, clean source
  and dependency evidence, exact executed commands, test/result counts, hashes,
  all finding dispositions, remaining risks, and untested regimes.
- [ ] Obtain a second brief review pass after all fixes and reruns. Confirm no
  Critical/Important finding remains and no fix introduced a new design or
  evidence mismatch.
- [ ] Only then sign the independent final gate and accept—or reject—the final
  comparison and single recommendation.
