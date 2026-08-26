# Task 5R Evidence-Remediation Code Review

## Decision

**Implementation review: PASS. Retained evidence: PENDING.**

The review covered the approved amendment at `61ae8249af`, paired benchmark
and analyzer implementation at `7d1252e2e9`, transaction tooling at
`4902313371`, and review fixes at `d55ef90e46`. It found and resolved four
acceptance-relevant evidence defects. The second focused pass found no
remaining Critical or Important implementation issue.

This decision does not accept a performance result. No fully paired 135-shard
run, final manifest, Task 6 calculation, break-even estimate, hybrid result, or
method recommendation existed during the review.

External research: not used — the approved design, locked repository code,
schemas, and executable behavior fully determine this review.

## Scope and invariants

- Production EAMxx/MAM4xx physics, numerical methods, solver arithmetic,
  boundaries, limiters, interfaces, and gitlinks were not changed.
- The varying axis is nonzero-state occupancy. Every method executes 25
  aerosol pairs and 51 transport fields; width-scaling claims are prohibited.
- IMEX exposes one supported combined coefficient-copy/factorization/solve
  diagnostic phase. Separate IMEX factorization and solve values are blank.
- Phase and fence observations are diagnostic only. Comparative functions use
  only primary uninstrumented complete-advance timings.
- A retained shard does not exist unless its paired artifacts, hashes, and
  completion record all validate.

## Findings and fixes

### 1. Even-sample medians could differ by one ULP

Python used the library median expression `(lower + upper) / 2`, while the C++
benchmark deliberately uses the overflow-safe expression
`lower + (upper - lower) / 2`. A hand-derived subnormal fixture demonstrated a
one-ULP difference and failed before the fix.

The analyzer now uses the exact C++ expression for primary, instrumented,
phase/fence, and bootstrap medians. The fixture passes and paired validation
requires exact reconstructed medians.

### 2. Summary and sidecar aliases could bypass string comparison

The benchmark originally rejected only identical path strings. A red CTest
proved that `aliased_output.csv` and `./aliased_output.csv` could name the same
file and still be accepted.

The parser now compares normalized filesystem paths before opening either
stream. Expected-failure CTests also require the intended diagnostic text, so
a crash or unrelated nonzero exit can no longer count as a pass.

### 3. Final validation could depend on disposable scratch inputs

An identity manifest could bind a benchmark executable and dependency record
that existed only below scratch. Later validation would fail if those inputs
were cleaned, despite the result directory otherwise being complete.

Identity creation now copies the exact executable and dependency-provenance
record exclusively into `<run>/inputs`, sets read-only/executable modes, hashes
the retained copies, and binds every shard command to the retained executable.
Sealing rejects any noncanonical path or unexpected run-directory entry.

### 4. Timing stopping-policy fields were not fully cross-checked

The analyzer accepted valid-looking rows whose `stop_reason` contradicted the
recorded CI. It also did not require a retained `maximum_repetitions` row to
reach the manifest's configured maximum.

Validation now requires `ci_target` only for relative half-width below `0.05`,
requires `maximum_repetitions` only when that target was not met, and, for a
retained run, requires the configured maximum repetition count. Tests cover
wrong relative width, instrumented CI, nonfinite values, overhead, stop reason,
and invalid-row timing fields.

## Executed evidence

All builds were outside the synced source tree under
`/Users/burr114/e3sm-scratch/ndrop-screening/build/`, with build parallelism
capped at four.

| Gate | Result |
|---|---|
| Focused Python suite | 35/35 passed |
| Debug ASan/UBSan registered suite | 13/13 passed in 452.66 s; no sanitizer finding |
| Optimized registered suite | 13/13 passed in 9.70 s |
| Explicit optimized numerical gate | 5/5 passed in 7.63 s |
| Canonical path and diagnostic failure tests | 3/3 focused optimized tests passed |

Commands:

```text
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s docs/analysis/evidence/ndrop_integrator -p 'test_*.py' -v

cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/task5r-final-debug --parallel 4
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/task5r-final-debug --output-on-failure

cmake --build /Users/burr114/e3sm-scratch/ndrop-screening/build/task5r-final-release --parallel 4
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/task5r-final-release --output-on-failure
ctest --test-dir /Users/burr114/e3sm-scratch/ndrop-screening/build/task5r-final-release \
  -R '^ndrop_(forward_euler|imex|backward_euler|cross_scheme|mam4xx_exact_forward_euler)$' \
  --output-on-failure
```

Validation executable hashes:

- optimized: `8bc3c27b03a8a44d962893164ab514d1b7b89092e6a2285988e0f9bee919f6a7`
- debug ASan/UBSan: `1fbea7a2929e5205d49dee3247e83344d9c1c1f715977f78e2ef4f54601573c3`

These hashes identify pre-retained validation executables only. The retained
run must use a newly sealed clean-detached identity, rerun the numerical gate
with its exact optimized binary, and record its own hashes and commands.

## Second focused pass

The second pass checked schema nullability, raw-sample correspondence, exact
median arithmetic, stopping rules, canonical paths, retained input copying,
exclusive writes, source/binary/provenance binding, completion-marker order,
deterministic merge order, and the primary-versus-diagnostic timing type
boundary. No unresolved Critical or Important issue remains in the Task 5R
implementation.

## Remaining gates and limitations

- Final clean-detached dependency and source validation is pending.
- The exact retained timing binary and dependency-provenance manifest do not
  yet exist.
- Storage preflight, all 135 aligned shards, merged artifact hashes, and the
  final manifest are pending.
- Actual matched-error, break-even, FE/IMEX hybrid, performance tables, and
  the single method recommendation remain blocked.
- GPU performance, production-width variability, whole-model behavior, and
  supported-machine EAMxx integration remain out of scope.
