# EAMxx-MAM4xx `ndrop` Boundary and Disabled-Mixing Investigation

Date: 2026-07-14

Status: Share-ready technical brief; no production physics was modified.

## Executive Result

| Concern | Classification | Decisive evidence |
|---|---|---|
| Can `top_lev=0` produce an invalid index? | **CONFIRMED DEFECT** | The default L128 production configuration supplies zero, and the exact 128-level C++ function produced a Kokkos bounds violation at index `-1`. |
| Can `enable_aero_vertical_mix=false` apply paired updates inconsistently? | **CONFIRMED DEFECT** | The switch is an accepted EAMxx parameter. In deterministic one-, two-, and three-substep runs, all 25 aerosol pair ledgers failed while enabled controls conserved them; cloud-borne aerosol changed while interstitial aerosol and droplet number did not. |

These classifications follow the requested standard: each state is reachable
through an allowed production-facing interface, and each concern has an
executed bounds or invariant violation. Neither classification rests on source
inspection alone.

The level-zero boundary defect is tracked publicly in
[MAM4xx issue 564](https://github.com/eagles-project/mam4xx/issues/564).

## Inspected Revisions

| Repository | Revision |
|---|---|
| E3SM | `12e44bb469efb0d6c0cb94c7b796a86967883ff5` |
| `externals/mam4xx` | `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db` |
| E3SM MAM4xx gitlink | `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db` |

The MAM4xx submodule was clean during the investigation. Existing unrelated
working-tree changes were preserved.

On 2026-07-15, the current checkout was E3SM revision
`781c3f5c4cc07f5bc8b5a8d748b0aa923f9b8522`, still using MAM4xx revision
`1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`. The relevant EAMxx ACI
interface, wrapper, MAM4xx defaults, and component-configuration files have no
changes between the evidence revision and this checkout. The classifications
below therefore remain applicable to the current source.

Executable evidence used AppleClang on arm64 macOS, double-precision
`haero::Real`, Kokkos 4.5.1 Serial, `KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK`, AddressSanitizer, and
UndefinedBehaviorSanitizer. The evidence driver configured MAM4xx for 128
vertical levels and compiled the inspected `ndrop.hpp` directly.

## Production Call Path

```text
MAMAci::run_impl
  -> MAMAci::call_function_dropmixnuc
    -> mam4::ndrop::dropmixnuc
      -> mam4::ndrop::update_from_explmix
        -> mam4::ndrop::explmix
```

Relevant source anchors:

| Behavior | Source |
|---|---|
| EAMxx parameter presence checks | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:41-53` |
| Parameter values loaded without range checks | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:323-327` |
| Production call and argument forwarding | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:499-573` |
| EAMxx wrapper and `dropmixnuc` invocation | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:199-231,347-436` |
| Cloud-borne state replacement | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:440-455` |
| Interstitial tendency application | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:458-500` |
| Explicit mixing implementation | `externals/mam4xx/src/mam4xx/ndrop.hpp:1167-1436` |

## Concern 1: Boundary Indexing

### Reachability

EAMxx/MAM4xx defaults select 128 vertical levels in
`components/eamxx/cime_config/config_component.xml:25-33`. For 128 levels, the
default `top_level_mam4xx` is zero in
`components/eamxx/cime_config/namelist_defaults_eamxx.xml:283-289`.

The constructor requires that `top_level_mam4xx` exist but does not validate
its value. `MAMAci::run_impl` forwards it unchanged into `dropmixnuc`.
Therefore, `top_lev=0` is not merely a synthetic input: it is the current
default L128 production input.

### Invalid derived index

The active loop is zero-based and spans `[top_lev,pver)`. At each active level:

```cpp
const int km1 = max(k - 1, top_lev);
eddy_diff_km(k) =
    zn(k) * eddy_diff(k - 1) * csbot(k - 1) * zs(km1);
```

The clamped `km1` protects `zs` and state-neighbor accesses, but it does not
protect `eddy_diff(k-1)` or `csbot(k-1)`. At `k=top_lev=0`, both indices are
`-1`. See `externals/mam4xx/src/mam4xx/ndrop.hpp:1219-1240`.

The legacy Fortran implementation declares an extra interface boundary and
explicitly assigns zero above the first active level before evaluating the
corresponding coefficient:

```fortran
ekk(top_lev-1) = 0
ekkm(k) = zn(k) * ekk(k-1) * zs(km1)
```

See `components/eam/src/physics/cam/ndrop.F90:891-908`. No equivalent C++
boundary initialization is present.

### Executed boundary matrix

Each value ran in a separate bounds-checked 128-level process:

| `top_lev` | Result |
|---:|---|
| `-1` | Exit 134; Kokkos view index `-1`, extent 128 |
| `0` | Exit 134; Kokkos view index `-1`, extent 128 |
| `1` | Exit 0; 25/25 controlled ledgers passed |
| `6` | Exit 0; 25/25 controlled ledgers passed |
| `127` | Exit 0; no bounds or ledger violation |
| `128` | Exit 0 with no active-level work; invalid upper input was not rejected |

The retained driver deliberately skips its host-side substep calculation for
`top_lev<=0`. A `LAUNCH_UPDATE top_lev=0` marker was emitted before the Kokkos
index error, establishing that the abort came from the exact MAM4xx function
rather than from replicated harness arithmetic.

For `top_lev=1` and `6`, changing only `eddy_diff(top_lev-1)` from zero to a
nonzero sentinel changed subcycling and the final state fingerprint. This
separately demonstrates that `dropmixnuc` consumes an in-bounds but locally
uninitialized upper-interface value when `top_lev>0`.

### Classification

**CONFIRMED DEFECT**. A current default production value reaches the exact
index expression, and the 128-level bounds-checked execution violated the view
extent. A GPU build without bounds checks may not abort; it can instead read
invalid device memory or use an arbitrary value. That changes the symptom, not
the invalid index.

### Why production runs can complete

The negative index is evaluated on each L128 MAM4xx ACI call when the active
loop reaches `k=0`, but optimized production builds normally do not enable
Kokkos view bounds checking. An address immediately outside a view can remain
inside mapped GPU memory and therefore need not raise an illegal-address
exception.

The current allocation pattern can further mask the symptom. `eddy_diff` and
`csbot` are separate row-major `(ncol,nlev)` scratch views, and `dropmixnuc`
receives one row subview per column. In a typical unchecked address
calculation, level `-1` for column `c>0` aliases level `nlev-1` of column
`c-1`. Because `dropmixnuc` explicitly sets `eddy_diff(nlev-1)=0`, the
erroneous product often evaluates to zero and accidentally resembles the
intended upper no-flux boundary. Column zero still addresses memory before the
allocation, and no valid interface contract permits cross-column aliasing.

### Anticipated propagation

The invalid values directly affect only `eddy_diff_km(0)`. The effective
column `max_tinv`, `dtmix`, and `nsubmix` are then calculated from that
coefficient; they are not themselves read out of bounds. Because `nsubmix` is
a column scalar, a corrupt level-zero coefficient can change the timestep used
by every active-level update in the explicit loop.

For a finite coefficient, the direct upper-boundary transport term normally
vanishes: `km1=0`, `overlapm(0)=1`, and the upper-neighbor state difference is
zero. Expected finite-value outcomes are:

| Invalid value | Anticipated consequence |
|---|---|
| Zero, or a value that does not control the column reduction | No visible numerical or performance effect |
| Spuriously positive and limiting | Excess `nsubmix`, approximately proportional additional explicit-loop work, GPU tail imbalance, and a different finite-step solution throughout the column |
| Negative and masking the true level-zero loss rate | Too few substeps, possible forward-Euler instability, positivity clipping, artificial sources, and aerosol or droplet-number errors |
| NaN, infinity, or unmapped address | NaN propagation, pathological integer `nsubmix`, excessive looping, or device failure |

These propagation modes follow from the executed index, equations, and memory
layout. Their frequency and magnitude in production GPU simulations have not
been measured. Successful EAMxx runs therefore do not disprove the defect;
they indicate that the unchecked read has usually been masked or non-limiting
in those configurations.

## Concern 2: Disabled Vertical Mixing

### Reachability and affected fields

`enable_aero_vertical_mix` is an accepted Boolean EAMxx parameter and is read
without a restriction that it remain true. Its default is true, so the false
state is configurable rather than the current default.

The explicit loop evolves 25 paired aerosol constituents:

- four modal-number mixing ratios; and
- 21 modal species-mass mixing ratios.

Each has interstitial `I` and cloud-borne `C` representations. Cloud-droplet
number `Nc` is coupled through the same activation sources.

The current switch placement is:

| Explicit update | When switch is false |
|---|---|
| Droplet number | Skipped |
| Interstitial aerosol number | Skipped |
| Interstitial aerosol mass | Skipped |
| Cloud-borne aerosol number | Still applied |
| Cloud-borne aerosol mass | Still applied |
| No-cloud cloud-borne-to-interstitial transfer | Still applied |

See `externals/mam4xx/src/mam4xx/ndrop.hpp:1321-1408` and `:1415-1435`.

### Cloud-fraction-driven repartitioning outside the switch

The phrase “cloud-fraction transfer” refers to same-level partition changes,
not transport between vertical levels. Before the explicit substep loop,
`update_from_newcld` and `update_from_cldn_profile` respond to cloud evolution:

- cloud growth activates an amount `Delta` and applies
  `I <- I-Delta`, `C <- C+Delta`;
- cloud shrinkage returns a fraction of cloud-borne aerosol to interstitial;
  and
- cloud disappearance applies `I <- I+C`, `C <- 0` and removes cloud-droplet
  number at that level.

The no-cloud conversion is also repeated after subcycling. These operations
are local phase repartitioning and, when both paired updates are applied,
preserve aerosol number or mass in `I+C`. They are not the behavior classified
as defective. The defect is that with the switch false, the explicit substep
still advances cloud-borne activation while skipping its paired interstitial
and droplet-number updates, and it can expose an unwritten destination buffer.

### Intended paired ledger

For each aerosol constituent `j`, activation supplies opposite source terms:

```text
dI_j/dt = -S_j + interstitial mixing
dC_j/dt = +S_j + cloud-borne mixing
```

With no-flux boundaries and before positivity clipping, the intended column
ledger is

```text
L_j = sum_k (pdel(k)/g) * (I_j(k) + C_j(k)).
```

Activation cancels in `I+C`, and paired interface fluxes telescope. The
approved test contract interprets the switch as disabling turbulent transport,
not as permitting one side of activation to proceed alone.

### Comparison with Fortran EAMv3

The switch-induced imbalance is not present in the E3SM v3.0.0 Fortran
implementation inspected at revision
`399d4301138617088dd93214123d6c025e061302`. The Fortran code has no equivalent
`enable_aero_vertical_mix` option. For every modal number and species-mass
field, it calls `explmix` once for cloud-borne aerosol with `is_unact=.false.`
and once for interstitial aerosol with `is_unact=.true.`. The first call applies
`+source`; the second applies `-source`.

Thus the activation contribution cancels in the combined interstitial plus
cloud-borne amount, provided neither update is changed by positivity clipping
and no relevant boundary flux is present. This comparison rules out the same
asymmetric option-path defect in Fortran EAMv3; it does not establish exact
conservation for every operation or physical state in that code.

### Ping-pong state exposure

`dropmixnuc` initializes only aerosol time buffer zero
(`externals/mam4xx/src/mam4xx/ndrop.hpp:1586-1607`). Subcycling starts with
`nsav=0`, `nnew=1` and swaps them after the first substep (`:1671-1681` and
`:1286-1299`).

When the switch is false, cloud-borne destination buffer `nnew` is written but
interstitial destination buffer `nnew` is not. Later source calculations and
final tendencies can read that unwritten buffer. Existing
`update_from_explmix` validation initializes both buffers and always passes
`true`, so it masks this state (`externals/mam4xx/src/validation/ndrop/update_from_explmix.cpp:16-31,118-149,168-176`).

### Deterministic test construction

The primary 128-level test isolated switch behavior from the boundary defect by
using `top_lev=1` and explicit zero-flux boundaries. It used:

```text
dtmicro = 1 s
cldn = 1
zn = zs = csbot = 1
I_j(k) = 1 + 0.01*j + 0.001*(k+1)
C_j(k) = 0.2 + 0.002*j + 0.0002*(k+1)
qcld(k) = 10 + 0.01*k
nact = 0.2*D
mact = 0.1*D
D = 0.01, 0.5, 1.0
```

Those three `D` values produced `nsubmix=1`, `2`, and `3`. Destination buffer
one was initialized equal to the physical initial interstitial state in the
primary test, so its result does not depend on poison data. Cloud fraction one
also disables the post-loop no-cloud transfer. All post-state values remained
finite and positive, so positivity clipping did not cause the ledger errors.

Acceptance was

```text
abs(Lpost-Lpre) <= 100*epsilon*max(abs(Lpre),1)
```

for every constituent.

### Executed ledger matrix

| `nsubmix` | Switch | Max relative ledger error | Failed ledgers | Interstitial L1 change | Cloud-borne L1 change | Droplet L1 change | Droplet-source L1 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | enabled | `1.64e-16` | 0/25 | `4.324257` | `4.324357` | `1.180064` | `1.179864` |
| 1 | disabled | `1.65e-3` | 25/25 | `0` | `4.324357` | `0` | `1.179864` |
| 2 | enabled | `2.93e-16` | 0/25 | `212.79845` | `212.79845` | `57.5294` | `56.0656` |
| 2 | disabled | `8.27e-2` | 25/25 | `0` | `216.21285` | `0` | `58.9932` |
| 3 | enabled | `3.12e-16` | 0/25 | `414.502575` | `414.502575` | `110.351675` | `102.889159` |
| 3 | disabled | `1.65e-1` | 25/25 | `0` | `432.4257` | `0` | `117.9864` |

The enabled control conserved every pair to roundoff. Every disabled case
violated every pair ledger. The disabled code diagnosed activation sources and
changed cloud-borne aerosol, but did not apply the paired interstitial change or
the coupled droplet-number change.

### Destination-buffer sensitivity

Only the initial contents of interstitial destination buffer one were then
scaled by 0, 1, and 7:

- Enabled fingerprints were bit-for-bit identical for each `nsubmix`.
- Disabled fingerprints changed at every `nsubmix`.
- At `nsubmix=1`, disabled fingerprints ranged from `30656.8` to `386305.3`.
- At `nsubmix=3`, disabled fingerprints ranged from `41992.4` to `437920.0`.

This proves an unwritten-buffer dependency independently of the paired-source
violation. The scale-one matrix already proves the ledger defect without this
dependency.

### EAMxx application and classification

After `dropmixnuc`, EAMxx installs cloud-borne aerosol by direct replacement and
applies interstitial changes as tendencies over `dt`. See
`eamxx_mam_aci_process_interface.cpp:595-616` and
`eamxx_mam_aci_functions.hpp:440-500`. The inconsistent pair therefore reaches
the EAMxx aerosol state; no later pairing operation in this call path repairs
the ledger.

**CONFIRMED DEFECT**. `false` is a valid production-facing parameter value, and
the exact function executed all 25 paired-ledger violations with deterministic,
positive inputs. The separate buffer-initialization sweep also demonstrated
state dependence on a destination that production does not initialize.

## Additional Subcycling Finding

The investigation corrected an earlier cost-analysis assumption. The source
assigns `dtmin=dtmicro` before a Kokkos minimum reduction, but
`Kokkos::Min<Real>::init` replaces the accumulator with the maximum
representable `Real`. The initial `dtmicro` is not included in the reduction.
For at least one `tinv>1e-6 s^-1`, the as-executed relationship is

```text
nsubmix = floor(dtmicro * max_tinv / 0.9) + 1.
```

Thus `nsubmix` can be one under weak mixing. This matters directly to the
disabled-switch defect because odd/even buffer parity changes which unwritten
state is consumed or returned. It also changes the low end of previous cost
estimates. This is a separate numerical issue and was not assigned one of the
two requested classifications.

## Existing Test Baseline

The configured stock target was built and run before the focused matrix:

```text
/opt/homebrew/bin/cmake --build /private/tmp/esulf/mb \
  --target mam4_ndrop_unit_tests -j2
/private/tmp/esulf/mb/src/tests/mam4_ndrop_unit_tests
```

Result: 51 assertions passed in five test cases. The stock tests cover scalar
helpers but not the `top_lev` boundary, disabled switch, full pair ledger, or
production ping-pong initialization.

Focused cases were invoked as separate processes, for example:

```text
ndrop_evidence boundary 0
ndrop_evidence boundary 1
ndrop_evidence boundary 127
ndrop_evidence boundary 128
ndrop_evidence ledger 1 0.01 1
ndrop_evidence ledger 0 0.01 1
ndrop_evidence ledger 1 0.5 1
ndrop_evidence ledger 0 0.5 1
ndrop_evidence ledger 1 1.0 1
ndrop_evidence ledger 0 1.0 1
```

The last ledger argument is the initial destination-buffer scale.

The exact standalone source and configuration instructions are retained in
`docs/analysis/evidence/ndrop_evidence_driver.cpp` and
`docs/analysis/evidence/README.md`.

Final machine checks against that attached source reported:

```text
EXACT_BOUNDARY_AUDIT: 6 cases, 0 failures
ATTACHED_LEDGER_AUDIT: 18 cases, 0 failures
ATTACHED_UPPER_INTERFACE_AUDIT: 4 cases, 0 failures
```

## High Droplet-Activation Attribution Finding

Output-field availability and the activation-input path were additionally
checked on 2026-07-15 against upstream E3SM `master` revision
`276942b7b4b504e3b4d59d3691f68df59300f8ae` and its MAM4xx gitlink
`b154d7f23b2e910242abbbab5dc72250d7b950c5`.

The SHOC field consumed directly by EAMxx MAM4xx ACI is `w_variance`, not the
SHOC prognostic field `tke`. EAMxx interpolates `w_variance` from midpoints to
interfaces, constructs

```text
tke_int = 1.5 * interp_mid_to_int(w_variance)
wsub_raw(k) = sqrt(0.5*(tke_int(k)+tke_int(k+1))*(2/3))
wsub(k) = max(wsubmin, wsub_raw(k))
```

and passes `wsub` into `dropmixnuc`. MAM4xx then constructs

```text
wtke(k) = max(wsub(k), 0.1 m/s)
```

and passes `wtke` to `activate_modal`. Thus all raw activation velocities below
`0.1 m/s` receive the same activation velocity. The configurable default
`wsubmin` is `0.001 m/s`, so the second hard-coded floor is normally the
controlling lower bound. The interpolation, velocity floors, and
`w_variance` time-level treatment are priority checks for explaining high
droplet activation in EAMxx relative to EAMv3.

Current output fields provide the major environmental context: droplet-number
response, six CCN spectra, `w_variance`, SHOC diffusivity and cloud fields,
modal aerosol size and partitioning, vertical geometry and density, cloud
fractions, and P3 state and process response. They do not directly expose
`wsub_raw`, `wsub`, `wtke`, maximum supersaturation, mode-resolved activation
fractions and rates, activated number/mass fluxes, explicit pre-ACI aerosol
copies, or the actual adaptive mixing work. Existing output therefore supports
correlation and partial offline reconstruction, but not definitive causal
attribution.

The complete field lists, focused and extended cadences, averaging choices,
direct-output availability notes, and missing-attribution diagnostics are
maintained in
[`docs/design/eamxx_mam4xx_ndrop.md`](../design/eamxx_mam4xx_ndrop.md#comprehensive-existing-output-recommendation-for-activation-attribution).

## Limitations

- The focused driver executes the exact MAM4xx inner function with the
  production level count and buffer contract, but it is a deterministic
  standalone column rather than a full EAMxx atmosphere run.
- A full EAMxx/CIME case was not created because the local macOS host is not a
  supported E3SM machine under the repository's case-creation rules.
- Production reachability was established separately from EAMxx configuration
  and call-path source; it was not inferred from synthetic inputs.
- The boundary test used a CPU bounds-checked Kokkos backend. GPU symptoms can
  differ, but the zero-based `k-1` index expression is backend-independent.
- No fix was implemented or tested, and no claim is made about corrected
  answers or bit-for-bit behavior.
- No production performance profile was collected. The cost report remains an
  analytic estimate, now corrected for one-substep behavior.

## Recommended Actions

1. Restore an explicit no-flux upper-interface boundary in C++ and validate
   `0 <= top_lev < pver` before launching the device kernel. Retain
   `top_lev=0` as valid for the current L128 configuration.
2. Add bounds tests for `0`, `1`, `6`, `pver-1`, `-1`, and `pver` on a
   bounds-checked host backend and at least one GPU backend.
3. Define the disabled-switch contract in the interface. Under the approved
   contract, suppress turbulent transport terms while retaining symmetric
   `-S/+S` activation transfers and coupled droplet-number evolution.
4. Write every ping-pong destination on every substep. A disabled operation
   should copy through its old state rather than leave the destination
   untouched.
5. Add per-constituent `I+C` ledger tests for one, even, and odd substep counts,
   plus a full EAMxx application test covering cloud-borne replacement and
   interstitial tendency application.
6. Audit whether `dtmin=dtmicro` was intended to cap the Kokkos reduction. If
   so, implement the cap explicitly and add weak/strong-mixing substep tests.
7. Before interpreting the planned five-day `ne256` subcycle matrix, correct
   the level-zero boundary or demonstrate an equivalent safe boundary. Then
   record actual `nsubmix`, a fenced `dropmixnuc` region time, and preferably
   `max_tinv` with its limiting level. Keep `enable_aero_vertical_mix=true`
   until the disabled-switch defect is corrected and tested.

## Related Documentation

- `docs/design/eamxx_mam4xx_ndrop.md`: reconstructed as-built design, equation
  map, interface contract, executed evidence, and validation plan.
- `docs/analysis/eamxx_ndrop_aci_cost_report.md`: broader ACI cost, resolution,
  subcycling, `wsub`, and SHOC/CLUBB comparison.
- `docs/analysis/evidence/`: standalone evidence driver and reproduction
  instructions.

## Revision History

| Date | Author / Agent | Change Summary |
|---|---|---|
| 2026-07-21 | User / Codex | Recorded MAM4xx issue 564 as the public tracker for the level-zero boundary defect. |
| 2026-07-21 | User / Codex | Verified against E3SM v3.0.0 that Fortran applies both signs of every activation transfer and does not contain the EAMxx disabled-mixing switch defect. |
| 2026-07-21 | User / Codex | Clarified the conservative same-level cloud-fraction-driven repartitioning that remains outside the vertical-mixing switch and distinguished it from the disabled-substep defect. |
| 2026-07-15 | User / Codex | Confirmed applicability to the current checkout, explained why unchecked production runs can complete, documented anticipated subcycling impacts, and connected the defects to the planned five-day performance matrix. |
| 2026-07-15 | Codex | Added high-activation attribution findings and current-versus-missing diagnostic inventory. |
| 2026-07-14 | User / Codex | Recorded executed boundary and paired-ledger evidence and assigned both final confirmed-defect classifications. |
