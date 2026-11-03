# EAMxx-MAM4xx `ndrop` Activation and Explicit Mixing

> **[RECONSTRUCTED FROM EXISTING CODE — reviewed and approved 2026-07-14]**
>
> This document describes the implementation inspected at E3SM revision
> `12e44bb469efb0d6c0cb94c7b796a86967883ff5` and MAM4xx revision
> `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`. It is an approved as-built design
> reconstruction, not a production redesign. Two suspected behaviors are called
> out as investigation targets; neither was classified as a confirmed defect
> until executable evidence was available. The local evidence uses Kokkos 4.5.1
> with its Serial backend, debug view-bounds checking, ASan, and UBSan.

The output-diagnostic availability and activation-input trace documented below
were additionally checked against upstream E3SM `master` revision
`276942b7b4b504e3b4d59d3691f68df59300f8ae` and its MAM4xx gitlink
`b154d7f23b2e910242abbbab5dc72250d7b950c5` on 2026-07-15.

The diagnostic implementation was rebased on 2026-08-17 at E3SM commit
`e2c578e40f`, based on upstream E3SM revision
`fe750f50e98187afe02103349d8363acd7745a9a`. It uses MAM4xx commit
`0d076d8ebb004a18f05f51f97b08edd9def6826f`, based on upstream MAM4xx revision
`8fbf7750bafa8b14533d72397aeb027184d0ec37`, and EKAT revision
`3a1485aecfed7a036e4fa30ce845498fd37975e2`. This newer MAM4xx revision always
performs aerosol vertical mixing and no longer implements a disabled-mixing
path. It also includes MAM4xx PR #567 (`0b100fb3d`), which replaces the unsafe
`eddy_diff(k-1)` access at `top_lev=0` with the clipped `eddy_diff(km1)` index.

## Scope and Status

This document covers the EAMxx call to MAM4xx `dropmixnuc`, with emphasis on:

- the vertical indexing contract and the `top_lev` boundary;
- activation and vertical-mixing subcycling;
- interstitial and cloud-borne aerosol transfer semantics;
- the behavior of `enable_aero_vertical_mix`;
- the state and tendency updates returned to EAMxx; and
- the validation used to classify the two suspected behaviors.

No production-physics modification is authorized by this document. The
validation plan was approved before investigation-specific test code was
written or executed.

The user authorized a diagnostic-only pull request on 2026-08-11 to expose the
actual explicit mixing-loop count. That authorization covers the three output
fields described in the instrumentation section below and their tests. It does
not authorize a change to the `ndrop` equations, loop bounds, activation,
mixing, clipping, or model coupling.

The proposed isolated Forward Euler, IMEX, and Backward Euler comparison below
is a new numerical-method design that was **approved on 2026-07-23**. This
approval authorizes the staged standalone test and validation work described
below; it does not authorize production implementation.

> **[VALIDATION-GATE REMEDIATION AMENDMENT — approved 2026-07-30]**
>
> An independent review reran and examined the Forward Euler and IMEX gates
> and assigned both a `FAIL`. The numerical formulations were judged sound in
> their tested normal-magnitude regime, but the standalone state transition
> omitted the current MAM4xx zero-cloud resuspension/reset, finite-input
> overflow could produce nonfinite diagnostics without failure, and parts of
> the already-approved convergence and per-field clipping test contract were
> not implemented. The user approved the remediation contract on 2026-07-30.
>
> The corrected Forward Euler and IMEX gates passed independent extra-high
> review on 2026-08-06 at candidate commit `2c4427e`. The current locked-stack
> Backward Euler method-local gate passed on 2026-08-24 at commits
> `2f43c8c1e4` and `c6f7f18a5d`. Cross-method validation and production
> implementation remain separate, unauthorized stages.

### Phase 0 design gate

This document is Phase 0 / the design gate of the **Isolated `ndrop`
Integrator Validation and Performance Plan**. Once approved, it is intended to
enable the later harness and Forward Euler phase, IMEX gate, Backward Euler
gate, cross-scheme validation, performance evaluation, and final report. Phase
0 supplies the canonical operator equations, boundary and indexing contracts,
column invariants and ledgers, and the classification contract for the two
originally suspected behaviors: `top_lev` boundary indexing and disabled
vertical mixing.

This Phase 0 approval would authorize only the separate standalone test plan.
Any future production implementation of IMEX or Backward Euler must reference
this design and the Isolated `ndrop` Integrator Validation and Performance Plan
and must pass a separate production design and review gate before entering an
EAMxx or MAM4xx production path.

The focused results and recommendations are summarized in
`docs/analysis/eamxx_mam4xx_ndrop_defect_brief.md`.

## Scientific Basis

### Process representation

`ndrop` diagnoses cloud-droplet activation and transfers modal aerosol number
and species mass between interstitial and cloud-borne representations. It also
performs an explicit vertical mixing calculation that couples cloud droplet
number, interstitial aerosol, and cloud-borne aerosol.

The activation parameterization follows the Abdul-Razzak and Ghan family of
parameterizations. The implementation also cites Ghan et al. for droplet
nucleation and cloud-fraction treatment.

### Prognostic quantities

For each vertical level `k` and aerosol constituent `j`, define:

- \(I_{j,k}\): interstitial aerosol number or mass mixing ratio;
- \(C_{j,k}\): cloud-borne aerosol number or mass mixing ratio;
- \(N_{c,k}\): cloud-droplet number mixing ratio, represented locally as
  `qcld`;
- \(S_{j,k}\): activation transfer rate from interstitial to cloud-borne
  aerosol; and
- \(S_{N,k}\): corresponding source of cloud-droplet number.

The aerosol state consists of 25 paired constituents: four modal-number
mixing ratios and 21 modal species-mass mixing ratios. Each constituent has an
interstitial and a cloud-borne representation. The droplet-number field is not
one member of those aerosol pairs.

The intended local activation ledger for each aerosol constituent is

\[
  \left.\frac{dI_{j,k}}{dt}\right|_{\mathrm{act}}=-S_{j,k},
  \qquad
  \left.\frac{dC_{j,k}}{dt}\right|_{\mathrm{act}}=+S_{j,k}.
\]

Thus activation alone preserves \(I_{j,k}+C_{j,k}\). Cloud decay or absence
can return cloud-borne material to the interstitial representation by the
inverse paired transfer.

### Explicit mixing equations

For one explicit substep of length \(\Delta t_m\), the implementation updates
droplet number as

\[
N_{c,k}^{n+1}=\max\left(0, N_{c,k}^{n}+\Delta t_m\left[
S_{N,k}+D_k^+\left(O_k^+N_{c,k+1}^{n}-N_{c,k}^{n}\right)
+D_k^-\left(O_k^-N_{c,k-1}^{n}-N_{c,k}^{n}\right)
\right]\right).
\]

For an interstitial aerosol constituent,

\[
I_{j,k}^{n+1}=\max\left(0, I_{j,k}^{n}+\Delta t_m\left[
-S_{j,k}+D_k^+\left(I_{j,k+1}^{n}-I_{j,k}^{n}
+C_{j,k+1}^{n}(1-O_k^+)\right)
+D_k^-\left(I_{j,k-1}^{n}-I_{j,k}^{n}
+C_{j,k-1}^{n}(1-O_k^-)\right)
\right]\right),
\]

and for its cloud-borne pair,

\[
C_{j,k}^{n+1}=\max\left(0, C_{j,k}^{n}+\Delta t_m\left[
S_{j,k}+D_k^+\left(O_k^+C_{j,k+1}^{n}-C_{j,k}^{n}\right)
+D_k^-\left(O_k^-C_{j,k-1}^{n}-C_{j,k}^{n}\right)
\right]\right).
\]

Here \(D_k^+\) and \(D_k^-\) are the mixing rates across the lower and upper
interfaces of cell `k`, and \(O_k^+\) and \(O_k^-\) are cloud-overlap factors.
Before positivity clipping, summing the paired aerosol equations cancels the
activation source. It also resuspends cloud-borne aerosol mixed into the
non-cloudy fraction into the interstitial representation.

The pressure-thickness-weighted column ledger for constituent `j` is

\[
  L_j=\sum_k \frac{\Delta p_k}{g}\left(I_{j,k}+C_{j,k}\right).
\]

With paired interface fluxes and no-flux top and bottom boundaries, vertical
mixing should preserve \(L_j\) before clipping. Separate positivity clipping
of \(I\) and \(C\) can break exact conservation when a controlled test enters
the clipping regime; conservation tests therefore must also report whether a
limiter was active.

### Assumptions and limitations

- C++ vertical indices are zero-based.
- The column has `pver = mam4::nlev` levels with valid indices
  `[0, pver-1]`.
- `top_lev` is intended to identify the first active MAM4xx level.
- The explicit mixing formulation assumes valid no-flux boundary
  coefficients at the first and last active levels.
- The current public interface does not state or enforce a numeric range for
  `top_lev`.
- The physical meaning of disabling aerosol vertical mixing is not documented
  precisely enough to determine whether activation transfers should remain
  active. This is a design ambiguity requiring review.

## Numerical Approach

### Host-model coupling and process order

For each EAMxx column, `MAMAci::run_impl` gathers thermodynamic, cloud,
turbulence, aerosol, and geometry fields and invokes
`MAMAci::call_function_dropmixnuc`. That wrapper builds MAM4xx column views and
calls `mam4::ndrop::dropmixnuc` in a Kokkos team kernel.

Within `dropmixnuc`, the relevant process order is:

1. initialize work arrays and load interstitial and cloud-borne aerosol into
   first of the two alternating work arrays;
2. diagnose activation and transfer associated with new cloud;
3. diagnose activation rates for existing cloud and cloud base;
4. calculate a stable explicit mixing substep and repeat the explicit update;
5. transfer cloud-borne aerosol back to interstitial aerosol where cloud is
   absent; and
6. return cloud-borne state directly and interstitial changes as tendencies.

In EAMxx, heterogeneous freezing consumes the old cloud-borne aerosol state
before the returned `ndrop` cloud-borne state is installed. EAMxx then replaces
cloud-borne aerosol with the returned state and applies the returned
interstitial tendencies over the process time step.

### ACI subgrid velocity and the velocity consumed by `ndrop`

The direct turbulence input to EAMxx MAM4xx ACI is the SHOC midpoint field
`w_variance`, not the SHOC prognostic field `tke`. EAMxx linearly interpolates
`w_variance` from midpoints to interfaces using layer thickness `dz`, with the
top and bottom boundary interface values set equal to the corresponding
midpoint values. It then defines the ACI-local interface work field

\[
  \mathrm{tke}_{int}(i)=\frac{3}{2}
  \operatorname{interp}_{mid\rightarrow int}
  \left(w_{\mathrm{variance}}\right)_i .
\]

For active midpoint level `k`, EAMxx calculates

\[
  w_{\mathrm{sub,raw},k}=
  \sqrt{\frac{1}{2}\left[
    \mathrm{tke}_{int,k}+\mathrm{tke}_{int,k+1}
  \right]\frac{2}{3}},
\]

and applies the configurable lower bound

\[
  w_{\mathrm{sub},k}=\max(w_{\mathrm{submin}},
  w_{\mathrm{sub,raw},k}).
\]

The current default is `wsubmin=0.001 m s-1`. Levels above `top_lev` are set
directly to that value. The `wsub` profile is passed to `dropmixnuc`, which
applies a second, hard-coded lower bound:

\[
  w_{\mathrm{tke},k}=\max(w_{\mathrm{sub},k},0.1\ \mathrm{m\ s^{-1}}).
\]

`wtke`, rather than `w_variance`, SHOC `tke`, or the pre-floor `wsub`, is the
velocity passed to `activate_modal`. Consequently, every raw activation
velocity below `0.1 m s-1` is promoted to `0.1 m s-1`. This floor, the
midpoint-to-interface interpolation, and the time level of `w_variance` are
first-order comparison targets when diagnosing high droplet activation in
EAMxx relative to EAMv3. The separately calculated `wsig` and `wsubice` fields
are bounded for other uses but are not the liquid-droplet activation velocity
passed to `activate_modal`.

### Subcycling

The source assigns `dtmin = dtmicro` before invoking a Kokkos minimum
reduction. That assignment is not an input to the reduction: the
`Kokkos::Min<Real>` reducer initializes its accumulator to
`reduction_identity<Real>::min()` (the maximum representable `Real`) and
overwrites `dtmin` with the reduction result. The as-executed calculation is
therefore

\[
  \Delta t_{\mathrm{red}}=\min_k\left(\frac{1}{t_{\mathrm{inv},k}}
  \;\text{for}\;t_{\mathrm{inv},k}>10^{-6}\ \mathrm{s}^{-1}\right),
\]

not the minimum of that value and \(\Delta t_{\mathrm{micro}}\). If no level
enters the conditional reduction, the reduction identity is retained.

\[
  \Delta t_{\mathrm{candidate}}=0.9\Delta t_{\mathrm{red}},\qquad
  n_{\mathrm{submix}}=\left\lfloor
  \frac{\Delta t_{\mathrm{micro}}}{\Delta t_{\mathrm{candidate}}}
  \right\rfloor+1,
\]

\[
  \Delta t_m=\frac{\Delta t_{\mathrm{micro}}}{n_{\mathrm{submix}}}.
\]

Consequently, weak mixing can produce `nsubmix=1`; the previously reconstructed
minimum of two was incorrect. There is no hard upper cap on `nsubmix` in the
inspected C++ implementation. The activation rates are limited using the
interface mixing rate before subcycling. Whether the ineffective
`dtmin=dtmicro` preassignment is intentional or a porting error is outside the
two primary classifications, but the test matrix must cover one, even, and odd
substep counts because buffer exposure depends on parity.

### Fortran EAMv3 explicit mixing subloop

The Fortran EAMv3 implementation also uses an adaptive explicit mixing
subloop. For each column, `dropmixnuc`:

1. reads the turbulence diffusivity `kvh` and forms the interface diffusivity
   `ekd`, normally limited to `[0.01,100] m2 s-1`;
2. replaces `ekd` at diagnosed cloud bases with `wbar/zs`;
3. forms the density-weighted interface coefficients `ekk`, with explicit zero
   values above and below the active column;
4. calculates `tinv=ekkp+ekkm`, `dtmix`, and `nsubmix`; and
5. executes `do n=1,nsubmix`, updating droplet number and every cloud-borne and
   interstitial aerosol number and mass field with `explmix`.

The Fortran call receives `dtmicro=deltatin` from `microp_aero`. Its
`nsubmix_bnd=min(nsubmix,100)` value only selects a diagnostic histogram bin;
the actual loop still runs the full `nsubmix` count. The substeps remain
sequential because each substep uses the aerosol state produced by the
previous one. Vertical levels are expressed as Fortran array operations and
loops rather than EAMxx Kokkos team/vector ranges.

Thus explicit mixing subcycling was inherited from Fortran; it was not
introduced by the C++ port. At the older evidence revision, important port
differences included the C++ boundary representation, Kokkos execution
structure, and the EAMxx `enable_aero_vertical_mix` option described below.

### Vertical indexing and boundaries

Most active-level loops use the half-open Kokkos range
`[top_lev, pver)`. Neighbor indices are usually clamped as

```text
kp1 = min(k + 1, pver - 1)
km1 = max(k - 1, top_lev)
```

At the current diagnostic revision, the upper-interface mixing coefficient is
formed as

```text
eddy_diff_km(k) = zn(k) * eddy_diff(km1) * csbot(km1) * zs(km1)
```

This is the repair from MAM4xx PR #567. In particular, `top_lev=0` no longer
causes a `-1` access in this expression.

#### Historical evidence-revision behavior (boundary indexing)

The earlier evidence revision formed the coefficient as

```text
eddy_diff_km(k) = zn(k) * eddy_diff(k - 1) * csbot(k - 1) * zs(km1)
```

Consequences at that revision were:

- when `top_lev=0`, the first active iteration reads `eddy_diff(-1)` and
  `csbot(-1)`;
- clamping `km1` does not protect those two reads;
- when `top_lev>0`, the same expression reads level `top_lev-1`, but
  `dropmixnuc` initializes `eddy_diff` and `csbot` only from `top_lev` through
  the surface, so the upper-boundary value is not established locally; and
- no range check or assertion rejects `top_lev<0`, `top_lev>=pver`, or the
  production value `top_lev=0`.

The legacy Fortran implementation uses one-based active levels and declares
the interface array with an additional zero boundary index. It explicitly sets
the interface above the top active level to zero before evaluating the same
form of mixing rate. The corresponding C++ boundary initialization was not
evident in the older call path. The executed evidence for this repaired defect
remains documented below so results from affected revisions stay reproducible.

### `enable_aero_vertical_mix` behavior

At the current diagnostic revision, EAMxx still passes this Boolean argument
for interface compatibility, but MAM4xx does not use its value. The explicit
mixing loop always updates droplet number, interstitial aerosol, and
cloud-borne aerosol. Current results therefore cannot be described as an
enabled-versus-disabled comparison.

#### Historical evidence-revision behavior (disabled vertical mixing)

At MAM4xx revision `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`, the switch
controlled these operations:

| Update | Switch condition at the evidence revision |
|---|---|
| Mix/update cloud-droplet number `qcld` | Conditional |
| Mix/update interstitial aerosol number | Conditional |
| Mix/update interstitial aerosol mass | Conditional |
| Mix/update cloud-borne aerosol number | Unconditional |
| Mix/update cloud-borne aerosol mass | Unconditional |
| No-cloud cloud-borne-to-interstitial transfer after subcycling | Unconditional |

The activation source has opposite signs in the interstitial and cloud-borne
equations. Consequently, disabling only the interstitial update while retaining
the cloud-borne update does not preserve cancellation between the matching
interstitial loss and cloud-borne gain in the equations actually advanced.

The aerosol arrays used two time buffers. At entry, only buffer zero was loaded
from the EAMxx state. The source and destination indices swapped after each
substep. When aerosol vertical mixing was disabled, the interstitial
destination array was not written by the conditional `explmix` calls, but
subsequent substeps and final tendency construction could read the swapped
array. Validation fixtures at that revision initialized both arrays, which did
not reproduce the production initialization pattern.

Cloud-fraction-driven same-level repartitioning occurred outside the explicit
calls disabled by the switch. This is not vertical transport. For each paired
interstitial/cloud-borne constituent `(I,C)` at one model level:

- when cloud fraction grows by more than `0.01`, activation computes a positive
  amount `Delta` and applies `I <- I-Delta`, `C <- C+Delta`;
- when cloud fraction shrinks, a fraction of `C` is returned with the same
  paired form but `Delta<0`; and
- when cloud is absent, the code applies `I <- I+C`, `C <- 0` and removes the
  corresponding cloud-droplet number.

Before positivity clipping or the disabled-array defect, these local paired
operations preserved `I+C`. They remained active when
`enable_aero_vertical_mix=false` because the switch disabled selected explicit
neighbor-mixing updates, not all activation and cloud-phase repartitioning in
`dropmixnuc`. Retaining conservative local repartitioning while disabling
turbulent transport is physically coherent; the confirmed historical defect
was that the substep loop retained the cloud-borne `+S` update while
suppressing the corresponding interstitial `-S` and droplet-number updates and
could read an unwritten alternating work array.

### Comparison with Fortran EAMv3

The same switch-induced defect is not present in the E3SM v3.0.0 Fortran
implementation inspected at revision
`399d4301138617088dd93214123d6c025e061302`. That implementation has no
`enable_aero_vertical_mix` option. For every modal number and species-mass
field, it calls `explmix` for both aerosol forms using the same activation
amount:

- the cloud-borne call uses `is_unact=.false.` and adds the activation amount;
- the interstitial call uses `is_unact=.true.` and subtracts the activation
  amount.

The activation terms therefore cancel in `I+C` when both updates remain
positive and there are no relevant boundary fluxes. Positivity clipping can
still alter exact conservation if an update would otherwise become negative,
but EAMv3 does not contain the EAMxx option path that applies only the
cloud-borne increase.

## Module and Interface Design

### Public call path

The inspected production path is:

```text
MAMAci::run_impl
  -> MAMAci::call_function_dropmixnuc
    -> mam4::ndrop::dropmixnuc
      -> mam4::ndrop::update_from_explmix
        -> mam4::ndrop::explmix
```

The primary locations are:

- `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp`
- `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp`
- `externals/mam4xx/src/mam4xx/ndrop.hpp`

At the current diagnostic revision, the `dropmixnuc` signature begins near
`externals/mam4xx/src/mam4xx/ndrop.hpp:1279`, and the EAMxx wrapper signature
begins near
`components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:230`.

### Source evidence map

The current diagnostic path is anchored by function name as follows:

| Current evidence | Function or location |
|---|---|
| Repaired `km1` indexing, timestep selection, and explicit loop | `mam4::ndrop::update_from_explmix` in `externals/mam4xx/src/mam4xx/ndrop.hpp` |
| Compatibility Boolean accepted but unused | `mam4::ndrop::dropmixnuc` signature in the same file |
| Exact `nsubmix` returned to EAMxx | `MAMAci::call_function_dropmixnuc` in `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp` |
| Exact-count output fields and precision bounds | `MAMAci::create_requests`, `MAMAci::initialize_impl`, and `detail::accumulate_ndrop_nsubmix` in the EAMxx MAM ACI interface and functions files |

The detailed line anchors below apply to the older evidence revisions. They
are retained to make the executed defect investigation reproducible and should
not be used as current-main source locations:

| Evidence | File and lines |
|---|---|
| Parameter presence checks, but no value/range checks | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:41-53` |
| Parameter values loaded | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:323-327` |
| Production forwarding into `call_function_dropmixnuc` | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:499-573` |
| Heterogeneous freezing before aerosol state installation | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:595-616` |
| EAMxx construction of interstitial/cloud-borne inputs and `dropmixnuc` call | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:347-436` |
| `w_variance` midpoint-to-interface interpolation and ACI-local `tke` construction | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:40-67` |
| `wsub`, `wsig`, and `wsubice` construction and bounds | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:69-97` |
| `wsubmin` default and ACI configuration | `components/eamxx/cime_config/namelist_defaults_eamxx.xml:283-289` |
| Hard-coded `wmixmin=0.1 m s-1`, `wtke` construction, and activation calls | `externals/mam4xx/src/mam4xx/ndrop.hpp:1521-1665` |
| Direct cloud-borne replacement | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:440-455` |
| Interstitial tendency application | `components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:458-500` |
| Explicit update equations and positivity limiter | `externals/mam4xx/src/mam4xx/ndrop.hpp:1092-1164` |
| Derived boundary indices and `k-1` reads | `externals/mam4xx/src/mam4xx/ndrop.hpp:1219-1240` |
| Explicit-step selection and alternating-array swap | `externals/mam4xx/src/mam4xx/ndrop.hpp:1264-1299` |
| Conditional droplet and interstitial updates; unconditional cloud-borne updates | `externals/mam4xx/src/mam4xx/ndrop.hpp:1321-1408` |
| Unconditional no-cloud resuspension | `externals/mam4xx/src/mam4xx/ndrop.hpp:1415-1435` |
| Boundary-array initialization and buffer-zero state load | `externals/mam4xx/src/mam4xx/ndrop.hpp:1544-1607` |
| Subcycle call and final aerosol outputs | `externals/mam4xx/src/mam4xx/ndrop.hpp:1671-1757` |
| Legacy zero-flux upper boundary | `components/eam/src/physics/cam/ndrop.F90:891-908` |
| Existing validation hardcodes `top_lev=6`, enables mixing, and initializes both buffers | `externals/mam4xx/src/validation/ndrop/update_from_explmix.cpp:16-31,118-149,168-176` |
| Existing scalar `explmix` smoke test | `externals/mam4xx/src/tests/mam4_ndrop_unit_tests.cpp:184-220` |
| L72/L128 `top_level_mam4xx` defaults | `components/eamxx/cime_config/namelist_defaults_eamxx.xml:283-289` |
| Default MAM4xx build selects 128 levels | `components/eamxx/cime_config/config_component.xml:25-33` |
| EAMxx level count forwarded to MAM4xx | `components/eamxx/src/physics/mam/CMakeLists.txt:23-27` and `externals/mam4xx/src/mam4xx/aero_config.hpp.in:22-24` |
| `mac_aero_mic` membership and resolution-dependent subcycle defaults | `components/eamxx/cime_config/namelist_defaults_eamxx.xml:558-573` |
| Generic process subcycle loop and `dt_sub=dt/number_of_subcycles` | `components/eamxx/src/share/atm_process/atmosphere_process.cpp:150-221` |
| Existing ACI input/output field registrations | `components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp:91-174` |
| Explicit mixing loop used for work-count estimates | `externals/mam4xx/src/mam4xx/ndrop.hpp:1286-1413` |
| P3 history fields explicitly noted as not averaged over subcycles | `components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp:147-150` |

### Configuration inputs

| Input | Role | Units | Dimensions | Valid range/current contract |
|---|---|---:|---|---|
| `top_level_mam4xx` / `top_lev` | First active vertical level | index | scalar | Zero-based in C++; no enforced range |
| `enable_aero_vertical_mix` | Retained EAMxx compatibility input | boolean | scalar | Current MAM4xx ignores the value and always performs mixing; older revisions used the defective partial-disable behavior documented above |
| `dt` / `dtmicro` | Interval advanced by one `dropmixnuc` call | s | scalar | Positive and finite assumed; not checked locally |
| `temp` | Midlevel temperature | K | `pver` | Positive and finite assumed |
| `pmid`, `pint` | Midlevel and interface pressure | Pa | `pver`, host interface extent | Positive, monotone atmospheric profile assumed |
| `pdel`, `rpdel` | Layer pressure thickness and reciprocal | Pa, Pa-1 | `pver` | Positive, finite, mutually consistent assumed |
| `zm` | Midlevel geopotential height | m | `pver` | Must decrease with increasing `k`; asserted in `dropmixnuc` |
| `state_q` | Host state including interstitial modal number and mass | # kg-1 or kg kg-1 | `pver x pcnst` | Nonnegative aerosol entries assumed |
| `ncldwtr` | Initial cloud-droplet number mixing ratio | # kg-1 | `pver` | Nonnegative and finite assumed |
| `v_diffusivity` | Interface turbulent diffusivity input | m2 s-1 | host interface extent | Finite assumed; values used by `ndrop` are clipped to `[0.01,100]` m2 s-1 except the zero-flux surface |
| `cldn`, `cldo` | Current and previous liquid-cloud fraction | 1 | `pver` each | Intended `[0,1]`; no local assertion |
| `wsub` | Subgrid vertical velocity used by activation | m s-1 | `pver` | Finite assumed; locally floored at 0.1 m s-1 into `wtke` |
| `qqcw_fld` | Cloud-borne modal number and mass, input and replacement output | # kg-1 or kg kg-1 | `25 x pver` | Nonnegative assumed; output clipped nonnegative |
| `lspectype_amode`, `lmassptr_amode`, `numptr_amode`, `nspec_amode`, `mam_idx`, `mam_cnst_idx` | Aerosol species/mode mapping | index/count | compile-time mode/species extents | Must be consistent, in bounds, and use the expected legacy one-based mapping values where the C++ code subtracts one |
| `specdens_amode`, `spechygro`, `voltonumbhi_amode`, `voltonumblo_amode`, `exp45logsig`, `alogsig`, `aten` | Aerosol physical and distribution parameters | mixed | compile-time mode/species extents or scalar | Finite, physically admissible MAM4 values assumed |

The current namelist defaults set `top_level_mam4xx=6` for 72 levels and
`top_level_mam4xx=0` for 128 levels. The EAMxx component default uses 128
levels, so `top_lev=0` is reachable from the current default production
configuration. The process constructor checks that the parameters exist but
does not validate their values.

### State and outputs

| Quantity | Entry representation | Local evolution | EAMxx application |
|---|---|---|---|
| Interstitial aerosol number/mass | `progs.q_aero_i` | Paired transfer and explicit mixing in `raercol` | Returned change divided by `dt`, then applied as a tendency |
| Cloud-borne aerosol number/mass | `progs.q_aero_c` | Paired transfer and explicit mixing in `raercol_cw` | Returned as a replacement state after heterogeneous freezing |
| Cloud-droplet number | `qcld` work field | Activation source and explicit mixing | Returned through the ACI process output path |
| Activation rates | Temporary `nact`/`mact` | Diagnosed and limited before explicit substeps | Diagnostic/work data, not a separate aerosol state |

The 25 aerosol pairs comprise:

| Mode | Number fields | Species-mass fields | Total paired constituents |
|---|---:|---:|---:|
| Accumulation | 1 | 7 | 8 |
| Aitken | 1 | 4 | 5 |
| Coarse | 1 | 7 | 8 |
| Primary carbon | 1 | 3 | 4 |
| **Total** | **4** | **21** | **25** |

`dropmixnuc` does not directly update water-vapor mixing ratio, cloud-liquid
mass, temperature, or pressure.

The remaining outputs and internal work fields are:

| Name(s) | Role | Units | Dimensions/interpretation |
|---|---|---:|---|
| `tendnd` | Total droplet-number tendency consumed by cloud microphysics | # kg-1 s-1 | `pver` |
| `factnum` | Modal number activation fraction, also consumed by heterogeneous freezing | 1 | modes x `pver` |
| `ndropcol` | Layer contribution to column droplet number | # m-2 | `pver` diagnostic |
| `ndropmix` | Droplet-number mixing contribution | # kg-1 s-1 | `pver` diagnostic |
| `nsource` | Droplet-number activation source | # kg-1 s-1 | `pver` diagnostic |
| `wtke` | Activation velocity after the local lower bound | m s-1 | `pver` diagnostic/work field |
| `ccn` | CCN concentrations at prescribed supersaturations | # m-3 | `pver x supersaturation bins` |
| `coltend`, `coltend_cw` | Pressure-weighted interstitial/cloud-borne column tendency contributions | # m-2 s-1 or kg m-2 s-1 | `25 x pver` diagnostics |
| `raercol`, `raercol_cw` | Two alternating interstitial/cloud-borne work states | # kg-1 or kg kg-1 | `pver x 2 x 25` |
| `nact`, `mact` | Modal number/mass activation rates | s-1 | `pver x modes` |
| `eddy_diff`, `zn`, `csbot`, `zs`, `eddy_diff_kp`, `eddy_diff_km` | Geometry and explicit mixing coefficients | mixed; final coefficients s-1 | `pver` work fields |
| `overlapp`, `overlapm` | Cloud overlap below/above | 1 | `pver`, expected `[0,1]` |
| `qncld`, `srcn`, `source` | Saved droplet number and explicit source work fields | # kg-1 and # kg-1 s-1 or kg kg-1 s-1 | `pver` |
| `dz`, `csbot_cscen`, `raertend`, `qqcwtend` | Layer geometry, density ratio, and local tendency work fields | mixed | `pver` |

### Error and edge-case handling

The inspected implementation relies on caller-provided dimensions and values.
It does not fail loudly for an invalid `top_lev`, and it applies positivity
clipping after each explicit update. There is no conservation correction after
clipping. Extremely large inverse time scales can produce a large `nsubmix`
because no explicit cap is present.

## Reachability Assessment

The first two subsections record source-level reachability separately from
synthetic execution. The later subsections add the executed boundary and ledger
evidence; final classifications remain in the validation section.

### `top_lev` states

| State | Public interface | Current production path | Standalone test | Documented invalid? |
|---|---|---|---|---|
| `top_lev=0` | Accepted | Yes: default 128-level configuration | Yes | No |
| `top_lev=1` | Accepted | Configurable; no default identified | Yes | No |
| `top_lev=6` | Accepted | Yes: default 72-level configuration | Yes; existing validation uses it | No |
| `top_lev=pver-1` | Accepted | Configurable; no default identified | Yes | No |
| `top_lev<0` | Accepted as an `int`; not checked | No valid production configuration identified | Yes | No enforced contract |
| `top_lev>=pver` | Accepted as an `int`; not checked | No valid production configuration identified | Yes | No enforced contract |

### Disabled-mixing state

`enable_aero_vertical_mix=false` is still accepted by the EAMxx parameter
interface and is reachable through configuration. At the current MAM4xx
revision, however, the argument is unused and both Boolean values execute the
same always-mixing path. The deterministic enabled/disabled comparison below
applies only to the older evidence revision, where `false` reached the
defective partial-disable path.

### Historical executed boundary evidence

The current diagnostic branch contains the PR #567 index repair. The results
in this subsection deliberately reproduce the older, affected revision and are
not current-main test results.

On 2026-07-14, an external 128-level evidence driver compiled the inspected
`update_from_explmix` header without modifying MAM4xx production source. It used
Kokkos 4.5.1 Serial with `KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK`, AddressSanitizer,
and UndefinedBehaviorSanitizer. Each boundary value ran in a separate process.

| `top_lev` | Exit/result | Interpretation |
|---:|---|---|
| `-1` | Exit 134: Kokkos view access at index `-1`, extent 128 | Invalid input is not rejected before entering the kernel |
| `0` | Exit 134: Kokkos view access at index `-1`, extent 128 | Executed bounds violation for the default L128 production input |
| `1` | Exit 0; all 25 controlled ledgers within tolerance | Next boundary runs when the upper-interface coefficient is explicitly initialized |
| `6` | Exit 0; all 25 controlled ledgers within tolerance | Legacy-compatible active-level boundary runs under explicit zero-flux initialization |
| `127` | Exit 0; no bounds or ledger violation | Last physical level is accepted as the first active level |
| `128` | Exit 0 with no active-level work | Out-of-range upper input is silently accepted rather than rejected |

The final attached driver skips its replicated host-side substep calculation
for `top_lev<=0` and emits a launch marker immediately before calling
`update_from_explmix`. The marker precedes the `top_lev=0` Kokkos error, so the
recorded bounds violation originates in the exact MAM4xx function.

For `top_lev=1` and `6`, changing only `eddy_diff(top_lev-1)` from 0 to
2 s-1-equivalent units changed the reduction result and final state
fingerprint. This confirms that an in-bounds value above the active range is
consumed. `dropmixnuc` does not initialize that value for these boundaries.

The stock `mam4_ndrop_unit_tests` also passed 51 assertions in five test cases.
Those tests do not exercise `update_from_explmix` boundaries or paired column
ledgers, so they neither contradict nor establish the boundary result.

#### Why ordinary EAMxx runs can complete

The invalid index is evaluated whenever the L128 MAM4xx ACI path reaches the
`k=0` iteration, but an invalid access and a detected runtime error are not the
same event. Production GPU builds normally do not enable
`KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK`, so `View::operator()` does not reject the
negative index before performing the load. GPU virtual-memory protection is at
a much coarser granularity than one array element, and an adjacent address can
remain mapped even though it is outside the Kokkos view. AddressSanitizer and
the executed Kokkos bounds-checked host test intentionally make this otherwise
silent error fail loudly.

The current EAMxx allocation pattern can also mask the numerical symptom. The
`eddy_diff` and `csbot` scratch arrays are separate row-major
`(ncol,nlev)` Kokkos views, and `dropmixnuc` receives a row subview for each
column. In a typical unchecked address calculation, level `-1` for column
`c>0` aliases level `nlev-1` of column `c-1`. The code explicitly sets the
surface value `eddy_diff(nlev-1)=0`, so the erroneous product often evaluates
to zero, accidentally resembling the intended upper zero-flux boundary. This
is not a valid boundary implementation: column zero still addresses memory
before the allocation, cross-column initialization is not synchronized, and
allocator, backend, or compiler changes can alter the value or produce a
device fault.

Therefore the expected production symptom is not necessarily a crash. It can
be no visible effect, nondeterministic or configuration-dependent changes in
`eddy_diff_km`, `max_tinv`, and `nsubmix`, or an illegal-memory failure. The
negative-index expression is exercised by default L128 MAM4xx ACI runs; it is
not exercised by every EAMxx configuration, such as L72 with `top_lev=6` or a
configuration without MAM4xx ACI.

This paragraph describes the older evidence revision. Current MAM4xx main
contains PR #567 (`0b100fb3d`), which uses the already clipped `km1` index for
all factors in `eddy_diff_km`. The diagnostic branch described above includes
that repair and does not execute this negative-index expression.

#### Propagation and anticipated impacts

The out-of-bounds values are read only while constructing
`eddy_diff_km(0)`. The effective column maximum inverse timescale and
`nsubmix` are then computed from that coefficient; they are not themselves
read out of bounds. Because `nsubmix` and `dtmix` are column scalars, a corrupt
top-level coefficient can change the explicit timestep used at every active
level in that column.

For a finite `eddy_diff_km(0)`, the direct upper-boundary transport term is
normally algebraically zero. At `k=0`, the code sets `km1=0` and computes
`overlapm(0)=1`, so the upper-neighbor differences passed to both `explmix`
overloads vanish. The most likely finite-value effects are therefore:

- no effect when the aliased value is zero or when the resulting top-level
  inverse timescale does not control the column reduction;
- unnecessary substeps and approximately proportional additional loop cost
  when a spuriously positive value makes level zero limiting; and
- too few substeps when a negative value hides the true level-zero loss rate,
  with possible forward-Euler instability, positivity clipping, artificial
  sources, and aerosol or droplet-number errors.

Even when extra substeps preserve stability, changing `nsubmix` changes the
finite-step numerical solution throughout the column because all 51 explicit
updates use the resulting `dtmix`. If either invalid read produces a NaN or
infinity, the zero boundary difference does not protect the update: NaN times
zero remains NaN, while a very large or infinite inverse timescale can also
make the integer `nsubmix` calculation pathological. The unchecked read can
therefore cause NaNs, an excessive loop count, an illegal-memory failure, or
dependence on memory allocation, column ordering, backend, and execution
timing. These are anticipated propagation modes from the equations and memory
layout; their production frequency and magnitude have not yet been measured.

### Historical executed disabled-mixing evidence

Current MAM4xx has removed the operational disabled path. The same instrumented
128-level driver exercised the exact older
`update_from_explmix` implementation at `top_lev=1`, isolating switch behavior
from the separately confirmed `top_lev=0` defect. The controlled column used:

- `dtmicro=1 s`, unit `zn`, `zs`, and `csbot`, cloud fraction one, and explicit
  zero-flux upper and lower boundaries;
- nonzero, level-varying interstitial and cloud-borne values for all 25 aerosol
  constituents;
- nonzero number and mass activation rates below their rate limit;
- interior mixing coefficients selected to produce `nsubmix=1`, `2`, and `3`;
- destination buffer one initialized equal to the physical initial
  interstitial state for the primary ledger comparison; and
- no positivity clipping in the primary comparison.

The pressure-weighted ledger used the approved invariant
`sum_k (pdel/g) (I+C)`. With unit `zn=g/pdel`, its test weights are one.

| `nsubmix` | Switch | Maximum relative ledger error | Failed constituent ledgers | Interstitial L1 change | Cloud-borne L1 change | Droplet L1 change | Diagnosed droplet-source L1 |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | enabled | `1.64e-16` | 0/25 | `4.324257` | `4.324357` | `1.180064` | `1.179864` |
| 1 | disabled | `1.65e-3` | 25/25 | `0` | `4.324357` | `0` | `1.179864` |
| 2 | enabled | `2.93e-16` | 0/25 | `212.79845` | `212.79845` | `57.5294` | `56.0656` |
| 2 | disabled | `8.27e-2` | 25/25 | `0` | `216.21285` | `0` | `58.9932` |
| 3 | enabled | `3.12e-16` | 0/25 | `414.502575` | `414.502575` | `110.351675` | `102.889159` |
| 3 | disabled | `1.65e-1` | 25/25 | `0` | `432.4257` | `0` | `117.9864` |

All primary-case post-state values were finite and positive. Thus the disabled
ledger errors are not caused by positivity clipping. The disabled cases
diagnose a positive droplet source while leaving droplet number unchanged, and
change cloud-borne aerosol while leaving its interstitial pair unchanged.

A second matrix varied only interstitial destination-buffer initialization by
scales 0, 1, and 7. Enabled results were bit-for-bit identical at each substep
count. Disabled state fingerprints changed substantially at all three substep
counts; for `nsubmix=1` they ranged from `30656.8` to `386305.3`, and for
`nsubmix=3` from `41992.4` to `437920.0`. This establishes a separate
unwritten-buffer dependency. The scale-one primary matrix is sufficient to
establish the mass/number conservation violation without relying on poisoned or
uninitialized data.

## Verification and Validation Plan

Validation is STRICT because the result will guide production atmospheric
physics development. The reconstructed design and the test contract that
disables turbulent transport while retaining paired activation transfers were
approved before the investigation-specific tests were run.

### Boundary-index matrix

Exercise `top_lev` at:

- `0`, matching the default 128-level production configuration;
- `1`, the next boundary value;
- `6`, matching the 72-level production default and existing validation;
- `pver-1`, the upper accepted boundary implied by current loop construction;
- `-1` and `pver`, to characterize current fail-loud behavior for values that
  should be invalid under a reviewed contract.

For each case:

1. run a bounds-checked or sanitizer-enabled host configuration where locally
   available;
2. initialize all work arrays with deterministic values and poison or sentinel
   values outside the intended initialized range;
3. verify zero negative or out-of-range accesses;
4. verify that each intended active level is updated exactly once and inactive
   levels remain unchanged;
5. verify an explicit zero-flux top boundary; and
6. verify that invalid inputs fail deterministically before entering a device
   kernel once an input contract is approved.

Acceptance criteria are zero bounds violations, zero use of poisoned boundary
data, exact active-level coverage, and deterministic rejection of inputs
outside the approved range.

### Historical disabled vertical-mixing matrix

This matrix is expected-failure evidence for the older defective revision; it
is not a successful validation case and does not apply to current MAM4xx, where
the Boolean argument is unused. Use deterministic, nonzero interstitial and
cloud-borne number and mass fields.
Choose cloud fractions and activation rates that exercise paired transfers but
avoid positivity clipping in the principal conservation case. Compare enabled
and disabled runs while varying the mixing-rate profile enough to force one,
even, and odd `nsubmix` values. The executed 128-level evidence harness used
unit geometry, a one-second process interval, and representative interior
mixing coefficients of 0.01, 0.5, and 1.0 s-1-equivalent values, respectively.
It reported the value implied by the actual Kokkos reduction rather than an
independently assumed `dtmicro` cap.

For every aerosol constituent, calculate:

\[
L_j^{\mathrm{pre}}=\sum_k\frac{\Delta p_k}{g}(I_{j,k}^{\mathrm{pre}}+C_{j,k}^{\mathrm{pre}}),
\quad
L_j^{\mathrm{post}}=\sum_k\frac{\Delta p_k}{g}(I_{j,k}^{\mathrm{post}}+C_{j,k}^{\mathrm{post}}).
\]

Check:

- conservation of each mass and number ledger;
- paired-field source cancellation;
- nonnegativity and finite values;
- whether one member changes while its paired member does not;
- independence from initial contents of the destination work array;
- consistency across even and odd `nsubmix`; and
- agreement with the historical test definition of
  `enable_aero_vertical_mix=false`.

For a controlled no-clipping case, require

\[
|L_j^{\mathrm{post}}-L_j^{\mathrm{pre}}|
\le 100\epsilon\max\left(|L_j^{\mathrm{pre}}|,1\right),
\]

where \(\epsilon\) is the machine epsilon of the configured `Real` type and
the ledger is expressed in consistently scaled test units. Also require no
negative value below \(-10\epsilon\) times the field scale and no NaN or
infinity. A separate limiter-activation case should document, but not conceal,
any non-conservation introduced by positivity clipping.

### Production-path evidence

A synthetic kernel test can establish an interface hazard or an executed
bounds/ledger violation, but production reachability must be established
separately. For the older `top_lev=0` defect, retain the source trace from the
default L128 namelist through `MAMAci::run_impl` into `dropmixnuc`. For the
historical disabled switch, retain a source trace showing that a legal EAMxx
parameter list could supply `false`; explicitly distinguish this from current
MAM4xx, where that value no longer changes the executed path.

### Classification rule

Each concern will receive exactly one of these classifications after tests:

- **CONFIRMED DEFECT** only for a valid reachable call path with an executed
  invariant or bounds violation;
- **CONFIRMED SAFE** only when a demonstrated contract, guard, loop bound, or
  executed test prevents the suspected behavior;
- **INTERFACE HAZARD** for an allowed public input that is unsafe even if the
  production defaults do not exercise it;
- **CODE-READING RISK** when suspicious code is present but reachability or
  consequences are not established; or
- **UNRESOLVED** when local evidence cannot distinguish the possibilities.

Final classifications for the two requested concerns are:

- Boundary indexing is **CONFIRMED DEFECT**: `top_lev=0` is supplied by the
  default L128 production configuration and the exact 128-level C++ function
  produced an executed Kokkos bounds violation at index `-1`.
- Disabled vertical mixing is **CONFIRMED DEFECT**: `false` is accepted by the
  EAMxx parameter interface and forwarded to the production function; the
  executed one/even/odd substep matrix violated all 25 paired ledgers while the
  enabled control conserved them, and it separately demonstrated dependence on
  the unwritten destination buffer.

## Proposed Isolated Time-Integrator Study

> **[STANDALONE NUMERICAL EXPERIMENT - APPROVED 2026-07-23]**

This first phase is a simple standalone screening study that can run locally on
the development Mac. It will compare the current adaptive Forward Euler
implementation with an IMEX method and Backward Euler for the
`update_from_explmix` state advance and supporting coefficient assembly. It
does not change production source, activation physics, cloud-fraction
evolution outside this routine, or EAMxx coupling.

The purpose is initial evidence: verify the methods, establish first-order
behavior, measure positivity clipping and its conservation error, and obtain
CPU timing and operation-count estimates. This phase will not claim GPU or
whole-model performance. GPU implementation and profiling are deferred until
the local results justify that work.

The original study baseline was E3SM revision
`e12e7d61dba924717e8f35fe8288247c1766ee30` with MAM4xx revision
`1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`. That stack and its later
remediation evidence remain historical. The current locked stack and expanded
screening contract are defined below. Generated dependency checkouts and builds
must use task-specific paths below `/Users/burr114/e3sm-scratch/`, while
retained results and provenance must use task-specific paths below
`/Users/burr114/e3sm-results/`. No build or run tree may be created in the
OneDrive source tree.

### Cross-method screening execution contract

> **[CURRENT CROSS-METHOD SCREENING AMENDMENT - APPROVED 2026-08-18]**

This amendment authorizes only standalone test, validation, analysis, and
local CPU timing code. It does not authorize a change to production EAMxx or
MAM4xx equations, interfaces, gitlinks, or runtime behavior. Forward Euler
(FE), IMEX Euler, and Backward Euler (BE) receive independent numerical gates;
a method that fails its gate is retained as evidence but excluded from timing
and cross-method conclusions.

The current comparison stack is locked to these exact revisions:

| Repository | Current screening revision |
|---|---|
| E3SM source baseline | `e2c578e40f6602a108eff577091de86a453624ab` |
| MAM4xx | `0d076d8ebb004a18f05f51f97b08edd9def6826f` |
| EKAT | `3a1485aecfed7a036e4fa30ce845498fd37975e2` |
| Kokkos | `23745d8494203711a6090d0a201c3ae3bf57853b` |

The archived 2026-08-06 gate stack remains E3SM
`3f8d0812602730b7bcca2b02e5a4ecae42fd62d8`, MAM4xx
`1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`, EKAT
`65cda3ff61abe22add622b56cee6234e29fc51c4`, Kokkos
`6a19eb9454da490dbf785443d8f021028baeb161`, and Haero
`8f9b61aed4252461d76eb7c206927e4a005930f5`. It is historical evidence, not a
substitute for the current gate. The checked-in `dependency_lock.json` is the
machine-readable authority for both stacks. The preparation script must fail
closed on any missing commit, gitlink mismatch, attached HEAD, or dirty
dependency checkout and must record the source worktree separately from the
clean detached dependency checkouts.

#### First-active-level controller versus physical boundary

Let `t` be the first active level. The physical operator has zero flux through
the top interface and therefore uses

\[
E_t=0,\qquad D_t^{\uparrow}=0
\]

in FE tendencies and every IMEX or BE matrix. Current repaired MAM4xx may also
form

\[
D_{top,controller}(t)=
z_n(t)\,eddy\_diff(t)\,csbot(t)\,z_s(t)
\]

for the FE stability reduction at the first active level. This is a timestep
controller, not a physical exchange coefficient. It may increase `nsubmix`
but must never enter a tendency, ledger, IMEX matrix, BE matrix, or implicit
residual. The current-MAM4xx adapter must compare the returned `nsubmix`
directly for `t=0,1,6`, including a profile for which this controller alone
sets the stable FE step. Legacy MAM4xx `1b02f272...` at `t=0` remains an
expected failure and cannot satisfy the current gate.

#### Signed no-flux operator and implicit rows

Number active levels from `0` through `n-1` after removing levels above `t`.
Let `w_k=\Delta p_k/g`, let nonnegative `E_k` be the shared exchange at the
interface above level `k`, and define

\[
d_k^{\uparrow}=E_k/w_k,\qquad
d_k^{\downarrow}=E_{k+1}/w_k,
\]

with `E_0=E_n=0`. Let `o_k^{up}` and `o_k^{down}` be the cloud-overlap factors
toward levels `k-1` and `k+1`. Thus
`d_0^{up}=d_{n-1}^{down}=0`; no clamped neighbor is a physical boundary flux.
The signed scalar cloud/droplet matrix `A_C=I-h L_C` has rows

\[
(A_C)_{k,k}=1+h(d_k^{\uparrow}+d_k^{\downarrow}),
\]

\[
(A_C)_{k,k-1}=-h d_k^{\uparrow}o_k^{up},\qquad
(A_C)_{k,k+1}=-h d_k^{\downarrow}o_k^{down}.
\]

The scalar interstitial matrix `A_I=I-h L_I` has the same diagonal and

\[
(A_I)_{k,k-1}=-h d_k^{\uparrow},\qquad
(A_I)_{k,k+1}=-h d_k^{\downarrow}.
\]

Missing off-diagonal entries in the top and bottom rows are exactly zero. The
IMEX interstitial right-hand side receives the positive return terms

\[
h d_k^{\uparrow}(1-o_k^{up})C_{k-1}^{n+1}
+h d_k^{\downarrow}(1-o_k^{down})C_{k+1}^{n+1}.
\]

For BE, order each level block as `x_k=(C_k,I_k)^T`, write
`lambda_k=1+h(d_k^{up}+d_k^{down})`, and, away from the bottom source branch,
use

\[
L_k=-h d_k^{up}
\begin{bmatrix}o_k^{up}&0\\1-o_k^{up}&1\end{bmatrix},
\qquad B_k=\lambda_k I_2,
\]

\[
U_k=
\begin{bmatrix}
-h d_k^{down}o_k^{down}&-h a_k\\
-h d_k^{down}(1-o_k^{down})&-h d_k^{down}+h a_k
\end{bmatrix}.
\]

Here the approved non-bottom activation coupling is exactly
`S_k=a_k I_(k+1)`. Consequently `L_0=0`. At a no-flux bottom `U_(n-1)=0`
except that the retained current-C++ bottom source
`S_(n-1)=a_(n-1)(I_(n-1)+C_(n-1))` is incorporated into

\[
B_{n-1}=
\begin{bmatrix}
\lambda-ha&-ha\\ha&\lambda+ha
\end{bmatrix}.
\]

This bottom row is valid only while the current `max(0,tmpa)` branch remains
the same. If a candidate state would change that branch, the BE step must fail
with an outside-linear-branch diagnostic; it must not silently solve a
different equation. Number and mass use distinct `a_k` matrix families when
their activation rates differ. After the aerosol blocks, the droplet scalar
system uses `A_C` with the current droplet source independently assembled from
the solved modal-number state. Its source signs and bottom contribution must
be checked independently of the aerosol block bands.

The host-only BE implementation uses a 2-by-2 block Thomas factorization and
solve. A zero or unsafe scalar or 2-by-2 pivot fails deterministically when its
magnitude is at most

\[
100\,\epsilon\,\lVert A\rVert_{\infty}.
\]

The dense partial-pivoting oracle must assemble the physical equations
directly into a `2*n` by `2*n` matrix. It may not reuse block bands, BE
assembly helpers, or the block Thomas factorization. The raw BE solution is
checked for positivity before the common zero-cloud projection; the projection
is then tested separately for isolated, contiguous, and all-zero masks.

#### Matrix-family reuse

Factorization reuse is allowed once per distinct matrix family per column.
IMEX has one cloud/droplet scalar family and one interstitial scalar family.
BE has separate number and mass 2-by-2 block families per activation mode when
their coefficients differ, plus the droplet scalar family. Identical matrices
may share a factorization across right-hand sides, but no matrix, factorization,
or right-hand side may be reused across columns. Repeated matched-error steps
may reuse only within the same column, family, timestep, coefficient set, and
bottom-branch state.

#### Diagnostic equations and classifications

For each aerosol pair `j`, every method reports

\[
\Delta L_{j,total}=L_{j,post}-L_{j,pre},
\]

\[
\Delta L_{j,clip}=\sum_s\sum_k w_k
[(I_{j,k}^{clip}+C_{j,k}^{clip})-(I_{j,k}^{raw}+C_{j,k}^{raw})],
\]

\[
\Delta L_{j,operator}=\Delta L_{j,total}-\Delta L_{j,clip}.
\]

IMEX and BE set `Delta L_clip=0` because clipping is prohibited. A zero initial
ledger uses an absolute, not relative, conservation tolerance. FE additionally
forms at every substep and for every field the physical residual

\[
r_s=q_{raw}^{s+1}-q_{clip}^{s}-\Delta t_m F(q_{clip}^{s}).
\]

IMEX and BE form

\[
r=q^{n+1}-q^n-h(F_{explicit}+F_{implicit})
\]

directly from the physical interface coefficients and activation equations.
Residual normalizations may contain only like-unit terms from that field; mass,
aerosol number, and droplet number are never combined. All diagnostic paths
are exercised after multiplying state magnitudes by `1e-12`, `1`, and `1e12`.
Any nonfinite state, matrix, right-hand side, solution, ledger, residual, or
normalization fails closed.

All raw negative values are counted. A negative is material when it is below
`-10*epsilon*field_scale`, where the scale is the larger of the field's maximum
initial and maximum raw absolute value. Roundoff negatives are reported but do
not fail a method. Any required admissible IMEX or BE case with a material
negative fails that method. FE retains production per-substep clipping and
reports its frequency, magnitude, and ledger change. Any optional implicit
fallback clipping is a separately labeled experiment and cannot repair a
method gate.

#### Task 2 as-built synchronization

The standalone FE/IMEX requalification was implemented at E3SM harness commits
`e67097b6bc29e919329a729d48776f866a32038e` and
`07c786d82d7571fd01527b26cc3c21de09e09acd`. `OperatorData` and
`PreparedOperator` now carry `first_level_controller_rate` in s^-1. The value is
validated as finite and nonnegative, participates only in the FE
`production_nsubmix` reduction and its reported maximum inverse timescale, and
is absent from `apply_cloud_operator`, `apply_interstitial_operator`, FE
tendencies, IMEX matrices, residuals, and all aerosol ledgers.

`AdvanceDiagnostics` now exposes signed `total_ledger_delta`,
`clipping_ledger_delta`, and `operator_ledger_delta` arrays for all 25 aerosol
pairs, global and field-specific roundoff/material negative counts, and the
controller rate. FE assembles an independent like-unit physical residual for
each raw field at every substep before clipping. IMEX retains separate matrix
and independently assembled physical residuals. Both paths fail closed on
nonfinite residual or ledger arithmetic; zero initial aerosol ledgers use the
approved absolute normalization. No equation, source sign, boundary row,
limiter, bottom-source convention, or zero-cloud projection differs from the
approved contract.

The current MAM4xx adapter now builds against the locked rank-3-view API with
no archived Haero dependency. It directly compares the returned `nsubmix`,
final alternating buffer, droplets, and all aerosol reservoirs at
`top_lev=0,1,6`. In its controller-limiting profile, exact current MAM4xx and
the canonical harness both return two substeps while the same physical
no-flux operator without the controller returns one. The archived
`1b02f272...` failure remains a separate historical artifact and is not
registered in the current CTest gate.

#### Task 3 as-built synchronization

The standalone Backward Euler gate was implemented at E3SM harness commits
`2f43c8c1e4ebd6b97152414ab02dca760bc883e6` and
`c6f7f18a5d58e118cf4146a790fe69895a1c4adb`. It orders each aerosol level as
`(C,I)`, assembles the approved signed `L`, `B`, and `U` blocks, factors each
distinct activation matrix family once with a host 2-by-2 block Thomas solve,
and then advances droplets with one scalar cloud matrix using sources from the
solved modal-number state. The first-level FE controller is reported but does
not enter any BE matrix, source, residual, or ledger.

Unsafe scalar pivots fail at `100*epsilon*||A||_infinity`. The solver reports
matrix and independently assembled physical residuals, per-field raw minima,
roundoff/material negative classifications, exact factorization and right-hand-
side counts, and signed total/clipping/operator ledger changes. Positivity is
checked before the common zero-cloud projection; clipping is never applied.
Nonfinite inputs, matrices, right-hand sides, solutions, branch diagnostics,
residuals, and ledger arithmetic fail at the detected stage.

The committed dense oracle independently assembles a `2*n` matrix directly
from the physical equations and solves it with scalar partial pivoting. It does
not reuse BE bands, factorization, source-mapping helpers, or the block solver.
The canonical no-flux limiter makes the admissible bottom activation rate zero,
so a test-only inspection separately verifies the retained nonzero bottom block
`[[lambda-ha,-ha],[ha,lambda+ha]]` and the positive bottom droplet source.
Changing that source branch remains a deterministic failure. Debug/ASan/UBSan
and optimized gates pass; no production source, interface, or gitlink changed.

#### Task 4 as-built synchronization

The common cross-method gate was implemented at E3SM harness commit
`7807aedfdecafe5e30b783e98e81ac4b723f56e3`. Its adaptive Dormand--Prince
5(4) oracle independently reconstructs the complete physical right-hand side,
including shared-interface rates, overlap, activation limiting, separate
number/mass maps, the retained current-C++ bottom sources, and all 51 state
fields. Analytic two-level diffusion and activation oracles and an `rtol`
`1e-12` versus `1e-13` reference check independently gate the reference.

The exact Cartesian matrix contains 1,875 cases and 5,625 method records. Each
method performs three bitwise-identical serial advances; every record retains
all-field weighted-L1/Linf errors, physical residuals, raw positivity,
clipping, and signed total/operator/clipping aerosol ledgers. Error norms and
their reductions fail closed on nonfinite inputs or intermediate arithmetic,
and the no-material-negative rule applies to FE, IMEX, and BE alike.

Clean debug/ASan/UBSan and optimized Kokkos Serial runs each pass 11,880 checks
with identical machine-readable results. All three methods pass all 1,875
cases, all-field first-order slopes and `R^2`, and pairwise Richardson
common-solution bounds. The maximum reference two-tolerance change is
`2.0331196e-13`; maximum residual and ledger errors are respectively below
`4.1e-16` and `2.5e-15`, with no material negative or clipping case. This
authorizes Task 5 local CPU timing only. Coarse one-call Linf errors reach
`0.1129882`, so Task 4 does not establish full-interval production accuracy,
matched-error cost, GPU performance, or a production recommendation.

#### Required numerical matrix and gates

The accuracy gate covers level counts `1,2,3,127,128`; overlap values
`0,0.25,0.5,0.75,1`; activation strengths `0,0.1,0.5,0.9` times the limiter
and limiter-active strength `2.0`; active aerosol-pair counts `1,4,8,16,25`
(`3,9,17,33,51` transported fields); zero, typical, and strong mixing; and
deterministic seed `0x4E44524F`.

Analytic exact solutions are used where available. Otherwise, an independently
assembled adaptive Dormand-Prince 5(4) oracle evaluates the physical right-hand
side with `rtol=1e-12` and field-scaled absolute tolerances. Repeating the
reference with tolerances ten times tighter must change every field by no more
than `2e-11`. Each method must have observed first-order slope in
`[0.85,1.15]`, `R^2>=0.98`, double-precision residual at most `1e-11`, ledger
error at most `max(1e-11,5*nlev*epsilon)`, no material negative, and three
serial runs with bitwise-identical outputs. Only methods passing their complete
gate enter performance comparisons.

The Task 4 oracle independently reconstructs pressure-normalized interface
rates, overlap ratios, activation limiting, constituent-to-mode maps, the
current-C++ bottom source, and all 51 physical right-hand sides directly from
`OperatorData`; it does not call `prepare_operator`, an integrator matrix or
band assembler, or an integrator tendency/source helper. It uses the standard
Dormand-Prince 5(4) embedded coefficients and an infinity norm of
`|q_5-q_4|/(atol_j+rtol*max(|q_n|,|q_5|))`. For each transported field `j`,
`atol_j=rtol*max(max_k|q_j(0)|,1e-30)`. Accepted-step growth is limited to a
factor of five, rejected-step shrinkage to at least a factor of 0.2, and both
use a 0.9 safety factor. A run fails on a nonfinite stage, error norm, state,
step underflow, or more than 100,000 attempted steps. The canonical zero-cloud
state transition is applied once after the continuous reference advance, so a
one-call candidate comparison has the same projection boundary; smooth order
fits use strictly positive clouds and repeated candidate steps without a
projection discontinuity.

The deterministic common matrix uses ACI interval `0.5 s`; zero, typical, and
strong interface-rate scales `0`, `0.05`, and `0.5 s^-1`; pressure weights in
`[0.8,1.2] kg m^-2`; alternating cloud fractions `0.8` and `0.8*r` for the
declared overlap ratio `r`; and the approved seed `0x4E44524F` for bounded
state/profile perturbations. Number activation uses mode multipliers
`0.7,0.8,0.9,1.0`; mass activation uses `0.55,0.65,0.75,0.85`; each multiplier
is applied to the declared activation strength times the local below rate
before canonical limiting. Inactive aerosol pairs are exactly zero. The
controller-only first-level rate is one quarter of the selected mixing scale
and enters only FE substep selection. These dimensionless-rate products keep
the common matrix within the method-local positivity domain while still
exercising the limiter-active category. Task 4 writes one deterministic CSV
record per case and method plus a Markdown summary; retained outputs belong
under `/Users/burr114/e3sm-results/ndrop-screening/task4-cross-scheme`.

#### Required local CPU screen

Performance cases cover levels `1,2,3,16,32,64,72,127,128`; active aerosol
pairs `1,4,8,16,25`; column batches `1,16,64,256,1024`; typical and strong
mixing; ACI intervals `33.33,50,75,100,150,300` seconds; and FE `nsubmix`
`1,2,3,5,8,16,21`. The timed full advance includes assembly, limiter and
controller logic, matrix and right-hand-side work, copies, factorizations,
solves, positivity checks, FE clipping, state update, and required fences.
Separable subtimers and matched diagnostics-off/on runs report instrumentation
overhead.

Each case uses 10 warmups and 30 through 100 measured repetitions. Sampling
stops only when a fixed-seed 10,000-resample bootstrap median 95% confidence
interval has less than 5% half-width or 100 repetitions are reached. Methods
are interleaved and output checksums are consumed. Runtime-versus-error targets
are `1e-2,3e-3,1e-3,3e-4,1e-4`; extrapolation is prohibited. Break-even is the
smallest measured FE `nsubmix` whose median runtime is at least the candidate
median, with interpolation only between adjacent measured FE values and with
bootstrap uncertainty carried through. Points with FE clipping or implicit
material negativity are excluded. Hybrid results are threshold and packing-
overhead estimates only, and local Mac measurements cannot support a GPU
speedup claim.

### Numerical contract: operator and integrator

The **operator** is the fixed set of vertical mixing coefficients, cloud-
overlap factors, activation sources, and explicit no-flux boundaries embodied
in the approved explicit equations above. The **integrator** is the method used
to advance that operator through time; changing an integrator in the isolated
comparison must not silently change the operator.

The canonical current-MAM4xx **state transition** also includes a common
post-integration zero-cloud operation that is separate from the continuous
operator and from the choice of integrator. After the final Forward Euler
substep or implicit solve, every level whose current cloud fraction is exactly
zero must:

\[
N_{c,k}\leftarrow 0,\qquad
I_{j,k}\leftarrow I_{j,k}+C_{j,k},\qquad
C_{j,k}\leftarrow 0
\quad\text{for all }j.
\]

This operation reproduces `update_from_explmix` and preserves every aerosol
ledger \(L_j\). It must be applied identically after Forward Euler, IMEX, and
Backward Euler. Exact-current comparisons must include columns containing
isolated zero-cloud levels and contiguous zero-cloud regions. Analytic
convergence cases must use strictly positive cloud fractions so this
non-smooth projection is not mistaken for time-integration error; zero-cloud
behavior is validated separately as an exact state-transition and ledger
property. IMEX and Backward Euler positivity checks apply to the raw solved
state before this operation so resuspension or reset cannot hide a materially
negative candidate solution.

- **Forward Euler** treats vertical mixing, cloud-to-interstitial return, and
  activation explicitly at the saved substep state and uses the production
  adaptive subcycling.
- **IMEX** treats activation explicitly at the old state while treating
  vertical mixing and cloud-to-interstitial return implicitly.
- **Backward Euler** treats vertical mixing, cloud-to-interstitial return, and
  activation implicitly at the new state.

### Scope and canonical test operator

The primary comparison will reproduce the current C++ semi-discrete equations,
including behavior that may later be corrected independently. This isolates
the time integrator as the only intended algorithmic difference.

Let \(L_C\) be the cloud-overlap-weighted mixing operator used for cloud-borne
aerosol and droplet number, \(L_I\) the full-overlap interstitial mixing
operator, and \(R\) the return of cloud-borne material mixed into non-cloudy
air. For aerosol constituent \(j\), let \(S_j(I_j,C_j)\) be the current C++
activation transfer and let \(S_N(I,C)\) be the droplet-number source. With
mixing coefficients, overlap, and activation rates frozen over one ACI call,
the equations are

\[
\frac{dC_j}{dt}=L_C C_j+S_j(I_j,C_j),
\]

\[
\frac{dI_j}{dt}=L_I I_j+R C_j-S_j(I_j,C_j),
\]

\[
\frac{dN_c}{dt}=L_C N_c+S_N(I,C).
\]

For all non-bottom active levels, the current source has the form
\(S_{j,k}=a_{j,k}I_{j,k+1}\). Modal number uses `nact`; species mass uses
`mact`. At the bottom level, the current C++ implementation uses `nact` for
both number and mass and forms a source proportional to \(I_{j,k}+C_{j,k}\).
The inspected Fortran implementation instead uses `mact` for mass and includes
an additional cloud-borne loss term. The primary experiment must retain the
current C++ formula. A separate, non-gating `mact`-at-bottom sensitivity may be
reported, but it must not be described as complete Fortran equivalence.

The canonical test operator will receive one nonnegative exchange coefficient
per physical interface and explicit zero-flux top and bottom boundaries. For
pressure weight \(w_k=\Delta p_k/g\), the cell rates must be derived from the
shared interface coefficient so that

\[
w_kD_k^+=w_{k+1}D_{k+1}^-.
\]

The cloud-borne, interstitial, and cloud-to-interstitial return operators must
all be assembled from these same interface data and overlap factors. This
prevents a standalone test from accidentally constructing a nonconservative
operator. An input that supplies incompatible cell rates must fail validation.

The harness must not reproduce the historical `eddy_diff(-1)` access. Current
exact-function comparison uses the locked repaired MAM4xx revision at
`top_lev=0,1,6`; only the archived affected revision remains restricted to a
safe active range such as `top_lev=1`. Passing repaired-current tests must not
be reported as retroactive validation of the archived defective path.

### Test-only interface

No new production API is proposed. The isolated harness will expose the
conceptual interface

```text
advance(method, state, operator_data, dt) -> advance_result
```

where:

| Item | Contents and contract |
|---|---|
| `method` | `forward_euler`, `imex`, or `backward_euler` |
| `state` | `N_c[nlev]`, `I[25][nlev]`, and `C[25][nlev]`; finite and nonnegative on input |
| `operator_data` | Shared interface exchange coefficients, overlap, pressure weights, capped `nact`/`mact`, mode map, and explicit boundary conditions |
| `dt` | Positive ACI interval in seconds |
| `advance_result` | New state, status, executed substeps, limiter count, per-field raw minima and clipping diagnostics, normalized residual, and normalized ledger error |

Invalid dimensions, nonfinite values, negative mixing or activation rates,
incompatible interface rates, invalid boundaries, and nonpositive `dt` must
return a deterministic failure before a solve. Finite inputs whose derived
coefficients, sources, matrix entries, right-hand sides, solutions, weighted
ledgers, residual norms, or normalizations overflow or otherwise become
nonfinite must return a deterministic failure at the point of detection.
Success therefore means that the complete state and every reported diagnostic
are finite. The implementation must not use comparisons or maximum reductions
that silently convert a NaN diagnostic into an apparently acceptable zero.

### Positivity and clipping

In summary, Forward Euler reproduces production subcycling and clipping,
whereas test-only IMEX and Backward Euler do not clip and report material
negative values as positivity failures.

Nonnegative physical state is a strict requirement. Forward Euler must
reproduce the production subcycling and componentwise `max(0,qnew)` clipping
semantics exactly, including the existing activation-rate limiter and clipping
after every explicit substep. It must retain the raw pre-clipping result for
diagnosis.

IMEX and Backward Euler are test-only schemes and **must not silently clip**.
Negative values consistent with the declared roundoff tolerance are acceptable
and remain visible in diagnostics; larger negative values are reported as
positivity failures for that method and case.

For Forward Euler, record the raw minimum, count of clipped values, and
weighted negative deficit for every field. For aerosol constituent \(j\),
record the ledger change introduced only by clipping,

\[
\Delta L_{j,clip}=\sum_k w_k\left[(I_{j,k}^{clip}+C_{j,k}^{clip})
 -(I_{j,k}^{raw}+C_{j,k}^{raw})\right],
\]

as well as the total pre-to-post ledger change. When clipping is inactive,
activation plus no-flux mixing must preserve the pressure-weighted column
ledger \(L_j\). When production Forward Euler clipping is active, its
contribution is reported separately and is not classified as an operator or
paired-transfer defect. No downstream or global fixer is assumed in this
standalone phase. Whether EAMxx later masks or redistributes these errors is a
separate source-tracing and coupled-model TODO.

### Forward Euler reference

Forward Euler will retain the production calculation of `nsubmix`, use
\(\Delta t_m=\Delta t/n_{submix}\), recompute the state-dependent activation
source each substep, alternate old/new buffers, and apply componentwise
`max(0,qnew)` clipping. A separate unclipped fixed-step Forward Euler path may
be used only to measure convergence in smooth cases known to remain positive.

The canonical Forward Euler implementation must agree with the exact current
C++ function for safe boundaries, one, even, and odd substep counts, and both
limiter-active and limiter-inactive cases before either candidate method is
implemented.

### IMEX definition

IMEX evaluates activation from the old state and treats vertical mixing and
cloud-to-interstitial return implicitly:

\[
(I-\Delta t L_C)C_j^{n+1}=C_j^n+\Delta t S_j(I_j^n,C_j^n),
\]

\[
(I-\Delta t L_I)I_j^{n+1}=I_j^n-\Delta t S_j(I_j^n,C_j^n)
 +\Delta t R C_j^{n+1},
\]

\[
(I-\Delta t L_C)N_c^{n+1}=N_c^n+\Delta t S_N(I^n,C^n).
\]

The cloud-borne and droplet systems share one scalar tridiagonal matrix; the
interstitial fields share a second. The implementation must factor or reduce
each shared matrix once per column and apply it to all relevant right-hand
sides. Local execution will use EKAT Thomas. A later GPU phase may evaluate
EKAT cyclic reduction after the standalone screen is complete.

This IMEX method is not assumed to produce a nonnegative raw state because its
activation transfer is explicit. A material negative result is a reported
positivity failure; it is not replaced by a clipped test result.

### Backward Euler definition

Backward Euler treats mixing, cloud-to-interstitial return, and activation at
the new time level:

\[
C_j^{n+1}=C_j^n+\Delta t\left[L_C C_j^{n+1}
 +S_j(I_j^{n+1},C_j^{n+1})\right],
\]

\[
I_j^{n+1}=I_j^n+\Delta t\left[L_I I_j^{n+1}+R C_j^{n+1}
 -S_j(I_j^{n+1},C_j^{n+1})\right].
\]

Each aerosol pair therefore forms a 2-by-2 block-tridiagonal system. Modal
number fields share matrices by `nact`; species fields share matrices by mode
and `mact`, except for the retained current-C++ bottom source. After solving
the four aerosol-number pairs, droplet number is advanced with

\[
(I-\Delta t L_C)N_c^{n+1}=N_c^n+\Delta t S_N(I^{n+1},C^{n+1}).
\]

A host block-Thomas implementation is the local correctness and timing
baseline. It must agree with an independently coded dense partial-pivoting
solve for small systems. If the implicit bottom source changes the active
branch of the current `max(0,tmpa)` expression or produces a negative combined
`I+C` state, the case must be reported as outside the simple linear solve or as
a positivity failure rather than silently changing the equation. A block
cyclic-reduction GPU solver is deferred.

### Linearity and algebraic structure

The full activation parameterization is nonlinear when supersaturation,
particle size, activation fractions, diffusivity, or cloud geometry are
recomputed from the evolving model state. The isolated integrator problem in
this design is narrower: mixing coefficients, overlap, `nact`, `mact`, and all
geometry are frozen over an ACI call. Under that contract, the mixing operator
is linear and each activation source has a frozen coefficient multiplying the
tracer state. Neither the approved IMEX method nor the proposed Backward Euler
method requires a nonlinear iteration or Newton solve.

The complete IMEX aerosol update is not one scalar tridiagonal system. For each
interstitial/cloud-borne pair it has the block-lower-triangular form

\[
\begin{bmatrix} A_C & 0 \\ -\Delta t R & A_I \end{bmatrix}
\begin{bmatrix} C^{n+1} \\ I^{n+1} \end{bmatrix} =
\begin{bmatrix} b_C \\ b_I \end{bmatrix},
\qquad
A_C=I-\Delta tL_C,\quad A_I=I-\Delta tL_I.
\]

It is therefore solved in sequence: solve the cloud-borne and droplet scalar
tridiagonal systems, form the cloud-to-interstitial return from the new
cloud-borne state, and then solve the interstitial scalar tridiagonal systems.
This structure permits the existing scalar EKAT solver to be used with many
right-hand sides. Backward Euler couples activation between `C` and `I` at the
new time level, so each aerosol pair is block tridiagonal, but it remains a
linear block solve while the coefficients remain frozen.

The problem becomes genuinely nonlinear if a future method recomputes
activation rates, sizes, supersaturation, diffusivity, overlap, or branch
conditions from the unknown new state. Production Forward Euler clipping is
also nonlinear and nonsmooth. The candidate IMEX and Backward Euler methods do
not include that clipping; zero-cloud resuspension/reset is applied after the
advance for a fixed cloud mask. The bottom-level `max(0,tmpa)` branch remains
an explicitly reported excluded case whenever a candidate solve would change
its active branch.

### Linear-solver and integration-library boundary

EKAT is the **E3SM Kokkos Application Toolkit** in `externals/ekat`. Its
header-only `ekat_tridiag.hpp` supports diagonally dominant scalar
tridiagonal systems with one matrix and one or many right-hand sides, or many
independent matrices. It provides serial and team Thomas algorithms for CPU
execution and Kokkos-team cyclic reduction intended for GPU execution. The
arrays use `LayoutRight`, and the implementation is designed to solve many
independent atmospheric columns concurrently.

EKAT supplies the algebraic solve, not the `ndrop` integration method. An IMEX
implementation still needs a thin `ndrop`-specific layer to evaluate the
explicit activation source, assemble matrices and right-hand sides, solve the
cloud-borne fields before the interstitial fields, apply the common zero-cloud
state transition, and report positivity and ledger diagnostics. The preferred
first production prototype, if separately approved, is therefore a small
`ndrop` IMEX-Euler implementation using EKAT's multiple-right-hand-side
tridiagonal solver rather than a new tridiagonal solver.

General integration packages remain useful reference or later-stage options.
SUNDIALS ARKODE supplies additive Runge-Kutta IMEX methods, and PETSc supplies
`TSARKIMEX`, but either package would still require `ndrop`-specific explicit
and implicit callbacks, state adapters, matrix/Jacobian integration, the
zero-cloud transition, positivity checks, and GPU column batching. That is
disproportionate integration work for the current fixed-step, first-order
screen. SUNDIALS could be valuable as an independent standalone reference, for
higher-order or adaptive IMEX experiments, or if future physics makes the
implicit problem nonlinear. Before production use, it would need a dependency,
Kokkos-version, device-memory, and batched-performance gate; the SUNDIALS
installation documentation inspected for this review reports Kokkos testing
only through 4.3.01, whereas the local E3SM evidence used Kokkos 4.5.1. Neither
SUNDIALS nor PETSc was used by the current standalone harness.

### Other implicit candidates (not authorized)

Backward Euler remains the required fully implicit baseline because it is
simple, strongly damps stiff modes, and exposes the coupled block system
directly. If that baseline shows unacceptable first-order damping or accuracy,
the next method should be **one** two-stage, second-order, L-stable scheme:
SDIRK2/ESDIRK2 or TR-BDF2. Both can reuse the same basic matrix structure over
their stages and are more relevant than RK4 to this stiff problem, but they
roughly double the solve work and do not guarantee positivity.

| Candidate | Potential value | Critical limitation for `ndrop` | Priority |
|---|---|---|---|
| SDIRK2/ESDIRK2 | Second order and L-stable; repeated diagonal coefficient can permit matrix reuse | Two implicit stages; positivity is not guaranteed | First higher-order candidate |
| TR-BDF2 | Second order and L-stable with favorable stiff damping | Two-stage implementation and validation; positivity is not guaranteed | Comparable alternative to SDIRK2 |
| BDF2 | Second order with one solve after startup | Requires reliable process-local history and startup; other split physics modifies aerosol between ACI calls | Lower priority |
| Crank-Nicolson or implicit midpoint | Second order with one solve | A-stable but not L-stable; stiff modes can oscillate and produce negative values | Not recommended for the first screen |
| Rosenbrock-W/ROS2 | Linearly implicit and well matched to a frozen linear operator | Two solves and additional stage machinery; not a conventional fully implicit step | Secondary option |
| Higher-order Radau/implicit Runge-Kutta | High-order stiff accuracy | Multiple coupled stages are too costly for an initial 51-field column solve | Defer |
| Modified Patankar methods | Can target simultaneous positivity and conservation for production-destruction systems | Requires a new derivation for diffusion, overlap, and paired transfer; not a drop-in implicit solve | Separate research track |

Higher formal order inside `ndrop` may not increase whole-model order if
frozen coefficients or outer ACI/process splitting dominate the error. The
decision sequence is therefore: complete the common cross-method accuracy
matrix, compare BE damping and matched-error cost with FE and IMEX, and add
SDIRK2/ESDIRK2 or TR-BDF2 only if the evidence shows that a higher-order
implicit method is needed. RK4 remains outside the planned screening suite.

### Algebraic invariants

For pressure weights \(w_k=\Delta p_k/g=1/z_{n,k}\), coefficient assembly must
demonstrate the weighted cancellation

\[
\sum_k w_k\left[(L_C C)_k+(L_I I)_k+(R C)_k\right]=0.
\]

Because the same \(S_j\) enters the paired equations with opposite signs, each
raw, unclipped method must preserve

\[
\sum_k w_k(I_{j,k}+C_{j,k})
\]

for all 25 aerosol constituents under zero-flux boundaries. The proof must be
mapped to matrix coefficients and supplemented by randomized tests constructed
from compatible shared-interface coefficients. Production Forward Euler
clipping may break this invariant, but its contribution must be isolated with
\(\Delta L_{j,clip}\); IMEX and Backward Euler must report a positivity failure
instead of modifying the raw ledger through clipping.

### Validation and acceptance (high level)

Validation will exercise identity and zero operators, uniform-state
preservation, analytic no-flux diffusion modes, and simple activation-transfer
cases. Forward Euler, IMEX, and Backward Euler are each expected to show
first-order convergence under timestep refinement and to satisfy tight solver-
residual and column-ledger tolerances. Detailed cases and numerical thresholds
are defined in the separate validation plan.

These criteria apply equally to all three methods. A method that has not passed
its own validation gate cannot be used either to validate another method or to
reject another method on accuracy or performance grounds.

### Known defects and excluded cases

**The older `top_lev=0` bounds failure and partial disabled-mixing path are
excluded from successful numerical validation.**

At the evidence revisions, the production-reachable `top_lev=0` invalid
boundary access and the defective `enable_aero_vertical_mix=false` path were
known defects, not successful validation cases. Tests may reproduce them as
expected-failure evidence, but excluded paths must never be counted as passing
numerical validation. Current MAM4xx contains the PR #567 index repair and has
removed the operational partial-disable path; current validation must exercise
those repaired semantics.

These cases and their classifications are documented in
[`docs/analysis/eamxx_mam4xx_ndrop_defect_brief.md`](../analysis/eamxx_mam4xx_ndrop_defect_brief.md).
This design must remain consistent with that brief as either artifact is
updated.

### Validity domain and regimes

The operator and ledger definitions are intended to apply from small isolated
columns through large level counts including typical L128-like columns, across
cloud-overlap values from zero through fractional overlap to full overlap, and
from typical turbulence to strong-mixing and cloud-base regimes. Validation
must demonstrate that the equations, explicit no-flux boundaries, and ledger
accounting remain meaningful throughout these qualitative regimes; known
defective interfaces remain excluded as stated above.

### Staged implementation and review gates

Implementation may advance only through these gates:

1. **Design approval:** approve this section before any prototype code.
2. **Harness and Forward Euler gate:** establish exact-current-function
   agreement under safe boundaries and independently review signs, indexing,
   units, source selection, clipping, and buffer parity.
3. **IMEX gate:** complete code review and pass the entire numerical validation
   suite before beginning Backward Euler.
4. **Backward Euler gate:** validate block Thomas against the independent dense
   solve, review the coupled matrices, and pass the suite before cross-scheme
   comparison.
5. **Cross-scheme freeze:** all three schemes must pass from one revision in
   debug/sanitizer and optimized local configurations before CPU timing.
6. **Local screening report:** complete accuracy, clipping, conservation, and
   CPU timing comparisons before deciding whether to begin a GPU phase.

When a gate fails, subsequent implementation stops until the defect is fixed,
the design is updated if necessary, and the failed gate is rerun. Optimization
after a gate requires rerunning that scheme's complete validation suite.

The original local runs reported the Forward Euler and IMEX gates as passed on
2026-07-28. An independent review on 2026-07-29 reran both suites and confirmed
the reported normal-regime assertion counts and convergence results, but
assigned both gates a `FAIL`. Reviewer-created tests demonstrated that the
standalone state transition disagreed with current MAM4xx at zero-cloud levels
and that finite, extreme-magnitude inputs could yield nonfinite ledger or
residual diagnostics without a failure. The review also confirmed that the
Forward Euler artifact lacked its required analytic and coupled convergence
suite, the IMEX activation oracle did not exercise all aerosol fields, neither
gate implemented the required coupled all-field convergence case, and Forward
Euler clipping diagnostics were not fully per-field. The remediation suite
passed in debug/sanitizer and optimized local configurations on 2026-07-30,
including exact-current zero-cloud comparisons. A second independent review on
2026-08-06 found no numerical-method defect but retained `FAIL` for both gates
because the committed activation oracles were insensitive to mode
permutations, Forward Euler clipping and finite-stage regression coverage was
incomplete, and IMEX lacked a committed limiter-active oracle and focused
finite-checker tests. The corrected candidate includes those cases and passed
all three registered executables in sanitizer/debug and optimized builds. A
fresh independent extra-high review reran the suites, exercised mode and matrix
mutations plus independent hand calculations, and assigned strict `PASS`
decisions to both gates. The current locked-stack Backward Euler method-local
gate subsequently passed in debug/sanitizer and optimized configurations on
2026-08-24. The common Task 4 Cartesian gate then passed for FE, IMEX, and BE
in both configurations from one committed revision; the approved local CPU
timing phase may proceed. Executed evidence is recorded in
[`docs/analysis/evidence/ndrop_integrator/README.md`](../analysis/evidence/ndrop_integrator/README.md).

### Accuracy and convergence validation

Validation uses all 25 aerosol pairs and droplet number and includes:

- zero-operator identity and uniform-state preservation with full cloud
  overlap; fractional-overlap cases instead require the combined aerosol
  ledger invariant;
- exact-current zero-cloud resuspension/reset for isolated zero-cloud levels,
  contiguous zero-cloud regions, and zero-overlap profiles;
- an analytic discrete no-flux diffusion mode applied independently to all 51
  fields;
- analytic activation-transfer cases that exercise droplet number, all four
  modal-number pairs, and every species-mass pair through the appropriate
  mode-specific activation coefficient;
- a representative coupled self-convergence comparison in which all methods
  approach the same timestep-refined solution, with nonuniform mixing,
  fractional but strictly positive cloud overlap, activation in every mode,
  and nonzero state in all 51 fields;
- overlap values zero, fractional, and one;
- level counts 1, 2, 3, 127, and 128;
- documented typical and strong/cloud-base mixing profiles;
- one, even, and odd Forward Euler substep counts; and
- invalid inputs, solver failure, positivity failure, and limiter-active cases.

Finite-value failure tests must separately cover nonfinite input rejection and
overflow or contamination in derived coefficients, activation sources,
matrices, right-hand sides, solutions where constructible, ledgers, residuals,
and normalizations. If the approved M-matrix contract makes a nonfinite
solution impossible from finite matrix and right-hand-side values, a focused
test of the solution checker may be used instead of fabricating an invalid
physical operator. Every such case must return a deterministic non-success
status and a diagnostic identifying the failed stage. At least one IMEX
residual test must construct the residual directly from the shared physical
interface coefficients rather than reusing the tridiagonal coefficients
assembled for the solver.

RK4 is not a candidate method and is not required for this initial phase. The
analytic diffusion and activation cases establish formal order. The coupled
case uses three-level self-convergence,

\[
p=\log_2\left(\frac{\lVert q_h-q_{h/2}\rVert}
 {\lVert q_{h/2}-q_{h/4}\rVert}\right),
\]

plus agreement among the methods at the finest timestep. An independent
higher-order reference may be added later only if the analytic and
self-convergence evidence disagree or cannot identify the asymptotic range.

Use at least four timestep halvings and fit the three finest points that remain
above roundoff. In no-clipping cases, both normalized weighted-L1 and
L-infinity errors must demonstrate first-order convergence with fitted slope
in `[0.85,1.15]` and `R^2 >= 0.98`. Formal order is not required after clipping
activates because the non-smooth projection can change it; those cases instead
report raw error, clipped error, and ledger impact.

Errors must be normalized separately for every field; dimensional aerosol
mass, aerosol number, and droplet number must never be summed into one norm.
For field \(j\), define

\[
E_{1,j}=\frac{\sum_k w_k|q_{j,k}-q_{j,k}^{ref}|}
 {\max(\sum_k w_k|q_{j,k}^{ref}|,\sum_k w_k|q_{j,k}^{0}|)},
\]

\[
E_{\infty,j}=\frac{\max_k|q_{j,k}-q_{j,k}^{ref}|}
 {\max(\max_k|q_{j,k}^{ref}|,\max_k|q_{j,k}^{0}|)},
\]

Convergence cases must initialize every assessed field with a nonzero norm so
these denominators are defined. Tests whose exact state is identically zero
use the separately declared absolute roundoff tolerance rather than a relative
error. Report every field, category maxima, and the maximum across fields.
Pass/fail and matched-error screening use the maximum field-relative error.

For analytic cases, \(q^{ref}\) is the exact solution. For the coupled case,
do not designate one candidate method as truth. After first-order
self-convergence and finest-step cross-method agreement are established, use
the first-order Richardson estimate

\[
\widehat E_{h,j}=\lVert q_{h,j}-q_{h/2,j}\rVert/(2^1-1)
\]

as the screening error estimate for runtime-versus-error curves. Do not report
a matched-error comparison for a Forward Euler timestep at which clipping
occurs or for an IMEX or Backward Euler timestep that has a positivity failure.
For those cases, report the Forward Euler raw/clipped difference and ledger
change or the candidate raw negative state and failure status, respectively.

For each field in the finest no-clipping coupled case, the difference between
any two methods must decrease under refinement and must not exceed twice the
larger of their Richardson error estimates. Otherwise the methods have not
demonstrated convergence to a common solution and the cross-scheme gate fails.

Normalized linear-solver residuals must be at most `1e-11` in double precision
and `5e-6` in single precision. With `eps` equal to the configured `Real`
machine epsilon, normalized ledger error must not exceed
`max(ledger_floor, 5*nlev*eps)`, where `ledger_floor` is `1e-11` for double and
`1e-6` for single precision. A raw negative value is significant when its
magnitude exceeds `10*eps` times the larger of the field's maximum initial and
maximum raw absolute values. All raw negatives are counted in the positivity
diagnostics; Forward Euler additionally records the amount changed by clipping.
Forward Euler must retain separate raw minima, clipping counts, and weighted
deficits for droplet number and for each interstitial and cloud-borne aerosol
field. Clipping-active tests must verify those per-field values and the
resulting clipping-only ledger change. The suite must seek an admissible
production-adaptive clipping case rather than relying only on a forced unstable
step. If the approved rate and subcycling constraints make adaptive clipping
unreachable, the evidence must state and test that property and retain the
forced-step case only as a diagnostic-unit test.

Required configurations are a double-precision CPU debug/sanitizer build and
an optimized local Mac build using the available Kokkos backend. The optimized
executable used for timing must rerun the numerical acceptance suite.
Gate evidence must record exact source and dependency revisions, compiler and
flags, commands, pass/fail output, and hashes of retained result files. The
design, harness sources, tests, and evidence record must be committed together
as one coherent revision before either remediated gate is declared passed.

### Local performance-screening contract

The primary local timed region includes every operation needed to replace the
explicit update: coefficient and source assembly, activation-rate limiting,
matrix assembly, copies required because a solver overwrites coefficients,
factorization, all right-hand-side solves, method-required positivity handling,
state updates, and required Kokkos fences. Detailed Forward Euler clipping and
ledger diagnostics are collected in a separate instrumented pass. A matched
instrumentation-off/on run must report their overhead. Subtimers may separate
assembly, solve, and positivity-handling cost, but only the complete
uninstrumented advance supports a local CPU speedup estimate.

Benchmark one column and batches of 16, 64, 256, and 1024 columns; Forward
Euler `nsubmix` values 1, 2, 3, 5, 8, 16, and 21; and ACI intervals 33.33, 50,
75, 100, 150, and 300 s. Use identical active levels, inputs, and production-
feasible layouts for all methods. Homogeneous batches may control the workload
but must not reuse one column's matrix or factorization across columns. Any
layout conversion or host/device staging required by one method belongs inside
its timed region.

Use at least 10 untimed warmups and 30 measured repetitions. Continue to at
most 100 repetitions until the bootstrap 95% confidence-interval half-width
for the median is below 5%. Report complete-region time per column-call,
assembly, solve and positivity-handling time, memory/workspace, operation
counts, result status, and speedup relative to Forward Euler. Candidate
positivity checks must be timed; candidate clipping is prohibited.

For equal-interval timing, assemble coefficients once and advance every method
over the same frozen ACI interval. For timestep-refined or matched-error
timing, retain the same frozen operator, assemble it once, and include every
repeated integration step and its method-required positivity handling in the
timed region. This is an internal-integrator comparison and must not be
represented as changing the outer MAC/ACI subcycle count, which would re-run
other physics and rebuild the operator.

Synthetic L128 typical and strong profiles are sufficient for this screening
phase. Local batching, right-hand-side concurrency, operation counts, workspace
size, and sequential vertical depth can indicate whether a method is worth a
GPU implementation, but CPU timing cannot emulate GPU occupancy, warp behavior,
memory bandwidth, or synchronization cost. No result from this phase may be
labeled GPU performance. Production relevance requires sampled EAMxx columns
and later profiling on a representative GPU.

#### Standalone matched-error interval study (approved 2026-10-07)

This small, separately identified study measures local Kokkos Serial CPU cost
for the _same frozen initial state, operator, and total interval_ `T` at a
declared error tolerance. It does not change the earlier Task 5R matrix or its
`ndrop_summary_v2`/`ndrop_phase_samples_v1` interpretation. It makes no claim
about a production population, GPU speedup, whole ACI, or method selection in
EAMxx. The approved production-state replay gates remain independent.

The four core fixtures are the Cartesian product of `nlev={72,128}` and
`(mixing_scale,T)={(0.05,50 s),(0.5,150 s)}` using the existing deterministic
nonuniform benchmark generator, all 25 aerosol pairs/51 fields initialized
nonzero, strictly positive cloud fractions, and a one-column timing batch.
These are dimensionless synthetic moderate/stiff cases, not calibrated
atmospheric columns. Record the complete generator inputs and dimensionless
`T*max(rate_above+rate_below)` and `T*max(limited_activation_rate)` for each
fixture. Test isolated, contiguous, and all-zero cloud masks separately as
correctness cases; their zero-cloud transition is applied only after `T`.
The matched-error artifact orders field errors as droplet number, followed by
interstitial and cloud-borne fields for constituent 0, then constituent 1,
through constituent 24; it retains both norms against both references.

For each method use a positive integer internal count `N`, step `h=T/N`, and
the same frozen physical operator over all `N` steps. FE uses its existing
forced-substep full-interval entry point. New standalone IMEX and BE interval
entry points accept `(initial, operator_data, T, N, control)`; existing one-step
entry points retain their signatures and bitwise behavior at `N=1`. Prepare
the operator and immutable `h`-dependent matrices once per column-call;
refresh state-dependent activation sources and right-hand sides every step.
Retain the existing solver arithmetic, branch/positivity rules, numerical
checks, and update order. Apply common zero-cloud resuspension/reset once,
after the last internal step. The primary uninstrumented timer includes the
operator/matrix preparation, all allocations and copies required by the
method, every solve and positivity operation, final projection, output state,
and required fence. No setup or factor is shared across columns or timing
repetitions. Cross-step factorization reuse is permitted only by the existing
matrix-family rule; it is not required. The IMEX public combined Thomas
operation is retained. Detailed diagnostics and reference work remain outside
the primary timing measurement.

For these coupled fixtures, the primary accuracy measure is the maximum of all
51 separately normalized weighted-L1 and L-infinity field errors of the
_timed `N`-step state_ against each of two independently assembled DOPRI54
reference runs, using the field normalization defined above. The reference
runs use the existing field-scaled absolute tolerances, relative tolerances
`1e-12` and `1e-13`, and must both succeed with a maximum fieldwise
two-reference disagreement at most `2e-11`. Both sets of candidate-reference
errors must be at most the target. The two-tolerance check is an empirical
consistency test, not a rigorous bound. Richardson differences, observed
orders, and the earlier first-order/cross-method gates remain diagnostic and
must be reported; the existing `|q_h-q_(h/2)|` quantity estimates the _fine_
state's leading first-order error and must not be presented as the error of
the timed coarse state. The older Task 5R screening column retains its
historical meaning and is never silently reinterpreted.

The targets are `1e-2`, `3e-3`, `1e-3`, `3e-4`, and `1e-4`. Before timing,
evaluate accuracy at `N=1,2,4,...,8192` independently for each method and
fixture; FE may start at the first ladder point no smaller than its production
controller count. Retain failures and their first reasons. Clipped FE,
material-negative implicit, nonfinite, branch, residual, ledger, or unresolved
reference points are ineligible. For each target select the smallest
_tested_ qualifying `N` without consulting runtime; freeze and hash the
selection before an independent timing pass. If none qualifies, report
`not_attained_within_cap`, without extrapolation. Deduplicate identical
selected configurations and report the achieved error and any tolerance
overshoot; no global minimum-cost claim follows from a dyadic ladder.

The confirmation pass uses at least ten warmups and 30--100 complete,
method-interleaved timing blocks per fixture, the existing fixed-seed
10,000-resample bootstrap and exact benchmark median convention, and the
existing 5% median-interval stopping target. Stop only at complete blocks.
For a same-fixture, same-target eligible method pair, resample shared block
indices and compute the ratio of the two bootstrap median times. Report the
ratio, 95% interval, sample count, and an inconclusive result when the interval
includes one. A missing eligible pair is unavailable and is never reweighted
or called a tie. This study uses integer seed `0x4e44524f` with Python's
`random.Random` integer-seed stream, 10,000 index-vector replicates, the exact
benchmark median (`lower+(upper-lower)/2` for even samples), and sorted
zero-based replicate indices 249 and 9749 for the 95% endpoints. These rules,
Python version, and analysis-file hash are bound in the artifact; no default
quantile/interpolation choice is permitted. Instrumented phase observations
remain diagnostic-only.

Use a distinct versioned matched-error artifact with fixture/state/operator
identity, fixed `T`, method and `N`, source/dependency/binary identity,
reference self-check, per-field errors and gate outcomes, frozen selections,
raw paired timing blocks, setup/factorization policy, completion hashes, and
explicit failure/missing statuses. The exploratory fixed-count results cannot
be promoted into this study. Run focused interval and reference tests plus the
existing optimized and debug ASan/UBSan numerical suites before a timed pass;
test `N=1` equivalence, odd/even step composition, state-dependent source
refresh, one-time projection, diagnostics-on/off state equality, failure
aggregation, and artifact corruption. A retained result additionally needs a
clean exact source/binary/dependency chain. The initial bounded study uses a
single serial worker and a two-hour benchmark-stage cap; hitting the cap
produces explicit unfinished outcomes and does not relax a numerical gate.

#### Production-controller-matched standalone extension (approved 2026-10-08)

This extension changes the primary standalone comparison from a common
arbitrary absolute-error target to the accuracy actually delivered by Forward
Euler at the existing production controller count. It remains an offline
Kokkos Serial study of frozen synthetic operators. It does not authorize or
represent production-state replay, climatological weighting, GPU or whole-ACI
performance, or a production method recommendation. The preceding
absolute-target study and its `ndrop_matched_*_v1` artifacts remain historical
exploratory evidence and are not reinterpreted or merged.

For each case, obtain `N_FE=production_nsubmix(prepared_operator,T)` without
rounding it to the dyadic ladder. Execute Forward Euler over the complete
frozen interval `T` with exactly `N_FE`. Both independently assembled DOPRI54
references must pass the existing two-tolerance self-check. Define the
case-specific target

\[
 E_{\mathrm{FE}}=\max(E_{\mathrm{FE,normal}},
                        E_{\mathrm{FE,tight}}),
\]

where each argument is the existing maximum across the 51 separately
normalized weighted-L1 and L-infinity field errors. Retain both complete
per-field vectors against both references. Exact finite comparison
`E_method <= E_FE` is the primary matching rule; no runtime, confidence
interval, or later method result may affect selection. Report per-field
candidate/FE comparisons as diagnostics, but do not add an unapproved
fieldwise-dominance gate.

A case is ineligible when `N_FE` is nonpositive, either reference is
unresolved, or the controller-count FE result fails any existing clipping,
material-negativity, finite-value, branch, residual, ledger, or state gate.
For IMEX and Backward Euler independently, test positive integer counts in
strict ascending order starting at one. Stop at the first numerically eligible
count whose primary error is no greater than `E_FE`; therefore every smaller
integer count is present in the accuracy artifact with its measured error or
first failure reason. The count cap is 8192. Reaching the cap or the bounded
accuracy-stage resource limit without a match yields an explicit
`not_attained_within_cap` or `accuracy_budget_exhausted` outcome and no timed
pair. After selection, prescribed adjacent/refined diagnostic points may be
evaluated only in a separately labeled confirmation pass; they cannot revise
the frozen selection without creating and approving a new selection artifact.

The initial continuity fixtures remain the existing L72/L128 moderate and
strong synthetic operators. A second, separately labeled synthetic
controller-stratification matrix targets actual returned `N_FE` values
`{1,2,4,8,16,32,64,128,256,512}` at L72 and L128. Fixture construction may
vary only the existing synthetic mixing-scale input; every fixture records
the requested bin and actual controller result and is rejected if they do not
match. These unweighted fixtures diagnose workload dependence only and must
not be called typical, climatological, or production-frequency evidence.
Adding production states or weights remains subject to the independent
production-state replay approval sequence.

Use separate closed artifact schemas
`ndrop_controller_matched_accuracy_v1`,
`ndrop_controller_matched_selection_v1`,
`ndrop_controller_matched_timing_v1`, and
`ndrop_controller_matched_analysis_v1`. Bind the exact case/operator identity,
controller inputs and returned `N_FE`, controller-count FE state and reference
errors, every tested implicit count and first failure reason, selected counts,
achieved errors and candidate-to-FE error ratios, state checksums, source and
dependency identities, exact binary and commands, completion records, and
final manifest. The timing validator must require each timed state checksum to
equal its accuracy-qualified checksum.

Only after the complete accuracy artifact validates and the immutable
selection is frozen may a separately approved timing pass run. Time FE at
exact `N_FE` and each implicit method at its selected count using the same
complete uninstrumented region and final-only projection already defined
above. Retain ten warmups, 30--100 complete interleaved blocks, the 5% median
precision rule, the fixed 10,000-resample paired bootstrap, and the exact RNG,
median, and endpoint rules. Report FE/IMEX, FE/BE, and IMEX/BE time ratios;
missing matches remain unavailable and are not reweighted.

Before that timing approval, focused tests must prove exact non-dyadic
controller-count FE execution, reconstruction of `E_FE`, strict ascending
integer selection, retention of all lower failures, candidates below/equal to/
above `N_FE`, cap and budget outcomes, reference and FE failures, per-field
ordering, diagnostics-on/off state identity, timing-to-accuracy checksum
equality, one-time zero-cloud projection, and artifact corruption rejection.
The optimized numerical/support suite and debug ASan/UBSan suite must pass.
No operator equation, numerical arithmetic or reduction order, solver,
limiter, clipping rule, boundary, coupling, production API, dependency, or
gitlink changes are authorized.

#### Physical-column common-accuracy extension (approved 2026-10-08)

The user approved the in-chat physical-column design and local implementation.
This separate standalone extension leaves every prior fixture and artifact
unchanged. It authorizes local offline implementation and untimed synthetic/
repository-fixture validation only; another timed pass requires separate
approval after the selections and numerical evidence are frozen. It grants no
production capture/replay, directional-rate amendment, Task 5 access,
production-source/API, dependency/gitlink, or production-recommendation authority.

**Fixture assumptions and equation map.** Four fixtures 28--31 are respectively
L72 moderate/strong and L128 moderate/strong, all at `T=150 s`. Read the
repository `components/homme/test/vcoord/acme-72{i,m}.ascii` and
`scream-128{i,m}.ascii` A/B tables, retaining their bytes/hashes; do not
substitute an evenly spaced or inferred grid. These are repository coordinate
families, not proof of equality to a particular production initial-condition
file. Set `p0=100000 Pa`, `ps=101780 Pa`, `Rd=287 J kg^-1 K^-1`,
`g=9.80665 m s^-2`, and fixed `Tair=289 K`. The isothermal hydrostatic
geometry approximation is explicitly an idealization, not a reconstruction
of the complete DYCOMS thermodynamic sounding:
`p=A*p0+B*ps`, `z=(Rd*Tair/g)*log(ps/p)`,
`dz=(Rd*Tair/g)*log(p_bottom/p_top)`, `rho=p/(Rd*Tair)`,
and cell mass `w=(p_bottom-p_top)/g` [kg m^-2]. Preserve the full grids and
top-to-bottom ordering; require positive monotone pressures, midpoint pressures
strictly inside interfaces, positive thicknesses, and finite derived values.

Cloud base/top are 600/840 m, anchored to the
[DYCOMS-II RF01 case specification](https://atmos.uw.edu/~breth/GCSS/DYCOMS_RF01.pdf).
Cell liquid-cloud fraction is `0.8` times the geometric fraction of the cell
intersecting that cloud layer; it is zero elsewhere. Background interface
diffusivity is 3/30 m2 s^-1 below 840 m for moderate/strong, and 0.01 m2 s^-1
above. Prescribed cloud-base `wtke` is 0.3/1.0 m s^-1. Those turbulence values
are sensitivity assumptions, not observations or SHOC closure outputs.
At each diagnosed cloud-base interface (`cldn[k]>0.01` and
`cldn[k]-cldn[k+1]>0.01`), use the locked source's dimensional relation
`K=wtke*(z[k]-z[k+1])`. The standalone shared finite-volume exchange is
`E_face=rho_face*K/(z_above-z_below)` [kg m^-2 s^-1], with exact zero no-flux
endpoints. Ordinary preparation divides this by cell mass and applies the
existing overlap/activation limiter. The controller-only first-active-level
duplicate is the derived first downward rate, never a new physical tendency.
This is a physically scaled standalone operator, not an assertion that its
shared exchange represents every locked production directional-rate pair.

Use prescribed activation fractions `{.35,.4,.45,.5}` for number and
`{.25,.3,.35,.4}` for mass, with pre-limiter rates
`fraction*wtke/sqrt(2*pi)*(rho_face/rho_cell)*(cldn[k]-cldn[k+1])/dz[k]`
only at those cloud bases. These are idealized Gaussian-flux fractions, not
outputs of a production activation calculation; retain limiter outcomes.
Initial modal total number concentrations are `{100,200,1,10}` cm^-3,
dry geometric mean diameters `{.12,.04,1,.10}` micrometers, geometric standard
deviations `{1.8,1.6,1.8,1.6}`, and assumed effective dry density
`1500 kg m^-3`. Total mode mass concentration is
`N*(pi/6)*Dg^3*exp(4.5*log(sigma)^2)*1500` [kg m^-3], not an independent
arbitrary mass amplitude. Divide number/mass concentrations by cell density
for native mixing ratios. Distribute mode mass equally among its existing
7/4/7/3 species as an explicitly idealized composition, without changing the
canonical 25-constituent mapping. Partition fraction `0.1*cldn` cloud-borne;
initialize droplet number to `80 cm^-3*cldn/rho`. Moderate/strong initial
states, cloud masks, geometry, and composition are identical within a grid
pair. All 25 pairs are represented; clear-cell cloud-borne/droplet zeros are
physical occupancy, never reduced executed width.

**Accuracy and execution.** Record an actual adaptive FE baseline using
`production_nsubmix` unchanged. Separately search every integer FE count from
that count upward, and every integer IMEX/BE count from one upward, for common
error ceilings `1e-2` and `1e-3`, through the existing 8192 cap. Stop a method
at its first eligible tight-target match; retain every lower failure and error
and all explicit cap/budget outcomes. FE at the controller count uses the
adaptive entry point; larger counts use the existing forced-count test entry
point and are labeled accuracy-refined FE, never production-selected work.
Use the existing maximum of the 51 fieldwise normalized weighted-L1 and
L-infinity errors against both reference tolerances, self-check `2e-11`, and
all existing positivity, clipping, finite, residual, ledger and boundary gates.
Error matching is to the same ceiling, not identical achieved errors.
For each deduplicated selection retain `N,2N,4N` confirmation against both
references; missing/failed/non-decreasing confirmation withholds timing, not
selection, without relaxing a gate. Baseline and common-target costs are
separate reporting questions. Retain controller and physical stiffness,
cloud-base/limiting levels and spacings, actual/selected counts, both-reference
per-field errors, and explicit exclusions. Do not guarantee small counts or
calibrate inputs to a requested FE bin.

Use separate `ndrop_physical_column_<kind>_v1` artifacts, with closed fixture
identity, copied grid/executable/source/build hashes, completion records and
accuracy-only frozen selections. The new CSV header and validator must retain
the baseline/measurement/confirmation distinction and blank inapplicable
failure metrics. Future timing uses only the complete uninstrumented advance,
ten warmups, 30--100 complete interleaved blocks, existing 5% stopping and
fixed-seed 10,000 paired-bootstrap conventions. Detailed/reference work stays
outside primary timing; unavailable pairs are not reweighted. No timings or
method preferences follow from untimed accuracy tests.

**Implementation/testing scope.** Modify only this document and the existing
standalone `benchmark_ndrop_matched_error.cpp`, `ndrop_matched_artifacts.py`,
`test_ndrop_matched_artifacts.py`, `CMakeLists.txt`, and `README.md` under
`docs/analysis/evidence/ndrop_integrator/`. Numerical implementations are
unchanged. Test grid/header/extent/order corruption, dimensional exchange and
modal mass/number consistency, paired-state identity, zero-cloud handling,
controller-selected versus accuracy-refined counts, reference/failure metric
applicability, complete integer search, common-ceiling selection, confirmation
and missing-pair rejection, copied-grid/source/completion corruption, and
diagnostics-enabled/disabled state equality. Run the unchanged numerical and
benchmark-support unit suites in debug ASan/UBSan and optimized builds, but
select tests explicitly to exclude timed benchmark smoke tests. Record
baseline/index/outside-scope preservation and a focused static review/recheck.
This idealization can support local algorithm screening; atmospheric
representativeness and production fidelity remain unestablished.

#### Physically scaled controller-matched extension (approved 2026-10-08)

The user's in-chat approval to implement the proposed extension authorizes
local implementation, offline unit tests, and untimed accuracy/selection for
the unchanged physical fixtures 28--31 above. It does not authorize a timed
run, production replay/capture, performance analysis, or a production claim.
The fixed 1%/0.1% study remains separate and unchanged.

Apply the existing production-controller-matched selection contract to these
four idealizations: obtain the actual natural `N_FE` from the unchanged
adaptive FE entry point and set `E_FE` to its maximum normalized error across
all 51 fields, both norms and both certified references. Never refine the
primary FE workload to satisfy a fixed ceiling. Search IMEX and BE from one
in consecutive integers for their first eligible `E_method <= E_FE`. Preserve
every lower failure/error and the existing cap, resource, reference,
positivity, clipping, finite, residual, ledger, and boundary gates. `N,2N,4N`
confirmation can withhold timing but cannot reselect a configuration. Matching
the maximum does not imply fieldwise dominance or establish that FE accuracy
is scientifically sufficient.

Use the separate `ndrop_physical_controller_<kind>_v1` family for identity,
accuracy, completion, selection, timing admission, timing, analysis, and final
manifest. The closed accuracy CSV extends the existing controller-matched
header with `turbulence`, `cloud_base_wtke_m_s`, `controller_rate_t`,
`cloud_base_level`, `cloud_base_spacing_m`, `limiting_level`,
`limiting_spacing_m`, and `activation_limiter_count`; `mixing_scale` denotes
the fixture's background diffusivity in m2 s-1, not a dimensionless multiplier.
All physical characteristics remain bound across every case row. The limiter
count is blank for candidate-free terminals, otherwise a nonnegative count.
Identity binds `fixture_policy=dycoms_inspired_isothermal_repository_grid_v1`,
`comparison_policy=natural_fe_controller_error_v1`, exact four-case coverage,
copied executable and four grid hashes, source/build/dependency hashes,
commands and unchanged statistical protocols. Reject mixed families, stale
or missing completions, altered grids, targets not equal to FE error, refined
FE primary counts, incomplete integer searches, and unqualified timing states.
The as-built offline entry point is
`ndrop_controller_matched_artifacts.py run-accuracy --physical-columns`, then
`freeze`. It invokes `--physical-controller-accuracy` with the copied grids;
future separately approved timing invokes `--physical-controller-timing`.

Only the complete diagnostics-disabled advance and required fence supply
future primary timing. Keep references, error checks, serialization, and
bootstrap work outside it. Retain ten warmups, 30--100 complete interleaved
blocks, the existing 5% stopping rule and fixed 10,000-resample paired
bootstrap; missing eligible pairs remain unavailable, without reweighting.
Separate explicit timing approval is required after sealed accuracy,
selection and timing-admission evidence is reviewed.

Change only this document, `benchmark_ndrop_matched_error.cpp`,
`ndrop_controller_matched_artifacts.py`, its existing Python test file,
`README.md`, and necessary `CMakeLists.txt` integration in the standalone
directory. Reuse the existing physical-fixture diagnostics-on/off and
unchanged numerical unit gates in debug ASan/UBSan and optimized builds.
Add artifact corruption/failure and untimed local-binary tests proving
controller-only FE workload, first implicit matches, target reconstruction,
confirmation separation, copied-grid/provenance integrity, and future timing
checksum/pairing admission using explicitly fabricated timing fixtures only.
No numerical implementation changes are authorized. This remains an unweighted
standalone shared-exchange idealization, not established locked directional-rate
production fidelity; the deferred representation issue is not resolved here.

#### Paired physical high-work geometry sensitivity (approved 2026-10-08)

The user approved the proposed six moderate/strong pairs (12 cases, including
the four existing native cases) for local implementation, offline tests, and
untimed controller-matched accuracy/selection. Separate post-selection timing
approval remains required. This approval adds no production capture/replay,
directional-rate, retained Task 5 access, production-source/API, dependency,
gitlink, performance-analysis, publication, staging, or commit authority.

**Frozen fixtures.** Retain native IDs 28--31 and their numerical inputs
bit-for-bit. Add IDs 32/33 = L72-derived 2x moderate/strong; 34/35 = L72-derived
4x; 36/37 = L128-derived 2x; and 38/39 = L128-derived 4x. Subdivide each whole
native cell whose open height interval intersects 400--1000 m; no other cell
is changed. Split its log-pressure interval equally into r=2 or 4 children,
with `p_j=p_top*exp(j*log(p_bottom/p_top)/r)` for interior interfaces and
original endpoints copied exactly. Child midpoint pressure is the geometric
mean of its interfaces; unmodified native midpoints remain unchanged. Use the
same hydrostatic equations, cloud slab 600--840 m, prescribed moderate/strong
K=3/30 m2 s^-1 and wtke=0.3/1.0 m s^-1, 150 s interval, activation construction,
and no-flux endpoints as the physical-column extension. Recompute geometric
cloud fractions and face exchange on each mesh, not production calculations.
The actual level counts are 72/79/93 and 128/141/167; never label derived
meshes as native L72/L128, pad/truncate them, or calibrate inputs to FE counts.
The L128 counts correct the proposal's derived 140/164 labels: the exact grid
parser finds 13 intersecting parent cells rather than 12. The approved band
predicate, subdivision, and every physical input are unchanged.

**Conservative initialization.** Subdivide the native common old state, not
an evolved state. For a parent's children with dry masses w_i and geometric
cloud fractions c_i, define `c_bar=sum(w_i*c_i)/sum(w_i)`. Keep each child's
interstitial mixing ratio equal to its parent's; set its cloud-borne and
droplet mixing ratios to their respective parent value times `c_i/c_bar`.
When `c_bar=0`, require the corresponding parent cloud-borne/droplet values to
be zero and copy zeros. Copy unsplit states verbatim. This preserves each of
the 51 parent-cell dry-mass-weighted field integrals, modal mass/number ratios,
nonnegativity, and clear-child zeros, up to bounded floating-point roundoff.
The initialized child concentrations are conservative representations of the
native state, not newly evaluated constant-volume concentrations. Require
parent and column integral agreement within 128 binary64 epsilon relative
error (exact zero agreement for zero integrals). Initial states/geometry are
bitwise identical within each moderate/strong pair. This fixture-only remap
changes no integrator arithmetic, reduction order, limiter, or projection.

**Admission, artifacts, and reporting.** Reuse the natural-controller FE
error target, strict ascending integer IMEX/BE search, 8192 cap, 7200 s accuracy
budget, two certified references, all existing correctness gates, and N/2N/4N
confirmation without relaxation or reselection. Keep the previous common
ceilings and all previous artifact families separate. Introduce the closed
`ndrop_physical_geometry_controller_<kind>_v1` family for identity,
accuracy, completion, selection, timing-admission, timing, analysis, and final evidence, with
`fixture_policy=dycoms_inspired_conservative_cloud_band_refinement_v1` and
`comparison_policy=natural_fe_controller_error_v1`. Bind exact 12-case coverage,
source/design/build/dependency/executable and copied native-grid hashes. Extend
the physical-controller CSV by `parent_grid_levels`, `geometry_subdivision`,
`refinement_lower_m`, `refinement_upper_m`, `limiting_height_m`,
`limiting_above_diffusivity_m2_s`, `limiting_below_diffusivity_m2_s`, and
`limiting_cloud_base_faces`. The last field counts the limiting cell's two
adjacent faces using the existing cloud-base predicate (0--2), distinguishing
cloud-base involvement from background diffusion without claiming a unique
mechanism for a mixed-rate limiter. Native subdivision is 1, band bounds
400/1000 m are policy metadata even for unmodified native controls, and
actual `levels` never means parent-grid size. All identity metadata are fixed
across a fixture's rows. Preserve all lower failures, achieved errors, both
reference field vectors, actual FE counts, limiter activity, limiting levels/
spacings and exclusions. Missing pairs are not reweighted. Numerical admission
does not itself authorize timing. The separately approved bounded timing stage
below supersedes the initial execution prohibition only within its exact scope.
The as-built untimed entry points are
`benchmark_ndrop_matched_error --check-physical-geometry-fixtures`,
`ndrop_controller_matched_artifacts.py run-accuracy --physical-geometry`, and
`freeze`; the driver invokes `--physical-geometry-controller-accuracy` against
the copied native grids. All outputs are exploratory local scratch evidence,
not a remediated retained Task 5 dataset or production replay admission.

**Verification and limits.** Modify only this document and existing standalone
`benchmark_ndrop_matched_error.cpp`, `ndrop_controller_matched_artifacts.py`,
`test_ndrop_controller_matched_artifacts.py`, `CMakeLists.txt`, and `README.md`.
Prove full-grid/endpoints/unsplit preservation, child extents and monotonicity,
dry-mass and all-field integral conservation, paired identity, lognormal modal
ratios, clear-cell handling, unchanged native outputs, diagnostic-on/off state
identity, natural controller counts, exact metadata/coverage, corruption and
completion rejection, strict selection and confirmation separation. Run only
explicitly selected untimed numerical/support/fixture tests in optimized and
ASan/UBSan builds, then bounded untimed accuracy and frozen selection. These
are unweighted prescribed-turbulence resolution sensitivities, not the current
EAM/EAMxx grid, a full LES/SHOC state, observed tail frequencies, or validated
locked directional-rate production fidelity. Native strong limiting cells
need not lie at cloud base; high work and implicit benefit are not guaranteed.
Retain naturally emerging counts and exclusions rather than tuning K, wtke,
dt, geometry, or tolerances after seeing results. Future timing retains the
complete uninstrumented advance, 10 warmups, 30--100 paired blocks, existing
5% stopping and fixed 10,000-resample bootstrap, under separate approval.

**Bounded timing stage (approved 2026-10-12).** The user's in-chat approval
authorizes narrow local timing/analysis software enablement and one bounded
standalone timing pass for these 12 frozen cases and their 36 configurations.
Use the existing same-file offline allowlist above; change no fixture, numerical
implementation, tolerance, gate, controller or statistical rule. Bind the prior
accuracy SHA-256 `7bf06451e74ad1b9930e8fd8df00e85ffbcbcb7772798ef5b87144b501d92cb6`
and selection SHA-256 `c553ce935b85a3ec675f9099ca900c3cfddcf4c2a975c3f088a16cbdd36ebdf7`
in the local approval record. Because the timing flag changes the executable,
build/test in new unsynced scratch, rerun full untimed accuracy with that exact
optimized timing binary, and require byte-identical accuracy CSV and unchanged
selected configurations/checksums before sealing its fresh identity/selection.
Do not update or merge the earlier dataset. Renew scope approval if those
numerical results change; provenance-only identity hashes may change normally.

The geometry driver uses `--physical-geometry-controller-timing` and requires
an immutable local `timing_approval.json` with exactly `artifact_schema`
(`ndrop_physical_geometry_controller_timing_approval_v1`), `approval_reference`
(`user_in_thread_2026-10-12_geometry_timing`), `authority`
(`bounded_local_12_case_geometry_timing_and_analysis`), `identity_sha256`,
`accuracy_sha256`, `selection_sha256`, `timing_admission_sha256`,
`prior_accuracy_sha256`, `prior_selection_sha256`, and `command` (the exact
retained-executable argument vector, including copied grids and fixed timing
parameters). Require current accuracy to equal the approved prior accuracy;
all current hashes and the reconstructed command must match. This ordinary
JSON/hash record binds already granted human authority; it is not an approval
service, automatic grant, or proof of production authority. Missing/stale records
fail before timing. Existing timing outputs/partials/completions/execution
records prevent another timed process in the same directory. A failed or
additional pass requires renewed bounded approval and a new output location.

Measure only the complete diagnostics-disabled advance and required fence:
coefficient/source preparation, required allocations/copies, limiting, all
factorizations/solves, positivity handling and state updates remain included.
Fixture construction, certified reference/error checks, checksums, serialization,
stopping bootstrap and paired-ratio analysis stay outside primary timing.
Retain ten warmups per configuration, 30--100 complete interleaved blocks,
existing 5% median-precision stopping, and the 7200 s budget. This is one local
Serial process; builds use at most two workers. Expected generated storage is
under 1 GB, with no retained Task 5 dataset. Record actual counts and stop reasons.

After validated timing and completion, write immutable `timing_execution.json`
with exactly `artifact_schema` (`ndrop_physical_geometry_controller_execution_v1`),
`identity_sha256`, `selection_sha256`, `timing_approval_sha256`, `timing_sha256`,
`command`, `process_exit_code` (0), `timing_rows`, and `warmups_by_configuration`
(ordered records with exactly `fixture`, `method`, `steps`, `warmups`, all
warmups=10). Counts are source/command-bound execution evidence, not fabricated
measurements. Validators reconstruct them and require approval/execution records
for analysis and final admission; final manifests hash both. Geometry analysis
retains all physical and refinement metadata, achieved errors, paired raw samples,
time ratios/intervals, exclusions, exact command/provenance and the distinct
`physical_geometry_controller_idealization_unweighted` evidence scope. Use the
unchanged ratio-of-medians paired bootstrap and missing-pair/no-reweighting rule.

Test corruption of approval/execution/hash bindings, geometry parser coverage,
qualified timing checksums and complete blocks, median/ratio reconstruction,
and unchanged original artifact families with explicitly fabricated samples
only. Re-run relevant optimized and ASan/UBSan numerical/support/fixture suites
and Python offline tests before timing, then perform static review and a brief
recheck. The only new conclusions permitted are bounded local method-time
comparisons on these prescribed-turbulence idealizations. Production source,
capture/replay, directional-rate resolution, dependencies/gitlinks, retained
Task 5 access, remote/HPC, production/model-wide claims, publication, staging,
commits and pushes remain excluded.

#### Controlled moderate/strong mixing fixtures (pending approval)

**Status: PENDING EXPLICIT APPROVAL; DOCUMENT-ONLY PROPOSAL.** The user's
request to prepare this extension authorizes this document revision only.
It does not approve the fixture contract, implementation, or another timed
run. This is a four-case addition to the production-controller-matched
standalone study, not a replacement or reinterpretation of its existing
fixtures, artifacts, or executed evidence. All production-state replay
approvals and the unresolved directional-rate/shared-exchange issue remain
separate and unchanged.

**Question and controlled inputs.** Compare FE, IMEX, and BE under moderate
versus strong synthetic mixing while holding the total interval, initial
state, grid, cloud fractions, pressure weights, and prescribed activation-rate
inputs fixed within each level-count pair. The older continuity fixtures
change both `T` (50 versus 150 s) and mixing scale, and the existing generator
also scales activation with mixing; they are not a mixing-only contrast.
Here "moderate" and "strong" label the chosen synthetic scales, not observed
atmospheric frequencies, diffusivities, TKE, or calibrated production columns.

| Fixture ID | Fixture family | Levels | Mixing label | `mixing_scale` | `activation_mixing_scale` | `T` |
|---:|---|---:|---|---:|---:|---:|
| 24 | `controlled_mixing` | 72 | moderate | 0.05 | 0.05 | 150 s |
| 25 | `controlled_mixing` | 72 | strong | 0.5 | 0.05 | 150 s |
| 26 | `controlled_mixing` | 128 | moderate | 0.05 | 0.05 | 150 s |
| 27 | `controlled_mixing` | 128 | strong | 0.5 | 0.05 | 150 s |

Retain all 25 nonzero aerosol pairs/51 fields and the exact deterministic
`initial_state` and nonuniform geometry in
`benchmark_ndrop_matched_error.cpp`. For each level count, form the ordinary
existing generator's operator at mixing scale `0.05` as the activation-input
baseline. Copy its four-mode number and mass activation-rate arrays, including
their zero bottom entries, without modification into **both** new cases.
Their coefficients remain `{0.35,0.4,0.45,0.5}` and `{0.25,0.3,0.35,0.4}`,
respectively, times the baseline downward mixing rate. Thus
`activation_mixing_scale` identifies the generator input used to construct
the fixed activation arrays; it is not a new integrator input or a change to
the activation equation. Prepared mixing and activation rates retain their
existing `s^-1` meaning; `T*rate` diagnostics are dimensionless.

Use the existing interface-exchange formula separately at each case's
`mixing_scale`, with unchanged pressure weights and overlaps. Retain the
existing controller-only first-level term `0.25*mixing_scale`; it follows
mixing strength and never enters the physical tendency. Invoke the ordinary
operator preparation and limiter, without bypassing either. The fixed
activation inputs are below the moderate case's downward-rate limits and
therefore also below the strong case's limits; require the prepared activation
arrays to remain identical and the activation-limiter count to be zero.
A mismatch invalidates the fixture rather than being repaired by altering
rates, limits, or arithmetic. State-dependent activation sources and transfers
are still refreshed normally and may differ as the evolving states differ.
The fixed quantity is the prescribed rate input, not an artificial fixed
source or flux.

**Selection, numerical admission, and timing.** Use the unchanged controller
to obtain and execute each case's actual `N_FE`; do not calibrate the new
cases to a requested bin or force their counts. `requested_nsubmix` is blank.
Use that case's certified FE error as `E_FE`, retaining the existing
independent DOPRI54 references at `rtol=1e-12` and `1e-13`, field-scaled
absolute tolerances, and self-consistency limit `2e-11`. Search each implicit
method at every positive integer count in ascending order through the existing
8192 cap, stopping at the first numerically eligible `E_method <= E_FE`.
Preserve all existing finite-value, clipping, positivity, branch, residual,
ledger, state, cap, budget, and missing-pair rules. Before timing, retain
eligible `N,2N,4N` refinement evidence against both references for all three
methods and the pilot's separately hash-bound timing admission: both aggregate
norms must strictly decrease against both references. Missing or failed
refinements withhold timing without changing selection; observed orders
remain diagnostic and the existing formal convergence gates remain required.

After a separate timing approval, measure only the complete uninstrumented
advance with unchanged setup/factor reuse policy, source refresh, allocations,
solves, positivity/state handling, final-only zero-cloud reset, and fence.
Retain ten warmups, 30--100 complete interleaved method blocks, the existing
5% per-method median-precision stopping rule, and the existing paired-ratio
Python RNG/seed, 10,000 resamples, even-sample median, and endpoint conventions.
Detailed and phase diagnostics remain outside primary timing. Stop only at
complete blocks and retain stopping counts/reasons and every exclusion.

**Artifact identity and reporting.** Leave fixtures 0--23 and all existing
`ndrop_controller_matched_*_v1` definitions unchanged. Isolate this four-case
extension in `ndrop_controlled_mixing_<kind>_v1` artifacts, where `<kind>` is
exactly `identity`, `accuracy`, `selection`, `timing_admission`, `timing`,
`completion`, `analysis`, or `final`. Reuse the corresponding existing
encoding, fields, numerical rules, and hash/completion transaction mechanics,
with only these input-identity additions: identity records
`fixture_policy="controlled_mixing_v1"`, `activation_mixing_scale=0.05`, and
the exact fixture list `[24,25,26,27]`; the accuracy CSV inserts
`activation_mixing_scale` immediately after `mixing_scale`; each analysis
fixture entry includes that same value. The fixture family and four exact
table rows are closed definitions, not caller-supplied populations. Reject
mixed policies/schemas, missing or changed activation identity, and any
coverage other than these four fixtures. Bind the exact generator/source,
dependencies, executable, commands, selections, numerical/refinement evidence,
timed-state checksums, and outputs as in the existing standalone protocol.
Do not merge or substitute prior artifacts for this extension.

For each grid show moderate and strong cases side by side: actual mixing and
activation stiffness diagnostics, reference self-check and tolerances,
controller count, selected implicit counts, achieved errors against both
references, per-field candidate/FE diagnostic comparisons, raw paired times,
FE/IMEX, FE/BE, and IMEX/BE ratios/intervals, and exclusions. Error matching is
**within each case** against its own FE target. Targets can differ between
mixing regimes; side-by-side costs are not a common-error cross-regime speedup
or a fieldwise-dominance result. No cross-regime paired confidence interval,
new common error target, population weighting, or production preference is
authorized by this extension.

**Permitted future work and approval boundary.** A separate explicit approval
of this section would permit only fixture construction, corresponding closed
artifact/validator/report support, documentation, and local offline tests in
`benchmark_ndrop_matched_error.cpp`, `ndrop_controller_matched_artifacts.py`,
`test_ndrop_controller_matched_artifacts.py`, and `README.md` under
`docs/analysis/evidence/ndrop_integrator/`; `CMakeLists.txt` may change only if
needed to register an untimed fixture/integration check. No integrator,
reference implementation, production adapter, numerical method, production
source/API, dependency, or gitlink edit is included. That approval would also
authorize bounded **untimed** accuracy/refinement and frozen-selection
preparation for these four cases: one Serial worker, the existing 7200 s
accuracy-stage cap, and a new task-specific unsynced directory below
`/Users/burr114/e3sm-scratch/`. It does not authorize timed execution or
performance analysis. Present the exact timed commands, new output paths,
resource bound, and passed admission evidence for separate approval afterward.

Future tests must demonstrate exact common initial/geometry/cloud/activation
inputs within each grid pair, correct exchange/controller inputs, unchanged
prepared activation and zero limiting, actual adaptive FE execution, all-field
reference errors, ascending integer matching, finite/numerical/refinement
admission, unchanged diagnostics-on/off states, complete paired blocks and
timing checksums, closed identity/coverage and corruption rejection. Verify
that legacy fixture generation and its v1 artifacts remain unchanged. Run
the existing debug ASan/UBSan and optimized numerical/support suites before
timing. All such tests and four-case results are **proposed**, not executed
evidence. No production capture/replay, retained Task 5 access, remote/HPC
activity, publication, staging, commit, or production/model-wide claim is
authorized. Any requirement needing another source file or a numerical change
must be reported for renewed approval rather than silently substituted.

#### Task 5 evidence-remediation amendment

**Status: APPROVED AND IMPLEMENTED; RETAINED EXECUTION PENDING.** The user
approved this amendment at commit
`61ae8249af0eb1c7e2fbc5d2da98a80dd1cff0a7` on 2026-08-25. The standalone
benchmark, analyzer, and transaction tooling implement it through review-fix
commit `d55ef90e46` and subsequent artifact-promotion hardening commit
`19c663a61ef73f99eb0b001d2947fd24b06e5e7b`; no remediated retained
performance run has yet passed the paired-artifact contract.

This amendment repairs three evidence-contract defects found by static review
of Task 5 commit `dba6d00598d86a239a2bf2743d95325821c9378c`. It changes
benchmark labels, timing diagnostics, serialization, and provenance only. It
does not change an operator equation, numerical method, solver arithmetic,
state transition, limiter, clipping rule, positivity threshold, boundary,
coupling, or production EAMxx/MAM4xx interface.

##### Fixed executed width and nonzero-state occupancy

All methods continue to execute droplet number plus 25 interstitial/cloud-borne
aerosol pairs, or 51 transported fields. The matrix input previously named
`active_pairs` changes only which leading aerosol pairs are initialized
nonzero. It must therefore be renamed `nonzero_aerosol_pairs`, and
`transported_fields` must be replaced by `nonzero_state_fields`. Every summary
row must satisfy

\[
n_{\mathrm{nonzero\ state}}=1+2n_{\mathrm{nonzero\ aerosol\ pairs}},
\qquad
n_{\mathrm{executed\ aerosol\ pairs}}=25,
\qquad
n_{\mathrm{executed\ transport\ fields}}=51.
\]

Every summary row declares `artifact_schema=ndrop_summary_v2` and contains `run_id`,
`nonzero_aerosol_pairs`, `nonzero_state_fields`,
`executed_aerosol_pairs=25`, and `executed_transport_fields=51`. Neither the
analysis nor the report may call the occupancy input an executed width or use
it to claim transported-width scaling. Occupancy sensitivity may be reported
only under its explicit nonzero-state label.

##### Linear-solver phase meanings

The primary uninstrumented complete advance remains the only timing admitted
to CPU speedup, matched-error, break-even, hybrid, or recommendation
calculations. Instrumented phase and fence samples are diagnostic evidence
only and may not substitute for primary time or support cross-method phase
ratios.

Each summary row has `linear_phase_mode` with exact method mapping:

| Method | `linear_phase_mode` | Separate factorization | Separate solve | Combined factorization-solve |
|---|---|---|---|---|
| Forward Euler | `none` | unavailable | unavailable | unavailable |
| IMEX Euler | `combined` | unavailable | unavailable | measured |
| Backward Euler | `separate` | measured | measured | unavailable |

The exact summary timing columns are
`factorization_seconds_per_column`, `solve_seconds_per_column`, and
`combined_factorization_solve_seconds_per_column`. Unavailable CSV values are
blank, never numeric zero. The one IMEX combined accumulator receives both
existing scalar-system observations: each begins before allocation and copying
of that system's tridiagonal coefficient views and ends after its supported
public serial `ekat::tridiag::thomas` factorization-and-solve call. Its
cumulative value therefore contains the cloud and interstitial coefficient-
view allocation/copy, both Thomas operations, and all right-hand-side solves.
The implementation must not call EKAT `impl` helpers or introduce another
solver merely to split this phase. Backward Euler retains its existing
separate factorization and solve boundaries.

##### Raw instrumented-sample sidecar

Every benchmark invocation requires a run identifier matching
`[A-Za-z0-9._-]{1,128}` and a phase-sample output path. The summary and sidecar
carry the same run identifier. The `ndrop_phase_samples_v1` sidecar contains
one row for every zero-based instrumented repetition of every summary row with
`validation_ok=true`, including timed rows later excluded by eligibility.
Rows with `validation_ok=false` have zero repetitions and no sidecar rows.

Each sidecar row records `artifact_schema`, `run_id`, `case_id`, `method`,
`linear_phase_mode`, `instrumented_repetition`, and seconds per column for the
instrumented complete advance, setup, separate factorization, separate solve,
combined factorization-solve, update, diagnostics, and the final fence.
The exact timing names are `instrumented_seconds_per_column`,
`setup_seconds_per_column`, `factorization_seconds_per_column`,
`solve_seconds_per_column`,
`combined_factorization_solve_seconds_per_column`,
`update_seconds_per_column`, `diagnostics_seconds_per_column`, and
`fence_seconds_per_column`.
Method-inapplicable linear fields follow the same blank-value contract as the
summary. Serialization occurs only after measurement and is outside every
timed region. The sidecar instrumented total must equal the corresponding
summary raw instrumented sample element by element; recomputing each phase and
fence median with the benchmark's exact median algorithm must reproduce its
summary value.

##### Transaction, provenance, and retained-run rules

A clean detached source commit, binary hash, locked-dependency provenance
hash, compiler/build/machine identity, exact commands, schemas, run identifier,
and output paths are bound in an immutable pre-run identity manifest. The
retained run identifier is
`task5-remediated-<UTC-YYYYMMDDTHHMMSSZ>-<12-character-source-SHA>`.
Before the manifest is sealed, the exact executable and dependency-provenance
record are copied exclusively into the new retained run's `inputs` directory;
the manifest and every shard command bind those retained copies. This prevents
later validation from depending on disposable scratch build products.
Summary and sidecar shards are written as separate partial files. After a
successful process exit, both are closed, validated, and hashed, then renamed.
An authoritative per-shard completion JSON containing both paths, hashes,
counts, run identifier, and identity-manifest hash is written and renamed
last. A shard is nonexistent for analysis unless that completion record and
both artifact hashes validate.

The full 24,300-row run uses 180-case shards aligned to the nine-row method
groups, producing 135 shard pairs. After all shard completion records pass,
the summary is merged in case-ID order and the sidecar in
`(case_id,instrumented_repetition)` order. A final immutable manifest records
the merged hashes and complete coverage. Any Task 5 result produced before
this amendment is historical pre-remediation evidence and may not be merged
into or substituted for the remediated dataset.

##### Amendment validation and failure behavior

Tests must prove exact occupancy sets and formulas, constant executed width,
unique composite keys and exact matrix marginals, method/phase-mode mapping,
blank-value rules, required CLI arguments, write/flush failures, invalid-row
sample absence, timed-but-ineligible sample retention, unique contiguous
sidecar repetition keys, elementwise summary/sidecar equality, and exact phase
median reconstruction. Old, unversioned, incomplete, mixed-run, stale-hash, or
completion-marker-free artifacts fail closed before analysis.

Fresh debug ASan/UBSan and optimized complete suites must pass, followed by a
full optimized numerical rerun with the exact timing binary, peer-style review,
and a second focused pass. Temporary worktrees and builds remain below
`/Users/burr114/e3sm-scratch/`; retained artifacts and provenance remain below
`/Users/burr114/e3sm-results/`. Actual matched-error, break-even, hybrid, and
recommendation calculations remain blocked until the remediated retained run
passes this contract.

The pre-retained implementation review on 2026-08-25 passed 35 focused Python
tests, all 13 registered debug ASan/UBSan tests with no sanitizer finding, all
13 optimized tests, and the explicit five-test optimized numerical gate. The
review also proved expected diagnostic text for every fail-closed CLI test,
canonical output-path separation, exact C++ median reconstruction, stopping-
policy consistency, and retained executable/provenance copying. These are
implementation-validation results only. A clean detached build, exact timing-
binary numerical rerun, fully paired 135-shard run, final manifest, and Task 6
analysis remain pending.

### Open issues for the proposed study

| Issue | Required decision or evidence | Interim treatment | Blocks |
|---|---|---|---|
| Numerical-method design | Approved by the user on 2026-07-23 | Follow the staged review gates below | Resolved for standalone prototype code |
| Bottom mass source | MAM4xx physics-owner decision on C++ versus intended formulation | Current C++ primary; note as active TODO | Production adoption only |
| Historical `top_lev=0` defect | Current PR #567 repair plus bounds and conservation regression evidence | Use the repaired `km1` expression | Supported-machine production integration |
| Historical disabled-mixing defect | Confirm current MAM4xx always-mixing behavior | Use current MAM4xx, where the partial-disable path was removed | Historical affected revisions only |
| Downstream fixer behavior | Trace whether EAMxx later masks or redistributes Forward Euler clipping errors | Report standalone Forward Euler raw and clipped ledgers with no fixer | Coupled scientific interpretation only |
| GPU implementation | Production-representative backend, solver, compiler, precision, and GPU | Defer until local screening supports it | GPU performance conclusion |
| Runtime distributions | Sampled coefficients and `nsubmix` from EAMxx | Synthetic typical/strong profiles | Production representativeness |
| Task 5 evidence remediation | Approved amendment and reviewed implementation are complete; a new retained run must still pass | Preserve existing Task 5 artifacts as historical; make no performance conclusion | Task 5 final evidence and Tasks 6--8 |

## Proposed Diagnostics for Substepping Experiments

This section defines a proposed diagnostics contract for the planned five-day
`ne256` experiments with `dt_phys=600 s` and
`mac_aero_mic:number_of_subcycles = 2, 4, 6, 12`. The corresponding ACI call
intervals are 300, 150, 100, and 50 s. This is a documentation and experiment
design update only; production instrumentation requires separate approval.

### Cost model and inner-loop work

Changing `mac_aero_mic:number_of_subcycles` changes the call frequency of every
process in the group. In a MAM4xx `pg2` configuration these are TMS, SHOC,
cloud fraction, MAM4xx ACI, and P3. The performance comparison must therefore
separate fixed per-call group cost from the adaptive `ndrop` mixing cost.

For physics step `p`, ACI call `r`, and column `c`, define

\[
W_{\mathrm{mix},p} =
\sum_{r=1}^{N_{\mathrm{MAC},p}}\sum_c
n_{\mathrm{submix},c,r}\,n_{\mathrm{active},c},
\qquad
n_{\mathrm{active},c}=pver-top\_lev.
\]

With four modes, 21 aerosol species-mass fields, and
`enable_aero_vertical_mix=true`, each active level in each explicit substep
executes 51 `explmix` updates: one droplet-number update, 25 cloud-borne aerosol
updates, and 25 interstitial aerosol updates. The source contains 56 nested
`TeamVectorRange` traversals and 24 explicit team barriers per substep. A useful
derived work proxy is therefore

\[
W_{\mathrm{explmix},p}=51W_{\mathrm{mix},p}.
\]

This proxy is not a replacement for measured GPU kernel time, but it separates
changes in adaptive work from fixed costs and is cheap to derive from actual
`nsubmix` values. The 56 `TeamVectorRange` traversals are nested inside the
per-column Kokkos kernel; they are not 56 independent GPU kernel launches.
Source-line or program-counter profiling is needed to apportion time within
that kernel.

### Provisional time-integrator performance model

This subsection compares the present Forward Euler subcycling with two
candidate implicit formulations. It is a source-based estimate, not measured
performance. It assumes the level-zero boundary defect is corrected, mixing is
enabled, and the estimated `max_tinv` ranges above represent the runtime
distribution.

The principal alternatives and their present status are:

| Method | Main advantage | Main disadvantage or risk | GPU work pattern | Current status |
|---|---|---|---|---|
| Current subcycled Forward Euler | Exact production behavior; very low overhead when `nsubmix` is one; repeated short steps can resolve the frozen operator accurately | Stability requires many sequential vertical passes in strong mixing; work varies by column; clipping can alter aerosol ledgers | Many fields can run in parallel, but levels and substeps are sequential and high-`nsubmix` columns create a long tail | Standalone gate passed independent review at the tested boundaries; no new production or GPU timing |
| IMEX-Euler | Removes the explicit mixing stability restriction; uses two shared scalar tridiagonal matrix families with many right-hand sides; fixed work per call | Explicit activation is lagged over the ACI interval; raw positivity can fail; solve overhead can exceed FE for `nsubmix` of one or two; first-order | EKAT cyclic reduction can expose parallelism within each column and across fields/columns | Standalone gate passed independent review; no production implementation or GPU timing |
| Fully implicit Backward Euler | Removes the mixing restriction and evaluates activation transfer at the new state; likely more robust when explicit-source lag matters | Requires custom 2-by-2 block-tridiagonal solves; more flops and workspace; first-order stiff damping can reduce accuracy; positivity is not guaranteed | Fixed work, but the coupled block solver is not supplied by the current EKAT scalar interface | Standalone host method-local gate passed; common Task 4 comparison, production implementation, and GPU timing remain separate work |
| Hybrid FE/IMEX | Retains the cheap FE path for low-`nsubmix` columns and targets IMEX at the expensive tail | Requires classification/batching or compaction, a threshold policy, two maintained paths, and validation at the switch boundary | Potentially best production GPU balance if high-`nsubmix` columns can be batched efficiently | Provisional concept only; requires profiling and a separate design gate |

The provisional GPU ranking is hybrid FE/IMEX if high-`nsubmix` columns can be
batched efficiently, otherwise IMEX for strong-mixing columns and current FE
for low-`nsubmix` columns. Fully implicit Backward Euler becomes preferable
only if evidence shows that IMEX source lag or positivity failures outweigh
its extra block-solve cost. This ranking is a source-based hypothesis, not a
measured GPU result.

Removing a stability restriction does not establish adequate accuracy. Several
short Forward Euler substeps can be more accurate than one first-order IMEX or
Backward Euler step over the full ACI interval. Performance decisions must
therefore compare runtime at matched error, not only runtime for equal outer
intervals, and must reject candidate timesteps that fail the raw-positivity or
ledger criteria.

All listed default horizontal resolutions use the same L128 grid. Approximate
layer thicknesses are 24-52 m in the lowest kilometre and 54-127 m from 1-2 km.
The stiffness therefore scales approximately as `K/dz^2`, or as `wtke/dz` on
the cloud-base override path, but does not increase directly with horizontal
resolution. Horizontal resolution changes the ACI interval and total number of
columns.

There are 51 vertically mixed fields: cloud-droplet number, 25 cloud-borne
aerosol number/mass fields, and the corresponding 25 interstitial fields. Many
fields share matrix coefficients. EKAT already provides GPU cyclic-reduction
and CPU Thomas solvers for a scalar tridiagonal matrix with multiple right-hand
sides (`externals/ekat/src/algorithm/ekat_tridiag.hpp`), and SHOC already uses
that facility. An optimized IMEX implementation can therefore factor or reduce
each shared matrix once per column and apply it to many aerosol fields.

For provisional estimates, one optimized IMEX call is assigned a cost of
**1.5-3 Forward Euler substeps**. A fully coupled Backward Euler call is assigned
**3-5 Forward Euler substeps** because the paired activation terms make the
system more complicated than EKAT's existing scalar multi-right-hand-side
case. The system is still linear under the frozen-operator contract; the added
cost comes from its 2-by-2 block structure, not a nonlinear iteration. These
ranges include matrix/RHS assembly and GPU synchronization, but not development
cost or changes elsewhere in ACI. A field-by-field implementation that does
not reuse matrices can be substantially slower.

The tables report estimated vertical-mixing-loop runtime relative to current
Forward Euler runtime; values below one are faster. They use a representative
`nsubmix` from each previously estimated range. The Forward Euler work column
gives total loop passes per physics step as range (representative value).

#### Typical turbulence

| Resolution | ACI interval (s) | `nsubmix` per call, range (representative) | Forward Euler passes per physics step | IMEX runtime / Forward Euler | Backward Euler runtime / Forward Euler | Expected result |
|---|---:|---:|---:|---:|---:|---|
| `ne4np4` | 150 | 1-6 (3) | 24-144 (72) | 0.50-1.00 | 1.00-1.67 | IMEX is near break-even to about 2x faster; full Backward Euler is unlikely to help. |
| `ne30np4` | 150 | 1-6 (3) | 12-72 (36) | 0.50-1.00 | 1.00-1.67 | Same per-call result as `ne4`; absolute work is spread over many more columns. |
| `ne120np4` | 75 | 1-3 (2) | 6-18 (12) | 0.75-1.50 | 1.50-2.50 | IMEX is borderline; full Backward Euler is likely slower. |
| `ne256np4` | 33.33 | 1-2 (1) | 6-12 (6) | 1.50-3.00 | 3.00-5.00 | The current one-step explicit path is likely faster. |
| `ne512np4` | 50 | 1-2 (1) | 2-4 (2) | 1.50-3.00 | 3.00-5.00 | The current one-step explicit path is likely faster. |
| `ne1024np4` | 50 | 1-2 (1) | 1-2 (1) | 1.50-3.00 | 3.00-5.00 | The current one-step explicit path is likely faster. |

#### Strong turbulence or cloud-base control

| Resolution | ACI interval (s) | `nsubmix` per call, range (representative) | Forward Euler passes per physics step | IMEX runtime / Forward Euler | Backward Euler runtime / Forward Euler | Expected result |
|---|---:|---:|---:|---:|---:|---|
| `ne4np4` | 150 | 9-21 (15) | 216-504 (360) | 0.10-0.20 | 0.20-0.33 | IMEX is about 5-10x faster; full Backward Euler is about 3-5x faster. |
| `ne30np4` | 150 | 9-21 (15) | 108-252 (180) | 0.10-0.20 | 0.20-0.33 | Same per-call speedup as `ne4`, with much larger domain-wide savings. |
| `ne120np4` | 75 | 5-11 (8) | 30-66 (48) | 0.19-0.38 | 0.38-0.63 | IMEX is about 2.7-5.3x faster; full Backward Euler is about 1.6-2.7x faster. |
| `ne256np4` | 33.33 | 2-5 (4) | 12-30 (24) | 0.38-0.75 | 0.75-1.25 | IMEX likely helps; full Backward Euler is borderline. |
| `ne512np4` | 50 | 3-7 (5) | 6-14 (10) | 0.30-0.60 | 0.60-1.00 | IMEX is about 1.7-3.3x faster; full Backward Euler is break-even to about 1.7x faster. |
| `ne1024np4` | 50 | 3-7 (5) | 3-7 (5) | 0.30-0.60 | 0.60-1.00 | Same per-call result as `ne512`; the total saving scales with the much larger column count. |

These are loop-level estimates, not whole-model speedups. If the measured
vertical-mixing loop fraction of baseline model runtime is `f` and a candidate
has table ratio `r`, the maximum whole-model speedup implied by this model is
`1 / ((1-f) + f*r)`. Implicit integration also replaces a variable-length
column loop with fixed solve work, which may reduce GPU tail imbalance in
heterogeneous strong-turbulence cases. Conversely, it imposes unnecessary
fixed work on the many columns for which `nsubmix` is one. A hybrid policy that
retains Forward Euler for `nsubmix <= 2` and batches larger-`nsubmix` columns for
IMEX is therefore a plausible performance target, but it requires profiling
and a separate numerical design review.

### Diagnostics currently present in the code

These diagnostics can be used without modifying production code. Some require
enabling an existing option or selecting the field in an output stream. Run
metadata and timers should be collected for all five-day performance runs;
three-dimensional scientific fields are better suited to a shorter companion
run because their I/O can materially change GPU performance.

| Category | Currently available diagnostic | Recommended use |
|---|---|---|
| Run provenance | `dt_phys`; configured `mac_aero_mic:number_of_subcycles`; the derived ACI interval `dt_phys/number_of_subcycles`; E3SM and MAM4xx revisions; compiler; Kokkos backend; GPU model; MPI task and thread layout | Record once per run. These are configuration/build metadata rather than history fields. |
| Existing performance timers | Generic EAMxx GPTL `::run` timers for `mac_aero_mic` and its children: TMS, SHOC, cloud fraction, MAM4xx ACI, and P3 | Retain total time and call count, then derive time per physics step and per child call. These timers separate process costs but do not isolate `dropmixnuc` or its explicit mixing loop. |
| ACI droplet response | `nc`, `nc_inp_to_aci`, `nc_nuceat_tend`, `nsource`, `ndropmix` | Compare end-state droplet number and activation/mixing response across subcycle counts. |
| ACI aerosol activation response | `ccn_0p02`, `ccn_0p05`, `ccn_0p1`, `ccn_0p2`, `ccn_0p5`, `ccn_1p0` | Retain in the diagnostic companion run to identify changes in activation environment. |
| Direct ACI turbulence input | `w_variance` | Retain profiles or selected levels. This, rather than SHOC prognostic `tke`, is the direct turbulence input used to construct ACI-local `wsub`. |
| SHOC and cloud context | `eddy_diff_heat`, `eddy_diff_mom`, `tke`, `pbl_height`, `cldfrac_liq`, `cldfrac_liq_prev`, `qc`, `p_mid`, `pseudo_density` | Diagnose turbulence, cloud base, liquid cloud, density, and vertical structure. `tke` is useful context but is not directly consumed by ACI. |
| Optional SHOC diagnostics | `brunt`, `shoc_mix`, `isotropy`, `shoc_cond`, `shoc_evap`, and the registered interface flux/moment fields | Available when `extra_shoc_diags=true`; use only in a short companion run. |
| Aerosol partitioning | Four interstitial modal-number tracers, four cloud-borne modal-number fields, and the paired interstitial/cloud-borne species-mass fields | Compute offline column burdens and interstitial-plus-cloud-borne (`I+C`) ledgers at output times. This cannot expose within-call clipping or transient imbalance. |
| Cloud-fraction response | `cldfrac_liq`, `cldfrac_ice`, `cldfrac_tot`, `cldfrac_ice_for_analysis`, `cldfrac_tot_for_analysis` | Compare cloud occurrence and phase partitioning across cases. |
| P3 state and bulk response | `T_mid`, `qv`, `qc`, `nc`, `qr`, `nr`, `qi`, `ni`, `qm`, `bm`, `precip_liq_surf_mass`, `precip_ice_surf_mass`, `precip_total_tend`, `nevapr`, effective radii, and `diag_equiv_reflectivity` | Diagnose downstream cloud and precipitation response. The two surface precipitation masses are accumulated fields. |
| Optional P3 process diagnostics | `qr2qv_evap`, `qi2qv_sublim`, `qc2qr_accret`, `qc2qr_autoconv`, `qv2qi_vapdep`, `qc2qi_berg`, collection/freezing/melting terms, and sedimentation terms | Available when `extra_p3_diags=true`; use in a short companion run to attribute response changes. |

The existing fields have three important limitations. First, normal output after
the process group generally sees the final MAC subcycle value, not a
physics-step integral or average. This is particularly important for
tendency-like ACI and P3 fields. Second, `nc_inp_to_aci` is registered with
tendency units even though it stores a number mixing ratio. Third, the six CCN
fields are registered with volume rather than inverse-volume units despite
concentration comments in the implementation. These metadata problems must be
documented or corrected before quantitative use.

At the baseline revisions, no registered field exposes actual `nsubmix`,
`max_tinv`, `dtmix`, the limiting level, ACI-local `wsub`/`wtke`, inner-loop
limiter activity, or a timer scoped specifically to `dropmixnuc`. The draft
diagnostic change now adds the actual count, its sum over a physics step, and
its maximum over a physics step. Those fields have not yet been collected in
an EAMxx simulation. The other quantities in this list remain unavailable.

#### Comprehensive existing-output recommendation for activation attribution

The following table is the recommended existing-output contract for a focused
EAMxx-versus-EAMv3 high-activation investigation. The focused stream should
cover approximately 6--24 simulated hours around the event. `1 nsteps`
samples the end of each outer atmosphere step; it does not expose individual
`mac_aero_mic` or internal `ndrop` subcycles. Configuration records and timers
are not history variables.

| Diagnostic set | Variables or records | Focused cadence | Focused averaging | Extended-run cadence and averaging |
|---|---|---:|---|---|
| Configuration and provenance | `dt_phys`, `mac_aero_mic::number_of_subcycles`, derived ACI timestep, E3SM revision, MAM4xx revision, compiler, Kokkos backend, GPU, MPI layout, and thread layout | Once per run | N/A | Once per run |
| Performance timers | `EAMxx::mac_aero_mic::run`, `EAMxx::tms::run`, `EAMxx::shoc::run`, `EAMxx::cld_fraction::run`, `EAMxx::mam4_aci::run`, and `EAMxx::p3::run` | Once per run segment | N/A | Once per run or simulation segment |
| ACI droplet response | `nc`, `nc_inp_to_aci`, `nc_nuceat_tend`, `nsource`, `ndropmix` | `1 nsteps` | `instant` | `nc`: hourly `instant`; tendencies: hourly `average` |
| CCN | `ccn_0p02`, `ccn_0p05`, `ccn_0p1`, `ccn_0p2`, `ccn_0p5`, `ccn_1p0` | `1 nsteps` | `instant` | Hourly `average` |
| Direct ACI turbulence input | `w_variance` | `1 nsteps` | `instant` | Hourly `average` |
| SHOC and cloud context | `eddy_diff_heat`, `eddy_diff_mom`, `tke`, `pbl_height`, `cldfrac_liq`, `cldfrac_liq_prev`, `qc`, `p_mid`, `pseudo_density` | `1 nsteps` | `instant` | Hourly `average`; optionally hourly `instant` state snapshots |
| ACI thermodynamic and vertical geometry | `T_mid`, `p_mid`, `p_int`, `pseudo_density`, `dz`, `z_mid`, `AtmosphereDensity` | `1 nsteps` | `instant` | Hourly `instant`; optionally hourly averages for context |
| Optional SHOC | `brunt`, `shoc_mix`, `isotropy`, `shoc_cond`, `shoc_evap`, `wthl_sec`, `wqw_sec`, `qw_sec`, `uw_sec`, `vw_sec`, `w3` | `1 nsteps` | `instant` | Hourly `average` |
| Aerosol modal size | `dgnum` | `1 nsteps` | `instant` | Hourly `average` |
| Interstitial aerosol number | `num_a1`, `num_a2`, `num_a3`, `num_a4` | `1 nsteps` | `instant` | Hourly `average` |
| Interstitial aerosol mass | `so4_a1`, `pom_a1`, `soa_a1`, `bc_a1`, `dst_a1`, `nacl_a1`, `mom_a1`, `so4_a2`, `soa_a2`, `nacl_a2`, `mom_a2`, `dst_a3`, `nacl_a3`, `so4_a3`, `bc_a3`, `pom_a3`, `soa_a3`, `mom_a3`, `pom_a4`, `bc_a4`, `mom_a4` | `1 nsteps` | `instant` | Hourly `average` |
| Cloud-borne aerosol number | `num_c1`, `num_c2`, `num_c3`, `num_c4` | `1 nsteps` | `instant` | Hourly `average` |
| Cloud-borne aerosol mass | `so4_c1`, `pom_c1`, `soa_c1`, `bc_c1`, `dst_c1`, `nacl_c1`, `mom_c1`, `so4_c2`, `soa_c2`, `nacl_c2`, `mom_c2`, `dst_c3`, `nacl_c3`, `so4_c3`, `bc_c3`, `pom_c3`, `soa_c3`, `mom_c3`, `pom_c4`, `bc_c4`, `mom_c4` | `1 nsteps` | `instant` | Hourly `average` |
| Cloud fraction | `cldfrac_liq`, `cldfrac_liq_prev`, `cldfrac_ice`, `cldfrac_tot`, `cldfrac_ice_for_analysis`, `cldfrac_tot_for_analysis` | `1 nsteps` | `instant` | Hourly `average` |
| P3 thermodynamic and hydrometeor state | `T_mid`, `qv`, `qc`, `qr`, `qi`, `qm`, `nc`, `nr`, `ni`, `bm` | `1 nsteps` | `instant` | Hourly `instant`; optionally hourly averages |
| P3 precipitation accumulators | `precip_liq_surf_mass`, `precip_ice_surf_mass` | `1 nsteps` | `instant` | Hourly `instant`; difference successive records for interval precipitation |
| P3 response rates | `precip_total_tend`, `nevapr` | `1 nsteps` | `instant` | Hourly `average` |
| P3 diagnostic properties | `eff_radius_qc`, `eff_radius_qi`, `eff_radius_qr`, `diag_equiv_reflectivity` | `1 nsteps` | `instant` | Hourly `average`; optionally hourly `instant` for joint distributions |
| Optional P3 processes | `qr2qv_evap`, `qi2qv_sublim`, `qc2qr_accret`, `qc2qr_autoconv`, `qv2qi_vapdep`, `qc2qi_berg`, `qc2qr_ice_shed`, `qc2qi_collect`, `qr2qi_collect`, `qc2qi_hetero_freeze`, `qr2qi_immers_freeze`, `qi2qr_melt`, `qr_sed`, `qc_sed`, `qi_sed` | `1 nsteps` | `instant` | Hourly `average` |

`p_int` and `dgnum` are persistent registered fields. `dz`, `z_mid`, and
`AtmosphereDensity` are case-sensitive derived diagnostics computed by the
EAMxx output system. The aerosol lists above contain only the configured
species-mode combinations; for example, primary-carbon mode 4 contains only
`pom`, `bc`, and `mom`. The spelling is `nacl_*`, not `ncl_*`.

#### Missing diagnostics required for causal attribution

The existing-output table can establish whether large `nc_nuceat_tend`
coincides with high CCN, high turbulence, rapid cloud growth, or altered
aerosol partitioning. It cannot by itself attribute an EAMxx--EAMv3 difference
to the activation-velocity calculation, modal activation, or internal mixing.
The following quantities are not directly available from current history
output.

| Missing diagnostic | Role in attribution | Current status |
|---|---|---|
| `wsub_raw` | Velocity derived from interpolated `w_variance` before the configurable floor | Internal value only |
| `wsub` | Profile passed from EAMxx into `dropmixnuc` | Internal work field only |
| `wtke` | Velocity actually passed into `activate_modal`, after the hard-coded `0.1 m s-1` floor | Internal work field only |
| `wsub_floor_active`, `wtke_floor_active` | Identify where either lower bound controls activation | Not implemented |
| `w_variance_int` or ACI-local `tke_int` | Test midpoint-to-interface interpolation and time-level differences | Internal work fields only |
| Mode-resolved `factnum` and activated-mass fraction | Separate available aerosol from the fraction activated in each mode | `factnum` internal; mass fraction local only |
| Mode-resolved `nact` and `mact` | Quantify activation rates supplied to explicit mixing | Internal work fields only |
| Maximum supersaturation | Directly distinguish dynamical forcing from aerosol limitation | Local variable in `activate_modal` |
| Mode-resolved `fluxn`, `fluxm`, and `flux_fullact` | Explain conversion of activated fractions into number and mass fluxes | Local variables only |
| `cloud_frac_used_by_ndrop` and `delta_cloud_frac` | Expose the cloud-water gate and the `0.01` cloud-growth thresholds actually traversed | Reconstructable, but not directly output |
| Explicit pre-ACI aerosol number, mass, and cloud-borne copies | Ensure comparison uses the state consumed by activation rather than a post-ACI state | Not implemented |
| Actual `nsubmix`, `max_tinv`, `dtmix`, and limiting level | Separate activation behavior from adaptive explicit-mixing work and timestep sensitivity | Draft fields expose actual `nsubmix`, its physics-step sum, and its physics-step maximum; the remaining quantities are internal only. |

Until these quantities are instrumented, an every-step existing-field stream
supports correlation and reconstruction but not definitive causal attribution.
In particular, offline reconstruction of `wsub` from `w_variance` and `dz`
cannot reveal all inner MAC calls, and it does not directly record where the
`wsubmin` and `wmixmin` floors were active.

#### Performance-only subset of existing diagnostics

For the current fixed-grid `ne256` subcycle comparison, use the following
reduced set before adding scientific-response diagnostics. The third row is
the minimum existing-field context while actual `nsubmix` and `max_tinv`
remain unavailable. All currently planned runs are five days; there is no
separate extended-run stream. Keep output identical across cases and retain
the I/O timers so its cost remains visible in the timing comparison.

| Diagnostic set | Variables/records | Five-day-run cadence | Averaging/summary |
|---|---|---:|---|
| Configuration and provenance | `dt_phys`, `mac_aero_mic::number_of_subcycles`, derived ACI timestep, E3SM and MAM4xx revisions, compiler, Kokkos backend, GPU and node counts, MPI/task/thread layout, output-stream configuration, and segment duration | Once per run | N/A |
| Application, process, and I/O timers | `EAMxx::run`, `EAMxx::mac_aero_mic::run`, child `::run` timers for TMS, SHOC, cloud fraction, MAM4xx ACI, and P3, plus `EAMxx::IO::standard` and named output-stream timers | End of run or timing segment | Retain total time, call count, time per call or physics step, and rank variability when available |
| Minimum existing `ndrop` driver context | `eddy_diff_heat`, `w_variance`; recommended cloud-regime context: `cldfrac_liq`, `cldfrac_liq_prev` | Hourly | `instant` |

`eddy_diff_heat` supplies the diffusivity that EAMxx interpolates to ACI
interfaces. `w_variance` is the direct turbulence input used to construct
ACI-local `wsub`, while current and previous liquid-cloud fraction control
cloud-base diffusivity replacement and data-dependent activation paths. These
fields still cannot reconstruct actual loop work: normal output sees the final
MAC-subcycle state, and temporal averages obscure the maxima and threshold
crossings that control `nsubmix`.

The explicit timestep-driving path is

```text
SHOC eddy_diff_heat(midpoint)
  -> EAMxx linear interpolation to kvh_int(interface)
  -> ndrop eddy_diff(k) = bound(kvh_int(k+1), 0.01, 100 m2/s)
  -> eddy_diff_kp(k), eddy_diff_km(k)
  -> tinv(k) = eddy_diff_kp(k) + eddy_diff_km(k)
  -> column reduction -> dtmix and nsubmix.
```

At diagnosed cloud base, `update_from_cldn_profile` replaces the bounded SHOC
value with `eddy_diff=wtke/zs`. Consequently, `eddy_diff_heat` is a genuine
upstream driver of the mixing timestep in ordinary layers but is not always
the final diffusivity controlling the column maximum. Geometry, pressure,
temperature, cloud-base state, and the level-zero indexing defect can also
alter the resulting rate.

EAMxx also supports `MAX` output streams, but these update a tally every model
step and maxima of separate fields need not be co-located or simultaneous.
They are not recommended for the current five-day timing matrix. Hourly
`instant` output provides limited context with lower tally overhead and must
not be used to reconstruct `max_tinv`.

Five simulated days is short, but hourly native-grid 3-D output at `ne256` can
still be large enough to affect runtime and storage. If necessary, use the same
coarsened remap or supported horizontal sampling for these four context fields
in every case, or reduce their common cadence. Preserve the native-grid process
timers and record the output configuration as run provenance.

Remove ACI droplet and CCN response, optional SHOC diagnostics, aerosol modal
size and partitioning, P3 state and precipitation, and the broader
thermodynamic/geometry set from the performance-only recommendation. For these
fixed-resolution experiments, geometry and active-level count are invariant
run metadata; their full 3-D output adds I/O without explaining differences
among subcycle cases. Restore those fields later for scientific-response and
conservation analysis.

### Diagnostics requiring code instrumentation

The draft pull request implements the smallest useful part of this section:
the exact `nsubmix` count and two physics-step summaries. The user approved
that diagnostic-only work on 2026-08-11. The additional fields and timers in
the tables below remain proposals and require their own implementation review.
GPU implementations should accumulate values on the device and write them
infrequently; copying every call to the host or adding an MPI reduction would
perturb the measurement.

#### Implemented in the draft diagnostic pull request

| Output field | Meaning |
|---|---|
| `ndrop_nsubmix` | Exact number of explicit mixing iterations used by the most recent MAM4xx ACI call in each column. |
| `ndrop_sum_nsubmix` | Sum of the exact counts over all MAM4xx ACI calls in the enclosing physics step, for each column. |
| `ndrop_max_nsubmix` | Largest exact count among those calls, for each column. |

MAM4xx returns the same integer that controls the production loop; EAMxx does
not reconstruct the value from saved inputs. EAMxx updates the three fields on
the device inside the existing column kernel. On the first MAM4xx ACI call of a
physics step, that same kernel starts a new sum and maximum; no separate reset
kernel is launched. The summaries remain available after the final call. No
per-call host copy or MPI reduction is added. The fields use the model's `Real`
type in memory, while normal history files default to single precision. To
preserve exact integer values in both places, an always-active device check
limits each count and the accumulated sum to `16,777,216` (`2^24`) in every
build. The unit test exercises the limit and rejection beyond it. EAMxx
initializes all three fields to zero. Before the first ACI call, zero means
"not yet calculated," not an executed zero-iteration loop.

The rebased MAM4xx diagnostic commit passed its focused `ndrop` unit test and
the complete 80-test `ndrop` comparison suite. Two existing maximum-
supersaturation comparisons were reported as skipped by their configured
tolerance rule; there were no numerical failures. The comparison scripts
require NumPy. One attempted rerun selected the system Python without NumPy
and stopped during test startup; the corrected rerun used the existing bundled
Python environment and passed. The EAMxx wrapper change has not been compiled
locally because this Mac is not a supported EAMxx test machine. Before merge,
EAMxx CI or a supported machine must exercise the device write, field output,
and reset behavior across multiple local and enclosing process-group subcycles
and across two physics steps. The same path must run on a GPU backend.

#### Minimum performance and adaptive-work set

| Proposed diagnostic | Units and required aggregation | Purpose |
|---|---|---|
| `MAMAci::call_function_dropmixnuc` profiling region | Total time, call count, mean and maximum time per call | Isolate the wrapper-plus-`dropmixnuc` region from ice nucleation, heterogeneous freezing, and state installation. Place the region around the existing wrapper call and its post-call fence; packing and interface interpolation remain included. |
| `ndrop_nsubmix` | Per column and ACI call; physics-step sum, maximum, mean, and histogram | Directly measure explicit-loop trip count. Do not reconstruct it from the intended algorithm. |
| `ndrop_sum_nsubmix` | `sum_r nsubmix(c,r)` per column per physics step, plus global sum | Measure total explicit-loop calls as ACI state and `max_tinv` evolve between MAC subcycles. |
| `ndrop_max_tinv` and `ndrop_no_tinv_above_threshold` | Maximum s^-1 among levels with `tinv>1e-6 s^-1`, or zero plus a true empty-reduction flag; physics-step maximum and histogram | Record the rate controlling `nsubmix` without misrepresenting an empty Kokkos reduction. |
| `ndrop_dtmix` | s per column and call; physics-step minimum, mean, and maximum | Record the actual explicit mixing time step after integer subdivision. |
| `ndrop_limiting_level`, pressure, and height | Zero-based index, Pa, and m at `argmax(tinv)`; physics-step extremum | Determine whether cost is controlled by cloud base, a boundary level, or another regime. |
| Runtime `ndrop_active_levels`, `top_level_mam4xx`, and `enable_aero_vertical_mix` record/assertion | Count, index, and Boolean | Verify actual loop extent and switch state rather than relying only on namelist provenance. |
| `ndrop_work_units` and `ndrop_explmix_updates` | `sum(nsubmix*nactive)` and `51` times that sum per physics step | Regress measured region time against expected explicit-loop work. |
| GPU profile metrics | Kernel/region time, launch count, occupancy, achieved memory bandwidth, and source-line or program-counter samples over a short no-I/O segment | Determine whether cost is controlled by memory traffic, barriers, occupancy, the explicit loop, or surrounding ACI work. Requires Kokkos profiling regions plus the platform profiler. |

Recommended `nsubmix` histogram bins are `1`, `2`, `3`, `4`, `5-8`, `9-16`,
`17-32`, and `33+`.

#### Criticality ranking for performance interpretation

For the fixed-grid `ne256` substepping experiments, the minimum defensible
performance instrumentation is smaller than the complete table above.

| Rank | Instrumentation | Decision enabled |
|---|---|---|
| Essential | A fenced profiling/timer region around `MAMAci::call_function_dropmixnuc`, retaining total time and call count | Determines how much runtime is actually in the wrapper-plus-`dropmixnuc` region rather than elsewhere in MAM4xx ACI or the MAC group. |
| Essential | Device-aggregated actual `nsubmix`: count, sum, maximum, mean, and the recommended histogram over all columns and ACI calls | Determines whether runtime changes are explained by adaptive explicit-loop trip counts and exposes a small expensive tail that a mean would hide. |
| Essential derived quantity | `Wmix=sum(nsubmix*nactive)` per output interval; `51*Wmix` may be reported as the expected number of `explmix` updates | Supplies the work variable for regression against measured region time. At fixed L128 and fixed `top_lev`, this is directly proportional to `sum(nsubmix)` and need not be stored independently. |
| Essential validity measurement | Matched instrumentation-off/on timing | Establishes whether counters, reductions, or output materially perturb GPU performance. |
| High value | `max_tinv`, empty-reduction flag, and limiting-level index, reduced to a maximum and histogram | Explains why `nsubmix` changed and detects whether the invalid level-zero boundary is controlling the result. This becomes explanatory rather than validity-critical after that defect is fixed. |
| High value, short run | GPU profiler metrics for the fenced region: elapsed kernel time, occupancy, achieved memory bandwidth, and source/PC samples | Explains departures from work-proportional timing, including memory limitation, barriers, poor occupancy, and long-tail column imbalance. Use a warmed-up no-I/O segment rather than the full five-day run. |
| Secondary | Limiting-level `eddy_diff_kp`, `eddy_diff_km`, diffusivity override/bound counters, and ACI-local `wsub`/`wtke` | Explains the physical and numerical cause of large `max_tinv`; collect only if the essential diagnostics show unexpectedly large or changing `nsubmix`. |
| Physics/numerics rather than performance | Activation-limiter counts, positivity-clipping amounts, paired ledgers, and physics-step tendency accumulation | Needed to assess scientific integrity and timestep sensitivity, but not required to attribute the first-order runtime cost. |

`ndrop_dtmix` is redundant for performance accounting because it is exactly
`dt_aci/nsubmix` after integer subdivision. A separate stored value is useful
only as a consistency check. Similarly, per-column/per-call history output is
not required for the main five-day runs: device-resident interval sums,
maxima, counts, and histograms preserve the relevant performance information
with much lower overhead.

#### Internal drivers and numerical-behavior set

The primary timing stream should retain reduced values at the limiting level.
Full three-dimensional fields belong in a shorter diagnostic companion run.

| Proposed diagnostic | Units and aggregation | Interpretation |
|---|---|---|
| `eddy_diff_kp`, `eddy_diff_km`, and their sum at the limiting level | s^-1 at `argmax(tinv)`, with optional sampled profiles | Separate lower- and upper-interface contributions to the stability limit. |
| `eddy_diff` before and after the cloud-base override | m2 s^-1 at the limiting level, plus override count | Distinguish SHOC diffusivity from `eddy_diff=wtke/zs` imposed at diagnosed cloud base. |
| `eddy_diff` low/high bound counts | Counts at 0.01 and 100 m2 s^-1 | Identify whether hard diffusivity bounds control the result. |
| `zn`, `zs`, `csbot`, `pdel`, and `dz` at the limiting level | m2 kg^-1, m^-1, kg m^-3, Pa, and m | Diagnose vertical-resolution and density scaling in the interface rates. Although pressure thickness and geometry exist elsewhere, the reduced values at the `ndrop` limiting level require instrumentation. |
| ACI-local `wsub` and `wtke` | m s^-1; sampled profiles and limiting-level values | Connect existing `w_variance` to activation velocity and the cloud-base diffusivity override. |
| `overlapp` and `overlapm` | Nondimensional; sampled profiles or limiting-level values | Identify overlap regimes that alter activation and mixing. |
| `nact` and `mact` maxima and limiter counts | s^-1 by mode; maximum and count clipped to `eddy_diff_kp` | Determine whether activation-rate limiting changes as ACI timestep changes. |
| Positivity-limiter count and clipped amount | Count and # kg^-1 or kg kg^-1 by droplet, interstitial, and cloud-borne category | Detect clipping before `max(qnew,0)` changes conservation or scientific response. |
| Paired aerosol ledger residual | Per constituent, or maximum relative change in `sum pdel/g*(I+C)` | Detect within-call transfer inconsistency and clipping-related nonconservation. |
| Physics-step accumulated ACI/P3 tendencies | Time-weighted mean over all MAC subcycles | Make tendency comparisons independent of which subcycle happened to run last. This requires accumulation support even though the underlying instantaneous fields already exist. |
| Instrumentation overhead | Matched instrumentation-off/on region and application timing | Quantify whether the diagnostics themselves perturb GPU performance. |

### Aggregation and output cadence

Normal history output after `mac_aero_mic` sees the state or computed field from
the final inner call. For a subcycling experiment, a final-call tendency is not
the physics-step-integrated effect. For any tendency-like quantity `X`, retain
the time-weighted physics-step mean

\[
\overline{X}_p = \frac{1}{\Delta t_{phys}}
\sum_{r=1}^{N_{MAC}} \Delta t_{aci,r} X_r.
\]

Use a sum for work and timer counts, a maximum for stability limits and limiter
extrema, and an end-of-physics-step snapshot for prognostic state. The current
experiment phase consists only of the matched five-day performance runs. Each
case should use identical hourly `instant` output for `eddy_diff_heat`,
`w_variance`, `cldfrac_liq`, and `cldfrac_liq_prev`, plus end-of-run timer
summaries. Separate extended runs, every-step companion output, focused
inner-call windows, and full-run profiling are deferred.

After the draft count fields pass EAMxx testing, their per-column sum and
maximum can be sampled with hourly `instant` output in the same five-day
matrix. Each hourly record represents the completed physics step at the output
time; it does not accumulate counts over the preceding hour. Histograms and
work counters are still unimplemented. A short warmed-up profiler segment or
matched instrumentation-off/on segment remains the recommended follow-up for
quantifying kernel behavior and diagnostic overhead, but it is not part of the
current run matrix.

### Validity gates for the planned runs

- Use a MAM4xx revision containing PR #567 and verify the `top_lev=0` boundary
  test before interpreting `ndrop` scaling. The diagnostic branch satisfies
  the revision requirement; supported-machine EAMxx testing remains pending.
- Record actual `nsubmix`. The inspected Kokkos reduction overwrites the
  preassigned `dtmicro` value, so weak mixing can produce one explicit substep.
- Do not use the older `enable_aero_vertical_mix=false` path as a validation
  case. Current MAM4xx has removed that operational switch and always performs
  the mixing update.
- Keep output fields, cadence, remapping, compression, MPI/GPU layout, and
  profiler settings identical within each timing comparison.
- Exclude initialization and profiler warm-up from steady-state timing, and use
  repeated timing samples to report run-to-run variability.

## Production-State Replay and Method-Selection Amendment

> **[APPROVED - 2026-09-14]** The user approved the immutable document
> snapshot with SHA-256
> `d79a14253ca3acb5108d54f37b3ddcc9aaa2a9d9c4ecef391721c27868df2138`
> for the transition `AMENDMENT_PENDING` to `AMENDMENT_APPROVED`. That
> approval authorizes only preparation of the serialized-schema and
> C++-interface addendum below.

The document has three distinct status classes:

| Status | Scope and authority |
|---|---|
| **Historical/executed evidence** | The defect investigations, analytic checks, 1,875-case synthetic Cartesian gate, method-local gates, and recorded completed validation runs remain immutable evidence of what was executed on their recorded stacks. They are regression evidence, not production method-selection evidence, and grant no new execution authority. |
| **Previously approved standalone work** | The existing test-only FE, IMEX, BE, reference, benchmark, and evidence-transaction work remains approved within its recorded local standalone scope. It does not authorize production instrumentation, a retained Task 5 run, supported-machine work, or a production recommendation. |
| **Approved production-state replay amendment** | The immutable approved snapshot supersedes earlier synthetic performance-screening and diagnostic plans where they conflict and freezes this architecture. Its approval authorizes preparation of the exact interface/schema addendum only. After that addendum receives its required follow-on approval, it may authorize local implementation and testing of the offline capture-schema validator, exact locked-source replay adapter, census/selection logic, replay benchmark, and analysis rules described here. Amendment approval does not authorize a production-source edit, production capture, access to prior retained Task 5 results, a supported-machine or scheduler command, a GPU run, publication, or a production-method recommendation. Each remains blocked by the approvals listed below. |

Nothing in this amendment changes an operator equation, limiter, boundary,
source sign, state transition, candidate method, or numerical acceptance gate.
The analytic solutions and 1,875-case gate remain mandatory regression gates,
but cannot select a production method.

Validation remains **STRICT** because the work compares numerical methods for
production aerosol--cloud microphysics. A method preference requires executed
accuracy, robustness, and timing evidence from the same provenance-locked
production-state replay population. Until physics owners provide a production
accuracy budget, results are Pareto/status evidence rather than a claim that
one method is the "fastest acceptable" method.

### Evidence hierarchy and locked stack

Evidence is interpreted in this order:

1. analytic and independently assembled numerical references establish
   correctness and convergence;
2. the existing 1,875-case synthetic matrix establishes deterministic
   regression coverage across controlled regimes;
3. exact production-boundary captures establish state frequency, operator
   frequency, Forward Euler work, and same-state method behavior;
4. primary uninstrumented timings establish standalone cost on that replay
   population; and
5. a whole-ACI production-layout test is required for any production or GPU
   recommendation.

The locked physics target is MAM4xx
`0d076d8ebb004a18f05f51f97b08edd9def6826f`, with EKAT
`3a1485aecfed7a036e4fa30ce845498fd37975e2` and the Kokkos revision recorded in
`dependency_lock.json`. Capture-enabled evidence must separately bind the exact
instrumentation-only descendant or canonical patch hash, the actual
post-implementation E3SM source commit and tree, the approved changed-file
allowlist, the complete canonical base-to-instrumented diff hash, clean source
trees, the executable and build identity, and unchanged dependency gitlink
names, modes, and object IDs. Provenance validation must prove ancestry from
the locked target or clean application of the canonical patch to that target,
and must fail closed on an extra changed path, dirty tree, diff or patch
mismatch, source mismatch, or changed gitlink. Separate otherwise-valid Git or
manifest provenance does not establish observational equivalence.

The uninitialized historical MAM4xx gitlink in this worktree must not be
treated as the production source and must not be mixed into a capture, replay,
build, or conclusion. No gitlink may be changed to construct an instrumented
stack. Temporary dependencies and builds belong below
`/Users/burr114/e3sm-scratch/`; retained evidence belongs below
`/Users/burr114/e3sm-results/`.

### Target-model declaration

The benchmark targets the locked current MAM4xx **C implementation** at
`0d076d8ebb004a18f05f51f97b08edd9def6826f`, including its current bottom
activation/source branch and its exact EAMxx coupling inputs. Any difference
between that C branch and the EAMv3 Fortran bottom-source formulation is a
separately labeled diagnostic only. It must not be corrected, blended into the
reference target, or used to prefer a candidate integrator unless a separate
atmospheric-physics design change is documented and explicitly approved.

For this replay target, an admissible record satisfies
`0 <= top_lev < pver`. Physical exchange uses the locked no-flux top and bottom
boundaries. The first-active-level controller contribution may affect only the
production Forward Euler `nsubmix` selection; it remains absent from physical
tendencies, ledgers, and implicit matrices. Adding or changing production
configuration validation, a production API range guard, or boundary physics is
outside this amendment and remains blocked by a separate approved design
change.

### Exact production capture boundary

Capture is a two-phase transaction. Phase 1 occurs inside `dropmixnuc`,
immediately before its first call to `update_from_explmix`. At that point all
of the following production work must already have completed for the current
ACI call:

- activation diagnosis;
- cloud-growth, cloud-decay, and cloud-disappearance handling;
- cloud-base diffusivity overrides;
- cloud-overlap calculation; and
- number and mass activation-rate limiting.

Phase 1 must make a capture-owned, immutable deep copy of the common old
droplet buffer, both active interstitial and cloud-borne aerosol buffers, and
all required inputs before any `update_from_explmix` mutation. The copy must
complete using backend-appropriate synchronization or an existing equivalent
completion guarantee before production mutates the source storage. A retained
alias, pointer, or view into production storage that may subsequently be reused
or overwritten is prohibited. The capture-disabled path may not acquire a new
copy, fence, event, or synchronization requirement.

The Phase 1 payload contains every value consumed by the production update and
substep controller: `dtmicro` [s], `top_lev` [zero-based level index], current
liquid-cloud fraction `cldn` [1] with full `pver` extent, `eddy_diff`
[m2 s-1], `zn` [m2 kg-1], `zs` [m-1], `csbot` [kg m-3], both directional
cloud overlaps [1], post-limiter `nact` and `mact` [s-1], `pdel` [Pa], `dz`
[m], pressure [Pa], height [m], and source identifiers. `cldn` must have the
recorded vertical orientation, contain only finite values in `[0,1]`, and be
covered by the payload checksum. Source identifiers
must resolve the case/run, ensemble member if any, physics step, MAC/ACI
subcycle, column/global column, capture sequence, grid, level count, and the
exact E3SM/MAM4xx/EKAT/Kokkos revisions and configuration hashes. The schema
must also record buffer identity/parity, scalar type, array extents, vertical
index orientation, units, and checksums.

Phase 2 occurs at the caller only after the complete production `dropmixnuc`
call and all of its substeps return. Return alone is insufficient when relevant
backend work remains queued: every producer of the actual returned production
`nsubmix`, final buffer parity, completion metadata, and any final-state payload
required by the approved contract must have completed through
backend-appropriate synchronization or an existing equivalent completion
guarantee before the value is read or copied. Phase 2 values must be preserved
in capture-owned storage before their production storage can be reused or
overwritten, and must be attached to the same immutable source identity and
initial parity recorded by Phase 1.

Validation, hashing, atomic promotion, and sealing occur only after both phases
match and all completion, ownership, extent, orientation, unit, finite-value,
range, and checksum requirements pass. Missing, nonfinite, inconsistent,
out-of-range, partial, mismatched, or failed records are `UNAVAILABLE`, not
replay cases.

Production captures use the realistic full L128 nonuniform pressure/height
geometry. L72 records are separately labeled compatibility evidence. An
arbitrary prefix, suffix, or truncation to another `nlev` is not a physical
grid and cannot enter production conclusions.

### Exact operator reconstruction gate

Replay must reconstruct the actual production formulas from the captured raw
inputs. It must not replace them with `K/dz^2`, assume `zs=1/dz`, impose a
uniform grid, or substitute the standalone harness's simplified exchange
coefficient. For active level `k`, let the production implementation define
the directional rates

\[
r_k^{\uparrow}=R^{\uparrow}(eddy\_diff, z_n,z_s,csbot,pdel,dz,p,z)_k,
\qquad
r_k^{\downarrow}=R^{\downarrow}(eddy\_diff,z_n,z_s,csbot,pdel,dz,p,z)_k,
\]

where `R` denotes the exact locked-source expression, including cloud-base
overrides and the no-flux boundary. These rates have units s\(^{-1}\). The
captured overlaps multiply only the production terms to which the locked
source applies them. The post-limiter `nact` and `mact` rates must retain their
number/mass and mode mappings; they may not be recomputed from an idealized
activation fraction.

Before any candidate comparison, reconstructed values must agree with
production for every captured level and field in:

- both directional physical rates and overlap factors;
- each post-limiter number and mass activation rate;
- the identity and value of the limiting level;
- the controller-only first-active-level contribution; and
- the exact integer `nsubmix`.

Agreement is bitwise where the replay uses the same scalar type and operation
order. Otherwise the test must set a documented roundoff-only tolerance from
an operation-order analysis and demonstrate that it cannot alter the limiting
level or `nsubmix`. Any mismatch is `UNAVAILABLE` for comparison and fails the
operator-replay gate. The controller contribution remains absent from all
physical tendencies and implicit matrices.

For diagnostic separation, define

\[
\chi_{mix}=\Delta t_{aci}\max_k(r_k^{\uparrow}+r_k^{\downarrow}),
\]

\[
\chi_{act}=\Delta t_{aci}\max_{k,m}
\{nact_{m,k},mact_{m,k}\},
\]

and

\[
\chi_{FE}=\Delta t_{aci}\max_k
\left(r_k^{\uparrow}+r_k^{\downarrow}
+\mathbf{1}_{k=t}r_{top,controller}\right).
\]

All three are dimensionless and must be retained separately. For the current
production selector, the replay must independently recover

\[
nsubmix=\left\lfloor\frac{\chi_{FE}}{0.9}\right\rfloor+1
\]

subject to the locked source's threshold and representability branches. The
separate diagnostics prevent activation stiffness or the controller from being
hidden inside a mixing-only label.

### Frozen-state comparison contract

Production-count FE, one-step IMEX Euler, and one-step BE start from the same
captured active buffers, use the same frozen operator inputs, and advance the
same complete ACI interval `dtmicro`. FE includes its production `nsubmix`
selection and all production substeps. IMEX and BE each take exactly one step.
All methods retain the previously approved boundaries, source signs, limiter
semantics, zero-cloud transition, state-field maps, residuals, and aerosol
ledger diagnostics. Direct IMEX or BE factorization, pivot, solve, finite,
positivity, branch, or residual failure remains fail-closed; no clipping or
fallback may convert it into a passing direct-method result.

The independently assembled Dormand--Prince 5(4) physical-RHS reference is
run at `rtol=1e-12` and `rtol=1e-13`, with the existing field-scaled absolute
tolerances. The maximum fieldwise normalized disagreement must be no greater
than `2e-11`. A case failing this self-consistency check is
`REFERENCE_UNRESOLVED`; none of its candidate errors may enter a preference.
The existing analytic and synthetic gates remain mandatory before production
replay, but passing them cannot substitute for this per-state reference gate.

The provisional repository-documented ne256 MAC/ACI intervals are 50, 100,
150, and 300 s. They may label proposed replay/sensitivity strata, but must not
be described as executed until an authoritative sensitivity-run manifest binds
each interval to an actual run. Use 25 s only for selected convergence checks,
not as a production-frequency stratum.

### Census, stability pressure, and censoring

Before census capture, its approved population manifest must freeze the cases,
ensemble members, time windows, spin-up exclusions, timesteps, columns,
MAC/ACI calls, cadence, expected counts, and exclusion rules. Adaptive thinning
is prohibited. The sealed census must account for every expected record, with
zero unaccounted records. If an unavailable record cannot be assigned its
source stratum or weight, method-selection aggregation remains blocked until an
approved amendment resolves the population accounting.

The source census must account for every eligible production capture before
sampling or benchmarking. For each source stratum it records case frequency,
the histogram of exact production `nsubmix`, and cumulative estimated FE work
(at minimum the sum of production substeps, with any richer work weight defined
and frozen before analysis). FE stability pressure is classified solely by
computed production `nsubmix`:

| Class | Production `nsubmix` |
|---|---:|
| weak | 1 |
| moderate | 2--4 |
| strong | 5--16 |
| near-stiff | 17--64 |
| extreme tail | greater than 64 |

Every class is cross-tabulated with `chi_mix`, `chi_act`, `chi_FE`, activation
limiter activity, source frequency, and cumulative FE work. Core benchmarking
covers exact production counts through 64. The extended tail covers 65--512.
Counts above 512 are `RESOURCE_CENSORED`, not FE numerical failures; they
receive no FE/candidate timing or error ratio. No aggregate recommendation is
permitted if censored records exceed 1% of eligible source cases or 1% of
estimated source FE work. The extreme-tail materiality decision is therefore
an output of the complete census, not an assumption.

### Pilot selection and preregistered expansion

Replay selection operates only after the complete source census and record
validation. A stratum is the Cartesian key
`timestep x FE-regime x activation-regime`, where timestep is the captured
`dtmicro`, FE regime is the `nsubmix` classification above, and activation
regime is exactly one of:

- `inactive`: `chi_act=0`;
- `active_unlimited`: `chi_act>0` and no number or mass activation limiter is
  active; or
- `limiter_active`: at least one number or mass activation limiter is active.

The pilot contains one deterministic full-state production replay from every
populated replay-eligible stratum. The candidate pool contains only available
production records with `1 <= nsubmix <= 512`; a stratum populated only by
records above 512 remains census evidence and has no replay case. Within a
stratum, selection is the valid record with the smallest cryptographic hash of
`selection_seed || immutable_source_id`; the
hash algorithm, seed, source-ID serialization, census manifest hash, and
tie-break rule are frozen in a preregistration artifact before any FE, IMEX, or
BE candidate result is inspected.

Every populated observed extreme-tail stratum with `nsubmix` from 65 through
512 is represented under the same rule. A production tail record is never
replaced by a synthetic profile. At most two synthetic or perturbed stress
sentinels may be added only when an operator boundary or numerical failure
branch lacks an observed case. Sentinels are labeled `STRESS_SENTINEL`, carry
zero census weight, are excluded from source-frequency, FE-work, timing-ratio,
and recommendation aggregates, and serve only robustness diagnosis.

After pilot results are available, a populated production stratum may expand
to at most three total replay cases. Additional records are selected by the
same next-smallest frozen hash, never by favorable candidate behavior, and
only when at least one preregistered trigger fires:

1. **Within-stratum variation:** the census range or preregistered quantile
   spread of `chi_mix`, `chi_act`, `chi_FE`, operator norm, cloud-overlap
   occupancy, or state magnitude crosses its frozen material-variation
   threshold.
2. **Reference certification:** the pilot is `REFERENCE_UNRESOLVED`, or a
   second certified record is required to distinguish an isolated reference
   failure from a stratum-wide limitation. An unresolved record is retained
   and is never replaced silently.
3. **Source/work materiality:** the stratum crosses a frozen share threshold
   for source-case frequency or estimated FE work in the complete census.
4. **Decision sensitivity:** a preregistered leave-one-stratum-out or
   worst-case perturbation calculation shows that the stratum can change a
   reported pairwise status or Pareto conclusion.

The exact numeric trigger thresholds, quantile definitions, materiality
weights, and decision-sensitivity calculation must be frozen with the seed and
census hash before candidate results are inspected. Expansion decisions and
their triggering values are recorded even when no additional eligible record
exists. No discretionary expansion, post-result threshold change, or more
than three production cases per populated stratum is permitted without a new
approved amendment.

Before candidate execution, the approved preregistration must bind the exact
sealed census-manifest hash and exact source/tree and executable hashes for the
replay harness, FE, IMEX, BE, reference implementation, timing driver, and
analysis implementation. Any change to one of those identities invalidates the
preregistration and requires renewed approval before candidate execution or
inspection.

### Accuracy, robustness, timing, and decision statuses

Each method/case must first pass the approved finite-value, raw positivity,
conservation/ledger, physical-residual, boundary, branch, and reference gates.
Report both normalized pressure-weighted L1 and Linf error. Timing comparisons
use only primary uninstrumented complete-advance samples; phase-instrumented
samples diagnose cost composition and instrumentation overhead but cannot
enter a speed ratio or preference.

Timing uses a specified paired procedure. One timing block contains one
primary uninstrumented complete-advance observation for every compared method
on the same replay case, with interleaved method order. Only a replay for which
both methods, the operator gate, and the reference gate are eligible enters a
pair. For methods `A` and `B`, bootstrap replicate `b` resamples complete
matched timing blocks within each replay case and computes

\[
\rho_{i,b}=\frac{\operatorname{median}T_{A,i,b}}
                  {\operatorname{median}T_{B,i,b}}.
\]

It then resamples selected production cases as complete matched units within
each stratum and combines case ratios through the weighted log ratio

\[
\rho_b=\exp\left(
\frac{\sum_i w_i\log\rho_{i,b}}{\sum_i w_i}\right).
\]

Use exactly 10,000 hierarchical paired-bootstrap replicates, a preregistered
seed, and the percentile 95% interval. Preserve pairing, strata, and all method
values when resampling; never resample method columns independently. Report
source-frequency weights and estimated-FE-work weights as two separate
aggregates. Within each aggregate, every selected production case receives its
stratum's frozen census weight divided by the number of selected production
cases in that stratum. A timing advantage for `A` requires the entire interval
for `rho=T_A/T_B` to be below one; an interval containing one is timing-
indeterminate. Report per-case ratios as descriptive values. Ratios are never
formed for a missing, failed, unresolved-reference, sentinel, or censored
case.

The approved interface/schema addendum must freeze the exact bootstrap random
number generator and deterministic stream-use order, seed representation and
interpretation, with-replacement draw counts at both resampling levels, the
even-sample median, percentile endpoint/interpolation rule, weighted-quantile
ordering and interpolation convention, and the handling of one-case strata and
missing or insufficient eligible pairs. The census-bound preregistration must
instantiate those rules with the exact seeds and must freeze the numerical tie
tolerance. No omitted convention may be supplied by an implementation default.

Numerical errors are deterministic outcomes of each certified replay and do
not receive bootstrap confidence intervals. For each method and for both
normalized pressure-weighted L1 and Linf, report the complete per-case values,
the source-frequency-weighted and FE-work-weighted arithmetic means, weighted
median, weighted 90th and 99th percentiles, and maximum, separately by stratum
and over the eligible population. Weights come only from the frozen complete
census; stress sentinels have zero weight. Zero errors remain zero and are not
regularized to manufacture ratios.

For deterministic pairwise accuracy classification, `A` is _no worse_ than
`B` only when both error norms for every matched eligible replay satisfy the
preregistered numerical tie tolerance relative to `B`. It is accuracy-
favorable when it is no worse and is strictly better in at least one norm on at
least one production replay; it is an accuracy tie when both methods are no
worse; otherwise the errors are crossed. This casewise rule prevents a weighted
mean from hiding a tail regression. Weighted summaries describe population
magnitude and materiality but cannot overturn the casewise classification.

A method `A` `DOMINATES` `B` only when `A` passes all robustness gates, is
accuracy-favorable or tied, and has a favorable paired-bootstrap timing ratio
under both source-frequency and FE-work weighting. Opposing speed and accuracy
directions, or conflicting source/work-weighted timing conclusions, are
`TRADEOFF`. An interval containing one, insufficient eligible pairs, or a
casewise relation that establishes neither dominance nor a directional
tradeoff is `INDETERMINATE`. Gate and availability precedence remains
`UNAVAILABLE`, `NUMERICAL_FAILURE`, `REFERENCE_UNRESOLVED`, or
`RESOURCE_CENSORED` as applicable. Same-timestep ratios remain descriptive.
Without a physics-owner accuracy budget, aggregation may report Pareto fronts
and status frequencies but not "fastest acceptable method."

### Instrumentation tiers

| Tier | Required content | Use and exclusion |
|---|---|---|
| **1. Capture/replay correctness and census** | Immutable source identifiers and stack/configuration hashes; capture role; capture-owned initial and final buffer identity/parity and all droplet/aerosol state required by the approved contract; `dtmicro`, `top_lev`, current `cldn` [1] with full `pver` extent, recorded orientation, finite `[0,1]` validation, and checksum coverage; actual returned production `nsubmix`; `eddy_diff`, `zn`, `zs`, `csbot`, overlaps, post-limiter `nact`/`mact`, `pdel`, `dz`, pressure, and height; reconstructed directional rates, observed `dtmin`/`nsubmix`, reconstructed selector maximum, candidate count and complete minimizing-level tie set, controller contribution, `chi_mix`, `chi_act`, `chi_FE`; schema/extents/units/checksums; source frequency and cumulative FE-work counters. | Required to admit a capture, reconstruct the exact operator, form strata, freeze selection, and complete the census. Phase 1 comprises the approved P1a/P1b/P1c subcheckpoints; Phase 2 follows the completed production call and owns returned values before reuse. Missing or inconsistent content fails closed as `UNAVAILABLE`. |
| **2. Timing validity** | Primary instrumentation-off complete-advance samples; matched case and timing-block identifiers; interleaved method order; warm-up/sample counts; fences required by the measured algorithm; checksums preventing dead-code removal; executable/build/compiler/machine identity; paired-bootstrap seed and weights; a matched instrumentation-on overhead run kept in a separate artifact. | Required for any timing ratio. Only primary uninstrumented samples enter ratios, bootstrap intervals, Pareto status, or recommendations. Instrumented phase values and overhead runs are diagnostic only. |
| **3. Deferred explanatory diagnostics** | Per-phase setup/factorization/solve/update/fence timings; GPU kernel and memory profiling; packing, compaction, dispatch, workspace, occupancy, and launch diagnostics; broader SHOC/cloud/ACI output fields; diffusivity-bound and override profiles; expanded causal-attribution tendencies and profiler traces. | Deferred until Tier 1 census/replay and Tier 2 timing justify them, or until the whole-ACI GPU gate is separately approved. They may explain behavior but cannot repair a correctness gate or enter primary speed ratios. |

The contract-only approval of document SHA-256
`a829ffc01d09a9960e9bc5df04a15199a39833fd6942cb696a8c942bb1292c38`
replaced the former scalar “limiting level” requirement with the exact OI-1
evidence contract shown above. This approval changes no other instrumentation
tier and grants no implementation or execution authority.

Every capture artifact and completion manifest has one immutable checksummed
role, exactly `qualification` or `census`. Roles cannot be relabeled or
promoted. A dedicated qualification validator accepts only `qualification`
artifacts and may compute only observational-equivalence and
instrumentation-overhead evidence. Census, replay, selection, primary timing,
error aggregation, and method-selection consumers must reject every
`qualification` artifact. Qualification records have zero census weight and
are inadmissible for timing, error, Pareto, or method-selection conclusions.

### Whole-ACI production-layout and GPU gate

Standalone CPU or GPU replay may select a prototype for further testing only.
A production recommendation additionally requires a whole-ACI,
production-layout GPU experiment that includes production packing and
unpacking, heterogeneous columns and `nsubmix`, method dispatch/compaction if
used, required fences, workspace allocation/reuse, and instrumentation
overhead. It must demonstrate unchanged downstream scientific state and
tendencies at the approved comparison tolerance, not merely kernel agreement.
The primary comparison remains an uninstrumented run; a matched instrumented
run measures overhead. Supported-machine resources, exact scheduler commands,
project/account, queue or QOS, resources, walltime, paths, and explicit user
approval are required before capture or GPU execution.

### Interfaces, tests, and failure behavior

The proposed implementation introduces four test/evidence interfaces, whose
exact serialized schema and C++ signatures must be presented for renewed
approval before source changes:

1. a production-boundary capture writer receiving read-only production views
   and scalar/source metadata;
2. a validator/reader producing one immutable L128 or labeled L72 replay
   record; and
3. an exact locked-source operator adapter producing physical rates, limiting
   metadata, controller contribution, and production `nsubmix` for the existing
   FE/IMEX/BE harness; and
4. a dedicated qualification validator that accepts only `qualification`
   artifacts and emits only equivalence and instrumentation-overhead evidence.

Capture must be observational: it may not alter view ownership, buffer parity,
fence placement in the uninstrumented control, branch decisions, or downstream
science. Qualification requires three supported-backend arms with identical
inputs, compiler configuration, backend, layout, and seeds: the locked base
build without instrumentation; the instrumentation descendant with capture
disabled; and the same descendant with capture enabled. Base versus disabled
must establish that the source edit is inert. Disabled versus enabled must
check buffer parity, branch decisions, actual `nsubmix`, state, tendencies,
downstream science, synchronization, and the approved instrumentation-overhead
bound.

The observational-equivalence approval record must bind the exact tested
instrumentation descendant or canonical patch hash, complete diff hash and
approved changed-file allowlist, tested executable/build identity, relevant
backend, compiler configuration, layout, seeds, capture settings, qualification
manifest hash, and three-arm evidence hashes. The addendum must specify the
backend-specific Phase 1 and Phase 2 completion mechanism and freeze the exact
changed-file allowlist before the production capture-writer edit is approved.

Required tests include schema round trips and corruption failures;
unit/shape/range/nonfinite guards; exact buffer and source-identifier binding;
locked-source operator equality including cloud-base overrides and limiting
levels; full L128 nonuniform geometry plus separately labeled L72
compatibility; FE count and buffer-parity replay; IMEX/BE fail-closed probes;
analytic and 1,875-case regression gates; two-tolerance reference checks;
deterministic weighted-L1/Linf aggregation and casewise-status fixtures;
hierarchical paired-bootstrap ratio fixtures; pilot/expansion selection and
complete census accounting; 64/65/512/513 boundary cases; both 1% censoring
thresholds; timing provenance and instrumented-data exclusion; and a whole-ACI
downstream-science control.

No missing case may be silently dropped. Schema/provenance/operator mismatch is
`UNAVAILABLE`; a direct candidate failure is `NUMERICAL_FAILURE`; a failed
reference self-check is `REFERENCE_UNRESOLVED`; and an exact count above 512 is
`RESOURCE_CENSORED`. Executed evidence, proposed tests, and untested regimes
must remain explicitly separated in every report.

### Authorization boundary and blocking approvals

Approval of this amendment freezes the scientific, numerical, sampling,
statistical, target-model, and instrumentation contracts in this section. It
authorizes the next offline design step only: specify the exact serialized
schemas and C++ interfaces in a follow-on addendum. Only explicit approval of
that addendum may authorize implementation and local testing of the validator,
locked-source replay adapter, census/selection code, replay benchmark, and
analysis in the confined worktree and approved scratch/result roots. This
amendment does not authorize changing production MAM4xx/EAMxx source,
inserting capture instrumentation into a production build, accessing prior
retained Task 5 results, connecting to an HPC system, issuing scheduler or
transfer commands, running production capture or GPU jobs, publishing results,
changing physics, or recommending a production method.

Approvals use ordinary immutable JSON or document records identified by
SHA-256. Each record names the permitted transition and the exact referenced
document or manifest hash. A status label alone grants no authority. This is a
fixed approval sequence, not a generalized workflow, signing, or external
approval system:

| State | Required transition evidence | Authority granted and exclusions |
|---|---|---|
| `AMENDMENT_PENDING` | Current unapproved amendment | Read-only review only. |
| `AMENDMENT_APPROVED` | Explicit approval referencing the amendment hash | Authorizes only preparation of the exact interface/schema addendum. |
| `ADDENDUM_APPROVED` | Explicit approval referencing the reviewed addendum hash | Authorizes only local offline implementation and testing of the validator, locked-source replay adapter, census/selection code, replay benchmark, analysis software, schemas, manifests, and associated offline tests. It does not authorize a production-source edit or production-state execution. |
| `CAPTURE_WRITER_APPROVED` | Separate approval referencing the addendum, exact production files and signatures, and changed-file allowlist | Authorizes only the named observational production capture-writer edit. It does not authorize a capture run. |
| `QUALIFICATION_CAPTURE_APPROVED` | Instrumentation-only provenance proof plus a bounded qualification run/resource manifest and command-specific HPC approvals | Authorizes only the three-arm qualification capture. This authority is non-transferable and does not authorize census capture. Qualification artifacts are inadmissible for census, timing, error, replay, or method-selection conclusions. |
| `OBSERVATIONAL_EQUIVALENCE_APPROVED` | Explicit approval referencing the successful qualification manifest and three-arm evidence, exact tested instrumentation implementation and executable/build identity, and the equivalence-relevant backend, compiler configuration, layout, seeds, and capture settings | Qualifies that exact capture-enabled implementation and configuration for a separately approved census. It does not itself authorize census capture. |
| `CENSUS_CAPTURE_APPROVED` | Observational-equivalence approval plus a census-specific frozen population/resource manifest, instrumentation-only proof, and command-specific approvals | Authorizes only the bounded census capture after verifying that its producer matches the qualified capture-enabled implementation, executable/build identity, backend, compiler configuration, layout, seeds, and equivalence-relevant capture settings. Separately valid provenance is not a substitute. A change to a bound identity or setting requires renewed three-arm qualification. The census population and immutable `census` role may intentionally differ from the qualification population and `qualification` role; that exception authorizes no other difference. |
| `CENSUS_SEALED` | A complete validated census manifest with every expected record accounted for and all completion records and hashes valid | Authorizes preparation of census-bound selection/statistics preregistration only. Candidate execution and result inspection remain prohibited. |
| `PREREGISTRATION_APPROVED` | Explicit approval referencing the sealed census hash, exact candidate/reference/timing/analysis identities, selection rules, statistical conventions, seeds, weights, and analysis version | Authorizes materializing the frozen selection and preparing a bounded replay execution manifest. Any bound change invalidates this approval. |
| `REPLAY_EXECUTION_APPROVED` | Preregistration plus an exact replay-specific run/resource manifest and command approvals | Authorizes only the bounded pilot or expansion candidate replay and inspection named by that manifest. Each expansion requires renewed bounded approval. It does not authorize GPU execution, publication, physics changes, or a production recommendation. |

Run/resource approval must identify the exact supported machine, locked stack
and clean checkout, cases and timesteps, input-data authorization,
output/retention paths, capacity and inode checks, and expected CPU/GPU,
memory, I/O, network, and storage effects. Command-specific approval requires
the applicable live facility-policy review, exact shell and scheduler scripts,
project/account, non-priority queue or QOS unless separately authorized,
nodes/tasks/CPUs/GPUs/memory, walltime, allocation charge, and human approval
for every command. Approval of this design is not job-submission approval.

Later production recommendation remains separately blocked by a physics-owner
accuracy budget, completed census with acceptable censoring, certified replay
evidence, and the approved whole-ACI production-layout GPU gate.

## Production-State Replay Serialized-Schema and C++-Interface Addendum

> **[ADDENDUM APPROVED - 2026-09-24]** The user explicitly approved the exact
> addendum content at document SHA-256
> `ca30c32c85451ff2fa14666119cfa7da8f665c9d5f2eb662c6a40574b513f470`
> for the transition `AMENDMENT_APPROVED` to `ADDENDUM_APPROVED`. This
> authorizes only local offline implementation and testing of the validator,
> locked-source replay adapter, census and selection code, replay benchmark,
> analysis software, schemas, manifests, and associated offline tests in the
> approved allowlist. It does not authorize any production MAM4xx or EAMxx
> source edit (including the capture writer), physics or numerical change,
> dependency or gitlink change, capture, production-state replay or benchmark
> execution, performance analysis, method-selection conclusion, retained Task
> 5 result access, remote or HPC activity, publication, staging, committing, or
> pushing. This status revision records the approval without substantively
> changing the approved schemas, interfaces, or contracts.
>
> The addendum was prepared under the approved amendment snapshot
> `d79a14253ca3acb5108d54f37b3ddcc9aaa2a9d9c4ecef391721c27868df2138`,
> whose bytes were recovered and verified from the local Git object database.
> The user also approved document snapshot
> `a829ffc01d09a9960e9bc5df04a15199a39833fd6942cb696a8c942bb1292c38`
> for the contract-only transition `AMENDMENT_APPROVED` to
> `AMENDMENT_APPROVED`. That approval resolved only OI-1 and OI-2 exactly as
> documented: the P1a/P1b/P1c and post-reset capture/evidence contract, and the
> conditional replay-eligible `R`/`C`/`U` estimand and denominators. Every
> unlisted contract and later approval remains unchanged.

This addendum makes the amendment's serialization and interface requirements
implementation-exact. It does not change the target model, equations,
candidate methods, limiters, clipping, boundaries, coupling, acceptance gates,
sampling, censoring, timing eligibility, or decision rules. A later
implementation must use these names, field sets, orderings, types, and failure
states exactly; an incompatible change requires an approved document
amendment before code changes.

The **Offline diagnostic-interface amendment (OI-5)** below is approved at
document SHA-256 `2300c9ef1448fb161042038e0075d8646cdfc733624ba2c0f7160e21b2f2fbe0`
for `ADDENDUM_APPROVED` to `ADDENDUM_APPROVED`. It authorizes only its four
diagnostic-only offline file additions, exact diagnostic exposure, two explicit
representation corrections, reference wrapper classifier, and associated local
offline implementation/testing. The original addendum and all later gates
remain in force. Approval is not evidence of completed implementation or tests.

### Common serialization contract

All retained artifacts in this addendum use UTF-8 JSON without a byte-order
mark and with one trailing line feed. Duplicate keys, unknown keys, non-UTF-8
text, JSON floating-point numbers, `NaN`, and infinity are invalid. Canonical
JSON sorts object keys lexicographically by their UTF-8 byte sequence, emits no
insignificant whitespace, escapes quotation mark, reverse solidus, and U+0000
through U+001F only, uses the two-character escapes `\b`, `\t`, `\n`, `\f`, and
`\r` for those five control characters and lowercase `\u00xx` for the others,
emits every other Unicode scalar value directly as UTF-8, preserves array
order, and emits integers in minimal base-10 form.
Booleans and `null` use the JSON spellings. Scientific real values present in
JSON are strings called `real_hex`. Unless a field explicitly names the binary
payload's scalar encoding, `real_hex` is canonical IEEE-754 binary64 encoded
from its 64-bit representation as follows: normal values use
`[-]0x1.` followed by exactly 13 lowercase hexadecimal fraction digits and a
minimal signed decimal exponent; subnormal values use `[-]0x0.` followed by
exactly 13 fraction digits and exponent `p-1022`; positive and negative zero
are respectively `0x0.0000000000000p+0` and
`-0x0.0000000000000p+0`. Thus one is `0x1.0000000000000p+0`, one half is
`0x1.0000000000000p-1`, and the least positive subnormal is
`0x0.0000000000001p-1022`. Spellings such as `0x1p+0`, `0x1.0p+0`, or
`0x2p-1` are invalid even when numerically equivalent. Nonfinite values are
invalid except where a schema represents failure through an enum and `null`.

Canonical parsing must recover the exact binary64 bits, including the sign of
zero and every subnormal. Hashing and uniqueness compare canonical bytes.
Numeric sorting uses IEEE-754 `totalOrder`; therefore `-0` precedes `+0`, and
bit-distinct values are never silently normalized. A source-stratum timestep
must additionally be finite and positive, so neither signed zero is valid
there. Binary payloads retain their declared binary32 or binary64 bit patterns
directly and are not converted through JSON. Future tests must include
round-trip, signed-zero, least-subnormal, alternative-spelling rejection,
unique-timestep, total-order, and byte-identical hash fixtures.

Every field named `*_sha256` is exactly 64 lowercase hexadecimal characters.
A hash adjacent to a `*_path` field hashes the byte-for-byte final artifact at
that path. A standalone content hash such as `source_identity_sha256`,
`population_sha256`, `run_resource_sha256`, `identity_sha256`,
`census_sha256`, `sentinel_content_sha256`, `preregistration_sha256`,
`replay_run_resource_sha256`, `execution_manifest_sha256`, or
`final_manifest_sha256` hashes the canonical content specified by its schema,
with the self-hash field omitted
where stated; it is not silently substituted for a file hash.

Artifact `*_path` fields are normalized POSIX paths relative to the containing
immutable run root; absolute paths, `..`, symbolic-link escape, and path aliases
are rejected. The `working_directory`, members of `input_paths`, and
`output_root` in `ndrop_run_resource_v1` and
`ndrop_replay_run_resource_v1` are the only stored host-location
exceptions: they are normalized absolute paths, contain no `..`, are bound to
the named machine/facility, and are never resolved as retained-artifact paths.
Every `run_id` matches
`[A-Za-z0-9._-]{1,128}`. Timestamps use UTC `YYYY-MM-DDTHH:MM:SSZ`. Schema
objects are closed: a missing required field or an extra field is invalid.

The following schema identifiers are reserved exactly:

| Artifact | `artifact_schema` |
|---|---|
| Approval record | `ndrop_approval_v1` |
| Capture population | `ndrop_capture_population_v1` |
| Run/resource manifest | `ndrop_run_resource_v1` |
| Capture run identity | `ndrop_capture_identity_v1` |
| Complete capture metadata | `ndrop_capture_record_v1` |
| Per-record completion | `ndrop_capture_completion_v1` |
| Qualification timing observation | `ndrop_capture_qualification_timing_v1` |
| Qualification field inventory | `ndrop_capture_qualification_field_inventory_v1` |
| Qualification production-sequence observation | `ndrop_capture_qualification_sequence_v1` |
| Qualification per-call trace | `ndrop_capture_qualification_trace_v1` |
| Qualification arm output | `ndrop_capture_qualification_arm_v1` |
| Qualification report | `ndrop_capture_qualification_v1` |
| Census diagnostic cross-tabulation | `ndrop_capture_census_diagnostics_v1` |
| Sealed census manifest | `ndrop_capture_census_v1` |
| Replay preregistration | `ndrop_replay_preregistration_v1` |
| Stress-sentinel definition | `ndrop_replay_sentinel_v1` |
| Expansion decision evidence | `ndrop_replay_expansion_decision_v1` |
| Replay run/resource manifest | `ndrop_replay_run_resource_v1` |
| Per-method replay result | `ndrop_replay_result_v1` |
| Candidate state payload | `ndrop_replay_state_v1` |
| Numerical replay diagnostics | `ndrop_replay_diagnostics_v1` |
| Primary timing observation | `ndrop_replay_timing_v1` |
| Tier-2 instrumentation overhead | `ndrop_replay_overhead_timing_v1` |
| Diagnostic timing observation | `ndrop_replay_diagnostic_timing_v1` |
| Replay execution completion | `ndrop_replay_execution_manifest_v1` |
| Final replay/analysis manifest | `ndrop_replay_final_manifest_v1` |

Unless a field definition below is stricter, schema primitives are fixed as
follows. A “nonempty string” contains 1--256 UTF-8 bytes and no NUL; a name or
ID has the stated regex or this bound. An unsigned integer must fit its stated
width; an unqualified count is unsigned 64-bit. A path/hash pair is either both
non-null or both null, and parallel path/hash arrays have equal length. Sets
are serialized as duplicate-free arrays in canonical byte order; sequences
retain their declared semantic order. `case_ids`, `member_ids`, column IDs,
and command-approval references are unique; case/member arrays are canonical
byte sorted, column IDs numeric sorted, and commands remain execution ordered.
Every approval or data-authorization reference is a closed object containing
exactly `path` and `sha256`; descriptive authority text remains only in the
referenced approval record. Command fields are strings; command-list fields
are ordered arrays of nonempty strings and never shell fragments supplied by a
default. A nullable field is nullable only where this addendum says so.

Every `source_stratum` is either null where explicitly allowed or a closed
object containing exactly `dtmicro`, `fe_regime`, and `activation_regime`.
`fe_regime` is `weak`, `moderate`, `strong`, `near_stiff`, or `extreme_tail`
with the approved `nsubmix` boundaries; `activation_regime` is `inactive`,
`active_unlimited`, or `limiter_active`. Canonical stratum ordering is numeric
`dtmicro`, then that FE enum order, then that activation enum order. All schema
validators reject a nested unknown key just as they reject a root unknown key.
These rules do not supply a missing threshold, seed, population, resource, or
approval value.

An approval record contains exactly `artifact_schema`, `created_utc`,
`approver`, `transition_from`, `transition_to`, `referenced_artifact_path`,
`referenced_artifact_sha256`, and `authority_and_exclusions`. The last field is
the complete approved scope text, not a status label. Each state transition
must reference the exact required document or manifest hash from the amendment
state table. No approval record is inferred from a filename, directory, branch,
tag, commit message, or state string.

A contract-only amendment that does not advance the state uses
`transition_from=AMENDMENT_APPROVED` and
`transition_to=AMENDMENT_APPROVED`; its authority text must enumerate the
exact requirements replaced and state that every unlisted contract and every
later transition remains unchanged. This is an ordinary hash-bound approval
record, not a new state or approval framework.

For the approved OI-5 revision of the already-approved addendum, the
record instead uses `transition_from=ADDENDUM_APPROVED` and
`transition_to=ADDENDUM_APPROVED`, with the same exact-hash and enumerated-scope
requirements. The user's explicit approval binds document SHA-256
`2300c9ef1448fb161042038e0075d8646cdfc733624ba2c0f7160e21b2f2fbe0`;
this status bookkeeping does not expand that authority.

Files are written to exclusive `.partial` paths, flushed, closed, reopened,
validated, and SHA-256 hashed before atomic promotion in the same filesystem.
The per-record completion JSON is written, flushed, validated, hashed by its
run manifest, and promoted last. A payload or metadata file without its valid
completion record is nonexistent to every consumer. A failed write, flush,
close, reopen, validation, hash, or rename produces no available record and
must be represented by a valid `UNAVAILABLE` completion record when that final
record can itself be written; otherwise the expected record remains
unaccounted and blocks census sealing.

### Immutable source identity

The `source_identity` object is common to all per-record capture, census,
selection, replay-result, and replay-timing artifacts and contains exactly the
fields below. Qualification arm-level timing instead binds the qualification
identity because it is not a production record:

| Field | Type and constraint |
|---|---|
| `case_id` | nonempty string, maximum 256 UTF-8 bytes |
| `member_id` | nonempty string; use literal `none` when no ensemble member exists |
| `physics_step` | unsigned 64-bit integer |
| `mac_subcycle` | unsigned 32-bit integer |
| `aci_call` | unsigned 32-bit integer within that MAC subcycle |
| `global_column_id` | unsigned 64-bit integer |
| `capture_sequence` | unsigned 64-bit integer, unique and contiguous from zero within a run |
| `grid_id` | nonempty string, maximum 128 UTF-8 bytes |
| `pver` | unsigned 32-bit integer; production conclusions require `128`, while compatibility records use `72` |

`source_identity_sha256` is the SHA-256 of the canonical JSON encoding of this
object. It is the record key and filename stem. A duplicate key, gap in the
expected sequence, or disagreement between an embedded identity and its hash
is fail-closed. Level index zero is model top and increasing level index is
downward. Constituent order is the locked MAM4xx `mam_idx` order: accumulation
number then seven masses, Aitken number then four masses, coarse number then
seven masses, and primary-carbon number then three masses. Mode order is
accumulation, Aitken, coarse, primary carbon. The identity manifest stores the
canonical mapping object and its hash; a record may not supply a replacement
mapping.

### Population and run/resource schemas

`ndrop_capture_population_v1` contains exactly `artifact_schema`,
`created_utc`, `population_id`, `capture_role`, `grid_role`, `case_ids`, `member_ids`,
`time_windows`, `spinup_exclusions`, `dtmicro_values`, `column_ids`,
`mac_aci_calls`, `cadence`, `exclusion_rules`, `adaptive_thinning`,
`expected_source_identities`, `expected_record_count`,
`input_data_authorization`, and `population_sha256`. `capture_role` is exactly
`qualification` or `census`. `grid_role` is exactly `production_l128` or
`compatibility_l72`; every expected source identity must respectively have
`pver=128` or `pver=72`. `time_windows` is an ordered array of closed
objects containing exactly `case_id`, `member_id`, `start_model_time`,
`end_model_time`, and `calendar`; the start is inclusive and the end is
exclusive. `spinup_exclusions` uses closed objects containing exactly
`case_id`, `member_id`, and `end_model_time`. `mac_aci_calls` is an ordered
array of closed objects containing exactly `mac_subcycle` and `aci_call`.
`cadence` is a closed object containing exactly `first_eligible_call`,
`call_stride`, and `column_stride`, all nonnegative integers with positive
strides. `dtmicro_values` are unique positive canonical `real_hex` strings,
sorted by the numeric total-order rule above, and `column_ids` are unique
unsigned 64-bit integers. `exclusion_rules` is the complete ordered array of
closed objects containing exactly `rule_id`, `scope`, `predicate`, `field`,
`comparison`, `value_kind`, `real_value`, `integer_value`, `string_values`,
and `reason`. `rule_id` is unique and matches `[A-Za-z0-9._-]{1,128}`.
`scope` is `source_identity`, `model_time`, or `capture_input`; `predicate` is
`comparison` or `set_membership`; `comparison` is `eq`, `ne`, `lt`, `le`,
`gt`, `ge`, `in`, or `not_in` and must match the predicate; and `value_kind`
is `real_hex`, `integer`, or `string_set`. Exactly the corresponding value
field is non-null. `string_values` is sorted, unique, and nonempty when used.
Rules are evaluated in array order against the pre-capture value of the named,
schema-defined field, and the first matching rule supplies the retained
exclusion reason. A field without a declared type and unit is invalid.
Thresholds and set members are fixed by the later population approval, not by
implementation defaults. A rule added or changed at runtime is invalid.
`adaptive_thinning` must be false.

`expected_source_identities` is the authoritative ordered list of the complete
source-identity objects, with `capture_sequence` contiguous from zero. All
descriptive case, member, window, timestep, column, call, cadence, spin-up, and
exclusion fields must reproduce exactly that list. `expected_record_count`
equals its length. `population_sha256` hashes canonical JSON with that field
omitted. Qualification and census population objects may differ, but their
roles and identities cannot be relabeled.

`ndrop_run_resource_v1` contains exactly `artifact_schema`, `created_utc`,
`run_id`, `capture_role`, `population_manifest_path`,
`population_manifest_sha256`, `machine`, `facility`, `project_account`,
`queue_or_qos`, `priority_authorized`, `nodes`, `tasks`, `cpus_per_task`,
`gpus_per_node`, `memory_bytes`, `walltime_seconds`,
`expected_allocation_cost`, `working_directory`, `input_paths`,
`input_sha256s`, `output_root`,
`expected_output_bytes`, `expected_output_inodes`, `available_bytes`,
`available_inodes`, `network_effects`, `exact_environment_commands`,
`exact_build_command`, `exact_run_commands`, `retention_policy`,
`data_authorization`, `facility_policy_review_path`,
`facility_policy_review_sha256`, `command_approval_paths`,
`command_approval_sha256s`, and `run_resource_sha256`. Local-only entries use
literal `local` for `facility` and `null` only for scheduler fields that do not
exist. Every command is a complete string in execution order.
`priority_authorized` must be false unless a separate command-specific approval
explicitly names the priority class. Capacity values are recorded immediately
before the authorized run. `input_paths` and `input_sha256s` have equal length
and bind every external run input file byte-for-byte; a directory input is
represented by a retained canonical inventory manifest and its file hash.
`run_resource_sha256` hashes
canonical JSON with that field omitted.

### Capture-run identity schema

`ndrop_capture_identity_v1` is sealed before a qualification or census run and
contains exactly these fields:

| Field | Required value or object |
|---|---|
| `artifact_schema`, `created_utc`, `run_id` | Common values above |
| `capture_role` | exactly `qualification` or `census`; immutable for the run |
| `grid_role` | exact population value `production_l128` or `compatibility_l72` |
| `approval_records` | ordered closed objects containing exactly `path` and `sha256` for every approval authorizing this run |
| `population_manifest_path`, `population_manifest_sha256` | exact qualification or census population definition |
| `run_resource_manifest_path`, `run_resource_manifest_sha256` | exact approved bounded resources and commands |
| `locked_mam4xx_commit` | `0d076d8ebb004a18f05f51f97b08edd9def6826f` |
| `instrumentation_source_kind` | exactly `descendant_commit` or `canonical_patch` |
| `instrumented_mam4xx_identity` | exact 40-character lowercase descendant commit when `instrumentation_source_kind=descendant_commit`, or exact canonical patch SHA-256 otherwise |
| `instrumented_mam4xx_tree` | exact clean tree produced by the descendant or canonical patch application |
| `base_e3sm_commit`, `base_e3sm_tree` | exact clean E3SM endpoint before the instrumentation-only edit |
| `instrumented_e3sm_commit`, `instrumented_e3sm_tree` | actual clean post-implementation E3SM commit and tree |
| `base_to_instrumented_diff_path`, `base_to_instrumented_diff_sha256` | complete canonical diff |
| `changed_file_allowlist` | exact ordered path array from this addendum |
| `dependency_gitlinks` | ordered closed objects containing exactly `path`, `mode`, and `object_id`; each must equal the locked E3SM tree |
| `dependency_provenance_path`, `dependency_provenance_sha256` | current locked E3SM/MAM4xx/EKAT/Kokkos provenance; historical Haero provenance is prohibited |
| `executable_path`, `executable_sha256`, `build_manifest_path`, `build_manifest_sha256` | exact retained executable and build identity |
| `compiler`, `compiler_flags`, `backend`, `execution_space`, `layout` | complete qualification-relevant build/runtime strings |
| `capture_settings` | closed object containing exactly `enabled`, `record_cadence`, `record_capacity`, `overflow_policy`, and `output_root`; overflow policy is `unavailable` |
| `seed_manifest_path`, `seed_manifest_sha256` | exact model and qualification seeds; no inferred default |
| `model_configuration_path`, `model_configuration_sha256` | exact retained runtime configuration consumed by the capture producer |
| `grid_definition_path`, `grid_definition_sha256` | exact retained grid and vertical-coordinate definition used by the captured population |
| `scalar_encoding` | exactly `ieee754_binary32_le` or `ieee754_binary64_le` |
| `source_mapping`, `source_mapping_sha256` | closed object containing exactly ordered `constituents` and `modes` arrays, plus its canonical hash |
| `record_schema`, `completion_schema` | exactly `ndrop_capture_record_v1` and `ndrop_capture_completion_v1` |
| `expected_record_count` | positive unsigned 64-bit integer fixed by the population manifest |
| `record_directory` | normalized relative directory |
| `identity_sha256` | SHA-256 of canonical JSON with this field omitted |

An identity supplied through `EAMXX_NDROP_CAPTURE_IDENTITY` must have
`capture_settings.enabled=true`; false is valid only as a rejected/negative
test fixture and cannot construct a `CaptureWriter`.

The `source_mapping.modes` array has exactly four closed entries containing
`index`, `name`, and `species_count`: `(0,accumulation,7)`,
`(1,aitken,4)`, `(2,coarse,7)`, and `(3,primary_carbon,3)`. The
`source_mapping.constituents` array has exactly 25 closed entries containing
`index`, `mode_index`, `within_mode_index`, `name`, `quantity_kind`, and
`unit`. In increasing `index` it is exactly:
`accumulation_number`, `accumulation_so4`, `accumulation_pom`,
`accumulation_soa`, `accumulation_bc`, `accumulation_dst`,
`accumulation_nacl`, `accumulation_mom`, `aitken_number`, `aitken_so4`,
`aitken_soa`, `aitken_nacl`, `aitken_mom`, `coarse_number`, `coarse_dst`,
`coarse_nacl`, `coarse_so4`, `coarse_bc`, `coarse_pom`, `coarse_soa`,
`coarse_mom`, `primary_carbon_number`, `primary_carbon_pom`,
`primary_carbon_bc`, and `primary_carbon_mom`. Number entries have
`within_mode_index=0`, `quantity_kind=number`, and `unit=number_per_kg`;
successive mass entries use locked species order with
`within_mode_index=1..species_count`, `quantity_kind=mass`, and
`unit=kg_per_kg`. This is the locked `mam_idx` order derived from
`src/mam4xx/aero_modes.hpp:127-154`; no run may rename, reorder, extend, or
sparsify it.

The current dependency-provenance artifact contains exactly an ordered array
of closed entries with `name`, `path`, `gitlink_object_id`, `source_commit`,
`source_tree`, and `applicability`, followed by its own `provenance_sha256`.
The current names are `e3sm`, `mam4xx`, `ekat`, and `kokkos`; applicability is
`target`, `adapter`, `utility`, or `execution_backend`. A null gitlink is
allowed only for the top-level E3SM entry or a dependency nested through an
already bound dependency. The locked current adapter has no Haero dependency,
as recorded by the implemented Task 2 adapter and confirmed by the locked E3SM
gitlinks. A `haero` entry is therefore a role/provenance mismatch, not a
current dependency. Historical Haero hashes may remain in explicitly labeled
historical evidence only and cannot qualify a capture or replay stack.

The manifest validator recomputes the tree, diff, executable, build,
dependency, gitlink, mapping, approval, and population hashes before accepting
a run. An extra changed path, dirty repository, unlisted path, changed gitlink,
different executable, or mismatch between an approval and its referenced hash
invalidates the identity. The role cannot be changed by copying, renaming, or
embedding the artifact in another manifest.

The canonical complete-diff artifact is UTF-8 with LF line endings and
contains, in order, exactly one `MAM4xx` section and one `E3SM` section. Each
section has this exact byte framing, with angle-bracket terms replaced by
their ASCII values and no blank lines:

```text
NDROP-COMPLETE-DIFF-V1\n
repository=<MAM4xx-or-E3SM>\n
base_commit=<40-lowercase-hex>\n
base_tree=<40-lowercase-hex>\n
instrumented_commit=<40-lowercase-hex-or-64-lowercase-patch-sha256>\n
instrumented_tree=<40-lowercase-hex>\n
diff_bytes=<minimal-unsigned-decimal>\n
<exactly-diff_bytes-bytes>
```

The first line appears once at offset zero. The remaining six metadata lines
and diff body appear first for `MAM4xx`, then immediately for `E3SM`; each diff
body must already end in one LF and no separator byte is added. Diff bytes are
the exact output of `git -c core.quotePath=true diff --binary --full-index
--no-ext-diff --no-color <base> <instrumented> --` under `LC_ALL=C`, with
external diff drivers disabled. Empty output is represented by
`diff_bytes=0\n` followed immediately by the next section or end of file. The
MAM4xx base is the locked commit above; the E3SM base is the clean parent
selected for the instrumentation descendant and is recorded rather than
inferred. The validator regenerates both sections, verifies the metadata and
changed paths against the addendum allowlist, and hashes the complete framed
artifact. A merge commit, non-parent base, or canonical patch application must
state both exact endpoints in the source identity and reproduce the same bytes;
otherwise provenance fails.

### Locked-source verification and approved OI-1 capture contract

Read-only inspection of the exact locked sources confirmed the following. In
MAM4xx commit `0d076d8ebb004a18f05f51f97b08edd9def6826f`, `ColumnView` is an
unmanaged rank-1 device view and `View2D`/`View3D` are device views. The locked
`update_from_explmix` accepts rank-3 `raercol` and `raercol_cw` views at
`src/mam4xx/ndrop.hpp:1053-1080`; its team reduction first produces overlaps
and directional rates at lines 1101-1122, applies the in-place activation-rate
limiters at lines 1131-1133, completes that reduction at line 1151, and only
then produces `nsubmix` at lines 1153-1159. The final zero-cloud
resuspension/reset is performed at lines 1254-1275. The caller creates the raw
activation rates at lines 1468-1497 and calls `update_from_explmix` at lines
1503-1506. The locked source stores neither a limiting-level index nor a
separate `max_tinv` value.

The exact EAMxx caller at
`components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:230-267` uses
rank-4 work arrays across columns; each column is passed to the locked MAM4xx
call as rank-3 views at lines 328-329 and 395-422. These facts confirm an
addendum defect, not a demonstrated compilation or runtime failure: the
array-of-`View1D` declarations previously proposed below did not match the
locked API and are corrected in this revision.

They also exposed a material conflict in the original approved amendment.
Before the OI-1 contract-only approval, its single Phase 1 was fixed immediately
before the call to `update_from_explmix`, yet it required that overlaps,
directional rates, post-limiter rates, limiter activity, selector metadata, and
limiting-level evidence already exist. They do not exist at that checkpoint in
the locked source. Relabeling offline reconstruction as observed production
data or moving production calculations remains prohibited.

The **approved OI-1 contract resolution** was accepted through the contract-only
approval of document SHA-256
`a829ffc01d09a9960e9bc5df04a15199a39833fd6942cb696a8c942bb1292c38`.
It amends only (a) the single-checkpoint meaning of capture Phase 1 and (b) the
requirement for a directly observed scalar limiting level and separate
`max_tinv`. Transaction Phase 1 comprises these ordered, enabled-path-only
subcheckpoints:

1. **P1a, pre-update snapshot.** After the existing barrier at locked
   `src/mam4xx/ndrop.hpp:1497` and before the call at lines 1503--1506, copy the
   common old `qcld`/interstitial/cloud-borne state, current `cldn`, required
   raw inputs, initial buffer parity, and stable buffer identities into
   capture-owned storage. This remains the immutable replay initial condition.
2. **P1b, limiter observation.** At each existing assignment at lines
   1132--1133, and only in the capture-enabled overload, retain the exact
   pre-assignment operand, `eddy_diff_kp` cap, stored post-assignment value, and
   direct activity flag for qualification; retain only the post-limiter value
   and direct flag in the census payload. For finite nonnegative operands the
   flag is true exactly when `pre_rate > cap`; equality is not limiter activity.
   The original `mam4::min` assignment remains the sole producer of the stored
   rate. No comparison result feeds production.
3. **P1c, selector observation.** After the existing reduction barrier at line
   1151 and the existing `nsubmix` calculation at lines 1153--1159, but before
   the first substep begins at line 1169, retain completed overlaps,
   directional rates, post-limiter rates, the actual reduction output
   `dtmin`, and the actual `nsubmix`. Use the existing barriers at lines 1151
   and 1168; the amendment authorizes no new reduction, altered reduction
   order, or blanket synchronization.

The locked source does not produce a unique limiting level, an empty-reduction
flag, a separate selector maximum, or the controller-only FE contribution.
Those values therefore become explicitly **offline reconstructions** rather
than observed production fields. The adapter forms, in increasing production
level order, the finite candidate set
`C={(k,tinv_k): top_lev<=k<pver and tinv_k>1e-6}` using the locked scalar width
and operation order. It retains `selector_candidate_count=|C|` and every level
whose `1/tinv_k` has the same scalar bit pattern as observed `dtmin` in the
ordered array `selector_minimizer_levels`. No unique limiting level is
fabricated: an empty array is valid only when `C` is empty, multiple levels
remain a tie set, and the smallest member may be displayed only as a labeled
reporting projection. A nonempty `C` with no bitwise match, or a reconstructed
`nsubmix` that differs from the observed integer, fails operator certification.
The controller-only FE rate and selector maximum are reconstructed and labeled
as such; neither is serialized as a production observation.

After `update_from_explmix` returns, the existing caller barrier at line 1508
is the first source-anchored point at which the post-substep zero-cloud
resuspension/reset at lines 1254--1275 is complete for the team. The enabled
overload copies final state and final parity to capture-owned device storage at
that point, before the work views or local `nnew` can be reused. The wrapper
then retains the actual returned `nsubmix`. Transaction Phase 2 remains at the
caller after the complete production wrapper call and the existing EAMxx
completion fence; it copies the completed device slots to capture-owned host
storage, validates the phase identity, hashes, and seals. A backend may replace
that mapping only with an already existing equivalent completion guarantee
bound by qualification. No production alias may outlive its checkpoint.

The minimum operations now specified by the approved contract are capture-enabled
stores into preallocated slots at P1a, P1b, P1c, and the post-reset point;
qualification-only pre/cap/post trace stores; the already required enabled
device-to-host copy; and validation. The original non-capture function remains
the disabled path, so capture-disabled execution performs no capture load,
store, allocation, copy, barrier, fence, atomic, or file operation. This
contract changes no calculation, branch input, arithmetic, reduction, state
transition, limiter, or return value. The OI-1 contract-only approval alone
granted no implementation authority. `ADDENDUM_APPROVED` now authorizes the
specified offline interfaces and tests, but no production capture operation
may be added or run until `CAPTURE_WRITER_APPROVED`; records lacking the
required P1b/P1c or post-reset evidence are `UNAVAILABLE`.

The source-anchored producer/lifetime map is:

| Quantity | Evidence class and producer | First valid | Mutation/reuse boundary |
|---|---|---|---|
| `dtmicro`, `top_lev`, source identity, `cldn`, `pdel`, `pressure`, `height` | Raw caller inputs | On entry to locked `dropmixnuc`, lines 1279-1297 | Caller-owned; capture requires the pre-call owned copy before launch storage may be reused |
| `qcld`, active `raercol`, active `raercol_cw` common old state | Observed production state after cloud growth/decay/disappearance | After barriers at lines 1467 and 1497 | Mutated by `update_from_explmix` beginning with the limiter/rate work and first substep at lines 1131-1178; must be copied before the call at lines 1503-1506 |
| `eddy_diff`, `zn`, `zs`, `csbot`, `dz` | Observed setup/raw work values | Produced at lines 1382-1417 and cloud-profile work through 1497 | `eddy_diff` may receive the locked cloud-profile override through line 1493; stable for explicit integration after line 1497 |
| Pre-limiter `nact`, `mact` | Observed raw activation outputs; retained only in the qualification trace under the approved OI-1 contract | After `update_from_cldn_profile` and barrier, lines 1481-1497 | Overwritten in place by the limiters at lines 1131-1133; P1b must capture the exact operands before each assignment |
| `overlapp`, `overlapm`, `eddy_diff_kp`, `eddy_diff_km` | Observed production work values | Produced inside the reduction at lines 1101-1122 and complete after line 1151 | Stable through the substep loop, but nonexistent at the former single pre-call checkpoint |
| Post-limiter `nact`, `mact` and direct limiter flags | Rates are observed after lines 1131-1133; direct flags require P1b observation of the exact pre/cap operands at those assignments | After the reduction barrier at line 1151 | P1c retains post rates; P1b retains flags and qualification-only pre/cap/post evidence before storage reuse |
| Observed `dtmin` and `nsubmix` | Direct selector outputs of the locked reduction and arithmetic | `dtmin` completes at line 1151; `nsubmix` follows at lines 1153-1159 | P1c must preserve both before the substep loop or local reuse |
| Controller contribution, selector maximum, candidate count, and minimizing-level tie set | Offline reconstructions from retained observed inputs; the locked source stores none of them separately | Available only after offline reconstruction; the set is compared to observed `dtmin` and `nsubmix` | Must never be relabeled as observed production evidence; no scalar level is substituted for a tie set |
| Returned `nsubmix` and final parity | `nsubmix` is observed `int`; parity is the local `nnew` after swaps at lines 1169-1175 | After the complete substep loop | Final parity is preserved after the barrier at line 1508; the wrapper retains returned `nsubmix`; caller Phase 2 attaches both only after backend completion |
| Final `qcld`, `raercol`, `raercol_cw` | Observed final production state | Only after zero-cloud resuspension/reset at lines 1254-1275 and the enclosing barrier at line 1508 | Must be copied to capture-owned storage at that point before later output work or scratch reuse |
| EAMxx wrapper column/slot identity | Observed wrapper mapping | Column `icol=team.league_rank()` at wrapper lines 324-329 | Valid for that launch only; immutable host mapping must outlive completion and serialization |

“Offline reconstruction” is reserved for values recomputed by the locked
adapter from retained raw/observed inputs. It never changes the evidence class
to “observed production.” The interface, payload, and backend tables below
encode the approved OI-1 resolution. Their local offline implementation and
testing is authorized by `ADDENDUM_APPROVED`; their production-source
integration remains subject to `CAPTURE_WRITER_APPROVED`.

### Two-phase capture record and binary payload

One successful record consists of
`<source_identity_sha256>.capture.json`,
`<source_identity_sha256>.capture.bin`, and
`<source_identity_sha256>.completion.json`. The metadata JSON contains exactly:

| Field | Type and constraint |
|---|---|
| `artifact_schema` | `ndrop_capture_record_v1` |
| `run_id`, `capture_role` | exact identity-manifest values |
| `identity_manifest_path`, `identity_manifest_sha256` | exact sealed run identity |
| `source_identity`, `source_identity_sha256` | common immutable identity |
| `scalar_encoding` | exact identity-manifest value |
| `pver`, `top_lev` | unsigned integers with `0 <= top_lev < pver` |
| `initial_read_buffer`, `initial_write_buffer` | distinct integers in `{0,1}` |
| `final_read_buffer`, `final_write_buffer` | distinct integers in `{0,1}` |
| `droplet_buffer_id`, `interstitial_buffer_pair_id`, `cloudborne_buffer_pair_id` | nonempty run-local opaque identities, maximum 128 bytes each; memory addresses are prohibited |
| `production_nsubmix` | positive signed 64-bit integer returned by the completed production call |
| `phase1_completion` | closed object containing exactly `pre_update_snapshot`, `limiter_observation`, and `selector_observation`; each is a closed object containing exactly `source_anchor`, `producer_scope=team`, `mechanism=existing_team_completion`, and `complete=true` |
| `phase2_completion` | closed object with `producer_scope=execution_space`, exact backend mechanism string, and `complete=true` |
| `ownership_complete` | Boolean; must be true |
| `payload_path`, `payload_sha256`, `payload_bytes` | exact binary payload identity |
| `payload_real_count`, `payload_flag_count` | exactly `123*pver+2` and `8*pver` |

The three `source_anchor` values are exactly
`ndrop.hpp:1497-before-1503`, `ndrop.hpp:1132-1133`, and
`ndrop.hpp:1151-and-1153-1159-before-1169`, respectively. All three completion
objects must be present and true; a partial Phase 1 is not serializable.

The payload starts with the eight ASCII bytes `NDRPCAP1`, followed by one byte
whose value is `4` or `8` for the scalar width, one byte with value `1` for
little endian, and six zero reserved bytes. It then stores all real values in
the manifest's IEEE scalar width and all flags as unsigned bytes. Multidimensional
arrays are field-major: constituent or mode first, then level; level is always
the fastest-varying index. The real section order is exact:

| Order | Field | Units | Extent and validation |
|---:|---|---:|---|
| 1 | `dtmicro` | s | scalar, finite and positive |
| 2 | `observed_dtmin` | s | scalar, exact locked-reduction output; finite and positive for an available record |
| 3 | `qcld_initial` | # kg-1 | `pver`, finite and nonnegative |
| 4 | `interstitial_initial` | mixed # kg-1 or kg kg-1 | `25*pver`, finite and nonnegative |
| 5 | `cloudborne_initial` | mixed # kg-1 or kg kg-1 | `25*pver`, finite and nonnegative |
| 6 | `cldn` | 1 | `pver`, finite and in `[0,1]` |
| 7 | `eddy_diff` | m2 s-1 | `pver`, finite |
| 8 | `zn` | m2 kg-1 | `pver`, finite |
| 9 | `zs` | m-1 | `pver`, finite |
| 10 | `csbot` | kg m-3 | `pver`, finite |
| 11 | `overlapp` | 1 | `pver`, finite and in `[0,1]` |
| 12 | `overlapm` | 1 | `pver`, finite and in `[0,1]` |
| 13 | `nact` | s-1 | `4*pver`, post-limiter, finite and nonnegative |
| 14 | `mact` | s-1 | `4*pver`, post-limiter, finite and nonnegative |
| 15 | `pdel` | Pa | `pver`, finite and positive |
| 16 | `dz` | m | `pver`, finite and positive |
| 17 | `pressure` | Pa | production `pmid`, `pver`, finite, positive, increasing downward |
| 18 | `height` | m | production `zm`, `pver`, finite, decreasing downward |
| 19 | `production_eddy_diff_kp` | s-1 | `pver`, finite and nonnegative |
| 20 | `production_eddy_diff_km` | s-1 | `pver`, finite and nonnegative |
| 21 | `qcld_final` | # kg-1 | `pver`, finite and nonnegative |
| 22 | `interstitial_final` | mixed # kg-1 or kg kg-1 | `25*pver`, finite and nonnegative |
| 23 | `cloudborne_final` | mixed # kg-1 or kg kg-1 | `25*pver`, finite and nonnegative |

The flag section stores `number_limiter_active[mode][level]` followed by
`mass_limiter_active[mode][level]`, each exactly `4*pver` bytes and each byte
exactly zero or one. These flags record whether the locked production limiter
changed the corresponding pre-limiter rate; they do not change or recompute a
rate. They are **directly observed production flags**, not values independently
reconstructible from the retained post-limiter rate alone. Certification
requires the approved OI-1 checkpoint/evidence contract and a later approved
source edit to retain the exact pre-assignment operand, cap, stored result, and strict
`pre_rate > cap` flag at locked-source lines 1132--1133, with no
different-order recomputation. The qualification-only `branch_trace` must
retain those four values for every mode and level. Focused future
qualification fixtures must exercise below, exactly-at, and above-cap inputs
for number and mass in every mode, require equality to produce a false flag,
and require the trace relation and direct flag to agree exactly. Those
qualification trace values are not census payload fields and have zero census
weight. Independent reconstruction from a census record's post-limiter values
alone is explicitly prohibited; source review, the qualified trace evidence,
the approved unchanged instrumentation identity, and direct observation are
the certification chain. If any link is absent, the flag is uncertified and
the record is `UNAVAILABLE`. This flag evidence remains required for
activation-regime classification and sampling.

`observed_dtmin` and `production_nsubmix` are the only direct selector outputs.
The offline adapter retains reconstructed controller rate, selector maximum,
candidate count, and the complete minimizing-level tie set separately so the
physical operator cannot silently include the controller-only FE contribution
or relabel a reconstruction as observed. The payload byte count is therefore
`16 + scalar_width*(123*pver+2) + 8*pver` for the approved OI-1 field set.

P1a fills the common old state, raw inputs, initial parity, and buffer
identities before the first `update_from_explmix` mutation. P1b fills the flags
at the limiter assignments, and P1c fills fields 2 and 11--14 and 19--20 after
their producers complete. These requirements are contract- and
addendum-approved and may be represented and tested offline, but the actual
production fill operations remain unimplemented and require
`CAPTURE_WRITER_APPROVED`. Phase 2 attaches fields 21--23, actual
returned `production_nsubmix`, final parity, and completion metadata only at the caller
after the complete `dropmixnuc` call, including the zero-cloud
resuspension/reset at locked-source lines 1254-1275, and after every producer
has completed. The same source identity, role, scalar encoding, `pver`, and
buffer identities must match across phases. Phase 2 copies every value to
capture-owned storage before production storage can be reused. No alias or
view into production storage is serializable.

The full payload and its metadata are checksum-covered. Inactive levels remain
present; a validator must not truncate, fill, or normalize them. The replay
adapter operates on production levels `[top_lev,pver)` and retains the full
payload for provenance. A byte-count, order, unit, orientation, range, parity,
phase, or checksum mismatch is `UNAVAILABLE`.

### Per-record completion schema

`ndrop_capture_completion_v1` contains exactly:

| Field | Type and constraint |
|---|---|
| `artifact_schema`, `created_utc`, `run_id`, `capture_role` | common exact values |
| `identity_manifest_path`, `identity_manifest_sha256` | sealed run identity |
| `source_identity`, `source_identity_sha256` | expected record key |
| `status` | exactly `AVAILABLE` or `UNAVAILABLE` |
| `unavailable_reason` | `null` for `AVAILABLE`; otherwise exactly one of `capacity_exhausted`, `phase1_copy_failed`, `limiter_observation_failed`, `selector_observation_failed`, `post_reset_preservation_failed`, `phase1_completion_failed`, `phase2_producer_incomplete`, `phase2_copy_failed`, `phase_mismatch`, `ownership_failed`, `validation_failed`, `serialization_failed`, `hash_failed`, or `promotion_failed` |
| `capture_metadata_path`, `capture_metadata_sha256`, `capture_payload_path`, `capture_payload_sha256`, `payload_bytes` | required non-null exact values for `AVAILABLE`; all `null` for `UNAVAILABLE` |
| `phase1_complete`, `phase2_complete`, `ownership_complete` | all true for `AVAILABLE`; exact observed booleans for `UNAVAILABLE` |
| `source_stratum` | closed object containing `dtmicro`, FE regime, and activation regime when determinable, otherwise `null` |
| `source_frequency_weight`, `estimated_fe_work_weight` | nonnegative integers when determinable, otherwise `null` |

An `AVAILABLE` completion is admitted only after reopening and reproducing the
metadata hash, payload hash, byte count, all validations, and both phase
completion guarantees. An `UNAVAILABLE` completion accounts for an expected
source identity but carries no replay payload. If its stratum or either weight
is `null`, the census may be sealed for accounting only, but all
method-selection aggregation remains blocked pending an approved amendment.
For a `census` record, `source_frequency_weight` is exactly one and an
available record's `estimated_fe_work_weight` is exactly its captured
`production_nsubmix`; an unavailable record carries those same values only
when they are independently known from the completed capture phases. For a
`qualification` record both weights are exactly zero.

### Qualification and census schemas

Each `ndrop_capture_qualification_timing_v1` JSON Lines object contains exactly
`artifact_schema`, `run_id`, `qualification_identity_path`,
`qualification_identity_sha256`, `timing_block`, `interleaved_position`, `arm`,
`warmup_block_count`, `measured_block_count`, `minimum_measured_blocks`,
`maximum_measured_blocks`, `stopping_metric`, `stopping_interval_lower`,
`stopping_interval_upper`, `stopping_relative_half_width`, `stopping_reason`,
`seconds_per_complete_call`, `state_checksum`, `executable_sha256`, `compiler`,
`machine`, `backend`, and `execution_space`. `arm` is exactly `disabled` or
`enabled`; every block contains one matched observation from each arm with
interleaved order. The integer, Boolean, positive-`real_hex`, checksum, and
identity constraints are the same as primary timing below. It is diagnostic
qualification evidence only and is rejected by every census, replay, primary
timing, error, and method-selection consumer.
It uses matched disabled/enabled complete-call blocks with exactly 10 warmup
blocks and the same 30--100-block deterministic stopping rule specified below;
it has no phase fields. A partial block or fewer than 30 complete blocks cannot
support the overhead approval. Thus a Boolean completion label cannot
substitute for the retained counts and stopping evidence.

Each `ndrop_capture_qualification_field_inventory_v1` canonical JSON object
contains exactly `artifact_schema`, `created_utc`, `e3sm_source_commit`,
`e3sm_source_tree`, `process_dag_path`, `process_dag_sha256`, `scope`,
`processes`, and `inventory_sha256`. `scope` is exactly
`mamaci_call_boundary` or `downstream_science`. `processes` is ordered by the
approved production process DAG and contains closed objects with exactly
`process_name`, `process_instance`, `subcycle_count`, and `fields`.
`process_name` and `process_instance` are nonempty strings;
`subcycle_count` is a positive integer and must reproduce the approved
subcycle-configuration artifact. `fields` is ordered by `order_index` and
contains closed objects with exactly `field_name`, `role`, `data_type`,
`units`, `grid_name`, `layout_names`, `extents`, `diagnostic_hash_name`, and
`order_index`. `role` is `input`, `output`, or `internal` for
`mamaci_call_boundary`, and `downstream` for `downstream_science`;
`data_type` is `real` or `integer`; `layout_names` is a nonempty ordered array
of nonempty dimension-name strings; `extents` is an equal-length ordered array
of positive integers; `order_index` is contiguous from zero within a process;
and all remaining strings are nonempty. The inventory must reproduce the
source-declared EAMxx field registration and the exact order printed by the
approved diagnostics path; an unlisted, missing, differently laid-out, or
differently ordered field fails qualification rather than inheriting an
implementation default. `inventory_sha256` hashes canonical JSON with that
field omitted. The call-boundary and downstream inventories are separate
objects, and neither may substitute for the other.

`ndrop_capture_qualification_sequence_v1` is a canonical JSON Lines extraction
of the existing EAMxx per-process diagnostics. Each row represents one complete
MAMAci subcycle and contains exactly `artifact_schema`, `run_id`, `arm`,
`sequence_id`, `call_ordinal`, `source_identity`, `source_identity_sha256`,
`atmosphere_step`, `process_subcycle`, `pre_model_time`, `post_model_time`,
`pre_field_hashes`, `post_field_hashes`, `raw_log_path`, `raw_log_sha256`,
`previous_row_sha256`, and `row_sha256`. `call_ordinal` is contiguous from zero
across the retained production sequence. `previous_row_sha256` is null only on
the first row and otherwise equals the preceding row's canonical self-hash.
Each field-hash array is ordered exactly as printed and contains closed objects
with `field_name`, `layout`, `grid_name`, `role`, and `eamxx_bfb_hash`; `role`
is `input`, `output`, or `internal`, and the EAMxx hash is exactly 16 lowercase
hexadecimal digits. Every entry, array order, role, layout, and grid must match
the bound `mamaci_call_boundary` field inventory. `row_sha256` hashes canonical
JSON with that field omitted. The raw log is retained and file-hashed;
extraction must reproduce every MAMAci `pre-sc-*` and `pst-sc-*` block with no
gap, duplicate, or unmatched boundary.

Each enabled-arm `ndrop_capture_qualification_trace_v1` JSON Lines object
contains exactly `artifact_schema`, `run_id`, `sequence_id`, `source_identity`,
`source_identity_sha256`, `call_ordinal`, `qualification_sequence_row_sha256`,
`input_state_sha256`, `capture_completion_path`, `capture_completion_sha256`,
`pre_update_state_sha256`, `post_zero_cloud_state_sha256`,
`limiter_trace_path`, `limiter_trace_sha256`, `branch_vector_sha256`,
`production_nsubmix`, `initial_read_buffer`, `initial_write_buffer`,
`final_read_buffer`, `final_write_buffer`, `standard_output_row_sha256`, and
`trace_row_sha256`. Rows are in the same contiguous call order as the sequence
artifact. The limiter trace has the pre/cap/post/direct-flag fields specified
above. `branch_vector_sha256` hashes the canonical ordered zero-cloud, limiter,
selector, and parity decisions. The trace is capture-enabled qualification
evidence only; base and disabled arms must not create it.

Read-only source inspection confirms that no new all-arm evidence interface is
needed. The existing `AtmosphereProcess::run` retains one process object across
subcycles and invokes `run_impl` repeatedly
(`components/eamxx/src/share/atm_process/atmosphere_process.cpp:142-213`). With
`internal_diagnostics_level=2`, it emits pre- and post-`run_impl` hashes for
every process subcycle at lines 176--187. The existing hash path separately
hashes every input, output, and internal field and fences each device reduction
(`atmosphere_process_hash.cpp:196-257`). The generic atmosphere-driver test
retains one `AtmosphereDriver` across repeated steps
(`components/eamxx/tests/eamxx_ad_test.cpp:61-76`), and the locked MAMAci
`run_impl` contains exactly one `call_function_dropmixnuc` invocation at
`eamxx_mam_aci_process_interface.cpp:583-595` before its downstream ACI work and
final fence at lines 597--642. The existing instant-output path can retain
end-of-step downstream fields every step, as demonstrated by
`components/eamxx/tests/single-process/mam/aci/output.yaml:1-67`.

Qualification must therefore run the exact approved production process DAG,
subcycle configuration, layout, inputs, seeds, and number of consecutive steps
in one driver instance without resetting state, scratch, parity, or
accumulators between calls. The per-subcycle diagnostic blocks establish every
MAMAci call boundary, including repeated-call storage reuse; the frequency-one
standard output and the corresponding diagnostics of later processes establish
the approved downstream-science field set. A final atmosphere-step output alone
is insufficient whenever more than one MAMAci subcycle occurred. The
qualification manifest freezes the exact ordered input/output/internal field
inventory and downstream field list; a missing field or diagnostic block
fails closed.

These existing diagnostics intentionally synchronize and therefore cannot be
used for overhead. The separate disabled/enabled overhead execution sets
`internal_diagnostics_level=0`, retains no sequence hashes in its timed region,
and uses only the diagnostic complete-call timing artifact already specified.
No capture-only tracing is added to the base or disabled arm.

Each `ndrop_capture_qualification_arm_v1` manifest contains exactly
`artifact_schema`, `created_utc`, `run_id`, `arm`,
`population_manifest_path`, `population_manifest_sha256`,
`run_resource_manifest_path`, `run_resource_manifest_sha256`,
`e3sm_source_commit`, `e3sm_source_tree`, `mam4xx_source_kind`,
`mam4xx_source_identity`, `mam4xx_source_tree`, `executable_path`,
`executable_sha256`,
`build_manifest_path`, `build_manifest_sha256`, `seed_manifest_path`,
`seed_manifest_sha256`, `process_dag_path`, `process_dag_sha256`,
`subcycle_configuration_path`, `subcycle_configuration_sha256`,
`field_inventory_path`, `field_inventory_sha256`,
`downstream_field_inventory_path`, `downstream_field_inventory_sha256`,
`raw_diagnostics_log_path`,
`raw_diagnostics_log_sha256`, `sequence_evidence_path`,
`sequence_evidence_sha256`, `branch_trace_path`, `branch_trace_sha256`,
`nsubmix_path`, `nsubmix_sha256`, `parity_path`, `parity_sha256`, `state_path`,
`state_sha256`, `tendency_path`, `tendency_sha256`,
`downstream_science_path`, `downstream_science_sha256`, `sequence_id`,
`record_count`,
`capture_identity_path`, `capture_identity_sha256`,
`capture_completion_paths`, `capture_completion_sha256s`, and
`arm_manifest_sha256`. `arm` is exactly `base`, `disabled`, or `enabled`.
For base and disabled, the three branch/`nsubmix`/parity path/hash pairs and
the capture-identity path/hash pair are null and the completion arrays are
empty; those arms may not create capture-only traces. They are all required
for enabled. Raw-diagnostics, sequence-evidence, state, tendency, and
downstream-science paths are required for all arms and use the existing
diagnostic and standard science-output paths. Every populated output artifact
has one deterministic row per qualification-population identity in contiguous
production-sequence order. Path/hash arrays have equal length, all identities
and counts reconcile, and the self-hash omits only `arm_manifest_sha256`.
State and tendency artifacts are canonical extractions from the exact paired
per-subcycle sequence rows; downstream-science artifacts are canonical
frequency-one step records with explicit mappings to every preceding MAMAci
call in that step. The enabled `branch_trace` must be a valid
`ndrop_capture_qualification_trace_v1` with the same sequence keys.

`ndrop_capture_qualification_v1` contains exactly `artifact_schema`,
`created_utc`, `qualification_identity_path`,
`qualification_identity_sha256`, `qualification_approval_path`,
`qualification_approval_sha256`, `population_manifest_path`,
`population_manifest_sha256`, `base_arm`, `disabled_arm`, `enabled_arm`,
`base_vs_disabled`, `disabled_vs_enabled`, `overhead`, `record_completions`,
and `qualification_passed`. Each arm object contains the exact source commit
and tree, MAM4xx commit or patch hash, executable and build hashes, compiler,
backend, execution space, layout, seeds, capture setting, and output-manifest
identity through exactly the keys `e3sm_source_commit`,
`e3sm_source_tree`, `mam4xx_source_kind`, `mam4xx_source_identity`,
`mam4xx_source_tree`, `executable_path`,
`executable_sha256`, `build_manifest_path`, `build_manifest_sha256`,
`compiler`, `compiler_flags`, `backend`, `execution_space`, `layout`,
`seed_manifest_path`, `seed_manifest_sha256`, `capture_setting`,
`sequence_id`, `process_dag_path`, `process_dag_sha256`,
`subcycle_configuration_path`, `subcycle_configuration_sha256`,
`field_inventory_path`, `field_inventory_sha256`,
`downstream_field_inventory_path`, `downstream_field_inventory_sha256`,
`output_manifest_path`, and `output_manifest_sha256`. Each output manifest must
be a valid `ndrop_capture_qualification_arm_v1`. The three arm objects
must match in inputs, compiler, compiler flags, backend, execution space,
layout, seeds, sequence ID, process DAG, subcycle configuration, and both field
inventories;
`capture_setting` is respectively `absent`, `disabled`, and `enabled`. Here
`absent` means the locked base source has no instrumentation;
`disabled` means the instrumented descendant is built but
`EAMXX_NDROP_CAPTURE_IDENTITY` is absent; and `enabled` means that same
descendant and executable receive the sealed identity path.
`mam4xx_source_kind` is `locked_base` for the base arm and exactly the capture
identity's `instrumentation_source_kind` for both descendant arms.

`base_vs_disabled` contains exactly `record_count`, `sequence_hash_equal`,
`state_hash_equal`, `tendency_hash_equal`, `downstream_science_hash_equal`,
`no_capture_artifact_created`, `no_capture_operation_executed`, and `passed`.
It passes only when every paired pre/post field hash at every MAMAci subcycle
and all standard science outputs are bytewise equal, sequence order and storage
reuse match, and the disabled path created no capture artifact or performed a
capture allocation, copy, store, barrier, fence, or file operation.

`disabled_vs_enabled` contains exactly `record_count`, `sequence_hash_equal`,
`state_hash_equal`, `tendency_hash_equal`, `downstream_science_hash_equal`,
`branch_reconstruction_equal`, `nsubmix_reconstruction_equal`,
`parity_reconstruction_equal`, `synchronization_passed`,
`capture_validation_passed`, and `passed`. It first requires bytewise-equal
pre/post field hashes at every corresponding MAMAci subcycle plus bytewise-
equal downstream science outputs. For each enabled capture it then reconstructs
the locked branch decisions, actual `nsubmix`, and final parity from the
captured common old state and inputs and requires exact equality with the
enabled production values.

Direct hidden-branch tracing on the disabled arm is prohibited. The exact proof
premises are instead cumulative: the complete approved source diff proves that
the disabled executable invokes the original non-capture path; base/disabled
per-call field hashes prove that edit inert across the retained sequence;
disabled/enabled pre/post and downstream hashes prove equal observable science
at every boundary; and the enabled direct trace plus independent replay proves
its limiter, selector, zero-cloud, parity, and `nsubmix` decisions. All arms use
identical inputs, compiler configuration, backend, layout, seeds, process DAG,
subcycling, and sequence length. Any broken premise makes the comparison fail;
no final-step-only induction, implicit tolerance, inferred missing input, or
disabled-path capture trace is allowed.
`overhead` contains exactly `approved_fractional_bound`,
`diagnostic_timing_path`, `diagnostic_timing_sha256`, `observation_count`,
`warmup_block_count`, `measured_block_count`, `stopping_reason`,
`all_blocks_complete`, `seed_uint64`, `ratio_median`, `ratio_ci95_lower`,
`ratio_ci95_upper`, and `passed`. It binds the approved upper bound, raw matched
disabled/enabled complete-call observations in a separate diagnostic timing
artifact, the enabled/disabled ratio and its 95% interval, and passes only when
the interval's upper endpoint minus one is no greater than the approved bound.
The bound and three ratio values are nonnegative `real_hex` strings, the
observation count is positive, and the seed is an unsigned 64-bit integer.
The warmup and measured counts must reproduce every timing row,
`all_blocks_complete` must be true, and fewer than 30 measured matched blocks
fails the overhead gate.
The overhead interval uses the same SplitMix64, 10,000-replicate, median, and
percentile rules specified below with the exact qualification seed bound in
the qualification identity; it is qualification evidence only.
The overhead commands must set `internal_diagnostics_level=0` for both arms;
the separate sequence-equivalence commands use level 2 and cannot contribute a
time sample. The manifest binds both command sets and rejects a timing row from
the synchronized diagnostics execution.
Only the dedicated qualification validator may create this report, and it must
reject any non-`qualification` identity or completion. `record_completions` is
ordered by capture sequence and each closed entry contains exactly
`source_identity_sha256`, `completion_path`, `completion_sha256`, and `status`;
it must reproduce the enabled arm's complete expected population.
`qualification_passed` is true only when all comparisons, synchronization,
capture validations, completion accounting, and the overhead test pass.

`ndrop_capture_census_v1` contains exactly `artifact_schema`, `created_utc`,
`grid_role`, `census_identity_path`, `census_identity_sha256`,
`census_capture_approval_path`, `census_capture_approval_sha256`,
`observational_equivalence_approval_path`,
`observational_equivalence_approval_sha256`, `population_manifest_path`,
`population_manifest_sha256`, `count_accounting`,
`source_frequency_accounting`, `estimated_fe_work_accounting`,
`aggregation_blocked`,
`record_completions`, `stratum_counts`, `nsubmix_histogram`,
`census_diagnostics_path`, `census_diagnostics_sha256`,
and `census_sha256`.
`grid_role` must reproduce the population and identity manifests.
`count_accounting` contains exactly `total_census`, `available`,
`replay_eligible`, `resource_censored`, `unavailable`,
`unavailable_with_known_stratum`, `unavailable_with_unknown_stratum`, and
`unaccounted`. `available` is the exact sum of `replay_eligible` and
`resource_censored`; replay-eligible means `AVAILABLE` with
`1 <= production_nsubmix <= 512`, and resource-censored means `AVAILABLE`
with `production_nsubmix > 512`. Each weight-accounting object contains exactly
`known_total`, `unknown_record_count`, `available`, `replay_eligible`,
`resource_censored`, and `unavailable_known`. All weights are nonnegative
integers, `known_total=available+unavailable_known`, and
`available=replay_eligible+resource_censored`. An unknown-weight record is
counted only in `unknown_record_count`; no zero is substituted.
`record_completions` is ordered by `capture_sequence` and contains every
expected source identity and a closed object containing exactly
`source_identity_sha256`, `completion_path`, `completion_sha256`, `status`,
`source_stratum`, `source_frequency_weight`, and `estimated_fe_work_weight`
exactly once. Each `stratum_counts` entry contains exactly `dtmicro`,
`fe_regime`, `activation_regime`, `count_accounting`,
`source_frequency_accounting`, and `estimated_fe_work_accounting`; the nested
objects have the same closed keys and reconciliation rules above, except
`total_census` counts only records assignable to that stratum. Each
`nsubmix_histogram` entry contains exactly `production_nsubmix`,
`record_count`, `source_frequency_weight`, and `estimated_fe_work_weight`,
ordered by increasing `production_nsubmix`. Census sealing requires
`count_accounting.unaccounted=0`, exact count and weight reconciliation, role
`census`, and valid completion hashes. `aggregation_blocked` is true whenever
an unavailable record lacks a stratum or either weight, any weight
reconciliation fails, or another required aggregation input or approval is
absent or invalid.
`census_sha256` hashes canonical JSON with that field omitted. Qualification
records are rejected rather than omitted.

`ndrop_capture_census_diagnostics_v1` retains the amendment-required census
cross-tabulations and contains exactly `artifact_schema`, `created_utc`,
`census_identity_path`, `census_identity_sha256`,
`population_manifest_path`, `population_manifest_sha256`,
`completion_paths`, `completion_sha256s`, `population_scope`,
`stratum_diagnostics`, `overall_diagnostics`, and `diagnostics_sha256`.
It is produced and hashed before the census manifest, so it does not contain a
census-manifest hash; the census manifest binds its path and file hash, avoiding
a circular hash. Completion arrays are in capture-sequence order and reproduce
the census input set. `population_scope` is exactly
`available_before_resource_cap_with_exclusions_reported`.

Each `stratum_diagnostics` entry contains exactly `dtmicro`, `fe_regime`,
`activation_regime`, `count_accounting`, `source_frequency_accounting`,
`estimated_fe_work_accounting`, `nsubmix_histogram`, `limiter_activity`, and
`metric_summaries`, and is ordered by canonical timestep, FE-regime enum order,
then activation-regime enum order. `limiter_activity` contains exactly
`inactive_record_count`, `number_active_record_count`,
`mass_active_record_count`, and `either_active_record_count`.
`metric_summaries` has exactly `chi_mix`, `chi_act`, `chi_fe`,
`operator_norm`, `cloud_overlap_occupancy`, and `state_magnitude`. Each summary
contains exactly `available_count`, `minimum`, `maximum`, `source_weighted_mean`,
`fe_work_weighted_mean`, `source_weighted_p50`, `source_weighted_p90`,
`source_weighted_p99`, `fe_work_weighted_p50`, `fe_work_weighted_p90`, and
`fe_work_weighted_p99`; values are finite canonical `real_hex` or null only
when the corresponding count or positive known weight is zero. The nested
accounting and histogram objects use the exact census grammar above.
`overall_diagnostics` has exactly `count_accounting`,
`source_frequency_accounting`, `estimated_fe_work_accounting`,
`nsubmix_histogram`, `limiter_activity`, and `metric_summaries` with the same
semantics. These diagnostics contain no candidate result and are admissible
before selection; aggregate counts alone cannot replace them.

A production preregistration accepts only a census with
`grid_role=production_l128` and only source identities with `pver=128`.
`compatibility_l72` records may be validated and reported in a separate
compatibility census, but they have zero production-frequency and
method-selection admission: they cannot populate a production stratum, enter
pilot or expansion selection, receive production weights, or enter a timing,
error, Pareto, or recommendation aggregate. A file rename or an otherwise
valid `AVAILABLE`/`1..512` status cannot change that exclusion. L72 is neither
padded nor truncated to L128, and the Cartesian production stratum definition
is unchanged.

#### Population estimand, exclusions, and mixed-tail strata

The **approved OI-2 estimand resolution** was accepted through the
contract-only approval of document SHA-256
`a829ffc01d09a9960e9bc5df04a15199a39833fd6942cb696a8c942bb1292c38`.
It fixes the denominator that was not uniquely selected by the original
amendment: the timing/error distribution is conditional on the replay-eligible
production population, with every excluded population weight reported rather
than redistributed.

For each weight basis `b` (`source_frequency` or `estimated_fe_work`) and
production stratum `s`, define `R_s,b` as the weight of `AVAILABLE` L128 records
with `1 <= production_nsubmix <= 512`, `C_s,b` as the weight of `AVAILABLE`
records with `production_nsubmix > 512`, and `U_s,b` as independently known
weight of `UNAVAILABLE` records. Source-frequency record weight is one.
Estimated-FE-work record weight is observed `production_nsubmix` for available
records and is present for an unavailable record only when independently known;
no zero or inferred value is substituted. Define
`R_b=sum_s R_s,b`, `C_b=sum_s C_s,b`, `U_b=sum_s U_s,b`,
`A_b=R_b+C_b`, and `T_b=A_b+U_b`. Records with unknown stratum or weight are
counted separately and are in none of these sums.

The report gives `R_b`, `C_b`, `U_b`, `A_b`, and `T_b` for both bases;
`R_b/T_b` as replay-eligible coverage of known weight; `C_b/T_b` and `U_b/T_b`
as excluded known-weight fractions; and the count of unknown records. It must
not describe known-weight coverage as total-population coverage when any weight
is unknown. The two approved 1% tests remain exactly `C_b/A_b` for each basis;
an unavailable record is excluded from that denominator, a zero `A_b` is
invalid, and failure of either test prohibits an aggregate recommendation.

Let `S_s` be the deterministic selected production cases in stratum `s` that
the approved pilot-plus-expansion rule names, and let `n_s=|S_s|`. For a
preregistered case statistic `y_s,i`--the paired timing statistic for a method
pair or the deterministic error metric for a method--the stratum statistic is
`mean_s(y)=sum_(i in S_s)y_s,i/n_s`. For basis `b`, the aggregate numerator is
`N_b=sum_s R_s,b*mean_s(y)`, its denominator is `D_b=R_b`, and the aggregate is
`N_b/D_b`. Equivalently, every selected case in `s` has effective weight
`R_s,b/n_s`. Weighted quantiles and the unchanged hierarchical paired bootstrap
use those same effective weights and denominator; stress sentinels retain zero
weight. No `C_s,b` or `U_s,b` enters a timing/error numerator or denominator.

Thus a mixed extreme-tail stratum may contain replayable `65--512` records in
`R_s,b` and censored `>512` records in `C_s,b`; only the former define the
conditional estimand, while both are visible in census coverage and the 1%
tests. Known unavailable weights remain in `U_s,b` and coverage only. If any
unavailable record lacks its stratum or either required weight, the census may
seal for accounting but all method-selection aggregation remains blocked
pending another approved amendment.

If any case already named in `S_s` has a missing or insufficient eligible
method pair, the corresponding pairwise aggregate is `INDETERMINATE`; the case,
stratum, and weight are not dropped or reassigned. The report records missing
and insufficient effective weight as the sum of `R_s,b/n_s` over affected
selected cases, while retaining the whole stratum in `D_b` to demonstrate why
no numeric aggregate is emitted. Deterministic selection, at most three cases
per populated replay-eligible stratum, both censoring tests, and every
missing-pair/no-reweighting rule remain unchanged. Until this conditional
estimand is implemented under an approved addendum, bound in an approved
preregistration, and supported by authorized replay evidence, aggregate timing,
error, and method-selection outputs remain blocked.

### Stress-sentinel and replay run/resource schemas

`ndrop_replay_sentinel_v1` defines only an optional, zero-weight perturbation
of one valid census record and contains exactly `artifact_schema`,
`created_utc`, `sentinel_id`, `census_manifest_path`,
`census_manifest_sha256`, `parent_completion_path`,
`parent_completion_sha256`, `parent_source_identity_sha256`, `mutations`,
`target_branch`, `zero_census_weight`, and `sentinel_content_sha256`.
`sentinel_id` matches `[A-Za-z0-9._-]{1,128}`, `zero_census_weight` is true,
and the parent must be an available record in that exact census. `mutations`
is a nonempty ordered array of closed objects containing exactly `field`,
`flat_index`, `value_kind`, `real_hex`, `integer`, and `boolean`.
`value_kind` is exactly `real`, `integer`, or `boolean`; exactly its
corresponding value is non-null.
`flat_index` is a nonnegative integer or is null only for a scalar.

The mutable field allowlist is exactly `top_lev`, `dtmicro`, `qcld_initial`,
`interstitial_initial`, `cloudborne_initial`, `cldn`, `eddy_diff`, `zn`, `zs`,
`csbot`, `overlapp`, `overlapm`, `nact`, `mact`, `pdel`, `dz`, `pressure`,
`height`, `number_limiter_active`, and `mass_limiter_active`. Flattening uses
the capture-payload order. Observed `dtmin`, production directional rates,
reconstructed controller/selector diagnostics, `nsubmix`, parity, and final
state cannot be mutated or copied as sentinel truth; they are recomputed or
treated as unavailable as required by the locked replay gates. Duplicate
mutation targets are invalid. The fully mutated input must pass the same
shape, orientation, finite-value, range, no-flux, and mapping guards as a
production replay. `target_branch` is a nonempty diagnostic justification and
confers no expected pass status. `sentinel_content_sha256` hashes canonical
JSON with that field omitted. This deliberately supports perturbed sentinels
only; a wholly synthetic sentinel requires an approved addendum revision.

`ndrop_replay_expansion_decision_v1` is the immutable evidence for one proposed
expansion stage and contains exactly `artifact_schema`, `created_utc`, `run_id`,
`preregistration_path`, `preregistration_sha256`,
`preregistration_approval_path`, `preregistration_approval_sha256`,
`prior_execution_manifest_paths`, `prior_execution_manifest_sha256s`,
`stage_index`, `stratum_decisions`, `selected_new_replay_case_sha256s`, and
`expansion_decision_sha256`. Prior-manifest arrays are nonempty, equal length,
ordered pilot then prior expansion stages, and validated before any trigger is
evaluated. The self-hash omits only `expansion_decision_sha256`.

Each `stratum_decisions` entry contains exactly `dtmicro`, `fe_regime`,
`activation_regime`, `prior_selected_case_sha256s`,
`remaining_eligible_case_sha256s`, `trigger_evaluations`, `decision`,
`selected_next_hashes`, `unfilled_slots`, and `unfilled_reason`. Entries use
canonical stratum order; hash arrays are in selection-hash order. Every
`trigger_evaluations` entry contains exactly `trigger_id`, `family`,
`applicable`, `observed_kind`, `observed_real`, `observed_boolean`,
`observed_status`, `threshold_kind`, `threshold_real`, `threshold_boolean`,
`comparison`, and `fired`. `family` is `variation`, `reference`,
`materiality`, or `decision_sensitivity`. For `variation` and `materiality`,
`observed_kind=real`, `threshold_kind=real`, only `observed_real` and
`threshold_real` are non-null, and `comparison=ge`. For `reference`,
`observed_kind=status`, `observed_status` is exactly `certified`,
`pilot_reference_unresolved`, or `second_record_required`,
`threshold_kind=none`, every threshold value is null, and `comparison=eq`;
`fired` is the exact preregistered condition match. For
`decision_sensitivity`, `observed_kind=boolean`,
`threshold_kind=boolean`, only `observed_boolean` and
`threshold_boolean=true` are non-null, and `comparison=eq`. A trigger that is
not applicable has `observed_kind=none`, every observed value null,
`applicable=false`, and `fired=false`, while its threshold fields still
reproduce the preregistered rule where that family has a threshold. No other
kind/field combination is valid. `decision` is `expand` or `do_not_expand`.
`selected_next_hashes` contains at most the next two hashes and never raises
the stratum total above three. `unfilled_reason` is null when
`unfilled_slots=0`, otherwise exactly `no_remaining_eligible`,
`case_cap_reached`, or `required_pair_unavailable`. Trigger decisions and
observed values are retained even when no slot can be filled.

This artifact may inspect only already authorized execution manifests and the
sealed census/preregistration content. It is produced after the preceding
stage completes and before a new expansion run/resource approval. It cannot
authorize execution. The next expansion manifest and approval bind its file
hash; candidate files outside completed prior manifests are prohibited inputs.

`ndrop_replay_run_resource_v1` contains exactly `artifact_schema`,
`created_utc`, `run_id`, `preregistration_path`,
`preregistration_sha256`, `preregistration_approval_path`,
`preregistration_approval_sha256`, `execution_stage`, `stage_index`,
`expansion_decision_path`, `expansion_decision_sha256`,
`selected_replay_case_sha256s`, `machine`, `facility`, `project_account`,
`queue_or_qos`, `priority_authorized`, `nodes`, `tasks`, `cpus_per_task`,
`gpus_per_node`, `memory_bytes`, `walltime_seconds`,
`expected_allocation_cost`, `working_directory`, `input_paths`,
`input_sha256s`, `output_root`,
`expected_output_bytes`, `expected_output_inodes`, `available_bytes`,
`available_inodes`, `network_effects`, `exact_environment_commands`,
`exact_build_command`, `exact_run_commands`, `retention_policy`,
`data_authorization`, `facility_policy_review_path`,
`facility_policy_review_sha256`, `command_approval_paths`,
`command_approval_sha256s`, and `replay_run_resource_sha256`.
`execution_stage` is exactly `pilot` or `expansion`; `stage_index` is zero for
the pilot and a positive, contiguous integer for each approved expansion. The
two expansion-decision fields are null for the pilot and required for an
expansion; they must validate as `ndrop_replay_expansion_decision_v1`, match
the stage and selected hashes, and reference all preceding execution manifests.
The
selected case hashes are unique and ordered; `run_id` must equal the
preregistration `run_id`, and
each stage uses an exclusive output subdirectory under its approved output
root. A production replay uses its
`source_identity_sha256`, while a sentinel uses its
`sentinel_content_sha256`.
Pilot selection must reproduce the approved preregistration directly.
Expansion selection must reproduce its preregistered trigger and next-hash
rules from already authorized pilot results; a new bounded replay-execution
approval is required for each expansion manifest. Resource, capacity,
facility, command, priority, path, and self-hash semantics are identical to
`ndrop_run_resource_v1`. The manifest cannot authorize a case absent from the
approved preregistration or a sentinel absent from its `stress_sentinels`
array.

### Replay preregistration and result schemas

`ndrop_replay_preregistration_v1` contains exactly:

| Field | Required content |
|---|---|
| `artifact_schema`, `created_utc`, `run_id` | common values |
| `census_manifest_path`, `census_manifest_sha256` | exact sealed, non-blocked `production_l128` census |
| `implementation_identities` | closed objects named `replay_harness`, `forward_euler`, `imex`, `backward_euler`, `reference`, `timing_driver`, and `analysis`, each containing exactly `source_commit`, `source_tree`, and `executable_sha256` |
| `selection_hash` | exactly `sha256` |
| `selection_seed_uint64` | unsigned 64-bit integer encoded as 16 lowercase hexadecimal digits; interpreted as big-endian numeric text, not host bytes |
| `source_identity_encoding` | exactly the canonical JSON bytes defined above |
| `population_estimand` | closed object containing exactly `scope=conditional_replay_eligible_production`, `eligible_predicate=census_and_production_l128_and_available_and_nsubmix_1_to_512`, `source_frequency_denominator=sum_replay_eligible_source_weight`, `estimated_fe_work_denominator=sum_replay_eligible_fe_work_weight`, `within_stratum_rule=equal_selected_case_weight`, `censored_weight_treatment=coverage_only`, `unavailable_weight_treatment=coverage_only`, and `excluded_weight_treatment=no_redistribution` |
| `strata` | complete ordered list of closed objects containing exactly `dtmicro`, `fe_regime`, `activation_regime`, `count_accounting`, `source_frequency_accounting`, `estimated_fe_work_accounting`, `aggregate_source_frequency_weight`, and `aggregate_estimated_fe_work_weight`; aggregate weights are exactly `R_s,source` and `R_s,work` |
| `pilot_rule` | closed object containing exactly `algorithm=minimum_sha256_per_populated_stratum`, `eligibility=AVAILABLE_and_nsubmix_1_to_512_and_production_l128`, and `production_records_per_stratum=1` |
| `tail_rule` | closed object containing exactly `minimum_nsubmix=65`, `maximum_nsubmix=512`, `maximum_stress_sentinels=2`, and `sentinel_weight=0` |
| `stress_sentinels` | ordered array of zero to two closed objects containing exactly `path` and `file_sha256` for sealed `ndrop_replay_sentinel_v1` artifacts included in this preregistration |
| `expansion_rule` | closed object containing exactly `maximum_production_records_per_stratum=3`, `selection=next_smallest_sha256`, `combination=any_true`, `variation_triggers`, `reference_triggers`, `materiality_triggers`, and `decision_sensitivity_triggers`; nested grammar is fixed below and no trigger has a default |
| `timing_rng` | closed object containing exactly `algorithm=splitmix64_v1`, `seed_uint64`, `state_wrap_bits=64`, `index_rule=modulo`, `stream_order=method_pair_weight_basis_replicate_stratum_all_case_blocks_then_case_resample`, `stopping_replicates=10000`, and `stopping_stream_order=case_hash_then_block_count_then_method_then_replicate` |
| `bootstrap_replicates` | exactly `10000` |
| `median_rule`, `percentile_rule`, `weighted_quantile_rule` | closed rule objects with exact keys and values below |
| `tie_tolerance` | closed object with exactly `weighted_l1` and `linf`; each contains exactly nonnegative `absolute` and `relative` `real_hex` values fixed before execution |
| `diagnostic_failure_reasons` | canonical array exactly `input_nonfinite`, `assembly_nonfinite`, `factorization_unsafe_pivot`, `factorization_nonfinite`, `solve_nonfinite`, `raw_material_negative`, `clipping_nonfinite`, `projection_nonfinite`, `matrix_residual_nonfinite`, `matrix_residual_exceeded`, `physical_residual_nonfinite`, `physical_residual_exceeded`, `ledger_nonfinite`, `ledger_mismatch`, `boundary_mismatch`, `branch_mismatch`, `diagnostic_incomplete`, in that order; it may not be extended after result inspection |
| `missing_pair_rule` | exactly `indeterminate_no_reweighting` |
| `preregistration_sha256` | SHA-256 of canonical JSON with this field omitted |

Each `variation_triggers` entry contains exactly `trigger_id`, `metric`,
`statistic`, `lower_quantile`, `upper_quantile`, `comparison`, `threshold`, and
`enabled`. IDs are unique. `metric` is `chi_mix`, `chi_act`, `chi_fe`,
`operator_norm`, `cloud_overlap_occupancy`, or `state_magnitude`; `statistic`
is `range` or `quantile_spread`; quantiles are null for `range` and otherwise
canonical binary64 `real_hex` values satisfying
`0<=lower_quantile<upper_quantile<=1`; `comparison` is `ge`; `threshold` is a
nonnegative `real_hex`; and `enabled` is Boolean. The statistic is computed
from the complete replay-eligible census stratum before candidate inspection.

Each `reference_triggers` entry contains exactly `trigger_id`, `condition`, and
`enabled`; `condition` is `pilot_reference_unresolved` or
`second_record_for_failure_isolation`. Each `materiality_triggers` entry
contains exactly `trigger_id`, `weight_basis`, `denominator`, `comparison`,
`threshold`, and `enabled`; `weight_basis` is `source_frequency` or
`estimated_fe_work`, `denominator` is `complete_known_census_weight`,
`comparison` is `ge`, and `threshold` is a canonical `real_hex` in `[0,1]`.
Each `decision_sensitivity_triggers` entry contains exactly `trigger_id`,
`calculation`, `outcome`, and `enabled`; `calculation` is
`leave_one_stratum_out` or `worst_case_perturbation`, and `outcome` is
`pairwise_status_change` or `pareto_membership_change`. All four arrays are
ordered by `trigger_id`, contain no duplicate ID, and are evaluated completely;
the expansion decision is true if any enabled trigger is true. Numeric values
are fixed in the approved preregistration and are never inferred from candidate
results or implementation defaults.

`median_rule` contains exactly
`algorithm=sorted_middle`, `even_rule=lower_plus_half_difference`,
`nonfinite=reject`, and `arithmetic=ieee754_binary64_round_to_nearest_ties_even`.
`percentile_rule` contains exactly `algorithm=nearest_observation_bounds`,
`lower_probability`, `upper_probability`, `lower_index=floor_p_times_n_minus_1`,
`upper_index=ceil_p_times_n_minus_1`, and `interpolation=none`; its probabilities
are the canonical `real_hex` encodings of 0.025 and 0.975. The
`weighted_quantile_rule` contains exactly
`algorithm=first_cumulative_weight_ge_q_total`,
`sort_key=value_then_source_identity_sha256`, `zero_weight=exclude`,
`zero_total=invalid`, `q0=minimum`, `q1=maximum`, and `interpolation=none`.
These objects, the tie-tolerance values, and the RNG object are hashed parts of
the preregistration; changing any key or value invalidates it.

The preregistration is sealed first and is then approved by a separate
`ndrop_approval_v1` record whose `referenced_artifact_sha256` is the
byte-for-byte preregistration file hash; validation also reproduces the
artifact's internal `preregistration_sha256`. The preregistration cannot
contain or hash its own approval. Selection ranks the
byte string consisting of the eight big-endian bytes of
`selection_seed_uint64` followed directly by the canonical source-identity
JSON bytes. The smallest 32-byte SHA-256 digest wins; equal digests are broken
by lexicographically smaller canonical source-identity bytes. Candidate output
must not exist before the preregistration approval record is final. Any change
to a bound source, tree, executable, census, seed, threshold, or convention
invalidates the preregistration. Pilot and next-hash expansion selection use
only the exact `pilot_rule.eligibility` pool; a censored record cannot displace
an eligible production replay.

The timing generator is SplitMix64 with unsigned 64-bit wraparound. Starting
from the preregistered unsigned integer state, one draw performs, in order,
`state += 0x9e3779b97f4a7c15`,
`z=(state^(state>>30))*0xbf58476d1ce4e5b9`,
`z=(z^(z>>27))*0x94d049bb133111eb`, and returns
`z^(z>>31)`. An index into `N` entries is `draw % N`. The final hierarchical
bootstrap and the sequential stopping calculation each initialize a separate
SplitMix64 state to the exact same preregistered unsigned seed; neither carries
state into the other. Within the stopping calculation, state carries forward
through evaluations in the declared `stopping_stream_order`; a rejected or
missing evaluation consumes no draws. Within the final bootstrap, stream
consumption is fixed by method pair in lexicographic order, weight basis in the order
`source_frequency`, `estimated_fe_work`, bootstrap replicate 0--9999, stratum
key order, source-identity hash order, and timing-block order. Within each
replicate and stratum, first iterate every selected case in source-identity
hash order and draw with replacement exactly its observed number of complete
matched timing blocks. After all case ratios in that stratum are computed,
draw with replacement exactly its selected production-case count from the
same hash-ordered case list. Each drawn case occurrence carries
the stratum's frozen aggregate weight divided by its selected production-case
count; duplicate draws repeat that occurrence and its weight. A one-case
stratum draws that case on every case draw but still resamples its timing
blocks. If any selected production case in a populated stratum lacks a matched
eligible observation for the method pair, that pairwise aggregate is
`INDETERMINATE`; neither the missing case nor its stratum weight is
redistributed. Eligible per-case ratios may still be reported descriptively.

For a stratum with one selected production case, block resampling measures only
the timing variation observed for that case. Repeatedly drawing that same case
does not measure unobserved between-case variation or establish uniform
performance over the source stratum. Every interval, table, and conclusion for
such a stratum must carry the label
`single_selected_case_no_between_case_inference`; it may support the approved
pilot comparison for that selected case but cannot be described as evidence
that all records in the stratum have the same performance. This is a reporting
limitation only and does not change the approved sampling, expansion triggers,
bootstrap draws, weights, or decision rule.

For `N` sorted finite values, the median is element `N/2` when `N` is odd and
`lower + (upper-lower)/2` when `N` is even, where `lower` and `upper` are the
two central values. For 10,000 sorted bootstrap ratios, the 95% percentile
endpoints are elements
`floor(0.025*(N-1))` and `ceil(0.975*(N-1))`; there is no interpolation. A
weighted quantile sorts by `(value,source_identity_sha256)` and returns the
first value whose cumulative nonnegative weight is at least `q` times total
positive weight; `q=0` returns the minimum, `q=1` the maximum, and there is no
interpolation. Zero-total weight is invalid.

For either norm, method A is no worse than B exactly when
`e_A <= e_B + atol + rtol*max(abs(e_A),abs(e_B))`. It is strictly better only
when `e_A + atol + rtol*max(abs(e_A),abs(e_B)) < e_B`. The preregistered
`atol` and `rtol` values are applied casewise and normwise; they may not be
estimated from candidate results.

Each `ndrop_replay_result_v1` JSON object contains exactly `artifact_schema`,
`run_id`, `preregistration_path`, `preregistration_sha256`,
`preregistration_approval_path`, `preregistration_approval_sha256`,
`replay_run_resource_path`, `replay_run_resource_sha256`,
`replay_execution_approval_path`, `replay_execution_approval_sha256`,
`replay_case_sha256`, `source_identity`, `source_identity_sha256`,
`sentinel_path`, `sentinel_file_sha256`, `sentinel_content_sha256`,
`record_role`, `method`, `eligibility_scope`,
`case_status`, `gate_results`, `production_nsubmix`, `candidate_nsubmix`,
`reference_self_check`, `weighted_l1_by_field`, `linf_by_field`,
`maximum_weighted_l1`, `maximum_linf`, `candidate_state_path`,
`candidate_state_sha256`, `candidate_state_schema`, `diagnostics_path`,
`diagnostics_sha256`, `diagnostics_schema`, `primary_timing_path`,
`primary_timing_sha256`, `overhead_timing_path`, `overhead_timing_sha256`,
`tier3_diagnostic_timing_path`, and `tier3_diagnostic_timing_sha256`.
`record_role` must be
`census` for a production result or `STRESS_SENTINEL` for a zero-weight
sentinel. For a serialized replay-result row, `case_status` is exactly
`ELIGIBLE`, `UNAVAILABLE`, `NUMERICAL_FAILURE`, or
`REFERENCE_UNRESOLVED`. `RESOURCE_CENSORED` is an exact census/final-report
disposition, not a replay result: a record above 512 is retained in census
completion, histogram, coverage, and censoring evidence but cannot be selected
or acquire an `ndrop_replay_result_v1` row.
`eligibility_scope` is exactly `method_selection`, `robustness_only`, or
`none`: only an eligible census result has `method_selection`, an eligible
sentinel has `robustness_only`, and every non-eligible result has `none`.
Inapplicable paths, hashes, counts, and error arrays are `null`, never numeric
zero. The 51-element error arrays use `real_hex` strings in the fixed droplet,
interstitial-constituent, cloudborne-constituent order.

For a census result, `replay_case_sha256` equals `source_identity_sha256` and
all three sentinel fields are null. For a sentinel result, `source_identity` is its
parent census identity, all three sentinel fields are required, and
`replay_case_sha256` equals `sentinel_content_sha256`. The file hash must match
the preregistration path/`file_sha256` entry. Thus two sentinels derived from
one parent cannot collide with each other or with the production replay.

`method` is exactly `forward_euler`, `imex`, or `backward_euler`.
`gate_results` is a closed object with the keys `input`, `operator`, `finite`,
`raw_positivity`, `ledger`, `physical_residual`, `boundary`, `branch`,
`reference`, and `timing_path_match`. Each value is a closed object containing
exactly `applicability`, `outcome`, `passed`, `metric`, `threshold`, `count`,
and `message`. `applicability` is `required` or `not_applicable`; `outcome` is
`not_attempted`, `passed`, or `failed`; `passed` is the corresponding Boolean
only after an attempted gate and otherwise `null`; `metric` and `threshold`
are `real_hex` or `null`; `count` is a nonnegative integer or `null`; and
`message` is a nonempty deterministic explanation. A failed gate is not
omitted and a not-attempted gate does not fabricate a zero metric.

`reference_self_check` is always a closed object containing exactly
`applicability`, `attempt_status`, `outcome`, `failure_stage`,
`failure_reason`, `primary_rtol`, `tight_rtol`, `primary_completed`,
`tight_completed`, `maximum_normalized_disagreement`, `threshold`, and
`passed`. `applicability` is `required` or `not_applicable`;
`attempt_status` is `not_attempted` or `attempted`; and `outcome` is
`not_attempted`, `completed`, or `failed`. `failure_stage` is `null` unless the
outcome is failed, when it is exactly `primary_setup`, `primary_advance`,
`tight_setup`, `tight_advance`, or `pair_comparison`. `failure_reason` is null
only for a completed outcome and otherwise exactly `candidate_failed_first`,
`record_not_replayable`, `sentinel_excluded`, `nonfinite_stage`,
`step_underflow`, `step_limit_exhausted`, `incomplete_reference_pair`, or
`comparison_nonfinite`. `primary_completed` and `tight_completed` are
Booleans. `maximum_normalized_disagreement` and `passed` are non-null only when
both references completed and were comparable; otherwise both are null.
`primary_rtol`, `tight_rtol`, and `threshold` are always the canonical
`real_hex` encodings of the approved `1e-12`, `1e-13`, and `2e-11` values and
are never replaced by observed values. No failed or incomplete reference
requires a fabricated finite disagreement. The allowed combinations are
closed: `not_applicable` requires `attempt_status=not_attempted`,
`outcome=not_attempted`, both completion flags false, a null stage, reason
`record_not_replayable` or `sentinel_excluded`, and null metric/pass fields; a
candidate failure before reference requires `applicability=required`, the same
not-attempted/false/null values, and reason `candidate_failed_first`; an
attempted failure sets `attempt_status=attempted`, `outcome=failed`, completion
flags to the stages actually completed, one non-null failure stage/reason, and
null metric/pass fields; a completed comparison requires both completion flags
true, null stage/reason, a finite metric, and Boolean `passed`. The
`gate_results.reference` entry must reproduce this outcome and cannot claim
`passed` when `reference_self_check.passed` is null or false.

The result applicability and status precedence are exact:

| Role/status | Candidate/gates | Reference | State/diagnostics | Error arrays | Timing |
|---|---|---|---|---|---|
| Census `UNAVAILABLE` | Only reached validation gates are populated; later gates are not attempted | Not applicable, `record_not_replayable` | Null | Null | Null |
| Census `RESOURCE_CENSORED` | Census-only disposition; no replay-result row and no candidate execution | No replay-result object | No replay-result object | No replay-result object | No replay-result object |
| Census `NUMERICAL_FAILURE` | Candidate was attempted and its exact failed stage is retained | Not attempted, `candidate_failed_first` | State is null unless a complete valid state exists; diagnostics retain every completed stage and null later stages | Null | Any retained observation is `eligible=false` and excluded |
| Census `REFERENCE_UNRESOLVED` | Candidate gates passed | Attempted; failed/incomplete or completed above threshold | Complete candidate state and diagnostics required | Null | Any retained observation is `eligible=false` and excluded |
| Census `ELIGIBLE` | Every required gate passed | Completed, comparable, and passed | Complete candidate state and diagnostics required | Complete 51-field arrays | Only valid matched primary rows may be eligible |
| Sentinel `UNAVAILABLE` or `NUMERICAL_FAILURE` | Same failure representation as above | Not applicable, `sentinel_excluded` | As reached | Always null | Always null |
| Sentinel `ELIGIBLE` | Required robustness gates passed | Not applicable, `sentinel_excluded` | Complete candidate state and diagnostics required | Always null | Only optional Tier-3 diagnostic timing, never primary or overhead |

Population disposition precedence is `UNAVAILABLE`, then pre-execution
`RESOURCE_CENSORED`; neither enters candidate execution. Within an admitted
replay, precedence is `UNAVAILABLE`, an attempted candidate's
`NUMERICAL_FAILURE`, `REFERENCE_UNRESOLVED`, then `ELIGIBLE`. Because candidate
gates execute before the reference, a direct candidate failure leaves the
reference not attempted; a reference failure cannot relabel it. Sentinels may
not use `RESOURCE_CENSORED` or `REFERENCE_UNRESOLVED`. These rules preserve
unresolved-reference exclusion from every method preference and aggregate.

`candidate_state_schema` is `ndrop_replay_state_v1` whenever a state path is
present and `null` otherwise. The state payload has this byte-exact 17-byte
header: offsets 0--8 are the nine ASCII bytes `NDRPSTAT1`
(`4e 44 52 50 53 54 41 54 31`); offset 9 is the scalar width byte, exactly
`04` or `08`; offset 10 is endian byte `01` for little endian; offsets 11--16
are six zero reserved bytes. Payload data begin at offset 17. They store
full-`pver` `qcld`, 25 interstitial constituents, and 25 cloudborne
constituents in that order, with level fastest. The exact size is
`17 + scalar_width*51*pver`: 14,705 and 29,393 bytes for binary32 and binary64
at L72, and 26,129 and 52,241 bytes respectively at L128. Inactive production
levels are copied unchanged from the capture. The state hash covers the
complete header and payload bytes.

#### Numerical diagnostics companion artifact

The approved OI-5 amendment specifies the collection interface and supersedes
only the residual-unit/zero-normalizer sentence and raw-nonfinite failure-code
restriction below with its two explicit representation corrections. The other
companion definitions remain unchanged.

Final states and generic gate messages cannot recover pre-clipping,
pre-projection, matrix, or physical-residual evidence. Every attempted
candidate therefore has one separately hashed canonical-JSON
`ndrop_replay_diagnostics_v1` object containing exactly `artifact_schema`,
`run_id`, `preregistration_path`, `preregistration_sha256`,
`preregistration_approval_path`, `preregistration_approval_sha256`,
`replay_run_resource_path`, `replay_run_resource_sha256`,
`replay_execution_approval_path`, `replay_execution_approval_sha256`,
`replay_case_sha256`, `source_identity_sha256`, `record_role`, `method`,
`case_status`, `field_order`, `pair_order`, `completion`,
`total_ledger_delta`, `clipping_ledger_delta`, `operator_ledger_delta`,
`raw_minimum_by_field`, `roundoff_negative_count_by_field`,
`material_negative_count_by_field`, `fe_clipping_count_by_field`,
`fe_weighted_deficit_by_field`, `matrix_residual`, `physical_residual`, and
`diagnostics_sha256`. The self-hash omits only `diagnostics_sha256`; the result
artifact binds both its path and file hash and reproduces the internal hash.

`field_order` is exactly the 51-entry locked order: `qcld`, then the 25
interstitial constituents in `source_mapping.constituents` order, then those
same 25 cloud-borne constituents. `pair_order` is exactly that 25-entry
constituent order. Each entry is a closed object containing exactly `name`,
`phase`, `constituent_index`, `quantity_kind`, and `native_unit`.
`quantity_kind` is `droplet_number`, `aerosol_number`, or `aerosol_mass`;
`native_unit` is respectively `number_per_kg` or `kg_per_kg`.
`constituent_index` is null only for `qcld`; `phase` is `droplet`,
`interstitial`, or `cloudborne`. These arrays must reproduce the identity
manifest's canonical source mapping and cannot be supplied by a result.

`completion` contains exactly `status`, `expected_substeps`,
`completed_substeps`, `failure_stage`, and `failure_reason`. `status` is
`not_started`, `partial`, or `complete`; substep counts are nonnegative and
`completed_substeps<=expected_substeps`; both failure fields are null only for
`complete`. `failure_stage` is exactly `input`, `assembly`, `factorization`,
`solve`, `raw_state`, `clipping`, `projection`, `matrix_residual`,
`physical_residual`, `ledger`, `boundary`, `branch`, or `diagnostics`.
`failure_reason` is null only when status is `complete`; otherwise it is
exactly one entry from the preregistration's fixed
`diagnostic_failure_reasons` array. The allowed stage-to-reason mapping is
exactly: `input` to `input_nonfinite`; `assembly` to `assembly_nonfinite`;
`factorization` to `factorization_unsafe_pivot` or
`factorization_nonfinite`; `solve` to `solve_nonfinite`; `raw_state` to
`raw_material_negative`; `clipping` to `clipping_nonfinite`; `projection` to
`projection_nonfinite`; each residual stage to its same-named `nonfinite` or
`exceeded` reason; `ledger` to `ledger_nonfinite` or `ledger_mismatch`;
`boundary` to `boundary_mismatch`; `branch` to `branch_mismatch`; and
`diagnostics` to `diagnostic_incomplete`. No other pairing, unclassified
nonfinite value, or implementation-defined extension is valid.

The three ledger arrays have exactly 25 entries. Each entry is a finite
`real_hex` or null if the corresponding arithmetic did not complete. Number
constituents use `number_per_m2` and mass constituents use `kg_per_m2`, as
fixed by `pair_order`. They are whole-advance signed sums. For FE,
`clipping_ledger_delta` sums every per-substep raw-to-clipped change and
`operator_ledger_delta=total_ledger_delta-clipping_ledger_delta`; for IMEX and
BE clipping is prohibited and every completed clipping-ledger entry is the
canonical positive zero, not null, while an uncompleted entry is null.

Each per-field array has exactly 51 entries. Raw minima are the minimum finite
values seen after each method update and before any FE clip or common
zero-cloud projection, aggregated over all attempted substeps. Roundoff and
material-negative counts sum occurrences over levels and substeps using the
approved `-10*epsilon*field_scale` split. FE clipping counts likewise sum
every changed level/substep, and FE weighted deficits are the nonnegative
pressure-weighted sums of `clipped-raw`, in `number_per_m2` or `kg_per_m2` by
field. For IMEX and BE the two FE-only arrays contain 51 null entries, not
numeric zeros. Each uncomputed or nonfinite entry in any otherwise fixed-size
array is null and must be explained by `completion`; consumers may not replace
it with zero.

`matrix_residual` and `physical_residual` each contain exactly
`applicability`, `aggregation`, `numerator_by_field`,
`normalizer_by_field`, `normalized_by_field`, `maximum_normalized`, and
`worst_substep_by_field`. `applicability` is `required` or
`not_applicable`; `aggregation` is
`maximum_over_substeps_then_fields`. Each `*_by_field` array has 51 finite
`real_hex` or null entries, and `worst_substep_by_field` has 51 nonnegative
integers or null entries. Numerators and normalizers use each field's native
mixing-ratio unit; a completed normalizer is finite and positive and contains
only like-unit terms. Normalized values and their maximum are dimensionless.
FE matrix residual is not applicable and all its value fields are null. IMEX
and BE matrix residuals are required for solved fields; all methods require an
independently assembled physical residual using the approved equations. For
FE the physical residual is evaluated at every substep before clipping and the
stored per-field value is the maximum normalized residual; IMEX and BE use
their one-step physical residual. The common zero-cloud projection occurs
after these diagnostics and cannot erase their evidence.

Diagnostics are correctness evidence collected outside the primary
uninstrumented timing observation. Their construction, serialization, and
hashing are excluded from primary timing; no diagnostics artifact can be used
as a timing source.

Each `ndrop_replay_timing_v1` object is one primary uninstrumented observation
and contains exactly `artifact_schema`, `run_id`, `preregistration_path`,
`preregistration_sha256`, `preregistration_approval_path`,
`preregistration_approval_sha256`, `replay_run_resource_path`,
`replay_run_resource_sha256`, `replay_execution_approval_path`,
`replay_execution_approval_sha256`,
`replay_case_sha256`, `source_identity_sha256`, `method`, `timing_block`,
`interleaved_position`, `warmup_block_count`, `measured_block_count`,
`minimum_measured_blocks`, `maximum_measured_blocks`, `stopping_metric`,
`stopping_interval_lower`, `stopping_interval_upper`,
`stopping_relative_half_width`, `stopping_reason`, `eligible`,
`ineligibility_reason`, `seconds_per_column`, `state_checksum`,
`executable_sha256`, `compiler`, `machine`, `backend`, and `execution_space`.
Only `eligible=true` rows enter a ratio. A stress-sentinel row must set
`eligible=false`; primary sentinel rows are otherwise prohibited.

The production-replay timing policy retains the applicable previously approved
local timing sufficiency rule: exactly 10 complete untimed warmup blocks,
at least 30 and at most 100 complete measured blocks, and a fixed-seed
10,000-resample median interval check after each complete block from 30 onward.
For each case and method, relative half-width is
`(upper-lower)/(2*median)` using the preregistered percentile rule; times and
the median are strictly positive. Stop only after an entire interleaved block
when every method is at or below 0.05, or after block 100. Partial blocks are
discarded for every method. `stopping_metric` is exactly
`bootstrap_median_relative_half_width`; `stopping_reason` is
`precision_reached`, `maximum_blocks_reached`, or `execution_failed`.
The count and stopping fields are identical on every retained row for the same
case. This rule is about timing sufficiency, not the later hierarchical
population bootstrap, whose approved procedure is unchanged.

A method pair is `missing` for a selected case when it has zero complete
matched eligible blocks and `insufficient` when it has 1--29. Either makes the
pairwise population aggregate `INDETERMINATE` without reweighting. Thirty
through 100 complete blocks are sufficient; reaching 100 without the 5%
criterion is retained as `maximum_blocks_reached`, reported explicitly, and
the wider final interval governs the timing decision. `ineligibility_reason`
is null for eligible rows and otherwise `case_status`, `missing_pair`,
`insufficient_blocks`, `checksum_mismatch`, or `binding_mismatch`.

Each mandatory Tier-2 `ndrop_replay_overhead_timing_v1` JSON Lines object
contains exactly `artifact_schema`, `run_id`, `preregistration_path`,
`preregistration_sha256`, `replay_run_resource_path`,
`replay_run_resource_sha256`, `replay_case_sha256`,
`source_identity_sha256`, `method`, `overhead_block`,
`interleaved_position`, `arm`, `instrumentation_scope`,
`warmup_block_count`, `measured_block_count`, `minimum_measured_blocks`,
`maximum_measured_blocks`, `stopping_metric`, `stopping_interval_lower`,
`stopping_interval_upper`, `stopping_relative_half_width`, `stopping_reason`,
`all_blocks_complete`,
`seconds_per_column`, `state_checksum`, and `executable_sha256`. `arm` is
`uninstrumented` or `correctness_instrumented`; `instrumentation_scope` is
exactly `required_correctness_diagnostics_only`. Each block contains a matched
pair of arms and uses the same 10/30--100 complete-block policy, stopping
metric, fixed-seed interval rule, and three stopping reasons as primary timing.
The count and stopping fields are identical on both rows and every retained
block for the same case/method; `all_blocks_complete` is true or the overhead
artifact is invalid. It records only complete-advance overhead; it has no
phase fields and cannot enter a primary method ratio or preference.

Each optional Tier-3 `ndrop_replay_diagnostic_timing_v1` object contains exactly
`artifact_schema`, `run_id`, `preregistration_path`,
`preregistration_sha256`, `preregistration_approval_path`,
`preregistration_approval_sha256`, `replay_run_resource_path`,
`replay_run_resource_sha256`, `replay_execution_approval_path`,
`replay_execution_approval_sha256`,
`replay_case_sha256`, `source_identity_sha256`, `method`, `timing_block`,
`tier3_approval_path`, `tier3_approval_sha256`, `diagnostic_kind`,
`phase_status`,
`instrumented_seconds_per_column`, `setup_seconds_per_column`,
`factorization_seconds_per_column`, `solve_seconds_per_column`,
`combined_factorization_solve_seconds_per_column`,
`update_seconds_per_column`, `diagnostics_seconds_per_column`,
`fence_seconds_per_column`, `state_checksum`, and `executable_sha256`.
`phase_status` contains exactly `setup`, `factorization`, `solve`,
`combined_factorization_solve`, `update`, `diagnostics`, and `fence`; each is
`measured`, `method_inapplicable`, or `deferred_unmeasured`. A measured phase
has a nonnegative `real_hex`; the other two statuses require null. The complete
instrumented time is positive. A Tier-3 artifact requires its own later
hash-bound approval; absent that approval, all Tier-3 path/hash/count fields in
mandatory manifests are null/zero. Mandatory Tier-2 evidence never invents a
phase value or creates a Tier-3 artifact.

Timing blocks and positions are nonnegative integers; primary and instrumented
complete-advance seconds are positive `real_hex`; and checksums and executable
hashes are 64 lowercase hexadecimal characters. Result, diagnostics, primary
timing, overhead timing, and optional Tier-3 timing artifacts use canonical
JSON/JSON Lines with a final line feed and deterministic replay-case/method/
block ordering. Serialization occurs after measurement and outside every
timed region.

`ndrop_replay_execution_manifest_v1` is the authoritative completion record
for one bounded pilot or expansion execution and is written last, after the
result, correctness diagnostics, primary timing, Tier-2 overhead, and any
separately approved Tier-3 artifact have been flushed, closed, reopened,
validated, hashed, and atomically promoted. It contains exactly
`artifact_schema`, `created_utc`, `run_id`, `census_manifest_path`,
`census_manifest_sha256`, `preregistration_path`,
`preregistration_sha256`, `preregistration_approval_path`,
`preregistration_approval_sha256`, `replay_run_resource_path`,
`replay_run_resource_sha256`, `replay_execution_approval_path`,
`replay_execution_approval_sha256`, `result_path`, `result_sha256`,
`result_count`, `diagnostics_path`, `diagnostics_sha256`, `diagnostics_count`,
`primary_timing_path`, `primary_timing_sha256`, `primary_timing_count`,
`overhead_timing_path`, `overhead_timing_sha256`, `overhead_timing_count`,
`tier3_diagnostic_timing_path`, `tier3_diagnostic_timing_sha256`,
`tier3_diagnostic_timing_count`, `expansion_decision_path`,
`expansion_decision_sha256`, `execution_stage`, `stage_index`,
`selected_replay_case_sha256s`, `selection_accounting`,
`expansion_accounting`, `timing_accounting`, `status_counts`, and
`execution_manifest_sha256`.
Every selected replay-case hash must be represented by exactly one result per
method, all row bindings must equal the manifest bindings, all keys and counts
must be exact, and no qualification artifact may appear. The self-hash omits
only `execution_manifest_sha256`. Expansion-decision fields are null for the
pilot and required for an expansion. Tier-3 path/hash fields are null and count
zero unless separately approved. Until this record validates, all replay
result, diagnostics, and timing files are nonexistent to analysis; an execution failure may
be reported operationally but cannot produce an analyzable partial manifest.

In both execution and final manifests, `selection_accounting` is a closed
object containing exactly `populated_strata`, `expected_pilot_cases`,
`selected_production_cases`, `selected_stress_sentinels`, and
`missing_selected_cases`; each is a nonnegative integer and the missing count
must be zero. `expansion_accounting` is a closed object containing exactly
`triggered_strata`, `new_selected_cases`, `strata_at_case_cap`, and
`unfilled_triggered_slots`, all nonnegative integers. `status_counts` is a
closed object containing exactly `ELIGIBLE`, `UNAVAILABLE`,
`NUMERICAL_FAILURE`, `REFERENCE_UNRESOLVED`, and `RESOURCE_CENSORED`, each a
nonnegative count. In an execution manifest these are result-row counts,
`RESOURCE_CENSORED` is exactly zero, and their sum equals `result_count`; the
resource-censored population count remains in census/population coverage, not
in a replay result. A final manifest retains the same result-row accounting
and reports the census disposition separately through `population_coverage`.
`timing_accounting` is a closed object containing exactly
`selected_production_cases`, `cases_with_primary_timing`,
`cases_with_tier2_overhead`, `missing_method_pairs`,
`insufficient_method_pairs`, `maximum_blocks_without_precision`,
`eligible_method_pairs`, and `tier3_artifact_present`. Counts are nonnegative,
the Boolean Tier-3 field matches the nullable artifact binding, and all values
reproduce the timing rows. Missing and insufficient pairs are separate and
neither is silently reweighted.

`ndrop_replay_final_manifest_v1` contains exactly `artifact_schema`,
`created_utc`, `run_id`, `census_manifest_path`, `census_manifest_sha256`,
`preregistration_path`, `preregistration_sha256`,
`preregistration_approval_path`, `preregistration_approval_sha256`,
`replay_execution_manifest_paths`, `replay_execution_manifest_sha256s`,
`expansion_decision_paths`, `expansion_decision_sha256s`, `result_path`,
`result_sha256`, `result_count`, `diagnostics_path`, `diagnostics_sha256`,
`diagnostics_count`, `primary_timing_path`, `primary_timing_sha256`,
`primary_timing_count`, `overhead_timing_path`, `overhead_timing_sha256`,
`overhead_timing_count`, `tier3_diagnostic_timing_path`,
`tier3_diagnostic_timing_sha256`, `tier3_diagnostic_timing_count`,
`selection_accounting`, `expansion_accounting`, `timing_accounting`,
`population_coverage`, `status_counts`,
`analysis_source_commit`, `analysis_source_tree`, `analysis_executable_path`,
`analysis_executable_sha256`, `report_paths`, `report_sha256s`, and
`final_manifest_sha256`. It is valid only when every input hash revalidates,
every selected replay-case hash is accounted for, no qualification artifact is
present, all result/timing keys are unique, path and hash arrays have equal
length, all expansion decisions are bound in stage order, and all counts
reproduce the artifacts. Tier-3 fields are null/null/zero when no later Tier-3
approval exists. `final_manifest_sha256` hashes canonical JSON with that field
omitted.

`population_coverage` contains exactly `count_accounting`,
`source_frequency_accounting`, `estimated_fe_work_accounting`,
`source_frequency_censor_fraction`, `estimated_fe_work_censor_fraction`,
`missing_pair_source_weight`, `missing_pair_fe_work_weight`,
`insufficient_pair_source_weight`, `insufficient_pair_fe_work_weight`,
`source_frequency_replay_eligible_coverage`,
`estimated_fe_work_replay_eligible_coverage`,
`source_frequency_censored_known_fraction`,
`estimated_fe_work_censored_known_fraction`,
`source_frequency_unavailable_known_fraction`,
`estimated_fe_work_unavailable_known_fraction`,
`single_selected_case_strata`, and `estimand`. It reproduces the census and
approved preregistration without redistribution. Fractions and known weights
are canonical `real_hex`; an unknown required denominator or weight is null and
blocks method-selection reporting. `single_selected_case_strata` is the
canonical ordered array of stratum keys carrying the required bootstrap
limitation label. Under the approved OI-2 estimand, the two
`*_censor_fraction` values are `C_b/A_b`, replay-eligible coverage is
`R_b/T_b`, censored-known fraction is `C_b/T_b`, and unavailable-known fraction
is `U_b/T_b`; a zero denominator is invalid.
`estimand` is a closed object containing exactly `population`,
`eligible_predicate`, `source_weight_denominator`,
`fe_work_weight_denominator`, `stratum_case_weight_rule`,
`censored_weight_role`, `unavailable_weight_role`, and
`missing_pair_rule`. Their exact values are respectively
`conditional_replay_eligible_production`,
`census_and_production_l128_and_available_and_nsubmix_1_to_512`, `R_source`,
`R_fe_work`, `equal_within_stratum_then_census_stratum_weight`,
`coverage_only`, `coverage_only`, and `indeterminate_no_reweighting`.

The analysis implementation accepts replay inputs only through one or more
valid `ndrop_replay_execution_manifest_v1` records, never through raw result or
timing paths. The path and hash arrays must have equal nonzero length and be
ordered pilot first and then by approved expansion stage. The final manifest
must reproduce every execution-manifest hash and count, reject duplicate replay
cases across stages, and bind the deterministically merged result/timing
artifacts and diagnostics before binding the analysis executable and reports.
The analysis
source/tree/executable identity must equal the preregistered `analysis`
identity, and `report_paths` and `report_sha256s` must have equal length.

### Exact device-capture interfaces

The instrumentation descendant adds only the declarations below in namespace
`mam4::ndrop::capture_v1`. They are device-side copies into storage allocated
by the EAMxx writer; they perform no allocation, filesystem access,
serialization, physics, reduction, or candidate calculation. In the locked
source, `View1D`, `View2D`, `View3D`, and `ConstView2D` are respectively
`DeviceType::view_1d<Real>`, `view_2d<Real>`, `view_3d<Real>`, and
`view_2d<const Real>` (`src/mam4xx/ndrop.hpp:18-21`). `ColumnView` and
`ConstColumnView` are unmanaged rank-1 default-device views
(`src/mam4xx/mam4_types.hpp:59-64`); the default memory space is CUDA space
when CUDA is enabled and host space otherwise (lines 23-35). Physical layout is
the exact `DeviceType` layout recorded by the build manifest and is never
inferred to be contiguous.

Logical indexing is fixed by the production accesses: `raercol` and
`raercol_cw` are rank-3 `[level][buffer][constituent]` views with extents
`[pver][2][25]`; `nact` and `mact` are rank-2 `[level][mode]` views with
extents `[pver][4]`; and the wrapper's rank-4 aerosol work views are
`[column][level][buffer][constituent]`, from which a rank-3 column subview is
passed. Payload serialization may be mode-major or constituent-major as
specified above, but the copy must explicitly index these logical dimensions;
no array-of-view conversion is authorized or assumed.

The complete locked `dropmixnuc` parameter list is the one at
`src/mam4xx/ndrop.hpp:1279-1306`: team; `dtmicro`; `temp`, `pmid`, `pint`,
`pdel`, `rpdel`, `zm`; `state_q`; `ncldwtr`, `v_diffusivity`, `cldn`; the four
mode arrays `voltonumbhi_amode`, `voltonumblo_amode`, `exp45logsig`, and
`alogsig`; `aten`; the compatibility Boolean; `qcld`, `wsub`, `cldo`;
`qqcw_fld`, `ptend_q`, `tendnd`, `factnum`, `ndropcol`, `ndropmix`, `nsource`;
`int& nsubmix`; `wtke`, `ccn`, `coltend`, `coltend_cw`; `top_lev`; and work
views `raercol_cw`, `raercol`, `nact`, `mact`, `eddy_diff`, `zn`, `csbot`,
`zs`, `overlapp`, `overlapm`, `eddy_diff_kp`, `eddy_diff_km`, `qncld`,
`srcn`, `source`, `dz`, `csbot_cscen`, `raertend`, and `qqcwtend`, in exactly
that order and with the const/reference qualifications in the locked header.
The complete EAMxx wrapper list at
`components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:230-267` is, in
order: `team_policy`, `dt`, `dry_atmosphere`, `rpdel`, `kvh_mid`, `kvh_int`,
`wsub`, `cloud_frac`, `cloud_frac_prev`, `dry_aero`, `nlev`, `top_lev`,
`enable_aero_vertical_mix`, `coltend`, `coltend_cw`, `qcld`, `ndropcol`,
`ndropmix`, `nsource`, `wtke`, `ccn`, `ndrop_nsubmix`,
`ndrop_sum_nsubmix`, `ndrop_max_nsubmix`,
`reset_ndrop_nsubmix_accumulators`, `qqcw_fld_work`, `ptend_q`, `factnum`,
`tendnd`, `raercol_cw`, `raercol`, `state_q_work`, `nact`, `mact`, and the
15-entry `dropmixnuc_scratch_mem`, with the exact const/reference and view
ranks in that header. Its exact MAM4xx call mapping is lines 395-422. A future compile-only compatibility
test must instantiate both unchanged signatures and both capture overloads
against those headers after `CAPTURE_WRITER_APPROVED` and before qualification;
this revision does not compile it.

The declarations below correct the rank mismatch and encode the approved OI-1
contract, but they remain proposed addendum interfaces with no implementation
authority. The existing non-capture `dropmixnuc` signature remains available
and unchanged.

```cpp
namespace mam4::ndrop::capture_v1 {

using ByteView = DeviceType::view_1d<std::uint8_t>;
using IntView = DeviceType::view_1d<std::int64_t>;

struct DeviceSlot final {
  View1D real_payload;
  View1D qualification_limiter_trace;
  ByteView limiter_flags;
  IntView integer_metadata;
  std::size_t real_offset;
  std::size_t qualification_trace_offset;
  std::size_t flag_offset;
  std::size_t integer_offset;
  bool qualification_trace_enabled;
  bool enabled;
};

static constexpr std::size_t slot_integer_count = 10;

KOKKOS_INLINE_FUNCTION
void record_limiter_activity(const DeviceSlot& slot, int level, int mode,
                             bool number_limiter, Real pre_rate, Real cap,
                             Real post_rate);

KOKKOS_INLINE_FUNCTION
void snapshot_pre_update(
    const ThreadTeam& team, const DeviceSlot& slot, Real dtmicro, int top_lev,
    int initial_read_buffer, int initial_write_buffer,
    const ColumnView& qcld, const View3D& raercol,
    const View3D& raercol_cw, const ColumnView& cldn,
    const ColumnView& eddy_diff, const ColumnView& zn, const ColumnView& zs,
    const ColumnView& csbot, const ColumnView& pdel, const ColumnView& dz,
    const ColumnView& pressure, const ColumnView& height);

KOKKOS_INLINE_FUNCTION
void snapshot_selector_observation(
    const ThreadTeam& team, const DeviceSlot& slot, Real observed_dtmin,
    int production_nsubmix, const ColumnView& overlapp,
    const ColumnView& overlapm, const View2D& nact, const View2D& mact,
    const ColumnView& eddy_diff_kp, const ColumnView& eddy_diff_km);

KOKKOS_INLINE_FUNCTION
void preserve_post_reset_producers(
    const ThreadTeam& team, const DeviceSlot& slot, int production_nsubmix,
    int final_read_buffer, int final_write_buffer, const ColumnView& qcld,
    const View3D& raercol, const View3D& raercol_cw);

} // namespace mam4::ndrop::capture_v1
```

The four views in each enabled slot refer only to writer-owned allocation
pools. Their minimum extents are respectively
`real_offset+(123*pver+2)`,
`qualification_trace_offset+(qualification_trace_enabled ? 24*pver : 0)`,
`flag_offset+8*pver`, and `integer_offset+slot_integer_count`; an out-of-range
offset is `capacity_exhausted`. The qualification trace stores pre, cap, and
post arrays for number then mass in mode-major, level-fastest order and is
absent for census slots; an enabled census slot has a zero-extent trace view,
`qualification_trace_offset=0`, and
`qualification_trace_enabled=false`. The ten signed-64-bit metadata entries at
`integer_offset+[0..9]` are exactly `pre_update_complete`,
`initial_read_buffer`, `initial_write_buffer`,
`limiter_observation_complete`, `selector_observation_complete`,
`phase1_complete`, `phase2_complete`, `production_nsubmix`,
`final_read_buffer`, and `final_write_buffer`. Boolean entries are zero or one
and parity entries are zero or one. Device code may set only its assigned span.

The writer's `device_slots()` extent equals the launched column count and slot
index equals `team.league_rank()`/wrapper `icol`. A host-only immutable table
binds that index to its `SourceIdentity`; no device integer substitutes for the
identity hash. Capacity-exhausted or otherwise unavailable columns receive
`enabled=false`, all offsets equal `std::numeric_limits<std::size_t>::max()`,
and no device access. Allocations and index bindings live from
`prepare_slot` through caller Phase 2 completion, host copy, validation, and
record promotion; no span is reused until its completion record is final.

The locked producer type of `nsubmix` remains `int&`; no production API or
selector arithmetic changes. Phase 2 first validates `nsubmix>0`, then widens
it losslessly to signed 64-bit metadata and serialization. A reader rejects a
serialized value greater than `INT_MAX` before a checked narrowing conversion
to the existing offline `int`; negative, zero, overflowed, or non-round-tripping
values are `UNAVAILABLE`.

The capture-enabled overload is named `dropmixnuc_with_capture` and has the
same parameters, order, types, and return behavior as the locked
`mam4::ndrop::dropmixnuc`, followed by exactly one final parameter:

```cpp
const mam4::ndrop::capture_v1::DeviceSlot& capture_slot
```

This rank-correct interface encodes the approved OI-1 checkpoint contract but
is not current implementation authority. `snapshot_pre_update` is P1a;
`record_limiter_activity` is P1b and writes the direct strict-comparison flag,
plus pre/cap/post only when `qualification_trace_enabled`; and
`snapshot_selector_observation` is P1c after the existing reduction barrier and
`nsubmix` calculation. It stores observed `dtmin` and `nsubmix`, not a fabricated
maximum or scalar limiting level. `preserve_post_reset_producers` runs after the
existing barrier at locked-source line 1508 and copies fields 21--23 and final
integer metadata only after all substeps and the zero-cloud transition.
Every array-copy interface writes the full declared `pver` extent in the payload
order; `top_lev` controls replay validation and active-level interpretation, not
capture truncation.
Only P1c, after the barrier has completed every per-level P1b store, may set
`limiter_observation_complete`, `selector_observation_complete`, and then
`phase1_complete`, in that order.

The original `dropmixnuc` does not delegate through a runtime-enabled capture
path. The instrumented EAMxx descendant calls the original function when
capture is disabled and the separately approved capture overload only when
capture is enabled. No captured value feeds a production branch, state update,
tendency, limiter, return value, or buffer choice. No capture function adds a
reduction or barrier; it relies only on the source-anchored existing completion
points above. When `slot.enabled` is false, every capture function returns
before a load, store, barrier, atomic, allocation, copy, or file operation.

### Exact EAMxx writer interface

The host-side writer is confined to namespace
`scream::mam_aci::ndrop_capture_v1` and has this exact public interface:

```cpp
enum class CaptureRole : std::uint8_t { qualification, census };

struct SourceIdentity final {
  std::string case_id;
  std::string member_id;
  std::uint64_t physics_step;
  std::uint32_t mac_subcycle;
  std::uint32_t aci_call;
  std::uint64_t global_column_id;
  std::uint64_t capture_sequence;
  std::string grid_id;
  std::uint32_t pver;
};

struct CaptureConfig final {
  CaptureRole role;
  std::string run_id;
  std::filesystem::path run_root;
  std::filesystem::path identity_manifest;
  std::size_t record_capacity;
};

enum class CompletionCode : std::uint8_t {
  available,
  capacity_exhausted,
  phase1_copy_failed,
  limiter_observation_failed,
  selector_observation_failed,
  post_reset_preservation_failed,
  phase1_completion_failed,
  phase2_producer_incomplete,
  phase2_copy_failed,
  phase_mismatch,
  ownership_failed,
  validation_failed,
  serialization_failed,
  hash_failed,
  promotion_failed
};

class CaptureWriter final {
 public:
  explicit CaptureWriter(const CaptureConfig& config);
  mam4::ndrop::capture_v1::DeviceSlot prepare_slot(
      const SourceIdentity& source_identity,
      std::string droplet_buffer_id,
      std::string interstitial_buffer_pair_id,
      std::string cloudborne_buffer_pair_id);
  Kokkos::View<const mam4::ndrop::capture_v1::DeviceSlot*, DefaultDevice>
  device_slots() const noexcept;
  void complete_phase2_after_dropmixnuc(
      const DefaultDevice::execution_space& execution_space);
  void mark_unavailable(const SourceIdentity& source_identity,
                        CompletionCode code) noexcept;
  void finalize_records();
};
```

Construction validates the sealed capture identity and exclusively allocates
the fixed-capacity capture-owned device and host buffers. `prepare_slot` is a
host operation performed before the column launch and never retains a
production view. Capacity exhaustion returns an unavailable slot and must not
overwrite or wrap. `complete_phase2_after_dropmixnuc` is called immediately
after the existing post-`call_function_dropmixnuc` completion point. It
verifies device completion, copies only enabled capture slots into
capture-owned host storage, waits for that enabled-path copy, validates both
phases, and performs the transaction above. `finalize_records` closes and
revalidates every expected per-record completion transaction and fails if any
is absent; it does not write a qualification report or census manifest. Those
role-level artifacts are created only by their dedicated offline validators.
Destruction is not a success path and never implicitly seals a run.

The only runtime activation interface is the environment variable
`EAMXX_NDROP_CAPTURE_IDENTITY`. When it is absent, capture is disabled before
any capture object, slot, allocation, copy, barrier, fence, or file operation
is created and the original wrapper is called. When present, its value is the
absolute path to the already sealed `ndrop_capture_identity_v1`; all role,
capacity, cadence, output-root, seed, and source-mapping settings come from
that validated manifest and cannot be overridden by other environment or
namelist values. The instrumentation descendant's CMake file may select its
separately retained clean MAM4xx source only through the cache path
`EAMXX_NDROP_CAPTURE_MAM4XX_SOURCE_DIR`; an empty value uses the unchanged
gitlink checkout. The nonempty path must validate against
`instrumented_mam4xx_identity` before configuration succeeds. These two names
are instrumentation interfaces only and do not change a production physics
API or namelist.

The approved addendum fixes the caller-side Phase 2 sequence below. It may be
implemented in production only if `CAPTURE_WRITER_APPROVED` later names the
exact source implementation. The existing global fence at
`eamxx_mam_aci_process_interface.cpp:595` is the proposed
producer-completion guarantee; `copy_execution_space` is the instance used for
the enabled device-to-host copy and its completion fence:

```cpp
call_function_dropmixnuc_with_capture(
    /* the unchanged call_function_dropmixnuc arguments */,
    capture_writer.device_slots());
Kokkos::fence();  // existing production completion point
capture_writer.complete_phase2_after_dropmixnuc(
    copy_execution_space);
```

`call_function_dropmixnuc_with_capture` has the same parameter list, order,
types, and return behavior as the locked EAMxx
`call_function_dropmixnuc`, followed by exactly one final parameter:

```cpp
Kokkos::View<const mam4::ndrop::capture_v1::DeviceSlot*, DefaultDevice>
    capture_slots
```

The disabled caller invokes the original `call_function_dropmixnuc`; it does
not construct a slot view or enter the capture-enabled overload.

For a backend on which the existing global fence is not a proven guarantee for
every execution-space instance that produced a slot, this interface is
unsupported and must return for addendum revision rather than adding an
unapproved blanket synchronization or wrapper parameter. The subsequent
device-to-host copy must be synchronous or followed by a
`copy_execution_space` fence before host validation. No additional fence,
event, allocation, copy, or file operation is permitted on the disabled
branch. A backend whose completion mapping cannot be proven by a focused test
is unsupported and cannot enter qualification.

### Exact offline validator and replay interfaces

The offline C++ implementation reuses `ndrop_test::State`,
`ndrop_test::OperatorData`, and the already approved FE, IMEX, BE, and reference
entry points. It adds these exact types and signatures in namespace
`ndrop_test::production_replay_v1`:

```cpp
enum class EvidenceRole : std::uint8_t { qualification, census };

struct SourceIdentity final {
  std::string case_id;
  std::string member_id;
  std::uint64_t physics_step;
  std::uint32_t mac_subcycle;
  std::uint32_t aci_call;
  std::uint64_t global_column_id;
  std::uint64_t capture_sequence;
  std::string grid_id;
  std::uint32_t pver;
};

enum class ValidationCode : std::uint8_t {
  ok,
  unavailable,
  schema_mismatch,
  role_mismatch,
  provenance_mismatch,
  completion_missing,
  hash_mismatch,
  shape_mismatch,
  unit_mismatch,
  orientation_mismatch,
  range_failure,
  nonfinite_value,
  phase_mismatch,
  ownership_failure,
  operator_mismatch,
  nsubmix_mismatch
};

struct ValidatedCapture final {
  SourceIdentity source_identity;
  EvidenceRole role;
  int pver;
  int top_lev;
  int production_nsubmix;
  int initial_read_buffer;
  int final_read_buffer;
  Real dtmicro;
  State full_initial_state;
  State full_production_final_state;
  std::vector<Real> cldn;
  std::vector<Real> eddy_diff;
  std::vector<Real> zn;
  std::vector<Real> zs;
  std::vector<Real> csbot;
  std::vector<Real> overlapp;
  std::vector<Real> overlapm;
  std::array<std::vector<Real>, num_modes> nact;
  std::array<std::vector<Real>, num_modes> mact;
  std::vector<Real> pdel;
  std::vector<Real> dz;
  std::vector<Real> pressure;
  std::vector<Real> height;
  std::vector<Real> production_eddy_diff_kp;
  std::vector<Real> production_eddy_diff_km;
  std::array<std::vector<std::uint8_t>, num_modes> number_limiter_active;
  std::array<std::vector<std::uint8_t>, num_modes> mass_limiter_active;
  Real observed_dtmin;
};

struct CaptureValidationResult final {
  ValidationCode code;
  std::string message;
  std::optional<ValidatedCapture> value;
};

struct CensusAdmissionResult;

class CensusReplayInput final {
 public:
  const ValidatedCapture& capture() const noexcept;
  const std::string& census_manifest_sha256() const noexcept;
  const std::string& completion_sha256() const noexcept;
  const std::string& identity_manifest_sha256() const noexcept;

 private:
  friend CensusAdmissionResult admit_census_replay_input(
      const std::filesystem::path&, const std::filesystem::path&,
      const std::filesystem::path&, const std::filesystem::path&);
  CensusReplayInput(ValidatedCapture, std::string, std::string, std::string);
  ValidatedCapture capture_;
  std::string census_manifest_sha256_;
  std::string completion_sha256_;
  std::string identity_manifest_sha256_;
};

struct CensusAdmissionResult final {
  ValidationCode code;
  std::string message;
  std::optional<CensusReplayInput> value;
};

struct ArtifactValidationResult final {
  ValidationCode code;
  std::string message;
};

struct ReplayCase final {
  std::string replay_case_sha256;
  std::optional<std::string> stress_sentinel_id;
  SourceIdentity source_identity;
  State initial_state;
  std::optional<State> production_final_state;
  OperatorData operator_data;
  Real dtmicro;
  int production_top_lev;
  std::optional<int> observed_production_nsubmix;
  int replay_nsubmix;
  std::optional<int> production_final_buffer;
  std::vector<std::uint8_t> limiter_activity;
};

struct OperatorGateResult final {
  ValidationCode code;
  std::string message;
  std::optional<ReplayCase> replay_case;
  Real chi_mix;
  Real chi_act;
  Real chi_fe;
  Real reconstructed_controller_rate;
  std::optional<Real> reconstructed_selector_max_tinv;
  int reconstructed_selector_candidate_count;
  std::vector<int> reconstructed_selector_minimizer_levels;
  int reconstructed_nsubmix;
};

CaptureValidationResult read_capture(
    const std::filesystem::path& completion_path,
    EvidenceRole required_role,
    const std::filesystem::path& identity_manifest);

CensusAdmissionResult admit_census_replay_input(
    const std::filesystem::path& completion_path,
    const std::filesystem::path& census_manifest,
    const std::filesystem::path& census_capture_approval,
    const std::filesystem::path& observational_equivalence_approval);

OperatorGateResult reconstruct_and_certify_locked_operator(
    const CensusReplayInput& input,
    const std::string& expected_mam4xx_commit,
    const std::string& expected_adapter_source_sha256);

OperatorGateResult construct_and_validate_stress_sentinel(
    const std::filesystem::path& sentinel_path,
    const CensusReplayInput& admitted_parent,
    const std::filesystem::path& preregistration,
    const std::string& expected_adapter_source_sha256);

ArtifactValidationResult validate_qualification_only(
    const std::filesystem::path& qualification_report,
    const std::filesystem::path& qualification_approval);

ArtifactValidationResult validate_census_only(
    const std::filesystem::path& census_manifest,
    const std::filesystem::path& census_capture_approval,
    const std::filesystem::path& observational_equivalence_approval);

ArtifactValidationResult validate_replay_preregistration(
    const std::filesystem::path& preregistration,
    const std::filesystem::path& census_manifest,
    const std::filesystem::path& preregistration_approval);

ArtifactValidationResult validate_replay_run_resource(
    const std::filesystem::path& replay_run_resource,
    const std::filesystem::path& preregistration,
    const std::filesystem::path& replay_execution_approval);

ArtifactValidationResult validate_replay_execution_manifest(
    const std::filesystem::path& replay_execution_manifest,
    const std::filesystem::path& replay_run_resource,
    const std::filesystem::path& replay_execution_approval);
```

`read_capture` returns a value only for a complete matching role, identity,
schema, payload, and completion transaction. The operator adapter maps
production level `top_lev+j` to replay level `j`; constructs no-flux interface
endpoints; maps `eddy_diff_km` to the upper rate, `eddy_diff_kp` to the lower
rate, `overlapm` to upper overlap, and `overlapp` to lower overlap; and uses the
existing pressure-weight, activation, state, zero-cloud, and controller
contracts unchanged. It verifies every reconstructed rate, direct limiter flag,
observed `dtmin`, controller value, selector candidate set, `chi` diagnostic,
parity, and integer `nsubmix` before returning a `ReplayCase`.
Same-type/same-order comparisons are bitwise. A separately approved
roundoff-only tolerance is required before any different-order comparison;
absent that approval the comparison fails. The minimizing-level vector is an
ordered offline reconstruction and may contain zero, one, or multiple levels;
no scalar production limiting level exists in `ValidatedCapture`. The
reconstructed selector maximum is null exactly when the candidate count is
zero; no numeric zero or finite placeholder is fabricated for an empty set.

`ValidatedCapture` is a schema-validation value only and cannot enter a
production replay directly. `admit_census_replay_input` is the sole constructor
of `CensusReplayInput`: it requires role `census`,
`grid_role=production_l128`, membership of the exact completion hash in the
sealed census, matching census and observational-equivalence approvals, and
revalidation of the identity, provenance, payload, and completion hashes. The
wrapper owns its validated data and exposes only const access. Qualification
records, compatibility-L72 records, caller-created aggregates, caller-mutated
copies, or a bare `ValidatedCapture` cannot satisfy the
`reconstruct_and_certify_locked_operator` signature. Stress sentinels are
constructed only from a separately admitted census parent and retain that
parent's bindings plus the sentinel hash.

`ReplayCase::initial_state` and a production replay's populated
`ReplayCase::production_final_state` contain exactly the active
`[top_lev,pver)` slices in the locked 51-field order;
`ValidatedCapture` retains the full-`pver` states without normalization or
truncation. `ReplayCase::limiter_activity` is the flattened number flags in
mode-major, level-fastest order followed by the mass flags in that same order.
For a production replay, `replay_case_sha256` is its source-identity hash and
`stress_sentinel_id` is empty, both production optionals are populated, and
`replay_nsubmix` equals the observed production count. The sentinel constructor
accepts only a preregistered `ndrop_replay_sentinel_v1`, applies its ordered
mutations to its validated census parent, recomputes all derived operator
values, assigns the sentinel content hash as `replay_case_sha256`, and returns
no production final
state, observed production count, final parity, or operator-equality
certification; `replay_nsubmix` is the locked-source count reconstructed for
the perturbed input and is never labeled as an observed production value. Such
a case can exercise only the unchanged analytic, finite, positivity, residual,
boundary, branch, and
robustness gates; it cannot enter census, reference-error, timing-ratio, or
method-selection aggregates.

`validate_qualification_only` rejects census artifacts and can emit only the
qualification report fields. `validate_census_only` rejects qualification
artifacts and verifies that census approval names the exact qualified enabled
implementation, executable/build identity, backend, compiler, layout, seeds,
and equivalence-relevant capture settings. Neither function may relabel an
artifact or repair a missing hash. The three replay validators require the
exact preceding approval and manifest hashes, reject qualification artifacts,
and do not accept a raw result or timing path in place of its execution
completion manifest.

### Offline diagnostic-interface amendment (OI-5; approved)

**APPROVED** at document SHA-256
`2300c9ef1448fb161042038e0075d8646cdfc733624ba2c0f7160e21b2f2fbe0`
for `ADDENDUM_APPROVED` to `ADDENDUM_APPROVED`. This is a narrow revision
of the addendum approved at `ca30c32c85451ff2fa14666119cfa7da8f665c9d5f2eb662c6a40574b513f470`.
The pre-edit document bytes match bookkeeping SHA-256
`8e4e84b45470af8c21fe568b54426817f0a0f1b7e84be0d7bc7e47ccadf3a23a`.
Read-only source evidence below is anchored to E3SM commit
`19c663a61ef73f99eb0b001d2947fd24b06e5e7b`; paths in this subsection are
relative to `docs/analysis/evidence/ndrop_integrator/`. No implementation or
test execution is recorded by this amendment.

#### Requirement-to-capability audit

| Approved requirement | Existing capability and source anchor | Demonstrated gap | Smallest mandatory change |
|---|---|---|---|
| Signed 25-pair total/clipping/operator ledgers; FE raw minima, clipping deficits/counts; all-method negative classifications | `AdvanceDiagnostics` and `FieldClipping`, `ndrop_integrator.hpp:60-106`; FE accumulation at `ndrop_integrator.cpp:643-843`; BE minima at `ndrop_backward_euler.cpp:432-462` | Successful-call values already exist; zero-initialized arrays do not identify completed entries after an early return | Reuse values; add validity/completion metadata only |
| All 51 pre-projection raw minima | IMEX `update_positivity`, `ndrop_imex.cpp:261-289`, retains a global solution minimum and per-field negative counts | IMEX per-field raw minima are not filled, and zero-cloud projection at lines 499-502 destroys the corresponding final-state evidence | Record field minima in the existing positivity traversal only when the new detail is requested |
| Per-field matrix/physical residual numerator, normalizer, normalized maximum and worst substep | FE `ndrop_integrator.cpp:115-193`; IMEX `ndrop_imex.cpp:94-240`; BE `ndrop_backward_euler.cpp:260-308,405-429,485-570` | Helpers retain only scalar maxima; FE raw states are overwritten each substep; implicit states are subsequently projected | Copy existing intermediate residual quantities into fixed-size per-field maxima before destruction; retain substep indices, not state trajectories |
| Complete and partial diagnostic artifacts | Existing early returns and `AdvanceResult::message`; `AdvanceDiagnostics` has no completion masks | Uncomputed zero/default entries can masquerade as evidence; existing schema lacks a reason for a nonfinite raw FE update | Add typed completion, masks, and the narrowly specified `raw_nonfinite` reason |
| Faithful residual normalization | FE/IMEX/BE physical residuals and IMEX matrix residuals use dimensionless row-scaled terms; BE matrix residuals use native-unit row terms; all helpers explicitly handle zero scale | The approved companion's universal native-unit/strictly-positive-normalizer sentence cannot describe these calculations exactly | Amend only that representation sentence as specified below; preserve every computed normalized value and gate |
| Complete reference-failure reporting | `ndrop_reference.cpp:285-440` returns `ok`, deterministic messages, and completed counters; `test_cross_scheme.cpp:934-989` already checks several exact messages | The addendum does not specify a wrapper mapping; lack of typed fields alone is not a source-interface blocker | Specify the exhaustive, source-bound wrapper mapping below; no reference source edit |

Optional structured reference statuses, general observers/callbacks, full
substep trajectories, and diagnostic refactoring are excluded. The retained
maxima/counts and failure metadata serve the already-approved auditability
requirement; no new scientific diagnostic or tolerance is introduced.

#### Minimal C++ extension and serialization correspondence

Append the following types in `ndrop_test` before `AdvanceDiagnostics` and
append the stated members to the existing structures. Existing function
signatures, member meanings, defaults, and entry points remain available.
Only `<cstdint>` and `<optional>` are additionally needed in this header.

```cpp
enum class ReplayDiagnosticStatus : std::uint8_t {
  not_started, partial, complete
};
enum class ReplayDiagnosticStage : std::uint8_t {
  input, assembly, factorization, solve, raw_state, clipping, projection,
  matrix_residual, physical_residual, ledger, boundary, branch, diagnostics
};
enum class ReplayDiagnosticReason : std::uint8_t {
  input_nonfinite, assembly_nonfinite, factorization_unsafe_pivot,
  factorization_nonfinite, solve_nonfinite, raw_material_negative,
  raw_nonfinite, clipping_nonfinite, projection_nonfinite,
  matrix_residual_nonfinite, matrix_residual_exceeded,
  physical_residual_nonfinite, physical_residual_exceeded, ledger_nonfinite,
  ledger_mismatch, boundary_mismatch, branch_mismatch, diagnostic_incomplete
};
struct ReplayResidualMaximum {
  Real numerator = 0;
  Real normalizer = 0;
  Real normalized = 0;
  std::uint32_t substep = 0;
};
struct ReplayDiagnosticDetail {
  ReplayDiagnosticStatus status = ReplayDiagnosticStatus::not_started;
  std::uint32_t expected_substeps = 0;
  std::uint32_t completed_substeps = 0;
  std::optional<ReplayDiagnosticStage> failure_stage;
  std::optional<ReplayDiagnosticReason> failure_reason;
  std::array<bool, 51> raw_minimum_valid{};
  std::array<bool, 51> negative_counts_valid{};
  std::array<bool, 51> fe_clipping_count_valid{};
  std::array<bool, 51> fe_weighted_deficit_valid{};
  std::array<bool, 25> total_ledger_valid{};
  std::array<bool, 25> clipping_ledger_valid{};
  std::array<bool, 25> operator_ledger_valid{};
  std::array<std::optional<ReplayResidualMaximum>, 51> matrix_residual{};
  std::array<std::optional<ReplayResidualMaximum>, 51> physical_residual{};
};
// Append to AdvanceDiagnostics:
// std::optional<ReplayDiagnosticDetail> replay_detail;
// Append to AdvanceControl, after the existing members:
// bool collect_replay_diagnostics = false;
```

Detail is engaged only when both `collect_detailed_diagnostics` and
`collect_replay_diagnostics` are true. All other combinations leave it
disengaged and execute no added collection work. The replay wrapper rejects a
request for a correctness artifact with either flag false before calling the
integrator. This uses the existing diagnostics control and one opt-in flag;
it does not introduce another advance, solver, callback, or observer API.
The returned detail owns only fixed-size scalar metadata. It retains no view,
pointer to solver storage, or full-state/substep history.

Indices are `0=qcld`, `1+j=interstitial[j]`, and
`26+j=cloudborne[j]`, for `j=0..24`. The serializer takes existing droplet
and `FieldClipping` members for minima, negative counts, FE clipping counts,
and weighted deficits, and existing ledger arrays for the three signed
ledgers. It emits a value only when the corresponding validity flag is true;
otherwise it emits null. `negative_counts_valid` covers both already-computed
roundoff/material counters together. IMEX's added minima use the same
finite raw solution and field traversal as its existing positivity check;
they never feed that check. IMEX/BE FE-only flags remain false and serialize
as 51 null entries. Completed implicit clipping-ledger entries are actual
zero contributions, with true validity, as already required.

Each residual optional maps directly to the corresponding four serialized
per-field arrays (`numerator`, `normalizer`, `normalized`, `worst_substep`).
An empty optional produces null in all four arrays. FE's matrix array remains
empty and `applicability=not_applicable`; implicit matrix arrays and all
physical arrays are required. A completed zero residual is an engaged tuple,
not an empty optional. `maximum_normalized` is the maximum of completed
per-field values, null if none completed. Complete artifacts require every
applicable entry; partial artifacts cannot certify a gate or timing row.

FE records the physical tuple at every existing substep residual call before
clipping; IMEX/BE record their one-step tuples before zero-cloud projection.
Keep the tuple with the largest normalized value, including the numerator
and normalizer from that same evaluation. Exact ties keep the earliest
zero-based substep; IMEX and BE therefore use index zero. IMEX cloud RHS
0 maps to field 0, cloud RHS `j+1` to `26+j`, and interstitial RHS `j` to
`1+j`. BE block component `x` maps to `26+j`, `y` to `1+j`, and the scalar
droplet solve to 0. The original matrix-family and solve ordering is unchanged.

The status/count fields map directly to serialized `completion`. Expected
substeps are the existing selected FE count or one for IMEX/BE, zero only
before that count is established; a wrapper with a known admitted FE count
must retain it even if setup fails. Increment completed substeps only after
all 51 field updates and required per-step clipping have completed, before
the separate zero-cloud projection. Failure during projection or final ledger
work may therefore have `completed_substeps=expected_substeps` but remains
`partial`. `not_started` means no candidate update/solve began; `complete`
requires all required diagnostics and gates to finish successfully. A failed
gate with fully computed metrics retains those metrics and has `partial`
status and its actual stage/reason. Every returned incomplete detail has
non-null stage/reason; enum spellings map exactly to the companion grammar.

Set validity only after the corresponding arithmetic completes, at its
existing producer. On an early return, retain finite aggregates of completed
evaluations and leave unreached entries null; these are explicitly partial
evidence, never whole-advance results. If an attempted contribution to an
aggregate is nonfinite or incomplete, invalidate that affected aggregate
rather than serialize its stale prefix as complete. Independent completed
fields/stages remain available. No infinity sentinel or initialized zero is
evidence of completion. Ledger masks are per pair and per decomposition term,
so a later pair failure does not erase completed earlier pairs. Helper-local
successes must be copied before a subsequent field, pair, or stage can return.
Invalid direct-call options/shapes, which the replay admission gate excludes,
use `diagnostics/diagnostic_incomplete` without inventing a scientific failure;
the existing `AdvanceResult::ok`, message, and return behavior remain unchanged.

Existing source failure sites populate the stage/reason corresponding to
their actual computation. Preserve the first failing stage and the existing
short-circuit order. If the existing final residual gate combines matrix and
physical maxima and both exceed its threshold, record `matrix_residual`
first, without changing the existing gate or message. Boundary/branch checks
performed by the approved wrapper may finalize its owned copy of the detail
with the already-specified failure code; they cannot overwrite an earlier
integrator failure. Unknown or unmappable diagnostic failure is
`diagnostics/diagnostic_incomplete` and remains fail-closed.

#### Approved diagnostic representation corrections

OI-5 approves exactly these two companion/preregistration corrections, in
addition to the interface and permissions above. They are effective under the
hash-bound OI-5 approval; the schema field sets and artifact names remain
unchanged.

1. Replace the companion's universal native-unit, positive-normalizer sentence
   with: FE/IMEX/BE **physical** residual tuples and IMEX **matrix** tuples
   retain the existing dimensionless `residual_max`, `scale_max`, and
   `normalized` values from their row-scaled helpers. No row scaling or
   reduction order changes. BE matrix tuples retain native-unit `residual`
   and `scale` from the row maximizing the existing row-normalized residual
   separately for each field; exact ties keep the smallest level index.
   A normalizer is finite and nonnegative. Preserve the existing expression
   `scale > 0 ? residual/scale : residual` exactly, including its zero-scale
   branch; never introduce a positive floor. Each tuple's units follow this
   fixed method/residual-kind mapping, with no new schema key. Only like-unit
   terms enter each original calculation. The per-field maxima must reproduce
   the existing aggregate normalized residual on a complete call.
2. Insert `raw_nonfinite` immediately after `raw_material_negative` in the
   preregistration's canonical `diagnostic_failure_reasons` array and permit
   `raw_state` to map to either reason. This represents the existing FE
   nonfinite-update returns at `ndrop_integrator.cpp:669,719`; it adds no gate,
   threshold, or failure behavior.

The first correction resolves a representation conflict, not a request to
redefine physical residuals. Row-dependent scaling means the two maxima need
not share a row or a common dimensional scale; multiplying both by an invented
scale would not recover native residual evidence. Requiring additional
unscaled numerators would add arithmetic, overflow cases, and a second
normalization unrelated to the existing gate. The approved minimum is to
retain exactly the values that established that gate. The OI-5 approval
expressly accepts this replacement of the native-unit/positive-normalizer
requirement. No gate, normalized value,
reference tolerance, or numerical equation is relaxed.

#### Reference failures: wrapper solution, no source extension

`ReferenceResult::ok`, `message`, state, and counters already support the
required result distinction when combined with wrapper knowledge of which
call was attempted. The source does not expose RK stage internals, but the
approved schema requires only setup/advance/pair-comparison stage and the
listed failure reasons. Adding a structured status to `ReferenceResult` would
duplicate the following finite source-bound mapping and require an unnecessary
edit to every reference return. Recommend the wrapper solution in the already
allowed `ndrop_production_replay.cpp`; leave `ReferenceResult` and
`ndrop_reference.cpp` unchanged.

For primary or tight invocation `X`, match exact messages as follows:

| Existing return/message | Serialized stage | Reason/outcome |
|---|---|---|
| `ok=true`, message `ok`, complete finite state | None | Mark only X completed; compare only after both calls complete |
| `reference input validation failed:` followed by a space and the exact message from the existing input validator, or `reference options must be finite and positive` | `X_setup` | `incomplete_reference_pair` |
| `reference operator contains a nonfinite coefficient` | `X_setup` | `nonfinite_stage` |
| `reference exceeded maximum attempted steps` | `X_advance` | `step_limit_exhausted` |
| `reference step size underflowed` | `X_advance` | `step_underflow` |
| Exact `reference stage N is nonfinite`, for N=1..7; `reference embedded solution is nonfinite`; `reference error norm is nonfinite`; `reference final state is nonfinite`; `reference zero-cloud postprocess overflowed` | `X_advance` | `nonfinite_stage` |

Primary failure prevents the tight call. Tight failure retains
`primary_completed=true`, `tight_completed=false`. Preserve the original
failure message in `gate_results.reference.message`; `failure_stage` expands
X to `primary` or `tight`. An input-validation suffix must reproduce a
read-only call to the same existing validator with the same inputs; no loose
substring or guessed stage classification is accepted. Unknown messages,
exceptions, missing returns, or inconsistent success/state/counter combinations
fail closed as `incomplete_reference_pair` at the active X advance (or X setup
if the invocation never began). They are explicitly uncertified interface
failures, never fabricated nonfinite/underflow diagnoses or completed results.
They require implementation review before admitting that source identity.
If the process cannot retain a result, the existing execution-completion gate
rejects the run; a wrapper cannot promise recovery from process termination.

A nonfinite comparison after two completed references is
`pair_comparison/comparison_nonfinite`. A finite comparison above `2e-11`
remains a completed comparison with `passed=false`. All existing nullability,
candidate-first precedence, sentinel exclusion, and `1e-12`/`1e-13` tolerances
remain exact. Never infer completed steps or function evaluations beyond the
existing counters. The wrapper source and reference source remain bound by
the existing source/tree/executable provenance; a changed return vocabulary
requires classifier review and its offline tests before renewed preregistration.

#### Approved file permissions and timing boundary

Only these four additions to the existing offline allowlist are necessary:

| Additional path under `docs/analysis/evidence/ndrop_integrator/` | Permission and necessity |
|---|---|
| `ndrop_integrator.hpp` | Append the exact diagnostic types, optional result detail, and opt-in control; preserve existing signatures and members |
| `ndrop_integrator.cpp` | Expose FE per-field/per-substep residual intermediates and completion/validity at existing producers/returns; do not replace FE or its helpers' calculations |
| `ndrop_imex.cpp` | Expose existing RHS residual intermediates, retain missing per-field raw minima, and record completion/validity; preserve the public EKAT Thomas calls and existing combined phase boundaries |
| `ndrop_backward_euler.cpp` | Expose existing block-component/scalar and physical residual intermediates and completion/validity; preserve families, pivots, factorization reuse, solves, branch decisions, and separate phase boundaries |

`ndrop_reference.cpp` is explicitly excluded: the complete wrapper solution
above needs no change there. Existing test files, `ndrop_phase_timer.hpp`,
production files, build dependencies, and gitlinks need no additional
permission. Tests and serialization changes belong in the already-approved
replay test/schema/wrapper files; build integration and documentation use the
already-approved CMake and README paths. Private helper signatures in the three
candidate implementation files may accept only the optional output needed for
this detail. Preserve their existing callers/probes, return values, arithmetic,
expression and reduction order, and failure/short-circuit behavior. No solver
or numerical advance is duplicated in the wrapper.

Default `collect_replay_diagnostics=false` preserves all existing callers.
Primary timing sets both diagnostic flags false and `phase_timings=nullptr`.
It performs no detail construction, array initialization, residual evaluation,
copy, allocation, clock call, or serialization introduced by OI-5; the optional
result member remains disengaged. Correctness collection sets both flags true
outside primary timing. Mandatory Tier-2 overhead may measure that separate
correctness-instrumented arm under its existing diagnostic-only rules. No
Tier-3 phase instrumentation is authorized here. Numerical advancement,
positivity checks, FE clipping, finite checks, and zero-cloud state transitions
remain on their existing paths in every mode. Existing diagnostics-only gates
retain their existing control semantics; an uninstrumented success never
substitutes for the separate correctness gate.

#### Future offline acceptance tests and approval boundary

Specify, but do not execute in this revision, these tests in the already
allowed `test_ndrop_production_replay.cpp`, `test_ndrop_replay_schema.cpp`,
and `test_ndrop_replay_artifacts.py`:

- Verify all 51/25 mappings, units, masks, implicit FE-only nulls, completed
  zeros, raw minima before reset, exact ledger reuse, and complete/partial
  companion round trips. Synthetic early-return fixtures must retain completed
  fields while leaving unreached/failed arithmetic null, including setup,
  partial FE substeps, projection, and final ledger failures.
- Verify FE streaming maxima across multiple substeps, distinct worst
  substeps by field, exact tie selection, implicit RHS/component mapping,
  zero-scale branches, and aggregate reproduction. Use repository/synthetic
  fixtures at the existing `1e-12`, `1`, `1e12` scales; verify the existing
  residual formulas independently without using the new collector as oracle.
- Compare pre-change source behavior with the amended build using identical
  local fixtures/compiler settings; require bitwise-identical numerical
  states and previously populated diagnostics in the same collection mode,
  and identical messages, gates, step counts, parity, limiter/clipping/negative counts,
  factorization/solve counts, and branch decisions. Compare new detail on/off
  with legacy detailed collection on (the formerly unpopulated IMEX per-field
  minima and new detail are the only added outputs), then detailed on/off on
  admitted valid fixtures. Diagnostics-only failure fixtures retain their existing difference
  from the lean path; they cannot acquire timing eligibility.
- Check that default/disabled collection leaves detail disengaged and enters
  no added collection operation, using source-path inspection and offline
  allocation/coverage checks, not timing benchmarks. Required numerical work
  must still execute. No production state or retained result is a test input.
- Exhaustively test every reference message mapping, both invocation roles,
  unknown/truncated messages, exceptions, inconsistent returns, candidate-first
  and sentinel not-attempted states, incomplete pairs, and nonfinite/above-gate
  comparisons. Reuse the reachable failure cases in
  `test_cross_scheme.cpp:934-989`; test otherwise unreachable classifications
  with fabricated `ReferenceResult` values in the wrapper tests, clearly
  distinguished from actual reference executions. Preserve all counters.
- Run the unchanged existing numerical and benchmark-support **unit tests**
  plus the new offline suites in debug ASan/UBSan and optimized builds after
  approval. Generated outputs use a new task-specific local scratch directory.
  These tests do not authorize a timed benchmark or production-state replay.

The numerical diagnostic interface is implementation-blocking; no assumed
completion, final-state reconstruction, or partial substitute resolves it.
The user's hash-bound OI-5 approval is recorded as `ADDENDUM_APPROVED` to
`ADDENDUM_APPROVED`, an ordinary revision at the existing state. It authorizes
only the four narrowly
scoped offline file additions, exact diagnostic exposure, two explicit
representation corrections, wrapper reference classifier, and associated
local offline implementation/testing above. Every other addendum requirement
and later approval remains unchanged. This bookkeeping records that authority
without expanding it. In particular it does not advance `CAPTURE_WRITER_APPROVED` or
authorize capture, qualification/census execution, production-state replay,
timed benchmarks, performance analysis, retained-result access, remote/HPC
activity, publication, staging, or committing.

### Backend completion and ownership mapping

Under the approved OI-1 contract, the completion mapping to be frozen by a
later `CAPTURE_WRITER_APPROVED` record is:

| Boundary | Serial/OpenMP execution | CUDA/HIP execution | Failure rule |
|---|---|---|---|
| P1a pre-update owned-state preservation | Common old state and raw inputs copy in the same team after the barrier at locked-source line 1497 and complete before `update_from_explmix` may mutate source storage | Same device-team rule | Alias retention, an early mutation, or absent completion is `UNAVAILABLE`. |
| P1b limiter observation | Unique per-level/mode stores at the existing assignments, with qualification-only pre/cap/post trace and census direct flags; no new barrier | Same | A missing trace in qualification, equality labeled active, or any captured value feeding production is `UNAVAILABLE`. |
| P1c selector preservation | Existing barrier at line 1151 completes array producers; observed `dtmin`, `nsubmix`, post rates, overlaps, and directional rates are copied before line 1169 without a new reduction or barrier | Same | A fabricated maximum/level, changed reduction order, or missing completion is `UNAVAILABLE`. |
| Post-reset preservation | Final values and parity copy to the preallocated capture slot only after the existing barrier at line 1508, which follows all substeps and the zero-cloud resuspension/reset | Same | No production alias may escape the team and no pre-reset state may be labeled final. |
| Caller Phase 2 | Existing post-call global `Kokkos::fence()`, then capture-owned host copy and copy-instance completion | Existing post-call global `Kokkos::fence()` proving all producers complete, then deep copy on the recorded copy instance and its fence before host read | Unproven queue or copy completion is `phase2_producer_incomplete`. |
| Serialization | Host only after both capture-owned copies validate | Host only after both capture-owned copies validate | Device or timed-region serialization is prohibited. |

The qualification manifest selects exactly one supported backend and records
which row mechanism was executed. Portability of the interface does not
qualify an untested backend. An existing equivalent completion guarantee may
replace a listed fence only when its source location and semantics are bound in
the approved qualification manifest and tested for that backend.

### Changed-file allowlists

`ADDENDUM_APPROVED` permits local offline implementation and tests only in this
exact allowlist:

```text
docs/analysis/evidence/ndrop_integrator/CMakeLists.txt
docs/analysis/evidence/ndrop_integrator/README.md
docs/analysis/evidence/ndrop_integrator/ndrop_replay_schema.hpp
docs/analysis/evidence/ndrop_integrator/ndrop_replay_schema.cpp
docs/analysis/evidence/ndrop_integrator/ndrop_production_replay.hpp
docs/analysis/evidence/ndrop_integrator/ndrop_production_replay.cpp
docs/analysis/evidence/ndrop_integrator/benchmark_ndrop_production_replay.cpp
docs/analysis/evidence/ndrop_integrator/ndrop_replay_artifacts.py
docs/analysis/evidence/ndrop_integrator/analyze_production_replay.py
docs/analysis/evidence/ndrop_integrator/test_ndrop_replay_schema.cpp
docs/analysis/evidence/ndrop_integrator/test_ndrop_production_replay.cpp
docs/analysis/evidence/ndrop_integrator/test_ndrop_replay_artifacts.py
docs/analysis/evidence/ndrop_integrator/test_analyze_production_replay.py
```

It does not permit edits to existing integrator equations. Implementation
edits are limited to the diagnostic-only OI-5 extension below.
If an implementation requires a file outside this list, work stops for an
approved addendum revision.

OI-5's four diagnostic-only additions, listed and justified in its approved
subsection, are part of the effective offline allowlist under its hash-bound
approval. `ndrop_reference.cpp` remains excluded.

A later `CAPTURE_WRITER_APPROVED` record may authorize only this separate exact
production/instrumentation allowlist:

```text
MAM4xx descendant: src/mam4xx/ndrop.hpp
MAM4xx descendant: src/mam4xx/ndrop_capture.hpp
MAM4xx descendant: src/tests/CMakeLists.txt
MAM4xx descendant: src/tests/mam4_ndrop_capture_unit_tests.cpp
E3SM: components/eamxx/src/physics/mam/CMakeLists.txt
E3SM: components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp
E3SM: components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.hpp
E3SM: components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp
E3SM: components/eamxx/src/physics/mam/eamxx_mam_ndrop_capture.hpp
E3SM: components/eamxx/src/physics/mam/eamxx_mam_ndrop_capture.cpp
```

The MAM4xx instrumentation source is a clean descendant of the locked target
or an exact clean application of the canonical patch in a separately retained
dependency source tree. The E3SM gitlink continues to name its locked object
ID; the build manifest separately binds the exact retained MAM4xx
instrumentation source used by the executable. No gitlink edit, dependency
substitution other than that exactly bound MAM4xx instrumentation source,
dirty tree, generated source change, or path outside the list is permitted.
The complete diff hash covers additions and deletions in both lists.

### Required offline tests before any capture-writer approval

The addendum implementation must add focused fail-closed tests for:

- canonical JSON, path, schema, unknown-field, duplicate-key, and SHA-256
  validation; canonical `real_hex` round trips for normals, subnormals,
  positive and negative zero; rejection of equivalent noncanonical spellings;
  IEEE total-order fixtures; unique-timestep and stratum-order fixtures; and
  stable canonical hashes;
- byte-exact `ndrop_replay_state_v1` fixtures for binary32 and binary64 at both
  L128 and labeled L72, proving the nine-byte `NDRPSTAT1` magic, offsets
  0--16, 17-byte header, exact field order, exact total sizes, and rejection
  of every corrupted header, width, endian, reserved-byte, extent, or length;
- exact canonical complete-diff framing, endpoint metadata, section lengths,
  concatenation order, and combined-artifact hash reproduction;
- exact current-`cldn` retention, zero-cloud transition inputs, and checksum;
- Phase 1/Phase 2 identity, ownership, completion, parity, and mismatch
  failures, including no acceptance of an alias-backed or completion-free
  record;
- all role/status applicability and precedence combinations, including
  reference not-attempted, completed, failed-before-a-pair, nonfinite-stage,
  step-underflow, step-limit, and incomplete-pair outcomes without fabricated
  metrics; all `UNAVAILABLE` reason codes; and blocking of an unaccounted or
  unstratifiable census record;
- immutable qualification/census roles, qualification-only consumption, and
  rejection by every census/replay/timing/error/method-selection consumer;
- exact active-level mapping, no-flux endpoints, production rates, limiter
  activity, observed `dtmin`, reconstructed selector-candidate/minimizer tie
  sets, controller-only FE contribution, `chi_mix`, `chi_act`, `chi_FE`, and
  returned `nsubmix`; limiter flags are tested as direct strict comparisons
  through the qualification trace, including equality as inactive, and are
  never inferred from equal post-limiter values;
- exact locked aliases, rank-3 shapes, wrapper rank-4-to-rank-3 subviews,
  parameter order, memory space, `DeviceSlot` metadata offsets/lifecycle,
  slot-to-column mapping, unavailable and capacity-exhausted slots, and
  checked `int`/signed-64-bit `nsubmix` conversions in a future compile-only
  compatibility target at the approved offline/capture-writer stages;
- full-state initial/final equality against the existing production-count FE
  adapter without changing FE, IMEX, BE, or reference arithmetic;
- exact `ndrop_replay_diagnostics_v1` ordering, units, shapes, substep
  aggregation, method applicability, status nullability, ledger balance,
  raw-minimum classifications, FE clipping deficits, matrix residuals, and
  independently assembled physical residuals;
- unique census identities, contiguous capture sequences, complete accounting,
  no adaptive thinning, exact census cross-tabulations, L128-production versus
  L72-compatibility exclusion, and 64/65/512/513 censoring boundaries;
- mixed replayable/censored records in one stratum, unavailable records with
  known strata and weights, unavailable records with unknown strata or
  weights, exact `R`/`C`/`U` numerators and denominators, conditional-estimand
  coverage, explicit excluded weight, both one-percent censoring tests, and no
  redistribution to selected records;
- deterministic pilot and next-hash expansion, one-case strata, missing-pair
  indeterminacy, stress-sentinel schema/mutation validation, unique replay-case
  keys, zero weighting, all preregistered trigger ceilings, and complete
  expansion-decision records when no additional eligible record exists;
- SplitMix64 stream fixtures, exact 10,000-draw hierarchical paired bootstrap,
  even median, percentile endpoints, weighted quantiles, two weight bases, and
  casewise numerical tie rules;
- required warm-up and measured counts, complete interleaved matched blocks,
  stopping calculations and reasons, minimum eligible-pair handling, and the
  distinction between missing and insufficient pairs;
- strict separation of primary uninstrumented timing, mandatory Tier-2-only
  overhead evidence, and optional separately approved Tier-3 phases, including
  `method_inapplicable` versus `deferred_unmeasured`;
- qualification production-sequence identity, contiguous process-subcycle
  rows, retained raw-log hashes, repeated-call storage reuse, exact pre/post
  per-field equality, downstream step mapping, rejection of final-step output
  as evidence for hidden inner calls, diagnostics-level exclusion from overhead
  timing, and base/disabled/enabled three-arm evidence without disabled-path
  capture instrumentation;
- rejection of historical Haero provenance as a current dependency and proof
  of the exact current E3SM/MAM4xx/EKAT/Kokkos dependency set;
- rejection of partial replay outputs, mixed execution bindings, and any
  analysis input lacking a valid replay-execution completion manifest;
- rejection of an altered census, adapter, candidate, reference, timing, or
  analysis source/tree/executable after preregistration; and
- rejection of an unapproved, mismatched, over-capacity, or wrong-stage replay
  run/resource manifest and of an expansion lacking renewed bounded execution
  approval.

These remain proposed tests, not executed evidence. The hash-bound addendum
approval authorizes their local offline implementation and execution only.
Passing them would not authorize production instrumentation or any capture,
replay, benchmark, performance-analysis, remote, or HPC run.

The OI-5 interface-completeness and unchanged-behavior tests are separately
specified above and authorized by its hash-bound approval. Existing
test expectations and historical executed evidence are not revised.

### Focused blocker-resolution decisions and open issues

The contract-only approval of document SHA-256
`a829ffc01d09a9960e9bc5df04a15199a39833fd6942cb696a8c942bb1292c38`
resolved OI-1 and OI-2 without advancing the approval state. The subsequent
approval of the exact addendum content at SHA-256
`ca30c32c85451ff2fa14666119cfa7da8f665c9d5f2eb662c6a40574b513f470`
advanced the state to `ADDENDUM_APPROVED` and authorizes only the enumerated
local offline implementation and tests. The four focused blockers now have the
following disposition:

| Blocker | Supporting local evidence | Approved resolution or disposition | Exact approval still needed |
|---|---|---|---|
| Capture checkpoints and limiter evidence | Locked MAM4xx `ndrop.hpp:1101-1159`, `1254-1275`, and `1467-1508` show the actual limiter, selector, reset, and lifetime boundaries; the source stores `dtmin` and `nsubmix`, but no maximum or limiting level. | **Approved as OI-1.** Phase 1 is P1a/P1b/P1c plus post-reset capture-owned preservation; observed `dtmin`/`nsubmix` and the exact offline candidate/minimizer tie set replace observed maximum/scalar-level claims. | Offline interface and test implementation is authorized by `ADDENDUM_APPROVED`; the production capture-writer/source edit still requires `CAPTURE_WRITER_APPROVED`. |
| Mixed-tail population estimand | Replay observations cannot exist for `>512` or unavailable records, while the approved census retains their weights and both 1% censoring tests. | **Approved as OI-2.** The estimand is conditional on the replay-eligible population with exact `R`/`C`/`U` accounting, denominators `R_b`, equal selected-case weight within stratum, explicit excluded coverage, and fail-closed missing-pair rules. | Offline census/selection/analysis implementation and testing is authorized; census capture, preregistration, replay execution, and evidence gates remain separate. |
| Qualification across a production sequence | Existing `AtmosphereProcess::run` hashes every process subcycle at level 2; the hash implementation fences and records every input/output/internal field; one driver object executes repeated steps; MAMAci contains one `dropmixnuc` wrapper call per `run_impl`. | Use the existing all-arm per-subcycle hash path for sequence equivalence, enabled-only direct capture trace for hidden decisions, and a separate level-0 matched timing execution for overhead. No new source interface is proposed. | The evidence schema and proof chain are included in `ADDENDUM_APPROVED`; actual qualification remains blocked by `CAPTURE_WRITER_APPROVED`, bounded qualification approval, and observational-equivalence approval. |
| Approved-baseline provenance | The local unreachable Git blob `35e6a2ef497613544bc6749591eae7068a98f74c` is 197,696 bytes and hashes exactly to the approved SHA-256 `d79a14253ca3acb5108d54f37b3ddcc9aaa2a9d9c4ecef391721c27868df2138`. | Treat those verified bytes as the authoritative approved snapshot; retain the approval statement and compare later work bytewise against that blob. Do not invent a commit association. | None. A prospective replacement baseline is unnecessary. |

The baseline search was confined to existing local evidence: all 15 reachable
revisions of this path, 60 reflog revisions, the empty stash set, exact hash
mentions in local documentation/logs, same-named local documents, 10
unreachable commits, and 42 unreachable Git blobs. The matching blob was not
found in a reachable or unreachable commit, so its commit/time attribution is unknown; that does not
alter the byte-exact match to the hash named in the historical approval.

#### Open Issues Log

| ID | Class | Human input or evidence required | Owner / next action | Unblock condition |
|---|---|---|---|---|
| OI-1 | Resolved by contract-only approval | None. The exact P1a/P1b/P1c, direct-limiter, observed-`dtmin`, and minimizer-tie-set contract is approved at snapshot `a829ffc01d09a9960e9bc5df04a15199a39833fd6942cb696a8c942bb1292c38`. | Offline implementer for the approved interface/tests; capture-writer owner for any production-source integration. | Contract and addendum gates are resolved; the source edit remains gated by `CAPTURE_WRITER_APPROVED`. |
| OI-2 | Resolved by contract-only approval | None. The exact conditional replay-eligible estimand and `R`/`C`/`U` denominators are approved at the same snapshot. | Offline implementer for the approved census/selection/analysis software and tests. | Contract and addendum gates are resolved; production evidence remains gated by census, preregistration, and replay approvals. |
| OI-3 | Deferred to qualification gate | Freeze the actual production-sequence population, process DAG, subcycling, per-field hash inventory, supported backend/compiler/layout, commands, and overhead bound. | Qualification owner at `QUALIFICATION_CAPTURE_APPROVED`. | Bounded three-arm manifest and command approvals name every value; this is not an addendum-design ambiguity. |
| OI-4 | Deferred to capture-writer gate | Prove the enabled checkpoint operations compile on the approved backend, use only the existing completion points, and leave the disabled path operation-free. | Capture-writer owner; prepare the next separately reviewed transition. | Compile-only compatibility and focused local tests pass under a hash-bound `CAPTURE_WRITER_APPROVED` record. |
| OI-5 | Resolved by same-state approval; implementation/test evidence pending | Approved at document SHA-256 `2300c9ef1448fb161042038e0075d8646cdfc733624ba2c0f7160e21b2f2fbe0`: four diagnostic-only file additions, exact exposure, two representation corrections, reference wrapper classifier, and associated local offline implementation/testing. | Offline implementer proceeds only within those permissions. | Contract permission is resolved; local offline conformance tests must still establish implementation completeness and unchanged behavior. All production and execution gates remain. |

Assumptions are limited to the source identities already bound by this
addendum: MAM4xx
`0d076d8ebb004a18f05f51f97b08edd9def6826f` and locally inspected E3SM
`e2c578e40f6602a108eff577091de86a453624ab`. No qualification command was run,
so existing diagnostics are source evidence of capability, not executed
equivalence evidence. If a later approved configuration fails to emit the
frozen field inventory or map every MAMAci subcycle, qualification fails closed
and returns for document review; final-step equality is not a fallback.

### Addendum authorization boundary and remaining approvals

This addendum is **`ADDENDUM_APPROVED`**. The user's explicit approval binds the
exact addendum content at SHA-256
`ca30c32c85451ff2fa14666119cfa7da8f665c9d5f2eb662c6a40574b513f470`
and authorizes only local offline implementation and testing of the validator,
locked-source replay adapter, census and selection code, replay benchmark,
analysis software, schemas, manifests, and associated offline tests in the
approved allowlist. This document-only status revision performs none of that
implementation. The approval grants no capture-writer/source edit,
production-state capture, replay or benchmark execution, performance analysis,
method-selection conclusion, retained-result access, remote or HPC activity,
publication, staging, commit, or push authority. The OI-5 diagnostic
interface revision is approved at SHA-256
`2300c9ef1448fb161042038e0075d8646cdfc733624ba2c0f7160e21b2f2fbe0`
for the same-state transition and only its enumerated offline scope. The next state-advancing
transition remains `ADDENDUM_APPROVED` to `CAPTURE_WRITER_APPROVED`, after the
required offline evidence; OI-5 does not grant that authority.

The remaining transitions and gate are ordered and fail closed:

| Transition or gate | Required immutable evidence | Exact authority and exclusions |
|---|---|---|
| `ADDENDUM_APPROVED` to `CAPTURE_WRITER_APPROVED` | Separate approval record binding the approved addendum hash, exact source signatures and checkpoints, production changed-file allowlist, and canonical patch/diff identity | Only the named observational instrumentation/source edit and its approved local tests. It grants no capture run. |
| `CAPTURE_WRITER_APPROVED` to `QUALIFICATION_CAPTURE_APPROVED` | Instrumentation-only provenance plus an immutable bounded three-arm qualification population/resource manifest and exact command/facility approvals | Only the named qualification capture. Its authority is non-transferable and its artifacts cannot enter census, replay, timing, error, or method-selection evidence. |
| `QUALIFICATION_CAPTURE_APPROVED` to `OBSERVATIONAL_EQUIVALENCE_APPROVED` | Approval record binding the successful three-arm manifest/evidence hashes, exact tested implementation and executable/build, backend, compiler configuration, layout, seeds, capture settings, synchronization proof, and approved overhead bound | Qualifies only those exact identities and settings for a separately approved census. It does not authorize census capture. |
| `OBSERVATIONAL_EQUIVALENCE_APPROVED` to `CENSUS_CAPTURE_APPROVED` | Census-specific frozen population/resource manifest, command-specific approvals, and verification that its producer exactly matches the qualified capture-enabled implementation, executable/build, backend, compiler configuration, layout, seeds, and equivalence-relevant capture settings | Only that bounded census capture. Qualification authority is non-transferable. Qualification and census may differ only in population and immutable role; any other bound change requires renewed three-arm qualification. Independently valid provenance is not a substitute. |
| `CENSUS_CAPTURE_APPROVED` to `CENSUS_SEALED` | Complete census validation and accounting with every expected record classified, all completion/hash checks valid, zero unaccounted records, and the exact final census-manifest hash | This is a data-validation and sealing gate, not an approval state. Passing it permits preparation of census-bound preregistration only; candidate execution and inspection remain prohibited. |
| `CENSUS_SEALED` to `PREREGISTRATION_APPROVED` | Approval record binding the sealed census hash, population estimand, selection/expansion/statistical/tie rules and values, seeds, and exact replay-harness, FE, IMEX, BE, reference, timing-driver, and analysis source/tree/executable hashes | Permits materializing the frozen selection and preparing a bounded replay run/resource manifest. Any bound change invalidates preregistration and requires renewed approval. |
| `PREREGISTRATION_APPROVED` to `REPLAY_EXECUTION_APPROVED` | Exact bounded pilot or expansion run/resource manifest and command-specific approvals | Only the named candidate execution and result inspection. Every expansion requires a renewed bounded approval. It grants no GPU/whole-ACI execution, publication, physics/API change, or production recommendation. |

GPU or whole-ACI work, publication, any physics or production-API change, and
any production-method recommendation require their separately named later
approvals and prerequisites under the amendment. Access to retained Task 5
results remains separately prohibited; neither amendment approval nor
`ADDENDUM_APPROVED` grants it.

### Current limitations

- OI-5 is approved, with implementation/test evidence pending: the offline audit confirmed missing transient residual
  exposure, IMEX per-field minima, and partial-completion metadata. It also
  identified the residual-unit/zero-normalizer and raw-nonfinite representation
  gaps. Its approved four-file extension addresses these without a numerical
  change; reference failures can be classified through the existing interface.
  No implementation or test run has resolved this blocker yet.
- The exact bytes of approved amendment snapshot
  `d79a14253ca3acb5108d54f37b3ddcc9aaa2a9d9c4ecef391721c27868df2138`
  were recovered as local Git blob
  `35e6a2ef497613544bc6749591eae7068a98f74c` and verified byte-for-byte. The blob
  has no locally identified containing commit, so no commit/time attribution is
  claimed; baseline-to-current content comparison is nevertheless available.
- Authoritative executed sensitivity-test manifests are unavailable, so the
  proposed ne256 intervals are not executed evidence.
- Physics owners have not supplied a production accuracy budget, so only
  Pareto/status conclusions are allowed.
- The source census has not established extreme-tail materiality.
- Production capture and whole-ACI GPU work require a supported machine,
  resources, policy review, and explicit approval.
- The rank, memory-space, and locked production/wrapper parameter facts have
  been source-anchored, but the original amendment's single pre-update Phase 1 boundary
  cannot observe all required directional-rate, overlap, limiter-activity,
  and selector values. The exact P1a/P1b/P1c and post-reset resolution is now
  contract-approved through OI-1, and its offline interface and tests are
  authorized by `ADDENDUM_APPROVED` but remain unimplemented. The production
  capture-writer edit remains separately unapproved and unimplemented.
- Existing per-process diagnostics can expose every MAMAci subcycle in a
  repeated production sequence and the existing driver retains state across
  steps. This resolves the interface-design question, not the evidence gate:
  no three-arm sequence has been run, and a missing field or subcycle in the
  later frozen configuration fails qualification rather than falling back to
  final-step output.
- The approved amendment does not uniquely determine whether production
  timing/error aggregates estimate the replay-eligible population or the
  complete pre-cap census in mixed 65--512/>512 strata. OI-2 now resolves that
  ambiguity in favor of the exact conditional replay-eligible estimand. Its
  offline aggregation software and tests may now be implemented, but no
  production-data aggregate may be executed or reported before the census,
  preregistration, replay-execution, and evidence gates pass.
- The qualification backend/compiler/layout, exact qualification population
  and commands, overhead bound, census population manifest, and
  preregistration seeds, trigger thresholds, weights, and tie-tolerance values
  remain fail-closed inputs at their respective later approval states.
- Existing retained Task 5 results are outside this amendment's authorized
  inputs and have not been accessed or modified.

### Consolidated addendum-review disposition

This table records the B1--B15 review disposition without treating a proposed
test or unresolved proof obligation as executed evidence:

| Finding | Disposition | Document resolution or remaining blocker |
|---|---|---|
| B1 | Corrected | The candidate-state header is byte-exact and consistently 17 bytes; all four L128/L72 and binary32/binary64 future fixtures are specified. |
| B2 | Artifact and OI-5 interface amendment approved; implementation/test evidence pending | `ndrop_replay_diagnostics_v1` requires ordered ledgers, minima/counts/deficits, and residuals outside primary timing. OI-5's hash-bound approval authorizes the minimum diagnostic-only interface/permission corrections and two explicit representation corrections; conformance tests remain required. |
| B3 | Corrected | The closed applicability/precedence matrix and reference outcome object represent not-attempted, completed, and pre-pair failures without fabricated metrics. |
| B4 | Corrected | Closed nested grammars, canonical `real_hex`, IEEE ordering, and byte-framed complete-diff hashing are specified; run-specific values remain at their already assigned approval states. |
| B5 | Corrected | Versioned census cross-tabulations and per-stratum expansion-decision records are retained and hash-bound without circular dependencies. |
| B6 | Corrected | Production L128 and compatibility L72 admission, selection, weighting, and aggregation are explicitly separated. |
| B7 | Corrected | Applicable warm-up, measurement, stopping, complete-block, and eligible-pair evidence is explicit without reviving superseded synthetic-screening rules. |
| B8 | Contract and addendum resolved; production integration blocked | Actual ranks, aliases, dimensions, memory spaces, parameter lists, `DeviceSlot`, checked integer conversions, census admission, and future compile-only checking are specified. OI-1 approves the rank-correct P1a/P1b/P1c contract; offline implementation/testing is authorized, while production integration remains gated by `CAPTURE_WRITER_APPROVED`. |
| B9 | Contract and addendum resolved; production integration blocked | Limiter flags are direct strict `pre_rate > cap` observations, equality is inactive, and qualification retains exact pre/cap/post evidence. Census flags are never reconstructed from post-limiter rates. Offline implementation/testing is authorized; capture-writer integration and qualification remain separately gated. |
| B10 | Corrected; execution pending later gate | Existing EAMxx level-2 per-process diagnostics, repeated-driver execution, enabled-only trace, and separate level-0 overhead timing form the exact sequence proof without a new source interface or disabled-path capture instrumentation. No qualification evidence has been executed. |
| B11 | Corrected | Mandatory Tier-2 overhead evidence and optional, separately gated Tier-3 phase evidence have disjoint schemas and applicability states. |
| B12 | Corrected | Current E3SM/MAM4xx/EKAT/Kokkos provenance is segregated from historical Haero evidence; Haero cannot enter the current identity. |
| B13 | Contract and addendum resolved; production evidence blocked | Exact `R`/`C`/`U` counts, weights, coverage, one-percent tests, conditional numerators/denominators, mixed-tail treatment, and missing-pair indeterminacy are specified and approved through OI-2. Offline implementation/testing is authorized; census, preregistration, replay execution, and evidence gates remain. |
| B14 | Contract and addendum resolved; production integration blocked | The approved OI-1 source-anchored P1a/P1b/P1c and post-reset preservation map resolves producer/lifetime ambiguity without changing arithmetic or adding synchronization. Offline implementation/testing is authorized; production capture integration remains gated by `CAPTURE_WRITER_APPROVED`. |
| B15 | Corrected | Reports must state that a one-case-stratum timing bootstrap cannot measure unobserved between-case variation or establish stratum-wide uniformity. |

The exact approved baseline is verified, B10 has a source-supported resolution,
and the contract-only approval resolves B8/B9/B14 and B13 at the contract
level. The exact addendum content at SHA-256
`ca30c32c85451ff2fa14666119cfa7da8f665c9d5f2eb662c6a40574b513f470`
is now `ADDENDUM_APPROVED`; local offline implementation and testing in its
allowlist is authorized. No implementation occurs in this status-recording
revision, and production-source or execution work remains blocked by every
subsequent state-machine approval and evidence gate.

## Known Limitations and Open Decisions

- For the locked replay target, admissible records require
  `0 <= top_lev < pver` and use the locked no-flux physical top and bottom
  boundaries. A new production range guard, configuration rejection, public
  API contract, or boundary-physics change remains a separate open design
  decision requiring explicit approval.
- The earlier investigation interpreted `enable_aero_vertical_mix=false` as
  disabling turbulent transport while retaining conservative activation and
  cloud/interstitial phase transfers. Current MAM4xx instead ignores the
  compatibility argument and always performs mixing. EAMxx namelist text
  should be updated so users are not offered a setting that has no effect.
- GPU bounds behavior may differ from a bounds-checked host build. Host
  sanitizer evidence should be paired with deterministic invariant tests that
  are portable to GPU builds.
- Validation at the older evidence revision initialized both alternating work
  arrays and therefore cannot by itself establish the safety of that historical
  partial-disable path.
- Tasks 2 through 4 now exercise repaired current-MAM4xx boundaries at
  `top_lev=0,1,6`, verify paired aerosol ledger decompositions, gate FE, IMEX,
  and host Backward Euler independently, and pass the complete common
  cross-method parameter matrix. Local CPU timing, matched-error analysis, and
  whole-model boundary evidence remain pending.
- The `Kokkos::Min` reducer overwrites, rather than incorporates, the
  preassigned `dtmin=dtmicro`. This permits `nsubmix=1` under weak mixing and
  differs from the apparent source-level intent; the numerical and performance
  implications require separate follow-up.

## Scientific References

- Abdul-Razzak, H., and S. J. Ghan (1998), A parameterization of aerosol
  activation: 1. Single aerosol type, _Journal of Geophysical Research_,
  103(D6), 6123-6131.
- Abdul-Razzak, H., and S. J. Ghan (2000), A parameterization of aerosol
  activation: 2. Multiple aerosol types, _Journal of Geophysical Research_,
  105(D5), 6837-6844.
- Ghan, S. J., C. C. Chuang, and J. E. Penner (1993), A parameterization of
  cloud droplet nucleation, _Atmospheric Research_, 30, 198-221.
- Alexander, R. (1977), [Diagonally implicit Runge-Kutta methods for stiff
  O.D.E.'s](https://epubs.siam.org/doi/10.1137/0714068), _SIAM Journal on
  Numerical Analysis_, 14, 1006-1021.
- [SUNDIALS ARKODE user documentation](https://sundials.readthedocs.io/en/latest/arkode/Usage/General.html)
  and [Kokkos installation notes](https://sundials.readthedocs.io/en/latest/sundials/Install_link.html).
- [PETSc `TSARKIMEX` documentation](https://petsc.org/release/manualpages/TS/TSType/).
- Kopecz, S., and A. Meister (2018), [On order conditions for modified
  Patankar-Runge-Kutta schemes](https://www.sciencedirect.com/science/article/pii/S0168927417301861),
  _Applied Numerical Mathematics_, 123, 159-179.

## Revision History

| Date | Author / Agent | Change Summary |
|---|---|---|
| 2026-10-12 | User / Codex | Approved the bounded local timing stage for the 12 frozen physical geometry cases: narrow same-file timing enablement, unchanged statistical/numerical rules, exact timing-binary accuracy requalification with unchanged CSV/selections, ordinary hash-bound approval/execution records, and one local paired timing/analysis pass. No production, capture/replay, directional-rate, dependency/gitlink, retained Task 5, remote/HPC, publication or Git-write authority is added. |
| 2026-10-08 | User / Codex | Approved and recorded the paired physical geometry sensitivity before implementation: four unchanged native controls plus eight conservative 2x/4x cloud-band subdivisions, identical moderate/strong old states, emergent production-controller counts and unchanged matching/reference/convergence gates; separate closed untimed identities and separate future timing approval. No numerical or production implementation, capture/replay, performance-analysis, dependency/gitlink, or retained Task 5 authority is added. |
| 2026-10-08 | User / Codex | Approved and recorded the narrow physical-controller extension before implementation: unchanged physical fixtures 28--31, natural adaptive FE error as the implicit accuracy target, independent artifact family, unchanged numerical/convergence/timing gates, offline tests and untimed selection only; separate timing approval and all production gates remain required. |
| 2026-10-08 | User / Codex | Approved the in-chat redesign and local implementation of a separate four-case physically scaled standalone column study: repository L72/L128 coordinate families, explicitly idealized DYCOMS cloud geometry/modal aerosol/turbulence inputs, natural controller FE baseline versus common 1%/0.1% accuracy-refined configurations, certified references, unchanged numerical/timing gates and separate post-selection timing approval. Recorded the fixture equation/units and assumptions before implementation; no production or directional-rate authority is added. |
| 2026-10-08 | Codex | Prepared a pending four-case controlled moderate/strong synthetic mixing extension at L72/L128 and a common 150 s interval, fixing activation inputs to the existing moderate generator baseline. Specified separate fixture/artifact identities, unchanged controller-selected matched-error and timing rules, future conformance tests, narrow offline file scope, and separate implementation/untimed-selection versus timed-run approvals. Existing approved fixtures, numerical/scientific contracts, evidence, and production gates are unchanged; this document-only preparation grants no implementation or execution authority. |
| 2026-10-08 | User / Codex | Approved the production-controller-matched standalone extension: exact production-selected FE count and its certified reference error define each synthetic case target; IMEX/BE use strict ascending integer accuracy-only selection; prior absolute-target artifacts remain separate; exact controller-stratified synthetic fixtures, new closed artifacts, existing complete timing/bootstrap rules, and a separate post-selection timing approval are required. Production replay, physics, production source/API, dependencies, and production conclusions remain unchanged. |
| 2026-10-07 | User / Codex | Approved a bounded standalone matched-error interval study for frozen synthetic operators: method-specific internal refinement, reference-checked error of the timed state, narrow IMEX/BE interval entry points, final-only cloud projection, complete uninstrumented timing, frozen accuracy-first selection, paired timing blocks, and separate artifacts. Production replay, physics, and Task 5R evidence contracts are unchanged. |
| 2026-09-28 | User / Codex | Recorded the user's OI-5 approval at document SHA-256 `2300c9ef1448fb161042038e0075d8646cdfc733624ba2c0f7160e21b2f2fbe0` for `ADDENDUM_APPROVED` to `ADDENDUM_APPROVED`. Reconciled only status, authority, allowlist cross-references, Open Issues, and limitations. The approval covers only the four diagnostic-only additions, exact diagnostic exposure, two representation corrections, reference wrapper classifier, and associated local offline implementation/testing. OI-5 content is not substantively changed; implementation/test evidence remains separate, and every capture/production/execution/publication gate remains unchanged. |
| 2026-09-24 | Codex | Prepared pending OI-5 after verifying that the pre-edit document matches bookkeeping SHA-256 `8e4e84b45470af8c21fe568b54426817f0a0f1b7e84be0d7bc7e47ccadf3a23a`. Audited the approved diagnostics against current offline source; specified opt-in transient residual/completion exposure, missing IMEX field minima, four justified offline allowlist additions, explicit residual-unit/zero-normalizer and raw-nonfinite representation corrections, and future conformance tests. Chose an exhaustive source-bound reference wrapper classifier instead of adding `ndrop_reference.cpp` permission. Updated only directly affected status/cross-references; OI-5 remains pending separate same-state approval and grants no implementation or execution authority. |
| 2026-09-24 | User / Codex | Recorded the user's hash-bound approval of the exact addendum content at document SHA-256 `ca30c32c85451ff2fa14666119cfa7da8f665c9d5f2eb662c6a40574b513f470` for `AMENDMENT_APPROVED` to `ADDENDUM_APPROVED`. This authorizes only the enumerated local offline implementation and testing of the validator, locked-source replay adapter, census/selection code, replay benchmark, analysis software, schemas, manifests, and associated offline tests. Reconciled status, authority, Open Issues, limitations, dispositions, and remaining gates without substantively changing the approved addendum. This document-only revision performs no implementation or execution and grants no production-source edit, capture, replay, benchmark, performance-analysis, retained-result, remote/HPC, publication, staging, commit, or push authority. |
| 2026-09-23 | User / Codex | Recorded the user's contract-only approval of document snapshot `a829ffc01d09a9960e9bc5df04a15199a39833fd6942cb696a8c942bb1292c38` for `AMENDMENT_APPROVED` to `AMENDMENT_APPROVED`. OI-1 now replaces the former single pre-update/scalar-limiting-level requirement with the exact P1a/P1b/P1c, post-reset, direct-limiter, observed-`dtmin`/`nsubmix`, and offline minimizer-tie-set contract; OI-2 now fixes the conditional replay-eligible `R`/`C`/`U` estimand and denominators. Reconciled addendum schemas, interfaces, limitations, dispositions, and authority text. The addendum remains unapproved but is ready for a separate hash-bound `AMENDMENT_APPROVED` to `ADDENDUM_APPROVED` approval; this document-only revision grants no implementation, production-source edit, capture, replay, benchmark, analysis, retained-result, remote/HPC, staging, or commit authority. |
| 2026-09-22 | User / Codex | Performed the focused document-only blocker-resolution pass. Specified, pending separate approval, the exact P1a/P1b/P1c and post-reset capture checkpoints, direct limiter tie semantics, observed `dtmin`/`nsubmix`, and offline selector tie-set evidence; specified the conditional replay-eligible mixed-tail estimand with exact `R`/`C`/`U` accounting; replaced independently reset one-call qualification with source-anchored repeated production-sequence evidence using existing EAMxx per-subcycle hashes and separate unsynchronized overhead timing; and recovered the approved `d79a142...` bytes as verified local Git blob `35e6a2ef...`. Added an open-issues log and exact approval boundary. OI-1/OI-2 remain unapproved; this revision grants no implementation, source-edit, qualification, capture, replay, benchmark, analysis, remote/HPC, staging, or commit authority. |
| 2026-09-22 | User / Codex | Adjudicated consolidated addendum findings B1--B15 against the pending document and exact locally available locked source. Corrected binary framing, reference/status applicability, closed encodings, numerical diagnostics, census/expansion evidence, L128/L72 admission, timing and instrumentation tiers, source API/type facts, dependency provenance, and reporting limits; added source-anchored producer/lifetime and qualification proof requirements. The review also confirmed blockers requiring a separate amendment for the in-call capture checkpoint/limiter evidence and mixed-tail estimand, plus local verification of the exact approved-baseline bytes and the one-call qualification path. The addendum remains pending and not approval-ready; this document-only revision grants no implementation or execution authority. |
| 2026-09-14 | User / Codex | Recorded explicit approval of the production-state replay amendment snapshot `d79a14253ca3acb5108d54f37b3ddcc9aaa2a9d9c4ecef391721c27868df2138` for `AMENDMENT_PENDING` to `AMENDMENT_APPROVED`, then prepared the pending serialized-schema and C++-interface addendum. The addendum fixes canonical artifacts, two-phase capture payload and ownership, qualification/census/preregistration/sentinel/replay-run/result schemas, exact device/writer/validator interfaces, backend completion mapping, changed-file allowlists, statistical conventions, fail-closed tests, and remaining approvals without changing any scientific or numerical contract. This document-only revision grants no implementation or execution authority. |
| 2026-09-14 | User / Codex | Hardened the pending production-state replay amendment without changing its scientific or numerical contract: defined capture-owned synchronized pre-update and post-`dropmixnuc` transaction phases; added current `cldn` to Tier 1; separated locked physics from instrumentation-only provenance; froze census eligibility and statistical implementation bindings; replaced the circular approval list with a SHA-256-bound state sequence including non-transferable qualification and explicit census approval; required three-arm observational qualification and consumer-enforced evidence roles; reconciled replay-only `top_lev`/no-flux semantics; and updated Task 5R provenance. These document-only changes grant no implementation or execution authority. |
| 2026-09-11 | User / Codex | Made the pending production-state replay amendment approval-ready without changing its scientific or numerical contract: separated historical, previously approved standalone, and pending authority; added deterministic pilot-plus-preregistered expansion; specified hierarchical matched-case timing-ratio bootstrap and deterministic error aggregation; declared the locked current MAM4xx C target; tiered correctness, timing, and deferred explanatory instrumentation; and enumerated all approvals that still block production capture. |
| 2026-09-11 | Codex | Proposed a pending-approval production-state replay redesign: exact `dropmixnuc` pre-update capture, locked-source operator equality, separate `chi_mix`/`chi_act`/`chi_FE`, production-count FE versus one-step IMEX/BE, L128 census and censoring, accuracy/robustness/timing dominance statuses, provisional ne256 intervals, two-tolerance reference failure behavior, and the whole-ACI production-layout GPU gate. Existing analytic and 1,875-case gates remain regression evidence only; no source, numerical method, capture, run, or retained Task 5 result was changed. |
| 2026-08-25 | Codex | Synchronized the approved Task 5R amendment as implemented through `d55ef90e46`, with subsequent artifact-promotion hardening at `19c663a61ef73f99eb0b001d2947fd24b06e5e7b`: retained `ndrop_summary_v2` and `ndrop_phase_samples_v1` semantics, fixed executed-width/occupancy labels, combined IMEX timing, paired completion manifests, diagnostic-only phase typing, exact C++ median and stopping-policy validation, canonical output separation, and retained copies of the exact executable and dependency record. Pre-retained Python, ASan/UBSan, optimized, and optimized numerical gates pass; the clean detached retained run and all performance conclusions remain pending. |
| 2026-08-25 | Codex | Added the pending Task 5 evidence-remediation amendment after static review: fixed the executed-width contract at 25 aerosol pairs/51 fields while relabeling the varying input as nonzero-state occupancy; specified method-explicit linear phase modes with supported combined IMEX coefficient-copy/factorization/solve timing; required a repetition-level raw phase/fence sidecar; bound paired artifacts through immutable identity and completion manifests; and restricted all comparisons to primary uninstrumented time. Equations, solver arithmetic, and production code remain unchanged; explicit approval is required before implementation. |
| 2026-08-24 | Codex | Synchronized the Task 4 common cross-method gate as built: independent adaptive Dormand--Prince 5(4) physical-RHS oracle, analytic and two-tolerance reference checks, complete 1,875-case Cartesian matrix, per-field errors and diagnostic ledgers, all-method raw-positivity and fail-closed norm gates, three-run bitwise determinism, all-field first-order fits, pairwise Richardson common-solution checks, deterministic CSV/table outputs, and clean debug/sanitizer plus optimized evidence. FE, IMEX, and BE pass Task 4; only the local CPU timing phase is authorized next. |
| 2026-08-24 | Codex | Synchronized the Task 3 host Backward Euler implementation as built: signed 2-by-2 block rows, exact matrix-family reuse, unsafe-pivot failure, scalar droplet solve after modal-number blocks, independent dense and physical-equation oracles, raw positivity before projection, fail-closed diagnostics, bottom-stencil inspection, analytic first-order convergence, and clean locked-stack debug/sanitizer plus optimized evidence. The method-local gate passes; production physics remains unchanged. |
| 2026-08-24 | Codex | Synchronized the Task 2 FE/IMEX implementation as built: controller-only substep selection, signed total/clipping/operator ledgers, like-unit independent residuals, explicit roundoff/material classifications, fail-closed diagnostic arithmetic, scaled and zero-ledger tests, three-run bitwise state checks, and the locked current-MAM4xx rank-3 adapter at `top_lev=0,1,6`. Debug/sanitizer and optimized standalone gates pass; production physics remains unchanged. |
| 2026-08-18 | User / Codex | Replaced disposable `/private/tmp` screening storage with the approved persistent local split: temporary worktrees, dependency checkouts, and builds below `/Users/burr114/e3sm-scratch/`; retained results and provenance below `/Users/burr114/e3sm-results/`. Numerical and production-physics contracts are unchanged. |
| 2026-08-18 | User / Codex | Approved and froze the standalone cross-method screening contract: current and archived dependency revisions, controller-only first-active-level FE rate, signed no-flux IMEX and BE rows, block-Thomas and independent dense-oracle rules, diagnostic decompositions, matrix-family reuse, complete numerical and performance grids, acceptance thresholds, and test-only scope. Production physics remains unchanged. |
| 2026-08-17 | User / Codex | Documented the approved diagnostic-only implementation of the exact `nsubmix` count, per-physics-step sum, and per-physics-step maximum; rebased it onto current E3SM and MAM4xx; distinguished the repaired current paths from the two older defects; and kept all operator equations unchanged. |
| 2026-08-06 | Independent reviewer / Codex | Recorded strict independent `PASS` decisions for the corrected Forward Euler and IMEX gates at candidate `2c4427e` after fresh sanitizer/debug and optimized reruns, oracle-sensitivity mutations, and independent clipping and fractional-overlap hand checks. The Backward Euler gate may proceed. |
| 2026-08-06 | Codex | Remediated the second independent FE/IMEX gate review with mode-sensitive independent activation oracles, complete exact clipping diagnostics, an IMEX limiter-active direct oracle, focused finite-stage failure probes, a reproducible CMake/CTest driver, and fresh sanitizer/debug plus optimized evidence. Both gates remain pending independent review. |
| 2026-08-03 | Codex | Clarified that the frozen isolated operator is linear; documented the sequential scalar-tridiagonal IMEX structure, the block-tridiagonal Backward Euler structure, the EKAT versus integration-library boundary, current FE/IMEX/BE/hybrid trade-offs, and deferred higher-order implicit candidates. No production implementation was authorized. |
| 2026-07-30 | Codex | Synchronized the remediation as built: added exact zero-cloud resuspension/reset, deterministic finite-value failures, scaled and independently assembled residuals, per-field Forward Euler clipping diagnostics, all-field analytic convergence, coupled FE/IMEX convergence, and debug/sanitizer plus optimized evidence. The automated suites pass; independent gate re-review remains required. |
| 2026-07-30 | User / Codex | Recorded explicit approval of the validation-gate remediation amendment, authorizing standalone harness and test remediation while keeping Backward Euler blocked until both remediated gates pass. |
| 2026-07-29 | Codex | Added the pending remediation amendment after independent review failed both completed gates: made zero-cloud resuspension/reset part of the canonical state transition, required fail-loud handling of nonfinite derived quantities and diagnostics, expanded all-field activation and coupled convergence coverage, restored per-field Forward Euler clipping requirements, required an independently assembled residual check, and required coherent committed gate provenance. |
| 2026-07-28 | Codex | Recorded successful standalone Forward Euler and IMEX gates, including per-field convergence, residual, ledger, positivity, randomized-regime, sanitizer, optimized-build, and exact-current regression evidence; Backward Euler remains unstarted. |
| 2026-07-23 | User / Codex | Recorded explicit approval of the standalone numerical-method design, clearing Phase 0 for staged test and validation work but not production implementation. |
| 2026-07-21 | User / Codex | Linked the document to Phase 0 of the isolated-integrator plan; clarified operator versus integrator, high-level acceptance, validity regimes, excluded defects, and non-authorization; restored production clipping only for Forward Euler and fail-loud positivity handling for test-only IMEX and Backward Euler. |
| 2026-07-21 | User / Codex | Narrowed the proposed study to Mac-local standalone screening, replaced required RK4 with analytic and self-convergence evidence, required clipping with explicit conservation-error accounting, and deferred GPU implementation, bottom-source resolution, and downstream-fixer tracing. |
| 2026-07-21 | User / Codex | Drafted the pending-review isolated Forward Euler, IMEX, and Backward Euler numerical design, staged validation gates, convergence criteria, and complete GPU timing contract; no test or production code was changed. |
| 2026-07-21 | User / Codex | Added a provisional GPU-aware Forward Euler, IMEX, and fully implicit Backward Euler performance model across default resolutions and turbulence regimes; no production code was changed. |
| 2026-07-21 | User / Codex | Documented the adaptive explicit vertical-mixing subloop in Fortran EAMv3 and its main differences from the EAMxx Kokkos implementation. |
| 2026-07-21 | User / Codex | Verified that E3SM v3.0.0 Fortran applies both signs of each activation transfer and does not contain the EAMxx disabled-mixing switch defect. |
| 2026-07-21 | User / Codex | Clarified that “cloud-fraction transfer” means same-level interstitial/cloud-borne repartitioning during cloud growth, decay, or disappearance, not vertical transport. |
| 2026-07-15 | Codex | Documented the explicit `eddy_diff_heat` to interface diffusivity to `tinv`/`nsubmix` path and the cloud-base override that can supersede it. |
| 2026-07-15 | User / Codex | Restricted the current experiment plan to matched five-day runs, removed separate extended and focused-output streams, and set existing `ndrop` context to identical hourly instantaneous output. |
| 2026-07-15 | Codex | Documented the source-traced `w_variance` to `wsub` to `wtke` calculation, the hard-coded activation-velocity floor, the diagnostics missing for EAMxx-versus-EAMv3 causal attribution, and the comprehensive existing-output field/cadence contract; no implementation code was changed. |
| 2026-07-15 | Codex | Added and prioritized the GPU-aware diagnostics contract, reduced existing diagnostics to a performance-only subset with required `eddy_diff_heat`/`w_variance` context, and documented how the `top_lev=0` invalid access can be masked and propagate; no production instrumentation was implemented. |
| 2026-07-14 | Codex | Corrected the accumulation/Aitken mode labels during the second peer-review pass; the constituent counts and executed test mapping already followed the source `{7, 4, 7, 3}` order. |
| 2026-07-14 | User / Codex | Recorded explicit approval of the reconstructed design, corrected Kokkos subcycling contract, and focused validation plan. |
| 2026-07-14 | Codex | Recorded the executed one/even/odd disabled-mixing matrix, mass/number conservation and destination-array violations, and final confirmed-defect classification. |
| 2026-07-14 | Codex | Recorded Kokkos-bounds-checked 128-level boundary evidence and classified the production-reachable `top_lev=0` access as a confirmed defect; retained disabled mixing as a code-reading risk pending ledger tests. |
| 2026-07-14 | Codex | Corrected the as-executed subcycle selection after instrumented boundary work exposed that the Kokkos minimum reducer overwrites the preassigned `dtmicro` cap; expanded validation to one, even, and odd substep counts. |
| 2026-07-14 | Codex | Reconstructed the existing EAMxx-MAM4xx `ndrop` indexing, subcycling, paired-field updates, coupling, and validation requirements from source; marked for review before investigation-specific tests. |
