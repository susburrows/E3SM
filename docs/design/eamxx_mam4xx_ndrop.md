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
> review on 2026-08-06 at candidate commit `2c4427e`. The standalone Backward
> Euler gate may now begin; production implementation remains unauthorized.

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

#### Historical evidence-revision behavior

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

#### Historical evidence-revision behavior

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

The harness must not reproduce the known `eddy_diff(-1)` access.
Exact-current-function comparison is restricted to a safe active range such
as `top_lev=1`; successful isolated tests must not be reported as validation of
the defective production `top_lev=0` path.

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
\begin{bmatrix} C^{n+1} \\ I^{n+1} \end{bmatrix}
=
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
decision sequence is therefore: finish and validate Backward Euler, compare
its damping and matched-error cost with FE and IMEX, and add SDIRK2/ESDIRK2 or
TR-BDF2 only if the evidence shows that a higher-order implicit method is
needed. RK4 remains outside the planned screening suite.

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
decisions to both gates. Backward Euler has not begun. Executed evidence is
recorded in
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
| Fully implicit Backward Euler | Removes the mixing restriction and evaluates activation transfer at the new state; likely more robust when explicit-source lag matters | Requires custom 2-by-2 block-tridiagonal solves; more flops and workspace; first-order stiff damping can reduce accuracy; positivity is not guaranteed | Fixed work, but the coupled block solver is not supplied by the current EKAT scalar interface | Initial design gate approved; standalone implementation and validation remain separate work |
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

## Known Limitations and Open Decisions

- The intended valid range of `top_lev` must be approved. The current L128
  default requires that zero either be valid or be rejected earlier by a
  redesigned configuration.
- The intended top-interface boundary condition must be stated explicitly.
  The legacy implementation indicates zero flux, but the C++ contract does not.
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
- The existing tests do not span the boundary matrix or verify paired aerosol
  ledgers.
- The `Kokkos::Min` reducer overwrites, rather than incorporates, the
  preassigned `dtmin=dtmicro`. This permits `nsubmix=1` under weak mixing and
  differs from the apparent source-level intent; the numerical and performance
  implications require separate follow-up.

## Scientific References

- Abdul-Razzak, H., and S. J. Ghan (1998), A parameterization of aerosol
  activation: 1. Single aerosol type, *Journal of Geophysical Research*,
  103(D6), 6123-6131.
- Abdul-Razzak, H., and S. J. Ghan (2000), A parameterization of aerosol
  activation: 2. Multiple aerosol types, *Journal of Geophysical Research*,
  105(D5), 6837-6844.
- Ghan, S. J., C. C. Chuang, and J. E. Penner (1993), A parameterization of
  cloud droplet nucleation, *Atmospheric Research*, 30, 198-221.
- Alexander, R. (1977), [Diagonally implicit Runge-Kutta methods for stiff
  O.D.E.'s](https://epubs.siam.org/doi/10.1137/0714068), *SIAM Journal on
  Numerical Analysis*, 14, 1006-1021.
- [SUNDIALS ARKODE user documentation](https://sundials.readthedocs.io/en/latest/arkode/Usage/General.html)
  and [Kokkos installation notes](https://sundials.readthedocs.io/en/latest/sundials/Install_link.html).
- [PETSc `TSARKIMEX` documentation](https://petsc.org/release/manualpages/TS/TSType/).
- Kopecz, S., and A. Meister (2018), [On order conditions for modified
  Patankar-Runge-Kutta schemes](https://www.sciencedirect.com/science/article/pii/S0168927417301861),
  *Applied Numerical Mathematics*, 123, 159-179.

## Revision History

| Date | Author / Agent | Change Summary |
|---|---|---|
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
