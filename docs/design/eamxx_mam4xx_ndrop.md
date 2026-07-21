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
   ping-pong buffer zero;
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
introduced by the C++ port. Important port differences include the C++
boundary representation, Kokkos execution structure, and the EAMxx
`enable_aero_vertical_mix` option.

### Vertical indexing and boundaries

Most active-level loops use the half-open Kokkos range
`[top_lev, pver)`. Neighbor indices are usually clamped as

```text
kp1 = min(k + 1, pver - 1)
km1 = max(k - 1, top_lev)
```

The upper-interface mixing coefficient is currently formed as

```text
eddy_diff_km(k) = zn(k) * eddy_diff(k - 1) * csbot(k - 1) * zs(km1)
```

Consequences of the as-built indexing are:

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
form of mixing rate. The corresponding C++ boundary initialization or guard is
not evident in the inspected call path.

### `enable_aero_vertical_mix` behavior

Inside the explicit substep loop, the switch controls these operations:

| Update | Switch condition in current code |
|---|---|
| Mix/update cloud-droplet number `qcld` | Conditional |
| Mix/update interstitial aerosol number | Conditional |
| Mix/update interstitial aerosol mass | Conditional |
| Mix/update cloud-borne aerosol number | Unconditional |
| Mix/update cloud-borne aerosol mass | Unconditional |
| No-cloud cloud-borne-to-interstitial transfer after subcycling | Unconditional |

The activation source has opposite signs in the interstitial and cloud-borne
equations. Consequently, disabling only the interstitial update while retaining
the cloud-borne update does not preserve the paired source cancellation in the
equations actually advanced.

The aerosol arrays use two time buffers. At entry, only buffer zero is loaded
from the EAMxx state. The source and destination buffer indices swap after each
substep. When aerosol vertical mixing is disabled, the interstitial destination
buffer is not written by the conditional `explmix` calls, but subsequent
substeps and final tendency construction can read the swapped buffer. Existing
validation fixtures initialize both buffers, which does not reproduce this
production initialization pattern.

Cloud-fraction-driven same-level repartitioning occurs outside the
switch-guarded explicit calls. This is not vertical transport. For each paired
interstitial/cloud-borne constituent `(I,C)` at one model level:

- when cloud fraction grows by more than `0.01`, activation computes a positive
  amount `Delta` and applies `I <- I-Delta`, `C <- C+Delta`;
- when cloud fraction shrinks, a fraction of `C` is returned with the same
  paired form but `Delta<0`; and
- when cloud is absent, the code applies `I <- I+C`, `C <- 0` and removes the
  corresponding cloud-droplet number.

Before positivity clipping or the disabled-buffer defect, these local paired
operations preserve `I+C`. They remain active when
`enable_aero_vertical_mix=false` because the switch guards selected explicit
neighbor-mixing updates, not all activation and cloud-phase repartitioning in
`dropmixnuc`. Retaining conservative local repartitioning while disabling
turbulent transport is physically coherent; the confirmed defect is that the
substep loop retains the cloud-borne `+S` update while suppressing its paired
interstitial `-S` and droplet-number updates and can read an unwritten ping-pong
buffer.

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

The current `dropmixnuc` signature is declared in
`externals/mam4xx/src/mam4xx/ndrop.hpp:1438-1477`. The EAMxx wrapper signature
is in
`components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp:199-231`.

### Source evidence map

The following anchors apply to the inspected revisions and are intended to make
the investigation reproducible:

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
| Explicit-step selection and ping-pong swap | `externals/mam4xx/src/mam4xx/ndrop.hpp:1264-1299` |
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
| `enable_aero_vertical_mix` | Controls selected explicit updates | boolean | scalar | Default `true`; precise disabled-state semantics undocumented |
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
| Interstitial aerosol number/mass | `progs.q_aero_i` | Paired transfer and optional explicit mixing in `raercol` | Returned change divided by `dt`, then applied as a tendency |
| Cloud-borne aerosol number/mass | `progs.q_aero_c` | Paired transfer and unconditional explicit cloud-borne update in `raercol_cw` | Returned as a replacement state after heterogeneous freezing |
| Cloud-droplet number | `qcld` work field | Activation source and optional explicit mixing | Returned through the ACI process output path |
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
| `raercol`, `raercol_cw` | Interstitial/cloud-borne ping-pong state | # kg-1 or kg kg-1 | `pver x 2 x 25` |
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

`enable_aero_vertical_mix=false` is accepted by the EAMxx parameter interface
and is reachable through configuration. No current default using `false` has
been identified. A deterministic standalone test can exercise the state. The
namelist description says that the switch disables vertical mixing of
interstitial aerosol and liquid number during activation, which is consistent
with some conditional calls but does not specify the required behavior of the
cloud-borne pair, paired activation sources, or ping-pong buffers.

### Executed boundary evidence

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

### Executed disabled-mixing evidence

The same instrumented 128-level driver exercised the exact
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
establish the paired-ledger violation without relying on poisoned or
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

### Disabled vertical-mixing matrix

Use deterministic, nonzero interstitial and cloud-borne number and mass fields.
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
- independence from initial contents of the destination ping-pong buffer;
- consistency across even and odd `nsubmix`; and
- agreement with the reviewed meaning of `enable_aero_vertical_mix=false`.

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
separately. For `top_lev=0`, retain the source trace from the default L128
namelist through `MAMAci::run_impl` into `dropmixnuc`. For the disabled switch,
retain a source trace showing that a legal EAMxx parameter list can supply
`false`; explicitly distinguish this from a default production configuration.

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

No currently registered field exposes actual `nsubmix`, `max_tinv`, `dtmix`,
the limiting level, ACI-local `wsub`/`wtke`, inner-loop limiter activity, or a
timer scoped specifically to `dropmixnuc`.

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
| Actual `nsubmix`, `max_tinv`, `dtmix`, and limiting level | Separate activation behavior from adaptive explicit-mixing work and timestep sensitivity | Internal only |

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

The following quantities are not currently available as EAMxx output fields or
appropriately scoped timers. Production instrumentation remains subject to a
separate design approval. GPU implementations should use device-resident
accumulators and infrequent output; per-call host copies or MPI reductions
would perturb the measurement.

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

If production instrumentation is later approved, reduced device-accumulated
`nsubmix` and work histograms can be written hourly in the same five-day
matrix. A short warmed-up profiler segment or matched instrumentation-off/on
segment remains the recommended follow-up for quantifying kernel behavior and
diagnostic overhead, but it is not part of the current run matrix.

### Validity gates for the planned runs

- The current L128 `top_lev=0` boundary-index defect can directly corrupt
  `eddy_diff_km`, `max_tinv`, and `nsubmix`. The defect must be corrected, or a
  demonstrably equivalent safe boundary must be used, before interpreting
  `ndrop` scaling from these experiments.
- Record actual `nsubmix`. The inspected Kokkos reduction overwrites the
  preassigned `dtmicro` value, so weak mixing can produce one explicit substep.
- Keep `enable_aero_vertical_mix` fixed and true unless the separately confirmed
  disabled-switch defect has been corrected and tested.
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
- The approved investigation contract interprets
  `enable_aero_vertical_mix=false` as disabling turbulent transport while
  retaining conservative paired activation and cloud/interstitial phase
  transfers. Production naming and namelist documentation should be clarified
  when the defect is fixed.
- GPU bounds behavior may differ from a bounds-checked host build. Host
  sanitizer evidence should be paired with deterministic invariant tests that
  are portable to GPU builds.
- Existing MAM4xx validation initializes both ping-pong buffers and therefore
  cannot by itself establish the safety of production initialization when the
  switch is false.
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

## Revision History

| Date | Author / Agent | Change Summary |
|---|---|---|
| 2026-07-21 | User / Codex | Documented the adaptive explicit vertical-mixing subloop in Fortran EAMv3 and its main differences from the EAMxx Kokkos implementation. |
| 2026-07-21 | User / Codex | Verified that E3SM v3.0.0 Fortran applies both signs of each activation transfer and does not contain the EAMxx disabled-mixing switch defect. |
| 2026-07-21 | User / Codex | Clarified that “cloud-fraction transfer” means same-level interstitial/cloud-borne repartitioning during cloud growth, decay, or disappearance, not vertical transport. |
| 2026-07-15 | Codex | Documented the explicit `eddy_diff_heat` to interface diffusivity to `tinv`/`nsubmix` path and the cloud-base override that can supersede it. |
| 2026-07-15 | User / Codex | Restricted the current experiment plan to matched five-day runs, removed separate extended and focused-output streams, and set existing `ndrop` context to identical hourly instantaneous output. |
| 2026-07-15 | Codex | Documented the source-traced `w_variance` to `wsub` to `wtke` calculation, the hard-coded activation-velocity floor, the diagnostics missing for EAMxx-versus-EAMv3 causal attribution, and the comprehensive existing-output field/cadence contract; no implementation code was changed. |
| 2026-07-15 | Codex | Added and prioritized the GPU-aware diagnostics contract, reduced existing diagnostics to a performance-only subset with required `eddy_diff_heat`/`w_variance` context, and documented how the `top_lev=0` invalid access can be masked and propagate; no production instrumentation was implemented. |
| 2026-07-14 | Codex | Corrected the accumulation/Aitken mode labels during the second peer-review pass; the constituent counts and executed test mapping already followed the source `{7, 4, 7, 3}` order. |
| 2026-07-14 | User / Codex | Recorded explicit approval of the reconstructed design, corrected Kokkos subcycling contract, and focused validation plan. |
| 2026-07-14 | Codex | Recorded the executed one/even/odd disabled-mixing matrix, paired-ledger and destination-buffer violations, and final confirmed-defect classification. |
| 2026-07-14 | Codex | Recorded Kokkos-bounds-checked 128-level boundary evidence and classified the production-reachable `top_lev=0` access as a confirmed defect; retained disabled mixing as a code-reading risk pending ledger tests. |
| 2026-07-14 | Codex | Corrected the as-executed subcycle selection after instrumented boundary work exposed that the Kokkos minimum reducer overwrites the preassigned `dtmicro` cap; expanded validation to one, even, and odd substep counts. |
| 2026-07-14 | Codex | Reconstructed the existing EAMxx-MAM4xx `ndrop` indexing, subcycling, paired-field updates, coupling, and validation requirements from source; marked for review before investigation-specific tests. |
