# EAMxx MAM4xx ndrop ACI Cost Analysis

Date: 2026-06-16

Last updated: 2026-08-03

## Scope

This report summarizes the findings from a code and configuration review of
EAMxx MAM4xx aerosol-cloud interaction (ACI), focusing on the `ndrop`
activation/mixing path, its explicit mixing subcycling, and the turbulence
quantities that drive its cost. This is an analysis artifact only; no
implementation code was changed.

The focused bounds and paired-ledger investigation is reported in
[`eamxx_mam4xx_ndrop_defect_brief.md`](eamxx_mam4xx_ndrop_defect_brief.md).
The reconstructed as-built equations, interfaces, validation contract, and
diagnostic design are maintained in
[`../design/eamxx_mam4xx_ndrop.md`](../design/eamxx_mam4xx_ndrop.md).

The main code path reviewed is:

```text
components/eamxx/src/physics/mam/eamxx_mam_aci_process_interface.cpp
components/eamxx/src/physics/mam/eamxx_mam_aci_functions.hpp
externals/mam4xx/src/mam4xx/ndrop.hpp
```

For comparison with E3SMv3, the Fortran EAM path reviewed is:

```text
components/eam/src/physics/cam/microp_aero.F90
components/eam/src/physics/cam/ndrop.F90
```

## Current Investigation Snapshot

The executable boundary and paired-ledger evidence was collected against:

| Repository | Evidence revision |
|---|---|
| E3SM | `12e44bb469efb0d6c0cb94c7b796a86967883ff5` |
| `externals/mam4xx` | `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db` |

On 2026-07-15, the current checkout was E3SM revision
`781c3f5c4cc07f5bc8b5a8d748b0aa923f9b8522` with the same MAM4xx gitlink.
The relevant EAMxx ACI interface, wrapper, defaults, and component-configuration
files have no changes between the evidence revision and this checkout. The
source reachability conclusions and executed MAM4xx evidence therefore remain
applicable.

The current experiment matrix is limited to matched five-day `ne256` runs with
an experiment-specific fixed physics timestep of 600 s:

| Case | MAC/ACI subcycles | ACI timestep |
|---|---:|---:|
| `SUB-02` | 2 | 300 s |
| `SUB-04` | 4 | 150 s |
| `SUB-06` | 6 | 100 s |
| `SUB-12` | 12 | 50 s |

The 600 s value is a test setting, not the default `ne256` physics timestep
used in the resolution-estimate table below. No separate extended runs are
currently planned. No production GPU profile or actual `nsubmix` distribution
has yet been collected.

The minimum existing-output contract for every case is:

- complete configuration, revision, compiler/backend, GPU, and MPI/thread
  provenance once per run;
- end-of-run `EAMxx::run`, `mac_aero_mic`, child-process, and I/O timer
  summaries; and
- identical hourly instantaneous `eddy_diff_heat` and `w_variance`, with
  `cldfrac_liq` and `cldfrac_liq_prev` recommended to identify cloud-base
  regimes.

At `ne256`, even these four native-grid 3-D fields can create substantial I/O.
Use identical cadence, remapping or supported sampling, compression, and I/O
layout in every case, and retain the I/O timers. These fields provide context
but cannot reconstruct actual inner-loop work because normal history output
sees only the final MAC subcycle.

## Main Conclusions

`ndrop` in EAMxx/MAM4xx does not contain an obvious nonlinear convergence
iteration. The most obvious expensive calculation is an explicit vertical
mixing subcycle loop in `mam4::ndrop::update_from_explmix`.

The number of explicit mixing substeps, `nsubmix`, is controlled by the fastest
local vertical exchange rate in the column:

```text
tinv(k) = eddy_diff_kp(k) + eddy_diff_km(k)
max_tinv = max_k tinv(k)
```

where `tinv` has units of `s^-1`. Physically, `tinv(k)` is the first-order loss
rate from layer `k` by turbulent exchange with the neighboring layers above and
below. The inverse, `1/tinv`, is the local mixing timescale.

The as-executed subcycling logic is:

```text
dtmin = Kokkos minimum-reduction identity (maximum representable Real)
if tinv(k) > 1e-6:
    dtmin = min(dtmin, 1/tinv(k))
dtmix = 0.9 * dtmin
nsubmix = int(dtmicro / dtmix) + 1
dtmix = dtmicro / nsubmix
```

Although the source assigns `dtmin = dtmicro` before the reduction, the
`Kokkos::Min<Real>` reducer initializes its accumulator independently and
overwrites that scalar. Thus `dtmicro` is not included as an initial minimum.
This behavior was confirmed against Kokkos 4.5.1 and with an instrumented
128-level driver.

If the maximum exchange rate dominates,

```text
nsubmix = floor(dtmicro * max_tinv / 0.9) + 1
```

`nsubmix` can therefore be 1 under weak mixing. If no level satisfies
`tinv > 1e-6 s^-1`, the reduction identity is retained and the same one-step
result follows. There is no hard cap on `nsubmix` in the C++ EAMxx/MAM4xx
implementation.

Two correctness findings are now established by executed focused tests:

1. **Boundary indexing is a confirmed defect.** The default L128 input
   `top_lev=0` reaches `eddy_diff(-1)` and `csbot(-1)`. A Kokkos
   bounds-checked 128-level execution aborts at index `-1`.
2. **Disabled vertical mixing is a confirmed defect.** The public
   `enable_aero_vertical_mix=false` setting advances cloud-borne aerosol while
   skipping paired interstitial and droplet-number updates. Deterministic one-,
   two-, and three-substep tests failed all 25 aerosol pair ledgers, while the
   enabled controls conserved them to roundoff.

The level-zero invalid read need not crash an optimized GPU run. Bounds checks
are normally disabled, and the row-major scratch layout can make a negative
row-subview index alias the preceding column's surface value, where
`eddy_diff` is explicitly zero. For a finite bad coefficient, the direct
upper-boundary flux term is also algebraically zero because `km1=0` and
`overlapm(0)=1`. The main finite-value propagation is therefore through the
column-wide `max_tinv`, `nsubmix`, and `dtmix`: a positive value can add
unnecessary work, while a negative value can under-subcycle and trigger
forward-Euler clipping. NaN, infinity, unmapped memory, or execution-order
changes can produce more severe or nondeterministic symptoms.

## ACI Time Step and Subcycling

`dtmicro` in `dropmixnuc` is the timestep passed to the `mam4_aci` process. It
is not the internal P3 microphysics substepping timestep. If `mam4_aci` is
inside a subcycled atmosphere process group, `dtmicro` is the already-subcycled
group timestep.

For MAM4xx compsets, `mam4_aci` is normally inside the `mac_aero_mic` group:

```text
tms, shoc, cld_fraction, mam4_aci, p3
```

for PG2 MAM4xx configurations. Default `mac_aero_mic` subcycles found in
`components/eamxx/cime_config/namelist_defaults_eamxx.xml` are:

| Resolution | MAC/ACI subcycles |
|---|---:|
| `ne4np4` | 24 |
| `ne30np4` | 12 |
| `ne120np4` | 6 |
| `ne256np4` | 6 |
| `ne512np4` | 2 |
| `ne1024np4` | 1 |

Thus, under the default MAM4xx process layout, `dropmixnuc` is called once per
MAC/ACI subcycle, and each call may itself run `nsubmix` explicit mixing
iterations.

### Fortran EAMv3 also has the explicit mixing subloop

The Fortran implementation calculates `dtmix` and `nsubmix` from its vertical
mixing rates and then executes `do n=1,nsubmix`. Within every substep it updates
droplet number and calls `explmix` for each cloud-borne and interstitial aerosol
number and mass field. The loop is therefore not new to EAMxx.

The principal implementation differences are that Fortran uses ordinary loops
and array operations instead of Kokkos team/vector ranges, initializes explicit
zero-flux interface cells at both column boundaries, and always advances both
aerosol forms. Its apparent limit of 100 applies only to a diagnostic histogram
index; it does not cap the actual `nsubmix` loop count.

## Inputs to nsubmix

The immediate inputs are:

```text
dtmicro
eddy_diff_kp(k)
eddy_diff_km(k)
```

The mixing rates are computed approximately as:

```text
eddy_diff_kp(k) = zn(k) * eddy_diff(k)     * csbot(k)     * zs(k)
eddy_diff_km(k) = zn(k) * eddy_diff(k - 1) * csbot(k - 1) * zs(km1)
```

where:

```text
zn = gravity * rpdel
rpdel = 1 / pdel
zs = inverse distance between midlevels
csbot = interface-density-related factor
```

Away from cloud-base activation, `eddy_diff` comes from `eddy_diff_heat` after
interpolation to interfaces and is bounded to `[0.01, 100] m2/s` in
`dropmixnuc`.

At cloud-base activation points, `update_from_cldn_profile` overwrites:

```text
eddy_diff = wtke / zs
```

This path bypasses the `[0.01, 100] m2/s` limiter. Since `zs` is inverse length,
this is equivalent to an effective diffusivity proportional to
`wtke * vertical_spacing`. It makes cloud-base exchange rates especially
sensitive to the diagnosed subgrid vertical velocity `wtke`.

## Fields Updated by ndrop Vertical Mixing

With the default `enable_aero_vertical_mix = true`, the explicit mixing loop in
`mam4::ndrop::update_from_explmix` updates three classes of state:

1. Cloud-droplet number mixing ratio, represented internally by `qcld` and
   returned to EAMxx as part of `nc_nuceat_tend`.
2. Interstitial modal aerosol number and species mass mixing ratios.
3. Cloud-borne modal aerosol number and species mass mixing ratios.

The complete MAM4 modal field mapping is:

| Mode | Interstitial fields | Cloud-borne fields |
|---|---|---|
| 1, accumulation | `num_a1`; `so4_a1`, `pom_a1`, `soa_a1`, `bc_a1`, `dst_a1`, `nacl_a1`, `mom_a1` | `num_c1`; `so4_c1`, `pom_c1`, `soa_c1`, `bc_c1`, `dst_c1`, `nacl_c1`, `mom_c1` |
| 2, Aitken | `num_a2`; `so4_a2`, `soa_a2`, `nacl_a2`, `mom_a2` | `num_c2`; `so4_c2`, `soa_c2`, `nacl_c2`, `mom_c2` |
| 3, coarse | `num_a3`; `dst_a3`, `nacl_a3`, `so4_a3`, `bc_a3`, `pom_a3`, `soa_a3`, `mom_a3` | `num_c3`; `dst_c3`, `nacl_c3`, `so4_c3`, `bc_c3`, `pom_c3`, `soa_c3`, `mom_c3` |
| 4, primary carbon | `num_a4`; `pom_a4`, `bc_a4`, `mom_a4` | `num_c4`; `pom_c4`, `bc_c4`, `mom_c4` |

This is 25 aerosol constituents in each reservoir: four modal number fields and
21 modal species-mass fields. The mode/species mapping is defined by
`mam4::mode_aero_species` in `externals/mam4xx/src/mam4xx/aero_modes.hpp`; the
EAMxx field-name mapping is in
`components/eamxx/src/physics/mam/mam_coupling.hpp`.

## Turning Off Vertical Mixing with `enable_aero_vertical_mix`

The default value is `true`. Although `false` is an accepted configuration,
the current implementation does not consistently disable vertical mixing and
does not preserve the aerosol state correctly. The disabled path is therefore
a confirmed defect and should not be used in production.

### What happens in the explicit substep loop

For each of the 25 aerosol number or mass quantities, let `I` denote its
interstitial amount, `C` its cloud-borne amount, and `S` the activation rate.
With the switch enabled, the relevant source terms have opposite signs:

```text
interstitial aerosol:  dI/dt includes -S
cloud-borne aerosol:   dC/dt includes +S
cloud-droplet number:  dNc/dt includes the activated number source
```

The decrease in `I` balances the increase in `C`, so activation transfers
aerosol between the two forms rather than creating aerosol.

Setting `enable_aero_vertical_mix=false` changes the loop as follows:

| Quantity | Behavior when the switch is `false` |
|---|---|
| Cloud-droplet number `qcld` | Skips both its vertical-mixing update and the droplet-number source calculated in the substep. |
| Interstitial aerosol `raercol` | Skips its vertical-mixing update, the `-S` activation change for aerosol number and mass, and terms that return cloud-borne aerosol from non-overlapping cloud regions. |
| Cloud-borne aerosol `raercol_cw` | Still applies vertical mixing and the `+S` activation change for aerosol number and mass. |
| Cloud-borne aerosol where cloud fraction is zero | Still enters the code that returns it to interstitial aerosol after subcycling, although that calculation can already contain an invalid interstitial work value. |

Consequently, the option does not actually turn off all aerosol vertical
mixing. More importantly, it retains the positive cloud-borne activation term
while removing the equal interstitial loss. Ignoring other terms, one substep
therefore produces

```text
change in (I + C) = dtmix * S
```

instead of zero. This is a spurious source of aerosol number and mass.
The focused conservation test used complete cloud overlap, so this activation
imbalance was demonstrated independently of any additional error caused by
skipping the cloud-overlap return terms.

The calculation of `nsubmix`, the substep loop, the aerosol-source calculations,
and all 25 cloud-borne `explmix` updates still run. Setting the option to `false`
therefore neither provides a correct physical state nor removes most of the
adaptive loop structure for a performance experiment.

### A separate work-array failure

The subcycling code uses two work arrays for the old and new aerosol values,
but `dropmixnuc` initializes only array zero from the EAMxx state. Subcycling
starts with array zero as the old state and array one as the destination.

With the switch off, the first substep writes new cloud-borne values into array
one but does not write new interstitial values there. If `nsubmix=1`, the final
interstitial tendencies are calculated from that unwritten destination. If
there are additional substeps, array one becomes the old interstitial state
used by the next activation-source calculations. Because neither destination
array receives an interstitial substep update, later output can then use stale
or unwritten values depending on the substep count. Thus the failure is not
limited to even or odd substep counts.

Focused tests separated this work-array problem from the uncancelled activation
source:

- with the second array initialized to a valid copy of the initial state, all
  25 aerosol number and mass conservation checks failed for `nsubmix=1`, `2`,
  and `3` when the switch was off;
- otherwise identical enabled cases conserved interstitial plus cloud-borne
  aerosol to roundoff; and
- changing only the initial contents of the normally unwritten interstitial
  array changed the disabled result at every tested substep count.

These results establish both an aerosol-conservation error and a dependence on
an unwritten work array. Keep the switch at its default value, `true`, until the
disabled path is corrected and tested. The executed evidence is documented in
[`eamxx_mam4xx_ndrop_defect_brief.md`](eamxx_mam4xx_ndrop_defect_brief.md).

### Changes made before and after the explicit vertical mixing loop

Separate code within `ndrop` responds to changes in cloud fraction. These
changes occur within one model level; they are not vertical transport between
neighboring levels:

- **Growing cloud:** activated aerosol moves from interstitial to cloud-borne,
  and cloud-droplet number increases.
- **Shrinking cloud:** some cloud-borne aerosol returns to interstitial, and
  cloud-droplet number decreases.
- **No cloud:** cloud-borne aerosol returns to interstitial, and cloud-droplet
  number is set to zero.

For an aerosol number or mass quantity, these local changes use the form

```text
I <- I - Delta
C <- C + Delta
```

where `Delta` is positive for transfer into cloud and negative for return to
interstitial aerosol. These two changes preserve `I+C` when applied to valid
inputs. They are not the source of the demonstrated conservation defect. The
problem occurs inside the explicit substep loop, followed by possible use of
the unwritten interstitial destination in the post-loop no-cloud calculation.

### Comparison with Fortran EAMv3

The E3SM v3.0.0 Fortran `ndrop` code does not have this disabled-path defect.
It has no equivalent `enable_aero_vertical_mix` switch. For each aerosol number
or mass field, the Fortran code always applies both sides of activation: it
adds the activated amount to cloud-borne aerosol and subtracts the same amount
from interstitial aerosol. The two changes therefore cancel in their sum unless
the non-negativity correction or a boundary flux intervenes.

This conclusion was checked directly at E3SM v3.0.0 revision
`399d4301138617088dd93214123d6c025e061302`. It is specific to the asymmetric
disabled path in EAMxx; it is not a claim that every Fortran `ndrop` operation
is exactly conservative in all circumstances.

### How the results reach the EAMxx state

| Internal result | How EAMxx applies it |
|---|---|
| Final `qcld` | Converted to `nc_nuceat_tend`; P3 applies the tendency to cloud-droplet number `nc` where liquid cloud water is present. |
| Final interstitial aerosol `raercol` | Converted to `ptend_q`; MAM4xx ACI applies that tendency to interstitial aerosol tracers. |
| Final cloud-borne aerosol `raercol_cw` | Copied through `qqcw_fld_work` and directly replaces the cloud-borne fields. |

Consequently, `ptend_q` is not a vertical-mixing-only tendency. It combines
activation, same-level transfer caused by cloud-fraction changes, and explicit
vertical mixing. Likewise, `nc_nuceat_tend` is the net `ndrop` droplet-number
tendency. Its diagnostic decomposition is approximately

```text
nc_nuceat_tend = nsource + ndropmix
```

apart from the `max(qcld,1e-6)` floor used in the total tendency. `nsource` and
`ndropmix` are diagnostics; EAMxx does not apply them separately.

The loop does not update cloud-water mass `qc`, water vapor `qv`, temperature,
rain, ice, or aerosol gases such as `SO2`, `H2SO4`, `DMS`, and `SOAG`. SHOC
fields such as `eddy_diff_heat`, `w_variance`, and `tke` drive the calculation
but are not themselves updated by `ndrop`.

### Separate top-boundary issue

The `top_lev` indexing defect is independent of the switch behavior. In the
default L128 configuration, `top_lev=0` causes the first iteration to read
`eddy_diff(-1)` and `csbot(-1)`. In L72, `top_lev=6` keeps the index in bounds,
but `dropmixnuc` does not locally initialize the value at `top_lev-1`.

## Resolution Dependence

Vertical resolution can directly increase `nsubmix`, because the mixing rate
scales approximately like:

```text
tinv ~ K / dz^2
```

for diffusivity-limited mixing, or like:

```text
tinv ~ wtke / dz
```

in the cloud-base overwrite path.

Horizontal resolution does not directly enter `nsubmix` except through the
default timestep, process subcycling, and changes to resolved meteorology or
turbulence. It does, however, directly increase total cost by increasing the
number of columns.

## Vertical Grid Considerations

For MAM4xx EAMxx configurations in this checkout, the default vertical grid is
L128:

```text
components/eamxx/cime_config/config_component.xml
SCREAM_NUM_VERTICAL_LEV 128
```

The L128 vertical coordinate file inspected locally was:

```text
inputdata/atm/scream/init/vertical_coordinates_L128_20220927.nc
```

Under a standard-column approximation (`ps = 100000 Pa`, scale height
`8000 m`), approximate layer thicknesses in cloud-forming tropospheric regions
are:

| Region | Approximate levels | Approximate dz range |
|---|---:|---:|
| 0-1 km, near-surface/low boundary layer | 21 | 24-52 m |
| 1-2 km, low cloud | 13 | 54-127 m |
| 2-3 km, lower troposphere | 6 | 143-203 m |
| 3-4 km | 5 | 212-260 m |
| 3-5 km | 10 | 212-310 m |

This vertical grid is the same across the listed default horizontal MAM4xx
resolutions, so default horizontal refinement does not by itself refine the
vertical grid in cloud-forming layers.

## Estimated nsubmix and Mixing Work

The following table is an estimate, not an observed diagnostic histogram.

Assumptions used:

```text
Vertical grid: MAM4xx L128
dt_phys: working default convention dt_phys = 6 * se_tstep
dt_ndrop: dt_phys / MAC_subcycles
typical max_tinv: 0.005-0.03 s^-1
strong/cloud-base max_tinv: 0.05-0.12 s^-1
```

The local `ne4` case confirms `dt_phys = 3600 s = 6 * 600 s`, but `dt_phys` is
case-specific through `ATM_NCPL`. Created cases should be checked with:

```text
./xmlquery ATM_NCPL
```

| Resolution | Vertical levels | Dynamics time step (s) | Physics time step (s) | MAC/ACI subcycles | `dt_ndrop` per ACI call (s) | Estimated `max_tinv` typical (s^-1) | Estimated `max_tinv` strong/cloud-base (s^-1) | `nsubmix` typical | `nsubmix` strong/cloud-base | Mixing time step typical (s) | Mixing time step strong/cloud-base (s) | Explicit mixing loop iterations per physics step, typical | Explicit mixing loop iterations per physics step, strong/cloud-base |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ne4np4` | 128 | 600 | 3600 | 24 | 150 | 0.005-0.03 | 0.05-0.12 | 1-6 | 9-21 | 25-150 | 7-17 | 24-144 | 216-504 |
| `ne30np4` | 128 | 300 | 1800 | 12 | 150 | 0.005-0.03 | 0.05-0.12 | 1-6 | 9-21 | 25-150 | 7-17 | 12-72 | 108-252 |
| `ne120np4` | 128 | 75 | 450 | 6 | 75 | 0.005-0.03 | 0.05-0.12 | 1-3 | 5-11 | 25-75 | 7-15 | 6-18 | 30-66 |
| `ne256np4` | 128 | 33.33 | 200 | 6 | 33.33 | 0.005-0.03 | 0.05-0.12 | 1-2 | 2-5 | 17-33 | 7-17 | 6-12 | 12-30 |
| `ne512np4` | 128 | 16.67 | 100 | 2 | 50 | 0.005-0.03 | 0.05-0.12 | 1-2 | 3-7 | 25-50 | 7-17 | 2-4 | 6-14 |
| `ne1024np4` | 128 | 8.33 | 50 | 1 | 50 | 0.005-0.03 | 0.05-0.12 | 1-2 | 3-7 | 25-50 | 7-17 | 1-2 | 3-7 |

The high-resolution rows do not explode in this estimate because the default
physics and/or process timesteps are shorter. However, unusually large
`max_tinv` values can still raise `nsubmix`, and the C++ implementation does
not cap it.

## wtke and Why It Could Be High

In EAMxx/MAM4xx `ndrop`, `wtke` is not prognosed by `ndrop`. It is derived from
the MAM ACI input `wsub`:

```text
wtke(k) = max(wsub(k), 0.1 m/s)
```

That hardcoded `0.1 m/s` floor is separate from the MAM ACI parameter
`wsubmin`, whose default is `0.001 m/s`. Therefore, widespread `wtke = 0.1 m/s`
does not indicate strong turbulence; it indicates the `ndrop` floor.

There are three distinct ranges to keep separate:

| Quantity | Range in current default EAMxx SHOC-to-MAM path | Where the bound comes from |
|---|---:|---|
| Raw `wsub` computed by MAM ACI | about `0.016-5.77 m/s` | Indirectly from SHOC `tke` bounds and default `w2tune = 1` |
| Raw `wsub` allowed by MAM ACI itself | `wsubmin` to no explicit upper cap | MAM ACI applies only `max(wsubmin, wsub)` to liquid-activation `wsub` |
| Effective `wtke` used inside `ndrop` | about `0.1-5.77 m/s` in the current default SHOC path | `ndrop` applies `wtke = max(wsub, 0.1)` |

The finite `5.77 m/s` upper value is not a MAM ACI cap. It follows from the
current SHOC producer. SHOC caps `tke` at `50 m2/s2`, and with default
`w2tune = 1` diagnoses:

```text
w_variance = (2/3) * tke
wsub ~= sqrt(w_variance)
max wsub ~= sqrt((2/3) * 50) = 5.77 m/s
```

Similarly, SHOC's minimum `tke = 0.0004 m2/s2` gives:

```text
min wsub ~= sqrt((2/3) * 0.0004) = 0.016 m/s
```

which is above the default MAM ACI `wsubmin = 0.001 m/s`. If `shoc_1p5tke =
true`, SHOC sets `w_variance = 0`, MAM ACI raw `wsub` falls to `wsubmin`, and
`ndrop` still uses `wtke = 0.1 m/s`.

For EAMxx/SHOC, the path is:

```text
SHOC tke -> SHOC w_variance -> MAM wsub -> ndrop wtke
```

SHOC diagnoses:

```text
w_variance = w2tune * (2/3) * tke
```

MAM then interpolates this to interfaces and computes approximately:

```text
wsub ~= sqrt(averaged w_variance)
```

Possible causes of systematically high `wtke` include:

1. The `ndrop` hard floor of `0.1 m/s`.
2. High SHOC TKE from strong shear, buoyancy production, cloud-top cooling, or
   boundary-layer turbulence.
3. Increased `w2tune`, which scales `w_variance` linearly and `wsub` like
   `sqrt(w2tune)`.
4. The SHOC TKE cap is large enough to allow high `wsub`; with `tke = 50 m2/s2`
   and `w2tune = 1`, `wsub` can approach `sqrt(33) ~= 5.8 m/s`.
5. Incorrect field semantics or units for `w_variance`.
6. Process ordering or stale field values if `mam4_aci` is moved relative to
   `shoc`.
7. Very loose ACI interval checks for `w_variance` in the current interface.

## Comparison With E3SMv3 CLUBB MAM wsub

E3SMv3 default EAM physics uses:

```text
-nlev 80 -clubb_sgs -microphys p3
```

In Fortran EAM/MAM, `microp_aero_wsub_scheme = 1` and
`microp_aero_wsubmin = 0.001` by default.

For CLUBB, `microp_aero.F90` reads CLUBB vertical velocity variance:

```text
WP2_nadv
```

then computes:

```text
tke = (3/2) * WP2_nadv
wsub(k) = sqrt(0.5 * (tke(k) + tke(k+1)) * (2/3))
wsub(k) = min(wsub(k), 10 m/s)
wsub(k) = max(wsubmin, wsub(k))
```

Since the factors `3/2` and `2/3` cancel, this is effectively:

```text
wsub(k) = sqrt(0.5 * [WP2_nadv(k) + WP2_nadv(k+1)])
```

Thus, EAMxx/SHOC and E3SMv3/CLUBB use similar MAM-level formulas: both compute
`wsub` as the square root of a vertical velocity variance. The main difference
is the turbulence scheme providing that variance:

| Aspect | E3SMv3/CLUBB MAM | EAMxx/SHOC MAM |
|---|---|---|
| Source field | `WP2_nadv` | `w_variance` |
| Source meaning | CLUBB vertical velocity variance | SHOC diagnosed vertical velocity variance |
| Vertical location | Interface-like `pverp` field | Midpoint field interpolated to interfaces |
| Formula | `sqrt(mean adjacent WP2)` | `sqrt(mean interpolated w_variance)` |
| Pre-`ndrop` cap | `wsub <= 10 m/s` in Fortran EAM MAM | no direct MAM ACI `wsub` cap; current default SHOC producer indirectly limits raw `wsub` to about `5.77 m/s` |
| Pre-`ndrop` floor | `0.001 m/s` default | `0.001 m/s` default |
| Effective `ndrop` range | `0.1-10 m/s` | about `0.1-5.77 m/s` for current default SHOC input |

Therefore, a systematic difference between EAMxx and E3SMv3 is more likely to
come from SHOC versus CLUBB vertical velocity variance than from the MAM
square-root formula itself.

## SHOC TKE Calculation and Tuning Dependence

SHOC carries `tke` as an updated tracer. The core local SGS TKE update is:

```text
P_buoy  = (g / T0) * wthv_sec           default SHOC path
P_shear = tk * shear_term
D       = Cee / shoc_mix * tke^(3/2)

tke_new = tke_old + dt * (max(0, P_shear + P_buoy) - D)
```

Then it is bounded:

```text
0.0004 m2/s2 <= tke <= 50 m2/s2
```

The most directly relevant SHOC tuning variables are:

| Variable | Default | Main role |
|---|---:|---|
| `length_fac` | 0.5 | Controls `shoc_mix`, which directly changes dissipation `D ~ tke^(3/2) / shoc_mix`. |
| `coeff_km` | 0.1 | Controls momentum diffusivity `tk`, which affects future shear production. |
| `coeff_kh` | 0.1 | Controls heat diffusivity `tkh`, affecting thermodynamic mixing and buoyancy flux indirectly. |
| `lambda_low` | 0.001 | Stability correction lower bound for isotropic timescale. |
| `lambda_high` | 0.08 | Stability correction upper bound for isotropic timescale. |
| `lambda_slope` | 2.65 | Stability correction slope. |
| `lambda_thresh` | 0.02 | Stability threshold for damping. |
| `w2tune` | 1.0 | Does not primarily control TKE, but directly controls diagnosed `w_variance` consumed by MAM ACI. |
| `shoc_1p5tke` | false | Switches to a different 1.5-TKE closure path. |

For the `wtke`/MAM question, `w2tune` is especially important because it is
directly in the mapping from SHOC TKE to the MAM activation velocity variance.

## History of coeff_kh and coeff_km

The local checkout is shallow, so upstream history was checked using a filtered
clone of:

```text
https://github.com/E3SM-Project/E3SM
```

Corrected history:

| Commit | Date | Finding |
|---|---:|---|
| `00fafb066f9` | 2023-10-03 | Introduced runtime options `Ckh` and `Ckm` for SHOC, both defaulting to `0.1`. This replaced hardcoded local values `const Scalar Ckh = 0.1` and `const Scalar Ckm = 0.1` in `shoc_eddy_diffusivities_impl.hpp`. |
| `6636dabafd` | 2025-02-20 | Did not introduce `Ckh`/`Ckm`; it renamed `namelist_defaults_scream.xml` to `namelist_defaults_eamxx.xml` as a 100% file rename. The parent already had `Ckh=0.1` and `Ckm=0.1`. |
| `a5ed75748f` | 2025-04-07 | Renamed namelist keys from `Ckh`/`Ckm` to `coeff_kh`/`coeff_km`; numerical values stayed `0.1`. |
| `48debfbb08` | 2025-04-09 | Added an alternate 1.5-TKE-style eddy diffusivity path using `Ckh*shoc_mix*sqrt(tke)` and `Ckm*shoc_mix*sqrt(tke)`. |
| `36d8330f01` | 2026-03-18 | Reverted eddy diffusivities back to the default SHOC form using `Ckh*isotropy*tke` and `Ckm*isotropy*tke`, retaining the coefficient values. |

No evidence was found that the default numerical values of `coeff_kh` or
`coeff_km` changed from `0.1`; the main changes were introduction as runtime
options, later rename, and changes to the surrounding eddy-diffusivity closure
logic.

## Candidate Integration Strategies

The isolated numerical work changes the integration of the frozen `ndrop`
mixing/activation operator, not the activation physics itself. With mixing,
overlap, geometry, and activation coefficients held fixed, the test problem is
linear. The broader activation calculation is nonlinear when those quantities
are recomputed, and production Forward Euler clipping is a nonlinear,
nonsmooth operation.

| Strategy | Main performance opportunity | Main accuracy or implementation risk | Current evidence |
|---|---|---|---|
| Current subcycled Forward Euler | Cheapest when `nsubmix` is one; no solve setup | Strong mixing requires many sequential substeps and vertical passes; clipping can change the aerosol ledger | Exact-current standalone comparisons pass at safe boundaries; automated remediation suite passes pending independent re-review |
| IMEX-Euler | Replaces stability-driven substeps with fixed scalar tridiagonal solves over many right-hand sides | Explicit activation is lagged for the ACI interval; raw positivity can fail; overhead can lose when `nsubmix` is one or two | CPU correctness and first-order convergence suite passes pending independent re-review; no production or GPU timing |
| Fully implicit Backward Euler | Treats activation and mixing at the new state and removes the explicit mixing restriction | Needs a custom 2-by-2 block-tridiagonal solver; more work and storage; first-order damping; positivity is not guaranteed | Designed but not implemented; blocked behind FE and IMEX gate re-review |
| Hybrid FE/IMEX | Keeps low-cost FE columns and applies IMEX only to the expensive high-`nsubmix` tail | Requires device-side classification/batching, a tested threshold, and two maintained numerical paths | Provisional GPU-oriented concept only |

The provisional GPU preference is hybrid FE/IMEX if expensive columns can be
batched efficiently. Otherwise, use IMEX as the leading candidate for the
strong-mixing tail and retain Forward Euler where `nsubmix` is one or two.
Backward Euler is primarily a fallback if IMEX source lag or raw-positivity
failures prove material. This ordering is not based on measured GPU timing.

For IMEX, the complete aerosol pair is block lower triangular, not one scalar
tridiagonal matrix. It can nevertheless be solved as two sequential scalar
tridiagonal families: cloud-borne aerosol first, then interstitial aerosol
after adding the cloud-to-interstitial return. EKAT, the E3SM Kokkos
Application Toolkit, already supplies CPU Thomas and GPU-oriented cyclic-
reduction solvers for scalar tridiagonal matrices with multiple right-hand
sides. An `ndrop` implementation would still need custom source evaluation,
matrix/right-hand-side assembly, solve ordering, zero-cloud handling,
positivity checks, and ledger diagnostics.

SUNDIALS ARKODE and PETSc `TSARKIMEX` are possible reference or later-stage
integration frameworks, but neither removes that `ndrop`-specific work. For the
current fixed-step first-order method, a small custom IMEX update using EKAT is
the lower-complexity production concept. SUNDIALS becomes more attractive for
an independent oracle, adaptive or higher-order IMEX, or future nonlinear
implicit physics; dependency, Kokkos-version, device-memory, and column-batch
performance would need separate qualification.

Stability and accuracy must be evaluated separately. Several short Forward
Euler steps can be more accurate than one first-order implicit step over the
same ACI interval. Any performance comparison used for a method decision must
therefore include matched-error results as well as equal-interval timings.
Backward Euler is the required implicit baseline; if its damping or first-order
accuracy is inadequate, the first higher-order candidates are a two-stage
second-order L-stable SDIRK2/ESDIRK2 or TR-BDF2 method. BDF2 is less convenient
because other split physics modifies aerosol between ACI calls, and Crank-
Nicolson is a poor first choice because stiff modes can oscillate. Modified
Patankar methods are a separate positivity-and-conservation research track.
RK4 is not part of the current screening plan.

## Recommended Diagnostics

### Existing diagnostics for the current five-day runs

| Priority | Diagnostic | Cadence or summary |
|---|---|---|
| Required | Configuration and provenance listed in the current investigation snapshot | Once per run |
| Required | `EAMxx::run`, `EAMxx::mac_aero_mic::run`, TMS, SHOC, cloud-fraction, MAM4xx-ACI, P3, and EAMxx I/O timers | End of run; retain total, calls, time per call or physics step, and rank variability when available |
| Required context | `eddy_diff_heat`, `w_variance` | Hourly `instant`, identical stream in every case |
| Recommended context | `cldfrac_liq`, `cldfrac_liq_prev` | Hourly `instant`, identical stream in every case |

`eddy_diff_heat` is a real upstream timestep driver:

```text
eddy_diff_heat(midpoint)
  -> interpolated interface kvh_int
  -> bounded ndrop eddy_diff
  -> eddy_diff_kp/km
  -> tinv
  -> dtmix and nsubmix
```

At diagnosed cloud base, `eddy_diff=wtke/zs` replaces the bounded SHOC value,
which is why `w_variance` and liquid cloud fraction are useful additional
context. Hourly snapshots remain correlational; they do not expose all ACI
calls or the limiting level.

### Highest-priority instrumentation

| Priority | Proposed measurement | Performance question answered |
|---|---|---|
| Essential | Fenced timer/profiling region around `MAMAci::call_function_dropmixnuc` | How much ACI time is in the wrapper plus `dropmixnuc`? |
| Essential | Device-aggregated actual `nsubmix`: count, sum, maximum, mean, and histogram | Does runtime scale with explicit-loop trip count, and is a high-cost tail present? |
| Essential derived | `Wmix=sum(nsubmix*nactive)` and optionally `51*Wmix` | Does measured region time scale with expected explicit mixing work? |
| Essential validity | Matched instrumentation-off/on timing | Do the counters and reductions perturb GPU performance? |
| High value | `max_tinv`, empty-reduction flag, and limiting-level index | Why did `nsubmix` change, and is level zero incorrectly controlling it? |
| High value, short segment | GPU occupancy, achieved bandwidth, kernel time, and source/program-counter samples | Are barriers, memory traffic, occupancy, or column imbalance limiting performance? |

`dtmix` need not be stored for performance accounting because the final value
is exactly `dt_aci/nsubmix`. Internal `wsub`, `wtke`, interface coefficients,
limiter counts, and paired ledgers remain important second-stage physics and
numerics diagnostics, but they are not required for first-order runtime
attribution.

## Remaining Caveats

The table values are analytic estimates, not run output. They assume default
MAM4xx L128 vertical resolution and a working default timestep convention
`dt_phys = 6 * se_tstep`. Created cases can override physics coupling frequency
through `ATM_NCPL`, so case-specific timing should be verified from the case
XML before making performance claims.

No conclusion here should be interpreted as a validated performance profile.
The immediate next step is the matched five-day `ne256` matrix using the
existing timer/context contract above. A production-quality attribution still
requires actual `nsubmix`, the scoped `dropmixnuc` timer, and preferably
`max_tinv` with its limiting level. The level-zero boundary defect must be
fixed, or an equivalent safe boundary demonstrated, before interpreting those
adaptive-work measurements.

## Revision History

| Date | Author / Agent | Change Summary |
|---|---|---|
| 2026-08-03 | Codex | Added the current FE, IMEX, Backward Euler, and hybrid comparison; clarified the frozen problem's linear structure, EKAT and external-library roles, matched-error requirement, GPU implications, and deferred higher-order candidates. |
| 2026-07-21 | User / Codex | Corrected the disabled-mixing description to distinguish vertical transport from activation, cover aerosol number as well as mass, and explain the unwritten interstitial destination for every tested substep count. |
| 2026-07-21 | User / Codex | Clarified that Fortran EAMv3 contains the same adaptive explicit mixing subloop and summarized its execution and boundary differences from EAMxx. |
| 2026-07-21 | User / Codex | Verified that E3SM v3.0.0 Fortran always applies both sides of aerosol activation and does not contain the EAMxx disabled-mixing switch defect. |
| 2026-07-21 | User / Codex | Rewrote the `enable_aero_vertical_mix` section in plain language and defined the same-level interstitial/cloud-borne transfers during cloud growth, decay, and disappearance. |
| 2026-07-15 | User / Codex | Consolidated thread findings, current revisions, confirmed defects, masked production symptoms, the five-day `ne256` run matrix, existing-output contract, and prioritized GPU instrumentation. |
| 2026-07-14 | User / Codex | Corrected the as-executed Kokkos reduction behavior, one-substep estimates, field-update map, and confirmed-defect references. |
| 2026-06-16 | Codex | Initial `ndrop` ACI cost analysis. |
