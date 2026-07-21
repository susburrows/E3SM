# `ndrop` Evidence Driver

`ndrop_evidence_driver.cpp` is the standalone driver used for the 2026-07-14
boundary and disabled-mixing investigation. It compiles the repository's
unmodified `externals/mam4xx/src/mam4xx/ndrop.hpp` and is intentionally not part
of the E3SM or MAM4xx production build.

The reported run used:

- E3SM `12e44bb469efb0d6c0cb94c7b796a86967883ff5`;
- MAM4xx `1b02f2720a50d7a2d66a6975f7fefe6c7aebd2db`;
- double-precision `haero::Real`;
- Kokkos 4.5.1 Serial with debug mode;
- AppleClang on arm64 macOS; and
- `KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK`, ASan, and UBSan.

## Configuration

Use a short, unsynced build path. Given existing EKAT and HAERO installation
prefixes:

```bash
E3SM_ROOT=/path/to/E3SM
BUILD=/private/tmp/ndrop-ev
EKAT_PREFIX=/path/to/ekat/install
HAERO_PREFIX=/path/to/haero/install

cmake -S "$E3SM_ROOT/externals/mam4xx" -B "$BUILD/mb128" \
  -DNUM_VERTICAL_LEVELS=128 \
  -DENABLE_TESTS=OFF \
  -DENABLE_SKYWALKER=OFF \
  -DMAM4XX_HAERO_DIR="$HAERO_PREFIX" \
  -DEkat_DIR="$EKAT_PREFIX/lib/cmake/ekat" \
  -DKokkos_DIR="$EKAT_PREFIX/lib/cmake/Kokkos" \
  -Dhaero_DIR="$HAERO_PREFIX/lib/cmake/haero" \
  -DCMAKE_BUILD_TYPE=Debug
```

Compile the driver with the same compiler and architecture used by the EKAT
and HAERO libraries. The command used locally was equivalent to:

```bash
clang++ -std=gnu++17 -g -O0 \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -DEKAT_DEFAULT_BFB -DEKAT_HAS_KOKKOS -DKOKKOS_DEPENDENCE \
  -DKOKKOS_ENABLE_DEBUG_BOUNDS_CHECK \
  -I"$BUILD/mb128" \
  -I"$E3SM_ROOT/externals/mam4xx/src" \
  -I"$BUILD/mb128/src" \
  -isystem "$HAERO_PREFIX/include" \
  -isystem "$EKAT_PREFIX/include/ekat" \
  -isystem "$EKAT_PREFIX/include/kokkos" \
  "$E3SM_ROOT/docs/analysis/evidence/ndrop_evidence_driver.cpp" \
  "$HAERO_PREFIX/lib/libhaero.a" \
  "$EKAT_PREFIX/lib/libekat_kokkosutils.a" \
  "$EKAT_PREFIX/lib/libekat_core.a" \
  "$EKAT_PREFIX/lib/libkokkoscontainers.a" \
  "$EKAT_PREFIX/lib/libkokkoscore.a" \
  "$EKAT_PREFIX/lib/libkokkossimd.a" \
  -o "$BUILD/ndrop_evidence"
```

Library names and additional link dependencies can vary with the local EKAT
configuration. Reusing the flags and link libraries from a configured
`mam4_ndrop_unit_tests` target is the most reliable approach.

## Boundary Matrix

Run each expected abort in a separate process:

```bash
"$BUILD/ndrop_evidence" boundary -1
"$BUILD/ndrop_evidence" boundary 0
"$BUILD/ndrop_evidence" boundary 1
"$BUILD/ndrop_evidence" boundary 6
"$BUILD/ndrop_evidence" boundary 127
"$BUILD/ndrop_evidence" boundary 128
```

The investigated revision aborts with a Kokkos index `-1` for `-1` and `0`.
Values `1`, `6`, and `127` return, while `128` returns without rejecting the
invalid upper input.

The driver's host-side substep reporter intentionally returns early for
`top_lev<=0`. Confirm that `LAUNCH_UPDATE top_lev=0` appears before the bounds
error; this ensures the exact MAM4xx function, not replicated reporter
arithmetic, produced the violation.

To test dependence on the interface above the active range:

```bash
"$BUILD/ndrop_evidence" boundary 1 0
"$BUILD/ndrop_evidence" boundary 1 2
"$BUILD/ndrop_evidence" boundary 6 0
"$BUILD/ndrop_evidence" boundary 6 2
```

## Ledger Matrix

The arguments are `ENABLED`, interior mixing coefficient, and destination
buffer scale:

```bash
"$BUILD/ndrop_evidence" ledger 1 0.01 1
"$BUILD/ndrop_evidence" ledger 0 0.01 1
"$BUILD/ndrop_evidence" ledger 1 0.5 1
"$BUILD/ndrop_evidence" ledger 0 0.5 1
"$BUILD/ndrop_evidence" ledger 1 1.0 1
"$BUILD/ndrop_evidence" ledger 0 1.0 1
```

These coefficients produce one, two, and three explicit substeps. Repeat with
destination scales 0 and 7 to test dependence on the unwritten ping-pong
buffer. The focused brief records the expected invariant results and numerical
values.
