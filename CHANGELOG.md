# Changelog

All notable changes to this project are documented here.
Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added

### Changed

### Fixed

## [1.0.0] - 2026-09-08

### Added

- Rebuilt as **numsolve**, an explainable numerical-methods CLI, on top of
  the original version solver (preserved unmodified under
  `archive/original/`).
- Structured elimination trace (`PivotStep`/`RowSwapStep`/
  `RowReplacementStep`) rendered in standard textbook row-operation
  notation, instead of only printing before/after matrix snapshots.
- PLU factorization derived from the exact same elimination multipliers
  as the trace (`P*A = L*U`), with a `||P*A - L*U||_inf` verification,
  framed as one process viewed two ways rather than two independent
  methods.
- Hager (1984) / Higham (1988) 1-norm condition number estimator
  (`estimated kappa_1`, `rcond_1`), implemented via triangular solves
  against the factorization -- no BLAS/LAPACK dependency, no explicit
  matrix inversion.
- Pivot diagnostics (row-swap count, smallest pivot used, "worst
  candidate ratio" showing when partial pivoting rescued a near-zero
  divide) and element growth factor.
- Solution residual (`||Ax-b||_inf` and a scaled relative residual).
- Replaced the original fixed `EPSILON = 1e-9` tolerance with a
  scale-aware threshold (`machine_epsilon * n * ||A||_inf`).
- Plain-text augmented-matrix parser, canonical JSON parser, and stdin
  input, replacing the original's interactive-prompt-only input (which
  made it impossible to script or test).
- 28 Catch2 tests, including a scale-invariance regression test and
  condition-estimator tests against hand-computable reference values.

### Fixed

- The scale-aware tolerance initially had a defensive floor of one
  machine epsilon "just in case division by zero"; that floor silently
  reintroduced the exact bug it was meant to fix for matrices at a
  smaller scale (e.g. `1e-20 * I`), caught by the scale-invariance
  regression test before this was ever shipped. Removed the floor; a
  genuinely zero matrix still correctly yields tolerance 0.
