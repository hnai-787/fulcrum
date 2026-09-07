# Project Notes

## Source

Migrated from `air-university-cybersecurity-projects/projects/gaussian-elimination-solver`
into this workspace as an independent project on 2026-09-07.

## Cleanup decisions

- The compiled `gaussian-elimination-solver.exe` in the source folder was
  **not** copied (build output, excluded per this workspace's project
  `.gitignore` convention). Rebuild with `g++ src/gaussian_elimination_solver.cpp -o solver`.

## Assumptions

None beyond what's stated in the README.

## Remaining work

None identified.

## 2026-09-08: Rebuilt as numsolve (explainable numerical-methods CLI)

### What changed and why

Raw Gaussian elimination is a solved problem (numpy/LAPACK exist); the
sharper direction agreed on for this project was an *explainable*
numerical-methods tool: full elimination trace in textbook notation, PLU
factorization presented as the same process rather than a second
algorithm, and real numerical diagnostics (conditioning, pivot quality,
element growth, residual) instead of a bare "unique/infinite/none"
verdict. Grounded in a research pass on Hager/Higham condition
estimation, GEPP stability literature, and standard row-operation
notation conventions.

The original interactive-only, fixed-epsilon program was preserved
unmodified under `archive/academic-original/`.

### Key engineering decisions and why

- **Scale-aware tolerance, not a fixed epsilon.** The original
  `EPSILON = 1e-9` is provably wrong for small-scale-but-well-conditioned
  matrices (e.g. `1e-20 * I`, condition number 1, every pivot smaller
  than the epsilon). Replaced with `machine_epsilon * n * ||A||_inf`.
- **One elimination process, two views, not two methods.** PLU
  factorization is derived directly from the same recorded multipliers
  the elimination trace uses -- `L`'s entries are literally the
  multipliers already computed. Verified via `||PA - LU||_inf` and
  `||Ax-b||_inf` as two genuine mathematical invariants, rather than
  misleadingly presenting "Gaussian elimination" and "LU decomposition"
  as independently-agreeing methods (they aren't independent under
  partial pivoting).
- **Hager/Higham 1-norm condition estimator, not an exact SVD-based
  kappa.** Avoids a BLAS/LAPACK dependency entirely; reuses the
  triangular solves already needed elsewhere. Always reported and
  labeled as an estimate.
- **Dropped interactive-prompt-only input.** The original program could
  only be driven via a sequence of `cin >>` prompts, which made it
  impossible to script or unit-test. numsolve reads a file or stdin
  instead -- a deliberate CLI-shape change, not an accidental regression
  (interactive mode was intentionally not preserved; see README
  Limitations).
- **Reused the CMake + vcpkg (nlohmann-json, Catch2) toolchain** already
  set up for the sibling `fwlint` project in this same workspace, rather
  than re-provisioning dependencies from scratch.

### A real bug found and fixed via the test suite

The scale-aware tolerance function's first version included a defensive
floor -- `max(tau, machine_epsilon)` -- added out of an unfounded worry
about tolerance ever being exactly zero. That floor (~2.22e-16) is far
larger than the pivots of a `1e-20 * I` matrix, so it silently
reintroduced the *exact* bug the scale-aware tolerance was built to fix:
the scale-invariance regression test caught `fac.rank == 0` (expected 3)
immediately. Removed the floor entirely -- a genuinely all-zero matrix
still correctly yields tolerance 0 and reports no pivot, no floor needed.
This is now a permanent regression test in `tests/test_elimination.cpp`.

### Verification performed

`cmake --build` and the full test suite (28 tests / 70 assertions, all
passing) were actually run. The CLI was run against: the classic textbook
3×3 unique-solution system (verified against the known hand solution
x=(2,3,-1)), the original coursework's no-solution and
infinite-solutions cases, a 5×5 Hilbert matrix constructed so the true
solution is the all-ones vector (demonstrating, with real numbers, that a
tiny residual does not imply an accurate solution when the matrix is
ill-conditioned), and the rectangular-matrix rejection path. Every number
quoted in the README's worked examples is real, captured CLI output.

### Remaining work / honest limitations

See README "Limitations" and "Future Enhancements" -- notably: square
systems only, no LaTeX export yet, no CI pipeline run (not pushed to
GitHub in this task).
