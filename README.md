# Gaussian Elimination Solver

## Course Information

| Field | Details |
|---|---|
| Course | Linear Algebra |
| Semester | Semester 2 — Spring 2024 |
| University | Air University, Islamabad |
| Student | Hussain Ali (232095) |

## Overview

A C++ program that solves systems of linear equations (`Ax = b`) using
Gaussian elimination with partial pivoting, and correctly distinguishes
unique, infinite, and no-solution cases.

## Problem Statement

Implement Gaussian elimination robustly enough to handle all three solution
classes a linear system can have, including the numerical edge case of
near-zero pivots.

## Objectives

- Implement partial pivoting to improve numerical stability.
- Detect and report unique-solution, infinite-solution, and no-solution cases.
- Show the augmented matrix before and after elimination for transparency.

## Tools and Technologies

- C++ (`<cmath>`, `<iomanip>`, `<iostream>`, `<vector>`)

## Features

- Partial pivoting for numerical stability.
- `EPSILON = 1e-9` tolerance for near-zero comparisons.
- Clear before/after matrix printout.
- Three sample datasets covering each solution class.

## Methodology

1. Implement elimination with partial pivoting in `src/gaussian_elimination_solver.cpp`.
2. Classify the result (unique / infinite / none) based on the reduced matrix.
3. Validate against three hand-built systems, one per solution class.
4. Keep earlier iterations in `original/` for comparison.

## Repository Structure

```text
gaussian-elimination-solver/
  README.md
  PROJECT_NOTES.md
  src/
    gaussian_elimination_solver.cpp
  original/
    general_gaussian_elimination_solver.cpp
    fixed_6x6_consistency_checker.cpp
  sample-data/
    unique-solution-input.txt
    infinite-solutions-input.txt
    no-solution-input.txt
  screenshots/
  project.yaml
```

## Setup Instructions

```bash
g++ src/gaussian_elimination_solver.cpp -o solver
```

## Usage

```bash
./solver < sample-data/unique-solution-input.txt
```

## How to Review

1. Start with this README.
2. Read `sample-data/*.txt` for the three input systems.
3. Run the solver against each and compare to `screenshots/` for expected output.

## Screenshots

See `screenshots/` — 3 screenshots, one per solution class.

## Results

`sample-data/unique-solution-input.txt` is a real 3×3 system
`[[2,1,-1],[-3,-1,2],[-2,1,2]] x = [8,-11,-3]`; the other two files cover
the infinite- and no-solution cases. All three are verified against the
program's screenshotted output.

## Limitations

- Dense matrices only, no sparse-matrix optimization.
- No complex-number support.

## Future Enhancements

- LU decomposition as an alternative solving strategy.
- A small test harness instead of manual sample-data verification.

## Safety and Privacy

No secrets, credentials, or private data are involved — this is a pure
numerical-methods exercise.

## Ethical Notice

Academic coursework; no ethical concerns apply.
