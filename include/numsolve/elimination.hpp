#ifndef NUMSOLVE_ELIMINATION_HPP
#define NUMSOLVE_ELIMINATION_HPP

#include "numsolve/matrix.hpp"
#include "numsolve/trace.hpp"

namespace numsolve {

// A scale-aware "numerical zero" threshold: tau = machineEps * n * ||A||_inf.
// A fixed absolute epsilon (e.g. 1e-9) is wrong in general -- it misjudges
// a perfectly well-conditioned but small-scale matrix (like 1e-20 * I) as
// singular, and can be too permissive for a huge-scale matrix. See README
// "Design decisions" for the worked example and citation.
double scaleAwareTolerance(const Matrix& a);

// The result of factoring A (not the augmented [A|b] system) into
// P*A = L*U via Gaussian elimination with partial pivoting. L and U are
// produced directly from the same multipliers recorded in `trace` --
// there is exactly one elimination process here, not two independent
// "methods" that happen to agree (see README "Design decisions").
struct FactorizationResult {
    Matrix L;
    Matrix U;
    std::vector<std::size_t> permutation;  // permutation[i] = original row now at position i
    EliminationTrace trace;
    std::size_t rank = 0;
    std::size_t n = 0;
    double tolerance = 0.0;
    double maxElementOriginal = 0.0;
    double maxElementDuringElimination = 0.0;
};

// Factors a square matrix A. Throws std::invalid_argument if A is not
// square (rectangular systems are still supported for solving/
// classification via eliminateAugmented, just not for this PLU view --
// see README "Limitations").
FactorizationResult factorize(const Matrix& a, double toleranceOverride = -1.0);

// Applies the exact sequence of row swaps and row replacements recorded
// during factorize() to a right-hand-side vector, exactly as if it had
// been carried along as an extra augmented column throughout elimination.
Vector applyEliminationToVector(const FactorizationResult& fac, const Vector& b);

enum class SolutionStatus { Unique, Infinite, None };

struct ClassificationResult {
    SolutionStatus status;
    std::size_t rankA;
    std::size_t rankAugmented;
    std::size_t variables;
};

// Classifies a system from its factorization and eliminated right-hand
// side using the standard rank comparison:
//   rank(A) < rank([A|b])          -> no solution
//   rank(A) == rank([A|b]) < n     -> infinite solutions
//   rank(A) == rank([A|b]) == n    -> unique solution
ClassificationResult classify(const FactorizationResult& fac, const Vector& eliminatedB, double tolerance);

// Back-substitution using U and the eliminated right-hand side. Only
// meaningful when classify() reports Unique.
Vector backSubstitute(const FactorizationResult& fac, const Vector& eliminatedB);

}  // namespace numsolve

#endif
