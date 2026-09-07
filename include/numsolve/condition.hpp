#ifndef NUMSOLVE_CONDITION_HPP
#define NUMSOLVE_CONDITION_HPP

#include "numsolve/elimination.hpp"

namespace numsolve {

// Estimated (not exact) 1-norm condition number, via the classic
// Hager (1984) / Higham (1988) 1-norm estimator: ||A^-1||_1 is estimated
// through a handful of triangular solves against the already-computed
// P*A=L*U factorization, without ever forming A^-1 explicitly. See
// README "Design decisions" for why an estimate (and why the 1-norm)
// rather than an exact kappa_2 via SVD.
struct ConditionEstimate {
    bool singular = false;
    double matrixNorm1 = 0.0;
    double inverseNorm1Estimate = 0.0;
    double conditionEstimate = 0.0;   // kappa_1(A) estimate
    double reciprocalEstimate = 0.0;  // rcond_1(A) estimate = 1 / kappa_1
    std::size_t iterations = 0;
};

ConditionEstimate estimateCondition(const Matrix& a, const FactorizationResult& fac);

}  // namespace numsolve

#endif
