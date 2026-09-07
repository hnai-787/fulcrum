#ifndef NUMSOLVE_DIAGNOSTICS_HPP
#define NUMSOLVE_DIAGNOSTICS_HPP

#include <optional>

#include "numsolve/condition.hpp"
#include "numsolve/elimination.hpp"

namespace numsolve {

struct PivotDiagnostics {
    std::size_t rowSwaps = 0;
    double smallestPivotMagnitude = 0.0;
    // min over every pivot step of |candidate before swap| / |pivot actually
    // used|; 1.0 if no step ever needed a smaller-than-selected candidate.
    // A very small ratio means partial pivoting rescued the computation
    // from an almost-zero divide -- see README "Design decisions".
    double worstCandidateRatio = 1.0;
};

struct GrowthDiagnostics {
    double maxOriginal = 0.0;
    double maxDuringElimination = 0.0;
    double growthFactor = 0.0;
};

struct ResidualDiagnostics {
    bool available = false;
    double absoluteResidualInf = 0.0;
    double relativeResidual = 0.0;
};

struct FactorizationCheck {
    double normOfPaMinusLu = 0.0;
};

struct Diagnostics {
    PivotDiagnostics pivot;
    GrowthDiagnostics growth;
    ConditionEstimate condition;
    ResidualDiagnostics residual;
    FactorizationCheck factorizationCheck;
};

PivotDiagnostics computePivotDiagnostics(const FactorizationResult& fac);
GrowthDiagnostics computeGrowthDiagnostics(const FactorizationResult& fac);
FactorizationCheck verifyFactorization(const Matrix& a, const FactorizationResult& fac);
ResidualDiagnostics computeResidual(const Matrix& a, const Vector& b, const Vector& x);

}  // namespace numsolve

#endif
