#include "numsolve/diagnostics.hpp"

#include <algorithm>
#include <cmath>

namespace numsolve {

PivotDiagnostics computePivotDiagnostics(const FactorizationResult& fac) {
    PivotDiagnostics diag;
    diag.smallestPivotMagnitude = -1.0;
    diag.worstCandidateRatio = 1.0;

    for (const TraceStep& step : fac.trace) {
        if (step.kind == StepKind::Swap) {
            ++diag.rowSwaps;
        } else if (step.kind == StepKind::Pivot && step.pivot.pivotFound) {
            double usedMagnitude = 0.0;
            for (const PivotCandidate& c : step.pivot.candidates) {
                if (c.row == step.pivot.selectedRow) {
                    usedMagnitude = std::fabs(c.value);
                    break;
                }
            }
            if (diag.smallestPivotMagnitude < 0.0 || usedMagnitude < diag.smallestPivotMagnitude) {
                diag.smallestPivotMagnitude = usedMagnitude;
            }
            if (usedMagnitude > 0.0) {
                double ratio = std::fabs(step.pivot.candidateBeforeSwap) / usedMagnitude;
                diag.worstCandidateRatio = std::min(diag.worstCandidateRatio, ratio);
            }
        }
    }

    if (diag.smallestPivotMagnitude < 0.0) diag.smallestPivotMagnitude = 0.0;
    return diag;
}

GrowthDiagnostics computeGrowthDiagnostics(const FactorizationResult& fac) {
    GrowthDiagnostics diag;
    diag.maxOriginal = fac.maxElementOriginal;
    diag.maxDuringElimination = fac.maxElementDuringElimination;
    diag.growthFactor = (diag.maxOriginal > 0.0) ? diag.maxDuringElimination / diag.maxOriginal : 0.0;
    return diag;
}

FactorizationCheck verifyFactorization(const Matrix& a, const FactorizationResult& fac) {
    Matrix permutedA(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) permutedA(i, j) = a(fac.permutation[i], j);
    }
    Matrix lu = multiply(fac.L, fac.U);
    FactorizationCheck check;
    check.normOfPaMinusLu = normInf(subtract(permutedA, lu));
    return check;
}

ResidualDiagnostics computeResidual(const Matrix& a, const Vector& b, const Vector& x) {
    ResidualDiagnostics diag;
    Vector ax = multiply(a, x);
    Vector residual(ax.size());
    for (std::size_t i = 0; i < ax.size(); ++i) residual[i] = b[i] - ax[i];

    diag.absoluteResidualInf = normInf(residual);
    double denom = normInf(a) * normInf(x) + normInf(b);
    diag.relativeResidual = (denom > 0.0) ? diag.absoluteResidualInf / denom : diag.absoluteResidualInf;
    diag.available = true;
    return diag;
}

}  // namespace numsolve
