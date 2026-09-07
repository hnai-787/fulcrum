#include "numsolve/elimination.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace numsolve {

double scaleAwareTolerance(const Matrix& a) {
    // No absolute floor here on purpose: flooring tau at machine epsilon
    // would silently reintroduce the exact bug this function exists to
    // avoid (a well-conditioned but small-scale matrix like 1e-20*I would
    // have every pivot misjudged as "zero"). If the matrix is genuinely
    // all-zero, scale is 0, tau is 0, and the strict `> tolerance` pivot
    // check still correctly finds no pivot -- no floor is needed for that
    // case either.
    std::size_t n = std::max(a.rows(), a.cols());
    double scale = normInf(a);
    return std::numeric_limits<double>::epsilon() * static_cast<double>(std::max<std::size_t>(n, 1)) * scale;
}

FactorizationResult factorize(const Matrix& a, double toleranceOverride) {
    if (a.rows() != a.cols()) {
        throw std::invalid_argument("factorize: matrix must be square (rectangular systems use eliminateAugmented-style classification only)");
    }

    std::size_t n = a.rows();
    FactorizationResult result;
    result.n = n;
    result.U = a;
    result.L = Matrix::identity(n);
    result.permutation.resize(n);
    for (std::size_t i = 0; i < n; ++i) result.permutation[i] = i;

    result.maxElementOriginal = maxAbsElement(a);
    result.maxElementDuringElimination = result.maxElementOriginal;
    result.tolerance = toleranceOverride >= 0.0 ? toleranceOverride : scaleAwareTolerance(a);

    std::size_t pivotRow = 0;
    for (std::size_t col = 0; col < n && pivotRow < n; ++col) {
        PivotStep pivotStep;
        pivotStep.column = col;
        pivotStep.pivotRow = pivotRow;
        pivotStep.candidateBeforeSwap = result.U(pivotRow, col);

        std::size_t bestRow = pivotRow;
        for (std::size_t r = pivotRow; r < n; ++r) {
            pivotStep.candidates.push_back({r, result.U(r, col)});
            if (std::fabs(result.U(r, col)) > std::fabs(result.U(bestRow, col))) bestRow = r;
        }
        pivotStep.selectedRow = bestRow;
        pivotStep.pivotFound = std::fabs(result.U(bestRow, col)) > result.tolerance;
        result.trace.push_back({StepKind::Pivot, pivotStep, {}, {}});

        if (!pivotStep.pivotFound) {
            continue;  // no usable pivot in this column; try the next column at the same pivotRow
        }

        if (bestRow != pivotRow) {
            result.U.swapRows(pivotRow, bestRow);
            std::swap(result.permutation[pivotRow], result.permutation[bestRow]);
            for (std::size_t c = 0; c < pivotRow; ++c) {
                std::swap(result.L(pivotRow, c), result.L(bestRow, c));
            }
            result.trace.push_back({StepKind::Swap, {}, {pivotRow, bestRow}, {}});
        }

        double pivotVal = result.U(pivotRow, col);
        for (std::size_t r = pivotRow + 1; r < n; ++r) {
            double m = result.U(r, col) / pivotVal;
            if (m != 0.0) {
                for (std::size_t c = col; c < n; ++c) {
                    result.U(r, c) -= m * result.U(pivotRow, c);
                }
                result.L(r, pivotRow) = m;
                result.trace.push_back({StepKind::Replace, {}, {}, {r, pivotRow, m}});
            }
            for (std::size_t c = col; c < n; ++c) {
                result.maxElementDuringElimination = std::max(result.maxElementDuringElimination, std::fabs(result.U(r, c)));
            }
        }

        ++pivotRow;
    }

    result.rank = pivotRow;
    return result;
}

Vector applyEliminationToVector(const FactorizationResult& fac, const Vector& b) {
    Vector v = b;
    for (const TraceStep& step : fac.trace) {
        if (step.kind == StepKind::Swap) {
            std::swap(v[step.swapStep.row1], v[step.swapStep.row2]);
        } else if (step.kind == StepKind::Replace) {
            v[step.replace.target] -= step.replace.multiplier * v[step.replace.source];
        }
    }
    return v;
}

ClassificationResult classify(const FactorizationResult& fac, const Vector& eliminatedB, double tolerance) {
    ClassificationResult result;
    result.variables = fac.n;
    result.rankA = fac.rank;
    result.rankAugmented = fac.rank;

    for (std::size_t r = fac.rank; r < fac.n; ++r) {
        if (std::fabs(eliminatedB[r]) > tolerance) {
            result.rankAugmented = fac.rank + 1;
            break;
        }
    }

    if (result.rankAugmented > result.rankA) {
        result.status = SolutionStatus::None;
    } else if (result.rankA < result.variables) {
        result.status = SolutionStatus::Infinite;
    } else {
        result.status = SolutionStatus::Unique;
    }
    return result;
}

Vector backSubstitute(const FactorizationResult& fac, const Vector& eliminatedB) {
    std::size_t n = fac.n;
    Vector x(n, 0.0);
    for (std::size_t ii = 0; ii < n; ++ii) {
        std::size_t r = n - 1 - ii;
        double value = eliminatedB[r];
        for (std::size_t c = r + 1; c < n; ++c) value -= fac.U(r, c) * x[c];
        x[r] = value / fac.U(r, r);
    }
    return x;
}

}  // namespace numsolve
