#include "numsolve/condition.hpp"

#include <cmath>

#include "numsolve/linalg.hpp"

namespace numsolve {

ConditionEstimate estimateCondition(const Matrix& a, const FactorizationResult& fac) {
    ConditionEstimate result;
    result.matrixNorm1 = norm1(a);

    if (fac.rank < fac.n) {
        result.singular = true;
        return result;
    }

    std::size_t n = fac.n;
    Vector x(n, 1.0 / static_cast<double>(n));
    double estimate = 0.0;
    long prevJ = -1;
    constexpr std::size_t maxIterations = 5;

    for (std::size_t iter = 0; iter < maxIterations; ++iter) {
        Vector y = solveWithFactorization(fac, x);
        estimate = norm1(y);

        Vector s(n);
        for (std::size_t i = 0; i < n; ++i) s[i] = (y[i] >= 0.0) ? 1.0 : -1.0;

        Vector z = solveTransposeWithFactorization(fac, s);

        std::size_t j = 0;
        double best = std::fabs(z[0]);
        for (std::size_t i = 1; i < n; ++i) {
            if (std::fabs(z[i]) > best) {
                best = std::fabs(z[i]);
                j = i;
            }
        }

        result.iterations = iter + 1;
        if (prevJ >= 0 && best <= std::fabs(z[static_cast<std::size_t>(prevJ)])) {
            break;
        }

        x.assign(n, 0.0);
        x[j] = 1.0;
        prevJ = static_cast<long>(j);
    }

    result.inverseNorm1Estimate = estimate;
    result.conditionEstimate = result.matrixNorm1 * result.inverseNorm1Estimate;
    result.reciprocalEstimate = (result.conditionEstimate > 0.0) ? 1.0 / result.conditionEstimate : 0.0;
    return result;
}

}  // namespace numsolve
