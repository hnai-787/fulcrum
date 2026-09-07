#include "numsolve/linalg.hpp"

#include <stdexcept>

namespace numsolve {

Vector solveLower(const Matrix& l, const Vector& rhs) {
    std::size_t n = l.rows();
    Vector x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double value = rhs[i];
        for (std::size_t j = 0; j < i; ++j) value -= l(i, j) * x[j];
        x[i] = value / l(i, i);
    }
    return x;
}

Vector solveUpper(const Matrix& u, const Vector& rhs) {
    std::size_t n = u.rows();
    Vector x(n, 0.0);
    for (std::size_t ii = 0; ii < n; ++ii) {
        std::size_t i = n - 1 - ii;
        double value = rhs[i];
        for (std::size_t j = i + 1; j < n; ++j) value -= u(i, j) * x[j];
        x[i] = value / u(i, i);
    }
    return x;
}

Vector permute(const Vector& v, const std::vector<std::size_t>& permutation) {
    Vector result(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) result[i] = v[permutation[i]];
    return result;
}

Vector unpermute(const Vector& v, const std::vector<std::size_t>& permutation) {
    Vector result(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) result[permutation[i]] = v[i];
    return result;
}

Vector solveWithFactorization(const FactorizationResult& fac, const Vector& rhs) {
    Vector pb = permute(rhs, fac.permutation);
    Vector y = solveLower(fac.L, pb);
    return solveUpper(fac.U, y);
}

Vector solveTransposeWithFactorization(const FactorizationResult& fac, const Vector& rhs) {
    Matrix ut = transpose(fac.U);
    Matrix lt = transpose(fac.L);
    Vector y1 = solveLower(ut, rhs);
    Vector y2 = solveUpper(lt, y1);
    return unpermute(y2, fac.permutation);
}

}  // namespace numsolve
