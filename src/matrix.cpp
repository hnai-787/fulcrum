#include "numsolve/matrix.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace numsolve {

Matrix multiply(const Matrix& a, const Matrix& b) {
    if (a.cols() != b.rows()) throw std::invalid_argument("multiply: dimension mismatch");
    Matrix result(a.rows(), b.cols());
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t k = 0; k < a.cols(); ++k) {
            double aik = a(i, k);
            if (aik == 0.0) continue;
            for (std::size_t j = 0; j < b.cols(); ++j) {
                result(i, j) += aik * b(k, j);
            }
        }
    }
    return result;
}

Vector multiply(const Matrix& a, const Vector& x) {
    if (a.cols() != x.size()) throw std::invalid_argument("multiply: dimension mismatch");
    Vector result(a.rows(), 0.0);
    for (std::size_t i = 0; i < a.rows(); ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < a.cols(); ++j) sum += a(i, j) * x[j];
        result[i] = sum;
    }
    return result;
}

double normInf(const Matrix& a) {
    double best = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < a.cols(); ++j) sum += std::fabs(a(i, j));
        best = std::max(best, sum);
    }
    return best;
}

double norm1(const Matrix& a) {
    double best = 0.0;
    for (std::size_t j = 0; j < a.cols(); ++j) {
        double sum = 0.0;
        for (std::size_t i = 0; i < a.rows(); ++i) sum += std::fabs(a(i, j));
        best = std::max(best, sum);
    }
    return best;
}

double normInf(const Vector& v) {
    double best = 0.0;
    for (double x : v) best = std::max(best, std::fabs(x));
    return best;
}

double norm1(const Vector& v) {
    double sum = 0.0;
    for (double x : v) sum += std::fabs(x);
    return sum;
}

Matrix subtract(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) throw std::invalid_argument("subtract: dimension mismatch");
    Matrix result(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) result(i, j) = a(i, j) - b(i, j);
    }
    return result;
}

double maxAbsElement(const Matrix& a) {
    double best = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) best = std::max(best, std::fabs(a(i, j)));
    }
    return best;
}

Matrix transpose(const Matrix& a) {
    Matrix result(a.cols(), a.rows());
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) result(j, i) = a(i, j);
    }
    return result;
}

}  // namespace numsolve
