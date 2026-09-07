#ifndef NUMSOLVE_MATRIX_HPP
#define NUMSOLVE_MATRIX_HPP

#include <cstddef>
#include <vector>

namespace numsolve {

using Vector = std::vector<double>;

class Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, double fill = 0.0)
        : rows_(rows), cols_(cols), data_(rows, Vector(cols, fill)) {}

    static Matrix identity(std::size_t n) {
        Matrix m(n, n);
        for (std::size_t i = 0; i < n; ++i) m(i, i) = 1.0;
        return m;
    }

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }

    double& operator()(std::size_t r, std::size_t c) { return data_[r][c]; }
    double operator()(std::size_t r, std::size_t c) const { return data_[r][c]; }

    Vector& row(std::size_t r) { return data_[r]; }
    const Vector& row(std::size_t r) const { return data_[r]; }

    void swapRows(std::size_t a, std::size_t b) { std::swap(data_[a], data_[b]); }

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<Vector> data_;
};

// Elementwise A*B (standard matrix product).
Matrix multiply(const Matrix& a, const Matrix& b);
Vector multiply(const Matrix& a, const Vector& x);

// Infinity norm: max absolute row sum.
double normInf(const Matrix& a);
// 1-norm: max absolute column sum.
double norm1(const Matrix& a);
// Infinity norm of a vector (max absolute entry).
double normInf(const Vector& v);
// 1-norm of a vector (sum of absolute entries).
double norm1(const Vector& v);

Matrix subtract(const Matrix& a, const Matrix& b);
Matrix transpose(const Matrix& a);
double maxAbsElement(const Matrix& a);

}  // namespace numsolve

#endif
