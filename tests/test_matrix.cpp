#include <catch2/catch_test_macros.hpp>

#include "numsolve/matrix.hpp"

using namespace numsolve;

TEST_CASE("identity matrix", "[matrix]") {
    Matrix i = Matrix::identity(3);
    for (std::size_t r = 0; r < 3; ++r)
        for (std::size_t c = 0; c < 3; ++c) REQUIRE(i(r, c) == (r == c ? 1.0 : 0.0));
}

TEST_CASE("multiply matrix by vector", "[matrix]") {
    Matrix a(2, 2);
    a(0, 0) = 1; a(0, 1) = 2;
    a(1, 0) = 3; a(1, 1) = 4;
    Vector x = {1, 1};
    Vector y = multiply(a, x);
    REQUIRE(y[0] == 3.0);
    REQUIRE(y[1] == 7.0);
}

TEST_CASE("multiply matrix by matrix", "[matrix]") {
    Matrix a = Matrix::identity(2);
    Matrix b(2, 2);
    b(0, 0) = 5; b(0, 1) = 6; b(1, 0) = 7; b(1, 1) = 8;
    Matrix c = multiply(a, b);
    REQUIRE(c(0, 0) == 5.0);
    REQUIRE(c(1, 1) == 8.0);
}

TEST_CASE("norms", "[matrix]") {
    Matrix a(2, 2);
    a(0, 0) = 1; a(0, 1) = -2;
    a(1, 0) = -3; a(1, 1) = 4;
    // 1-norm: max absolute column sum = max(1+3, 2+4) = 6
    REQUIRE(norm1(a) == 6.0);
    // inf-norm: max absolute row sum = max(1+2, 3+4) = 7
    REQUIRE(normInf(a) == 7.0);
    REQUIRE(maxAbsElement(a) == 4.0);
}

TEST_CASE("transpose", "[matrix]") {
    Matrix a(2, 3);
    a(0, 0) = 1; a(0, 1) = 2; a(0, 2) = 3;
    a(1, 0) = 4; a(1, 1) = 5; a(1, 2) = 6;
    Matrix t = transpose(a);
    REQUIRE(t.rows() == 3);
    REQUIRE(t.cols() == 2);
    REQUIRE(t(2, 0) == 3.0);
    REQUIRE(t(0, 1) == 4.0);
}
