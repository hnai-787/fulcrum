#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "numsolve/condition.hpp"
#include "numsolve/elimination.hpp"

using namespace numsolve;

TEST_CASE("condition estimate is exact for a diagonal matrix", "[condition]") {
    // A = diag(1, 1000). ||A||_1 = 1000, ||A^-1||_1 = max(1, 0.001) = 1,
    // so kappa_1(A) = 1000 exactly -- a known, hand-computable reference
    // value the Hager/Higham estimator should reproduce closely.
    Matrix a(2, 2);
    a(0, 0) = 1.0;
    a(1, 1) = 1000.0;

    FactorizationResult fac = factorize(a);
    ConditionEstimate est = estimateCondition(a, fac);

    REQUIRE_FALSE(est.singular);
    REQUIRE(std::fabs(est.conditionEstimate - 1000.0) / 1000.0 < 0.01);
}

TEST_CASE("condition estimate is 1.0 for the identity matrix", "[condition]") {
    Matrix a = Matrix::identity(4);
    FactorizationResult fac = factorize(a);
    ConditionEstimate est = estimateCondition(a, fac);
    REQUIRE(std::fabs(est.conditionEstimate - 1.0) < 1e-9);
}

TEST_CASE("singular matrix is reported as singular, not given a finite estimate", "[condition]") {
    Matrix a(2, 2);
    a(0, 0) = 1; a(0, 1) = 1;
    a(1, 0) = 2; a(1, 1) = 2;  // row 2 = 2 * row 1 -> singular

    FactorizationResult fac = factorize(a);
    ConditionEstimate est = estimateCondition(a, fac);
    REQUIRE(est.singular);
}

TEST_CASE("Hilbert matrix H4 is reported as severely ill-conditioned", "[condition]") {
    // The classic textbook example of a small, severely ill-conditioned
    // matrix (true kappa_1(H4) is on the order of 1e4-1e5).
    constexpr int n = 4;
    Matrix h(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) h(i, j) = 1.0 / (i + j + 1);

    FactorizationResult fac = factorize(h);
    ConditionEstimate est = estimateCondition(h, fac);
    REQUIRE_FALSE(est.singular);
    REQUIRE(est.conditionEstimate > 1000.0);
}

TEST_CASE("rcond is the reciprocal of the condition estimate", "[condition]") {
    Matrix a(2, 2);
    a(0, 0) = 1.0;
    a(1, 1) = 100.0;
    FactorizationResult fac = factorize(a);
    ConditionEstimate est = estimateCondition(a, fac);
    REQUIRE(std::fabs(est.reciprocalEstimate * est.conditionEstimate - 1.0) < 1e-9);
}
