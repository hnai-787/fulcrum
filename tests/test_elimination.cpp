#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "numsolve/diagnostics.hpp"
#include "numsolve/elimination.hpp"
#include "numsolve/linalg.hpp"

using namespace numsolve;

namespace {

Matrix mat(std::initializer_list<std::initializer_list<double>> rows) {
    std::size_t r = rows.size();
    std::size_t c = rows.begin()->size();
    Matrix m(r, c);
    std::size_t i = 0;
    for (const auto& row : rows) {
        std::size_t j = 0;
        for (double v : row) m(i, j++) = v;
        ++i;
    }
    return m;
}

}  // namespace

TEST_CASE("unique solution: textbook 3x3 system", "[elimination]") {
    // 2x + y - z = 8, -3x - y + 2z = -11, -2x + y + 2z = -3
    // Known solution: x=2, y=3, z=-1
    Matrix a = mat({{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}});
    Vector b = {8, -11, -3};

    FactorizationResult fac = factorize(a);
    Vector eb = applyEliminationToVector(fac, b);
    ClassificationResult c = classify(fac, eb, fac.tolerance);
    REQUIRE(c.status == SolutionStatus::Unique);

    Vector x = backSubstitute(fac, eb);
    REQUIRE(std::fabs(x[0] - 2.0) < 1e-9);
    REQUIRE(std::fabs(x[1] - 3.0) < 1e-9);
    REQUIRE(std::fabs(x[2] - (-1.0)) < 1e-9);
}

TEST_CASE("requires pivoting: zero in the natural pivot position", "[elimination]") {
    // Without row exchange the first pivot (0,0) is zero.
    Matrix a = mat({{0, 1}, {1, 1}});
    Vector b = {1, 2};

    FactorizationResult fac = factorize(a);
    Vector eb = applyEliminationToVector(fac, b);
    ClassificationResult c = classify(fac, eb, fac.tolerance);
    REQUIRE(c.status == SolutionStatus::Unique);

    Vector x = backSubstitute(fac, eb);
    // x + y = 2, y = 1 -> x = 1, y = 1
    REQUIRE(std::fabs(x[0] - 1.0) < 1e-9);
    REQUIRE(std::fabs(x[1] - 1.0) < 1e-9);

    // Confirm a swap actually happened.
    bool sawSwap = false;
    for (const auto& step : fac.trace) if (step.kind == StepKind::Swap) sawSwap = true;
    REQUIRE(sawSwap);
}

TEST_CASE("scale-invariance regression: tiny-but-well-conditioned matrix is not misjudged as singular",
          "[elimination]") {
    // A fixed absolute epsilon like the original coursework's 1e-9 would
    // wrongly call every pivot of 1e-20*I "zero". A scale-aware tolerance
    // must not.
    Matrix a = Matrix::identity(3);
    for (std::size_t i = 0; i < 3; ++i) a(i, i) = 1e-20;
    Vector b = {1e-20, 2e-20, 3e-20};

    FactorizationResult fac = factorize(a);
    REQUIRE(fac.rank == 3);
    Vector eb = applyEliminationToVector(fac, b);
    ClassificationResult c = classify(fac, eb, fac.tolerance);
    REQUIRE(c.status == SolutionStatus::Unique);

    Vector x = backSubstitute(fac, eb);
    REQUIRE(std::fabs(x[0] - 1.0) < 1e-6);
    REQUIRE(std::fabs(x[1] - 2.0) < 1e-6);
    REQUIRE(std::fabs(x[2] - 3.0) < 1e-6);
}

TEST_CASE("inconsistent system reports no solution", "[elimination]") {
    // x + y = 1, x + y = 2 -- contradictory
    Matrix a = mat({{1, 1}, {1, 1}});
    Vector b = {1, 2};

    FactorizationResult fac = factorize(a);
    Vector eb = applyEliminationToVector(fac, b);
    ClassificationResult c = classify(fac, eb, fac.tolerance);
    REQUIRE(c.status == SolutionStatus::None);
}

TEST_CASE("rank-deficient consistent system reports infinite solutions", "[elimination]") {
    // x + y = 2, 2x + 2y = 4 -- same equation twice
    Matrix a = mat({{1, 1}, {2, 2}});
    Vector b = {2, 4};

    FactorizationResult fac = factorize(a);
    REQUIRE(fac.rank == 1);
    Vector eb = applyEliminationToVector(fac, b);
    ClassificationResult c = classify(fac, eb, fac.tolerance);
    REQUIRE(c.status == SolutionStatus::Infinite);
}

TEST_CASE("PA = LU factorization identity holds", "[elimination]") {
    Matrix a = mat({{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}});
    FactorizationResult fac = factorize(a);
    FactorizationCheck check = verifyFactorization(a, fac);
    REQUIRE(check.normOfPaMinusLu < 1e-9);
}

TEST_CASE("L is unit lower triangular and U is upper triangular", "[elimination]") {
    Matrix a = mat({{0, 1}, {1, 1}});  // forces a swap
    FactorizationResult fac = factorize(a);
    for (std::size_t i = 0; i < 2; ++i) {
        REQUIRE(fac.L(i, i) == 1.0);
        for (std::size_t j = i + 1; j < 2; ++j) REQUIRE(fac.L(i, j) == 0.0);
        for (std::size_t j = 0; j < i; ++j) REQUIRE(fac.U(i, j) == 0.0);
    }
}

TEST_CASE("residual is tiny for a well-conditioned solved system", "[elimination]") {
    Matrix a = mat({{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}});
    Vector b = {8, -11, -3};
    FactorizationResult fac = factorize(a);
    Vector eb = applyEliminationToVector(fac, b);
    Vector x = backSubstitute(fac, eb);
    ResidualDiagnostics r = computeResidual(a, b, x);
    REQUIRE(r.absoluteResidualInf < 1e-9);
}

TEST_CASE("element growth factor is 1.0 when no growth occurs", "[elimination]") {
    Matrix a = Matrix::identity(3);
    FactorizationResult fac = factorize(a);
    GrowthDiagnostics g = computeGrowthDiagnostics(fac);
    REQUIRE(std::fabs(g.growthFactor - 1.0) < 1e-9);
}

TEST_CASE("pivot diagnostics report a rescue when partial pivoting avoids a tiny candidate", "[elimination]") {
    // Natural (0,0) candidate is tiny; row 2 offers a much larger one.
    Matrix a = mat({{1e-15, 1}, {1, 1}});
    Vector b = {1, 2};
    FactorizationResult fac = factorize(a);
    PivotDiagnostics diag = computePivotDiagnostics(fac);
    REQUIRE(diag.rowSwaps == 1);
    REQUIRE(diag.worstCandidateRatio < 1e-10);
}
