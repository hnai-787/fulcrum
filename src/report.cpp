#include "numsolve/report.hpp"

#include <cmath>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <sstream>

namespace numsolve {

using nlohmann::json;

namespace {

std::string fmt(double v, int precision = 4) {
    std::ostringstream oss;
    if (std::fabs(v) >= 1e5 || (std::fabs(v) > 0 && std::fabs(v) < 1e-4)) {
        oss << std::scientific << std::setprecision(std::max(2, precision - 2));
    } else {
        oss << std::fixed << std::setprecision(precision);
    }
    oss << v;
    return oss.str();
}

std::string signedTerm(double m) { return (m >= 0.0) ? ("- " + fmt(m)) : ("+ " + fmt(-m)); }

}  // namespace

std::string renderMatrix(const Matrix& m, const std::string& label) {
    std::ostringstream oss;
    if (!label.empty()) oss << label << "\n";
    for (std::size_t i = 0; i < m.rows(); ++i) {
        oss << "[ ";
        for (std::size_t j = 0; j < m.cols(); ++j) oss << std::setw(10) << fmt(m(i, j)) << " ";
        oss << "]\n";
    }
    return oss.str();
}

std::string renderVector(const Vector& v, const std::string& label) {
    std::ostringstream oss;
    if (!label.empty()) oss << label << "\n";
    for (double x : v) oss << "[ " << std::setw(10) << fmt(x) << " ]\n";
    return oss.str();
}

std::string renderEliminationTrace(const LinearSystem& system, const FactorizationResult& fac,
                                    const Vector& /*eliminatedB*/) {
    std::ostringstream oss;
    Matrix work = system.a;
    Vector wb = system.b;
    bool pendingSnapshot = false;
    int stepNumber = 0;

    auto printSnapshotIfPending = [&]() {
        if (!pendingSnapshot) return;
        Matrix augmented(work.rows(), work.cols() + 1);
        for (std::size_t r = 0; r < work.rows(); ++r) {
            for (std::size_t c = 0; c < work.cols(); ++c) augmented(r, c) = work(r, c);
            augmented(r, work.cols()) = wb[r];
        }
        oss << renderMatrix(augmented) << "\n";
        pendingSnapshot = false;
    };

    for (const TraceStep& step : fac.trace) {
        if (step.kind == StepKind::Pivot) {
            printSnapshotIfPending();
            ++stepNumber;
            const PivotStep& p = step.pivot;
            oss << "Step " << stepNumber << " -- choose pivot in column " << (p.column + 1) << "\n";
            oss << "Candidates:\n";
            for (const PivotCandidate& c : p.candidates) {
                oss << "  R" << (c.row + 1) << "  " << fmt(c.value);
                if (c.row == p.selectedRow) oss << "  <- largest magnitude";
                oss << "\n";
            }
            if (!p.pivotFound) {
                oss << "No candidate exceeds the numerical tolerance in this column -- "
                       "no pivot here (matrix is rank-deficient in this column).\n\n";
            } else if (p.selectedRow != p.pivotRow) {
                oss << "Partial pivoting selects R" << (p.selectedRow + 1) << ".\n";
                oss << "R" << (p.pivotRow + 1) << " <-> R" << (p.selectedRow + 1) << "\n";
                std::swap(work.row(p.pivotRow), work.row(p.selectedRow));
                std::swap(wb[p.pivotRow], wb[p.selectedRow]);
                pendingSnapshot = true;
            }
        } else if (step.kind == StepKind::Replace) {
            const RowReplacementStep& r = step.replace;
            oss << "m = " << fmt(r.multiplier) << "\n";
            oss << "R" << (r.target + 1) << " <- R" << (r.target + 1) << " " << signedTerm(r.multiplier)
                << "·R" << (r.source + 1) << "\n";
            for (std::size_t c = 0; c < work.cols(); ++c) work(r.target, c) -= r.multiplier * work(r.source, c);
            wb[r.target] -= r.multiplier * wb[r.source];
            pendingSnapshot = true;
        }
    }
    printSnapshotIfPending();

    return oss.str();
}

std::string renderClassification(const ClassificationResult& c) {
    std::ostringstream oss;
    oss << "System classification\n";
    oss << "rank(A)        = " << c.rankA << "\n";
    oss << "rank([A|b])    = " << c.rankAugmented << "\n";
    oss << "variables      = " << c.variables << "\n\n";

    switch (c.status) {
        case SolutionStatus::None:
            oss << "Result: NO SOLUTION -- rank(A) < rank([A|b]) means the equations contradict each other.\n";
            break;
        case SolutionStatus::Infinite:
            oss << "Result: INFINITE SOLUTIONS -- consistent (rank(A) = rank([A|b])) but rank < variables, "
                << "so " << (c.variables - c.rankA) << " variable(s) are free.\n";
            break;
        case SolutionStatus::Unique:
            oss << "Result: UNIQUE SOLUTION -- rank(A) = rank([A|b]) = variables.\n";
            break;
    }
    return oss.str();
}

std::string renderPluView(const FactorizationResult& fac, const FactorizationCheck& check) {
    std::ostringstream oss;
    oss << "PLU interpretation\n";
    oss << "==================\n";
    oss << "The elimination multipliers recorded above populate L directly "
           "(this is the same elimination, not a second algorithm):\n\n";
    oss << renderMatrix(fac.L, "L =");
    oss << "\n" << renderMatrix(fac.U, "U =");
    oss << "\nPermutation (entry i is which original row now sits at position i):\n";
    oss << "[ ";
    for (std::size_t i = 0; i < fac.permutation.size(); ++i) oss << (fac.permutation[i] + 1) << " ";
    oss << "]\n\n";
    oss << "Verification: ||P*A - L*U||_inf = " << fmt(check.normOfPaMinusLu, 3) << "\n";
    return oss.str();
}

std::string renderDiagnostics(const Diagnostics& d) {
    std::ostringstream oss;
    oss << "Numerical diagnostics\n";
    oss << "======================\n";
    oss << "Conditioning (Hager/Higham 1-norm estimate)\n";
    if (d.condition.singular) {
        oss << "  matrix is singular -- condition number is effectively infinite\n";
    } else {
        oss << "  ||A||_1                 " << fmt(d.condition.matrixNorm1, 6) << "\n";
        oss << "  estimated ||A^-1||_1     " << fmt(d.condition.inverseNorm1Estimate, 6) << "\n";
        oss << "  estimated kappa_1(A)     " << fmt(d.condition.conditionEstimate, 6) << "\n";
        oss << "  estimated rcond_1(A)     " << fmt(d.condition.reciprocalEstimate, 6) << "\n";
        oss << "  estimator iterations     " << d.condition.iterations << "\n";
    }
    oss << "\nPivoting\n";
    oss << "  row swaps                " << d.pivot.rowSwaps << "\n";
    oss << "  smallest pivot magnitude " << fmt(d.pivot.smallestPivotMagnitude, 6) << "\n";
    oss << "  worst candidate ratio    " << fmt(d.pivot.worstCandidateRatio, 6)
        << (d.pivot.worstCandidateRatio < 1e-6 ? "  <- partial pivoting rescued a near-zero divide" : "") << "\n";
    oss << "\nElement growth\n";
    oss << "  max |A_original|         " << fmt(d.growth.maxOriginal, 6) << "\n";
    oss << "  max |A| during elimination " << fmt(d.growth.maxDuringElimination, 6) << "\n";
    oss << "  growth factor            " << fmt(d.growth.growthFactor, 6) << "\n";
    if (d.residual.available) {
        oss << "\nSolution verification\n";
        oss << "  ||Ax-b||_inf             " << fmt(d.residual.absoluteResidualInf, 6) << "\n";
        oss << "  relative residual        " << fmt(d.residual.relativeResidual, 6) << "\n";
    }
    return oss.str();
}

std::string renderJson(const LinearSystem& system, const FactorizationResult& fac,
                        const ClassificationResult& classification, const Vector* solution,
                        const Diagnostics& diagnostics) {
    json root;
    root["classification"] = {
        {"status", classification.status == SolutionStatus::Unique     ? "unique"
                   : classification.status == SolutionStatus::Infinite ? "infinite"
                                                                        : "none"},
        {"rank_a", classification.rankA},
        {"rank_augmented", classification.rankAugmented},
        {"variables", classification.variables},
    };
    if (solution) root["solution"] = *solution;

    root["plu"] = {
        {"rank", fac.rank},
        {"factorization_check_inf_norm", diagnostics.factorizationCheck.normOfPaMinusLu},
    };

    root["diagnostics"] = {
        {"condition",
         {
             {"singular", diagnostics.condition.singular},
             {"norm_1_A", diagnostics.condition.matrixNorm1},
             {"estimated_norm_1_A_inverse", diagnostics.condition.inverseNorm1Estimate},
             {"estimated_kappa_1", diagnostics.condition.conditionEstimate},
             {"estimated_rcond_1", diagnostics.condition.reciprocalEstimate},
             {"estimator", "hager-higham-1norm"},
             {"iterations", diagnostics.condition.iterations},
         }},
        {"pivoting",
         {
             {"row_swaps", diagnostics.pivot.rowSwaps},
             {"smallest_pivot_magnitude", diagnostics.pivot.smallestPivotMagnitude},
             {"worst_candidate_ratio", diagnostics.pivot.worstCandidateRatio},
         }},
        {"growth",
         {
             {"max_original", diagnostics.growth.maxOriginal},
             {"max_during_elimination", diagnostics.growth.maxDuringElimination},
             {"growth_factor", diagnostics.growth.growthFactor},
         }},
    };
    if (diagnostics.residual.available) {
        root["diagnostics"]["residual"] = {
            {"absolute_inf_norm", diagnostics.residual.absoluteResidualInf},
            {"relative", diagnostics.residual.relativeResidual},
        };
    }

    (void)system;
    return root.dump(2);
}

}  // namespace numsolve
