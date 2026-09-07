#ifndef NUMSOLVE_REPORT_HPP
#define NUMSOLVE_REPORT_HPP

#include <string>

#include "numsolve/diagnostics.hpp"
#include "numsolve/elimination.hpp"
#include "numsolve/parser.hpp"

namespace numsolve {

std::string renderMatrix(const Matrix& m, const std::string& label = "");
std::string renderVector(const Vector& v, const std::string& label = "");

// The step-by-step "what a student would write by hand" view: every
// pivot selection and row operation, in Ri <-> Rj / Ri <- Ri - m*Rj
// notation, ending in back substitution.
std::string renderEliminationTrace(const LinearSystem& system, const FactorizationResult& fac,
                                    const Vector& eliminatedB);

std::string renderClassification(const ClassificationResult& classification);

// The "same elimination, viewed as a factorization" section: L, U, the
// permutation, and the P*A ≈ L*U identity check.
std::string renderPluView(const FactorizationResult& fac, const FactorizationCheck& check);

std::string renderDiagnostics(const Diagnostics& diagnostics);

std::string renderJson(const LinearSystem& system, const FactorizationResult& fac,
                        const ClassificationResult& classification, const Vector* solution,
                        const Diagnostics& diagnostics);

}  // namespace numsolve

#endif
