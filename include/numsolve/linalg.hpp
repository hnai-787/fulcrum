#ifndef NUMSOLVE_LINALG_HPP
#define NUMSOLVE_LINALG_HPP

#include "numsolve/elimination.hpp"
#include "numsolve/matrix.hpp"

namespace numsolve {

// General forward substitution: solves L*x = rhs where L is lower
// triangular (diagonal need not be 1 -- used both for the factorization's
// own unit-diagonal L and, transposed, for upper-triangular U).
Vector solveLower(const Matrix& l, const Vector& rhs);

// General back substitution: solves U*x = rhs where U is upper triangular.
Vector solveUpper(const Matrix& u, const Vector& rhs);

Vector permute(const Vector& v, const std::vector<std::size_t>& permutation);
Vector unpermute(const Vector& v, const std::vector<std::size_t>& permutation);

// Solves A*x = rhs using an already-computed P*A = L*U factorization:
// forward-substitute L*y = P*rhs, then back-substitute U*x = y.
Vector solveWithFactorization(const FactorizationResult& fac, const Vector& rhs);

// Solves A^T*z = rhs using the same factorization, via
// A^T = U^T * L^T * P  =>  solve U^T y1 = rhs, L^T y2 = y1, z = P^T y2.
// U^T and L^T are triangular the other way around, so this reuses
// solveLower/solveUpper on transposed (small, cheap-to-copy) matrices
// rather than needing separate "transposed solve" implementations.
Vector solveTransposeWithFactorization(const FactorizationResult& fac, const Vector& rhs);

}  // namespace numsolve

#endif
