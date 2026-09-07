#ifndef NUMSOLVE_TRACE_HPP
#define NUMSOLVE_TRACE_HPP

#include <cstddef>
#include <vector>

namespace numsolve {

// Every trace type below mirrors a step a student would perform by hand,
// using the standard textbook elementary-row-operation notation:
//   Ri <-> Rj              (RowSwapStep)
//   Ri <- Ri - m*Rj         (RowReplacementStep)
// Computation never prints directly; it only records these steps, and a
// renderer turns them into text/JSON later (see report.hpp). This keeps
// the engine testable independently of any presentation format.

struct PivotCandidate {
    std::size_t row;
    double value;
};

// Describes the search for a pivot in one column: which rows were
// candidates, which one was selected (largest magnitude = partial
// pivoting), and what the *unpivoted* candidate would have been -- used
// later to explain "partial pivoting avoided a near-zero divide" instead
// of just labeling a swap as dangerous.
struct PivotStep {
    std::size_t column;
    std::size_t pivotRow;
    std::vector<PivotCandidate> candidates;
    std::size_t selectedRow;
    double candidateBeforeSwap;
    bool pivotFound;
};

struct RowSwapStep {
    std::size_t row1;
    std::size_t row2;
};

struct RowReplacementStep {
    std::size_t target;
    std::size_t source;
    double multiplier;  // target <- target - multiplier * source
};

enum class StepKind { Pivot, Swap, Replace };

struct TraceStep {
    StepKind kind;
    PivotStep pivot{};
    RowSwapStep swapStep{};
    RowReplacementStep replace{};
};

using EliminationTrace = std::vector<TraceStep>;

}  // namespace numsolve

#endif
