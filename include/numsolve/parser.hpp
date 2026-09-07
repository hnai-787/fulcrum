#ifndef NUMSOLVE_PARSER_HPP
#define NUMSOLVE_PARSER_HPP

#include <string>

#include "numsolve/matrix.hpp"

namespace numsolve {

struct LinearSystem {
    Matrix a;
    Vector b;
};

struct ParseResult {
    bool ok = false;
    LinearSystem system;
    std::string error;
};

// Plain-text augmented-matrix format:
//   line 1: "n" (square) or "rows cols" (cols = number of variables)
//   next `rows` lines: `cols` coefficients, an optional "|" (ignored), then one b value
// Example (3 variables):
//   3
//   2  1 -1 |  8
//   -3 -1 2 | -11
//   -2  1 2 | -3
ParseResult parseAugmentedText(const std::string& source);

// Canonical JSON: {"A": [[...], ...], "b": [...]}
ParseResult parseJson(const std::string& source);

ParseResult parseFile(const std::string& path);

}  // namespace numsolve

#endif
