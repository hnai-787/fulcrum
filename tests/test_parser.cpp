#include <catch2/catch_test_macros.hpp>

#include "numsolve/parser.hpp"

using namespace numsolve;

TEST_CASE("parses a square augmented-matrix text system", "[parser]") {
    std::string src =
        "3\n"
        "2  1 -1 |  8\n"
        "-3 -1 2 | -11\n"
        "-2  1 2 | -3\n";
    ParseResult result = parseAugmentedText(src);
    REQUIRE(result.ok);
    REQUIRE(result.system.a.rows() == 3);
    REQUIRE(result.system.a.cols() == 3);
    REQUIRE(result.system.a(0, 0) == 2.0);
    REQUIRE(result.system.b[0] == 8.0);
    REQUIRE(result.system.b[2] == -3.0);
}

TEST_CASE("parses without the optional pipe separator", "[parser]") {
    std::string src = "2\n1 0 5\n0 1 6\n";
    ParseResult result = parseAugmentedText(src);
    REQUIRE(result.ok);
    REQUIRE(result.system.b[0] == 5.0);
}

TEST_CASE("rows cols header form is accepted", "[parser]") {
    std::string src = "2 2\n1 0 5\n0 1 6\n";
    ParseResult result = parseAugmentedText(src);
    REQUIRE(result.ok);
    REQUIRE(result.system.a.rows() == 2);
}

TEST_CASE("truncated input fails closed with a clear error", "[parser]") {
    std::string src = "3\n1 2 3 | 4\n";  // declares 3 rows, only gives 1
    ParseResult result = parseAugmentedText(src);
    REQUIRE_FALSE(result.ok);
    REQUIRE_FALSE(result.error.empty());
}

TEST_CASE("non-numeric token fails closed", "[parser]") {
    ParseResult result = parseAugmentedText("2\n1 x | 3\n0 1 | 4\n");
    REQUIRE_FALSE(result.ok);
}

TEST_CASE("parses canonical JSON", "[parser]") {
    std::string src = R"({"A": [[2, 1], [1, -1]], "b": [3, 0]})";
    ParseResult result = parseJson(src);
    REQUIRE(result.ok);
    REQUIRE(result.system.a(0, 0) == 2.0);
    REQUIRE(result.system.b[1] == 0.0);
}

TEST_CASE("malformed JSON fails closed", "[parser]") {
    ParseResult result = parseJson("{ not json");
    REQUIRE_FALSE(result.ok);
}

TEST_CASE("JSON with mismatched row lengths fails closed", "[parser]") {
    std::string src = R"({"A": [[1, 2], [3]], "b": [1, 2]})";
    ParseResult result = parseJson(src);
    REQUIRE_FALSE(result.ok);
}
