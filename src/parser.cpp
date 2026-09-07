#include "numsolve/parser.hpp"

#include <cctype>
#include <fstream>
#include <nlohmann/json.hpp>
#include <sstream>

namespace numsolve {

using nlohmann::json;

ParseResult parseAugmentedText(const std::string& source) {
    ParseResult result;
    std::istringstream stream(source);
    std::string firstLine;

    if (!std::getline(stream, firstLine)) {
        result.error = "empty input";
        return result;
    }

    std::istringstream headerStream(firstLine);
    std::vector<long long> headerValues;
    long long value;
    while (headerStream >> value) headerValues.push_back(value);

    std::size_t rows, cols;
    if (headerValues.size() == 1) {
        rows = cols = static_cast<std::size_t>(headerValues[0]);
    } else if (headerValues.size() == 2) {
        rows = static_cast<std::size_t>(headerValues[0]);
        cols = static_cast<std::size_t>(headerValues[1]);
    } else {
        result.error = "first line must contain 1 integer (n) or 2 integers (rows cols)";
        return result;
    }

    if (rows == 0 || cols == 0) {
        result.error = "matrix dimensions must be positive";
        return result;
    }
    if (rows > 200 || cols > 200) {
        result.error = "matrix dimensions too large for this educational solver (limit 200)";
        return result;
    }

    Matrix a(rows, cols);
    Vector b(rows);

    // Read the remaining content as a flat token stream, skipping "|"
    // separators wherever they appear, and pull (cols + 1) numbers per row.
    std::string rest((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
    std::istringstream tokenStream(rest);
    std::string token;

    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            do {
                if (!(tokenStream >> token)) {
                    result.error = "unexpected end of input while reading coefficients";
                    return result;
                }
            } while (token == "|");
            try {
                a(r, c) = std::stod(token);
            } catch (const std::exception&) {
                result.error = "invalid numeric value: \"" + token + "\"";
                return result;
            }
        }
        do {
            if (!(tokenStream >> token)) {
                result.error = "unexpected end of input while reading constants vector";
                return result;
            }
        } while (token == "|");
        try {
            b[r] = std::stod(token);
        } catch (const std::exception&) {
            result.error = "invalid numeric value: \"" + token + "\"";
            return result;
        }
    }

    result.ok = true;
    result.system = {a, b};
    return result;
}

ParseResult parseJson(const std::string& source) {
    ParseResult result;
    json root;
    try {
        root = json::parse(source);
    } catch (const json::parse_error& e) {
        result.error = std::string("JSON syntax error: ") + e.what();
        return result;
    }

    if (!root.contains("A") || !root.contains("b") || !root["A"].is_array() || !root["b"].is_array()) {
        result.error = "JSON must have array fields \"A\" (matrix of rows) and \"b\" (vector)";
        return result;
    }

    std::size_t rows = root["A"].size();
    if (rows == 0 || !root["A"][0].is_array()) {
        result.error = "\"A\" must be a non-empty array of rows";
        return result;
    }
    std::size_t cols = root["A"][0].size();

    Matrix a(rows, cols);
    for (std::size_t r = 0; r < rows; ++r) {
        if (!root["A"][r].is_array() || root["A"][r].size() != cols) {
            result.error = "every row of \"A\" must have the same length";
            return result;
        }
        for (std::size_t c = 0; c < cols; ++c) {
            if (!root["A"][r][c].is_number()) {
                result.error = "non-numeric entry in \"A\"";
                return result;
            }
            a(r, c) = root["A"][r][c].get<double>();
        }
    }

    if (root["b"].size() != rows) {
        result.error = "\"b\" must have the same number of entries as \"A\" has rows";
        return result;
    }
    Vector b(rows);
    for (std::size_t r = 0; r < rows; ++r) {
        if (!root["b"][r].is_number()) {
            result.error = "non-numeric entry in \"b\"";
            return result;
        }
        b[r] = root["b"][r].get<double>();
    }

    result.ok = true;
    result.system = {a, b};
    return result;
}

ParseResult parseFile(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        ParseResult result;
        result.error = "could not open file: " + path;
        return result;
    }
    std::ostringstream ss;
    ss << file.rdbuf();
    std::string content = ss.str();

    bool looksLikeJson = path.size() >= 5 && path.compare(path.size() - 5, 5, ".json") == 0;
    return looksLikeJson ? parseJson(content) : parseAugmentedText(content);
}

}  // namespace numsolve
