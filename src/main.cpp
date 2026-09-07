#include <iostream>
#include <sstream>

#include "numsolve/condition.hpp"
#include "numsolve/diagnostics.hpp"
#include "numsolve/elimination.hpp"
#include "numsolve/parser.hpp"
#include "numsolve/report.hpp"

using namespace numsolve;

namespace {

void printHelp() {
    std::cout <<
        "numsolve -- explainable Gaussian elimination / PLU solver\n\n"
        "USAGE:\n"
        "  numsolve solve <file|-> [--steps] [--view elimination|plu] [--diagnostics]\n"
        "                          [--format text|json] [--all]\n"
        "  numsolve --help\n\n"
        "  <file> is a plain augmented-matrix text file or a .json file (see README).\n"
        "  Use \"-\" to read from stdin.\n"
        "  --all shows the elimination trace, the PLU view, and diagnostics together.\n";
}

std::string readAll(std::istream& in) {
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        printHelp();
        return argc < 2 ? 2 : 0;
    }

    std::string command = argv[1];
    if (command != "solve") {
        std::cerr << "Error: unrecognized command \"" << command << "\"\n\n";
        printHelp();
        return 2;
    }

    if (argc < 3) {
        std::cerr << "Error: solve requires a <file|-> argument\n";
        return 2;
    }

    std::string path = argv[2];
    bool showSteps = false, showPlu = false, showDiagnostics = false, jsonOutput = false;
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--steps") showSteps = true;
        else if (arg == "--view") {
            if (i + 1 >= argc) { std::cerr << "Error: --view requires a value\n"; return 2; }
            std::string view = argv[++i];
            if (view == "elimination") showSteps = true;
            else if (view == "plu") showPlu = true;
            else { std::cerr << "Error: --view must be \"elimination\" or \"plu\"\n"; return 2; }
        } else if (arg == "--diagnostics") showDiagnostics = true;
        else if (arg == "--all") { showSteps = showPlu = showDiagnostics = true; }
        else if (arg == "--format") {
            if (i + 1 >= argc) { std::cerr << "Error: --format requires a value\n"; return 2; }
            std::string fmt = argv[++i];
            if (fmt == "json") jsonOutput = true;
            else if (fmt != "text") { std::cerr << "Error: --format must be \"text\" or \"json\"\n"; return 2; }
        } else {
            std::cerr << "Error: unrecognized option \"" << arg << "\"\n";
            return 2;
        }
    }

    ParseResult parsed = (path == "-") ? [] {
        std::string content = readAll(std::cin);
        // stdin has no filename to sniff; try JSON first, fall back to text.
        ParseResult asJson = parseJson(content);
        return asJson.ok ? asJson : parseAugmentedText(content);
    }() : parseFile(path);

    if (!parsed.ok) {
        std::cerr << "Error: could not parse input -- " << parsed.error << "\n";
        return 2;
    }

    const LinearSystem& system = parsed.system;
    if (system.a.rows() != system.a.cols()) {
        std::cerr << "Error: this solver currently requires a square system (got "
                  << system.a.rows() << "x" << system.a.cols() << "). See README Limitations.\n";
        return 2;
    }

    FactorizationResult fac;
    try {
        fac = factorize(system.a);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    }

    Vector eliminatedB = applyEliminationToVector(fac, system.b);
    ClassificationResult classification = classify(fac, eliminatedB, fac.tolerance);

    Vector solution;
    const Vector* solutionPtr = nullptr;
    if (classification.status == SolutionStatus::Unique) {
        solution = backSubstitute(fac, eliminatedB);
        solutionPtr = &solution;
    }

    Diagnostics diagnostics;
    diagnostics.pivot = computePivotDiagnostics(fac);
    diagnostics.growth = computeGrowthDiagnostics(fac);
    diagnostics.condition = estimateCondition(system.a, fac);
    diagnostics.factorizationCheck = verifyFactorization(system.a, fac);
    if (solutionPtr) diagnostics.residual = computeResidual(system.a, system.b, solution);

    if (jsonOutput) {
        std::cout << renderJson(system, fac, classification, solutionPtr, diagnostics) << "\n";
        return 0;
    }

    std::cout << renderMatrix(system.a, "A =") << "\n" << renderVector(system.b, "b =") << "\n";

    if (showSteps) {
        std::cout << "Elimination trace (Gaussian elimination with partial pivoting)\n";
        std::cout << "================================================================\n";
        std::cout << renderEliminationTrace(system, fac, eliminatedB) << "\n";
    }

    std::cout << renderClassification(classification) << "\n";

    if (solutionPtr) {
        std::cout << "Solution:\n";
        for (std::size_t i = 0; i < solution.size(); ++i) {
            std::cout << "  x" << (i + 1) << " = " << solution[i] << "\n";
        }
        std::cout << "\n";
    }

    if (showPlu) {
        std::cout << renderPluView(fac, diagnostics.factorizationCheck) << "\n";
    }
    if (showDiagnostics) {
        std::cout << renderDiagnostics(diagnostics) << "\n";
    }

    return 0;
}
