/*
    Linear Algebra Project
    Gaussian Elimination Solver in C++

    Features:
    - Solves Ax = b using Gaussian Elimination with partial pivoting
    - Supports any n x n system
    - Detects:
        1. Unique solution
        2. Infinite solutions
        3. No solution
    - Prints augmented matrix before and after elimination
    - Includes safer numerical comparison using EPSILON
*/

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const double EPSILON = 1e-9;

enum class SolutionStatus {
    Unique,
    Infinite,
    None
};

struct SolverResult {
    SolutionStatus status;
    vector<double> solution;
    vector<vector<double>> rowEchelonMatrix;
};

bool isZero(double value) {
    return fabs(value) < EPSILON;
}

void printAugmentedMatrix(const vector<vector<double>>& matrix) {
    cout << fixed << setprecision(3);

    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); j++) {
            if (j == row.size() - 1) {
                cout << " | ";
            }
            cout << setw(10) << row[j] << " ";
        }
        cout << '\n';
    }
}

SolverResult solveGaussianElimination(vector<vector<double>> augmented) {
    int rows = static_cast<int>(augmented.size());
    int variables = static_cast<int>(augmented[0].size()) - 1;

    int pivotRow = 0;
    vector<int> pivotColumnForRow(rows, -1);

    for (int col = 0; col < variables && pivotRow < rows; col++) {
        int bestRow = pivotRow;

        for (int r = pivotRow + 1; r < rows; r++) {
            if (fabs(augmented[r][col]) > fabs(augmented[bestRow][col])) {
                bestRow = r;
            }
        }

        if (isZero(augmented[bestRow][col])) {
            continue;
        }

        swap(augmented[pivotRow], augmented[bestRow]);
        pivotColumnForRow[pivotRow] = col;

        double pivot = augmented[pivotRow][col];

        for (int c = col; c <= variables; c++) {
            augmented[pivotRow][c] /= pivot;
        }

        for (int r = pivotRow + 1; r < rows; r++) {
            double factor = augmented[r][col];

            for (int c = col; c <= variables; c++) {
                augmented[r][c] -= factor * augmented[pivotRow][c];
            }
        }

        pivotRow++;
    }

    // Check for inconsistency: 0x + 0y + ... = non-zero
    for (int r = 0; r < rows; r++) {
        bool allCoefficientZero = true;

        for (int c = 0; c < variables; c++) {
            if (!isZero(augmented[r][c])) {
                allCoefficientZero = false;
                break;
            }
        }

        if (allCoefficientZero && !isZero(augmented[r][variables])) {
            return {SolutionStatus::None, {}, augmented};
        }
    }

    int rank = 0;
    for (int r = 0; r < rows; r++) {
        if (pivotColumnForRow[r] != -1) {
            rank++;
        }
    }

    if (rank < variables) {
        return {SolutionStatus::Infinite, {}, augmented};
    }

    vector<double> solution(variables, 0.0);

    for (int r = rank - 1; r >= 0; r--) {
        int pivotCol = pivotColumnForRow[r];

        if (pivotCol == -1) {
            continue;
        }

        double value = augmented[r][variables];

        for (int c = pivotCol + 1; c < variables; c++) {
            value -= augmented[r][c] * solution[c];
        }

        solution[pivotCol] = value;
    }

    return {SolutionStatus::Unique, solution, augmented};
}

vector<vector<double>> inputSystem() {
    int n;

    cout << "Enter the number of unknowns: ";
    if (!(cin >> n)) {
        throw invalid_argument("Number of unknowns must be a valid integer.");
    }

    if (n <= 0) {
        throw invalid_argument("Number of unknowns must be greater than zero.");
    }

    if (n > 50) {
        throw invalid_argument("Number of unknowns is too large for this educational solver.");
    }

    vector<vector<double>> augmented(n, vector<double>(n + 1));

    cout << "\nEnter the coefficients of matrix A:\n";
    for (int i = 0; i < n; i++) {
        cout << "Equation " << i + 1 << " coefficients:\n";
        for (int j = 0; j < n; j++) {
            cout << "A[" << i + 1 << "][" << j + 1 << "] = ";
            if (!(cin >> augmented[i][j]) || !isfinite(augmented[i][j])) {
                throw invalid_argument("Matrix coefficients must be valid finite numbers.");
            }
        }
    }

    cout << "\nEnter the constants vector b:\n";
    for (int i = 0; i < n; i++) {
        cout << "b[" << i + 1 << "] = ";
        if (!(cin >> augmented[i][n]) || !isfinite(augmented[i][n])) {
            throw invalid_argument("Constants vector values must be valid finite numbers.");
        }
    }

    return augmented;
}

void printResult(const SolverResult& result) {
    if (result.status == SolutionStatus::None) {
        cout << "\nResult: The system has no solution.\n";
        return;
    }

    if (result.status == SolutionStatus::Infinite) {
        cout << "\nResult: The system has infinite solutions.\n";
        return;
    }

    cout << "\nResult: The system has a unique solution:\n";
    cout << fixed << setprecision(6);

    for (size_t i = 0; i < result.solution.size(); i++) {
        cout << "x" << i + 1 << " = " << result.solution[i] << '\n';
    }
}

int main() {
    try {
        vector<vector<double>> augmented = inputSystem();

        cout << "\nAugmented Matrix Before Gaussian Elimination:\n";
        printAugmentedMatrix(augmented);

        SolverResult result = solveGaussianElimination(augmented);

        cout << "\nAugmented Matrix After Gaussian Elimination:\n";
        printAugmentedMatrix(result.rowEchelonMatrix);

        printResult(result);
    }
    catch (const exception& ex) {
        cerr << "Error: " << ex.what() << '\n';
        return 1;
    }

    return 0;
}
