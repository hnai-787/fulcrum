// This code is for the second part : Solve 6x6 system of equations with unknown coefficients and tell whether the system is consistent�or�not
#include <iostream>
#include <cmath>
#include <limits> // Include the <limits> header for clearing the input buffer

using namespace std;

const int N = 6; // Number of equations (fixed to 6)

// Function to perform Gaussian elimination
bool gaussianElimination(double mat[N][N + 1]) {
    for (int i = 0; i < N; i++) {
        // Find pivot for column i
        int max_row = i;
        for (int j = i + 1; j < N; j++) {
            if (abs(mat[j][i]) > abs(mat[max_row][i])) {
                max_row = j;
            }
        }

        // Swap the maximum row with current row (column by column)
        for (int k = i; k <= N; k++) {
            double temp = mat[i][k];
            mat[i][k] = mat[max_row][k];
            mat[max_row][k] = temp;
        }

        // Make all rows below this one 0 in current column
        for (int j = i + 1; j < N; j++) {
            double factor = mat[j][i] / mat[i][i];
            for (int k = i; k <= N; k++) {
                mat[j][k] -= factor * mat[i][k];
            }
        }
    }

    // Check if the system is consistent
    for (int i = 0; i < N; i++) {
        bool all_zero = true;
        for (int j = 0; j < N; j++) {
            if (mat[i][j] != 0) {
                all_zero = false;
                break;
            }
        }
        if (all_zero && mat[i][N] != 0) {
            return false; // Inconsistent system
        }
    }

    return true; // Consistent system
}

int main() {
    double equations[N][N + 1];

    cout << "Enter the coefficients of the 6 equations row by row (6 coefficients per equation):" << endl;
    for (int i = 0; i < N; i++) {
        cout << "Equation " << (i + 1) << ": ";
        int count = 0; // Counter for entered coefficients
        while (count < N) {
            cin >> equations[i][count];
            count++;
        }
        // Clear input buffer
        cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    if (gaussianElimination(equations)) {
        cout << "The system of equations is consistent." << endl;
        cout << "The solution is:" << endl;
        for (int i = 0; i < N; i++) {
            cout << "x" << (i + 1) << " = " << equations[i][N] / equations[i][i] << endl;
        }
    }
    else {
        cout << "The system of equations is inconsistent." << endl;
    }

    return 0;
}
