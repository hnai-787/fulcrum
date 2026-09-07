/*
Project Statement:
Write a C++ program to solve a system of linear equations using Gaussian Elimination method.
The program should input the number of unknowns (size of the system), the coefficients matrix A, and the constants vector b of the system Ax = b.
The program should form the augmented matrix [A|b] and perform Gaussian Elimination to solve the system of linear equations.
The program should output the result based on the status of the solution (unique solution, infinite solutions, or no solution).
The program should also output the augmented matrix before and after Gaussian Elimination.
*/
#include <iostream> // For standard input/output
#include <vector> // For vector container
#include <cmath> // For fabs function
#include <iomanip> // For setprecision function
using namespace std;

// Epsilon value for comparing double values
const double EPSILON = 1e-9;

void printMatrix(const vector<vector<double>>& matrix) {
    // Get the dimensions of the matrix
    int m = matrix.size();
    int n = matrix[0].size();

    // Print the matrix
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // Print the value with 2 decimal places
            cout << fixed << setprecision(1) << matrix[i][j];
            
            // Print the separator
            if (j == n - 2) {
                cout << "  :  ";
            } else {
                cout << "  ";
            }
        }
        cout << endl;
    }
}

// Function to perform Gaussian Elimination
pair<vector<double>, string> gaussianElimination(vector<vector<double>>& matrix, int n) {
    // Initialize the rank of the matrix
    int rank = 0;

    // Perform Gaussian Elimination
    for (int col = 0; col < n; ++col) {
        // Find the pivot row
        int pivotRow = col;

        // Find the row with the maximum value in the current column
        for (int i = col + 1; i < n; ++i) {
            // Check if the current value is greater than the pivot value
            if (fabs(matrix[i][col]) > fabs(matrix[pivotRow][col])) {
                pivotRow = i; // Update the pivot row
            }
        }

        // If the maximum value is less than EPSILON, continue to the next column
        if (fabs(matrix[pivotRow][col]) < EPSILON) { 
            continue; // Skip this column
        }

        // Swap the rows to make the pivot element the leading one in the current column
        swap(matrix[col], matrix[pivotRow]);
        // Increment the rank of the matrix
        rank++;

        // Perform row operations to make the current column 0 below the pivot
        for (int i = col + 1; i < n; ++i) {
            // Calculate the ratio to make the current element 0 in the current column
            double ratio = matrix[i][col] / matrix[col][col];
            
            // Perform row operation to make the current element 0 in the current column
            for (int j = col; j <= n; ++j) {
                // Subtract the product of the ratio and the pivot row from the current row
                matrix[i][j] -= ratio * matrix[col][j];
            }
        }
    }

    // Back substitution to find the solution
    vector<double> solution(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        // Check for inconsistent or dependent systems
        if (fabs(matrix[i][i]) < EPSILON) {
            // Check if the constant term is non-zero
            if (fabs(matrix[i][n]) > EPSILON) {
                // Inconsistent system
                return {{}, "no solution"}; // Return no solution
            } 
            else {
                // Dependent system
                continue; // Skip to the next row
            }
        }
        // Calculate the solution for the current variable
        solution[i] = matrix[i][n] / matrix[i][i];

        // Perform back substitution to update the remaining variables
        for (int j = i - 1; j >= 0; --j) {
            // Subtract the product of the solution and the coefficient from the constant term
            matrix[j][n] -= matrix[j][i] * solution[i];
        }
    }

    // Determine the type of solution based on the rank of the matrix
    if (rank < n) {
        // Infinite solutions
        return {solution, "infinite solutions"};
    } 
    else {
        // Unique solution
        return {solution, "unique solution"};
    }
}

// Main function
int main() {
    // Input the number of unknowns (size of the system)
    int n;
    cout << "Enter the number of unknowns: ";
    cin >> n;

    // Input the coefficients matrix A and the constants vector b of the system Ax = b
    // Initialize the coefficients matrix A
    vector<vector<double>> A(n, vector<double>(n));

    // Input the coefficients matrix A 
    cout << "Enter the coefficients matrix A (" << n << "x" << n << "):\n";
    
    // Loop to input the elements of the coefficients matrix A
    for (int i = 0; i < n; ++i) {
        // Loop to input the elements of each row of the coefficients matrix A
        for (int j = 0; j < n; ++j) {
            // Input the element A[i][j]
            cout << "A[" << i + 1 << "][" << j + 1 << "] = ";
            cin >> A[i][j];
        }
    }

    // Input the constants vector b
    // Initialize the constants vector b
    vector<double> b(n);

    // Input the constants vector b
    cout << "Enter the constants vector b (" << n << "x1):\n";
    // Loop to input the elements of the constants vector b
    for (int i = 0; i < n; ++i) {
        // Input the element b[i]
        cout << "b[" << i + 1 << "] = ";
        cin >> b[i];
    }

    // Form the augmented matrix
    // Initialize the augmented matrix with the coefficients matrix A and the constants vector b
    vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));

    // Loop to form the augmented matrix
    for (int i = 0; i < n; ++i) {
        // Copy the coefficients matrix A and the constants vector b to the augmented matrix
        for (int j = 0; j < n; ++j) {
            // Copy the element A[i][j] to the augmented matrix
            augmentedMatrix[i][j] = A[i][j];
        }
        // Copy the element b[i] to the augmented matrix
        augmentedMatrix[i][n] = b[i];
    }

    // Print the augmented matrix before Gaussian Elimination
    cout << "\nAugmented Matrix before Gaussian Elimination:\n";
    printMatrix(augmentedMatrix);

    // Perform Gaussian Elimination to solve the system of linear equations
    auto result = gaussianElimination(augmentedMatrix, n);
    vector<double> solution = result.first;
    string status = result.second;

    // Output the result based on the status of the solution
    if (status == "no solution") {
        cout << "\nThe system has no solution." << endl;
    } 
    else if (status == "infinite solutions") {
        cout << "\nThe system has infinite solutions." << endl;
    } 
    else {
        // Print the solution
        cout << "\nThe system has a unique solution:\n";
        // Loop to output the values of the unknowns
        for (int i = 0; i < n; ++i) {
            cout << "x" << i + 1 << " = " << fixed << setprecision(2) << solution[i] << endl;
        }
    }

    // Print the augmented matrix after Gaussian Elimination
    cout << "\nAugmented Matrix after Gaussian Elimination:\n";
    printMatrix(augmentedMatrix);

    // Pause the console before exiting
    cout << endl << endl;
    system("pause");

    return 0;
}
