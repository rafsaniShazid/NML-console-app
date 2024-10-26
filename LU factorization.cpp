#include <iostream>
#include <vector>
using namespace std;

// Function to perform LU Decomposition
void LU_Decomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, int n) {
    for (int i = 0; i < n; i++) {
        // Upper Triangular Matrix (U)
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += (L[i][j] * U[j][k]);
            }
            U[i][k] = A[i][k] - sum;
        }

        // Lower Triangular Matrix (L)
        for (int k = i; k < n; k++) {
            if (i == k) {
                L[i][i] = 1;  // Diagonal as 1
            } else {
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += (L[k][j] * U[j][i]);
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Function to solve Ly = b using forward substitution
vector<double> Forward_Substitution(vector<vector<double>>& L, vector<double>& b, int n) {
    vector<double> y(n, 0);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }
    return y;
}

// Function to solve Ux = y using backward substitution
vector<double> Backward_Substitution(vector<vector<double>>& U, vector<double>& y, int n) {
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
    return x;
}

int main() {
    int n;

    // Input the number of variables (size of matrix)
    cout << "Enter the number of variables: ";
    cin >> n;

    // Input the coefficients of the matrix A
    vector<vector<double>> A(n, vector<double>(n));
    cout << "Enter the coefficients of the matrix row-wise:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    // Input the constants on the right-hand side vector b
    vector<double> b(n);
    cout << "Enter the constants (b vector):\n";
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    // Initialize L and U matrices
    vector<vector<double>> L(n, vector<double>(n, 0));  // Lower triangular matrix
    vector<vector<double>> U(n, vector<double>(n, 0));  // Upper triangular matrix

    // Perform LU Decomposition
    LU_Decomposition(A, L, U, n);

    // Print Lower Triangular Matrix L
    cout << "Lower Triangular Matrix (L):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << L[i][j] << "\t";
        }
        cout << endl;
    }

    // Print Upper Triangular Matrix U
    cout << "\nUpper Triangular Matrix (U):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << U[i][j] << "\t";
        }
        cout << endl;
    }

    // Step 1: Solve Ly = b using forward substitution
    vector<double> y = Forward_Substitution(L, b, n);

    // Step 2: Solve Ux = y using backward substitution
    vector<double> x = Backward_Substitution(U, y, n);

    // Print the solution vector x
    cout << "\nSolution of the system (x vector):\n";
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }

    return 0;
}