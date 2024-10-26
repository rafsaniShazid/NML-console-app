#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

bool invertMatrix(vector<vector<double>>& matrix, vector<vector<double>>& inverse) {
    int n = matrix.size();
    
    // Initializing the inverse matrix as the identity matrix
    inverse.assign(n, vector<double>(n, 0.0));
    for (int i=0; i<n; ++i) {
        inverse[i][i] = 1.0;
    }
    
    for (int i=0; i<n; ++i) {
        double pivot=matrix[i][i];
        int pivotRow=i;
        
        for (int j=i+1;j<n;++j) {
            if (fabs(matrix[j][i]) > fabs(pivot)) {
                pivot=matrix[j][i];
                pivotRow=j;
            }
        }

        // If the pivot is zero, the matrix is singular and not invertible
        if (fabs(pivot) < 1e-9) {
            cout << "Matrix is singular and cannot be inverted." << endl;
            return false;
        }

        // Swap rows in both matrix and inverse
        if (pivotRow !=i) {
            swap(matrix[i], matrix[pivotRow]);
            swap(inverse[i], inverse[pivotRow]);
        }

        // Scale the pivot row to make pivot element 1
        pivot = matrix[i][i];
        for (int j=0;j<n;++j) {
            matrix[i][j] /= pivot;
            inverse[i][j] /= pivot;
        }

        // Eliminate the current column in other rows
        for (int j=0; j<n; ++j) {
            if (j != i) {
                double factor=matrix[j][i];
                for (int k=0;k<n; ++k) {
                    matrix[j][k] -= factor*matrix[i][k];
                    inverse[j][k] -= factor*inverse[i][k];
                }
            }
        }
    }
    return true;
}

int main() {
    int n;
    cout << "Enter the size of the matrix: ";
    cin >> n;

    vector<vector<double>> matrix(n, vector<double>(n));
    vector<vector<double>> inverse;

    cout << "Enter the elements of the matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> matrix[i][j];
        }
    }

    if (invertMatrix(matrix, inverse)) {
        cout << "The inverse matrix is:" << endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << fixed << setprecision(5) << inverse[i][j] << " ";
            }
            cout << endl;
        }
    }

    return 0;
}
