#include <bits/stdc++.h>
typedef long long ll;
#define forr(i, j, k) for(ll i = j; i < k; i++)
using namespace std;

void gauss_jordan_elimination() {
    //  Enter the number of equations
    int n;
    cout << "Enter the number of equations: ";
    cin >> n;

    // Taking Input of the augmented matrix
    vector<vector<double>> mat(n, vector<double>(n + 1));
    cout << "Enter the augmented matrix (each row with n+1 values):\n";
    forr(i, 0, n) {
        forr(j, 0, n + 1) {
            cin >> mat[i][j];
        }
    }

    //Performing Gauss-Jordan Elimination
    forr(i, 0, n) {
        // Finding the pivot element for column i
        double max_element = abs(mat[i][i]);
        int maxRow = i;
        forr(k, i + 1, n) {
            if (abs(mat[k][i]) > max_element) {
                max_element = abs(mat[k][i]);
                maxRow = k;
            }
        }

        // Swaping the current row with the row containing the highest element
        swap(mat[maxRow], mat[i]);

        // Normalizing the pivot row
        double pivot = mat[i][i];
        forr(j, i, n + 1) {
            mat[i][j] /= pivot;
        }

        // Making all other rows 0 in the current column
        forr(k, 0, n) {
            if (k != i) {
                double factor = mat[k][i];
                forr(j, i, n + 1) {
                    mat[k][j] -= factor * mat[i][j];
                }
            }
        }
    }

    //Displayong the results
    cout << "Solution:\n";
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << mat[i][n] << "\n";
    }
}

int main() {
    // Calling the gauss_jordan_elimination function
    gauss_jordan_elimination();
    return 0;
}
