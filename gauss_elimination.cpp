#include <bits/stdc++.h>
typedef long long ll;
#define forr(i, j, k) for(ll i = j; i < k; i++)
using namespace std;

void gauss_elimination() {
    // Taking the input of number of equations
    int n;
    cout << "Enter the number of equations: ";
    cin >> n;

    // Taking the input of the augmented matrix
    vector<vector<double>> mat(n, vector<double>(n + 1));
    cout << "Enter the augmented matrix (each row with n+1 values):\n";
    forr(i, 0, n) {
        forr(j, 0, n + 1) {
            cin >> mat[i][j];
        }
    }

    //Performing Gaussian Elimination
    forr(i, 0, n) {
        double max_element = abs(mat[i][i]);
        int maxRow = i;
        // Loop to find the max element of a column
        forr(k, i + 1, n) {
            if (abs(mat[k][i]) > max_element) {
                max_element = abs(mat[k][i]);
                maxRow = k;
            }
        }
        // Swaping the current row with the row containing the highest element
        swap(mat[maxRow], mat[i]);

        // Making all rows below this one 0 in the current column
        forr(k, i + 1, n) {
            double factor = mat[k][i] / mat[i][i];
            forr(j, i, n + 1) {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }

    //Performing Back Substitution
    vector<double> result(n);
    for (int i = n - 1; i >= 0; --i) {
        result[i] = mat[i][n];
        for (int j = i + 1; j < n; ++j) {
            result[i] -= mat[i][j] * result[j];
        }
        result[i] /= mat[i][i];
    }

    //Displaying the results
    cout << "Solution:\n";
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << result[i] << "\n";
    }
}

int main() {
    // Calling the gauss_elimination function
    gauss_elimination();
    return 0;
}
