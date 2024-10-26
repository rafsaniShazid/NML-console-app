#include <bits/stdc++.h>
using namespace std;
void Jacobi_iteration(){
    cout << "Enter the number of equations: ";
    int n;
    cin >> n;

    vector<vector<float>> coeff(n, vector<float>(n));
    vector<float> rhs(n);
    vector<float> x(n, 0.0), x_old(n, 0.0); 
    float precision;
    
    cout << "Enter the coefficients of each equation followed by the constant term:\n";
    for (int i = 0; i < n; ++i) {
        cout << "Equation " << i + 1 << ": ";
        for (int j = 0; j < n; ++j) {
            cin >> coeff[i][j];
        }
        cin >> rhs[i];
    }

    cout << "Enter precision: ";
    cin >> precision;

    cout << setprecision(6) << fixed;

    int step = 1;
    bool stop;
    do {
        stop = true;
        cout << "Step " << step << ":\n";
        
        for (int i = 0; i < n; ++i) {
            float sum = rhs[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= coeff[i][j] * x[j];
                }
            }
            // storing the old values
            x_old[i] = x[i]; 
            x[i] = sum / coeff[i][i];

            cout << "x[" << i + 1 << "] = " << x[i] << " ";
        }
        cout << "\n";

        // Checking the precisions
        for (int i = 0; i < n; ++i) {
            if (fabs(x[i] - x_old[i]) > precision) {
                stop = false;
                break;
            }
        }

        step++;
    } while (!stop);

    cout << "\nSolution:\n";
    for (int i = 0; i < n; ++i) {
        cout << "x[" << i + 1 << "] = " << x[i] << "\n";
    }

}

int main() {
    Jacobi_iteration();
    return 0;
}
