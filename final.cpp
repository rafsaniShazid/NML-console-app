#include<bits/stdc++.h>
typedef long long ll;
#define forr(i, j, k) for(ll i = j; i < k; i++)
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

void gauss_siedel(){

    int n;
    cout << "Enter the number of variables/equations: ";
    cin >> n;

    // Vectors to hold coefficients and constants
    vector<vector<float>> a(n, vector<float>(n)); // Coefficients
    vector<float> r(n);                           // Constants
    vector<float> x(n, 0);                        // Solution vector (initial guesses)
    vector<float> x_new(n);                       // Updated solution
    float e;                                      // Precision

    
    cout << "Enter the coefficients for each equation (row-wise):\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> a[i][j]; 
        }
        cin >> r[i]; 
    }

    cout << setprecision(6) << fixed;
    cout << "Enter precision: ";
    cin >> e;

    int step = 1;
    float max_error; 
    do
    {
        
        for (int i = 0; i < n; i++)
        {
            float sum = r[i];
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    sum -= a[i][j] * x[j];
                }
            }
            x_new[i] = sum / a[i][i]; 
        }

        
        cout << step << "\t";
        for (int i = 0; i < n; i++)
        {
            cout << "x" << i + 1 << "= " << x_new[i] << ", ";
        }
        cout << endl;

        
        max_error = 0; 
        for (int i = 0; i < n; i++)
        {
            max_error = max(max_error, fabs(x_new[i] - x[i])); 
        }

        
        x = x_new; 

        step++;
    } while (max_error > e); 
}

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

void LU_factorization() {
    int n;

    //Taking Input of the number of variables (size of matrix)
    cout << "Enter the number of variables: ";
    cin >> n;

    //Taking Input of the coefficients of the matrix A
    vector<vector<double>> A(n, vector<double>(n));
    cout << "Enter the coefficients of the matrix row-wise:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    //Taking Input of the constants on the right-hand side vector b
    vector<double> b(n);
    cout << "Enter the constants (b vector):\n";
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    // Initializing L and U matrices
    vector<vector<double>> L(n, vector<double>(n, 0));  // Lower triangular matrix
    vector<vector<double>> U(n, vector<double>(n, 0));  // Upper triangular matrix

    // Performing LU Decomposition
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

    // Printing Lower Triangular Matrix L
    cout << "Lower Triangular Matrix (L):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << L[i][j] << "\t";
        }
        cout << endl;
    }

    // Printing Upper Triangular Matrix U
    cout << "\nUpper Triangular Matrix (U):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << U[i][j] << "\t";
        }
        cout << endl;
    }

    // Solving Ly = b using forward substitution
    vector<double> y(n, 0);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }

    //Solving Ux = y using backward substitution
    vector<double> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    // Printing the solution vector x
    cout << "\nSolution of the system (x vector):\n";
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}

double function1(vector<double> &coeff, double x)
{
    double result = 0;
    for (int i = 0; i < coeff.size(); i++)
    {
        result += coeff[i] * pow(x, coeff.size() - i - 1);
    }
    return result;
}

double bisectionMethod()
{
    int deg;
    cout << "Enter the degree of the polynomial: ";
    cin >> deg;

    vector<double> coeff(deg + 1);
    cout << "Enter the coefficients (from highest to lowest degree): ";
    for (int i = 0; i <= deg; i++)
    {
        cin >> coeff[i];
    }

    double a, b, error_accept;
    cout << "Enter the interval [a, b]: ";
    cin >> a >> b;

    cout << "Enter the acceptable error: ";
    cin >> error_accept;

    if (function1(coeff, a) * function1(coeff, b) >= 0)
    {
        cout << "Bisection method won't work with the given range." << endl;
        return -1;
    }

    double c;
    double error = abs(b - a);

    while (error > error_accept)
    {
        c = (a + b) / 2;
        if (function1(coeff, c) == 0.0)
            break;

        if (function1(coeff, c) * function1(coeff, a) < 0)
            b = c;
        else
            a = c;

        error = abs(b - a);
    }

    cout << "Root found at: " << c << " with error: " << error << endl;
    return c;
}

double function2(const vector<double> &coeff, double x) {
    double result = 0;
    for (int i = 0; i < coeff.size(); i++) {
        result += coeff[i] * pow(x, coeff.size() - i - 1);
    }
    return result;
}

double falsePositionMethod() {
    int deg;
    cout << "Enter the degree of the polynomial: ";
    cin >> deg;

    vector<double> coeff(deg + 1);
    cout << "Enter the coefficients (from highest to lowest degree): ";
    for (int i = 0; i <= deg; i++) {
        cin >> coeff[i];
    }

    double a, b, error_accept;
    cout << "Enter the interval [a, b]: ";
    cin >> a >> b;

    cout << "Enter the acceptable error: ";
    cin >> error_accept;

    if (function2(coeff, a) * function2(coeff, b) >= 0) {
        cout << "False Position method is invalid in this range." << endl;
        return -1;
    }

    double c_before = a;
    double c;

    do {
        c = (a * function2(coeff, b) - b * function2(coeff, a)) /
            (function2(coeff, b) - function2(coeff, a));

        if (function2(coeff, c) * function2(coeff, a) < 0)
            b = c;
        else
            a = c;

        double error = abs(c - c_before);
        c_before = c;

        cout << "Current approximation: " << c << ", Error: " << error << endl;

    } while (abs(c - c_before) > error_accept);

    cout << "Root found at: " << c << " with error: " << abs(c - c_before) << endl;
    return c;
}

void invertMatrix() {
    int n;
    cout << "Enter the size of the matrix: ";
    cin >> n;

    vector<vector<double>> matrix(n, vector<double>(n));
    vector<vector<double>> inverse(n, vector<double>(n, 0.0));

    cout << "Enter the elements of the matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> matrix[i][j];
        }
        inverse[i][i] = 1.0;
    }

    for (int i = 0; i < n; ++i) {
        double pivot = matrix[i][i];
        int pivotRow = i;

        for (int j = i + 1; j < n; ++j) {
            if (fabs(matrix[j][i]) > fabs(pivot)) {
                pivot = matrix[j][i];
                pivotRow = j;
            }
        }

        if (fabs(pivot) < 1e-9) {
            cout << "Matrix is singular and cannot be inverted." << endl;
            return;
        }

        if (pivotRow != i) {
            swap(matrix[i], matrix[pivotRow]);
            swap(inverse[i], inverse[pivotRow]);
        }

        pivot = matrix[i][i];
        for (int j = 0; j < n; ++j) {
            matrix[i][j] /= pivot;
            inverse[i][j] /= pivot;
        }

        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double factor = matrix[j][i];
                for (int k = 0; k < n; ++k) {
                    matrix[j][k] -= factor * matrix[i][k];
                    inverse[j][k] -= factor * inverse[i][k];
                }
            }
        }
    }

    cout << "The inverse matrix is:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << fixed << setprecision(5) << inverse[i][j] << " ";
        }
        cout << endl;
    }
}

void runge_kutta() {
    double x0 = 0;
    double y0 = 1;
    double x;
    double h;

    cout << "Enter the value of x and h: ";
    cin >> x >> h;

    double y = y0;
    int n = static_cast<int>((x - x0) / h);

    for (int i = 0; i < n; i++) {
        double k1 = h * (x0 + y);
        double k2 = h * (x0 + h / 2.0 + y + k1 / 2.0);
        double k3 = h * (x0 + h / 2.0 + y + k2 / 2.0);
        double k4 = h * (x0 + h + y + k3);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        x0 += h;
    }

    cout << "The value of y at x=" << x << " is: " << y << endl;
}


int main(){
    
    return 0;
}
