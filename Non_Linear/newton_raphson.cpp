#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double function(const vector<double> &coeff, double x)
{
    double result = 0;
    for (int i = 0; i < coeff.size(); i++)
    {
        result += coeff[i] * pow(x, coeff.size() - i - 1);
    }
    return result;
}

double newtonRaphsonMethod()
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

    vector<double> deriv(deg);
    for (int i = 0; i < deg; i++)
    {
        deriv[i] = coeff[i] * (deg - i);
    }

    double x0;
    cout << "Enter the initial guess x0: ";
    cin >> x0;

    int n;
    cout << "Enter the number of iterations: ";
    cin >> n;

    double x1;

    for (int i = 0; i < n; i++)
    {
        double f_x = function(coeff, x0);
        double df_x = function(deriv, x0);

        if (fabs(df_x) < 1e-10)
        {
            cout << "Derivative is too close to zero; Newton-Raphson method fails." << endl;
            return -1;
        }

        x1 = x0 - f_x / df_x;

        if (fabs(x1 - x0) < 1e-10)
        {
            cout << "Root found at: " << x1 <<  endl;
            return x1;
        }

        x0 = x1;
    }

    cout << "Root found at: " << x1  << endl;
    return x1;
}

int main()
{
    newtonRaphsonMethod();
    return 0;
}
