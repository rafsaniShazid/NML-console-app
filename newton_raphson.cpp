#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double function(vector<double> &coeff, double x)
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
    cout << "Enter the coefficients of the derivative (from highest to lowest degree): ";
    for (int i = 0; i < deg; i++)
    {
        cin >> deriv[i];
    }

    double x0;
    cout << "Enter the initial guess x0: ";
    cin >> x0;

    int n;
    cout << "Enter the number of iterations: ";
    cin >> n;

    double x1;

    for (int i = 1; i <= n; i++)
    {
        double df_x = function(deriv, x0);
        if (df_x == 0)
        {
            cout << "Derivative is zero; Newton-Raphson method fails." << endl;
            return -1;
        }

        x1 = x0 - function(coeff, x0) / df_x;
        x0 = x1;
    }

    cout << "Root found at: " << x1 << " after " << n << " iterations." << endl;
    return x1;
}

int main()
{
    newtonRaphsonMethod();

    return 0;
}