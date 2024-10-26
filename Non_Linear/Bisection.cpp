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

    if (function(coeff, a) * function(coeff, b) >= 0)
    {
        cout << "Bisection method won't work with the given range." << endl;
        return -1;
    }

    double c;
    double error = abs(b - a);

    while (error > error_accept)
    {
        c = (a + b) / 2;
        if (function(coeff, c) == 0.0)
            break;

        if (function(coeff, c) * function(coeff, a) < 0)
            b = c;
        else
            a = c;

        error = abs(b - a);
    }

    cout << "Root found at: " << c << " with error: " << error << endl;
    return c;
}

int main()
{
    bisectionMethod();
    return 0;
}
