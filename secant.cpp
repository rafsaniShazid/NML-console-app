#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double function( vector<double>& coeff, double x) {
    double result = 0;
    for (int i = 0; i < coeff.size(); i++) {
        result += coeff[i] * pow(x, coeff.size() - i - 1);
    }
    return result;
}
double secantMethod() {
    int deg;
    cout << "Enter the degree of the polynomial: ";
    cin >> deg;

    vector<double> coeff(deg + 1);
    cout << "Enter the coefficients (from highest to lowest degree): ";
    for (int i = 0; i <= deg; i++) {
        cin >> coeff[i];
    }

    double x0, x1;
    cout << "Enter the initial guesses x0 and x1: ";
    cin >> x0 >> x1;

    int n;
    cout << "Enter the number of iterations: ";
    cin >> n;

    double xi = x1;

    for (int i = 1; i <= n; i++) {
        double fx0 = function(coeff, x0);
        double fx1 = function(coeff, x1);

        if (fx1 == fx0) {
            cout << "Division by zero error in Secant method." << endl;
            return -1;
        }

        xi = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        x0 = x1;
        x1 = xi;
    }

    cout << "Root found at: " << xi << " after " << n << " iterations." << endl;
    return xi;
}
int main() {

    secantMethod();


    return 0;
}