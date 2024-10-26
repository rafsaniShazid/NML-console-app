#include <iostream>
#include <cmath>
using namespace std;

// Define the function f(x, y) representing dy/dx = f(x, y)
double f(double x, double y) {
    return x + y;   //dy/dx = x + y
}

// Runge-Kutta 4th order method
double rungeKutta(double x0, double y0, double x, double h) {
    // Initializing values of x and y
    double y = y0;

    // Number of steps to reach from x0 to x
    int n = static_cast<int>((x - x0) / h);

    // Iterate to find the value of y at each step
    for (int i = 0; i < n; i++) {
        double k1 = h * f(x0, y);
        double k2 = h * f(x0+h/2.0, y+k1/2.0);
        double k3 = h * f(x0+h/2.0, y+k2/2.0);
        double k4 = h * f(x0+h, y+k3);

        // Updating y using the RK4 formula
        y += (k1+2*k2+2*k3+k4)/6.0;
        x0 += h;
    }

    return y;
}

int main() {
    double x0 = 0;      // Initial x value
    double y0 = 1;      // Initial y value, e.g., y(0) = 1
    double x;       // x value at which we need the solution
    double h;     // Step size
    cout<<"Enter the value of x and h"<<endl;
    cin>>x>>h;

    double result = rungeKutta(x0, y0, x, h);
    cout << "The value of y at x=" << x << " is: " << result << endl;

    return 0;
}
