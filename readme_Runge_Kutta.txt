# Runge-Kutta 4th Order (RK4) Method in C++

This C++ program implements the Runge-Kutta 4th order method (RK4) to solve ordinary differential equations (ODEs) of the form \( y' = f(x, y) \). The RK4 method is a numerical technique that approximates the solution of an ODE given an initial value and a step size.

## Overview

The Runge-Kutta 4th Order (RK4) method is widely used for solving ODEs because of its high accuracy and efficiency. The formula calculates \( y_{n+1} \) based on the current value \( y_n \), the derivative \( f(x, y) \), and intermediate steps to estimate the next value.

## Formula

Given an initial value problem:


y' = f(x, y), \quad y(x_0) = y_0


The RK4 formula to compute \( y_{n+1} \) from \( y_n \) is:


y_{n+1} = y_n + \frac{1}{6} \left( k_1 + 2k_2 + 2k_3 + k_4 \right)


Where:
- \( k_1 = h \f(x_n, y_n) \)
- \( k_2 = h \f(x_n + \frac{h}{2}, y_n + \frac{k_1}{2}) \)
- \( k_3 = h \f(x_n + \frac{h}{2}, y_n + \frac{k_2}{2}) \)
- \( k_4 = h \f(x_n + h, y_n + k_3) \)

## Code Explanation

The program defines:
- **`f(x, y)`**: The differential equation \( y' = f(x, y) \). Replace the function body with your specific equation.
- **`rungeKutta(x0, y0, x, h)`**: A function that computes \( y \) at a specified \( x \) using RK4, given the initial values \( x_0 \), \( y_0 \), target \( x \), and step size \( h \).

