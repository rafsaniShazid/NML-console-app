1. Runge-Kutta 4th Order (RK4) Method in C++

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

##Code Explanation

The program defines:
- **`f(x, y)`**: The differential equation \( y' = f(x, y) \). Replace the function body with your specific equation.
- **`rungeKutta(x0, y0, x, h)`**: A function that computes \( y \) at a specified \( x \) using RK4, given the initial values \( x_0 \), \( y_0 \), target \( x \), and step size \( h \).

Matrix Inversion in C++
Overview
This program performs matrix inversion on an n×n matrix using Gaussian elimination with partial pivoting. The program takes the matrix size and elements as input from the console, computes the inverse (if it exists), and outputs the result.

Requirements
Compiler: C++11 or later (to support vector, iostream, and other standard libraries)
OS Compatibility: Works on Windows, MacOS, and Linux.
How It Works
The program follows these steps:

Initialize Identity Matrix:

The inverse matrix starts as an identity matrix of the same size as the input matrix.
Gaussian Elimination with Partial Pivoting:

The program iteratively scales each row to make the pivot element equal to 1.
The rows are swapped to ensure numerical stability using partial pivoting.
Each pivot element in the matrix row is used to clear out elements above and below it in the same column across other rows.
Invertibility Check:

If a zero pivot is encountered, the matrix is singular (non-invertible), and the program displays a message accordingly.
Output:

If the matrix is invertible, the inverse matrix is displayed with each element formatted to 5 decimal places for precision.

Usage:

Compilation:To compile the program, use a C++ compiler.


The program expects:
An integer n, the size of the matrix.
The elements of the matrix, given row by row.

Example Input:
Enter the size of the matrix: 3
Enter the elements of the matrix:
2 5 7
6 3 4
5 -2 -3

Sample Output:
For the above input, the program will output:

The inverse matrix is:
1.00000 -1.00000 1.00000 
-38.00000 41.00000 -34.00000 
27.00000 -29.00000 24.00000 

Edge Cases:
Non-invertible matrix: If the matrix is singular, the program outputs "Matrix is singular and cannot be inverted."
Precision: The program rounds each element of the inverse matrix to 5 decimal places.

Notes:
Floating-Point Precision: Small inaccuracies may occur in floating-point calculations; however, they are minimal for most use cases.
Matrix Size Limit: No explicit size limit is set, but performance may decrease for very large matrices due to computational complexity.


3.LU Decomposition in C++
Overview:
This program uses LU Decomposition to solve a system of linear equations represented by  A⋅x=b.  The LU decomposition factorizes matrix A into a lower triangular matrix L and an upper triangular matrix U, where A=L⋅U. The program then uses forward and backward substitution to solve for vector x.

Requirements:
Compiler: C++11 or later (supports vector and other standard libraries)
Compatibility: Works on Windows, MacOS, and Linux.

How It Works
The program proceeds in three main steps:

LU Decomposition:

The matrix A is decomposed into two matrices: L (lower triangular matrix) and U (upper triangular matrix).L has 1s on its diagonal, and U is an upper triangular matrix.
Forward Substitution:

Solves L⋅y=b for y using forward substitution.
Backward Substitution:

Solves 
U⋅x=y for x using backward substitution. The result is the solution vector x that satisfies the original equation A⋅x=b.

Usage

Compile the program using any C++ compiler.

Run the compiled program from the terminal:

Input Format
The program prompts for:

An integer n, the number of variables (or the size of the matrix).
The n×n matrix A, with coefficients entered row by row.
The vector b, containing constants on the right-hand side.

Example Input:
Enter the number of variables: 3
Enter the coefficients of the matrix row-wise:
2 1 1
4 -6 0
-2 7 2
Enter the constants (b vector):
5 -2 9

Sample Output:
For the above input, the output will be:

Lower Triangular Matrix (L):
1.00000     0.00000     0.00000    
2.00000     1.00000     0.00000    
-1.00000    -0.50000    1.00000    

Upper Triangular Matrix (U):
2.00000     1.00000     1.00000    
0.00000     -8.00000    -2.00000    
0.00000     0.00000     1.00000    

Solution of the system (x vector):
x1 = 1.00000
x2 = -1.00000
x3 = 2.00000

Edge Cases:
Singular Matrix: The program does not handle singular matrices and assumes an invertible matrix.

Precision: The output is rounded to 5 decimal places for readability.
