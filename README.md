
# Numerical Methods in C++

This repository contains implementations of several numerical methods in C++ for solving systems of linear equations and ordinary differential equations (ODEs). Each method includes an explanation, input format, example, algorithm, and edge cases.

## Table of Contents
1. [Jacobi Iteration Method](#1-jacobi-iteration-method)
2. [Gauss-Seidel Iteration Method](#2-gauss-seidel-iteration-method)
3. [Gaussian Elimination Method](#3-gaussian-elimination-method)
4. [Gauss-Jordan Elimination Method](#4-gauss-jordan-elimination-method)
5. [Runge-Kutta 4th Order Method](#5-runge-kutta-4th-order-method)
6. [Matrix Inversion](#6-matrix-inversion)
7. [LU Decomposition](#7-lu-decomposition)

---

## 1. Jacobi Iteration Method

### Overview
This program solves a system of linear equations using the **Jacobi Iteration Method**. The Jacobi method is an iterative algorithm suitable for diagonally dominant systems, expressed as \( Ax = b \).

### Requirements
- **Compiler**: C++11 or later
- **OS Compatibility**: Windows, macOS, Linux

### Input Format
1. Enter the number of equations (integer).
2. Enter the coefficients of each equation followed by the constant term.
3. Specify the desired precision for the iterative solution.

#### Example Input
```
Enter the number of equations: 3
Enter the coefficients of each equation followed by the constant term:
Equation 1: 10 -1 2 6
Equation 2: -1 11 -1 25
Equation 3: 2 -1 10 -11
Enter precision: 0.001
```

#### Sample Output
```
Step 1:
x[1] = 0.600000 x[2] = 2.272727 x[3] = -1.100000 
Step 2:
x[1] = 1.137273 x[2] = 2.009091 x[3] = -0.790909
...
```

---

## 2. Gauss-Seidel Iteration Method

### Overview
This program uses the **Gauss-Seidel Iteration Method** to solve a system of linear equations. It’s an iterative method similar to the Jacobi method but usually converges faster.

### Requirements
- **Compiler**: C++11 or later
- **OS Compatibility**: Windows, macOS, Linux

### Input Format
1. Enter the number of equations (integer).
2. Enter the coefficients of each equation followed by the constant term.
3. Specify the desired precision for the iterative solution.

#### Example Input
```
Enter the number of equations: 3
Enter the coefficients of each equation followed by the constant term:
Equation 1: 4 -1 0 3
Equation 2: -1 4 -1 5
Equation 3: 0 -1 4 6
Enter precision: 0.001
```

#### Sample Output
```
Iteration 1:
x[1] = 0.750000 x[2] = 1.437500 x[3] = 1.859375 
...
```

---

## 3. Gaussian Elimination Method

### Overview
The **Gaussian Elimination Method** is a direct method for solving a system of linear equations by reducing the matrix to row echelon form.

### Requirements
- **Compiler**: C++11 or later
- **OS Compatibility**: Windows, macOS, Linux

### Input Format
1. Enter the number of equations (integer).
2. Enter the augmented matrix coefficients row by row.

#### Example Input
```
Enter the number of equations: 3
Enter the augmented matrix coefficients:
Row 1: 2 1 -1 8
Row 2: -3 -1 2 -11
Row 3: -2 1 2 -3
```

#### Sample Output
```
Solution: x = 2, y = 3, z = -1
```

---

## 4. Gauss-Jordan Elimination Method

### Overview
The **Gauss-Jordan Elimination Method** is an extension of Gaussian elimination, reducing the matrix to its reduced row echelon form.

### Requirements
- **Compiler**: C++11 or later
- **OS Compatibility**: Windows, macOS, Linux

### Input Format
1. Enter the number of equations (integer).
2. Enter the augmented matrix coefficients row by row.

#### Example Input
```
Enter the number of equations: 2
Enter the augmented matrix coefficients:
Row 1: 1 1 2
Row 2: 2 4 8
```

#### Sample Output
```
Solution: x = 2, y = 0
```

---

## 5. Runge-Kutta 4th Order Method

### Overview
The **Runge-Kutta 4th Order Method** is a numerical technique to solve ordinary differential equations (ODEs). It is a popular method due to its accuracy.

### Requirements
- **Compiler**: C++11 or later
- **OS Compatibility**: Windows, macOS, Linux

### Input Format
1. Enter the differential equation.
2. Enter the initial values, step size, and the number of steps.

#### Example Input
```
Enter the initial conditions (x0, y0): 0 1
Enter step size: 0.1
Enter the number of steps: 10
```

#### Sample Output
```
After 1 step: x = 0.100000, y = 1.105170
After 2 steps: x = 0.200000, y = 1.221402
...
```

---

## 6. Matrix Inversion

### Overview
This method performs **Matrix Inversion** to solve a system of equations. It is based on finding the inverse of the coefficient matrix.

### Requirements
- **Compiler**: C++11 or later
- **OS Compatibility**: Windows, macOS, Linux

### Input Format
1. Enter the size of the matrix (n x n).
2. Enter the matrix elements.

#### Example Input
```
Enter the size of the matrix: 3
Enter the matrix elements:
1 2 3
0 1 4
5 6 0
```

#### Sample Output
```
Inverse Matrix:
-24 18 5
20 -15 -4
-5 4 1
```

---

## 7. LU Decomposition

### Overview
The **LU Decomposition** method decomposes a matrix into a lower triangular matrix (L) and an upper triangular matrix (U) to solve a system of equations.

### Requirements
- **Compiler**: C++11 or later
- **OS Compatibility**: Windows, macOS, Linux

### Input Format
1. Enter the size of the matrix (n x n).
2. Enter the matrix elements.

#### Example Input
```
Enter the size of the matrix: 3
Enter the matrix elements:
2 3 1
4 8 3
2 6 3
```

#### Sample Output
```
Lower Matrix L:
1 0 0
2 1 0
1 0.5 1

Upper Matrix U:
2 3 1
0 2 -1
0 0 2
```

---

## License
This project is open-source and free to use.

## Contributing
Contributions are welcome! Please feel free to submit a pull request or report any issues.

## Contact
For any questions or feedback, please feel free to contact with rafsaniShazid(github username).
