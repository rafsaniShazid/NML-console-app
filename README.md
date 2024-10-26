# Numerical Methods Solver in C++

This repository provides a set of programs for solving various mathematical problems, including linear equations, non-linear equations, and differential equations using different numerical methods in C++.

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Supported Methods](#supported-methods)
   - [Solution of Linear Equations](#solution-of-linear-equations)
     - [Jacobi Iteration Method](#jacobi-iteration-method)
     - [Gauss-Seidel Method](#gauss-seidel-method)
     - [Gaussian Elimination Method](#gaussian-elimination-method)
     - [Gauss-Jordan Elimination Method](#gauss-jordan-elimination-method)
     - [LU Factorization](#lu-factorization)
   - [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
     - [Bisection Method](#bisection-method)
     - [False Position Method](#false-position-method)
     - [Secant Method](#secant-method)
     - [Newton-Raphson Method](#newton-raphson-method)
   - [Solution of Differential Equations](#solution-of-differential-equations)
     - [Runge-Kutta Method](#runge-kutta-method)
5. [License](#license)

## Overview
This project aims to provide solutions for a variety of numerical problems using popular numerical methods. Each method is implemented in C++ with user-friendly input prompts.

## Installation
1. **Requirements**:
   - C++11 or later.
   - A C++ compiler such as GCC or Clang.
   
2. **Compilation**:
   Use any standard C++ compiler to compile the code. For example:
   
   ```bash
   g++ main.cpp -o numerical_solver
