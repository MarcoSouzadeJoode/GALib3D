# GALib3D

A small C library for 3D Geometric Algebra (Euclidean signature +,+,+), designed for personal experimentation and calculations with multivectors.

This library provides basic operations for multivectors, vectors, bivectors, and trivectors in 3D space. It is lightweight, human-written, and includes several duplicate or alternative functions for verification and experimentation purposes.

 # Features

Multivector Arithmetic: Addition, subtraction, geometric product, wedge product, duals, and reversals.

Vector Operations: Dot and cross products for quick sanity checks.

Blade Detection: Functions to check if a multivector is a scalar, vector, bivector, or trivector.

Geometric Algebra Utilities:

dotpluswedge: Computes the geometric product of two vectors (scalar + bivector).

pure_blade_grade: Returns the grade of a pure blade or -1 if mixed.

Placeholders for Advanced Operations: Functions for multivector exponentials and bivector normalization are prepared for future implementation.
