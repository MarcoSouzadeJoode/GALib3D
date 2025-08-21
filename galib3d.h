/*
A little collection of functions & structs to deal with  3D Multivectors
in Geometric Algebra, with Euclidean signature (+,+,+) (Vanilla/Vector Geometric Algebra; VGA)
for my personal use, contains many "duplicate" functions (can be useful for double checking)

// code writen by a human
Marco Souza de Joode, summer 2025


*/
#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

// 3D multivector, i.e., scalar + 3D vector + 3D bivector + pseudoscalar
typedef struct {
    // using the following basis:
    // v = a + b1e1 + b2e2 + b3e3 + c1e23 + c2e31 + c3e12 + de123
    double a;  // scalar
    double b1, b2, b3; // vector components
    double c1, c2, c3; // bi-vector = pseudovector components
    double d; // trivector = pseudoscalar component
} Multivector;

typedef struct {
    double x, y, z;
} Vec3;

double dot(Vec3 u, Vec3 v);
Vec3 cross(Vec3 u, Vec3 v);

Multivector add(Multivector u, Multivector v);
Multivector sub(Multivector u, Multivector v);
Multivector dual(Multivector v);
Multivector dotpluswedge(Vec3 u, Vec3 v);
Multivector wedge(Multivector u, Multivector v);
Multivector gp(Multivector u, Multivector v);


Multivector reverse(Multivector M);
Multivector bivector_normalize(Multivector B);



Multivector rotor(Multivector B, double theta);
Multivector rotate_vector(Multivector v, Multivector B, double theta);


int is_scalar(Multivector M);
int is_vec(Multivector M);
int is_bivec(Multivector M);
int is_trivec(Multivector M);
int pure_blade_grade(Multivector M);


#endif
