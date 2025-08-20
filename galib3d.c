/*
A little collection of functions & structs to deal with  3D Multivectors
in Geometric Algebra, with Euclidean signature (+,+,+) (Vanilla/Vector Geometric Algebra; VGA)
for my personal use, contains many "duplicate" functions (can be useful for double checking)

// code writen by a human
Marco Souza de Joode, summer 2025


*/

#include "galib3d.h"


Multivector add(Multivector u, Multivector v) {
    Multivector w;
    w.a  = u.a  + v.a;
    w.b1 = u.b1 + v.b1;
    w.b2 = u.b2 + v.b2;
    w.b3 = u.b3 + v.b3;
    w.c1 = u.c1 + v.c1;
    w.c2 = u.c2 + v.c2;
    w.c3 = u.c3 + v.c3;
    w.d  = u.d  + v.d;
    return w;
}

Multivector sub(Multivector u, Multivector v) {
    Multivector w;
    w.a  = u.a  - v.a;
    w.b1 = u.b1 - v.b1;
    w.b2 = u.b2 - v.b2;
    w.b3 = u.b3 - v.b3;
    w.c1 = u.c1 - v.c1;
    w.c2 = u.c2 - v.c2;
    w.c3 = u.c3 - v.c3;
    w.d  = u.d  - v.d;
    return w;
}

Multivector dual(Multivector v) {
    Multivector D;
    D.a = v.d;
    D.d = v.a; // scalar <-> pseudoscalar parts exchanged
    D.b1 = v.c1;
    D.b2 = v.c2;
    D.b3 = v.c3;

    D.c1 = v.b1;
    D.c2 = v.b2;
    D.c3 = v.b3;

    return D;
}

double dot(Vec3 u, Vec3 v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

Vec3 cross(Vec3 u, Vec3 v) {
    // vector cross product
    // useful for sanity checks
    Vec3 w;
    w.x = u.y * v.z - u.z * v.y;
    w.y = u.z * v.x - u.x * v.z;
    w.z = u.x * v.y - u.y * v.x;
    return w;
}


Multivector dotpluswedge(Vec3 u, Vec3 v) {
    // calculating the geometric product of two vectors = dot product + wedge product
    // e.g., of a term with lower and higher grade (scalar + bivector)

    Multivector w;

    // scalar part
    w.a = dot(u, v);

    // geom. product of two vectors does not contain a vectorial part
    w.b1 = 0;
    w.b2 = 0;
    w.b3 = 0;

    // bivector part: plane perpendicular to cross product direction = dual
    w.c1 = u.y * v.z - u.z * v.y;
    w.c2 = u.z * v.x - u.x * v.z;
    w.c3 = u.x * v.y - u.y * v.x;


    // no trivector part
    w.d = 0;
    return w;
}


Multivector wedge(Multivector u, Multivector v) {
    // the wedge product (=join)
    Multivector w;

    // the scalar part is the only scalar contribution
    w.a = u.a * v.a;

    // to the vectorial part, scalar * vector is the only contribution
    w.b1 = u.a * v.b1 + v.a * u.b1;
    w.b2 = u.a * v.b2 + v.a * u.b2;
    w.b3 = u.a * v.b3 + v.a * u.b3;

    // to the bivector part, we have a vector ^ vector contribution + a scalar * bivector contribution
    w.c1 = u.b2 * v.b3 - u.b3 * v.b2 + u.a * v.c1 + v.a * u.c1;   
    w.c2 = u.b3 * v.b1 - u.b1 * v.b3 + u.a * v.c2 + v.a * u.c2;   
    w.c3 = u.b1 * v.b2 - u.b2 * v.b1 + u.a * v.c3 + v.a * u.c3;
    
    // trivector part (oriented volume) = (scalar * trivector terms) + (vector ^ bivector terms)
    w.d = (u.a * v.d + v.a * u.d) + (u.b1 * v.c1 + u.b2 * v.c2 + u.b3 * v.c3);

    return w;

}

Multivector gp(Multivector u, Multivector v) {
    // geometric product of two 3D = 8 component multivectors
    Multivector w;
    w.a = u.a * v.a + (u.b1 * v.b1 + u.b2 * v.b2 + u.b3 * v.b3) - u.c1 * v.c1 - u.c2 * v.c2 - u.c3 * v.c3 - (u.d * v.d);

    w.b1 = u.a*v.b1 + u.b1*v.a + (u.c3*v.b2 - u.c2*v.b3) + (u.b3*v.c2 - u.b2*v.c3) - (u.d * v.c1 + u.c1 * v.d);
    w.b2 = u.a*v.b2 + u.b2*v.a + (u.c1*v.b3 - u.c3*v.b1) + (u.b1*v.c3 - u.b3*v.c1) - (u.d * v.c2 + u.c2 * v.d);
    w.b3 = u.a*v.b3 + u.b3*v.a + (u.c2*v.b1 - u.c1*v.b2) + (u.b2*v.c1 - u.b1*v.c2) - (u.d * v.c3 + u.c3 * v.d);

    w.c1 = u.a*v.c1 + u.c1*v.a + (u.b2*v.b3 - u.b3*v.b2) + (u.b1*v.d + u.d*v.b1) + ( -u.c2*v.c3 + u.c3*v.c2);
    w.c2 = u.a*v.c2 + u.c2*v.a + (u.b3*v.b1 - u.b1*v.b3) + (u.b2*v.d + u.d*v.b2) + ( -u.c3*v.c1 + u.c1*v.c3);
    w.c3 = u.a*v.c3 + u.c3*v.a + (u.b1*v.b2 - u.b2*v.b1) + (u.b3*v.d + u.d*v.b3) + ( -u.c1*v.c2 + u.c2*v.c1);

    w.d =  u.a*v.d + u.d*v.a + (u.b1*v.c1 + u.b2*v.c2 + u.b3*v.c3) + (u.c1*v.b1 + u.c2*v.b2 + u.c3*v.b3);

    return w;
}