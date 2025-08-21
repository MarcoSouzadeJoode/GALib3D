# GALib3D
A few C functions to manipulate multivectors in 3D geometric algebra ($Cl_{3,0}( {\rm I\!R} $) Clifford algebra of 3-dimensional Euclidean space over the real numbers). Computing geometric products, wedges, duals of multivectors.

Represent multivectors in 3D (scalar + vector + bivector + trivector).

I am not including any constructors from other containers ... the chance of me implementing your use case are low, and the struct is very simple.


# Contents
Compute basic operations on multivectors:

-Addition (add)
-Subtraction (sub)
-Geometric product (gp)
-Wedge product (wedge)
-Dual (dual)

Utility functions for vectors - mostly for tests of the previous functions:

-Dot product (dot)
-Cross product (cross)

Specialized vector operations:

- dotpluswedge computes the geometric product of two vectors as scalar + bivector.
(Should be $uv = u\cdot v + v \wedge w $)


```c
typedef struct {
    double a;   // scalar
    double b1;  // vector x-component
    double b2;  // vector y-component
    double b3;  // vector z-component
    double c1;  // bivector e23
    double c2;  // bivector e31
    double c3;  // bivector e12
    double d;   // trivector e123
} Multivector;
```


The crucial snippet is this (takes some time to figure out signs by hand)
```c
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
```

# Example Use
```c
#include "galib3d.h"

int main() {
    Multivector e1 = {0,1,0,0,0,0,0,0};
    Multivector e2 = {0,0,1,0,0,0,0,0};
    Multivector e12 = wedge(e1, e2);  // bivector e12

    Multivector e123 = wedge(e1, e12); // trivector

    
    Multivector u = {1, 2.5, 3, 4, 0, 0,0, -1};
    Multivector u = {1, 2.5, 3, -4, 0, 0,0, +1};
    Multivector w;
    Multivector w_dual;



    w = gp(u, v);
    w_dual = dual(w);
    return 0;
}
```

# TODO
- further tests of wedges
- rotations, reflections, connection to quaterion scripts
- regressive product
- PGA ?
  
