#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "../galib3d.h"

int equal_mv(Multivector u, Multivector v) {
    //approx equality tester
    double eps = 1e-12;
    return (fabs(u.a - v.a) < eps) &&
           (fabs(u.b1 - v.b1) < eps) &&
           (fabs(u.b2 - v.b2) < eps) &&
           (fabs(u.b3 - v.b3) < eps) &&
           (fabs(u.c1 - v.c1) < eps) &&
           (fabs(u.c2 - v.c2) < eps) &&
           (fabs(u.c3 - v.c3) < eps) &&
           (fabs(u.d - v.d) < eps);
}

void print_multivector(Multivector m) {
    printf("a=%f b=(%f,%f,%f) c=(%f,%f,%f) d=%f\n", m.a, m.b1,m.b2,m.b3,m.c1,m.c2,m.c3,m.d);
}


void test_dotpluswedge_vs_gp() {
    Vec3 u = {2.0, 3.0, -1.0};
    Vec3 v = {-5.0, 1.0, 2.0};

    Multivector dw = dotpluswedge(u, v);
    Multivector u_mv = {0, u.x, u.y, u.z, 0,0,0,0};
    Multivector v_mv = {0, v.x, v.y, v.z, 0,0,0,0};
    Multivector gp_uv = gp(u_mv, v_mv);

    // zero out vector/trivector parts in dotpluswedge
    Multivector gp_check = {gp_uv.a, 0,0,0, gp_uv.c1, gp_uv.c2, gp_uv.c3, 0};
    print_multivector(gp_uv);
    print_multivector(dw);
    assert(equal_mv(dw, gp_uv));
    printf("dotpluswedge vs gp test passed!\n");
}



int main() {
    // Scalars
    Multivector s1 = {2,0,0,0,0,0,0,0};
    Multivector s2 = {2,0,0,0,0,0,0,0};

    // Vectors
    Multivector e1 = {0,1,0,0,0,0,0,0};
    Multivector e2 = {0,0,1,0,0,0,0,0};
    Multivector e3 = {0,0,0,1,0,0,0,0};

    // Bivectors
    Multivector e23 = {0,0,0,0,1,0,0,0};
    Multivector e31 = {0,0,0,0,0,1,0,0};
    Multivector e12 = {0,0,0,0,0,0,1,0};

    // Trivector
    Multivector e123 = {0,0,0,0,0,0,0,1};


    //testing scalar * multivec
    print_multivector(gp(s1, e31)); 
    print_multivector(gp(s1, e12)); 
    print_multivector(gp(s1, e1)); 
    print_multivector(gp(s1, e123)); 


    // Test geometric product
    print_multivector(gp(e1,e1)); // should be scalar 1
    print_multivector(gp(e1,e2)); // should be e12
    print_multivector(gp(e1,e23)); // should be e123
    print_multivector(gp(e123,e1)); // should be bivector
    print_multivector(gp(e123,e123)); // should be scalar -1

    // more on bivs:
    printf("more bivector tests\n");
    print_multivector(gp(e12, e23)); // e13 c : (0 -1 0)
    print_multivector(gp(e23, e12)); // =-e12 e23 = -e13 = e31 c : (0 1 0)
    print_multivector(gp(e12, e31)); // =-e31 e12 = -e32 = +e23 (1, 0, 0)
    print_multivector(gp(e23, e31)); // = e21 = -e12 = (0, 0, -1)


    printf("\n");

    test_dotpluswedge_vs_gp();

    printf("wedges\n");

    e12 = wedge(e1, e2);  // bivector e12

    e123 = wedge(e1, e23); // trivector

    print_multivector(e12);
    print_multivector(e123);


    return 0;
}
