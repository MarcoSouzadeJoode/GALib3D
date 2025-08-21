#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../galib3d.h"

Multivector pseudoscalar() {
    Multivector I = {0,0,0,0,0,0,0,1}; // only trivector component
    return I;
}

int main() {
    srand(time(NULL));
    Multivector I = pseudoscalar();
    double eps = 1e-10;

    for (int t = 0; t < 500; t++) {
        // random vectors
        Vec3 va = {
            ((double)rand()/RAND_MAX)*2.0 - 1.0,
            ((double)rand()/RAND_MAX)*2.0 - 1.0,
            ((double)rand()/RAND_MAX)*2.0 - 1.0
        };
        Vec3 vb = {
            ((double)rand()/RAND_MAX)*2.0 - 1.0,
            ((double)rand()/RAND_MAX)*2.0 - 1.0,
            ((double)rand()/RAND_MAX)*2.0 - 1.0
        };

        Multivector A = {0, va.x, va.y, va.z, 0,0,0,0};
        Multivector B = {0, vb.x, vb.y, vb.z, 0,0,0,0};

        // wedge product in GA
        Multivector W = wedge(A,B);

        // classical cross product
        Vec3 c = cross(va, vb);

        Multivector C = (Multivector){0, c.x, c.y, c.z, 0,0,0,0};
        Multivector DualCross = gp(I, C);

        double dx = W.c1 - DualCross.c1;
        double dy = W.c2 - DualCross.c2;
        double dz = W.c3 - DualCross.c3;
        double err = sqrt(dx*dx + dy*dy + dz*dz);

        printf("Test %d: wedge vs dual(cross), error = %.3e\n", t, err);
        if (err > eps) {
            printf("  wedge: (c1=%.6f, c2=%.6f, c3=%.6f)\n", W.c1, W.c2, W.c3);
            printf("  dualX: (c1=%.6f, c2=%.6f, c3=%.6f)\n", DualCross.c1, DualCross.c2, DualCross.c3);
        }
    }

    return 0;
}
