/*
Test: rotating 10^6 vectors around 10^6 axes using the quaternion module
and the GALib3D module. They seem to match

Marco Souza de Joode, 2025
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "galib3d.h"
#include "./include/quaternion_math.h"

#define N 1000000
#define EPS 1e-6
#define PI 3.14159265358979


void normalize(struct vector *v) {
    double norm = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
    if (norm > 0) {
        v->x /= norm;
        v->y /= norm;
        v->z /= norm;
    }
}

int main() {
    srand(time(NULL));

    int match_count = 0;
    int mismatch_count = 0;

    for (int i = 0; i < N; i++) {
        // Random vector
        struct vector v = {
            ((double) rand() / RAND_MAX) * 2.0 - 1.0,
            ((double) rand() / RAND_MAX) * 2.0 - 1.0,
            ((double) rand() / RAND_MAX) * 2.0 - 1.0
        };

        // Random axis
        struct vector axis = {
            ((double) rand() / RAND_MAX) * 2.0 - 1.0,
            ((double) rand() / RAND_MAX) * 2.0 - 1.0,
            ((double) rand() / RAND_MAX) * 2.0 - 1.0
        };
        normalize(&axis);

        // Random angle
        double theta = ((double) rand() / RAND_MAX) * 2.0 * PI;

        // --- GA rotation ---
        Multivector mv = {0, v.x, v.y, v.z, 0, 0, 0, 0};
        Multivector B = {0, 0, 0, 0, axis.x, axis.y, axis.z, 0};
        Multivector v_rot_ga = rotate_vector(mv, B, theta);

        // --- Quaternion rotation ---
        struct quaternion q = q_versor(&axis, theta);
        q_normalize(&q);

        struct quaternion v_quat = {0.0, v.x, v.y, v.z};
        struct quaternion q_inv = q_inverse(&q);
        struct quaternion temp = q_mult(&q, &v_quat);
        struct quaternion v_rot_q = q_mult(&temp, &q_inv);

        // Compare results
        double dx = v_rot_ga.b1 - v_rot_q.x;
        double dy = v_rot_ga.b2 - v_rot_q.y;
        double dz = v_rot_ga.b3 - v_rot_q.z;
        double err = sqrt(dx*dx + dy*dy + dz*dz);

        if (err < EPS) {
            match_count++;
        } else {
            mismatch_count++;
            if (mismatch_count < 10) { 
                printf("Mismatch example:\n");
                printf("  GA: (%f, %f, %f)\n", v_rot_ga.b1, v_rot_ga.b2, v_rot_ga.b3);
                printf("  Q : (%f, %f, %f)\n", v_rot_q.x, v_rot_q.y, v_rot_q.z);
                printf("  error = %e\n", err);
            }
        }
    }

    printf("\nTotal tests: %d\n", N);
    printf("Matches: %d\n", match_count);
    printf("Mismatches: %d\n", mismatch_count);

    return 0;
}
