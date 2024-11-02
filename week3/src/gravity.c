#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "vec_op.h"

#define PI 3.14159265358979
#define N 3
#define dim 3

typedef struct simulation_parameters {
    double G;
    double h;
    double t;
    double total_T;
    char filename[64];
} sim_param;

typedef struct simulation {
    particle particles[N];
    sim_param parameters;
} simulation;

void init (simulation* S);

void calculate_accelerations (simulation* S) {
    for (int i = 0; i < N; i++) {
    V3 acc = {0.0, 0.0, 0.0};
        for (int j = 0; j < N; j++) {
                
            if (i != j) {
                V3 diff, acc_contribution;
                double dist, prefactor;

                diff = v_subtract(S->particles[j].r, S->particles[i].r);
                dist = v_mag(diff);

                assert(dist > 1e14);

                prefactor = S->parameters.G * S->particles[j].m  / (dist * dist * dist);
                acc_contribution = v_scale(diff, prefactor / S->particles[i].m);
                v_add_ref(&acc, &acc_contribution);
            }
        }
        S->particles[i].a = acc;
    }
}


int main() {
    simulation S;
    init(&S);

    return 0;
};

void init (simulation* S) {

    S->parameters.G = 1.0;
    strcpy(S->parameters.filename, "../runs/tmp");
    S->parameters.t = 0.0;
    S->parameters.total_T = 7.0;
    S->parameters.h = 1e-4;

    S->particles[0].r = (V3){0.0, 4.0, 0.0};
    S->particles[1].r = (V3){-3.0, 0.0, 0.0};
    S->particles[2].r = (V3){0.0, 0.0, 0.0};

    for (int i = 0; i < N; i++) {
        S->particles[i].v = (V3){0.0, 0.0, 0.0};
        S->particles[i].a = (V3){0.0, 0.0, 0.0};
    };



}

