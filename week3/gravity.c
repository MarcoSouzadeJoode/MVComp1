#include <stdio.h>
#include <math.h>
#include <string.h>

#include "vec_op.h"

#define PI 3.14159265358979
#define N 3
#define dim 3

typedef struct V3 {
    double x;
    double y;
    double z;
} V3;

typedef struct particle {
    V3 r;
    V3 v;
    V3 a;
    double m;
} particle;

/*
V3 v_add(V3 a, V3 b);
V3 v_subtract(V3 a, V3 b);
V3 v_scale(V3 a, double k);
double v_dot(V3 a, V3 b);
double v_mag(V3 a);
*/

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


int main() {
    return 0;
};

void init (simulation * S) {
    V3 v_null = {0.0, 0.0, 0.0};

    S->parameters.G = 1.0;
    strcpy(S->parameters.filename, "../runs/tmp");
    S->parameters.t = 0.0;
    S->parameters.total_T = 7.0;
    S->parameters.h = 1e-4;



    S->particles[0].r = (V3){0.0, 4.0, 0.0};
    S->particles[1].r = (V3){-3.0, 0.0, 0.0};
    S->particles[2].r = (V3){0.0, 0.0, 0.0};

    for (int i = 0; i < dim; i++) {
        S->particles[i].v = v_null;
        S->particles[0].a = v_null;
    }



}

