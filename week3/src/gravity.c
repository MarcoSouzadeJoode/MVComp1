#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "vec_op.h"

#define PI 3.14159265358979
#define N 3
#define dim 3

typedef struct particle {
    V3 r;
    V3 v;
    V3 a;
    double m;
} particle;

typedef struct simulation_parameters {
    double G;
    double h;
    double t;
    double total_T;
    double kinetic_energy;
    double potential_energy;
    double angular_momentum;
    char filename[64];
} sim_param;

typedef struct simulation {
    particle particles[N];
    sim_param parameters;
} simulation;

void init (simulation* S);
void pr_header(FILE *file);
void pr_step(simulation *S, FILE *file);

void calculate_accelerations (simulation* S) {
    for (int i = 0; i < N; i++) {
    V3 acc = {0.0, 0.0, 0.0};
        for (int j = 0; j < N; j++) {
                
            if (i != j) {
                V3 diff, acc_contribution;
                double dist, prefactor;

                diff = v_subtract(S->particles[i].r, S->particles[j].r);
                dist = v_mag(diff);

                assert(dist > 1e-14);

                prefactor = - S->parameters.G * S->particles[j].m  / (dist * dist * dist);
                acc_contribution = v_scale(diff, prefactor / S->particles[i].m);
                //printf("%f\n", acc_contribution.x);
                acc = v_add(acc, acc_contribution);
            }
        }
        S->particles[i].a = acc;
    }
}




void init (simulation* S) {

    S->parameters.G = 1.0;
    strcpy(S->parameters.filename, "runs/tmp");
    S->parameters.t = 0.0;
    S->parameters.total_T = 7.0;
    S->parameters.h = 1e-4;

    S->particles[0].r = (V3){0.0, 4.0, 0.0};
    S->particles[1].r = (V3){-3.0, 0.0, 0.0};
    S->particles[2].r = (V3){0.0, 0.0, 0.0};

    S->particles[0].m = 3.0;
    S->particles[1].m = 4.0;
    S->particles[2].m = 5.0;



    for (int i = 0; i < N; i++) {
        S->particles[i].v = (V3){0.0, 0.0, 0.0};
        S->particles[i].a = (V3){0.0, 0.0, 0.0};
    };



}

void velocity_verlet(simulation *S) {
    particle tmp[N];

    
    for (int i = 0; i < N; i++) {
        tmp[i] = S->particles[i];

        S->particles[i].r = v_add(S->particles[i].r,
                                v_add(v_scale(S->particles[i].v,S->parameters.h),
                                v_scale(S->particles[i].a,S->parameters.h * S->parameters.h * 0.5)));
    }

    calculate_accelerations(S);

    for (int i = 0; i < N; i++) {
        S->particles[i].v = v_add(tmp[i].v, v_scale(v_add(tmp[i].a, S->particles[i].a), 0.5 * S->parameters.h));    
    }
    

}

void RK2_step(simulation* S) {
    particle k1h[N], k2h[N], tmp[N];

    calculate_accelerations(S);


    // storing for later restoration
    for (int i = 0; i < N; i++) {
    tmp[i].r = S->particles[i].r;
    tmp[i].v = S->particles[i].v;
    }


    // first substep
    for (int i = 0; i < N; i++) {
    k1h[i].a = v_scale(S->particles[i].a, S->parameters.h);
    k1h[i].v = v_scale(S->particles[i].v, S->parameters.h);

    S->particles[i].r = v_add(S->particles[i].r, k1h[i].v);
    S->particles[i].v = v_add(S->particles[i].v, k1h[i].a);

    }

    calculate_accelerations(S);

    // second substep: based on new accelerations
    for (int i = 0; i < N; i++) {
    k2h[i].a = v_scale(S->particles[i].a, S->parameters.h);
    k2h[i].v = v_scale(S->particles[i].v, S->parameters.h);
    }

    // re-storing original positions & velocities
    for (int i = 0; i < N; i++) {
    S->particles[i].r = v_add(tmp[i].r, v_scale(v_add(k1h[i].v, k2h[i].v), 0.5));
    S->particles[i].v = v_add(tmp[i].v, v_scale(v_add(k1h[i].a, k2h[i].a), 0.5));
    }

}

void midpoint_step(simulation* S) {
    particle tmp[N];
    double h = S->parameters.h;

    calculate_accelerations(S);

    // save initial state
    for (int i = 0; i < N; i++) {
        tmp[i] = S->particles[i];
    }

    for (int i = 0; i < N; i++) {
        S->particles[i].v = v_add(S->particles[i].v, v_scale(S->particles[i].a, 0.5 * h));
        S->particles[i].r = v_add(S->particles[i].r, v_scale(S->particles[i].v, 0.5 * h));
    }

    calculate_accelerations(S);

    // update r, v using midpoint derivative values

    for (int i = 0; i < N; i++) {
        S->particles[i].v = v_add(tmp[i].v, v_scale(S->particles[i].a, h));
        S->particles[i].r = v_add(tmp[i].r, v_scale(S->particles[i].v, h));
    }


}


void RK4_step(simulation* S) {
    particle k1[N], k2[N], k3[N], k4[N], tmp[N];
    particle two_k2[N], two_k3[N];

    double h = S->parameters.h;
    // storing for later restoration
    for (int i = 0; i < N; i++) {
    tmp[i].r = S->particles[i].r;
    tmp[i].v = S->particles[i].v;
    }

    calculate_accelerations(S);

    // first substep
    for (int i = 0; i < N; i++) {
    k1[i].a = S->particles[i].a;
    k1[i].v = S->particles[i].v;

    S->particles[i].r = v_add(tmp[i].r, v_scale(k1[i].v, 0.5*h));
    S->particles[i].v = v_add(tmp[i].v, v_scale(k1[i].a, 0.5*h));

    }

    calculate_accelerations(S);

    // second substep: based on new accelerations
    for (int i = 0; i < N; i++) {
    k2[i].a = S->particles[i].a;
    k2[i].v = S->particles[i].v;

    S->particles[i].r = v_add(tmp[i].r, v_scale(k2[i].v, 0.5*h));
    S->particles[i].v = v_add(tmp[i].v, v_scale(k2[i].a, 0.5*h));
    }

    calculate_accelerations(S);

    // third substep: based on new accelerations
    for (int i = 0; i < N; i++) {
    k3[i].a = S->particles[i].a;
    k3[i].v = S->particles[i].v;

    S->particles[i].r = v_add(tmp[i].r, v_scale(k3[i].v, h));
    S->particles[i].v = v_add(tmp[i].v, v_scale(k3[i].a, h));
    }

    calculate_accelerations(S);


    for (int i = 0; i < N; i++) {
    // calculating last slope k4:
    k4[i].a = S->particles[i].a;
    k4[i].v = S->particles[i].v;

    two_k2[i].v = v_scale(k2[i].v, 2.0);
    two_k2[i].a = v_scale(k2[i].a, 2.0);

    two_k3[i].v = v_scale(k3[i].v, 2.0);
    two_k3[i].a = v_scale(k3[i].a, 2.0);


    S->particles[i].r = v_add(tmp[i].r, 
    v_scale(v_add(k1[i].v, v_add(two_k2[i].v, v_add(two_k3[i].v, k4[i].v)))
    , h / 6.0));
    
    S->particles[i].v = v_add(tmp[i].v, 
    v_scale(v_add(k1[i].a, v_add(two_k2[i].a, v_add(two_k3[i].a, k4[i].a)))
    , h / 6.0));

    //printf("updated\n");

    }

}

void energy(simulation* S) {
    double EK = 0, EP = 0;

    for (int i = 0; i < N; i++) {
        EK += 0.5 * S->particles[i].m * v_dot(S->particles[i].v, S->particles[i].v);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j > i) {
                EP += S->parameters.G * S->particles[i].m * S->particles[j].m /
                 (v_mag(v_subtract(S->particles[i].r, S->particles[j].r)));
            }
        }
    }

    S->parameters.kinetic_energy = EK;
    S->parameters.potential_energy = EP;
}

void angular_momentum(simulation* S) {
    double L = 0;
    for (int i = 0; i < N; i++) {
        L+= S->particles[i].m * v_mag( v_cross (S->particles[i].r, S->particles[i].v));
    }
    S->parameters.angular_momentum = L;
}


void simulate(simulation* S) {

    FILE *file = fopen(S->parameters.filename, "w");

    pr_header(file);
    
    while (S->parameters.t < S->parameters.total_T) {
        pr_step(S, file);
        midpoint_step(S);
        energy(S);
        angular_momentum(S);
        S->parameters.t += S->parameters.h;
    }

    fclose(file);

}

void pr_step(simulation *S, FILE *file) {
    fprintf(file, "%+e\t", S->parameters.t);
    for (int i = 0; i < N; i++) {
        fprintf(file, "%+e\t%+e\t%+e\t", S->particles[i].r.x, S->particles[i].r.y, S->particles[i].r.z);
        fprintf(file, "%+e\t%+e\t%+e\t", S->particles[i].v.x, S->particles[i].v.y, S->particles[i].v.z);
        fprintf(file, "%+e\t%+e\t%+e\t", S->particles[i].a.x, S->particles[i].a.y, S->particles[i].a.z);
    }
    fprintf(file, "%+e\t%+e\t%+e\n", S->parameters.kinetic_energy, S->parameters.potential_energy, S->parameters.angular_momentum);
}

void pr_header(FILE *file) {
    fprintf(file, "t\t");
    for (int i = 0; i < N; i++) {
        fprintf(file, "x%d\ty%d\tz%d\t", i, i, i);
        fprintf(file, "vx%d\tvy%d\tvz%d\t", i, i, i);
        fprintf(file, "ax%d\tay%d\taz%d\t", i, i, i);
    }
    fprintf(file, "EK\tEP\tL\n");
}

int main() {
    simulation S;
    init(&S);
    simulate(&S);

    return 0;
};