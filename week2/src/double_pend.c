#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979

struct sim_param {
    double g;
    double h;
    double t;
    double total_T;
    char filename[64];
};

struct pendulum {
    double m1;
    double m2;
    double L1;
    double L2;
    double F;
};

struct state_vector {
    double phi1;
    double phi2;
    double q1;
    double q2;
};

struct state_change {
    double phi1_dot;
    double phi2_dot;
    double q1_dot;
    double q2_dot;
};


void init(struct sim_param* sim, struct pendulum* pend, struct state_vector* X) {
    strncpy(sim->filename, "../runs/tmp", sizeof(sim->filename) - 1);
    sim->filename[sizeof(sim->filename) - 1] = '\0';

    sim->g = 1.0;
    sim->t = 0.0;
    sim->h = 0.05;
    sim->total_T = 100;

    pend->L1 = 2.0;
    pend->L2 = 1.0;
    pend->m1 = 0.5;
    pend->m2 = 1.0;

    //constant a term often used
    pend->F = pend->m1 * pend->L1*pend->L1 + pend->m2 * pend->L2*pend->L2;

    X->phi1 = 50.0 * PI / 180.0;
    X->phi2 = -120.0 * PI / 180.0;

    X->q1 = 0.0;
    X->q2 = 0.0;
}



void slope_of_state_vector(struct sim_param* sim, struct pendulum* p, struct state_vector* X, struct state_change *K) {
    double phi1_dot, phi2_dot;
    double cos_dif, sin_dif, detM;

    cos_dif = cos(X->phi1 - X->phi2);
    sin_dif = sin(X->phi1 - X->phi2);

    /*
    detM = p->L1 * p->L2 * p->m2 * (p->m1+p->m2)
        - (p->m2 * p->L1 * p->L2 * cos_dif)*(p->m2 * p->L1 * p->L2 * cos_dif);
    */

    detM = p->L1 * p->L1 * p->L2 * p->L2 * p->m2 * (p->m1 + p->m2) 
        - (p->m2 * p->L1 * p->L2 * cos_dif) * (p->m2 * p->L1 * p->L2 * cos_dif);


    //K->phi1_dot = (X->q1*p->m2*p->L2*p->L2 - X->q2 * p->m2 * p->L1 * p->L2 * cos_dif) / detM;
    //K->phi2_dot = ((- X->q1 * p->m2 * p->L1 * p->L2 * cos_dif) + X->q2 * p->L1*p->L1* (p->m1*p->m2)) / detM;

    K->phi1_dot = (X->q1 * p->m2 * p->L2 * p->L2 - X->q2 * p->m2 * p->L1 * p->L2 * cos_dif) / detM;
    K->phi2_dot = (-X->q1 * p->m2 * p->L1 * p->L2 * cos_dif + X->q2 * p->L1 * p->L1 * (p->m1 + p->m2)) / detM;

    K->q1_dot = - p->m2 * p->L1 * p->L2 * K->phi1_dot * K->phi2_dot*sin_dif 
        - sim->g * p->L1 * sin(X->phi1) * (p->m1 + p->m2);
    K->q2_dot = +  p->m2 * p->L1 * p->L2 * K->phi1_dot * K->phi2_dot*sin_dif
        - sim->g * p->L2 * sin(X->phi2) * (p->m1 + p->m2);

}



void slope_of_state_vector_messed_up(struct sim_param* sim, struct pendulum* p, struct state_vector* X, struct state_change *K) {
    double phi1_dot, phi2_dot;
    double cos_dif, sin_dif, detM, A, B, pdpd;

    cos_dif = cos(X->phi1 - X->phi2);
    sin_dif = sin(X->phi1 - X->phi2);
    A = p->m2 * p->L1 * p->L2 * cos_dif;
    B = p->m2 * p->L2 * p->L2;

    /*
    detM = p->L1 * p->L2 * p->m2 * (p->m1+p->m2)
        - (p->m2 * p->L1 * p->L2 * cos_dif)*(p->m2 * p->L1 * p->L2 * cos_dif);
    */

    detM =  B * p->F - A * A;
    printf("%f\n",detM);

    // phi_1 dot times phi_2 dot
    pdpd = (A*B * X->q1 * X->q1 - A*p->F * X->q2 * X->q2 + (A*A - B*p->F)*X->q1*X->q2) / detM / detM;

    K->phi1_dot = X->q1 / p->F - (A/B) * X->q2;
    K->phi2_dot = (-A * X->q1 + p->F * X->q2)/detM;

    //pdpd = K->phi1_dot * K->phi2_dot;
    
    K->q1_dot = -2 * p->m2 * p->L1 * p->L2 * pdpd * sin_dif 
        - sim->g * p->L1 * sin(X->phi1) * (p->m1 + p->m2);
    
    K->q2_dot = +2 * p->m2 * p->L1 * p->L2 * pdpd * sin_dif
        - p->m2 * sim->g * p->L2 * sin(X->phi2);

}


void midpoint(struct sim_param* sim, struct pendulum* pend, struct state_vector *X, struct state_change *K) {
    struct state_vector tmp;
    struct state_change k1, k2;

    slope_of_state_vector(sim, pend, X, &k1);

    // temp = X + kh/2
    tmp.phi1 = X->phi1 + 0.5 * sim->h * k1.phi1_dot;
    tmp.phi2 = X->phi2 + 0.5 * sim->h * k1.phi2_dot;
    tmp.q1 = X->q1 + 0.5 * sim->h * k1.q1_dot;
    tmp.q2 = X->q2 + 0.5 * sim->h * k1.q2_dot;

    slope_of_state_vector(sim, pend, &tmp, &k2);

    K->phi1_dot = k2.phi1_dot;
    K->phi2_dot = k2.phi2_dot;
    K->q1_dot = k2.q1_dot;
    K->q2_dot = k2.q2_dot;
}


void RK2(struct sim_param* sim, struct pendulum* pend, struct state_vector *X, struct state_change *K) {
    struct state_vector tmp;
    struct state_change k1, k2;

    slope_of_state_vector(sim, pend, X, &k1);

    // temp = X + kh
    tmp.phi1 = X->phi1 +  sim->h * k1.phi1_dot;
    tmp.phi2 = X->phi2 +  sim->h * k1.phi2_dot;
    tmp.q1 = X->q1 +  sim->h * k1.q1_dot;
    tmp.q2 = X->q2 +   sim->h * k1.q2_dot;

    slope_of_state_vector(sim, pend, &tmp, &k2);

    K->phi1_dot = 0.5 * (k1.phi1_dot + k2.phi1_dot);
    K->phi2_dot = 0.5 * (k1.phi2_dot + k2.phi2_dot);
    K->q1_dot = 0.5 * (k1.q1_dot + k2.q1_dot);
    K->q2_dot = 0.5 * (k1.q2_dot + k2.q2_dot);
}


void RK4(struct sim_param* sim, struct pendulum* pend, struct state_vector *X, struct state_change *K) {
    struct state_vector tmp, tmp2, tmp3;
    struct state_change k1, k2, k3, k4;

    slope_of_state_vector(sim, pend, X, &k1);

    tmp.phi1 = X->phi1 + 0.5 * sim->h * k1.phi1_dot;
    tmp.phi2 = X->phi2 + 0.5 * sim->h * k1.phi2_dot;
    tmp.q1 = X->q1 + 0.5 * sim->h * k1.q1_dot;
    tmp.q2 = X->q2 + 0.5 * sim->h * k1.q2_dot;

    slope_of_state_vector(sim, pend, &tmp, &k2);

    tmp2.phi1 = X->phi1 + 0.5 * sim->h * k2.phi1_dot;
    tmp2.phi2 = X->phi2 + 0.5 * sim->h * k2.phi2_dot;
    tmp2.q1 = X->q1 + 0.5 * sim->h * k2.q1_dot;
    tmp2.q2 = X->q2 + 0.5 * sim->h * k2.q2_dot;

    slope_of_state_vector(sim, pend, &tmp2, &k3);

    tmp3.phi1 = X->phi1 + sim->h * k3.phi1_dot;
    tmp3.phi2 = X->phi2 +  sim->h * k3.phi2_dot;
    tmp3.q1 = X->q1 +  sim->h * k3.q1_dot;
    tmp3.q2 = X->q2 + sim->h * k3.q2_dot;

    slope_of_state_vector(sim, pend, &tmp3, &k4);

    K->phi1_dot = 1.0/6.0 * (k1.phi1_dot+2*k2.phi1_dot + 2*k3.phi1_dot + k4.phi1_dot);
    K->phi2_dot = 1.0/6.0 * (k1.phi2_dot+2*k2.phi2_dot + 2*k3.phi2_dot + k4.phi2_dot);
    K->q1_dot =  1.0/6.0 * (k1.q1_dot+2*k2.q1_dot + 2*k3.q1_dot + k4.q1_dot);
    K->q2_dot =  1.0/6.0 * (k1.q2_dot+2*k2.q2_dot + 2*k3.q2_dot + k4.q2_dot);
}


double energy(struct sim_param* sim, struct pendulum* p,struct state_vector *X, struct state_change *K) {
    double T, V;

    T = 0.5 * p->m1 * (p->L1 * K->phi1_dot)*(p->L1 * K->phi1_dot)
        + 0.5 * p->m2 * ((p->L1 * K->phi1_dot)*(p->L1 * K->phi1_dot) + (p->L2 * K->phi2_dot)*(p->L2 * K->phi2_dot)
        + 2 * p->L1 * p->L2 * K->phi1_dot * K->phi2_dot * cos(X->phi1 - X->phi2));
    
    //V = p->m1 * sim->g * p->L1 * (1 - cos(X->phi1))
    //    + p->m2 * sim->g*(p->L1 *(1 - cos(X->phi1)) + p->L2*(1 - cos(X->phi2)) );

    V = p->m1 * sim->g * p->L1 * (1 - cos(X->phi1))
    + p->m2 * sim->g * (p->L1 * (1 - cos(X->phi1)) + p->L2 * (1 - cos(X->phi2)));

    
    return T + V;
}

void pr_file(struct sim_param* sim, struct state_vector *X, double E, FILE *file) {
    fprintf(file, "%+e\t%+e\t%+e\t%+e\t%+e\t%+e\n",sim->t,X->phi1,X->phi2,X->q1,X->q2,E);
}

void simulation(struct sim_param* sim, struct pendulum* pend, struct state_vector *X) {
    struct state_change K;
    double E;

    FILE *file = fopen(sim->filename, "w");
    fprintf(file, "t\tphi1\tphi2\tq1\tq2\tE\n");

    while (sim->t < sim->total_T) {
        midpoint(sim, pend, X, &K);
        E = energy(sim, pend, X, &K);
        pr_file(sim, X, E, file);
        printf("E  = %f\n", E);

        X->phi1 += sim->h * K.phi1_dot;
        X->phi2 += sim->h * K.phi2_dot;
        X->q1 += sim->h * K.q1_dot;
        X->q2 += sim->h * K.q2_dot;


        sim->t += sim->h;
    }

    fclose(file);
}


int main (){
    struct sim_param sim;
    struct pendulum pend;
    struct state_vector X;

    init(&sim, &pend, &X);
    simulation(&sim, &pend, &X);
}