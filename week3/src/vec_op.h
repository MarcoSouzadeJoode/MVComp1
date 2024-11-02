#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979



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


V3 v_add(V3 a, V3 b);
V3 v_subtract(V3 a, V3 b);
V3 v_scale(V3 a, double k);
double v_dot(V3 a, V3 b);
double v_mag(V3 a);


V3 v_null = {0.0, 0.0, 0.0};


V3 v_add(V3 a, V3 b) {
    return (V3){a.x + b.x, a.y + b.y, a.z + b.z};
}

V3 v_add_ref(V3 *a, V3* b) {
    return (V3){a->x + b->x, a->y + b->y, a->z + b->z};
}

V3 v_subtract(V3 a, V3 b) {
    return (V3){a.x - b.x, a.y - b.y, a.z - b.z};
}

V3 v_scale(V3 a, double k) {
    return (V3){k * a.x, k*a.y, k*a.z};
}

double v_dot(V3 a, V3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double v_mag(V3 a) {
    return v_dot(a, a);
}



