#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979



typedef struct V3 {
    double x;
    double y;
    double z;
} V3;




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

V3 v_cross(V3 a, V3 b) {
    V3 res;
    res.x = a.y * b.z - a.z * b.y;
    res.y = a.z * b.x - a.x * b.z;
    res.z = a.x * b.y - a.y * b.x;
    return res;
}

double v_dot(V3 a, V3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double v_mag(V3 a) {
    return sqrt(v_dot(a, a));
}





