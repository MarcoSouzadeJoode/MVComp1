#include <stdio.h>


int main() {

// x is done using int division
// y is done float division: i/2 = 3
int i = 7;
float x = 2*(i/2);
float y = 2*(i/2.);
printf("%e %e \n", x,y);



// very big numbers!
// we should add / subtract numbers of similar magnitudes,
// otherwise numerical precision is lost!!

// xx is correct
// yy is wrong!

// in python (see attached .py script) the behaviour
// is correct in both cases. Integer division must be said explicitly 
// using the // operator.

float a, b, c;
a = 1.0e17;
b = -1.0e17;
c = 1.0;

float xx, yy;
xx = (a + b) + c;
yy = a + (b + c);

printf("xx, yy: %e %e", xx, yy);


// we overflow the float type (1e40 > 1e38)
float d = 1e20;
float e;
e = d*d;
printf("d, e : %e %e\n", d,e);

//we can use doubles then, and the issue is solved:
// this is also what python does.

printf("Using double precision\n");
double f = 1e20;
double g;
g = f*f;
printf("f, g : %e %e\n", f, g);


}