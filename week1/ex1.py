#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 08:51:05 2024

@author: marco

int i = 7;
float x = 2*(i/2);
float y = 2*(i/2.);
printf("%e %e \n", y,z);
"""
import numpy as np
import matplotlib.pyplot as plt
i = 7
x = 2*(i/2)
y = 2*(i/2.)
print(x, y)
# prints 7.0 7.0

# equivalent to C behaviour now, performing int division:
    
i = 7
x = 2*(i//2)
y = 2*(i/2.)
print(x, y)


a = 1.0e17;
b = -1.0e17;
c = 1.0;
xx = (a + b) + c;
yy = a + (b + c);

print(xx, yy)

# trying using int arithmetic

a = 1 * 10**17
b = -1 * 10**17
c= 1
xx = (a+b) + c
yy = a + (b+c)

print(xx, yy)


d = 1e20
e = d*d
print("e, d: ",e,d)



### EXCERSIZE: Accumulation of error (14 points)
def f(x):
    return (np.exp(-x) + x -1) / x**2


X_lin = np.linspace(0.1, 0.00001, 1000)
Y_lin = [f(x) for x in X_lin]

plt.plot(X_lin, Y_lin)
plt.xlabel("X")
plt.ylabel("f(X)")
plt.show()


X_log = np.logspace(-3, -6, 10000)

Y_log = [f(x) for x in X_log]
plt.plot(X_log, Y_log, linewidth=0.1, c="k", label="f(X)")
plt.hlines(0.5, 1e-6, 1e-3, color="r", label="1/2", linewidth=1, zorder=-1)
#plt.yscale("log")
plt.legend()
plt.xscale("log")
plt.xlabel("X")
plt.ylabel("f(X)")
plt.show()


## The issue starts occurring around X = 1e-5
## The problem is due to the fact that 1.000000 is much larger than the values
## of exp(-x) in this numerical range, and so information is lost in the 
## finite mantissa of the double precision floats. 

## This can be circumvented by using the numpy np.expm1 (exp - 1) function,
## and we have also written a function, which evaluates the function using
## Taylor expansion - pushing the limit up to the ~ 1e-300 boundary of doubles.

## 


def f_smarter(x):
    return (np.expm1(-x) + x) / x**2

X_log = np.logspace(-5, -20, 10000)
Y_log_sm = [f_smarter(x) for x in X_log]
plt.plot(X_log, Y_log_sm, linewidth=2, c="k", label="f(X)")
plt.hlines(0.5, 1e-5, 1e-9, color="r", label="1/2", linewidth=1, zorder=-1)
#plt.yscale("log")
plt.legend()
plt.xscale("log")
plt.xlabel("X")
plt.ylabel("f(X)")
plt.show()


## Expanding the exp, subtracting and dividing analytically:

def f_taylor(x):
    if x > 1e-2:
        return f_smarter(x)
    else:
        return (0.5 - (1/6) * x - (1/24)*x**2 - (1/120)*x**3 + (1/720)*x**4)

X_log = np.logspace(100, -300, 10000)
Y_log_taylor = [f_taylor(x) for x in X_log]
plt.plot(X_log, Y_log_taylor, linewidth=1, c="k", label="f(X) using taylor")

plt.hlines(0.5, 1e100, 1e-300, color="r", label="1/2", linewidth=2, zorder=-1)
#plt.yscale("log")
plt.legend()
plt.xscale("log")
plt.xlabel("X")
plt.ylabel("f(X)")
plt.show()

## This appears to have solved the problem, as for very low values of x,
## the value 0.5 is returned. 

## However, it is also possible to do this:
def f_sanitized(x):
    if (x > 1e-5) and (x < 1e100):
        return f_smarter(x)
    elif (x >= 1e100):
        return 0
    else:
        return 0.5


X_log_big = np.logspace(280, -300, 10000)
Y_log_san = [f_sanitized(x) for x in X_log_big]
plt.plot(X_log, Y_log_san, linewidth=1, c="k", label="f(X) using an escamotage")
plt.legend(fontsize=8)
plt.xscale("log")
plt.xlabel("X")
plt.ylabel("f(X)")
plt.show()