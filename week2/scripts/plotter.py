#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 20:30:57 2024

@author: marco
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df_rk4 = pd.read_csv("../runs/RUN_RK2", delimiter="\t")
df_rk2 = pd.read_csv("../runs/RUN_RK4", delimiter="\t")
df_m = pd.read_csv("../runs/RUN_midpoint", delimiter="\t")




L1 = 2.0
L2 = 1.0  


x1_m = L1 * np.sin(df_m["phi1"])
y1_m = -L1 * np.cos(df_m["phi1"])  
x2_m = x1_m + L2 * np.sin(df_m["phi2"])
y2_m = y1_m - L2 * np.cos(df_m["phi2"])


x1_rk2 = L1 * np.sin(df_rk2["phi1"])
y1_rk2 = -L1 * np.cos(df_rk2["phi1"])
x2_rk2 = x1_rk2 + L2 * np.sin(df_rk2["phi2"])
y2_rk2 = y1_rk2 - L2 * np.cos(df_rk2["phi2"])


x1_rk4 = L1 * np.sin(df_rk4["phi1"])
y1_rk4 = -L1 * np.cos(df_rk4["phi1"])
x2_rk4 = x1_rk4 + L2 * np.sin(df_rk4["phi2"])
y2_rk4 = y1_rk4 - L2 * np.cos(df_rk4["phi2"])

#plt.ylim(-2*np.pi, 2*np.pi)
plt.plot(df_m["t"], df_m["phi1"], label="phi1 mid")
plt.plot(df_rk2["t"], df_rk2["phi1"], label="phi1 RK2")
plt.plot(df_rk4["t"], df_rk4["phi1"], label="phi1 RK4")
plt.legend()
plt.show()

# ENERGY
plt.plot(df_m["t"], (df_m["E"] - df_m["E"][0] )/ df_m["E"][0] , label="E mid")
plt.plot(df_rk2["t"], (df_rk2["E"] - df_rk2["E"][0] )/ df_rk2["E"][0], label="E RK2")
plt.plot(df_rk4["t"],(df_rk4["E"] - df_rk4["E"][0] )/ df_rk4["E"][0], label="E RK4")
plt.legend()
plt.show()

plt.plot(df_m["t"], df_m["phi2"],label="phi2 mid")
plt.plot(df_rk2["t"], df_rk2["phi2"],label="phi2 RK2")
plt.plot(df_rk4["t"], df_rk4["phi2"],label="phi2 RK4")
plt.legend()
plt.show()

plt.plot(df_m["t"], df_m["q1"], label="q1 mid")
plt.plot(df_rk2["t"], df_rk2["q1"], label="q1 RK2")
plt.plot(df_rk4["t"], df_rk4["q1"], label="q1 RK4")
plt.legend()
plt.show()


plt.plot(df_m["t"], df_m["q2"],label="q2 mid")
plt.plot(df_rk2["t"], df_rk2["q2"],label="q2 RK2")
plt.plot(df_rk4["t"], df_rk4["q2"],label="q2 RK4")
plt.legend()
plt.show()



plt.plot(df_rk4["phi2"], df_rk4["q2"],label="q2 RK4")
plt.show()

"""


"""




#plt.plot(x1_m, y1_m, label="PENDULUM 1", c="b")
plt.plot(x2_m, y2_m, c="b", label="mid")

#plt.plot(x1_rk2, y1_rk2, label="PENDULUM 1 RK2", c="k")
plt.plot(x2_rk2, y2_rk2, c="k",label="rk2")

#plt.plot(x1_rk4, y1_rk4, label="PENDULUM 1", c="g")
plt.plot(x2_rk4, y2_rk4, c="r", label="rk4")
plt.legend()

plt.show()


def lines(x1, x2, y1, y2, N):
    k = (y2-y1)/(x2-x1)
    q = y1 - k*x1
    
    xsp = np.linspace(x1, x2, N)
    ysp = k*xsp + q
    return (xsp, ysp)


"""
for i, _ in enumerate(x1[:]):
    plt.xlim(-3.5, 3.5)
    plt.ylim(-3.5, 3.5)
    
    
    plt.scatter(0, 0, c="k")
    plt.scatter(x1[i], y1[i], s=50,c="r",label="PENDULUM 1",zorder=10)
    plt.scatter(x2[i], y2[i], s=150,c="b",label="PENDULUM 2",zorder=10)
    
    xspA, yspA = lines(0, x1[i], 0, y1[i], 100)
    plt.plot(xspA, yspA, c="k")
    
    xspB, yspB = lines(x1[i], x2[i], y1[i], y2[i], 100)
    plt.plot(xspB, yspB, c="blue",zorder=-1)
    plt.plot(xspA, yspA, c="k",zorder=-5)
    
    plt.savefig(f"anim/frame_{i:05}.jpg", dpi=100)
    plt.clf()
"""
    