#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 12:52:41 2024

@author: marco
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Loading data
df = pd.read_csv("../runs/tmp", delimiter="\t")

plt.plot(df["x0"] - df["x1"], df["y0"] - df["y1"])
plt.show()

plt.plot(df["x0"], df["y0"])
plt.plot(df["x1"], df["y1"])

plt.show()

E0 = df["EK"][2]+ df["EP"][2]
E = df["EK"]+ df["EP"]


plt.ylim(-100, 100)
plt.plot(df["t"], df["EP"], label="m")
plt.plot(df["t"], df["EK"], label="m")
plt.show()


plt.plot(df["t"], df["EP"]  + df["EK"], label="m")
plt.show()


plt.plot(df["t"], df["L"], label="m")
plt.show()


EK = 0.5 * (3*(  df["vx0"]**2 +  df["vy0"]**2 + df["vz0"]**2)
            +4*(  df["vx1"]**2 +  df["vy1"]**2 + df["vz1"]**2)
            +5*(  df["vx2"]**2 +  df["vy2"]**2 + df["vz2"]**2))


plt.ylim(0,5)
plt.plot(df["t"], EK, label="m")
plt.plot(df["t"], df["EK"], label="m")

plt.show()