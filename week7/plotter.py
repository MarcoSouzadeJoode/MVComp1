#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 23:01:57 2024

@author: marco
"""

# Load the data from the uploaded file and plot it using Python
import matplotlib.pyplot as plt
import numpy as np

# Load the data file path
file_path = 'gravity_100.dat'
file_path2 = 'gravity_10000.dat'
file_path3 = 'LU_100.dat'
file_path4 = 'LU_10000.dat'

# Load the data, skipping the header
data = np.loadtxt(file_path, skiprows=1)
data2 = np.loadtxt(file_path2, skiprows=1)
data3 = np.loadtxt(file_path3, skiprows=1)
data4 = np.loadtxt(file_path4, skiprows=1)

# Extract columns
x = data[:, 0]
phi = data[:, 1]

x2 = data2[:, 0]
phi2 = data2[:, 1]

x3 = data3[:, 0]
phi3 = data3[:, 1]


x4 = data4[:, 0]
phi4 = data4[:, 1]

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(x, phi, label=r'$\Phi(x)$, 100', color='blue', linewidth=2)
plt.plot(x3, phi3, label=r'$\Phi(x)$, LU decomp, 100', color='blue', linewidth=5, alpha=0.3)
plt.plot(x2, phi2, label=r'$\Phi(x)$, 10000', color='red', linewidth=2)
plt.plot(x4, phi4, label=r'$\Phi(x)$, LU decomp, 10000', color='red', linewidth=5, alpha=0.3)


# Label the axes and the plot
plt.title("1D Gravity Potential", fontsize=16)
plt.xlabel("X (Position)", fontsize=14)
plt.ylabel(r"$\Phi$ (Potential)", fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
plt.legend(fontsize=12)
plt.tight_layout()

# Show the plot
plt.show()
