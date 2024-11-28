#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 22:21:21 2024

@author: marco
"""

import numpy as np
import matplotlib.pyplot as plt

# Data
theta_crit = np.array([0.2, 0.4, 0.6, 0.8])
tree_time = np.array([5.798673391342163, 2.172332286834717, 1.021904468536377, 0.6720137596130371])
summation_time = np.array([2.13323974609375, 2.1668589115142822, 2.1559505462646484, 2.2780115604400635])
mean_error = np.array([0.0004, 0.0031, 0.0096, 0.0222])

# Plotting
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot times
ax1.plot(theta_crit, tree_time, 'o-', label='Tree Time (s)', color='blue')
ax1.set_xlabel('Theta Crit')
ax1.set_ylabel('Execution Time (s)', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.legend(loc='upper left')
ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

# Add a secondary y-axis for mean error
ax2 = ax1.twinx()
ax2.plot(theta_crit, mean_error, 's-', label='Mean Error', color='red')
ax2.set_ylabel('Mean Error', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax2.legend(loc='upper right')

# Title and layout
plt.title('Execution Times and Mean Error vs Theta Crit')
plt.tight_layout()
plt.show()
