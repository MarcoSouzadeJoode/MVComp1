import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Data
npart = np.array([500, 1000, 2000, 4000])
tree_time = np.array([0.4194779396057129, 1.072082281112671, 2.6269755363464355, 6.39818549156189])
summation_time = np.array([0.567838191986084, 2.2728068828582764, 9.066545009613037, 37.03518319129944])
mean_error = np.array([0.0126, 0.0096, 0.0078, 0.0070])

# Logarithmic transformation for regression
log_npart = np.log10(npart)
log_tree_time = np.log10(tree_time)
log_summation_time = np.log10(summation_time)

# Fit regression lines
slope_tree, intercept_tree, _, _, _ = linregress(log_npart, log_tree_time)
slope_sum, intercept_sum, _, _, _ = linregress(log_npart, log_summation_time)

# Generate fitted lines
fitted_tree = 10**(slope_tree * log_npart + intercept_tree)
fitted_sum = 10**(slope_sum * log_npart + intercept_sum)

# Plot execution times with log-log axes
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))

ax1.loglog(npart, tree_time, 'o-', label='Tree Time')
ax1.loglog(npart, summation_time, 'o-', label='Direct Summation Time')
ax1.loglog(npart, fitted_tree, '--', label=f'Tree Fit (slope={slope_tree:.2f})')
ax1.loglog(npart, fitted_sum, '--', label=f'Summation Fit (slope={slope_sum:.2f})')

ax1.set_xlabel('Number of Particles (NPART)')
ax1.set_ylabel('Execution Time (s)')
ax1.set_title('Execution Time and Mean Error vs NPART (Log-Log Scale)')
ax1.legend(loc='upper left')
ax1.grid(True, which="both", linestyle='--', linewidth=0.5)

# Add a secondary y-axis for the mean error
ax2.plot(npart, mean_error, 's-', color='red', label='Mean Error')
ax2.set_ylabel('mean Error')
ax2.set_xscale("log")
ax2.set_yscale("log")

ax2.legend()
ax2.grid(True, which="both", linestyle='--', linewidth=0.5)

# Tight layout
plt.tight_layout()
plt.show()