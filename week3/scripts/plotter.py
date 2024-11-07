import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Loading data
df_m = pd.read_csv("../runs/midpoint", delimiter="\t")
df_v = pd.read_csv("../runs/verlet", delimiter="\t")
df_2 = pd.read_csv("../runs/rk2", delimiter="\t")
df_4 = pd.read_csv("../runs/rk4", delimiter="\t")

# Define color map for consistency across plots
colors = {
    'midpoint': '#1f77b4',
    'verlet': '#aec7e8',
    'rk2': '#ff7f0e',
    'rk4': '#2ca02c'
}

# Position plots
def plot_positions(df, title, color):
    plt.figure()
    plt.plot(df["x0"], df["y0"], label="Particle 0", color=color, linewidth=1)
    plt.plot(df["x1"], df["y1"], label="Particle 1", color=color, linestyle="--", linewidth=1)
    plt.plot(df["x2"], df["y2"], label="Particle 2", color=color, linestyle=":", linewidth=1)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Positions - {title}")
    plt.legend()
    plt.grid()
    plt.show()

# Plotting positions for each method
plot_positions(df_m, "Midpoint", colors['midpoint'])
plot_positions(df_v, "Verlet", colors['verlet'])
plot_positions(df_2, "RK2", colors['rk2'])
plot_positions(df_4, "RK4", colors['rk4'])

# Energy plot (all methods on one plot)
plt.figure()
plt.plot(df_m["t"], df_m["EK"] + df_m["EP"], label="Midpoint", color=colors['midpoint'])
plt.plot(df_v["t"], df_v["EK"] + df_v["EP"], label="Verlet", color=colors['verlet'])
plt.plot(df_2["t"], df_2["EK"] + df_2["EP"], label="RK2", color=colors['rk2'])
plt.plot(df_4["t"], df_4["EK"] + df_4["EP"], label="RK4", color=colors['rk4'])
plt.xlabel("Time")
plt.ylabel("Total Energy (EK + EP)")
plt.title("Total Energy over Time for Different Methods")
plt.legend()
plt.grid()
plt.show()

# Angular Momentum plot (all methods on one plot)
plt.figure()
plt.plot(df_m["t"], df_m["L"], label="Midpoint", color=colors['midpoint'])
plt.plot(df_v["t"], df_v["L"], label="Verlet", color=colors['verlet'])
plt.plot(df_2["t"], df_2["L"], label="RK2", color=colors['rk2'])
plt.plot(df_4["t"], df_4["L"], label="RK4", color=colors['rk4'])
plt.xlabel("Time")
plt.ylabel("Angular Momentum (L)")
plt.title("Angular Momentum over Time for Different Methods")
plt.legend()
plt.grid()
plt.show()
