## Week 3

task 1 

(a) see src, C code 

(b) see src, C code

(c) figures in plots folder

(d) figures in plots folder

(e) figures in plots folder

(f) The velocity verlet method seems to be performing the best. Then rk4, we should have plotted separatly to see it a bit better. 

Expected the velocity verlet method to perform best since it is a symplectic method. 

task 2 

code in adaptive.jl

Adaptive method based on an energy criterion abs(E-E0)/E0 > epsilon. Different values for epsilon where tested. see folder plots_adaptive, figures ending with +red.png start with h = 0.05. The number before +red is the epsilon value, we see that with higher epsilon value it takes a bit longer before the stepsize decreases, but all methods get stuck with the minimum stepsize of 1e-4 (0.0001). For comparison not_adaptive_0.05+_h.png and not_adaptive_0.0001+_h.png are plots for non adaptive with const stepsize 0.05 and 1e-4. The one starting out with 1e-4 is significantly better and we get stuck in the lower stepsize it is computationally more expensive. If we start with the low stepsize (files with name ending +lowstart.png) it is as if we just fix the step size since our creterion keeps us at the minimum stepsize, so also computational expensive. A compromise is that we change the criterion to abs(E-Eprevious)/Eprevious > epsilon. See plots with ending +prev.png. Since the energy diverges from the initial energy in a chaotic system like this comparing to the previous energy instead can give a better balance between performance and accuracy since the we do not necassarily get stuck with the lower timestep after a while, but still lower the time step at points where the energy difference becomes big. 