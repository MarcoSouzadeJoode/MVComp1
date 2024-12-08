# README week7
# Marco Souza de Joode, Jan Straub

7.1.1 : Utilizing the liberty offered by the footnote (choosing Fortran):
solved using

subroutine solve_tridiagonal(a, b, c, r, N, y)

in file week7.f90

Tested in

subroutine tridiagonal_times_vector(a, b, c, y, N, r)


7.1.2a : Solved in week7.f90, in subroutine

subroutine gravity_problem_1d(xmax, N, phi)

utilizing the previous subroutine.

7.1.2b : Overplotted, see figure potentials_1d.png


7.1.2c : See same figure. the potential is very similar, the BC is obeyed
even closer.


7.2.1,2 : Solved in subroutine LU_solver(xmax, N, phi). Instead of scipy,
the standard LAPACK implementation has been used, i.e., calling the subroutines

DGETRF(..)
DGETRS(..)

Documentation e.g., here: 
https://netlib.org/lapack/explore-html-3.6.1/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html

It takes noticeably longer.

7.3: see Jupyter notebooks ex731 and ex732, and the included figures.



