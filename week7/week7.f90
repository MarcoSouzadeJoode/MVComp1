! Marco Souza de Joode, Jan Straub, 2024
! week 7 excersize


program week7
    implicit none
    integer :: N
    real(8), allocatable :: a(:), b(:), c(:)
    real(8), allocatable :: r(:), y(:), r_test(:)
    real(8), allocatable :: x(:),phi(:)
    real(8) :: xmax
    character(len=100) :: filename


    N = 10

    allocate(a(N))
    allocate(b(N))
    allocate(c(N))
    allocate(r(N))
    allocate(y(N))
    allocate(r_test(N))





    call random_number(a)
    call random_number(b)
    call random_number(c)
    call random_number(r)

    a(1) = 0
    c(N) = 0

    ! part (7.1.1)
    call print_vector(a, N)
    call solve_tridiagonal(a, b, c, r, N, y)
    call tridiagonal_times_vector(a, b, c, y, N, r_test)
    print *, "r_test"
    call print_vector(r_test, N)
    call vector_equality(r, r_test, N)

    ! part (7.1.2)

    xmax = 3.d0
    N = 100
    allocate(phi(N))
    allocate(x(N))
    filename = "_.dat"
    call gravity_problem_1d(xmax, N, phi)
    call linspace(-xmax, xmax, N, x)
    call print_potential_to_file(x, phi, N, filename)

    deallocate(phi)
    deallocate(x)


    N = 10000
    allocate(phi(N))
    allocate(x(N))
    filename = "_.dat"
    call gravity_problem_1d(xmax, N, phi)
    call linspace(-xmax, xmax, N, x)
    call print_potential_to_file(x, phi, N, filename)

    deallocate(phi)
    deallocate(x)

    N = 1000000
    allocate(phi(N))
    allocate(x(N))
    filename = "_.dat"
    call gravity_problem_1d(xmax, N, phi)
    call linspace(-xmax, xmax, N, x)
    call print_potential_to_file(x, phi, N, filename)

    deallocate(phi)
    deallocate(x)

    N = 100
    allocate(phi(N))
    allocate(x(N))
    filename = "LU_100.dat"
    call LU_solver(xmax, N, phi)
    call linspace(-xmax, xmax, N, x)
    call print_potential_to_file(x, phi, N, filename)

    deallocate(phi)
    deallocate(x)

    N = 10000
    allocate(phi(N))
    allocate(x(N))
    filename = "LU_10000.dat"
    call LU_solver(xmax, N, phi)
    call linspace(-xmax, xmax, N, x)
    call print_potential_to_file(x, phi, N, filename)

    deallocate(phi)
    deallocate(x)


    ! https://www.netlib.org/lapack/explore-3.1.1-html/dgetrf.f.html

contains



subroutine solve_tridiagonal(a, b, c, r, N, y)
    implicit none
    real(8), intent(in) :: a(:), b(:), c(:), r(:)
    integer, intent(in) :: N
    integer :: i
    real(8), intent(out) :: y(:)

    ! temporary arrays
    real(8), allocatable :: at(:), bt(:), rt(:)

    allocate(at(N))
    allocate(bt(N))
    allocate(rt(N))

    ! FORWARD ELIMINATION
    bt(1) = b(1)
    rt(1) = r(1)

    do i = 2, n
        bt(i) = b(i) - (c(i-1) * a(i))/bt(i-1)
        rt(i) = r(i) - rt(i-1) * (a(i)/bt(i-1))
    end do

    ! BACKWARDS SUBSTITUTION
    y(N) = rt(N) / bt(N)

    do i = N-1, 1, -1
        y(i) = (rt(i) - c(i) * y(i+1)) / bt(i)
    end do

end subroutine solve_tridiagonal


subroutine tridiagonal_times_vector(a, b, c, y, N, r)
implicit none
    real(8), intent(in) :: a(:), b(:), c(:), y(:)
    integer, intent(in) :: N
    integer :: i
    real(8), intent(out) :: r(:)

    r(1) = y(1) * b(1) + y(2) * c(1)

    do i = 2, N - 1
        r(i) = y(i-1) * a(i) + y(i) * b(i) + y(i+1) * c(i)
    end do

    r(N) = y(N-1) * a(N) + y(N) * b(N)

end subroutine tridiagonal_times_vector


subroutine print_vector(x, dim)
implicit none
    real(8), intent(in) :: x(:)
    integer, intent(in) :: dim
    integer :: i

    do i = 1, dim
        write(*,'(E12.4)') x(i)
    end do
end subroutine print_vector

subroutine vector_equality(x, y, dim)
implicit none
    real(8), intent(in) :: x(:), y(:)
    integer, intent(in) :: dim
    integer :: i, err
    real(8) :: eps

    err = 0
    eps = 1e-10

    do i = 1, dim
        if (abs(x(i) - y(i)) > eps) then
            err = err + 1
        end if   
    end do

    if (err > 0) then
        print *, "VECTORS NOT EQUAL"
    else
        print "(A, ES10.2)", "VECTORS EQUAL WITH PRECISION ", eps
    end if

end subroutine vector_equality

subroutine gravity_problem_1d(xmax, N, phi)
implicit none
    REAL(8), INTENT(IN) :: XMAX
    real(8), allocatable :: rho(:), a(:), b(:), c(:), x(:)
    real(8), intent(out) ::  phi(:)
    integer :: i
    integer, intent(in) :: N
    real(8) :: dx, q

    allocate(x(N))
    allocate(rho(N))
    allocate(a(N))
    allocate(b(N))
    allocate(c(N))

    dx = 2*xmax / N

    call linspace(-xmax, xmax, N, x)
    
    q = 1.d0 / dx / dx

    do i = 1, N
        a(i) = q
        b(i) = -2.d0 * q
        c(i) = q

        if (abs(x(i)) <= 1) then
            rho(i) = 1
        else
            rho(i) = 0
        end if
    end do

    call solve_tridiagonal(a, b, c, rho, N, phi)

end subroutine gravity_problem_1d



subroutine LU_solver(xmax, N, phi)
implicit none
    integer, intent(in) :: N
    real(8) ::  rho(n)
    real(8), intent(out) :: phi(:)
    real(8) :: a(n), b(n), c(n), xr(n)
    real(8), intent(in) :: xmax
    real(8) :: X(n,n)
    real(8) :: dx, q

    integer :: i, j, INFO
    integer :: ipiv(n)


    X = 0.d0

    ! populating matrix A

    dx = 2*xmax / N

    call linspace(-xmax, xmax, N, xr)
    
    q = 1.d0 / dx / dx

    do i = 1, N
        a(i) = q
        b(i) = -2.d0 * q
        c(i) = q

        if (abs(xr(i)) <= 1) then
            rho(i) = 1
        else
            rho(i) = 0
        end if
    end do

    do i = 1, n
        X(i, i) = b(i) 
        if (i > 1) then
            X(i, i-1) = a(i-1)
        endif
        if (i < n) then
            X(i, i+1) = c(i) 
        endif
    end do

    print *, "got here"

    ! LAPACK has to be linked .. see makefile

    call DGETRF(N, N, X, N, ipiv, INFO)

    if (INFO /= 0) then
    print *, "DGETRF failed with INFO =", INFO
    stop
    end if

    print *, "got even here"

    ! X is throughput, i.e., populated w L and U now.
 
    phi = rho
    ! initializing phi with rho.
    ! the next call overwrites it with the solution to the system, i.e., the potential

    call DGETRS('N', n, 1, X, n, ipiv, phi, n, info)
    if (info /= 0) then
        print *, "Error: Solution failed with INFO =", info
        stop
    end if

end subroutine LU_solver


subroutine linspace(min, max, N, arr)
implicit none
    real(8), intent(out) :: arr(:)
    real(8), intent(in) :: min, max
    integer, intent(in) :: N
    integer :: i
    real(8) :: range, step

    range = max - min
    step = range / (N-1)

    do i = 1, N
        arr(i) = min + (i-1) * step
    end do
end subroutine linspace

subroutine print_potential_to_file(x, phi, N, filename)
    real(8), intent(in) :: x(:), phi(:)
    integer, intent(in) :: N
    character(len=*), intent(in) :: filename

    integer :: i

    open(1, file=filename, status='replace', action='write')
        write(1, "(A)") "X" // char(9) // "PHI"
        do i = 1, N
            write(1, "(F8.3,A, F8.3)") x(i), char(9), phi(i)
        end do
    close(1)

    print *, "written to ", filename
end subroutine print_potential_to_file



end program week7