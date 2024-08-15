program pe_simulation
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer :: nz, nr, ir, i, j
    real(dp) :: zs, f, omega, eps, c0, d, h, h2, k0, rmax, deltar, fac
    real(dp), allocatable :: z(:), r(:), c(:)
    complex(dp), allocatable :: psi(:,:)
    complex(dp), allocatable :: A(:,:), B(:,:), C(:,:), L(:,:), U(:,:), y1(:), y(:)
    real(dp) :: pi = 4.0_dp * atan(1.0_dp)
    complex(dp) :: ii = (0.0_dp, 1.0_dp)

    ! *** set up Munk SSP ***
    zs = 1000.0_dp
    f = 50.0_dp
    omega = 2.0_dp * pi * f
    eps = 0.00737_dp
    c0 = 1500.0_dp
    d = 5000.0_dp
    nz = 250
    h = d / nz
    h2 = h * h
    allocate(z(nz))
    do i = 1, nz
        z(i) = (i - 1) * h
    end do

    allocate(c(nz))
    do i = 1, nz
        c(i) = c0 * (1.0_dp + eps * ((2.0_dp * (z(i) - 1300.0_dp) / 1300.0_dp) - 1.0_dp + exp(-(2.0_dp * (z(i) - 1300.0_dp) / 1300.0_dp))))
    end do

    ! *** Gaussian starter ***
    k0 = omega / c0
    zs = 1000.0_dp
    fac = 10.0_dp
    nr = 2000
    rmax = 100000.0_dp
    deltar = rmax / nr
    allocate(r(nr))
    do i = 1, nr
        r(i) = (i - 1) * deltar
    end do

    allocate(psi(nz, nr))
    psi(:, 1) = sqrt(k0 / fac) * exp(-(k0 / fac)**2 * (z - zs)**2)

    ! *** Form marching matrix ***
    allocate(A(nz, nz), B(nz, nz), C(nz, nz), L(nz, nz), U(nz, nz), y1(nz), y(nz))
    A = 0.0_dp
    B = 0.0_dp
    C = 0.0_dp
    L = 0.0_dp
    U = 0.0_dp
    y1 = 0.0_dp
    y = 0.0_dp

    do i = 1, nz
        do j = 1, nz
            if (i == j) then
                A(i, j) = -2.0_dp / h2 + k0**2 * ((c0 / c(i))**2 - 1.0_dp)
            else if (abs(i - j) == 1) then
                A(i, j) = 1.0_dp / h2
            end if
        end do
    end do

    B = 2.0_dp * ii * k0 / deltar * eye(nz) - A / 2.0_dp
    C = 2.0_dp * ii * k0 / deltar * eye(nz) + A / 2.0_dp

    ! *** March out in range ***
    call lu_decomp(C, L, U)

    do ir = 1, nr - 1
        y1 = matmul(B, psi(:, ir))
        y = solve_lower(L, y1)
        psi(:, ir + 1) = solve_upper(U, y)
    end do

    ! *** Write data to file ***
    open(unit=10, file='output_r.txt', status='unknown')
    do i = 1, nr
        write(10, *) r(i)
    end do
    close(10)

    open(unit=11, file='output_z.txt', status='unknown')
    do i = 1, nz
        write(11, *) z(i)
    end do
    close(11)

    open(unit=12, file='output_psi.txt', status='unknown')
    do i = 1, nz
        do j = 1, nr
            write(12, *) real(psi(i, j)), aimag(psi(i, j))
        end do
    end do
    close(12)

    deallocate(z, r, c, psi, A, B, C, L, U, y1, y)

contains

    function eye(n) result(I)
        integer, intent(in) :: n
        complex(dp) :: I(n, n)
        integer :: i
        I = 0.0_dp
        do i = 1, n
            I(i, i) = 1.0_dp
        end do
    end function eye

    subroutine lu_decomp(A, L, U)
        complex(dp), intent(in) :: A(:,:)
        complex(dp), intent(out) :: L(:,:), U(:,:)
        integer :: n, i, j, k
        n = size(A, 1)
        L = 0.0_dp
        U = 0.0_dp
        do k = 1, n
            U(k, k) = A(k, k)
            do i = k + 1, n
                L(i, k) = A(i, k) / U(k, k)
                U(k, i) = A(k, i)
            end do
            do i = k + 1, n
                do j = k + 1, n
                    A(i, j) = A(i, j) - L(i, k) * U(k, j)
                end do
            end do
        end do
    end subroutine lu_decomp

    function solve_lower(L, b) result(x)
        complex(dp), intent(in) :: L(:,:), b(:)
        complex(dp) :: x(size(b))
        integer :: n, i, j
        n = size(b)
        x = 0.0_dp
        do i = 1, n
            x(i) = b(i)
            do j = 1, i - 1
                x(i) = x(i) - L(i, j) * x(j)
            end do
            x(i) = x(i) / L(i, i)
        end do
    end function solve_lower

    function solve_upper(U, b) result(x)
        complex(dp), intent(in) :: U(:,:), b(:)
        complex(dp) :: x(size(b))
        integer :: n, i, j
        n = size(b)
        x = 0.0_dp
        do i = n, 1, -1
            x(i) = b(i)
            do j = i + 1, n
                x(i) = x(i) - U(i, j) * x(j)
            end do
            x(i) = x(i) / U(i, i)
        end do
    end function solve_upper

end program pe_simulation
