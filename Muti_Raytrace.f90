program ray_calculation
    implicit none
    ! 定义常量
    real, parameter :: depth = 5000.0           ! 水深
    real, parameter :: zu = 0.0                 ! 水面位置
    real, parameter :: c0 = 1492.0              ! 初始声速，1492m/s
    real, parameter :: eps = 0.0074             ! eps 是一个常数，表示声速的相对变化率
    real, parameter :: z0 = 1300                ! 声轴位置
    real, parameter :: pi = 3.14159             ! pi值
    real, parameter :: theta0 = -10 * pi / 180.0 ! 声线出射角
    integer, parameter :: iterationNum = 600000  ! 迭代次数
    real :: h = 1                               ! 迭代步长
    integer, parameter :: rayID = 2             ! 声线数量

    ! 定义数组
    real :: z(rayID, iterationNum)              ! 初值距离
    real :: r(rayID, iterationNum)              ! 水平距离
    real :: c(rayID, iterationNum)              ! 给定点的声速
    real :: ksi(rayID, iterationNum)            ! 水平方向切向量
    real :: zeta(rayID, iterationNum)           ! 深度方向切向量
    real :: Eatc(rayID, iterationNum)           ! 声速的修正因子

    ! 定义其他变量
    integer :: j

    ! 设置初始条件
    do j = 1, rayID
        c(j, 1) = c0
        r(j, 1) = 0.0
        z(j, 1) = 10.0
        ksi(j, 1) = cos(theta0) / c0
        zeta(j, 1) = sin(theta0) / c0
        Eatc(j, 1) = 0.0
    end do

    ! 从源点以一定角度发射多条声线
    do j = 1, rayID
        call ray_tracing_sub(z, r, c, ksi, zeta, Eatc, h, 1, iterationNum, depth, zu, z0, c0, eps, j)
    end do

contains

    subroutine ray_tracing_sub(z, r, c, ksi, zeta, Eatc, h, reflect, iterationNum, depth, zu, z0, c0, eps, j)
        implicit none
        integer, intent(in) :: iterationNum, reflect, j
        real, intent(inout) :: z(:,:), r(:,:), c(:,:), ksi(:,:), zeta(:,:), Eatc(:,:)
        real, intent(inout) :: h
        real, intent(in) :: depth, zu, z0, c0, eps
        integer :: i

        ! 调试输出：初始条件
        write (*,*) 'Calculating ray ', j

        do i = 1, iterationNum - 1
            z(j, i+1) = z(j, i) + h * c(j, i) * zeta(j, i)
            r(j, i+1) = r(j, i) + h * c(j, i) * ksi(j, i)

            ! MUNK声速剖面
            Eatc(j, i+1) = 2 * (z(j, i+1) - z0) / z0
            c(j, i+1) = c0 * (1 + eps * (Eatc(j, i+1) + exp(-Eatc(j, i+1)) - 1))

            ! 更新水平方向切向量
            ksi(j, i+1) = ksi(j, i) ! - h * (1.0 / c(j, i)**2) * 0.1

            ! 更新深度方向切向量
            zeta(j, i+1) = reflect * (zeta(j, i) + h * 1.2 / c(j, i)**2 * (1.0 + 2.4 / c0 * z(j, i)) ** (-1.5))

            ! 海底反射，且修正步长
            if (z(j, i+1) >= depth) then
                zeta(j, i+1) = -zeta(j, i+1)
                z(j, i+1) = 2 * depth - z(j, i+1)
                h = (depth - z(j, i)) / (c(j, i) * zeta(j, i))
                ! 调试输出：海底反射
                write (*,*) 'Bottom reflection at iteration ', i, ' for ray ', j
                write (*,*) 'New zeta = ', zeta(j, i+1), ' New z = ', z(j, i+1), ' New h = ', h
            end if

            ! 海面反射，且修正步长
            if (z(j, i+1) <= zu) then
                zeta(j, i+1) = -zeta(j, i+1)
                z(j, i+1) = 2 * zu - z(j, i+1)
                h = (zu - z(j, i)) / (c(j, i) * zeta(j, i))
                ! 调试输出：海面反射
                write (*,*) 'Surface reflection at iteration ', i, ' for ray ', j
                write (*,*) 'New zeta = ', zeta(j, i+1), ' New z = ', z(j, i+1), ' New h = ', h
            end if

            ! 调试输出：每1000次迭代输出一次
            if (mod(i, 1000) == 0) then
                write (*,*) 'Iteration ', i, ' for ray ', j
                write (*,*) 'z = ', z(j, i+1), ' r = ', r(j, i+1), ' c = ', c(j, i+1)
                write (*,*) 'ksi = ', ksi(j, i+1), ' zeta = ', zeta(j, i+1)
            end if
        end do
    end subroutine ray_tracing_sub

end program ray_calculation
