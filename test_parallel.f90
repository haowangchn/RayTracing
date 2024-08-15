! 需要配置openmp并行环境或者更MPI环境
! 命令行并行命令：
!   编译：gfortran -fopenmp -o test_parallel.f90
!   运行：./test_parallel
! 

program ray_tracing
    use omp_lib
    implicit none

    ! 定义常量
    real, parameter :: depth = 5000.0         ! 水深
    real, parameter :: zu = 0.0               ! 水面位置
    real, parameter :: c0 = 1500.0            ! 初始声速，1500m/s
    real, parameter :: eps = 0.00737          ! eps 是一个常数，表示声速的相对变化率
    real, parameter :: pi = 3.14159           ! pi值
    real, parameter :: theta0 = -10 * pi / 180.0 ! 声线出射角
    integer, parameter :: n = 1000000           ! 迭代次数
    real, parameter :: h = 0.1                ! 迭代步长
    real, parameter :: f   = 50.0
    real :: omega = 2 * pi * f                ! 频率 Hz
    integer :: reflect = 1                    ! 声线垂向切向量符号辅助改变---声线轨迹反转
    real :: temph                             ! 辅助迭代步长---声线轨迹反转

    ! 定义数组
    real :: zp(n)                             ! 深度变量 zp 计算声速剖面 profile---声速剖面计算
    real :: cp(n)                             ! 声速变量 cp，仅是深度方向依赖的函数---声速剖面计算
    real :: c(n), r(n), ksi(n), z(n), zeta(n) ! 声线轨迹数据---轨迹计算
    real :: x(n)                              ! 辅助声速计算变量---声速剖面
    real :: s(n)                              ! 轨迹---幅值计算
    real :: A0(n)                             ! 幅值---幅值计算
    complex :: p(n)                           ! 声压---声压场计算
    complex :: p_auxi(n)                        ! 声压指数项---声压场计算

    ! 定义其他变量
    integer :: i

    ! 声明要写入数据的文件名，设置文件名
    character(len=100) :: RayTrace
    character(len=100) :: SoundProfile
    ! 设置文件名
    RayTrace = 'test_RayTrace.txt'
    SoundProfile = 'SoundProfile.txt'

    ! 定义初始条件---声线轨迹数据
    c(1) = 1550.0              ! 声速数值中的初值
    r(1) = 0.0                 ! 初始水平位置
    z(1) = 10.0                ! 源点深度位置，75m
    ksi(1) = cos(theta0) / c0  ! 初始的切向量 (ksi0,  )
    zeta(1) = sin(theta0) / c0 ! 初始的切向量 ( , zeta0)
    s(1) = 0                   ! 初始轨迹---幅值计算

    ! 并行化主循环
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n - 1
        z(i+1) = z(i) + h * c(i) * zeta(i)
        r(i+1) = r(i) + h * c(i) * ksi(i)

        ! 计算声速剖面
        x(i+1) = 2 * ( z(i+1) - 1300 ) / 1300
        c(i+1) = c0 * ( 1 + eps * ( x(i+1) - 1 + exp( -x(i+1) ) ) )

        ! 更新水平方向切向量
        ksi(i+1) = ksi(i)  - h * 1.0 / c(i)**2 * 0.1
        ! 更新深度方向切向量
        zeta(i+1) = reflect * ( zeta(i) + h * 1.2 / c(i)**2 * (1.0 + 2.4 / c0 * z(i)) ** (-1.5) )

        ! 处理反射
        if (z(i+1) >= depth) then
            zeta(i+1) = -zeta(i+1)
            z(i+1) = 2*depth - z(i+1)
        end if

        if (z(i+1) <= zu) then
            zeta(i+1) = -zeta(i+1)
            z(i+1) = 2*zu - z(i+1)
        end if

        ! 输出检测
        ! write (*,*) z(i), r(i), c(i), ksi(i), zeta(i)
    end do
    !$OMP END PARALLEL DO

    ! 并行化幅值计算
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n-1
        s(i+1) = s(i) + ((z(i+1) - z(i))**2 + (r(i+1) - r(i))**2)**(0.5)
        A0(i+1) = 1 / (4 * pi) * (c(i+1) / c0) / ((s(i+1)**2))
    end do 
    !$OMP END PARALLEL DO

    ! 并行化声压场计算
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n-1
        p_auxi(i) = cmplx(0.0, omega * log(c(i) / c(1)))
        p(i) = A0(i) * exp( p_auxi(i) )
    end do
    !$OMP END PARALLEL DO

    ! 写入文件
    write (*,*) 'Writing to ', RayTrace
    open(unit=1, file=RayTrace, status='replace', action='write')
    write(1, '(A)') '     Depth     Horizontal Length Sound Speed       The tangential vector         A0 &
        &            p(i)'

    do i = 1, n-1
        ! 写入数据到文件
        write(1, '(8(F14.6, 1X))') z(i), r(i), c(i), ksi(i), zeta(i), A0(i), p(i)
    end do 
    ! 关闭文件
    close(1)
    write (*,*) 'Results written to ', RayTrace

end program ray_tracing
