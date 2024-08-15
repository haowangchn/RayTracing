! 修改1 采用MUNK声速剖面
! 声速剖面是对的
! % *** set up Munk SSP ***
! z   = linspace( 0, d, nz );    % grid coordinates
! x = 2 * ( z - 1300 ) / 1300;  
! c = c0 * ( 1 + eps * ( x - 1 + exp( -x ) ) );

program ray_tracing
    use omp_lib
    implicit none

    ! 定义常量
    real, parameter :: depth = 5000.0          ! 水深
    real, parameter :: zu = 0.0                ! 水面位置
    real, parameter :: c0 = 1492.0             ! 初始声速，m/s
    real, parameter :: eps = 0.0074            ! eps 是一个常数，表示声速的相对变化率
    real, parameter :: z0 = 1300               ! 声轴位置
    real, parameter :: pi = 3.14159            ! pi值
    real, parameter :: theta0 = 60 * pi / 180.0 ! 声线出射角
    real, parameter :: theta = 30 * pi / 180.0 !接收处的角度
    integer, parameter :: n = 600000           ! 迭代次数
    real :: h = 0.1                              ! 迭代步长
    real, parameter :: f   = 50.0
    real :: omega = 2 * pi * f                 ! 频率 Hz
    integer :: reflect = 1                     ! 声线垂向切向量符号辅助改变---声线轨迹反转
    real :: temph                              ! 辅助迭代步长---声线轨迹反转
    integer :: m                               ! 声线数量
    real :: dtheta0 = 1 * pi / 180.0           ! dtheta, 变化的角度


    ! 定义数组
    real :: zp(n)                              ! 深度变量 zp 计算声速剖面 profile---声速剖面计算
    real :: cp(n)                              ! 声速变量 cp，仅是深度方向依赖的函数---声速剖面计算
    real :: c(n), r(n), ksi(n), z(n), zeta(n)  ! 声线轨迹数据---轨迹计算
    real :: Eatc(n)                            ! 声速的修正因子---声速剖面
    real :: s(n)                               ! 轨迹---幅值计算
    real :: A0(n)                              ! 幅值---幅值计算
    real :: J(n)                               ! 幅值---fuzhi
    complex :: p(n)                            ! 声压---声压场计算
    complex :: p_auxi(n)                       ! 声压指数项---声压场计算

    ! 定义其他变量
    integer :: i

    ! 声明要写入数据的文件名，设置文件名
    character(len=100) :: RayTrace
    character(len=100) :: SoundProfile
    ! 设置文件名
    RayTrace = 'test_RayTrace.txt'
    SoundProfile = 'SoundProfile.txt'


    ! 定义初始条件---声线轨迹数据
    c(1) = 1492.0                              ! 声速数值中的初值
    r(1) = 0.0                                 ! 初始水平位置
    z(1) = 1300                                ! 源点深度位置，75m
    ksi(1) = cos(theta0) / c0                  ! 初始的切向量 (ksi0,  )
    zeta(1) = sin(theta0) / c0                 ! 初始的切向量 ( , zeta0)
    s(1) = 0                                   ! 初始轨迹---幅值计算
    A0(1) = 0.25 * pi


    ! print *, 'Depth :', 'Horizontal Length :', 'Sound Speed: ','The tangential vector :'
    !! TODO 从源点以一定角度发射多条声线
    do i = 1, n - 1

        z(i+1) = z(i) + h * c(i) * zeta(i)
        r(i+1) = r(i) + h * c(i) * ksi(i)

        ! 更新水平方向切向量
        ksi(i+1) = ksi(i)
        ! 更新深度方向切向量
        zeta(i+1) = reflect * ( zeta(i) + h / (c(i)**2) * c0 *eps / 650 * (1 - exp((1300 - z(i)) / 650)))
        ! zeta(i+1) = reflect * ( zeta(i) + h * 1.2 / c(i)**2 * (1.0 + 2.4 / c0 * z(i)) ** (-1.5) )

        ! 计算声速剖面
        Eatc(i+1) = 2 * ( z(i+1) - z0 ) / z0
        c(i+1) = c0 * ( 1 + eps * ( Eatc(i+1) + exp( -Eatc(i+1)) - 1  ) )

        !! XXX 考虑反射，应当引入一个变量控制反射方向
        ! 海底反射，且修正步长
        if (z(i+1) >= depth) then
            ! 计算到达海底的精确步长
            temph = (depth - z(i)) / (c(i) * zeta(i))
            ! 更新z和r到达海底的精确位置
            z(i+1) = z(i) + temph * c(i) * zeta(i)
            r(i+1) = r(i) + temph * c(i) * ksi(i)
            ! 反转zeta方向
            zeta(i+1) = -zeta(i)

            ! 使用小步长进行平滑调整
            temph = h - temph
            do while (temph > 0)
                ! 计算下一步的步长，取较小的步长进行调整
                if (temph > h) then
                    temph = h
                end if
                z(i+1) = z(i+1) + temph * c(i) * zeta(i+1)
                r(i+1) = r(i+1) + temph * c(i) * ksi(i+1)
                temph = temph - h
            end do
        endif

        ! 海面反射，且修正步长
        if (z(i+1) <= zu) then
            ! 计算到达海面的精确步长
            temph = (zu - z(i)) / (c(i) * zeta(i))
            ! 更新z和r到达海面的精确位置
            z(i+1) = z(i) + temph * c(i) * zeta(i)
            r(i+1) = r(i) + temph * c(i) * ksi(i)
            ! 反转zeta方向
            zeta(i+1) = -zeta(i)

            ! 使用小步长进行平滑调整
            temph = h - temph
            do while (temph > 0)
                ! 计算下一步的步长，取较小的步长进行调整
                if (temph > h) then
                    temph = h
                end if
                z(i+1) = z(i+1) + temph * c(i) * zeta(i+1)
                r(i+1) = r(i+1) + temph * c(i) * ksi(i+1)
                temph = temph - h
            end do
        end if


        ! 输出检测
        ! write (*,*) z(i), r(i), c(i), ksi(i), zeta(i)
    end do



    ! -----------------------计算幅值------------------------------
    !! TODO: 计算幅值
    do i = 1, n-1
        J(i) = r(i) / sin(theta) * (r(i+1) - r(i)) / dtheta0
        A0(i+1) = A0(1) * (c(i) * J(1) / c0 / J(i))
        ! s(i+1) = s(i) + ((z(i+1) - z(i))**2 + (r(i+1) - r(i))**2)**(0.5)
        ! ! write(*,*) s(i)
        ! A0(i+1) = A0(1) * (c(i+1) / c0) / ((s(i+1)**2))
        ! ! write(*,*) A0(i)
    end do 
    write (*,*) "J(1)", J(1)
    write (*,*) A0(1)
    write (*,*) "J(2)", J(2)
    write (*,*) A0(2)
    write (*,*) A0(3)
    write (*,*) A0(4)
    write (*,*) A0(5)
    ! -----------------------幅值计算结束---------------------------




    ! -----------------------计算声压场-----------------------------
    !! TODO计算声压，写成一个子程序
    do i = 1, n-1
        p_auxi(i) = cmplx(0.0, omega * log(c(i) / c(1)))
        p(i) = A0(i) * exp( p_auxi(i) )
        ! write (*,*) "Wrinting pressure : ", p(i)

    end do
    ! ---------------------------声压场结束--------------------------




    ! ----------------------------写入文件---------------------------
    write (*,*) 'Writing to ......', RayTrace
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
    ! ----------------------------写入文件结束------------------------

    ! 调用子程序计算声速剖面并写入文件
    ! call calculate_and_write_sound_profile(depth, zp, cp, c0, SoundProfile)

end program ray_tracing





    ! ! 如果透过海底，修正步长，然后直接反射
    ! if (z(n+1) >= depth) then
    !     h = (depth - z(n)) / c(n) / zeta(n)
    !     zeta(n+1) = zeta(n) - 1 / c(n)**2 * (-1.2 * (1 + 2.4 / c0 * z(n)))
    ! end if


    ! TODO 考虑海面如何反射
    ! ! 如果到达海面，修正步长直接反射
    ! if (z(n+1) <= zu) then
    !     zeta(n+1) = zeta(n) - 1 / c(n)**2 * (-1.2 * (1 + 2.4 / c0 * z(n)))
    ! end if