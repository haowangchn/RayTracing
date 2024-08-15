! 水深65米，声速1500m/s，密度1025kg/m^3
! 海底，无界，声速1620m/s，密度2600kg/m^3
! 空气，声速无限小、密度无限小

! 源点：水深30m，
! 场点1：水平距离100米，水深30米
! 场点2：水平距离1000米，水深30米

program ray_tracing
    use omp_lib
    implicit none
    real :: start_time, end_time

!---------------------------------------变量声明-------------------------------------------
    integer, parameter :: n_p = 40                        ! 迭代次数--声速剖面
    real :: h_p = 5                                          ! 迭代步长--声速剖面
    real, allocatable :: z_p(:), c_p(:), Eatc_p(:), zeta_p(:) ! 声速剖面计算--声速剖面
    real, parameter :: z0 = 1300                              ! 声轴位置--声速剖面
    real, parameter :: c0 = 1500                            ! 初始声速，m/s--声速剖面@@@@@@@@@@@@@@@

    real, parameter :: depth = 5000.0                         ! 水深
    real, parameter :: zu = 0.0                               ! 水面位置
    real, parameter :: eps = 0.0074                           ! eps 是一个常数，表示声速的相对变化率
    real, parameter :: pi = 3.14159                           ! pi值
    real, parameter :: theta0 = 0 * pi / 180.0               ! 声线出射角
    real, parameter :: dtheta0 = 0.1 * pi / 180.0             ! dtheta--声束计算
    real, parameter :: dphi = 1 * pi / 180.0                  ! dphi--声线束步长---特征声线计算
    real, parameter :: theta = 30 * pi / 180.0                ! 场点接收处的角度  !! FIXME

    real, parameter :: rsp = 0.0                              ! 源点水平坐标---源点 r or z values of source point
    real, parameter :: zsp = 30.0                           ! 源点深度坐标---源点
    ! real, parameter :: rfp = 0.19                             ! 场点水平坐标---场点 r or z values of field point
    ! real, parameter :: zfp = 1000.112                         ! 场点深度坐标---场点
    real, parameter :: rfp = 100.0
    real, parameter :: zfp = 30.0

    integer, parameter :: n = 50                              ! 迭代次数--声线轨迹及相应声速剖面
    integer, parameter :: m = 80                              ! 声线数量--声线组!! NOTE @@@@@@
    real, parameter :: ds = 5                                  ! 迭代步长--声线轨迹及相应声速剖面
    integer, parameter :: reflect = 1                         ! 声线垂向切向量符号辅助改变---声线轨迹反转
    real :: temph                                             ! 辅助迭代步长---声线轨迹反转 !! NOTE @@@@@@
    real :: phi_a, phi_b                                      ! 跨越介质的时候折射角变化
    real :: csp = c0                                          ! 发射点初始声速!! NOTE @@@@@@

    real :: delta_r = 0                                       ! ---辅助计算特征声线
    real :: delta_z = 0                                       ! ---辅助计算特征声线


    integer :: eigenray_num = 0                               ! 记录特征声线编号---计算场点声压
    integer :: eigenray_itern = 0                             ! 记录特征声线迭代次数---计算场点声压
    integer :: i, k, ip, jj                                   ! 定义迭代变量                        
    character(len=100) :: RayTrace
    character(len=100) :: SoundProfile

    ! 数组
    real, allocatable :: c(:,:), r(:,:), ksi(:,:), z(:,:), zeta(:,:)
    real, allocatable :: Eatc(:,:), s(:,:)                    ! 轨迹和幅值数组
    ! 单条声束
    real, allocatable :: phi(:)
    real, allocatable :: r_a(:,:), z_a(:,:), ksi_a(:,:), zeta_a(:,:), Eatc_a(:,:), c_a(:,:)    ! 单条声束下界
    real, allocatable :: r_b(:,:), z_b(:,:), ksi_b(:,:), zeta_b(:,:), Eatc_b(:,:), c_b(:,:)    ! 单条声束上界
    real, allocatable :: r_c(:,:), z_c(:,:), c_c(:,:)                                          ! 单条声束中线
    
    ! 声明------ 计算声压和幅值 ------
    ! real :: f   = 250                                       ! 频率 Hz
    integer :: f = 0
    integer :: f_min = 100
    integer :: f_max = 2000
    integer :: f_step = 100
    real :: s_i, s_j
    real :: omega
    real :: AOJ = 10000                                       ! 随便赋的初值
    real, allocatable :: A0(:,:), J(:,:)                      ! 幅值和J数组
    complex, allocatable :: p(:,:), p_auxi(:,:)               ! 声压和辅助声压数组---计算场点声压                                         
    complex :: p_fp      ! 场点声压
    real :: p_ref       ! 参考声压（Pa）
    real :: spl_real, spl_imag         ! 声压级（dB）
    real :: p_fp_mag    ! 测得声压的幅值（实数）
    real :: p_fp_real, p_fp_imag


    ! 声明------ 捕捉场点的矩形 ------
    real :: r1, r2, z1, z2                                    ! 矩形四个角点
    real :: r_min, r_max, z_min, z_max                        ! 矩形的边界
    real :: min_distance = 1.0e+30                            ! 距离场点最近的矩形的距离，初始化最小距离为一个较大的值
    real :: current_distance = 0.0

    real :: c_sea = 1500                 ! 海水中的声速
    real :: c_haid = 1620               ! 海底中的声速
    real :: c_air = 15                   ! 空气中的声速


!---------------------------------------变量初始化-------------------------------------------
    allocate( c_p(n_p), z_p(n_p), Eatc_p(n_p),zeta_p(n_p) )                                ! 声速剖面计算

    allocate( c(m, n), r(m, n), ksi(m, n), z(m, n), zeta(m, n) )
    allocate( Eatc(m, n), s(m, n) )

    allocate( phi(m))
    allocate( r_a(m, n), z_a(m, n), ksi_a(m, n), zeta_a(m, n), Eatc_a(m, n), c_a(m, n) )   ! 单条声束下界
    allocate( r_b(m, n), z_b(m, n), ksi_b(m, n), zeta_b(m, n), Eatc_b(m, n), c_b(m, n) )   ! 单条声束上界
    allocate( r_c(m, n), z_c(m, n), c_c(m, n) )                                            ! 单条声束中线
 
    ! 初始化------ 计算声压和幅值 ------
    allocate(A0(m, n), J(m, n))

    ! J(k, 1) = 
    allocate(p(m, n), p_auxi(m, n))

    ! 设置文件名
    RayTrace = 'fp_RayTrace.txt'
    SoundProfile = 'fp_SoundProfile.txt'

    ! 数组赋初值
    c_p(1) = c_sea                                                ! 声速数值中的初值---声速剖面计算
    z_p(1) = 0                                                 ! z方向初值---声速剖面计算
    zeta_p(1) = sin(theta0) / c0                               ! 初始纵向切向量 ( , zeta0)

    ! 捕捉场点的矩形初始化
    r1 = 0.0                                                   ! 矩形四个角点
    r2 = 0.0 
    z1 = 0.0 
    z2 = 0.0                                         
    r_min = 0.0                                                ! 矩形的边界
    r_max = 0.0 
    z_min = 0.0 
    z_max = 0.0                                     

    ! 计算声压
    p_ref = 20.0e-6

    ! 获取开始时间
    call cpu_time(start_time)

!------------------------------------------------程序计算部分-------------------------------------------------
    !--------------------------声速剖面计算--------------------------验证正确
    do ip = 1, n_p-1
        if (z_p(ip) < 0) then
            c_p(ip) = 15
        else if (z_p(ip) >= 0 .and. z_p(ip) <= 65) then
            c_p(ip) = 1500
        else if (z_p(ip) >= 65) then
            c_p(ip) = 1620
        end if

        z_p(ip+1) = z_p(ip) + h_p
    end do
    !--------------------------声速剖面计算结束--------------------------

    !--------------------------声线轨迹计算--------------------------
    ! 计算m条声线束的轨迹
    do k = 1, m ! 声线数量--声线组
        ! 定义初始条件---声线轨迹数据
        phi(k) = theta0 + dphi * (k - m / 2)                ! 以出射角为中心，角度步长为1°，发出一个扇面的声束，声束数量是m 
        c_a(k, 1) = c_sea 
        r_a(k, 1) = rsp                                     ! 初始水平位置---源点
        z_a(k, 1) = zsp                                     ! 初始深度位置---源点
        ksi_a(k, 1) = cos(phi(k)) / c_sea                   ! 初始的水平切向量 (ksi0,  )---声束下界
        zeta_a(k, 1) = sin(phi(k)) / c_sea                  ! 初始的垂直切向量

        c_b(k, 1) = c_sea 
        r_b(k, 1) = rsp                                     ! 初始水平位置---源点
        z_b(k, 1) = zsp                                     ! 初始深度位置---源点
        ksi_b(k, 1) = cos(phi(k) + dtheta0) / c_sea         ! 初始的水平切向量 (ksi0,  )---声束上界
        zeta_b(k, 1) = sin(phi(k) + dtheta0) / c_sea        ! 初始的垂直切向量

        ! 计算一组声线的声束，_a，_b分别为每条声束的下界和上界。
        do i = 1, n-1  ! n为迭代次数。此处为n-1，否则为n时，会出现数组越界的异常值
            !------声束下边界------
            if (z_a(k, i) < 0) then
                ! 考虑反射，更新切向量
                ksi_a(k, i+1) = ksi_a(k, 1)
                zeta_a(k, i+1) = -1 * zeta_a(k, 1)             

                c_a(k, i) = 15.0

                r_a(k, i+1) = r_a(k, i) + ds * ksi_a(k, i) * c_a(k, i)    ! 更新水平距离
                z_a(k, i+1) = z_a(k, i) + ds * zeta_a(k, i) * c_a(k, i)   ! 更新深度

            else if (z_a(k, i) >= 0 .and. z_a(k, i) <= 65) then
                c_a(k, i) = 1500

                ksi_a(k, i+1) = ksi_a(k, i)
                zeta_a(k, i+1) = zeta_a(k, i)

                r_a(k, i+1) = r_a(k, i) + ds * ksi_a(k, i) * c_a(k, i)    ! 更新水平距离
                z_a(k, i+1) = z_a(k, i) + ds * zeta_a(k, i) * c_a(k, i)   ! 更新深度

            else if (z_a(k, i) > 65) then
                ! 考虑反射，更新切向量

                ksi_a(k, i+1) = ksi_a(k, 1)
                zeta_a(k, i+1) = -1 * zeta_a(k, 1)

                c_a(k, i) = 1620.0

                r_a(k, i+1) = r_a(k, i) + ds * ksi_a(k, i) * c_a(k, i)    ! 更新水平距离
                z_a(k, i+1) = z_a(k, i) + ds * zeta_a(k, i) * c_a(k, i)   ! 更新深度

            end if

            !------声束上边界------
            if (z_b(k, i) < 0) then
                ! 计算折射角并更新切向量 (空气中的声速极低，假设全反射)
                ksi_b(k, i+1) = ksi_b(k, 1)
                zeta_b(k, i+1) = -1 * zeta_b(k, 1)                      !

                r_b(k, i+1) = r_b(k, i) + ds * ksi_b(k, i+1) * c_b(k, i)    ! 更新水平距离
                z_b(k, i+1) = z_b(k, i) + ds * zeta_b(k, i+1) * c_b(k, i)   ! 更新深度

            else if (z_b(k, i) >= 0 .and. z_b(k, i) <= 65) then
                c_b(k, i) = 1500

                ksi_b(k, i+1) = ksi_b(k, 1)
                zeta_b(k, i+1) = zeta_b(k, 1)

                r_b(k, i+1) = r_b(k, i) + ds * ksi_b(k, i+1) * c_b(k, i)    ! 更新水平距离
                z_b(k, i+1) = z_b(k, i) + ds * zeta_b(k, i+1) * c_b(k, i)   ! 更新深度

            else if (z_b(k, i) > 65) then
                ! 更新切向量
                ksi_b(k, i+1) = ksi_b(k, 1)
                zeta_b(k, i+1) = -1 * zeta_b(k, 1)

                c_b(k, i) = 1620.0

                r_b(k, i+1) = r_b(k, i) + ds * ksi_b(k, i+1) * c_b(k, i)    ! 更新水平距离
                z_b(k, i+1) = z_b(k, i) + ds * zeta_b(k, i+1) * c_b(k, i)   ! 更新深度

            end if

            ! !------声束中线------
            ! r_c(k, i) = 0.5 * (r_a(k, i) + r_b(k, i))
            ! z_c(k, i) = 0.5 * (z_a(k, i) + z_b(k, i))
            ! c_c(k, i) = 0.5 * (c_a(k, i) + c_b(k, i))
        end do
    end do



    !--------------------------声线轨迹计算结束--------------------------   
    write (*,*) "-------------------------------find the eigenray----------------------------------"
    ! 第一种好像不太靠谱的捕捉场点的方式：
    ! 根据场点数据查找场点被哪两条声线包围，类似于捞网---求特征声线
    ! 初始化最小距离为一个较大的值
    do i = 1, n-1
        do k = 1, m-1
            ! 获取当前矩形的四个角点
            r1 = r_a(k, i)
            z1 = z_a(k, i)
            r2 = r_a(k + 1, i + 1)
            z2 = z_a(k + 1, i + 1)

            ! 计算矩形的边界
            r_min = min(r1, r2)
            r_max = max(r1, r2)
            z_min = min(z1, z2)
            z_max = max(z1, z2)

            ! 当矩形的边界紧密相邻时，场点可能会恰好位于多个矩形的交界处，可能会出现返回多个值的现象
            ! 判断场点是否在矩形内部
            if (rfp >= r_min .AND. rfp <= r_max .AND. zfp >= z_min .AND. zfp <= z_max) then
                ! 计算场点到矩形中心的距离
                current_distance = sqrt((rfp - (r_min + r_max) / 2)**2 + (zfp - (z_min + z_max) / 2)**2)
                
                ! 如果当前距离比最小距离更小，则更新
                if (current_distance < min_distance) then
                    min_distance = current_distance
                    eigenray_num = k
                    eigenray_itern = i
                end if
            end if
        end do
    end do

    ! ! 根据特征声线数据，对场点数据进行插值。输出场点声压。
    ! !! TODO 测试场点位于海面、海底、源点很近距离时的极端情况。

    ! -----------------------计算幅值、声压------------------------------

    if (eigenray_num > 0 .AND. eigenray_num <= 100 .AND. eigenray_itern > 0 .AND. eigenray_itern <= 100) then
        k = eigenray_num
        i = eigenray_itern

        s_i = ds * i
        J(k, i) = - s_i**2 * cos(theta0)  ! 式子(3.54)
        A0(k, i) = 0.25 / pi * sqrt(abs(c_a(k, i) * cos(theta0) / c0 / J(k, i)))    ! 式子(3.56)

        s_j = ds * (i + 1)
        J(k+1, i) = - s_j**2 * cos(theta0)  ! 式子(3.54)
        A0(k+1, i) = 0.25 / pi * sqrt(abs(c_a(k+1, i) * cos(theta0) / c0 / J(k+1, i))) 
        
        write (*,*) A0(k,i)  ! 式子(3.56)
        write (*,*) A0(k+1,i)  ! 式子(3.56)  

        ! 打开文件进行写入
        open(unit=3, file='frequency_pressure.txt', status='replace', action='write')
        ! 写入标题行
        write(3, '(A)') 'Frequency Real_Pressure_Level Imaginary_Pressure_Level'

        ! 计算频率-压强曲线
        do f = f_min, f_max, f_step
            omega = 2 * pi * f
            p_auxi(k, i) = cmplx(0.0, omega * 1 / (c_a(k, i)) * s_i)               ! 书中式子（3.57）
            p(k, i) = A0(k, i) * exp(p_auxi(k, i))                                 ! 书中式子（3.57）

            p_auxi(k+1, i) = cmplx(0.0, omega * 1 / (c_a(k+1, i)) * s_j)           ! 书中式子（3.57）
            p(k+1, i) = A0(k+1, i) * exp(p_auxi(k+1, i))                           ! 书中式子（3.57'])

            ! 计算声压幅值
            p_fp = 0.5 * (p(k, i) + p(k+1, i))

            ! 提取实部和虚部
            p_fp_real = real(p_fp)
            p_fp_imag = aimag(p_fp)

            ! 将频率和声压级写入文件
            write(3, '(I10, SP, F10.5, SP, F10.5)') f, p_fp_real, p_fp_imag

        end do 

        ! 关闭文件
        close(3)
    else
        write(*,*) "Error: eigenray_num or eigenray_itern is out of range"
    end if
    ! -----------------------幅值、声压计算结束---------------------------


!------------------------------------------------程序计算部分结束-------------------------------------------------

    ! ! ----------------------------写入声速剖面---------------------------
    ! ! write (*,*) 'Writing to ......', SoundProfile
    ! open(unit=1, file=SoundProfile, status='replace', action='write')
    ! write(1, '(A)') '     Depth     SoundProfile'

    ! do ip = 1, n_p
    !     write(1, '(F14.1, 1X, F14.2)') z_p(ip), c_p(ip)
    !     ! write(1, '(F14.1, 1X, F14.2, 1X, A)') z_p(ip), c_p(ip), '/'
    ! end do 

    ! close(1)
    ! write (*,*) 'Results written to ', SoundProfile
    ! ! ----------------------------写入声速剖面------------------------


    ! ----------------------------写入声线轨迹---------------------------
    write (*,*) 'Writing to ......', RayTrace
    write (*,*) 'rsp:', rsp
    write (*,*) 'zsp:', zsp
    write (*,*) 'rfp:', rfp
    write (*,*) 'zfp:', zfp
    ! 打开文件进行写入
    open(unit=2, file=RayTrace, status='replace', action='write')

    ! 写入源点信息
    write(2, '(A)') 'The source point is : '
    write(2, '(2F14.6)') rsp, zsp

    ! 写入场点信息
    write(2, '(A)') 'The field point is : '
    write(2, '(2F14.6)') rfp, zfp

    write(2, '(A)') "The closest ray number k and k+1: "
    write(2, '(2I14)') eigenray_num, eigenray_num + 1


    ! 写入声线数据
    do k = 1, m ! m为声线束数量
        ! 在每条声线束的数据前添加一个标题行
        write(2, '(A, I0)') 'Ray Bundle ', k
        do i = 1, n  ! n为迭代次数
            write(2, '(2F14.6)') z_a(k, i), r_a(k, i)
        end do
    end do 
    ! 关闭文件
    close(2)
    write (*,*) 'Results written to ', RayTrace
    ! ----------------------------声线轨迹写入结束------------------------


    ! 获取结束时间
    call cpu_time(end_time)
    print *, 'Total running time: ', end_time - start_time, 's'
end program ray_tracing