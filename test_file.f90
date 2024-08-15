program ray_tracing
    use omp_lib
    implicit none
    real :: start_time, end_time

!---------------------------------------变量声明-------------------------------------------
    integer, parameter :: n_p = 500000                        ! 迭代次数--声速剖面
    real :: h_p = 10                                          ! 迭代步长--声速剖面
    real, allocatable :: z_p(:), c_p(:), Eatc_p(:), zeta_p(:) ! 声速剖面计算--声速剖面
    real, parameter :: z0 = 1300                              ! 声轴位置--声速剖面
    real, parameter :: c0 = 1492.0                            ! 初始声速，m/s--声速剖面


    real, parameter :: depth = 5000.0                         ! 水深
    real, parameter :: zu = 0.0                               ! 水面位置
    real, parameter :: eps = 0.0074                           ! eps 是一个常数，表示声速的相对变化率
    real, parameter :: pi = 3.14159                           ! pi值
    real, parameter :: theta0 = 60 * pi / 180.0               ! 声线出射角
    real, parameter :: dtheta0 = 0.1 * pi / 180.0             ! dtheta--声束计算
    real, parameter :: dphi = 1 * pi / 180.0                  ! dphi--声线束步长---特征声线计算
    real, parameter :: theta = 30 * pi / 180.0                ! 场点接收处的角度  !! FIXME

    real, parameter :: rsp = 0.0                              ! 源点水平坐标---源点 r or z values of source point
    real, parameter :: zsp = 1000.0                           ! 源点深度坐标---源点
    ! real, parameter :: rfp = 0.19                             ! 场点水平坐标---场点 r or z values of field point
    ! real, parameter :: zfp = 1000.112                         ! 场点深度坐标---场点
    real, parameter :: rfp = 1.0
    real, parameter :: zfp = 1000.0

    integer, parameter :: n = 5000                           ! 迭代次数--声线轨迹及相应声速剖面
    integer, parameter :: m = 25                              ! 声线数量--声线组!! NOTE @@@@@@
    real, parameter :: h = 10                                ! 迭代步长--声线轨迹及相应声速剖面
    real, parameter :: f   = 100.0
    real :: omega = 2 * pi * f                                ! 频率 Hz
    integer, parameter :: reflect = 1                         ! 声线垂向切向量符号辅助改变---声线轨迹反转
    real :: temph                                             ! 辅助迭代步长---声线轨迹反转 !! NOTE @@@@@@
    real :: csp = c0                                          ! 发射点初始声速!! NOTE @@@@@@

    real :: delta_r = 0                                       ! ---辅助计算特征声线
    real :: delta_z = 0                                       ! ---辅助计算特征声线


    integer :: eigenray_num = 0                               ! 记录特征声线编号---计算场点声压
    integer :: eigenray_itern = 0                             ! 记录特征声线迭代次数---计算场点声压
    integer :: i, k, ip, jj                                       ! 定义迭代变量                        
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
    real, allocatable :: A0(:,:), J(:,:)                      ! 幅值和J数组
    complex, allocatable :: p(:,:), p_auxi(:,:)               ! 声压和辅助声压数组---计算场点声压
    complex :: p_fp                                           ! 场点声压

    ! 声明------ 捕捉场点的矩形 ------
    real :: r1, r2, z1, z2                                    ! 矩形四个角点
    real :: r_min, r_max, z_min, z_max                        ! 矩形的边界
    real :: min_distance = 1.0e+30                            ! 距离场点最近的矩形的距离，初始化最小距离为一个较大的值
    real :: current_distance = 0.0

    ! ! 声明最小距离和相应的声线编号
    ! real :: r_min, r_max, z_min, z_max
    ! real :: distance, min_distance1, min_distance2
    ! integer :: closest_ray1, closest_ray2
    ! integer :: iter1, iter2
    ! logical :: found_ray1, found_ray2
    ! real :: r1, z1, r2, z2, z_current, z_next, z_min_current, z_max_current
    ! real, dimension(:), allocatable :: z_distances
    ! integer, dimension(:), allocatable :: z_indices
    ! real :: temp_distance
    ! integer :: temp_index
    ! integer :: l


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
    RayTrace = 'test_RayTrace.txt'
    SoundProfile = 'test_SoundProfile.txt'

    ! 数组赋初值
    c_p(1) = c0                                                ! 声速数值中的初值---声速剖面计算
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

    ! ! 初始化最小距离和相应的声线编号
    ! min_distance1 = 1.0e+30
    ! min_distance2 = 1.0e+30
    ! closest_ray1 = -1
    ! closest_ray2 = -1
    ! found_ray1 = .false.
    ! found_ray2 = .false.  
    ! allocate(z_distances(m-1))
    ! allocate(z_indices(m-1))                                

    ! 获取开始时间
    call cpu_time(start_time)
!------------------------------------------------程序计算部分-------------------------------------------------
    ! 计算声速剖面--验证成功
    do ip = 1, n_p
        ! 更新深度方向切向量
        zeta_p(ip+1) = zeta_p(ip) - h_p / (c_p(ip)**2) * c0 * eps / 650 * (1 - exp((1300 - z_p(ip)) / 650))
        z_p(ip+1) = z_p(ip) + h_p * c_p(ip) * zeta_p(ip)

        Eatc_p(ip+1) = 2 * ( z_p(ip+1) - z0 ) / z0
        c_p(ip+1) = c0 * ( 1 + eps * ( Eatc_p(ip+1) + exp( -Eatc_p(ip+1)) - 1  ) )
    end do 

    ! 计算m条声线束的轨迹
    do k = 1, m ! 声线数量--声线组
        ! 定义初始条件---声线轨迹数据
        phi(k) = theta0 + dphi * (k - m / 2)                ! 以出射角为中心，角度步长为1°，发出一个扇面的声束，声束数量是m            
        
        c_a(k, 1) = csp 
        r_a(k, 1) = rsp                                     ! 初始水平位置---源点
        z_a(k, 1) = zsp                                     ! 初始深度位置---源点
        ksi_a(k, 1) = cos(phi(k)) / c0                      ! 初始的水平切向量 (ksi0,  )---声束下界
        zeta_a(k, 1) = sin(phi(k)) / c0

        c_b(k, 1) = csp 
        r_b(k, 1) = rsp                                     ! 初始水平位置---源点
        z_b(k, 1) = zsp                                     ! 初始深度位置---源点
        ksi_b(k, 1) = cos(phi(k) + dtheta0) / c0            ! 初始的水平切向量 (ksi0,  )---声束上界
        zeta_b(k, 1) = sin(phi(k)+ dtheta0) / c0

        ! 计算一组声线的声束，_a，_b分别为每条声束的下界和上界。
        do i = 1, n  ! n为迭代次数。
            !------声束下边界------
            r_a(k, i+1) = r_a(k, i) + h * c_a(k, i) * ksi_a(k, i)     ! 更新水平距离，对应切向量不变，r只随声速变化。
            z_a(k, i+1) = z_a(k, i) + h * c_a(k, i) * zeta_a(k, i)    ! 更新深度

            ksi_a(k, i+1) = ksi_a(k, i)
            zeta_a(k, i+1) = zeta_a(k, i) - h / (c_a(k, i)**2) * c0 * eps / 650 * (1 - exp((1300 - z_a(k, i)) / 650))  ! 更新深度距离
            
            Eatc_a(k, i+1) = 2 * ( z_a(k, i+1) - z0 ) / z0            ! 辅助变量---更新声速
            c_a(k, i+1) = c0 * ( 1 + eps * ( Eatc_a(k, i+1) + exp( -Eatc_a(k, i+1)) - 1  ) )                           ! 更新声速


            !------声束上边界------
            r_b(k, i+1) = r_b(k, i) + h * c_b(k, i) * ksi_b(k, i)     ! 更新水平距离
            z_b(k, i+1) = z_b(k, i) + h * c_b(k, i) * zeta_b(k, i)    ! 更新深度

            ksi_b(k, i+1) = ksi_b(k, i)
            zeta_b(k, i+1) = zeta_b(k, i) - h / (c_b(k, i)**2) * c0 * eps / 650 * (1 - exp((1300 - z_b(k, i)) / 650))  ! 更新深度距离

            Eatc_b(k, i+1) = 2 * ( z_b(k, i+1) - z0 ) / z0            ! 辅助变量---更新声速
            c_b(k, i+1) = c0 * ( 1 + eps * ( Eatc_b(k, i+1) + exp( -Eatc_b(k, i+1)) - 1  ) )                           ! 更新声速
            
            !------声束中线------
            r_c(k, i) = 0.5 * (r_a(k, i) + r_b(k, i))
            z_c(k, i) = 0.5 * (z_a(k, i) + z_b(k, i))
            c_c(k, i) = 0.5 * (c_a(k, i) + c_b(k, i))
        end do
    end do
    
    
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

    ! 输出最接近的矩形的信息
    write (*,*) "Congratulations! The field point is found!"
    write (*,*) "The closest ray number k and k+1: ", eigenray_num, eigenray_num + 1
    write (*,*) "The closest iter count: ", eigenray_itern, eigenray_itern + 1

    ! ! 第二种捕捉场点的方式：
    ! ! 步骤 1: 检测水平方向上的声线是否到达场点的水平位置
    ! do i = 1, n-1
    !     do k = 1, m-1
    !         ! 获取当前矩形的四个角点
    !         r1 = r_a(k, i)
    !         r2 = r_a(k + 1, i + 1)

    !         ! 计算矩形的水平边界
    !         r_min = min(r1, r2)
    !         r_max = max(r1, r2)

    !         ! 判断场点是否在矩形的水平范围内
    !         if (rfp >= r_min .AND. rfp <= r_max) then
    !             ! 初始化距离数组
    !             do jj = 1, m-1
    !                 z_current = z_a(k, jj)
    !                 z_next = z_a(k + 1, jj + 1)
    !                 z_distances(jj) = min(abs(zfp - z_current), abs(zfp - z_next))
    !                 z_indices(jj) = jj
    !             end do

    !             ! 对距离进行排序
    !             do jj = 1, m-2
                    
    !                 do l = jj + 1, m-1
    !                     if (z_distances(jj) > z_distances(l)) then
                            
    !                         temp_distance = z_distances(jj)
    !                         z_distances(jj) = z_distances(l)
    !                         z_distances(l) = temp_distance
    !                         temp_index = z_indices(jj)
    !                         z_indices(jj) = z_indices(l)
    !                         z_indices(l) = temp_index
    !                     end if
    !                 end do
    !             end do

    !             ! 获取最接近的两条声线
    !             iter1 = z_indices(1)
    !             iter2 = z_indices(2)
    !             closest_ray1 = k
    !             closest_ray2 = k + 1
    !             found_ray1 = .true.
    !             found_ray2 = .true.
    !             exit
    !         end if
    !     end do
    ! end do

    ! ! 输出最接近的两条声线的信息
    ! if (found_ray1 .and. found_ray2) then
    !     write (*,*) "The closest two rays to the field point are:"
    !     write (*,*) "Ray 1: Number = ", closest_ray1, ", Iteration = ", iter1
    !     write (*,*) "Ray 2: Number = ", closest_ray2, ", Iteration = ", iter2
    ! else
    !     write (*,*) "Unable to find two closest rays to the field point."
    ! end if
    ! eigenray_num = closest_ray1
    ! eigenray_itern = iter1



    ! 根据特征声线数据，对场点数据进行插值。输出场点声压。
    !! TODO 测试场点位于海面、海底、源点很近距离时的极端情况。

    ! -----------------------计算幅值、声压------------------------------
    !! TODO: 计算幅值
    if (eigenray_num > 0 .AND. eigenray_num <= m .AND. eigenray_itern > 0 .AND. eigenray_itern <= n) then
        do k = eigenray_num, eigenray_num + 1
            A0(k, 1) = 0.25 * pi
            r_a(k, 1) = 0.01 !! XXX
            J(k, 1) = r_a(k, 1) / sin(theta) * (r_b(k, 1) - r_a(k, 1)) / dtheta0
            ! write (*,*) " J(k, 1) is ", J(k, 1)
            ! write (*,*) " r(k, 1) is ", r_a(k, 1)
            ! write (*,*) " A0 ", A0(k, 1)
            do i = 2, eigenray_itern

                J(k, i) = r_a(k, i) / sin(theta) * (r_b(k, i) - r_a(k, i)) / dtheta0
                ! write (*,*) " J(k, 1) is ", J(k, 1)
                A0(k, i + 1) = A0(k, 1) * (c_c(k, i) * J(k, 1) / c0 / J(k, i))
                ! write (*,*) " The A0 (", k, i, ") is", A0(k, i)
                p_auxi(k, i) = cmplx(0.0, omega * log(c_c(k, i) / c_c(k, 1)))
                ! write (*,*) " The p_auxi (", k, i, ") is", p_auxi(k, i)
                p(k, i) = A0(k, i) * exp(p_auxi(k, i))
                ! write(*,*) "Writing pressure: ", p(k, i)
                ! write (*,*) " ------------k i", k, i, "------------ "
                ! write (*,*) " The speed is ", c_c(k, i)
                ! write (*,*) " J(k, 1) is ", J(k, 1)
                ! write (*,*) " J(k, i) is ", J(k, i)
                ! write(*,*) "Writing pressure: ", p(k, i)
            end do
        end do 
        p_fp = 0.5 * (p(eigenray_num, eigenray_itern) + p(eigenray_num+1, eigenray_itern))
        write (*, '(A, F0.5, SP, F0.5, "i")') &
            &"The pressure of the field of point is : ", real(p_fp), aimag(p_fp)
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

    ! ! ----------------------------随声线轨迹变化的声速剖面---------------------------
    ! write (*,*) 'Writing to ......', SoundProfile
    ! open(unit=1, file=SoundProfile, status='replace', action='write')
    ! write(1, '(A)') '     Depth     SoundProfile'

    ! do i = 1, n  ! n为迭代次数
    !     write(1, '(8(F14.6, 1X))') z_a(1, i), c_a(1, i)
    ! end do 
    ! close(1)
    ! write (*,*) 'Results written to ', SoundProfile
    ! ! ----------------------------随声线轨迹变化的声速剖面------------------------


    ! ! ----------------------------写入声线轨迹---------------------------
    ! write (*,*) 'Writing to ......', RayTrace
    ! write (*,*) 'rsp:', rsp
    ! write (*,*) 'zsp:', zsp
    ! write (*,*) 'rfp:', rfp
    ! write (*,*) 'zfp:', zfp
    ! ! 打开文件进行写入
    ! open(unit=2, file=RayTrace, status='replace', action='write')

    ! ! 写入源点信息
    ! write(2, '(A)') 'The source point is : '
    ! write(2, '(2F14.6)') rsp, zsp

    ! ! 写入场点信息
    ! write(2, '(A)') 'The field point is : '
    ! write(2, '(2F14.6)') rfp, zfp

    ! write(2, '(A)') "The closest ray number k and k+1: "
    ! write(2, '(2I14)') eigenray_num, eigenray_num + 1


    ! ! 写入声线数据
    ! do k = 1, m ! m为声线束数量
    !     ! 在每条声线束的数据前添加一个标题行
    !     write(2, '(A, I0)') 'Ray Bundle ', k
    !     do i = 1, n  ! n为迭代次数
    !         write(2, '(2F14.6)') z_a(k, i), r_a(k, i)
    !     end do
    ! end do 
    ! ! 关闭文件
    ! close(2)
    ! write (*,*) 'Results written to ', RayTrace
    ! ! ----------------------------声线轨迹写入结束------------------------


    ! 调用子程序计算声速剖面并写入文件
    ! call calculate_and_write_sound_profile(depth, zp, cp, c0, SoundProfile)
    ! 获取结束时间
    call cpu_time(end_time)
    print *, 'Total running time: ', end_time - start_time, 's'
end program ray_tracing