    ! ----------------------------写入声速剖面---------------------------
    write (*,*) 'Writing to ......', SoundProfile
    open(unit=1, file=SoundProfile, status='replace', action='write')
    write(1, '(A)') '     Depth     SoundProfile'

    do ip = 1, n_p
        write(1, '(8(F14.6, 1X))') z_p(ip), c_p(ip)
    end do 

    close(1)
    write (*,*) 'Results written to ', SoundProfile
    ! ----------------------------写入声速剖面------------------------