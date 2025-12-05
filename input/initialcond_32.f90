subroutine list_fill(nx, nz, scalefac, data_in, data_out)
    implicit none
    integer, intent(in) :: nx, nz
    real(8), intent(in) :: scalefac        ! was integer
    integer :: i, j, nxout, nzout, nxcopy, jmax
    real(8), intent(in)  :: data_in(*)
    real(8), intent(out) :: data_out(*)

    ! nx already includes +2 pad; scale base (nx-2) and add pad back
    nxout = int( (nx - 2) * scalefac ) + 2
    nzout = int(  nz       * scalefac )

    do i = 1, nxout * nzout
        data_out(i) = 0.0d0
    end do

    ! copy only the overlapping region
    nxcopy = min(nx, nxout)
    jmax   = min(nz, nzout)
    do j = 1, jmax
        data_out( (j-1)*nxout + 1 : (j-1)*nxout + nxcopy ) = &
        data_in ( (j-1)*nx    + 1 : (j-1)*nx    + nxcopy )
    end do
end subroutine


subroutine read_write(filein, fileout, NPL, isv, isp, scalefac)
    implicit none

    character(len=*), intent(in) :: filein, fileout
    integer, intent(in) :: NPL, isv, isp
    real(8), intent(in) :: scalefac

    real(8) :: t, Re, alp, bet, mpgx, dthetai, dummRe
    integer :: nband, iter
    integer :: N(4,0:4)
    integer, allocatable :: dummint(:)

    ! plane buffers
    real(8), allocatable :: plane_in(:), plane_out(:)

    ! per-plane header fields
    integer :: j, jin, dummI, nxin, nzin
    integer :: nxout_hdr, nzout_hdr
    integer :: start

    ! y-grid arrays (sizes preserved, we’re only scaling x & z)
    real(8), allocatable :: yu(:), dthdyu(:)

    write(6,*) "starting readwrite", trim(fileout)

    allocate(dummint(88))

    ! ---- read input file header ----
    open(40, file=filein, form='unformatted')
    read(40) t, Re, alp, bet, mpgx, nband, iter, dummint
    read(40) N

    ! grid arrays on y (keep as in source file)
    if (isv == 0) then
        allocate(yu    ( N(4,0):N(4,3)+1 ))
        allocate(dthdyu( N(4,0):N(4,3)+1 ))
    else
        allocate(yu    ( N(3,0):N(3,3)+1 ))   ! (leave as in your code)
        allocate(dthdyu( N(3,0):N(3,3)+1 ))
    end if
    read(40) yu, dthetai, dthdyu

    ! input plane buffer (original sizes)
    allocate(plane_in( (N(1,3)+2) * N(2,3) ))

    ! output plane buffer (scaled sizes)
    allocate(plane_out( ( int(N(1,3)*scalefac) + 2 ) * int(N(2,3)*scalefac) ))

    ! ---- update N for output (base sizes only; +2 is implicit per plane header) ----
    N(1,1) = int(N(1,1) * scalefac)
    N(1,2) = int(N(1,2) * scalefac)
    N(1,3) = int(N(1,3) * scalefac)

    N(2,1) = int(N(2,1) * scalefac)
    N(2,2) = int(N(2,2) * scalefac)
    N(2,3) = int(N(2,3) * scalefac)

    ! ---- open output and write header ----
    open(50, file=fileout, form='unformatted', status='replace')
    write(50) t, Re, alp, bet, mpgx, nband, iter, dummint
    write(50) N
    write(50) yu, dthetai, dthdyu

    if (isp == 0) then
        start = 1
    else
        start = 2
    end if

    ! ---- planes 1 .. NPL/2 ----
    do j = 1, NPL/2
        read(40) jin, dummI, nxin, nzin, dummRe, plane_in

        call list_fill(nxin, nzin, scalefac, plane_in, plane_out)

        nxout_hdr = int( (nxin - 2) * scalefac ) + 2
        nzout_hdr = int(  nzin       * scalefac )

        write(50) jin, dummI, nxout_hdr, nzout_hdr, dummRe, plane_out
    end do

    write(6,*) "nxin=", nxin, "nzin=", nzin

    ! ---- planes NPL/2+1 .. NPL ----
    do j = NPL/2 + 1, NPL
        read(40) jin, dummI, nxin, nzin, dummRe, plane_in

        call list_fill(nxin, nzin, scalefac, plane_in, plane_out)

        nxout_hdr = int( (nxin - 2) * scalefac ) + 2
        nzout_hdr = int(  nzin       * scalefac )

        write(50) jin, dummI, nxout_hdr, nzout_hdr, dummRe, plane_out
    end do

    write(6,*) "nxin=", nxin, "nzin=", nzin
    write(6,*) "N=", N

    close(40)
    close(50)

    write(6,*) "exiting readwrite", trim(fileout)

end subroutine read_write

program initial_conditions
    implicit none

    ! allocate variables
    character(len=:), allocatable :: input
    character(len=:), allocatable :: output
    character(len=:), allocatable :: dirin
    character(len=:), allocatable :: dirout
    integer :: planes
    real(8) :: scalefac

    ! ---- choose scaling ----
    scalefac = 0.5d0     ! downscale by 2 in x & z (keeps low-k content)

    dirin  = '.'
    dirout = '.'

    ! ---- u1 ----
planes = 154
input  = dirin  // '/u1_0064x0064x0153_t14299.dat'
output = dirout // '/u1_0032x0032x0153_t00000.dat'
call read_write(input, output, planes, 0, 0, scalefac)
write(6,*) "u1 done!"

! ---- u2 (v-grid metadata in header section; leave isv=1) ----
planes = 153
input  = dirin  // '/u2_0064x0064x0153_t14299.dat'
output = dirout // '/u2_0032x0032x0153_t00000.dat'
call read_write(input, output, planes, 1, 0, scalefac)
write(6,*) "u2 done!"

! ---- u3 ----
planes = 154
input  = dirin  // '/u3_0064x0064x0153_t14299.dat'
output = dirout // '/u3_0032x0032x0153_t00000.dat'
call read_write(input, output, planes, 0, 0, scalefac)
write(6,*) "u3 done!"

! ---- p (isp=1) ----
planes = 152
input  = dirin  // '/p_0064x0064x0153_t14299.dat'
output = dirout // '/p_0032x0032x0153_t00000.dat'
call read_write(input, output, planes, 0, 1, scalefac)
write(6,*) "p done!"

end program initial_conditions
