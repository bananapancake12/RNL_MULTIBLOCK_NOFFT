include 'rft_buff.f08'
include 'cft_buff.f08'

program four_to_phys
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!  FOUR TO PHYS  !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
  integer nx,nz,ny,nx2,nz2,nxold,nzold,nband,iter
  integer j,jd,iband,id
  real(8), allocatable:: buffSR1(:,:),buffSR2(:,:),buffSR1OLD(:,:)
  real(8), allocatable:: buffR(:),buffC(:)
  integer, allocatable:: N(:,:),Ngal(:,:),dummint(:)
  real(8) t,Re,alp,bet,mpgx,dumm_y
  real(8), allocatable:: yu(:),dthdyu(:)
  real(8), allocatable:: yv(:),dthdyv(:)
  real(8) dthetai
  character*50 fnameima,fnameimb,filout
  character*4 ext1,ext2,ext3
  character*5 ext4

  open(40,file='input_FtP.txt',form='formatted')
  read(40,20) nx
  read(40,20) nz
  read(40,20) ny
  read(40,10) t

  write(ext1,'(i4.4)') nx
  write(ext2,'(i4.4)') nz
  write(ext3,'(i4.4)') ny
  write(ext4,'(i5.5)') int(10d0*t)


!!!!!!!!!!!!!!!!!!!!!!!! u1 !!!!!!!!!!!!!!!!!!!!!!!!!!!
  filout='output/u1_'
  fnameima=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  filout='output/u1_phys_'
  fnameimb=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  open(10,file=fnameima,form='unformatted')
  open(20,file=fnameimb,form='unformatted')
  allocate(dummint(88))
  read (10) t,Re,alp,bet,mpgx,nband,iter,dummint
  write(20) t,Re,alp,bet,mpgx,nband,iter,dummint
!  deallocate(dummint)
  allocate(N   (4,0:nband+1))
  allocate(Ngal(4,0:nband+1))
  read (10) N
  ! Defining Ngal from N
  Ngal = N
  do iband = 1,nband
    Ngal(1,iband) = N(1,iband)!*3/2
    Ngal(2,iband) = N(2,iband)!*3/2
  end do
  write(20) Ngal
  allocate(yu    (N(4,0):N(4,nband)+1))
  allocate(dthdyu(N(4,0):N(4,nband)+1))
  read (10) yu,dthetai,dthdyu
  write(20) yu,dthetai,dthdyu

  do iband=1,nband
    ! Initialisation
    nx =N(1,iband)+2
    nz =N(2,iband)
    nx2=Ngal(1,iband)+2
    nz2=Ngal(2,iband)
    allocate(buffSR1(nx ,nz ))                                       ! Buffer to store u 
    allocate(buffSR2(nx2,nz2))                                       ! Buffer to store u 
    allocate(buffR(2*Ngal(1,iband)+18),buffC(10*Ngal(2,iband)+19))   ! Buffer for the fft's
    call rfti(Ngal(1,iband),buffR)
    call cfti(Ngal(2,iband),buffC)
    j=N(4,iband-1)
    if (iband==1) then  ! Band 1 (wall)
      read (10) jd,id,nx,nz,dumm_y,buffSR1       ! Read from file. Stored in SR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    else                ! Bands 2 & 3 (interface)
      dumm_y=yu(N(4,iband-1))
      call resize(buffSR1,buffSR1OLD,nx,nz,nxold,nzold)
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) j ,1 ,nx2,nz2,dumm_y,buffSR2
      deallocate(buffSR1OLD)
    end if
    do j=N(4,iband-1)+1,N(4,iband)  ! The rest
      read (10) jd,id,nx,nz,dumm_y,buffSR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    end do
    j=N(4,iband)+1
    if (iband==nband) then  ! Band 3 (wall)
      read (10) jd,id,nx,nz,dumm_y,buffSR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    else                    ! Bands 1 & 2 (interface). Saves the next value and the size, and this is used in the next band with the old number of points. Used for instance in y-x plots
      nxold=nx
      nzold=nz
      allocate(buffSR1OLD(nxold,nzold))
!      call rft(buffSR1,nx,nz,-1,buffR)
!      call cft(buffSR1,nx,2,nx/2,-1,buffC)
      buffSR1OLD=buffSR1
    end if
    deallocate(buffSR1)
    deallocate(buffSR2)
    deallocate(buffR,buffC)
  end do
  write(*,*) 'Writing ',fnameimb

!!!!!!!!!!!!!!!!!!!!!!!! u2 !!!!!!!!!!!!!!!!!!!!!!!!!!!
  filout='output/u2_'
  fnameima=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  filout='output/u2_phys_'
  fnameimb=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  open(10,file=fnameima,form='unformatted')
  open(20,file=fnameimb,form='unformatted')
  read (10)
  write(20) t,Re,alp,bet,mpgx,nband,iter,dummint
  read (10)
  write(20) Ngal
  allocate(yv    (N(3,0):N(3,nband)+1))
  allocate(dthdyv(N(3,0):N(3,nband)+1))
  read (10) yv,dthetai,dthdyv
  write(20) yv,dthetai,dthdyv

  do iband=1,nband
    nx =N(1,iband)+2
    nz =N(2,iband)
    nx2=Ngal(1,iband)+2
    nz2=Ngal(2,iband)
    allocate(buffSR1(nx ,nz ))
    allocate(buffSR2(nx2,nz2))
    allocate(buffR(2*Ngal(1,iband)+18),buffC(10*Ngal(2,iband)+19))
    call rfti(Ngal(1,iband),buffR)
    call cfti(Ngal(2,iband),buffC)
    j=N(3,iband-1)
    if (iband==1) then
      read (10) jd,id,nx,nz,dumm_y,buffSR1       ! Read from file. Stored in SR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    else
      dumm_y=yv(N(3,iband-1))
      call resize(buffSR1,buffSR1OLD,nx,nz,nxold,nzold)
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
      deallocate(buffSR1OLD)
    end if
    do j=N(3,iband-1)+1,N(3,iband)
      read (10) jd,id,nx,nz,dumm_y,buffSR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    end do
    j=N(3,iband)+1
    if (iband==nband) then
      read (10) jd,id,nx,nz,dumm_y,buffSR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    else
      nxold=nx
      nzold=nz
      allocate(buffSR1OLD(nxold,nzold))
!      call rft(buffSR1,nx,nz,-1,buffR)
!      call cft(buffSR1,nx,2,nx/2,-1,buffC)
      buffSR1OLD=buffSR1
    end if
    deallocate(buffSR1)
    deallocate(buffSR2)
    deallocate(buffR,buffC)
  end do
  write(*,*) 'Writing ',fnameimb

!!!!!!!!!!!!!!!!!!!!!!!! u3 !!!!!!!!!!!!!!!!!!!!!!!!!!!
  filout='output/u3_'
  fnameima=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  filout='output/u3_phys_'
  fnameimb=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  open(10,file=fnameima,form='unformatted')
  open(20,file=fnameimb,form='unformatted')
  read (10)
  write(20) t,Re,alp,bet,mpgx,nband,iter,dummint
  read (10)
  write(20) Ngal
  read (10)
  write(20) yu,dthetai,dthdyu

  do iband=1,nband
    nx =N(1,iband)+2
    nz =N(2,iband)
    nx2=Ngal(1,iband)+2
    nz2=Ngal(2,iband)
    allocate(buffSR1(nx ,nz ))
    allocate(buffSR2(nx2,nz2))
    allocate(buffR(2*Ngal(1,iband)+18),buffC(10*Ngal(2,iband)+19))
    call rfti(Ngal(1,iband),buffR)
    call cfti(Ngal(2,iband),buffC)
    j=N(4,iband-1)
    if (iband==1) then  ! Band 1 (wall)
      read (10) jd,id,nx,nz,dumm_y,buffSR1       ! Read from file. Stored in SR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    else                ! Bands 2 & 3 (interface)
      dumm_y=yu(N(4,iband-1))
      call resize(buffSR1,buffSR1OLD,nx,nz,nxold,nzold)
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) j ,1 ,nx2,nz2,dumm_y,buffSR2
      deallocate(buffSR1OLD)
    end if
    do j=N(4,iband-1)+1,N(4,iband)  ! The rest
      read (10) jd,id,nx,nz,dumm_y,buffSR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    end do
    j=N(4,iband)+1
    if (iband==nband) then  ! Band 3 (wall)
      read (10) jd,id,nx,nz,dumm_y,buffSR1
      call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
      call cft(buffSR2,nx2,2,nx2/2,1,buffC)
      call rft(buffSR2,nx2,  nz2,  1,buffR)
      write(20) jd,id,nx2,nz2,dumm_y,buffSR2
    else                    ! Bands 1 & 2 (interface). Saves the next value and the size, and this is used in the next band with the old number of points. Used for instance in y-x plots
      nxold=nx
      nzold=nz
      allocate(buffSR1OLD(nxold,nzold))
!      call rft(buffSR1,nx,nz,-1,buffR)
!      call cft(buffSR1,nx,2,nx/2,-1,buffC)
      buffSR1OLD=buffSR1
    end if
    deallocate(buffSR1)
    deallocate(buffSR2)
    deallocate(buffR,buffC)
  end do
  write(*,*) 'Writing ',fnameimb

!!!!!!!!!!!!!!!!!!!!!!!! p  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  filout='output/p_'
  fnameima=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  filout='output/p_phys_'
  fnameimb=filout(1:index(filout,' ')-1)//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
  open(10,file=fnameima,form='unformatted')
  open(20,file=fnameimb,form='unformatted')
  read (10)
  write(20) t,Re,alp,bet,mpgx,nband,iter,dummint
  deallocate(dummint)
  read (10)
  write(20) Ngal
  read (10)
  write(20) yu,dthetai,dthdyu

  do iband=1,nband
    nx =N(1,iband)+2
    nz =N(2,iband)
    nx2=Ngal(1,iband)+2
    nz2=Ngal(2,iband)
    allocate(buffSR1(nx ,nz ))
    allocate(buffSR2(nx2,nz2))
    allocate(buffR(2*Ngal(1,iband)+18),buffC(10*Ngal(2,iband)+19))
    call rfti(Ngal(1,iband),buffR)
    call cfti(Ngal(2,iband),buffC)
    if (iband==1) then  ! Band 1 (wall) !C! No pressure point for 1st point
      do j=N(4,0)+1,N(4,1)-1  ! The rest
        read (10) jd,id,nx,nz,dumm_y,buffSR1
        call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
        call cft(buffSR2,nx2,2,nx2/2,1,buffC)
        call rft(buffSR2,nx2,  nz2,  1,buffR)
        write(20) jd,id,nx2,nz2,dumm_y,buffSR2
      end do
    elseif(iband==3)then!else                ! Bands 2 & 3 (interface)
      do j=N(4,2)+2,N(4,3)  ! The rest
        read (10) jd,id,nx,nz,dumm_y,buffSR1
        call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
        call cft(buffSR2,nx2,2,nx2/2,1,buffC)
        call rft(buffSR2,nx2,  nz2,  1,buffR)
        write(20) jd,id,nx2,nz2,dumm_y,buffSR2
      end do
    else ! iband = 2
      do j=N(4,1),N(4,2)+1  ! The rest
        read (10) jd,id,nx,nz,dumm_y,buffSR1
        call resize(buffSR2,buffSR1,nx2,nz2,nx,nz) ! Add zeros N->Ngal. Stored in SR2
        call cft(buffSR2,nx2,2,nx2/2,1,buffC)
        call rft(buffSR2,nx2,  nz2,  1,buffR)
        write(20) jd,id,nx2,nz2,dumm_y,buffSR2
      end do
    end if
    deallocate(buffSR1)
    deallocate(buffSR2)
    deallocate(buffR,buffC)
  end do
  write(*,*) 'Writing ',fnameimb

10 FORMAT(7X,D10.1)
20 FORMAT(7X,I10)

end program

subroutine resize(buffSR,buffSROLD,nx,nz,nxold,nzold)

  implicit none
  integer i,k,dk,nx,nz,nxold,nzold
  real(8) buffSR(nx,nz),buffSROLD(nxold,nzold)

  dk=nzold-nz

  buffSR=0d0
  do k=1,min(nz/2,nzold/2)
    do i=1,min(nx,nxold)
      buffSR(i,k)=buffSROLD(i,k)
    end do
  end do
!   do k=nz-min(nz/2,nzold/2)+2,nz
  do k=nz-min(nz/2,nzold/2)+1,nz
    do i=1,min(nx,nxold)
      buffSR(i,k)=buffSROLD(i,k+dk)
    end do
  end do

end subroutine

