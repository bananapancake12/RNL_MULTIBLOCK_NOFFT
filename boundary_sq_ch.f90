subroutine boundary
! This function creates the ribbed geometry.
! It defines lists with the points at the immersed boundaries (list_ib1 and list_ib2),
!  it also defines the weights of the interpolation (A_ib1 and A_ib2).

! ATTENTION, it raises an error if the number of points in z (Ngal(2,1)) is not a multiple of the number of ribs (nribs)

! It is called by 'getbounds' in 'start.f90'

!                         | ------------------ Periodic unity --------------------|
!   0        *     *      *      *      *                           *      *      *      *      *                     
!  -1        *     *      *      *      *                           *      *      *      *      *                     
!  -2        *     *      *      *      *                           *      *      *      *      *                     
!  -3        *     *      *      *      *                           *      *      *      *      *                      
!  -4        *     *      *      *      *                           *      *      *      *      *                     
!       *    *     *      *      *      *      *      *      *      *      *      *      *      *      *     *  
!                         1      2      3                           7      8     9(1)  10(2)  11(3)


!   ------------------------ ny22
!   \  /\  /\  /\  /\  /\  /
!    \/  \/  \/  \/  \/  \/
!   ------------------------ ny12
!
!
!
!   ------------------------ ny21
!    /\  /\  /\  /\  /\  /\
!   /  \/  \/  \/  \/  \/  \
!   ------------------------ ny11



  use declaration
  implicit none
  integer k,j,j1,k1,ilist,ic,dny
  !integer nlist_ib1,nlist_ib2,nyu11,nyu21,nyu12,nyu22,nyv11,nyv21,nyv12,nyv22
  integer nlist_ib1,nlist_ib2
  integer, allocatable:: list_ib1(:,:),list_ib2(:,:)
  real(8), allocatable:: A_ib1(:),A_ib2(:)

  ! Check that the number of grid points in x and z is a multiple of the number of tiles 
  if (ntilex.eq.0) then
    write(*,*) 'ERROR: ntilex equals 0'
    stop
  end if
  if (ntilez.eq.0) then
    write(*,*) 'ERROR: ntilez equals 0'
    stop
  end if
  if (mod(Ngal(1,1),ntilez)/=0) then
    write(*,*) 'ERROR: nx not multiple of ntilex'
    stop
  end if
  if (mod(Ngal(2,1),ntilez)/=0) then
    write(*,*) 'ERROR: nz not multiple of ntilez'
    stop
  end if

  dnx = Ngal(1,1)/ntilex  ! Width:  Number of points per tile in streamwise direction
  dnz = Ngal(2,1)/ntilez  ! Width:  Number of points per tile in spanwise   direction
  dny = -N(3,0)           ! Height: Number of points from the valley to the tip of a riblet.

! LOWER WALL

  nyu11=N(4,0)+1                  ! Lower band. Lower boundary. First point with fluid points (if no dummy planes)
  nyu21=0                         ! Lower band. Riblets' tips.  Last  point with immersed boundaries.
  nyv11=N(3,0)+1                  ! Lower band. Lower boundary. First point with fluid points (if no dummy planes)
  nyv21=0                         ! Lower band. Riblets' tips.  Last  point with immersed boundaries.
  !nlist_ib1=nribs*dny*npeak      ! Number of solid points.     Riblets x height x width
  nlist_ib1=0
  allocate(list_ib1(6,nlist_ib1))
  allocate(A_ib1(nlist_ib1))
!   
!   ! Iterates over the 'solid' points
!   ! Stores in 'list_ib1' the coordinates j and k of solid points in the lower band,
!   !  and also the coordinates of the points involved in computing its state.
!   do ic=1,nribs
!     do j=1,dny
!       j1=N(3,0)+j
!       ! Left half of the riblet
!       do k=1,(npeak+1)/2
!         ilist=(ic-1)*(dny*npeak)+(j-1)*npeak+k
!         k1=(ic-1)*dnz+k
!         list_ib1(1,ilist)=k1     ! k coordinate of the solid point
!         list_ib1(2,ilist)=j1     ! j coordinate of the solid point
!         list_ib1(3,ilist)=k1     ! ---- 
!         list_ib1(4,ilist)=j1     ! ---- Coordinates of the points involved in computing the 'solid' point
!         list_ib1(5,ilist)=k1     ! ---- There are only two point because the problem is solved in a 2D approach
!         list_ib1(6,ilist)=j1     ! ----
!       end do
!       ! Right half of the riblet
!       do k=(npeak+1)/2+1,npeak
!         ilist=(ic-1)*(dny*npeak)+(j-1)*npeak+k
!         k1=ic*dnz-(npeak-k)
!         list_ib1(1,ilist)=k1
!         list_ib1(2,ilist)=j1
!         list_ib1(3,ilist)=k1
!         list_ib1(4,ilist)=j1
!         list_ib1(5,ilist)=k1
!         list_ib1(6,ilist)=j1
!       end do
!     end do
!   end do
  A_ib1=0d0    ! Weights of the immersed boundaries. Set to 0 because the points fall on top of the boundaries.
               !  No need of interpolation.

! UPPER WALL
! Similar code to the upper part

  nyu12=nn+2                      ! Upper band. Riblets' tips.
  nyu22=N(4,nband)                ! Upper band. Upper boundary.
  nyv12=nn+1                      ! Upper band. Riblets' tips.
  nyv22=N(3,nband)                ! Upper band. Upper boundary.
  !nlist_ib2=nribs*dny*npeak
  nlist_ib2=0
  allocate(list_ib2(6,nlist_ib2))
  allocate(A_ib2(nlist_ib2))
!   do ic=1,nribs
!     do j=1,dny
!       j1=nn+j
!       do k=1,(npeak+1)/2
!         ilist=(ic-1)*(dny*npeak)+(j-1)*npeak+k
!         k1=(ic-1)*dnz+k
!         list_ib2(1,ilist)=k1
!         list_ib2(2,ilist)=j1
!         list_ib2(3,ilist)=k1
!         list_ib2(4,ilist)=j1
!         list_ib2(5,ilist)=k1
!         list_ib2(6,ilist)=j1
!       end do
!       do k=(npeak+1)/2+1,npeak
!         ilist=(ic-1)*(dny*npeak)+(j-1)*npeak+k
!         k1=ic*dnz-(npeak-k)
!         list_ib2(1,ilist)=k1
!         list_ib2(2,ilist)=j1
!         list_ib2(3,ilist)=k1
!         list_ib2(4,ilist)=j1
!         list_ib2(5,ilist)=k1
!         list_ib2(6,ilist)=j1
!       end do
!     end do
!   end do
  A_ib2=0d0

  open(10,file='output/boundary_'//ext1//'x'//ext2//'x'//ext3//'.dat',form='unformatted')
  write(10) Lx,Ly,Lz,nband
  write(10) Ngal,nlist_ib1,nlist_ib2,nyu11,nyu21,nyu12,nyu22,nyv11,nyv21,nyv12,nyv22
  write(10) list_ib1,list_ib2,A_ib1,A_ib2
  close(10)

  deallocate(list_ib1,list_ib2,A_ib1,A_ib2)

end subroutine
