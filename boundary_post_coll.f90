subroutine boundary
! This function creates the geometry.
! It defines lists with the points at the immersed boundaries (list_ib1 and list_ib2),
!   it also defines the weights of the interpolation (A_ib1 and A_ib2).

! ATTENTION, it raises an error if the number of points in x (Ngal(1,1)) and z (Ngal(2,1))
!   is not a multiple of the number of tiles (ntilex, ntilez)

! It is called by 'getbounds' in 'start.f90'

!                           | ------------------ Periodic unity --------------------|
!     0        *     *      *      *      *                           *      *      *      *      *                     
!    -1        *     *      *      *      *                           *      *      *      *      *                     
!    -2        *     *      *      *      *                           *      *      *      *      *                     
! Y  -3        *     *      *      *      *                           *      *      *      *      *                      
!    -4        *     *      *      *      *                           *      *      *      *      *                     
!   wall  *    *     *      *      *      *      *      *      *      *      *      *      *      *      *     *  
!                           1      2      3                           7      8     9(1)  10(2)  11(3)
!                                                      X,Z


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
  integer i,k,j,i1,j1,k1,ilist,ix,iz,dny,shift,grid
  integer nlist_ib1,nlist_ib2,ny11,ny21,ny12,ny22
  integer points_tile
  integer, allocatable:: list_ib1(:,:,:),list_ib2(:,:,:)
  real(8), allocatable:: A_ib1(:,:,:),A_ib2(:,:,:)

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

  dnx = Ngal(1,1)/ntilex  ! Width:  Number of points per tile the in streamwise direction
  dnz = Ngal(2,1)/ntilez  ! Width:  Number of points per tile the in spanwise   direction
  dny = -N(3,0)           ! Height: Number of points from the valley to the tip of a riblet.

! CREATE THE PATTERN. 
!  This part can be generalise calling different functions
!  Here we could implement a switch to use different geometries
! if (geometry_flag = 17) then
!   call boundary_blabla
!            ____________
! Post      |3 |      | 4| 
!           |__|      |__|
!           |            |
!           |__        __|
!           |1 |      | 2|
!           |__|______|__|

  ! First tile
  points_tile = npeakx*npeakz*dny   ! Total number of immersed boundary points in the pattern tile

  ! LOWER WALL

  nyu11 = N(4,0)+1                  ! Lower band. Lower boundary. First point with fluid points, one point above the wall
  nyu21 = 0                         ! Lower band. Roughness tips. Last  point with immersed boundaries.
  nyv11 = N(3,0)+1                  ! Lower band. Lower boundary. First point with fluid points, one point above the wall
  nyv21 = 0                         ! Lower band. Roughness tips. Last  point with immersed boundaries. Roughness tips.
  nlist_ib1 = ntilex*ntilez*points_tile    ! Number of solid points.     Riblets x height x width
  allocate(list_ib1(6,nlist_ib1,2))
  allocate(A_ib1(1,nlist_ib1,2))

  A_ib1 = 0d0    ! Weights of the immersed boundaries. Set to 0 because the most points don't need any interpolation

  ilist = 0
  do j = 1,dny-1 ! The last plane, j=dny, is calculated separately
    j1 = N(3,0)+j
    do i = 1,(npeakx+1)/2
      i1 = i
  ! Bottom-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1     ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
      end do
  ! Bottom-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1     ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
      end do
    end do
    do i=(npeakx+1)/2+1,npeakx
      i1 = dnx-npeakx+i
  ! Top-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1     ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
      end do
  ! Top-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1     ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
      end do
    end do
  end do
  ! Top plane in ugrid is calculated separately because of the weighting functions.
  do j = dny,dny
    j1 = N(3,0)+j
    do i = 1,(npeakx+1)/2
      i1 = i
  ! Bottom-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
        A_ib1   ( 1,ilist,ugrid) = (yu(j1)-yv(j1))/(yu(j1+1)-yv(j1))  ! So that u=0 at yv=0
      end do
  ! Bottom-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
        A_ib1   ( 1,ilist,ugrid) = (yu(j1)-yv(j1))/(yu(j1+1)-yv(j1))  ! So that u=0 at yv=0
      end do
    end do
    do i=(npeakx+1)/2+1,npeakx
      i1 = dnx-npeakx+i
  ! Top-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
        A_ib1   ( 1,ilist,ugrid) = (yu(j1)-yv(j1))/(yu(j1+1)-yv(j1))  ! So that u=0 at yv=0
      end do
  ! Top-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib1( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,ugrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,ugrid) = k1     ! ---- 
        list_ib1( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib1( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib1( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib1( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib1( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib1( 5,ilist,vgrid) = k1     ! ---- 
        list_ib1( 6,ilist,vgrid) = j1     ! ---- 
        A_ib1   ( 1,ilist,ugrid) = (yu(j1)-yv(j1))/(yu(j1+1)-yv(j1))  ! So that u=0 at yv=0
      end do
    end do
  end do

  if (ilist/=points_tile) then
    write(*,*) 'ERROR: ilist is not equal to points_tile'
    stop
  end if

  ! UPPER WALL

  nyu12 = nn+2                      ! Upper band. Roughness tips. First point with immersed boundaries.
  nyu22 = N(4,nband)                ! Upper band. Upper boundary. Last  point with fluid points, one point bellow the wall
  nyv12 = nn+1                      ! Upper band. Roughness tips. First point with immersed boundaries. Roughness tips
  nyv22 = N(3,nband)                ! Upper band. Upper boundary. Last  point with fluid points, one point bellow the wall
  nlist_ib2 = ntilex*ntilez*points_tile    ! Number of solid points.    Points per tile x tiles in x x tiles in z
  allocate(list_ib2(6,nlist_ib2,2))
  allocate(A_ib2   (1,nlist_ib2,2))
  A_ib2 = 0d0    ! Weights of the immersed boundaries. Set to 0 because the most points don't need any interpolation

  ! Bottom plane in ugrid is calculated separately because of the weighting functions.
  ilist = 0
  do j = 1,1
    j1 = nn+j
    do i = 1,(npeakx+1)/2
      i1 = i
  ! Bottom-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1     ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
        A_ib2   ( 1,ilist,ugrid) = (yu(j1+1)-yv(j1))/(yu(j1)-yv(j1))  ! So that u=0 at yv=0
      end do
  ! Bottom-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1     ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
        A_ib2   ( 1,ilist,ugrid) = (yu(j1+1)-yv(j1))/(yu(j1)-yv(j1))  ! So that u=0 at yv=0
      end do
    end do
    do i=(npeakx+1)/2+1,npeakx
      i1 = dnx-npeakx+i
  ! Top-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1     ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
        A_ib2   ( 1,ilist,ugrid) = (yu(j1+1)-yv(j1))/(yu(j1)-yv(j1))  ! So that u=0 at yv=0
      end do
  ! Top-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1     ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
        A_ib2   ( 1,ilist,ugrid) = (yu(j1+1)-yv(j1))/(yu(j1)-yv(j1))  ! So that u=0 at yv=0
      end do
    end do
  end do

  do j = 2,dny
    j1 = nn+j
    do i = 1,(npeakx+1)/2
      i1 = i
  ! Bottom-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
      end do
  ! Bottom-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
      end do
    end do
    do i=(npeakx+1)/2+1,npeakx
      i1 = dnx-npeakx+i
  ! Top-left corner of the post 
      do k = 1,(npeakz+1)/2
        ilist = ilist + 1
        k1 = k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
      end do
  ! Top-right corner of the post 
      do k = (npeakz+1)/2+1,npeakz
        ilist = ilist + 1
        k1 = dnz-npeakz+k
        list_ib2( 1,ilist,ugrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,ugrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,ugrid) = j1+1   ! j coordinate of the solid point 
        list_ib2( 4,ilist,ugrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,ugrid) = k1     ! ---- 
        list_ib2( 6,ilist,ugrid) = j1+1   ! ---- 
        list_ib2( 1,ilist,vgrid) = i1     ! i coordinate of the solid point
        list_ib2( 2,ilist,vgrid) = k1     ! k coordinate of the solid point
        list_ib2( 3,ilist,vgrid) = j1     ! j coordinate of the solid point 
        list_ib2( 4,ilist,vgrid) = i1     ! ---- Coordinates of the points involved in computing the 'solid' point
        list_ib2( 5,ilist,vgrid) = k1     ! ---- 
        list_ib2( 6,ilist,vgrid) = j1     ! ---- 
      end do
    end do
  end do

  if (ilist/=points_tile) then
    write(*,*) 'ERROR: ilist is not equal to points_tile'
    stop
  end if

! REPLICATE THE PATTERN
!   This section should be common for all geometries
  do grid = 1,2
    do ix = 1,ntilex
      do iz = 1,ntilez
        shift = points_tile*(iz-1) + points_tile*ntilez*(ix-1) 
        do ilist = 1,points_tile
          list_ib1(1,ilist+shift,grid) = list_ib1(1,ilist,grid) + dnx*(ix-1)     ! i coordinate
          list_ib1(2,ilist+shift,grid) = list_ib1(2,ilist,grid) + dnz*(iz-1)     ! k coordinate
          list_ib1(3,ilist+shift,grid) = list_ib1(3,ilist,grid)                  ! j coordinate
          list_ib1(4,ilist+shift,grid) = list_ib1(4,ilist,grid) + dnx*(ix-1)
          list_ib1(5,ilist+shift,grid) = list_ib1(5,ilist,grid) + dnz*(iz-1)
          list_ib1(6,ilist+shift,grid) = list_ib1(6,ilist,grid)
          list_ib2(1,ilist+shift,grid) = list_ib2(1,ilist,grid) + dnx*(ix-1)
          list_ib2(2,ilist+shift,grid) = list_ib2(2,ilist,grid) + dnz*(iz-1)
          list_ib2(3,ilist+shift,grid) = list_ib2(3,ilist,grid)
          list_ib2(4,ilist+shift,grid) = list_ib2(4,ilist,grid) + dnx*(ix-1)
          list_ib2(5,ilist+shift,grid) = list_ib2(5,ilist,grid) + dnz*(iz-1)
          list_ib2(6,ilist+shift,grid) = list_ib2(6,ilist,grid)
          A_ib1   (1,ilist+shift,grid) = A_ib1   (1,ilist,grid)
          A_ib2   (1,ilist+shift,grid) = A_ib2   (1,ilist,grid)
        end do
      end do
    end do
  end do


! Save the lists into a file
  open(10,file=trim(dirout)//'boundary_'//ext1//'x'//ext2//'x'//ext3//'.dat',form='unformatted')
  write(10) Lx,Ly,Lz
  write(10) Ngal,nlist_ib1,nlist_ib2,nyu11,nyu21,nyu12,nyu22,nyv11,nyv21,nyv12,nyv22
  write(10) list_ib1,list_ib2,A_ib1,A_ib2
  close(10)

  deallocate(list_ib1,list_ib2,A_ib1,A_ib2)

end subroutine 
