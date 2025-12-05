subroutine init_sl_stats(myid)

  use declaration
  implicit none
  
  integer myid

  allocate(du1dy_PL(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(du3dy_PL(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  allocate(u1_f_PL(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  allocate(u3_f_PL(N(1,bandPL(myid))+2,N(2,bandPL(myid)),jgal(vgrid,1)-1:jgal(vgrid,2)+1))
  
  allocate(bslip_u1_M(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  allocate(bslip_du1dy_M(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  allocate(bslip_u3_M(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  allocate(bslip_du3dy_M(N(1,bandPL(myid))+2,N(2,bandPL(myid))))

  bslip_u1_M=0d0
  bslip_du1dy_M=0d0
  bslip_u3_M=0d0
  bslip_du3dy_M=0d0
  
  istat_sl=0d0
  
   if (myid==0) then
     write(ext4,'(i5.5)') int(t)
     fnameimc = 'sl_stats_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
   end if
  
endsubroutine

subroutine sl_stats(u1,u3,myid,status,ierr)

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer      :: status(MPI_STATUS_SIZE),ierr,myid
  type(cfield) :: u1(sband:eband),u3(sband:eband)
  
  integer :: i,k,kk,j,column,iband
  
     
  !real(8) :: bslip_u1_f_bot(N(1,bandPL(myid))+2,N(2,bandPL(myid)))
  !real(8) :: bslip_u1_f_top(N(1,bandPL(myid))+2,N(2,bandPL(myid)))
  !real(8) :: bslip_u3_f_bot(N(1,bandPL(myid))+2,N(2,bandPL(myid)))
  !real(8) :: bslip_u3_f_top(N(1,bandPL(myid))+2,N(2,bandPL(myid)))
  
  !complex(8) :: u1_f   (N(1,bandPL(myid))/2+1,N(2,bandPL(myid)))
  !complex(8) :: du1dy_f(N(1,bandPL(myid))/2+1,N(2,bandPL(myid)))
  !complex(8) :: u3_f   (N(1,bandPL(myid))/2+1,N(2,bandPL(myid)))
  !complex(8) :: du3dy_f(N(1,bandPL(myid))/2+1,N(2,bandPL(myid)))
  
  
  !real(8),pointer :: bslip_u1_f_bot(:,:)
  !real(8),pointer :: bslip_u1_f_top(:,:)
  !real(8),pointer :: bslip_u3_f_bot(:,:)
  !real(8),pointer :: bslip_u3_f_top(:,:)
  
  complex(8),pointer :: u1_f   (:,:)
  complex(8),pointer :: du1dy_f   (:,:)
  complex(8),pointer :: u3_f   (:,:)
  complex(8),pointer :: du3dy_f   (:,:)
  
  !allocate(bslip_u1_f_bot(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  !allocate(bslip_u1_f_top(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  !allocate(bslip_u3_f_bot(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  !allocate(bslip_u3_f_top(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  
  allocate(u1_f   (N(1,bandPL(myid))/2+1,N(2,bandPL(myid))))
  allocate(du1dy_f(N(1,bandPL(myid))/2+1,N(2,bandPL(myid))))
  allocate(u3_f   (N(1,bandPL(myid))/2+1,N(2,bandPL(myid))))
  allocate(du3dy_f(N(1,bandPL(myid))/2+1,N(2,bandPL(myid))))

istat_sl=istat_sl+1
  
  do iband=1,2
    do column = 1,columns_num(iband,myid)
      j = jlim(1,vgrid,2)
      du1dy_columns(iband)%f(j,column)=(u1(iband)%f(j+1,column)-u1(iband)%f(j,column))*dthdyv(j)*ddthetavi
      du3dy_columns(iband)%f(j,column)=(u3(iband)%f(j+1,column)-u3(iband)%f(j,column))*dthdyv(j)*ddthetavi
    enddo
  enddo
  
  do iband=2,3
    do column = 1,columns_num(iband,myid)
      j = jlim(2,vgrid,2)
      du1dy_columns(iband)%f(j,column)=-(u1(iband)%f(j+1,column)-u1(iband)%f(j,column))*dthdyv(j)*ddthetavi
      du3dy_columns(iband)%f(j,column)=-(u3(iband)%f(j+1,column)-u3(iband)%f(j,column))*dthdyv(j)*ddthetavi
    enddo
  enddo
  
  du1dy_PL = 0d0
  du3dy_PL = 0d0
  u1_f_PL = 0d0
  u3_f_PL = 0d0
  
  call modes_to_planes_phys_lims_2(du1dy_PL,du1dy_columns,jlim(1,vgrid,2),jlim(1,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  call modes_to_planes_phys_lims_2(du3dy_PL,du3dy_columns,jlim(1,vgrid,2),jlim(1,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  call modes_to_planes_phys_lims_2(du1dy_PL,du1dy_columns,jlim(2,vgrid,2),jlim(2,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  call modes_to_planes_phys_lims_2(du3dy_PL,du3dy_columns,jlim(2,vgrid,2),jlim(2,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  
  call modes_to_planes_phys_lims_2(u1_f_PL,u1_itp,jlim(1,vgrid,2),jlim(1,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  call modes_to_planes_phys_lims_2(u3_f_PL,u3_itp,jlim(1,vgrid,2),jlim(1,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  call modes_to_planes_phys_lims_2(u1_f_PL,u1_itp,jlim(2,vgrid,2),jlim(2,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  call modes_to_planes_phys_lims_2(u3_f_PL,u3_itp,jlim(2,vgrid,2),jlim(2,vgrid,2),vgrid,myid,bandPL(myid),status,ierr)
  
    
  if(myid==0)then
    !bslip_u1_f_bot=u1_f_PL(:,:,jlim(1,vgrid,2))/du1dy_PL(:,:,jlim(1,vgrid,2))
    !bslip_u3_f_bot=u3_f_PL(:,:,jlim(1,vgrid,2))/du3dy_PL(:,:,jlim(1,vgrid,2))
 
    
  
  do i = 1,N(1,bandPL(myid))+2
      do k = 1,N(2,bandPL(myid))

        bslip_u1_M(i,k)    = bslip_u1_M(i,k) + u1_f_PL(i,k,jlim(1,vgrid,2))
        bslip_du1dy_M(i,k) = bslip_du1dy_M(i,k) + du1dy_PL(i,k,jlim(1,vgrid,2))
        bslip_u3_M(i,k)    = bslip_u3_M(i,k) + u3_f_PL(i,k,jlim(1,vgrid,2))
        bslip_du3dy_M(i,k) = bslip_du3dy_M(i,k) + du3dy_PL(i,k,jlim(1,vgrid,2))
       
      enddo
    enddo
    
    
  elseif(myid==np-1)then
    !bslip_u1_f_top=u1_f_PL(:,:,jlim(2,vgrid,2))/du1dy_PL(:,:,jlim(2,vgrid,2))
    !bslip_u3_f_top=u3_f_PL(:,:,jlim(2,vgrid,2))/du3dy_PL(:,:,jlim(2,vgrid,2))
 
    
  
  do i = 1,N(1,bandPL(myid))+2
      do k = 1,N(2,bandPL(myid))

        bslip_u1_M(i,k)    = bslip_u1_M(i,k) + u1_f_PL(i,k,jlim(2,vgrid,2))
        bslip_du1dy_M(i,k) = bslip_du1dy_M(i,k) + du1dy_PL(i,k,jlim(2,vgrid,2))
        bslip_u3_M(i,k)    = bslip_u3_M(i,k) + u3_f_PL(i,k,jlim(2,vgrid,2))
        bslip_du3dy_M(i,k) = bslip_du3dy_M(i,k) + du3dy_PL(i,k,jlim(2,vgrid,2))
       
      enddo
    enddo
    
    
    
  
  endif
  
  
  !deallocate(bslip_u1_f_bot)
  !deallocate(bslip_u1_f_top)
  !deallocate(bslip_u3_f_bot)
  !deallocate(bslip_u3_f_top)
  
  deallocate(u1_f)
  deallocate(du1dy_f)
  deallocate(u3_f)
  deallocate(du3dy_f)

endsubroutine

subroutine write_sl_stats(myid,status,ierr)

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr,myid
  
  !real(8) :: xreal(N(1,bandPL(myid))+2,N(2,bandPL(myid)))
  !real(8) :: xcomp(N(1,bandPL(myid))+2,N(2,bandPL(myid)))
  integer, allocatable :: dummint(:)
  
  real(8),pointer :: xreal(:,:)
  real(8),pointer :: xcomp(:,:)
     
  allocate(xreal(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
  allocate(xcomp(N(1,bandPL(myid))+2,N(2,bandPL(myid))))
 
  if(myid==np-1)then !send

    call MPI_SEND(bslip_u1_M(1,1),(N(1,bandPL(myid))+2)*N(2,bandPL(myid))   ,MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(bslip_du1dy_M(1,1),(N(1,bandPL(myid))+2)*N(2,bandPL(myid)),MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(bslip_u3_M(1,1),(N(1,bandPL(myid))+2)*N(2,bandPL(myid))   ,MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(bslip_du3dy_M(1,1),(N(1,bandPL(myid))+2)*N(2,bandPL(myid)),MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)

  elseif(myid==0)then !recieve
    fnameimc = trim(dirout)//'sl_stats_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
    open(10,file=fnameimc,form='unformatted')
    allocate(dummint(88))
    dummint = 0
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    deallocate(dummint)
    write(10) N
    write(10) istat_sl

    call recvbreal(xreal,bslip_u1_M(1,1),myid,status,ierr)   
    call recvbreal(xreal,bslip_du1dy_M(1,1),myid,status,ierr)
    call recvbreal(xreal,bslip_u3_M(1,1),myid,status,ierr)    
    call recvbreal(xreal,bslip_du3dy_M(1,1),myid,status,ierr)

  endif
  
  
  bslip_u1_M = 0d0
  bslip_du1dy_M = 0d0
  bslip_u3_M = 0d0
  bslip_du3dy_M = 0d0
  
  istat_sl=0d0

  deallocate(xreal)
  deallocate(xcomp)  

endsubroutine  


subroutine recvbreal(xread,xlocal,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    recvX    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr

  integer j,iproc,myid
  integer msize
  real(8) xread(N(1,bandPL(myid))+2,N(2,bandPL(myid)))
  real(8) xlocal(N(1,bandPL(myid))+2,N(2,bandPL(myid)))

  
  call MPI_RECV(xread         ,(N(1,bandPL(myid))+2)*N(2,bandPL(myid)),MPI_REAL8,np-1,np-1,MPI_COMM_WORLD,status,ierr)
  xlocal=xlocal+xread
  write(10) xlocal

end subroutine

subroutine recvbcomp(xread,xlocal,myid,status,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    recvX    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer status(MPI_STATUS_SIZE),ierr

  integer j,iproc,myid
  integer msize
  complex(8) xread(N(1,bandPL(myid))/2+1,N(2,bandPL(myid))/2)
  complex(8) xlocal(N(1,bandPL(myid))/2+1,N(2,bandPL(myid))/2)


  call MPI_RECV(xread,(N(1,bandPL(myid))/2+1)*N(2,bandPL(myid))/2,MPI_COMPLEX8,np-1,np-1,MPI_COMM_WORLD,status,ierr)
  xlocal=xlocal+xread
  write(10) xlocal

end subroutine
