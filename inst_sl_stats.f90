subroutine inst_sl_stats(u1,u3,myid,status,ierr)

  use declaration
  implicit none

  include 'mpif.h'             ! MPI variables
  integer      :: status(MPI_STATUS_SIZE),ierr,myid
  type(cfield) :: u1(sband:eband),u3(sband:eband)
  
  integer :: i,k,kk,j,column,iband
  integer, allocatable :: dummint(:)
  
  !real(8) :: xreal(N(1,bandPL(myid))+2,N(2,bandPL(myid)))

  real(8),pointer :: xreal(:,:)
  allocate(xreal(N(1,bandPL(myid))+2,N(2,bandPL(myid))))

  
   if (myid==0) then
     write(ext4,'(i5.5)') int(t)
     fnameimc = 'sl_stats_inst_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
   end if
    
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
  
  if(myid==np-1)then !send
    call MPI_SEND(u1_f_PL(1,1,jlim(2,vgrid,2)),(N(1,bandPL(myid))+2)*N(2,bandPL(myid))   ,MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(du1dy_PL(1,1,jlim(2,vgrid,2)),(N(1,bandPL(myid))+2)*N(2,bandPL(myid)),MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(u3_f_PL(1,1,jlim(2,vgrid,2)),(N(1,bandPL(myid))+2)*N(2,bandPL(myid))   ,MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)
    call MPI_SEND(du3dy_PL(1,1,jlim(2,vgrid,2)),(N(1,bandPL(myid))+2)*N(2,bandPL(myid)),MPI_REAL8   ,0,myid,MPI_COMM_WORLD,ierr)

  elseif(myid==0)then !recieve
    write(ext4,'(i5.5)') int(10d0*t)!int(t)!
    fnameimc = trim(dirout)//'sl_stats_inst_'//ext1//'x'//ext2//'x'//ext3//'_'//ext4//'.dat'
    open(10,file=fnameimc,form='unformatted')
    allocate(dummint(88))
    dummint = 0
    write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
    deallocate(dummint)
    write(10) N
    
    call recvbreal_inst(xreal,u1_f_PL(1,1,jlim(1,vgrid,2)),myid,status,ierr)   
    call recvbreal_inst(xreal,du1dy_PL(1,1,jlim(1,vgrid,2)),myid,status,ierr)
    call recvbreal_inst(xreal,u3_f_PL(1,1,jlim(1,vgrid,2)),myid,status,ierr)    
    call recvbreal_inst(xreal,du3dy_PL(1,1,jlim(1,vgrid,2)),myid,status,ierr)

  close(10)
  
  endif

  deallocate(xreal)

endsubroutine

subroutine recvbreal_inst(xread,xlocal,myid,status,ierr)
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

  
  call MPI_RECV(xread,(N(1,bandPL(myid))+2)*N(2,bandPL(myid)),MPI_REAL8,np-1,np-1,MPI_COMM_WORLD,status,ierr)

  write(10) xlocal
  write(10) xread

end subroutine