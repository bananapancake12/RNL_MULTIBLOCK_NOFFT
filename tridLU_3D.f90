! !  ! ****** package of utilities for tridiagonal LU factorisation *****
! !  !                                   RGM ETSI Aeronauticos 2008
! !  ! *************************************************************


!  ! ****** package of utilities for tridiagonal LU factorisation *****
!  !                                   RGM ETSI Aeronauticos 2008
!  ! *************************************************************


! *************************************************************************************************************** !
!
! Last modified by C 20/05/2016
!   - Set up for SHS
!   - Tridiag set up for FFT in middle of bands 1 and 3
!   - Split SolU and SolV as V BC can be solved in Fourier space
!
! TODO
!   - More efficient FFT
!   - DG needs to contatin option for non-zero v velocity at wall
!
! *************************************************************************************************************** !
!
! Contains: LUsolU    - Called from littleharsh for solveU : calls LU_build and LU_dec
!                       - Solves the (1-L) matrix for the U and W velocities
!                       - Decomposed into (1-Lxz)(1-y)
!                       - Contains FFT to solve wall in physical space and interface in Fourier space
!           LUsolV    - Called from littleharsh for solveV : calls LU_build and LU_dec
!                       - Solves the (1-L) matrix for the V velocities
!                       - Decomposed into (1-Lxz)(1-y)
!                       - No FFT all Fourier
!           LUsolP    - Called from littleharsh for solveP : calls LU_buildP and LU_decP
!                       - Solves the lplacian of the pressure (divergance of the gradient)
!           LUsol0    - Called in flowrate_corr (part of const. flow rate) : calls LU_build0 and LU_dec0
!                       - Solves mean mode to calculate pressure gradient for constant mass flow rate
!           LU_build  - Called from Lusol
!                       - Builds the (1-L) matrices. Different for u/w and v
!           LU_buildP - Called from start (with allocation)
!                       - Builds the DG matrix for the pressure
!                       - Only called once as no dt dependence
!           LU_build0 - Called from Lusol0
!                       - Builds LHS for LUsol0
!           LU_dec    - Called from LUsol
!                       - LU Decomposition
!           LU_decP   - Called from LUsolP
!                       - LU Decomposition
!                       - Only called once as no dt dependence
!           LU_dec0   - Called from LUsol0
!                       - LU Decomposition
!
! *************************************************************************************************************** !



subroutine LUsolU_W(u,rhsu,w,rhsw,a,grid,iband,myid)
   !----------------------------------------------------------------------*
   !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
   !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
   ! INPUT m_j number of unknowns
   !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
   !       a(2,:,i,k) inverse of the diagonal of L_j
   !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
   !       (note that the diagonal of U_j is 1)
   !       x(nystart:nyend,i,k)   rhs of the problem
   ! OUTPUT
   !       all untouched but the solution 'x'
   !----------------------------------------------------------------------*

   use declaration
   implicit none
  
   include 'mpif.h'             ! MPI variables
   integer status(MPI_STATUS_SIZE),ierr

   integer    :: i,k,j,column,iband,grid
   integer    :: iopt,myid
   type(cfield)  u(sband:eband)
   type(cfield)  w(sband:eband)
   type(cfield)  rhsu(sband:eband)
   type(cfield)  rhsw(sband:eband)
    
   real(8)    :: a(3,jlim(1,grid,iband):jlim(2,grid,iband))


   ! -Original
      do column = 1,columns_num(iband,myid)
         call LU_buildU(jlim(1,grid,iband),jlim(2,grid,iband),column,myid,iband,a)
         call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),myid,a)

        u(iband)%f(jlim(1,grid,iband),column)=rhsu(iband)%f(jlim(1,grid,iband),column)*a(2,jlim(1,grid,iband))
        w(iband)%f(jlim(1,grid,iband),column)=rhsw(iband)%f(jlim(1,grid,iband),column)*a(2,jlim(1,grid,iband))
        
        do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
          u(iband)%f(j,column) = (rhsu(iband)%f(j,column)-a(1,j)*u(iband)%f(j-1,column))*a(2,j)
        end do

        do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
          u(iband)%f(j,column) =  u(iband)%f(j,column)-a(3,j)*u(iband)%f(j+1,column)
        end do

        do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
         w(iband)%f(j,column) = (rhsw(iband)%f(j,column)-a(1,j)*w(iband)%f(j-1,column))*a(2,j)
       end do

       do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
         w(iband)%f(j,column) =  w(iband)%f(j,column)-a(3,j)*w(iband)%f(j+1,column)
       end do
      end do

      if(myid==0) then 
        write(6,*) "solving u"
      end if 

      ! if(myid==0) then
      !   do column = 1,10
      !     write(6,*) u(2)%f(7,column)
      !   end do
      ! end if 


      if(myid==0) then 
        write(6,*) "solving w"
      end if 

      ! if(myid==0) then
      !   do column = 1,10
      !     write(6,*) w(2)%f(7,column)
      !   end do
      ! end if 

  
end subroutine


subroutine LUsolV(u,rhsu,a,grid,iband,myid)
   !----------------------------------------------------------------------*
   !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
   !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
   ! INPUT m_j number of unknowns
   !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
   !       a(2,:,i,k) inverse of the diagonal of L_j
   !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
   !       (note that the diagonal of U_j is 1)
   !       x(nystart:nyend,i,k)   rhs of the problem
   ! OUTPUT
   !       all untouched but the solution 'x'
   !----------------------------------------------------------------------*

   use declaration
   implicit none
  
   include 'mpif.h'             ! MPI variables
   integer status(MPI_STATUS_SIZE),ierr

   integer    :: i,k,j,column,iband,grid
   integer    :: iopt,myid
   type(cfield) u(sband:eband)
   type(cfield) rhsu(sband:eband)
    
   real(8)    :: a(3,jlim(1,grid,iband):jlim(2,grid,iband))

   ! -Original
      do column = 1,columns_num(iband,myid)
         call LU_buildV(jlim(1,grid,iband),jlim(2,grid,iband),column,myid,iband,a)
         call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),myid,a)

        u(iband)%f(jlim(1,grid,iband),column)=rhsu(iband)%f(jlim(1,grid,iband),column)*a(2,jlim(1,grid,iband))

        do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
          u(iband)%f(j,column) = (rhsu(iband)%f(j,column)-a(1,j)*u(iband)%f(j-1,column))*a(2,j)
        end do

        do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
          u(iband)%f(j,column) =  u(iband)%f(j,column)-a(3,j)*u(iband)%f(j+1,column)
        end do
      end do

      if (myid ==0) then
        write(6,*) "solving v"
      end if 

      ! if(myid==0) then
      !   do column = 1,10
      !     write(6,*) u(2)%f(7,column)
      !   end do
      ! end if 
  

end subroutine


subroutine LUsolP(x,myid,iband,nystart,nyend)
   !----------------------------------------------------------------------*
   !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
   !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
   ! INPUT m_j number of unknowns
   !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
   !       a(2,:,i,k) inverse of the diagonal of L_j
   !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
   !       (note that the diagonal of U_j is 1)
   !       x(nystart:nyend,i,k)   rhs of the problem
   ! OUTPUT
   !       all untouched but the solution 'x'
   !----------------------------------------------------------------------*

   use declaration
   implicit none
   integer i,k,j,column,myid,iband
   integer nystart,nyend
   complex(8) x(nystart:nyend,columns_num(iband,myid))
   
   do column = 1,columns_num(iband,myid)
      i = columns_i(column,iband,myid)
      k = columns_k(column,iband,myid)
      x(nystart,column) = x(nystart,column)*DG(iband)%f_dg(2,nystart,column)
      do j = nystart+1,nyend
         x(j,column) = (x(j,column)-DG(iband)%f_dg(1,j,column)*x(j-1,column))*DG(iband)%f_dg(2,j,column)
      end do
      do j = nyend-1,nystart,-1
         x(j,column) =  x(j,column)-DG(iband)%f_dg(3,j,column)*x(j+1,column)
      end do

   end do


end subroutine

subroutine LUsol0(x,nystart,nyend)
   !----------------------------------------------------------------------*
   !                  LUsol for real (0,0) mode
   !----------------------------------------------------------------------*

   use declaration
   implicit none
   integer j
   integer nystart,nyend
   real(8) x(nystart:nyend)
   real(8), allocatable:: a(:,:)
   allocate(a(3,nystart:nyend))

   call LU_build0(nystart,nyend,a)
   call LU_dec0  (nystart,nyend,a)

   x(nystart) = x(nystart)*a(2,nystart)
   do j = nystart+1,nyend
      x(j) = (x(j)-a(1,j)*x(j-1))*a(2,j)
   end do
   do j = nyend-1,nystart,-1
      x(j) = x(j)-a(3,j)*x(j+1)
   end do

   deallocate(a)

end subroutine


subroutine LU_buildU(nystart,nyend,column,myid,iband,a)
   !-------------------------------------------------------!
   !       specifies original values of a(1:3,j,i,k)       !
   !-------------------------------------------------------!

   use declaration
   implicit none
   integer i,k,j,grid,myid,column,iband
   integer nystart,nyend
   real(8) k2x,k2z,beta
   real(8) a(3,nystart:nyend)

   beta = bRK(kRK)*dt/Re
  

   !!!!!!!!!!!!!!!!!!!     u velocity:       !!!!!!!!!!!!!!!!!!!
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid)
        k2x = k2F_x(i)
        k2z = k2F_z(k)

          a(1,nystart) = 0d0
          a(2,nystart) = 1d0  
          a(3,nystart) = gridweighting(iband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0

        do j = nystart+1,nyend-1
          a(1,j) =    -beta* dyu2i(1,j)
          a(2,j) = 1d0-beta*(dyu2i(2,j)+k2x+k2z)
          a(3,j) =    -beta* dyu2i(3,j)
        end do

          a(1,nyend) = gridweighting(iband,2)!0d0 !Free-shear -1 no-slip 0
          a(2,nyend) = 1d0
          a(3,nyend) = 0d0
end subroutine


subroutine LU_buildV(nystart,nyend,column,myid,iband,a)
   !-------------------------------------------------------!
   !       specifies original values of a(1:3,j,i,k)       !
   !-------------------------------------------------------!

   use declaration
   implicit none
   integer i,k,j,grid,myid,column,iband
   integer nystart,nyend
   real(8) k2x,k2z,beta
   real(8) a(3,nystart:nyend)

   beta = bRK(kRK)*dt/Re
  


    !!!!!!!!!!!!!!!!!!!     v velocity:       !!!!!!!!!!!!!!!!!!!  
        i = columns_i(column,iband,myid)
        k = columns_k(column,iband,myid)
        k2x = k2F_x(i)
        k2z = k2F_z(k)

          a(1,nystart) = 0d0
          a(2,nystart) = 1d0  
          a(3,nystart) = 0d0

        do j = nystart+1,nyend-1
          a(1,j) =    -beta* dyv2i(1,j)
          a(2,j) = 1d0-beta*(dyv2i(2,j)+k2x+k2z)
          a(3,j) =    -beta* dyv2i(3,j)
        end do

          a(1,nyend) = 0d0
          a(2,nyend) = 1d0
          a(3,nyend) = 0d0
end subroutine


!  
subroutine LU_buildP(nystart,nyend,myid,iband,a)
   !-------------------------------------------------------!
   !       specifies original values of a(1:3,j,i,k)       !
   !-------------------------------------------------------!

   use declaration
   implicit none
   type(rfield_dg)  a(sband:eband)
   integer column,i,k,j,myid,iband
   integer nystart,nyend
   real(8) k2x,k2z

   real(8) D_vec_b(nystart-1:nyend+1)
   real(8) G_vec_b(nystart-1:nyend+1)
   real(8) D_vec_t(nystart-1:nyend+1)
   real(8) G_vec_t(nystart-1:nyend+1)

   do j = nystart,nyend
      D_vec_b(j-1) = -ddthetavi*dthdyu(j)
      G_vec_b(j  ) = -ddthetavi*dthdyv(j)
      D_vec_t(j  ) =  ddthetavi*dthdyu(j)
      G_vec_t(j+1) =  ddthetavi*dthdyv(j)
   end do
  
  
  
  
   do column = 1,columns_num(iband,myid)
      i = columns_i(column,iband,myid)
      k = columns_k(column,iband,myid)
    
      !For exact wavenumbers
      k2x = k2F_x(i)
      k2z = k2F_z(k)
    
      !    !For 2nd order centrered difference wavenumbers
      !    k2x = k1F_x(i)*k1F_x(i)
      !    k2z = k1F_z(k)*k1F_z(k)
   
      a(iband)%f_dg(2,nystart,column) = (D_vec_t(nystart)*G_vec_b(nystart  )) + k2x + k2z
      a(iband)%f_dg(3,nystart,column) = (D_vec_t(nystart)*G_vec_t(nystart+1))
      a(iband)%f_dg(2,nyend,  column) = (D_vec_b(nyend-1)*G_vec_t(nyend    )) + k2x + k2z
      a(iband)%f_dg(1,nyend,  column) = (D_vec_b(nyend-1)*G_vec_b(nyend-1  ))
      do j = nystart+1,nyend-1
         a(iband)%f_dg(2,j,column) = (D_vec_b(j-1)*G_vec_t(j  )) + (D_vec_t(j)*G_vec_b(j)) + k2x + k2z
         a(iband)%f_dg(1,j,column) = (D_vec_b(j-1)*G_vec_b(j-1))
         a(iband)%f_dg(3,j,column) = (D_vec_t(j  )*G_vec_t(j+1))
      end do

   end do
  
   if(myid==0)then
      if(iband==midband)then
         a(iband)%f_dg(1,nystart,1) = 0d0
         a(iband)%f_dg(2,nystart,1) = 1d0/((yu(1+1)-yu(1)))**2 !C!
         a(iband)%f_dg(3,nystart,1) = 0d0
      end if
   end if

   ! For modified wavenumbers, need BC on last pressure mode (as k2x=k2z=0)
   !  do column = 1,columns_num(iband,myid)
   !  i = columns_i(column,iband,myid)
   !  k = columns_k(column,iband,myid)
   !  if(abs(k1F_x(i)*k1F_x(i))<10e-10.and.abs(k1F_z(k)*k1F_z(k))<10e-10)then
   !  a(iband)%f_dg(1,nystart,column) = 0d0
   !  a(iband)%f_dg(2,nystart,column) = 1d0/((yu(1+1)-yu(1)))**2 !C!
   !  a(iband)%f_dg(3,nystart,column) = 0d0
   !  endif
   !  enddo
  

   call LU_decP(nystart,nyend,columns_num(iband,myid),a(iband)%f_dg(:,nystart:nyend,1:columns_num(iband,myid)),iband,myid)

  
end subroutine

subroutine LU_build0(nystart,nyend,a)
   !----------------------------------------------------------------------*
   !                  LU_build for real (0,0) mode
   !----------------------------------------------------------------------*

   use declaration
   implicit none
   integer j
   integer nystart,nyend
   real(8) a(3,nystart:nyend)
   real(8) beta

   !!!!!!!!!!!!!!!!!!     velocity     !!!!!!!!!!!!!!!!!!
   beta = bRK(kRK)*dt/Re
   a(1,nystart) = 0d0
   a(2,nystart) = 1d0
   a(3,nystart) = gridweighting(midband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
   do j = nystart+1,nyend-1
      a(1,j) =    -beta*dyu2i(1,j)
      a(2,j) = 1d0-beta*dyu2i(2,j)
      a(3,j) =    -beta*dyu2i(3,j)
   end do
   a(1,nyend) = gridweighting(midband,2)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
   a(2,nyend) = 1d0
   a(3,nyend) = 0d0

end subroutine

subroutine LU_dec(nystart,nyend,myid,a)
   !----------------------------------------------------------------------*
   !         GIVEN 'a_j' FOR j=nystart:nyend
   !         PERFORMS a_j=L_j*U_j
   !         WITH 'a_j' allmost TRIDIAGONAL
   ! INPUT  nystart:nyend matrix leading dimension
   !        iopt selects original value of 'a'
   !        a_j is given as a tridiagonal with
   !        a(1,nystart:nyend-1) lower diagonal,
   !        a(2,:) main diagonal,
   !        a(3,nystart+1:nyend) upper diagonal
   !        the matrix 'a_j' is assembeld as:
   !  a_j(2,1) a_j(3,1)    0        0         0      .............       0
   !  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
   !     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
   !     0     ...................................................       0
   !     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
   !     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
   !
   ! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
   !         a(2,nystart:nyend) inverse of the diagonal of L_j
   !         a(3,nystart:nyend-1) upper diagonal of U_j
   !         (notice that the diagonal of U_j is 1)
   !----------------------------------------------------------------------*

   use declaration
   implicit none
   integer j,iband,myid
   integer nystart,nyend
   real(8) a(3,nystart:nyend)
 
   a(2,nystart) = 1d0/a(2,nystart)
   do j = nystart+1,nyend
     a(3,j-1) = a(3,j-1)*a(2,j-1)
     a(2,j) = 1d0/(a(2,j)-a(1,j)*a(3,j-1))
   end do

end subroutine

!  
subroutine LU_decP(nystart,nyend,columns,a,iband,myid)
   !----------------------------------------------------------------------*
   !         GIVEN 'a_j' FOR j=nystart:nyend
   !         PERFORMS a_j=L_j*U_j
   !         WITH 'a_j' allmost TRIDIAGONAL
   ! INPUT  nystart:nyend matrix leading dimension
   !        iopt selects original value of 'a'
   !        a_j is given as a tridiagonal with
   !        a(1,nystart:nyend-1) lower diagonal,
   !        a(2,:) main diagonal,
   !        a(3,nystart+1:nyend) upper diagonal
   !        the matrix 'a_j' is assembeld as:
   !  a_j(2,1) a_j(3,1)    0        0         0      .............       0
   !  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
   !     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
   !     0     ...................................................       0
   !     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
   !     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
   !
   ! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
   !         a(2,nystart:nyend) inverse of the diagonal of L_j
   !         a(3,nystart:nyend-1) upper diagonal of U_j
   !         (notice that the diagonal of U_j is 1)
   !----------------------------------------------------------------------*

   use declaration
   implicit none

   integer j,column
   integer columns,iband,myid
   integer nystart,nyend
   real(8) a(3,nystart:nyend,columns)
  
   do column = 1,columns
      a(2,nystart,column) = 1d0/a(2,nystart,column)
      do j = nystart+1,nyend
         a(3,j-1,column) = a(3,j-1,column)*a(2,j-1,column)
         a(2,j  ,column) = 1d0/(a(2,j,column) - a(1,j,column)*a(3,j-1,column))
      end do
   end do

end subroutine

subroutine LU_dec0(nystart,nyend,a)
   !----------------------------------------------------------------------*
   !                  LU_dec for real (0,0) mode
   !----------------------------------------------------------------------*

   use declaration
   implicit none
   integer j
   integer nystart,nyend
   real(8) a(3,nystart:nyend)

   a(2,nystart) = 1d0/a(2,nystart)
   do j = nystart+1,nyend
      a(3,j-1) = a(3,j-1)*a(2,j-1)
      a(2,j  ) = 1d0/(a(2,j)-a(1,j)*a(3,j-1))
   end do

end subroutine











! ! *************************************************************************************************************** !
! !
! ! Last modified by C 20/05/2016
! !   - Set up for SHS
! !   - Tridiag set up for FFT in middle of bands 1 and 3
! !   - Split SolU and SolV as V BC can be solved in Fourier space
! !
! ! TODO
! !   - More efficient FFT
! !   - DG needs to contatin option for non-zero v velocity at wall
! !
! ! *************************************************************************************************************** !
! !
! ! Contains: LUsolU    - Called from littleharsh for solveU : calls LU_build and LU_dec
! !                       - Solves the (1-L) matrix for the U and W velocities
! !                       - Decomposed into (1-Lxz)(1-y)
! !                       - Contains FFT to solve wall in physical space and interface in Fourier space
! !           LUsolV    - Called from littleharsh for solveV : calls LU_build and LU_dec
! !                       - Solves the (1-L) matrix for the V velocities
! !                       - Decomposed into (1-Lxz)(1-y)
! !                       - No FFT all Fourier
! !           LUsolP    - Called from littleharsh for solveP : calls LU_buildP and LU_decP
! !                       - Solves the lplacian of the pressure (divergance of the gradient)
! !           LUsol0    - Called in flowrate_corr (part of const. flow rate) : calls LU_build0 and LU_dec0
! !                       - Solves mean mode to calculate pressure gradient for constant mass flow rate
! !           LU_build  - Called from Lusol
! !                       - Builds the (1-L) matrices. Different for u/w and v
! !           LU_buildP - Called from start (with allocation)
! !                       - Builds the DG matrix for the pressure
! !                       - Only called once as no dt dependence
! !           LU_build0 - Called from Lusol0
! !                       - Builds LHS for LUsol0
! !           LU_dec    - Called from LUsol
! !                       - LU Decomposition
! !           LU_decP   - Called from LUsolP
! !                       - LU Decomposition
! !                       - Only called once as no dt dependence
! !           LU_dec0   - Called from LUsol0
! !                       - LU Decomposition
! !
! ! *************************************************************************************************************** !


! subroutine LUsolU_W(u,rhsu,w,rhsw,a,grid,iband,myid)
!    !----------------------------------------------------------------------*
!    !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
!    !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
!    ! INPUT m_j number of unknowns
!    !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
!    !       a(2,:,i,k) inverse of the diagonal of L_j
!    !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
!    !       (note that the diagonal of U_j is 1)
!    !       x(nystart:nyend,i,k)   rhs of the problem
!    ! OUTPUT
!    !       all untouched but the solution 'x'
!    !----------------------------------------------------------------------*

!    use declaration
!    implicit none
  
!    include 'mpif.h'             ! MPI variables
!    integer status(MPI_STATUS_SIZE),ierr

!    integer    :: i,k,j,column,iband,grid
!    integer    :: iopt,myid
!    type(cfield)  u(sband:eband)
!    type(cfield)  w(sband:eband)
!    type(cfield)  rhsu(sband:eband)
!    type(cfield)  rhsw(sband:eband)
    
!    real(8)    :: a(3,jlim(1,grid,iband):jlim(2,grid,iband))


!    ! -Original
!       do column = 1,columns_num(iband,myid)
!          call LU_buildU(jlim(1,grid,iband),jlim(2,grid,iband),column,myid,iband,a)
!          call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),myid,a)

!         u(iband)%f(jlim(1,grid,iband),column)=rhsu(iband)%f(jlim(1,grid,iband),column)*a(2,jlim(1,grid,iband))
!         w(iband)%f(jlim(1,grid,iband),column)=rhsw(iband)%f(jlim(1,grid,iband),column)*a(2,jlim(1,grid,iband))
        
!         do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
!           u(iband)%f(j,column) = (rhsu(iband)%f(j,column)-a(1,j)*u(iband)%f(j-1,column))*a(2,j)
!         end do

!         do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
!           u(iband)%f(j,column) =  u(iband)%f(j,column)-a(3,j)*u(iband)%f(j+1,column)
!         end do

!         do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
!          w(iband)%f(j,column) = (rhsw(iband)%f(j,column)-a(1,j)*w(iband)%f(j-1,column))*a(2,j)
!        end do

!        do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
!          w(iband)%f(j,column) =  w(iband)%f(j,column)-a(3,j)*w(iband)%f(j+1,column)
!        end do
!       end do

  
! end subroutine


! subroutine LUsolV(u,rhsu,a,grid,iband,myid)
!    !----------------------------------------------------------------------*
!    !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
!    !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
!    ! INPUT m_j number of unknowns
!    !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
!    !       a(2,:,i,k) inverse of the diagonal of L_j
!    !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
!    !       (note that the diagonal of U_j is 1)
!    !       x(nystart:nyend,i,k)   rhs of the problem
!    ! OUTPUT
!    !       all untouched but the solution 'x'
!    !----------------------------------------------------------------------*

!    use declaration
!    implicit none
  
!    include 'mpif.h'             ! MPI variables
!    integer status(MPI_STATUS_SIZE),ierr

!    integer    :: i,k,j,column,iband,grid
!    integer    :: iopt,myid
!    type(cfield) u(sband:eband)
!    type(cfield) rhsu(sband:eband)
    
!    real(8)    :: a(3,jlim(1,grid,iband):jlim(2,grid,iband))

!    ! -Original
!       do column = 1,columns_num(iband,myid)
!          call LU_buildV(jlim(1,grid,iband),jlim(2,grid,iband),column,myid,iband,a)
!          call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),myid,a)

!         u(iband)%f(jlim(1,grid,iband),column)=rhsu(iband)%f(jlim(1,grid,iband),column)*a(2,jlim(1,grid,iband))

!         do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
!           u(iband)%f(j,column) = (rhsu(iband)%f(j,column)-a(1,j)*u(iband)%f(j-1,column))*a(2,j)
!         end do

!         do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
!           u(iband)%f(j,column) =  u(iband)%f(j,column)-a(3,j)*u(iband)%f(j+1,column)
!         end do
!       end do
  

! end subroutine


! ! subroutine LUsolU(x,u,grid,myid)
! ! !----------------------------------------------------------------------*
! ! !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
! ! !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
! ! ! INPUT m_j number of unknowns
! ! !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
! ! !       a(2,:,i,k) inverse of the diagonal of L_j
! ! !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
! ! !       (note that the diagonal of U_j is 1)
! ! !       x(nystart:nyend,i,k)   rhs of the problem
! ! ! OUTPUT
! ! !       all untouched but the solution 'x'
! ! !----------------------------------------------------------------------*

! !   use declaration
! !   implicit none
  
! !   include 'mpif.h'             ! MPI variables
! !   integer status(MPI_STATUS_SIZE),ierr

! !   integer    :: i,k,j,column,iband,grid
! !   integer    :: iopt,myid
! !   type(cfield) x(sband:eband)
! !   type(cfield) u(sband:eband)
    
! !   real(8) , allocatable:: a(:,:,:,:)

  
! !   !a(diag,column,j,iband)
! !   allocate(a(3,maxval(columns_num),jlim(1,grid,2):jlim(2,grid,2),nband))

! !   do iband = sband,eband
! !     call LU_build(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,a)
! !     call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,a(:,:,:,:))
! !   enddo
  

! ! ! -Original   
! !   do iband = sband,eband
! !     do column = 1,columns_num(iband,myid)
! !       x(iband)%f(jlim(1,grid,iband),column)=x(iband)%f(jlim(1,grid,iband),column)*a(2,column,jlim(1,grid,iband),iband)
! !       do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
! !         x(iband)%f(j,column) = (x(iband)%f(j,column)-a(1,column,j,iband)*x(iband)%f(j-1,column))*a(2,column,j,iband)
! !       end do
! !       do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
! !         x(iband)%f(j,column) =  x(iband)%f(j,column)-a(3,column,j,iband)*x(iband)%f(j+1,column)
! !       end do
! !     end do
! !   enddo

! !   deallocate(a)
  
! ! end subroutine


! subroutine LUsolU(x,grid,myid)
! !----------------------------------------------------------------------*
! !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
! !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
! ! INPUT m_j number of unknowns
! !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
! !       a(2,:,i,k) inverse of the diagonal of L_j
! !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
! !       (note that the diagonal of U_j is 1)
! !       x(nystart:nyend,i,k)   rhs of the problem
! ! OUTPUT
! !       all untouched but the solution 'x'
! !----------------------------------------------------------------------*

!   use declaration
!   implicit none
  
!   include 'mpif.h'             ! MPI variables
!   integer status(MPI_STATUS_SIZE),ierr

!   integer    :: i,k,j,column,iband,grid
!   integer    :: iopt,myid
!   type(cfield) x(sband:eband)
!   real(8), allocatable:: xPL(:,:,:)
    
!   real(8) , allocatable:: ay(:,:,:,:)
!   real(8) , allocatable:: axz(:,:)

  
!   !physlim_bot = 7   !? Last physical point - ~ 5 points?
!   !physlim_top = 144 !NO WALL! 145 !? Last physical point

!   !ay(BC,diag,j,iband)
!   allocate(ay(2,3,jlim(1,grid,2):jlim(2,grid,2),nband))
!   allocate(axz(maxval(columns_num),nband))

!   do iband = sband,eband
!     call LU_build(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,axz,ay)
!     call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,ay(1,:,:,:))
!     call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,ay(2,:,:,:))
!   enddo
  
! ! ! FFT CHECK
! !   if(myid==0.and.grid==ugrid)then
! !   do j=jlim(1,ugrid,2),25
! !     print *, j, yu(j), abs(ay(1,1,j,2)-ay(2,1,j,2)),abs(ay(1,2,j,2)-ay(2,2,j,2)),abs(ay(1,3,j,2)-ay(2,3,j,2))
! !   enddo
! !   endif
! !   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! !   stop
  
 
! !  ay(2,:,:,:)=ay(1,:,:,:) !NS - Wall check
 
!    do iband = sband,eband
!     do column = 1,columns_num(iband,myid)
!       do j=jlim(1,grid,iband)+1,jlim(2,grid,iband)-1 !C! Don't include first and last as BC's
!         x(iband)%f(j,column)=x(iband)%f(j,column)/axz(column,iband)
        
!       enddo
!     enddo
!   enddo
  
!   allocate(xPL(N(1,bandPL_FFT(myid))+2,N(2,bandPL_FFT(myid)),limPL_FFT(grid,1,myid):limPL_FFT(grid,2,myid)))

  
! !----------------------------------------------------------------------*
! !  
! !  ____________________________________________
! ! |                   ¦                        |
! ! |   ^       |       ¦      ^        |        |
! ! |   | (5)   | (6)   ¦      | (5)    | (6)    |
! ! |   |       v       ¦      |        v        |
! ! |-------------------¦------------------------|   FFT Interface
! ! |   ^       |       ¦      ^        |        |
! ! |   |       |       ¦      | (4)    | (7)    |
! ! |   |       |       ¦      |        v        |
! ! |_ _|_ _ _ _|_ _ _ _¦________________________|   Band interface
! ! |   |       |       |
! ! |   |       |       |
! ! |   |       |       |
! ! |   | (3)   | (8)   |
! ! |   |       |       |
! ! |   |       |       | 
! ! |_ _|_ _ _ _|_ _ _ _|________________________    Band interface
! ! |   |       |       ¦                        |
! ! |   |       |       ¦      ^        |        |
! ! |   |       |       ¦      | (2)    | (9)    |
! ! |   |       v       ¦      |        v        |
! ! |-------------------¦------------------------|   FFT Interface
! ! |   ^       |       ¦      ^        |        |
! ! |   | (1)   | (10)  ¦      | (1)    | (10)   |
! ! |   |       v       ¦      |        v        |
! ! |___________________¦________________________|
! !
! !
! !----------------------------------------------------------------------*
  
!   xPL=0d0
!   !FFT physlim_bot
!   call modes_to_planes_phys_lims(xPL,x,jlim(1,grid,2),physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)
!   call modes_to_planes_phys_lims(xPL,x,physlim_top+1,jlim(2,grid,2),grid,myid,bandPL_FFT(myid),status,ierr)
!   do j=limPL_FFT(grid,1,myid),limPL_FFT(grid,2,myid)
!     if(j<=physlim_bot)then
!       call four_to_phys_N(xPL(1,1,j),bandPL_FFT(myid))
!     endif
!     if(j>=physlim_top+1)then
!       call four_to_phys_N(xPL(1,1,j),bandPL_FFT(myid))
!     endif
!   enddo
!   call planes_to_modes_phys_lims(x,xPL,jlim(1,grid,2),physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)
!   call planes_to_modes_phys_lims(x,xPL,physlim_top+1,jlim(2,grid,2),grid,myid,bandPL_FFT(myid),status,ierr)
   
!   !First point in band
!   !Physical
!   do iband = 1,2 ! (1)
!     do column = 1,columns_num(iband,myid)
!       !NS or FS
! x(iband)%f(jlim(1,grid,iband),column) = real(x(iband)%f(jlim(1,grid,iband),column)*&
! &ay(planeBC(2*columns_i(column,iband,myid)+1,columns_k(column,iband,myid)),2,jlim(1,grid,iband),iband))&
! &+im*aimag(x(iband)%f(jlim(1,grid,iband),column)*&
! &ay(planeBC(2*columns_i(column,iband,myid)+2,columns_k(column,iband,myid)),2,jlim(1,grid,iband),iband))
!     enddo
!   enddo
   
!    !First point in band
!    !Fourier
!    iband = 3 ! (4)
!    do column = 1,columns_num(iband,myid)
!      !Taking NS
!      x(iband)%f(jlim(1,grid,iband),column) = x(iband)%f(jlim(1,grid,iband),column)*ay(1,2,jlim(1,grid,iband),iband)
!    enddo
  
!  !Heading UP
!  !Physical
!   do iband = 1,2 ! (1)
!     do j = jlim(1,grid,iband)+1,physlim_bot
!       do column = 1,columns_num(iband,myid)     
!         !NS or FS   
! x(iband)%f(j,column) = real((x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+1,columns_k(column,iband,myid)),1,j,iband)*&
! &x(iband)%f(j-1,column))*ay(planeBC(2*columns_i(column,iband,myid)+1,columns_k(column,iband,myid)),2,j,iband))+&
! &im*aimag((x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+2,columns_k(column,iband,myid)),1,j,iband)*&
! &x(iband)%f(j-1,column))*ay(planeBC(2*columns_i(column,iband,myid)+2,columns_k(column,iband,myid)),2,j,iband))
!       enddo
!     enddo
!   enddo
  
!   xPL=0d0
!   !FFT physlim_bot
!   call modes_to_planes_phys_lims(xPL,x,physlim_bot,physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)
!   do j=limPL_FFT(grid,1,myid),limPL_FFT(grid,2,myid)
!     if(j==physlim_bot)then
!       call phys_to_four_N(xPL(1,1,j),bandPL_FFT(myid))
!     endif
!   enddo
!   call planes_to_modes_phys_lims(x,xPL,physlim_bot,physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)

!   !Block 1 upto phys lim
!   !Fourier
!   iband = 1 ! (2)
!   do j = physlim_bot+1,jlim(2,grid,iband)
!     do column = 1,columns_num(iband,myid)
!       !Taking NS
!       x(iband)%f(j,column) = (x(iband)%f(j,column)-ay(1,1,j,iband)*x(iband)%f(j-1,column))*ay(1,2,j,iband)
!     enddo
!   enddo
   
!   !Block 2 upto phys lim
!   !Fourier
!   iband = 2 ! (3)
!   do j = physlim_bot+1,physlim_top
!     do column = 1,columns_num(iband,myid)
!       !Taking NS
!       x(iband)%f(j,column) = (x(iband)%f(j,column)-ay(1,1,j,iband)*x(iband)%f(j-1,column))*ay(1,2,j,iband)
!     enddo
!   enddo
  
!   !Block 3 upto phys lim
!   !Fourier
!   iband = 3 ! (4)
!   do j = jlim(1,grid,iband)+1,physlim_top
!     do column = 1,columns_num(iband,myid)
!       !Taking NS
!       x(iband)%f(j,column) = (x(iband)%f(j,column)-ay(1,1,j,iband)*x(iband)%f(j-1,column))*ay(1,2,j,iband)
!     enddo
!   enddo
  
!   xPL=0d0
!   !FFT physlim_top
!   call modes_to_planes_phys_lims(xPL,x,physlim_top,physlim_top,grid,myid,bandPL_FFT(myid),status,ierr)
!   do j=limPL_FFT(grid,1,myid),limPL_FFT(grid,2,myid)
!     if(j==physlim_top)then
!       call four_to_phys_N(xPL(1,1,j),bandPL_FFT(myid))
!     endif
!   enddo
!   call planes_to_modes_phys_lims(x,xPL,physlim_top,physlim_top,grid,myid,bandPL_FFT(myid),status,ierr)
   
!   !Blocks 2 & 3 upto top
!   !Physical
!   do iband = 2,3 ! (5)
!     do j = physlim_top+1,jlim(2,grid,iband)
!       do column = 1,columns_num(iband,myid)
!         !NS or FS
! x(iband)%f(j,column) = real((x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+1,columns_k(column,iband,myid)),1,j,iband)*&
! &x(iband)%f(j-1,column))*ay(planeBC(2*columns_i(column,iband,myid)+1,columns_k(column,iband,myid)),2,j,iband))&
! &+im*aimag((x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+2,columns_k(column,iband,myid)),1,j,iband)*&
! &x(iband)%f(j-1,column))*ay(planeBC(2*columns_i(column,iband,myid)+2,columns_k(column,iband,myid)),2,j,iband))
!       enddo
!     enddo
!   enddo
  
!   !Heading DOWN
!   !Blocks 2 & 3 down to physlim_top
!   !Physical
!   do iband = 2,3 ! (6)
!     do j = jlim(2,grid,iband)-1,physlim_top,-1
!       do column = 1,columns_num(iband,myid)
!         !NS or FS
! x(iband)%f(j,column) =  real(x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+1,columns_k(column,iband,myid)),3,j,iband)*x(iband)%f(j+1,column))&
! &+im*aimag(x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+2,columns_k(column,iband,myid)),3,j,iband)*x(iband)%f(j+1,column))
!       enddo
!     enddo
!   enddo
  
!   xPL=0d0
!   !FFT physlim_top
!   call modes_to_planes_phys_lims(xPL,x,physlim_top,physlim_top,grid,myid,bandPL_FFT(myid),status,ierr)
!   do j=limPL_FFT(grid,1,myid),limPL_FFT(grid,2,myid)
!     if(j==physlim_top)then
!       call phys_to_four_N(xPL(1,1,j),bandPL_FFT(myid))
!     endif
!   enddo
!   call planes_to_modes_phys_lims(x,xPL,physlim_top,physlim_top,grid,myid,bandPL_FFT(myid),status,ierr)
  
!   !Block 3 down to end
!   !Fourier
!   iband = 3 ! (7)
!   do j = physlim_top-1,jlim(1,grid,iband),-1
!     do column = 1,columns_num(iband,myid)
!       !Taking NS
!       x(iband)%f(j,column) =  x(iband)%f(j,column)-ay(1,3,j,iband)*x(iband)%f(j+1,column)
!     enddo
!   enddo
  
!   !Block 2 downto physlim_bot
!   !Fourier
!   iband = 2 ! (8)
!   do j = physlim_top-1,physlim_bot,-1
!     do column = 1,columns_num(iband,myid)
!       x(iband)%f(j,column) =  x(iband)%f(j,column)-ay(1,3,j,iband)*x(iband)%f(j+1,column)
!     enddo
!   enddo
  
!   !Block 1 down to physlim_bot
!   !Fourier
!   iband = 1 ! (9)
!   do j = jlim(2,grid,iband)-1,physlim_bot,-1
!     do column = 1,columns_num(iband,myid)
!       x(iband)%f(j,column) =  x(iband)%f(j,column)-ay(1,3,j,iband)*x(iband)%f(j+1,column)
!     enddo
!   enddo
  
!   xPL=0d0
!   !FFT physlim_bot
!   call modes_to_planes_phys_lims(xPL,x,physlim_bot,physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)
!   do j=limPL_FFT(grid,1,myid),limPL_FFT(grid,2,myid)
!     if(j==physlim_bot)then
!       call four_to_phys_N(xPL(1,1,j),bandPL_FFT(myid))
!     endif
!   enddo
!   call planes_to_modes_phys_lims(x,xPL,physlim_bot,physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)
   
!   !Blocks 1 & 2 down to bottom
!   !Physical
!   do iband = 1,2 ! (10)
!     do j = physlim_bot-1,jlim(1,grid,iband),-1
!       do column = 1,columns_num(iband,myid)
!          !NS or FS
! x(iband)%f(j,column) =  real(x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+1,columns_k(column,iband,myid)),3,j,iband)*x(iband)%f(j+1,column))&
! +im*aimag(x(iband)%f(j,column)-&
! &ay(planeBC(2*columns_i(column,iband,myid)+2,columns_k(column,iband,myid)),3,j,iband)*x(iband)%f(j+1,column))
!       enddo
!     enddo
!   enddo

   
!    xPL=0d0
!    !FFT remaining planes to Fourier
!    call modes_to_planes_phys_lims(xPL,x,jlim(1,grid,2),physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)
!    call modes_to_planes_phys_lims(xPL,x,physlim_top+1,jlim(2,grid,2),grid,myid,bandPL_FFT(myid),status,ierr)
!    do j=limPL_FFT(grid,1,myid),limPL_FFT(grid,2,myid)
!      if(j<=physlim_bot)then
!        call phys_to_four_N(xPL(1,1,j),bandPL_FFT(myid))
!      endif
!      if(j>=physlim_top+1)then
!        call phys_to_four_N(xPL(1,1,j),bandPL_FFT(myid))
!      endif
!    enddo
!    call planes_to_modes_phys_lims(x,xPL,jlim(1,grid,2),physlim_bot,grid,myid,bandPL_FFT(myid),status,ierr)
!    call planes_to_modes_phys_lims(x,xPL,physlim_top+1,jlim(2,grid,2),grid,myid,bandPL_FFT(myid),status,ierr)

! ! -Original   
! !   do iband = sband,eband
! !     do column = 1,columns_num(iband,myid)
! !       x(iband)%f(jlim(1,grid,iband),column)=x(iband)%f(jlim(1,grid,iband),column)*ay(1,2,jlim(1,grid,iband),iband)
! !       do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
! !         x(iband)%f(j,column) = (x(iband)%f(j,column)-ay(1,1,j,iband)*x(iband)%f(j-1,column))*ay(1,2,j,iband)
! !       end do
! !       do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
! !         x(iband)%f(j,column) =  x(iband)%f(j,column)-ay(1,3,j,iband)*x(iband)%f(j+1,column)
! !       end do
! !     end do
! !   enddo

!   deallocate(xPL)
!   deallocate(axz)
!   deallocate(ay)

! end subroutine

! ! subroutine LUsolV(x,grid,ufield,myid)
! ! !----------------------------------------------------------------------*
! ! !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
! ! !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
! ! ! INPUT m_j number of unknowns
! ! !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
! ! !       a(2,:,i,k) inverse of the diagonal of L_j
! ! !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
! ! !       (note that the diagonal of U_j is 1)
! ! !       x(nystart:nyend,i,k)   rhs of the problem
! ! ! OUTPUT
! ! !       all untouched but the solution 'x'
! ! !----------------------------------------------------------------------*

! !   use declaration
! !   implicit none
  
! !   include 'mpif.h'             ! MPI variables
! !   integer status(MPI_STATUS_SIZE),ierr

! !   integer    :: i,k,j,column,iband,grid,ufield
! !   integer    :: iopt,myid
! !   type(cfield) x(sband:eband)
! !   real(8), allocatable:: xPL(:,:,:)
    
! !   real(8) , allocatable:: ay(:,:,:,:)
! !   real(8) , allocatable:: axz(:,:)

! !   !ay(BC,diag,j,iband)
! !   allocate(ay(2,3,jlim(1,grid,2):jlim(2,grid,2),nband))
! !   allocate(axz(maxval(columns_num),nband))

! !   do iband = sband,eband
! !     call LU_build(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,axz,ay,ufield)
! !     call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,ay(1,:,:,:))
! !     !call LU_dec(jlim(1,grid,iband),jlim(2,grid,iband),grid,myid,iband,ay(2,:,:,:))
! !   enddo
  
! !    do iband = sband,eband
! !     do column = 1,columns_num(iband,myid)
! !       do j=jlim(1,grid,iband)+1,jlim(2,grid,iband)-1 !C! Don't include first and last as BC's
! !         x(iband)%f(j,column)=x(iband)%f(j,column)/axz(column,iband)
        
! !       enddo
! !     enddo
! !   enddo
  
! ! ! -Original   
! !   do iband = sband,eband
! !     do column = 1,columns_num(iband,myid)
! !       x(iband)%f(jlim(1,grid,iband),column)=x(iband)%f(jlim(1,grid,iband),column)*ay(1,2,jlim(1,grid,iband),iband)
! !       do j = jlim(1,grid,iband)+1,jlim(2,grid,iband)
! !         x(iband)%f(j,column) = (x(iband)%f(j,column)-ay(1,1,j,iband)*x(iband)%f(j-1,column))*ay(1,2,j,iband)
! !       end do
! !       do j = jlim(2,grid,iband)-1,jlim(1,grid,iband),-1
! !         x(iband)%f(j,column) =  x(iband)%f(j,column)-ay(1,3,j,iband)*x(iband)%f(j+1,column)
! !       end do
! !     end do
! !   enddo

! !   deallocate(axz)
! !   deallocate(ay)

! ! end subroutine

! subroutine LUsolP(x,myid,iband,nystart,nyend)
! !----------------------------------------------------------------------*
! !      GIVEN a_j=L*U AND a_j*x_j=f_j FOR k=nstart:nend
! !      FIRST SOLVE L_j*y_j=f_j AND THEN U_j*x_j=y_j
! ! INPUT m_j number of unknowns
! !       a(1,nystart:nyend-1,i,k) lower diagonal of L_j
! !       a(2,:,i,k) inverse of the diagonal of L_j
! !       a(3,nystart+1:nyend,i,k) upper diagonal of U_j
! !       (note that the diagonal of U_j is 1)
! !       x(nystart:nyend,i,k)   rhs of the problem
! ! OUTPUT
! !       all untouched but the solution 'x'
! !----------------------------------------------------------------------*

!   use declaration
!   implicit none
!   integer i,k,j,column,myid,iband
!   integer nystart,nyend
!   complex(8) x(nystart:nyend,columns_num(iband,myid))
   
!   do column = 1,columns_num(iband,myid) 
!     i = columns_i(column,iband,myid)
!     k = columns_k(column,iband,myid)
!     x(nystart,column) = x(nystart,column)*DG(iband)%f_dg(2,nystart,column)
!     do j = nystart+1,nyend
!       x(j,column) = (x(j,column)-DG(iband)%f_dg(1,j,column)*x(j-1,column))*DG(iband)%f_dg(2,j,column)
!     end do
!     do j = nyend-1,nystart,-1
!       x(j,column) =  x(j,column)-DG(iband)%f_dg(3,j,column)*x(j+1,column)
!     end do

!   end do


! end subroutine

! subroutine LUsol0(x,nystart,nyend)
! !----------------------------------------------------------------------*
! !                  LUsol for real (0,0) mode
! !----------------------------------------------------------------------*

!   use declaration
!   implicit none
!   integer j
!   integer nystart,nyend
!   real(8) x(nystart:nyend)
!   real(8), allocatable:: a(:,:)
!   allocate(a(3,nystart:nyend))

!   call LU_build0(nystart,nyend,a)
!   call LU_dec0  (nystart,nyend,a)

!   x(nystart) = x(nystart)*a(2,nystart)
!   do j = nystart+1,nyend
!     x(j) = (x(j)-a(1,j)*x(j-1))*a(2,j)
!   end do
!   do j = nyend-1,nystart,-1
!     x(j) = x(j)-a(3,j)*x(j+1)
!   end do

!   deallocate(a)

! end subroutine


! subroutine LU_buildU(nystart,nyend,column,myid,iband,a)
!    !-------------------------------------------------------!
!    !       specifies original values of a(1:3,j,i,k)       !
!    !-------------------------------------------------------!

!    use declaration
!    implicit none
!    integer i,k,j,grid,myid,column,iband
!    integer nystart,nyend
!    real(8) k2x,k2z,beta
!    real(8) a(3,nystart:nyend)

!    beta = bRK(kRK)*dt/Re
  

!    !!!!!!!!!!!!!!!!!!!     u velocity:       !!!!!!!!!!!!!!!!!!!
!         i = columns_i(column,iband,myid)
!         k = columns_k(column,iband,myid)
!         k2x = k2F_x(i)
!         k2z = k2F_z(k)

!           a(1,nystart) = 0d0
!           a(2,nystart) = 1d0  
!           a(3,nystart) = gridweighting(iband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0

!         do j = nystart+1,nyend-1
!           a(1,j) =    -beta* dyu2i(1,j)
!           a(2,j) = 1d0-beta*(dyu2i(2,j)+k2x+k2z)
!           a(3,j) =    -beta* dyu2i(3,j)
!         end do

!           a(1,nyend) = gridweighting(iband,2)!0d0 !Free-shear -1 no-slip 0
!           a(2,nyend) = 1d0
!           a(3,nyend) = 0d0
! end subroutine


! subroutine LU_buildV(nystart,nyend,column,myid,iband,a)
!    !-------------------------------------------------------!
!    !       specifies original values of a(1:3,j,i,k)       !
!    !-------------------------------------------------------!

!    use declaration
!    implicit none
!    integer i,k,j,grid,myid,column,iband
!    integer nystart,nyend
!    real(8) k2x,k2z,beta
!    real(8) a(3,nystart:nyend)

!    beta = bRK(kRK)*dt/Re
  


!     !!!!!!!!!!!!!!!!!!!     v velocity:       !!!!!!!!!!!!!!!!!!!  
!         i = columns_i(column,iband,myid)
!         k = columns_k(column,iband,myid)
!         k2x = k2F_x(i)
!         k2z = k2F_z(k)

!           a(1,nystart) = 0d0
!           a(2,nystart) = 1d0  
!           a(3,nystart) = 0d0

!         do j = nystart+1,nyend-1
!           a(1,j) =    -beta* dyv2i(1,j)
!           a(2,j) = 1d0-beta*(dyv2i(2,j)+k2x+k2z)
!           a(3,j) =    -beta* dyv2i(3,j)
!         end do

!           a(1,nyend) = 0d0
!           a(2,nyend) = 1d0
!           a(3,nyend) = 0d0
! end subroutine

! subroutine LU_build(nystart,nyend,grid,myid,iband,axz,ay,ufield)
! !-------------------------------------------------------!
! !       specifies original values of a(1:3,j,i,k)       !
! !-------------------------------------------------------!

!   use declaration
!   implicit none
!   integer i,k,j,grid,myid,column,iband,ufield
!   integer nystart,nyend
!   real(8) k2x,k2z,beta
 
!   real(8) axz(maxval(columns_num),nband)
  
!   !ay(BC,diag,j,iband)
!   real(8) ay(2,3,jlim(1,grid,iband):jlim(2,grid,iband),nband)

!   beta = bRK(kRK)*dt/Re
  
!   !!!!!!!!!!!!!!!!!!!     u velocity:       !!!!!!!!!!!!!!!!!!!
!   if (grid==ugrid) then
!     if(ufield==1)then
!         ay(1,1,nystart,iband) = 0d0
!         ay(1,2,nystart,iband) = 1d0  
!         ay(1,3,nystart,iband) = -gridweighting_bc_u1 !gridweighting(iband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
!     elseif(ufield==3)then
!         ay(1,1,nystart,iband) = 0d0
!         ay(1,2,nystart,iband) = 1d0  
!         ay(1,3,nystart,iband) = -gridweighting_bc_u3 !gridweighting(iband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
!     endif
!       do j = nystart+1,nyend-1
!         ay(1,1,j      ,iband) =    -beta* dyu2i(1,j)
!         ay(1,2,j      ,iband) = 1d0-beta*(dyu2i(2,j))
!         ay(1,3,j      ,iband) =    -beta* dyu2i(3,j)
!       end do
!     if(ufield==1)then
!         ay(1,1,nyend  ,iband) = -gridweighting_bc_u1!gridweighting(iband,2)!0d0 !Free-shear -1 no-slip 0
!         ay(1,2,nyend  ,iband) = 1d0
!         ay(1,3,nyend  ,iband) = 0d0
!     elseif(ufield==3)then
! 	ay(1,1,nyend  ,iband) = -gridweighting_bc_u3!gridweighting(iband,2)!0d0 !Free-shear -1 no-slip 0
!         ay(1,2,nyend  ,iband) = 1d0
!         ay(1,3,nyend  ,iband) = 0d0
!     endif
!         ay(2,:,:,iband)=ay(1,:,:,iband)
!         ay(2,3,nystart,iband) = -1d0
!         ay(2,1,nyend  ,iband) = -1d0
    
!   !!!!!!!!!!!!!!!!!!!     v velocity:       !!!!!!!!!!!!!!!!!!!  
!   else if (grid==vgrid) then
!         ay(1,1,nystart,iband) = 0d0
!         ay(1,2,nystart,iband) = 1d0  
!         ay(1,3,nystart,iband) = 0d0
!       do j = nystart+1,nyend-1
!         ay(1,1,j      ,iband) =    -beta* dyv2i(1,j)
!         ay(1,2,j      ,iband) = 1d0-beta*(dyv2i(2,j))
!         ay(1,3,j      ,iband) =    -beta* dyv2i(3,j)
!       end do
!         ay(1,1,nyend  ,iband) = 0d0
!         ay(1,2,nyend  ,iband) = 1d0
!         ay(1,3,nyend  ,iband) = 0d0
        
!         ay(2,:,:,iband)=ay(1,:,:,iband)
    
!   end if
        
!    do column = 1,columns_num(iband,myid) 
!      i = columns_i(column,iband,myid)
!      k = columns_k(column,iband,myid)
!      k2x = k2F_x(i)
!      k2z = k2F_z(k)
!      axz(column,iband)=1-beta*(k2x+k2z)

!    enddo

! end subroutine
! !  
! subroutine LU_buildP(nystart,nyend,myid,iband,a)
! !-------------------------------------------------------!
! !       specifies original values of a(1:3,j,i,k)       !
! !-------------------------------------------------------!

!   use declaration
!   implicit none
!   type(rfield_dg)  a(sband:eband)
!   integer column,i,k,j,myid,iband
!   integer nystart,nyend
!   real(8) k2x,k2z

!   real(8) D_vec_b(nystart-1:nyend+1)
!   real(8) G_vec_b(nystart-1:nyend+1)
!   real(8) D_vec_t(nystart-1:nyend+1)
!   real(8) G_vec_t(nystart-1:nyend+1)

!   do j = nystart,nyend
!       D_vec_b(j-1) = -ddthetavi*dthdyu(j)
!       G_vec_b(j  ) = -ddthetavi*dthdyv(j)
!       D_vec_t(j  ) =  ddthetavi*dthdyu(j)
!       G_vec_t(j+1) =  ddthetavi*dthdyv(j)
!   end do
  
  
!   do column = 1,columns_num(iband,myid)
!     i = columns_i(column,iband,myid)
!     k = columns_k(column,iband,myid)
    
!     !For exact wavenumbers
!      k2x = k2F_x(i) 
!      k2z = k2F_z(k)
    
!     !For 2nd order centrered difference wavenumbers
! !    k2x = k1F_x(i)*k1F_x(i) 
! !    k2z = k1F_z(k)*k1F_z(k)
   
!     a(iband)%f_dg(2,nystart,column) = (D_vec_t(nystart)*G_vec_b(nystart  )) + k2x + k2z
!     a(iband)%f_dg(3,nystart,column) = (D_vec_t(nystart)*G_vec_t(nystart+1))
!     a(iband)%f_dg(2,nyend,  column) = (D_vec_b(nyend-1)*G_vec_t(nyend    )) + k2x + k2z
!     a(iband)%f_dg(1,nyend,  column) = (D_vec_b(nyend-1)*G_vec_b(nyend-1  ))
!     do j = nystart+1,nyend-1
!       a(iband)%f_dg(2,j,column) = (D_vec_b(j-1)*G_vec_t(j  )) + (D_vec_t(j)*G_vec_b(j)) + k2x + k2z
!       a(iband)%f_dg(1,j,column) = (D_vec_b(j-1)*G_vec_b(j-1))
!       a(iband)%f_dg(3,j,column) = (D_vec_t(j  )*G_vec_t(j+1))
!     end do

!   end do
  
!   if(myid==0)then
!     if(iband==midband)then 
!       a(iband)%f_dg(1,nystart,1) = 0d0
!       a(iband)%f_dg(2,nystart,1) = 1d0/((yu(1+1)-yu(1)))**2 !C!
!       a(iband)%f_dg(3,nystart,1) = 0d0
!     end if
!   end if

! ! For modified wavenumbers, need BC on last pressure mode (as k2x=k2z=0)
! !  do column = 1,columns_num(iband,myid)
! !  i = columns_i(column,iband,myid)
! !  k = columns_k(column,iband,myid)
! !  if(abs(k1F_x(i)*k1F_x(i))<10e-10.and.abs(k1F_z(k)*k1F_z(k))<10e-10)then
! !  a(iband)%f_dg(1,nystart,column) = 0d0
! !  a(iband)%f_dg(2,nystart,column) = 1d0/((yu(1+1)-yu(1)))**2 !C!
! !  a(iband)%f_dg(3,nystart,column) = 0d0
! !  endif
! !  enddo
  

!   call LU_decP(nystart,nyend,columns_num(iband,myid),a(iband)%f_dg(:,nystart:nyend,1:columns_num(iband,myid)),iband,myid)

  
! end subroutine

! subroutine LU_build0(nystart,nyend,a)
! !----------------------------------------------------------------------*
! !                  LU_build for real (0,0) mode
! !----------------------------------------------------------------------*

!   use declaration
!   implicit none
!   integer j
!   integer nystart,nyend
!   real(8) a(3,nystart:nyend)
!   real(8) beta

!   !!!!!!!!!!!!!!!!!!     velocity     !!!!!!!!!!!!!!!!!!
!   beta = bRK(kRK)*dt/Re
!   a(1,nystart) = 0d0
!   a(2,nystart) = 1d0
!   a(3,nystart) = gridweighting(midband,1)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
!   do j = nystart+1,nyend-1
!     a(1,j) =    -beta*dyu2i(1,j)
!     a(2,j) = 1d0-beta*dyu2i(2,j)
!     a(3,j) =    -beta*dyu2i(3,j)
!   end do
!   a(1,nyend) = gridweighting(midband,2)!0d0 !ay(3,nystart)=0d0 !Free-shear -1 no-slip 0
!   a(2,nyend) = 1d0
!   a(3,nyend) = 0d0

! end subroutine

! subroutine LU_dec(nystart,nyend,grid,myid,iband,ay)
! !----------------------------------------------------------------------*
! !         GIVEN 'a_j' FOR j=nystart:nyend
! !         PERFORMS a_j=L_j*U_j
! !         WITH 'a_j' allmost TRIDIAGONAL
! ! INPUT  nystart:nyend matrix leading dimension
! !        iopt selects original value of 'a'
! !        a_j is given as a tridiagonal with
! !        a(1,nystart:nyend-1) lower diagonal,
! !        a(2,:) main diagonal,
! !        a(3,nystart+1:nyend) upper diagonal
! !        the matrix 'a_j' is assembeld as:
! !  a_j(2,1) a_j(3,1)    0        0         0      .............       0
! !  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
! !     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
! !     0     ...................................................       0
! !     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
! !     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
! !
! ! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
! !         a(2,nystart:nyend) inverse of the diagonal of L_j
! !         a(3,nystart:nyend-1) upper diagonal of U_j
! !         (notice that the diagonal of U_j is 1)
! !----------------------------------------------------------------------*

!   use declaration
!   implicit none
!   integer j,column,iband,myid,grid
!   integer nystart,nyend
!   real(8) ay(3,jlim(1,grid,2):jlim(2,grid,2),nband)

!     ay(2,nystart,iband) = 1d0/ay(2,nystart,iband)
!     do j = nystart+1,nyend
!       ay(3,j-1,iband) = ay(3,j-1,iband)*ay(2,j-1,iband)
!       ay(2,j  ,iband) = 1d0/(ay(2,j,iband)-ay(1,j,iband)*ay(3,j-1,iband))
!     end do

  

! end subroutine
! !  
! subroutine LU_decP(nystart,nyend,columns,a,iband,myid)
! !----------------------------------------------------------------------*
! !         GIVEN 'a_j' FOR j=nystart:nyend
! !         PERFORMS a_j=L_j*U_j
! !         WITH 'a_j' allmost TRIDIAGONAL
! ! INPUT  nystart:nyend matrix leading dimension
! !        iopt selects original value of 'a'
! !        a_j is given as a tridiagonal with
! !        a(1,nystart:nyend-1) lower diagonal,
! !        a(2,:) main diagonal,
! !        a(3,nystart+1:nyend) upper diagonal
! !        the matrix 'a_j' is assembeld as:
! !  a_j(2,1) a_j(3,1)    0        0         0      .............       0
! !  a_j(1,2) a_j(2,2) a_j(3,2)    0         0      .............       0
! !     0     a_j(1,3) a_j(2,3) a_j(3,3)     0      .............       0
! !     0     ...................................................       0
! !     0     ...................... a_j(1,nyend-1) a_j(2,nyend-1)  a_j(3,nyend-1)
! !     0     ..........................     0      a_j(1,nyend  )  a_j(2,nyend  )
! !
! ! OUTPUT  a(1,nystart+1:nyend) lower diagonal of L_j
! !         a(2,nystart:nyend) inverse of the diagonal of L_j
! !         a(3,nystart:nyend-1) upper diagonal of U_j
! !         (notice that the diagonal of U_j is 1)
! !----------------------------------------------------------------------*

!   use declaration
!   implicit none

!   integer j,column
!   integer columns,iband,myid
!   integer nystart,nyend
!   real(8) a(3,nystart:nyend,columns)
  
!   do column = 1,columns
!     a(2,nystart,column) = 1d0/a(2,nystart,column)
!     do j = nystart+1,nyend
!       a(3,j-1,column) = a(3,j-1,column)*a(2,j-1,column)     
!       a(2,j  ,column) = 1d0/(a(2,j,column) - a(1,j,column)*a(3,j-1,column))
!     end do
!   end do 

! end subroutine

! subroutine LU_dec0(nystart,nyend,a)
! !----------------------------------------------------------------------*
! !                  LU_dec for real (0,0) mode
! !----------------------------------------------------------------------*

!   use declaration
!   implicit none
!   integer j
!   integer nystart,nyend
!   real(8) a(3,nystart:nyend)

!   a(2,nystart) = 1d0/a(2,nystart)
!   do j = nystart+1,nyend
!     a(3,j-1) = a(3,j-1)*a(2,j-1)
!     a(2,j  ) = 1d0/(a(2,j)-a(1,j)*a(3,j-1))
!   end do

! end subroutine
