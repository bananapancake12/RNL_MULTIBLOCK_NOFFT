        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  3 15:05:31 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RECVXC__genmod
          INTERFACE 
            SUBROUTINE RECVXC(XC,XCM,BUFFC,MX,MZ,GRID,STATUS,IERR)
              INTEGER(KIND=4) :: GRID
              INTEGER(KIND=4) :: MZ
              INTEGER(KIND=4) :: MX
              REAL(KIND=8) :: XC(MX,MZ,0:NY((GRID,1))-NY((GRID,0)))
              REAL(KIND=8) :: XCM(MX,MZ,LIMPL_INCW((GRID,1,0))-NY((GRID,&
     &0)):LIMPL_INCW((GRID,2,0))-NY((GRID,0)))
              REAL(KIND=8) :: BUFFC(MX,MZ,0:NY((GRID,1))-NY((GRID,0)))
              INTEGER(KIND=4) :: STATUS(5)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE RECVXC
          END INTERFACE 
        END MODULE RECVXC__genmod
