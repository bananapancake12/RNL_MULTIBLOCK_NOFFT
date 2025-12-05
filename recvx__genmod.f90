        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  3 15:05:31 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RECVX__genmod
          INTERFACE 
            SUBROUTINE RECVX(X,XM,GRID,STATUS,IERR)
              INTEGER(KIND=4) :: GRID
              REAL(KIND=8) :: X(NY((GRID,0)):NY((GRID,3))+1)
              REAL(KIND=8) :: XM(LIMPL_INCW((GRID,1,0)):LIMPL_INCW((GRID&
     &,2,0)))
              INTEGER(KIND=4) :: STATUS(5)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE RECVX
          END INTERFACE 
        END MODULE RECVX__genmod
