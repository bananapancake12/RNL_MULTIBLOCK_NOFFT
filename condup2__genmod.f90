        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  3 15:05:31 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CONDUP2__genmod
          INTERFACE 
            SUBROUTINE CONDUP2(YMC,X1,X2,XSIGN1,XSIGN2,GRID,MYID)
              USE DECLARATION
              INTEGER(KIND=4) :: MYID
              INTEGER(KIND=4) :: GRID
              REAL(KIND=8) :: YMC(DNX,DNZ,-LIMPL_INCW((GRID,2,MYID))+NY(&
     &(GRID,3))+1:-LIMPL_INCW((GRID,1,MYID))+NY((GRID,3))+1)
              REAL(KIND=8) :: X1(IGAL,KGAL,LIMPL_INCW((GRID,1,MYID)):   &
     &LIMPL_INCW((GRID,2,MYID)))
              REAL(KIND=8) :: X2(IGAL,KGAL,LIMPL_INCW((GRID,1,MYID)):   &
     &LIMPL_INCW((GRID,2,MYID)))
              REAL(KIND=8) :: XSIGN1
              REAL(KIND=8) :: XSIGN2
            END SUBROUTINE CONDUP2
          END INTERFACE 
        END MODULE CONDUP2__genmod
