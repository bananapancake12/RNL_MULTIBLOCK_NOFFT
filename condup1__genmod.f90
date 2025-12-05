        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  3 15:05:31 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CONDUP1__genmod
          INTERFACE 
            SUBROUTINE CONDUP1(YMC,Y2MC,X,XSIGN,GRID,MYID)
              USE DECLARATION
              INTEGER(KIND=4) :: MYID
              INTEGER(KIND=4) :: GRID
              REAL(KIND=8) :: YMC(DNX,DNZ,-LIMPL_INCW((GRID,2,MYID))+NY(&
     &(GRID,3))+1:-LIMPL_INCW((GRID,1,MYID))+NY((GRID,3))+1)
              REAL(KIND=8) :: Y2MC(DNX,DNZ,-LIMPL_INCW((GRID,2,MYID))+NY&
     &((GRID,3))+1:-LIMPL_INCW((GRID,1,MYID))+NY((GRID,3))+1)
              REAL(KIND=8) :: X(IGAL,KGAL,LIMPL_INCW((GRID,1,MYID)):    &
     &LIMPL_INCW((GRID,2,MYID)))
              REAL(KIND=8) :: XSIGN
            END SUBROUTINE CONDUP1
          END INTERFACE 
        END MODULE CONDUP1__genmod
