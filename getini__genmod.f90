        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  3 15:05:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETINI__genmod
          INTERFACE 
            SUBROUTINE GETINI(U1,U2,U3,P,DIV,MYID,STATUS,IERR)
              USE DECLARATION
              TYPE (CFIELD) :: U1(1:3)
              TYPE (CFIELD) :: U2(1:3)
              TYPE (CFIELD) :: U3(1:3)
              TYPE (CFIELD) :: P(1:3)
              TYPE (CFIELD) :: DIV(1:3)
              INTEGER(KIND=4) :: MYID
              INTEGER(KIND=4) :: STATUS(5)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE GETINI
          END INTERFACE 
        END MODULE GETINI__genmod
