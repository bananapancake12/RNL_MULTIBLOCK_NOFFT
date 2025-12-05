        !COMPILER-GENERATED INTERFACE MODULE: Wed Dec  3 15:05:30 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MBLOCK_INI__genmod
          INTERFACE 
            SUBROUTINE MBLOCK_INI(U1,U2,U3,P,MYID,STATUS,IERR)
              USE DECLARATION
              TYPE (CFIELD) :: U1(1:3)
              TYPE (CFIELD) :: U2(1:3)
              TYPE (CFIELD) :: U3(1:3)
              TYPE (CFIELD) :: P(1:3)
              INTEGER(KIND=4) :: MYID
              INTEGER(KIND=4) :: STATUS(5)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE MBLOCK_INI
          END INTERFACE 
        END MODULE MBLOCK_INI__genmod
