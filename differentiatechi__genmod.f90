        !COMPILER-GENERATED INTERFACE MODULE: Mon Jun 24 16:59:15 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIFFERENTIATECHI__genmod
          INTERFACE 
            SUBROUTINE DIFFERENTIATECHI(HOLDER,NODES,NODE,I,SZ,XSUM,    &
     &WEIGHTS)
              INTEGER(KIND=4) :: SZ
              REAL(KIND=8) :: HOLDER(SZ)
              REAL(KIND=8) :: NODES(SZ)
              REAL(KIND=8) :: NODE
              INTEGER(KIND=4) :: I
              REAL(KIND=8) :: XSUM
              REAL(KIND=8) :: WEIGHTS(SZ)
            END SUBROUTINE DIFFERENTIATECHI
          END INTERFACE 
        END MODULE DIFFERENTIATECHI__genmod
