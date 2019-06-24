        !COMPILER-GENERATED INTERFACE MODULE: Mon Jun 24 16:59:15 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE XTOF__genmod
          INTERFACE 
            SUBROUTINE XTOF(SZ,FS,CHI,F)
              INTEGER(KIND=4) :: FS
              INTEGER(KIND=4) :: SZ
              REAL(KIND=8) :: CHI(SZ,SZ+2)
              REAL(KIND=8) :: F(FS,SZ+2)
            END SUBROUTINE XTOF
          END INTERFACE 
        END MODULE XTOF__genmod
