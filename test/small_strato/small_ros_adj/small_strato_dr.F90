    PROGRAM small_strato_dr
    USE ROS_ADJ_f90_Integrator
    USE User_Parameters
    IMPLICIT NONE
    INTEGER,PARAMETER :: NADJ = 5
    CHARACTER(LEN=12), PARAMETER, DIMENSION(7) :: SPC_NAMES = (/ &
      'O1D         ','O           ','O3          ', &
      'NO          ','NO2         ','M           ', &
      'O2          ' /)
    INTEGER, PARAMETER, DIMENSION(6) :: MONITOR = (/ &
        1,  2,  3,  4,  5,  7 /)
    CHARACTER(LEN=12), PARAMETER, DIMENSION(1) :: SMASS = (/ &
      'N           ' /)

    INTEGER :: i, j, ind_1 = 1, ind_2 = 2, NMONITOR=6, NMASS=1, NP=0

    INTEGER, PARAMETER :: NSPEC = 7
    DOUBLE PRECISION, DIMENSION(NVAR,NADJ) :: Y_adj, ATOL_adj, RTOL_adj

!~~>  Control (in) and status (out) arguments for the integration
    DOUBLE PRECISION, DIMENSION(20) :: RCNTRL, RSTATUS
    INTEGER,       DIMENSION(20) :: ICNTRL, ISTATUS
    DOUBLE PRECISION :: CFACTOR, x
    DOUBLE PRECISION :: ATOL(NVAR), RTOL(NVAR), VAR(NVAR), &
                       DVAL(NSPEC)

    DOUBLE PRECISION:: TSTART, TEND, T

!~~~> Tolerances for calculating concentrations       
    DO i=1,NVAR
      RTOL(i) = 1.0d-4
      ATOL(i) = 1.0d-3
    END DO

    DO j=1,NADJ
      DO i=1,NVAR
        RTOL_adj(i,j) = 1.0d-4
        ATOL_adj(i,j) = 1.0d-10
      END DO
    END DO

!      Initialize
    CFACTOR = 1.000000e+00
    x = (0.)*CFACTOR
    DO i = 1, NVAR
      VAR(i) = x
    END DO
    x = (0.)*CFACTOR
    DO i = 1, NFIX
      FIX(i) = x
    END DO
    VAR(1) = (9.906E+01)*CFACTOR
    VAR(2) = (6.624E+08)*CFACTOR
    VAR(3) = (5.326E+11)*CFACTOR
    VAR(4) = (8.725E+08)*CFACTOR
    VAR(5) = (2.240E+08)*CFACTOR
    FIX(1) = (8.120E+16)*CFACTOR
    FIX(2) = (1.697E+16)*CFACTOR
    TSTART = (12*3600)
    TEND = TSTART + (3*24*3600)

!~~~> Default control options
    ICNTRL(1:20) = 0
    RCNTRL(1:20) = 0.0d0

!~~~> Begin time loop

    T = TSTART
!  get mass
    DVAL(1) = VAR(4)+VAR(5)+FIX(1)+FIX(1)


    CALL INTEGRATE_ADJ(TIN=T,TOUT=TEND,NVAR=NVAR,NNZERO=NNZ,Lambda=Y_adj,Y=VAR, &
       NP=NP,RTOL=RTOL,ATOL=ATOL,NADJ=NADJ,ATOL_adj=ATOL_adj,RTOL_adj=RTOL_adj,&
FUN=FUN_U, JAC=JAC_U,ADJINIT=adjinit,HESSTR_VEC= HESSTR_VEC_U, &
      RSTATUS_U=RSTATUS, ISTATUS_U=ISTATUS,ICNTRL_U=ICNTRL)
    DVAL(1) = VAR(4)+VAR(5)+FIX(1)+FIX(1)


!~~~> End time loop ~~~~~~~~~~
    print *,VAR
   
    OPEN(20, FILE='rkadj_ADJ_results.m')
    WRITE(6,*) '**************************************************'
    WRITE(6,*) ' Concentrations and Sensitivities at final time '
    WRITE(6,*) ' were written in the file rkadj_ADJ_results.m'
    WRITE(6,*) '**************************************************'
    DO j=1,NADJ
      WRITE(20,993) ( Y_adj(i,j), i=1,NVAR )
    END DO

    WRITE(6,995) TRIM(SPC_NAMES(ind_1)),TRIM(SPC_NAMES(ind_1)), &
                   Y_adj(ind_1,1)
    WRITE(6,995) TRIM(SPC_NAMES(ind_2)),TRIM(SPC_NAMES(ind_2)), &
                   Y_adj(ind_2,2)
    WRITE(6,995) TRIM(SPC_NAMES(ind_2)),TRIM(SPC_NAMES(ind_1)), &
                   Y_adj(ind_1,2)
    WRITE(6,995) TRIM(SPC_NAMES(ind_1)),TRIM(SPC_NAMES(ind_2)), &
                   Y_adj(ind_2,1)


993 FORMAT(1000(E24.16,2X))
995 FORMAT('ADJ: d[',A,'](tf)/d[',A,'](t0)=',E14.7)

      !~~~> The entire matrix of sensitivities
    WRITE(6,996) ( 'd ',TRIM(SPC_NAMES(i)), i=1,NVAR ) 
    DO j=1,NADJ
      WRITE(6,997) TRIM(SPC_NAMES(j)),( Y_adj(i,j), i=1,NVAR )
    END DO
996 FORMAT(12X,100('  ',A2,A6,4X))
997 FORMAT('d/d',A6,' = ',100(E12.5,2X))

    END PROGRAM small_strato_dr


