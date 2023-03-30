    PROGRAM small_strato_dr
    USE RK_f90_Integrator
    USE User_Parameters
    IMPLICIT NONE
    CHARACTER(LEN=12), PARAMETER, DIMENSION(7) :: SPC_NAMES = (/ &
      'O1D         ','O           ','O3          ', &
      'NO          ','NO2         ','M           ', &
      'O2          ' /)
    INTEGER, PARAMETER, DIMENSION(6) :: MONITOR = (/ &
        1,  2,  3,  4,  5,  7 /)
    CHARACTER(LEN=12), PARAMETER, DIMENSION(1) :: SMASS = (/ &
      'N           ' /)

    INTEGER :: i, ind_1 = 1, ind_2 = 2, NMONITOR=6, NMASS=1

    INTEGER, PARAMETER :: NSPEC = 7
!~~>  Control (in) and status (out) arguments for the integration
    DOUBLE PRECISION, DIMENSION(20) :: RCNTRL, RSTATUS
    INTEGER,       DIMENSION(20) :: ICNTRL, ISTATUS
    DOUBLE PRECISION :: CFACTOR, x
    DOUBLE PRECISION :: ATOL(NVAR), RTOL(NVAR), VAR(NVAR)
    DOUBLE PRECISION:: TSTART, TEND, T, l2,l2tmp,error
    DOUBLE PRECISION,DIMENSION(5):: ANS = (/154.2258448,1034529328.3727823,& 
832159397391.5885010,905729267.2488701,190770732.7509280/)
    INTEGER :: iter,hz,clock0,clock1,cpu
    DOUBLE PRECISION :: time
!~~~> Tolerances for calculating concentrations       
    DO i=1,NVAR
      RTOL(i) = 1.0d-2
      ATOL(i) = 1.0d-1
    END DO

    DO iter = 1, 5
!~~~> Tolerances for calculating concentrations       
    DO i=1,NVAR
      RTOL(i) = 0.1*RTOL(i)
      ATOL(i) = 0.1*ATOL(i)
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
    TEND = TSTART + (30*24*3600)

!~~~> Default control options
    ICNTRL(1:20) = 0
    RCNTRL(1:20) = 0.0d0

!~~~> Begin time loop

    T = TSTART

    call system_clock(count_rate=hz)
    call system_clock(count=clock0)
    CALL INTEGRATE(TIN=T,TOUT=TEND,N=NVAR,NNZERO=NNZ,VAR=VAR,RTOL=RTOL,&
     ATOL=ATOL,FUN=FUN_U,JAC=JAC_U,RSTATUS_U=RSTATUS,ISTATUS_U=ISTATUS,&
         ICNTRL_U=ICNTRL)
    call system_clock(count=clock1)
    cpu = clock1-clock0
    time = real(cpu)/(real(hz))
!~~~> End time loop ~~~~~~~~~~
!    print *,VAR
    l2=0d0
    l2tmp=0d0
    do i=1,NVAR
         l2=l2+(VAR(i)-ANS(i))*(VAR(i)-ANS(i))
         l2tmp =l2tmp+ANS(i)*ANS(i)
    end do
    error=sqrt(l2/l2tmp)
    write(6,919) error,time
919 FORMAT(2x,'L2 norm:',F12.9,'     time elapsed:',F12.9)
 

     WRITE(6,991) RSTATUS(1),     &
               ( TRIM(SPC_NAMES(MONITOR(i))),           &
                 VAR(MONITOR(i))/CFACTOR, i=1,NVAR )

991   FORMAT(' AT TIME',E10.3,2X,200(A,'=',E11.4,'; '))
    end do
    END PROGRAM small_strato_dr


