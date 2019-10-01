! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! INTEGRATE - Integrator routine
!   Arguments :
!      TIN       - Start Time for Integration
!      TOUT      - End Time for Integration
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  ERK - Explicit Runge-Kutta methods                                     !
!            * RK2(3): 3 stages, order 2                                  !
!            * RK3(2): 4 stages, order 3                                  !
!            * RK4(3): 5 stages, order 4                                  !
!            * Dopri5: 7 stages, order 5                                  !
!            * Verner: 8 stages, order 6                                  !
!            * Dopri853: 12 stages, order 8                               !
!                                                                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE ERK_f90_Integrator
     IMPLICIT NONE
     ! Arguments 
      INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4,  &
               Nrej=5, Ndec=6, Nsol=7, Nsng=8,               &
               Ntexit=1, Nhexit=2, Nhnew=3
CONTAINS
SUBROUTINE INTEGRATE( TIN, TOUT, NVAR, VAR, RTOL, ATOL, FUN, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U )

   IMPLICIT NONE
!~~~>  Statistics on the work performed by the ERK method
   INTEGER, INTENT(IN) :: NVAR
   DOUBLE PRECISION, INTENT(INOUT) :: VAR(NVAR)
   DOUBLE PRECISION, INTENT(IN) :: RTOL(NVAR)
   DOUBLE PRECISION, INTENT(IN) :: ATOL(NVAR)
   DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: Ierr_U

!   INTEGER, SAVE :: Ntotal = 0
   DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
   INTEGER       :: ICNTRL(20), ISTATUS(20), Ierr

   EXTERNAL :: FUN
   
   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0


   !~~~> fine-tune the integrator:
   ICNTRL(2) = 0	! 0 - vector tolerances, 1 - scalar tolerances
   ICNTRL(6) = 0	! starting values of Newton iterations: interpolated (0), zero (1)
  ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   
   T1 = TIN; T2 = TOUT
   CALL ERK( NVAR,T1,T2,VAR,RTOL,ATOL,FUN,  &
               RCNTRL,ICNTRL,RSTATUS,ISTATUS,Ierr )

   !~~~> Debug option: print number of steps
!   Ntotal = Ntotal + ISTATUS(Nstp)
!   PRINT*,'NSTEPS=',ISTATUS(Nstp), '(',Ntotal,')'
   

   IF (Ierr < 0) THEN
        PRINT *,'ERK: Unsuccessful exit at T=',TIN,' (Ierr=',Ierr,')'
   ENDIF
  
   !for Debug
!   WRITE(UNIT=6,FMT=3844) TOUT,VAR(1:2),ISTATUS(3:5)
!   3844 FORMAT('Final T =', f5.2,' Value=', 2E12.5,' Nsteps=',I6, &
!    ' NAsteps=',I6,' NRsteps=', I6)
 
   ! if optional parameters are given for output they to return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(Ierr_U))    Ierr_U       = Ierr

   END SUBROUTINE INTEGRATE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERK(NVAR, Tinitial, Tfinal, Y, RelTol, AbsTol, FUN, &
                       RCNTRL, ICNTRL, RSTATUS, ISTATUS, Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using an Explicit
!    Runge-Kutta (ERK) method.
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!-     NVAR       = dimension of input system 
!-     Y(NVAR)    = vector of initial conditions (at T=Tinitial)
!-    [Tinitial,Tfinal]  = time range of integration
!     (if Tinitial>Tfinal the integration is performed backwards in time)
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE Fun( NVAR, T, Y, P ) = ODE function,
!                       returns Ydot = Y' = F(T,Y)
!-    ICNTRL(1:20)    = integer inputs parameters
!-    RCNTRL(1:20)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!-    Y(NVAR)         -> vector of final states (at T->Tfinal)
!-    ISTATUS(1:20)   -> integer output parameters
!-    RSTATUS(1:20)   -> real output parameters
!-    Ierr            -> job status upon return
!                        success (positive value) or
!                        failure (negative value)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    Note: For input parameters equal to zero the default values of the
!          corresponding variables are used.
!~~~>  
!    ICNTRL(1) = not used
!
!    ICNTRL(2) = 0: AbsTol, RelTol are NVAR-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3) = Method
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0 the default value of 10000 is used
!
!
!~~~>  Real parameters
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!                  It is strongly recommended to keep Hmin = ZERO
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                 (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly generaler
!                  than the predicted value  (default=0.9)
!    RCNTRL(10) -> Qmin
!    RCNTRL(11) -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
!                  step size is kept constant and the LU factorization
!                  reused (default Qmin=1, Qmax=1.2)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call adds the current no. of fcn calls
!      to previous value of ISTATUS(1), and similar for the other params.
!      Set ISTATUS(1:10) = 0 before call to avoid this accumulation.
!
!    ISTATUS(1) = No. of function calls
!    ISTATUS(2) = No. of jacobian calls
!    ISTATUS(3) = No. of steps
!    ISTATUS(4) = No. of accepted steps
!    ISTATUS(5) = No. of rejected steps (except at the beginning)
!    ISTATUS(6) = No. of LU decompositions
!    ISTATUS(7) = No. of forward/backward substitutions
!    ISTATUS(8) = No. of singular matrix decompositions
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                   computed Y upon return
!    RSTATUS(2)  -> Hexit,last accepted step before return
!    RSTATUS(3)  -> Hnew, last predicted step before return
!        For multiple restarts, use Hnew as Hstart in the following run
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      INTEGER, INTENT(IN)          :: NVAR, ICNTRL(20)
      DOUBLE PRECISION, INTENT(IN)    :: Tinitial, Tfinal, &
                    RelTol(NVAR), AbsTol(NVAR), RCNTRL(20)
      DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)
      INTEGER, INTENT(OUT)         :: Ierr
      INTEGER, INTENT(INOUT)       :: ISTATUS(20) 
      DOUBLE PRECISION, INTENT(OUT)   :: RSTATUS(20)
       
!~~~>  ERK method coefficients, up to 12 stages
      INTEGER, PARAMETER :: Smax = 12
      INTEGER, PARAMETER :: RK2=1, RK3=2, RK4=3, RK5=4, RK6=5, RK8=6
      DOUBLE PRECISION :: rkA(Smax,Smax), rkB(Smax), rkC(Smax), &
                       rkE(Smax),rkELO 
      INTEGER :: erkMethod, rkS ! The number of stages
! Local variables      
      DOUBLE PRECISION :: Hmin, Hmax, Hstart, Roundoff,    &
                       FacMin, Facmax, FacSafe, FacRej, &
                       Qmin, Qmax
      INTEGER       :: ITOL, Max_no_steps, i
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
      DOUBLE PRECISION DLAMCH
      EXTERNAL FUN

      rkA = ZERO
      rkB = ZERO
      rkC = ZERO
      rkE = ZERO
      ISTATUS(1:20) = 0
      RSTATUS(1:20) = ZERO
      Ierr = 0

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:NVAR) and RelTol(1:NVAR)
      IF (ICNTRL(2) == 0) THEN
         ITOL = 1
      ELSE
         ITOL = 0
      END IF

!~~~> ICNTRL(3) - method selection       
      SELECT CASE (ICNTRL(3))
      CASE (1)
         CALL Erk23
      CASE (2)
         CALL Erk3_Heun
      CASE (3)
         CALL Erk43
      CASE (0,4)
         CALL Dopri5
      CASE (5)
         CALL Verner
      CASE (6)
         CALL Dopri853
      CASE DEFAULT
         CALL Dopri5
      END SELECT
      
!~~~>   The maximum number of time steps admitted
      IF (ICNTRL(4) == 0) THEN
         Max_no_steps = 10000
      ELSEIF (ICNTRL(4) > 0) THEN
         Max_no_steps=ICNTRL(4)
      ELSE
         PRINT * ,'User-selected ICNTRL(4)=',ICNTRL(4)
         CALL ERK_ErrorMsg(-1,Tinitial,ZERO,Ierr)
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
      Roundoff = DLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
      IF (RCNTRL(1) == ZERO) THEN
         Hmin = ZERO
      ELSEIF (RCNTRL(1) > ZERO) THEN
         Hmin = RCNTRL(1)
      ELSE
         PRINT * , 'User-selected RCNTRL(1)=', RCNTRL(1)
         CALL ERK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Upper bound on the step size: (positive value)
      IF (RCNTRL(2) == ZERO) THEN
         Hmax = ABS(Tfinal-Tinitial)
      ELSEIF (RCNTRL(2) > ZERO) THEN
         Hmax = MIN(ABS(RCNTRL(2)),ABS(Tfinal-Tinitial))
      ELSE
         PRINT * , 'User-selected RCNTRL(2)=', RCNTRL(2)
         CALL ERK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Starting step size: (positive value)
      IF (RCNTRL(3) == ZERO) THEN
         Hstart = MAX(Hmin,Roundoff)
      ELSEIF (RCNTRL(3) > ZERO) THEN
         Hstart = MIN(ABS(RCNTRL(3)),ABS(Tfinal-Tinitial))
      ELSE
         PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
         CALL ERK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
  
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax
      IF (RCNTRL(4) == ZERO) THEN
         FacMin = 0.2
      ELSEIF (RCNTRL(4) > ZERO) THEN
         FacMin = RCNTRL(4)
      ELSE
         PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
         CALL ERK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
      IF (RCNTRL(5) == ZERO) THEN
         FacMax = 10.0
      ELSEIF (RCNTRL(5) > ZERO) THEN
         FacMax = RCNTRL(5)
      ELSE
         PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
         CALL ERK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
      IF (RCNTRL(6) == ZERO) THEN
         FacRej = 0.1
      ELSEIF (RCNTRL(6) > ZERO) THEN
         FacRej = RCNTRL(6)
      ELSE
         PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
         CALL ERK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
      IF (RCNTRL(7) == ZERO) THEN
         FacSafe = 0.9
      ELSEIF (RCNTRL(7) > ZERO) THEN
         FacSafe = RCNTRL(7)
      ELSE
         PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
         CALL ERK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF

!~~~> Qmin, Qmax: IF Qmin < Hnew/Hold < Qmax, STEP SIZE = CONST.
      IF(RCNTRL(10) == ZERO)THEN
         Qmin=ONE
      ELSE
         Qmin=RCNTRL(10)
      END IF
      IF(RCNTRL(11) == ZERO)THEN
         Qmax=1.2D0
      ELSE
         Qmax=RCNTRL(11)
      END IF

!~~~>  Check if tolerances are reasonable
      IF (ITOL == 0) THEN
         IF (AbsTol(1) <= ZERO .OR. RelTol(1) <= 10.D0*Roundoff) THEN
            PRINT * , ' Scalar AbsTol = ',AbsTol(1)
            PRINT * , ' Scalar RelTol = ',RelTol(1)
            CALL ERK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
         END IF
      ELSE
         DO i=1,NVAR
            IF (AbsTol(i) <= 0.D0.OR.RelTol(i) <= 10.D0*Roundoff) THEN
              PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
              PRINT * , ' RelTol(',i,') = ',RelTol(i)
              CALL ERK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
            END IF
         END DO
      END IF
    
    IF (Ierr < 0) RETURN

    CALL ERK_Integrator( NVAR,Tinitial,Tfinal,Y,FUN,Ierr )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CONTAINS  !  PROCEDURES INTERNAL TO ERK
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERK_Integrator( NVAR,Tinitial,Tfinal,Y,FUN,Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      IMPLICIT NONE      
      INTEGER, INTENT(IN) :: NVAR
      DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)
      DOUBLE PRECISION, INTENT(IN) :: Tinitial, Tfinal
      INTEGER, INTENT(OUT) :: Ierr
      
!~~~> Local variables:      
      DOUBLE PRECISION :: K(NVAR,rkS), TMP(NVAR), SCAL(NVAR), &
                       T, H, Hratio, Err, Fac, Hnew,Tdirection      
      INTEGER :: j, istage  
      LOGICAL :: Reject, FirstStep
          
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
      EXTERNAL  FUN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initializations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      T = Tinitial
      Tdirection = SIGN(ONE,Tfinal-Tinitial)
      H = MAX(ABS(Hmin),ABS(Hstart))
      IF (ABS(H) <= 10.D0*Roundoff) H=1.0D-6
      H=MIN(ABS(H),Hmax)
      H=SIGN(H,Tdirection)
      Reject=.FALSE.
      FirstStep=.TRUE.

      CALL ERK_ErrorScale(NVAR, ITOL, AbsTol, RelTol, Y, SCAL)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( (Tfinal-T)*Tdirection - Roundoff > ZERO )

      IF (ISTATUS(Nstp) > Max_no_steps) THEN
             CALL ERK_ErrorMsg(-6,T,H,Ierr); RETURN
      END IF   
      IF ( (T+0.1d0*H == T) .OR. (ABS(H) <= Roundoff) ) THEN
             CALL ERK_ErrorMsg(-7,T,H,Ierr); RETURN
      END IF   

stages:DO istage = 1, rkS

       CALL Set2zero(NVAR,K(1,istage))
   
       CALL Set2zero(NVAR,TMP)

       IF(istage > 1) THEN
           DO j = 1, istage-1
               CALL DAXPY(NVAR,H*rkA(istage,j),K(1,j),1,TMP,1)
           END DO
       END IF
        
       CALL DAXPY(NVAR,ONE,Y,1,TMP,1)

       CALL FUN(NVAR,T+rkC(istage)*H,TMP,K(1,istage))

       ISTATUS(Nfun) = ISTATUS(Nfun) + 1
 
  END DO stages

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Error estimation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      CALL Set2zero(NVAR,TMP)
      DO i = 1,rkS
         IF (rkE(i)/=ZERO) CALL DAXPY(NVAR,H*rkE(i),K(1,i),1,TMP,1)
      END DO  

      CALL ERK_ErrorNorm(NVAR, TMP, SCAL, Err)

!~~~> Computation of new step size Hnew
      Fac  = FacSafe*(Err)**(-ONE/rkELO)
      Fac  = MAX(FacMin,MIN(FacMax,Fac))
      Hnew = H*Fac

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Accept/Reject step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
accept: IF ( Err < ONE ) THEN !~~~> Step is accepted

         FirstStep=.FALSE.
         ISTATUS(Nacc) = ISTATUS(Nacc) + 1

!~~~> Update time and solution
         
!for debug
!         print *,'Accepted step. Time= ',T,'; Stepsize= ',H         
         
         T  =  T + H
         ! Y(:) <-- Y(:) + Sum_j rkE(j)*K_j(:)
         DO i = 1,rkS 
            IF (rkB(i)/=ZERO) CALL DAXPY(NVAR,H*rkB(i),K(1,i),1,Y,1)
         END DO  
       
!~~~> Update scaling coefficients
         CALL ERK_ErrorScale(NVAR, ITOL, AbsTol, RelTol, Y, SCAL)

!~~~> Next time step
         Hnew = Tdirection*MIN(ABS(Hnew),Hmax)
         ! Last T and H
         RSTATUS(Ntexit) = T
         RSTATUS(Nhexit) = H
         RSTATUS(Nhnew)  = Hnew
         ! No step increase after a rejection
         IF (Reject) Hnew = Tdirection*MIN(ABS(Hnew),ABS(H))
         Reject = .FALSE.
         IF ((T+Hnew/Qmin-Tfinal)*Tdirection > ZERO) THEN
            H = Tfinal-T
         ELSE
            Hratio=Hnew/H
            H = Hnew
         END IF
	 
      ELSE accept !~~~> Step is rejected

         IF (FirstStep .OR. Reject) THEN
             H = FacRej*H
         ELSE
             H = Hnew
         END IF
         Reject  = .TRUE.
         IF (ISTATUS(Nacc) >= 1) ISTATUS(Nrej) = ISTATUS(Nrej) + 1

         !for debug
!         print *,'Rejected step. Time= ',T,'; Stepsize= ',H

      END IF accept
      
      END DO Tloop

      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE ERK_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_ErrorScale(NVAR, ITOL, AbsTol, RelTol, Y, SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER :: i, NVAR, ITOL
      DOUBLE PRECISION :: AbsTol(NVAR), RelTol(NVAR), &
                       Y(NVAR), SCAL(NVAR)
      IF (ITOL == 0) THEN
        DO i=1,NVAR
          SCAL(i) = ONE / ( AbsTol(1)+RelTol(1)*ABS(Y(i)) )
        END DO
      ELSE
        DO i=1,NVAR
          SCAL(i) = ONE / ( AbsTol(i)+RelTol(i)*ABS(Y(i)) )
        END DO
      END IF
      END SUBROUTINE ERK_ErrorScale
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_ErrorNorm(NVAR, Y, SCAL, Err)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      
      INTEGER :: i, NVAR
      DOUBLE PRECISION :: Y(NVAR), SCAL(NVAR), Err      
      Err = ZERO
      DO i=1,NVAR
           Err = Err+(Y(i)*SCAL(i))**2
      END DO
      Err = MAX( SQRT(Err/DBLE(NVAR)), 1.0d-10 )
!
      END SUBROUTINE ERK_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_ErrorMsg(Code,T,H,Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DOUBLE PRECISION, INTENT(IN) :: T, H
      INTEGER, INTENT(IN)  :: Code
      INTEGER, INTENT(OUT) :: Ierr

      Ierr = Code
      PRINT * , &
        'Forced exit from ERK due to the following error:'

      SELECT CASE (Code)
        CASE (-1)
          PRINT * , '--> Improper value for maximal no of steps'
        CASE (-2)
          PRINT * , '--> Improper value for maximal no of Newton iterations'
        CASE (-3)
          PRINT * , '--> Hmin/Hmax/Hstart must be positive'
        CASE (-4)
          PRINT * , '--> FacMin/FacMax/FacRej must be positive'
        CASE (-5)
          PRINT * , '--> Improper tolerance values'
        CASE (-6)
          PRINT * , '--> No of steps exceeds maximum bound'
        CASE (-7)
          PRINT * , '--> Step size too general: T + 10*H = T', &
                    ' or H < Roundoff'
        CASE (-8)
          PRINT * , '--> Matrix is repeatedly singular'
        CASE DEFAULT
          PRINT *, 'Unknown Error code: ', Code
      END SELECT

      PRINT *, "T=", T, "and H=", H

      END SUBROUTINE ERK_ErrorMsg
      
      SUBROUTINE Erk23
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initializes coefficients (the Butcher tableau) for the RK2(3)
! Runge-Kutta method. This is a 2(3) Fehlberg-type method: second
! order accurate with a third order error control mechanism
!
! A = [ 0 0 0;
!       1 0 0;
!       0.25 0.25 0 ]
!
! b = [ 0.5 0.5 0 ]
!
! c = [ 0 1 0.5 ]
!
! e = [ 1/3 1/3 -2/3 ]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE

      rkS = 3
      rkELO = 4.d0
      erkMethod = RK2

      rkC(2) = 1.d0
      rkC(3) = 0.5d0

      rkA(2,1) = 1.d0
      rkA(3,1) = 0.25d0
      rkA(3,2) = 0.25d0

      rkB(1) = 0.5d0
      rkB(2) = 0.5d0

      rkE(1) = 1.d0/3.d0
      rkE(2) = 1.d0/3.d0
      rkE(3) = -2.d0/3.d0

      END SUBROUTINE Erk23

      SUBROUTINE Erk3_Heun
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initializes parameters (Butcher tableau) for the RK3(2)
! Runge-Kutta method built on top of Heun's third order method.
! It incorporates a second-order accurate error control mechanism.
!
! A = [ 0 0 0 0;
!       1/3 0 0 0;
!       0 2/3 0 0;
!       0.25 0 0.75 0]
!
! b = [ 0.25 0 0.75 0 ]
!
! c = [ 0 1/3 2/3 0 ]
!
! e = [ 0.25 -0.5 0.25 0 ]
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IMPLICIT NONE

      rkS = 4
      rkELO = 4.d0
      erkMethod = RK3

      rkC(2) = 1.d0 / 3.d0
      rkC(3) = 2.d0 / 3.d0

      rkA(2,1) = 1.d0 / 3.d0
      rkA(3,2) = 2.d0 / 3.d0
      rkA(4,1) = 0.25d0
      rkA(4,3) = 0.75d0

      rkB(1) = 0.25d0
      rkB(3) = 0.75d0

      rkE(1) = 0.25d0
      rkE(2) = -0.5d0
      rkE(3) = 0.25d0

      END SUBROUTINE Erk3_Heun

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Erk43
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Initializes coefficients (the Butcher tableau) for the RK4(3)
! Runge-Kutta method. This method is built by appending an
! additional stage to the well-known "3/8-rule", thus incorporating
! a 3rd order error control mechanism.
!
! A = [ 0 0 0 0 0;
!       1/3 0 0 0 0;
!       -1/3 1 0 0 0;
!       1 -1 -1 0 0;
!       1/8 3/8 3/8 1/8 0];
!
! b = [ 1/8 3/8 3/8 1/8 0]
!
! c = [ 0 1/3 2/3 1 1 ]
!
! e = [ 1/24 -1/8 1/8 1/8 -1/6 ]
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      erkMethod = RK4
      rkELO = 5.d0
      rkS = 5

      rkA(2,1) = 1.d0/3.d0
      rkA(3,1) = -1.d0/3.d0
      rkA(3,2) = 1.d0
      rkA(4,1) = 1.d0
      rkA(4,2) =-1.d0
      rkA(4,3) = 1.d0
      rkA(5,1) = 1.d0/8.d0
      rkA(5,2) = 3.d0/8.d0
      rkA(5,3) = 3.d0/8.d0
      rkA(5,4) = 1.d0/8.d0

      rkB(1)   =  1.d0/8.d0
      rkB(2)   =  3.d0/8.d0
      rkB(3)   =  3.d0/8.d0
      rkB(4)   =  1.d0/8.d0

      rkC(2)   = 1.d0/3.d0
      rkC(3)   = 2.d0/3.d0
      rkC(4)   = 1.d0
      rkC(5)   = 1.d0

      rkE(1)   =  1.d0/24.d0
      rkE(2)   =  -1.d0/8.d0
      rkE(3)   =  1.d0/8.d0
      rkE(4)   =  1.d0/8.d0
      rkE(5)   =  -1.d0/6.d0

      END SUBROUTINE Erk43


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Dopri5
!     7 stages, order 5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      erkMethod = RK5
      rkELO = 6.d0
      rkS      = 7

      rkA(2,1) = .2d0
      rkA(3,1) = 3.d0/40.d0
      rkA(3,2) = 9.d0/40.d0
      rkA(4,1) = 44.d0/45.d0
      rkA(4,2) = -56.d0/15.d0
      rkA(4,3) = 32.d0/9.d0
      rkA(5,1) = 19372.d0/6561.d0
      rkA(5,2) = -25360.d0/2187.d0
      rkA(5,3) = 64448.d0/6561.d0
      rkA(5,4) = -212.d0/729.d0
      rkA(6,1) = 9017.d0 / 3168.d0
      rkA(6,2) = -355.d0 / 33.d0
      rkA(6,3) = 46732.d0 / 5247.d0
      rkA(6,4) = 49.d0 / 176.d0
      rkA(6,5) = -5103.d0 / 18656.d0
      rkA(7,1) = 35.d0 / 384.d0
      rkA(7,3) = 500.d0 / 1113.d0
      rkA(7,4) = 125.d0 / 192.d0
      rkA(7,5) = -2187.d0 / 6784.d0
      rkA(7,6) = 11.d0 / 84.d0

      rkB(:)   = rkA(7,:)

      rkC(1)   = 0.d0
      rkC(2)   = .2d0
      rkC(3)   = .3d0
      rkC(4)   = .8d0
      rkC(5)   = 8.d0/9.d0
      rkC(6)   = 1.d0
      rkC(7)   = 1.d0

      rkE(1)   = 71.d0/57600.d0
      rkE(3)   = -71.d0/16695.d0
      rkE(4)   = 71.d0/1920.d0
      rkE(5)   = -17253.d0/339200.d0
      rkE(6)   = 22.d0/525.d0
      rkE(7)   = -1.d0/40.d0

      END SUBROUTINE Dopri5

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Verner
! 8 stages order 6
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      erkMethod = RK6
      rkELO = 7.d0
      rkS = 8

      rkA(2,1) = 1.d0/6.d0
      rkA(3,1) = 4.d0/75.d0
      rkA(3,2) = 16.d0/75.d0
      rkA(4,1) = 5.d0/6.d0
      rkA(4,2) = -8.d0/3.d0
      rkA(4,3) = 5.d0/2.d0
      rkA(5,1) = -165.d0/64.d0
      rkA(5,2) = 55.d0/6.d0
      rkA(5,3) = -425.d0/64.d0
      rkA(5,4) = 85.d0/96.d0
      rkA(6,1) = 12.d0/5.d0
      rkA(6,2) = -8.d0
      rkA(6,3) = 4015.d0/612.d0
      rkA(6,4) = -11.d0/36.d0
      rkA(6,5) = 88.d0/255.d0
      rkA(7,1) = -8263.d0/15000.d0
      rkA(7,2) = 124.d0/75.d0
      rkA(7,3) = -643.d0/680.d0
      rkA(7,4) = -81.d0/250.d0
      rkA(7,5) = 2484.d0/10625.d0
      rkA(7,6) = 0.d0
      rkA(8,1) = 3501.d0/1720.d0
      rkA(8,2) = -300.d0/43.d0
      rkA(8,3) = 297275.d0/52632.d0
      rkA(8,4) = -319.d0/2322.d0
      rkA(8,5) = 24068.d0/84065.d0
      rkA(8,6) = 0.d0
      rkA(8,7) = 3850.d0/26703.d0

      rkB(1) = 3.d0/40.d0
      rkB(2) = 0.d0
      rkB(3) = 875.d0/2244.d0
      rkB(4) = 23.d0/72.d0
      rkB(5) = 264.d0/1955.d0
      rkB(6) = 0.d0
      rkB(7) = 125.d0/11592.d0
      rkB(8) = 43.d0/616.d0
      
      rkC(2) = 1.d0/6.d0
      rkC(3) = 4.d0/15.d0
      rkC(4) = 2.d0/3.d0
      rkC(5) = 5.d0/6.d0
      rkC(6) = 1.0d0
      rkC(7) = 1.d0/15.d0
      rkC(8) = 1.d0

      rkE(1) = -1.d0/160.d0
      rkE(2) = 0.d0
      rkE(3) = 875.d0/2244.d0 - 2375.d0/5984.d0
      rkE(4) = 23.d0/72.d0 - 5.d0/16.d0
      rkE(5) = 264.d0/1955.d0 - 12.d0/85.d0
      rkE(6) = -3.d0/44.d0
      rkE(7) = 125.d0/11592.d0
      rkE(8) = 43.d0/616.d0

      END SUBROUTINE Verner


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Dopri853
!     12 stages order 8
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      erkMethod = RK8
      rkELO = 9.d0
      rkS = 12

      rkA(2,1) =    5.26001519587677318785587544488D-2
      rkA(3,1) =    1.97250569845378994544595329183D-2
      rkA(3,2) =    5.91751709536136983633785987549D-2
      rkA(4,1) =    2.95875854768068491816892993775D-2
      rkA(4,3) =    8.87627564304205475450678981324D-2
      rkA(5,1) =    2.41365134159266685502369798665D-1
      rkA(5,3) =   -8.84549479328286085344864962717D-1
      rkA(5,4) =    9.24834003261792003115737966543D-1
      rkA(6,1) =    3.7037037037037037037037037037D-2
      rkA(6,4) =    1.70828608729473871279604482173D-1
      rkA(6,5) =    1.25467687566822425016691814123D-1
      rkA(7,1) =    3.7109375D-2
      rkA(7,4) =    1.70252211019544039314978060272D-1
      rkA(7,5) =    6.02165389804559606850219397283D-2
      rkA(7,6) =   -1.7578125D-2

      rkA(8,1) =    3.70920001185047927108779319836D-2
      rkA(8,4) =    1.70383925712239993810214054705D-1
      rkA(8,5) =    1.07262030446373284651809199168D-1
      rkA(8,6) =   -1.53194377486244017527936158236D-2
      rkA(8,7) =    8.27378916381402288758473766002D-3
      rkA(9,1) =    6.24110958716075717114429577812D-1
      rkA(9,4) =   -3.36089262944694129406857109825D0
      rkA(9,5) =   -8.68219346841726006818189891453D-1
      rkA(9,6) =    2.75920996994467083049415600797D1
      rkA(9,7) =    2.01540675504778934086186788979D1
      rkA(9,8) =   -4.34898841810699588477366255144D1
      rkA(10,1) =   4.77662536438264365890433908527D-1
      rkA(10,4) =  -2.48811461997166764192642586468D0
      rkA(10,5) =  -5.90290826836842996371446475743D-1
      rkA(10,6) =   2.12300514481811942347288949897D1
      rkA(10,7) =   1.52792336328824235832596922938D1
      rkA(10,8) =  -3.32882109689848629194453265587D1
      rkA(10,9) =  -2.03312017085086261358222928593D-2

      rkA(11,1) =  -9.3714243008598732571704021658D-1
      rkA(11,4) =   5.18637242884406370830023853209D0
      rkA(11,5) =   1.09143734899672957818500254654D0
      rkA(11,6) =  -8.14978701074692612513997267357D0
      rkA(11,7) =  -1.85200656599969598641566180701D1
      rkA(11,8) =   2.27394870993505042818970056734D1
      rkA(11,9) =   2.49360555267965238987089396762D0
      rkA(11,10) = -3.0467644718982195003823669022D0
      rkA(12,1) =   2.27331014751653820792359768449D0
      rkA(12,4) =  -1.05344954667372501984066689879D1
      rkA(12,5) =  -2.00087205822486249909675718444D0
      rkA(12,6) =  -1.79589318631187989172765950534D1
      rkA(12,7) =   2.79488845294199600508499808837D1
      rkA(12,8) =  -2.85899827713502369474065508674D0
      rkA(12,9) =  -8.87285693353062954433549289258D0
      rkA(12,10) =  1.23605671757943030647266201528D1
      rkA(12,11) =  6.43392746015763530355970484046D-1

      rkB(1) =   5.42937341165687622380535766363D-2
      rkB(6) =   4.45031289275240888144113950566D0
      rkB(7) =   1.89151789931450038304281599044D0
      rkB(8) =  -5.8012039600105847814672114227D0
      rkB(9) =   3.1116436695781989440891606237D-1
      rkB(10) = -1.52160949662516078556178806805D-1
      rkB(11) =  2.01365400804030348374776537501D-1
      rkB(12) =  4.47106157277725905176885569043D-2

      rkC(2)  = 0.526001519587677318785587544488D-01
      rkC(3)  = 0.789002279381515978178381316732D-01
      rkC(4)  = 0.118350341907227396726757197510D+00
      rkC(5)  = 0.281649658092772603273242802490D+00
      rkC(6)  = 0.333333333333333333333333333333D+00
      rkC(7)  = 0.25D+00
      rkC(8)  = 0.307692307692307692307692307692D+00
      rkC(9)  = 0.651282051282051282051282051282D+00
      rkC(10) = 0.6D+00
      rkC(11) = 0.857142857142857142857142857142D+00
      rkC(12) = 1.d0

      rkE(1)  =  0.1312004499419488073250102996D-01
      rkE(6)  = -0.1225156446376204440720569753D+01
      rkE(7)  = -0.4957589496572501915214079952D+00
      rkE(8)  =  0.1664377182454986536961530415D+01
      rkE(9)  = -0.3503288487499736816886487290D+00
      rkE(10) =  0.3341791187130174790297318841D+00
      rkE(11) =  0.8192320648511571246570742613D-01
      rkE(12) = -0.2235530786388629525884427845D-01

      END SUBROUTINE Dopri853


     SUBROUTINE SET2ZERO(N,Y)
!--------------------------------------------------------------
!     copies zeros into the vector y:  y <- 0
!     after BLAS
!--------------------------------------------------------------

      INTEGER ::  i,M,MP1,N
      DOUBLE PRECISION ::  Y(N)
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0
      IF (N.LE.0) RETURN

      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = ZERO
        END DO
        IF( N .LT. 8 ) RETURN
      END IF
      MP1 = M+1
      DO i = MP1,N,8
        Y(i)     = ZERO
        Y(i + 1) = ZERO
        Y(i + 2) = ZERO
        Y(i + 3) = ZERO
        Y(i + 4) = ZERO
        Y(i + 5) = ZERO
        Y(i + 6) = ZERO
        Y(i + 7) = ZERO
      END DO
      END SUBROUTINE SET2ZERO

      SUBROUTINE WADD(N,X,Y,Z)
      INTEGER :: i, M, MP1, N
      DOUBLE PRECISION :: X(N),Y(N),Z(N)

      IF (N.LE.0) RETURN

      M = MOD(N,5)
      IF( M /= 0 ) THEN
         DO i = 1,M
            Z(i) = X(i) + Y(i)
         END DO
         IF( N < 5 ) RETURN
      END IF
      MP1 = M+1
      DO i = MP1,N,5
         Z(i)     = X(i)     + Y(i)
         Z(i + 1) = X(i + 1) + Y(i + 1)
         Z(i + 2) = X(i + 2) + Y(i + 2)
         Z(i + 3) = X(i + 3) + Y(i + 3)
         Z(i + 4) = X(i + 4) + Y(i + 4)
      END DO
      END SUBROUTINE WADD
   END SUBROUTINE ERK ! AND ALL ITS INTERNAL PROCEDURES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! End of INTEGRATE function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END MODULE ERK_f90_Integrator
