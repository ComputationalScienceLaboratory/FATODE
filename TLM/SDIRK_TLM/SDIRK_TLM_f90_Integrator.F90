! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! INTEGRATE - Integrator routine
!   Arguments :
!      TIN       - Start Time for Integration
!      TOUT      - End Time for Integration
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  SDIRK-TLM - Tangent Linear Model of Singly-Diagonally-Implicit RK      !
!            * Sdirk 2a, 2b: L-stable, 2 stages, order 2                  !
!            * Sdirk 3a:     L-stable, 3 stages, order 2, adj-invariant   !
!            * Sdirk 4a, 4b: L-stable, 5 stages, order 4                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE SDIRK_TLM_f90_Integrator

  IMPLICIT NONE
  PUBLIC
  SAVE
 
!~~~>  Statistics on the work performed by the SDIRK method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4,  &
           Nrej=5, Ndec=6, Nsol=7, Nsng=8,               &
           Ntexit=1, Nhexit=2, Nhnew=3
                 
CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INTEGRATE_TLM(N, NTLM, NNZERO, Y, Y_TLM, TIN, TOUT, ATOL_TLM,RTOL_TLM, &
 ATOL,RTOL,FUN, JAC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    IMPLICIT NONE
    INTEGER :: N,NNZERO
!~~~> Y - Concentrations
    DOUBLE PRECISION :: Y(N)
!~~~> NTLM - No. of sensitivity coefficients
    INTEGER NTLM
!~~~> Y_TLM - Sensitivities of concentrations
!     Note: Y_TLM (1:N,j) contains sensitivities of
!               Y(1:N) w.r.t. the j-th parameter, j=1...NTLM
    DOUBLE PRECISION :: Y_TLM(N,NTLM)
    DOUBLE PRECISION :: TIN  ! TIN - Start Time
    DOUBLE PRECISION :: TOUT ! TOUT - End Time
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RTOL_TLM(N,NTLM),ATOL_TLM(N,NTLM),&
                                           ATOL(N),RTOL(N)
    INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
    INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
    INTEGER,       INTENT(OUT), OPTIONAL :: IERR_U

    INTEGER :: IERR

    DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
    INTEGER       :: ICNTRL(20), ISTATUS(20)
    EXTERNAL :: FUN,JAC
    ICNTRL(1:20) = 0
    RCNTRL(1:20) = 0.0

    !~~~> fine-tune the integrator:
    ICNTRL(2)  = 0   ! 0=vector tolerances, 1=scalar tolerances
    ICNTRL(5)  = 8   ! Max no. of Newton iterations
    ICNTRL(6)  = 0   ! Starting values for Newton are interpolated (0) or zero (1)
    ICNTRL(7)  = 0   ! How to solve TLM: 0=modified Newton, 1=direct
    ICNTRL(9)  = 0   ! TLM Newton Iterations influence
    ICNTRL(12) = 0   ! TLM Truncation Error influence

    !~~~> if optional parameters are given, and if they are >0,
    !     then use them to overwrite default settings
    IF (PRESENT(ICNTRL_U)) THEN
      WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
    END IF
    IF (PRESENT(RCNTRL_U)) THEN
      WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
    END IF


    T1 = TIN; T2 = TOUT
    CALL SdirkTLM(  N, NTLM, NNZERO, T1, T2, Y, Y_TLM, RTOL, ATOL, &
                  RTOL_TLM, ATOL_TLM, FUN,JAC,RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)


    ! if optional parameters are given for output
    ! use them to store information in them
    IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
    IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
    IF (PRESENT(IERR_U)) IERR_U = IERR

    IF (IERR < 0) THEN
      PRINT *,'SDIRK-TLM: Unsuccessful exit at T=', TIN,' (IERR=',IERR,')'
    ENDIF

  END SUBROUTINE INTEGRATE_TLM



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SdirkTLM(N, NTLM, NNZERO, Tinitial, Tfinal, Y, Y_TLM, RelTol, AbsTol, &
                          RelTol_TLM, AbsTol_TLM,FUN, JAC,  &
			  RCNTRL, ICNTRL, RSTATUS, ISTATUS, Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using a Singly-Diagonally-Implicit
!    Runge-Kutta (SDIRK) method.
!
!    This implementation is based on the book and the code Sdirk4:
!
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.
!    This code is based on the SDIRK4 routine in the above book.
!
!    Methods:
!            * Sdirk 2a, 2b: L-stable, 2 stages, order 2                  
!            * Sdirk 3a:     L-stable, 3 stages, order 2, adjoint-invariant   
!            * Sdirk 4a, 4b: L-stable, 5 stages, order 4                  
!
!    (C)  Adrian Sandu, July 2005
!    Virginia Polytechnic Institute and State University
!    Contact: sandu@cs.vt.edu
!    Revised by Philipp Miehe and Adrian Sandu, May 2006
!    Revised by Hong Zhang and Adrian Sandu, Feb 2011
!    This implementation is part of FATODE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!-     Y(N)    = vector of initial conditions (at T=Tinitial)
!-    [Tinitial,Tfinal]  = time range of integration
!     (if Tinitial>Tfinal the integration is performed backwards in time)
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE ode_Fun( T, Y, Ydot ) = ODE function,
!                       returns Ydot = Y' = F(T,Y)
!- SUBROUTINE ode_Fun( T, Y, Ydot ) = Jacobian of the ODE function,
!                       returns Jcb = dF/dY
!-    ICNTRL(1:20)    = integer inputs parameters
!-    RCNTRL(1:20)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!-    Y(N)         -> vector of final states (at T->Tfinal)
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
!    ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3) = Method
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0 the default value of 100000 is used
!
!    ICNTRL(5)  -> maximum number of Newton iterations
!        For ICNTRL(5)=0 the default value of 8 is used
!
!    ICNTRL(6)  -> starting values of Newton iterations:
!        ICNTRL(6)=0 : starting values are interpolated (the default)
!        ICNTRL(6)=1 : starting values are zero
!
!    ICNTRL(7)  -> method to solve TLM equations:
!        ICNTRL(7)=0 : modified Newton re-using LU (the default)
!        ICNTRL(7)=1 : direct solution (additional one LU factorization per stage)
!
!    ICNTRL(9) -> switch for TLM Newton iteration error estimation strategy
!		ICNTRL(9) = 0: base number of iterations as forward solution
!		ICNTRL(9) = 1: use RTOL_TLM and ATOL_TLM to calculate
!				error estimation for TLM at Newton stages
!
!    ICNTRL(12) -> switch for TLM truncation error control
!		ICNTRL(12) = 0: TLM error is not used
!		ICNTRL(12) = 1: TLM error is computed and used
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
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!                  than the predicted value  (default=0.9)
!    RCNTRL(8)  -> ThetaMin. If Newton convergence rate smaller
!                  than ThetaMin the Jacobian is not recomputed;
!                  (default=0.001)
!    RCNTRL(9)  -> NewtonTol, stopping criterion for Newton's method
!                  (default=0.03)
!    RCNTRL(10) -> Qmin
!    RCNTRL(11) -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
!                  step size is kept constant and the LU factorization
!                  reused (default Qmin=1, Qmax=1.2)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call to Rosenbrock adds the current no. of fcn calls
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
      USE LS_Solver
      IMPLICIT NONE

! Arguments      
      INTEGER, INTENT(IN)          :: N, NTLM, NNZERO, ICNTRL(20)
      DOUBLE PRECISION, INTENT(IN)    :: Tinitial, Tfinal, &
                    RelTol(N), AbsTol(N), RCNTRL(20), &
		    RelTol_TLM(N,NTLM), AbsTol_TLM(N,NTLM)
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N), Y_TLM(N,NTLM)
      INTEGER, INTENT(OUT)         :: Ierr
      INTEGER, INTENT(INOUT)       :: ISTATUS(20) 
      DOUBLE PRECISION, INTENT(OUT)   :: RSTATUS(20)
       
!~~~>  SDIRK method coefficients, up to 5 stages
      INTEGER, PARAMETER :: Smax = 5
      INTEGER, PARAMETER :: S2A=1, S2B=2, S3A=3, S4A=4, S4B=5
      DOUBLE PRECISION :: rkGamma, rkA(Smax,Smax), rkB(Smax), rkC(Smax), &
                       rkD(Smax),  rkE(Smax), rkBhat(Smax), rkELO,    &
                       rkAlpha(Smax,Smax), rkTheta(Smax,Smax),DLAMCH
      INTEGER :: sdMethod, rkS ! The number of stages

! Local variables      
      DOUBLE PRECISION :: Hmin, Hmax, Hstart, Roundoff,    &
                       FacMin, Facmax, FacSafe, FacRej, &
                       ThetaMin, NewtonTol, Qmin, Qmax
      INTEGER       :: ITOL, NewtonMaxit, Max_no_steps, i
      LOGICAL       :: StartNewton, DirectTLM, TLMNewtonEst, TLMtruncErr
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
      EXTERNAL :: FUN,JAC

      rkA = ZERO
      rkB = ZERO
      rkC = ZERO
      rkD = ZERO
      rkE = ZERO
      rkBhat = ZERO
      rkAlpha = ZERO
      rkTheta = ZERO      
      Ierr = 0
      ISTATUS(1:20) = 0
      RSTATUS(1:20) = ZERO

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
      IF (ICNTRL(2) == 0) THEN
         ITOL = 1
      ELSE
         ITOL = 0
      END IF

!~~~> ICNTRL(3) - method selection       
      SELECT CASE (ICNTRL(3))
      CASE (1)
         CALL Sdirk2a
      CASE (2)
         CALL Sdirk2b
      CASE (3)
         CALL Sdirk3a
      CASE (0,4)
         CALL Sdirk4a
      CASE (5)
         CALL Sdirk4b
      CASE DEFAULT
         CALL Sdirk2a
      END SELECT
      
!~~~>   The maximum number of time steps admitted
      IF (ICNTRL(4) == 0) THEN
         Max_no_steps = 200000
      ELSEIF (ICNTRL(4) > 0) THEN
         Max_no_steps=ICNTRL(4)
      ELSE
         PRINT * ,'User-selected ICNTRL(4)=',ICNTRL(4)
         CALL SDIRK_ErrorMsg(-1,Tinitial,ZERO,Ierr)
   END IF
   !~~~> The maximum number of Newton iterations admitted
      IF(ICNTRL(5) == 0)THEN
         NewtonMaxit=8
      ELSE
         NewtonMaxit=ICNTRL(5)
         IF(NewtonMaxit <= 0)THEN
             PRINT * ,'User-selected ICNTRL(5)=',ICNTRL(5)
             CALL SDIRK_ErrorMsg(-2,Tinitial,ZERO,Ierr)
         END IF
      END IF
!~~~> StartNewton:  Use extrapolation for starting values of Newton iterations
      IF (ICNTRL(6) == 0) THEN
         StartNewton = .TRUE.
      ELSE
         StartNewton = .FALSE.
      END IF      
!~~~> Solve TLM equations directly or by Newton iterations 
      DirectTLM = (ICNTRL(7) == 1)
!~~~> Newton iteration error control selection  
      IF (ICNTRL(9) == 0) THEN 
         TLMNewtonEst = .FALSE.
      ELSE
         TLMNewtonEst = .TRUE.
      END IF      
!~~~> TLM truncation error control selection  
      IF (ICNTRL(12) == 0) THEN 
         TLMtruncErr = .FALSE.
      ELSE
         TLMtruncErr = .TRUE.
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
         CALL SDIRK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Upper bound on the step size: (positive value)
      IF (RCNTRL(2) == ZERO) THEN
         Hmax = ABS(Tfinal-Tinitial)
      ELSEIF (RCNTRL(2) > ZERO) THEN
         Hmax = MIN(ABS(RCNTRL(2)),ABS(Tfinal-Tinitial))
      ELSE
         PRINT * , 'User-selected RCNTRL(2)=', RCNTRL(2)
         CALL SDIRK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Starting step size: (positive value)
      IF (RCNTRL(3) == ZERO) THEN
         Hstart = MAX(Hmin,Roundoff)
      ELSEIF (RCNTRL(3) > ZERO) THEN
         Hstart = MIN(ABS(RCNTRL(3)),ABS(Tfinal-Tinitial))
      ELSE
         PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
         CALL SDIRK_ErrorMsg(-3,Tinitial,ZERO,Ierr)
      END IF
   
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax
      IF (RCNTRL(4) == ZERO) THEN
         FacMin = 0.2
      ELSEIF (RCNTRL(4) > ZERO) THEN
         FacMin = RCNTRL(4)
      ELSE
         PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
      IF (RCNTRL(5) == ZERO) THEN
         FacMax = 10.0
      ELSEIF (RCNTRL(5) > ZERO) THEN
         FacMax = RCNTRL(5)
      ELSE
         PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
      IF (RCNTRL(6) == ZERO) THEN
         FacRej = 0.1
      ELSEIF (RCNTRL(6) > ZERO) THEN
         FacRej = RCNTRL(6)
      ELSE
         PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
      IF (RCNTRL(7) == ZERO) THEN
         FacSafe = 0.9
      ELSEIF (RCNTRL(7) > ZERO) THEN
         FacSafe = RCNTRL(7)
      ELSE
         PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
         CALL SDIRK_ErrorMsg(-4,Tinitial,ZERO,Ierr)
      END IF

!~~~> ThetaMin: decides whether the Jacobian should be recomputed
      IF(RCNTRL(8) == 0.D0)THEN
         ThetaMin = 1.0d-3
      ELSE
         ThetaMin = RCNTRL(8)
      END IF

!~~~> Stopping criterion for Newton's method
      IF(RCNTRL(9) == ZERO)THEN
         NewtonTol = 3.0d-2
      ELSE
         NewtonTol = RCNTRL(9)
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
            CALL SDIRK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
         END IF
      ELSE
         DO i=1,N
            IF (AbsTol(i) <= 0.D0.OR.RelTol(i) <= 10.D0*Roundoff) THEN
              PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
              PRINT * , ' RelTol(',i,') = ',RelTol(i)
              CALL SDIRK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
            END IF
         END DO
      END IF
    
    IF (Ierr < 0) RETURN
    CALL LSS_Init(N,NNZERO)
    CALL SDIRK_IntegratorTLM( N,NTLM,Tinitial,Tfinal,Y,Y_TLM,Ierr )

    PRINT*,'Tangent Linear Model  STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)
    CALL LSS_Free

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CONTAINS  !  PROCEDURES INTERNAL TO SDIRK
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE SDIRK_IntegratorTLM( N,NTLM,Tinitial,Tfinal,Y,Y_TLM,Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE

!~~~> Arguments:      
      INTEGER, INTENT(IN) :: N, NTLM
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N), Y_TLM(N,NTLM)
      DOUBLE PRECISION, INTENT(IN) :: Tinitial, Tfinal
      INTEGER, INTENT(OUT) :: Ierr
      
!~~~> Local variables:      
      DOUBLE PRECISION :: Z(N,rkS), G(N), TMP(N),         &   
                       NewtonRate, SCAL(N), DZ(N),        &
                       T, H, Theta, Hratio, NewtonPredictedErr, &
                       Qnewton, Err, Fac, Hnew, Tdirection,     &
                       NewtonIncrement, NewtonIncrementOld, &
		       SCAL_TLM(N), Yerr(N), Yerr_TLM(N,NTLM), ThetaTLM
      DOUBLE PRECISION :: Z_TLM(N,rkS,NTLM)
      INTEGER :: itlm, j, IER, istage, NewtonIter, saveNiter, NewtonIterTLM
      LOGICAL :: Reject, FirstStep, SkipJac, SkipLU, NewtonDone, Transp
      DOUBLE PRECISION :: HGammaInv
      INTEGER :: ConsecutiveSng
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initializations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      T = Tinitial
      Tdirection = SIGN(ONE,Tfinal-Tinitial)
      H = MAX(ABS(Hmin),ABS(Hstart))
      IF (ABS(H) <= 10.D0*Roundoff) H=1.0D-6
      H=MIN(ABS(H),Hmax)
      H=SIGN(H,Tdirection)
      SkipLU  = .FALSE.
      SkipJac = .FALSE.
      Reject =  .FALSE.
      FirstStep=.TRUE.
      Transp = .FALSE.
      CALL SDIRK_ErrorScale(N, ITOL, AbsTol, RelTol, Y, SCAL)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( (Tfinal-T)*Tdirection - Roundoff > ZERO )


!~~~>  Compute E = 1/(h*gamma)-Jac and its LU decomposition
      IF ( .NOT.SkipLU ) THEN ! This time around skip the Jac update and LU
         ConsecutiveSng = 0
         IER = 1
Hloop1:   DO WHILE (IER /= 0)
           HGammaInv = ONE/(H*rkGamma)
!~~~>  Compute the Jacobian
           IF (.NOT. SkipJac) THEN
             CALL LSS_Jac(T,Y,JAC)
             ISTATUS(Njac) = ISTATUS(Njac) + 1
           END IF

           CALL LSS_Decomp(HGammaInv,IER)
           ISTATUS(Ndec) = ISTATUS(Ndec) + 1

           IF (IER /= 0) THEN
             WRITE (6,*) ' MATRIX IS SINGULAR, ISING=',IER,' T=',T,' H=',H
             ISTATUS(Nsng) = ISTATUS(Nsng) + 1;
             ConsecutiveSng = ConsecutiveSng + 1;
             IF (ConsecutiveSng >= 6) RETURN ! Failure
             H = 0.5d0*H
             SkipJac = .TRUE.
             SkipLU  = .FALSE.
             Reject  = .TRUE.
           END IF
      END DO Hloop1

!         CALL SDIRK_PrepareMatrix ( H, T, Y, FJAC, &
!                   SkipJac, SkipLU, E, IP, Reject, IER )
         IF (IER /= 0) THEN
             CALL SDIRK_ErrorMsg(-8,T,H,Ierr); RETURN
         END IF
      END IF      

      IF (ISTATUS(Nstp) > Max_no_steps) THEN
             CALL SDIRK_ErrorMsg(-6,T,H,Ierr); RETURN
      END IF   
      IF ( (T+0.1d0*H == T) .OR. (ABS(H) <= Roundoff) ) THEN
             CALL SDIRK_ErrorMsg(-7,T,H,Ierr); RETURN
      END IF   

stages:DO istage = 1, rkS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Simplified Newton iterations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~>  Starting values for Newton iterations
       CALL SET2ZERO(N,Z(1,istage))
       
!~~~>   Prepare the loop-independent part of the right-hand side
       CALL SET2ZERO(N,G)
       IF (istage > 1) THEN
           DO j = 1, istage-1
               ! Gj(:) = sum_j Theta(i,j)*Zj(:) = H * sum_j A(i,j)*Fun(Zj(:))
               CALL DAXPY(N,rkTheta(istage,j),Z(1,j),1,G,1)
               ! Zi(:) = sum_j Alpha(i,j)*Zj(:)
	       IF (StartNewton) THEN
                 CALL DAXPY(N,rkAlpha(istage,j),Z(1,j),1,Z(1,istage),1)
	       END IF
           END DO
       END IF

       !~~~>  Initializations for Newton iteration
       NewtonDone = .FALSE.
       Fac = 0.5d0 ! Step reduction factor if too many iterations
            
NewtonLoop:DO NewtonIter = 1, NewtonMaxit

!~~~>   Prepare the loop-dependent part of the right-hand side
 	    CALL WADD(N,Y,Z(1,istage),TMP)         	! TMP <- Y + Zi
            CALL FUN(N,T+rkC(istage)*H,TMP,DZ)	! DZ <- Fun(Y+Zi)
            ISTATUS(Nfun) = ISTATUS(Nfun) + 1
!            DZ(1:N) = G(1:N) - Z(1:N,istage) + (H*rkGamma)*DZ(1:N)
	    CALL DSCAL(N, H*rkGamma, DZ, 1)
	    CALL DAXPY (N, -ONE, Z(1,istage), 1, DZ, 1)
            CALL DAXPY (N, ONE, G,1, DZ,1)

!~~~>   Solve the linear system
            HGammaInv = ONE/(H*rkGamma)
            CALL DSCAL(N,HGammaInv,DZ,1)
            CALL LSS_Solve(Transp,DZ)
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1

!            CALL SDIRK_Solve ( H, N, E, IP, IER, DZ )
            
!~~~>   Check convergence of Newton iterations
            CALL SDIRK_ErrorNorm(N, DZ, SCAL, NewtonIncrement)
            IF ( NewtonIter == 1 ) THEN
                Theta      = ABS(ThetaMin)
                NewtonRate = 2.0d0 
            ELSE
                Theta = NewtonIncrement/NewtonIncrementOld
                IF (Theta < 0.99d0) THEN
                    NewtonRate = Theta/(ONE-Theta)
                    ! Predict error at the end of Newton process 
                    NewtonPredictedErr = NewtonIncrement &
                               *Theta**(NewtonMaxit-NewtonIter)/(ONE-Theta)
                    IF (NewtonPredictedErr >= NewtonTol) THEN
                      ! Non-convergence of Newton: predicted error too large
                      Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                      Fac = 0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                      EXIT NewtonLoop
                    END IF
                ELSE ! Non-convergence of Newton: Theta too large
                    EXIT NewtonLoop
                END IF
            END IF
            NewtonIncrementOld = NewtonIncrement
            ! Update solution: Z(:) <-- Z(:)+DZ(:)
            CALL DAXPY(N,ONE,DZ,1,Z(1,istage),1) 
            
            ! Check error in Newton iterations
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
            IF (NewtonDone) THEN
    	       ! Tune error in TLM variables by defining the minimal number of Newton iterations.
	       saveNiter = NewtonIter+1
	       EXIT NewtonLoop
	    END IF
            
            END DO NewtonLoop
            
            IF (.NOT.NewtonDone) THEN
                 !CALL RK_ErrorMsg(-12,T,H,Ierr);
                 H = Fac*H; Reject=.TRUE.
                 SkipJac = .TRUE.; SkipLU = .FALSE.
                 
                 CYCLE Tloop
            END IF

!~~~>  End of simplified Newton iterations for forward variables

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Solve for TLM variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~>  Direct solution for TLM variables
DirTLM:IF (DirectTLM) THEN 

         TMP(1:N) = Y(1:N) + Z(1:N,istage)
         SkipJac = .FALSE.
         ConsecutiveSng = 0
         IER = 1
Hloop2:   DO WHILE (IER /= 0)
           HGammaInv = ONE/(H*rkGamma)
!~~~>  Compute the Jacobian
           IF (.NOT. SkipJac) THEN
             CALL LSS_Jac1(T+rkC(istage)*H,TMP,JAC)
             ISTATUS(Njac) = ISTATUS(Njac) + 1
           END IF
           CALL LSS_Decomp_TLM(HGammaInv,IER)
           ISTATUS(Ndec) = ISTATUS(Ndec) + 1
           IF (IER /= 0) THEN
             WRITE (6,*) ' MATRIX IS SINGULAR, ISING=',IER,' T=',T,' H=',H
             ISTATUS(Nsng) = ISTATUS(Nsng) + 1;
             ConsecutiveSng = ConsecutiveSng + 1;
             IF (ConsecutiveSng >= 6) RETURN ! Failure
             H = 0.5d0*H
             SkipJac = .TRUE.
             SkipLU  = .FALSE.
             Reject  = .TRUE.
           END IF
      END DO Hloop2

!         CALL SDIRK_PrepareMatrix ( H, T+rkC(istage)*H, TMP, Jac, &
!                   SkipJac, SkipLU, E_TLM, IP_TLM, Reject, IER )
         IF (IER /= 0) CYCLE TLoop

TlmL:    DO itlm = 1, NTLM
            G(1:N) = Y_TLM(1:N,itlm)
            IF (istage > 1) THEN
              ! Gj(:) = sum_j Theta(i,j)*Zj_TLM(:) 
              !       = H * sum_j A(i,j)*Jac(Zj(:))*(Yj_TLM+Zj_TLM)
              DO j = 1, istage-1
                  CALL DAXPY(N,rkTheta(istage,j),Z_TLM(1,j,itlm),1,G,1)
              END DO
            END IF
!~~~>   Solve the linear system
            HGammaInv = ONE/(H*rkGamma)
            CALL DSCAL(N,HGammaInv,G,1)
            CALL LSS_Solve_TLM(Transp,G)
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1
            Z_TLM(1:N,istage,itlm) = G(1:N) - Y_TLM(1:N,itlm)
         END DO TlmL

       ELSE DirTLM

!~~~>  Jacobian of the current stage solution
       TMP(1:N) = Y(1:N) + Z(1:N,istage)
       CALL LSS_Jac1(T+rkC(istage)*H,TMP,JAC)
       ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Simplified Newton iterations for TLM variables
TlmLoop:DO itlm = 1,NTLM
       NewtonRate = MAX(NewtonRate,Roundoff)**0.8d0

!~~~>  Starting values for Newton iterations
       CALL SET2ZERO(N,Z_TLM(1,istage,itlm))
       
!~~~>   Prepare the loop-independent part of the right-hand side
       CALL LSS_Mul_Jac1(DZ,Y_TLM(1,itlm))
!       DZ = MATMUL(Jac,Y_TLM(1:N,itlm))    ! DZ <- Jac(Y+Z)*Y_TLM
       G(1:N) = (H*rkGamma)*DZ(1:N)
       IF (istage > 1) THEN
           ! Gj(:) = sum_j Theta(i,j)*Zj_TLM(:) 
           !       = H * sum_j A(i,j)*Jac(Zj(:))*(Yj_TLM+Zj_TLM)
           DO j = 1, istage-1
               CALL DAXPY(N,rkTheta(istage,j),Z_TLM(1,j,itlm),1,G,1)
           END DO
       END IF
       
       
       !~~~>  Initializations for Newton iteration
       IF (TLMNewtonEst) THEN
          NewtonDone = .FALSE.
          Fac = 0.5d0 ! Step reduction factor if too many iterations

          CALL SDIRK_ErrorScale(N,ITOL,AbsTol_TLM(1,itlm),RelTol_TLM(1,itlm), &
               Y_TLM(1,itlm),SCAL_TLM)
       END IF
            
NewtonLoopTLM:DO NewtonIterTLM = 1, NewtonMaxit

!~~~>   Prepare the loop-dependent part of the right-hand side
            CALL LSS_Mul_Jac1(DZ,Z_TLM(1,istage,itlm))
!            DZ = MATMUL(Jac,Z_TLM(1:N,istage,itlm))    ! DZ <- Jac(Y+Z)*Z_TLM
            DZ(1:N) = (H*rkGamma)*DZ(1:N)+G(1:N)-Z_TLM(1:N,istage,itlm)
!~~~>   Solve the linear system
            HGammaInv = ONE/(H*rkGamma)
            CALL DSCAL(N,HGammaInv,DZ,1)
            CALL LSS_Solve(Transp,DZ)
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1
            
	    IF (TLMNewtonEst) THEN
!~~~>   Check convergence of Newton iterations
            CALL SDIRK_ErrorNorm(N, DZ, SCAL_TLM, NewtonIncrement)
            IF ( NewtonIterTLM <= 1 ) THEN
                ThetaTLM      = ABS(ThetaMin)
                NewtonRate = 2.0d0 
            ELSE
                ThetaTLM = NewtonIncrement/NewtonIncrementOld
                IF (ThetaTLM < 0.99d0) THEN
                    NewtonRate = ThetaTLM/(ONE-ThetaTLM)
                    ! Predict error at the end of Newton process 
                    NewtonPredictedErr = NewtonIncrement &
                               *ThetaTLM**(NewtonMaxit-NewtonIterTLM)/(ONE-ThetaTLM)
                    IF (NewtonPredictedErr >= NewtonTol) THEN
                      ! Non-convergence of Newton: predicted error too large
                      Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                      Fac = 0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIterTLM))
                      EXIT NewtonLoopTLM
                    END IF
                ELSE ! Non-convergence of Newton: Theta too large
                    EXIT NewtonLoopTLM
                END IF
            END IF
            NewtonIncrementOld = NewtonIncrement
	    END IF !(TLMNewtonEst)
	    
            ! Update solution: Z_TLM(:) <-- Z_TLM(:)+DZ(:)
            CALL DAXPY(N,ONE,DZ,1,Z_TLM(1,istage,itlm),1) 
            
            ! Check error in Newton iterations
            IF (TLMNewtonEst) THEN
               NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
               IF (NewtonDone) EXIT NewtonLoopTLM
	    ELSE
  	       ! Minimum number of iterations same as FWD iterations
               IF (NewtonIterTLM>=saveNiter) EXIT NewtonLoopTLM
	    END IF
            
            END DO NewtonLoopTLM
            
            IF ((TLMNewtonEst) .AND. (.NOT.NewtonDone)) THEN
                 !CALL RK_ErrorMsg(-12,T,H,Ierr);
                 H = Fac*H; Reject=.TRUE.
                 SkipJac = .TRUE.; SkipLU = .FALSE.      
                 CYCLE Tloop
            END IF

      END DO TlmLoop
!~~~> End of simplified Newton iterations for TLM

    END IF DirTLM

   END DO stages

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Error estimation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      CALL SET2ZERO(N,Yerr)
      DO i = 1,rkS
         IF (rkE(i)/=ZERO) CALL DAXPY(N,rkE(i),Z(1,i),1,Yerr,1)
      END DO  
!~~~>   Solve the linear system
      HGammaInv = ONE/(H*rkGamma)
      CALL DSCAL(N,HGammaInv,Yerr,1)
      CALL LSS_Solve(Transp,Yerr)
      ISTATUS(Nsol) = ISTATUS(Nsol) + 1

      CALL SDIRK_ErrorNorm(N, Yerr, SCAL, Err)

      IF (TLMtruncErr) THEN
        CALL SET2ZERO(N*NTLM,Yerr_TLM)
        DO itlm=1,NTLM
          DO j=1,rkS  
            IF (rkE(j) /= ZERO) CALL DAXPY(N,rkE(j),Z_TLM(1,j,itlm),1,Yerr_TLM(1,itlm),1)
          END DO 
!~~~>   Solve the linear system
          CALL DSCAL(N,HGammaInv,Yerr_TLM(1,itlm),1)
          CALL LSS_Solve(Transp,Yerr_TLM(1,itlm))
          ISTATUS(Nsol) = ISTATUS(Nsol) + 1
        END DO
        CALL SDIRK_ErrorNorm_TLM(N,NTLM, Yerr_TLM, Err)
      END IF

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
         T  =  T + H
         ! Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
         DO i = 1,rkS 
            IF (rkD(i)/=ZERO) THEN
                 CALL DAXPY(N,rkD(i),Z(1,i),1,Y,1)
                 DO itlm = 1, NTLM
                   CALL DAXPY(N,rkD(i),Z_TLM(1,i,itlm),1,Y_TLM(1,itlm),1)
                 END DO  
            END IF     
         END DO  
       
!~~~> Update scaling coefficients
         CALL SDIRK_ErrorScale(N, ITOL, AbsTol, RelTol, Y, SCAL)

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
            ! If step not changed too much keep Jacobian and reuse LU
            SkipLU = ( (Theta <= ThetaMin) .AND. (Hratio >= Qmin) &
                                     .AND. (Hratio <= Qmax) )
            ! For TLM: do not skip LU (decrease TLM error)
   	    SkipLU = .FALSE.
            IF (.NOT.SkipLU) H = Hnew
         END IF
         ! If convergence is fast enough, do not update Jacobian
!         SkipJac = (Theta <= ThetaMin)
         SkipJac = .FALSE.

      ELSE accept !~~~> Step is rejected

         IF (FirstStep .OR. Reject) THEN
             H = FacRej*H
         ELSE
             H = Hnew
         END IF
         Reject = .TRUE.
         SkipJac = .TRUE.
         SkipLU  = .FALSE. 
         IF (ISTATUS(Nacc) >= 1) ISTATUS(Nrej) = ISTATUS(Nrej) + 1
    
      END IF accept
      END DO Tloop

      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE SDIRK_IntegratorTLM


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_ErrorScale(N, ITOL, AbsTol, RelTol, Y, SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER :: N, i, ITOL
      DOUBLE PRECISION :: AbsTol(N), RelTol(N), &
                       Y(N), SCAL(N)
      IF (ITOL == 0) THEN
        DO i=1,N
          SCAL(i) = ONE / ( AbsTol(1)+RelTol(1)*ABS(Y(i)) )
        END DO
      ELSE
        DO i=1,N
          SCAL(i) = ONE / ( AbsTol(i)+RelTol(i)*ABS(Y(i)) )
        END DO
      END IF
      END SUBROUTINE SDIRK_ErrorScale
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_ErrorNorm(N, Y, SCAL, Err)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      
      INTEGER :: i, N
      DOUBLE PRECISION :: Y(N), SCAL(N), Err      
      Err = ZERO
      DO i=1,N
           Err = Err+(Y(i)*SCAL(i))**2
      END DO
      Err = MAX( SQRT(Err/DBLE(N)), 1.0d-10 )
!
      END SUBROUTINE SDIRK_ErrorNorm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE SDIRK_ErrorNorm_TLM(N,NTLM, Y_TLM, FWD_Err)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      
      INTEGER :: itlm, NTLM, N
      DOUBLE PRECISION :: Y_TLM(N,NTLM), SCAL_TLM(N), FWD_Err, Err  
      
      DO itlm=1,NTLM
        CALL SDIRK_ErrorScale(N,ITOL,AbsTol_TLM(1,itlm),RelTol_TLM(1,itlm), &
               Y_TLM(1,itlm),SCAL_TLM)
	CALL SDIRK_ErrorNorm(N, Y_TLM(1,itlm), SCAL_TLM, Err)
	FWD_Err = MAX(FWD_Err, Err)
      END DO
!
      END SUBROUTINE SDIRK_ErrorNorm_TLM


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE SDIRK_ErrorMsg(Code,T,H,Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   DOUBLE PRECISION, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: Ierr

   Ierr = Code
   PRINT * , &
     'Forced exit from SDIRK due to the following error:'

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
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE SDIRK_ErrorMsg
      

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk4a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S4A
! Number of stages
      rkS = 5

! Method coefficients
      rkGamma = .2666666666666666666666666666666667d0

      rkA(1,1) = .2666666666666666666666666666666667d0
      rkA(2,1) = .5000000000000000000000000000000000d0
      rkA(2,2) = .2666666666666666666666666666666667d0
      rkA(3,1) = .3541539528432732316227461858529820d0
      rkA(3,2) = -.5415395284327323162274618585298197d-1
      rkA(3,3) = .2666666666666666666666666666666667d0
      rkA(4,1) = .8515494131138652076337791881433756d-1
      rkA(4,2) = -.6484332287891555171683963466229754d-1
      rkA(4,3) = .7915325296404206392428857585141242d-1
      rkA(4,4) = .2666666666666666666666666666666667d0
      rkA(5,1) = 2.100115700566932777970612055999074d0
      rkA(5,2) = -.7677800284445976813343102185062276d0
      rkA(5,3) = 2.399816361080026398094746205273880d0
      rkA(5,4) = -2.998818699869028161397714709433394d0
      rkA(5,5) = .2666666666666666666666666666666667d0

      rkB(1)   = 2.100115700566932777970612055999074d0
      rkB(2)   = -.7677800284445976813343102185062276d0
      rkB(3)   = 2.399816361080026398094746205273880d0
      rkB(4)   = -2.998818699869028161397714709433394d0
      rkB(5)   = .2666666666666666666666666666666667d0

      rkBhat(1)= 2.885264204387193942183851612883390d0
      rkBhat(2)= -.1458793482962771337341223443218041d0
      rkBhat(3)= 2.390008682465139866479830743628554d0
      rkBhat(4)= -4.129393538556056674929560012190140d0
      rkBhat(5)= 0.d0

      rkC(1)   = .2666666666666666666666666666666667d0
      rkC(2)   = .7666666666666666666666666666666667d0
      rkC(3)   = .5666666666666666666666666666666667d0
      rkC(4)   = .3661315380631796996374935266701191d0
      rkC(5)   = 1.d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.d0
      rkD(2)   = 0.d0
      rkD(3)   = 0.d0
      rkD(4)   = 0.d0
      rkD(5)   = 1.d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   = -.6804000050475287124787034884002302d0
      rkE(2)   = 1.558961944525217193393931795738823d0
      rkE(3)   = -13.55893003128907927748632408763868d0
      rkE(4)   = 15.48522576958521253098585004571302d0
      rkE(5)   = 1.d0

! Local order of Err estimate
      rkElo    = 4

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = 1.875000000000000000000000000000000d0
      rkTheta(3,1) = 1.708847304091539528432732316227462d0
      rkTheta(3,2) = -.2030773231622746185852981969486824d0
      rkTheta(4,1) = .2680325578937783958847157206823118d0
      rkTheta(4,2) = -.1828840955527181631794050728644549d0
      rkTheta(4,3) = .2968246986151577397160821594427966d0
      rkTheta(5,1) = .9096171815241460655379433581446771d0
      rkTheta(5,2) = -3.108254967778352416114774430509465d0
      rkTheta(5,3) = 12.33727431701306195581826123274001d0
      rkTheta(5,4) = -11.24557012450885560524143016037523d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = 2.875000000000000000000000000000000d0
      rkAlpha(3,1) = .8500000000000000000000000000000000d0
      rkAlpha(3,2) = .4434782608695652173913043478260870d0
      rkAlpha(4,1) = .7352046091658870564637910527807370d0
      rkAlpha(4,2) = -.9525565003057343527941920657462074d-1
      rkAlpha(4,3) = .4290111305453813852259481840631738d0
      rkAlpha(5,1) = -16.10898993405067684831655675112808d0
      rkAlpha(5,2) = 6.559571569643355712998131800797873d0
      rkAlpha(5,3) = -15.90772144271326504260996815012482d0
      rkAlpha(5,4) = 25.34908987169226073668861694892683d0
               
!~~~> Coefficients for continuous solution
!          rkD(1,1)= 24.74416644927758d0
!          rkD(1,2)= -4.325375951824688d0
!          rkD(1,3)= 41.39683763286316d0
!          rkD(1,4)= -61.04144619901784d0
!          rkD(1,5)= -3.391332232917013d0
!          rkD(2,1)= -51.98245719616925d0
!          rkD(2,2)= 10.52501981094525d0
!          rkD(2,3)= -154.2067922191855d0
!          rkD(2,4)= 214.3082125319825d0
!          rkD(2,5)= 14.71166018088679d0
!          rkD(3,1)= 33.14347947522142d0
!          rkD(3,2)= -19.72986789558523d0
!          rkD(3,3)= 230.4878502285804d0
!          rkD(3,4)= -287.6629744338197d0
!          rkD(3,5)= -18.99932366302254d0
!          rkD(4,1)= -5.905188728329743d0
!          rkD(4,2)= 13.53022403646467d0
!          rkD(4,3)= -117.6778956422581d0
!          rkD(4,4)= 134.3962081008550d0
!          rkD(4,5)= 8.678995715052762d0
!
!         DO i=1,4  ! CONTi <-- Sum_j rkD(i,j)*Zj
!           CALL SET2ZERO(N,CONT(1,i))
!           DO j = 1,rkS
!             CALL DAXPY(N,rkD(i,j),Z(1,j),1,CONT(1,i),1)
!           END DO
!         END DO
          
          rkELO = 4.0d0
          
      END SUBROUTINE Sdirk4a

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk4b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S4B
! Number of stages
      rkS = 5

! Method coefficients
      rkGamma = .25d0

      rkA(1,1) = 0.25d0
      rkA(2,1) = 0.5d00
      rkA(2,2) = 0.25d0
      rkA(3,1) = 0.34d0
      rkA(3,2) =-0.40d-1
      rkA(3,3) = 0.25d0
      rkA(4,1) = 0.2727941176470588235294117647058824d0
      rkA(4,2) =-0.5036764705882352941176470588235294d-1
      rkA(4,3) = 0.2757352941176470588235294117647059d-1
      rkA(4,4) = 0.25d0
      rkA(5,1) = 1.041666666666666666666666666666667d0
      rkA(5,2) =-1.020833333333333333333333333333333d0
      rkA(5,3) = 7.812500000000000000000000000000000d0
      rkA(5,4) =-7.083333333333333333333333333333333d0
      rkA(5,5) = 0.25d0

      rkB(1)   =  1.041666666666666666666666666666667d0
      rkB(2)   = -1.020833333333333333333333333333333d0
      rkB(3)   =  7.812500000000000000000000000000000d0
      rkB(4)   = -7.083333333333333333333333333333333d0
      rkB(5)   =  0.250000000000000000000000000000000d0

      rkBhat(1)=  1.069791666666666666666666666666667d0
      rkBhat(2)= -0.894270833333333333333333333333333d0
      rkBhat(3)=  7.695312500000000000000000000000000d0
      rkBhat(4)= -7.083333333333333333333333333333333d0
      rkBhat(5)=  0.212500000000000000000000000000000d0

      rkC(1)   = 0.25d0
      rkC(2)   = 0.75d0
      rkC(3)   = 0.55d0
      rkC(4)   = 0.50d0
      rkC(5)   = 1.00d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.0d0
      rkD(2)   = 0.0d0
      rkD(3)   = 0.0d0
      rkD(4)   = 0.0d0
      rkD(5)   = 1.0d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   =  0.5750d0
      rkE(2)   =  0.2125d0
      rkE(3)   = -4.6875d0
      rkE(4)   =  4.2500d0
      rkE(5)   =  0.1500d0

! Local order of Err estimate
      rkElo    = 4

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = 2.d0
      rkTheta(3,1) = 1.680000000000000000000000000000000d0
      rkTheta(3,2) = -.1600000000000000000000000000000000d0
      rkTheta(4,1) = 1.308823529411764705882352941176471d0
      rkTheta(4,2) = -.1838235294117647058823529411764706d0
      rkTheta(4,3) = 0.1102941176470588235294117647058824d0
      rkTheta(5,1) = -3.083333333333333333333333333333333d0
      rkTheta(5,2) = -4.291666666666666666666666666666667d0
      rkTheta(5,3) =  34.37500000000000000000000000000000d0
      rkTheta(5,4) = -28.33333333333333333333333333333333d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = 3.
      rkAlpha(3,1) = .8800000000000000000000000000000000d0
      rkAlpha(3,2) = .4400000000000000000000000000000000d0
      rkAlpha(4,1) = .1666666666666666666666666666666667d0
      rkAlpha(4,2) = -.8333333333333333333333333333333333d-1
      rkAlpha(4,3) = .9469696969696969696969696969696970d0
      rkAlpha(5,1) = -6.d0
      rkAlpha(5,2) = 9.d0
      rkAlpha(5,3) = -56.81818181818181818181818181818182d0
      rkAlpha(5,4) = 54.d0
      
      END SUBROUTINE Sdirk4b

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk2a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S2A
! Number of stages
      rkS = 2

! Method coefficients
      rkGamma = .2928932188134524755991556378951510d0

      rkA(1,1) = .2928932188134524755991556378951510d0
      rkA(2,1) = .7071067811865475244008443621048490d0
      rkA(2,2) = .2928932188134524755991556378951510d0

      rkB(1)   = .7071067811865475244008443621048490d0
      rkB(2)   = .2928932188134524755991556378951510d0

      rkBhat(1)= .6666666666666666666666666666666667d0
      rkBhat(2)= .3333333333333333333333333333333333d0

      rkC(1)   = 0.292893218813452475599155637895151d0
      rkC(2)   = 1.0d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.0d0
      rkD(2)   = 1.0d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   =  0.4714045207910316829338962414032326d0
      rkE(2)   = -0.1380711874576983496005629080698993d0

! Local order of Err estimate
      rkElo    = 2

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = 2.414213562373095048801688724209698d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = 3.414213562373095048801688724209698d0
          
      END SUBROUTINE Sdirk2a

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk2b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S2B
! Number of stages
      rkS      = 2

! Method coefficients
      rkGamma  = 1.707106781186547524400844362104849d0

      rkA(1,1) = 1.707106781186547524400844362104849d0
      rkA(2,1) = -.707106781186547524400844362104849d0
      rkA(2,2) = 1.707106781186547524400844362104849d0

      rkB(1)   = -.707106781186547524400844362104849d0
      rkB(2)   = 1.707106781186547524400844362104849d0

      rkBhat(1)= .6666666666666666666666666666666667d0
      rkBhat(2)= .3333333333333333333333333333333333d0

      rkC(1)   = 1.707106781186547524400844362104849d0
      rkC(2)   = 1.0d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.0d0
      rkD(2)   = 1.0d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   = -.4714045207910316829338962414032326d0
      rkE(2)   =  .8047378541243650162672295747365659d0

! Local order of Err estimate
      rkElo    = 2

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) = -.414213562373095048801688724209698d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) = .5857864376269049511983112757903019d0
      
      END SUBROUTINE Sdirk2b


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Sdirk3a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sdMethod = S3A
! Number of stages
      rkS = 3

! Method coefficients
      rkGamma = .2113248654051871177454256097490213d0

      rkA(1,1) = .2113248654051871177454256097490213d0
      rkA(2,1) = .2113248654051871177454256097490213d0
      rkA(2,2) = .2113248654051871177454256097490213d0
      rkA(3,1) = .2113248654051871177454256097490213d0
      rkA(3,2) = .5773502691896257645091487805019573d0
      rkA(3,3) = .2113248654051871177454256097490213d0

      rkB(1)   = .2113248654051871177454256097490213d0
      rkB(2)   = .5773502691896257645091487805019573d0
      rkB(3)   = .2113248654051871177454256097490213d0

      rkBhat(1)= .2113248654051871177454256097490213d0
      rkBhat(2)= .6477918909913548037576239837516312d0
      rkBhat(3)= .1408832436034580784969504064993475d0

      rkC(1)   = .2113248654051871177454256097490213d0
      rkC(2)   = .4226497308103742354908512194980427d0
      rkC(3)   = 1.d0

! Ynew = Yold + h*Sum_i {rkB_i*k_i} = Yold + Sum_i {rkD_i*Z_i}
      rkD(1)   = 0.d0
      rkD(2)   = 0.d0
      rkD(3)   = 1.d0

! Err = h * Sum_i {(rkB_i-rkBhat_i)*k_i} = Sum_i {rkE_i*Z_i}
      rkE(1)   =  0.9106836025229590978424821138352906d0
      rkE(2)   = -1.244016935856292431175815447168624d0
      rkE(3)   =  0.3333333333333333333333333333333333d0

! Local order of Err estimate
      rkElo    = 2

! h*Sum_j {rkA_ij*k_j} = Sum_j {rkTheta_ij*Z_j}
      rkTheta(2,1) =  1.0d0
      rkTheta(3,1) = -1.732050807568877293527446341505872d0
      rkTheta(3,2) =  2.732050807568877293527446341505872d0

! Starting value for Newton iterations: Z_i^0 = Sum_j {rkAlpha_ij*Z_j}
      rkAlpha(2,1) =   2.0d0
      rkAlpha(3,1) = -12.92820323027550917410978536602349d0
      rkAlpha(3,2) =   8.83012701892219323381861585376468d0

      END SUBROUTINE Sdirk3a
!--------------------------------------------------------------
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
!--------------------------------------------------------------
      SUBROUTINE WADD(N,X,Y,Z)
!--------------------------------------------------------------
!     adds two vectors: z <- x + y
!     BLAS - like
!--------------------------------------------------------------
!     USE rkadj_Precision

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


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   END SUBROUTINE SdirkTLM ! AND ALL ITS INTERNAL PROCEDURES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


END MODULE SDIRK_TLM_f90_Integrator
! End of INTEGRATE function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


