! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! INTEGRATE - Integrator routine
!   Arguments :
!      TIN       - Start Time for Integration
!      TOUT      - End Time for Integration
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  ERK_ADJ - Adjoint Model of Explicit Runge-Kutta methods                !
!            * RK2(3): 3 stages, order 2                                  !
!            * RK3(2): 4 stages, order 3                                  !
!            * RK4(3): 5 stages, order 4                                  !
!            * Dopri5: 7 stages, order 5                                  !
!            * Verner: 8 stages, order 6                                  !
!            * Dopri853: 12 stages, order 8                               !
!                                                                         !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE ERK_adj_f90_Integrator

  IMPLICIT NONE
  PUBLIC
  SAVE
  
!~~~>  Statistics on the work performed by the ERK method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, Nrej=5, Ndec=6, &
                        Nsol=7, Nsng=8, Ntexit=1, Nhexit=2, Nhnew=3
                 
CONTAINS

  SUBROUTINE INTEGRATE_ADJ( NVAR, NP, NADJ, NNZ, Y, Lambda, TIN, TOUT, &
           ATOL, RTOL, FUN, JAC, AdjInit,  ICNTRL_U, RCNTRL_U, ISTATUS_U, &
           RSTATUS_U, Ierr_U, Mu, JACP, DRDY, DRDP, QFUN, Q )

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NVAR,NP,NADJ
   INTEGER, INTENT(IN), OPTIONAL :: NNZ

   DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)
!  Quadrature term
   DOUBLE PRECISION, INTENT(INOUT),OPTIONAL :: Q(NADJ)
!~~~> NADJ - No. of cost functionals for which adjoints
!                are evaluated simultaneously
!            If single cost functional is considered (like in
!                most applications) simply set NADJ = 1      
!~~~> Lambda - Sensitivities w.r.t. concentrations
!     Note: Lambda (1:N,j) contains sensitivities of
!           the j-th cost functional w.r.t. Y(1:N), j=1...NADJ
   DOUBLE PRECISION,INTENT(INOUT)  :: Lambda(NVAR,NADJ)
   DOUBLE PRECISION, INTENT(INOUT), OPTIONAL ::  Mu(NP, NADJ)
!~~~> Tolerances for adjoint calculations
!     (used for full continuous adjoint, and for controlling
!      iterations when used to solve the discrete adjoint)   
   DOUBLE PRECISION, INTENT(IN)  ::  ATOL(NVAR), RTOL(NVAR)
   DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time
   ! Optional input parameters and statistics
   INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,       INTENT(OUT), OPTIONAL :: Ierr_U

   DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
   INTEGER       :: ICNTRL(20), ISTATUS(20), Ierr
   EXTERNAL      FUN,JAC,AdjInit
   EXTERNAL      QFUN,JACP,DRDY,DRDP
   OPTIONAL      QFUN,JACP,DRDY,DRDP
   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0

   ! If optional parameters are given, and if they are >0, 
   ! then they overwrite default settings. 
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
      
   T1 = TIN; T2 = TOUT

! Evaluating sensitivities w.r.t parameters requires NP>0, functions MU and JACP
! are provided.
   IF(NP>0 .AND. PRESENT(Mu) .AND. PRESENT(JACP)) THEN
     IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDY) .AND. PRESENT(DRDP)) THEN
!~~~> This is the case that cost function contains a quadrature term NADJ 
!     should be 1; Q is defined to store the value of the quadrature
!     term at the last step and functions QFUN,DRDY,DRDP must be provided.
!     cost function = g(y,p) + Integrate{r(y,p)}
       CALL ERKADJ1( N=NVAR, NP=NP, NADJ=NADJ, Tinitial=T1, Tfinal=T2, Y=Y,   &
                     Lambda=Lambda, Mu=Mu, FUN=FUN, JAC=JAC, RelTol=RTOL,     &
                     AbsTol=ATOL, JACP=JACP, DRDY=DRDY, DRDP=DRDP, QFUN=QFUN, &
                     AdjInit=AdjInit, RCNTRL=RCNTRL, ICNTRL=ICNTRL,           &
                     RSTATUS=RSTATUS, ISTATUS=ISTATUS, Ierr=Ierr, Q=Q )      
     ELSE
!~~~> No quadrature term is involved. cost function = g(y,p)
       CALL ERKADJ1( N=NVAR, NP=NP, NADJ=NADJ, Tinitial=T1, Tfinal=T2, Y=Y,   &
                     Lambda=Lambda, Mu=Mu, FUN=FUN, JAC=JAC, RelTol=RTOL,     &
                     AbsTol=ATOL, JACP=JACP,                                  &
                     AdjInit=AdjInit, RCNTRL=RCNTRL, ICNTRL=ICNTRL,           &
                     RSTATUS=RSTATUS, ISTATUS=ISTATUS, Ierr=Ierr)

     END IF
   ELSE
! Evaluating sensitivites w.r.t only initial conditions
     IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDY) )THEN
!~~~> This is the case that cost function contains a quadrature term
! NADJ should be 1; Q is defined to store the value of the quadrature
! term at the last step and functions QFUN,DRDY must be provided. 
! cost function = g(y) + Integrate{r(y)}   
       CALL ERKADJ2( N=NVAR, NP=NP, NADJ=NADJ, Tinitial=T1, Tfinal=T2, Y=Y,   &
                     Lambda=Lambda, FUN=FUN, JAC=JAC, RelTol=RTOL,            &
                     AbsTol=ATOL, DRDY=DRDY, QFUN=QFUN, Q=Q,                  &
                     AdjInit=AdjInit, RCNTRL=RCNTRL, ICNTRL=ICNTRL,           &
                     RSTATUS=RSTATUS, ISTATUS=ISTATUS, Ierr=Ierr)
     ELSE
!~~~> No quadrature term is involved. cost function = g(y)   
       CALL ERKADJ2( N=NVAR, NP=NP, NADJ=NADJ, Tinitial=T1, Tfinal=T2, Y=Y,   &
                     Lambda=Lambda, FUN=FUN, JAC=JAC, RelTol=RTOL,            &
                     AbsTol=ATOL, AdjInit=AdjInit, RCNTRL=RCNTRL,            &
                     ICNTRL=ICNTRL, RSTATUS=RSTATUS, ISTATUS=ISTATUS, Ierr=Ierr)
     END IF
   END IF

   IF (Ierr < 0) THEN
        PRINT *,'ERK: Unsuccessful exit at T=',TIN,' (Ierr=',Ierr,')'
   ENDIF
   
   ! if optional parameters are given for output they to return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(Ierr_U))    Ierr_U       = Ierr

   END SUBROUTINE INTEGRATE_ADJ


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERKADJ1(N, NP, NADJ, Tinitial, Tfinal, Y, Lambda, Mu, FUN, JAC, &
                      JACP, RelTol, AbsTol, AdjInit, RCNTRL, ICNTRL, RSTATUS, &
                      ISTATUS, Ierr, DRDY, DRDP, QFUN, Q)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using an Explicit
!    Runge-Kutta (ERK) method.
!
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
!        For ICNTRL(4)=0 the default value of 10000 is used
!        Note: use a conservative estimate, since the checkpoint
!              buffers are allocated to hold Max_no_steps
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
      IMPLICIT NONE

! Arguments      
      INTEGER, INTENT(IN)          :: N, NP,NADJ,ICNTRL(20)
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N), Lambda(N,NADJ),Mu(NP,NADJ)
      DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: Q(NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tinitial, Tfinal, &
                    RelTol(N), AbsTol(N), RCNTRL(20)
      INTEGER, INTENT(OUT)         :: Ierr
      INTEGER, INTENT(INOUT)       :: ISTATUS(20) 
      DOUBLE PRECISION, INTENT(OUT)   :: RSTATUS(20)
       
!~~~>  ERK method coefficients, up to 5 stages
      INTEGER, PARAMETER :: Smax = 12
      INTEGER, PARAMETER :: RK2=1, RK3=2, RK4=3, RK5=4, RK6=5, RK8=6
      DOUBLE PRECISION :: rkA(Smax,Smax), rkB(Smax), rkC(Smax),rkE(Smax),rkELO
      INTEGER :: erkMethod, rkS ! The number of stages

!~~~>  Checkpoints in memory buffers
      INTEGER :: stack_ptr = 0 ! last written entry in checkpoint
      DOUBLE PRECISION, DIMENSION(:),     POINTER :: chk_H, chk_T
      DOUBLE PRECISION, DIMENSION(:,:),   POINTER :: chk_Y
      DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: chk_Z
! Local variables      
      DOUBLE PRECISION :: Hmin, Hmax, Hstart, Roundoff,    &
                       FacMin, Facmax, FacSafe, FacRej, &
                       Qmin, Qmax

      INTEGER       :: ITOL, Max_no_steps, i
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
      DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
      LOGICAL :: GetQuad
      EXTERNAL      QFUN,FUN,JAC,JACP,DRDY,DRDP,AdjInit
      OPTIONAL      QFUN,DRDY,DRDP

      DOUBLE PRECISION :: DLAMCH
      
      rkA = ZERO
      rkB = ZERO
      rkC = ZERO
      rkE = ZERO
      GetQuad = .FALSE.
!~~~> Compute the quadrature term if required sourcs are provided
      IF( PRESENT(Q) .AND. PRESENT(QFUN) .AND. &
        PRESENT(DRDP) .AND. PRESENT(DRDY) )THEN
          GetQuad = .TRUE.
      END IF

!~~~>  Initialize statistics
      ISTATUS(1:20) = 0
      RSTATUS(1:20) = ZERO
      Ierr          = 0

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
         Max_no_steps = ICNTRL(4)
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
         DO i=1,N
            IF (AbsTol(i) <= 0.D0.OR.RelTol(i) <= 10.D0*Roundoff) THEN
              PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
              PRINT * , ' RelTol(',i,') = ',RelTol(i)
              CALL ERK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
            END IF
         END DO
      END IF
    
    IF (Ierr < 0) RETURN
!~~~>  Allocate memory buffers    
    CALL ERK_AllocBuffers

!~~~>  Call forward integration    
    CALL ERK_FwdInt( N, NADJ, Tinitial, Tfinal, Y, GetQuad, Ierr )

    PRINT*,'FORWARD STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

!~~~> Initialize adjoint variables
    CALL AdjInit(N, NP, NADJ, Tfinal, Y, Lambda, Mu)

!~~~>  Call adjoint integration    
    CALL ERK_DadjInt( N, NP, NADJ, Lambda, Mu, GetQuad, Ierr )

    PRINT*,'ADJOINT STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)
!~~~>  Free memory buffers    
    CALL ERK_FreeBuffers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CONTAINS  !  Procedures internal to ERKADJ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERK_FwdInt( NVAR,NADJ,Tinitial,Tfinal,Y,GetQuad,Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE

!~~~> Arguments:      
      INTEGER,INTENT(IN) :: NVAR,NADJ
      DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)
      DOUBLE PRECISION, INTENT(IN) :: Tinitial, Tfinal
      INTEGER, INTENT(OUT) :: Ierr
      LOGICAL, INTENT(IN) :: GetQuad 
!~~~> Local variables:      
      DOUBLE PRECISION :: Z(NVAR,rkS), K(NVAR,rkS), TMP(NVAR), SCAL(NVAR),    &
                          G(NVAR), R(NADJ), T, H, Hratio, Rerr, Fac, Hnew, Tdirection
      INTEGER :: j, istage
      LOGICAL :: Reject, FirstStep
      
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

      CALL ERK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)

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

       CALL Set2zero(NVAR,Z(1,istage))

       IF(istage > 1) THEN
           DO j = 1, istage-1
               CALL DAXPY(NVAR,H*rkA(istage,j),K(1,j),1,Z(1,istage),1)
           END DO
       END IF

       CALL DAXPY(NVAR,ONE,Y,1,Z(1,istage),1)

       CALL FUN(NVAR,T+rkC(istage)*H,Z(1,istage),K(1,istage))

       ISTATUS(Nfun) = ISTATUS(Nfun) + 1

  END DO stages

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Error estimation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      CALL SET2ZERO(N,TMP)
      DO i = 1,rkS
         IF (rkE(i)/=ZERO) CALL DAXPY(N,H*rkE(i),K(1,i),1,TMP,1)
      END DO  

      CALL ERK_ErrorNorm(N, TMP, SCAL, Rerr)

!~~~> Computation of new step size Hnew
      Fac  = FacSafe*(Rerr)**(-ONE/rkELO)
      Fac  = MAX(FacMin,MIN(FacMax,Fac))
      Hnew = H*Fac

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Accept/Reject step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
accept: IF ( Rerr < ONE ) THEN !~~~> Step is accepted

         FirstStep=.FALSE.
         ISTATUS(Nacc) = ISTATUS(Nacc) + 1

!~~~> Update the results for the quadrature term
         IF(GetQuad) Then
           CALL SET2ZERO(N,TMP)
           CALL SET2ZERO(NADJ,R)
           DO i = 1, rkS
             CALL WADD(N,Y,Z(1,i),TMP)
             CALL QFUN(N,NADJ,T+rkC(i)*H,TMP,R)
             CALL DAXPY(NADJ,H*rkB(i),R,1,Q,1)
           END DO
         END IF

!~~~> Checkpoint solution
         CALL ERK_Push( T, H, Y, Z )

!~~~> Update time and solution
         T  =  T + H
         ! Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
         DO i = 1,rkS 
            IF (rkB(i)/=ZERO) CALL DAXPY(N,H*rkB(i),K(1,i),1,Y,1)
         END DO  
       
!~~~> Update scaling coefficients
         CALL ERK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)

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
         
      END IF accept
      END DO Tloop

      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE ERK_FwdInt


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERK_DadjInt( N, NP, NADJ, Lambda, Mu, GetQuad, Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE

!~~~> Arguments:      
      INTEGER, INTENT(IN) :: N, NP, NADJ
      DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ), Mu(NP,NADJ)
      INTEGER, INTENT(OUT) :: Ierr
      LOGICAL,INTENT(IN) :: GetQuad
!~~~> Local variables:      
      DOUBLE PRECISION :: Alpha,Beta
      DOUBLE PRECISION :: Y(N)
      DOUBLE PRECISION :: Z(N,rkS),TMP(N), T, H
      DOUBLE PRECISION,ALLOCATABLE :: WY(:,:),WP(:,:),FPJAC(:,:),FJAC(:,:), U(:,:,:), V(:,:,:)
      INTEGER :: j, istage, iadj

      IF(GetQuad) THEN
        ALLOCATE(WY(N,NADJ),WP(NP,NADJ),STAT=j)
        IF(j .NE. 0) STOP 'allocation error for WY,WP'
        WY(:,:) = 0.0d0
        WP(:,:) = 0.0d0
      END IF
!      initiliazation
      ALLOCATE(FPJAC(N,NP),FJAC(N,N),U(N,NADJ,rkS),V(NP,NADJ,rkS),STAT=j)
      IF(j .NE. 0) STOP 'allocation error for FPJAC,FJAC,U,V'
      FPJAC(:,:) = 0.0d0
      
      Alpha = 1.0d0
      Beta = 0.0d0     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( stack_ptr > 0 )
        
   !~~~>  Recover checkpoints for stage values and vectors
      CALL ERK_Pop( T, H, Y, Z)

      ISTATUS(Nstp) = ISTATUS(Nstp) + 1

stages:DO istage = rkS, 1, -1

!~~~>  Jacobian of the current stage solution
      CALL Jac(N,T+rkC(istage)*H,Z(1,istage),FJAC)
      ISTATUS(Njac) = ISTATUS(Njac) + 1

      IF(GetQuad) CALL DRDY(NADJ,N,N,T+rkC(istage)*H,Z(1,istage),WY)
      CALL JACP(N,NP,T+rkC(istage)*H,Z(1,istage),FPJAC)
      IF(GetQuad) CALL DRDP(NADJ,N,NP,T+rkC(istage)*H,Z(1,istage),WP)
       
adj:   DO iadj = 1, NADJ
 

        CALL SET2ZERO(N,TMP)
        IF (istage < rkS) THEN
          DO j = istage+1, rkS
            CALL DAXPY(N,H*rkA(j,istage),U(1,iadj,j),1,TMP,1) ! TMP < h*A*U
          END DO
        END IF
        CALL DAXPY(N,H*rkB(istage),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B*Lambda

!~~~~>f_y^T * (b_i * Lambda_{n+1} + \sum_{j=i+1}^s a_{j,i} u_j)

        CALL DGEMV('T',N,N,ALPHA,FJAC,N,TMP,1,BETA,U(1,iadj,istage),1)  ! U = Jac^T * TMP

        CALL DGEMV('T',N,NP,Alpha,FPJAC,N,TMP,1,Beta,V(1,iadj,istage),1) ! V = Jacp^T *TMP
       
        IF(GetQuad) CALL DAXPY(N,H*rkB(istage),WY(1,iadj),1,U(1,iadj,istage),1) ! U = U + h*B*WY
        IF(GetQuad) CALL DAXPY(NP,H*rkB(istage),WP(1,iadj),1,V(1,iadj,istage),1) ! V = V + h*B*WP      
        
        END DO adj
 
      END DO stages

!~~~> Update adjoint solution
         ! Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
         DO istage = 1,rkS 
            DO iadj = 1,NADJ
                CALL DAXPY(N,ONE,U(1,iadj,istage),1,Lambda(1,iadj),1)
                !Lambda(1:N,iadj) = Lambda(1:N,iadj) + U(1:N,iadj,istage)
                CALL DAXPY(NP,ONE,V(1,iadj,istage),1,Mu(1,iadj),1) 
                !Mu(1:NP,iadj) = Mu(1:NP,iadj) + V(1:NP,iadj,istage)
            END DO
         END DO  

      END DO Tloop
      DEALLOCATE(FPJAC,FJAC,U,V,STAT=i)
      IF(i.NE.0) STOP 'deallocation error for FPJAC,FJAC, U, V'
      IF(GetQuad) THEN
        DEALLOCATE(WY,WP,STAT=i)
        IF(i .NE. 0) STOP 'deallocation error for WY,WP'
      END IF
      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE ERK_DadjInt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_AllocBuffers
!~~~>  Allocate buffer space for checkpointing
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       INTEGER :: i
   
       ALLOCATE( chk_H(Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer H'; STOP
       END IF   
       ALLOCATE( chk_T(Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer T'; STOP
       END IF   
       ALLOCATE( chk_Y(N,Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer Y'; STOP
       END IF   
       ALLOCATE( chk_Z(N,rkS,Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer K'; STOP
       END IF    
 
     END SUBROUTINE ERK_AllocBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE ERK_FreeBuffers
!~~~>  Dallocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       INTEGER :: i
   
       DEALLOCATE( chk_H, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer H'; STOP
       END IF   
       DEALLOCATE( chk_T, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer T'; STOP
       END IF   
       DEALLOCATE( chk_Y, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer Y'; STOP
       END IF   
       DEALLOCATE( chk_Z, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer K'; STOP
       END IF    
 
     END SUBROUTINE ERK_FreeBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE ERK_Push( T, H, Y, Z)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       DOUBLE PRECISION :: T, H, Y(N), Z(N,Smax)
!       INTEGER       :: P(N)
!       DOUBLE PRECISION :: E(N,N)
   
       stack_ptr = stack_ptr + 1
       IF ( stack_ptr > Max_no_steps ) THEN
         PRINT*,'Push failed: buffer overflow'
         STOP
       END IF  
       chk_H( stack_ptr ) = H
       chk_T( stack_ptr ) = T
       chk_Y(1:N,stack_ptr) = Y(1:N)
       chk_Z(1:N,1:rkS,stack_ptr) = Z(1:N,1:rkS)
  
      END SUBROUTINE ERK_Push

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_Pop( T, H, Y, Z)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       DOUBLE PRECISION :: T, H, Y(N), Z(N,Smax)
!       INTEGER       :: P(N)
!       DOUBLE PRECISION :: E(N,N)
   
       IF ( stack_ptr <= 0 ) THEN
         PRINT*,'Pop failed: empty buffer'
         STOP
       END IF  
       H = chk_H( stack_ptr )
       T = chk_T( stack_ptr )
       Y(1:N) = chk_Y(1:N,stack_ptr)
       Z(1:N,1:rkS) = chk_Z(1:N,1:rkS,stack_ptr)

       stack_ptr = stack_ptr - 1
  
      END SUBROUTINE ERK_Pop

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER :: i, ITOL
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
      END SUBROUTINE ERK_ErrorScale
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_ErrorNorm(N, Y, SCAL, Rerr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      
      INTEGER :: i, N
      DOUBLE PRECISION :: Y(N), SCAL(N), Rerr
      Rerr = ZERO
      DO i=1,N
           Rerr = Rerr+(Y(i)*SCAL(i))**2
      END DO
      Rerr = MAX( SQRT(Rerr/DBLE(N)), 1.0d-10 )
!
      END SUBROUTINE ERK_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      SUBROUTINE ERK_ErrorMsg(Code,T,H,Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
          PRINT * , '--> No of steps exceeds maximum bound', max_no_steps
        CASE (-7)
          PRINT * , '--> Step size too small: T + 10*H = T', &
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


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   END SUBROUTINE ERKADJ1 ! AND ALL ITS INTERNAL PROCEDURES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERKADJ2( N, NP, NADJ, Tinitial, Tfinal, Y, Lambda, FUN, JAC,    &
                       AdjInit, RelTol, AbsTol, RCNTRL, ICNTRL, RSTATUS,      &
                       ISTATUS, Ierr, DRDY, QFUN, Q)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using an Explicit
!    Runge-Kutta (ERK) method.
!
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
!        For ICNTRL(4)=0 the default value of 1500 is used
!        Note: use a conservative estimate, since the checkpoint
!              buffers are allocated to hold Max_no_steps
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
      IMPLICIT NONE

! Arguments      
      INTEGER, INTENT(IN)          :: N, NP, NADJ, ICNTRL(20)
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N), Lambda(N,NADJ)
      DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: Q(NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tinitial, Tfinal, &
                    RelTol(N), AbsTol(N), RCNTRL(20)
      INTEGER, INTENT(OUT)         :: Ierr
      INTEGER, INTENT(INOUT)       :: ISTATUS(20) 
      DOUBLE PRECISION, INTENT(OUT)   :: RSTATUS(20)
       
!~~~>  ERK method coefficients, up to 5 stages
      INTEGER, PARAMETER :: Smax = 12
      INTEGER, PARAMETER :: RK2=1, RK3=2, RK4=3, RK5=4, RK6=5, RK8=6
      DOUBLE PRECISION :: rkA(Smax,Smax), rkB(Smax), rkC(Smax),rkE(Smax),rkELO
      INTEGER :: erkMethod, rkS ! The number of stages

!~~~>  Checkpoints in memory buffers
      INTEGER :: stack_ptr = 0 ! last written entry in checkpoint
      DOUBLE PRECISION, DIMENSION(:),     POINTER :: chk_H, chk_T
      DOUBLE PRECISION, DIMENSION(:,:),   POINTER :: chk_Y
      DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: chk_Z
! Local variables      
      DOUBLE PRECISION :: Hmin, Hmax, Hstart, Roundoff,    &
                       FacMin, Facmax, FacSafe, FacRej, &
                       Qmin, Qmax

      INTEGER       :: ITOL, Max_no_steps, i
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
      DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
      LOGICAL  :: GetQuad
      EXTERNAL       QFUN,FUN,JAC,DRDY,AdjInit
      OPTIONAL       QFUN,DRDY
      DOUBLE PRECISION :: DLAMCH

      rkA = ZERO
      rkB = ZERO
      rkC = ZERO
      rkE = ZERO
      GetQuad = .FALSE.
      IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDY)) GetQuad = .TRUE.

!~~~>  Initialize statistics
      ISTATUS(1:20) = 0
      RSTATUS(1:20) = ZERO
      Ierr          = 0

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
         Max_no_steps = ICNTRL(4)
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
         DO i=1,N
            IF (AbsTol(i) <= 0.D0.OR.RelTol(i) <= 10.D0*Roundoff) THEN
              PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
              PRINT * , ' RelTol(',i,') = ',RelTol(i)
              CALL ERK_ErrorMsg(-5,Tinitial,ZERO,Ierr)
            END IF
         END DO
      END IF
    
    IF (Ierr < 0) RETURN
!~~~>  Allocate memory buffers    
    CALL ERK_AllocBuffers

!~~~>  Call forward integration    
    CALL ERK_FwdInt( N, NADJ, Tinitial, Tfinal, Y, GetQuad, Ierr )

    PRINT*,'FORWARD STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

!~~~> Initialize adjoint variables
    CALL AdjInit(N, NP, NADJ, Tfinal, Y, Lambda)

!~~~>  Call adjoint integration    
    CALL ERK_DadjInt( N, NADJ, Lambda, GetQuad, Ierr )

    PRINT*,'ADJOINT STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)
!~~~>  Free memory buffers    
    CALL ERK_FreeBuffers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CONTAINS  !  Procedures internal to ERKADJ2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERK_FwdInt( NVAR,NADJ,Tinitial,Tfinal,Y,GetQuad,Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE

!~~~> Arguments:      
      INTEGER,INTENT(IN) :: NVAR,NADJ
      DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)
      DOUBLE PRECISION, INTENT(IN) :: Tinitial, Tfinal
      INTEGER, INTENT(OUT) :: Ierr
      LOGICAL :: GetQuad
!~~~> Local variables:      
      DOUBLE PRECISION :: Z(NVAR,rkS),K(NVAR,rkS), TMP(NVAR), SCAL(NVAR), &
                          G(NVAR), R(NADJ), T, H, Hratio, Rerr, Fac, Hnew, Tdirection
      INTEGER :: j, istage
      LOGICAL :: Reject, FirstStep
      
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

      CALL ERK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)

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

       CALL Set2zero(NVAR,Z(1,istage))

       IF(istage > 1) THEN
           DO j = 1, istage-1
               CALL DAXPY(NVAR,H*rkA(istage,j),K(1,j),1,Z(1,istage),1)
           END DO
       END IF

       CALL DAXPY(NVAR,ONE,Y,1,Z(1,istage),1)

       CALL FUN(NVAR,T+rkC(istage)*H,Z(1,istage),K(1,istage))

       ISTATUS(Nfun) = ISTATUS(Nfun) + 1

  END DO stages

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Error estimation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      CALL SET2ZERO(N,TMP)
      DO i = 1,rkS
         IF (rkE(i)/=ZERO) CALL DAXPY(N,H*rkE(i),K(1,i),1,TMP,1)
      END DO  

      CALL ERK_ErrorNorm(N, TMP, SCAL, Rerr)

!~~~> Computation of new step size Hnew
      Fac  = FacSafe*(Rerr)**(-ONE/rkELO)
      Fac  = MAX(FacMin,MIN(FacMax,Fac))
      Hnew = H*Fac

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Accept/Reject step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
accept: IF ( Rerr < ONE ) THEN !~~~> Step is accepted

         FirstStep=.FALSE.
         ISTATUS(Nacc) = ISTATUS(Nacc) + 1

!~~~> Update the results for the quadrature term
         IF(GetQuad) Then
           CALL SET2ZERO(N,TMP)
           CALL SET2ZERO(NADJ,R)
           DO i = 1, rkS
             CALL WADD(N,Y,Z(1,i),TMP)
             CALL QFUN(N,NADJ,T+rkC(i)*H,TMP,R)
             CALL DAXPY(NADJ,H*rkB(i),R,1,Q,1)
           END DO
         END IF

!~~~> Checkpoint solution
         CALL ERK_Push( T, H, Y, Z )

!~~~> Update time and solution
         T  =  T + H
         ! Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
         DO i = 1,rkS 
            IF (rkB(i)/=ZERO) CALL DAXPY(N,H*rkB(i),K(1,i),1,Y,1)
         END DO  
       
!~~~> Update scaling coefficients
         CALL ERK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)

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
         
      END IF accept
      END DO Tloop

      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE ERK_FwdInt



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE ERK_DadjInt( N, NADJ, Lambda, GetQuad, Ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE

!~~~> Arguments:      
      INTEGER, INTENT(IN) :: N, NADJ
      DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ)
      INTEGER, INTENT(OUT) :: Ierr
      LOGICAL :: GetQuad
!~~~> Local variables:      
      DOUBLE PRECISION :: Alpha,Beta
      DOUBLE PRECISION :: Y(N)
      DOUBLE PRECISION :: Z(N,rkS), TMP(N), T, H
      DOUBLE PRECISION, ALLOCATABLE :: FJAC(:,:), U(:,:,:), WY(:,:)
      INTEGER :: j, istage, iadj

!      initiliazation

      IF(GetQuad) THEN
        ALLOCATE(WY(N,NADJ),STAT=j)
        IF(j .NE. 0) STOP 'allocation error for WY'
        WY(:,:) = 0.0d0
      END IF

      ALLOCATE(FJAC(N,N),U(N,NADJ,rkS),STAT=j)
      IF(j .NE. 0) STOP 'allocation error for FJAC,U'
      Alpha = 1
      Beta = 0     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( stack_ptr > 0 )
        
   !~~~>  Recover checkpoints for stage values and vectors
      CALL ERK_Pop( T, H, Y, Z)

      ISTATUS(Nstp) = ISTATUS(Nstp) + 1

stages:DO istage = rkS, 1, -1

!~~~>  Jacobian of the current stage solution
      CALL Jac(N,T+rkC(istage)*H,Z(1,istage),FJAC)
      ISTATUS(Njac) = ISTATUS(Njac) + 1

      IF(GetQuad) CALL DRDY(NADJ,N,N,T+rkC(istage)*H,Z(1,istage),WY)
       
adj:   DO iadj = 1, NADJ

        CALL SET2ZERO(N,TMP)
        IF (istage < rkS) THEN
          DO j = istage+1, rkS
            CALL DAXPY(N,H*rkA(j,istage),U(1,iadj,j),1,TMP,1) ! TMP < h*A*U
          END DO
        END IF
        CALL DAXPY(N,H*rkB(istage),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B*Lambda

!~~~~>f_y^T * (b_i * Lambda_{n+1} + \sum_{j=i+1}^s a_{j,i} u_j)

        CALL DGEMV('T',N,N,ALPHA,FJAC,N,TMP,1,BETA,U(1,iadj,istage),1)  ! U = Jac^T * TMP

        IF(GetQuad) CALL DAXPY(N,H*rkB(istage),WY(1,iadj),1,U(1,iadj,istage),1)
! U = U + h*B*WY
        
        END DO adj
 
      END DO stages

!~~~> Update adjoint solution
         ! Y(:) <-- Y(:) + Sum_j rkD(j)*Z_j(:)
         DO istage = 1,rkS 
            DO iadj = 1,NADJ
                CALL DAXPY(N,ONE,U(1,iadj,istage),1,Lambda(1,iadj),1)
                !Lambda(1:N,iadj) = Lambda(1:N,iadj) + U(1:N,iadj,istage)
            END DO
         END DO  

      END DO Tloop

      DEALLOCATE(FJAC,U,STAT=j)
      IF(j .NE. 0) STOP 'deallocate error for FJAC or U'
      IF(GetQuad) THEN
        DEALLOCATE(WY,STAT=i)
        IF(i .NE. 0) STOP 'deallocation error for WY'
      END IF

      ! Successful return
      Ierr  = 1
  
      END SUBROUTINE ERK_DadjInt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_AllocBuffers
!~~~>  Allocate buffer space for checkpointing
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       INTEGER :: i
   
       ALLOCATE( chk_H(Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer H'; STOP
       END IF   
       ALLOCATE( chk_T(Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer T'; STOP
       END IF   
       ALLOCATE( chk_Y(N,Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer Y'; STOP
       END IF   
       ALLOCATE( chk_Z(N,rkS,Max_no_steps), STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed allocation of buffer K'; STOP
       END IF     
 
     END SUBROUTINE ERK_AllocBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE ERK_FreeBuffers
!~~~>  Dallocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       INTEGER :: i
   
       DEALLOCATE( chk_H, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer H'; STOP
       END IF   
       DEALLOCATE( chk_T, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer T'; STOP
       END IF   
       DEALLOCATE( chk_Y, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer Y'; STOP
       END IF   
       DEALLOCATE( chk_Z, STAT=i )
       IF (i/=0) THEN
          PRINT*,'Failed deallocation of buffer K'; STOP
       END IF   
 
     END SUBROUTINE ERK_FreeBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE ERK_Push( T, H, Y, Z)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       DOUBLE PRECISION :: T, H, Y(N), Z(N,Smax)
!       INTEGER       :: P(N)
!       DOUBLE PRECISION :: E(N,N)
   
       stack_ptr = stack_ptr + 1
       IF ( stack_ptr > Max_no_steps ) THEN
         PRINT*,'Push failed: buffer overflow'
         STOP
       END IF  
       chk_H( stack_ptr ) = H
       chk_T( stack_ptr ) = T
       chk_Y(1:N,stack_ptr) = Y(1:N)
       chk_Z(1:N,1:rkS,stack_ptr) = Z(1:N,1:rkS)
  
      END SUBROUTINE ERK_Push

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_Pop( T, H, Y, Z)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       DOUBLE PRECISION :: T, H, Y(N), Z(N,Smax)
!       INTEGER       :: P(N)
!       DOUBLE PRECISION :: E(N,N)
   
       IF ( stack_ptr <= 0 ) THEN
         PRINT*,'Pop failed: empty buffer'
         STOP
       END IF  
       H = chk_H( stack_ptr )
       T = chk_T( stack_ptr )
       Y(1:N) = chk_Y(1:N,stack_ptr)
       Z(1:N,1:rkS) = chk_Z(1:N,1:rkS,stack_ptr)

       stack_ptr = stack_ptr - 1
  
      END SUBROUTINE ERK_Pop

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_ErrorScale(ITOL, AbsTol, RelTol, Y, SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER :: i, ITOL
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
      END SUBROUTINE ERK_ErrorScale
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE ERK_ErrorNorm(N, Y, SCAL, Rerr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      
      INTEGER :: i, N
      DOUBLE PRECISION :: Y(N), SCAL(N), Rerr
      Rerr = ZERO
      DO i=1,N
           Rerr = Rerr+(Y(i)*SCAL(i))**2
      END DO
      Rerr = MAX( SQRT(Rerr/DBLE(N)), 1.0d-10 )
!
      END SUBROUTINE ERK_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      SUBROUTINE ERK_ErrorMsg(Code,T,H,Ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
          PRINT * , '--> No of steps exceeds maximum bound', max_no_steps
        CASE (-7)
          PRINT * , '--> Step size too small: T + 10*H = T', &
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
      rkELO = 3
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
      rkELO = 3.d0
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
      rkELO = 5
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
      rkELO = 6
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
      rkELO = 7
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
      rkELO = 9
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


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   END SUBROUTINE ERKADJ2 ! AND ALL ITS INTERNAL PROCEDURES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
END MODULE ERK_adj_f90_Integrator


! End of INTEGRATE function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


