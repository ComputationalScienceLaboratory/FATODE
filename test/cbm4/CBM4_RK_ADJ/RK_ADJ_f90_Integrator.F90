!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!RungeKuttaADJ - Adjoint Model of Fully Implicit 3-stage Runge-Kutta:         !
!          * Radau-2A   quadrature (order 5)                                  !
!          * Radau-1A   quadrature (order 5)                                  !
!          * Lobatto-3C quadrature (order 4)                                  !
!          * Gauss      quadrature (order 6)                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE RK_adj_f90_Integrator
  IMPLICIT NONE
  PUBLIC
  SAVE

!~~~>  Statistics on the work performed by the Runge-Kutta method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, Nrej=5, Ndec=6,       &
                        Nsol=7, Nsng=8, Ntexit=1, Nhacc=2, Nhnew=3  
CONTAINS

!******************************************************************************
  SUBROUTINE INTEGRATE_ADJ(NVAR, NP, NADJ, NNZERO, Y, Lambda, TIN, TOUT,  &
               ATOL_adj, RTOL_adj, ATOL, RTOL, FUN, JAC, AdjInit, ICNTRL_U, &
               RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, Mu, JACP, DRDY, DRDP, &
               QFUN, Q )
    IMPLICIT NONE
    INTEGER :: NVAR,NNZERO,NP

!~~~> Y - ODE state vector, initial conditions as input, final solution as
!output
    DOUBLE PRECISION, INTENT(INOUT)  :: Y(NVAR)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL  :: Q(NADJ)

!~~~> NADJ - No. of cost functionals for which adjoints are evaluated simultaneously. 
!   If single cost functional is considered (like in most applications) simply set NADJ = 1      
    INTEGER, INTENT(IN) :: NADJ

!~~~> Lambda - Sensitivities
!     Note: Lambda (1:NVAR,j) contains sensitivities of the j-th cost functional w.r.t. Y(1:NVAR), j=1...NADJ
!           Mu(1:NP,j) contains sensitivities of the j-th cost functional w.r.t parameters 
    DOUBLE PRECISION, INTENT(INOUT) :: Lambda(NVAR,NADJ)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL ::  Mu(NP,NADJ)

!~~~> Tolerances for adjoint calculations
!   (used for full continuous adjoint, and for controlling the convergence of iterations for solving the discrete adjoint)
    DOUBLE PRECISION, INTENT(IN) :: ATOL_adj(NVAR,NADJ), RTOL_adj(NVAR,NADJ), &
                                    ATOL(NVAR), RTOL(NVAR)
    DOUBLE PRECISION :: TIN  ! TIN - Start Time
    DOUBLE PRECISION :: TOUT ! TOUT - End Time
    INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
    INTEGER,       INTENT(INOUT), OPTIONAL :: ISTATUS_U(20)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: RSTATUS_U(20)
    INTEGER,       INTENT(INOUT), OPTIONAL :: IERR_U
    INTEGER :: IERR
    DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
    INTEGER       :: ICNTRL(20), ISTATUS(20)
    EXTERNAL :: FUN,JAC,AdjInit
    EXTERNAL :: QFUN,DRDY,DRDP,JACP
    OPTIONAL :: QFUN,DRDY,DRDP,JACP
    ICNTRL(1:20) = 0
    RCNTRL(1:20) = 0.0

!~~~> fine-tune the integrator:
    ICNTRL(2)  = 0   ! 0=vector tolerances, 1=scalar tolerances
    ICNTRL(5)  = 8   ! Max no. of Newton iterations
    ICNTRL(6)  = 0   ! Starting values for Newton are: 0=interpolated, 1=zero
    ICNTRL(7)  = 1   ! Adj. system solved by: 1=iteration, 2=direct, 3=adaptive
    ICNTRL(8)  = 0   ! Adj. LU decomp: 0=compute, 1=save from fwd
    ICNTRL(9)  = 2   ! Adjoint: 1=none, 2=discrete, 3=full continuous, 4=simplified continuous
    ICNTRL(10) = 1   ! Error estimator: 0=classic, 1=SDIRK ???
    ICNTRL(11) = 1   ! Step controller: 1=Gustaffson, 2=classic

!~~~> if optional parameters are given, and if they are >0, then use them to overwrite default settings
    IF (PRESENT(ICNTRL_U)) THEN
      WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
    END IF
    IF (PRESENT(RCNTRL_U)) THEN
      WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
    END IF
    T1 = TIN; T2 = TOUT   

! Evaluating sensitivities w.r.t parameters requires NP>0, functions MU and JACP are provided.
    IF(NP>0 .AND. PRESENT(Mu) .AND. PRESENT(JACP) ) THEN
      IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND.          &
        PRESENT(DRDY) .AND. PRESENT(DRDP)) THEN
        print *,'mode1'
!~~~> This is the case that cost function contains a quadrature term NADJ 
!     should be 1; Q is defined to store the value of the quadrature
!     term at the last step and functions QFUN,DRDY,DRDP must be provided.
!     cost function = g(y,p) + Integrate{r(y,p)}
        CALL RungeKuttaADJ1( N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, Q=Q,          &
               NADJ=NADJ, Lambda=Lambda, Mu=Mu, Tstart=T1, Tend=T2,           &
               RelTol=RTOL, AbsTol=ATOL, RelTol_adj=RTOL_adj,                 &
               AbsTol_adj=ATOL_adj, JACP=JACP, DRDY=DRDY, DRDP=DRDP,          &
               AdjInit=AdjInit, RCNTRL=RCNTRL, ICNTRL=ICNTRL, RSTATUS=RSTATUS,&
               ISTATUS=ISTATUS, IERR=IERR, FUN=FUN, QFUN=QFUN, JAC=JAC ) 
      ELSE

!~~~> No quadrature term is involved. cost function = g(y,p)
        CALL RungeKuttaADJ1( N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, NADJ=NADJ,    &
               Lambda=Lambda, Mu=Mu, Tstart=T1, Tend=T2, RelTol=RTOL,         &
               AbsTol=Atol, RelTol_adj=RTOL_adj, AbsTol_adj=ATOL_adj,         &
               JACP=JACP, AdjInit=AdjInit, RCNTRL=RCNTRL, ICNTRL=ICNTRL,      &
               RSTATUS=RSTATUS, ISTATUS=ISTATUS, IERR=IERR, FUN=FUN, JAC=JAC )
      END IF
        
    ELSE

! Evaluating sensitivites w.r.t only initial conditions
      IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDY))THEN
!~~~> This is the case that cost function contains a quadrature term
! NADJ should be 1; Q is defined to store the value of the quadrature
! term at the last step and functions QFUN,DRDY must be provided. 
! cost function = g(y) + Integrate{r(y)}            
        CALL RungeKuttaADJ2( N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, Q=Q,          &
               NADJ=NADJ, Lambda=Lambda, Tstart=T1, Tend=T2, RelTol=RTOL,     &
               AbsTol=ATOL, RelTol_adj=RTOL_adj, AbsTol_adj=ATOL_adj,         &
               AdjInit=AdjInit, RCNTRL=RCNTRL, ICNTRL=ICNTRL, RSTATUS=RSTATUS,&
               ISTATUS=ISTATUS, IERR=IERR, FUN=FUN, JAC=JAC, QFUN=QFUN,       &
               DRDY=DRDY )
      ELSE

          !~~~> No quadrature term is involved. cost function = g(y)        
        CALL RungeKuttaADJ2( N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, NADJ=NADJ,    &
               Lambda=Lambda, Tstart=T1, Tend=T2, RelTol=RTOL, AbsTol=ATOL,   &
               RelTol_adj=RTOL_adj, AbsTol_adj=ATOL_adj, AdjInit=AdjInit,     &
               RCNTRL=RCNTRL, ICNTRL=ICNTRL, RSTATUS=RSTATUS, ISTATUS=ISTATUS,&
               IERR=IERR, FUN=FUN, JAC=JAC )
      END IF          
    END IF

    ! if optional parameters are given for output use them to store information in them
    IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
    IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
    IF (PRESENT(IERR_U)) IERR_U = IERR

    IF (IERR < 0) THEN
      PRINT *,'Runge-Kutta-ADJ: Unsuccessful exit at T=',TIN,' (IERR=',IERR,')'
    ENDIF

  END SUBROUTINE INTEGRATE_ADJ


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE RungeKuttaADJ1( N, NP, NNZERO, Y, NADJ, Lambda, Mu, Tstart, Tend,&
               RelTol, AbsTol, RelTol_adj, AbsTol_adj, FUN, JAC, JACP,        &
               AdjInit, RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR, DRDY, DRDP,   &
               QFUN, Q )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  This implementation is based on the book and the code Radau5:
!
!         E. HAIRER AND G. WANNER
!         "SOLVING ORDINARY DIFFERENTIAL EQUATIONS II. 
!              STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS."
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
!         SPRINGER-VERLAG (1991)
!
!         UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!         CH-1211 GENEVE 24, SWITZERLAND
!         E-MAIL:  HAIRER@DIVSUN.UNIGE.CH,  WANNER@DIVSUN.UNIGE.CH
!
!   Methods:
!          * Radau-2A   quadrature (order 5)                              
!          * Radau-1A   quadrature (order 5)                              
!          * Lobatto-3C quadrature (order 4)                              
!          * Gauss      quadrature (order 6)                              
!                                                                         
!   (C)  Adrian Sandu, August 2005                                       
!   Virginia Polytechnic Institute and State University                  
!   Contact: sandu@cs.vt.edu                                             
!   Revised by Philipp Miehe and Adrian Sandu, May 2006                  
!   Revised by Hong Zhang and Adrian Sandu, Feb 2011
!   This implementation is part of FATODE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!       ----------------
!
!    Note: For input parameters equal to zero the default values of the
!          corresponding variables are used.
!
!     N           Dimension of the system
!     NP          Number of parameters of interest
!     NADJ        Number of adjoints
!     T           Initial time value
!
!     Tend        Final T value (Tend-T may be positive or negative)
!
!     Y(N)        Initial values for Y
!     
!     Q(NADJ)     Initial values for Q
!     RelTol,AbsTol   Relative and absolute error tolerances. 
!          for ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
!                        = 1: AbsTol, RelTol are scalars
!
!     Lambda,Mu   Initial values for adjoint variables
!~~~>  Integer input parameters:
!  
!    ICNTRL(1) = not used
!
!    ICNTRL(2) = 0: AbsTol, RelTol are NVAR-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3) = RK method selection       
!              = 1:  Radau-2A    (the default)
!              = 2:  Lobatto-3C
!              = 3:  Gauss
!              = 4:  Radau-1A
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0 the default value of 10000 is used
!
!    ICNTRL(5)  -> maximum number of Newton iterations
!        For ICNTRL(5)=0 the default value of 8 is used
!
!    ICNTRL(6)  -> starting values of Newton iterations:
!        ICNTRL(6)=0 : starting values are obtained from 
!                      the extrapolated collocation solution
!                      (the default)
!        ICNTRL(6)=1 : starting values are zero
!
!    ICNTRL(7)  -> method to solve the linear ADJ equations:
!        ICNTRL(7)=0,1 : modified Newton re-using LU (the default)
!                        with a fixed number of iterations
!        ICNTRL(7)=2 :   direct solution (additional one LU factorization
!                        of 3Nx3N matrix per step); good for debugging
!        ICNTRL(7)=3 :   adaptive solution (if Newton does not converge
!                        switch to direct)
!
!    ICNTRL(8)  -> checkpointing the LU factorization at each step:
!        ICNTRL(8)=0 : do *not* save LU factorization (the default)
!        ICNTRL(8)=1 : save LU factorization
!        Note: if ICNTRL(7)=1 the LU factorization is *not* saved
!
!    ICNTRL(9) -> Type of adjoint algorithm
!         = 0 : default is discrete adjoint ( of method ICNTRL(3) )
!         = 1 : no adjoint       
!         = 2 : discrete adjoint ( of method ICNTRL(3) )
!         = 3 : fully adaptive continuous adjoint ( with method ICNTRL(6) )
!         = 4 : simplified continuous adjoint ( with method ICNTRL(6) )
!
!    ICNTRL(10) -> switch for error estimation strategy
!		ICNTRL(10) = 0: one additional stage at c=0, 
!				see Hairer (default)
!		ICNTRL(10) = 1: two additional stages at c=0 
!				and SDIRK at c=1, stiffly accurate
!
!    ICNTRL(11) -> switch for step size strategy
!              ICNTRL(11)=1:  mod. predictive controller (Gustafsson, default)
!              ICNTRL(11)=2:  classical step size control
!              the choice 1 seems to produce safer results;
!              for simple problems, the choice 2 produces
!              often slightly faster runs
!
!~~~>  Real input parameters:
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!                  (highly recommended to keep Hmin = ZERO, the default)
!
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!
!    RCNTRL(3)  -> Hstart, the starting step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                 (default=0.1)
!
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!                  than the predicted value  (default=0.9)
!
!    RCNTRL(8)  -> ThetaMin. If Newton convergence rate smaller
!                  than ThetaMin the Jacobian is not recomputed;
!                  (default=0.001)
!
!    RCNTRL(9)  -> NewtonTol, stopping criterion for Newton's method
!                  (default=0.03)
!
!    RCNTRL(10) -> Qmin
!
!    RCNTRL(11) -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
!                  step size is kept constant and the LU factorization
!                  reused (default Qmin=1, Qmax=1.2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!    OUTPUT ARGUMENTS:
!    -----------------
!
!    T           -> T value for which the solution has been computed
!                     (after successful return T=Tend).
!
!    Y(N)        -> Numerical solution at T
!
!    Q(NADJ)     -> Numerical solution at T
!
!    IERR        -> Reports on successfulness upon return:
!                    = 1 for success
!                    < 0 for error (value equals error code)
!
!    ISTATUS(1)  -> No. of function calls
!    ISTATUS(2)  -> No. of Jacobian calls
!    ISTATUS(3)  -> No. of steps
!    ISTATUS(4)  -> No. of accepted steps
!    ISTATUS(5)  -> No. of rejected steps (except at very beginning)
!    ISTATUS(6)  -> No. of LU decompositions
!    ISTATUS(7)  -> No. of forward/backward substitutions
!    ISTATUS(8)  -> No. of singular matrix decompositions
!
!    RS, the time corresponding to the
!                     computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!                   For multiple restarts, use Hnew as Hstart 
!                     in the subsequent run
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USE LS_Solver
    IMPLICIT NONE      
    INTEGER :: N,NP,NNZERO
    INTEGER, INTENT(IN)     :: NADJ
    INTEGER, INTENT(INOUT)  :: IERR
    DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ), Mu(NP,NADJ)
    DOUBLE PRECISION, INTENT(IN) :: AbsTol(N), RelTol(N), AbsTol_adj(N,NADJ), &
                                    RelTol_adj(N,NADJ)
    DOUBLE PRECISION, INTENT(INOUT) :: Y(N),RCNTRL(20),RSTATUS(20)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: Q(NADJ)
    INTEGER, INTENT(IN) :: ICNTRL(20)
    INTEGER, INTENT(INOUT) :: ISTATUS(20)
    LOGICAL :: StartNewton, Gustafsson, SdirkError, SaveLU, GetQuad
    INTEGER :: ITOL
    DOUBLE PRECISION, INTENT(IN):: Tstart,Tend
    DOUBLE PRECISION :: Texit

!~~~> Control arguments
    INTEGER :: Max_no_steps, NewtonMaxit, AdjointType, rkMethod
    DOUBLE PRECISION :: Hmin,Hmax,Hstart,Qmin,Qmax
    DOUBLE PRECISION :: Roundoff, ThetaMin, NewtonTol
    DOUBLE PRECISION :: FacSafe,FacMin,FacMax,FacRej

! Runge-Kutta method parameters
    INTEGER, PARAMETER :: R2A=1, R1A=2, L3C=3, GAU=4, L3A=5
    DOUBLE PRECISION :: rkT(3,3), rkTinv(3,3), rkTinvAinv(3,3), rkAinvT(3,3), &
                        rkA(3,3), rkB(3), rkC(3), rkD(0:3), rkE(0:3),         &
                        rkBgam(0:4), rkBhat(0:4), rkTheta(0:3), rkGamma,      &
                        rkAlpha, rkBeta, rkELO

! ADJ method parameters
    INTEGER, PARAMETER :: RK_adj_f90_none = 1, RK_adj_f90_discrete = 2,       &
                          RK_adj_f90_continuous = 3,                          &
                          RK_adj_f90_simple_continuous = 4
    INTEGER :: AdjointSolve                      
    INTEGER, PARAMETER :: Solve_direct = 1, Solve_fixed = 2, Solve_adaptive = 3
    INTEGER :: stack_ptr = 0 ! last written entry
    DOUBLE PRECISION, DIMENSION(:), POINTER :: chk_H, chk_T
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_Y, chk_Z
    INTEGER, DIMENSION(:), POINTER :: chk_NiT
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_dY, chk_d2Y

!~~~> Local variables
    INTEGER :: i 
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
    DOUBLE PRECISION :: DLAMCH
    EXTERNAL QFUN,FUN,JAC,JACP,DRDY,DRDP,AdjInit
    OPTIONAL QFUN,DRDY,DRDP

    GetQuad = .FALSE.
!~~~> Compute the quadrature term if required sourcs are provided
    IF( PRESENT(Q) .AND. PRESENT(QFUN) .AND. &
      PRESENT(DRDP) .AND. PRESENT(DRDY) )THEN
        GetQuad = .TRUE.
    END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        SETTING THE PARAMETERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IERR = 0
    ISTATUS(1:20) = 0
    RSTATUS(1:20) = ZERO
       
!~~~> ICNTRL(1) - autonomous system - not used       
!~~~> ITOL: 1 for vector and 0 for scalar AbsTol/RelTol
    IF (ICNTRL(2) == 0) THEN
      ITOL = 1
    ELSE
      ITOL = 0
    END IF

!~~~> Error control selection  
    IF (ICNTRL(10) == 0) THEN 
      SdirkError = .FALSE.
    ELSE
      SdirkError = .TRUE.
    END IF      

!~~~> Method selection  
    SELECT CASE (ICNTRL(3))     
    CASE (1)
      CALL Radau2A_Coefficients
    CASE (0,2)
      CALL Lobatto3C_Coefficients
    CASE (3)
      CALL Gauss_Coefficients
    CASE (4)
      CALL Radau1A_Coefficients
    CASE DEFAULT
      WRITE(6,*) 'ICNTRL(3)=',ICNTRL(3)
      CALL RK_ErrorMsg(-13,Tstart,ZERO,IERR)
    END SELECT

!~~~> Max_no_steps: the maximal number of time steps
    IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 10000
    ELSE
      Max_no_steps=ICNTRL(4)
      IF (Max_no_steps <= 0) THEN
        WRITE(6,*) 'ICNTRL(4)=',ICNTRL(4)
        CALL RK_ErrorMsg(-1,Tstart,ZERO,IERR)
      END IF
    END IF

!~~~> NewtonMaxit    maximal number of Newton iterations
    IF (ICNTRL(5) == 0) THEN
      NewtonMaxit = 8
    ELSE
      NewtonMaxit=ICNTRL(5)
      IF (NewtonMaxit <= 0) THEN
        WRITE(6,*) 'ICNTRL(5)=',ICNTRL(5)
        CALL RK_ErrorMsg(-2,Tstart,ZERO,IERR)
      END IF
    END IF

!~~~> StartNewton:  Use extrapolation for starting values of Newton iterations
    IF (ICNTRL(6) == 0) THEN
      StartNewton = .TRUE.
    ELSE
      StartNewton = .FALSE.
    END IF      

!~~~>  How to solve the linear adjoint system
    SELECT CASE (ICNTRL(7))     
    CASE (0,1)
      AdjointSolve = Solve_fixed
    CASE (2)
      AdjointSolve = Solve_direct
    CASE (3)
      AdjointSolve = Solve_adaptive
    CASE DEFAULT  
      PRINT * , 'User-selected adjoint solution: ICNTRL(7)=', ICNTRL(7)
      PRINT * , 'Implemented: =(0,1) (fixed), =2 (direct), =3 (adaptive)'
      CALL rk_ErrorMsg(-9,Tstart,ZERO,IERR)
      RETURN      
    END SELECT

!~~~>  Discrete or continuous adjoint formulation
    SELECT CASE (ICNTRL(9))     
    CASE (0,2)
      AdjointType = RK_adj_f90_discrete
    CASE (1)
      AdjointType = RK_adj_f90_none
    CASE DEFAULT  
      PRINT * , 'User-selected adjoint type: ICNTRL(9)=', ICNTRL(9)
      PRINT * , 'Implemented: =(0,2) (discrete adj) and =1 (no adj)'
      CALL rk_ErrorMsg(-9,Tstart,ZERO,IERR)
      RETURN      
    END SELECT

!~~~> Save or not the forward LU factorization
    SaveLU = (ICNTRL(8) /= 0)
    SaveLU = .FALSE. 
    IF (AdjointSolve == Solve_direct) SaveLU = .FALSE.

!~~~> Gustafsson: step size controller
    IF (ICNTRL(11) == 0) THEN
      Gustafsson = .TRUE.
    ELSE
      Gustafsson = .FALSE.
    END IF

!~~~> Roundoff: smallest number s.t. 1.0 + Roundoff > 1.0
    Roundoff = DLAMCH('E');

!~~~> Hmin = minimal step size
    IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
    ELSE
      Hmin = MIN(ABS(RCNTRL(1)),ABS(Tend-Tstart))
    END IF

!~~~> Hmax = maximal step size
    IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
    ELSE
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
    END IF

!~~~> Hstart = starting step size
    IF (RCNTRL(3) == ZERO) THEN
      Hstart = ZERO
    ELSE
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
    END IF

!~~~> FacMin: lower bound on step decrease factor
    IF(RCNTRL(4) == ZERO)THEN
      FacMin = 0.2d0
    ELSE
      FacMin = RCNTRL(4)
    END IF

!~~~> FacMax: upper bound on step increase factor
    IF(RCNTRL(5) == ZERO)THEN
      FacMax = 8.D0
    ELSE
      FacMax = RCNTRL(5)
    END IF

!~~~> FacRej: step decrease factor after 2 consecutive rejections
    IF(RCNTRL(6) == ZERO)THEN
      FacRej = 0.1d0
    ELSE
      FacRej = RCNTRL(6)
    END IF

!~~~> FacSafe:  by which the new step is slightly smaller
!               than the predicted value
    IF (RCNTRL(7) == ZERO) THEN
      FacSafe=0.9d0
    ELSE
      FacSafe=RCNTRL(7)
    END IF
    IF ( (FacMax < ONE) .OR. (FacMin > ONE) .OR. (FacSafe <= 1.0d-3)  &
         .OR. (FacSafe >= ONE) ) THEN
      WRITE(6,*)'RCNTRL(4:7)=',RCNTRL(4:7)
      CALL RK_ErrorMsg(-4,Tstart,ZERO,IERR)
    END IF

!~~~> ThetaMin:  decides whether the Jacobian should be recomputed
    IF (RCNTRL(8) == ZERO) THEN
      ThetaMin = 1.0d-3
    ELSE
      ThetaMin=RCNTRL(8)
      IF (ThetaMin <= 0.0d0 .OR. ThetaMin >= 1.0d0) THEN
        WRITE(6,*) 'RCNTRL(8)=', RCNTRL(8)
        CALL RK_ErrorMsg(-5,Tstart,ZERO,IERR)
      END IF
    END IF

!~~~> NewtonTol:  stopping crierion for Newton's method
    IF (RCNTRL(9) == ZERO) THEN
      NewtonTol = 3.0d-2
    ELSE
      NewtonTol = RCNTRL(9)
      IF (NewtonTol <= Roundoff) THEN
        WRITE(6,*) 'RCNTRL(9)=',RCNTRL(9)
        CALL RK_ErrorMsg(-6,Tstart,ZERO,IERR)
      END IF
    END IF

!~~~> Qmin AND Qmax: IF Qmin < Hnew/Hold < Qmax then step size = const.
    IF (RCNTRL(10) == ZERO) THEN
      Qmin=1.D0
    ELSE
      Qmin=RCNTRL(10)
    END IF
    IF (RCNTRL(11) == ZERO) THEN
      Qmax=1.2D0
    ELSE
      Qmax=RCNTRL(11)
    END IF
    IF (Qmin > ONE .OR. Qmax < ONE) THEN
      WRITE(6,*) 'RCNTRL(10:11)=',Qmin,Qmax
      CALL RK_ErrorMsg(-7,Tstart,ZERO,IERR)
    END IF

!~~~> Check if tolerances are reasonable
    IF (ITOL == 0) THEN
      IF (AbsTol(1) <= ZERO.OR.RelTol(1) <= 10.d0*Roundoff) THEN
        WRITE (6,*) 'AbsTol/RelTol=',AbsTol,RelTol 
        CALL RK_ErrorMsg(-8,Tstart,ZERO,IERR)
      END IF
    ELSE
      DO i=1,N
        IF (AbsTol(i) <= ZERO.OR.RelTol(i) <= 10.d0*Roundoff) THEN
          WRITE (6,*) 'AbsTol/RelTol(',i,')=',AbsTol(i),RelTol(i)
          CALL RK_ErrorMsg(-8,Tstart,ZERO,IERR)
        END IF
      END DO
    END IF

!~~~> Parameters are wrong
    IF (IERR < 0) RETURN

!~~~>  Allocate checkpoint space or open checkpoint files
    IF (AdjointType == RK_adj_f90_discrete) THEN
      CALL rk_AllocateDBuffers()
    ELSEIF ( (AdjointType == RK_adj_f90_continuous).OR.                       &
           (AdjointType == RK_adj_f90_simple_continuous) ) THEN
      CALL rk_AllocateCBuffers
    END IF

    CALL LSS_Init(N,NNZERO)

!~~~> Call the core method
    CALL RK_FwdIntegrator( N,NADJ,Tstart,Tend,Y,AdjointType,GetQuad,IERR )
    PRINT*,'FORWARD STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc), ' Rej=',ISTATUS(Nrej),&
           ' Singular=',Nsng,' Fun=', ISTATUS(Nfun),' Jac=',ISTATUS(Njac),    &
           ' Sol=', ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

!~~~>  If Forward integration failed return   
    IF (IERR<0) RETURN

!~~~> Initialize adjoint variables
    CALL AdjInit(N,NP,NADJ,Tend,Y,Lambda,Mu)

    SELECT CASE (AdjointType)   
    CASE (RK_adj_f90_discrete)   
      CALL rk_DadjInt( N,NADJ,NP,Lambda,Mu,Tstart,Tend,Texit,GetQuad,IERR )
    CASE (RK_adj_f90_continuous) 
      CALL rk_CadjInt( NADJ, Lambda, Tend, Tstart, Texit, IERR )
    CASE (RK_adj_f90_simple_continuous)
      CALL rk_SimpleCadjInt( NADJ, Lambda, Tstart, Tend, Texit, IERR )
    END SELECT ! AdjointType
    PRINT*,'ADJOINT STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc), ' Rej=',ISTATUS(Nrej),&
           ' Singular=', Nsng, ' Fun=',ISTATUS(Nfun), ' Jac=',ISTATUS(Njac),  &
           ' Sol=', ISTATUS(Nsol), ' Dec=', ISTATUS(Ndec)

!~~~>  Free checkpoint space or close checkpoint files
    IF (AdjointType == RK_adj_f90_discrete) THEN
      CALL rk_FreeDBuffers
    ELSEIF ( (AdjointType == RK_adj_f90_continuous) .OR. &
           (AdjointType == RK_adj_f90_simple_continuous) ) THEN
      CALL rk_FreeCBuffers
    END IF
    CALL LSS_Free
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS ! Internal procedures to RungeKuttaADJ1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_AllocateDBuffers()
!~~~>  Allocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
      ALLOCATE( chk_Z(N*3,Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer Z'; STOP
      END IF   
      ALLOCATE( chk_NiT(Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer NiT'; STOP
      END IF   
 
    END SUBROUTINE rk_AllocateDBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_FreeDBuffers()
!~~~>  Dallocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        PRINT*,'Failed deallocation of buffer Z'; STOP
      END IF   
      DEALLOCATE( chk_NiT, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer NiT'; STOP
      END IF   

    END SUBROUTINE rk_FreeDBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_AllocateCBuffers()
!~~~>  Allocate buffer space for continuous adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
      ALLOCATE( chk_dY(N,Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer dY'; STOP
      END IF   
      ALLOCATE( chk_d2Y(N,Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer d2Y'; STOP
      END IF   
 
    END SUBROUTINE rk_AllocateCBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_FreeCBuffers()
!~~~>  Dallocate buffer space for continuous adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: i
      print*,'cbuffers deallocate???'
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
      DEALLOCATE( chk_dY, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer dY'; STOP
      END IF   
      DEALLOCATE( chk_d2Y, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer d2Y'; STOP
      END IF   
 
    END SUBROUTINE rk_FreeCBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_DPush( T, H, Y, Zstage, NewIt)!, Jcb )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION :: T, H, Y(N), Zstage(N*3)
      INTEGER :: NewIt
      stack_ptr = stack_ptr + 1
      IF ( stack_ptr > Max_no_steps ) THEN
        PRINT*,'Push failed: buffer overflow'
        STOP
      END IF  
      chk_H( stack_ptr ) = H
      chk_T( stack_ptr ) = T
      ! CALL DCOPY(NVAR,Y,1,chk_Y(1,stack_ptr),1)
      ! CALL DCOPY(NVAR*3,Zstage,1,chk_Z(1,stack_ptr),1)
      chk_Y(1:N,stack_ptr) = Y(1:N)
      chk_Z(1:3*N,stack_ptr) = Zstage(1:3*N)
      chk_NiT( stack_ptr ) = NewIt
  
    END SUBROUTINE rk_DPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_DPop( T, H, Y, Zstage, NewIt) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION :: T, H, Y(N), Zstage(N*3) ! , Jcb(LU_NONZERO)
      INTEGER :: NewIt   
      IF ( stack_ptr <= 0 ) THEN
        PRINT*,'Pop failed: empty buffer'
        STOP
      END IF  
      H = chk_H( stack_ptr )
      T = chk_T( stack_ptr )
      ! CALL DCOPY(NVAR,chk_Y(1,stack_ptr),1,Y,1)
      Y(1:N) = chk_Y(1:N,stack_ptr)
      ! CALL DCOPY(NVAR*3,chk_Z(1,stack_ptr),1,Zstage,1)
      Zstage(1:3*N) = chk_Z(1:3*N,stack_ptr)
      NewIt = chk_NiT( stack_ptr )
      stack_ptr = stack_ptr - 1
  
    END SUBROUTINE rk_DPop
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_CPush(T, H, Y, dY, d2Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION :: T, H, Y(N), dY(N), d2Y(N)   
      stack_ptr = stack_ptr + 1
      IF ( stack_ptr > Max_no_steps ) THEN
        PRINT*,'Push failed: buffer overflow'
        STOP
      END IF  
      chk_H( stack_ptr ) = H
      chk_T( stack_ptr ) = T
      ! CALL DCOPY(NVAR,Y,1,chk_Y(1,stack_ptr),1)
      ! CALL DCOPY(NVAR,dY,1,chk_dY(1,stack_ptr),1)
      ! CALL DCOPY(NVAR,d2Y,1,chk_d2Y(1,stack_ptr),1)
      chk_Y(1:N,stack_ptr)  =  Y(1:N)
      chk_dY(1:N,stack_ptr) = dY(1:N)
      chk_d2Y(1:N,stack_ptr)= d2Y(1:N)
  
    END SUBROUTINE rk_CPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_CPop( T, H, Y, dY, d2Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION :: T, H, Y(N), dY(N), d2Y(N)   
      IF ( stack_ptr <= 0 ) THEN
        PRINT*,'Pop failed: empty buffer'
        STOP
      END IF  
      H = chk_H( stack_ptr )
      T = chk_T( stack_ptr )
      ! CALL DCOPY(NVAR,chk_Y(1,stack_ptr),1,Y,1)
      ! CALL DCOPY(NVAR,chk_dY(1,stack_ptr),1,dY,1)
      ! CALL DCOPY(NVAR,chk_d2Y(1,stack_ptr),1,d2Y,1)
      Y(1:N) = chk_Y(1:N,stack_ptr)
      dY(1:N) = chk_dY(1:N,stack_ptr)
      d2Y(1:N) = chk_d2Y(1:N,stack_ptr)
      stack_ptr = stack_ptr - 1
  
    END SUBROUTINE rk_CPop


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_FwdIntegrator( N,NADJ,Tstart,Tend,Y,AdjointType,GetQuad,IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
!~~~> Arguments
      INTEGER,  INTENT(IN)         :: N, NADJ
      DOUBLE PRECISION, INTENT(IN)    :: Tend, Tstart
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
      INTEGER,  INTENT(INOUT)        :: IERR
      INTEGER,  INTENT(IN)         :: AdjointType
      LOGICAL, INTENT(IN)          :: GetQuad
!~~~> Local variables
!      DOUBLE PRECISION    :: FJAC(N,N), E1(N,N)
!      COMPLEX(kind=dp) :: E2(N,N)   
      DOUBLE PRECISION, DIMENSION(N) :: Z1, Z2, Z3, Z4, SCAL, DZ1, DZ2, DZ3,  &
                                        DZ4, G, TMP, FCN0
      DOUBLE PRECISION  :: CONT(N,3), Tdirection, H, T, Hacc, Hnew, Hold, Fac,&
                           FacGus, Theta, Rerr, ErrOld, NewtonRate,           &
                           NewtonIncrement, Hratio, Qnewton,                  &
                           NewtonPredictedErr, NewtonIncrementOld, ThetaSD
      INTEGER :: NewtonIter, ISING, Nconsecutive, NewIt
      LOGICAL :: Reject, FirstStep, SkipJac, NewtonDone, SkipLU, Transp      
      DOUBLE PRECISION, DIMENSION(:), POINTER :: Zstage
            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initial setting
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (AdjointType == RK_adj_f90_discrete) THEN ! Save stage solution
        ALLOCATE(Zstage(N*3), STAT=i)
        IF (i/=0) THEN
          PRINT*,'Allocation of Zstage failed'
          STOP
        END IF
      END IF   
      T=Tstart

      Tdirection = SIGN(ONE,Tend-Tstart)
      H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , Hmax )
      IF (ABS(H) <= 10.d0*Roundoff) H = 1.0d-6
      H = SIGN(H,Tdirection)
      Hold      = H
      Reject    = .FALSE.
      FirstStep = .TRUE.
      SkipJac   = .FALSE.
      SkipLU    = .FALSE.
      Transp    = .FALSE.
      IF ((T+H*1.0001D0-Tend)*Tdirection >= ZERO) THEN
        H = Tend-T
      END IF
      Nconsecutive = 0
      CALL RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( (Tend-T)*Tdirection - Roundoff > ZERO )

        IF ( .NOT.SkipLU ) THEN 
        ! This time around skip the Jac update and LU

!~~~> Compute the Jacobian matrix
          IF ( .NOT.SkipJac ) THEN
            CALL LSS_Jac(T,Y,JAC)
            ISTATUS(Njac) = ISTATUS(Njac) + 1
          END IF

!~~~> Compute the matrices E1 and E2 and their decompositions
          CALL RK_Decomp(H,ISING)
          IF (ISING /= 0) THEN
            ISTATUS(Nsng) = ISTATUS(Nsng) + 1
            Nconsecutive = Nconsecutive + 1
            IF (Nconsecutive >= 5) THEN
              CALL RK_ErrorMsg(-12,T,H,IERR); RETURN
            END IF
            H=H*0.5d0
            Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
            CYCLE Tloop
          ELSE
            Nconsecutive = 0    
          END IF   
        END IF ! SkipLU
   
        ISTATUS(Nstp) = ISTATUS(Nstp) + 1
        IF (ISTATUS(Nstp) > Max_no_steps) THEN
          PRINT*,'Max number of time steps is ',Max_no_steps
          CALL RK_ErrorMsg(-9,T,H,IERR); RETURN
        END IF
        IF (0.1D0*ABS(H) <= ABS(T)*Roundoff)  THEN
          CALL RK_ErrorMsg(-10,T,H,IERR); RETURN
        END IF
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Loop for the simplified Newton iterations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
!~~~>  Starting values for Newton iteration
        IF ( FirstStep .OR. (.NOT.StartNewton) ) THEN
          CALL Set2zero(N,Z1)
          CALL Set2zero(N,Z2)
          CALL Set2zero(N,Z3)
        ELSE
          ! Evaluate quadratic polynomial
          CALL RK_Interpolate('eval',N,H,Hold,Z1,Z2,Z3,CONT)
        END IF
      
!~~~>  Initializations for Newton iteration
        NewtonDone = .FALSE.
        Fac = 0.5d0 ! Step reduction if too many iterations
      
NewtonLoop:DO  NewtonIter = 1, NewtonMaxit  
 
!~~~> Prepare the right-hand side
          CALL RK_PrepareRHS(N,T,H,Y,Z1,Z2,Z3,DZ1,DZ2,DZ3)
            
!~~~> Solve the linear systems
          CALL RK_Solve( N,H,DZ1,DZ2,DZ3,ISING )

          NewtonIncrement = SQRT( ( RK_ErrorNorm(N,SCAL,DZ1)**2 +     &
                          RK_ErrorNorm(N,SCAL,DZ2)**2 +               &
                          RK_ErrorNorm(N,SCAL,DZ3)**2 )/3.0d0 )            
          IF ( NewtonIter == 1 ) THEN
            Theta      = ABS(ThetaMin)
            NewtonRate = 2.0d0 
          ELSE
            Theta = NewtonIncrement/NewtonIncrementOld
            IF (Theta < 0.99d0) THEN
              NewtonRate = Theta/(ONE-Theta)
            ELSE ! Non-convergence of Newton: Theta too large
              EXIT NewtonLoop
            END IF
            IF ( NewtonIter < NewtonMaxit ) THEN 
              ! Predict error at the end of Newton process 
              NewtonPredictedErr = NewtonIncrement*                           &
                                   Theta**(NewtonMaxit-NewtonIter)/(ONE-Theta)
              IF (NewtonPredictedErr >= NewtonTol) THEN
              ! Non-convergence of Newton: predicted error too large
                Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                Fac=0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                EXIT NewtonLoop
              END IF
            END IF
          END IF
          NewtonIncrementOld = MAX(NewtonIncrement,Roundoff) 
          ! Update solution
          CALL DAXPY(N,-ONE,DZ1,1,Z1,1) ! Z1 <- Z1 - DZ1
          CALL DAXPY(N,-ONE,DZ2,1,Z2,1) ! Z2 <- Z2 - DZ2
          CALL DAXPY(N,-ONE,DZ3,1,Z3,1) ! Z3 <- Z3 - DZ3    
          ! Check error in Newton iterations
          NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
          IF (NewtonDone) THEN
            NewIt = NewtonIter
            EXIT NewtonLoop
          END IF  
          IF (NewtonIter == NewtonMaxit) THEN
            PRINT*, 'Slow or no convergence in Newton Iteration:',    &
                       ' Max no. of Newton iterations reached'
          END IF

        END DO NewtonLoop
         
        IF (.NOT.NewtonDone) THEN
         !CALL RK_ErrorMsg(-12,T,H,IERR);
          H = Fac*H; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
!???       ISTATUS(Nrej) = ISTATUS(Nrej) + 1
          CYCLE Tloop
        END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> SDIRK Stage
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF (SdirkError) THEN

!~~~>  Starting values for Newton iterations
          Z4(1:N) = Z3(1:N)
! ???       
!~~~>   Prepare the loop-independent part of the right-hand side
          CALL FUN(N,T,Y,DZ4)
          ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       
! G = H*rkBgam(0)*DZ4 + rkTheta(1)*Z1 + rkTheta(2)*Z2 + rkTheta(3)*Z3
          CALL Set2Zero(N, G)
          CALL DAXPY(N,rkBgam(0)*H, DZ4,1,G,1) 
          CALL DAXPY(N,rkTheta(1),Z1,1,G,1)
          CALL DAXPY(N,rkTheta(2),Z2,1,G,1)
          CALL DAXPY(N,rkTheta(3),Z3,1,G,1)

          !~~~>  Initializations for Newton iteration
          NewtonDone = .FALSE.
          Fac = 0.5d0 ! Step reduction factor if too many iterations
            
SDNewtonLoop:DO NewtonIter = 1, NewtonMaxit

!~~~>   Prepare the loop-dependent part of the right-hand side
            CALL WADD(N,Y,Z4,TMP)         ! TMP <- Y + Z4
            CALL FUN(N,T+H,TMP,DZ4)       ! DZ4 <- Fun(Y+Z4)
            ISTATUS(Nfun) = ISTATUS(Nfun) + 1
!         DZ4(1:N) = (G(1:N)-Z4(1:N))*(rkGamma/H) + DZ4(1:N)
            CALL DAXPY (N, -ONE*rkGamma/H, Z4, 1, DZ4, 1)
            CALL DAXPY (N, rkGamma/H, G,1, DZ4,1)

!~~~>   Solve the linear system
            CALL LSS_Solve(Transp,DZ4)
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1 
!~~~>   Check convergence of Newton iterations
            NewtonIncrement = RK_ErrorNorm(N,SCAL,DZ4)
            IF ( NewtonIter == 1 ) THEN
              ThetaSD      = ABS(ThetaMin)
              NewtonRate = 2.0d0 
            ELSE
              ThetaSD = NewtonIncrement/NewtonIncrementOld
              IF (ThetaSD < 0.99d0) THEN
                NewtonRate = ThetaSD/(ONE-ThetaSD)
                ! Predict error at the end of Newton process 
                NewtonPredictedErr = NewtonIncrement*                         &
                               ThetaSD**(NewtonMaxit-NewtonIter)/(ONE-ThetaSD)
                IF (NewtonPredictedErr >= NewtonTol) THEN
                  !Non-convergence of Newton: predicted error too large
                  !print*,'Error too large: ', NewtonPredictedErr
                  Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                  Fac=0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                  EXIT SDNewtonLoop
                END IF
              ELSE ! Non-convergence of Newton: ThetaSD too large
                !print*,'Theta too large: ',ThetaSD
                EXIT SDNewtonLoop
              END IF
            END IF
            NewtonIncrementOld = NewtonIncrement
            ! Update solution: Z4 <-- Z4 + DZ4
            CALL DAXPY(N,ONE,DZ4,1,Z4,1) 
            
            ! Check error in Newton iterations
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
            IF (NewtonDone) EXIT SDNewtonLoop
            
          END DO SDNewtonLoop
           
          IF (.NOT.NewtonDone) THEN
            H = Fac*H;Reject=.TRUE.;SkipJac = .TRUE.;SkipLU = .FALSE.
            CYCLE Tloop
          END IF
        END IF
!~~~>  End of implified SDIRK Newton iterations
            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Error estimation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF (SdirkError) THEN
!         DZ4(1:N) =  rkD(1)*Z1 + rkD(2)*Z2 + rkD(3)*Z3 - Z4    
          CALL Set2Zero(N, DZ4)
          IF (rkD(1) /= ZERO) CALL DAXPY(N, rkD(1), Z1, 1, DZ4, 1)
          IF (rkD(2) /= ZERO) CALL DAXPY(N, rkD(2), Z2, 1, DZ4, 1)
          IF (rkD(3) /= ZERO) CALL DAXPY(N, rkD(3), Z3, 1, DZ4, 1)
          CALL DAXPY(N, -ONE, Z4, 1, DZ4, 1)
          Rerr = RK_ErrorNorm(N,SCAL,DZ4)    
        ELSE
          CALL RK_ErrorEstimate(N,H,Y,T,Z1,Z2,Z3,SCAL,Rerr,FirstStep,Reject)
        END IF
!~~~> Computation of new step size Hnew
        Fac  = Rerr**(-ONE/rkELO)*                                    &
            MIN(FacSafe,(ONE+2*NewtonMaxit)/(NewtonIter+2*NewtonMaxit))
        Fac  = MIN(FacMax,MAX(FacMin,Fac))
        Hnew = Fac*H

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Accept/reject step 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
accept:IF (Rerr < ONE) THEN !~~~> STEP IS ACCEPTED
          IF (AdjointType == RK_adj_f90_discrete) THEN 
            ! Save stage solution
            ! CALL DCOPY(N,Z1,1,Zstage(1),1)
            ! CALL DCOPY(N,Z2,1,Zstage(N+1),1)
            ! CALL DCOPY(N,Z3,1,Zstage(2*N+1),1)
            Zstage(1:N)       = Z1(1:N) 
            Zstage(N+1:2*N)   = Z2(1:N)
            Zstage(2*N+1:3*N) = Z3(1:N)
            ! Push old Y - Y at the beginning of the stage
            CALL rk_DPush(T, H, Y, Zstage, NewIt)
          END IF

!~~~> Update the results for the quadrature term
          IF(GetQuad) CALL RK_UpdateQuad(N,NADJ,H,T,Y,Z1,Z2,Z3,Q)

          FirstStep=.FALSE.
          ISTATUS(Nacc) = ISTATUS(Nacc) + 1
          IF (Gustafsson) THEN
          !~~~> Predictive controller of Gustafsson
            IF (ISTATUS(Nacc) > 1) THEN
              FacGus=FacSafe*(H/Hacc)*(Rerr**2/ErrOld)**(-0.25d0)
              FacGus=MIN(FacMax,MAX(FacMin,FacGus))
              Fac=MIN(Fac,FacGus)
              Hnew = Fac*H
            END IF
            Hacc=H
            ErrOld=MAX(1.0d-2,Rerr)
          END IF
          Hold = H
          T = T+H 

          ! Update solution: Y <- Y + sum (d_i Z_i)
          IF (rkD(1) /= ZERO) CALL DAXPY(N,rkD(1),Z1,1,Y,1)
          IF (rkD(2) /= ZERO) CALL DAXPY(N,rkD(2),Z2,1,Y,1)
          IF (rkD(3) /= ZERO) CALL DAXPY(N,rkD(3),Z3,1,Y,1)

          ! Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3
          IF (StartNewton) THEN
            CALL RK_Interpolate('make',N,H,Hold,Z1,Z2,Z3,CONT)
          ENDIF
          CALL RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)
          RSTATUS(Ntexit) = T
          RSTATUS(Nhnew)  = Hnew
          RSTATUS(Nhacc)  = H
          Hnew = Tdirection*MIN( MAX(ABS(Hnew),Hmin) , Hmax )
          IF (Reject) Hnew = Tdirection*MIN(ABS(Hnew),ABS(H))
            Reject = .FALSE.
            IF ((T+Hnew/Qmin-Tend)*Tdirection >=  ZERO) THEN
              H = Tend-T
            ELSE
            Hratio = Hnew/H
            ! Reuse the LU decomposition
            SkipLU = (Theta<=ThetaMin).AND.(Hratio>=Qmin).AND.(Hratio<=Qmax)
!???            SkipLU = .false.
            IF (.NOT.SkipLU) H=Hnew
          END IF
          ! If convergence is fast enough, do not update Jacobian
          ! SkipJac = (Theta <= ThetaMin)
          SkipJac = .FALSE.

        ELSE accept !~~~> Step is rejected
          IF (FirstStep .OR. Reject) THEN
            H = FacRej*H
          ELSE
            H = Hnew
          END IF
          Reject   = .TRUE.
          SkipJac  = .TRUE.
          SkipLU   = .FALSE. 
          IF (ISTATUS(Nacc) >= 1) ISTATUS(Nrej) = ISTATUS(Nrej) + 1
        END IF accept
      
      END DO Tloop

!~~~> Save last state: only needed for continuous adjoint
      IF ( (AdjointType == RK_adj_f90_continuous) .OR.                &
           (AdjointType == RK_adj_f90_simple_continuous) ) THEN
        CALL Fun(T,Y,Fcn0)
        CALL rk_CPush( T, H, Y, Fcn0, DZ3)
!~~~> Deallocate stage buffer: only needed for discrete adjoint
      ELSEIF (AdjointType == RK_adj_f90_discrete) THEN 
        DEALLOCATE(Zstage, STAT=i)
        IF (i/=0) THEN
          PRINT*,'Deallocation of Zstage failed'
          STOP
        END IF
      END IF   
   
      ! Successful exit
      IERR = 1  

    END SUBROUTINE RK_FwdIntegrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_DadjInt( N,NADJ,NP,Lambda,Mu,Tstart,Tend,T,GetQuad,IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
!~~~> Arguments
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)     :: N, NADJ, NP
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ), Mu(NP,NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tstart, Tend
      DOUBLE PRECISION, INTENT(INOUT) :: T
      LOGICAL, INTENT(IN)             :: GetQuad
      INTEGER,  INTENT(INOUT)         :: IERR

!~~~> Local variables
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WY1, WY2, WY3, WP1,WP2,&
                                                   WP3, FPJAC1, FPJAC2, FPJAC3
      DOUBLE PRECISION, DIMENSION(N)    :: U1, U2, U3
      DOUBLE PRECISION, DIMENSION(NP)   :: V1, V2, V3
      DOUBLE PRECISION, DIMENSION(N) :: SCAL, DU1, DU2, DU3
      DOUBLE PRECISION :: Y(N), Zstage(3*N), X(3*N)
      DOUBLE PRECISION  :: H, NewtonRate, NewtonIncrement, NewtonIncrementOld
      INTEGER :: NewtonIter, ISING, iadj, NewIt
      LOGICAL :: Reject, NewtonDone, NewtonConverge, SkipJac, Transp
      DOUBLE PRECISION :: T1, Theta, Alpha,Beta
      DOUBLE PRECISION, DIMENSION(N) :: TMP, G1, G2, G3
            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initial setting
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> If there is a quadrature term
      IF (GetQuad) THEN
        ALLOCATE( WY1(N,NADJ), WY2(N,NADJ), WY3(N,NADJ), WP1(NP,NADJ),        &
                 WP2(NP,NADJ), WP3(NP,NADJ),STAT=ISING )
      
        IF(ISING .NE. 0) STOP 'allocation error for WYs,WPs'
        WY1(:,:) = 0.0d0
        WY2(:,:) = 0.0d0
        WY3(:,:) = 0.0d0
        WP1(:,:) = 0.0d0
        WP2(:,:) = 0.0d0
        WP3(:,:) = 0.0d0
      END IF

      ALLOCATE(FPJAC1(N,NP),FPJAC2(N,NP),FPJAC3(N,NP),STAT=ISING)
      IF(ISING .NE. 0) STOP 'allocation error for FPJACs'
      FPJAC1(:,:) = 0.0d0
      FPJAC2(:,:) = 0.0d0
      FPJAC3(:,:) = 0.0d0
      T1 = Tend      
!      Tdirection = SIGN(ONE,Tend-Tstart)
      NewtonConverge = .TRUE.
      Reject = .FALSE.
      SkipJac = .FALSE. 
      Transp = .FALSE.
!~~~> Time loop begins below 
TimeLoop:DO WHILE ( stack_ptr > 0 )

        IF (.not.Reject) THEN
   
!~~~>  Recover checkpoints for stage values and vectors
          CALL rk_DPop(T, H, Y, Zstage, NewIt)
!~~~>  Compute LU decomposition 
          IF (.NOT.SaveLU) THEN
!~~~> Compute the Jacobian matrix
            CALL LSS_Jac(T,Y,JAC)
            ISTATUS(Njac) = ISTATUS(Njac) + 1
!~~~> Compute the matrices E1 and E2 and their decompositions
            CALL RK_Decomp(H,ISING)
          END IF
          ISTATUS(Njac) = ISTATUS(Njac) + 3  
!~~~>   Jacobian values at stage vectors
          CALL WADD(N,Y,Zstage(1),TMP)       ! TMP  <- Y + Z1
          CALL LSS_Jac1(T+rkC(1)*H,TMP,JAC)
          CALL JACP(N,NP,T+rkC(1)*H,TMP,FPJAC1)
          IF(GetQuad) Then
            CALL DRDY(NADJ,N,N,T+rkC(1)*H,TMP,WY1)
            CALL DRDP(NADJ,N,NP,T+rkC(1)*H,TMP,WP1)
          END IF

          CALL WADD(N,Y,Zstage(1+N),TMP)     ! TMP  <- Y + Z2
          CALL LSS_Jac2(T+rkC(2)*H,TMP,JAC)
          CALL JACP(N,NP,T+rkC(2)*H,TMP,FPJAC2)
          IF(GetQuad) Then
            CALL DRDY(NADJ,N,N,T+rkC(2)*H,TMP,WY2)
            CALL DRDP(NADJ,N,NP,T+rkC(2)*H,TMP,WP2)
          END IF

          CALL WADD(N,Y,Zstage(1+2*N),TMP)   ! TMP  <- Y + Z3
          CALL LSS_Jac3(T+rkC(3)*H,TMP,JAC)
          CALL JACP(N,NP,T+rkC(3)*H,TMP,FPJAC3)
          IF(GetQuad) Then
            CALL DRDY(NADJ,N,N,T+rkC(3)*H,TMP,WY3)  
            CALL DRDP(NADJ,N,NP,T+rkC(3)*H,TMP,WP3)  
          END IF
        END IF ! .not.Reject

111     CONTINUE

        IF ( (AdjointSolve == Solve_adaptive .and. .not.NewtonConverge&
          ) .or. (AdjointSolve == Solve_direct) ) THEN
          CALL LSS_Decomp_Big(H,rkA,ISING) 
          ISTATUS(Ndec) = ISTATUS(Ndec) + 1
          !  Use full big algebra:    
        END IF
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Loop for the simplified Newton iterations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Adj:DO iadj = 1, NADJ      
!~~~>  Starting values for Newton iteration
          ! CALL DCOPY(N,Lambda(1,iadj),1,U1,1)
          ! CALL DCOPY(N,Lambda(1,iadj),1,U2,1)
          ! CALL DCOPY(N,Lambda(1,iadj),1,U3,1)
          CALL Set2Zero(N,U1)
          CALL Set2Zero(N,U2)
          CALL Set2Zero(N,U3)
      
!~~~>  Initializations for Newton iteration
          NewtonDone = .FALSE.

!~~~>    Right Hand Side - part G for all Newton iterations
          CALL SET2ZERO(N,G1)
          CALL SET2ZERO(N,G2)
          CALL SET2ZERO(N,G3)
          CALL LSS_Mul_Jactr1(TMP,Lambda(1,iadj))
          CALL DAXPY(N,-H*rkB(1),TMP,1,G1,1) ! R1 <- R1 - h*B_1*F1
          IF(GetQuad) CALL DAXPY(N,-H*rkB(1),WY1(1,iadj),1,G1,1)
          CALL LSS_Mul_Jactr2(TMP,Lambda(1,iadj))
          CALL DAXPY(N,-H*rkB(2),TMP,1,G2,1) ! R2 <- R2 - h*B_2*F2
          IF(GetQuad) CALL DAXPY(N,-H*rkB(2),WY2(1,iadj),1,G2,1)
          CALL LSS_Mul_Jactr3(TMP,Lambda(1,iadj))
          CALL DAXPY(N,-H*rkB(3),TMP,1,G3,1) ! R3 <- R3 - h*B_3*F3
          IF(GetQuad) CALL DAXPY(N,-H*rkB(3),WY3(1,iadj),1,G3,1)

          IF ( (AdjointSolve == Solve_adaptive .and. NewtonConverge)  &
            .or. (AdjointSolve == Solve_fixed) ) THEN

NewtonLoopAdj:DO  NewtonIter = 1, NewtonMaxit  

              !~~~> Prepare the right-hand side
              CALL RK_PrepareRHS_Adj( N,H,Lambda(1,iadj),U1,U2,U3,G1,G2,G3,   &
                                      DU1,DU2,DU3)

              !~~~> Solve the linear systems
              CALL RK_SolveTR( N,H,DU1,DU2,DU3,ISING )

!~~~> The following code performs an adaptive number of Newton
!     iterations for solving adjoint system
              IF (AdjointSolve == Solve_adaptive) THEN

                CALL RK_ErrorScale( N,ITOL,AbsTol_adj(1:N,iadj),      &
                        RelTol_adj(1:N,iadj),Lambda(1:N,iadj),SCAL )

                ! SCAL(1:N) = 1.0d0
                NewtonIncrement = SQRT( ( RK_ErrorNorm(N,SCAL,DU1)**2 +       &
                                RK_ErrorNorm(N,SCAL,DU2)**2 +                 &
                                RK_ErrorNorm(N,SCAL,DU3)**2 )/3.0d0 )
            
                IF ( NewtonIter == 1 ) THEN
                  Theta      = ABS(ThetaMin)
                  NewtonRate = 2.0d0 
                ELSE
                  Theta = NewtonIncrement/NewtonIncrementOld
                  IF (Theta < 0.99d0) THEN
                    NewtonRate = Theta/(ONE-Theta)
                  ELSE ! Non-convergence of Newton: Theta too large
                    Reject = .TRUE.
                    NewtonDone = .FALSE.
                    EXIT NewtonLoopAdj
                  END IF
                END IF
  
                NewtonIncrementOld = MAX(NewtonIncrement,Roundoff) 

              END IF ! (AdjointSolve == Solve_adaptive)
            
              ! Update solution
              CALL DAXPY(N,-ONE,DU1,1,U1,1) ! U1 <- U1 - DU1
              CALL DAXPY(N,-ONE,DU2,1,U2,1) ! U2 <- U2 - DU2
              CALL DAXPY(N,-ONE,DU3,1,U3,1) ! U3 <- U3 - DU3

              IF (AdjointSolve == Solve_adaptive) THEN
              ! When performing an adaptive number of iterations
              !       check the error in Newton iterations
                NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
                IF ((NewtonDone).and.(NewtonIter>NewIt+1)) THEN
                  EXIT NewtonLoopAdj
                END IF
              ELSE IF (AdjointSolve == Solve_fixed) THEN
                IF (NewtonIter>MAX(NewIt+1,6)) EXIT NewtonLoopAdj
              END IF   
            
            END DO NewtonLoopAdj
      
            IF ((AdjointSolve==Solve_adaptive).AND.(.NOT.NewtonDone)) THEN
            ! print*,'Newton iterations do not converge, switching to full system.'
              NewtonConverge = .FALSE.
              Reject = .TRUE.
              GOTO 111
            END IF
            ! compute adjoint varibales for Mu 
            Alpha = 1
            Beta  = 0
            CALL SET2ZERO(N,TMP)
            CALL DAXPY(N,H*rkA(1,1),U1,1,TMP,1) ! TMP < h*A_11*U1
            CALL DAXPY(N,H*rkA(2,1),U2,1,TMP,1) ! TMP < TMP + h*A_21*U2
            CALL DAXPY(N,H*rkA(3,1),U3,1,TMP,1) ! TMP < TMP + h*A_31*U3
            CALL DAXPY(N,H*rkB(1),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B_1*Lambda
            CALL DGEMV('T',N,NP,Alpha,FPJAC1,N,TMP,1,Beta,V1,1) ! V1 = Jacp^T *TMP
            IF(GetQuad) CALL DAXPY(NP,H*rkB(1),WP1(1,iadj),1,V1,1) ! V1 = V1+ h*B_1 * WP1      

            CALL SET2ZERO(N,TMP)
            CALL DAXPY(N,H*rkA(1,2),U1,1,TMP,1) ! TMP < h*A_12*U1
            CALL DAXPY(N,H*rkA(2,2),U2,1,TMP,1) ! TMP < TMP + h*A_22*U2
            CALL DAXPY(N,H*rkA(3,2),U3,1,TMP,1) ! TMP < TMP + h*A_32*U3
            CALL DAXPY(N,H*rkB(2),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B_2*Lambda
            CALL DGEMV('T',N,NP,Alpha,FPJAC2,N,TMP,1,Beta,V2,1) ! V2 = Jacp^T *TMP
            IF(GetQuad) CALL DAXPY(NP,H*rkB(2),WP2(1,iadj),1,V2,1) ! V2 = V2+ h*B_2 * WP2 

            CALL SET2ZERO(N,TMP)
            CALL DAXPY(N,H*rkA(1,3),U1,1,TMP,1) ! TMP < h*A_13*U1
            CALL DAXPY(N,H*rkA(2,3),U2,1,TMP,1) ! TMP < TMP + h*A_23*U2
            CALL DAXPY(N,H*rkA(3,3),U3,1,TMP,1) ! TMP < TMP + h*A_33*U3
            CALL DAXPY(N,H*rkB(3),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B_3*Lambda
            CALL DGEMV('T',N,NP,Alpha,FPJAC3,N,TMP,1,Beta,V3,1)! V3 = Jacp^T *TMP
            IF(GetQuad) CALL DAXPY(NP,H*rkB(3),WP3(1,iadj),1,V3,1) ! V3 = V3+ h*B_3 * WP3 

            ! Update adjoint solution for Lambda: Y_adj <- Y_adj + sum (U_i)
            CALL DAXPY(N,ONE,U1,1,Lambda(1,iadj),1)
            CALL DAXPY(N,ONE,U2,1,Lambda(1,iadj),1)
            CALL DAXPY(N,ONE,U3,1,Lambda(1,iadj),1)

            ! Update adjoint solution for Mu
            CALL DAXPY(NP,ONE,V1,1,Mu(1,iadj),1)
            CALL DAXPY(NP,ONE,V2,1,Mu(1,iadj),1)
            CALL DAXPY(NP,ONE,V3,1,Mu(1,iadj),1)

          ELSE ! NewtonConverge = .false.

            X(1:N)       = -G1(1:N)
            X(N+1:2*N)   = -G2(1:N)
            X(2*N+1:3*N) = -G3(1:N)
            Transp = .TRUE.
            CALL LSS_Solve_Big(Transp, X)
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1
            ! CALL WGESL('T',3*N,Jbig,IPbig,X)
            CALL SET2ZERO(N,TMP)
            CALL DAXPY(N,H*rkA(1,1),X(1:N),1,TMP,1) ! TMP < h*A_11*U1
            CALL DAXPY(N,H*rkA(2,1),X(N+1:2*N),1,TMP,1) ! TMP < TMP + h*A_21*U2
            CALL DAXPY(N,H*rkA(3,1),X(2*N+1:3*N),1,TMP,1) ! TMP < TMP + h*A_31*U3
            CALL DAXPY(N,H*rkB(1),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B_1*Lambda
            CALL DGEMV('T',N,NP,Alpha,FPJAC1,N,TMP,1,Beta,V1,1) ! V1 = Jacp^T *TMP
            IF(GetQuad) CALL DAXPY(NP,H*rkB(1),WP1(1,iadj),1,V1,1) ! V1 = V1+ h*B_1 * WP1      
            CALL SET2ZERO(N,TMP)
            CALL DAXPY(N,H*rkA(1,2),X(1:N),1,TMP,1) ! TMP < h*A_12*U1
            CALL DAXPY(N,H*rkA(2,2),X(N+1:2*N),1,TMP,1) ! TMP < TMP + h*A_22*U2
            CALL DAXPY(N,H*rkA(3,2),X(2*N+1:3*N),1,TMP,1) ! TMP < TMP + h*A_32*U3
            CALL DAXPY(N,H*rkB(2),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B_2*Lambda
            CALL DGEMV('T',N,NP,Alpha,FPJAC2,N,TMP,1,Beta,V2,1) ! V2 = Jacp^T *TMP
            IF(GetQuad) CALL DAXPY(NP,H*rkB(2),WP2(1,iadj),1,V2,1) ! V2 = V2+ h*B_2 * WP2 

            CALL SET2ZERO(N,TMP)
            CALL DAXPY(N,H*rkA(1,3),X(1:N),1,TMP,1) ! TMP < h*A_13*U1
            CALL DAXPY(N,H*rkA(2,3),X(N+1:2*N),1,TMP,1) ! TMP < TMP + h*A_23*U2
            CALL DAXPY(N,H*rkA(3,3),X(2*N+1:3*N),1,TMP,1) ! TMP < TMP + h*A_33*U3
            CALL DAXPY(N,H*rkB(3),Lambda(1,iadj),1,TMP,1) ! TMP < TMP + h*B_3*Lambda
            CALL DGEMV('T',N,NP,Alpha,FPJAC3,N,TMP,1,Beta,V3,1) ! V3 = Jacp^T *TMP
            IF(GetQuad) CALL DAXPY(NP,H*rkB(3),WP3(1,iadj),1,V3,1) ! V3 = V3+ h*B_3 * WP3 

            Lambda(1:N,iadj) = Lambda(1:N,iadj)+X(1:N)+X(N+1:2*N)+X(2*N+1:3*N)
            Mu(1:NP,iadj) = Mu(1:NP,iadj)+V1(1:NP)+V2(1:NP)+V3(1:NP)

            IF ((AdjointSolve==Solve_adaptive).AND.(iadj>=NADJ)) THEN
              NewtonConverge = .TRUE.
              Reject = .FALSE.
            END IF
     
          END IF ! NewtonConverge
 
        END DO Adj
 
        T1 = T1-H
        ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      END DO TimeLoop

      DEALLOCATE(FPJAC1,FPJAC2,FPJAC3,STAT=ISING)
      IF(ISING .NE. 0) STOP 'deallocation error for FPJACs'
      IF(GetQuad) THEN
        DEALLOCATE(WY1,WY2,WY3,WP1,WP2,WP3,STAT=ISING)
        IF(ISING .NE. 0) STOP 'deallocation error for WYs,WPs'
      END IF
      ! Successful exit
      IERR = 1  

    END SUBROUTINE RK_DadjInt


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_CadjInt ( NADJ, Y, Tstart, Tend, T, IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
!~~~> Arguments
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)     :: NADJ
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tstart, Tend
      DOUBLE PRECISION, INTENT(INOUT) :: T
      INTEGER,  INTENT(INOUT)        :: IERR

    END SUBROUTINE rk_CadjInt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_SimpleCadjInt ( NADJ, Y, Tstart, Tend, T, IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
!~~~> Arguments
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)     :: NADJ
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tstart, Tend
      DOUBLE PRECISION, INTENT(INOUT) :: T
      INTEGER,  INTENT(INOUT)        :: IERR

    END SUBROUTINE rk_SimpleCadjInt


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: T, H
      INTEGER, INTENT(IN)  :: Code
      INTEGER, INTENT(INOUT) :: IERR
      IERR = Code
      PRINT * , &
           'Forced exit from RungeKutta due to the following error:'
      SELECT CASE (Code)
      CASE (-1)
        PRINT * , '--> Improper value for maximal no of steps'
      CASE (-2)
        PRINT * , '--> Improper value for maximal no of Newton iterations'
      CASE (-3)
        PRINT * , '--> Hmin/Hmax/Hstart must be positive'
      CASE (-4)
        PRINT * , '--> Improper values for FacMin/FacMax/FacSafe/FacRej'
      CASE (-5)
        PRINT * , '--> Improper value for ThetaMin'
      CASE (-6)
        PRINT * , '--> Newton stopping tolerance too small'
      CASE (-7)
        PRINT * , '--> Improper values for Qmin, Qmax'
      CASE (-8)
        PRINT * , '--> Tolerances are too small'
      CASE (-9)
        PRINT * , '--> No of steps exceeds maximum bound'
      CASE (-10)
        PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
      CASE (-11)
        PRINT * , '--> Matrix is repeatedly singular'
      CASE (-12)
        PRINT * , '--> Non-convergence of Newton iterations'
      CASE (-13)
        PRINT * , '--> Requested RK method not implemented'
      CASE DEFAULT
        PRINT *, 'Unknown Error code: ', Code
      END SELECT

      WRITE(6,FMT="(5X,'T=',E12.5,'  H=',E12.5)") T, H 

    END SUBROUTINE RK_ErrorMsg


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N, ITOL
      DOUBLE PRECISION, INTENT(IN) :: AbsTol(*), RelTol(*), Y(N)
      DOUBLE PRECISION, INTENT(INOUT) :: SCAL(N)
      INTEGER :: i
   
      IF (ITOL==0) THEN
        DO i=1,N
          SCAL(i)= ONE/(AbsTol(1)+RelTol(1)*ABS(Y(i)))
        END DO
      ELSE
        DO i=1,N
          SCAL(i)=ONE/(AbsTol(i)+RelTol(i)*ABS(Y(i)))
        END DO
      END IF
      
    END SUBROUTINE RK_ErrorScale


!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SUBROUTINE RK_Transform(N,Tr,Z1,Z2,Z3,W1,W2,W3)
!!~~~>                 W <-- Tr x Z
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      IMPLICIT NONE
!      INTEGER :: N, i
!      DOUBLE PRECISION :: Tr(3,3),Z1(N),Z2(N),Z3(N),W1(N),W2(N),W3(N)
!      DOUBLE PRECISION :: x1, x2, x3
!      DO i=1,N
!          x1 = Z1(i); x2 = Z2(i); x3 = Z3(i)
!          W1(i) = Tr(1,1)*x1 + Tr(1,2)*x2 + Tr(1,3)*x3
!          W2(i) = Tr(2,1)*x1 + Tr(2,2)*x2 + Tr(2,3)*x3
!          W3(i) = Tr(3,1)*x1 + Tr(3,2)*x2 + Tr(3,3)*x3
!      END DO
!  END SUBROUTINE RK_Transform
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_Interpolate(action,N,H,Hold,Z1,Z2,Z3,CONT)
!~~~>   Constructs or evaluates a quadratic polynomial
!         that interpolates the Z solution at current step
!         and provides starting values for the next step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: N, i
      DOUBLE PRECISION :: H,Hold,Z1(N),Z2(N),Z3(N),CONT(N,3)
      DOUBLE PRECISION :: r, x1, x2, x3, den
      CHARACTER(LEN=4) :: action
       
      SELECT CASE (action) 
      CASE ('make')
         ! Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3
         den = (rkC(3)-rkC(2))*(rkC(2)-rkC(1))*(rkC(1)-rkC(3))
         DO i=1,N
             CONT(i,1)=(-rkC(3)**2*rkC(2)*Z1(i)+Z3(i)*rkC(2)*rkC(1)**2 &
                        +rkC(2)**2*rkC(3)*Z1(i)-rkC(2)**2*rkC(1)*Z3(i) &
                        +rkC(3)**2*rkC(1)*Z2(i)-Z2(i)*rkC(3)*rkC(1)**2)&
                        /den-Z3(i)
             CONT(i,2)= -( rkC(1)**2*(Z3(i)-Z2(i)) + rkC(2)**2*(Z1(i)  &
                          -Z3(i)) +rkC(3)**2*(Z2(i)-Z1(i)) )/den
             CONT(i,3)= ( rkC(1)*(Z3(i)-Z2(i)) + rkC(2)*(Z1(i)-Z3(i))  &
                           +rkC(3)*(Z2(i)-Z1(i)) )/den
         END DO
      CASE ('eval')
          ! Evaluate quadratic polynomial
          r = H/Hold
         x1 = ONE + rkC(1)*r
         x2 = ONE + rkC(2)*r
         x3 = ONE + rkC(3)*r
         DO i=1,N
            Z1(i) = CONT(i,1)+x1*(CONT(i,2)+x1*CONT(i,3))
            Z2(i) = CONT(i,1)+x2*(CONT(i,2)+x2*CONT(i,3))
            Z3(i) = CONT(i,1)+x3*(CONT(i,2)+x3*CONT(i,3))
         END DO
      END SELECT   
    END SUBROUTINE RK_Interpolate


    SUBROUTINE RK_UpdateQuad(N,NADJ,H,T,Y,Z1,Z2,Z3,Q)
      INTEGER :: N,NADJ
      DOUBLE PRECISION :: T,H
      DOUBLE PRECISION :: Q(NADJ),Y(N),Z1(N),Z2(N),Z3(N)
!~~~> local variables
      DOUBLE PRECISION :: F(NADJ),TMP(N)
      
      TMP(:) = 0.d0
      F(:) = 0.d0
      CALL WADD(N,Y,Z1,TMP)
      CALL QFUN(N,NADJ,T+rkC(1)*H,TMP,F)
      CALL DAXPY(NADJ,H*rkB(1),F,1,Q,1)
      CALL WADD(N,Y,Z2,TMP)
      CALL QFUN(N,NADJ,T+rkC(2)*H,TMP,F)
      CALL DAXPY(NADJ,H*rkB(2),F,1,Q,1)
      CALL WADD(N,Y,Z3,TMP)
      CALL QFUN(N,NADJ,T+rkC(3)*H,TMP,F)
      CALL DAXPY(NADJ,H*rkB(3),F,1,Q,1)
    END SUBROUTINE RK_UpdateQuad

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_PrepareRHS(N,T,H,Y,Z1,Z2,Z3,R1,R2,R3)
!~~~> Prepare the right-hand side for Newton iterations
!     R = Z - hA x F
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER :: N
      DOUBLE PRECISION :: T,H
      DOUBLE PRECISION, DIMENSION(N) :: Y,Z1,Z2,Z3,F,R1,R2,R3,TMP

      CALL DCOPY(N,Z1,1,R1,1) ! R1 <- Z1
      CALL DCOPY(N,Z2,1,R2,1) ! R2 <- Z2
      CALL DCOPY(N,Z3,1,R3,1) ! R3 <- Z3

      CALL WADD(N,Y,Z1,TMP)              ! TMP <- Y + Z1
      CALL FUN(N,T+rkC(1)*H,TMP,F)    ! F1 <- Fun(Y+Z1)
      CALL DAXPY(N,-H*rkA(1,1),F,1,R1,1) ! R1 <- R1 - h*A_11*F1
      CALL DAXPY(N,-H*rkA(2,1),F,1,R2,1) ! R2 <- R2 - h*A_21*F1
      CALL DAXPY(N,-H*rkA(3,1),F,1,R3,1) ! R3 <- R3 - h*A_31*F1

      CALL WADD(N,Y,Z2,TMP)              ! TMP <- Y + Z2
      CALL FUN(N,T+rkC(2)*H,TMP,F)    ! F2 <- Fun(Y+Z2)
      CALL DAXPY(N,-H*rkA(1,2),F,1,R1,1) ! R1 <- R1 - h*A_12*F2
      CALL DAXPY(N,-H*rkA(2,2),F,1,R2,1) ! R2 <- R2 - h*A_22*F2
      CALL DAXPY(N,-H*rkA(3,2),F,1,R3,1) ! R3 <- R3 - h*A_32*F2

      CALL WADD(N,Y,Z3,TMP)              ! TMP <- Y + Z3
      CALL FUN(N,T+rkC(3)*H,TMP,F)    ! F3 <- Fun(Y+Z3)
      CALL DAXPY(N,-H*rkA(1,3),F,1,R1,1) ! R1 <- R1 - h*A_13*F3
      CALL DAXPY(N,-H*rkA(2,3),F,1,R2,1) ! R2 <- R2 - h*A_23*F3
      CALL DAXPY(N,-H*rkA(3,3),F,1,R3,1) ! R3 <- R3 - h*A_33*F3
            
    END SUBROUTINE RK_PrepareRHS
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_PrepareRHS_Adj(N,H,Lambda,U1,U2,U3,G1,G2,G3,R1,R2,R3)
!~~~> Prepare the right-hand side for Newton iterations
!     R = Z_adj - hA x Jac*Z_adj - h J^t b lambda
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(IN) :: H
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: U1,U2,U3,Lambda,G1,G2,G3
      DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: R1,R2,R3
      DOUBLE PRECISION, DIMENSION(N) ::  F,TMP


      CALL DCOPY(N,G1,1,R1,1) ! R1 <- G1
      CALL DCOPY(N,G2,1,R2,1) ! R2 <- G2
      CALL DCOPY(N,G3,1,R3,1) ! R3 <- G3

      CALL SET2ZERO(N,F)
      CALL DAXPY(N,-H*rkA(1,1),U1,1,F,1) ! F1 <- -h*A_11*U1
      CALL DAXPY(N,-H*rkA(2,1),U2,1,F,1) ! F1 <- F1 - h*A_21*U2
      CALL DAXPY(N,-H*rkA(3,1),U3,1,F,1) ! F1 <- F1 - h*A_31*U3
      CALL LSS_Mul_Jactr1(TMP,F)
      CALL DAXPY(N,ONE,U1,1,TMP,1) ! R1 <- U1 -Jac(Y+Z1)^t*h*sum(A_j1*U_j)
      CALL DAXPY(N,ONE,TMP,1,R1,1) ! R1 <- U1 -Jac(Y+Z1)^t*h*sum(A_j1*U_j)

      CALL SET2ZERO(N,F)
      CALL DAXPY(N,-H*rkA(1,2),U1,1,F,1) ! F2 <- -h*A_12*U1
      CALL DAXPY(N,-H*rkA(2,2),U2,1,F,1) ! F2 <- F2 - h*A_22*U2
      CALL DAXPY(N,-H*rkA(3,2),U3,1,F,1) ! F2 <- F2 - h*A_32*U3
      CALL LSS_Mul_Jactr2(TMP,F)
      CALL DAXPY(N,ONE,U2,1,TMP,1) ! R2 <- U2 -Jac(Y+Z2)^t*h*sum(A_j2*U_j)
      CALL DAXPY(N,ONE,TMP,1,R2,1) ! R2 <- U2 -Jac(Y+Z2)^t*h*sum(A_j2*U_j)

      CALL SET2ZERO(N,F)
      CALL DAXPY(N,-H*rkA(1,3),U1,1,F,1) ! F3 <- -h*A_13*U1
      CALL DAXPY(N,-H*rkA(2,3),U2,1,F,1) ! F3 <- F3 - h*A_23*U2
      CALL DAXPY(N,-H*rkA(3,3),U3,1,F,1) ! F3 <- F3 - h*A_33*U3
      CALL LSS_Mul_Jactr3(TMP,F)
      CALL DAXPY(N,ONE,U3,1,TMP,1) ! R3 <- U3 -Jac(Y+Z3)^t*h*sum(A_j3*U_j)
      CALL DAXPY(N,ONE,TMP,1,R3,1) ! R3 <- U3 -Jac(Y+Z3)^t*h*sum(A_j3*U_j)

    END SUBROUTINE RK_PrepareRHS_Adj


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_Decomp(H,ISING)
   !~~~> Compute the matrices E1 and E2 and their decompositions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      
      INTEGER :: ISING
      DOUBLE PRECISION    :: H, Alpha, Beta, Gamma
!      DOUBLE PRECISION    :: FJAC(N,N),E1(N,N)
!      COMPLEX(kind=dp) :: E2(N,N)
      
      Gamma = rkGamma/H
      Alpha = rkAlpha/H
      Beta  = rkBeta /H

      CALL LSS_Decomp(Gamma, ISING)
      
      IF (ISING /= 0) THEN
         ISTATUS(Ndec) = ISTATUS(Ndec) + 1
         RETURN
      END IF
     
      CALL LSS_Decomp_Cmp(Alpha,Beta,ISING)
      ISTATUS(Ndec) = ISTATUS(Ndec) + 1
      
    END SUBROUTINE RK_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_Solve(N,H,R1,R2,R3,ISING)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      INTEGER :: N,ISING
      DOUBLE PRECISION    :: R1(N),R2(N),R3(N)
      DOUBLE PRECISION    :: H, x1, x2, x3
      INTEGER :: i
      LOGICAL :: Transp
!      
     ! Z <- h^{-1) T^{-1) A^{-1) x Z
      DO i=1,N
          x1 = R1(i)/H; x2 = R2(i)/H; x3 = R3(i)/H
          R1(i) = rkTinvAinv(1,1)*x1 + rkTinvAinv(1,2)*x2 + rkTinvAinv(1,3)*x3
          R2(i) = rkTinvAinv(2,1)*x1 + rkTinvAinv(2,2)*x2 + rkTinvAinv(2,3)*x3
          R3(i) = rkTinvAinv(3,1)*x1 + rkTinvAinv(3,2)*x2 + rkTinvAinv(3,3)*x3
      END DO
      Transp = .FALSE.
      CALL LSS_Solve(Transp,R1)
      CALL LSS_Solve_CMP(Transp, R2, R3)

      ! Z <- T x Z
      DO i=1,N
          x1 = R1(i); x2 = R2(i); x3 = R3(i)
          R1(i) = rkT(1,1)*x1 + rkT(1,2)*x2 + rkT(1,3)*x3
          R2(i) = rkT(2,1)*x1 + rkT(2,2)*x2 + rkT(2,3)*x3
          R3(i) = rkT(3,1)*x1 + rkT(3,2)*x2 + rkT(3,3)*x3
      END DO

      ISTATUS(Nsol) = ISTATUS(Nsol) + 2

    END SUBROUTINE RK_Solve


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_SolveTR(N,H,R1,R2,R3,ISING)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(INOUT) :: ISING
      
      DOUBLE PRECISION, INTENT(INOUT) :: R1(N),R2(N),R3(N)
      DOUBLE PRECISION    :: H, x1, x2, x3
      INTEGER :: i
      LOGICAL :: Transp
!      
     ! RHS <- h^{-1) (A^{-1) T^{-1))^t x RHS
      DO i=1,N
          x1 = R1(i)/H; x2 = R2(i)/H; x3 = R3(i)/H
          R1(i) = rkAinvT(1,1)*x1 + rkAinvT(2,1)*x2 + rkAinvT(3,1)*x3
          R2(i) = rkAinvT(1,2)*x1 + rkAinvT(2,2)*x2 + rkAinvT(3,2)*x3
          R3(i) = rkAinvT(1,3)*x1 + rkAinvT(2,3)*x2 + rkAinvT(3,3)*x3
      END DO
      Transp = .TRUE.
      CALL LSS_Solve(Transp,R1)
      
      R3(:) = -R3(:)
      CALL LSS_Solve_CMP(Transp,R2,R3)      
      R3(:) = -R3(:)      

      ! X <- (T^{-1})^t x X
      DO i=1,N
          x1 = R1(i); x2 = R2(i); x3 = R3(i)
          R1(i) = rkTinv(1,1)*x1 + rkTinv(2,1)*x2 + rkTinv(3,1)*x3
          R2(i) = rkTinv(1,2)*x1 + rkTinv(2,2)*x2 + rkTinv(3,2)*x3
          R3(i) = rkTinv(1,3)*x1 + rkTinv(2,3)*x2 + rkTinv(3,3)*x3
      END DO

      ISTATUS(Nsol) = ISTATUS(Nsol) + 2

    END SUBROUTINE RK_SolveTR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_ErrorEstimate(N,H,Y,T,Z1,Z2,Z3,SCAL,Err,FirstStep,Reject)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      
      INTEGER :: N
      DOUBLE PRECISION :: SCAL(N), Z1(N), Z2(N), Z3(N), F1(N), F2(N), F0(N),  &
                          Y(N), TMP(N), T, H
      INTEGER :: i
      LOGICAL FirstStep,Reject,Transp
      DOUBLE PRECISION :: HEE1,HEE2,HEE3,Err

      HEE1  = rkE(1)/H
      HEE2  = rkE(2)/H
      HEE3  = rkE(3)/H

      CALL FUN(N,T,Y,F0)
      ISTATUS(Nfun) = ISTATUS(Nfun) + 1

      DO  i=1,N
         F2(i)  = HEE1*Z1(i)+HEE2*Z2(i)+HEE3*Z3(i)
         TMP(i) = rkE(0)*F0(i) + F2(i)
      END DO
      Transp = .FALSE.
      CALL LSS_Solve(Transp,TMP)
      ISTATUS(Nsol) = ISTATUS(Nsol) + 1
      IF ((rkMethod==R1A).OR.(rkMethod==GAU).OR.(rkMethod==L3A)) then
       CALL LSS_Solve(Transp,TMP)
      END IF
      IF (rkMethod==GAU) CALL LSS_Solve(Transp,TMP)

      Err = RK_ErrorNorm(N,SCAL,TMP)
!
      IF (Err < ONE) RETURN
fNrej:IF (FirstStep.OR.Reject) THEN
          DO i=1,N
             TMP(i)=Y(i)+TMP(i)
          END DO
          CALL FUN(N,T,TMP,F1)
          ISTATUS(Nfun) = ISTATUS(Nfun) + 1
          DO i=1,N
             TMP(i)=F1(i)+F2(i)
          END DO

          CALL LSS_Solve(Transp,TMP)     
          ISTATUS(Nsol) = ISTATUS(Nsol) + 1
          Err = RK_ErrorNorm(N,SCAL,TMP)
       END IF fNrej
 
    END SUBROUTINE RK_ErrorEstimate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOUBLE PRECISION FUNCTION RK_ErrorNorm(N,SCAL,DY)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER :: N
      DOUBLE PRECISION :: SCAL(N),DY(N)
      INTEGER :: i

      RK_ErrorNorm = ZERO
      DO i=1,N
          RK_ErrorNorm = RK_ErrorNorm + (DY(i)*SCAL(i))**2
      END DO
      RK_ErrorNorm = MAX( SQRT(RK_ErrorNorm/N), 1.0d-10 )
 
    END FUNCTION RK_ErrorNorm
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Radau2A_Coefficients
!    The coefficients of the 3-stage Radau-2A method
!    (given to ~30 accurate digits)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
! The coefficients of the Radau2A method
      DOUBLE PRECISION :: b0

!      b0 = 1.0d0
      IF (SdirkError) THEN
        b0 = 0.2d-1
      ELSE
        b0 = 0.5d-1
      END IF

! The coefficients of the Radau2A method
      rkMethod = R2A

      rkA(1,1) =  1.968154772236604258683861429918299d-1
      rkA(1,2) = -6.55354258501983881085227825696087d-2
      rkA(1,3) =  2.377097434822015242040823210718965d-2
      rkA(2,1) =  3.944243147390872769974116714584975d-1
      rkA(2,2) =  2.920734116652284630205027458970589d-1
      rkA(2,3) = -4.154875212599793019818600988496743d-2
      rkA(3,1) =  3.764030627004672750500754423692808d-1
      rkA(3,2) =  5.124858261884216138388134465196080d-1
      rkA(3,3) =  1.111111111111111111111111111111111d-1

      rkB(1) = 3.764030627004672750500754423692808d-1
      rkB(2) = 5.124858261884216138388134465196080d-1
      rkB(3) = 1.111111111111111111111111111111111d-1

      rkC(1) = 1.550510257216821901802715925294109d-1
      rkC(2) = 6.449489742783178098197284074705891d-1
      rkC(3) = 1.0d0
      
      ! New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
      rkD(1) = 0.0d0
      rkD(2) = 0.0d0
      rkD(3) = 1.0d0

      ! Classical error estimator: 
      ! H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) = 1.0d0*b0
      rkE(1) = -10.04880939982741556246032950764708d0*b0
      rkE(2) = 1.382142733160748895793662840980412d0*b0
      rkE(3) = -.3333333333333333333333333333333333d0*b0

      ! Sdirk error estimator
      rkBgam(0) = b0
      rkBgam(1) = .3764030627004672750500754423692807d0-1.558078204724922382431975370686279d0*b0
      rkBgam(2) = .8914115380582557157653087040196118d0*b0+.5124858261884216138388134465196077d0
      rkBgam(3) = -.1637777184845662566367174924883037d0-.3333333333333333333333333333333333d0*b0
      rkBgam(4) = .2748888295956773677478286035994148d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -1.520677486405081647234271944611547d0-10.04880939982741556246032950764708d0*b0
      rkTheta(2) = 2.070455145596436382729929151810376d0+1.382142733160748895793662840980413d0*b0
      rkTheta(3) = -.3333333333333333333333333333333333d0*b0-.3744441479783868387391430179970741d0

      ! Local order of error estimator 
      IF (b0==0.0d0) THEN
        rkELO  = 6.0d0
      ELSE
        rkELO  = 4.0d0
      END IF

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 3.637834252744495732208418513577775d0
      rkAlpha = 2.681082873627752133895790743211112d0
      rkBeta  = 3.050430199247410569426377624787569d0

      rkT(1,1) =  9.443876248897524148749007950641664d-2
      rkT(1,2) = -1.412552950209542084279903838077973d-1
      rkT(1,3) = -3.00291941051474244918611170890539d-2
      rkT(2,1) =  2.502131229653333113765090675125018d-1
      rkT(2,2) =  2.041293522937999319959908102983381d-1
      rkT(2,3) =  3.829421127572619377954382335998733d-1
      rkT(3,1) =  1.0d0
      rkT(3,2) =  1.0d0
      rkT(3,3) =  0.0d0

      rkTinv(1,1) =  4.178718591551904727346462658512057d0
      rkTinv(1,2) =  3.27682820761062387082533272429617d-1
      rkTinv(1,3) =  5.233764454994495480399309159089876d-1
      rkTinv(2,1) = -4.178718591551904727346462658512057d0
      rkTinv(2,2) = -3.27682820761062387082533272429617d-1
      rkTinv(2,3) =  4.766235545005504519600690840910124d-1
      rkTinv(3,1) = -5.02872634945786875951247343139544d-1
      rkTinv(3,2) =  2.571926949855605429186785353601676d0
      rkTinv(3,3) = -5.960392048282249249688219110993024d-1

      rkTinvAinv(1,1) =  1.520148562492775501049204957366528d+1
      rkTinvAinv(1,2) =  1.192055789400527921212348994770778d0
      rkTinvAinv(1,3) =  1.903956760517560343018332287285119d0
      rkTinvAinv(2,1) = -9.669512977505946748632625374449567d0
      rkTinvAinv(2,2) = -8.724028436822336183071773193986487d0
      rkTinvAinv(2,3) =  3.096043239482439656981667712714881d0
      rkTinvAinv(3,1) = -1.409513259499574544876303981551774d+1
      rkTinvAinv(3,2) =  5.895975725255405108079130152868952d0
      rkTinvAinv(3,3) = -1.441236197545344702389881889085515d-1

      rkAinvT(1,1) = .3435525649691961614912493915818282d0
      rkAinvT(1,2) = -.4703191128473198422370558694426832d0
      rkAinvT(1,3) = .3503786597113668965366406634269080d0
      rkAinvT(2,1) = .9102338692094599309122768354288852d0
      rkAinvT(2,2) = 1.715425895757991796035292755937326d0
      rkAinvT(2,3) = .4040171993145015239277111187301784d0
      rkAinvT(3,1) = 3.637834252744495732208418513577775d0
      rkAinvT(3,2) = 2.681082873627752133895790743211112d0
      rkAinvT(3,3) = -3.050430199247410569426377624787569d0

    END SUBROUTINE Radau2A_Coefficients

    

    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Lobatto3C_Coefficients
!    The coefficients of the 3-stage Lobatto-3C method
!    (given to ~30 accurate digits)
!    The parameter b0 can be chosen to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
      DOUBLE PRECISION :: b0

      rkMethod = L3C

!      b0 = 1.0d0
      IF (SdirkError) THEN
        b0 = 0.2d0
      ELSE
        b0 = 0.5d0
      END IF
! The coefficients of the Lobatto3C method

      rkA(1,1) =  .1666666666666666666666666666666667d0
      rkA(1,2) = -.3333333333333333333333333333333333d0
      rkA(1,3) =  .1666666666666666666666666666666667d0
      rkA(2,1) =  .1666666666666666666666666666666667d0
      rkA(2,2) =  .4166666666666666666666666666666667d0
      rkA(2,3) = -.8333333333333333333333333333333333d-1
      rkA(3,1) =  .1666666666666666666666666666666667d0
      rkA(3,2) =  .6666666666666666666666666666666667d0
      rkA(3,3) =  .1666666666666666666666666666666667d0

      rkB(1) = .1666666666666666666666666666666667d0
      rkB(2) = .6666666666666666666666666666666667d0
      rkB(3) = .1666666666666666666666666666666667d0

      rkC(1) = 0.0d0
      rkC(2) = 0.5d0
      rkC(3) = 1.0d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) = .16666666666666666666666666666666667d0-b0
      rkBhat(2) = .66666666666666666666666666666666667d0
      rkBhat(3) = .16666666666666666666666666666666667d0
      
      ! New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
      rkD(1) = 0.0d0
      rkD(2) = 0.0d0
      rkD(3) = 1.0d0

      ! Classical error estimator: 
      !   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) =   .3808338772072650364017425226487022*b0
      rkE(1) = -1.142501631621795109205227567946107*b0
      rkE(2) = -1.523335508829060145606970090594809*b0
      rkE(3) =   .3808338772072650364017425226487022*b0

      ! Sdirk error estimator
      rkBgam(0) = b0
      rkBgam(1) = .1666666666666666666666666666666667d0-1.d0*b0
      rkBgam(2) = .6666666666666666666666666666666667d0
      rkBgam(3) = -.2141672105405983697350758559820354d0
      rkBgam(4) = .3808338772072650364017425226487021d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -3.d0*b0-.3808338772072650364017425226487021d0
      rkTheta(2) = -4.d0*b0+1.523335508829060145606970090594808d0
      rkTheta(3) = -.142501631621795109205227567946106d0+b0

      ! Local order of error estimator 
      IF (b0==0.0d0) THEN
        rkELO  = 5.0d0
      ELSE
        rkELO  = 4.0d0
      END IF

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 2.625816818958466716011888933765284d0
      rkAlpha = 1.687091590520766641994055533117359d0
      rkBeta  = 2.508731754924880510838743672432351d0

      rkT(1,1) = 1.d0
      rkT(1,2) = 1.d0
      rkT(1,3) = 0.d0
      rkT(2,1) = .4554100411010284672111720348287483d0
      rkT(2,2) = -.6027050205505142336055860174143743d0
      rkT(2,3) = -.4309321229203225731070721341350346d0
      rkT(3,1) = 2.195823345445647152832799205549709d0
      rkT(3,2) = -1.097911672722823576416399602774855d0
      rkT(3,3) = .7850032632435902184104551358922130d0

      rkTinv(1,1) = .4205559181381766909344950150991349d0
      rkTinv(1,2) = .3488903392193734304046467270632057d0
      rkTinv(1,3) = .1915253879645878102698098373933487d0
      rkTinv(2,1) = .5794440818618233090655049849008650d0
      rkTinv(2,2) = -.3488903392193734304046467270632057d0
      rkTinv(2,3) = -.1915253879645878102698098373933487d0
      rkTinv(3,1) = -.3659705575742745254721332009249516d0
      rkTinv(3,2) = -1.463882230297098101888532803699806d0
      rkTinv(3,3) = .4702733607340189781407813565524989d0

      rkTinvAinv(1,1) = 1.104302803159744452668648155627548d0
      rkTinvAinv(1,2) = .916122120694355522658740710823143d0
      rkTinvAinv(1,3) = .5029105849749601702795812241441172d0
      rkTinvAinv(2,1) = 1.895697196840255547331351844372453d0
      rkTinvAinv(2,2) = 3.083877879305644477341259289176857d0
      rkTinvAinv(2,3) = -1.502910584974960170279581224144117d0
      rkTinvAinv(3,1) = .8362439183082935036129145574774502d0
      rkTinvAinv(3,2) = -3.344975673233174014451658229909802d0
      rkTinvAinv(3,3) = .312908409479233358005944466882642d0

      rkAinvT(1,1) = 2.625816818958466716011888933765282d0
      rkAinvT(1,2) = 1.687091590520766641994055533117358d0
      rkAinvT(1,3) = -2.508731754924880510838743672432351d0
      rkAinvT(2,1) = 1.195823345445647152832799205549710d0
      rkAinvT(2,2) = -2.097911672722823576416399602774855d0
      rkAinvT(2,3) = .7850032632435902184104551358922130d0
      rkAinvT(3,1) = 5.765829871932827589653709477334136d0
      rkAinvT(3,2) = .1170850640335862051731452613329320d0
      rkAinvT(3,3) = 4.078738281412060947659653944216779d0

    END SUBROUTINE Lobatto3C_Coefficients

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Gauss_Coefficients
!    The coefficients of the 3-stage Gauss method
!    (given to ~30 accurate digits)
!    The parameter b3 can be chosen by the user
!    to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
      DOUBLE PRECISION :: b0
! The coefficients of the Gauss method


      rkMethod = GAU
      
!      b0 = 4.0d0
      b0 = 0.1d0
      
! The coefficients of the Gauss method

      rkA(1,1) =  .1388888888888888888888888888888889d0
      rkA(1,2) = -.359766675249389034563954710966045d-1
      rkA(1,3) =  .97894440153083260495800422294756d-2
      rkA(2,1) =  .3002631949808645924380249472131556d0
      rkA(2,2) =  .2222222222222222222222222222222222d0
      rkA(2,3) = -.224854172030868146602471694353778d-1
      rkA(3,1) =  .2679883337624694517281977355483022d0
      rkA(3,2) =  .4804211119693833479008399155410489d0
      rkA(3,3) =  .1388888888888888888888888888888889d0

      rkB(1) = .2777777777777777777777777777777778d0
      rkB(2) = .4444444444444444444444444444444444d0
      rkB(3) = .2777777777777777777777777777777778d0

      rkC(1) = .1127016653792583114820734600217600d0
      rkC(2) = .5000000000000000000000000000000000d0
      rkC(3) = .8872983346207416885179265399782400d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) =-1.4788305577012361475298775666303999d0*b0 &
                  +.27777777777777777777777777777777778d0
      rkBhat(2) =  .44444444444444444444444444444444444d0 &
                  +.66666666666666666666666666666666667d0*b0
      rkBhat(3) = -.18783610896543051913678910003626672d0*b0 &
                  +.27777777777777777777777777777777778d0

      ! New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
      rkD(1) = .1666666666666666666666666666666667d1
      rkD(2) = -.1333333333333333333333333333333333d1
      rkD(3) = .1666666666666666666666666666666667d1

      ! Classical error estimator: 
      !   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) = .2153144231161121782447335303806954d0*b0
      rkE(1) = -2.825278112319014084275808340593191d0*b0
      rkE(2) = .2870858974881495709929780405075939d0*b0
      rkE(3) = -.4558086256248162565397206448274867d-1*b0

      ! Sdirk error estimator
      rkBgam(0) = 0.d0
      rkBgam(1) = .2373339543355109188382583162660537d0
      rkBgam(2) = .5879873931885192299409334646982414d0
      rkBgam(3) = -.4063577064014232702392531134499046d-1
      rkBgam(4) = .2153144231161121782447335303806955d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -2.594040933093095272574031876464493d0
      rkTheta(2) = 1.824611539036311947589425112250199d0
      rkTheta(3) = .1856563166634371860478043996459493d0

      ! ELO = local order of classical error estimator 
      rkELO = 4.0d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 4.644370709252171185822941421408064d0
      rkAlpha = 3.677814645373914407088529289295970d0
      rkBeta  = 3.508761919567443321903661209182446d0

      rkT(1,1) =  .7215185205520017032081769924397664d-1
      rkT(1,2) = -.8224123057363067064866206597516454d-1
      rkT(1,3) = -.6012073861930850173085948921439054d-1
      rkT(2,1) =  .1188325787412778070708888193730294d0
      rkT(2,2) =  .5306509074206139504614411373957448d-1
      rkT(2,3) =  .3162050511322915732224862926182701d0
      rkT(3,1) = 1.d0
      rkT(3,2) = 1.d0
      rkT(3,3) = 0.d0

      rkTinv(1,1) =  5.991698084937800775649580743981285d0
      rkTinv(1,2) =  1.139214295155735444567002236934009d0
      rkTinv(1,3) =   .4323121137838583855696375901180497d0
      rkTinv(2,1) = -5.991698084937800775649580743981285d0
      rkTinv(2,2) = -1.139214295155735444567002236934009d0
      rkTinv(2,3) =   .5676878862161416144303624098819503d0
      rkTinv(3,1) = -1.246213273586231410815571640493082d0
      rkTinv(3,2) =  2.925559646192313662599230367054972d0
      rkTinv(3,3) =  -.2577352012734324923468722836888244d0

      rkTinvAinv(1,1) =  27.82766708436744962047620566703329d0
      rkTinvAinv(1,2) =   5.290933503982655311815946575100597d0
      rkTinvAinv(1,3) =   2.007817718512643701322151051660114d0
      rkTinvAinv(2,1) = -17.66368928942422710690385180065675d0
      rkTinvAinv(2,2) = -14.45491129892587782538830044147713d0
      rkTinvAinv(2,3) =   2.992182281487356298677848948339886d0
      rkTinvAinv(3,1) = -25.60678350282974256072419392007303d0
      rkTinvAinv(3,2) =   6.762434375611708328910623303779923d0
      rkTinvAinv(3,3) =   1.043979339483109825041215970036771d0
      
      rkAinvT(1,1) = .3350999483034677402618981153470483d0
      rkAinvT(1,2) = -.5134173605009692329246186488441294d0
      rkAinvT(1,3) = .6745196507033116204327635673208923d-1
      rkAinvT(2,1) = .5519025480108928886873752035738885d0
      rkAinvT(2,2) = 1.304651810077110066076640761092008d0
      rkAinvT(2,3) = .9767507983414134987545585703726984d0
      rkAinvT(3,1) = 4.644370709252171185822941421408064d0
      rkAinvT(3,2) = 3.677814645373914407088529289295970d0
      rkAinvT(3,3) = -3.508761919567443321903661209182446d0
      
    END SUBROUTINE Gauss_Coefficients



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Radau1A_Coefficients
!    The coefficients of the 3-stage Gauss method
!    (given to ~30 accurate digits)
!    The parameter b3 can be chosen by the user
!    to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
!      DOUBLE PRECISION :: b0 = 0.3d0
      DOUBLE PRECISION :: b0 = 0.1d0

! The coefficients of the Radau1A method

      rkMethod = R1A

      rkA(1,1) =  .1111111111111111111111111111111111d0
      rkA(1,2) = -.1916383190435098943442935597058829d0
      rkA(1,3) =  .8052720793239878323318244859477174d-1
      rkA(2,1) =  .1111111111111111111111111111111111d0
      rkA(2,2) =  .2920734116652284630205027458970589d0
      rkA(2,3) = -.481334970546573839513422644787591d-1
      rkA(3,1) =  .1111111111111111111111111111111111d0
      rkA(3,2) =  .5370223859435462728402311533676479d0
      rkA(3,3) =  .1968154772236604258683861429918299d0

      rkB(1) = .1111111111111111111111111111111111d0
      rkB(2) = .5124858261884216138388134465196080d0
      rkB(3) = .3764030627004672750500754423692808d0

      rkC(1) = 0.d0
      rkC(2) = .3550510257216821901802715925294109d0
      rkC(3) = .8449489742783178098197284074705891d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) = .11111111111111111111111111111111111d0-b0
      rkBhat(2) = .51248582618842161383881344651960810d0
      rkBhat(3) = .37640306270046727505007544236928079d0

      ! New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
      rkD(1) = .3333333333333333333333333333333333d0
      rkD(2) = -.8914115380582557157653087040196127d0
      rkD(3) = .1558078204724922382431975370686279d1

      ! Classical error estimator: 
      ! H* Sum (b_j-bhat_j) f(Z_j) = H*E(0)*F(0) + Sum E_j Z_j
      rkE(0) =   .2748888295956773677478286035994148d0*b0
      rkE(1) = -1.374444147978386838739143017997074d0*b0
      rkE(2) = -1.335337922441686804550326197041126d0*b0
      rkE(3) =   .235782604058977333559011782643466d0*b0

      ! Sdirk error estimator
      rkBgam(0) = 0.0d0
      rkBgam(1) = .1948150124588532186183490991130616d-1
      rkBgam(2) = .7575249005733381398986810981093584d0
      rkBgam(3) = -.518952314149008295083446116200793d-1
      rkBgam(4) = .2748888295956773677478286035994148d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -1.224370034375505083904362087063351d0
      rkTheta(2) = .9340045331532641409047527962010133d0
      rkTheta(3) = .4656990124352088397561234800640929d0

      ! ELO = local order of classical error estimator 
      rkELO = 4.0d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 3.637834252744495732208418513577775d0
      rkAlpha = 2.681082873627752133895790743211112d0
      rkBeta  = 3.050430199247410569426377624787569d0

      rkT(1,1) =  .424293819848497965354371036408369d0
      rkT(1,2) = -.3235571519651980681202894497035503d0
      rkT(1,3) = -.522137786846287839586599927945048d0
      rkT(2,1) =  .57594609499806128896291585429339d-1
      rkT(2,2) =  .3148663231849760131614374283783d-2
      rkT(2,3) =  .452429247674359778577728510381731d0
      rkT(3,1) = 1.d0
      rkT(3,2) = 1.d0
      rkT(3,3) = 0.d0

      rkTinv(1,1) = 1.233523612685027760114769983066164d0
      rkTinv(1,2) = 1.423580134265707095505388133369554d0
      rkTinv(1,3) = .3946330125758354736049045150429624d0
      rkTinv(2,1) = -1.233523612685027760114769983066164d0
      rkTinv(2,2) = -1.423580134265707095505388133369554d0
      rkTinv(2,3) = .6053669874241645263950954849570376d0
      rkTinv(3,1) = -.1484438963257383124456490049673414d0
      rkTinv(3,2) = 2.038974794939896109682070471785315d0
      rkTinv(3,3) = -.544501292892686735299355831692542d-1

      rkTinvAinv(1,1) =  4.487354449794728738538663081025420d0
      rkTinvAinv(1,2) =  5.178748573958397475446442544234494d0
      rkTinvAinv(1,3) =  1.435609490412123627047824222335563d0
      rkTinvAinv(2,1) = -2.854361287939276673073807031221493d0
      rkTinvAinv(2,2) = -1.003648660720543859000994063139137d+1
      rkTinvAinv(2,3) =  1.789135380979465422050817815017383d0
      rkTinvAinv(3,1) = -4.160768067752685525282947313530352d0
      rkTinvAinv(3,2) =  1.124128569859216916690209918405860d0
      rkTinvAinv(3,3) =  1.700644430961823796581896350418417d0

      rkAinvT(1,1) = 1.543510591072668287198054583233180d0
      rkAinvT(1,2) = -2.460228411937788329157493833295004d0
      rkAinvT(1,3) = -.412906170450356277003910443520499d0
      rkAinvT(2,1) = .209519643211838264029272585946993d0
      rkAinvT(2,2) = 1.388545667194387164417459732995766d0
      rkAinvT(2,3) = 1.20339553005832004974976023130002d0
      rkAinvT(3,1) = 3.637834252744495732208418513577775d0
      rkAinvT(3,2) = 2.681082873627752133895790743211112d0
      rkAinvT(3,3) = -3.050430199247410569426377624787569d0

    END SUBROUTINE Radau1A_Coefficients

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE RungeKuttaADJ1 ! and all its internal procedures
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE RungeKuttaADJ2( N, NP, NNZERO, Y, NADJ, Lambda, Tstart, Tend,    &
               RelTol, AbsTol, RelTol_adj, AbsTol_adj, FUN, JAC, AdjInit,     &
               RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR, DRDY, QFUN, Q )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  This implementation is based on the book and the code Radau5:
!
!         E. HAIRER AND G. WANNER
!         "SOLVING ORDINARY DIFFERENTIAL EQUATIONS II. 
!              STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS."
!         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14,
!         SPRINGER-VERLAG (1991)
!
!         UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
!         CH-1211 GENEVE 24, SWITZERLAND
!         E-MAIL:  HAIRER@DIVSUN.UNIGE.CH,  WANNER@DIVSUN.UNIGE.CH
!
!   Methods:
!          * Radau-2A   quadrature (order 5)                              
!          * Radau-1A   quadrature (order 5)                              
!          * Lobatto-3C quadrature (order 4)                              
!          * Gauss      quadrature (order 6)                              
!                                                                         
!   (C)  Adrian Sandu, August 2005                                       
!   Virginia Polytechnic Institute and State University                  
!   Contact: sandu@cs.vt.edu                                             
!   Revised by Philipp Miehe and Adrian Sandu, May 2006                  
!   Revised by Hong Zhang and Adrian Sandu, Feb 2011
!   This implementation is part of FATODE
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!       ----------------
!
!    Note: For input parameters equal to zero the default values of the
!          corresponding variables are used.
!
!     N           Dimension of the system
!     T           Initial time value
!
!     Tend        Final T value (Tend-T may be positive or negative)
!
!     Y(N)        Initial values for Y
!
!     Q(NADJ)     Initial values for Q
!
!     RelTol,AbsTol   Relative and absolute error tolerances. 
!          for ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
!                        = 1: AbsTol, RelTol are scalars
!
!~~~>  Integer input parameters:
!  
!    ICNTRL(1) = not used
!
!    ICNTRL(2) = 0: AbsTol, RelTol are NVAR-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3) = RK method selection       
!              = 1:  Radau-2A    (the default)
!              = 2:  Lobatto-3C
!              = 3:  Gauss
!              = 4:  Radau-1A
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0 the default value of 10000 is used
!
!    ICNTRL(5)  -> maximum number of Newton iterations
!        For ICNTRL(5)=0 the default value of 8 is used
!
!    ICNTRL(6)  -> starting values of Newton iterations:
!        ICNTRL(6)=0 : starting values are obtained from 
!                      the extrapolated collocation solution
!                      (the default)
!        ICNTRL(6)=1 : starting values are zero
!
!    ICNTRL(7)  -> method to solve the linear ADJ equations:
!        ICNTRL(7)=0,1 : modified Newton re-using LU (the default)
!                        with a fixed number of iterations
!        ICNTRL(7)=2 :   direct solution (additional one LU factorization
!                        of 3Nx3N matrix per step); good for debugging
!        ICNTRL(7)=3 :   adaptive solution (if Newton does not converge
!                        switch to direct)
!
!    ICNTRL(8)  -> checkpointing the LU factorization at each step:
!        ICNTRL(8)=0 : do *not* save LU factorization (the default)
!        ICNTRL(8)=1 : save LU factorization
!        Note: if ICNTRL(7)=1 the LU factorization is *not* saved
!
!    ICNTRL(9) -> Type of adjoint algorithm
!         = 0 : default is discrete adjoint ( of method ICNTRL(3) )
!         = 1 : no adjoint       
!         = 2 : discrete adjoint ( of method ICNTRL(3) )
!         = 3 : fully adaptive continuous adjoint ( with method ICNTRL(6) )
!         = 4 : simplified continuous adjoint ( with method ICNTRL(6) )
!
!    ICNTRL(10) -> switch for error estimation strategy
!               ICNTRL(10) = 0: one additional stage at c=0, 
!                               see Hairer (default)
!               ICNTRL(10) = 1: two additional stages at c=0 
!                               and SDIRK at c=1, stiffly accurate
!
!    ICNTRL(11) -> switch for step size strategy
!              ICNTRL(11)=1:  mod. predictive controller (Gustafsson, default)
!              ICNTRL(11)=2:  classical step size control
!              the choice 1 seems to produce safer results;
!              for simple problems, the choice 2 produces
!              often slightly faster runs
!
!~~~>  Real input parameters:
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!                  (highly recommended to keep Hmin = ZERO, the default)
!
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!
!    RCNTRL(3)  -> Hstart, the starting step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                 (default=0.1)
!
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!                  than the predicted value  (default=0.9)
!
!    RCNTRL(8)  -> ThetaMin. If Newton convergence rate smaller
!                  than ThetaMin the Jacobian is not recomputed;
!                  (default=0.001)
!
!    RCNTRL(9)  -> NewtonTol, stopping criterion for Newton's method
!                  (default=0.03)
!
!    RCNTRL(10) -> Qmin
!
!    RCNTRL(11) -> Qmax. If Qmin < Hnew/Hold < Qmax, then the
!                  step size is kept constant and the LU factorization
!                  reused (default Qmin=1, Qmax=1.2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!    OUTPUT ARGUMENTS:
!    -----------------
!
!    T           -> T value for which the solution has been computed
!                     (after successful return T=Tend).
!
!    Y(N)        -> Numerical solution at T
!
!    Q(NADJ)     -> Numerical solution of the quadrature term at T
!    IERR        -> Reports on successfulness upon return:
!                    = 1 for success
!                    < 0 for error (value equals error code)
!
!    ISTATUS(1)  -> No. of function calls
!    ISTATUS(2)  -> No. of Jacobian calls
!    ISTATUS(3)  -> No. of steps
!    ISTATUS(4)  -> No. of accepted steps
!    ISTATUS(5)  -> No. of rejected steps (except at very beginning)
!    ISTATUS(6)  -> No. of LU decompositions
!    ISTATUS(7)  -> No. of forward/backward substitutions
!    ISTATUS(8)  -> No. of singular matrix decompositions
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                     computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!                   For multiple restarts, use Hnew as Hstart 
!                     in the subsequent run
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USE LS_Solver
    IMPLICIT NONE
      
    INTEGER, INTENT(IN)     :: N,NNZERO,NP,NADJ
    DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ)
    DOUBLE PRECISION, INTENT(IN) :: AbsTol(N),RelTol(N), AbsTol_adj(N,NADJ),  &
                                      RelTol_adj(N,NADJ)
    DOUBLE PRECISION, INTENT(INOUT) :: Y(N),RCNTRL(20),RSTATUS(20)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: Q(NADJ)
    INTEGER :: ICNTRL(20), ISTATUS(20)
    LOGICAL :: StartNewton, Gustafsson, SdirkError, SaveLU, GetQuad
    INTEGER :: IERR, ITOL
    DOUBLE PRECISION, INTENT(IN):: Tstart,Tend
    DOUBLE PRECISION :: Texit

    !~~~> Control arguments
    INTEGER :: Max_no_steps, NewtonMaxit, AdjointType, rkMethod
    DOUBLE PRECISION :: Hmin,Hmax,Hstart,Qmin,Qmax
    DOUBLE PRECISION :: Roundoff, ThetaMin, NewtonTol
    DOUBLE PRECISION :: FacSafe,FacMin,FacMax,FacRej
    ! Runge-Kutta method parameters
    INTEGER, PARAMETER :: R2A=1, R1A=2, L3C=3, GAU=4, L3A=5
    DOUBLE PRECISION :: rkT(3,3),rkTinv(3,3),rkTinvAinv(3,3), rkAinvT(3,3),   &
                        rkA(3,3), rkB(3),  rkC(3), rkD(0:3), rkE(0:3),      &
                        rkBgam(0:4), rkBhat(0:4), rkTheta(0:3),             &
                       rkGamma,  rkAlpha, rkBeta, rkELO
    ! ADJ method parameters
    INTEGER, PARAMETER :: RK_adj_f90_none = 1, RK_adj_f90_discrete = 2,       &
                     RK_adj_f90_continuous = 3,RK_adj_f90_simple_continuous = 4
    INTEGER :: AdjointSolve                      
    INTEGER, PARAMETER :: Solve_direct = 1, Solve_fixed = 2, Solve_adaptive = 3
    INTEGER :: stack_ptr = 0 ! last written entry
    DOUBLE PRECISION, DIMENSION(:), POINTER :: chk_H, chk_T
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_Y, chk_Z
    INTEGER, DIMENSION(:), POINTER :: chk_NiT
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_dY, chk_d2Y
    !~~~> Local variables
    INTEGER :: i 
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE = 1.0d0
    DOUBLE PRECISION :: DLAMCH
    EXTERNAL QFUN,FUN,JAC,AdjInit,DRDY
    OPTIONAL QFUN,DRDY

    GetQuad = .FALSE.
    IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDY)) GetQuad = .TRUE.
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        SETTING THE PARAMETERS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IERR = 0
    ISTATUS(1:20) = 0
    RSTATUS(1:20) = ZERO
!~~~> ICNTRL(1) - autonomous system - not used       
!~~~> ITOL: 1 for vector and 0 for scalar AbsTol/RelTol
    IF (ICNTRL(2) == 0) THEN
      ITOL = 1
    ELSE
      ITOL = 0
    END IF
!~~~> Error control selection  
    IF (ICNTRL(10) == 0) THEN 
      SdirkError = .FALSE.
    ELSE
      SdirkError = .TRUE.
    END IF      
!~~~> Method selection  
    SELECT CASE (ICNTRL(3))     
    CASE (1)
      CALL Radau2A_Coefficients
    CASE (0,2)
      CALL Lobatto3C_Coefficients
    CASE (3)
      CALL Gauss_Coefficients
    CASE (4)
      CALL Radau1A_Coefficients
    CASE DEFAULT
      WRITE(6,*) 'ICNTRL(3)=',ICNTRL(3)
      CALL RK_ErrorMsg(-13,Tstart,ZERO,IERR)
    END SELECT
!~~~> Max_no_steps: the maximal number of time steps
    IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 10000
    ELSE
      Max_no_steps=ICNTRL(4)
      IF (Max_no_steps <= 0) THEN
        WRITE(6,*) 'ICNTRL(4)=',ICNTRL(4)
        CALL RK_ErrorMsg(-1,Tstart,ZERO,IERR)
      END IF
    END IF
!~~~> NewtonMaxit    maximal number of Newton iterations
    IF (ICNTRL(5) == 0) THEN
      NewtonMaxit = 8
    ELSE
      NewtonMaxit=ICNTRL(5)
      IF (NewtonMaxit <= 0) THEN
        WRITE(6,*) 'ICNTRL(5)=',ICNTRL(5)
        CALL RK_ErrorMsg(-2,Tstart,ZERO,IERR)
      END IF
    END IF
!~~~> StartNewton:  Use extrapolation for starting values of Newton iterations
    IF (ICNTRL(6) == 0) THEN
       StartNewton = .TRUE.
    ELSE
       StartNewton = .FALSE.
    END IF      
!~~~>  How to solve the linear adjoint system
    SELECT CASE (ICNTRL(7))     
    CASE (0,1)
      AdjointSolve = Solve_fixed
    CASE (2)
      AdjointSolve = Solve_direct
    CASE (3)
      AdjointSolve = Solve_adaptive
    CASE DEFAULT  
      PRINT * , 'User-selected adjoint solution: ICNTRL(7)=', ICNTRL(7)
      PRINT * , 'Implemented: =(0,1) (fixed), =2 (direct), =3 (adaptive)'
      CALL rk_ErrorMsg(-9,Tstart,ZERO,IERR)
      RETURN      
    END SELECT
!~~~>  Discrete or continuous adjoint formulation
    SELECT CASE (ICNTRL(9))     
    CASE (0,2)
      AdjointType = RK_adj_f90_discrete
    CASE (1)
      AdjointType = RK_adj_f90_none
    CASE DEFAULT  
      PRINT * , 'User-selected adjoint type: ICNTRL(9)=', ICNTRL(9)
      PRINT * , 'Implemented: =(0,2) (discrete adj) and =1 (no adj)'
      CALL rk_ErrorMsg(-9,Tstart,ZERO,IERR)
      RETURN      
    END SELECT
!~~~> Save or not the forward LU factorization
    SaveLU = (ICNTRL(8) /= 0)
    SaveLU = .FALSE. 
    IF (AdjointSolve == Solve_direct) SaveLU = .FALSE.
!~~~> Gustafsson: step size controller
    IF (ICNTRL(11) == 0) THEN
      Gustafsson = .TRUE.
    ELSE
      Gustafsson = .FALSE.
    END IF

!~~~> Roundoff: smallest number s.t. 1.0 + Roundoff > 1.0
    Roundoff = DLAMCH('E');

!~~~> Hmin = minimal step size
    IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
    ELSE
      Hmin = MIN(ABS(RCNTRL(1)),ABS(Tend-Tstart))
    END IF
!~~~> Hmax = maximal step size
    IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
    ELSE
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
    END IF
!~~~> Hstart = starting step size
    IF (RCNTRL(3) == ZERO) THEN
      Hstart = ZERO
    ELSE
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
    END IF
!~~~> FacMin: lower bound on step decrease factor
    IF(RCNTRL(4) == ZERO)THEN
      FacMin = 0.2d0
    ELSE
      FacMin = RCNTRL(4)
    END IF
!~~~> FacMax: upper bound on step increase factor
    IF(RCNTRL(5) == ZERO)THEN
       FacMax = 8.D0
    ELSE
       FacMax = RCNTRL(5)
    END IF
!~~~> FacRej: step decrease factor after 2 consecutive rejections
    IF(RCNTRL(6) == ZERO)THEN
      FacRej = 0.1d0
    ELSE
      FacRej = RCNTRL(6)
    END IF
!~~~> FacSafe:  by which the new step is slightly smaller than the predicted value
    IF (RCNTRL(7) == ZERO) THEN
      FacSafe=0.9d0
    ELSE
      FacSafe=RCNTRL(7)
    END IF
    IF ( (FacMax < ONE) .OR. (FacMin > ONE) .OR. &
           (FacSafe <= 1.0d-3) .OR. (FacSafe >= ONE) ) THEN
      WRITE(6,*)'RCNTRL(4:7)=',RCNTRL(4:7)
      CALL RK_ErrorMsg(-4,Tstart,ZERO,IERR)
    END IF

!~~~> ThetaMin:  decides whether the Jacobian should be recomputed
    IF (RCNTRL(8) == ZERO) THEN
      ThetaMin = 1.0d-3
    ELSE
      ThetaMin=RCNTRL(8)
      IF (ThetaMin <= 0.0d0 .OR. ThetaMin >= 1.0d0) THEN
        WRITE(6,*) 'RCNTRL(8)=', RCNTRL(8)
        CALL RK_ErrorMsg(-5,Tstart,ZERO,IERR)
      END IF
    END IF
!~~~> NewtonTol:  stopping crierion for Newton's method
    IF (RCNTRL(9) == ZERO) THEN
      NewtonTol = 3.0d-2
    ELSE
      NewtonTol = RCNTRL(9)
      IF (NewtonTol <= Roundoff) THEN
        WRITE(6,*) 'RCNTRL(9)=',RCNTRL(9)
        CALL RK_ErrorMsg(-6,Tstart,ZERO,IERR)
      END IF
    END IF
!~~~> Qmin AND Qmax: IF Qmin < Hnew/Hold < Qmax then step size = const.
    IF (RCNTRL(10) == ZERO) THEN
      Qmin=1.D0
    ELSE
      Qmin=RCNTRL(10)
    END IF
    IF (RCNTRL(11) == ZERO) THEN
      Qmax=1.2D0
    ELSE
      Qmax=RCNTRL(11)
    END IF
    IF (Qmin > ONE .OR. Qmax < ONE) THEN
      WRITE(6,*) 'RCNTRL(10:11)=',Qmin,Qmax
      CALL RK_ErrorMsg(-7,Tstart,ZERO,IERR)
    END IF
!~~~> Check if tolerances are reasonable
    IF (ITOL == 0) THEN
      IF (AbsTol(1) <= ZERO.OR.RelTol(1) <= 10.d0*Roundoff) THEN
        WRITE (6,*) 'AbsTol/RelTol=',AbsTol,RelTol 
        CALL RK_ErrorMsg(-8,Tstart,ZERO,IERR)
      END IF
    ELSE
      DO i=1,N
        IF (AbsTol(i) <= ZERO.OR.RelTol(i) <= 10.d0*Roundoff) THEN
          WRITE (6,*) 'AbsTol/RelTol(',i,')=',AbsTol(i),RelTol(i)
          CALL RK_ErrorMsg(-8,Tstart,ZERO,IERR)
        END IF
      END DO
    END IF

!~~~> Parameters are wrong
    IF (IERR < 0) RETURN

!~~~>  Allocate checkpoint space or open checkpoint files
    IF (AdjointType == RK_adj_f90_discrete) THEN
      CALL rk_AllocateDBuffers()
    ELSEIF ( (AdjointType == RK_adj_f90_continuous).OR. &
           (AdjointType == RK_adj_f90_simple_continuous) ) THEN
      CALL rk_AllocateCBuffers
    END IF

    CALL LSS_Init(N,NNZERO)
!~~~> Call the core method
    CALL RK_FwdIntegrator( N,NADJ,Tstart,Tend,Y,AdjointType,GetQuad,IERR )
    PRINT*,'FORWARD STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),' Rej=',ISTATUS(Nrej), &
           ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=', ISTATUS(Njac),    &
           ' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)
!~~~>  If Forward integration failed return   
    IF (IERR<0) RETURN

!~~~> Initialize adjoint variables
    CALL AdjInit(N,NP,NADJ,Tend,Y,Lambda)

    SELECT CASE (AdjointType)   
    CASE (RK_adj_f90_discrete)   
      CALL rk_DadjInt ( N,NADJ,Lambda,Tstart,Tend,Texit,GetQuad,IERR )
    CASE (RK_adj_f90_continuous) 
      CALL rk_CadjInt ( NADJ, Lambda, Tend, Tstart, Texit, IERR )
    CASE (RK_adj_f90_simple_continuous)
      CALL rk_SimpleCadjInt ( NADJ, Lambda, Tstart, Tend, Texit, IERR )
    END SELECT ! AdjointType
    PRINT*,'ADJOINT STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
         ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

!~~~>  Free checkpoint space or close checkpoint files
    IF (AdjointType == RK_adj_f90_discrete) THEN
      CALL rk_FreeDBuffers
    ELSEIF ( (AdjointType == RK_adj_f90_continuous) .OR. &
           (AdjointType == RK_adj_f90_simple_continuous) ) THEN
      CALL rk_FreeCBuffers
    END IF
    CALL LSS_Free
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS ! Internal procedures to RungeKuttaADJ2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_AllocateDBuffers()
!~~~>  Allocate buffer space for discrete adjoint
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
      ALLOCATE( chk_Z(N*3,Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer Z'; STOP
      END IF   
      ALLOCATE( chk_NiT(Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer NiT'; STOP
      END IF   
 
    END SUBROUTINE rk_AllocateDBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_FreeDBuffers()
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
        PRINT*,'Failed deallocation of buffer Z'; STOP
      END IF   
      DEALLOCATE( chk_NiT, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer NiT'; STOP
      END IF   

    END SUBROUTINE rk_FreeDBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_AllocateCBuffers()
!~~~>  Allocate buffer space for continuous adjoint
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
      ALLOCATE( chk_dY(N,Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer dY'; STOP
      END IF   
      ALLOCATE( chk_d2Y(N,Max_no_steps), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer d2Y'; STOP
      END IF   
 
    END SUBROUTINE rk_AllocateCBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_FreeCBuffers()
!~~~>  Dallocate buffer space for continuous adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: i
      print*,'cbuffers deallocate???'
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
      DEALLOCATE( chk_dY, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer dY'; STOP
      END IF   
      DEALLOCATE( chk_d2Y, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer d2Y'; STOP
      END IF   
 
    END SUBROUTINE rk_FreeCBuffers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_DPush( T, H, Y, Zstage, NewIt)!, Jcb )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION :: T, H, Y(N), Zstage(N*3)
      INTEGER :: NewIt

      stack_ptr = stack_ptr + 1
      IF ( stack_ptr > Max_no_steps ) THEN
        PRINT*,'Push failed: buffer overflow'
        STOP
      END IF  
      chk_H( stack_ptr ) = H
      chk_T( stack_ptr ) = T
      ! CALL DCOPY(NVAR,Y,1,chk_Y(1,stack_ptr),1)
      ! CALL DCOPY(NVAR*3,Zstage,1,chk_Z(1,stack_ptr),1)
      chk_Y(1:N,stack_ptr) = Y(1:N)
      chk_Z(1:3*N,stack_ptr) = Zstage(1:3*N)
      chk_NiT( stack_ptr ) = NewIt
  
    END SUBROUTINE rk_DPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_DPop( T, H, Y, Zstage, NewIt) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DOUBLE PRECISION :: T, H, Y(N), Zstage(N*3) ! , Jcb(LU_NONZERO)
      INTEGER :: NewIt
   
      IF ( stack_ptr <= 0 ) THEN
        PRINT*,'Pop failed: empty buffer'
        STOP
      END IF  
      H = chk_H( stack_ptr )
      T = chk_T( stack_ptr )
      ! CALL DCOPY(NVAR,chk_Y(1,stack_ptr),1,Y,1)
      Y(1:N) = chk_Y(1:N,stack_ptr)
      ! CALL DCOPY(NVAR*3,chk_Z(1,stack_ptr),1,Zstage,1)
      Zstage(1:3*N) = chk_Z(1:3*N,stack_ptr)
      NewIt = chk_NiT( stack_ptr )

      stack_ptr = stack_ptr - 1
  
    END SUBROUTINE rk_DPop
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_CPush(T, H, Y, dY, d2Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DOUBLE PRECISION :: T, H, Y(N), dY(N), d2Y(N)
      stack_ptr = stack_ptr + 1
      IF ( stack_ptr > Max_no_steps ) THEN
        PRINT*,'Push failed: buffer overflow'
        STOP
      END IF  
      chk_H( stack_ptr ) = H
      chk_T( stack_ptr ) = T
      ! CALL DCOPY(NVAR,Y,1,chk_Y(1,stack_ptr),1)
      ! CALL DCOPY(NVAR,dY,1,chk_dY(1,stack_ptr),1)
      ! CALL DCOPY(NVAR,d2Y,1,chk_d2Y(1,stack_ptr),1)
      chk_Y(1:N,stack_ptr)  =  Y(1:N)
      chk_dY(1:N,stack_ptr) = dY(1:N)
      chk_d2Y(1:N,stack_ptr)= d2Y(1:N)
  
    END SUBROUTINE rk_CPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_CPop( T, H, Y, dY, d2Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DOUBLE PRECISION :: T, H, Y(N), dY(N), d2Y(N)   
      IF ( stack_ptr <= 0 ) THEN
        PRINT*,'Pop failed: empty buffer'
        STOP
      END IF  
      H = chk_H( stack_ptr )
      T = chk_T( stack_ptr )
      ! CALL DCOPY(NVAR,chk_Y(1,stack_ptr),1,Y,1)
      ! CALL DCOPY(NVAR,chk_dY(1,stack_ptr),1,dY,1)
      ! CALL DCOPY(NVAR,chk_d2Y(1,stack_ptr),1,d2Y,1)
      Y(1:N) = chk_Y(1:N,stack_ptr)
      dY(1:N) = chk_dY(1:N,stack_ptr)
      d2Y(1:N) = chk_d2Y(1:N,stack_ptr)

      stack_ptr = stack_ptr - 1
  
    END SUBROUTINE rk_CPop


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_FwdIntegrator( N,NADJ,Tstart,Tend,Y,AdjointType,GetQuad,IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
!~~~> Arguments
      INTEGER,  INTENT(IN)            :: N,NADJ
      DOUBLE PRECISION, INTENT(IN)    :: Tend, Tstart
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
      INTEGER,  INTENT(INOUT)         :: IERR
      INTEGER,  INTENT(IN)            :: AdjointType
      LOGICAL, INTENT(IN)             :: GetQuad
!~~~> Local variables
      DOUBLE PRECISION, DIMENSION(N) :: Z1,Z2,Z3,Z4,SCAL,DZ1,DZ2,DZ3,DZ4,G,TMP,FCN0
      DOUBLE PRECISION  :: CONT(N,3), Tdirection,  H, T, Hacc, Hnew, Hold, Fac,&
                 FacGus, Theta, Err, ErrOld, NewtonRate, NewtonIncrement, &
                 Hratio, Qnewton, NewtonPredictedErr,NewtonIncrementOld, ThetaSD
      INTEGER :: NewtonIter, ISING, Nconsecutive, NewIt
      LOGICAL :: Reject, FirstStep, SkipJac, NewtonDone, SkipLU, Transp
      
      DOUBLE PRECISION, DIMENSION(:), POINTER :: Zstage

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initial setting
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (AdjointType == RK_adj_f90_discrete) THEN ! Save stage solution
        ALLOCATE(Zstage(N*3), STAT=i)
        IF (i/=0) THEN
          PRINT*,'Allocation of Zstage failed'
          STOP
        END IF
      END IF   
      T=Tstart

      Tdirection = SIGN(ONE,Tend-Tstart)
      H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , Hmax )
      IF (ABS(H) <= 10.d0*Roundoff) H = 1.0d-6
      H = SIGN(H,Tdirection)
      Hold      = H
      Reject    = .FALSE.
      FirstStep = .TRUE.
      SkipJac   = .FALSE.
      SkipLU    = .FALSE.
      Transp    = .FALSE.
      IF ((T+H*1.0001D0-Tend)*Tdirection >= ZERO) THEN
         H = Tend-T
      END IF
      Nconsecutive = 0
      CALL RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tloop: DO WHILE ( (Tend-T)*Tdirection - Roundoff > ZERO )

        IF ( .NOT.SkipLU ) THEN ! This time around skip the Jac update and LU
!~~~> Compute the Jacobian matrix
          IF ( .NOT.SkipJac ) THEN
            CALL LSS_Jac(T,Y,JAC)
            ISTATUS(Njac) = ISTATUS(Njac) + 1
          END IF
!~~~> Compute the matrices E1 and E2 and their decompositions
          CALL RK_Decomp(H,ISING)
          IF (ISING /= 0) THEN
            ISTATUS(Nsng) = ISTATUS(Nsng) + 1; Nconsecutive = Nconsecutive + 1
            IF (Nconsecutive >= 5) THEN
              CALL RK_ErrorMsg(-12,T,H,IERR); RETURN
            END IF
            H=H*0.5d0; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
            CYCLE Tloop
          ELSE
            Nconsecutive = 0    
          END IF   
        END IF ! SkipLU
   
        ISTATUS(Nstp) = ISTATUS(Nstp) + 1
        IF (ISTATUS(Nstp) > Max_no_steps) THEN
          PRINT*,'Max number of time steps is ',Max_no_steps
          CALL RK_ErrorMsg(-9,T,H,IERR); RETURN
        END IF
        IF (0.1D0*ABS(H) <= ABS(T)*Roundoff)  THEN
          CALL RK_ErrorMsg(-10,T,H,IERR); RETURN
        END IF
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Loop for the simplified Newton iterations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
      !~~~>  Starting values for Newton iteration
        IF ( FirstStep .OR. (.NOT.StartNewton) ) THEN
          CALL Set2zero(N,Z1)
          CALL Set2zero(N,Z2)
          CALL Set2zero(N,Z3)
        ELSE
          ! Evaluate quadratic polynomial
          CALL RK_Interpolate('eval',N,H,Hold,Z1,Z2,Z3,CONT)
        END IF
      
        !~~~>  Initializations for Newton iteration
        NewtonDone = .FALSE.
        Fac = 0.5d0 ! Step reduction if too many iterations
      
NewtonLoop:DO  NewtonIter = 1, NewtonMaxit  
          !~~~> Prepare the right-hand side
          CALL RK_PrepareRHS(N,T,H,Y,Z1,Z2,Z3,DZ1,DZ2,DZ3)
          !~~~> Solve the linear systems
          CALL RK_Solve( N,H,DZ1,DZ2,DZ3,ISING )
          NewtonIncrement = SQRT( ( RK_ErrorNorm(N,SCAL,DZ1)**2 + &
                               RK_ErrorNorm(N,SCAL,DZ2)**2 + &
                                RK_ErrorNorm(N,SCAL,DZ3)**2 )/3.0d0 )
            
          IF ( NewtonIter == 1 ) THEN
            Theta      = ABS(ThetaMin)
            NewtonRate = 2.0d0 
          ELSE
            Theta = NewtonIncrement/NewtonIncrementOld
            IF (Theta < 0.99d0) THEN
              NewtonRate = Theta/(ONE-Theta)
            ELSE ! Non-convergence of Newton: Theta too large
              EXIT NewtonLoop
            END IF
            IF ( NewtonIter < NewtonMaxit ) THEN 
            ! Predict error at the end of Newton process 
               NewtonPredictedErr = NewtonIncrement &
                             *Theta**(NewtonMaxit-NewtonIter)/(ONE-Theta)
               IF (NewtonPredictedErr >= NewtonTol) THEN
                 ! Non-convergence of Newton: predicted error too large
                 Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                 Fac=0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                 EXIT NewtonLoop
               END IF
            END IF
          END IF

          NewtonIncrementOld = MAX(NewtonIncrement,Roundoff) 
          ! Update solution
          CALL DAXPY(N,-ONE,DZ1,1,Z1,1) ! Z1 <- Z1 - DZ1
          CALL DAXPY(N,-ONE,DZ2,1,Z2,1) ! Z2 <- Z2 - DZ2
          CALL DAXPY(N,-ONE,DZ3,1,Z3,1) ! Z3 <- Z3 - DZ3
            
          ! Check error in Newton iterations
          NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
          IF (NewtonDone) THEN
            NewIt = NewtonIter
            EXIT NewtonLoop
          END IF  
          IF (NewtonIter == NewtonMaxit) THEN
            PRINT*, 'Slow or no convergence in Newton Iteration: Max no. of', &
                        'Newton iterations reached'
          END IF

        END DO NewtonLoop

            
        IF (.NOT.NewtonDone) THEN
          !CALL RK_ErrorMsg(-12,T,H,IERR);
          H = Fac*H; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
!???       ISTATUS(Nrej) = ISTATUS(Nrej) + 1
          CYCLE Tloop
        END IF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> SDIRK Stage
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF (SdirkError) THEN
!~~~>  Starting values for Newton iterations
          Z4(1:N) = Z3(1:N)
! ???       
!~~~>   Prepare the loop-independent part of the right-hand side
          CALL FUN(N,T,Y,DZ4)
          ISTATUS(Nfun) = ISTATUS(Nfun) + 1
      
!       G = H*rkBgam(0)*DZ4 + rkTheta(1)*Z1 + rkTheta(2)*Z2 + rkTheta(3)*Z3
          CALL Set2Zero(N, G)
          CALL DAXPY(N,rkBgam(0)*H, DZ4,1,G,1) 
          CALL DAXPY(N,rkTheta(1),Z1,1,G,1)
          CALL DAXPY(N,rkTheta(2),Z2,1,G,1)
          CALL DAXPY(N,rkTheta(3),Z3,1,G,1)

       !~~~>  Initializations for Newton iteration
          NewtonDone = .FALSE.
          Fac = 0.5d0 ! Step reduction factor if too many iterations
            
SDNewtonLoop:DO NewtonIter = 1, NewtonMaxit

!~~~>   Prepare the loop-dependent part of the right-hand side
            CALL WADD(N,Y,Z4,TMP)         ! TMP <- Y + Z4
            CALL FUN(N,T+H,TMP,DZ4)       ! DZ4 <- Fun(Y+Z4)
            ISTATUS(Nfun) = ISTATUS(Nfun) + 1
!            DZ4(1:N) = (G(1:N)-Z4(1:N))*(rkGamma/H) + DZ4(1:N)
            CALL DAXPY (N, -ONE*rkGamma/H, Z4, 1, DZ4, 1)
            CALL DAXPY (N, rkGamma/H, G,1, DZ4,1)

!~~~>   Solve the linear system
            CALL LSS_Solve(Transp,DZ4)
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1 
!~~~>   Check convergence of Newton iterations
            NewtonIncrement = RK_ErrorNorm(N,SCAL,DZ4)
            IF ( NewtonIter == 1 ) THEN
              ThetaSD      = ABS(ThetaMin)
              NewtonRate = 2.0d0 
            ELSE
              ThetaSD = NewtonIncrement/NewtonIncrementOld
              IF (ThetaSD < 0.99d0) THEN
                 NewtonRate = ThetaSD/(ONE-ThetaSD)
                 ! Predict error at the end of Newton process 
                 NewtonPredictedErr = NewtonIncrement &
                            *ThetaSD**(NewtonMaxit-NewtonIter)/(ONE-ThetaSD)
                 IF (NewtonPredictedErr >= NewtonTol) THEN
                      ! Non-convergence of Newton: predicted error too large
                      !print*,'Error too large: ', NewtonPredictedErr
                      Qnewton = MIN(10.0d0,NewtonPredictedErr/NewtonTol)
                      Fac = 0.8d0*Qnewton**(-ONE/(1+NewtonMaxit-NewtonIter))
                      EXIT SDNewtonLoop
                 END IF
              ELSE ! Non-convergence of Newton: ThetaSD too large
                    !print*,'Theta too large: ',ThetaSD
                    EXIT SDNewtonLoop
              END IF
            END IF
            NewtonIncrementOld = NewtonIncrement
            ! Update solution: Z4 <-- Z4 + DZ4
            CALL DAXPY(N,ONE,DZ4,1,Z4,1) 
            
            ! Check error in Newton iterations
            NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
            IF (NewtonDone) EXIT SDNewtonLoop
            
          END DO SDNewtonLoop
            
          IF (.NOT.NewtonDone) THEN
             H = Fac*H; Reject=.TRUE.; SkipJac = .TRUE.;  SkipLU = .FALSE.
             CYCLE Tloop
          END IF
        END IF
!~~~>  End of implified SDIRK Newton iterations

            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Error estimation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IF (SdirkError) THEN
!         DZ4(1:N) =  rkD(1)*Z1 + rkD(2)*Z2 + rkD(3)*Z3 - Z4    
          CALL Set2Zero(N, DZ4)
          IF (rkD(1) /= ZERO) CALL DAXPY(N, rkD(1), Z1, 1, DZ4, 1)
          IF (rkD(2) /= ZERO) CALL DAXPY(N, rkD(2), Z2, 1, DZ4, 1)
          IF (rkD(3) /= ZERO) CALL DAXPY(N, rkD(3), Z3, 1, DZ4, 1)
          CALL DAXPY(N, -ONE, Z4, 1, DZ4, 1)
          Err = RK_ErrorNorm(N,SCAL,DZ4)    
        ELSE
          CALL  RK_ErrorEstimate(N,H,Y,T, &
               Z1,Z2,Z3,SCAL,Err,FirstStep,Reject)
        END IF
!~~~> Computation of new step size Hnew
        Fac  = Err**(-ONE/rkELO)*          &
             MIN(FacSafe,(ONE+2*NewtonMaxit)/(NewtonIter+2*NewtonMaxit))
        Fac  = MIN(FacMax,MAX(FacMin,Fac))
        Hnew = Fac*H

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Accept/reject step 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
accept:IF (Err < ONE) THEN !~~~> STEP IS ACCEPTED
          IF (AdjointType == RK_adj_f90_discrete) THEN ! Save stage solution
            ! CALL DCOPY(N,Z1,1,Zstage(1),1)
            ! CALL DCOPY(N,Z2,1,Zstage(N+1),1)
            ! CALL DCOPY(N,Z3,1,Zstage(2*N+1),1)
            Zstage(1:N)       = Z1(1:N) 
            Zstage(N+1:2*N)   = Z2(1:N)
            Zstage(2*N+1:3*N) = Z3(1:N)
            ! Push old Y - Y at the beginning of the stage
            CALL rk_DPush(T, H, Y, Zstage, NewIt)
          END IF

          IF(GetQuad) CALL RK_UpdateQuad(N,NADJ,H,T,Y,Z1,Z2,Z3,Q)
          FirstStep=.FALSE.
          ISTATUS(Nacc) = ISTATUS(Nacc) + 1
          IF (Gustafsson) THEN
            !~~~> Predictive controller of Gustafsson
            IF (ISTATUS(Nacc) > 1) THEN
               FacGus=FacSafe*(H/Hacc)*(Err**2/ErrOld)**(-0.25d0)
               FacGus=MIN(FacMax,MAX(FacMin,FacGus))
               Fac=MIN(Fac,FacGus)
               Hnew = Fac*H
            END IF
            Hacc=H
            ErrOld=MAX(1.0d-2,Err)
          END IF
          Hold = H
          T = T+H 
         ! Update solution: Y <- Y + sum (d_i Z_i)
          IF (rkD(1) /= ZERO) CALL DAXPY(N,rkD(1),Z1,1,Y,1)
          IF (rkD(2) /= ZERO) CALL DAXPY(N,rkD(2),Z2,1,Y,1)
          IF (rkD(3) /= ZERO) CALL DAXPY(N,rkD(3),Z3,1,Y,1)
          ! Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3
          IF (StartNewton) CALL RK_Interpolate('make',N,H,Hold,Z1,Z2,Z3,CONT)
          CALL RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)
          RSTATUS(Ntexit) = T
          RSTATUS(Nhnew)  = Hnew
          RSTATUS(Nhacc)  = H
          Hnew = Tdirection*MIN( MAX(ABS(Hnew),Hmin) , Hmax )
          IF (Reject) Hnew = Tdirection*MIN(ABS(Hnew),ABS(H))
          Reject = .FALSE.
          IF ((T+Hnew/Qmin-Tend)*Tdirection >=  ZERO) THEN
            H = Tend-T
          ELSE
            Hratio=Hnew/H
            ! Reuse the LU decomposition
            SkipLU = (Theta<=ThetaMin) .AND. (Hratio>=Qmin) .AND. (Hratio<=Qmax)
!???            SkipLU = .false.
            IF (.NOT.SkipLU) H=Hnew
          END IF
          ! If convergence is fast enough, do not update Jacobian
          ! SkipJac = (Theta <= ThetaMin)
          SkipJac = .FALSE.

        ELSE accept !~~~> Step is rejected
          IF (FirstStep .OR. Reject) THEN
             H = FacRej*H
          ELSE
             H = Hnew
          END IF
          Reject   = .TRUE.
          SkipJac  = .TRUE.
          SkipLU   = .FALSE. 
          IF (ISTATUS(Nacc) >= 1) ISTATUS(Nrej) = ISTATUS(Nrej) + 1
        END IF accept
      
      END DO Tloop

!~~~> Save last state: only needed for continuous adjoint
      IF ( (AdjointType == RK_adj_f90_continuous) .OR. &
       (AdjointType == RK_adj_f90_simple_continuous) ) THEN
       CALL Fun(T,Y,Fcn0)
       CALL rk_CPush( T, H, Y, Fcn0, DZ3)
!~~~> Deallocate stage buffer: only needed for discrete adjoint
      ELSEIF (AdjointType == RK_adj_f90_discrete) THEN 
        DEALLOCATE(Zstage, STAT=i)
        IF (i/=0) THEN
          PRINT*,'Deallocation of Zstage failed'
          STOP
        END IF
      END IF   
   
      ! Successful exit
      IERR = 1  

    END SUBROUTINE RK_FwdIntegrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_DadjInt( N,NADJ,Lambda,Tstart,Tend,T,GetQuad,IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
!~~~> Arguments
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)     :: N, NADJ
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tstart, Tend
      DOUBLE PRECISION, INTENT(INOUT) :: T
      INTEGER,  INTENT(INOUT)         :: IERR
      LOGICAL, INTENT(IN)             :: GetQuad
!~~~> Local variables
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WY1,WY2,WY3
      DOUBLE PRECISION, DIMENSION(N) :: SCAL,U1,U2,U3,DU1,DU2,DU3
      DOUBLE PRECISION :: Y(N), Zstage(3*N), X(3*N)
      DOUBLE PRECISION  :: H, NewtonRate, NewtonIncrement, NewtonIncrementOld
      INTEGER :: NewtonIter, ISING, iadj, NewIt
      LOGICAL :: Reject, NewtonDone, NewtonConverge,SkipJac,Transp
      DOUBLE PRECISION :: T1, Theta
      DOUBLE PRECISION, DIMENSION(N) :: TMP, G1, G2, G3

!~~~> If there is a quadrature term
      IF (GetQuad) THEN
        ALLOCATE( WY1(N,NADJ), WY2(N,NADJ), WY3(N,NADJ), STAT=ISING )
        IF(ISING .NE. 0) STOP 'allocation error for WYs'
        WY1(:,:) = 0
        WY2(:,:) = 0
        WY3(:,:) = 0
      END IF
            
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Initial setting
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      T1 = Tend      
!      Tdirection = SIGN(ONE,Tend-Tstart)
      NewtonConverge = .TRUE.
      Reject = .FALSE.
      SkipJac = .FALSE. 
      Transp = .FALSE.
!~~~> Time loop begins below 
TimeLoop:DO WHILE ( stack_ptr > 0 )
        IF (.not.Reject) THEN  
          !~~~>  Recover checkpoints for stage values and vectors
          CALL rk_DPop(T, H, Y, Zstage, NewIt)
          !~~~>  Compute LU decomposition 
          IF (.NOT.SaveLU) THEN
          !~~~> Compute the Jacobian matrix
            CALL LSS_Jac(T,Y,JAC)
            ISTATUS(Njac) = ISTATUS(Njac) + 1
          !~~~> Compute the matrices E1 and E2 and their decompositions
            CALL RK_Decomp(H,ISING)
          END IF
          ISTATUS(Njac) = ISTATUS(Njac) + 3  
          !~~~>   Jacobian values at stage vectors
          CALL WADD(N,Y,Zstage(1),TMP)       ! TMP  <- Y + Z1
          CALL LSS_Jac1(T+rkC(1)*H,TMP,JAC)
          IF(GetQuad) CALL DRDY(NADJ,N,N,T+rkC(1)*H,TMP,WY1)
          CALL WADD(N,Y,Zstage(1+N),TMP)     ! TMP  <- Y + Z2
          CALL LSS_Jac2(T+rkC(2)*H,TMP,JAC)
          IF(GetQuad) CALL DRDY(NADJ,N,N,T+rkC(2)*H,TMP,WY2)
          CALL WADD(N,Y,Zstage(1+2*N),TMP)   ! TMP  <- Y + Z3
          CALL LSS_Jac3(T+rkC(3)*H,TMP,JAC)
          IF(GetQuad) CALL DRDY(NADJ,N,N,T+rkC(3)*H,TMP,WY3)
        END IF ! .not.Reject

111 CONTINUE

        IF ( (AdjointSolve == Solve_adaptive .and. .not.NewtonConverge) &
         .or. (AdjointSolve == Solve_direct) ) THEN

       
          CALL LSS_Decomp_Big(H,rkA,ISING) 
          ISTATUS(Ndec) = ISTATUS(Ndec) + 1
!  Use full big algebra:    
        
        END IF
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~>  Loop for the simplified Newton iterations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Adj:DO iadj = 1, NADJ      
      !~~~>  Starting values for Newton iteration
      ! CALL DCOPY(N,Lambda(1,iadj),1,U1,1)
      ! CALL DCOPY(N,Lambda(1,iadj),1,U2,1)
      ! CALL DCOPY(N,Lambda(1,iadj),1,U3,1)
          CALL Set2Zero(N,U1)
          CALL Set2Zero(N,U2)
          CALL Set2Zero(N,U3)
      
!~~~>  Initializations for Newton iteration
          NewtonDone = .FALSE.
!~~~>    Right Hand Side - part G for all Newton iterations
          CALL SET2ZERO(N,G1) 
          CALL SET2ZERO(N,G2)
          CALL SET2ZERO(N,G3)
          CALL LSS_Mul_Jactr1(TMP,Lambda(1,iadj))
          CALL DAXPY(N,-H*rkB(1),TMP,1,G1,1) ! R1 <- R1 - h*B_1*F1
          IF(GetQuad) CALL DAXPY(N,-H*rkB(1),WY1(1,iadj),1,G1,1)
          CALL LSS_Mul_Jactr2(TMP,Lambda(1,iadj))
          CALL DAXPY(N,-H*rkB(2),TMP,1,G2,1) ! R2 <- R2 - h*B_2*F2
          IF(GetQuad) CALL DAXPY(N,-H*rkB(2),WY2(1,iadj),1,G2,1)
          CALL LSS_Mul_Jactr3(TMP,Lambda(1,iadj))
          CALL DAXPY(N,-H*rkB(3),TMP,1,G3,1) ! R3 <- R3 - h*B_3*F3
          IF(GetQuad) CALL DAXPY(N,-H*rkB(3),WY3(1,iadj),1,G3,1)

!          CALL RK_PrepareRHS_G(N,H,Lambda(1,iadj),G1,G2,G3)
                                       
          IF ( (AdjointSolve == Solve_adaptive .and. NewtonConverge)      &
            .or. (AdjointSolve == Solve_fixed) ) THEN

NewtonLoopAdj:DO  NewtonIter = 1, NewtonMaxit  

            !~~~> Prepare the right-hand side
              CALL RK_PrepareRHS_Adj(N,H,Lambda(1,iadj), &
                        U1,U2,U3,             &
                        G1, G2, G3, DU1,DU2,DU3)

            !~~~> Solve the linear systems
              CALL RK_SolveTR( N,H,DU1,DU2,DU3,ISING )

!~~~> The following code performs an adaptive number of Newton
!     iterations for solving adjoint system
              IF (AdjointSolve == Solve_adaptive) THEN

                CALL RK_ErrorScale(N,ITOL,                      &
                 AbsTol_adj(1:N,iadj),RelTol_adj(1:N,iadj), &
                 Lambda(1:N,iadj),SCAL)

              ! SCAL(1:N) = 1.0d0
                NewtonIncrement = SQRT( ( RK_ErrorNorm(N,SCAL,DU1)**2 +   &
                                RK_ErrorNorm(N,SCAL,DU2)**2 +         &
                                RK_ErrorNorm(N,SCAL,DU3)**2 )/3.0d0 )
            
        
                IF ( NewtonIter == 1 ) THEN
                  Theta      = ABS(ThetaMin)
                  NewtonRate = 2.0d0 
                ELSE
                  Theta = NewtonIncrement/NewtonIncrementOld
                  IF (Theta < 0.99d0) THEN
                      NewtonRate = Theta/(ONE-Theta)
                  ELSE ! Non-convergence of Newton: Theta too large
                      Reject = .TRUE.
                      NewtonDone = .FALSE.
                      EXIT NewtonLoopAdj
                  END IF
                END IF
          
                NewtonIncrementOld = MAX(NewtonIncrement,Roundoff) 

              END IF ! (AdjointSolve == Solve_adaptive)
            
            ! Update solution
              CALL DAXPY(N,-ONE,DU1,1,U1,1) ! U1 <- U1 - DU1
              CALL DAXPY(N,-ONE,DU2,1,U2,1) ! U2 <- U2 - DU2
              CALL DAXPY(N,-ONE,DU3,1,U3,1) ! U3 <- U3 - DU3

              IF (AdjointSolve == Solve_adaptive) THEN
            ! When performing an adaptive number of iterations
            !       check the error in Newton iterations
                NewtonDone = (NewtonRate*NewtonIncrement <= NewtonTol)
                IF ((NewtonDone).and.(NewtonIter>NewIt+1)) EXIT NewtonLoopAdj
              ELSE IF (AdjointSolve == Solve_fixed) THEN
                IF (NewtonIter>MAX(NewIt+1,6)) EXIT NewtonLoopAdj
              END IF   
            
            END DO NewtonLoopAdj
      
            IF ((AdjointSolve==Solve_adaptive).AND.(.NOT.NewtonDone)) THEN
            ! print*,'Newton iterations do not converge, switching to full system.'
              NewtonConverge = .FALSE.
              Reject = .TRUE.
              GOTO 111
            END IF
      
            ! Update adjoint solution: Y_adj <- Y_adj + sum (U_i)
            CALL DAXPY(N,ONE,U1,1,Lambda(1,iadj),1)
            CALL DAXPY(N,ONE,U2,1,Lambda(1,iadj),1)
            CALL DAXPY(N,ONE,U3,1,Lambda(1,iadj),1)

          ELSE ! NewtonConverge = .false.

            X(1:N)       = -G1(1:N)
            X(N+1:2*N)   = -G2(1:N)
            X(2*N+1:3*N) = -G3(1:N)

            Transp = .TRUE.
            CALL LSS_Solve_Big(Transp, X)
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1
            ! CALL WGESL('T',3*N,Jbig,IPbig,X)
            Lambda(1:N,iadj) = Lambda(1:N,iadj)+X(1:N)+X(N+1:2*N)+X(2*N+1:3*N)
            IF ((AdjointSolve==Solve_adaptive).AND.(iadj>=NADJ)) THEN
              NewtonConverge = .TRUE.
              Reject = .FALSE.
            END IF
     
          END IF ! NewtonConverge
 
        END DO Adj
 
        T1 = T1-H
        ISTATUS(Nstp) = ISTATUS(Nstp) + 1
      END DO TimeLoop

      IF(GetQuad) DEALLOCATE(WY1,WY2,WY3,STAT=ISING)
      IF(ISING .NE. 0) STOP 'deallocation error for WYs'
      ! Successful exit
      IERR = 1  

    END SUBROUTINE RK_DadjInt


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_CadjInt (                        &
        NADJ, Y,                                &
        Tstart, Tend, T, IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
!~~~> Arguments
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)     :: NADJ
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tstart, Tend
      DOUBLE PRECISION, INTENT(INOUT) :: T
      INTEGER,  INTENT(INOUT)        :: IERR

    END SUBROUTINE rk_CadjInt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE rk_SimpleCadjInt (                  &
        NADJ, Y,                                &
        Tstart, Tend, T,                        &
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
!~~~> Arguments
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)     :: NADJ
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
      DOUBLE PRECISION, INTENT(IN)    :: Tstart, Tend
      DOUBLE PRECISION, INTENT(INOUT) :: T
      INTEGER,  INTENT(INOUT)        :: IERR

    END SUBROUTINE rk_SimpleCadjInt


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: T, H
      INTEGER, INTENT(IN)  :: Code
      INTEGER, INTENT(INOUT) :: IERR

      IERR = Code
      PRINT * , &
               'Forced exit from RungeKutta due to the following error:'

      SELECT CASE (Code)
      CASE (-1)
        PRINT * , '--> Improper value for maximal no of steps'
      CASE (-2)
        PRINT * , '--> Improper value for maximal no of Newton iterations'
      CASE (-3)
        PRINT * , '--> Hmin/Hmax/Hstart must be positive'
      CASE (-4)
        PRINT * , '--> Improper values for FacMin/FacMax/FacSafe/FacRej'
      CASE (-5)
        PRINT * , '--> Improper value for ThetaMin'
      CASE (-6)
        PRINT * , '--> Newton stopping tolerance too small'
      CASE (-7)
        PRINT * , '--> Improper values for Qmin, Qmax'
      CASE (-8)
        PRINT * , '--> Tolerances are too small'
      CASE (-9)
        PRINT * , '--> No of steps exceeds maximum bound'
      CASE (-10)
        PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
      CASE (-11)
        PRINT * , '--> Matrix is repeatedly singular'
      CASE (-12)
        PRINT * , '--> Non-convergence of Newton iterations'
      CASE (-13)
        PRINT * , '--> Requested RK method not implemented'
      CASE DEFAULT
        PRINT *, 'Unknown Error code: ', Code
      END SELECT

      WRITE(6,FMT="(5X,'T=',E12.5,'  H=',E12.5)") T, H 

    END SUBROUTINE RK_ErrorMsg


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_ErrorScale(N,ITOL,AbsTol,RelTol,Y,SCAL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N, ITOL
      DOUBLE PRECISION, INTENT(IN) :: AbsTol(*), RelTol(*), Y(N)
      DOUBLE PRECISION, INTENT(INOUT) :: SCAL(N)
      INTEGER :: i
   
      IF (ITOL==0) THEN
        DO i=1,N
          SCAL(i)= ONE/(AbsTol(1)+RelTol(1)*ABS(Y(i)))
        END DO
      ELSE
        DO i=1,N
          SCAL(i)=ONE/(AbsTol(i)+RelTol(i)*ABS(Y(i)))
        END DO
      END IF
      
    END SUBROUTINE RK_ErrorScale


!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SUBROUTINE RK_Transform(N,Tr,Z1,Z2,Z3,W1,W2,W3)
!!~~~>                 W <-- Tr x Z
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      IMPLICIT NONE
!      INTEGER :: N, i
!      DOUBLE PRECISION :: Tr(3,3),Z1(N),Z2(N),Z3(N),W1(N),W2(N),W3(N)
!      DOUBLE PRECISION :: x1, x2, x3
!      DO i=1,N
!          x1 = Z1(i); x2 = Z2(i); x3 = Z3(i)
!          W1(i) = Tr(1,1)*x1 + Tr(1,2)*x2 + Tr(1,3)*x3
!          W2(i) = Tr(2,1)*x1 + Tr(2,2)*x2 + Tr(2,3)*x3
!          W3(i) = Tr(3,1)*x1 + Tr(3,2)*x2 + Tr(3,3)*x3
!      END DO
!  END SUBROUTINE RK_Transform
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_Interpolate(action,N,H,Hold,Z1,Z2,Z3,CONT)
!~~~>   Constructs or evaluates a quadratic polynomial
!         that interpolates the Z solution at current step
!         and provides starting values for the next step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: N, i
      DOUBLE PRECISION :: H,Hold,Z1(N),Z2(N),Z3(N),CONT(N,3)
      DOUBLE PRECISION :: r, x1, x2, x3, den
      CHARACTER(LEN=4) :: action
       
      SELECT CASE (action) 
      CASE ('make')
         ! Construct the solution quadratic interpolant Q(c_i) = Z_i, i=1:3
         den = (rkC(3)-rkC(2))*(rkC(2)-rkC(1))*(rkC(1)-rkC(3))
         DO i=1,N
             CONT(i,1)=(-rkC(3)**2*rkC(2)*Z1(i)+Z3(i)*rkC(2)*rkC(1)**2 &
                        +rkC(2)**2*rkC(3)*Z1(i)-rkC(2)**2*rkC(1)*Z3(i) &
                        +rkC(3)**2*rkC(1)*Z2(i)-Z2(i)*rkC(3)*rkC(1)**2)&
                        /den-Z3(i)
             CONT(i,2)= -( rkC(1)**2*(Z3(i)-Z2(i)) + rkC(2)**2*(Z1(i)  &
                          -Z3(i)) +rkC(3)**2*(Z2(i)-Z1(i)) )/den
             CONT(i,3)= ( rkC(1)*(Z3(i)-Z2(i)) + rkC(2)*(Z1(i)-Z3(i))  &
                           +rkC(3)*(Z2(i)-Z1(i)) )/den
         END DO
      CASE ('eval')
          ! Evaluate quadratic polynomial
          r = H/Hold
         x1 = ONE + rkC(1)*r
         x2 = ONE + rkC(2)*r
         x3 = ONE + rkC(3)*r
         DO i=1,N
            Z1(i) = CONT(i,1)+x1*(CONT(i,2)+x1*CONT(i,3))
            Z2(i) = CONT(i,1)+x2*(CONT(i,2)+x2*CONT(i,3))
            Z3(i) = CONT(i,1)+x3*(CONT(i,2)+x3*CONT(i,3))
         END DO
       END SELECT   
    END SUBROUTINE RK_Interpolate


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_UpdateQuad(N,NADJ,H,T,Y,Z1,Z2,Z3,Q)
!~~~> Update the quadrature term
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: N, NADJ
      DOUBLE PRECISION :: T,H
      DOUBLE PRECISION :: Q(NADJ),Y(N),Z1(N),Z2(N),Z3(N)
!~~~> local variables
      DOUBLE PRECISION :: F(NADJ),TMP(N)
      TMP(:) = 0.d0
      F(:) = 0.d0
      CALL WADD(N,Y,Z1,TMP)
      CALL QFUN(N,NADJ,T+rkC(1)*H,TMP,F)
      CALL DAXPY(NADJ,H*rkB(1),F,1,Q,1)
      CALL WADD(N,Y,Z2,TMP)
      CALL QFUN(N,NADJ,T+rkC(2)*H,TMP,F)
      CALL DAXPY(NADJ,H*rkB(2),F,1,Q,1)
      CALL WADD(N,Y,Z3,TMP)
      CALL QFUN(N,NADJ,T+rkC(3)*H,TMP,F)
      CALL DAXPY(NADJ,H*rkB(3),F,1,Q,1)

    END SUBROUTINE RK_UpdateQuad


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_PrepareRHS(N,T,H,Y,Z1,Z2,Z3,R1,R2,R3)
!~~~> Prepare the right-hand side for Newton iterations
!     R = Z - hA x F
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER :: N
      DOUBLE PRECISION :: T,H
      DOUBLE PRECISION, DIMENSION(N) :: Y,Z1,Z2,Z3,F,R1,R2,R3,TMP

      CALL DCOPY(N,Z1,1,R1,1) ! R1 <- Z1
      CALL DCOPY(N,Z2,1,R2,1) ! R2 <- Z2
      CALL DCOPY(N,Z3,1,R3,1) ! R3 <- Z3

      CALL WADD(N,Y,Z1,TMP)              ! TMP <- Y + Z1
      CALL FUN(N,T+rkC(1)*H,TMP,F)    ! F1 <- Fun(Y+Z1)
      CALL DAXPY(N,-H*rkA(1,1),F,1,R1,1) ! R1 <- R1 - h*A_11*F1
      CALL DAXPY(N,-H*rkA(2,1),F,1,R2,1) ! R2 <- R2 - h*A_21*F1
      CALL DAXPY(N,-H*rkA(3,1),F,1,R3,1) ! R3 <- R3 - h*A_31*F1

      CALL WADD(N,Y,Z2,TMP)              ! TMP <- Y + Z2
      CALL FUN(N,T+rkC(2)*H,TMP,F)    ! F2 <- Fun(Y+Z2)
      CALL DAXPY(N,-H*rkA(1,2),F,1,R1,1) ! R1 <- R1 - h*A_12*F2
      CALL DAXPY(N,-H*rkA(2,2),F,1,R2,1) ! R2 <- R2 - h*A_22*F2
      CALL DAXPY(N,-H*rkA(3,2),F,1,R3,1) ! R3 <- R3 - h*A_32*F2

      CALL WADD(N,Y,Z3,TMP)              ! TMP <- Y + Z3
      CALL FUN(N,T+rkC(3)*H,TMP,F)    ! F3 <- Fun(Y+Z3)
      CALL DAXPY(N,-H*rkA(1,3),F,1,R1,1) ! R1 <- R1 - h*A_13*F3
      CALL DAXPY(N,-H*rkA(2,3),F,1,R2,1) ! R2 <- R2 - h*A_23*F3
      CALL DAXPY(N,-H*rkA(3,3),F,1,R3,1) ! R3 <- R3 - h*A_33*F3
            
   END SUBROUTINE RK_PrepareRHS
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   SUBROUTINE RK_PrepareRHS_Adj(N,H,Lambda, &
                      U1,U2,U3,G1,G2,G3,R1,R2,R3)
!~~~> Prepare the right-hand side for Newton iterations
!     R = Z_adj - hA x Jac*Z_adj - h J^t b lambda
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(IN) :: H
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: U1,U2,U3,Lambda,G1,G2,G3
      DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: R1,R2,R3
      DOUBLE PRECISION, DIMENSION(N) ::  F,TMP


      CALL DCOPY(N,G1,1,R1,1) ! R1 <- G1
      CALL DCOPY(N,G2,1,R2,1) ! R2 <- G2
      CALL DCOPY(N,G3,1,R3,1) ! R3 <- G3

      CALL SET2ZERO(N,F)
      CALL DAXPY(N,-H*rkA(1,1),U1,1,F,1) ! F1 <- -h*A_11*U1
      CALL DAXPY(N,-H*rkA(2,1),U2,1,F,1) ! F1 <- F1 - h*A_21*U2
      CALL DAXPY(N,-H*rkA(3,1),U3,1,F,1) ! F1 <- F1 - h*A_31*U3
      CALL LSS_Mul_Jactr1(TMP,F)
      CALL DAXPY(N,ONE,U1,1,TMP,1) ! R1 <- U1 -Jac(Y+Z1)^t*h*sum(A_j1*U_j)
      CALL DAXPY(N,ONE,TMP,1,R1,1) ! R1 <- U1 -Jac(Y+Z1)^t*h*sum(A_j1*U_j)

      CALL SET2ZERO(N,F)
      CALL DAXPY(N,-H*rkA(1,2),U1,1,F,1) ! F2 <- -h*A_11*U1
      CALL DAXPY(N,-H*rkA(2,2),U2,1,F,1) ! F2 <- F2 - h*A_21*U2
      CALL DAXPY(N,-H*rkA(3,2),U3,1,F,1) ! F2 <- F2 - h*A_31*U3
      CALL LSS_Mul_Jactr2(TMP,F)
      CALL DAXPY(N,ONE,U2,1,TMP,1) ! R2 <- U2 -Jac(Y+Z2)^t*h*sum(A_j2*U_j)
      CALL DAXPY(N,ONE,TMP,1,R2,1) ! R2 <- U2 -Jac(Y+Z2)^t*h*sum(A_j2*U_j)

      CALL SET2ZERO(N,F)
      CALL DAXPY(N,-H*rkA(1,3),U1,1,F,1) ! F3 <- -h*A_11*U1
      CALL DAXPY(N,-H*rkA(2,3),U2,1,F,1) ! F3 <- F3 - h*A_21*U2
      CALL DAXPY(N,-H*rkA(3,3),U3,1,F,1) ! F3 <- F3 - h*A_31*U3
      CALL LSS_Mul_Jactr3(TMP,F)
      CALL DAXPY(N,ONE,U3,1,TMP,1) ! R3 <- U3 -Jac(Y+Z3)^t*h*sum(A_j3*U_j)
      CALL DAXPY(N,ONE,TMP,1,R3,1) ! R3 <- U3 -Jac(Y+Z3)^t*h*sum(A_j3*U_j)

    END SUBROUTINE RK_PrepareRHS_Adj

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_Decomp(H,ISING)
   !~~~> Compute the matrices E1 and E2 and their decompositions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      
      INTEGER :: ISING
      DOUBLE PRECISION    :: H, Alpha, Beta, Gamma
!      DOUBLE PRECISION    :: FJAC(N,N),E1(N,N)
!      COMPLEX(kind=dp) :: E2(N,N)
      
      Gamma = rkGamma/H
      Alpha = rkAlpha/H
      Beta  = rkBeta /H

      CALL LSS_Decomp(Gamma, ISING)
      
      IF (ISING /= 0) THEN
         ISTATUS(Ndec) = ISTATUS(Ndec) + 1
         RETURN
      END IF
     
      CALL LSS_Decomp_Cmp(Alpha,Beta,ISING)
      ISTATUS(Ndec) = ISTATUS(Ndec) + 1
      
    END SUBROUTINE RK_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_Solve(N,H,R1,R2,R3,ISING)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      INTEGER :: N,ISING
      DOUBLE PRECISION    :: R1(N),R2(N),R3(N)
      DOUBLE PRECISION    :: H, x1, x2, x3
      INTEGER :: i
      LOGICAL :: Transp
!      
     ! Z <- h^{-1) T^{-1) A^{-1) x Z
      DO i=1,N
          x1 = R1(i)/H; x2 = R2(i)/H; x3 = R3(i)/H
          R1(i) = rkTinvAinv(1,1)*x1 + rkTinvAinv(1,2)*x2 + rkTinvAinv(1,3)*x3
          R2(i) = rkTinvAinv(2,1)*x1 + rkTinvAinv(2,2)*x2 + rkTinvAinv(2,3)*x3
          R3(i) = rkTinvAinv(3,1)*x1 + rkTinvAinv(3,2)*x2 + rkTinvAinv(3,3)*x3
      END DO
      Transp = .FALSE.
      CALL LSS_Solve(Transp,R1)
      CALL LSS_Solve_CMP(Transp, R2, R3)

      ! Z <- T x Z
      DO i=1,N
          x1 = R1(i); x2 = R2(i); x3 = R3(i)
          R1(i) = rkT(1,1)*x1 + rkT(1,2)*x2 + rkT(1,3)*x3
          R2(i) = rkT(2,1)*x1 + rkT(2,2)*x2 + rkT(2,3)*x3
          R3(i) = rkT(3,1)*x1 + rkT(3,2)*x2 + rkT(3,3)*x3
      END DO

      ISTATUS(Nsol) = ISTATUS(Nsol) + 2

    END SUBROUTINE RK_Solve


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_SolveTR(N,H,R1,R2,R3,ISING)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(INOUT) :: ISING
      
      DOUBLE PRECISION, INTENT(INOUT) :: R1(N),R2(N),R3(N)
      DOUBLE PRECISION    :: H, x1, x2, x3
      INTEGER :: i
      LOGICAL :: Transp
!      
     ! RHS <- h^{-1) (A^{-1) T^{-1))^t x RHS
      DO i=1,N
          x1 = R1(i)/H; x2 = R2(i)/H; x3 = R3(i)/H
          R1(i) = rkAinvT(1,1)*x1 + rkAinvT(2,1)*x2 + rkAinvT(3,1)*x3
          R2(i) = rkAinvT(1,2)*x1 + rkAinvT(2,2)*x2 + rkAinvT(3,2)*x3
          R3(i) = rkAinvT(1,3)*x1 + rkAinvT(2,3)*x2 + rkAinvT(3,3)*x3
      END DO
      Transp = .TRUE.
      CALL LSS_Solve(Transp,R1)
      
      R3(:) = -R3(:)
      CALL LSS_Solve_CMP(Transp,R2,R3)      
      R3(:) = -R3(:)      

      ! X <- (T^{-1})^t x X
      DO i=1,N
          x1 = R1(i); x2 = R2(i); x3 = R3(i)
          R1(i) = rkTinv(1,1)*x1 + rkTinv(2,1)*x2 + rkTinv(3,1)*x3
          R2(i) = rkTinv(1,2)*x1 + rkTinv(2,2)*x2 + rkTinv(3,2)*x3
          R3(i) = rkTinv(1,3)*x1 + rkTinv(2,3)*x2 + rkTinv(3,3)*x3
      END DO

      ISTATUS(Nsol) = ISTATUS(Nsol) + 2

    END SUBROUTINE RK_SolveTR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE RK_ErrorEstimate(N,H,Y,T, &
               Z1,Z2,Z3,SCAL,Err,     &
               FirstStep,Reject)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      
      INTEGER :: N
      DOUBLE PRECISION :: SCAL(N),Z1(N),Z2(N),Z3(N),F1(N),F2(N), &
                        F0(N),Y(N),TMP(N),T,H
      INTEGER :: i
      LOGICAL FirstStep,Reject,Transp
      DOUBLE PRECISION :: HEE1,HEE2,HEE3,Err

      HEE1  = rkE(1)/H
      HEE2  = rkE(2)/H
      HEE3  = rkE(3)/H

      CALL FUN(N,T,Y,F0)
      ISTATUS(Nfun) = ISTATUS(Nfun) + 1

      DO  i=1,N
         F2(i)  = HEE1*Z1(i)+HEE2*Z2(i)+HEE3*Z3(i)
         TMP(i) = rkE(0)*F0(i) + F2(i)
      END DO
      Transp = .FALSE.
      CALL LSS_Solve(Transp,TMP)
      ISTATUS(Nsol) = ISTATUS(Nsol) + 1
      IF ((rkMethod==R1A).OR.(rkMethod==GAU).OR.(rkMethod==L3A)) then
       CALL LSS_Solve(Transp,TMP)
      END IF
      IF (rkMethod==GAU) CALL LSS_Solve(Transp,TMP)

      Err = RK_ErrorNorm(N,SCAL,TMP)
!
      IF (Err < ONE) RETURN
fNrej:IF (FirstStep.OR.Reject) THEN
          DO i=1,N
             TMP(i)=Y(i)+TMP(i)
          END DO
          CALL FUN(N,T,TMP,F1)
          ISTATUS(Nfun) = ISTATUS(Nfun) + 1
          DO i=1,N
             TMP(i)=F1(i)+F2(i)
          END DO

          CALL LSS_Solve(Transp,TMP)     
          ISTATUS(Nsol) = ISTATUS(Nsol) + 1
          Err = RK_ErrorNorm(N,SCAL,TMP)
       END IF fNrej
 
    END SUBROUTINE RK_ErrorEstimate

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOUBLE PRECISION FUNCTION RK_ErrorNorm(N,SCAL,DY)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      
      INTEGER :: N
      DOUBLE PRECISION :: SCAL(N),DY(N)
      INTEGER :: i

      RK_ErrorNorm = ZERO
      DO i=1,N
          RK_ErrorNorm = RK_ErrorNorm + (DY(i)*SCAL(i))**2
      END DO
      RK_ErrorNorm = MAX( SQRT(RK_ErrorNorm/N), 1.0d-10 )
 
    END FUNCTION RK_ErrorNorm
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Radau2A_Coefficients
!    The coefficients of the 3-stage Radau-2A method
!    (given to ~30 accurate digits)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
! The coefficients of the Radau2A method
      DOUBLE PRECISION :: b0

!      b0 = 1.0d0
      IF (SdirkError) THEN
        b0 = 0.2d-1
      ELSE
        b0 = 0.5d-1
      END IF

! The coefficients of the Radau2A method
      rkMethod = R2A

      rkA(1,1) =  1.968154772236604258683861429918299d-1
      rkA(1,2) = -6.55354258501983881085227825696087d-2
      rkA(1,3) =  2.377097434822015242040823210718965d-2
      rkA(2,1) =  3.944243147390872769974116714584975d-1
      rkA(2,2) =  2.920734116652284630205027458970589d-1
      rkA(2,3) = -4.154875212599793019818600988496743d-2
      rkA(3,1) =  3.764030627004672750500754423692808d-1
      rkA(3,2) =  5.124858261884216138388134465196080d-1
      rkA(3,3) =  1.111111111111111111111111111111111d-1

      rkB(1) = 3.764030627004672750500754423692808d-1
      rkB(2) = 5.124858261884216138388134465196080d-1
      rkB(3) = 1.111111111111111111111111111111111d-1

      rkC(1) = 1.550510257216821901802715925294109d-1
      rkC(2) = 6.449489742783178098197284074705891d-1
      rkC(3) = 1.0d0
      
      ! New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
      rkD(1) = 0.0d0
      rkD(2) = 0.0d0
      rkD(3) = 1.0d0

      ! Classical error estimator: 
      ! H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) = 1.0d0*b0
      rkE(1) = -10.04880939982741556246032950764708d0*b0
      rkE(2) = 1.382142733160748895793662840980412d0*b0
      rkE(3) = -.3333333333333333333333333333333333d0*b0

      ! Sdirk error estimator
      rkBgam(0) = b0
      rkBgam(1) = .3764030627004672750500754423692807d0-1.558078204724922382431975370686279d0*b0
      rkBgam(2) = .8914115380582557157653087040196118d0*b0+.5124858261884216138388134465196077d0
      rkBgam(3) = -.1637777184845662566367174924883037d0-.3333333333333333333333333333333333d0*b0
      rkBgam(4) = .2748888295956773677478286035994148d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -1.520677486405081647234271944611547d0-10.04880939982741556246032950764708d0*b0
      rkTheta(2) = 2.070455145596436382729929151810376d0+1.382142733160748895793662840980413d0*b0
      rkTheta(3) = -.3333333333333333333333333333333333d0*b0-.3744441479783868387391430179970741d0

      ! Local order of error estimator 
      IF (b0==0.0d0) THEN
        rkELO  = 6.0d0
      ELSE      
        rkELO  = 4.0d0
      END IF    

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 3.637834252744495732208418513577775d0
      rkAlpha = 2.681082873627752133895790743211112d0
      rkBeta  = 3.050430199247410569426377624787569d0

      rkT(1,1) =  9.443876248897524148749007950641664d-2
      rkT(1,2) = -1.412552950209542084279903838077973d-1
      rkT(1,3) = -3.00291941051474244918611170890539d-2
      rkT(2,1) =  2.502131229653333113765090675125018d-1
      rkT(2,2) =  2.041293522937999319959908102983381d-1
      rkT(2,3) =  3.829421127572619377954382335998733d-1
      rkT(3,1) =  1.0d0
      rkT(3,2) =  1.0d0
      rkT(3,3) =  0.0d0

      rkTinv(1,1) =  4.178718591551904727346462658512057d0
      rkTinv(1,2) =  3.27682820761062387082533272429617d-1
      rkTinv(1,3) =  5.233764454994495480399309159089876d-1
      rkTinv(2,1) = -4.178718591551904727346462658512057d0
      rkTinv(2,2) = -3.27682820761062387082533272429617d-1
      rkTinv(2,3) =  4.766235545005504519600690840910124d-1
      rkTinv(3,1) = -5.02872634945786875951247343139544d-1
      rkTinv(3,2) =  2.571926949855605429186785353601676d0
      rkTinv(3,3) = -5.960392048282249249688219110993024d-1

      rkTinvAinv(1,1) =  1.520148562492775501049204957366528d+1
      rkTinvAinv(1,2) =  1.192055789400527921212348994770778d0
      rkTinvAinv(1,3) =  1.903956760517560343018332287285119d0
      rkTinvAinv(2,1) = -9.669512977505946748632625374449567d0
      rkTinvAinv(2,2) = -8.724028436822336183071773193986487d0
      rkTinvAinv(2,3) =  3.096043239482439656981667712714881d0
      rkTinvAinv(3,1) = -1.409513259499574544876303981551774d+1
      rkTinvAinv(3,2) =  5.895975725255405108079130152868952d0
      rkTinvAinv(3,3) = -1.441236197545344702389881889085515d-1

      rkAinvT(1,1) = .3435525649691961614912493915818282d0
      rkAinvT(1,2) = -.4703191128473198422370558694426832d0
      rkAinvT(1,3) = .3503786597113668965366406634269080d0
      rkAinvT(2,1) = .9102338692094599309122768354288852d0
      rkAinvT(2,2) = 1.715425895757991796035292755937326d0
      rkAinvT(2,3) = .4040171993145015239277111187301784d0
      rkAinvT(3,1) = 3.637834252744495732208418513577775d0
      rkAinvT(3,2) = 2.681082873627752133895790743211112d0
      rkAinvT(3,3) = -3.050430199247410569426377624787569d0

    END SUBROUTINE Radau2A_Coefficients

    

    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Lobatto3C_Coefficients
!    The coefficients of the 3-stage Lobatto-3C method
!    (given to ~30 accurate digits)
!    The parameter b0 can be chosen to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
      DOUBLE PRECISION :: b0

      rkMethod = L3C

!      b0 = 1.0d0
      IF (SdirkError) THEN
        b0 = 0.2d0
      ELSE
        b0 = 0.5d0
      END IF
! The coefficients of the Lobatto3C method

      rkA(1,1) =  .1666666666666666666666666666666667d0
      rkA(1,2) = -.3333333333333333333333333333333333d0
      rkA(1,3) =  .1666666666666666666666666666666667d0
      rkA(2,1) =  .1666666666666666666666666666666667d0
      rkA(2,2) =  .4166666666666666666666666666666667d0
      rkA(2,3) = -.8333333333333333333333333333333333d-1
      rkA(3,1) =  .1666666666666666666666666666666667d0
      rkA(3,2) =  .6666666666666666666666666666666667d0
      rkA(3,3) =  .1666666666666666666666666666666667d0

      rkB(1) = .1666666666666666666666666666666667d0
      rkB(2) = .6666666666666666666666666666666667d0
      rkB(3) = .1666666666666666666666666666666667d0

      rkC(1) = 0.0d0
      rkC(2) = 0.5d0
      rkC(3) = 1.0d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) = .16666666666666666666666666666666667d0-b0
      rkBhat(2) = .66666666666666666666666666666666667d0
      rkBhat(3) = .16666666666666666666666666666666667d0
      
      ! New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
      rkD(1) = 0.0d0
      rkD(2) = 0.0d0
      rkD(3) = 1.0d0

      ! Classical error estimator: 
      !   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) =   .3808338772072650364017425226487022*b0
      rkE(1) = -1.142501631621795109205227567946107*b0
      rkE(2) = -1.523335508829060145606970090594809*b0
      rkE(3) =   .3808338772072650364017425226487022*b0

      ! Sdirk error estimator
      rkBgam(0) = b0
      rkBgam(1) = .1666666666666666666666666666666667d0-1.d0*b0
      rkBgam(2) = .6666666666666666666666666666666667d0
      rkBgam(3) = -.2141672105405983697350758559820354d0
      rkBgam(4) = .3808338772072650364017425226487021d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -3.d0*b0-.3808338772072650364017425226487021d0
      rkTheta(2) = -4.d0*b0+1.523335508829060145606970090594808d0
      rkTheta(3) = -.142501631621795109205227567946106d0+b0

      ! Local order of error estimator 
      IF (b0==0.0d0) THEN
        rkELO  = 5.0d0
      ELSE      
        rkELO  = 4.0d0
      END IF    

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 2.625816818958466716011888933765284d0
      rkAlpha = 1.687091590520766641994055533117359d0
      rkBeta  = 2.508731754924880510838743672432351d0

      rkT(1,1) = 1.d0
      rkT(1,2) = 1.d0
      rkT(1,3) = 0.d0
      rkT(2,1) = .4554100411010284672111720348287483d0
      rkT(2,2) = -.6027050205505142336055860174143743d0
      rkT(2,3) = -.4309321229203225731070721341350346d0
      rkT(3,1) = 2.195823345445647152832799205549709d0
      rkT(3,2) = -1.097911672722823576416399602774855d0
      rkT(3,3) = .7850032632435902184104551358922130d0

      rkTinv(1,1) = .4205559181381766909344950150991349d0
      rkTinv(1,2) = .3488903392193734304046467270632057d0
      rkTinv(1,3) = .1915253879645878102698098373933487d0
      rkTinv(2,1) = .5794440818618233090655049849008650d0
      rkTinv(2,2) = -.3488903392193734304046467270632057d0
      rkTinv(2,3) = -.1915253879645878102698098373933487d0
      rkTinv(3,1) = -.3659705575742745254721332009249516d0
      rkTinv(3,2) = -1.463882230297098101888532803699806d0
      rkTinv(3,3) = .4702733607340189781407813565524989d0

      rkTinvAinv(1,1) = 1.104302803159744452668648155627548d0
      rkTinvAinv(1,2) = .916122120694355522658740710823143d0
      rkTinvAinv(1,3) = .5029105849749601702795812241441172d0
      rkTinvAinv(2,1) = 1.895697196840255547331351844372453d0
      rkTinvAinv(2,2) = 3.083877879305644477341259289176857d0
      rkTinvAinv(2,3) = -1.502910584974960170279581224144117d0
      rkTinvAinv(3,1) = .8362439183082935036129145574774502d0
      rkTinvAinv(3,2) = -3.344975673233174014451658229909802d0
      rkTinvAinv(3,3) = .312908409479233358005944466882642d0

      rkAinvT(1,1) = 2.625816818958466716011888933765282d0
      rkAinvT(1,2) = 1.687091590520766641994055533117358d0
      rkAinvT(1,3) = -2.508731754924880510838743672432351d0
      rkAinvT(2,1) = 1.195823345445647152832799205549710d0
      rkAinvT(2,2) = -2.097911672722823576416399602774855d0
      rkAinvT(2,3) = .7850032632435902184104551358922130d0
      rkAinvT(3,1) = 5.765829871932827589653709477334136d0
      rkAinvT(3,2) = .1170850640335862051731452613329320d0
      rkAinvT(3,3) = 4.078738281412060947659653944216779d0

    END SUBROUTINE Lobatto3C_Coefficients

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Gauss_Coefficients
!    The coefficients of the 3-stage Gauss method
!    (given to ~30 accurate digits)
!    The parameter b3 can be chosen by the user
!    to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
      DOUBLE PRECISION :: b0
! The coefficients of the Gauss method


      rkMethod = GAU
      
!      b0 = 4.0d0
      b0 = 0.1d0
      
! The coefficients of the Gauss method

      rkA(1,1) =  .1388888888888888888888888888888889d0
      rkA(1,2) = -.359766675249389034563954710966045d-1
      rkA(1,3) =  .97894440153083260495800422294756d-2
      rkA(2,1) =  .3002631949808645924380249472131556d0
      rkA(2,2) =  .2222222222222222222222222222222222d0
      rkA(2,3) = -.224854172030868146602471694353778d-1
      rkA(3,1) =  .2679883337624694517281977355483022d0
      rkA(3,2) =  .4804211119693833479008399155410489d0
      rkA(3,3) =  .1388888888888888888888888888888889d0

      rkB(1) = .2777777777777777777777777777777778d0
      rkB(2) = .4444444444444444444444444444444444d0
      rkB(3) = .2777777777777777777777777777777778d0

      rkC(1) = .1127016653792583114820734600217600d0
      rkC(2) = .5000000000000000000000000000000000d0
      rkC(3) = .8872983346207416885179265399782400d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) =-1.4788305577012361475298775666303999d0*b0 &
                  +.27777777777777777777777777777777778d0
      rkBhat(2) =  .44444444444444444444444444444444444d0 &
                  +.66666666666666666666666666666666667d0*b0
      rkBhat(3) = -.18783610896543051913678910003626672d0*b0 &
                  +.27777777777777777777777777777777778d0

      ! New solution: h Sum_j b_j f(Z_j) = sum d_j Z_j
      rkD(1) = .1666666666666666666666666666666667d1
      rkD(2) = -.1333333333333333333333333333333333d1
      rkD(3) = .1666666666666666666666666666666667d1

      ! Classical error estimator: 
      !   H* Sum (B_j-Bhat_j)*f(Z_j) = H*E(0)*f(0) + Sum E_j*Z_j
      rkE(0) = .2153144231161121782447335303806954d0*b0
      rkE(1) = -2.825278112319014084275808340593191d0*b0
      rkE(2) = .2870858974881495709929780405075939d0*b0
      rkE(3) = -.4558086256248162565397206448274867d-1*b0

      ! Sdirk error estimator
      rkBgam(0) = 0.d0
      rkBgam(1) = .2373339543355109188382583162660537d0
      rkBgam(2) = .5879873931885192299409334646982414d0
      rkBgam(3) = -.4063577064014232702392531134499046d-1
      rkBgam(4) = .2153144231161121782447335303806955d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -2.594040933093095272574031876464493d0
      rkTheta(2) = 1.824611539036311947589425112250199d0
      rkTheta(3) = .1856563166634371860478043996459493d0

      ! ELO = local order of classical error estimator 
      rkELO = 4.0d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 4.644370709252171185822941421408064d0
      rkAlpha = 3.677814645373914407088529289295970d0
      rkBeta  = 3.508761919567443321903661209182446d0

      rkT(1,1) =  .7215185205520017032081769924397664d-1
      rkT(1,2) = -.8224123057363067064866206597516454d-1
      rkT(1,3) = -.6012073861930850173085948921439054d-1
      rkT(2,1) =  .1188325787412778070708888193730294d0
      rkT(2,2) =  .5306509074206139504614411373957448d-1
      rkT(2,3) =  .3162050511322915732224862926182701d0
      rkT(3,1) = 1.d0
      rkT(3,2) = 1.d0
      rkT(3,3) = 0.d0

      rkTinv(1,1) =  5.991698084937800775649580743981285d0
      rkTinv(1,2) =  1.139214295155735444567002236934009d0
      rkTinv(1,3) =   .4323121137838583855696375901180497d0
      rkTinv(2,1) = -5.991698084937800775649580743981285d0
      rkTinv(2,2) = -1.139214295155735444567002236934009d0
      rkTinv(2,3) =   .5676878862161416144303624098819503d0
      rkTinv(3,1) = -1.246213273586231410815571640493082d0
      rkTinv(3,2) =  2.925559646192313662599230367054972d0
      rkTinv(3,3) =  -.2577352012734324923468722836888244d0

      rkTinvAinv(1,1) =  27.82766708436744962047620566703329d0
      rkTinvAinv(1,2) =   5.290933503982655311815946575100597d0
      rkTinvAinv(1,3) =   2.007817718512643701322151051660114d0
      rkTinvAinv(2,1) = -17.66368928942422710690385180065675d0
      rkTinvAinv(2,2) = -14.45491129892587782538830044147713d0
      rkTinvAinv(2,3) =   2.992182281487356298677848948339886d0
      rkTinvAinv(3,1) = -25.60678350282974256072419392007303d0
      rkTinvAinv(3,2) =   6.762434375611708328910623303779923d0
      rkTinvAinv(3,3) =   1.043979339483109825041215970036771d0
      
      rkAinvT(1,1) = .3350999483034677402618981153470483d0
      rkAinvT(1,2) = -.5134173605009692329246186488441294d0
      rkAinvT(1,3) = .6745196507033116204327635673208923d-1
      rkAinvT(2,1) = .5519025480108928886873752035738885d0
      rkAinvT(2,2) = 1.304651810077110066076640761092008d0
      rkAinvT(2,3) = .9767507983414134987545585703726984d0
      rkAinvT(3,1) = 4.644370709252171185822941421408064d0
      rkAinvT(3,2) = 3.677814645373914407088529289295970d0
      rkAinvT(3,3) = -3.508761919567443321903661209182446d0
      
    END SUBROUTINE Gauss_Coefficients


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE Radau1A_Coefficients
!    The coefficients of the 3-stage Gauss method
!    (given to ~30 accurate digits)
!    The parameter b3 can be chosen by the user
!    to tune the error estimator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      IMPLICIT NONE
!      DOUBLE PRECISION :: b0 = 0.3d0
      DOUBLE PRECISION :: b0 = 0.1d0

! The coefficients of the Radau1A method

      rkMethod = R1A

      rkA(1,1) =  .1111111111111111111111111111111111d0
      rkA(1,2) = -.1916383190435098943442935597058829d0
      rkA(1,3) =  .8052720793239878323318244859477174d-1
      rkA(2,1) =  .1111111111111111111111111111111111d0
      rkA(2,2) =  .2920734116652284630205027458970589d0
      rkA(2,3) = -.481334970546573839513422644787591d-1
      rkA(3,1) =  .1111111111111111111111111111111111d0
      rkA(3,2) =  .5370223859435462728402311533676479d0
      rkA(3,3) =  .1968154772236604258683861429918299d0

      rkB(1) = .1111111111111111111111111111111111d0
      rkB(2) = .5124858261884216138388134465196080d0
      rkB(3) = .3764030627004672750500754423692808d0

      rkC(1) = 0.d0
      rkC(2) = .3550510257216821901802715925294109d0
      rkC(3) = .8449489742783178098197284074705891d0

      ! Classical error estimator, embedded solution: 
      rkBhat(0) = b0
      rkBhat(1) = .11111111111111111111111111111111111d0-b0
      rkBhat(2) = .51248582618842161383881344651960810d0
      rkBhat(3) = .37640306270046727505007544236928079d0

      ! New solution: H* Sum B_j*f(Z_j) = Sum D_j*Z_j
      rkD(1) = .3333333333333333333333333333333333d0
      rkD(2) = -.8914115380582557157653087040196127d0
      rkD(3) = .1558078204724922382431975370686279d1

      ! Classical error estimator: 
      ! H* Sum (b_j-bhat_j) f(Z_j) = H*E(0)*F(0) + Sum E_j Z_j
      rkE(0) =   .2748888295956773677478286035994148d0*b0
      rkE(1) = -1.374444147978386838739143017997074d0*b0
      rkE(2) = -1.335337922441686804550326197041126d0*b0
      rkE(3) =   .235782604058977333559011782643466d0*b0

      ! Sdirk error estimator
      rkBgam(0) = 0.0d0
      rkBgam(1) = .1948150124588532186183490991130616d-1
      rkBgam(2) = .7575249005733381398986810981093584d0
      rkBgam(3) = -.518952314149008295083446116200793d-1
      rkBgam(4) = .2748888295956773677478286035994148d0

      ! H* Sum Bgam_j*f(Z_j) = H*Bgam(0)*f(0) + Sum Theta_j*Z_j
      rkTheta(1) = -1.224370034375505083904362087063351d0
      rkTheta(2) = .9340045331532641409047527962010133d0
      rkTheta(3) = .4656990124352088397561234800640929d0

      ! ELO = local order of classical error estimator 
      rkELO = 4.0d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~> Diagonalize the RK matrix:               
      ! rkTinv * inv(rkA) * rkT =          
      !           |  rkGamma      0           0     |
      !           |      0      rkAlpha   -rkBeta   |
      !           |      0      rkBeta     rkAlpha  |
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      rkGamma = 3.637834252744495732208418513577775d0
      rkAlpha = 2.681082873627752133895790743211112d0
      rkBeta  = 3.050430199247410569426377624787569d0

      rkT(1,1) =  .424293819848497965354371036408369d0
      rkT(1,2) = -.3235571519651980681202894497035503d0
      rkT(1,3) = -.522137786846287839586599927945048d0
      rkT(2,1) =  .57594609499806128896291585429339d-1
      rkT(2,2) =  .3148663231849760131614374283783d-2
      rkT(2,3) =  .452429247674359778577728510381731d0
      rkT(3,1) = 1.d0
      rkT(3,2) = 1.d0
      rkT(3,3) = 0.d0

      rkTinv(1,1) = 1.233523612685027760114769983066164d0
      rkTinv(1,2) = 1.423580134265707095505388133369554d0
      rkTinv(1,3) = .3946330125758354736049045150429624d0
      rkTinv(2,1) = -1.233523612685027760114769983066164d0
      rkTinv(2,2) = -1.423580134265707095505388133369554d0
      rkTinv(2,3) = .6053669874241645263950954849570376d0
      rkTinv(3,1) = -.1484438963257383124456490049673414d0
      rkTinv(3,2) = 2.038974794939896109682070471785315d0
      rkTinv(3,3) = -.544501292892686735299355831692542d-1

      rkTinvAinv(1,1) =  4.487354449794728738538663081025420d0
      rkTinvAinv(1,2) =  5.178748573958397475446442544234494d0
      rkTinvAinv(1,3) =  1.435609490412123627047824222335563d0
      rkTinvAinv(2,1) = -2.854361287939276673073807031221493d0
      rkTinvAinv(2,2) = -1.003648660720543859000994063139137d+1
      rkTinvAinv(2,3) =  1.789135380979465422050817815017383d0
      rkTinvAinv(3,1) = -4.160768067752685525282947313530352d0
      rkTinvAinv(3,2) =  1.124128569859216916690209918405860d0
      rkTinvAinv(3,3) =  1.700644430961823796581896350418417d0

      rkAinvT(1,1) = 1.543510591072668287198054583233180d0
      rkAinvT(1,2) = -2.460228411937788329157493833295004d0
      rkAinvT(1,3) = -.412906170450356277003910443520499d0
      rkAinvT(2,1) = .209519643211838264029272585946993d0
      rkAinvT(2,2) = 1.388545667194387164417459732995766d0
      rkAinvT(2,3) = 1.20339553005832004974976023130002d0
      rkAinvT(3,1) = 3.637834252744495732208418513577775d0
      rkAinvT(3,2) = 2.681082873627752133895790743211112d0
      rkAinvT(3,3) = -3.050430199247410569426377624787569d0

    END SUBROUTINE Radau1A_Coefficients

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END SUBROUTINE RungeKuttaADJ2 ! and all its internal procedures
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE RK_adj_f90_Integrator
