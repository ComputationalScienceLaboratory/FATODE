! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! INTEGRATE - Integrator routine
!   Arguments :
!      TIN       - Start Time for Integration
!      TOUT      - End Time for Integration
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Discrete adjoints of Rosenbrock,                                       !
!           for several Rosenbrock methods:                               !
!               * Ros2                                                    !
!               * Ros3                                                    !
!               * Ros4                                                    !
!               * Rodas3                                                  !
!               * Rodas4                                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE ROS_adj_f90_Integrator

  IMPLICIT NONE
  PUBLIC
  SAVE

!~~~>  Statistics on the work performed by the Rosenbrock method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                      Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                      Ntexit=1, Nhexit=2, Nhnew =3, NHESS = 1

CONTAINS ! Routines in the module ROS_adj_f90_Integrator

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE INTEGRATE_ADJ( NVAR, NP, NADJ, NNZERO,  Y, Lambda, Mu, TIN, TOUT,  &
             ATOL_adj, RTOL_adj, ATOL, RTOL, FUN, JAC, ADJINIT, HESSTR_VEC,     &
             JACP, DRDY, DRDP, HESSTR_VEC_F_PY, HESSTR_VEC_R_PY, HESSTR_VEC_R,  &
             ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Q, QFUN ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IMPLICIT NONE        
    INTEGER, INTENT(IN) :: NVAR, NNZERO, NP
!~~~> Y - Concentrations
    DOUBLE PRECISION, INTENT(INOUT)  :: Y(NVAR)
!~~~> NADJ - No. of cost functionals for which adjoints
!                are evaluated simultaneously
!            If single cost functional is considered (like in
!                most applications) simply set NADJ = 1      
    INTEGER, INTENT(IN) :: NADJ
!~~~> Q - Quadrature term in the cost function
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: Q(NADJ)
!~~~> Lambda - Sensitivities w.r.t. concentrations
!     Note: Lambda (1:N,j) contains sensitivities of
!           the j-th cost functional w.r.t. Y(1:N), j=1...NADJ
    DOUBLE PRECISION, INTENT(INOUT)  :: Lambda(NVAR,NADJ)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: Mu(NP,NADJ)
    DOUBLE PRECISION, INTENT(IN)  :: TIN  ! TIN  - Start Time
    DOUBLE PRECISION, INTENT(IN)  :: TOUT ! TOUT - End Time
!~~~> Tolerances for adjoint calculations
!     (used only for full continuous adjoint)   
    DOUBLE PRECISION, INTENT(IN)  :: ATOL_adj(NVAR,NADJ), RTOL_adj(NVAR,NADJ),&
                                    ATOL(NVAR), RTOL(NVAR)
!~~~> Optional input parameters and statistics
    INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
    INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)

    DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20)
    INTEGER       :: ICNTRL(20), ISTATUS(20), IERR

    EXTERNAL FUN,JAC,ADJINIT,HESSTR_VEC
    EXTERNAL QFUN,JACP,DRDY,DRDP,HESSTR_VEC_F_PY,HESSTR_VEC_R_PY,HESSTR_VEC_R
    OPTIONAL QFUN,JACP,DRDY,DRDP,HESSTR_VEC_F_PY,HESSTR_VEC_R_PY,HESSTR_VEC_R

    ICNTRL(1:20)  = 0
    RCNTRL(1:20)  = 0.0
    ISTATUS(1:20) = 0
    RSTATUS(1:20) = 0.0
      
!~~~> fine-tune the integrator:
!   ICNTRL(1) = 0       ! 0 = non-autonomous, 1 = autonomous
!   ICNTRL(2) = 1       ! 0 = scalar, 1 = vector tolerances
!   RCNTRL(3) = STEPMIN ! starting step
!   ICNTRL(3) = 5       ! choice of the method for forward integration
!   ICNTRL(6) = 1       ! choice of the method for continuous adjoint
!   ICNTRL(7) = 2       ! 1=none, 2=discrete, 3=full continuous, 4=simplified continuous adjoint
!   ICNTRL(8) = 1       ! Save fwd LU factorization: 0 = *don't* save, 1 = save


    ! if optional parameters are given, and if they are >=0, then they overwrite default settings
    IF (PRESENT(ICNTRL_U)) THEN
      WHERE(ICNTRL_U(:) >= 0) ICNTRL(:) = ICNTRL_U(:)
    END IF
    IF (PRESENT(RCNTRL_U)) THEN
      WHERE(RCNTRL_U(:) >= 0) RCNTRL(:) = RCNTRL_U(:)
    END IF

! Evaluating sensitivities w.r.t parameters requires NP>0, functions MU and JACP
! are provided.
    IF(NP > 0 .AND. PRESENT(Mu) .AND. PRESENT(JACP) .AND. PRESENT(HESSTR_VEC_F_PY))THEN     
      IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDP) .AND. PRESENT(DRDY)&
        .AND. PRESENT(HESSTR_VEC_R_PY) .AND. PRESENT(HESSTR_VEC_R)) THEN
!~~~> This is the case that cost function contains a quadrature term NADJ 
!     should be 1; Q is defined to store the value of the quadrature
!     term at the last step and subroutines QFUN,DRDY,DRDP must be provided.
!     cost function = g(y,p) + Integrate{r(y,p)}
        CALL RosenbrockADJ1(N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, NADJ=NADJ,      &
              Lambda=Lambda, Mu=Mu, Tstart=TIN, Tend=TOUT, FUN=FUN, JAC=JAC,  &
              ADJINIT=ADJINIT, HESSTR_VEC=HESSTR_VEC, DRDP=DRDP, DRDY=DRDY,   &
              JACP=JACP, HESSTR_VEC_F_PY=HESSTR_VEC_F_PY,                     &
              HESSTR_VEC_R_PY=HESSTR_VEC_R_PY, HESSTR_VEC_R=HESSTR_VEC_R,     &
              AbsTol=ATOL, RelTol=RTOL, AbsTol_adj=ATOL_adj,                  &
              RelTol_adj=RTOL_adj, RCNTRL=RCNTRL, ICNTRL=ICNTRL,              &
              RSTATUS=RSTATUS, ISTATUS=ISTATUS, IERR=IERR, Q=Q, QFUN=QFUN)
      ELSE
!~~~> No quadrature term is involved. cost function = g(y,p)
        CALL RosenbrockADJ1( N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, NADJ=NADJ,     &
              Lambda=Lambda, Mu=Mu, Tstart=TIN, Tend=TOUT, RelTol=RTOL,       &
              AbsTol=Atol, RelTol_adj=RTOL_adj, AbsTol_adj=ATOL_adj,          &
              JACP=JACP, AdjInit=AdjInit, RCNTRL=RCNTRL, ICNTRL=ICNTRL,       &
              RSTATUS=RSTATUS, ISTATUS=ISTATUS, IERR=IERR, FUN=FUN, JAC=JAC,  &
              HESSTR_VEC=HESSTR_VEC, HESSTR_VEC_F_PY=HESSTR_VEC_F_PY )   
      END if
    ELSE
! Evaluating sensitivites w.r.t only initial conditions
      IF(PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDY) .AND.&
         PRESENT(HESSTR_VEC_R) )THEN
        CALL RosenbrockADJ2( N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, NADJ=NADJ,     &
               Lambda=Lambda, Tstart=TIN, Tend=TOUT, FUN=FUN, JAC=JAC,         &
               ADJINIT=ADJINIT, HESSTR_VEC=HESSTR_VEC, DRDY=DRDY,              &
               HESSTR_VEC_R=HESSTR_VEC_R, AbsTol=ATOL, RelTol=RTOL,            &
               AbsTol_adj=ATOL_adj, RelTol_adj=RTOL_adj, RCNTRL=RCNTRL,        &
               ICNTRL=ICNTRL, RSTATUS=RSTATUS, ISTATUS=ISTATUS, IERR=IERR, Q=Q,&
               QFUN=QFUN)
      ELSE
        PRINT *,'mode 4'
!~~~> No quadrature term is involved. cost function = g(y)        
        CALL RosenbrockADJ2( N=NVAR, NP=NP, NNZERO=NNZERO, Y=Y, NADJ=NADJ,     &
               Lambda=Lambda, Tstart=TIN, Tend=TOUT, FUN=FUN, JAC=JAC,         &
              ADJINIT=ADJINIT, HESSTR_VEC=HESSTR_VEC, AbsTol=ATOL, RelTol=RTOL,&
               AbsTol_adj=ATOL_adj, RelTol_adj=RTOL_adj, RCNTRL=RCNTRL,        &
               ICNTRL=ICNTRL, RSTATUS=RSTATUS, ISTATUS=ISTATUS, IERR=IERR )
      END IF
    END IF

    IF (IERR < 0) THEN
      print *,'RosenbrockADJ: Unsucessful step at T=', TIN,' (IERR=',IERR,')'
    END IF

    !   STEPMIN = RSTATUS(Nhexit)
    !   if optional parameters are given for output 
    !         copy to them to return information
    IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
    IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)

  END SUBROUTINE INTEGRATE_ADJ

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE RosenbrockADJ1(N, NP, NNZERO, Y, NADJ, Lambda, Tstart, Tend, FUN,  &
               JAC, ADJINIT, HESSTR_VEC, AbsTol, RelTol, AbsTol_adj, RelTol_adj,&
               RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR, Mu, JACP, Q, QFUN, DRDP, &
               DRDY, HESSTR_VEC_F_PY, HESSTR_VEC_R_PY, HESSTR_VEC_R)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   
!    ADJ = Adjoint of the Tangent Linear Model of a Rosenbrock Method
!
!    Solves the system y'=F(t,y) using a RosenbrockADJ method defined by:
!
!     G = 1/(H*gamma(1)) - Jac(t0,Y0)
!     T_i = t0 + Alpha(i)*H
!     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
!     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
!         gamma(i)*dF/dT(t0, Y0)
!     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j 
!
!    For details on RosenbrockADJ methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.  
!    The codes contained in the book inspired this implementation.       
!
!    (C)  Adrian Sandu, August 2004
!    Virginia Polytechnic Institute and State University    
!    Contact: sandu@cs.vt.edu
!    Revised by Hong Zhang and Adrian Sandu, Feb 2011          
!    This implementation is part of FATODE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    
!~~~>   INPUT ARGUMENTS: 
!    
!-     Y(N)    = vector of initial conditions (at T=Tstart)
!      NADJ       -> dimension of linearized system, 
!                   i.e. the number of sensitivity coefficients
!-     Lambda(N,NADJ) -> vector of initial sensitivity conditions (at T=Tstart)
!-    [Tstart,Tend]  = time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)  
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE Fun( T, Y, Ydot ) = ODE function, 
!                       returns Ydot = Y' = F(T,Y) 
!- SUBROUTINE Jac( T, Y, Jcb ) = Jacobian of the ODE function,
!                       returns Jcb = dF/dY 
!-    ICNTRL(1:10)    = integer inputs parameters
!-    RCNTRL(1:10)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:  
!     
!-    Y(N)    -> vector of final states (at T->Tend)
!-    Lambda(N,NADJ) -> vector of final sensitivities (at T=Tend)
!-    ICNTRL(11:20)   -> integer output parameters
!-    RCNTRL(11:20)   -> real output parameters
!-    IERR       -> job status upon return
!       - succes (positive value) or failure (negative value) -
!           =  1 : Success
!           = -1 : Improper value for maximal no of steps
!           = -2 : Selected RosenbrockADJ method not implemented
!           = -3 : Hmin/Hmax/Hstart must be positive
!           = -4 : FacMin/FacMax/FacRej must be positive
!           = -5 : Improper tolerance values
!           = -6 : No of steps exceeds maximum bound
!           = -7 : Step size too small
!           = -8 : Matrix is repeatedly singular
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1)   = 1: F = F(y)   Independent of T (AUTONOMOUS)
!              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)
!
!    ICNTRL(2)   = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1:  AbsTol, RelTol are scalars
!
!    ICNTRL(3)  -> selection of a particular Rosenbrock method
!        = 0 :  default method is Rodas3
!        = 1 :  method is  Ros2
!        = 2 :  method is  Ros3 
!        = 3 :  method is  Ros4 
!        = 4 :  method is  Rodas3
!        = 5:   method is  Rodas4
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0) the default value of BUFSIZE is used
!
!    ICNTRL(6)  -> selection of a particular Rosenbrock method for the
!                continuous adjoint integration - for cts adjoint it
!                can be different than the forward method ICNTRL(3)
!         Note 1: to avoid interpolation errors (which can be huge!) 
!                   it is recommended to use only ICNTRL(7) = 2 or 4
!         Note 2: the performance of the full continuous adjoint
!                   strongly depends on the forward solution accuracy Abs/RelTol
!
!    ICNTRL(7) -> Type of adjoint algorithm
!         = 0 : default is discrete adjoint ( of method ICNTRL(3) )
!         = 1 : no adjoint       
!         = 2 : discrete adjoint ( of method ICNTRL(3) )
!         = 3 : fully adaptive continuous adjoint ( with method ICNTRL(6) )
!         = 4 : simplified continuous adjoint ( with method ICNTRL(6) )
!
!    ICNTRL(8)  -> checkpointing the LU factorization at each step:
!        ICNTRL(8)=0 : do *not* save LU factorization (the default)
!        ICNTRL(8)=1 : save LU factorization
!        Note: if ICNTRL(7)=1 the LU factorization is *not* saved
!
!    ICNTRL(9)   = 1: integrand r(y)  Independent of T (AUTONOMOUS)
!              = 0: integrand r(t,y) Depends on T (NON-AUTONOMOUS)
!
!~~~>  Real input parameters:
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO 
!
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!          
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!            (default=0.1)
!
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller 
!         than the predicted value  (default=0.9)
!
!    RCNTRL(8)  -> ThetaMin. If Newton convergence rate smaller
!                  than ThetaMin the Jacobian is not recomputed;
!                  (default=0.001)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call to RosenbrockADJ adds the corrent no. of fcn calls
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
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    For multiple restarts, use Hexit as Hstart in the following run 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USE LS_Solver
    IMPLICIT NONE
   
!~~~>  Arguments
    INTEGER                      :: N, NP, NNZERO
    DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: Q(NADJ)
    INTEGER, INTENT(IN)          :: NADJ
    DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ), Mu(NP,NADJ)
    DOUBLE PRECISION, INTENT(IN)    :: Tstart,Tend
    DOUBLE PRECISION, INTENT(IN)    :: AbsTol(N),RelTol(N)
    DOUBLE PRECISION, INTENT(IN)    :: AbsTol_adj(N,NADJ), RelTol_adj(N,NADJ)
    INTEGER, INTENT(IN)          :: ICNTRL(20)
    DOUBLE PRECISION, INTENT(IN)    :: RCNTRL(20)
    INTEGER, INTENT(INOUT)       :: ISTATUS(20)
    DOUBLE PRECISION, INTENT(INOUT) :: RSTATUS(20)
    INTEGER, INTENT(OUT)         :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
    INTEGER ::  ros_S, rosMethod
    INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5
    DOUBLE PRECISION :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
    LOGICAL :: ros_NewF(6)
    CHARACTER(LEN=12) :: ros_Name
!~~~>  Types of Adjoints Implemented
    INTEGER, PARAMETER :: Adj_none = 1, Adj_discrete = 2,      &
                    Adj_continuous = 3, Adj_simple_continuous = 4
!~~~>  Checkpoints in memory
    INTEGER, PARAMETER :: bufsize = 200000
    INTEGER :: stack_ptr = 0 ! last written entry
    DOUBLE PRECISION, DIMENSION(:),   POINTER :: chk_H, chk_T
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_Y, chk_K
 
!   DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_J, chk_P
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_dY, chk_d2Y
!~~~>  Local variables     
    DOUBLE PRECISION :: Roundoff, FacMin, FacMax, FacRej, FacSafe
    DOUBLE PRECISION :: Hmin, Hmax, Hstart
    DOUBLE PRECISION :: Texit
    INTEGER :: i, UplimTol, Max_no_steps
    INTEGER :: AdjointType, CadjMethod 
    LOGICAL :: Autonomous, QuadAutonomous, VectorTol, SaveLU, GetQuad
!~~~>   Parameters
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0
    DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
    DOUBLE PRECISION :: DLAMCH
    EXTERNAL FUN,JAC,ADJINIT,JACP,HESSTR_VEC,HESSTR_VEC_F_PY
    EXTERNAL QFUN,DRDP,DRDY,HESSTR_VEC_R_PY,HESSTR_VEC_R
    OPTIONAL QFUN,DRDP,DRDY,HESSTR_VEC_R_PY,HESSTR_VEC_R

    ros_A = ZERO
    ros_C = ZERO
    ros_M = ZERO
    ros_E = ZERO
    ros_Alpha = ZERO
    ros_Gamma = ZERO
    GetQuad = .FALSE.
!~~~> Compute the quadrature term if required sources are provided
    IF( PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDP) .AND. PRESENT(DRDY) &
       .AND. PRESENT(HESSTR_VEC_R_PY) .AND.    &
       PRESENT(HESSTR_VEC_R) ) THEN
      GetQuad = .TRUE.
    END IF
  
!~~~>  Initialize statistics
    ISTATUS(1:20) = 0
    RSTATUS(1:20) = ZERO
   
!~~~>  Autonomous or time dependent ODE. Default is time dependent.
    Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
    IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
    ELSE 
      VectorTol = .FALSE.
      UplimTol  = 1
    END IF

!~~~>   Initialize the particular Rosenbrock method selected
    SELECT CASE (ICNTRL(3))
      CASE (1)
        CALL Ros2
      CASE (2)
        CALL Ros3
      CASE (3)
        CALL Ros4
      CASE (0,4)
        CALL Rodas3
      CASE (5)
        CALL Rodas4
      CASE DEFAULT
        PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
        CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
        RETURN
    END SELECT
   
!~~~>   The maximum number of steps admitted
    IF (ICNTRL(4) == 0) THEN
      Max_no_steps = bufsize - 1
    ELSEIF (Max_no_steps > 0) THEN
      Max_no_steps=ICNTRL(4)
    ELSE 
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN      
    END IF

!~~~>  The particular Rosenbrock method chosen for integrating the cts adjoint
    IF (ICNTRL(6) == 0) THEN
!      CadjMethod = 4
      CadjMethod = 3
    ELSEIF ( (ICNTRL(6) >= 1).AND.(ICNTRL(6) <= 5) ) THEN
      CadjMethod = ICNTRL(6)
    ELSE  
      PRINT * , 'Unknown CADJ Rosenbrock method: ICNTRL(6)=', CadjMethod
      CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN      
    END IF
 
!~~~>  Discrete or continuous adjoint formulation
    IF ( ICNTRL(7) == 0 ) THEN
      AdjointType = Adj_discrete
    ELSEIF ( (ICNTRL(7) >= 1).AND.(ICNTRL(7) <= 4) ) THEN
      AdjointType = ICNTRL(7)
    ELSE  
      PRINT * , 'User-selected adjoint type: ICNTRL(7)=', AdjointType
      CALL ros_ErrorMsg(-9,Tstart,ZERO,IERR)
      RETURN      
    END IF

!~~~> Save or not the forward LU factorization
    SaveLU = (ICNTRL(8) /= 0) 

!~~~> The quadrature term is autonomous or not 
    QuadAutonomous = .NOT.(ICNTRL(9) == 0)

!~~~>  Unit roundoff (1+Roundoff>1)  
    Roundoff = DLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
    IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
    ELSEIF (RCNTRL(1) > ZERO) THEN 
      Hmin = RCNTRL(1)
    ELSE  
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN      
    END IF
!~~~>  Upper bound on the step size: (positive value)
    IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
    ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
    ELSE  
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN      
    END IF
!~~~>  Starting step size: (positive value)
    IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
    ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
    ELSE  
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN      
    END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax 
    IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2d0
    ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
    ELSE  
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
    END IF
    IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0d0
    ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
    ELSE  
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
    END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
    IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1d0
    ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
    ELSE  
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
    END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
    IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9d0
    ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
    ELSE  
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
    END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.d0*Roundoff) &
         .OR. (RelTol(i) >= 1.0d0) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
     
!~~~>  Allocate checkpoint space or open checkpoint files
    IF (AdjointType == Adj_discrete) THEN
      CALL ros_AllocateDBuffers( ros_S )
    ELSEIF ( (AdjointType == Adj_continuous).OR.  &
           (AdjointType == Adj_simple_continuous) ) THEN
      CALL ros_AllocateCBuffers
    END IF
    CALL LSS_Init(N,NNZERO)
!~~~>  CALL Forward Rosenbrock method   
    CALL ros_FwdInt(N, NADJ, Y, Tstart, Tend, Texit, AbsTol, RelTol, GetQuad, & 
!  Error indicator
        IERR)

    PRINT*,'FORWARD STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

!~~~>  If Forward integration failed return   
    IF (IERR<0) RETURN

!~~~>   Initialize the particular Rosenbrock method for continuous adjoint
    IF ( (AdjointType == Adj_continuous).OR. &
           (AdjointType == Adj_simple_continuous) ) THEN
      SELECT CASE (CadjMethod)
        CASE (1)
          CALL Ros2
        CASE (2)
          CALL Ros3
        CASE (3)
          CALL Ros4
        CASE (4)
          CALL Rodas3
        CASE (5)
          CALL Rodas4
        CASE DEFAULT
          PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=', ICNTRL(3)
          CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR) 
          RETURN     
      END SELECT
    END IF
!~~~> Initialize adjoint variables
    CALL ADJINIT(N,NP,NADJ,Tend,Y,Lambda,Mu)

    SELECT CASE (AdjointType)   
      CASE (Adj_discrete)   
        CALL ros_DadjInt (                          &
           N, NADJ, NP, Lambda, Mu,                 &
           Tstart, Tend, Texit, GetQuad,            &
           IERR )
      CASE (Adj_continuous) 
        CALL ros_CadjInt (                          &
           NADJ, Lambda,                            &
           Tend, Tstart, Texit,                     &
           AbsTol_adj, RelTol_adj,                  &
           IERR )
      CASE (Adj_simple_continuous)
        CALL ros_SimpleCadjInt (                    &
           NADJ, Lambda,                            &
           Tstart, Tend, Texit,                     &
           IERR )
    END SELECT ! AdjointType

    PRINT*,'ADJOINT STATISTICS'
    PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

    CALL LSS_Free
!~~~>  Free checkpoint space or close checkpoint files
    IF (AdjointType == Adj_discrete) THEN
      CALL ros_FreeDBuffers
    ELSEIF ( (AdjointType == Adj_continuous) .OR. &
           (AdjointType == Adj_simple_continuous) ) THEN
      CALL ros_FreeCBuffers
    END IF
   

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CONTAINS !  Procedures internal to RosenbrockADJ1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_AllocateDBuffers( S )
!~~~>  Allocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: i, S
   
      ALLOCATE( chk_H(bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer H'; STOP
      END IF   
      ALLOCATE( chk_T(bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer T'; STOP
      END IF   
      ALLOCATE( chk_Y(N*S,bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer Y'; STOP
      END IF   
      ALLOCATE( chk_K(N*S,bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer K'; STOP
      END IF  

    END SUBROUTINE ros_AllocateDBuffers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_FreeDBuffers
!~~~>  Dallocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: i
   
      DEALLOCATE( chk_H, STAT=i )
      IF (i/=0) THEN
        PRINT*, 'Failed deallocation of buffer H'; STOP
      END IF   
      DEALLOCATE( chk_T, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer T'; STOP
      END IF   
      DEALLOCATE( chk_Y, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer Y'; STOP
      END IF   
      DEALLOCATE( chk_K, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer K'; STOP
      END IF   
!  SaveLU disabled
!   IF (SaveLU) THEN 
!     DEALLOCATE( chk_J, STAT=i )
!     IF (i/=0) THEN
!        PRINT*,'Failed deallocation of buffer J'; STOP
!     END IF
!     DEALLOCATE( chk_P, STAT=i )
!     IF (i/=0) THEN
!        PRINT*,'Fainled deallocation of buffer P'; STOP
!     END IF   
!   END IF   
 
    END SUBROUTINE ros_FreeDBuffers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_AllocateCBuffers
!~~~>  Allocate buffer space for continuous adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: i
      ALLOCATE( chk_H(bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer H'; STOP
      END IF   
      ALLOCATE( chk_T(bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer T'; STOP
      END IF   
      ALLOCATE( chk_Y(N,bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer Y'; STOP
      END IF   
      ALLOCATE( chk_dY(N,bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer dY'; STOP
      END IF   
      ALLOCATE( chk_d2Y(N,bufsize), STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed allocation of buffer d2Y'; STOP
      END IF   
    END SUBROUTINE ros_AllocateCBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_FreeCBuffers
!~~~>  Dallocate buffer space for continuous adjoint
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
      DEALLOCATE( chk_dY, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer dY'; STOP
      END IF   
      DEALLOCATE( chk_d2Y, STAT=i )
      IF (i/=0) THEN
        PRINT*,'Failed deallocation of buffer d2Y'; STOP
      END IF    
    END SUBROUTINE ros_FreeCBuffers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_DPush( S, T, H, Ystage, K, P )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER :: S ! no of stages
      DOUBLE PRECISION :: T, H, Ystage(N*S), K(N*S)
      INTEGER       :: P(N)
   
      stack_ptr = stack_ptr + 1
      IF ( stack_ptr > bufsize ) THEN
        PRINT*,'Push failed: buffer overflow'
        STOP
      END IF  
      chk_H( stack_ptr ) = H
      chk_T( stack_ptr ) = T
      !CALL DCOPY(N*S,Ystage,1,chk_Y(1,stack_ptr),1)
      !CALL DCOPY(N*S,K,1,chk_K(1,stack_ptr),1)
      chk_Y(1:N*S,stack_ptr) = Ystage(1:N*S)
      chk_K(1:N*S,stack_ptr) = K(1:N*S)
  
    END SUBROUTINE ros_DPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_DPop( S, T, H, Ystage, K )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      INTEGER :: S ! no of stages
      DOUBLE PRECISION :: T, H, Ystage(N*S), K(N*S)
      !   INTEGER       :: P(N)
   
      IF ( stack_ptr <= 0 ) THEN
        PRINT*,'Pop failed: empty buffer'
        STOP
      END IF  
      H = chk_H( stack_ptr )
      T = chk_T( stack_ptr )
      !CALL DCOPY(N*S,chk_Y(1,stack_ptr),1,Ystage,1)
      !CALL DCOPY(N*S,chk_K(1,stack_ptr),1,K,1)
      Ystage(1:N*S) = chk_Y(1:N*S,stack_ptr)
      K(1:N*S)      = chk_K(1:N*S,stack_ptr)
      !CALL DCOPY(LU_NONZERO,chk_J(1,stack_ptr),1,Jcb,1)

      stack_ptr = stack_ptr - 1
  
    END SUBROUTINE ros_DPop
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_CPush( T, H, Y, dY, d2Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DOUBLE PRECISION :: T, H, Y(N), dY(N), d2Y(N)
   
      stack_ptr = stack_ptr + 1
      IF ( stack_ptr > bufsize ) THEN
        PRINT*,'Push failed: buffer overflow'
        STOP
      END IF  
      chk_H( stack_ptr ) = H
      chk_T( stack_ptr ) = T
      !CALL DCOPY(N,Y,1,chk_Y(1,stack_ptr),1)
      !CALL DCOPY(N,dY,1,chk_dY(1,stack_ptr),1)
      !CALL DCOPY(N,d2Y,1,chk_d2Y(1,stack_ptr),1)
      chk_Y(1:N,stack_ptr)   = Y(1:N)
      chk_dY(1:N,stack_ptr)  = dY(1:N)
      chk_d2Y(1:N,stack_ptr) = d2Y(1:N)
    END SUBROUTINE ros_CPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_CPop( T, H, Y, dY, d2Y )
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
      !CALL DCOPY(N,chk_Y(1,stack_ptr),1,Y,1)
      !CALL DCOPY(N,chk_dY(1,stack_ptr),1,dY,1)
      !CALL DCOPY(N,chk_d2Y(1,stack_ptr),1,d2Y,1)
      Y(1:N)   = chk_Y(1:N,stack_ptr)
      dY(1:N)  = chk_dY(1:N,stack_ptr)
      d2Y(1:N) = chk_d2Y(1:N,stack_ptr)

      stack_ptr = stack_ptr - 1
  
    END SUBROUTINE ros_CPop

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
      DOUBLE PRECISION, INTENT(IN) :: T, H
      INTEGER, INTENT(IN)  :: Code
      INTEGER, INTENT(OUT) :: IERR
   
      IERR = Code
      PRINT * , &
      'Forced exit from RosenbrockADJ due to the following error:' 
     
      SELECT CASE (Code)
        CASE (-1)    
          PRINT * , '--> Improper value for maximal no of steps'
        CASE (-2)    
          PRINT * , '--> Selected RosenbrockADJ method not implemented'
        CASE (-3)    
          PRINT * , '--> Hmin/Hmax/Hstart must be positive'
        CASE (-4)    
          PRINT * , '--> FacMin/FacMax/FacRej must be positive'
        CASE (-5) 
          PRINT * , '--> Improper tolerance values'
        CASE (-6) 
          PRINT * , '--> No of steps exceeds maximum buffer bound'
        CASE (-7) 
          PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
        CASE (-8)    
          PRINT * , '--> Matrix is repeatedly singular'
        CASE (-9)    
          PRINT * , '--> Improper type of adjoint selected'
        CASE DEFAULT
          PRINT *, 'Unknown Error code: ', Code
        END SELECT
   
        PRINT *, "T=", T, "and H=", H
     
    END SUBROUTINE ros_ErrorMsg
   
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_FwdInt ( N, NADJ, Y,             &
           Tstart, Tend, T,                         &
           AbsTol, RelTol, GetQuad,                 &
!~~~> Error indicator
           IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockADJ method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NADJ
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval   
      DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
      DOUBLE PRECISION, INTENT(OUT) ::  T      
!~~~> Input: tolerances      
      DOUBLE PRECISION, INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Output: Error indicator
      INTEGER, INTENT(OUT) :: IERR
!~~~> Quadrature term indicator
      LOGICAL :: GetQuad
! ~~~~ Local variables        
      DOUBLE PRECISION :: Ynew(N), Fcn0(N), Fcn(N) 
      DOUBLE PRECISION :: K(N*ros_S), dFdT(N)
      DOUBLE PRECISION, DIMENSION(:), POINTER :: Ystage
      DOUBLE PRECISION :: H, Hnew, HC, HG, Fac, Tau 
      DOUBLE PRECISION :: Err, Yerr(N)
      INTEGER :: Pivot(N), Direction, ioffset, i, j, istage
      LOGICAL :: RejectLastH, RejectMoreH, Singular, Transp
!~~~~~ For Quadrature term
      DOUBLE PRECISION, ALLOCATABLE :: L(:),dRdT(:),R(:),RY(:,:) 
      DOUBLE PRECISION :: Delta
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IF(GetQuad) THEN
        ALLOCATE(L(NADJ*ros_S),dRdT(NADJ),R(NADJ),RY(N,NADJ),STAT=i)
        IF(i .NE. 0) STOP 'allocation error for L,dRdT,R,RY'
      END IF

!~~~>  Allocate stage vector buffer if needed
      IF (AdjointType == Adj_discrete) THEN
        ALLOCATE(Ystage(N*ros_S), STAT=i)
        ! Uninitialized Ystage may lead to NaN on some compilers
        Ystage = 0.0d0
        IF (i/=0) THEN
          PRINT*,'Allocation of Ystage failed'
          STOP
        END IF
      END IF   
   
!~~~>  Initial preparations
      T = Tstart
      RSTATUS(Nhexit) = ZERO
      H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
      IF (ABS(H) <= 10.0*Roundoff) H = DeltaMin

      IF (Tend  >=  Tstart) THEN
        Direction = +1
      ELSE
        Direction = -1
      END IF
      H = Direction*H

      RejectLastH=.FALSE.
      RejectMoreH=.FALSE.
   
!~~~> Time loop begins below 

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) ) 
      
        IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
          CALL ros_ErrorMsg(-6,T,H,IERR)
          RETURN
        END IF
        IF ( ((T+0.1d0*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
          CALL ros_ErrorMsg(-7,T,H,IERR)
          RETURN
        END IF
   
!~~~>  Limit H if necessary to avoid going beyond Tend   
        RSTATUS(Nhexit) = H
        H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
        CALL FUN(N,T,Y,Fcn0)
        ISTATUS(Nfun) = ISTATUS(Nfun) + 1 

!~~~>  Compute the function derivative with respect to T
        IF (.NOT.Autonomous) THEN
          CALL ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT )
        END IF
  
!~~~>   Compute the Jacobian at current time
        CALL LSS_Jac(T,Y,JAC)
        ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO  
   
          CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1),Singular)

          IF (Singular) THEN ! More than 5 consecutive failed decompositions
            CALL ros_ErrorMsg(-8,T,H,IERR)
            RETURN
          END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S
      
   ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
            ioffset = N*(istage-1)
      
   ! For the 1st istage the function has been computed previously
            IF ( istage == 1 ) THEN
              CALL DCOPY(N,Fcn0,1,Fcn,1)
              IF (AdjointType == Adj_discrete) THEN ! Save stage solution
              ! CALL DCOPY(N,Y,1,Ystage(1),1)
                Ystage(1:N) = Y(1:N)
                CALL DCOPY(N,Y,1,Ynew,1)
              END IF   
   ! istage>1 and a new function evaluation is needed at the current istage
            ELSEIF ( ros_NewF(istage) ) THEN
              CALL DCOPY(N,Y,1,Ynew,1)
              DO j = 1, istage-1
                CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
                K(N*(j-1)+1),1,Ynew,1) 
              END DO
              Tau = T + ros_Alpha(istage)*Direction*H
              CALL FUN(N,Tau,Ynew,Fcn)
              ISTATUS(Nfun) = ISTATUS(Nfun) + 1
            END IF ! if istage == 1 elseif ros_NewF(istage)
    ! save stage solution every time even if ynew is not updated
            IF ( ( istage > 1 ).AND.(AdjointType == Adj_discrete) ) THEN
    ! CALL DCOPY(N,Ynew,1,Ystage(ioffset+1),1)
              Ystage(ioffset+1:ioffset+N) = Ynew(1:N)
            END IF   
            CALL DCOPY(N,Fcn,1,K(ioffset+1),1)
            DO j = 1, istage-1
              HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
              CALL DAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
            END DO
            IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
              HG = Direction*H*ros_Gamma(istage)
              CALL DAXPY(N,HG,dFdT,1,K(ioffset+1),1)
            END IF
            Transp = .FALSE.
            CALL LSS_Solve(Transp, K(ioffset+1))
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1
      
          END DO Stage     
            
!~~~>  Compute the new solution 
          CALL DCOPY(N,Y,1,Ynew,1)
          DO j=1,ros_S
            CALL DAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
          END DO

!~~~>  Compute the error estimation 
          CALL DSCAL(N,ZERO,Yerr,1)
          DO j=1,ros_S     
            CALL DAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
          END DO 
          Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
          Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
          Hnew = H*Fac  

!~~~>  Check the error magnitude and adjust step size
          ISTATUS(Nstp) = ISTATUS(Nstp) + 1
          IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
  
            ISTATUS(Nacc) = ISTATUS(Nacc) + 1
            IF (AdjointType == Adj_discrete) THEN ! Save current state
              CALL ros_DPush( ros_S, T, H, Ystage, K,  Pivot )
            ELSEIF ( (AdjointType == Adj_continuous) .OR. &
            (AdjointType == Adj_simple_continuous) ) THEN
              CALL LSS_Mul_Jac(K,Fcn0)
              IF (.NOT. Autonomous) THEN
                CALL DAXPY(N,ONE,dFdT,1,K(1),1)
              END IF   
              CALL ros_CPush( T, H, Y, Fcn0, K(1) )
            END IF ! if(AdjointType 
!~~~> Update the results for the quadrature term
            IF(GetQuad) THEN
              CALL DRDY(NADJ,N,N,T,Y,RY)
              CALL QFUN(N,NADJ,T,Y,R)
              IF(.NOT.Autonomous) THEN
                Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
                CALL QFUN(N,NADJ,T+Delta,Y,dRdT)
                CALL DAXPY(NADJ,(-ONE),R,1,dRdT,1)
                CALL DSCAL(NADJ,(ONE/Delta),dRdT,1)
              END IF
!~~~>   Compute the stages for quadrature term
              DO istage = 1, ros_S
! Current istage offset. Current istage vector is L(ioffset+1:ioffset+NADJ)
                ioffset = NADJ*(istage-1)
! For the 1st istage the function has been computed previously
                IF ( istage == 1 ) THEN
                
! istage>1 and a new function evaluation is needed at the current istage
                ELSEIF ( ros_NewF(istage) ) THEN
                  Tau = T + ros_Alpha(istage)*Direction*H
                  CALL QFUN(N, NADJ,Tau,Ystage(N*(istage-1)+1),R)
                END IF ! if istage == 1 elseif ros_NewF(istage)
       
                CALL DCOPY(NADJ,R,1,L(ioffset+1),1)
                DO j = 1, istage-1
                  HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
                  CALL DAXPY(NADJ,HC,L(NADJ*(j-1)+1),1,L(ioffset+1),1)
                END DO

                IF ((.NOT. Autonomous) .AND. (ros_Gamma(istage).NE.ZERO)) THEN
                  HG = Direction*H*ros_Gamma(istage)
                  CALL DAXPY(NADJ,HG,dRdT,1,L(ioffset+1),1)
                END IF
! +r_y*k_i           
                DO i=1,NADJ
                  DO j=1,N
                    L(ioffset+i) = L(ioffset+i)+RY(j,i)*K(N*(istage-1)+j)
                  END DO
                END DO
                CALL DSCAL(NADJ,Direction*H*ros_Gamma(1),L(ioffset+1),1)
              END DO ! istage=1,ros_S 
!~~~>  Compute the new solution 
              DO j=1,ros_S
                CALL DAXPY(NADJ,ros_M(j),L(NADJ*(j-1)+1),1,Q,1)
              END DO
            END IF ! IF(GetQuad)
     
            CALL DCOPY(N,Ynew,1,Y,1)
            T = T + Direction*H
            Hnew = MAX(Hmin,MIN(Hnew,Hmax))
            IF (RejectLastH) THEN  ! No step size increase after a rejected step
              Hnew = MIN(Hnew,H) 
            END IF   
            RSTATUS(Nhexit) = H
            RSTATUS(Nhnew)  = Hnew
            RSTATUS(Ntexit) = T
            RejectLastH = .FALSE.  
            RejectMoreH = .FALSE.
            H = Hnew      
            EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
          ELSE           !~~~> Reject step
            IF (RejectMoreH) THEN
              Hnew = H*FacRej
            END IF   
            RejectMoreH = RejectLastH
            RejectLastH = .TRUE.
            H = Hnew
            IF (ISTATUS(Nacc) >= 1) THEN
              ISTATUS(Nrej) = ISTATUS(Nrej) + 1
            END IF    
          END IF ! Err <= 1

        END DO UntilAccepted 

      END DO TimeLoop 
   
!~~~> Save last state: only needed for continuous adjoint
      IF ( (AdjointType == Adj_continuous) .OR. &
      (AdjointType == Adj_simple_continuous) ) THEN
        CALL FUN(N,T,Y,Fcn0)
        ISTATUS(Nfun) = ISTATUS(Nfun) + 1
        CALL LSS_Jac(T,Y,JAC)
        ISTATUS(Njac) = ISTATUS(Njac) + 1

        CALL LSS_Mul_Jac(K, Fcn0)
        IF (.NOT. Autonomous) THEN
          CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
          CALL DAXPY(N,ONE,dFdT,1,K(1),1)
        END IF   
        CALL ros_CPush( T, H, Y, Fcn0, K(1) )
!~~~> Deallocate stage buffer: only needed for discrete adjoint
      ELSEIF (AdjointType == Adj_discrete) THEN 
        DEALLOCATE(Ystage, STAT=i)
        IF (i/=0) THEN
          PRINT*,'Deallocation of Ystage failed'
          STOP
        END IF
      END IF   

      IF(GetQuad) THEN
        DEALLOCATE(L,dRdT,R,RY,STAT=i)
        IF(i .NE. 0) STOP 'deallocation error for L, dRdT, R, RY'
      END IF   
!~~~> Succesful exit
      IERR = 1  !~~~> The integration was successful

    END SUBROUTINE ros_FwdInt
   
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_DadjInt (                        &
           N, NADJ, NP, Lambda, Mu,                 &
           Tstart, Tend, T, GetQuad,                &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockSOA method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)     :: N, NADJ, NP
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ),Mu(NP,NADJ)
!!~~~> Input: integration interval   
      DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
      DOUBLE PRECISION, INTENT(OUT) ::  T      
!~~~> Quadrature term indicator
      LOGICAL, INTENT(IN) :: GetQuad 
!~~~> Output: Error indicator
      INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables 
      DOUBLE PRECISION, ALLOCATABLE :: FPJAC0(:,:),FPJAC1(:,:),RP0(:,:),RP1(:,:),&
                                   RY0(:,:),RY1(:,:),U(:,:),V(:,:)        
      DOUBLE PRECISION :: Ystage(N*ros_S), K(N*ros_S), ros_W(ros_S)
      DOUBLE PRECISION :: Tmp(N), Tmp2(N), Tmp3(NP)
      DOUBLE PRECISION :: H, HC, HA, Tau, Alpha, Beta, GammaUhat
      INTEGER :: Direction
      INTEGER :: j, m, istage, istart, jstart
!~~~>  Local parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
      DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
      LOGICAL :: Transp
      DOUBLE PRECISION :: Delta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Initialization
      IF (GetQuad) THEN
        ALLOCATE(RP0(NP,NADJ),RP1(NP,NADJ),RY0(N,NADJ),RY1(N,NADJ),STAT=j)
        IF(j .NE. 0) STOP 'allocation error for RPs,RYs'
        RP0(:,:) = ZERO
        RP1(:,:) = ZERO
        RY0(:,:) = ZERO
        RY1(:,:) = ZERO
      END IF

      ALLOCATE(FPJAC0(N,NP),FPJAC1(N,NP),U(N*ros_S,NADJ), V(N*ros_S,NADJ),STAT=j)
      IF(j .NE. 0) STOP 'allocation error for FPJACs,U,V'
      FPJAC0(:,:) = ZERO
      FPJAC1(:,:) = ZERO

      Transp = .TRUE.
      IF (Tend  >=  Tstart) THEN
        Direction = +1
      ELSE
        Direction = -1
      END IF               

!~~~> calculate ros_W first
      IF(GetQuad) THEN
        GammaUhat = ZERO
        DO istage = ros_S, 1    
          ros_W(istage) = H*ros_Gamma(1)*ros_M(istage)
          DO j = istage+1,ros_S
            HC = ros_C((j-1)*(j-2)/2+istage)
            ros_W(istage) = ros_W(istage)+ HC*ros_Gamma(1)*ros_W(j)
          END DO
          GammaUhat = GammaUhat + ros_W(istage)*ros_Gamma(istage)
        END DO
      END IF
   
!~~~> Time loop begins below 
TimeLoop: DO WHILE ( stack_ptr > 0 )
        
!~~~>  Recover checkpoints for stage values and vectors
        CALL ros_DPop( ros_S, T, H, Ystage, K )

        ISTATUS(Nstp) = ISTATUS(Nstp) + 1

!~~~>    Compute LU decomposition 
        IF (.NOT.SaveLU) THEN
          CALL LSS_Jac(T,Ystage(1),JAC)
          ISTATUS(Njac) = ISTATUS(Njac) + 1
          Tau = ONE/(Direction*H*ros_Gamma(1))
          ISTATUS(Njac) = ISTATUS(Njac) + 1
     
          CALL LSS_Decomp(Tau,j)
          ISTATUS(Ndec) = ISTATUS(Ndec) + 1
        END IF
   
        IF(GetQuad)CALL DRDY(NADJ,N,N,T,Ystage(1),RY0)        

!~~~>   Compute Hessian at the beginning of the interval
!   CALL HESS(N, T,Ystage(1),Hes0)
        Tmp(:)  = ZERO
        Tmp3(:) = ZERO
!~~~>   Compute the stages

Stage: DO istage = ros_S, 1, -1
      
!~~~> Current istage first entry 
          istart = N*(istage-1) + 1
     
!~~~> Compute U 
          ! U = m_i * lambda_{n+1} + w_i drdy^T
          DO m = 1,NADJ
            CALL DCOPY(N,Lambda(1,m),1,U(istart,m),1)
            CALL DSCAL(N,ros_M(istage),U(istart,m),1)
            IF(GetQuad) CALL DAXPY(N,ros_W(istage),RY0,1,U(istart,m),1)
          END DO ! m=1:NADJ
          ! V=  sum a_{ji}(jac^T * u_j +w_j * drdy^T) stored in Tmp
          ! U = U+V+sum c_{ji}/h u_j
          DO j = istage+1,ros_S
            jstart = N*(j-1) + 1
            HA = ros_A((j-1)*(j-2)/2+istage)
            HC = ros_C((j-1)*(j-2)/2+istage)/(Direction*H)
            DO m = 1,NADJ
              CALL DAXPY(N,HA,V(jstart,m),1,U(istart,m),1) 
              CALL DAXPY(N,HC,U(jstart,m),1,U(istart,m),1)
            END DO ! m=1:NADJ            
          END DO 
          DO m = 1,NADJ
            CALL LSS_Solve(Transp,U(istart,m))
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1
          END DO ! m=1:NADJ
          
!~~~> Compute V 
         ! V=dfdy^t*u+drdy^T*w 
          Tau = T + ros_Alpha(istage)*Direction*H
!       CALL JACO(N, Tau,Ystage(istart),Jac)
!       ~~~~> jac0 jac1 exchanged
          CALL LSS_Jac(Tau,Ystage(istart),JAC)
          ISTATUS(Njac) = ISTATUS(Njac) + 1
          IF(GetQuad) CALL DRDY(NADJ,N,N,Tau,Ystage(istart),RY1)
          DO m = 1,NADJ
            CALL LSS_Mul_Jactr(V(istart,m),U(istart,m))
            IF(GetQuad)CALL DAXPY(N,ros_W(istage),RY1(1,m),1,V(istart,m),1)
          END DO ! m=1:NADJ
          
        END DO Stage
            
!~~~>  Update the sensitivity variables

!~~~>  Compute Lambda  
        DO istage=1,ros_S
          istart = N*(istage-1) + 1
          DO m = 1,NADJ
            ! Add V_i
            CALL DAXPY(N,ONE,V(istart,m),1,Lambda(1,m),1)
            ! Add (H0xK_i)^T * U_i
            CALL HESSTR_VEC(N,T,Ystage(1),U(istart,m),K(istart),Tmp)
!           CALL HessTR_Vec ( NHESS, Hes0, U(istart,m), K(istart), Tmp )
            CALL DAXPY(N,ONE,Tmp,1,Lambda(1,m),1) 
            ! Add (HP0xK_i)^T * U_i
            IF(GetQuad) THEN
              CALL HESSTR_VEC_R(m,N,T,Ystage(1),ros_W(istage),K(istart),Tmp)
              CALL DAXPY(N,ONE,Tmp,1,Lambda(1,m),1)
            END IF
          END DO ! m=1:NADJ
        END DO
     
!~~~>  Compute Mu
        Beta = 1.0d0
        DO istage = 1,ros_S
          istart = N*(istage-1) + 1
         ! Calculate v_bar(i) = (dfdp)^T * u(i) + ros_W(i) * (drdp)^T 
          Tau = T + ros_Alpha(istage)*Direction*H
         ! dfdp(T_i,Y_i) of dimension N x NP
          CALL JACP(N,NP,Tau,Ystage(istart),FPJAC1)
         ! drdp(T_i,Y_i) of dimension NP x NADJ
          IF(GetQuad) CALL DRDP(NADJ,N,NP,Tau,Ystage(istart),RP1)
          DO m = 1,NADJ
            Tmp3(1:NP) = ZERO
           ! Add dfdp^T u_i+drdp^T w_i
            CALL DGEMV('T',N,NP,ONE,FPJAC1,N,U(istart,m),1,Beta,Tmp3,1)
            CALL DAXPY(NP,ONE,Tmp3,1,Mu(1,m),1) 
            IF(GetQuad) THEN
              CALL DAXPY(NP,ros_W(istage),RP1(1,m),1,Tmp3,1)
              ! Add v_bar (Tmp3)
              ! CALL DAXPY(N,ONE,Tmp3,1,Mu(1,m),1) omg!!!!!!
              CALL DAXPY(NP,ONE,Tmp3,1,Mu(1,m),1)
            END IF
             ! Add dfdpdy term
            CALL HESSTR_VEC_F_PY(N,NP, T,Ystage(1),U(istart,m),K(istart),Tmp3)
            ! CALL DAXPY(N,ONE,Tmp3,1,Mu(1,m),1) omg
            CALL DAXPY(NP,ONE,Tmp3,1,Mu(1,m),1) 
            IF(GetQuad) THEN
             ! Add drdpdy term
              CALL HESSTR_VEC_R_PY(m,N,NP,T,Ystage(1),ros_W(istage),K(istart),Tmp3)
              CALL DAXPY(NP,ONE,Tmp3,1,Mu(1,m),1)
            END IF
          END DO ! m=1:NADJ
        END DO

!~~~>  NonAutonomous term
 
        ! Tmp stores sum gamma_i U_i
        IF (.NOT.Autonomous) THEN
          ! Approximate d(f_y)/dt 
          Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
          CALL LSS_Jac2(T+Delta,Ystage(1),JAC)
          ISTATUS(Njac) = ISTATUS(Njac) + 1
          CALL LSS_Jac_Time_Derivative(Delta)
 
          ! Approximate d(f_p)/dt
          CALL JACP(N,NP,T,Ystage(1),FPJAC0)
          CALL JACP(N,NP,T+Delta,Ystage(1),FPJAC1)
          Alpha = -1.0d0
          CALL DAXPY(N*NP,Alpha,FPJAC0,1,FPJAC1,1)
          Alpha = 1/Delta
          CALL DSCAL(N*NP,Alpha,FPJAC1,1)

          DO m = 1,NADJ
            Tmp(1:N) = ZERO
            DO istage = 1, ros_S
              istart = N*(istage-1) + 1
              CALL DAXPY(N,ros_Gamma(istage),U(istart,m),1,Tmp,1)
            END DO  
            ! Add H * df/(dpdt)_0^T * \sum(gamma_i U_i) to Mu
            Alpha = H
            CALL DGEMV('T',N,NP,Alpha,FPJAC1,N,Tmp,1,Beta,Mu(1,m),1)
            ! Add H * dJac_dT_0^T * \sum(gamma_i U_i) to Lambda
            CALL LSS_Mul_djdttr(Tmp2,Tmp)
            CALL DAXPY(N,H,Tmp2,1,Lambda(1,m),1)
          END DO ! m=1:NADJ

        END IF ! .NOT.Autonomous
        
        IF(GetQuad .AND. .NOT. QuadAutonomous) THEN
          ! Approximate d(r_y)/dt) and add it to Lambda 
          CALL DRDY(NADJ,N,N,T+Delta,Ystage(1),RY1)
          Alpha = -1.0d0
          CALL DAXPY(N*NADJ,Alpha,RY0,1,RY1,1)
          Alpha = 1/Delta
          CALL DSCAL(N*NADJ,Alpha,RY1,1)
          Alpha = H*GammaUhat
          CALL DAXPY(N*NADJ,Alpha,RY1,1,Lambda,1) 
    
          ! Approximate d(r_p)/dt and add it to Mu
          CALL DRDP(NADJ,N,NP,T,Ystage(1),RP0)
          CALL DRDP(NADJ,N,NP,T+Delta,Ystage(1),RP1)
          Alpha = -1.0d0
          CALL DAXPY(NADJ*NP,Alpha,RP0,1,RP1,1)
          Alpha = 1/Delta
          CALL DSCAL(NADJ*NP,Alpha,RP1,1)
          Alpha = H*GammaUhat
          CALL DAXPY(NP*NADJ,Alpha,RP1,1,Mu,1)
        END IF
  
      END DO TimeLoop 
   
      DEALLOCATE(FPJAC0,FPJAC1,U,V,STAT=j)
      IF(j .NE. 0) STOP 'deallocation error for FPJACs,RPs,RYs,U,V'
      IF(GetQuad) THEN
        DEALLOCATE(RP0,RP1,RY0,RY1,STAT=j)
        IF(j .NE. 0) STOP 'deallocation error for RPs,RYs'
      END IF
!~~~> Succesful exit
      IERR = 1  !~~~> The integration was successful

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END SUBROUTINE ros_DadjInt
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_CadjInt (                        &
           NADJ, Y,                                 &
           Tstart, Tend, T,                         &
           AbsTol_adj, RelTol_adj,                  &
!~~~> Error indicator
           IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockADJ method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
   INTEGER, INTENT(IN) :: NADJ
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
!~~~> Input: integration interval   
   DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend      
!~~~> Input: adjoint tolerances   
   DOUBLE PRECISION, INTENT(IN) :: AbsTol_adj(N,NADJ), RelTol_adj(N,NADJ)
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   DOUBLE PRECISION, INTENT(OUT) ::  T      
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables        
   DOUBLE PRECISION :: Y0(N)
   DOUBLE PRECISION :: Ynew(N,NADJ), Fcn0(N,NADJ), Fcn(N,NADJ) 
   DOUBLE PRECISION :: K(N*ros_S,NADJ), dFdT(N,NADJ)
   DOUBLE PRECISION :: H, Hnew, HC, HG, Fac, Tau 
   DOUBLE PRECISION :: Err, Yerr(N,NADJ)
   INTEGER :: Direction, ioffset, j, istage,iadj
   LOGICAL :: RejectLastH, RejectMoreH, Singular, Transp
!~~~>  Local parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
   DOUBLE PRECISION :: Delta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = 0.0
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0*Roundoff) H = DeltaMin
   
   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF               
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.
   
!~~~> Time loop begins below 

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) ) 
      
   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1d0*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   
!~~~>  Limit H if necessary to avoid going beyond Tend   
   RSTATUS(Nhexit) = H
   H = MIN(H,ABS(Tend-T))

!~~~>   Interpolate forward solution
   CALL ros_cadj_Y( T, Y0 )     
!~~~>   Compute the Jacobian at current time
   CALL LSS_Jac(T,Y0,JAC)
   ISTATUS(Njac) = ISTATUS(Njac) + 1
   
!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
       Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
       CALL LSS_Jac2(T+Delta,Y0,JAC)
       ISTATUS(Njac) = ISTATUS(Njac) + 1
       CALL LSS_Jac_Time_Derivative(Delta)
      DO iadj = 1, NADJ
        CALL LSS_Mul_djdttr(dFdT(1,iadj),Y(1,iadj))
        CALL DSCAL(N,(-ONE),dFdT(1,iadj),1)
      END DO
   END IF

   CALL LSS_Rev_Jac()
!~~~>  Ydot = -J^T*Y
  
   DO iadj = 1, NADJ
     CALL LSS_Mul_Jactr(Fcn0(1,iadj),Y(1,iadj))
   END DO
    
!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO  
   
   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S
      
      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)
      
      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         DO iadj = 1, NADJ
           CALL DCOPY(N,Fcn0(1,iadj),1,Fcn(1,iadj),1)
         END DO
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL DCOPY(N*NADJ,Y,1,Ynew,1)
         DO j = 1, istage-1
           DO iadj = 1, NADJ
             CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
                K(N*(j-1)+1,iadj),1,Ynew(1,iadj),1) 
           END DO       
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL ros_cadj_Y( Tau, Y0 )     
         CALL LSS_Jac1(Tau,Y0,JAC)
         ISTATUS(Njac) = ISTATUS(Njac) + 1
         CALL LSS_Rev_Jac1()
         DO iadj = 1, NADJ
             CALL LSS_Mul_Jactr1(Fcn(1,iadj),Ynew(1,iadj))
         END DO
       END IF ! if istage == 1 elseif ros_NewF(istage)

       DO iadj = 1, NADJ
          CALL DCOPY(N,Fcn(1,iadj),1,K(ioffset+1,iadj),1)
       END DO
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HC,K(N*(j-1)+1,iadj),1, &
                  K(ioffset+1,iadj),1)
         END DO
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HG,dFdT(1,iadj),1,K(ioffset+1,iadj),1)
         END DO
       END IF
       
       Transp = .TRUE.
       DO iadj = 1, NADJ
          CALL LSS_Solve(Transp,K(ioffset+1,iadj))
          ISTATUS(Nsol) = ISTATUS(Nsol) + 1
       END DO
      
   END DO Stage     
            

!~~~>  Compute the new solution 
   DO iadj = 1, NADJ
      CALL DCOPY(N,Y(1,iadj),1,Ynew(1,iadj),1)
      DO j=1,ros_S
         CALL DAXPY(N,ros_M(j),K(N*(j-1)+1,iadj),1,Ynew(1,iadj),1)
      END DO
   END DO

!~~~>  Compute the error estimation 
   CALL DSCAL(N*NADJ,ZERO,Yerr,1)
   DO j=1,ros_S     
       DO iadj = 1, NADJ
        CALL DAXPY(N,ros_E(j),K(N*(j-1)+1,iadj),1,Yerr(1,iadj),1)
       END DO
   END DO
!~~~> Max error among all adjoint components    
   iadj = 1
   Err = ros_ErrorNorm ( Y(1,iadj), Ynew(1,iadj), Yerr(1,iadj), &
              AbsTol_adj(1,iadj), RelTol_adj(1,iadj), VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac  

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      CALL DCOPY(N*NADJ,Ynew,1,Y,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H) 
      END IF   
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.  
      RejectMoreH = .FALSE.
      H = Hnew      
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF   
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1) THEN
         ISTATUS(Nrej) = ISTATUS(Nrej) + 1
      END IF    
   END IF ! Err <= 1

   END DO UntilAccepted 

   END DO TimeLoop 
      
!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_CadjInt
  
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_SimpleCadjInt (                  &
        NADJ, Y,                                 &
        Tstart, Tend, T,                         &
!~~~> Error indicator
        IERR )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockADJ method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE LS_Solver
  IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
   INTEGER, INTENT(IN) :: NADJ
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
!~~~> Input: integration interval   
   DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   DOUBLE PRECISION, INTENT(OUT) ::  T      
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables        
   DOUBLE PRECISION :: Y0(N)
   DOUBLE PRECISION :: Ynew(N,NADJ), Fcn0(N,NADJ), Fcn(N,NADJ) 
   DOUBLE PRECISION :: K(N*ros_S,NADJ), dFdT(N,NADJ)
   DOUBLE PRECISION :: H, HC, HG, Tau 
   DOUBLE PRECISION :: ghinv
   INTEGER :: Direction, ioffset, j, istage, iadj
   INTEGER :: istack
!~~~>  Local parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
   DOUBLE PRECISION :: Delta
   Logical :: Transp
!~~~>  Locally called functions
!    DOUBLE PRECISION WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
!~~~>  INITIAL PREPARATIONS
   
   IF (Tend  >=  Tstart) THEN
     Direction = -1
   ELSE
     Direction = +1
   END IF               
   
!~~~> Time loop begins below 
TimeLoop: DO istack = stack_ptr,2,-1
        
   T = chk_T(istack)
   H = chk_H(istack-1)
   !CALL DCOPY(N,chk_Y(1,istack),1,Y0,1)
   Y0(1:N) = chk_Y(1:N,istack)
   
!~~~>   Compute the Jacobian at current time
   CALL LSS_Jac(T,Y0,JAC)
   ISTATUS(Njac) = ISTATUS(Njac) + 1
   
!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T)) 
      CALL LSS_Jac2(T+Delta,Y0,JAC)
      ISTATUS(Njac) = ISTATUS(Njac) +1
      CALL LSS_Jac_Time_Derivative (Delta)
      DO iadj = 1, NADJ
        CALL LSS_Mul_djdttr(dFdT(1,iadj),Y(1,iadj))
        CALL DSCAL(N,(-ONE),dFdT(1,iadj),1)
      END DO
   END IF

!~~~>  Ydot = -J^T*Y
   CALL LSS_Rev_Jac()
   DO iadj = 1, NADJ
     CALL LSS_Mul_Jac(Fcn0(1,iadj),Y(1,iadj))
   END DO
   
!~~~>    Construct Ghimj = 1/(H*ham) - Jac0
     ghinv = ONE/(Direction*H*ros_Gamma(1))
     CALL LSS_Decomp(ghinv,j)
     ISTATUS(Ndec) = ISTATUS(Ndec) + 1
     IF (j /= 0) THEN
       CALL ros_ErrorMsg(-8,T,H,IERR)
       PRINT*,' The matrix is singular !'
       STOP
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S
      
      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)
      
      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         DO iadj = 1, NADJ
           CALL DCOPY(N,Fcn0(1,iadj),1,Fcn(1,iadj),1)
         END DO
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL DCOPY(N*NADJ,Y,1,Ynew,1)
         DO j = 1, istage-1
           DO iadj = 1, NADJ
             CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
                K(N*(j-1)+1,iadj),1,Ynew(1,iadj),1) 
           END DO       
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL ros_Hermite3( chk_T(istack-1), chk_T(istack), Tau, &
             chk_Y(1:N,istack-1), chk_Y(1:N,istack),       &
             chk_dY(1:N,istack-1), chk_dY(1:N,istack), Y0 )
         CALL LSS_Jac1(Tau,Y0,JAC)
         ISTATUS(Njac) = ISTATUS(Njac) + 1
         CALL LSS_Rev_Jac1()
         DO iadj = 1, NADJ
             CALL LSS_Mul_Jactr1(Fcn(1,iadj),Ynew(1,iadj))
         END DO
       END IF ! if istage == 1 elseif ros_NewF(istage)

       DO iadj = 1, NADJ
          CALL DCOPY(N,Fcn(1,iadj),1,K(ioffset+1,iadj),1)
       END DO
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HC,K(N*(j-1)+1,iadj),1, &
                  K(ioffset+1,iadj),1)
         END DO
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HG,dFdT(1,iadj),1,K(ioffset+1,iadj),1)
         END DO
       END IF
  
       Transp = .TRUE.
       DO iadj = 1, NADJ
         CALL LSS_Solve(Transp, K(ioffset+1,iadj))
         ISTATUS(Nsol) = ISTATUS(Nsol) + 1
       END DO
      
   END DO Stage     
            

!~~~>  Compute the new solution 
   DO iadj = 1, NADJ
      DO j=1,ros_S
         CALL DAXPY(N,ros_M(j),K(N*(j-1)+1,iadj),1,Y(1,iadj),1)
      END DO
   END DO

   END DO TimeLoop 
      
!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful
  END SUBROUTINE ros_SimpleCadjInt
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DOUBLE PRECISION FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, & 
               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

! Input arguments   
   DOUBLE PRECISION, INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables     
   DOUBLE PRECISION :: Err, Scale, Ymax   
   INTEGER  :: i
   
   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)
   
  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

!~~~> Input arguments   
   DOUBLE PRECISION, INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N) 
!~~~> Output arguments   
   DOUBLE PRECISION, INTENT(OUT) :: dFdT(N)   
!~~~> Local variables     
   DOUBLE PRECISION :: Delta  
   DOUBLE PRECISION, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
   
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FUN(N,T+Delta,Y,dFdT)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL DAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL DSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix (H,Direction,gam,Singular)
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*gam) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   USE LS_Solver
   IMPLICIT NONE
  
!~~~> Input arguments
!
   DOUBLE PRECISION, INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
   LOGICAL, INTENT(OUT) ::  Singular
!   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   DOUBLE PRECISION, INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: ising, Nconsecutive
   DOUBLE PRECISION :: ghinv
   DOUBLE PRECISION, PARAMETER :: ONE  = 1.0, HALF = 0.5

   Nconsecutive = 0
   Singular = .TRUE.

   
   DO WHILE (Singular)
      ghinv = ONE/(Direction*H*gam)

      CALL LSS_Decomp(ghinv,ISING)
      ISTATUS(Ndec) = ISTATUS(Ndec) + 1

     IF (ising == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ising .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ising

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE ros_cadj_Y( T, Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!  Finds the solution Y at T by interpolating the stored forward trajectory
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
   IMPLICIT NONE
!~~~> Input variables 
   DOUBLE PRECISION, INTENT(IN) :: T
!~~~> Output variables     
   DOUBLE PRECISION, INTENT(OUT) :: Y(N)
!~~~> Local variables     
   INTEGER     :: i
   DOUBLE PRECISION, PARAMETER  :: ONE = 1.0d0

!   chk_H, chk_T, chk_Y, chk_dY, chk_d2Y

   IF( (T < chk_T(1)).OR.(T> chk_T(stack_ptr)) ) THEN
      PRINT*,'Cannot locate solution at T = ',T
      PRINT*,'Stored trajectory is between Tstart = ',chk_T(1)
      PRINT*,'    and Tend = ',chk_T(stack_ptr)
      STOP
   END IF
   DO i = 1, stack_ptr-1
     IF( (T>= chk_T(i)).AND.(T<= chk_T(i+1)) ) EXIT
   END DO 


!   IF (.FALSE.) THEN
!
!   CALL ros_Hermite5( chk_T(i), chk_T(i+1), T, &
!                chk_Y(1,i),   chk_Y(1,i+1),     &
!                chk_dY(1,i),  chk_dY(1,i+1),    &
!                chk_d2Y(1,i), chk_d2Y(1,i+1), Y )
!   
!   ELSE
                
   CALL ros_Hermite3( chk_T(i), chk_T(i+1), T, &
                chk_Y(1:N,i),   chk_Y(1:N,i+1),     &
                chk_dY(1:N,i),  chk_dY(1:N,i+1),    &
                Y )
                        
!   
!   END IF       

  END SUBROUTINE ros_cadj_Y
  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE ros_Hermite3( a, b, T, Ya, Yb, Ja, Jb, Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!  Template for Hermite interpolation of order 5 on the interval [a,b]
! P = c(1) + c(2)*(x-a) + ... + c(4)*(x-a)^3
! P[a,b] = [Ya,Yb], P'[a,b] = [Ja,Jb]
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
   IMPLICIT NONE
!~~~> Input variables 
   DOUBLE PRECISION, INTENT(IN) :: a, b, T, Ya(N), Yb(N)
   DOUBLE PRECISION, INTENT(IN) :: Ja(N), Jb(N)
!~~~> Output variables     
   DOUBLE PRECISION, INTENT(OUT) :: Y(N)
!~~~> Local variables     
   DOUBLE PRECISION :: Tau, amb(3), C(N,4)
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0
   INTEGER :: i, j
   
   amb(1) = 1.0d0/(a-b)
   DO i=2,3
     amb(i) = amb(i-1)*amb(1)
   END DO
   
   
! c(1) = ya;
   CALL DCOPY(N,Ya,1,C(1,1),1)
! c(2) = ja;
   CALL DCOPY(N,Ja,1,C(1,2),1)
! c(3) = 2/(a-b)*ja + 1/(a-b)*jb - 3/(a - b)^2*ya + 3/(a - b)^2*yb  ;
   CALL DCOPY(N,Ya,1,C(1,3),1)
   CALL DSCAL(N,-3.0*amb(2),C(1,3),1)
   CALL DAXPY(N,3.0*amb(2),Yb,1,C(1,3),1)
   CALL DAXPY(N,2.0*amb(1),Ja,1,C(1,3),1)
   CALL DAXPY(N,amb(1),Jb,1,C(1,3),1)
! c(4) =  1/(a-b)^2*ja + 1/(a-b)^2*jb - 2/(a-b)^3*ya + 2/(a-b)^3*yb ;
   CALL DCOPY(N,Ya,1,C(1,4),1)
   CALL DSCAL(N,-2.0*amb(3),C(1,4),1)
   CALL DAXPY(N,2.0*amb(3),Yb,1,C(1,4),1)
   CALL DAXPY(N,amb(2),Ja,1,C(1,4),1)
   CALL DAXPY(N,amb(2),Jb,1,C(1,4),1)
   
   Tau = T - a
   CALL DCOPY(N,C(1,4),1,Y,1)
   CALL DSCAL(N,Tau**3,Y,1)
   DO j = 3,1,-1
     CALL DAXPY(N,TAU**(j-1),C(1,j),1,Y,1)
   END DO       

  END SUBROUTINE ros_Hermite3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE ros_Hermite5( a, b, T, Ya, Yb, Ja, Jb, Ha, Hb, Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!  Template for Hermite interpolation of order 5 on the interval [a,b]
! P = c(1) + c(2)*(x-a) + ... + c(6)*(x-a)^5
! P[a,b] = [Ya,Yb], P'[a,b] = [Ja,Jb], P"[a,b] = [Ha,Hb]
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
   IMPLICIT NONE
!~~~> Input variables 
   DOUBLE PRECISION, INTENT(IN) :: a, b, T, Ya(N), Yb(N)
   DOUBLE PRECISION, INTENT(IN) :: Ja(N), Jb(N), Ha(N), Hb(N)
!~~~> Output variables     
   DOUBLE PRECISION, INTENT(OUT) :: Y(N)
!~~~> Local variables     
   DOUBLE PRECISION :: Tau, amb(5), C(N,6)
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, HALF = 0.5d0
   INTEGER :: i, j
   
   amb(1) = 1.0d0/(a-b)
   DO i=2,5
     amb(i) = amb(i-1)*amb(1)
   END DO
     
! c(1) = ya;
   CALL DCOPY(N,Ya,1,C(1,1),1)
! c(2) = ja;
   CALL DCOPY(N,Ja,1,C(1,2),1)
! c(3) = ha/2;
   CALL DCOPY(N,Ha,1,C(1,3),1)
   CALL DSCAL(N,HALF,C(1,3),1)
   
! c(4) = 10*amb(3)*ya - 10*amb(3)*yb - 6*amb(2)*ja - 4*amb(2)*jb  + 1.5*amb(1)*ha - 0.5*amb(1)*hb ;
   CALL DCOPY(N,Ya,1,C(1,4),1)
   CALL DSCAL(N,10.0*amb(3),C(1,4),1)
   CALL DAXPY(N,-10.0*amb(3),Yb,1,C(1,4),1)
   CALL DAXPY(N,-6.0*amb(2),Ja,1,C(1,4),1)
   CALL DAXPY(N,-4.0*amb(2),Jb,1,C(1,4),1)
   CALL DAXPY(N, 1.5*amb(1),Ha,1,C(1,4),1)
   CALL DAXPY(N,-0.5*amb(1),Hb,1,C(1,4),1)

! c(5) =   15*amb(4)*ya - 15*amb(4)*yb - 8.*amb(3)*ja - 7*amb(3)*jb + 1.5*amb(2)*ha - 1*amb(2)*hb ;
   CALL DCOPY(N,Ya,1,C(1,5),1)
   CALL DSCAL(N, 15.0*amb(4),C(1,5),1)
   CALL DAXPY(N,-15.0*amb(4),Yb,1,C(1,5),1)
   CALL DAXPY(N,-8.0*amb(3),Ja,1,C(1,5),1)
   CALL DAXPY(N,-7.0*amb(3),Jb,1,C(1,5),1)
   CALL DAXPY(N,1.5*amb(2),Ha,1,C(1,5),1)
   CALL DAXPY(N,-amb(2),Hb,1,C(1,5),1)
   
! c(6) =   6*amb(5)*ya - 6*amb(5)*yb - 3.*amb(4)*ja - 3.*amb(4)*jb + 0.5*amb(3)*ha -0.5*amb(3)*hb ;
   CALL DCOPY(N,Ya,1,C(1,6),1)
   CALL DSCAL(N, 6.0*amb(5),C(1,6),1)
   CALL DAXPY(N,-6.0*amb(5),Yb,1,C(1,6),1)
   CALL DAXPY(N,-3.0*amb(4),Ja,1,C(1,6),1)
   CALL DAXPY(N,-3.0*amb(4),Jb,1,C(1,6),1)
   CALL DAXPY(N, 0.5*amb(3),Ha,1,C(1,6),1)
   CALL DAXPY(N,-0.5*amb(3),Hb,1,C(1,6),1)
   
   Tau = T - a
   CALL DCOPY(N,C(1,6),1,Y,1)
   DO j = 5,1,-1
     CALL DSCAL(N,Tau,Y,1)
     CALL DAXPY(N,ONE,C(1,j),1,Y,1)
   END DO       

  END SUBROUTINE ros_Hermite5
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

    IMPLICIT NONE
    DOUBLE PRECISION g
   
    g = 1.0d0 + 1.0d0/SQRT(2.0d0)
   
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'   
!~~~> Number of stages
    ros_S = 2
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
   
    ros_A(1) = (1.d0)/g
    ros_C(1) = (-2.d0)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.d0)/(2.d0*g)
    ros_M(2)= (1.d0)/(2.d0*g)
! E_i = Coefficients for error estimator    
    ros_E(1) = 1.d0/(2.d0*g)
    ros_E(2) = 1.d0/(2.d0*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0d0
    ros_Alpha(2) = 1.0d0 
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = g
    ros_Gamma(2) =-g
   
!    ros_B(1) = 0.5d0
!    ros_B(2) = 0.5d0
   
 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   IMPLICIT NONE
    
    rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'   
!~~~> Number of stages
   ros_S = 3
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
   ros_A(1)= 1.d0
   ros_A(2)= 1.d0
   ros_A(3)= 0.d0

   ros_C(1) = -0.10156171083877702091975600115545d+01
   ros_C(2) =  0.40759956452537699824805835358067d+01
   ros_C(3) =  0.92076794298330791242156818474003d+01
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1d+01
   ros_M(2) =  0.61697947043828245592553615689730d+01
   ros_M(3) = -0.42772256543218573326238373806514d+00
! E_i = Coefficients for error estimator    
   ros_E(1) =  0.5d+00
   ros_E(2) = -0.29079558716805469821718236208017d+01
   ros_E(3) =  0.22354069897811569627360909276199d+00
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0d+00
   ros_Alpha(2)= 0.43586652150845899941601945119356d+00
   ros_Alpha(3)= 0.43586652150845899941601945119356d+00
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1)= 0.43586652150845899941601945119356d+00
   ros_Gamma(2)= 0.24291996454816804366592249683314d+00
   ros_Gamma(3)= 0.21851380027664058511513169485832d+01

!   ros_B(1) = -0.75457412385404315829818998646589d+00
!   ros_B(2) = 1.94100407061964420292840123379419d+00
!   ros_B(3) = -0.18642994676560104463021124732829d+00
  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3 
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

   IMPLICIT NONE
   
    rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'   
!~~~> Number of stages
   ros_S = 4
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
   ros_A(1) = 0.2000000000000000d+01
   ros_A(2) = 0.1867943637803922d+01
   ros_A(3) = 0.2344449711399156d+00
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0D0

   ros_C(1) =-0.7137615036412310d+01
   ros_C(2) = 0.2580708087951457d+01
   ros_C(3) = 0.6515950076447975d+00
   ros_C(4) =-0.2137148994382534d+01
   ros_C(5) =-0.3214669691237626d+00
   ros_C(6) =-0.6949742501781779d+00
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735d+01
   ros_M(2) = 0.2870493262186792d+00
   ros_M(3) = 0.4353179431840180d+00
   ros_M(4) = 0.1093502252409163d+01
!~~~> E_i  = Coefficients for error estimator    
   ros_E(1) =-0.2815431932141155d+00
   ros_E(2) =-0.7276199124938920d-01
   ros_E(3) =-0.1082196201495311d+00
   ros_E(4) =-0.1093502252409163d+01
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.D0
   ros_Alpha(2) = 0.1145640000000000d+01
   ros_Alpha(3) = 0.6552168638155900d+00
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1) = 0.5728200000000000d+00
   ros_Gamma(2) =-0.1769193891319233d+01
   ros_Gamma(3) = 0.7592633437920482d+00
   ros_Gamma(4) =-0.1049021087100450d+00

!   ros_B(1) = 0.2255570073418735d+01
!   ros_B(2) = 0.2870493262186792d+00
!   ros_B(3) = 0.4353179431840180d+00
!   ros_B(4) = 0.1093502252409163d+01

  END SUBROUTINE Ros4
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

   IMPLICIT NONE
   
    rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'   
!~~~> Number of stages
   ros_S = 4
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
 
   ros_A(1) = 0.0d+00
   ros_A(2) = 2.0d+00
   ros_A(3) = 0.0d+00
   ros_A(4) = 2.0d+00
   ros_A(5) = 0.0d+00
   ros_A(6) = 1.0d+00

   ros_C(1) = 4.0d+00
   ros_C(2) = 1.0d+00
   ros_C(3) =-1.0d+00
   ros_C(4) = 1.0d+00
   ros_C(5) =-1.0d+00 
   ros_C(6) =-(8.0d+00/3.0d+00) 
         
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0d+00
   ros_M(2) = 0.0d+00
   ros_M(3) = 1.0d+00
   ros_M(4) = 1.0d+00
!~~~> E_i  = Coefficients for error estimator    
   ros_E(1) = 0.0d+00
   ros_E(2) = 0.0d+00
   ros_E(3) = 0.0d+00
   ros_E(4) = 1.0d+00
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0d+00    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0d+00
   ros_Alpha(2) = 0.0d+00
   ros_Alpha(3) = 1.0d+00
   ros_Alpha(4) = 1.0d+00
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1) = 0.5d+00
   ros_Gamma(2) = 1.5d+00
   ros_Gamma(3) = 0.0d+00
   ros_Gamma(4) = 0.0d+00

!   ros_B(1) = 5.0d0/6.0d0
!   ros_B(2) = -1.0d0/6.0d0
!   ros_B(3) = -1.0d0/6.0d0
!   ros_B(4) = 0.5d0

  END SUBROUTINE Rodas3
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'   
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000d0
    ros_Alpha(2) = 0.386d0
    ros_Alpha(3) = 0.210d0 
    ros_Alpha(4) = 0.630d0
    ros_Alpha(5) = 1.000d0
    ros_Alpha(6) = 1.000d0
        
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = 0.2500000000000000d+00
    ros_Gamma(2) =-0.1043000000000000d+00
    ros_Gamma(3) = 0.1035000000000000d+00
    ros_Gamma(4) =-0.3620000000000023d-01
    ros_Gamma(5) = 0.0d0
    ros_Gamma(6) = 0.0d0

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
    ros_A(1) = 0.1544000000000000d+01
    ros_A(2) = 0.9466785280815826d+00
    ros_A(3) = 0.2557011698983284d+00
    ros_A(4) = 0.3314825187068521d+01
    ros_A(5) = 0.2896124015972201d+01
    ros_A(6) = 0.9986419139977817d+00
    ros_A(7) = 0.1221224509226641d+01
    ros_A(8) = 0.6019134481288629d+01
    ros_A(9) = 0.1253708332932087d+02
    ros_A(10) =-0.6878860361058950d+00
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0d+00

    ros_C(1) =-0.5668800000000000d+01
    ros_C(2) =-0.2430093356833875d+01
    ros_C(3) =-0.2063599157091915d+00
    ros_C(4) =-0.1073529058151375d+00
    ros_C(5) =-0.9594562251023355d+01
    ros_C(6) =-0.2047028614809616d+02
    ros_C(7) = 0.7496443313967647d+01
    ros_C(8) =-0.1024680431464352d+02
    ros_C(9) =-0.3399990352819905d+02
    ros_C(10) = 0.1170890893206160d+02
    ros_C(11) = 0.8083246795921522d+01
    ros_C(12) =-0.7981132988064893d+01
    ros_C(13) =-0.3152159432874371d+02
    ros_C(14) = 0.1631930543123136d+02
    ros_C(15) =-0.6058818238834054d+01

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0d+00
    ros_M(6) = 1.0d+00

!~~~> E_i  = Coefficients for error estimator    
    ros_E(1) = 0.0d+00
    ros_E(2) = 0.0d+00
    ros_E(3) = 0.0d+00
    ros_E(4) = 0.0d+00
    ros_E(5) = 0.0d+00
    ros_E(6) = 1.0d+00

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.
     
!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0d0

  END SUBROUTINE Rodas4


END SUBROUTINE RosenbrockADJ1 ! and its internal procedures


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE RosenbrockADJ2( N, NP, NNZERO, Y, NADJ, Lambda, Tstart, Tend, FUN, &
             JAC, ADJINIT, HESSTR_VEC, AbsTol, RelTol, AbsTol_adj, RelTol_adj,&
             RCNTRL, ICNTRL, RSTATUS, ISTATUS, IERR, Q, QFUN, DRDY,           &
             HESSTR_VEC_R)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   
!    ADJ = Adjoint of the Tangent Linear Model of a Rosenbrock Method
!
!    Solves the system y'=F(t,y) using a RosenbrockADJ method defined by:
!
!     G = 1/(H*gamma(1)) - Jac(t0,Y0)
!     T_i = t0 + Alpha(i)*H
!     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
!     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
!         gamma(i)*dF/dT(t0, Y0)
!     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j 
!
!    For details on RosenbrockADJ methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.  
!    The codes contained in the book inspired this implementation.       
!
!    (C)  Adrian Sandu, August 2004
!    Virginia Polytechnic Institute and State University    
!    Contact: sandu@cs.vt.edu
!    Revised by Hong Zhang and Adrian Sandu, Feb 2011          
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    
!~~~>   INPUT ARGUMENTS: 
!    
!-     Y(N)    = vector of initial conditions (at T=Tstart)
!      NADJ       -> dimension of linearized system, 
!                   i.e. the number of sensitivity coefficients
!-     Lambda(N,NADJ) -> vector of initial sensitivity conditions (at T=Tstart)
!-    [Tstart,Tend]  = time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)  
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE Fun( T, Y, Ydot ) = ODE function, 
!                       returns Ydot = Y' = F(T,Y) 
!- SUBROUTINE Jac( T, Y, Jcb ) = Jacobian of the ODE function,
!                       returns Jcb = dF/dY 
!-    ICNTRL(1:10)    = integer inputs parameters
!-    RCNTRL(1:10)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:  
!     
!-    Y(N)    -> vector of final states (at T->Tend)
!-    Lambda(N,NADJ) -> vector of final sensitivities (at T=Tend)
!-    ICNTRL(11:20)   -> integer output parameters
!-    RCNTRL(11:20)   -> real output parameters
!-    IERR       -> job status upon return
!       - succes (positive value) or failure (negative value) -
!           =  1 : Success
!           = -1 : Improper value for maximal no of steps
!           = -2 : Selected RosenbrockADJ method not implemented
!           = -3 : Hmin/Hmax/Hstart must be positive
!           = -4 : FacMin/FacMax/FacRej must be positive
!           = -5 : Improper tolerance values
!           = -6 : No of steps exceeds maximum bound
!           = -7 : Step size too small
!           = -8 : Matrix is repeatedly singular
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1)   = 1: F = F(y)   Independent of T (AUTONOMOUS)
!              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)
!
!    ICNTRL(2)   = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1:  AbsTol, RelTol are scalars
!
!    ICNTRL(3)  -> selection of a particular Rosenbrock method
!        = 0 :  default method is Rodas3
!        = 1 :  method is  Ros2
!        = 2 :  method is  Ros3 
!        = 3 :  method is  Ros4 
!        = 4 :  method is  Rodas3
!        = 5:   method is  Rodas4
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0) the default value of BUFSIZE is used
!
!    ICNTRL(6)  -> selection of a particular Rosenbrock method for the
!                continuous adjoint integration - for cts adjoint it
!                can be different than the forward method ICNTRL(3)
!         Note 1: to avoid interpolation errors (which can be huge!) 
!                   it is recommended to use only ICNTRL(7) = 2 or 4
!         Note 2: the performance of the full continuous adjoint
!                   strongly depends on the forward solution accuracy Abs/RelTol
!
!    ICNTRL(7) -> Type of adjoint algorithm
!         = 0 : default is discrete adjoint ( of method ICNTRL(3) )
!         = 1 : no adjoint       
!         = 2 : discrete adjoint ( of method ICNTRL(3) )
!         = 3 : fully adaptive continuous adjoint ( with method ICNTRL(6) )
!         = 4 : simplified continuous adjoint ( with method ICNTRL(6) )
!
!    ICNTRL(8)  -> checkpointing the LU factorization at each step:
!        ICNTRL(8)=0 : do *not* save LU factorization (the default)
!        ICNTRL(8)=1 : save LU factorization
!        Note: if ICNTRL(7)=1 the LU factorization is *not* saved
!
!~~~>  Real input parameters:
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO 
!
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!          
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!            (default=0.1)
!
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller 
!         than the predicted value  (default=0.9)
!
!    RCNTRL(8)  -> ThetaMin. If Newton convergence rate smaller
!                  than ThetaMin the Jacobian is not recomputed;
!                  (default=0.001)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call to RosenbrockADJ adds the corrent no. of fcn calls
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
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    For multiple restarts, use Hexit as Hstart in the following run 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   USE LS_Solver
   IMPLICIT NONE
   
!~~~>  Arguments
   INTEGER                         :: N, NNZERO   
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
   INTEGER, INTENT(IN)             :: NP,NADJ
   DOUBLE PRECISION, INTENT(IN),OPTIONAL    :: Q(NADJ)
   DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ)
   DOUBLE PRECISION, INTENT(IN)    :: Tstart,Tend
   DOUBLE PRECISION, INTENT(IN)    :: AbsTol(N),RelTol(N)
   DOUBLE PRECISION, INTENT(IN)    :: AbsTol_adj(N,NADJ), RelTol_adj(N,NADJ)
   INTEGER, INTENT(IN)             :: ICNTRL(20)
   DOUBLE PRECISION, INTENT(IN)    :: RCNTRL(20)
   INTEGER, INTENT(INOUT)          :: ISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)            :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5
   DOUBLE PRECISION :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Types of Adjoints Implemented
   INTEGER, PARAMETER :: Adj_none = 1, Adj_discrete = 2,      &
                   Adj_continuous = 3, Adj_simple_continuous = 4
!~~~>  Checkpoints in memory
   INTEGER, PARAMETER :: bufsize = 200000
   INTEGER :: stack_ptr = 0 ! last written entry
   DOUBLE PRECISION, DIMENSION(:),   POINTER :: chk_H, chk_T
   DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_Y, chk_K
 
!   DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_J, chk_P
   DOUBLE PRECISION, DIMENSION(:,:), POINTER :: chk_dY, chk_d2Y
!~~~>  Local variables     
   DOUBLE PRECISION :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   DOUBLE PRECISION :: Hmin, Hmax, Hstart
   DOUBLE PRECISION :: Texit
   INTEGER :: i, UplimTol, Max_no_steps
   INTEGER :: AdjointType, CadjMethod 
   LOGICAL :: Autonomous, QuadAutonomous, VectorTol, SaveLU, GetQuad
!~~~>   Parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
   DOUBLE PRECISION :: DLAMCH
   EXTERNAL QFUN,FUN,JAC,ADJINIT,HESSTR_VEC,HESSTR_VEC_R,DRDY
   OPTIONAL QFUN,HESSTR_VEC_R,DRDY
 
   ros_A = ZERO
   ros_C = ZERO
   ros_M = ZERO
   ros_E = ZERO
   ros_Alpha = ZERO
   ros_Gamma = ZERO
!~~~>  Initialize statistics
   GetQuad = .FALSE.
!~~~> Compute the quadrature term if required sources are provided
   IF( PRESENT(Q) .AND. PRESENT(QFUN) .AND. PRESENT(DRDY) .AND. &
       PRESENT(HESSTR_VEC_R) ) THEN
        GetQuad = .TRUE.
   END IF

   ISTATUS(1:20) = 0
   RSTATUS(1:20) = ZERO
   
!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE 
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (ICNTRL(3))
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT
   
!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = bufsize - 1
   ELSEIF (Max_no_steps > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE 
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN      
   END IF

!~~~>  The particular Rosenbrock method chosen for integrating the cts adjoint
   IF (ICNTRL(6) == 0) THEN
!      CadjMethod = 4
      CadjMethod = 3
   ELSEIF ( (ICNTRL(6) >= 1).AND.(ICNTRL(6) <= 5) ) THEN
      CadjMethod = ICNTRL(6)
   ELSE  
      PRINT * , 'Unknown CADJ Rosenbrock method: ICNTRL(6)=', CadjMethod
      CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN      
   END IF
 
!~~~>  Discrete or continuous adjoint formulation
   IF ( ICNTRL(7) == 0 ) THEN
       AdjointType = Adj_discrete
   ELSEIF ( (ICNTRL(7) >= 1).AND.(ICNTRL(7) <= 4) ) THEN
       AdjointType = ICNTRL(7)
   ELSE  
      PRINT * , 'User-selected adjoint type: ICNTRL(7)=', AdjointType
      CALL ros_ErrorMsg(-9,Tstart,ZERO,IERR)
      RETURN      
   END IF

!~~~> Save or not the forward LU factorization
      SaveLU = (ICNTRL(8) /= 0) 

!~~~> The quadrature term is autonomous or not 
    QuadAutonomous = .NOT.(ICNTRL(9) == 0)
 
!~~~>  Unit roundoff (1+Roundoff>1)  
   Roundoff = DLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN 
      Hmin = RCNTRL(1)
   ELSE  
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE  
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE  
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax 
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2d0
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE  
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0d0
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE  
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1d0
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE  
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9d0
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE  
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN      
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.d0*Roundoff) &
         .OR. (RelTol(i) >= 1.0d0) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
     
!~~~>  Allocate checkpoint space or open checkpoint files
   IF (AdjointType == Adj_discrete) THEN
       CALL ros_AllocateDBuffers( ros_S )
   ELSEIF ( (AdjointType == Adj_continuous).OR.  &
           (AdjointType == Adj_simple_continuous) ) THEN
       CALL ros_AllocateCBuffers
   END IF
   CALL LSS_Init(N,NNZERO)
!~~~>  CALL Forward Rosenbrock method   
   CALL ros_FwdInt(N,NADJ,Y,Tstart,Tend,Texit,AbsTol,RelTol,GetQuad,          & 
!  Error indicator
        IERR)

   PRINT*,'FORWARD STATISTICS'
   PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

!~~~>  If Forward integration failed return   
   IF (IERR<0) RETURN
!~~~> Initialize adjoint variables
   CALL ADJINIT(N,NP,NADJ,Tend,Y,Lambda)
!~~~>   Initialize the particular Rosenbrock method for continuous adjoint
   IF ( (AdjointType == Adj_continuous).OR. &
           (AdjointType == Adj_simple_continuous) ) THEN
   SELECT CASE (CadjMethod)
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=', ICNTRL(3)
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR) 
       RETURN     
   END SELECT
   END IF

   SELECT CASE (AdjointType)   
   CASE (Adj_discrete)   
     CALL ros_DadjInt (                          &
        N, NADJ, Lambda,                         &
        Tstart, Tend, Texit, GetQuad,            &
        IERR )
   CASE (Adj_continuous) 
     CALL ros_CadjInt (                          &
        NADJ, Lambda,                            &
        Tend, Tstart, Texit,                     &
        AbsTol_adj, RelTol_adj,                  &
        IERR )
   CASE (Adj_simple_continuous)
     CALL ros_SimpleCadjInt (                    &
        NADJ, Lambda,                            &
        Tstart, Tend, Texit,                     &
        IERR )
   END SELECT ! AdjointType

   PRINT*,'ADJOINT STATISTICS'
   PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)

   CALL LSS_Free
!~~~>  Free checkpoint space or close checkpoint files
   IF (AdjointType == Adj_discrete) THEN
      CALL ros_FreeDBuffers
   ELSEIF ( (AdjointType == Adj_continuous) .OR. &
           (AdjointType == Adj_simple_continuous) ) THEN
      CALL ros_FreeCBuffers
   END IF
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  Procedures internal to RosenbrockADJ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_AllocateDBuffers( S )
!~~~>  Allocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INTEGER :: i, S
   
   ALLOCATE( chk_H(bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer H'; STOP
   END IF   
   ALLOCATE( chk_T(bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer T'; STOP
   END IF   
   ALLOCATE( chk_Y(N*S,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer Y'; STOP
   END IF   
   ALLOCATE( chk_K(N*S,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer K'; STOP
   END IF  

 END SUBROUTINE ros_AllocateDBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_FreeDBuffers
!~~~>  Dallocate buffer space for discrete adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INTEGER :: i
   
   DEALLOCATE( chk_H, STAT=i )
   IF (i/=0) THEN
      PRINT*, 'Failed deallocation of buffer H'; STOP
   END IF   
   DEALLOCATE( chk_T, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer T'; STOP
   END IF   
   DEALLOCATE( chk_Y, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer Y'; STOP
   END IF   
   DEALLOCATE( chk_K, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer K'; STOP
   END IF   
!  SaveLU disabled
!   IF (SaveLU) THEN 
!     DEALLOCATE( chk_J, STAT=i )
!     IF (i/=0) THEN
!        PRINT*,'Failed deallocation of buffer J'; STOP
!     END IF
!     DEALLOCATE( chk_P, STAT=i )
!     IF (i/=0) THEN
!        PRINT*,'Fainled deallocation of buffer P'; STOP
!     END IF   
!   END IF   
 
 END SUBROUTINE ros_FreeDBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_AllocateCBuffers
!~~~>  Allocate buffer space for continuous adjoint
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INTEGER :: i
   
   ALLOCATE( chk_H(bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer H'; STOP
   END IF   
   ALLOCATE( chk_T(bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer T'; STOP
   END IF   
   ALLOCATE( chk_Y(N,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer Y'; STOP
   END IF   
   ALLOCATE( chk_dY(N,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer dY'; STOP
   END IF   
   ALLOCATE( chk_d2Y(N,bufsize), STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed allocation of buffer d2Y'; STOP
   END IF   
 
 END SUBROUTINE ros_AllocateCBuffers


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_FreeCBuffers
!~~~>  Dallocate buffer space for continuous adjoint
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
   DEALLOCATE( chk_dY, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer dY'; STOP
   END IF   
   DEALLOCATE( chk_d2Y, STAT=i )
   IF (i/=0) THEN
      PRINT*,'Failed deallocation of buffer d2Y'; STOP
   END IF   
 
 END SUBROUTINE ros_FreeCBuffers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_DPush( S, T, H, Ystage, K, P )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INTEGER :: S ! no of stages
   DOUBLE PRECISION :: T, H, Ystage(N*S), K(N*S)
   INTEGER       :: P(N)
   
   stack_ptr = stack_ptr + 1
   IF ( stack_ptr > bufsize ) THEN
     PRINT*,'Push failed: buffer overflow'
     STOP
   END IF  
   chk_H( stack_ptr ) = H
   chk_T( stack_ptr ) = T
   !CALL DCOPY(N*S,Ystage,1,chk_Y(1,stack_ptr),1)
   !CALL DCOPY(N*S,K,1,chk_K(1,stack_ptr),1)
   chk_Y(1:N*S,stack_ptr) = Ystage(1:N*S)
   chk_K(1:N*S,stack_ptr) = K(1:N*S)
  
  END SUBROUTINE ros_DPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_DPop( S, T, H, Ystage, K )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Retrieves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   INTEGER :: S ! no of stages
   DOUBLE PRECISION :: T, H, Ystage(N*S), K(N*S)
!   INTEGER       :: P(N)
   
   IF ( stack_ptr <= 0 ) THEN
     PRINT*,'Pop failed: empty buffer'
     STOP
   END IF  
   H = chk_H( stack_ptr )
   T = chk_T( stack_ptr )
   !CALL DCOPY(N*S,chk_Y(1,stack_ptr),1,Ystage,1)
   !CALL DCOPY(N*S,chk_K(1,stack_ptr),1,K,1)
   Ystage(1:N*S) = chk_Y(1:N*S,stack_ptr)
   K(1:N*S)      = chk_K(1:N*S,stack_ptr)
   !CALL DCOPY(LU_NONZERO,chk_J(1,stack_ptr),1,Jcb,1)

   stack_ptr = stack_ptr - 1
  
  END SUBROUTINE ros_DPop
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_CPush( T, H, Y, dY, d2Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Saves the next trajectory snapshot for discrete adjoints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   DOUBLE PRECISION :: T, H, Y(N), dY(N), d2Y(N)
   
   stack_ptr = stack_ptr + 1
   IF ( stack_ptr > bufsize ) THEN
     PRINT*,'Push failed: buffer overflow'
     STOP
   END IF  
   chk_H( stack_ptr ) = H
   chk_T( stack_ptr ) = T
   !CALL DCOPY(N,Y,1,chk_Y(1,stack_ptr),1)
   !CALL DCOPY(N,dY,1,chk_dY(1,stack_ptr),1)
   !CALL DCOPY(N,d2Y,1,chk_d2Y(1,stack_ptr),1)
   chk_Y(1:N,stack_ptr)   = Y(1:N)
   chk_dY(1:N,stack_ptr)  = dY(1:N)
   chk_d2Y(1:N,stack_ptr) = d2Y(1:N)
  END SUBROUTINE ros_CPush
  
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_CPop( T, H, Y, dY, d2Y )
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
   !CALL DCOPY(N,chk_Y(1,stack_ptr),1,Y,1)
   !CALL DCOPY(N,chk_dY(1,stack_ptr),1,dY,1)
   !CALL DCOPY(N,chk_d2Y(1,stack_ptr),1,d2Y,1)
   Y(1:N)   = chk_Y(1:N,stack_ptr)
   dY(1:N)  = chk_dY(1:N,stack_ptr)
   d2Y(1:N) = chk_d2Y(1:N,stack_ptr)

   stack_ptr = stack_ptr - 1
  
  END SUBROUTINE ros_CPop



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   DOUBLE PRECISION, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from RosenbrockADJ due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)    
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected RosenbrockADJ method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum buffer bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE (-9)    
      PRINT * , '--> Improper type of adjoint selected'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
   PRINT *, "T=", T, "and H=", H
     
 END SUBROUTINE ros_ErrorMsg
   
     
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_FwdInt (N, NADJ, Y,  Tstart, Tend, T, AbsTol, RelTol, GetQuad,&
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockADJ method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   USE LS_Solver
   IMPLICIT NONE
   INTEGER, INTENT(IN)              :: N, NADJ
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
!~~~> Input: integration interval   
   DOUBLE PRECISION, INTENT(IN)    :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   DOUBLE PRECISION, INTENT(OUT)   ::  T      
!~~~> Input: tolerances      
   DOUBLE PRECISION, INTENT(IN)    ::  AbsTol(N), RelTol(N)
!~~~> Input: Quadrature term indicator
   LOGICAL, INTENT(IN)             :: GetQuad
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT)            :: IERR
! ~~~~ Local variables        
   DOUBLE PRECISION :: Ynew(N), Fcn0(N), Fcn(N) 
   DOUBLE PRECISION :: K(N*ros_S), dFdT(N)
   DOUBLE PRECISION, DIMENSION(:), POINTER :: Ystage
   DOUBLE PRECISION :: H, Hnew, HC, HG, Fac, Tau 
   DOUBLE PRECISION :: Err, Yerr(N)
   INTEGER :: Pivot(N), Direction, ioffset, i, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular, Transp
!~~~~~ For Quadrature term
   DOUBLE PRECISION, ALLOCATABLE :: L(:),dRdT(:),R(:),RY(:,:)
   DOUBLE PRECISION, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
   DOUBLE PRECISION :: Delta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IF(GetQuad) THEN
     ALLOCATE(L(NADJ*ros_S),dRdT(NADJ),R(NADJ),RY(N,NADJ),STAT=i)
     IF(i .NE. 0) STOP 'allocation error for L,dRdT,R,RY'
   END IF

!~~~>  Allocate stage vector buffer if needed
   IF (AdjointType == Adj_discrete) THEN
        ALLOCATE(Ystage(N*ros_S), STAT=i)
        ! Uninitialized Ystage may lead to NaN on some compilers
        Ystage = 0.0d0
        IF (i/=0) THEN
          PRINT*,'Allocation of Ystage failed'
          STOP
        END IF
   END IF   
   
!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.
   
!~~~> Time loop begins below 

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) ) 
      
   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1d0*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   
!~~~>  Limit H if necessary to avoid going beyond Tend   
   RSTATUS(Nhexit) = H
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL FUN(N,T,Y,Fcn0)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1 

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF
  
!~~~>   Compute the Jacobian at current time
   CALL LSS_Jac(T,Y,JAC)
   ISTATUS(Njac) = ISTATUS(Njac) + 1

!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO  
   
   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1),Singular)

   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S
      
      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)
      
      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         CALL DCOPY(N,Fcn0,1,Fcn,1)
         IF (AdjointType == Adj_discrete) THEN ! Save stage solution
            ! CALL DCOPY(N,Y,1,Ystage(1),1)
            Ystage(1:N) = Y(1:N)
            CALL DCOPY(N,Y,1,Ynew,1)
         END IF   
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL DCOPY(N,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
            K(N*(j-1)+1),1,Ynew,1) 
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FUN(N,Tau,Ynew,Fcn)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
      ! save stage solution every time even if ynew is not updated
       IF ( ( istage > 1 ).AND.(AdjointType == Adj_discrete) ) THEN
         ! CALL DCOPY(N,Ynew,1,Ystage(ioffset+1),1)
         Ystage(ioffset+1:ioffset+N) = Ynew(1:N)
       END IF   
       CALL DCOPY(N,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL DAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL DAXPY(N,HG,dFdT,1,K(ioffset+1),1)
       END IF
       Transp = .FALSE.
       CALL LSS_Solve(Transp, K(ioffset+1))
       ISTATUS(Nsol) = ISTATUS(Nsol) + 1
      
   END DO Stage     
            

!~~~>  Compute the new solution 
   CALL DCOPY(N,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL DAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation 
   CALL DSCAL(N,ZERO,Yerr,1)
   DO j=1,ros_S     
        CALL DAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO 
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac  

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      IF (AdjointType == Adj_discrete) THEN ! Save current state
          CALL ros_DPush( ros_S, T, H, Ystage, K,  Pivot )
      ELSEIF ( (AdjointType == Adj_continuous) .OR. &
           (AdjointType == Adj_simple_continuous) ) THEN
          CALL LSS_Mul_Jac(K,Fcn0)
          IF (.NOT. Autonomous) THEN
             CALL DAXPY(N,ONE,dFdT,1,K(1),1)
          END IF   
          CALL ros_CPush( T, H, Y, Fcn0, K(1) )
      END IF     

!~~~> Update the results for the quadrature term
      IF(GetQuad) Then
        CALL DRDY(NADJ,N,N,T,Y,RY)
        CALL QFUN(N,NADJ,T,Y,R)
        IF(.NOT.Autonomous) THEN
          Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
          CALL QFUN(N,NADJ,T+Delta,Y,dRdT)
          CALL DAXPY(NADJ,(-ONE),R,1,dRdT,1)
          CALL DSCAL(NADJ,(ONE/Delta),dRdT,1)
        END IF

!~~~>   Compute the stages for quadrature term
        DO istage = 1, ros_S
! Current istage offset. Current istage vector is L(ioffset+1:ioffset+NADJ)
          ioffset = NADJ*(istage-1)
! For the 1st istage the function has been computed previously
          IF ( istage == 1 ) THEN

! istage>1 and a new function evaluation is needed at the current istage
          ELSEIF ( ros_NewF(istage) ) THEN
            Tau = T + ros_Alpha(istage)*Direction*H
            CALL QFUN(N, NADJ,Tau,Ystage(N*(istage-1)+1),R)
          END IF ! if istage == 1 elseif ros_NewF(istage)

          CALL DCOPY(NADJ,R,1,L(ioffset+1),1)
          DO j = 1, istage-1
            HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
            CALL DAXPY(NADJ,HC,L(NADJ*(j-1)+1),1,L(ioffset+1),1)
          END DO

          IF (.NOT. Autonomous .AND. ros_Gamma(istage).NE.ZERO) THEN
            HG = Direction*H*ros_Gamma(istage)
            CALL DAXPY(NADJ,HG,dRdT,1,L(ioffset+1),1)
          END IF
! +r_y*k_i           
          DO i=1,NADJ
            DO j=1,N
              L(ioffset+i) = L(ioffset+i)+RY(j,i)*K(N*(istage-1)+j)
            END DO
          END DO
          CALL DSCAL(NADJ,Direction*H*ros_Gamma(1),L(ioffset+1),1)
        END DO 
!~~~>  Compute the new solution 
        DO j=1,ros_S
          CALL DAXPY(NADJ,ros_M(j),L(NADJ*(j-1)+1),1,Q,1)
        END DO
      END IF

      CALL DCOPY(N,Ynew,1,Y,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H) 
      END IF   
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.  
      RejectMoreH = .FALSE.
      H = Hnew      
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF   
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1) THEN
         ISTATUS(Nrej) = ISTATUS(Nrej) + 1
      END IF    
   END IF ! Err <= 1

   END DO UntilAccepted 

   END DO TimeLoop 
   
!~~~> Save last state: only needed for continuous adjoint
   IF ( (AdjointType == Adj_continuous) .OR. &
       (AdjointType == Adj_simple_continuous) ) THEN
       CALL FUN(N,T,Y,Fcn0)
       ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       CALL LSS_Jac(T,Y,JAC)
       ISTATUS(Njac) = ISTATUS(Njac) + 1

       CALL LSS_Mul_Jac(K, Fcn0)
       IF (.NOT. Autonomous) THEN
           CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
           CALL DAXPY(N,ONE,dFdT,1,K(1),1)
       END IF   
       CALL ros_CPush( T, H, Y, Fcn0, K(1) )
!~~~> Deallocate stage buffer: only needed for discrete adjoint
   ELSEIF (AdjointType == Adj_discrete) THEN 
        DEALLOCATE(Ystage, STAT=i)
        IF (i/=0) THEN
          PRINT*,'Deallocation of Ystage failed'
          STOP
        END IF
   END IF   

   IF(GetQuad) THEN
     DEALLOCATE(L,dRdT,R,RY,STAT=i)
     IF(i .NE. 0) STOP 'deallocation error for L, dRdT, R, RY'
   END IF
   
!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_FwdInt
   
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBROUTINE ros_DadjInt (                        &
           N, NADJ, Lambda,                         &
           Tstart, Tend, T, GetQuad,                &
!~~~> Error indicator
           IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockSOA method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE LS_Solver
      IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
      INTEGER, INTENT(IN)             :: N, NADJ
!~~~> First order adjoint   
      DOUBLE PRECISION, INTENT(INOUT) :: Lambda(N,NADJ)
!!~~~> Input: integration interval   
      DOUBLE PRECISION, INTENT(IN)    :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
      DOUBLE PRECISION, INTENT(OUT)   ::  T      
!~~~> INPUT: Quadrature term indicator
      LOGICAL :: GetQuad
!~~~> Output: Error indicator
      INTEGER, INTENT(OUT)            :: IERR
! ~~~~ Local variables        
      DOUBLE PRECISION :: Ystage(N*ros_S), K(N*ros_S), ros_W(ros_S)
      DOUBLE PRECISION, ALLOCATABLE :: U(:,:),V(:,:),RY0(:,:),RY1(:,:)
      DOUBLE PRECISION :: Tmp(N), Tmp2(N)
      DOUBLE PRECISION :: H, HC, HA, Tau, Alpha, GammaUhat 
      INTEGER :: Direction
      INTEGER :: j, m, istage, istart, jstart
!~~~>  Local parameters
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
      DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
      LOGICAL :: Transp
      DOUBLE PRECISION :: Delta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ALLOCATE(U(N*ros_S,NADJ), V(N*ros_S,NADJ),STAT=j)
      IF(j .NE. 0) STOP 'allocation error for U,V'

      IF(GetQuad) THEN
        ALLOCATE(RY0(N,NADJ),RY1(N,NADJ),STAT=j)
        IF(j .NE. 0) STOP 'allocation error for RYs'
        RY0(:,:) = ZERO
        RY1(:,:) = ZERO
      END IF

      Transp = .TRUE.
      IF (Tend  >=  Tstart) THEN
        Direction = +1
      ELSE
        Direction = -1
      END IF               

      !~~~> calculate ros_W first
      IF(GetQuad) THEN
        GammaUhat =ZERO
        DO istage=1,ros_S
          ros_W(istage) = H*ros_Gamma(1)*ros_M(istage)
          DO j = istage+1,ros_S
            HC = ros_C((j-1)*(j-2)/2+istage)
            ros_W(istage) = ros_W(istage)+ HC*ros_Gamma(1)*ros_W(j)
          END DO
          GammaUhat = GammaUhat + ros_W(istage)*ros_Gamma(istage)
        END DO
      END IF
!~~~> Time loop begins below 
TimeLoop: DO WHILE ( stack_ptr > 0 )
        
   !~~~>  Recover checkpoints for stage values and vectors
        CALL ros_DPop( ros_S, T, H, Ystage, K )

        ISTATUS(Nstp) = ISTATUS(Nstp) + 1

!~~~>    Compute LU decomposition 
        IF (.NOT.SaveLU) THEN
          CALL LSS_Jac(T,Ystage(1),JAC)
          ISTATUS(Njac) = ISTATUS(Njac) + 1
          Tau = ONE/(Direction*H*ros_Gamma(1))
          ISTATUS(Njac) = ISTATUS(Njac) + 1
     
          CALL LSS_Decomp(Tau,j)
          ISTATUS(Ndec) = ISTATUS(Ndec) + 1
        END IF
            
        IF(GetQuad) CALL DRDY(NADJ,N,N,T,Ystage(1),RY0)
!~~~>   Compute Hessian at the beginning of the interval
!   CALL HESS(N, T,Ystage(1),Hes0)
   
!~~~>   Compute the stages
Stage: DO istage = ros_S, 1, -1
      
!~~~> Current istage first entry 
          istart = N*(istage-1) + 1
!~~~> Compute U
          DO m = 1,NADJ
            CALL DCOPY(N,Lambda(1,m),1,U(istart,m),1)
            CALL DSCAL(N,ros_M(istage),U(istart,m),1)
            IF(GetQuad) CALL DAXPY(N,ros_W(istage),RY0,1,U(istart,m),1)
          END DO ! m=1:NADJ

          DO j = istage+1, ros_S
            jstart = N*(j-1) + 1
            HA = ros_A((j-1)*(j-2)/2+istage)
            HC = ros_C((j-1)*(j-2)/2+istage)/(Direction*H)
            DO m = 1,NADJ
              CALL DAXPY(N,HA,V(jstart,m),1,U(istart,m),1) 
              CALL DAXPY(N,HC,U(jstart,m),1,U(istart,m),1) 
            END DO ! m=1:NADJ
          END DO

          DO m = 1,NADJ
            CALL LSS_Solve(Transp,U(istart,m))
            ISTATUS(Nsol) = ISTATUS(Nsol) + 1
          END DO ! m=1:NADJ

!~~~> Compute V 
          Tau = T + ros_Alpha(istage)*Direction*H
!       CALL JACO(N, Tau,Ystage(istart),Jac)
!~~~~> jac0 jac1 exchanged
          CALL LSS_Jac(Tau,Ystage(istart),JAC)
          ISTATUS(Njac) = ISTATUS(Njac) + 1
          IF(GetQuad) CALL DRDY(NADJ,N,N,Tau,Ystage(istart),RY1)
          DO m = 1,NADJ
            CALL LSS_Mul_Jactr(V(istart,m),U(istart,m))
            IF(GetQuad) CALL DAXPY(N,ros_W(istage),RY1(1,m),1,V(istart,m),1)
          END DO ! m=1:NADJ

        END DO Stage     
!~~~> Update the sensitivity variable
   
!~~~>  Compute Lambda 
        DO istage=1,ros_S
          istart = N*(istage-1) + 1
          DO m = 1,NADJ
            ! Add V_i
            CALL DAXPY(N,ONE,V(istart,m),1,Lambda(1,m),1)
            ! Add (H0xK_i)^T * U_i
            CALL HESSTR_VEC(N,T,Ystage(1),U(istart,m),K(istart),TMP)
!           CALL HessTR_Vec ( NHESS, Hes0, U(istart,m), K(istart), Tmp )
            CALL DAXPY(N,ONE,Tmp,1,Lambda(1,m),1)
           ! Add (HP0xK_i)^T * U_i
            IF(GetQuad) THEN
              CALL HESSTR_VEC_R(m,N,T,Ystage(1),ros_W(istage),K(istart),TMP)
              CALL DAXPY(N,ONE,Tmp,1,Lambda(1,m),1)
            END IF
          END DO ! m=1:NADJ
        END DO

!~~~> NonAutonomous term
        IF (.NOT.Autonomous) THEN
          ! Approximate d(f_y)/dt
          Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
          CALL LSS_Jac2(T+Delta,Ystage(1),JAC)
          ISTATUS(Njac) = ISTATUS(Njac) + 1
          CALL LSS_Jac_Time_Derivative(Delta)
          DO m = 1,NADJ
            Tmp(1:N) = ZERO
            DO istage = 1, ros_S
              istart = N*(istage-1) + 1
              CALL DAXPY(N,ros_Gamma(istage),U(istart,m),1,Tmp,1)
            END DO  
            ! Add H * dJac_dT_0^T * \sum(gamma_i U_i) to Lambda
            CALL LSS_Mul_djdttr(Tmp2,Tmp)
            CALL DAXPY(N,H,Tmp2,1,Lambda(1,m),1)
          END DO ! m=1:NADJ
        END IF

        IF (GetQuad .AND. .NOT. QuadAutonomous) THEN
          CALL DRDY(NADJ,N,N,T+Delta,Ystage(1),RY1)
          Alpha = -1.0d0
          CALL DAXPY(N*NADJ,Alpha,RY0,1,RY1,1)
          Alpha = 1/Delta
          CALL DSCAL(N*NADJ,Alpha,RY1,1)
          Alpha = H*GammaUhat
          CALL DAXPY(N*NADJ,Alpha,RY1,1,Lambda,1)
        END IF
 
      END DO TimeLoop 
   
      DEALLOCATE(U, V,STAT=j)
      IF(j .NE. 0) STOP 'deallocation error for U,V'

      IF(GetQuad) THEN
        DEALLOCATE(RY0,RY1,STAT=j)
        IF(j .NE. 0) STOP 'deallocation error for RYs'
      END IF 
!~~~> Succesful exit
      IERR = 1  !~~~> The integration was successful

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END SUBROUTINE ros_DadjInt
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_CadjInt (                        &
        NADJ, Y,                                 &
        Tstart, Tend, T,                         &
        AbsTol_adj, RelTol_adj,                  &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockADJ method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE LS_Solver
  IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
   INTEGER, INTENT(IN) :: NADJ
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
!~~~> Input: integration interval   
   DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend      
!~~~> Input: adjoint tolerances   
   DOUBLE PRECISION, INTENT(IN) :: AbsTol_adj(N,NADJ), RelTol_adj(N,NADJ)
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   DOUBLE PRECISION, INTENT(OUT) ::  T      
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables        
   DOUBLE PRECISION :: Y0(N)
   DOUBLE PRECISION :: Ynew(N,NADJ), Fcn0(N,NADJ), Fcn(N,NADJ) 
   DOUBLE PRECISION :: K(N*ros_S,NADJ), dFdT(N,NADJ)
   DOUBLE PRECISION :: H, Hnew, HC, HG, Fac, Tau 
   DOUBLE PRECISION :: Err, Yerr(N,NADJ)
   INTEGER :: Direction, ioffset, j, istage,iadj
   LOGICAL :: RejectLastH, RejectMoreH, Singular, Transp
!~~~>  Local parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
   DOUBLE PRECISION :: Delta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = 0.0
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.0*Roundoff) H = DeltaMin
   
   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF               
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.
   
!~~~> Time loop begins below 

TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) ) 
      
   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1d0*H) == T).OR.(H <= Roundoff) ) THEN  ! Step size too small
      CALL ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   
!~~~>  Limit H if necessary to avoid going beyond Tend   
   RSTATUS(Nhexit) = H
   H = MIN(H,ABS(Tend-T))

!~~~>   Interpolate forward solution
   CALL ros_cadj_Y( T, Y0 )     
!~~~>   Compute the Jacobian at current time
   CALL LSS_Jac(T,Y0,JAC)
   ISTATUS(Njac) = ISTATUS(Njac) + 1
   
!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
       Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
       CALL LSS_Jac2(T+Delta,Y0,JAC)
       ISTATUS(Njac) = ISTATUS(Njac) + 1
       CALL LSS_Jac_Time_Derivative(Delta)
      DO iadj = 1, NADJ
        CALL LSS_Mul_djdttr(dFdT(1,iadj),Y(1,iadj))
        CALL DSCAL(N,(-ONE),dFdT(1,iadj),1)
      END DO
   END IF

   CALL LSS_Rev_Jac()
!~~~>  Ydot = -J^T*Y
  
   DO iadj = 1, NADJ
     CALL LSS_Mul_Jactr(Fcn0(1,iadj),Y(1,iadj))
   END DO
    
!~~~>  Repeat step calculation until current step accepted
UntilAccepted: DO  
   
   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S
      
      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)
      
      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         DO iadj = 1, NADJ
           CALL DCOPY(N,Fcn0(1,iadj),1,Fcn(1,iadj),1)
         END DO
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL DCOPY(N*NADJ,Y,1,Ynew,1)
         DO j = 1, istage-1
           DO iadj = 1, NADJ
             CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
                K(N*(j-1)+1,iadj),1,Ynew(1,iadj),1) 
           END DO       
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL ros_cadj_Y( Tau, Y0 )     
         CALL LSS_Jac1(Tau,Y0,JAC)
         ISTATUS(Njac) = ISTATUS(Njac) + 1
         CALL LSS_Rev_Jac1()
         DO iadj = 1, NADJ
             CALL LSS_Mul_Jactr1(Fcn(1,iadj),Ynew(1,iadj))
         END DO
       END IF ! if istage == 1 elseif ros_NewF(istage)

       DO iadj = 1, NADJ
          CALL DCOPY(N,Fcn(1,iadj),1,K(ioffset+1,iadj),1)
       END DO
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HC,K(N*(j-1)+1,iadj),1, &
                  K(ioffset+1,iadj),1)
         END DO
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HG,dFdT(1,iadj),1,K(ioffset+1,iadj),1)
         END DO
       END IF
       
       Transp = .TRUE.
       DO iadj = 1, NADJ
          CALL LSS_Solve(Transp,K(ioffset+1,iadj))
          ISTATUS(Nsol) = ISTATUS(Nsol) + 1
       END DO
      
   END DO Stage     
            

!~~~>  Compute the new solution 
   DO iadj = 1, NADJ
      CALL DCOPY(N,Y(1,iadj),1,Ynew(1,iadj),1)
      DO j=1,ros_S
         CALL DAXPY(N,ros_M(j),K(N*(j-1)+1,iadj),1,Ynew(1,iadj),1)
      END DO
   END DO

!~~~>  Compute the error estimation 
   CALL DSCAL(N*NADJ,ZERO,Yerr,1)
   DO j=1,ros_S     
       DO iadj = 1, NADJ
        CALL DAXPY(N,ros_E(j),K(N*(j-1)+1,iadj),1,Yerr(1,iadj),1)
       END DO
   END DO
!~~~> Max error among all adjoint components    
   iadj = 1
   Err = ros_ErrorNorm ( Y(1,iadj), Ynew(1,iadj), Yerr(1,iadj), &
              AbsTol_adj(1,iadj), RelTol_adj(1,iadj), VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac  

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      CALL DCOPY(N*NADJ,Ynew,1,Y,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  ! No step size increase after a rejected step
         Hnew = MIN(Hnew,H) 
      END IF   
      RSTATUS(Nhexit) = H
      RSTATUS(Nhnew)  = Hnew
      RSTATUS(Ntexit) = T
      RejectLastH = .FALSE.  
      RejectMoreH = .FALSE.
      H = Hnew      
      EXIT UntilAccepted ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
   ELSE           !~~~> Reject step
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF   
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (ISTATUS(Nacc) >= 1) THEN
         ISTATUS(Nrej) = ISTATUS(Nrej) + 1
      END IF    
   END IF ! Err <= 1

   END DO UntilAccepted 

   END DO TimeLoop 
      
!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful

  END SUBROUTINE ros_CadjInt
  
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_SimpleCadjInt (                  &
        NADJ, Y,                                 &
        Tstart, Tend, T,                         &
!~~~> Error indicator
        IERR )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic RosenbrockADJ method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE LS_Solver
  IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
   INTEGER, INTENT(IN) :: NADJ
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N,NADJ)
!~~~> Input: integration interval   
   DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   DOUBLE PRECISION, INTENT(OUT) ::  T      
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables        
   DOUBLE PRECISION :: Y0(N)
   DOUBLE PRECISION :: Ynew(N,NADJ), Fcn0(N,NADJ), Fcn(N,NADJ) 
   DOUBLE PRECISION :: K(N*ros_S,NADJ), dFdT(N,NADJ)
   DOUBLE PRECISION :: H, HC, HG, Tau 
   DOUBLE PRECISION :: ghinv
   INTEGER :: Direction, ioffset, j, istage, iadj
   INTEGER :: istack
!~~~>  Local parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0 
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
   DOUBLE PRECISION :: Delta
   Logical :: Transp
!~~~>  Locally called functions
!    DOUBLE PRECISION WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
!~~~>  INITIAL PREPARATIONS
   
   IF (Tend  >=  Tstart) THEN
     Direction = -1
   ELSE
     Direction = +1
   END IF               
   
!~~~> Time loop begins below 
TimeLoop: DO istack = stack_ptr,2,-1
        
   T = chk_T(istack)
   H = chk_H(istack-1)
   !CALL DCOPY(N,chk_Y(1,istack),1,Y0,1)
   Y0(1:N) = chk_Y(1:N,istack)
   
!~~~>   Compute the Jacobian at current time
   CALL LSS_Jac(T,Y0,JAC)
   ISTATUS(Njac) = ISTATUS(Njac) + 1
   
!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T)) 
      CALL LSS_Jac2(T+Delta,Y0,JAC)
      ISTATUS(Njac) = ISTATUS(Njac) +1
      CALL LSS_Jac_Time_Derivative (Delta)
      DO iadj = 1, NADJ
        CALL LSS_Mul_djdttr(dFdT(1,iadj),Y(1,iadj))
        CALL DSCAL(N,(-ONE),dFdT(1,iadj),1)
      END DO
   END IF

!~~~>  Ydot = -J^T*Y
   CALL LSS_Rev_Jac()
   DO iadj = 1, NADJ
     CALL LSS_Mul_Jac(Fcn0(1,iadj),Y(1,iadj))
   END DO
   
!~~~>    Construct Ghimj = 1/(H*ham) - Jac0
     ghinv = ONE/(Direction*H*ros_Gamma(1))
     CALL LSS_Decomp(ghinv,j)
     ISTATUS(Ndec) = ISTATUS(Ndec) + 1
     IF (j /= 0) THEN
       CALL ros_ErrorMsg(-8,T,H,IERR)
       PRINT*,' The matrix is singular !'
       STOP
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S
      
      ! Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)
      
      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         DO iadj = 1, NADJ
           CALL DCOPY(N,Fcn0(1,iadj),1,Fcn(1,iadj),1)
         END DO
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL DCOPY(N*NADJ,Y,1,Ynew,1)
         DO j = 1, istage-1
           DO iadj = 1, NADJ
             CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
                K(N*(j-1)+1,iadj),1,Ynew(1,iadj),1) 
           END DO       
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL ros_Hermite3( chk_T(istack-1), chk_T(istack), Tau, &
             chk_Y(1:N,istack-1), chk_Y(1:N,istack),       &
             chk_dY(1:N,istack-1), chk_dY(1:N,istack), Y0 )
         CALL LSS_Jac1(Tau,Y0,JAC)
         ISTATUS(Njac) = ISTATUS(Njac) + 1
         CALL LSS_Rev_Jac1()
         DO iadj = 1, NADJ
             CALL LSS_Mul_Jactr1(Fcn(1,iadj),Ynew(1,iadj))
         END DO
       END IF ! if istage == 1 elseif ros_NewF(istage)

       DO iadj = 1, NADJ
          CALL DCOPY(N,Fcn(1,iadj),1,K(ioffset+1,iadj),1)
       END DO
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HC,K(N*(j-1)+1,iadj),1, &
                  K(ioffset+1,iadj),1)
         END DO
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         DO iadj = 1, NADJ
           CALL DAXPY(N,HG,dFdT(1,iadj),1,K(ioffset+1,iadj),1)
         END DO
       END IF
  
       Transp = .TRUE.
       DO iadj = 1, NADJ
         CALL LSS_Solve(Transp, K(ioffset+1,iadj))
         ISTATUS(Nsol) = ISTATUS(Nsol) + 1
       END DO
      
   END DO Stage     
            

!~~~>  Compute the new solution 
   DO iadj = 1, NADJ
      DO j=1,ros_S
         CALL DAXPY(N,ros_M(j),K(N*(j-1)+1,iadj),1,Y(1,iadj),1)
      END DO
   END DO

   END DO TimeLoop 
      
!~~~> Succesful exit
   IERR = 1  !~~~> The integration was successful
  END SUBROUTINE ros_SimpleCadjInt
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DOUBLE PRECISION FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, & 
               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

! Input arguments   
   DOUBLE PRECISION, INTENT(IN) :: Y(N), Ynew(N), &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables     
   DOUBLE PRECISION :: Err, Scale, Ymax   
   INTEGER  :: i
   
   Err = ZERO
   DO i=1,N
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = MAX(Err,1.0d-10)
   
  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

!~~~> Input arguments   
   DOUBLE PRECISION, INTENT(IN) :: T, Roundoff, Y(N), Fcn0(N) 
!~~~> Output arguments   
   DOUBLE PRECISION, INTENT(OUT) :: dFdT(N)   
!~~~> Local variables     
   DOUBLE PRECISION :: Delta  
   DOUBLE PRECISION, PARAMETER :: ONE = 1.0d0, DeltaMin = 1.0d-6
   
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FUN(N,T+Delta,Y,dFdT)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL DAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL DSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix (H,Direction,gam,Singular)
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*gam) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   USE LS_Solver
   IMPLICIT NONE
  
!~~~> Input arguments
!
   DOUBLE PRECISION, INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
   LOGICAL, INTENT(OUT) ::  Singular
!   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   DOUBLE PRECISION, INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: ising, Nconsecutive
   DOUBLE PRECISION :: ghinv
   DOUBLE PRECISION, PARAMETER :: ONE  = 1.0, HALF = 0.5

   Nconsecutive = 0
   Singular = .TRUE.

   
   DO WHILE (Singular)
      ghinv = ONE/(Direction*H*gam)

      CALL LSS_Decomp(ghinv,ISING)
      ISTATUS(Ndec) = ISTATUS(Ndec) + 1

     IF (ising == 0) THEN
!~~~>    If successful done
        Singular = .FALSE.
     ELSE ! ising .ne. 0
!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
           H = H*HALF
        ELSE  ! More than 5 consecutive failed decompositions
           RETURN
        END IF  ! Nconsecutive
      END IF    ! ising

   END DO ! WHILE Singular

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE ros_cadj_Y( T, Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!  Finds the solution Y at T by interpolating the stored forward trajectory
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
   IMPLICIT NONE
!~~~> Input variables 
   DOUBLE PRECISION, INTENT(IN) :: T
!~~~> Output variables     
   DOUBLE PRECISION, INTENT(OUT) :: Y(N)
!~~~> Local variables     
   INTEGER     :: i
   DOUBLE PRECISION, PARAMETER  :: ONE = 1.0d0

!   chk_H, chk_T, chk_Y, chk_dY, chk_d2Y

   IF( (T < chk_T(1)).OR.(T> chk_T(stack_ptr)) ) THEN
      PRINT*,'Cannot locate solution at T = ',T
      PRINT*,'Stored trajectory is between Tstart = ',chk_T(1)
      PRINT*,'    and Tend = ',chk_T(stack_ptr)
      STOP
   END IF
   DO i = 1, stack_ptr-1
     IF( (T>= chk_T(i)).AND.(T<= chk_T(i+1)) ) EXIT
   END DO 


!   IF (.FALSE.) THEN
!
!   CALL ros_Hermite5( chk_T(i), chk_T(i+1), T, &
!                chk_Y(1,i),   chk_Y(1,i+1),     &
!                chk_dY(1,i),  chk_dY(1,i+1),    &
!                chk_d2Y(1,i), chk_d2Y(1,i+1), Y )
!   
!   ELSE
                
   CALL ros_Hermite3( chk_T(i), chk_T(i+1), T, &
                chk_Y(1:N,i),   chk_Y(1:N,i+1),     &
                chk_dY(1:N,i),  chk_dY(1:N,i+1),    &
                Y )
                        
!   
!   END IF       

  END SUBROUTINE ros_cadj_Y
  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE ros_Hermite3( a, b, T, Ya, Yb, Ja, Jb, Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!  Template for Hermite interpolation of order 5 on the interval [a,b]
! P = c(1) + c(2)*(x-a) + ... + c(4)*(x-a)^3
! P[a,b] = [Ya,Yb], P'[a,b] = [Ja,Jb]
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
   IMPLICIT NONE
!~~~> Input variables 
   DOUBLE PRECISION, INTENT(IN) :: a, b, T, Ya(N), Yb(N)
   DOUBLE PRECISION, INTENT(IN) :: Ja(N), Jb(N)
!~~~> Output variables     
   DOUBLE PRECISION, INTENT(OUT) :: Y(N)
!~~~> Local variables     
   DOUBLE PRECISION :: Tau, amb(3), C(N,4)
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0
   INTEGER :: i, j
   
   amb(1) = 1.0d0/(a-b)
   DO i=2,3
     amb(i) = amb(i-1)*amb(1)
   END DO
   
   
! c(1) = ya;
   CALL DCOPY(N,Ya,1,C(1,1),1)
! c(2) = ja;
   CALL DCOPY(N,Ja,1,C(1,2),1)
! c(3) = 2/(a-b)*ja + 1/(a-b)*jb - 3/(a - b)^2*ya + 3/(a - b)^2*yb  ;
   CALL DCOPY(N,Ya,1,C(1,3),1)
   CALL DSCAL(N,-3.0*amb(2),C(1,3),1)
   CALL DAXPY(N,3.0*amb(2),Yb,1,C(1,3),1)
   CALL DAXPY(N,2.0*amb(1),Ja,1,C(1,3),1)
   CALL DAXPY(N,amb(1),Jb,1,C(1,3),1)
! c(4) =  1/(a-b)^2*ja + 1/(a-b)^2*jb - 2/(a-b)^3*ya + 2/(a-b)^3*yb ;
   CALL DCOPY(N,Ya,1,C(1,4),1)
   CALL DSCAL(N,-2.0*amb(3),C(1,4),1)
   CALL DAXPY(N,2.0*amb(3),Yb,1,C(1,4),1)
   CALL DAXPY(N,amb(2),Ja,1,C(1,4),1)
   CALL DAXPY(N,amb(2),Jb,1,C(1,4),1)
   
   Tau = T - a
   CALL DCOPY(N,C(1,4),1,Y,1)
   CALL DSCAL(N,Tau**3,Y,1)
   DO j = 3,1,-1
     CALL DAXPY(N,TAU**(j-1),C(1,j),1,Y,1)
   END DO       

  END SUBROUTINE ros_Hermite3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE ros_Hermite5( a, b, T, Ya, Yb, Ja, Jb, Ha, Hb, Y )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!  Template for Hermite interpolation of order 5 on the interval [a,b]
! P = c(1) + c(2)*(x-a) + ... + c(6)*(x-a)^5
! P[a,b] = [Ya,Yb], P'[a,b] = [Ja,Jb], P"[a,b] = [Ha,Hb]
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
   IMPLICIT NONE
!~~~> Input variables 
   DOUBLE PRECISION, INTENT(IN) :: a, b, T, Ya(N), Yb(N)
   DOUBLE PRECISION, INTENT(IN) :: Ja(N), Jb(N), Ha(N), Hb(N)
!~~~> Output variables     
   DOUBLE PRECISION, INTENT(OUT) :: Y(N)
!~~~> Local variables     
   DOUBLE PRECISION :: Tau, amb(5), C(N,6)
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, HALF = 0.5d0
   INTEGER :: i, j
   
   amb(1) = 1.0d0/(a-b)
   DO i=2,5
     amb(i) = amb(i-1)*amb(1)
   END DO
     
! c(1) = ya;
   CALL DCOPY(N,Ya,1,C(1,1),1)
! c(2) = ja;
   CALL DCOPY(N,Ja,1,C(1,2),1)
! c(3) = ha/2;
   CALL DCOPY(N,Ha,1,C(1,3),1)
   CALL DSCAL(N,HALF,C(1,3),1)
   
! c(4) = 10*amb(3)*ya - 10*amb(3)*yb - 6*amb(2)*ja - 4*amb(2)*jb  + 1.5*amb(1)*ha - 0.5*amb(1)*hb ;
   CALL DCOPY(N,Ya,1,C(1,4),1)
   CALL DSCAL(N,10.0*amb(3),C(1,4),1)
   CALL DAXPY(N,-10.0*amb(3),Yb,1,C(1,4),1)
   CALL DAXPY(N,-6.0*amb(2),Ja,1,C(1,4),1)
   CALL DAXPY(N,-4.0*amb(2),Jb,1,C(1,4),1)
   CALL DAXPY(N, 1.5*amb(1),Ha,1,C(1,4),1)
   CALL DAXPY(N,-0.5*amb(1),Hb,1,C(1,4),1)

! c(5) =   15*amb(4)*ya - 15*amb(4)*yb - 8.*amb(3)*ja - 7*amb(3)*jb + 1.5*amb(2)*ha - 1*amb(2)*hb ;
   CALL DCOPY(N,Ya,1,C(1,5),1)
   CALL DSCAL(N, 15.0*amb(4),C(1,5),1)
   CALL DAXPY(N,-15.0*amb(4),Yb,1,C(1,5),1)
   CALL DAXPY(N,-8.0*amb(3),Ja,1,C(1,5),1)
   CALL DAXPY(N,-7.0*amb(3),Jb,1,C(1,5),1)
   CALL DAXPY(N,1.5*amb(2),Ha,1,C(1,5),1)
   CALL DAXPY(N,-amb(2),Hb,1,C(1,5),1)
   
! c(6) =   6*amb(5)*ya - 6*amb(5)*yb - 3.*amb(4)*ja - 3.*amb(4)*jb + 0.5*amb(3)*ha -0.5*amb(3)*hb ;
   CALL DCOPY(N,Ya,1,C(1,6),1)
   CALL DSCAL(N, 6.0*amb(5),C(1,6),1)
   CALL DAXPY(N,-6.0*amb(5),Yb,1,C(1,6),1)
   CALL DAXPY(N,-3.0*amb(4),Ja,1,C(1,6),1)
   CALL DAXPY(N,-3.0*amb(4),Jb,1,C(1,6),1)
   CALL DAXPY(N, 0.5*amb(3),Ha,1,C(1,6),1)
   CALL DAXPY(N,-0.5*amb(3),Hb,1,C(1,6),1)
   
   Tau = T - a
   CALL DCOPY(N,C(1,6),1,Y,1)
   DO j = 5,1,-1
     CALL DSCAL(N,Tau,Y,1)
     CALL DAXPY(N,ONE,C(1,j),1,Y,1)
   END DO       

  END SUBROUTINE ros_Hermite5
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

    IMPLICIT NONE
    DOUBLE PRECISION g
   
    g = 1.0d0 + 1.0d0/SQRT(2.0d0)
   
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'   
!~~~> Number of stages
    ros_S = 2
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
   
    ros_A(1) = (1.d0)/g
    ros_C(1) = (-2.d0)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.d0)/(2.d0*g)
    ros_M(2)= (1.d0)/(2.d0*g)
! E_i = Coefficients for error estimator    
    ros_E(1) = 1.d0/(2.d0*g)
    ros_E(2) = 1.d0/(2.d0*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0d0
    ros_Alpha(2) = 1.0d0 
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = g
    ros_Gamma(2) =-g
   
 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   IMPLICIT NONE
    
    rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'   
!~~~> Number of stages
   ros_S = 3
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
   ros_A(1)= 1.d0
   ros_A(2)= 1.d0
   ros_A(3)= 0.d0

   ros_C(1) = -0.10156171083877702091975600115545d+01
   ros_C(2) =  0.40759956452537699824805835358067d+01
   ros_C(3) =  0.92076794298330791242156818474003d+01
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1d+01
   ros_M(2) =  0.61697947043828245592553615689730d+01
   ros_M(3) = -0.42772256543218573326238373806514d+00
! E_i = Coefficients for error estimator    
   ros_E(1) =  0.5d+00
   ros_E(2) = -0.29079558716805469821718236208017d+01
   ros_E(3) =  0.22354069897811569627360909276199d+00
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0d+00
   ros_Alpha(2)= 0.43586652150845899941601945119356d+00
   ros_Alpha(3)= 0.43586652150845899941601945119356d+00
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1)= 0.43586652150845899941601945119356d+00
   ros_Gamma(2)= 0.24291996454816804366592249683314d+00
   ros_Gamma(3)= 0.21851380027664058511513169485832d+01

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3 
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

   IMPLICIT NONE
   
    rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'   
!~~~> Number of stages
   ros_S = 4
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
   ros_A(1) = 0.2000000000000000d+01
   ros_A(2) = 0.1867943637803922d+01
   ros_A(3) = 0.2344449711399156d+00
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0D0

   ros_C(1) =-0.7137615036412310d+01
   ros_C(2) = 0.2580708087951457d+01
   ros_C(3) = 0.6515950076447975d+00
   ros_C(4) =-0.2137148994382534d+01
   ros_C(5) =-0.3214669691237626d+00
   ros_C(6) =-0.6949742501781779d+00
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735d+01
   ros_M(2) = 0.2870493262186792d+00
   ros_M(3) = 0.4353179431840180d+00
   ros_M(4) = 0.1093502252409163d+01
!~~~> E_i  = Coefficients for error estimator    
   ros_E(1) =-0.2815431932141155d+00
   ros_E(2) =-0.7276199124938920d-01
   ros_E(3) =-0.1082196201495311d+00
   ros_E(4) =-0.1093502252409163d+01
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0d0    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.D0
   ros_Alpha(2) = 0.1145640000000000d+01
   ros_Alpha(3) = 0.6552168638155900d+00
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1) = 0.5728200000000000d+00
   ros_Gamma(2) =-0.1769193891319233d+01
   ros_Gamma(3) = 0.7592633437920482d+00
   ros_Gamma(4) =-0.1049021087100450d+00

  END SUBROUTINE Ros4
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

   IMPLICIT NONE
   
    rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'   
!~~~> Number of stages
   ros_S = 4
   
!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
 
   ros_A(1) = 0.0d+00
   ros_A(2) = 2.0d+00
   ros_A(3) = 0.0d+00
   ros_A(4) = 2.0d+00
   ros_A(5) = 0.0d+00
   ros_A(6) = 1.0d+00

   ros_C(1) = 4.0d+00
   ros_C(2) = 1.0d+00
   ros_C(3) =-1.0d+00
   ros_C(4) = 1.0d+00
   ros_C(5) =-1.0d+00 
   ros_C(6) =-(8.0d+00/3.0d+00) 
         
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0d+00
   ros_M(2) = 0.0d+00
   ros_M(3) = 1.0d+00
   ros_M(4) = 1.0d+00
!~~~> E_i  = Coefficients for error estimator    
   ros_E(1) = 0.0d+00
   ros_E(2) = 0.0d+00
   ros_E(3) = 0.0d+00
   ros_E(4) = 1.0d+00
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0d+00    
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0d+00
   ros_Alpha(2) = 0.0d+00
   ros_Alpha(3) = 1.0d+00
   ros_Alpha(4) = 1.0d+00
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
   ros_Gamma(1) = 0.5d+00
   ros_Gamma(2) = 1.5d+00
   ros_Gamma(3) = 0.0d+00
   ros_Gamma(4) = 0.0d+00

  END SUBROUTINE Rodas3
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'   
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000d0
    ros_Alpha(2) = 0.386d0
    ros_Alpha(3) = 0.210d0 
    ros_Alpha(4) = 0.630d0
    ros_Alpha(5) = 1.000d0
    ros_Alpha(6) = 1.000d0
        
!~~~> Gamma_i = \sum_j  gamma_{i,j}    
    ros_Gamma(1) = 0.2500000000000000d+00
    ros_Gamma(2) =-0.1043000000000000d+00
    ros_Gamma(3) = 0.1035000000000000d+00
    ros_Gamma(4) =-0.3620000000000023d-01
    ros_Gamma(5) = 0.0d0
    ros_Gamma(6) = 0.0d0

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )    
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )  
     
    ros_A(1) = 0.1544000000000000d+01
    ros_A(2) = 0.9466785280815826d+00
    ros_A(3) = 0.2557011698983284d+00
    ros_A(4) = 0.3314825187068521d+01
    ros_A(5) = 0.2896124015972201d+01
    ros_A(6) = 0.9986419139977817d+00
    ros_A(7) = 0.1221224509226641d+01
    ros_A(8) = 0.6019134481288629d+01
    ros_A(9) = 0.1253708332932087d+02
    ros_A(10) =-0.6878860361058950d+00
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0d+00

    ros_C(1) =-0.5668800000000000d+01
    ros_C(2) =-0.2430093356833875d+01
    ros_C(3) =-0.2063599157091915d+00
    ros_C(4) =-0.1073529058151375d+00
    ros_C(5) =-0.9594562251023355d+01
    ros_C(6) =-0.2047028614809616d+02
    ros_C(7) = 0.7496443313967647d+01
    ros_C(8) =-0.1024680431464352d+02
    ros_C(9) =-0.3399990352819905d+02
    ros_C(10) = 0.1170890893206160d+02
    ros_C(11) = 0.8083246795921522d+01
    ros_C(12) =-0.7981132988064893d+01
    ros_C(13) =-0.3152159432874371d+02
    ros_C(14) = 0.1631930543123136d+02
    ros_C(15) =-0.6058818238834054d+01

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0d+00
    ros_M(6) = 1.0d+00

!~~~> E_i  = Coefficients for error estimator    
    ros_E(1) = 0.0d+00
    ros_E(2) = 0.0d+00
    ros_E(3) = 0.0d+00
    ros_E(4) = 0.0d+00
    ros_E(5) = 0.0d+00
    ros_E(6) = 1.0d+00

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.
     
!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0d0
     
  END SUBROUTINE Rodas4

END SUBROUTINE RosenbrockADJ2 ! and its internal procedures

END MODULE ROS_adj_f90_Integrator




! End of INTEGRATE function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


