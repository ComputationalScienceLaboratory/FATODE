! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! INTEGRATE - Integrator routine
!   Arguments :
!      TIN       - Start Time for Integration
!      TOUT      - End Time for Integration
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Rosenbrock_TLM - Implementation of the Tangent Linear Model            !
!               for several Rosenbrock methods:                           !
!               * Ros2                                                    !
!               * Ros3                                                    !
!               * Ros4                                                    !
!               * Rodas3                                                  !
!               * Rodas4                                                  !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

MODULE ROS_tlm_f90_Integrator

   IMPLICIT NONE
   PUBLIC
   SAVE

!~~~>  Statistics on the work performed by the Rosenbrock method
  INTEGER, PARAMETER :: Nfun=1, Njac=2, Nstp=3, Nacc=4, &
                        Nrej=5, Ndec=6, Nsol=7, Nsng=8, &
                        Nhes=9, Ntexit=1, Nhexit=2, Nhnew = 3


CONTAINS ! Functions in the module ROS_tlm_f90_Integrator

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE INTEGRATE_TLM( N, NTLM, NNZERO, Y, Y_tlm, TIN, TOUT, ATOL_tlm, RTOL_tlm,&
       ATOL,RTOL,FUN,JAC,HESS_VEC,ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE        
   INTEGER ::N, NNZERO
!~~~> Y - Concentrations
   DOUBLE PRECISION :: Y(N)
!~~~> NTLM - No. of sensitivity coefficients
   INTEGER NTLM
!~~~> Y_tlm - Sensitivities of concentrations
!     Note: Y_tlm (1:N,j) contains sensitivities of
!               Y(1:N) w.r.t. the j-th parameter, j=1...NTLM
   DOUBLE PRECISION :: Y_tlm(N,NTLM)
   DOUBLE PRECISION, INTENT(IN)         :: TIN  ! TIN  - Start Time
   DOUBLE PRECISION, INTENT(IN)         :: TOUT ! TOUT - End Time
!~~~> Optional input parameters and statistics
   INTEGER,  INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RTOL_tlm(N,NTLM),ATOL_tlm(N,NTLM),&
                                        ATOL(N),RTOL(N)

   INTEGER, SAVE :: IERR
   DOUBLE PRECISION ::  RCNTRL(20), RSTATUS(20)
   INTEGER ::   ICNTRL(20), ISTATUS(20)
   DOUBLE PRECISION,SAVE :: STEPMIN= 0.0
   EXTERNAL :: FUN,JAC,HESS_VEC
   ICNTRL(1:20)  = 0
   RCNTRL(1:20)  = 0.0
   ISTATUS(1:20) = 0
   RSTATUS(1:20) = 0.0
   
   ICNTRL(1) = 0       ! non-autonomous
   ICNTRL(2) = 1       ! vector tolerances
   ICNTRL(3) = 5       ! choice of the method
   ICNTRL(12) = 1      ! 0 - fwd trunc error only, 1 - tlm trunc error
   RCNTRL(3) = STEPMIN ! starting step

   ! if optional parameters are given, and if they are >=0, then they overwrite default settings
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) >= 0) ICNTRL(:) = ICNTRL_U(:)
   ENDIF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) >= 0) RCNTRL(:) = RCNTRL_U(:)
   ENDIF
   
   CALL RosenbrockTLM(N, NNZERO, Y, NTLM, Y_tlm,      &
         TIN,TOUT,ATOL,RTOL,ATOL_tlm,RTOL_tlm,     &
         FUN,JAC,HESS_VEC,RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
    
   IF (IERR < 0) THEN
     print *,'Rosenbrock: Unsucessful step at T=', &
         TIN,' (IERR=',IERR,')'
   END IF
   
   STEPMIN = RSTATUS(Nhexit)
   ! if optional parameters are given for output they return information
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)

END SUBROUTINE INTEGRATE_TLM


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE RosenbrockTLM(N,NNZERO,Y,NTLM,Y_tlm,                      &
           Tstart,Tend,AbsTol,RelTol,AbsTol_tlm,RelTol_tlm,   &
           FUN,JAC,HESS_VEC,RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   
!    TLM = Tangent Linear Model of a Rosenbrock Method
!
!    Solves the system y'=F(t,y) using a Rosenbrock method defined by:
!
!     G = 1/(H*gamma(1)) - Jac(t0,Y0)
!     T_i = t0 + Alpha(i)*H
!     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
!     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
!         gamma(i)*dF/dT(t0, Y0)
!     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j 
!
!    For details on Rosenbrock methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.  
!    The codes contained in the book inspired this implementation.       
!
!    (C)  Adrian Sandu, August 2004
!    Virginia Polytechnic Institute and State University    
!    Contact: sandu@cs.vt.edu
!    Revised by Philipp Miehe and Adrian Sandu, May 2006   
!    Revided by Hong Zhang and Adrian Sandu, Feb 2011                     !
!    This implementation is part of FATODE                  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    
!~~~>   INPUT ARGUMENTS: 
!    
!-     Y(N)    -> vector of initial conditions (at T=Tstart)
!      NTLM       -> dimension of linearized system, 
!                   i.e. the number of sensitivity coefficients
!-     Y_tlm(N*NTLM) -> vector of initial sensitivity conditions (at T=Tstart)
!-    [Tstart,Tend]    -> time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)  
!-    RelTol, AbsTol -> user precribed accuracy
!- SUBROUTINE Fun( T, Y, Ydot ) -> ODE function, 
!                       returns Ydot = Y' = F(T,Y) 
!- SUBROUTINE Jac( T, Y, Jcb ) -> Jacobian of the ODE function,
!                       returns Jcb = dF/dY 
!-    ICNTRL(1:20)    -> integer inputs parameters
!-    RCNTRL(1:20)    -> real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:  
!     
!-    Y(N)         -> vector of final states (at T->Tend)
!-    Y_tlm(N*NTLM)-> vector of final sensitivities (at T=Tend)
!-    ISTATUS(1:20)   -> integer output parameters
!-    RSTATUS(:20)    -> real output parameters
!-    IERR            -> job status upon return
!       - succes (positive value) or failure (negative value) -
!           =  1 : Success
!           = -1 : Improper value for maximal no of steps
!           = -2 : Selected Rosenbrock method not implemented
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
!        = 5 :  method is  Rodas4
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0) the default value of 100000 is used
!
!    ICNTRL(12) -> switch for TLM truncation error control
!		ICNTRL(12) = 0: TLM error is not used
!		ICNTRL(12) = 1: TLM error is computed and used
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO 
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!          
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMin,upper bound on step increase factor (default=6)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                       (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller 
!         than the predicted value  (default=0.9)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
!~~~>     OUTPUT PARAMETERS:
!
!    Note: each call to Rosenbrock adds the corrent no. of fcn calls
!      to previous value of ISTATUS(1), and similar for the other params.
!      Set ISTATUS(1:10) = 0 before call to avoid this accumulation.
!
!    ISTATUS(1) = No. of function calls
!    ISTATUS(2) = No. of Jacobian calls
!    ISTATUS(3) = No. of steps
!    ISTATUS(4) = No. of accepted steps
!    ISTATUS(5) = No. of rejected steps (except at the beginning)
!    ISTATUS(6) = No. of LU decompositions
!    ISTATUS(7) = No. of forward/backward substitutions
!    ISTATUS(8) = No. of singular matrix decompositions
!    ISTATUS(9) = No. of Hessian calls
!
!    RSTATUS(1)  -> Texit, the time corresponding to the 
!                   computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    For multiple restarts, use Hexit as Hstart in the following run 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  USE LS_Solver
  IMPLICIT NONE
   
!~~~>  Arguments   
   INTEGER,       INTENT(IN)    :: N, NTLM, NNZERO
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
   DOUBLE PRECISION, INTENT(INOUT) :: Y_tlm(N,NTLM)
   DOUBLE PRECISION, INTENT(IN)    :: Tstart, Tend
   DOUBLE PRECISION, INTENT(IN)    :: AbsTol(N),RelTol(N)
   DOUBLE PRECISION, INTENT(IN)    :: AbsTol_tlm(N,NTLM),RelTol_tlm(N,NTLM)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   DOUBLE PRECISION, INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   DOUBLE PRECISION, INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5
   DOUBLE PRECISION :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables     
   DOUBLE PRECISION :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   DOUBLE PRECISION :: Hmin, Hmax, Hstart, Hexit
   DOUBLE PRECISION :: Texit
   INTEGER :: i, UplimTol, Max_no_steps 
   LOGICAL :: Autonomous, VectorTol, TLMtruncErr
!~~~>   Parameters
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0, ONE  = 1.0d0
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
   DOUBLE PRECISION :: DLAMCH
   EXTERNAL :: FUN,JAC,HESS_VEC
    
   ros_A = ZERO
   ros_C = ZERO
   ros_M = ZERO
   ros_E = ZERO
   ros_Alpha = ZERO
   ros_Gamma = ZERO
!~~~> Initialize the statistics
   IERR          = 0
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
     CASE (4)
       CALL Rodas3
     CASE (0,5)
       CALL Rodas4
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT
   
!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (Max_no_steps > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE 
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN      
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
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hexit < FacMax 
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
     
   CALL LSS_Init(N,NNZERO)
!~~~>  CALL Rosenbrock method   
   CALL ros_TLM_Int(Y, NTLM, Y_tlm,              &
        Tstart, Tend, Texit,                     & 
!  Error indicator
        IERR)
   CALL LSS_Free

   PRINT*,'Tangent Linear Model  STATISTICS'
   PRINT*,'Step=',ISTATUS(Nstp),' Acc=',ISTATUS(Nacc),   &
        ' Rej=',ISTATUS(Nrej), ' Singular=',Nsng,' Fun=',ISTATUS(Nfun),' Jac=',&
        ISTATUS(Njac),' Sol=',ISTATUS(Nsol),' Dec=',ISTATUS(Ndec)


CONTAINS ! Procedures internal to RosenbrockTLM
 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   DOUBLE PRECISION, INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)    
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
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
     
 END SUBROUTINE ros_ErrorMsg
   
     
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_TLM_Int (Y, NTLM, Y_tlm,         &
        Tstart, Tend, T,                         &
!~~~> Error indicator
        IERR )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method 
!      defined by ros_S (no of stages)  
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE
   
!~~~> Input: the initial condition at Tstart; Output: the solution at T   
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
!~~~> Input: Number of sensitivity coefficients
   INTEGER, INTENT(IN) :: NTLM 
!~~~> Input: the initial sensitivites at Tstart; Output: the sensitivities at T   
   DOUBLE PRECISION, INTENT(INOUT) :: Y_tlm(N,NTLM)
!~~~> Input: integration interval   
   DOUBLE PRECISION, INTENT(IN) :: Tstart,Tend      
!~~~> Output: time at which the solution is returned (T=Tend if success)   
   DOUBLE PRECISION, INTENT(OUT) ::  T      
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR
! ~~~~ Local variables        
   DOUBLE PRECISION :: Ynew(N), Fcn0(N), Fcn(N) 
   DOUBLE PRECISION :: K(N*ros_S)
   DOUBLE PRECISION :: Ynew_tlm(N,NTLM), Fcn0_tlm(N,NTLM), Fcn_tlm(N,NTLM) 
   DOUBLE PRECISION :: K_tlm(N*ros_S,NTLM)
   DOUBLE PRECISION :: dFdT(N),Tmp(N),Delta
   DOUBLE PRECISION :: H, Hnew, HC, HG, Fac, Tau 
   DOUBLE PRECISION :: Err, Err0, Err1, Yerr(N), Yerr_tlm(N,NTLM)
   INTEGER  :: Pivot(N), Direction, ioffset, j, istage, itlm
   LOGICAL  :: RejectLastH, RejectMoreH, Singular, Transp
!~~~>  Local parameters
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-5
!~~~>  Locally called functions
!   DOUBLE PRECISION WLAMCH
!   EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   
!~~~>  Initial preparations
   T = Tstart
   RSTATUS(Nhexit) = ZERO
   H = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   IF (ABS(H) <= 10.D0*Roundoff) H = DeltaMin
   
   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF               
   H = Direction*H

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.
   Transp = .FALSE.
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
   Hexit = H
   H = MIN(H,ABS(Tend-T))

!~~~>   Compute the function at current time
   CALL Fun(N,T,Y,Fcn0)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>   Compute the Jacobian at current time
   CALL LSS_Jac(T,Y,JAC)
   ISTATUS(Njac) = ISTATUS(Njac) + 1
  
!~~~>   Compute the Hessian at current time
!   CALL HESS(N,T,Y,Hes0)
!   ISTATUS(Nhes) = ISTATUS(Nhes) + 1
   
!~~~>   Compute the TLM function at current time
   DO itlm = 1, NTLM
      CALL LSS_Mul_Jac(Fcn0_tlm(1,itlm),Y_TLM(1,itlm))
   END DO  
   
!~~~>  Compute the function and Jacobian derivatives with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, Fcn0, dFdT )
      Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
      CALL LSS_Jac2(T+Delta,Y,JAC)
      ISTATUS(Njac) = ISTATUS(Njac) + 1 
      CALL LSS_Jac_Time_Derivative (Delta)
   END IF
 
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

      ! Initialize stage solution
       CALL DCOPY(N,Y,1,Ynew,1)
       CALL DCOPY(N*NTLM,Y_tlm,1,Ynew_tlm,1)
      
      ! For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         CALL DCOPY(N,Fcn0,1,Fcn,1)
         CALL DCOPY(N*NTLM,Fcn0_tlm,1,Fcn_tlm,1)
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
         DO j = 1, istage-1
           CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j),    &
                     K(N*(j-1)+1),1,Ynew,1) 
           DO itlm=1,NTLM                    
              CALL DAXPY(N,ros_A((istage-1)*(istage-2)/2+j), &
                     K_tlm(N*(j-1)+1,itlm),1,Ynew_tlm(1,itlm),1) 
           END DO                    
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FUN(N,Tau,Ynew,Fcn)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
         CALL LSS_Jac1(Tau,Ynew,JAC)
         ISTATUS(Njac) = ISTATUS(Njac) + 1
         DO itlm=1,NTLM 
           CALL LSS_Mul_Jac1(Fcn_tlm(1,itlm),Ynew_tlm(1,itlm))
         END DO              
       END IF ! if istage == 1 elseif ros_NewF(istage)
       CALL DCOPY(N,Fcn,1,K(ioffset+1),1)
       DO itlm=1,NTLM   
          CALL DCOPY(N,Fcn_tlm(1,itlm),1,K_tlm(ioffset+1,itlm),1)
       END DO                
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL DAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
         DO itlm=1,NTLM 
           CALL DAXPY(N,HC,K_tlm(N*(j-1)+1,itlm),1,K_tlm(ioffset+1,itlm),1)
         END DO              
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL DAXPY(N,HG,dFdT,1,K(ioffset+1),1)
         DO itlm=1,NTLM 
           CALL LSS_Mul_Jac2(Tmp,Ynew_tlm(1,itlm))
           CALL DAXPY(N,HG,Tmp,1,K_tlm(ioffset+1,itlm),1)
         END DO              
       END IF
       CALL LSS_Solve(Transp,K(ioffset+1))
       ISTATUS(Nsol) = ISTATUS(Nsol) + 1
       DO itlm=1,NTLM
         CALL Hess_Vec(N,T,Y,K(ioffset+1),Y_tlm(1,itlm),Tmp)
!         CALL Hess_Vec ( Hes0, K(ioffset+1), Y_tlm(1,itlm), Tmp )
         CALL DAXPY(N,ONE,Tmp,1,K_tlm(ioffset+1,itlm),1)
         CALL LSS_Solve(Transp,K_tlm(ioffset+1,itlm))
         ISTATUS(Nsol) = ISTATUS(Nsol) + 1
       END DO                
      
   END DO Stage     
            

!~~~>  Compute the new solution 
   CALL DCOPY(N,Y,1,Ynew,1)
   DO j=1,ros_S
      CALL DAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
   END DO
   DO itlm=1,NTLM       
     CALL DCOPY(N,Y_tlm(1,itlm),1,Ynew_tlm(1,itlm),1)
     DO j=1,ros_S
       CALL DAXPY(N,ros_M(j),K_tlm(N*(j-1)+1,itlm),1,Ynew_tlm(1,itlm),1)
     END DO
   END DO

!~~~>  Compute the error estimation 
   CALL Set2zero(N,Yerr)
   DO j=1,ros_S     
        CALL DAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
   END DO 
   Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
   IF (TLMtruncErr) THEN
     Err1 = 0.0d0
     CALL Set2zero(N*NTLM,Yerr_tlm)
     DO itlm=1,NTLM
       DO j=1,ros_S  
         CALL DAXPY(N,ros_E(j),K_tlm(N*(j-1)+1,itlm),1,Yerr_tlm(1,itlm),1)
       END DO 
     END DO
     Err = ros_ErrorNorm_tlm(Y_tlm,Ynew_tlm,Yerr_tlm,AbsTol_tlm,RelTol_tlm,Err,VectorTol)
   END IF
  
!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac  

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      CALL DCOPY(N,Ynew,1,Y,1)
      CALL DCOPY(N*NTLM,Ynew_tlm,1,Y_tlm,1)
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

  END SUBROUTINE ros_TLM_Int
   
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DOUBLE PRECISION FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, & 
               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

! Input arguments   
   DOUBLE PRECISION, INTENT(IN) :: Y(N), Ynew(N),    &
          Yerr(N), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables     
   DOUBLE PRECISION :: Err, Scale, Ymax   
   INTEGER  :: i
   DOUBLE PRECISION, PARAMETER :: ZERO = 0.0d0
   
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
  DOUBLE PRECISION FUNCTION ros_ErrorNorm_tlm ( Y_tlm, Ynew_tlm, Yerr_tlm, & 
               AbsTol_tlm, RelTol_tlm, Fwd_Err, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr_tlm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE         

! Input arguments   
   DOUBLE PRECISION, INTENT(IN) :: Y_tlm(N,NTLM), Ynew_tlm(N,NTLM),    &
          Yerr_tlm(N,NTLM), AbsTol_tlm(N,NTLM), RelTol_tlm(N,NTLM), Fwd_Err
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables     
   DOUBLE PRECISION :: TMP, Err
   INTEGER  :: itlm
   
   Err = FWD_Err
   DO itlm = 1,NTLM
     TMP = ros_ErrorNorm(Y_tlm(1,itlm), Ynew_tlm(1,itlm),Yerr_tlm(1,itlm), &
     		AbsTol_tlm(1,itlm), RelTol_tlm(1,itlm), VectorTol)
     Err = MAX(Err, TMP)
   END DO

   ros_ErrorNorm_tlm = MAX(Err,1.0d-10)
   
  END FUNCTION ros_ErrorNorm_tlm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
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
   DOUBLE PRECISION, PARAMETER :: DeltaMin = 1.0d-6
   
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL FUN(N,T+Delta,Y,dFdT)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   CALL DAXPY(N,(-ONE),Fcn0,1,dFdT,1)
   CALL DSCAL(N,(ONE/Delta),dFdT,1)

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   USE LS_Solver
   IMPLICIT NONE         
   
!~~~> Input arguments   
   DOUBLE PRECISION, INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments   
   LOGICAL, INTENT(OUT) ::  Singular
!~~~> Inout arguments   
   DOUBLE PRECISION, INTENT(INOUT) :: H   ! step size is decreased when LU fails
!~~~> Local variables     
   INTEGER  :: i, ising, Nconsecutive
   DOUBLE PRECISION ::  ghinv
   DOUBLE PRECISION, PARAMETER :: ONE  = 1.0d0, HALF = 0.5d0
   
   Nconsecutive = 0
   Singular = .TRUE.
   
   DO WHILE (Singular)
   
!~~~>    Construct Ghimj = 1/(H*ham) - Jac0
     ghinv = ONE/(Direction*H*gam)
!~~~>    Compute LU decomposition 
     CALL LSS_Decomp(ghinv,ising)
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
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

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


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3 
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

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
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

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
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
END SUBROUTINE RosenbrockTLM
!  and all its internal procedures
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

END MODULE ROS_tlm_f90_Integrator




! End of INTEGRATE function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


