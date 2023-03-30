MODULE User_Parameters
    INTEGER, PARAMETER :: NVAR = 5
    INTEGER, PARAMETER :: NNZ = 19
    INTEGER, PARAMETER :: NREACT=10 , NFIX=2
! RCONST - Rate constants (local)
    DOUBLE PRECISION :: RCONST(NREACT)
   ! FIX - Concentrations of fixed species 
    DOUBLE PRECISION :: FIX(NFIX)
CONTAINS
    SUBROUTINE FUN_U(NVAR,T, V, FCT)
    IMPLICIT NONE
    INTEGER :: NVAR
   ! V - Concentrations of variable species (local)
   ! FCT - Time derivative of variable species concentrations
    DOUBLE PRECISION :: A(NREACT)
    DOUBLE PRECISION :: V(NVAR), FCT(NVAR)
    DOUBLE PRECISION :: T, SUN
    CALL Update_SUN(T,SUN)
    CALL Update_RCONST(SUN)
  ! Computation of equation rates
    A(1) = RCONST(1)*FIX(2)
    A(2) = RCONST(2)*V(2)*FIX(2)
    A(3) = RCONST(3)*V(3)
    A(4) = RCONST(4)*V(2)*V(3)
    A(5) = RCONST(5)*V(3)
    A(6) = RCONST(6)*V(1)*FIX(1)
    A(7) = RCONST(7)*V(1)*V(3)
    A(8) = RCONST(8)*V(3)*V(4)
    A(9) = RCONST(9)*V(2)*V(5)
    A(10) = RCONST(10)*V(5)
  ! Aggregate function
    FCT(1) = A(5)-A(6)-A(7)
    FCT(2) = 2*A(1)-A(2)+A(3)-A(4)+A(6)-A(9)+A(10)
    FCT(3) = A(2)-A(3)-A(4)-A(5)-A(7)-A(8)
    FCT(4) = -A(8)+A(9)+A(10)
    FCT(5) = A(8)-A(9)-A(10)
    END SUBROUTINE FUN_U

    SUBROUTINE JAC_U(NVAR,T, V, JF)
    IMPLICIT NONE
    INTEGER :: NVAR
! V - Concentrations of variable species (local)
    DOUBLE PRECISION :: V(NVAR), T, SUN
    DOUBLE PRECISION :: JF(NVAR,NVAR)
! B - Temporary array
    DOUBLE PRECISION :: B(16)
    CALL Update_SUN(T, SUN)
    CALL Update_RCONST(SUN)
    JF(:,:) = 0.0d0

    B(2) = RCONST(2)*FIX(2)
    B(4) = RCONST(3)
    B(5) = RCONST(4)*V(3)
    B(6) = RCONST(4)*V(2)
    B(7) = RCONST(5)
    B(8) = RCONST(6)*FIX(1)
    B(10) = RCONST(7)*V(3)
    B(11) = RCONST(7)*V(1)
    B(12) = RCONST(8)*V(4)
    B(13) = RCONST(8)*V(3)
    B(14) = RCONST(9)*V(5)
    B(15) = RCONST(9)*V(2)
    B(16) = RCONST(10)
! Construct the Jacobian terms from B's
    JF(1,1) = -B(8)-B(10)
    JF(1,3) = B(7)-B(11)
    JF(2,1) = B(8)
    JF(2,2) = -B(2)-B(5)-B(14)
    JF(2,3) = B(4)-B(6)
    JF(2,5) = -B(15)+B(16)
    JF(3,1) = -B(10)
    JF(3,2) = B(2)-B(5)
    JF(3,3) = -B(4)-B(6)-B(7)-B(11)-B(12)
    JF(3,4) = -B(13)
    JF(3,5) = 0
    JF(4,2) = B(14)
    JF(4,3) = -B(12)
    JF(4,4) = -B(13)
    JF(4,5) = B(15)+B(16)
    JF(5,2) = -B(14)
    JF(5,3) = B(12)
    JF(5,4) = B(13)
    JF(5,5) = -B(15)-B(16)
    END SUBROUTINE JAC_U

    SUBROUTINE Update_SUN(TIME, SUN)
    IMPLICIT NONE
    DOUBLE PRECISION :: TIME,SUN
    DOUBLE PRECISION :: SunRise, SunSet
    DOUBLE PRECISION :: Thour, Tlocal, Ttmp
    ! PI - Value of pi
    DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979d0

    SunRise = 4.5
    SunSet  = 19.5
    Thour = TIME/3600.0
    Tlocal = Thour - (INT(Thour)/24)*24

    IF ((Tlocal>=SunRise).AND.(Tlocal<=SunSet)) THEN
       Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise)
       IF (Ttmp.GT.0) THEN
          Ttmp =  Ttmp*Ttmp
       ELSE
          Ttmp = -Ttmp*Ttmp
       END IF
       SUN = ( 1.0 + COS(PI*Ttmp) )/2.0
    ELSE
       SUN = 0.0
    END IF

    END SUBROUTINE Update_SUN

    SUBROUTINE Update_RCONST(SUN)
    IMPLICIT NONE
    DOUBLE PRECISION :: SUN
    RCONST(1) = ((2.643E-10)*SUN*SUN*SUN)
    RCONST(2) = ((8.018E-17))
    RCONST(3) = ((6.120E-04)*SUN)
    RCONST(4) = ((1.576E-15))
    RCONST(5) = ((1.070E-03)*SUN*SUN)
    RCONST(6) = ((7.110E-11))
    RCONST(7) = ((1.200E-10))
    RCONST(8) = ((6.062E-15))
    RCONST(9) = ((1.069E-11))
    RCONST(10) = ((1.289E-02)*SUN)
    END SUBROUTINE Update_RCONST

    SUBROUTINE Update_PHOTO (SUN)
    DOUBLE PRECISION :: SUN
    RCONST(1) = ((2.643E-10)*SUN*SUN*SUN)
    RCONST(3) = ((6.120E-04)*SUN)
    RCONST(5) = ((1.070E-03)*SUN*SUN)
    RCONST(10) = ((1.289E-02)*SUN)
    END SUBROUTINE Update_PHOTO

END MODULE User_Parameters
