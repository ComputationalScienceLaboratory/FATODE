!shallow water equations

module swe2Dxy_Parameters

  IMPLICIT NONE
  public

  save

  !gravitational acceleration [m/s^2]
  double precision, parameter :: gAcc = 9.80616d0

  !domain boundaries [unit square]
  double precision, parameter ::  xMin = -3.d0
  double precision, parameter ::  xMax = 3.d0

  double precision, parameter ::  yMin = -3.d0
  double precision, parameter ::  yMax = 3.d0

  double precision            :: hMax, speed

  !number of grid points along each direction
  integer, parameter ::  Mx = 40
  integer, parameter ::  My = 40

  integer, parameter ::  cells = Mx*My
  integer, parameter ::  ndim = 3*cells
  integer, parameter ::  nnz   = 63*cells

  !grid spacings
  double precision, parameter :: dx = (xMax - xMin) / (Mx - 1d0)
  double precision, parameter :: dy = (yMax - yMin) / (My - 1d0)
 
  !number of time steps
  integer :: NTS = 50

  !CFL number
  double precision :: CFL = 0.5d0

  !time step length
  double precision :: dt

  !~~> Text file where we store the solution
  !    The solution U will be written at every time step
  integer :: heightFile = 999   !write the height here
  integer :: uveloFile  = 998   !write the u velocity component here
  integer :: vveloFile  = 997   !write the v velocity component here

  !~~> Observation frequency (i.e. observe state every so many time steps)
  integer, parameter   :: obsFrequency = 5
  !~~> Size of the observation space (i.e. observe this many grid points)
  !integer, parameter   :: obsSpaceSize = ndim
  integer, parameter   :: obsSpaceSize = 300

  !--> Rinv - diagonal of the inverse of the observation error cov matrix
  double precision                      :: Rinv(obsSpaceSize)
  !--> B - background error covariance matrix, Lchol - its Cholesky factor
  double precision, allocatable         :: Bcov(:,:), Lchol(:,:)
  !--> Zeta - uniformly distributed perturbation factors
  double precision                      :: Zeta(30, ndim)
  !~~> Current initial state and background state
  double precision                      :: back(ndim), ref(ndim), ic(ndim)

  integer                               :: genfiles

  character*30        :: itercostfile, iterrmsefile, itertimefile
  character*30        :: runscostfile, runsrmsefile, runstimefile
  character           :: datamodestr
  character*30        :: AssimFile, runsgradfile, runshessfile

  double precision, allocatable         :: JF(:,:), JG(:,:)
! for convience add new definitions
  double precision  :: U_par(3,Mx+4,My+4)
  double precision, allocatable    :: F_par(:,:,:), FZ_par(:,:,:)
  double precision, allocatable    :: temp_par(:,:,:)
  double precision, allocatable :: feigvalues_par(:,:,:), geigvalues_par(:,:,:)

  integer              :: reclen

  !--> Change this to > 0 for debugging purposes
  integer, parameter   :: verbose = 1

end module swe2Dxy_Parameters


subroutine initializeGaussianGrid(grid, A)

  !set up 2D height grid with a Gaussian pulse of amplitude A centered at
  !(0,0)
  use swe2Dxy_parameters

  implicit none

  double precision, intent(inout) :: grid(3,Mx+4,My+4)
  double precision, intent(in) :: A

  integer i,j
  double precision :: x(Mx+4), y(My+4)

  !x = linspace(xMin,xMax,Mx)
  !y = linspace(yMin,yMax,My)
  x(:) = 0d0
  y(:) = 0d0

  do i=3,Mx+2
    x(i) = xMin + (i-2)*dx
  end do

  do j=3,My+2
    y(j) = yMin + (j-2)*dy
  end do

  !zero x- and y-velocities
  grid = 0d0
  do j=3, My+2
    do i=3, Mx+2
      grid(1,i,j) =  A * exp(-x(j)*x(j) - y(i)*y(i))
    end do
  end do

  !make sure height is positive everywhere
  grid(1,:,:) = grid(1,:,:) + 100d0
  grid(2,:,:) = 0d0
  grid(3,:,:) = 0d0

end subroutine initializeGaussianGrid

subroutine mapgr(ind, p, k, i, j)

  use swe2Dxy_Parameters
  implicit none
  integer :: ind
  integer :: p, k, i, j
  integer :: x, y, z

  call mapgrinv(p, z, x, y)

  z=z+k
  x=x+i
  if (x .EQ. 2) x=Mx+2
  if (x .EQ. 1) x=Mx+1
  if (x .EQ. (Mx+3)) x=3
  if (x .EQ. (Mx+4)) x=4
  y=y+j
  if (y .EQ. 2) y=My+2
  if (y .EQ. 1) y=My+1
  if (y .EQ. (My+3)) y=3
  if (y .EQ. (My+4)) y=4
  
  ind = (z-1)*Mx*My+(x-3)*Mx+y-2
  
end subroutine

subroutine mapgrinv(p, k, i, j)

  use swe2Dxy_Parameters
  implicit none
  integer  :: p
  integer :: k, i, j
  integer :: aux

  k = (p-1)/(Mx*My)
  aux = mod(p-1,Mx*My)

  i = aux/Mx
  j = mod(aux,Mx)

  k=k+1
  i=i+1
  j=j+1
  
! take boundaries into account
  i=i+2
  j=j+2

end subroutine


! (3,Nx+4,Ny+4) tensor to (3*Nx*Ny) vector mapping operator
subroutine grid2vec(grid_in, vec_out)

  use swe2Dxy_parameters
  implicit none

  double precision, intent(in)  ::  grid_in(3,Mx+4,My+4)
  double precision, intent(out) ::  vec_out(3*Mx*My)

  integer :: k, i, j
  integer :: indx

  vec_out(:) = 0d0
  indx = 1
  do k = 1, 3
    do i = 3, Mx+2
      do j = 3, My+2
           vec_out(indx) = grid_in(k,i,j)
           indx = indx + 1
      end do
    end do
  end do
  
end subroutine grid2vec

! (3*Nx*Ny) vector to (3,Nx+4,Ny+4) tensor mapping operator
subroutine vec2grid(vec_in, grid_out)

  use swe2Dxy_parameters
  implicit none

  double precision, intent(out) ::  grid_out(3,Mx+4,My+4)
  double precision, intent(in)  ::  vec_in(3*Mx*My)

  integer :: k, i, j
  integer :: indx

  grid_out = 0d0
  indx = 1
  do k = 1, 3
    do i = 3, Mx+2
      do j = 3, My+2
           grid_out(k,i,j) = vec_in(indx)
           indx = indx + 1
      end do
    end do
  end do

  
end subroutine vec2grid


subroutine compute_F(U, F, feigvalues, geigvalues)

  use swe2Dxy_Parameters
  implicit none
  double precision, intent(inout) :: U(3,Mx+4,My+4)
  double precision, intent(out)   :: F(3,Mx+4,My+4)
  double precision, intent(out)   :: feigvalues(3,Mx+4,My+4), geigvalues(3,Mx+4,My+4)

  !~~> local variables
  double precision  :: dUdxMinus(3), dUdxPlus(3)
  double precision  :: dUdyMinus(3), dUdyPlus(3)
  double precision  :: rhs(3)
  double precision  :: eigf(3), eigg(3)
  double precision  :: Fu1(3,3), Fu2(3,3), Fu3(3,3)
  double precision  :: Gu1(3,3), Gu2(3,3), Gu3(3,3)
  double precision  :: tmpF(3), tmpG(3), tmpMV(3)

  integer           :: i, j,  i1, j1 
  double precision  :: hc, uc, vc, g, isqrt


  !take care of boundary conditions
  !periodic boundaries: initialize halo cells (halo 2x2 corners not needed)
  !y-axis
  U(:,2,3:My+2)    =  U(:,Mx+2,3:My+2)
  U(:,1,3:My+2)    =  U(:,Mx+1,3:My+2)
  U(:,Mx+3,3:My+2) =  U(:,3,3:My+2)
  U(:,Mx+4,3:My+2) =  U(:,4,3:My+2)
  !x-axis
  U(:,3:Mx+2,2)    =  U(:,3:Mx+2,My+2)
  U(:,3:Mx+2,1)    =  U(:,3:Mx+2,My+1)
  U(:,3:Mx+2,My+3) =  U(:,3:Mx+2,3)
  U(:,3:Mx+2,My+4) =  U(:,3:Mx+2,4)

  feigvalues(:,:,:)   =  0d0
  geigvalues(:,:,:)   =  0d0

! space discretization for each gridcell

  do j = 3,My+2
     do i = 3,Mx+2

              !split the fluxes according to the signs of the eigenvalues
              !of the Jacobian matrices evaluated at U(:,i,j)
              eigf(1) = U(2,i,j)/U(1,i,j)                       !u
              eigf(2) = U(2,i,j)/U(1,i,j) - sqrt(gAcc*U(1,i,j)) !u - sqrt(gh)
              eigf(3) = U(2,i,j)/U(1,i,j) + sqrt(gAcc*U(1,i,j)) !u + sqrt(gh)
              feigvalues(:,i,j)=eigf

              hc = U(1,i,j)
              uc = U(2,i,j) / U(1,i,j)
              vc = U(3,i,j) / U(1,i,j)
              g = gAcc

              Fu1(:,:) = 0d0
              Fu1(3,1) = -vc
              Fu1(3,3) = 1d0

              isqrt = 1d0 / sqrt(g*hc)
              Fu2(:,:) = 0d0
              Fu2(1,1) = 0.5d0 * (1d0 + uc * isqrt )
              Fu2(1,2) = - 0.5d0 * isqrt
              Fu2(2,1) = 0.5d0 * (uc*uc - g*hc) * isqrt
              Fu2(2,2) = 0.5d0 - 0.5d0 * uc * isqrt
              Fu2(3,1) = 0.5d0 * (1d0 + uc * isqrt) * vc
              Fu2(3,2) = - 0.5d0 * vc * isqrt

              isqrt = 1d0 / sqrt(g*hc)
              Fu3(:,:) = 0d0
              Fu3(1,1) = 0.5d0 - 0.5d0 * uc * isqrt
              Fu3(1,2) = 0.5d0 * isqrt
              Fu3(2,1) = 0.5d0 * (g*hc - uc*uc) * isqrt
              Fu3(2,2) = 0.5d0 * (1d0 + uc * isqrt)
              Fu3(3,1) = 0.5d0 * (1d0 - uc * isqrt) * vc
              Fu3(3,2) = 0.5d0 * vc * isqrt

! use 3rd order upwind finite differences

              dUdxMinus = (-U(:,i+2,j) + 6d0*U(:,i+1,j) - 3d0*U(:,i,j) - 2d0*U(:,i-1,j)) / (6d0 * dx)
              dUdxPlus = (2d0*U(:,i+1,j) + 3d0*U(:,i,j) - 6d0*U(:,i-1,j) + U(:,i-2,j)) / (6.d0 * dx)
              !dUdxMinus = (-U(:,i+2,j) + 4d0*U(:,i+1,j) - 3d0*U(:,i,j)) / (2d0*dx)
              !dUdxPlus = (3d0*U(:,i,j) - 4d0*U(:,i-1,j) + U(:,i-2,j)) / (2d0*dx)
              !dUdxMinus = (U(:,i+1,j) - U(:,i,j)) / dx
              !dUdxPlus  = (U(:,i,j) - U(:,i-1,j)) / dx

              if (eigf(1) .ge. 0d0) then
                     if (verbose .EQ. 2) print *, "Compute F1+ for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Fu1(i1,j1) * dUdxPlus(j1)
                       end do
                     end do
                     tmpF = eigf(1) * tmpMV
              else
                     if (verbose .EQ. 2) print *, "Compute F1- for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Fu1(i1,j1) * dUdxMinus(j1)
                       end do
                     end do
                     tmpF = eigf(1) * tmpMV
              end if

              if (eigf(2) .ge. 0d0) then
                     if (verbose .EQ. 2) print *, "Compute F2+ for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Fu2(i1,j1) * dUdxPlus(j1)
                       end do
                     end do
                     tmpF = tmpF + eigf(2) * tmpMV
              else
                     if (verbose .EQ. 2) print *, "Compute F2- for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Fu2(i1,j1) * dUdxMinus(j1)
                       end do
                     end do
                     tmpF = tmpF + eigf(2) * tmpMV
              end if

              if (eigf(3) .ge. 0d0) then
                     if (verbose .EQ. 2) print *, "Compute F3+ for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Fu3(i1,j1) * dUdxPlus(j1)
                       end do
                     end do
                     tmpF = tmpF + eigf(3) * tmpMV
              else
                     if (verbose .EQ. 2) print *, "Compute F3- for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Fu3(i1,j1) * dUdxMinus(j1)
                       end do
                     end do
                     tmpF = tmpF + eigf(3) * tmpMV
              end if

              eigg(1) = U(3,i,j)/U(1,i,j)                       !v
              eigg(2) = U(3,i,j)/U(1,i,j) - sqrt(gAcc*U(1,i,j)) !v - sqrt(gh)
              eigg(3) = U(3,i,j)/U(1,i,j) + sqrt(gAcc*U(1,i,j)) !v + sqrt(gh)
              geigvalues(:,i,j)=eigg

              Gu1(:,:) = 0d0
              Gu1(2,1) = -uc
              Gu1(2,2) = 1d0

              isqrt = 1d0 / sqrt(g*hc)
              Gu2(:,:) = 0d0
              Gu2(1,1) = 0.5d0 * (1d0 + vc * isqrt)
              Gu2(1,3) = - 0.5d0 * isqrt
              Gu2(2,1) = 0.5d0 * uc * (1d0 + vc * isqrt)
              Gu2(2,3) = - 0.5d0 * uc * isqrt
              Gu2(3,1) = 0.5d0 * (vc*vc - g*hc) * isqrt
              Gu2(3,3) = 0.5d0 - 0.5d0 * vc * isqrt

              isqrt = 1d0 / sqrt(g*hc)
              Gu3(:,:) = 0d0
              Gu3(1,1) = 0.5d0 * (1d0 - vc * isqrt)
              Gu3(1,3) = 0.5d0 * isqrt
              Gu3(2,1) = 0.5d0 * uc * (1d0 - vc * isqrt)
              Gu3(2,3) = 0.5d0 * uc * isqrt
              Gu3(3,1) = 0.5d0 * (-vc*vc + g*hc) * isqrt
              Gu3(3,3) = 0.5d0 + 0.5d0 * vc  * isqrt

! use 3rd order upwind finite differences

              dUdyMinus = (-U(:,i,j+2) + 6d0*U(:,i,j+1) - 3d0*U(:,i,j) - 2d0*U(:,i,j-1)) / (6d0 * dy)
              dUdyPlus = (2d0*U(:,i,j+1) + 3d0*U(:,i,j) - 6d0*U(:,i,j-1) + U(:,i,j-2)) / (6.d0 * dy)
              !dUdyMinus = (-U(:,i,j+2) + 4d0*U(:,i,j+1) - 3d0*U(:,i,j)) / (2d0*dy)
              !dUdyPlus = (3d0*U(:,i,j) - 4d0*U(:,i,j-1) + U(:,i,j-2)) / (2d0*dy)
              !dUdyMinus = (U(:,i,j+1) - U(:,i,j)) / dy
              !dUdyPlus  = (U(:,i,j) - U(:,i,j-1)) / dy

              if (eigg(1) .ge. 0d0) then
                     if (verbose .EQ. 2) print *, "Compute G1+ for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Gu1(i1,j1) * dUdyPlus(j1)
                       end do
                     end do
                     tmpG = eigg(1) * tmpMV
              else
                     if (verbose .EQ. 2) print *, "Compute G1- for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Gu1(i1,j1) * dUdyMinus(j1)
                       end do
                     end do
                     tmpG = eigg(1) * tmpMV
              end if

              if (eigg(2) .ge. 0d0) then
                     if (verbose .EQ. 2) print *, "Compute G2+ for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Gu2(i1,j1) * dUdyPlus(j1)
                       end do
                     end do
                     tmpG = tmpG + eigg(2) * tmpMV
              else
                     if (verbose .EQ. 2) print *, "Compute G2- for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Gu2(i1,j1) * dUdyMinus(j1)
                       end do
                     end do
                     tmpG = tmpG + eigg(2) * tmpMV
              end if

              if (eigg(3) .ge. 0d0) then
                     if (verbose .EQ. 2) print *, "Compute G3+ for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Gu3(i1,j1) * dUdyPlus(j1)
                       end do
                     end do
                     tmpG = tmpG + eigg(3) * tmpMV
              else
                     if (verbose .EQ. 2) print *, "Compute G3- for ", i, j
                     tmpMV(:) = 0d0
                     do i1=1,3
                       do j1=1,3
                         tmpMV(i1) = tmpMV(i1) + Gu3(i1,j1) * dUdyMinus(j1)
                       end do
                     end do
                     tmpG = tmpG + eigg(3) * tmpMV
              end if

              rhs = -tmpF - tmpG

              F(:,i,j) = rhs

    end do
  end do


end subroutine


subroutine compute_JacF(U, feigvalues, geigvalues)

  use swe2Dxy_Parameters
  implicit none
  
  double precision, intent(in)  :: U(3, Mx+4, My+4)
  double precision, intent(in)  :: feigvalues(3,Mx+4,My+4), geigvalues(3,Mx+4,My+4)

  double precision  :: NN1, N1, C1, S1, SS1, WW1, W1, E1, EE1
  double precision  :: NN2, N2, C2, S2, SS2, WW2, W2, E2, EE2
  double precision  :: NN3, N3, C3, S3, SS3, WW3, W3, E3, EE3
  double precision  :: g, SqrtC1g
!redeclare
!  double precision  :: temp(ndim)

  integer           :: p, k, i, j, ind 


  g = gAcc

  JF(:,:)=0d0
  JG(:,:)=0d0

do p = 1, Mx*My

  call mapgrinv(p,k,i,j)

  NN1 = U(k,i-2,j)
  N1  = U(k,i-1,j)
  C1  = U(k,i,j)
  S1  = U(k,i+1,j)
  SS1 = U(k,i+2,j)

  NN2 = U(k+1,i-2,j)
  N2  = U(k+1,i-1,j)
  C2  = U(k+1,i,j)
  S2  = U(k+1,i+1,j)
  SS2 = U(k+1,i+2,j)

  NN3 = U(k+2,i-2,j)
  N3  = U(k+2,i-1,j)
  C3  = U(k+2,i,j)
  S3  = U(k+2,i+1,j)
  SS3 = U(k+2,i+2,j)

  SqrtC1g = Sqrt(C1*g)

  if (feigvalues(1,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JF1+ for ", i, j

    JF(p,:)=0d0

  else

    if (verbose .EQ. 2) print *, "Compute JF1- for ", i, j

    JF(p,:)=0d0

  end if

  if (feigvalues(2,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JF2+ for ", i, j

    call mapgr(ind,p,0,-2,0) !NN1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5+(0.5*C2)/(C1*SqrtC1g))*(C2/C1-SqrtC1g))/(6.*dx)
    end if
    call mapgr(ind,p,0,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(((0.5+(0.5*C2)/(C1*SqrtC1g))*(C2/C1-SqrtC1g))/dx)
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.5 + (0.5*C2)/(C1*SqrtC1g))/(2.*dx) + (((-0.25*C2*g)/(C1*(C1*g)**1.5) - &
      & (0.5*C2)/(C1**2*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) +&
      & (0.0416667*g*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*(C1*g)**1.5)) + (-(C2/C1**2) - g/(2.*SqrtC1g))*(((0.5 + &
      & (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) - (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*SqrtC1g))
    end if  
    call mapgr(ind,p,0,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5+(0.5*C2)/(C1*SqrtC1g))*(C2/C1-SqrtC1g))/(3.*dx)
    end if  
    call mapgr(ind,p,1,-2,0) !NN2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.0833333*(C2/C1-SqrtC1g))/(dx*SqrtC1g))
    end if  
    call mapgr(ind,p,1,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*(C2/C1-SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*(-0.25/(dx*SqrtC1g) + (0.0833333*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx*SqrtC1g)) + &
      & (((0.5 + (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) - (0.0833333*(3*C2 - 6*N2 + NN2 + &
      & 2*S2))/(dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.166667*(C2/C1-SqrtC1g))/(dx*SqrtC1g))
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JF2- for ", i, j

    call mapgr(ind,p,0,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(((0.5d0 + (0.5d0*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(3d0*dx))
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*(-(0.5 + (0.5*C2)/(C1*SqrtC1g))/(2.*dx) +&
      & (((-0.25*C2*g)/(C1*(C1*g)**1.5) - (0.5*C2)/(C1**2*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*dx) +&
      & (0.0416667*g*(-3*C2 - 2*N2 + 6*S2 - SS2))/(dx*(C1*g)**1.5)) + (-(C2/C1**2) -&
      & g/(2.*SqrtC1g))*(((0.5 + (0.5*C2)/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*dx) -&
      & (0.0833333*(-3*C2 - 2*N2 + 6*S2 - SS2))/(dx*SqrtC1g))
    end if  
    call mapgr(ind,p,0,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/dx
    end if  
    call mapgr(ind,p,0,2,0)  !SS1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(6.*dx)
    end if  
    call mapgr(ind,p,1,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*(C2/C1 - SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+ (C2/C1 - SqrtC1g)*(0.25/(dx*SqrtC1g) +&
      & (0.0833333*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx*SqrtC1g)) +&
      & (((0.5 + (0.5*C2)/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*dx) -&
      & (0.0833333*(-3*C2 - 2*N2 + 6*S2 - SS2))/(dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.5*(C2/C1 - SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,1,2,0)  !SS2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*(C2/C1 - SqrtC1g))/(dx*SqrtC1g)
    end if  

  end if

  if (feigvalues(3,i,j) .GE. 0d0) then
  
    if (verbose .EQ. 2) print *, "Compute JF3+ for ", i, j

    call mapgr(ind,p,0,-2,0) !NN1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(6.*dx)
    end if  
    call mapgr(ind,p,0,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/dx)
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.5 - (0.5*C2)/(C1*SqrtC1g))/(2.*dx) +&
      & (((0.25*C2*g)/(C1*(C1*g)**1.5) + (0.5*C2)/(C1**2*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) -&
      & (0.0416667*g*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*(C1*g)**1.5)) + (-(C2/C1**2) + g/(2.*SqrtC1g))*(((0.5 - &
      & (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) + (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*SqrtC1g))
    end if  
    call mapgr(ind,p,0,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(3.*dx)
    end if  
    call mapgr(ind,p,1,-2,0) !NN2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,1,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.5*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*(0.25/(dx*SqrtC1g) - (0.0833333*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx*SqrtC1g)) + (((0.5 &
      & - (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) + (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JF3- for ", i, j

    call mapgr(ind,p,0,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(3.*dx)
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*(-(0.5 - (0.5*C2)/(C1*SqrtC1g))/(2.*dx) + &
      & (((0.25*C2*g)/(C1*(C1*g)**1.5) + (0.5*C2)/(C1**2*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*dx) -&
      & (0.0416667*g*(-3*C2 - 2*N2 + 6*S2 - SS2))/(dx*(C1*g)**1.5)) + (-(C2/C1**2) + &
      & g/(2.*SqrtC1g))*(((0.5 - (0.5*C2)/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*dx) + (0.0833333*(-3*C2 - 2*N2 + 6*S2 - &
      & SS2))/(dx*SqrtC1g))
    end if  
    call mapgr(ind,p,0,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/dx
    end if  
    call mapgr(ind,p,0,2,0)  !SS1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(6.*dx)
    end if  
    call mapgr(ind,p,1,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.166667*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*(-0.25/(dx*SqrtC1g) - (0.0833333*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx*SqrtC1g)) +&
      & (((0.5 - (0.5*C2)/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*dx) + (0.0833333*(-3*C2 - 2*N2 + 6*S2 - &
      & SS2))/(dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,1,2,0)  !SS2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.0833333*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  

  end if

!goto 777
!776	continue
  
  WW1 = U(k,i,j-2)
  W1  = U(k,i,j-1)
  C1  = U(k,i,j)
  E1  = U(k,i,j+1)
  EE1 = U(k,i,j+2)

  WW2 = U(k+1,i,j-2)
  W2  = U(k+1,i,j-1)
  C2  = U(k+1,i,j)
  E2  = U(k+1,i,j+1)
  EE2 = U(k+1,i,j+2)  
  
  WW3 = U(k+2,i,j-2)
  W3  = U(k+2,i,j-1)
  C3  = U(k+2,i,j)
  E3  = U(k+2,i,j+1)
  EE3 = U(k+2,i,j+2)

  SqrtC1g = Sqrt(C1*g)
      
  JG(p,:)=0d0

  if (geigvalues(2,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JG2+ for ", i, j

    call mapgr(ind,p,0,0,-2) !WW1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(6.*dy)
    end if  
    call mapgr(ind,p,0,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-(((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/dy)
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.5 + (0.5*C3)/(C1*SqrtC1g))/(2.*dy) +&
      & (((-0.25*C3*g)/(C1*(C1*g)**1.5) - (0.5*C3)/(C1**2*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*dy) +&
      & (0.0416667*g*(3*C3 + 2*E3 - 6*W3 + WW3))/(dy*(C1*g)**1.5)) + (-(C3/C1**2) - g/(2.*SqrtC1g))*(((0.5 +&
      & (0.5*C3)/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*dy) - (0.0833333*(3*C3 + 2*E3 - 6*W3 + WW3))/(dy*SqrtC1g))
    end if  
    call mapgr(ind,p,0,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(3.*dy)
    end if  
    call mapgr(ind,p,2,0,-2) !WW3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,-1) !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,0)  !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*(-0.25/(dy*SqrtC1g) + (0.0833333*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy*SqrtC1g)) +&
      & (((0.5 + (0.5*C3)/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*dy) - (0.0833333*(3*C3 + 2*E3 - 6*W3 + &
      & WW3))/(dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,2,0,1)  !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JG2- for ", i, j

    call mapgr(ind,p,0,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(3.*dy)
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*(-(0.5 + (0.5*C3)/(C1*SqrtC1g))/(2.*dy) + (((-0.25*C3*g)/(C1*(C1*g)**1.5) -&
      & (0.5*C3)/(C1**2*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) + (0.0416667*g*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(dy*(C1*g)**1.5)) +&
      & (-(C3/C1**2) - g/(2.*SqrtC1g))*(((0.5 + (0.5*C3)/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) -&
      & (0.0833333*(-3*C3 + 6*E3 - EE3 - 2*W3))/(dy*SqrtC1g))
    end if  
    call mapgr(ind,p,0,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/dy
    end if  
    call mapgr(ind,p,0,0,2)  !EE1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(6.*dy)
    end if  
    call mapgr(ind,p,2,0,-1) !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,0)  !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*(0.25/(dy*SqrtC1g) + (0.0833333*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy*SqrtC1g)) +&
      & (((0.5 + (0.5*C3)/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) - (0.0833333*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,2,0,1)  !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,2)  !EE3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  

  end if  
  
  if (geigvalues(3,i,j) .GE. 0d0) then
  
    if (verbose .EQ. 2) print *, "Compute JG3+ for ", i, j
  
    call mapgr(ind,p,0,0,-2) !WW1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(6.*dy)
    end if  
    call mapgr(ind,p,0,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-(((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/dy)
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.5 - (0.5*C3)/(C1*SqrtC1g))/(2.*dy) + (((0.25*C3*g)/(C1*(C1*g)**1.5) +&
      & (0.5*C3)/(C1**2*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*dy) - (0.0416667*g*(3*C3 + 2*E3 - 6*W3 + &
      & WW3))/(dy*(C1*g)**1.5)) +&
      & (-(C3/C1**2) + g/(2.*SqrtC1g))*(((0.5 - (0.5*C3)/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*dy) +&
      & (0.0833333*(3*C3 + 2*E3 - 6*W3 + WW3))/(dy*SqrtC1g))
    end if  
    call mapgr(ind,p,0,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(3.*dy)
    end if  
    call mapgr(ind,p,2,0,-2) !WW3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,-1) !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,0)  !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*(0.25/(dy*SqrtC1g) - (0.0833333*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy*SqrtC1g)) +&
      & (((0.5 - (0.5*C3)/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*dy) + (0.0833333*(3*C3 + 2*E3 - 6*W3 + &
      & WW3))/(dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,2,0,1)  !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JG3- for ", i, j
  
    call mapgr(ind,p,0,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(3.*dy)
    end if  
    call mapgr(ind,p,0,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*(-(0.5 - (0.5*C3)/(C1*SqrtC1g))/(2.*dy) + (((0.25*C3*g)/(C1*(C1*g)**1.5) +&
      & (0.5*C3)/(C1**2*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) - (0.0416667*g*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(dy*(C1*g)**1.5)) +&
      & (-(C3/C1**2) + g/(2.*SqrtC1g))*(((0.5 - (0.5*C3)/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) + (0.0833333*(-3*C3 + &
      & 6*E3 - EE3 - 2*W3))/(dy*SqrtC1g))
    end if  
    call mapgr(ind,p,0,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/dy
    end if  
    call mapgr(ind,p,0,0,2)  !EE1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(6.*dy)
    end if  
    call mapgr(ind,p,2,0,-1) !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,0)  !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*(-0.25/(dy*SqrtC1g) - (0.0833333*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy*SqrtC1g)) +&
      & (((0.5 - (0.5*C3)/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) + (0.0833333*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,2,0,1)  !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,2,0,2)  !EE3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  

  end if

end do

!777	continue
!JF(:,:)=-JF(:,:)-JG(:,:)
!stop
!goto 800

do p = Mx*My+1, 2*Mx*My

  call mapgrinv(p,k,i,j)

  NN1 = U(k-1,i-2,j)
  N1  = U(k-1,i-1,j)
  C1  = U(k-1,i,j)
  S1  = U(k-1,i+1,j)
  SS1 = U(k-1,i+2,j)

  NN2 = U(k,i-2,j)
  N2  = U(k,i-1,j)
  C2  = U(k,i,j)
  S2  = U(k,i+1,j)
  SS2 = U(k,i+2,j)

  NN3 = U(k+1,i-2,j)
  N3  = U(k+1,i-1,j)
  C3  = U(k+1,i,j)
  S3  = U(k+1,i+1,j)
  SS3 = U(k+1,i+2,j)

  SqrtC1g = Sqrt(C1*g)

  if (feigvalues(1,i,j) .GE. 0d0) then
    if (verbose .EQ. 2) print *, "Compute JF1+ for ", i, j
  else
    if (verbose .EQ. 2) print *, "Compute JF1- for ", i, j
  end if
  JF(p,:)=0d0

! goto 786
    
  if (feigvalues(2,i,j) .GE. 0d0) then
  
    if (verbose .EQ. 2) print *, "Compute JF2+ for ", i, j
  
    call mapgr(ind,p,-1,-2,0) !NN1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*(C2**2/C1**2-C1*g)*(C2/C1-SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5*(C2**2/C1**2 - C1*g)*(C2/C1-SqrtC1g))/(dx*SqrtC1g))
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.25*(C2**2/C1**2 - C1*g))/(dx*SqrtC1g) + &
      & (0.0833333*((-2*C2**2)/C1**3 - g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g) - &
      & (0.0416667*g*(C2**2/C1**2 - C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*(C1*g)**1.5) + & 
      & (((0.25*C2*g)/(C1*(C1*g)**1.5) + (0.5*C2)/(C1**2*SqrtC1g))*(3*C2 - 6*N2 + NN2 + 2*S2))/(6.*dx)) + &
      & (-(C2/C1**2) - g/(2.*SqrtC1g))*((0.0833333*(C2**2/C1**2 - C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g) &
      & + ((0.5 - (0.5*C2)/(C1*SqrtC1g))*(3*C2 - 6*N2 + NN2 + 2*S2))/(6.*dx))
    end if  
    call mapgr(ind,p,-1,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*(C2**2/C1**2 - C1*g)*(C2/C1 - SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,-2,0) !NN2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(6.*dx)
    end if  
    call mapgr(ind,p,0,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/dx)
    end if  
    call mapgr(ind,p,0,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.5 - (0.5*C2)/(C1*SqrtC1g))/(2.*dx) + &
      & (0.166667*C2*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1**2*dx*SqrtC1g) - (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g)) + &
      & ((0.0833333*(C2**2/C1**2 - C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g) + & 
      & ((0.5 - (0.5*C2)/(C1*SqrtC1g))*(3*C2 - 6*N2 + NN2 + 2*S2))/(6.*dx))/C1
    end if  
    call mapgr(ind,p,0,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5-(0.5*C2)/(C1*SqrtC1g))*(C2/C1-SqrtC1g))/(3.*dx)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JF2- for ", i, j

    call mapgr(ind,p,-1,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.166667*(C2**2/C1**2 - C1*g)*(C2/C1 - SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((-0.25*(C2**2/C1**2 - C1*g))/(dx*SqrtC1g) +&
      & (0.0833333*((-2*C2**2)/C1**3 - g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g) -&
      & (0.0416667*g*(C2**2/C1**2 - C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*(C1*g)**1.5) +&
      & (((0.25*C2*g)/(C1*(C1*g)**1.5) + (0.5*C2)/(C1**2*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 - SS2))/(6.*dx)) +&
      & (-(C2/C1**2) - g/(2.*SqrtC1g))*((0.0833333*(C2**2/C1**2 - C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g) +&
      & ((0.5 - (0.5*C2)/(C1*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 - SS2))/(6.*dx))
    end if  
    call mapgr(ind,p,-1,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*(C2**2/C1**2 - C1*g)*(C2/C1 - SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,2,0)  !SS1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.0833333*(C2**2/C1**2 - C1*g)*(C2/C1 - SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,-1,0)  !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(3.*dx)
    end if  
    call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*(-(0.5 - (0.5*C2)/(C1*SqrtC1g))/(2.*dx) +&
      & (0.166667*C2*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1**2*dx*SqrtC1g) -&
      & (0.0833333*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*SqrtC1g)) +&
      & ((0.0833333*(C2**2/C1**2 - C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g) +&
      & ((0.5 - (0.5*C2)/(C1*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 - SS2))/(6.*dx))/C1
    end if  
    call mapgr(ind,p,0,1,0)   !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/dx
    end if  
    call mapgr(ind,p,0,2,0)   !SS2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5 - (0.5*C2)/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(6.*dx)
    end if  

  end if

  if (feigvalues(3,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JF3+ for ", i, j
  
    call mapgr(ind,p,-1,-2,0) !NN1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*(-(C2**2/C1**2) + C1*g)*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.5*(-(C2**2/C1**2) + C1*g)*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.25*(-(C2**2/C1**2) + C1*g))/(dx*SqrtC1g) +&
      & (0.0833333*((2*C2**2)/C1**3 + g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g) -&
      & (0.0416667*g*(-(C2**2/C1**2) + C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*(C1*g)**1.5) +&
      & (((-0.25*C2*g)/(C1*(C1*g)**1.5) - (0.5*C2)/(C1**2*SqrtC1g))*(3*C2 - 6*N2 + NN2 + 2*S2))/(6.*dx)) +&
      & (-(C2/C1**2) + g/(2.*SqrtC1g))*((0.0833333*(-(C2**2/C1**2) + C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g) +&
      & ((0.5 + (0.5*C2)/(C1*SqrtC1g))*(3*C2 - 6*N2 + NN2 + 2*S2))/(6.*dx))
    end if  
    call mapgr(ind,p,-1,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*(-(C2**2/C1**2) + C1*g)*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,-2,0)  !NN2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(6.*dx)
    end if  
    call mapgr(ind,p,0,-1,0)  !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/dx)
    end if  
    call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.5 + (0.5*C2)/(C1*SqrtC1g))/(2.*dx) -&
      & (0.166667*C2*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1**2*dx*SqrtC1g) +  (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g)) +&
      & ((0.0833333*(-(C2**2/C1**2) + C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g) + ((0.5 + (0.5*C2)/(C1*SqrtC1g))*(3*C2 - &
      & 6*N2 + NN2 + 2*S2))/(6.*dx))/C1
    end if  
    call mapgr(ind,p,0,1,0)   !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(3.*dx)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JF3- for ", i, j
  
    call mapgr(ind,p,-1,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.166667*(-(C2**2/C1**2) + C1*g)*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((-0.25*(-(C2**2/C1**2) + C1*g))/(dx*SqrtC1g) +&
      & (0.0833333*((2*C2**2)/C1**3 + g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g) -&
      & (0.0416667*g*(-(C2**2/C1**2) + C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*(C1*g)**1.5) + (((-0.25*C2*g)/(C1*(C1*g)**1.5) -&
      & (0.5*C2)/(C1**2*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 - SS2))/(6.*dx)) + (-(C2/C1**2) + &
      & g/(2.*SqrtC1g))*((0.0833333*(-(C2**2/C1**2) +&
      & C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g) + ((0.5 + (0.5*C2)/(C1*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 - SS2))/(6.*dx))
    end if  
    call mapgr(ind,p,-1,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*(-(C2**2/C1**2) + C1*g)*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,2,0)  !SS1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.0833333*(-(C2**2/C1**2) + C1*g)*(C2/C1 + SqrtC1g))/(dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(3.*dx)
    end if  
    call mapgr(ind,p,0,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*(-(0.5 + (0.5*C2)/(C1*SqrtC1g))/(2.*dx) -&
      & (0.166667*C2*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1**2*dx*SqrtC1g) + (0.0833333*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*SqrtC1g)) &
      & +&
      & ((0.0833333*(-(C2**2/C1**2) + C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g)+((0.5 + (0.5*C2)/(C1*SqrtC1g))*(-3*C2 - &
      & 2*N2 + 6*S2 - SS2))/(6.*dx))/C1
    end if  
    call mapgr(ind,p,0,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/dx
    end if  
    call mapgr(ind,p,0,2,0)  !SS2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((0.5 + (0.5*C2)/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(6.*dx)
    end if  

  end if

!786	continue
!goto 787

  WW1 = U(k-1,i,j-2)
  W1  = U(k-1,i,j-1)
  C1  = U(k-1,i,j)
  E1  = U(k-1,i,j+1)
  EE1 = U(k-1,i,j+2)

  WW2 = U(k,i,j-2)
  W2  = U(k,i,j-1)
  C2  = U(k,i,j)
  E2  = U(k,i,j+1)
  EE2 = U(k,i,j+2)

  WW3 = U(k+1,i,j-2)
  W3  = U(k+1,i,j-1)
  C3  = U(k+1,i,j)
  E3  = U(k+1,i,j+1)
  EE3 = U(k+1,i,j+2)

  SqrtC1g = Sqrt(C1*g)

  if (geigvalues(1,i,j) .GE. 0d0) then
  
    if (verbose .EQ. 2) print *, "Compute JG1+ for ", i, j
  
    call mapgr(ind,p,-1,0,-2) !WW1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-C2/(6.*C1*dy)
    end if
    call mapgr(ind,p,-1,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C2/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-C2/(2.*C1*dy) + (C2*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*C1**2*dy) - (C3*(3*C2 + 2*E2 - 6*W2 + WW2))/(6.*C1**2*dy)
    end if  
    call mapgr(ind,p,-1,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-C2/(3.*C1*dy)
    end if  
    call mapgr(ind,p,0,0,-2)  !WW2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C3/(6.*C1*dy)
    end if  
    call mapgr(ind,p,0,0,-1)  !W2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-(C3/(C1*dy))
    end if  
    call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C3/(2.*C1*dy) - (3*C1 + 2*E1 - 6*W1 + WW1)/(6.*C1*dy)
    end if
    call mapgr(ind,p,0,0,1)   !E2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C3/(3.*C1*dy)
    end if  
    call mapgr(ind,p,1,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=(3*C2 + 2*E2 - 6*W2 + WW2)/(6.*C1*dy)
    end if  

  else 

    if (verbose .EQ. 2) print *, "Compute JG1- for ", i, j

    call mapgr(ind,p,-1,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C2/(3.*C1*dy)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C2/(2.*C1*dy) + (C2*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*C1**2*dy) - (C3*(-3*C2 + 6*E2 - EE2 - 2*W2))/(6.*C1**2*dy)
    end if  
    call mapgr(ind,p,-1,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-(C2/(C1*dy))
    end if  
    call mapgr(ind,p,-1,0,2)  !EE1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C2/(6.*C1*dy)
    end if  
    call mapgr(ind,p,0,0,-1)  !W2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-C3/(3.*C1*dy)
    end if  
    call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-C3/(2.*C1*dy) - (-3*C1 + 6*E1 - EE1 - 2*W1)/(6.*C1*dy)
    end if  
    call mapgr(ind,p,0,0,1)   !E2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=C3/(C1*dy)
    end if  
    call mapgr(ind,p,0,0,2)   !EE2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=-C3/(6.*C1*dy)
    end if  
    call mapgr(ind,p,1,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=(-3*C2 + 6*E2 - EE2 - 2*W2)/(6.*C1*dy)
    end if  

  end if

  if (geigvalues(2,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JG2+ for ", i, j
  
    call mapgr(ind,p,-1,0,-2) !WW1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*C2*(1 + C3/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.25*C2*(1 + C3/(C1*SqrtC1g)))/(C1*dy) +&
      & (0.0833333*C2*(-(C3*g)/(2.*C1*(C1*g)**1.5) - C3/(C1**2*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) -&
      & (0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1**2*dy) + (0.0416667*C2*g*(3*C3 + 2*E3 - 6*W3 + &
      & WW3))/(C1*dy*(C1*g)**1.5) + &
      & (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1**2*dy*SqrtC1g)) + (-(C3/C1**2) - g/(2.*SqrtC1g))*((0.0833333*C2*(1 + &
      & C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) - (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,-1,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*C2*(1 + C3/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(C1*dy)
    end if  
	call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.0833333*(1 + C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) - &
      & (0.0833333*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,1,0,-2)  !WW3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*C2*(C3/C1 - SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*C2*(C3/C1 - SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((-0.25*C2)/(C1*dy*SqrtC1g) + (0.0833333*C2*(3*C1 + 2*E1 - 6*W1 + &
      & WW1))/(C1**2*dy*SqrtC1g)) +&
      & ((0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) - (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + &
      & WW3))/(C1*dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,0,1)   !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*C2*(C3/C1 - SqrtC1g))/(C1*dy*SqrtC1g)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JG2- for ", i, j
  
    call mapgr(ind,p,-1,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*C2*(1 + C3/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((-0.25*C2*(1 + C3/(C1*SqrtC1g)))/(C1*dy) + &
      & (0.0833333*C2*(-(C3*g)/(2.*C1*(C1*g)**1.5) -&
      & C3/(C1**2*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) - (0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - &
      & 2*W1))/(C1**2*dy) +&
      & (0.0416667*C2*g*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*(C1*g)**1.5) + (0.0833333*C2*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(C1**2*dy*SqrtC1g)) +&
      & (-(C3/C1**2) - g/(2.*SqrtC1g))*((0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) -&
      & (0.0833333*C2*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,-1,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*C2*(1 + C3/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,2)  !EE1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.0833333*(1 + C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) -&
      & (0.0833333*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,1,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*C2*(C3/C1 - SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.25*C2)/(C1*dy*SqrtC1g) + (0.0833333*C2*(-3*C1 + 6*E1 - EE1 - &
      & 2*W1))/(C1**2*dy*SqrtC1g)) +&
      & ((0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) - (0.0833333*C2*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(C1*dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,0,1)   !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*C2*(C3/C1 - SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,2)   !EE3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*C2*(C3/C1 - SqrtC1g))/(C1*dy*SqrtC1g)
    end if  

  end if    

  if (geigvalues(3,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JG3+ for ", i, j

    call mapgr(ind,p,-1,0,-2) !WW1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*C2*(1 - C3/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.25*C2*(1 - C3/(C1*SqrtC1g)))/(C1*dy) + (0.0833333*C2*((C3*g)/(2.*C1*(C1*g)**1.5) &
      & +&
      & C3/(C1**2*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) - (0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + &
      & WW1))/(C1**2*dy) -&
      & (0.0416667*C2*g*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*(C1*g)**1.5) - (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + &
      & WW3))/(C1**2*dy*SqrtC1g)) + (-(C3/C1**2) +&
      & g/(2.*SqrtC1g))*((0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) +&
      & (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,-1,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*C2*(1 - C3/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.0833333*(1 - C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) +&
      & (0.0833333*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,1,0,-2)  !WW3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*C2*(C3/C1 + SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*C2*(C3/C1 + SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.25*C2)/(C1*dy*SqrtC1g) - (0.0833333*C2*(3*C1 + 2*E1 - 6*W1 + &
      & WW1))/(C1**2*dy*SqrtC1g)) +&
      & ((0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) + (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + &
      & WW3))/(C1*dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,0,1)   !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*C2*(C3/C1 + SqrtC1g))/(C1*dy*SqrtC1g)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JG3- for ", i, j

    call mapgr(ind,p,-1,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*C2*(1 - C3/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((-0.25*C2*(1 - C3/(C1*SqrtC1g)))/(C1*dy) + (0.0833333*C2*((C3*g)/(2.*C1*(C1*g)**1.5) &
      & +&
      & C3/(C1**2*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) - (0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - &
      & 2*W1))/(C1**2*dy) -&
      & (0.0416667*C2*g*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*(C1*g)**1.5) - (0.0833333*C2*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(C1**2*dy*SqrtC1g)) + (-(C3/C1**2) +&
      & g/(2.*SqrtC1g))*((0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) + (0.0833333*C2*(-3*C3 + 6*E3 - &
      & EE3 - 2*W3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,-1,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*C2*(1 - C3/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,-1,0,2)  !EE1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(C1*dy)
    end if  
    call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.0833333*(1 - C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) + &
      & (0.0833333*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,1,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*C2*(C3/C1 + SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((-0.25*C2)/(C1*dy*SqrtC1g) - (0.0833333*C2*(-3*C1 + 6*E1 - EE1 - &
      & 2*W1))/(C1**2*dy*SqrtC1g)) +&
      & ((0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) + (0.0833333*C2*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(C1*dy*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,1,0,1)   !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*C2*(C3/C1 + SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,2)   !EE3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*C2*(C3/C1 + SqrtC1g))/(C1*dy*SqrtC1g)
    end if

  end if

787  continue

end do

!stop

do p=2*Mx*My+1, ndim

  call mapgrinv(p,k,i,j)

!print *, p, k, i, j

  NN1 = U(k-2,i-2,j)
  N1  = U(k-2,i-1,j)
  C1  = U(k-2,i,j)
  S1  = U(k-2,i+1,j)
  SS1 = U(k-2,i+2,j)

  NN2 = U(k-1,i-2,j)
  N2  = U(k-1,i-1,j)
  C2  = U(k-1,i,j)
  S2  = U(k-1,i+1,j)
  SS2 = U(k-1,i+2,j)
  
  NN3 = U(k,i-2,j)
  N3  = U(k,i-1,j)
  C3  = U(k,i,j)
  S3  = U(k,i+1,j)
  SS3 = U(k,i+2,j)

  SqrtC1g = Sqrt(C1*g)

!goto 796

  if (feigvalues(1,i,j) .GE. 0d0) then
  
    if (verbose .EQ. 2) print *, "Compute JF1+ for ", i, j
  
    call mapgr(ind,p,-2,-2,0) !NN1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((C2*C3)/(6.*C1**2*dx))
    end if
    call mapgr(ind,p,-2,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2*C3)/(C1**2*dx)
    end if
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2*(-C3/(2.*C1*dx) + (C3*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*C1**2*dx)))/C1 - (C2*(-(C3*(3*C1 - 6*N1 + NN1 &
      & + 2*S1))/(6.*C1*dx) + (3*C3 - 6*N3 + NN3 + 2*S3)/(6.*dx)))/C1**2
    end if  
    call mapgr(ind,p,-2,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((C2*C3)/(3*C1**2*dx))
    end if  
    call mapgr(ind,p,-1,0,0)  !C2 - ???
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-(C3*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*C1*dx) + (3*C3 - 6*N3 + NN3 + 2*S3)/(6.*dx))/C1
    end if  
    call mapgr(ind,p,0,-2,0)  !NN3
    if ((ind .GE. 0) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+C2/(6.*C1*dx)
    end if  
    call mapgr(ind,p,0,-1,0)  !N3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(C2/(C1*dx))
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2*(1/(2*dx)-(3*C1-6*N1+NN1+2*S1)/(6*C1*dx)))/C1
    end if  
    call mapgr(ind,p,0,1,0)   !S3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+C2/(3*C1*dx)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JF1- for ", i, j

    call mapgr(ind,p,-2,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2*C3)/(3.*C1**2*dx)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2*(C3/(2.*C1*dx) + (C3*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*C1**2*dx)))/C1 - (C2*(-(C3*(-3*C1 - 2*N1 + &
      & 6*S1 - SS1))/(6.*C1*dx) + (-3*C3 - 2*N3 + 6*S3 - SS3)/(6.*dx)))/C1**2
    end if  
    call mapgr(ind,p,-2,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-((C2*C3)/(C1**2*dx))
    end if  
    call mapgr(ind,p,-2,2,0)  !SS1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2*C3)/(6.*C1**2*dx)
    end if
    call mapgr(ind,p,-1,0,0)  !C2 ???
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-(C3*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*C1*dx) + (-3*C3 - 2*N3 + 6*S3 - SS3)/(6.*dx))/C1
    end if  
    call mapgr(ind,p,0,-1,0)  !N3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(C2/(3.*C1*dx))
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2*(-1/(2.*dx) - (-3*C1 - 2*N1 + 6*S1 - SS1)/(6.*C1*dx)))/C1
    end if  
    call mapgr(ind,p,0,1,0)   !S3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+C2/(C1*dx)
    end if  
    call mapgr(ind,p,0,2,0)   !SS3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(C2/(6.*C1*dx))
    end if

  end if

  if (feigvalues(2,i,j) .GE. 0d0) then
  
    if (verbose .EQ. 2) print *, "Compute JF2+ for ", i, j  
  
    call mapgr(ind,p,-2,-2,0) !NN1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)-(-0.5*C3*(1 + C2/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.25*C3*(1 + C2/(C1*SqrtC1g)))/(C1*dx) +&
      & (0.0833333*C3*(-(C2*g)/(2.*C1*(C1*g)**1.5) - C2/(C1**2*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) -&
      & (0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1**2*dx) +&
      & (0.0416667*C3*g*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*(C1*g)**1.5) +&
      & (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1**2*dx*SqrtC1g)) + (-(C2/C1**2) -&
      & g/(2.*SqrtC1g))*((0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) -&
      & (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))
    end if  
    call mapgr(ind,p,-2,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*C3*(1 + C2/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-1,-2,0) !NN2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.0833333*C3*(C2/C1 - SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*C3*(C2/C1 - SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((-0.25*C3)/(C1*dx*SqrtC1g) +&
      & (0.0833333*C3*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1**2*dx*SqrtC1g)) +&
      & ((0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) -&
      & (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,-1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.166667*C3*(C2/C1 - SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.0833333*(1 + C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) - &
      & (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))
    end if  

  else
  
    if (verbose .EQ. 2) print *, "Compute JF2- for ", i, j    
  
    call mapgr(ind,p,-2,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.166667*C3*(1 + C2/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((-0.25*C3*(1 + C2/(C1*SqrtC1g)))/(C1*dx) +&
      & (0.0833333*C3*(-(C2*g)/(2.*C1*(C1*g)**1.5) - C2/(C1**2*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) -&
      & (0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1**2*dx) +&
      & (0.0416667*C3*g*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*(C1*g)**1.5) +&
      & (0.0833333*C3*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1**2*dx*SqrtC1g)) +&
      & (-(C2/C1**2) - g/(2.*SqrtC1g))*((0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) -&
      & (0.0833333*C3*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*SqrtC1g))
    end if  
    call mapgr(ind,p,-2,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*C3*(1 + C2/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,2,0)  !SS1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(C2/C1 - SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-1,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*C3*(C2/C1 - SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.25*C3)/(C1*dx*SqrtC1g) +&
      & (0.0833333*C3*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1**2*dx*SqrtC1g)) +&
      & ((0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) - (0.0833333*C3*(-3*C2 - 2*N2 + 6*S2 - &
      & SS2))/(C1*dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,-1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.5*C3*(C2/C1 - SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,2,0)  !SS2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*C3*(C2/C1 - SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.0833333*(1 + C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) -&
      & (0.0833333*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*SqrtC1g))
    end if  

  end if

  if (feigvalues(3,i,j) .GE. 0d0) then
  
    if (verbose .EQ. 2) print *, "Compute JF3+ for ", i, j      
  
    call mapgr(ind,p,-2,-2,0) !NN1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.5*C3*(1 - C2/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.25*C3*(1 - C2/(C1*SqrtC1g)))/(C1*dx) +&
      & (0.0833333*C3*((C2*g)/(2.*C1*(C1*g)**1.5) + C2/(C1**2*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) -&
      & (0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1**2*dx) - (0.0416667*C3*g*(3*C2 - 6*N2 + NN2 + &
      & 2*S2))/(C1*dx*(C1*g)**1.5) -&
      & (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1**2*dx*SqrtC1g)) + (-(C2/C1**2) + g/(2.*SqrtC1g))*((0.0833333*C3*(1 - &
      & C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) + (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))
    end if  
    call mapgr(ind,p,-2,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*C3*(1 - C2/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-1,-2,0) !NN2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.0833333*C3*(C2/C1 + SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,-1,0) !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.5*C3*(C2/C1 + SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,0,0)  !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.25*C3)/(C1*dx*SqrtC1g) - (0.0833333*C3*(3*C1 - 6*N1 + NN1 + &
      & 2*S1))/(C1**2*dx*SqrtC1g)) +&
      & ((0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) + (0.0833333*C3*(3*C2 - 6*N2 + NN2 + &
      & 2*S2))/(C1*dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,-1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*C3*(C2/C1 + SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.0833333*(1 - C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) + &
      & (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JF3- for ", i, j      
  
    call mapgr(ind,p,-2,-1,0) !N1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.166667*C3*(1 - C2/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((-0.25*C3*(1 - C2/(C1*SqrtC1g)))/(C1*dx) +&
      & (0.0833333*C3*((C2*g)/(2.*C1*(C1*g)**1.5) + C2/(C1**2*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) -&
      & (0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1**2*dx) - (0.0416667*C3*g*(-3*C2 - 2*N2 + 6*S2 - &
      & SS2))/(C1*dx*(C1*g)**1.5) -&
      & (0.0833333*C3*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1**2*dx*SqrtC1g)) + (-(C2/C1**2) + g/(2.*SqrtC1g))*((0.0833333*C3*(1 - &
      & C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) + (0.0833333*C3*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*SqrtC1g))
    end if  
    call mapgr(ind,p,-2,1,0)  !S1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*C3*(1 - C2/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-2,2,0)  !SS1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(C2/C1 + SqrtC1g))/(C1*dx)
    end if  
    call mapgr(ind,p,-1,-1,0)  !N2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.166667*C3*(C2/C1 + SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,-1,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((-0.25*C3)/(C1*dx*SqrtC1g) -&
      & (0.0833333*C3*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1**2*dx*SqrtC1g)) +&
      & ((0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) + (0.0833333*C3*(-3*C2 - 2*N2 + 6*S2 - &
      & SS2))/(C1*dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,-1,1,0)   !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.5*C3*(C2/C1 + SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,-1,2,0)   !SS2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(-0.0833333*C3*(C2/C1 + SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.0833333*(1 - C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) + &
      & (0.0833333*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*SqrtC1g))
    end if  

  end if

!796	continue
!goto 797
  
  WW1 = U(k-2,i,j-2)
  W1  = U(k-2,i,j-1)
  C1  = U(k-2,i,j)
  E1  = U(k-2,i,j+1)
  EE1 = U(k-2,i,j+2)

  WW2 = U(k-1,i,j-2)
  W2  = U(k-1,i,j-1)
  C2  = U(k-1,i,j)
  E2  = U(k-1,i,j+1)
  EE2 = U(k-1,i,j+2)

  WW3 = U(k,i,j-2)
  W3  = U(k,i,j-1)
  C3  = U(k,i,j)    
  E3  = U(k,i,j+1)
  EE3 = U(k,i,j+2)

  SqrtC1g = Sqrt(C1*g)

  JG(p,:)=0d0
  
  if (geigvalues(2,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JG2+ for ", i, j      
  
    call mapgr(ind,p,-2,0,-2) !WW1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*(C3**2/C1**2 - C1*g)*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*(C3**2/C1**2 - C1*g)*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.25*(C3**2/C1**2 - C1*g))/(dy*SqrtC1g) + &
      &(0.0833333*((-2*C3**2)/C1**3 - g)*(3*C1 + 2*E1 - 6*W1 + WW1))/(dy*SqrtC1g) - &
      &(0.0416667*g*(C3**2/C1**2 - C1*g)*(3*C1 + 2*E1 - 6*W1 + WW1))/(dy*(C1*g)**1.5) + (((0.25*C3*g)/(C1*(C1*g)**1.5) +&
      & (0.5*C3)/(C1**2*SqrtC1g))*(3*C3 + 2*E3 - 6*W3 + WW3))/(6.*dy)) + (-(C3/C1**2) -&
      & g/(2.*SqrtC1g))*((0.0833333*(C3**2/C1**2 - C1*g)*(3*C1 + 2*E1 - 6*W1 + WW1))/(dy*SqrtC1g) + ((0.5 - &
      & (0.5*C3)/(C1*SqrtC1g))*(3*C3 + 2*E3 - 6*W3 + WW3))/(6.*dy))
    end if  
    call mapgr(ind,p,-2,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*(C3**2/C1**2 - C1*g)*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,-2)  !WW3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(6.*dy)
    end if  
    call mapgr(ind,p,0,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-(((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/dy)
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.5 - (0.5*C3)/(C1*SqrtC1g))/(2.*dy) +&
      & (0.166667*C3*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1**2*dy*SqrtC1g) -  (0.0833333*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*SqrtC1g)) +&
      & ((0.0833333*(C3**2/C1**2 - C1*g)*(3*C1 + 2*E1 - 6*W1 + WW1))/(dy*SqrtC1g) + ((0.5 - (0.5*C3)/(C1*SqrtC1g))*(3*C3 + 2*E3 - &
      & 6*W3 + WW3))/(6.*dy))/C1
    end if  
    call mapgr(ind,p,0,0,1)  !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(3.*dy)
    end if  

  else
  
    if (verbose .EQ. 2) print *, "Compute JG2- for ", i, j      
  
    call mapgr(ind,p,-2,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*(C3**2/C1**2 - C1*g)*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((-0.25*(C3**2/C1**2 - C1*g))/(dy*SqrtC1g) +&
      & (0.0833333*((-2*C3**2)/C1**3 - g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*SqrtC1g) -&
      & (0.0416667*g*(C3**2/C1**2 - C1*g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*(C1*g)**1.5) + &
      &(((0.25*C3*g)/(C1*(C1*g)**1.5) + (0.5*C3)/(C1**2*SqrtC1g))*(-3*C3 + 6*E3 - EE3 - 2*W3))/(6.*dy)) +&
      & (-(C3/C1**2) - g/(2.*SqrtC1g))*((0.0833333*(C3**2/C1**2 - C1*g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*SqrtC1g) +&
      & ((0.5 - (0.5*C3)/(C1*SqrtC1g))*(-3*C3 + 6*E3 - EE3 - 2*W3))/(6.*dy))
    end if  
    call mapgr(ind,p,-2,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*(C3**2/C1**2 - C1*g)*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,2)  !EE1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*(C3**2/C1**2 - C1*g)*(C3/C1 - SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(3.*dy)
    end if  
    call mapgr(ind,p,0,0,0)  !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*(-(0.5 - (0.5*C3)/(C1*SqrtC1g))/(2.*dy) +&
      & (0.166667*C3*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1**2*dy*SqrtC1g) -	(0.0833333*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*SqrtC1g)) &
      & +&
      & ((0.0833333*(C3**2/C1**2 - C1*g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*SqrtC1g) + ((0.5 - (0.5*C3)/(C1*SqrtC1g))*(-3*C3 + 6*E3 &
      & - EE3 - 2*W3))/(6.*dy))/C1
    end if  
    call mapgr(ind,p,0,0,1)  !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/dy
    end if  
    call mapgr(ind,p,0,0,2)  !EE3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 - (0.5*C3)/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(6.*dy)
    end if  

  end if
  
  if (geigvalues(3,i,j) .GE. 0d0) then

    if (verbose .EQ. 2) print *, "Compute JG3+ for ", i, j      
  
    call mapgr(ind,p,-2,0,-2) !WW1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.0833333*(-(C3**2/C1**2) + C1*g)*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.5*(-(C3**2/C1**2) + C1*g)*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.25*(-(C3**2/C1**2) + C1*g))/(dy*SqrtC1g) +&
      & (0.0833333*((2*C3**2)/C1**3 + g)*(3*C1 + 2*E1 - 6*W1 + WW1))/(dy*SqrtC1g) -&
      & (0.0416667*g*(-(C3**2/C1**2) + C1*g)*(3*C1 + 2*E1 - 6*W1 + WW1))/(dy*(C1*g)**1.5) + (((-0.25*C3*g)/(C1*(C1*g)**1.5) -&
      & (0.5*C3)/(C1**2*SqrtC1g))*(3*C3 + 2*E3 - 6*W3 + WW3))/(6.*dy)) + (-(C3/C1**2) + &
      & g/(2.*SqrtC1g))*((0.0833333*(-(C3**2/C1**2) + C1*g)*(3*C1 + 2*E1 - 6*W1 + WW1))/(dy*SqrtC1g) + ((0.5 + &
      & (0.5*C3)/(C1*SqrtC1g))*(3*C3 + 2*E3 - 6*W3 + WW3))/(6.*dy))
    end if  
    call mapgr(ind,p,-2,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*(-(C3**2/C1**2) + C1*g)*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,-2)  !WW3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(6.*dy)
    end if  
    call mapgr(ind,p,0,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-(((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/dy)
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.5 + (0.5*C3)/(C1*SqrtC1g))/(2.*dy) - (0.166667*C3*(3*C1 + 2*E1 - 6*W1 + &
      & WW1))/(C1**2*dy*SqrtC1g) + &
      & (0.0833333*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*SqrtC1g)) + ((0.0833333*(-(C3**2/C1**2) + C1*g)*(3*C1 + 2*E1 - 6*W1 + &
      & WW1))/(dy*SqrtC1g) + ((0.5 + (0.5*C3)/(C1*SqrtC1g))*(3*C3 + 2*E3 - 6*W3 + WW3))/(6.*dy))/C1
    end if  
    call mapgr(ind,p,0,0,1)   !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(3.*dy)
    end if  

  else

    if (verbose .EQ. 2) print *, "Compute JG3- for ", i, j      
  
    call mapgr(ind,p,-2,0,-1) !W1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*(-(C3**2/C1**2) + C1*g)*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,0)  !C1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((-0.25*(-(C3**2/C1**2) + C1*g))/(dy*SqrtC1g) +&
      & (0.0833333*((2*C3**2)/C1**3 + g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*SqrtC1g) - &
      & (0.0416667*g*(-(C3**2/C1**2) + C1*g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*(C1*g)**1.5) + (((-0.25*C3*g)/(C1*(C1*g)**1.5) -&
      & (0.5*C3)/(C1**2*SqrtC1g))*(-3*C3 + 6*E3 - EE3 - 2*W3))/(6.*dy)) + (-(C3/C1**2) +&
      & g/(2.*SqrtC1g))*((0.0833333*(-(C3**2/C1**2) + C1*g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*SqrtC1g) + ((0.5 + &
      & (0.5*C3)/(C1*SqrtC1g))*(-3*C3 + 6*E3 - EE3 - 2*W3))/(6.*dy)) 
    end if  
    call mapgr(ind,p,-2,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.5*(-(C3**2/C1**2) + C1*g)*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,-2,0,2)  !EE1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.0833333*(-(C3**2/C1**2) + C1*g)*(C3/C1 + SqrtC1g))/(dy*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,-1) !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(3.*dy)
    end if  
    call mapgr(ind,p,0,0,0)  !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*(-(0.5 + (0.5*C3)/(C1*SqrtC1g))/(2.*dy) -&
      & (0.166667*C3*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1**2*dy*SqrtC1g) + (0.0833333*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*SqrtC1g)) &
      & +&
      & ((0.0833333*(-(C3**2/C1**2) + C1*g)*(-3*C1 + 6*E1 - EE1 - 2*W1))/(dy*SqrtC1g) + ((0.5 + (0.5*C3)/(C1*SqrtC1g))*(-3*C3 + &
      & 6*E3 - EE3 - 2*W3))/(6.*dy))/C1
    end if  
    call mapgr(ind,p,0,0,1)  !E3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/dy
    end if  
    call mapgr(ind,p,0,0,2)  !EE3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)-((0.5 + (0.5*C3)/(C1*SqrtC1g))*(C3/C1 + SqrtC1g))/(6.*dy)
    end if  

  end if          

!797	continue

end do

!800	continue

  JF(:,:)=-JF(:,:)-JG(:,:)
    

end subroutine

subroutine g_adcompute_f( u, adf, g_u, g_adu )
!g_adu = d(Jac*adf)/dy * g_u = (H x g_u) *adf
!***************************************************************
!***************************************************************
!** This routine was generated by the                         **
!** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
!***************************************************************
!***************************************************************
!==============================================
! referencing used modules
!==============================================
use swe2dxy_parameters

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

!==============================================
! define arguments
!==============================================
double precision adf(3,mx+4,my+4)
double precision g_adu(3,mx+4,my+4)
double precision g_u(3,mx+4,my+4)
double precision u(3,mx+4,my+4)

!==============================================
! define local variables
!==============================================
double precision adeigf(3)
double precision adeigg(3)
double precision adfu1(3,3)
double precision adfu2(3,3)
double precision adfu3(3,3)
double precision adgu1(3,3)
double precision adgu2(3,3)
double precision adgu3(3,3)
double precision adisqrt
double precision adrhs(3)
double precision adtmpf(3)
double precision adtmpg(3)
double precision adtmpmv(3)
double precision aduc
double precision advc
double precision dudxminus(3)
double precision dudxplus(3)
double precision dudyminus(3)
double precision dudyplus(3)
double precision eigf(3)
double precision eigg(3)
double precision fu1(3,3)
double precision fu2(3,3)
double precision fu3(3,3)
double precision g
double precision g_addudxminus(3)
double precision g_addudxplus(3)
double precision g_addudyminus(3)
double precision g_addudyplus(3)
double precision g_adeigf(3)
double precision g_adeigg(3)
double precision g_adfu1(3,3)
double precision g_adfu2(3,3)
double precision g_adfu3(3,3)
double precision g_adgu1(3,3)
double precision g_adgu2(3,3)
double precision g_adgu3(3,3)
double precision g_adhc
double precision g_adisqrt
double precision g_adtmpmv(3)
double precision g_aduc
double precision g_advc
double precision g_dudxminus(3)
double precision g_dudxplus(3)
double precision g_dudyminus(3)
double precision g_dudyplus(3)
double precision g_eigf(3)
double precision g_eigg(3)
double precision g_fu1(3,3)
double precision g_fu2(3,3)
double precision g_fu3(3,3)
double precision g_gu1(3,3)
double precision g_gu2(3,3)
double precision g_gu3(3,3)
double precision g_hc
double precision g_isqrt
double precision g_tmpmv(3)
double precision g_uc
double precision g_vc
double precision gu1(3,3)
double precision gu2(3,3)
double precision gu3(3,3)
double precision hc
integer i
integer i1
double precision isqrt
integer j
integer j1
double precision tmpmv(3)
double precision uc
double precision vc

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------
g_addudxminus(:) = 0.d0
g_addudxplus(:) = 0.d0
g_addudyminus(:) = 0.d0
g_addudyplus(:) = 0.d0
g_adeigf(:) = 0.d0
adeigf(:) = 0.d0
g_adeigg(:) = 0.d0
adeigg(:) = 0.d0
g_adfu1(:,:) = 0.d0
adfu1(:,:) = 0.d0
g_adfu2(:,:) = 0.d0
adfu2(:,:) = 0.d0
g_adfu3(:,:) = 0.d0
adfu3(:,:) = 0.d0
g_adgu1(:,:) = 0.d0
adgu1(:,:) = 0.d0
g_adgu2(:,:) = 0.d0
adgu2(:,:) = 0.d0
g_adgu3(:,:) = 0.d0
adgu3(:,:) = 0.d0
g_adhc = 0.d0
g_adisqrt = 0.d0
adisqrt = 0.d0
adrhs(:) = 0.d0
adtmpf(:) = 0.d0
adtmpg(:) = 0.d0
g_adtmpmv(:) = 0.d0
adtmpmv(:) = 0.d0
g_aduc = 0.d0
aduc = 0.d0
g_advc = 0.d0
advc = 0.d0
g_u(:,2,3:my+2) = g_u(:,mx+2,3:my+2)
u(:,2,3:my+2) = u(:,mx+2,3:my+2)
g_u(:,1,3:my+2) = g_u(:,mx+1,3:my+2)
u(:,1,3:my+2) = u(:,mx+1,3:my+2)
g_u(:,mx+3,3:my+2) = g_u(:,3,3:my+2)
u(:,mx+3,3:my+2) = u(:,3,3:my+2)
g_u(:,mx+4,3:my+2) = g_u(:,4,3:my+2)
u(:,mx+4,3:my+2) = u(:,4,3:my+2)
g_u(:,3:mx+2,2) = g_u(:,3:mx+2,my+2)
u(:,3:mx+2,2) = u(:,3:mx+2,my+2)
g_u(:,3:mx+2,1) = g_u(:,3:mx+2,my+1)
u(:,3:mx+2,1) = u(:,3:mx+2,my+1)
g_u(:,3:mx+2,my+3) = g_u(:,3:mx+2,3)
u(:,3:mx+2,my+3) = u(:,3:mx+2,3)
g_u(:,3:mx+2,my+4) = g_u(:,3:mx+2,4)
u(:,3:mx+2,my+4) = u(:,3:mx+2,4)
do j = my+2, 3, -1
  do i = mx+2, 3, -1
    g_eigf(1) = g_u(2,i,j)/u(1,i,j)-g_u(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))
    eigf(1) = u(2,i,j)/u(1,i,j)
    g_eigf(2) = g_u(2,i,j)/u(1,i,j)-g_u(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j))&
                +1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigf(2) = u(2,i,j)/u(1,i,j)-sqrt(gacc*u(1,i,j))
    g_eigf(3) = g_u(2,i,j)/u(1,i,j)+g_u(1,i,j)*((-(u(2,i,j)/(u(1,i,j)*u(1,i,j))))&
                +1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigf(3) = u(2,i,j)/u(1,i,j)+sqrt(gacc*u(1,i,j))
    g_hc = g_u(1,i,j)
    hc = u(1,i,j)
    g_uc = g_u(2,i,j)/u(1,i,j)-g_u(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))
    uc = u(2,i,j)/u(1,i,j)
    g_vc = g_u(3,i,j)/u(1,i,j)-g_u(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))
    vc = u(3,i,j)/u(1,i,j)
    g = gacc
    g_fu1(:,:) = 0.d0
    fu1(:,:) = 0.d0
    g_fu1(3,1) = -g_vc
    fu1(3,1) = -vc
    g_fu1(3,3) = 0.d0
    fu1(3,3) = 1.d0
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_fu2(:,:) = 0.d0
    fu2(:,:) = 0.d0
    g_fu2(1,1) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    fu2(1,1) = 0.5d0*(1.d0+uc*isqrt)
    g_fu2(1,2) = (-0.5d0)*g_isqrt
    fu2(1,2) = -(0.5d0*isqrt)
    g_fu2(2,1) = (-(0.5d0*g_hc*g*isqrt))+0.5d0*g_isqrt*(uc*uc-g*hc)+g_uc*uc*isqrt
    fu2(2,1) = 0.5d0*(uc*uc-g*hc)*isqrt
    g_fu2(2,2) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    fu2(2,2) = 0.5d0-0.5d0*uc*isqrt
    g_fu2(3,1) = 0.5d0*g_isqrt*uc*vc+0.5d0*g_uc*isqrt*vc+0.5d0*g_vc*(1.d0+uc*isqrt)
    fu2(3,1) = 0.5d0*(1.d0+uc*isqrt)*vc
    g_fu2(3,2) = (-(0.5d0*g_isqrt*vc))-0.5d0*g_vc*isqrt
    fu2(3,2) = -(0.5d0*vc*isqrt)
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_fu3(:,:) = 0.d0
    fu3(:,:) = 0.d0
    g_fu3(1,1) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    fu3(1,1) = 0.5d0-0.5d0*uc*isqrt
    g_fu3(1,2) = 0.5d0*g_isqrt
    fu3(1,2) = 0.5d0*isqrt
    g_fu3(2,1) = 0.5d0*g_hc*g*isqrt+0.5d0*g_isqrt*(g*hc-uc*uc)-g_uc*uc*isqrt
    fu3(2,1) = 0.5d0*(g*hc-uc*uc)*isqrt
    g_fu3(2,2) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    fu3(2,2) = 0.5d0*(1.d0+uc*isqrt)
    g_fu3(3,1) = (-(0.5d0*g_isqrt*uc*vc))-0.5d0*g_uc*isqrt*vc+0.5d0*g_vc*(1.d0-uc*isqrt)
    fu3(3,1) = 0.5d0*(1.d0-uc*isqrt)*vc
    g_fu3(3,2) = 0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt
    fu3(3,2) = 0.5d0*vc*isqrt
    g_dudxminus = g_u(:,i-1,j)*((-2)/(6.d0*dx))-g_u(:,i+2,j)/(6.d0*dx)+g_u(:,i+1,j)*(6/(6.d0*dx))+g_u(:,i,j)*((-3)/(6.d0*dx))
    dudxminus = ((-u(:,i+2,j))+6.d0*u(:,i+1,j)-3.d0*u(:,i,j)-2.d0*u(:,i-1,j))/(6.d0*dx)
    g_dudxplus = g_u(:,i-2,j)/(6.d0*dx)+g_u(:,i-1,j)*((-6)/(6.d0*dx))+g_u(:,i+1,j)*(2/(6.d0*dx))+g_u(:,i,j)*(3/(6.d0*dx))
    dudxplus = (2.d0*u(:,i+1,j)+3.d0*u(:,i,j)-6.d0*u(:,i-1,j)+u(:,i-2,j))/(6.d0*dx)
    g_eigg(1) = g_u(3,i,j)/u(1,i,j)-g_u(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))
    eigg(1) = u(3,i,j)/u(1,i,j)
    g_eigg(2) = g_u(3,i,j)/u(1,i,j)-g_u(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigg(2) = u(3,i,j)/u(1,i,j)-sqrt(gacc*u(1,i,j))
    g_eigg(3) = g_u(3,i,j)/u(1,i,j)+g_u(1,i,j)*((-(u(3,i,j)/(u(1,i,j)*u(1,i,j))))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigg(3) = u(3,i,j)/u(1,i,j)+sqrt(gacc*u(1,i,j))
    g_gu1(:,:) = 0.d0
    gu1(:,:) = 0.d0
    g_gu1(2,1) = -g_uc
    gu1(2,1) = -uc
    g_gu1(2,2) = 0.d0
    gu1(2,2) = 1.d0
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_gu2(:,:) = 0.d0
    gu2(:,:) = 0.d0
    g_gu2(1,1) = 0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt
    gu2(1,1) = 0.5d0*(1.d0+vc*isqrt)
    g_gu2(1,3) = (-0.5d0)*g_isqrt
    gu2(1,3) = -(0.5d0*isqrt)
    g_gu2(2,1) = 0.5d0*g_isqrt*uc*vc+0.5d0*g_uc*(1.d0+vc*isqrt)+0.5d0*g_vc*uc*isqrt
    gu2(2,1) = 0.5d0*uc*(1.d0+vc*isqrt)
    g_gu2(2,3) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    gu2(2,3) = -(0.5d0*uc*isqrt)
    g_gu2(3,1) = (-(0.5d0*g_hc*g*isqrt))+0.5d0*g_isqrt*(vc*vc-g*hc)+g_vc*vc*isqrt
    gu2(3,1) = 0.5d0*(vc*vc-g*hc)*isqrt
    g_gu2(3,3) = (-(0.5d0*g_isqrt*vc))-0.5d0*g_vc*isqrt
    gu2(3,3) = 0.5d0-0.5d0*vc*isqrt
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_gu3(:,:) = 0.d0
    gu3(:,:) = 0.d0
    g_gu3(1,1) = (-(0.5d0*g_isqrt*vc))-0.5d0*g_vc*isqrt
    gu3(1,1) = 0.5d0*(1.d0-vc*isqrt)
    g_gu3(1,3) = 0.5d0*g_isqrt
    gu3(1,3) = 0.5d0*isqrt
    g_gu3(2,1) = (-(0.5d0*g_isqrt*uc*vc))+0.5d0*g_uc*(1.d0-vc*isqrt)-0.5d0*g_vc*uc*isqrt
    gu3(2,1) = 0.5d0*uc*(1.d0-vc*isqrt)
    g_gu3(2,3) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    gu3(2,3) = 0.5d0*uc*isqrt
    g_gu3(3,1) = 0.5d0*g_hc*g*isqrt+0.5d0*g_isqrt*((-(vc*vc))+g*hc)-g_vc*vc*isqrt
    gu3(3,1) = 0.5d0*((-(vc*vc))+g*hc)*isqrt
    g_gu3(3,3) = 0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt
    gu3(3,3) = 0.5d0+0.5d0*vc*isqrt
    g_dudyminus = g_u(:,i,j-1)*((-2)/(6.d0*dy))-g_u(:,i,j+2)/(6.d0*dy)+g_u(:,i,j+1)*(6/(6.d0*dy))+g_u(:,i,j)*((-3)/(6.d0*dy))
    dudyminus = ((-u(:,i,j+2))+6.d0*u(:,i,j+1)-3.d0*u(:,i,j)-2.d0*u(:,i,j-1))/(6.d0*dy)
    g_dudyplus = g_u(:,i,j-2)/(6.d0*dy)+g_u(:,i,j-1)*((-6)/(6.d0*dy))+g_u(:,i,j+1)*(2/(6.d0*dy))+g_u(:,i,j)*(3/(6.d0*dy))
    dudyplus = (2.d0*u(:,i,j+1)+3.d0*u(:,i,j)-6.d0*u(:,i,j-1)+u(:,i,j-2))/(6.d0*dy)
    adrhs = adrhs+adf(:,i,j)
    adf(:,i,j) = 0.d0
    adtmpf = adtmpf-adrhs
    adtmpg = adtmpg-adrhs
    adrhs = 0.d0
    if (eigg(3) .ge. 0.d0) then
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudyplus(j1)*gu3(i1,j1)+g_gu3(i1,j1)*dudyplus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+gu3(i1,j1)*dudyplus(j1)
        end do
      end do
      g_adeigg(3) = g_adeigg(3)+sum(adtmpg*g_tmpmv)
      adeigg(3) = adeigg(3)+sum(adtmpg*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigg(3)*adtmpg
      adtmpmv = adtmpmv+adtmpg*eigg(3)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudyplus(j1) = g_addudyplus(j1)+g_adtmpmv(i1)*gu3(i1,j1)+g_gu3(i1,j1)*adtmpmv(i1)
          g_adgu3(i1,j1) = g_adgu3(i1,j1)+g_adtmpmv(i1)*dudyplus(j1)+g_dudyplus(j1)*adtmpmv(i1)
          adgu3(i1,j1) = adgu3(i1,j1)+adtmpmv(i1)*dudyplus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    else
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudyminus(j1)*gu3(i1,j1)+g_gu3(i1,j1)*dudyminus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+gu3(i1,j1)*dudyminus(j1)
        end do
      end do
      g_adeigg(3) = g_adeigg(3)+sum(adtmpg*g_tmpmv)
      adeigg(3) = adeigg(3)+sum(adtmpg*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigg(3)*adtmpg
      adtmpmv = adtmpmv+adtmpg*eigg(3)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudyminus(j1) = g_addudyminus(j1)+g_adtmpmv(i1)*gu3(i1,j1)+g_gu3(i1,j1)*adtmpmv(i1)
          g_adgu3(i1,j1) = g_adgu3(i1,j1)+g_adtmpmv(i1)*dudyminus(j1)+g_dudyminus(j1)*adtmpmv(i1)
          adgu3(i1,j1) = adgu3(i1,j1)+adtmpmv(i1)*dudyminus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    endif
    if (eigg(2) .ge. 0.d0) then
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudyplus(j1)*gu2(i1,j1)+g_gu2(i1,j1)*dudyplus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+gu2(i1,j1)*dudyplus(j1)
        end do
      end do
      g_adeigg(2) = g_adeigg(2)+sum(adtmpg*g_tmpmv)
      adeigg(2) = adeigg(2)+sum(adtmpg*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigg(2)*adtmpg
      adtmpmv = adtmpmv+adtmpg*eigg(2)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudyplus(j1) = g_addudyplus(j1)+g_adtmpmv(i1)*gu2(i1,j1)+g_gu2(i1,j1)*adtmpmv(i1)
          g_adgu2(i1,j1) = g_adgu2(i1,j1)+g_adtmpmv(i1)*dudyplus(j1)+g_dudyplus(j1)*adtmpmv(i1)
          adgu2(i1,j1) = adgu2(i1,j1)+adtmpmv(i1)*dudyplus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    else
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudyminus(j1)*gu2(i1,j1)+g_gu2(i1,j1)*dudyminus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+gu2(i1,j1)*dudyminus(j1)
        end do
      end do
      g_adeigg(2) = g_adeigg(2)+sum(adtmpg*g_tmpmv)
      adeigg(2) = adeigg(2)+sum(adtmpg*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigg(2)*adtmpg
      adtmpmv = adtmpmv+adtmpg*eigg(2)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudyminus(j1) = g_addudyminus(j1)+g_adtmpmv(i1)*gu2(i1,j1)+g_gu2(i1,j1)*adtmpmv(i1)
          g_adgu2(i1,j1) = g_adgu2(i1,j1)+g_adtmpmv(i1)*dudyminus(j1)+g_dudyminus(j1)*adtmpmv(i1)
          adgu2(i1,j1) = adgu2(i1,j1)+adtmpmv(i1)*dudyminus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    endif
    if (eigg(1) .ge. 0.d0) then
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudyplus(j1)*gu1(i1,j1)+g_gu1(i1,j1)*dudyplus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+gu1(i1,j1)*dudyplus(j1)
        end do
      end do
      g_adeigg(1) = g_adeigg(1)+sum(adtmpg*g_tmpmv)
      adeigg(1) = adeigg(1)+sum(adtmpg*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigg(1)*adtmpg
      adtmpmv = adtmpmv+adtmpg*eigg(1)
      adtmpg = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudyplus(j1) = g_addudyplus(j1)+g_adtmpmv(i1)*gu1(i1,j1)+g_gu1(i1,j1)*adtmpmv(i1)
          g_adgu1(i1,j1) = g_adgu1(i1,j1)+g_adtmpmv(i1)*dudyplus(j1)+g_dudyplus(j1)*adtmpmv(i1)
          adgu1(i1,j1) = adgu1(i1,j1)+adtmpmv(i1)*dudyplus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    else
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudyminus(j1)*gu1(i1,j1)+g_gu1(i1,j1)*dudyminus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+gu1(i1,j1)*dudyminus(j1)
        end do
      end do
      g_adeigg(1) = g_adeigg(1)+sum(adtmpg*g_tmpmv)
      adeigg(1) = adeigg(1)+sum(adtmpg*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigg(1)*adtmpg
      adtmpmv = adtmpmv+adtmpg*eigg(1)
      adtmpg = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudyminus(j1) = g_addudyminus(j1)+g_adtmpmv(i1)*gu1(i1,j1)+g_gu1(i1,j1)*adtmpmv(i1)
          g_adgu1(i1,j1) = g_adgu1(i1,j1)+g_adtmpmv(i1)*dudyminus(j1)+g_dudyminus(j1)*adtmpmv(i1)
          adgu1(i1,j1) = adgu1(i1,j1)+adtmpmv(i1)*dudyminus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    endif
    g_adu(:,i,j-2) = g_addudyplus/(6.d0*dy)+g_adu(:,i,j-2)
    g_adu(:,i,j-1) = g_addudyplus*((-6)/(6.d0*dy))+g_adu(:,i,j-1)
    g_adu(:,i,j+1) = g_addudyplus*(2/(6.d0*dy))+g_adu(:,i,j+1)
    g_adu(:,i,j) = g_addudyplus*(3/(6.d0*dy))+g_adu(:,i,j)
    g_addudyplus = 0.d0
    g_adu(:,i,j-1) = g_addudyminus*((-2)/(6.d0*dy))+g_adu(:,i,j-1)
    g_adu(:,i,j+2) = (-(g_addudyminus/(6.d0*dy)))+g_adu(:,i,j+2)
    g_adu(:,i,j+1) = g_addudyminus*(6/(6.d0*dy))+g_adu(:,i,j+1)
    g_adu(:,i,j) = g_addudyminus*((-3)/(6.d0*dy))+g_adu(:,i,j)
    g_addudyminus = 0.d0
    g_adisqrt = 0.5d0*g_adgu3(3,3)*vc+g_adisqrt+0.5d0*g_vc*adgu3(3,3)
    adisqrt = adisqrt+0.5d0*adgu3(3,3)*vc
    g_advc = 0.5d0*g_adgu3(3,3)*isqrt+g_advc+0.5d0*g_isqrt*adgu3(3,3)
    advc = advc+0.5d0*adgu3(3,3)*isqrt
    g_adgu3(3,3) = 0.d0
    adgu3(3,3) = 0.d0
    g_adhc = 0.5d0*g_adgu3(3,1)*g*isqrt+g_adhc+0.5d0*g_isqrt*adgu3(3,1)*g
    g_adisqrt = 0.5d0*g_adgu3(3,1)*((-(vc*vc))+g*hc)+g_adisqrt+0.5d0*g_hc*adgu3(3,1)*g-g_vc*adgu3(3,1)*vc
    adisqrt = adisqrt+0.5d0*adgu3(3,1)*((-(vc*vc))+g*hc)
    g_advc = (-(g_adgu3(3,1)*vc*isqrt))+g_advc-g_isqrt*adgu3(3,1)*vc-g_vc*adgu3(3,1)*isqrt
    advc = advc-adgu3(3,1)*vc*isqrt
    g_adgu3(3,1) = 0.d0
    adgu3(3,1) = 0.d0
    g_adisqrt = 0.5d0*g_adgu3(2,3)*uc+g_adisqrt+0.5d0*g_uc*adgu3(2,3)
    adisqrt = adisqrt+0.5d0*adgu3(2,3)*uc
    g_aduc = 0.5d0*g_adgu3(2,3)*isqrt+g_aduc+0.5d0*g_isqrt*adgu3(2,3)
    aduc = aduc+0.5d0*adgu3(2,3)*isqrt
    g_adgu3(2,3) = 0.d0
    adgu3(2,3) = 0.d0
    g_adisqrt = (-(0.5d0*g_adgu3(2,1)*uc*vc))+g_adisqrt-0.5d0*g_uc*adgu3(2,1)*vc-0.5d0*g_vc*adgu3(2,1)*uc
    adisqrt = adisqrt-0.5d0*adgu3(2,1)*uc*vc
    g_aduc = 0.5d0*g_adgu3(2,1)*(1.d0-vc*isqrt)+g_aduc-0.5d0*g_isqrt*adgu3(2,1)*vc-0.5d0*g_vc*adgu3(2,1)*isqrt
    aduc = aduc+0.5d0*adgu3(2,1)*(1.d0-vc*isqrt)
    g_advc = (-(0.5d0*g_adgu3(2,1)*uc*isqrt))+g_advc-0.5d0*g_isqrt*adgu3(2,1)*uc-0.5d0*g_uc*adgu3(2,1)*isqrt
    advc = advc-0.5d0*adgu3(2,1)*uc*isqrt
    g_adgu3(2,1) = 0.d0
    adgu3(2,1) = 0.d0
    g_adisqrt = 0.5d0*g_adgu3(1,3)+g_adisqrt
    adisqrt = adisqrt+0.5d0*adgu3(1,3)
    g_adgu3(1,3) = 0.d0
    adgu3(1,3) = 0.d0
    g_adisqrt = (-(0.5d0*g_adgu3(1,1)*vc))+g_adisqrt-0.5d0*g_vc*adgu3(1,1)
    adisqrt = adisqrt-0.5d0*adgu3(1,1)*vc
    g_advc = (-(0.5d0*g_adgu3(1,1)*isqrt))+g_advc-0.5d0*g_isqrt*adgu3(1,1)
    advc = advc-0.5d0*adgu3(1,1)*isqrt
    g_adgu3(1,1) = 0.d0
    adgu3(1,1) = 0.d0
    g_adgu3(:,:) = 0.d0
    adgu3(:,:) = 0.d0
    g_adhc =g_adhc-g_adisqrt*(1/(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc)))+g_hc*adisqrt*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*&
&sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*&
&sqrt(g*hc)*sqrt(g*hc)))
    g_adisqrt = 0.d0
    adisqrt = 0.d0
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_gu2(:,:) = 0.d0
    gu2(:,:) = 0.d0
    g_gu2(1,1) = 0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt
    gu2(1,1) = 0.5d0*(1.d0+vc*isqrt)
    g_gu2(1,3) = (-0.5d0)*g_isqrt
    gu2(1,3) = -(0.5d0*isqrt)
    g_gu2(2,1) = 0.5d0*g_isqrt*uc*vc+0.5d0*g_uc*(1.d0+vc*isqrt)+0.5d0*g_vc*uc*isqrt
    gu2(2,1) = 0.5d0*uc*(1.d0+vc*isqrt)
    g_gu2(2,3) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    gu2(2,3) = -(0.5d0*uc*isqrt)
    g_gu2(3,1) = (-(0.5d0*g_hc*g*isqrt))+0.5d0*g_isqrt*(vc*vc-g*hc)+g_vc*vc*isqrt
    gu2(3,1) = 0.5d0*(vc*vc-g*hc)*isqrt
    g_adisqrt = (-(0.5d0*g_adgu2(3,3)*vc))+g_adisqrt-0.5d0*g_vc*adgu2(3,3)
    adisqrt = adisqrt-0.5d0*adgu2(3,3)*vc
    g_advc = (-(0.5d0*g_adgu2(3,3)*isqrt))+g_advc-0.5d0*g_isqrt*adgu2(3,3)
    advc = advc-0.5d0*adgu2(3,3)*isqrt
    g_adgu2(3,3) = 0.d0
    adgu2(3,3) = 0.d0
    g_adhc = (-(0.5d0*g_adgu2(3,1)*g*isqrt))+g_adhc-0.5d0*g_isqrt*adgu2(3,1)*g
    g_adisqrt = 0.5d0*g_adgu2(3,1)*(vc*vc-g*hc)+g_adisqrt-0.5d0*g_hc*adgu2(3,1)*g+g_vc*adgu2(3,1)*vc
    adisqrt = adisqrt+0.5d0*adgu2(3,1)*(vc*vc-g*hc)
    g_advc = g_adgu2(3,1)*vc*isqrt+g_advc+g_isqrt*adgu2(3,1)*vc+g_vc*adgu2(3,1)*isqrt
    advc = advc+adgu2(3,1)*vc*isqrt
    g_adgu2(3,1) = 0.d0
    adgu2(3,1) = 0.d0
    g_adisqrt = (-(0.5d0*g_adgu2(2,3)*uc))+g_adisqrt-0.5d0*g_uc*adgu2(2,3)
    adisqrt = adisqrt-0.5d0*adgu2(2,3)*uc
    g_aduc = (-(0.5d0*g_adgu2(2,3)*isqrt))+g_aduc-0.5d0*g_isqrt*adgu2(2,3)
    aduc = aduc-0.5d0*adgu2(2,3)*isqrt
    g_adgu2(2,3) = 0.d0
    adgu2(2,3) = 0.d0
    g_adisqrt = 0.5d0*g_adgu2(2,1)*uc*vc+g_adisqrt+0.5d0*g_uc*adgu2(2,1)*vc+0.5d0*g_vc*adgu2(2,1)*uc
    adisqrt = adisqrt+0.5d0*adgu2(2,1)*uc*vc
    g_aduc = 0.5d0*g_adgu2(2,1)*(1.d0+vc*isqrt)+g_aduc+0.5d0*g_isqrt*adgu2(2,1)*vc+0.5d0*g_vc*adgu2(2,1)*isqrt
    aduc = aduc+0.5d0*adgu2(2,1)*(1.d0+vc*isqrt)
    g_advc = 0.5d0*g_adgu2(2,1)*uc*isqrt+g_advc+0.5d0*g_isqrt*adgu2(2,1)*uc+0.5d0*g_uc*adgu2(2,1)*isqrt
    advc = advc+0.5d0*adgu2(2,1)*uc*isqrt
    g_adgu2(2,1) = 0.d0
    adgu2(2,1) = 0.d0
    g_adisqrt = (-0.5d0)*g_adgu2(1,3)+g_adisqrt
    adisqrt = adisqrt-0.5d0*adgu2(1,3)
    g_adgu2(1,3) = 0.d0
    adgu2(1,3) = 0.d0
    g_adisqrt = 0.5d0*g_adgu2(1,1)*vc+g_adisqrt+0.5d0*g_vc*adgu2(1,1)
    adisqrt = adisqrt+0.5d0*adgu2(1,1)*vc
    g_advc = 0.5d0*g_adgu2(1,1)*isqrt+g_advc+0.5d0*g_isqrt*adgu2(1,1)
    advc = advc+0.5d0*adgu2(1,1)*isqrt
    g_adgu2(1,1) = 0.d0
    adgu2(1,1) = 0.d0
    g_adgu2(:,:) = 0.d0
    adgu2(:,:) = 0.d0
    g_adhc = g_adhc-g_adisqrt*(1/(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc)))+g_hc*adisqrt*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*&
&sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*&
&sqrt(g*hc)*sqrt(g*hc)))
    g_adisqrt = 0.d0
    adisqrt = 0.d0
    g_adgu1(2,2) = 0.d0
    adgu1(2,2) = 0.d0
    g_aduc = (-g_adgu1(2,1))+g_aduc
    aduc = aduc-adgu1(2,1)
    g_adgu1(2,1) = 0.d0
    adgu1(2,1) = 0.d0
    g_adgu1(:,:) = 0.d0
    adgu1(:,:) = 0.d0
    g_adu(3,i,j) = g_adeigg(3)/u(1,i,j)+g_adu(3,i,j)-g_u(1,i,j)*(adeigg(3)/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = g_adeigg(3)*((-(u(3,i,j)/(u(1,i,j)*u(1,i,j))))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)+g_adu(1,i,j)-g_u(3,i,j)*&
&(adeigg(3)/(u(1,i,j)*u(1,i,j)))+g_u(1,i,j)*adeigg(3)*(2*u(3,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))-2*1./&
&(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*sqrt(gacc*u(1,i,j)))*gacc)
    g_adeigg(3) = 0.d0
    adeigg(3) = 0.d0
    g_adu(3,i,j) = g_adeigg(2)/u(1,i,j)+g_adu(3,i,j)-g_u(1,i,j)*(adeigg(2)/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = (-(g_adeigg(2)*(u(3,i,j)/(u(1,i,j)*u(1,i,j))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)))+g_adu(1,i,j)-g_u(3,i,j)*&
&(adeigg(2)/(u(1,i,j)*u(1,i,j)))+g_u(1,i,j)*adeigg(2)*(2*u(3,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))+2*1./&
&(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*sqrt(gacc*u(1,i,j)))*gacc)
    g_adeigg(2) = 0.d0
    adeigg(2) = 0.d0
    g_adu(3,i,j) = g_adeigg(1)/u(1,i,j)+g_adu(3,i,j)-g_u(1,i,j)*(adeigg(1)/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = &
      & (-(g_adeigg(1)*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))))+g_adu(1,i,j)-g_u(3,i,j)*(adeigg(1)/(u(1,i,j)*u(1,i,j)))+g_u(1,i,&
&j)*adeigg(1)*(2*u(3,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j)))
    g_adeigg(1) = 0.d0
    adeigg(1) = 0.d0
    if (eigf(3) .ge. 0.d0) then
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudxplus(j1)*fu3(i1,j1)+g_fu3(i1,j1)*dudxplus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+fu3(i1,j1)*dudxplus(j1)
        end do
      end do
      g_adeigf(3) = g_adeigf(3)+sum(adtmpf*g_tmpmv)
      adeigf(3) = adeigf(3)+sum(adtmpf*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigf(3)*adtmpf
      adtmpmv = adtmpmv+adtmpf*eigf(3)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudxplus(j1) = g_addudxplus(j1)+g_adtmpmv(i1)*fu3(i1,j1)+g_fu3(i1,j1)*adtmpmv(i1)
          g_adfu3(i1,j1) = g_adfu3(i1,j1)+g_adtmpmv(i1)*dudxplus(j1)+g_dudxplus(j1)*adtmpmv(i1)
          adfu3(i1,j1) = adfu3(i1,j1)+adtmpmv(i1)*dudxplus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    else
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudxminus(j1)*fu3(i1,j1)+g_fu3(i1,j1)*dudxminus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+fu3(i1,j1)*dudxminus(j1)
        end do
      end do
      g_adeigf(3) = g_adeigf(3)+sum(adtmpf*g_tmpmv)
      adeigf(3) = adeigf(3)+sum(adtmpf*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigf(3)*adtmpf
      adtmpmv = adtmpmv+adtmpf*eigf(3)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudxminus(j1) = g_addudxminus(j1)+g_adtmpmv(i1)*fu3(i1,j1)+g_fu3(i1,j1)*adtmpmv(i1)
          g_adfu3(i1,j1) = g_adfu3(i1,j1)+g_adtmpmv(i1)*dudxminus(j1)+g_dudxminus(j1)*adtmpmv(i1)
          adfu3(i1,j1) = adfu3(i1,j1)+adtmpmv(i1)*dudxminus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    endif
    if (eigf(2) .ge. 0.d0) then
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudxplus(j1)*fu2(i1,j1)+g_fu2(i1,j1)*dudxplus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+fu2(i1,j1)*dudxplus(j1)
        end do
      end do
      g_adeigf(2) = g_adeigf(2)+sum(adtmpf*g_tmpmv)
      adeigf(2) = adeigf(2)+sum(adtmpf*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigf(2)*adtmpf
      adtmpmv = adtmpmv+adtmpf*eigf(2)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudxplus(j1) = g_addudxplus(j1)+g_adtmpmv(i1)*fu2(i1,j1)+g_fu2(i1,j1)*adtmpmv(i1)
          g_adfu2(i1,j1) = g_adfu2(i1,j1)+g_adtmpmv(i1)*dudxplus(j1)+g_dudxplus(j1)*adtmpmv(i1)
          adfu2(i1,j1) = adfu2(i1,j1)+adtmpmv(i1)*dudxplus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    else
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudxminus(j1)*fu2(i1,j1)+g_fu2(i1,j1)*dudxminus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+fu2(i1,j1)*dudxminus(j1)
        end do
      end do
      g_adeigf(2) = g_adeigf(2)+sum(adtmpf*g_tmpmv)
      adeigf(2) = adeigf(2)+sum(adtmpf*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigf(2)*adtmpf
      adtmpmv = adtmpmv+adtmpf*eigf(2)
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudxminus(j1) = g_addudxminus(j1)+g_adtmpmv(i1)*fu2(i1,j1)+g_fu2(i1,j1)*adtmpmv(i1)
          g_adfu2(i1,j1) = g_adfu2(i1,j1)+g_adtmpmv(i1)*dudxminus(j1)+g_dudxminus(j1)*adtmpmv(i1)
          adfu2(i1,j1) = adfu2(i1,j1)+adtmpmv(i1)*dudxminus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    endif
    if (eigf(1) .ge. 0.d0) then
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudxplus(j1)*fu1(i1,j1)+g_fu1(i1,j1)*dudxplus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+fu1(i1,j1)*dudxplus(j1)
        end do
      end do
      g_adeigf(1) = g_adeigf(1)+sum(adtmpf*g_tmpmv)
      adeigf(1) = adeigf(1)+sum(adtmpf*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigf(1)*adtmpf
      adtmpmv = adtmpmv+adtmpf*eigf(1)
      adtmpf = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudxplus(j1) = g_addudxplus(j1)+g_adtmpmv(i1)*fu1(i1,j1)+g_fu1(i1,j1)*adtmpmv(i1)
          g_adfu1(i1,j1) = g_adfu1(i1,j1)+g_adtmpmv(i1)*dudxplus(j1)+g_dudxplus(j1)*adtmpmv(i1)
          adfu1(i1,j1) = adfu1(i1,j1)+adtmpmv(i1)*dudxplus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    else
      g_tmpmv(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_tmpmv(i1) = g_dudxminus(j1)*fu1(i1,j1)+g_fu1(i1,j1)*dudxminus(j1)+g_tmpmv(i1)
          tmpmv(i1) = tmpmv(i1)+fu1(i1,j1)*dudxminus(j1)
        end do
      end do
      g_adeigf(1) = g_adeigf(1)+sum(adtmpf*g_tmpmv)
      adeigf(1) = adeigf(1)+sum(adtmpf*tmpmv)
      g_adtmpmv = g_adtmpmv+g_eigf(1)*adtmpf
      adtmpmv = adtmpmv+adtmpf*eigf(1)
      adtmpf = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_addudxminus(j1) = g_addudxminus(j1)+g_adtmpmv(i1)*fu1(i1,j1)+g_fu1(i1,j1)*adtmpmv(i1)
          g_adfu1(i1,j1) = g_adfu1(i1,j1)+g_adtmpmv(i1)*dudxminus(j1)+g_dudxminus(j1)*adtmpmv(i1)
          adfu1(i1,j1) = adfu1(i1,j1)+adtmpmv(i1)*dudxminus(j1)
        end do
      end do
      g_adtmpmv(:) = 0.d0
      adtmpmv(:) = 0.d0
    endif
    g_adu(:,i-2,j) = g_addudxplus/(6.d0*dx)+g_adu(:,i-2,j)
    g_adu(:,i-1,j) = g_addudxplus*((-6)/(6.d0*dx))+g_adu(:,i-1,j)
    g_adu(:,i+1,j) = g_addudxplus*(2/(6.d0*dx))+g_adu(:,i+1,j)
    g_adu(:,i,j) = g_addudxplus*(3/(6.d0*dx))+g_adu(:,i,j)
    g_addudxplus = 0.d0
    g_adu(:,i-1,j) = g_addudxminus*((-2)/(6.d0*dx))+g_adu(:,i-1,j)
    g_adu(:,i+2,j) = (-(g_addudxminus/(6.d0*dx)))+g_adu(:,i+2,j)
    g_adu(:,i+1,j) = g_addudxminus*(6/(6.d0*dx))+g_adu(:,i+1,j)
    g_adu(:,i,j) = g_addudxminus*((-3)/(6.d0*dx))+g_adu(:,i,j)
    g_addudxminus = 0.d0
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_fu3(:,:) = 0.d0
    fu3(:,:) = 0.d0
    g_fu3(1,1) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    fu3(1,1) = 0.5d0-0.5d0*uc*isqrt
    g_fu3(1,2) = 0.5d0*g_isqrt
    fu3(1,2) = 0.5d0*isqrt
    g_fu3(2,1) = 0.5d0*g_hc*g*isqrt+0.5d0*g_isqrt*(g*hc-uc*uc)-g_uc*uc*isqrt
    fu3(2,1) = 0.5d0*(g*hc-uc*uc)*isqrt
    g_fu3(2,2) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    fu3(2,2) = 0.5d0*(1.d0+uc*isqrt)
    g_fu3(3,1) = (-(0.5d0*g_isqrt*uc*vc))-0.5d0*g_uc*isqrt*vc+0.5d0*g_vc*(1.d0-uc*isqrt)
    fu3(3,1) = 0.5d0*(1.d0-uc*isqrt)*vc
    g_adisqrt = 0.5d0*g_adfu3(3,2)*vc+g_adisqrt+0.5d0*g_vc*adfu3(3,2)
    adisqrt = adisqrt+0.5d0*adfu3(3,2)*vc
    g_advc = 0.5d0*g_adfu3(3,2)*isqrt+g_advc+0.5d0*g_isqrt*adfu3(3,2)
    advc = advc+0.5d0*adfu3(3,2)*isqrt
    g_adfu3(3,2) = 0.d0
    adfu3(3,2) = 0.d0
    g_adisqrt = (-(0.5d0*g_adfu3(3,1)*uc*vc))+g_adisqrt-0.5d0*g_uc*adfu3(3,1)*vc-0.5d0*g_vc*adfu3(3,1)*uc
    adisqrt = adisqrt-0.5d0*adfu3(3,1)*uc*vc
    g_aduc = (-(0.5d0*g_adfu3(3,1)*isqrt*vc))+g_aduc-0.5d0*g_isqrt*adfu3(3,1)*vc-0.5d0*g_vc*adfu3(3,1)*isqrt
    aduc = aduc-0.5d0*adfu3(3,1)*isqrt*vc
    g_advc = 0.5d0*g_adfu3(3,1)*(1.d0-uc*isqrt)+g_advc-0.5d0*g_isqrt*adfu3(3,1)*uc-0.5d0*g_uc*adfu3(3,1)*isqrt
    advc = advc+0.5d0*adfu3(3,1)*(1.d0-uc*isqrt)
    g_adfu3(3,1) = 0.d0
    adfu3(3,1) = 0.d0
    g_adisqrt = 0.5d0*g_adfu3(2,2)*uc+g_adisqrt+0.5d0*g_uc*adfu3(2,2)
    adisqrt = adisqrt+0.5d0*adfu3(2,2)*uc
    g_aduc = 0.5d0*g_adfu3(2,2)*isqrt+g_aduc+0.5d0*g_isqrt*adfu3(2,2)
    aduc = aduc+0.5d0*adfu3(2,2)*isqrt
    g_adfu3(2,2) = 0.d0
    adfu3(2,2) = 0.d0
    g_adhc = 0.5d0*g_adfu3(2,1)*g*isqrt+g_adhc+0.5d0*g_isqrt*adfu3(2,1)*g
    g_adisqrt = 0.5d0*g_adfu3(2,1)*(g*hc-uc*uc)+g_adisqrt+0.5d0*g_hc*adfu3(2,1)*g-g_uc*adfu3(2,1)*uc
    adisqrt = adisqrt+0.5d0*adfu3(2,1)*(g*hc-uc*uc)
    g_aduc = (-(g_adfu3(2,1)*uc*isqrt))+g_aduc-g_isqrt*adfu3(2,1)*uc-g_uc*adfu3(2,1)*isqrt
    aduc = aduc-adfu3(2,1)*uc*isqrt
    g_adfu3(2,1) = 0.d0
    adfu3(2,1) = 0.d0
    g_adisqrt = 0.5d0*g_adfu3(1,2)+g_adisqrt
    adisqrt = adisqrt+0.5d0*adfu3(1,2)
    g_adfu3(1,2) = 0.d0
    adfu3(1,2) = 0.d0
    g_adisqrt = (-(0.5d0*g_adfu3(1,1)*uc))+g_adisqrt-0.5d0*g_uc*adfu3(1,1)
    adisqrt = adisqrt-0.5d0*adfu3(1,1)*uc
    g_aduc = (-(0.5d0*g_adfu3(1,1)*isqrt))+g_aduc-0.5d0*g_isqrt*adfu3(1,1)
    aduc = aduc-0.5d0*adfu3(1,1)*isqrt
    g_adfu3(1,1) = 0.d0
    adfu3(1,1) = 0.d0
    g_adfu3(:,:) = 0.d0
    adfu3(:,:) = 0.d0
    g_adhc = g_adhc-g_adisqrt*(1/(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc)))+g_hc*adisqrt*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*&
&sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*&
&sqrt(g*hc)*sqrt(g*hc)))
    g_adisqrt = 0.d0
    adisqrt = 0.d0
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_fu2(:,:) = 0.d0
    fu2(:,:) = 0.d0
    g_fu2(1,1) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    fu2(1,1) = 0.5d0*(1.d0+uc*isqrt)
    g_fu2(1,2) = (-0.5d0)*g_isqrt
    fu2(1,2) = -(0.5d0*isqrt)
    g_fu2(2,1) = (-(0.5d0*g_hc*g*isqrt))+0.5d0*g_isqrt*(uc*uc-g*hc)+g_uc*uc*isqrt
    fu2(2,1) = 0.5d0*(uc*uc-g*hc)*isqrt
    g_fu2(2,2) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    fu2(2,2) = 0.5d0-0.5d0*uc*isqrt
    g_fu2(3,1) = 0.5d0*g_isqrt*uc*vc+0.5d0*g_uc*isqrt*vc+0.5d0*g_vc*(1.d0+uc*isqrt)
    fu2(3,1) = 0.5d0*(1.d0+uc*isqrt)*vc
    g_adisqrt = (-(0.5d0*g_adfu2(3,2)*vc))+g_adisqrt-0.5d0*g_vc*adfu2(3,2)
    adisqrt = adisqrt-0.5d0*adfu2(3,2)*vc
    g_advc = (-(0.5d0*g_adfu2(3,2)*isqrt))+g_advc-0.5d0*g_isqrt*adfu2(3,2)
    advc = advc-0.5d0*adfu2(3,2)*isqrt
    g_adfu2(3,2) = 0.d0
    adfu2(3,2) = 0.d0
    g_adisqrt = 0.5d0*g_adfu2(3,1)*uc*vc+g_adisqrt+0.5d0*g_uc*adfu2(3,1)*vc+0.5d0*g_vc*adfu2(3,1)*uc
    adisqrt = adisqrt+0.5d0*adfu2(3,1)*uc*vc
    g_aduc = 0.5d0*g_adfu2(3,1)*isqrt*vc+g_aduc+0.5d0*g_isqrt*adfu2(3,1)*vc+0.5d0*g_vc*adfu2(3,1)*isqrt
    aduc = aduc+0.5d0*adfu2(3,1)*isqrt*vc
    g_advc = 0.5d0*g_adfu2(3,1)*(1.d0+uc*isqrt)+g_advc+0.5d0*g_isqrt*adfu2(3,1)*uc+0.5d0*g_uc*adfu2(3,1)*isqrt
    advc = advc+0.5d0*adfu2(3,1)*(1.d0+uc*isqrt)
    g_adfu2(3,1) = 0.d0
    adfu2(3,1) = 0.d0
    g_adisqrt = (-(0.5d0*g_adfu2(2,2)*uc))+g_adisqrt-0.5d0*g_uc*adfu2(2,2)
    adisqrt = adisqrt-0.5d0*adfu2(2,2)*uc
    g_aduc = (-(0.5d0*g_adfu2(2,2)*isqrt))+g_aduc-0.5d0*g_isqrt*adfu2(2,2)
    aduc = aduc-0.5d0*adfu2(2,2)*isqrt
    g_adfu2(2,2) = 0.d0
    adfu2(2,2) = 0.d0
    g_adhc = (-(0.5d0*g_adfu2(2,1)*g*isqrt))+g_adhc-0.5d0*g_isqrt*adfu2(2,1)*g
    g_adisqrt = 0.5d0*g_adfu2(2,1)*(uc*uc-g*hc)+g_adisqrt-0.5d0*g_hc*adfu2(2,1)*g+g_uc*adfu2(2,1)*uc
    adisqrt = adisqrt+0.5d0*adfu2(2,1)*(uc*uc-g*hc)
    g_aduc = g_adfu2(2,1)*uc*isqrt+g_aduc+g_isqrt*adfu2(2,1)*uc+g_uc*adfu2(2,1)*isqrt
    aduc = aduc+adfu2(2,1)*uc*isqrt
    g_adfu2(2,1) = 0.d0
    adfu2(2,1) = 0.d0
    g_adisqrt = (-0.5d0)*g_adfu2(1,2)+g_adisqrt
    adisqrt = adisqrt-0.5d0*adfu2(1,2)
    g_adfu2(1,2) = 0.d0
    adfu2(1,2) = 0.d0
    g_adisqrt = 0.5d0*g_adfu2(1,1)*uc+g_adisqrt+0.5d0*g_uc*adfu2(1,1)
    adisqrt = adisqrt+0.5d0*adfu2(1,1)*uc
    g_aduc = 0.5d0*g_adfu2(1,1)*isqrt+g_aduc+0.5d0*g_isqrt*adfu2(1,1)
    aduc = aduc+0.5d0*adfu2(1,1)*isqrt
    g_adfu2(1,1) = 0.d0
    adfu2(1,1) = 0.d0
    g_adfu2(:,:) = 0.d0
    adfu2(:,:) = 0.d0
    g_adhc = g_adhc-g_adisqrt*(1/(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc)))+g_hc*adisqrt*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*&
&sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*&
&sqrt(g*hc)*sqrt(g*hc)))
    g_adisqrt = 0.d0
    adisqrt = 0.d0
    g_adfu1(3,3) = 0.d0
    adfu1(3,3) = 0.d0
    g_advc = (-g_adfu1(3,1))+g_advc
    advc = advc-adfu1(3,1)
    g_adfu1(3,1) = 0.d0
    adfu1(3,1) = 0.d0
    g_adfu1(:,:) = 0.d0
    adfu1(:,:) = 0.d0
    g_adu(3,i,j) = g_adu(3,i,j)+g_advc/u(1,i,j)-g_u(1,i,j)*(advc/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = &
      & g_adu(1,i,j)-g_advc*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))-g_u(3,i,j)*(advc/(u(1,i,j)*u(1,i,j)))+g_u(1,i,j)*advc*(2*u(3,&
&i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j)))
    g_advc = 0.d0
    advc = 0.d0
    g_adu(2,i,j) = g_adu(2,i,j)+g_aduc/u(1,i,j)-g_u(1,i,j)*(aduc/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = &
      & g_adu(1,i,j)-g_aduc*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))-g_u(2,i,j)*(aduc/(u(1,i,j)*u(1,i,j)))+g_u(1,i,j)*aduc*(2*u(2,&
&i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j)))
    g_aduc = 0.d0
    aduc = 0.d0
    g_adu(1,i,j) = g_adhc+g_adu(1,i,j)
    g_adhc = 0.d0
    g_adu(2,i,j) = g_adeigf(3)/u(1,i,j)+g_adu(2,i,j)-g_u(1,i,j)*(adeigf(3)/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = g_adeigf(3)*((-(u(2,i,j)/(u(1,i,j)*u(1,i,j))))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)+g_adu(1,i,j)-g_u(2,i,j)*&
&(adeigf(3)/(u(1,i,j)*u(1,i,j)))+g_u(1,i,j)*adeigf(3)*(2*u(2,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))-2*1./&
&(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*sqrt(gacc*u(1,i,j)))*gacc)
    g_adeigf(3) = 0.d0
    adeigf(3) = 0.d0
    g_adu(2,i,j) = g_adeigf(2)/u(1,i,j)+g_adu(2,i,j)-g_u(1,i,j)*(adeigf(2)/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = (-(g_adeigf(2)*(u(2,i,j)/(u(1,i,j)*u(1,i,j))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)))+g_adu(1,i,j)-g_u(2,i,j)*&
&(adeigf(2)/(u(1,i,j)*u(1,i,j)))+g_u(1,i,j)*adeigf(2)*(2*u(2,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))+2*1./&
&(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*sqrt(gacc*u(1,i,j)))*gacc)
    g_adeigf(2) = 0.d0
    adeigf(2) = 0.d0
    g_adu(2,i,j) = g_adeigf(1)/u(1,i,j)+g_adu(2,i,j)-g_u(1,i,j)*(adeigf(1)/(u(1,i,j)*u(1,i,j)))
    g_adu(1,i,j) = &
      & (-(g_adeigf(1)*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))))+g_adu(1,i,j)-g_u(2,i,j)*(adeigf(1)/(u(1,i,j)*u(1,i,j)))+g_u(1,i,&
&j)*adeigf(1)*(2*u(2,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j)))
    g_adeigf(1) = 0.d0
    adeigf(1) = 0.d0
  end do
end do
g_adu(:,3:mx+2,4) = g_adu(:,3:mx+2,my+4)+g_adu(:,3:mx+2,4)
g_adu(:,3:mx+2,my+4) = 0.d0
g_adu(:,3:mx+2,3) = g_adu(:,3:mx+2,my+3)+g_adu(:,3:mx+2,3)
g_adu(:,3:mx+2,my+3) = 0.d0
g_adu(:,3:mx+2,my+1) = g_adu(:,3:mx+2,my+1)+g_adu(:,3:mx+2,1)
g_adu(:,3:mx+2,1) = 0.d0
g_adu(:,3:mx+2,my+2) = g_adu(:,3:mx+2,my+2)+g_adu(:,3:mx+2,2)
g_adu(:,3:mx+2,2) = 0.d0
g_adu(:,4,3:my+2) = g_adu(:,mx+4,3:my+2)+g_adu(:,4,3:my+2)
g_adu(:,mx+4,3:my+2) = 0.d0
g_adu(:,3,3:my+2) = g_adu(:,mx+3,3:my+2)+g_adu(:,3,3:my+2)
g_adu(:,mx+3,3:my+2) = 0.d0
g_adu(:,mx+1,3:my+2) = g_adu(:,mx+1,3:my+2)+g_adu(:,1,3:my+2)
g_adu(:,1,3:my+2) = 0.d0
g_adu(:,mx+2,3:my+2) = g_adu(:,mx+2,3:my+2)+g_adu(:,2,3:my+2)
g_adu(:,2,3:my+2) = 0.d0

end subroutine g_adcompute_f

subroutine g_g_compute_f( u, g_u, g_v, g_g_f )
! g_g_f = d(Jac^T * g_u)* g_v = (H x g_v) * g_u
!***************************************************************
!***************************************************************
!** This routine was generated by the                         **
!** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
!***************************************************************
!***************************************************************
!==============================================
! referencing used modules
!==============================================
use swe2dxy_parameters

!==============================================
! all entries are defined explicitly
!==============================================
implicit none

!==============================================
! define arguments
!==============================================
double precision g_g_f(3,mx+4,my+4)
double precision g_u(3,mx+4,my+4)
double precision g_v(3,mx+4,my+4)
double precision u(3,mx+4,my+4)

!==============================================
! define local variables
!==============================================
double precision dudxminus(3)
double precision dudxplus(3)
double precision dudyminus(3)
double precision dudyplus(3)
double precision eigf(3)
double precision eigg(3)
double precision fu1(3,3)
double precision fu2(3,3)
double precision fu3(3,3)
double precision g
double precision g_dudxminus(3)
double precision g_dudxminut(3)
double precision g_dudxplus(3)
double precision g_dudxplut(3)
double precision g_dudyminus(3)
double precision g_dudyminut(3)
double precision g_dudyplus(3)
double precision g_dudyplut(3)
double precision g_eigf(3)
double precision g_eigg(3)
double precision g_eigh(3)
double precision g_eigi(3)
double precision g_fu1(3,3)
double precision g_fu2(3,3)
double precision g_fu3(3,3)
double precision g_fu4(3,3)
double precision g_fu5(3,3)
double precision g_fu6(3,3)
double precision g_g_eigf(3)
double precision g_g_eigg(3)
double precision g_g_fu1(3,3)
double precision g_g_fu2(3,3)
double precision g_g_fu3(3,3)
double precision g_g_gu1(3,3)
double precision g_g_gu2(3,3)
double precision g_g_gu3(3,3)
double precision g_g_isqrt
double precision g_g_rhs(3)
double precision g_g_tmpf(3)
double precision g_g_tmpg(3)
double precision g_g_tmpmv(3)
double precision g_g_uc
double precision g_g_vc
double precision g_gu1(3,3)
double precision g_gu2(3,3)
double precision g_gu3(3,3)
double precision g_gu4(3,3)
double precision g_gu5(3,3)
double precision g_gu6(3,3)
double precision g_hc
double precision g_hd
double precision g_isqrt
double precision g_isqru
double precision g_tmpmv(3)
double precision g_tmpmw(3)
double precision g_uc
double precision g_ud
double precision g_vc
double precision g_vd
double precision gu1(3,3)
double precision gu2(3,3)
double precision gu3(3,3)
double precision hc
integer i
integer i1
double precision isqrt
integer j
integer j1
double precision tmpmv(3)
double precision uc
double precision vc

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------
g_u(:,2,3:my+2) = g_u(:,mx+2,3:my+2)
g_v(:,2,3:my+2) = g_v(:,mx+2,3:my+2)
u(:,2,3:my+2) = u(:,mx+2,3:my+2)
g_u(:,1,3:my+2) = g_u(:,mx+1,3:my+2)
g_v(:,1,3:my+2) = g_v(:,mx+1,3:my+2)
u(:,1,3:my+2) = u(:,mx+1,3:my+2)
g_u(:,mx+3,3:my+2) = g_u(:,3,3:my+2)
g_v(:,mx+3,3:my+2) = g_v(:,3,3:my+2)
u(:,mx+3,3:my+2) = u(:,3,3:my+2)
g_u(:,mx+4,3:my+2) = g_u(:,4,3:my+2)
g_v(:,mx+4,3:my+2) = g_v(:,4,3:my+2)
u(:,mx+4,3:my+2) = u(:,4,3:my+2)
g_u(:,3:mx+2,2) = g_u(:,3:mx+2,my+2)
g_v(:,3:mx+2,2) = g_v(:,3:mx+2,my+2)
u(:,3:mx+2,2) = u(:,3:mx+2,my+2)
g_u(:,3:mx+2,1) = g_u(:,3:mx+2,my+1)
g_v(:,3:mx+2,1) = g_v(:,3:mx+2,my+1)
u(:,3:mx+2,1) = u(:,3:mx+2,my+1)
g_u(:,3:mx+2,my+3) = g_u(:,3:mx+2,3)
g_v(:,3:mx+2,my+3) = g_v(:,3:mx+2,3)
u(:,3:mx+2,my+3) = u(:,3:mx+2,3)
g_u(:,3:mx+2,my+4) = g_u(:,3:mx+2,4)
g_v(:,3:mx+2,my+4) = g_v(:,3:mx+2,4)
u(:,3:mx+2,my+4) = u(:,3:mx+2,4)
do j = 3, my+2
  do i = 3, mx+2
    g_g_eigf(1) = &
      & (-(g_v(2,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))+g_v(1,i,j)*((-(g_u(2,i,j)/(u(1,i,j)*u(1,i,j))))+g_u(1,i,j)*(2*&
&u(2,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))))
    g_eigf(1) = g_u(2,i,j)/u(1,i,j)-g_u(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))
    g_eigh(1) = g_v(2,i,j)/u(1,i,j)-g_v(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))
    eigf(1) = u(2,i,j)/u(1,i,j)
    g_g_eigf(2) = &
      & (-(g_v(2,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))-g_v(1,i,j)*(g_u(2,i,j)/(u(1,i,j)*u(1,i,j))-g_u(1,i,j)*(2*u(2,i,&
&j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))+2*1./(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*sqrt(gacc*u(1,i,&
&j)))*gacc))
    g_eigf(2) = g_u(2,i,j)/u(1,i,j)-g_u(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    g_eigh(2) = g_v(2,i,j)/u(1,i,j)-g_v(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigf(2) = u(2,i,j)/u(1,i,j)-sqrt(gacc*u(1,i,j))
    g_g_eigf(3) = &
      & (-(g_v(2,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))+g_v(1,i,j)*((-(g_u(2,i,j)/(u(1,i,j)*u(1,i,j))))+g_u(1,i,j)*(2*&
&u(2,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))-2*1./(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*&
&sqrt(gacc*u(1,i,j)))*gacc))
    g_eigf(3) = g_u(2,i,j)/u(1,i,j)+g_u(1,i,j)*((-(u(2,i,j)/(u(1,i,j)*u(1,i,j))))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    g_eigh(3) = g_v(2,i,j)/u(1,i,j)+g_v(1,i,j)*((-(u(2,i,j)/(u(1,i,j)*u(1,i,j))))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigf(3) = u(2,i,j)/u(1,i,j)+sqrt(gacc*u(1,i,j))
    g_hc = g_u(1,i,j)
    g_hd = g_v(1,i,j)
    hc = u(1,i,j)
    g_g_uc = &
      & (-(g_v(2,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))+g_v(1,i,j)*((-(g_u(2,i,j)/(u(1,i,j)*u(1,i,j))))+g_u(1,i,j)*(2*u(2,i,&
&j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))))
    g_uc = g_u(2,i,j)/u(1,i,j)-g_u(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))
    g_ud = g_v(2,i,j)/u(1,i,j)-g_v(1,i,j)*(u(2,i,j)/(u(1,i,j)*u(1,i,j)))
    uc = u(2,i,j)/u(1,i,j)
    g_g_vc = &
      & (-(g_v(3,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))+g_v(1,i,j)*((-(g_u(3,i,j)/(u(1,i,j)*u(1,i,j))))+g_u(1,i,j)*(2*u(3,i,&
&j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))))
    g_vc = g_u(3,i,j)/u(1,i,j)-g_u(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))
    g_vd = g_v(3,i,j)/u(1,i,j)-g_v(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))
    vc = u(3,i,j)/u(1,i,j)
    g = gacc
    g_g_fu1(:,:) = 0.d0
    g_fu1(:,:) = 0.d0
    g_fu4(:,:) = 0.d0
    fu1(:,:) = 0.d0
    g_g_fu1(3,1) = -g_g_vc
    g_fu1(3,1) = -g_vc
    g_fu4(3,1) = -g_vd
    fu1(3,1) = -vc
    g_g_fu1(3,3) = 0.d0
    g_fu1(3,3) = 0.d0
    g_fu4(3,3) = 0.d0
    fu1(3,3) = 1.d0
    g_g_isqrt = &
      & g_hd*g_hc*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*&
&sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)))
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    g_isqru = -(g_hd*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_g_fu2(:,:) = 0.d0
    g_fu2(:,:) = 0.d0
    g_fu5(:,:) = 0.d0
    fu2(:,:) = 0.d0
    g_g_fu2(1,1) = 0.5d0*g_g_isqrt*uc+0.5d0*g_g_uc*isqrt+0.5d0*g_isqru*g_uc+0.5d0*g_ud*g_isqrt
    g_fu2(1,1) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    g_fu5(1,1) = 0.5d0*g_isqru*uc+0.5d0*g_ud*isqrt
    fu2(1,1) = 0.5d0*(1.d0+uc*isqrt)
    g_g_fu2(1,2) = (-0.5d0)*g_g_isqrt
    g_fu2(1,2) = (-0.5d0)*g_isqrt
    g_fu5(1,2) = (-0.5d0)*g_isqru
    fu2(1,2) = -(0.5d0*isqrt)
    g_g_fu2(2,1) = 0.5d0*g_g_isqrt*(uc*uc-g*hc)+g_g_uc*uc*isqrt-0.5d0*g_hd*g_isqrt*g+g_isqru*((-(0.5d0*g_hc*g))+g_uc*uc)+g_ud*&
&(g_isqrt*uc+g_uc*isqrt)
    g_fu2(2,1) = (-(0.5d0*g_hc*g*isqrt))+0.5d0*g_isqrt*(uc*uc-g*hc)+g_uc*uc*isqrt
    g_fu5(2,1) = (-(0.5d0*g_hd*g*isqrt))+0.5d0*g_isqru*(uc*uc-g*hc)+g_ud*uc*isqrt
    fu2(2,1) = 0.5d0*(uc*uc-g*hc)*isqrt
    g_g_fu2(2,2) = (-(0.5d0*g_g_isqrt*uc))-0.5d0*g_g_uc*isqrt-0.5d0*g_isqru*g_uc-0.5d0*g_ud*g_isqrt
    g_fu2(2,2) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    g_fu5(2,2) = (-(0.5d0*g_isqru*uc))-0.5d0*g_ud*isqrt
    fu2(2,2) = 0.5d0-0.5d0*uc*isqrt
    g_g_fu2(3,1) = &
      & 0.5d0*g_g_isqrt*uc*vc+0.5d0*g_g_uc*isqrt*vc+0.5d0*g_g_vc*(1.d0+uc*isqrt)+g_isqru*(0.5d0*g_uc*vc+0.5d0*g_vc*uc)+&
&g_ud*(0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt)+g_vd*(0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt)
    g_fu2(3,1) = 0.5d0*g_isqrt*uc*vc+0.5d0*g_uc*isqrt*vc+0.5d0*g_vc*(1.d0+uc*isqrt)
    g_fu5(3,1) = 0.5d0*g_isqru*uc*vc+0.5d0*g_ud*isqrt*vc+0.5d0*g_vd*(1.d0+uc*isqrt)
    fu2(3,1) = 0.5d0*(1.d0+uc*isqrt)*vc
    g_g_fu2(3,2) = (-(0.5d0*g_g_isqrt*vc))-0.5d0*g_g_vc*isqrt-0.5d0*g_isqru*g_vc-0.5d0*g_vd*g_isqrt
    g_fu2(3,2) = (-(0.5d0*g_isqrt*vc))-0.5d0*g_vc*isqrt
    g_fu5(3,2) = (-(0.5d0*g_isqru*vc))-0.5d0*g_vd*isqrt
    fu2(3,2) = -(0.5d0*vc*isqrt)
    g_g_isqrt = &
      & g_hd*g_hc*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*&
&sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)))
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    g_isqru = -(g_hd*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_g_fu3(:,:) = 0.d0
    g_fu3(:,:) = 0.d0
    g_fu6(:,:) = 0.d0
    fu3(:,:) = 0.d0
    g_g_fu3(1,1) = (-(0.5d0*g_g_isqrt*uc))-0.5d0*g_g_uc*isqrt-0.5d0*g_isqru*g_uc-0.5d0*g_ud*g_isqrt
    g_fu3(1,1) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    g_fu6(1,1) = (-(0.5d0*g_isqru*uc))-0.5d0*g_ud*isqrt
    fu3(1,1) = 0.5d0-0.5d0*uc*isqrt
    g_g_fu3(1,2) = 0.5d0*g_g_isqrt
    g_fu3(1,2) = 0.5d0*g_isqrt
    g_fu6(1,2) = 0.5d0*g_isqru
    fu3(1,2) = 0.5d0*isqrt
    g_g_fu3(2,1) = &
      & 0.5d0*g_g_isqrt*(g*hc-uc*uc)-g_g_uc*uc*isqrt+0.5d0*g_hd*g_isqrt*g+g_isqru*(0.5d0*g_hc*g-g_uc*uc)-g_ud*(g_isqrt*&
&uc+g_uc*isqrt)
    g_fu3(2,1) = 0.5d0*g_hc*g*isqrt+0.5d0*g_isqrt*(g*hc-uc*uc)-g_uc*uc*isqrt
    g_fu6(2,1) = 0.5d0*g_hd*g*isqrt+0.5d0*g_isqru*(g*hc-uc*uc)-g_ud*uc*isqrt
    fu3(2,1) = 0.5d0*(g*hc-uc*uc)*isqrt
    g_g_fu3(2,2) = 0.5d0*g_g_isqrt*uc+0.5d0*g_g_uc*isqrt+0.5d0*g_isqru*g_uc+0.5d0*g_ud*g_isqrt
    g_fu3(2,2) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    g_fu6(2,2) = 0.5d0*g_isqru*uc+0.5d0*g_ud*isqrt
    fu3(2,2) = 0.5d0*(1.d0+uc*isqrt)
    g_g_fu3(3,1) = &
      & (-(0.5d0*g_g_isqrt*uc*vc))-0.5d0*g_g_uc*isqrt*vc+0.5d0*g_g_vc*(1.d0-uc*isqrt)-g_isqru*(0.5d0*g_uc*vc+0.5d0*g_vc*&
&uc)-g_ud*(0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt)-g_vd*(0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt)
    g_fu3(3,1) = (-(0.5d0*g_isqrt*uc*vc))-0.5d0*g_uc*isqrt*vc+0.5d0*g_vc*(1.d0-uc*isqrt)
    g_fu6(3,1) = (-(0.5d0*g_isqru*uc*vc))-0.5d0*g_ud*isqrt*vc+0.5d0*g_vd*(1.d0-uc*isqrt)
    fu3(3,1) = 0.5d0*(1.d0-uc*isqrt)*vc
    g_g_fu3(3,2) = 0.5d0*g_g_isqrt*vc+0.5d0*g_g_vc*isqrt+0.5d0*g_isqru*g_vc+0.5d0*g_vd*g_isqrt
    g_fu3(3,2) = 0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt
    g_fu6(3,2) = 0.5d0*g_isqru*vc+0.5d0*g_vd*isqrt
    fu3(3,2) = 0.5d0*vc*isqrt
    g_dudxminus = g_u(:,i-1,j)*((-2)/(6.d0*dx))-g_u(:,i+2,j)/(6.d0*dx)+g_u(:,i+1,j)*(6/(6.d0*dx))+g_u(:,i,j)*((-3)/(6.d0*dx))
    g_dudxminut = g_v(:,i-1,j)*((-2)/(6.d0*dx))-g_v(:,i+2,j)/(6.d0*dx)+g_v(:,i+1,j)*(6/(6.d0*dx))+g_v(:,i,j)*((-3)/(6.d0*dx))
    dudxminus = ((-u(:,i+2,j))+6.d0*u(:,i+1,j)-3.d0*u(:,i,j)-2.d0*u(:,i-1,j))/(6.d0*dx)
    g_dudxplus = g_u(:,i-2,j)/(6.d0*dx)+g_u(:,i-1,j)*((-6)/(6.d0*dx))+g_u(:,i+1,j)*(2/(6.d0*dx))+g_u(:,i,j)*(3/(6.d0*dx))
    g_dudxplut = g_v(:,i-2,j)/(6.d0*dx)+g_v(:,i-1,j)*((-6)/(6.d0*dx))+g_v(:,i+1,j)*(2/(6.d0*dx))+g_v(:,i,j)*(3/(6.d0*dx))
    dudxplus = (2.d0*u(:,i+1,j)+3.d0*u(:,i,j)-6.d0*u(:,i-1,j)+u(:,i-2,j))/(6.d0*dx)
    if (eigf(1) .ge. 0.d0) then
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudxplut(j1)*g_fu1(i1,j1)+g_fu4(i1,j1)*g_dudxplus(j1)+g_g_fu1(i1,j1)*dudxplus(j1)+g_g_tmpmv(i1)
          g_tmpmv(i1) = g_dudxplus(j1)*fu1(i1,j1)+g_fu1(i1,j1)*dudxplus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudxplut(j1)*fu1(i1,j1)+g_fu4(i1,j1)*dudxplus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+fu1(i1,j1)*dudxplus(j1)
        end do
      end do
      g_g_tmpf = g_eigh(1)*g_tmpmv+g_g_eigf(1)*tmpmv+g_g_tmpmv*eigf(1)+g_tmpmw*g_eigf(1)
    else
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudxminut(j1)*g_fu1(i1,j1)+g_fu4(i1,j1)*g_dudxminus(j1)+g_g_fu1(i1,j1)*dudxminus(j1)+g_g_tmpmv(i1)
          g_tmpmv(i1) = g_dudxminus(j1)*fu1(i1,j1)+g_fu1(i1,j1)*dudxminus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudxminut(j1)*fu1(i1,j1)+g_fu4(i1,j1)*dudxminus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+fu1(i1,j1)*dudxminus(j1)
        end do
      end do
      g_g_tmpf = g_eigh(1)*g_tmpmv+g_g_eigf(1)*tmpmv+g_g_tmpmv*eigf(1)+g_tmpmw*g_eigf(1)
    endif
    if (eigf(2) .ge. 0.d0) then
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudxplut(j1)*g_fu2(i1,j1)+g_fu5(i1,j1)*g_dudxplus(j1)+g_g_fu2(i1,j1)*dudxplus(j1)+g_g_tmpmv(i1)
          g_tmpmv(i1) = g_dudxplus(j1)*fu2(i1,j1)+g_fu2(i1,j1)*dudxplus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudxplut(j1)*fu2(i1,j1)+g_fu5(i1,j1)*dudxplus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+fu2(i1,j1)*dudxplus(j1)
        end do
      end do
      g_g_tmpf = g_eigh(2)*g_tmpmv+g_g_eigf(2)*tmpmv+g_g_tmpf+g_g_tmpmv*eigf(2)+g_tmpmw*g_eigf(2)
    else
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudxminut(j1)*g_fu2(i1,j1)+g_fu5(i1,j1)*g_dudxminus(j1)+g_g_fu2(i1,j1)*dudxminus(j1)+g_g_tmpmv(i1)
          g_tmpmv(i1) = g_dudxminus(j1)*fu2(i1,j1)+g_fu2(i1,j1)*dudxminus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudxminut(j1)*fu2(i1,j1)+g_fu5(i1,j1)*dudxminus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+fu2(i1,j1)*dudxminus(j1)
        end do
      end do
      g_g_tmpf = g_eigh(2)*g_tmpmv+g_g_eigf(2)*tmpmv+g_g_tmpf+g_g_tmpmv*eigf(2)+g_tmpmw*g_eigf(2)
    endif
    if (eigf(3) .ge. 0.d0) then
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudxplut(j1)*g_fu3(i1,j1)+g_fu6(i1,j1)*g_dudxplus(j1)+g_g_fu3(i1,j1)*dudxplus(j1)+g_g_tmpmv(i1)
          g_tmpmv(i1) = g_dudxplus(j1)*fu3(i1,j1)+g_fu3(i1,j1)*dudxplus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudxplut(j1)*fu3(i1,j1)+g_fu6(i1,j1)*dudxplus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+fu3(i1,j1)*dudxplus(j1)
        end do
      end do
      g_g_tmpf = g_eigh(3)*g_tmpmv+g_g_eigf(3)*tmpmv+g_g_tmpf+g_g_tmpmv*eigf(3)+g_tmpmw*g_eigf(3)
    else
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudxminut(j1)*g_fu3(i1,j1)+g_fu6(i1,j1)*g_dudxminus(j1)+g_g_fu3(i1,j1)*dudxminus(j1)+g_g_tmpmv(i1)
          g_tmpmv(i1) = g_dudxminus(j1)*fu3(i1,j1)+g_fu3(i1,j1)*dudxminus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudxminut(j1)*fu3(i1,j1)+g_fu6(i1,j1)*dudxminus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+fu3(i1,j1)*dudxminus(j1)
        end do
      end do
      g_g_tmpf = g_eigh(3)*g_tmpmv+g_g_eigf(3)*tmpmv+g_g_tmpf+g_g_tmpmv*eigf(3)+g_tmpmw*g_eigf(3)
    endif
    g_g_eigg(1) = &
      & (-(g_v(3,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))+g_v(1,i,j)*((-(g_u(3,i,j)/(u(1,i,j)*u(1,i,j))))+g_u(1,i,j)*(2*&
&u(3,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))))
    g_eigg(1) = g_u(3,i,j)/u(1,i,j)-g_u(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))
    g_eigi(1) = g_v(3,i,j)/u(1,i,j)-g_v(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j)))
    eigg(1) = u(3,i,j)/u(1,i,j)
    g_g_eigg(2) = &
      & (-(g_v(3,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))-g_v(1,i,j)*(g_u(3,i,j)/(u(1,i,j)*u(1,i,j))-g_u(1,i,j)*(2*u(3,i,&
&j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))+2*1./(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*sqrt(gacc*u(1,i,&
&j)))*gacc))
    g_eigg(2) = g_u(3,i,j)/u(1,i,j)-g_u(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    g_eigi(2) = g_v(3,i,j)/u(1,i,j)-g_v(1,i,j)*(u(3,i,j)/(u(1,i,j)*u(1,i,j))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigg(2) = u(3,i,j)/u(1,i,j)-sqrt(gacc*u(1,i,j))
    g_g_eigg(3) = &
      & (-(g_v(3,i,j)*(g_u(1,i,j)/(u(1,i,j)*u(1,i,j)))))+g_v(1,i,j)*((-(g_u(3,i,j)/(u(1,i,j)*u(1,i,j))))+g_u(1,i,j)*(2*&
&u(3,i,j)*u(1,i,j)/(u(1,i,j)*u(1,i,j)*u(1,i,j)*u(1,i,j))-2*1./(2.*sqrt(gacc*u(1,i,j)))*gacc/(4*sqrt(gacc*u(1,i,j))*&
&sqrt(gacc*u(1,i,j)))*gacc))
    g_eigg(3) = g_u(3,i,j)/u(1,i,j)+g_u(1,i,j)*((-(u(3,i,j)/(u(1,i,j)*u(1,i,j))))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    g_eigi(3) = g_v(3,i,j)/u(1,i,j)+g_v(1,i,j)*((-(u(3,i,j)/(u(1,i,j)*u(1,i,j))))+1./(2.*sqrt(gacc*u(1,i,j)))*gacc)
    eigg(3) = u(3,i,j)/u(1,i,j)+sqrt(gacc*u(1,i,j))
    g_g_gu1(:,:) = 0.d0
    g_gu1(:,:) = 0.d0
    g_gu4(:,:) = 0.d0
    gu1(:,:) = 0.d0
    g_g_gu1(2,1) = -g_g_uc
    g_gu1(2,1) = -g_uc
    g_gu4(2,1) = -g_ud
    gu1(2,1) = -uc
    g_g_gu1(2,2) = 0.d0
    g_gu1(2,2) = 0.d0
    g_gu4(2,2) = 0.d0
    gu1(2,2) = 1.d0
    g_g_isqrt = &
      & g_hd*g_hc*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*&
&sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)))
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    g_isqru = -(g_hd*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_g_gu2(:,:) = 0.d0
    g_gu2(:,:) = 0.d0
    g_gu5(:,:) = 0.d0
    gu2(:,:) = 0.d0
    g_g_gu2(1,1) = 0.5d0*g_g_isqrt*vc+0.5d0*g_g_vc*isqrt+0.5d0*g_isqru*g_vc+0.5d0*g_vd*g_isqrt
    g_gu2(1,1) = 0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt
    g_gu5(1,1) = 0.5d0*g_isqru*vc+0.5d0*g_vd*isqrt
    gu2(1,1) = 0.5d0*(1.d0+vc*isqrt)
    g_g_gu2(1,3) = (-0.5d0)*g_g_isqrt
    g_gu2(1,3) = (-0.5d0)*g_isqrt
    g_gu5(1,3) = (-0.5d0)*g_isqru
    gu2(1,3) = -(0.5d0*isqrt)
    g_g_gu2(2,1) = &
      & 0.5d0*g_g_isqrt*uc*vc+0.5d0*g_g_uc*(1.d0+vc*isqrt)+0.5d0*g_g_vc*uc*isqrt+g_isqru*(0.5d0*g_uc*vc+0.5d0*g_vc*uc)+&
&g_ud*(0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt)+g_vd*(0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt)
    g_gu2(2,1) = 0.5d0*g_isqrt*uc*vc+0.5d0*g_uc*(1.d0+vc*isqrt)+0.5d0*g_vc*uc*isqrt
    g_gu5(2,1) = 0.5d0*g_isqru*uc*vc+0.5d0*g_ud*(1.d0+vc*isqrt)+0.5d0*g_vd*uc*isqrt
    gu2(2,1) = 0.5d0*uc*(1.d0+vc*isqrt)
    g_g_gu2(2,3) = (-(0.5d0*g_g_isqrt*uc))-0.5d0*g_g_uc*isqrt-0.5d0*g_isqru*g_uc-0.5d0*g_ud*g_isqrt
    g_gu2(2,3) = (-(0.5d0*g_isqrt*uc))-0.5d0*g_uc*isqrt
    g_gu5(2,3) = (-(0.5d0*g_isqru*uc))-0.5d0*g_ud*isqrt
    gu2(2,3) = -(0.5d0*uc*isqrt)
    g_g_gu2(3,1) = 0.5d0*g_g_isqrt*(vc*vc-g*hc)+g_g_vc*vc*isqrt-0.5d0*g_hd*g_isqrt*g+g_isqru*((-(0.5d0*g_hc*g))+g_vc*vc)+g_vd*&
&(g_isqrt*vc+g_vc*isqrt)
    g_gu2(3,1) = (-(0.5d0*g_hc*g*isqrt))+0.5d0*g_isqrt*(vc*vc-g*hc)+g_vc*vc*isqrt
    g_gu5(3,1) = (-(0.5d0*g_hd*g*isqrt))+0.5d0*g_isqru*(vc*vc-g*hc)+g_vd*vc*isqrt
    gu2(3,1) = 0.5d0*(vc*vc-g*hc)*isqrt
    g_g_gu2(3,3) = (-(0.5d0*g_g_isqrt*vc))-0.5d0*g_g_vc*isqrt-0.5d0*g_isqru*g_vc-0.5d0*g_vd*g_isqrt
    g_gu2(3,3) = (-(0.5d0*g_isqrt*vc))-0.5d0*g_vc*isqrt
    g_gu5(3,3) = (-(0.5d0*g_isqru*vc))-0.5d0*g_vd*isqrt
    gu2(3,3) = 0.5d0-0.5d0*vc*isqrt
    g_g_isqrt = &
      & g_hd*g_hc*(2*1./(2.*sqrt(g*hc))*g/(4*sqrt(g*hc)*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))+2*1/(2.*sqrt(g*hc))*g*1./(2.*&
&sqrt(g*hc))*g*sqrt(g*hc)/(sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)*sqrt(g*hc)))
    g_isqrt = -(g_hc*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    g_isqru = -(g_hd*(1.d0*1./(2.*sqrt(g*hc))*g/(sqrt(g*hc)*sqrt(g*hc))))
    isqrt = 1.d0/sqrt(g*hc)
    g_g_gu3(:,:) = 0.d0
    g_gu3(:,:) = 0.d0
    g_gu6(:,:) = 0.d0
    gu3(:,:) = 0.d0
    g_g_gu3(1,1) = (-(0.5d0*g_g_isqrt*vc))-0.5d0*g_g_vc*isqrt-0.5d0*g_isqru*g_vc-0.5d0*g_vd*g_isqrt
    g_gu3(1,1) = (-(0.5d0*g_isqrt*vc))-0.5d0*g_vc*isqrt
    g_gu6(1,1) = (-(0.5d0*g_isqru*vc))-0.5d0*g_vd*isqrt
    gu3(1,1) = 0.5d0*(1.d0-vc*isqrt)
    g_g_gu3(1,3) = 0.5d0*g_g_isqrt
    g_gu3(1,3) = 0.5d0*g_isqrt
    g_gu6(1,3) = 0.5d0*g_isqru
    gu3(1,3) = 0.5d0*isqrt
    g_g_gu3(2,1) = &
      & (-(0.5d0*g_g_isqrt*uc*vc))+0.5d0*g_g_uc*(1.d0-vc*isqrt)-0.5d0*g_g_vc*uc*isqrt-g_isqru*(0.5d0*g_uc*vc+0.5d0*g_vc*&
&uc)-g_ud*(0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt)-g_vd*(0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt)
    g_gu3(2,1) = (-(0.5d0*g_isqrt*uc*vc))+0.5d0*g_uc*(1.d0-vc*isqrt)-0.5d0*g_vc*uc*isqrt
    g_gu6(2,1) = (-(0.5d0*g_isqru*uc*vc))+0.5d0*g_ud*(1.d0-vc*isqrt)-0.5d0*g_vd*uc*isqrt
    gu3(2,1) = 0.5d0*uc*(1.d0-vc*isqrt)
    g_g_gu3(2,3) = 0.5d0*g_g_isqrt*uc+0.5d0*g_g_uc*isqrt+0.5d0*g_isqru*g_uc+0.5d0*g_ud*g_isqrt
    g_gu3(2,3) = 0.5d0*g_isqrt*uc+0.5d0*g_uc*isqrt
    g_gu6(2,3) = 0.5d0*g_isqru*uc+0.5d0*g_ud*isqrt
    gu3(2,3) = 0.5d0*uc*isqrt
    g_g_gu3(3,1) = 0.5d0*g_g_isqrt*((-(vc*vc))+g*hc)-g_g_vc*vc*isqrt+0.5d0*g_hd*g_isqrt*g+g_isqru*(0.5d0*g_hc*g-g_vc*vc)-g_vd*&
&(g_isqrt*vc+g_vc*isqrt)
    g_gu3(3,1) = 0.5d0*g_hc*g*isqrt+0.5d0*g_isqrt*((-(vc*vc))+g*hc)-g_vc*vc*isqrt
    g_gu6(3,1) = 0.5d0*g_hd*g*isqrt+0.5d0*g_isqru*((-(vc*vc))+g*hc)-g_vd*vc*isqrt
    gu3(3,1) = 0.5d0*((-(vc*vc))+g*hc)*isqrt
    g_g_gu3(3,3) = 0.5d0*g_g_isqrt*vc+0.5d0*g_g_vc*isqrt+0.5d0*g_isqru*g_vc+0.5d0*g_vd*g_isqrt
    g_gu3(3,3) = 0.5d0*g_isqrt*vc+0.5d0*g_vc*isqrt
    g_gu6(3,3) = 0.5d0*g_isqru*vc+0.5d0*g_vd*isqrt
    gu3(3,3) = 0.5d0+0.5d0*vc*isqrt
    g_dudyminus = g_u(:,i,j-1)*((-2)/(6.d0*dy))-g_u(:,i,j+2)/(6.d0*dy)+g_u(:,i,j+1)*(6/(6.d0*dy))+g_u(:,i,j)*((-3)/(6.d0*dy))
    g_dudyminut = g_v(:,i,j-1)*((-2)/(6.d0*dy))-g_v(:,i,j+2)/(6.d0*dy)+g_v(:,i,j+1)*(6/(6.d0*dy))+g_v(:,i,j)*((-3)/(6.d0*dy))
    dudyminus = ((-u(:,i,j+2))+6.d0*u(:,i,j+1)-3.d0*u(:,i,j)-2.d0*u(:,i,j-1))/(6.d0*dy)
    g_dudyplus = g_u(:,i,j-2)/(6.d0*dy)+g_u(:,i,j-1)*((-6)/(6.d0*dy))+g_u(:,i,j+1)*(2/(6.d0*dy))+g_u(:,i,j)*(3/(6.d0*dy))
    g_dudyplut = g_v(:,i,j-2)/(6.d0*dy)+g_v(:,i,j-1)*((-6)/(6.d0*dy))+g_v(:,i,j+1)*(2/(6.d0*dy))+g_v(:,i,j)*(3/(6.d0*dy))
    dudyplus = (2.d0*u(:,i,j+1)+3.d0*u(:,i,j)-6.d0*u(:,i,j-1)+u(:,i,j-2))/(6.d0*dy)
    if (eigg(1) .ge. 0.d0) then
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudyplut(j1)*g_gu1(i1,j1)+g_g_gu1(i1,j1)*dudyplus(j1)+g_g_tmpmv(i1)+g_gu4(i1,j1)*g_dudyplus(j1)
          g_tmpmv(i1) = g_dudyplus(j1)*gu1(i1,j1)+g_gu1(i1,j1)*dudyplus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudyplut(j1)*gu1(i1,j1)+g_gu4(i1,j1)*dudyplus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+gu1(i1,j1)*dudyplus(j1)
        end do
      end do
      g_g_tmpg = g_eigi(1)*g_tmpmv+g_g_eigg(1)*tmpmv+g_g_tmpmv*eigg(1)+g_tmpmw*g_eigg(1)
    else
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudyminut(j1)*g_gu1(i1,j1)+g_g_gu1(i1,j1)*dudyminus(j1)+g_g_tmpmv(i1)+g_gu4(i1,j1)*g_dudyminus(j1)
          g_tmpmv(i1) = g_dudyminus(j1)*gu1(i1,j1)+g_gu1(i1,j1)*dudyminus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudyminut(j1)*gu1(i1,j1)+g_gu4(i1,j1)*dudyminus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+gu1(i1,j1)*dudyminus(j1)
        end do
      end do
      g_g_tmpg = g_eigi(1)*g_tmpmv+g_g_eigg(1)*tmpmv+g_g_tmpmv*eigg(1)+g_tmpmw*g_eigg(1)
    endif
    if (eigg(2) .ge. 0.d0) then
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudyplut(j1)*g_gu2(i1,j1)+g_g_gu2(i1,j1)*dudyplus(j1)+g_g_tmpmv(i1)+g_gu5(i1,j1)*g_dudyplus(j1)
          g_tmpmv(i1) = g_dudyplus(j1)*gu2(i1,j1)+g_gu2(i1,j1)*dudyplus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudyplut(j1)*gu2(i1,j1)+g_gu5(i1,j1)*dudyplus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+gu2(i1,j1)*dudyplus(j1)
        end do
      end do
      g_g_tmpg = g_eigi(2)*g_tmpmv+g_g_eigg(2)*tmpmv+g_g_tmpg+g_g_tmpmv*eigg(2)+g_tmpmw*g_eigg(2)
    else
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudyminut(j1)*g_gu2(i1,j1)+g_g_gu2(i1,j1)*dudyminus(j1)+g_g_tmpmv(i1)+g_gu5(i1,j1)*g_dudyminus(j1)
          g_tmpmv(i1) = g_dudyminus(j1)*gu2(i1,j1)+g_gu2(i1,j1)*dudyminus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudyminut(j1)*gu2(i1,j1)+g_gu5(i1,j1)*dudyminus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+gu2(i1,j1)*dudyminus(j1)
        end do
      end do
      g_g_tmpg = g_eigi(2)*g_tmpmv+g_g_eigg(2)*tmpmv+g_g_tmpg+g_g_tmpmv*eigg(2)+g_tmpmw*g_eigg(2)
    endif
    if (eigg(3) .ge. 0.d0) then
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudyplut(j1)*g_gu3(i1,j1)+g_g_gu3(i1,j1)*dudyplus(j1)+g_g_tmpmv(i1)+g_gu6(i1,j1)*g_dudyplus(j1)
          g_tmpmv(i1) = g_dudyplus(j1)*gu3(i1,j1)+g_gu3(i1,j1)*dudyplus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudyplut(j1)*gu3(i1,j1)+g_gu6(i1,j1)*dudyplus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+gu3(i1,j1)*dudyplus(j1)
        end do
      end do
      g_g_tmpg = g_eigi(3)*g_tmpmv+g_g_eigg(3)*tmpmv+g_g_tmpg+g_g_tmpmv*eigg(3)+g_tmpmw*g_eigg(3)
    else
      g_g_tmpmv(:) = 0.d0
      g_tmpmv(:) = 0.d0
      g_tmpmw(:) = 0.d0
      tmpmv(:) = 0.d0
      do i1 = 1, 3
        do j1 = 1, 3
          g_g_tmpmv(i1) = g_dudyminut(j1)*g_gu3(i1,j1)+g_g_gu3(i1,j1)*dudyminus(j1)+g_g_tmpmv(i1)+g_gu6(i1,j1)*g_dudyminus(j1)
          g_tmpmv(i1) = g_dudyminus(j1)*gu3(i1,j1)+g_gu3(i1,j1)*dudyminus(j1)+g_tmpmv(i1)
          g_tmpmw(i1) = g_dudyminut(j1)*gu3(i1,j1)+g_gu6(i1,j1)*dudyminus(j1)+g_tmpmw(i1)
          tmpmv(i1) = tmpmv(i1)+gu3(i1,j1)*dudyminus(j1)
        end do
      end do
      g_g_tmpg = g_eigi(3)*g_tmpmv+g_g_eigg(3)*tmpmv+g_g_tmpg+g_g_tmpmv*eigg(3)+g_tmpmw*g_eigg(3)
    endif
    g_g_rhs = (-g_g_tmpf)-g_g_tmpg
    g_g_f(:,i,j) = g_g_rhs
  end do
end do

end subroutine g_g_compute_f

