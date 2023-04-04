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


subroutine compute_F(U,F,feigvalues,geigvalues)

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

  integer           :: i, j, i1, j1
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


subroutine compute_JacF(U,feigvalues,geigvalues)

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
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.5 + (0.5*C2)/(C1*SqrtC1g))/(2.*dx)+(((-0.25*C2*g)/(C1*(C1*g)**1.5)-&
      & (0.5*C2)/(C1**2*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) +&
      & (0.0416667*g*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*(C1*g)**1.5)) + (-(C2/C1**2) - g/(2.*SqrtC1g))*(((0.5 + & 
      & (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) - (0.0833333*(3*C2 - 6*N2 + NN2 + & 
      & 2*S2))/(dx*SqrtC1g))
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
      & (((0.5 + (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) - &
      & (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*SqrtC1g))/C1
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
      & (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) + (0.0833333*(3*C2 - 6*N2 + NN2 + &
      & 2*S2))/(dx*SqrtC1g))
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
      & - (0.5*C2)/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*dx) + &
      & (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(dx*SqrtC1g))/C1
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
      & g/(2.*SqrtC1g))*(((0.5 - (0.5*C2)/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*dx) + (0.0833333*(-3*C2 -&
      & 2*N2 + 6*S2 - SS2))/(dx*SqrtC1g))
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
      &  SS2))/(dx*SqrtC1g))/C1
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
      & (0.5*C3)/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*dy) - (0.0833333*(3*C3 + 2*E3 - 6*W3 + &
      &  WW3))/(dy*SqrtC1g))
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
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*(-(0.5 + (0.5*C3)/(C1*SqrtC1g))/(2.*dy)+(((-0.25*C3*g)/(C1*(C1*g)**1.5) -&
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
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*(0.25/(dy*SqrtC1g) + (0.0833333*(-3*C1 + 6*E1 - EE1 - &
      & 2*W1))/(C1*dy*SqrtC1g)) +&
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
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*(0.25/(dy*SqrtC1g) - (0.0833333*(3*C1 + 2*E1 - 6*W1 + &
      & WW1))/(C1*dy*SqrtC1g)) +&
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
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*(-(0.5 - (0.5*C3)/(C1*SqrtC1g))/(2.*dy)+(((0.25*C3*g)/(C1*(C1*g)**1.5) +&
      & (0.5*C3)/(C1**2*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) - (0.0416667*g*(-3*C3 + 6*E3 - EE3 - &
      & 2*W3))/(dy*(C1*g)**1.5)) +&
      & (-(C3/C1**2) + g/(2.*SqrtC1g))*(((0.5 - (0.5*C3)/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) + &
      & (0.0833333*(-3*C3 + 6*E3 - EE3 - 2*W3))/(dy*SqrtC1g))
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
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*(-0.25/(dy*SqrtC1g) - (0.0833333*(-3*C1 + 6*E1 - EE1 - & 
      &2*W1))/(C1*dy*SqrtC1g)) +&
      & (((0.5 - (0.5*C3)/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*dy) + (0.0833333*(-3*C3 + 6*E3 - EE3 -&
      &  2*W3))/(dy*SqrtC1g))/C1
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
      & (-(C2/C1**2) + g/(2.*SqrtC1g))*((0.0833333*(-(C2**2/C1**2) + C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g)+&
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
      & (0.166667*C2*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1**2*dx*SqrtC1g) +  (0.0833333*(3*C2 - 6*N2 + NN2 + & 
      & 2*S2))/(C1*dx*SqrtC1g)) +&
      & ((0.0833333*(-(C2**2/C1**2) + C1*g)*(3*C1 - 6*N1 + NN1 + 2*S1))/(dx*SqrtC1g) + ((0.5 + &
      & (0.5*C2)/(C1*SqrtC1g))*(3*C2 - 6*N2 + NN2 + 2*S2))/(6.*dx))/C1
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
      & (0.0416667*g*(-(C2**2/C1**2) + C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*(C1*g)**1.5) + &
      & (((-0.25*C2*g)/(C1*(C1*g)**1.5) -&
      & (0.5*C2)/(C1**2*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 - SS2))/(6.*dx)) + (-(C2/C1**2) + &
      & g/(2.*SqrtC1g))*((0.0833333*(-(C2**2/C1**2) +&
      & C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g) + ((0.5 + (0.5*C2)/(C1*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 -&
      & SS2))/(6.*dx))
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
      & (0.166667*C2*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1**2*dx*SqrtC1g) + (0.0833333*(-3*C2 - 2*N2 + 6*S2 -&
      & SS2))/(C1*dx*SqrtC1g)) +&
      & ((0.0833333*(-(C2**2/C1**2) + C1*g)*(-3*C1 - 2*N1 + 6*S1 - SS1))/(dx*SqrtC1g)+((0.5 + &
      & (0.5*C2)/(C1*SqrtC1g))*(-3*C2 - 2*N2 + 6*S2 - SS2))/(6.*dx))/C1
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
      JG(p,ind)=-C2/(2.*C1*dy) + (C2*(3*C1 + 2*E1 - 6*W1 + WW1))/(6.*C1**2*dy) - (C3*(3*C2 + 2*E2 - 6*W2 + &
      & WW2))/(6.*C1**2*dy)
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
      JG(p,ind)=C2/(2.*C1*dy) + (C2*(-3*C1 + 6*E1 - EE1 - 2*W1))/(6.*C1**2*dy) - (C3*(-3*C2 + 6*E2 - EE2 -&
      & 2*W2))/(6.*C1**2*dy)
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
      & (0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1**2*dy) + (0.0416667*C2*g*(3*C3 + 2*E3 -&
      & 6*W3 + WW3))/(C1*dy*(C1*g)**1.5) + &
      & (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1**2*dy*SqrtC1g)) + (-(C3/C1**2) - &
      & g/(2.*SqrtC1g))*((0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) &
      & - (0.0833333*C2*(3*C3 + 2*E3 - 6*W3 + WW3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,-1,0,1)  !E1
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(0.166667*C2*(1 + C3/(C1*SqrtC1g))*(C3/C1 - SqrtC1g))/(C1*dy)
    end if  
	call mapgr(ind,p,0,0,0)   !C2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((0.0833333*(1 + C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) -&
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
      JG(p,ind)=JG(p,ind)+(C3/C1 - SqrtC1g)*((-0.25*C2)/(C1*dy*SqrtC1g) + (0.0833333*C2*(3*C1 + 2*E1 - 6*W1 +&
      & WW1))/(C1**2*dy*SqrtC1g)) +&
      & ((0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) - (0.0833333*C2*(3*C3 + 2*E3 &
      & - 6*W3 + WW3))/(C1*dy*SqrtC1g))/C1
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
      & C3/(C1**2*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) - (0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(-3*C1 &
      & + 6*E1 - EE1 - 2*W1))/(C1**2*dy) +&
      & (0.0416667*C2*g*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*(C1*g)**1.5) + (0.0833333*C2*(-3*C3 + 6*E3 - EE3 -& 
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
      & ((0.0833333*C2*(1 + C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) - (0.0833333*C2*(-3*C3 + 6*E3 & 
      & - EE3 - 2*W3))/(C1*dy*SqrtC1g))/C1
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
      & C3/(C1**2*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) - (0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(3*C1 &
      & + 2*E1 - 6*W1 + WW1))/(C1**2*dy) -&
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
      & ((0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(3*C1 + 2*E1 - 6*W1 + WW1))/(C1*dy) + (0.0833333*C2*(3*C3 + 2*E3 &
      & - 6*W3 + WW3))/(C1*dy*SqrtC1g))/C1
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
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((-0.25*C2*(1 - C3/(C1*SqrtC1g)))/(C1*dy) + &
      & (0.0833333*C2*((C3*g)/(2.*C1*(C1*g)**1.5) +&
      & C3/(C1**2*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) - (0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(-3*C1 &
      & + 6*E1 - EE1 - 2*W1))/(C1**2*dy) -&
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
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((0.0833333*(1 - C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))&
      & /(C1*dy) + (0.0833333*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*SqrtC1g))
    end if  
    call mapgr(ind,p,1,0,-1)  !W3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(-0.166667*C2*(C3/C1 + SqrtC1g))/(C1*dy*SqrtC1g)
    end if  
    call mapgr(ind,p,1,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JG(p,ind)=JG(p,ind)+(C3/C1 + SqrtC1g)*((-0.25*C2)/(C1*dy*SqrtC1g) - (0.0833333*C2*(-3*C1 + 6*E1 &
      & - EE1 - 2*W1))/(C1**2*dy*SqrtC1g)) +&
      & ((0.0833333*C2*(1 - C3/(C1*SqrtC1g))*(-3*C1 + 6*E1 - EE1 - 2*W1))/(C1*dy) &
      & + (0.0833333*C2*(-3*C3 + 6*E3 - EE3 - 2*W3))/(C1*dy*SqrtC1g))/C1
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
      JF(p,ind)=JF(p,ind)+(C2*(-C3/(2.*C1*dx) + (C3*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*C1**2*dx)))/C1 &
      & - (C2*(-(C3*(3*C1 - 6*N1 + NN1 + 2*S1))/(6.*C1*dx) + (3*C3 - 6*N3 + NN3 + 2*S3)/(6.*dx)))/C1**2
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
      JF(p,ind)=JF(p,ind)+(C2*(C3/(2.*C1*dx) + (C3*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*C1**2*dx)))/C1 &
      & - (C2*(-(C3*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*C1*dx) + (-3*C3 - 2*N3 + 6*S3 - SS3)/(6.*dx)))/C1**2
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
      JF(p,ind)=JF(p,ind)+(-(C3*(-3*C1 - 2*N1 + 6*S1 - SS1))/(6.*C1*dx) + (-3*C3 - 2*N3 + 6*S3 - SS3)&
      & /(6.*dx))/C1
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
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.0833333*(1 + C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 &
      & + 2*S1))/(C1*dx) - (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))
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
      & ((0.0833333*C3*(1 + C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))/(C1*dx) &
      & - (0.0833333*C3*(-3*C2 - 2*N2 + 6*S2 - SS2))/(C1*dx*SqrtC1g))/C1
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
      JF(p,ind)=JF(p,ind)+(C2/C1 - SqrtC1g)*((0.0833333*(1 + C2/(C1*SqrtC1g))*(-3*C1 - 2*N1 + 6*S1 - SS1))&
      & /(C1*dx) -&
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
      & (0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1**2*dx) - (0.0416667*C3*g*(3*C2 &
      & - 6*N2 + NN2 + 2*S2))/(C1*dx*(C1*g)**1.5) -&
      & (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1**2*dx*SqrtC1g)) + (-(C2/C1**2) +&
      & g/(2.*SqrtC1g))*((0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) &
      & + (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))
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
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.25*C3)/(C1*dx*SqrtC1g) - (0.0833333*C3*(3*C1 - 6*N1 + NN1 +&
      & 2*S1))/(C1**2*dx*SqrtC1g)) +&
      & ((0.0833333*C3*(1 - C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 + 2*S1))/(C1*dx) &
      & + (0.0833333*C3*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))/C1
    end if  
    call mapgr(ind,p,-1,1,0)  !S2
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(0.166667*C3*(C2/C1 + SqrtC1g))/(C1*dx*SqrtC1g)
    end if  
    call mapgr(ind,p,0,0,0)   !C3
    if ((ind .GE. 1) .AND. (ind .LE. ndim)) then
      JF(p,ind)=JF(p,ind)+(C2/C1 + SqrtC1g)*((0.0833333*(1 - C2/(C1*SqrtC1g))*(3*C1 - 6*N1 + NN1 & 
      & + 2*S1))/(C1*dx) + (0.0833333*(3*C2 - 6*N2 + NN2 + 2*S2))/(C1*dx*SqrtC1g))
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

