!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine DRDP(nadj,n,nrp,t,y,rp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer :: nadj,n,nrp
double precision :: t,y(n),rp(nrp,nadj)
rp(:,:) = 0d0
end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine DRDY(nadj,n,nry,t,y,ry)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer :: nadj,n,nry
  double precision :: t,y(n),ry(nry,nadj)
  ry(:,:) = 0d0
  ry(1,1) = 1.0d0
  ry(2,2) = 2.0d0
  ry(3,3) = 3.0d0 
end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine JACP(n,np,t,y,fpjac)
!       fpjac = df/dp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
integer :: n,np
double precision :: t,y(n),fpjac(n,np)
  fpjac(1,1) = -y(1)
  fpjac(2,1) = y(1)
  fpjac(1,2) = y(2)*y(3)
  fpjac(2,2) = -y(2)*y(3)
  fpjac(2,3) = -y(2)*y(2)
  fpjac(3,3) = y(2)*y(2)
end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine jac(nvar,t,y,fjac)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  implicit none
  integer :: nvar
  double precision :: c(3)
  double precision ::t
  double precision ::y(nvar),fjac(nvar,nvar)
  c(1) = 0.04e0
  c(2) = 1.0e4
  c(3) = 3.0e7

  fjac(1,1) = -c(1)
  fjac(1,2) = c(2)*y(3)
  fjac(1,3) = c(2)*y(2)

  fjac(2,1) = c(1)
  fjac(2,2) = -c(2)*y(3)-2*c(3)*y(2)
  fjac(2,3) = -c(2)*y(2)

  fjac(3,2) = 2*c(3)*y(2)
end subroutine jac

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine fun(nvar,t, y, p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  implicit none
  integer :: nvar
  double precision,parameter :: c(3) = (/0.04e0, 1.0e4, 3.0e7/)
  double precision :: t
  double precision :: y(nvar), p(nvar)

  p(1) = -c(1)*y(1) +c(2)*y(2)*y(3)
  p(3) = c(3)*y(2)*y(2)
  p(2) = -p(1) - p(3)
      
end subroutine fun

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine qfun(n,nr,t, y, r )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  implicit none
  integer :: n,nr
  double precision :: t
  double precision :: y(n), r(nr)

  r(1) = y(1)
  r(2) = y(2)
  r(3) = y(3)

end subroutine qfun

!~~~> Initialize adjoint variables
subroutine adjinit(n,np,nadj,t,y,lambda,mu)
  integer :: n,np,nadj,k
  double precision :: t,y(n),lambda(n,nadj)
  double precision,optional :: mu(np,nadj)
!~~~> if number of parameters is not zero, extra adjoint varaible mu should be
!defined
  if(NP>0 .and. .not. present(mu)) stop 'undefined argument mu'
!~~~>  the adjoint values at the final time
  lambda(1:n,1:nadj) = 0.0d0
  do k=1,nadj
    lambda(k,k) = y(k)
  end do
  mu(1:np,1:nadj) = 0.d0
  mu(1,1) = -y(1)*y(1)
  mu(1,2) = y(1)*y(2)*y(3)
  mu(2,1) = y(1)*y(2)
  mu(2,2) = -y(2)*y(2)*y(3)
  mu(2,3) = -y(2)*y(2)*y(2)
  mu(3,3) = y(2)*y(2)*y(3)
end subroutine

  program roberts_ros_adj_dr

    use sdirk_adj_f90_integrator
    implicit none
    external fun,qfun,jac, adjinit,DRDP,DRDY,JACP
    integer, parameter :: nvar = 3, nnz=0, nadj=3, np=3
    integer ::i,k
    double precision :: atol(nvar),rtol(nvar), var(nvar), q(nadj)
    double precision :: atol_adj(nvar,nadj),rtol_adj(nvar,nadj),y_adj(nvar,nadj),yp_adj(np,nadj)
    double precision:: tstart, tend, t1, t2, tcpu
    double precision:: rstate(20),rcntrl(20)
    integer ::istate(20),icntrl(20)


    ! Print problem description
    print *,'Adjoint Sensitivity Example for Chemical Kinetics'
    print *, '-------------------------------------------------'
    print *, 'ODE: dy1/dt = -p1*y1 + p2*y2*y3 '
    print *, '     dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2'
    print *, '     dy3/dt =  p3*(y2)^2 '
    print *, ' Find dG/dp for '
    print *, '     G = int_t0^tB0 g(t,p,y) dt '
    print *, '     g(t,p,y) = y(t,p) '


    do i=1,nvar
      rtol(i) = 1e-6
    end do

    atol(1) = 1e-8
    atol(2) = 1e-14
    atol(3) = 1e-6

    do k=1,nadj
      do i=1,nvar
        rtol_adj(i,k) = 1e-1
        atol_adj(i,k) = 1e-2
      end do
    end do

    ! initialize y
    var(1) = 1.0d0
    var(2) = 0.0d0
    var(3) = 0.0d0
    q(:) = 0.d0
    tstart = 0.0d0
    tend = 4.0e7
    istate(:)=0
    icntrl(:)=0
    rcntrl(:)=0
! change initial step size
    rcntrl(3)=1e-3
      
    call cpu_time(t1)
    call integrate_adj(tin=tstart, tout=tend, nvar=nvar, np=np, nnzero=nnz,&
        lambda=y_adj, mu= yp_adj, y=var, rtol=rtol, atol=atol, nadj=nadj,&
        atol_adj=atol_adj, rtol_adj=rtol_adj, fun=fun, jac=jac,adjinit=adjinit,&
        drdp =DRDP, drdy=DRDY, jacp=JACP,qfun=qfun,q=q,&
        rstatus_u=rstate,rcntrl_u=rcntrl, istatus_u=istate, icntrl_u=icntrl)
    call cpu_time(t2)
    tcpu=t2-t1

    write(6,2010) tstart,tend,istate(3:5),tcpu
2010    FORMAT (2x,'Time interval [',f12.2,'-',f12.2, ']',/,&
          2x,'Nsteps=',I5,' NAsteps=',I5,' NRsteps=',I5,/,&
               2x, 'CPUtime=', f14.3,/,2x, //)

    print *, '-------------------------------------'
    print *, 'dG/dp', yp_adj(:,1)
    print *, 'lambda',y_adj(:,1)
    print *, 'VAR',var(:)
    print *, 'Q',q(:)
end program roberts_ros_adj_dr
