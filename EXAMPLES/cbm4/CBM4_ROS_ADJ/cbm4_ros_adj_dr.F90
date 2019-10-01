!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine hesstr_vec(n,t,y,u,v,htv)
! htv = (f_yy * v)^T *u = (d(f_y^T * u)/dy) *v
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  use cbm4_parameters, only: fix,rconst,hess,update_sun,update_rconst,&
                             hessian,cbm4_hesstr_vec
                             
  integer, intent(in) :: n
  double precision, intent(in) :: t,y(n),u(n),v(n)
  double precision, intent(inout)  :: htv(n)
  
  call update_sun(t)
  call update_rconst()
  call hessian(y,fix,rconst,hess)
  call cbm4_hesstr_vec(hess,u,v,htv)
  
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine hesstr_vec_f_py(n,np,t,y,u,k,htvg)
! htv = (f_py * k)^T *u = (d(f_p^T * u)/dy) *k
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  use cbm4_parameters, only: fix,rconst,update_sun,update_rconst,&
                             djac_drcoeff
  integer, intent(in) :: n,np
  double precision, intent(in) :: t,y(n),u(n),k(n)
  double precision, intent(inout) :: htvg(np)

  double precision :: djdr(n,np)
  integer :: i,j
  integer :: jcoeff(29) = (/4,11,18,21,24,25,36,37,41,44,49,50,52,54,55,&
                            59,64,65,66,67,68,70,73,75,76,77,78,79,81 /)

  call update_sun(t)
  call update_rconst()
  call djac_drcoeff(y,fix,k,np,jcoeff,djdr)
!  u and k 
  do i=1,np
    htvg(i) = 0.d0
    do j=1,n
      htvg(i) = htvg(i) + djdr(j,i)*u(j)
    end do
  end do

  end subroutine


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine jacp(n, np, t, y, fpjac)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  use cbm4_parameters, only:fix,update_sun,update_rconst,dFun_dRcoeff
  implicit none
  integer, intent(in) :: n,np
  double precision, intent(in) :: t,y(n)
  double precision, intent(inout) :: fpjac(n,np)
   
  integer :: jcoeff(29) = (/4,11,18,21,24,25,36,37,41,44,49,50,52,54,55,&
                            59,64,65,66,67,68,70,73,75,76,77,78,79,81 /) 
  integer :: i 
  fpjac(:,:) = 0.0d0

  call update_sun(t)
  call update_rconst()
  call dFun_dRcoeff( y, fix, np, jcoeff, fpjac)

  end subroutine jacp

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine jac(n, t, y, fjac)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  use cbm4_parameters, only:fix,rconst,update_sun,update_rconst,cbm4_jac
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) ::t,y(n)
  double precision, intent(inout) :: fjac(n,n)
  fjac(:,:) = 0.0d0
  call update_sun(t)
  call update_rconst()
  call cbm4_jac(y,fix,rconst,fjac)
  end subroutine jac
 
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine fun(n, t, y, p )
 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  use cbm4_parameters, only:fix,rconst,update_sun,update_rconst,cbm4_fun
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: t,y(n)
  double precision, intent(inout) :: p(n)
  
  call update_sun(t)
  call update_rconst()
  call cbm4_fun(y,fix,rconst,p)

  ! add emisions 
  p(31) = p(31) + 2.0d3
  p(26) = p(26) + 2.0d3
  end subroutine fun
  
!~~~> initialize adjoint variables
  subroutine adjinit(n,np,nadj,t,y,lambda,mu)
  implicit none
  integer, intent(in) :: n,np,nadj
  double precision, intent(in) :: t,y(n)
  double precision, intent(inout) :: lambda(n,nadj)
  double precision, intent(inout), optional :: mu(np,nadj)
  
  integer :: k

  if(np>0 .and. .not. present(mu)) stop 'underfined argument mu'
  lambda(1:n,1:nadj) = 0.0d0

  do k=1,nadj
    lambda(k,k) = 1.0d0
  end do
  if(np>0 .and. present(mu)) mu(1:np,1:nadj) = 0.0d0
  end subroutine

  subroutine output_statics(var,tstart,tend, tcpu, istate)
  use cbm4_parameters
  implicit none
  double precision :: tstart, tend, tcpu
  double precision :: var(nvar)
  integer :: istate(20)
  double precision :: tmp(nvar), delta(nvar)
  double precision :: l2, l2tmp, error, atol
  integer :: i
 
  write(6,2010) tstart,tend,istate(3:5),tcpu
  2010   format (2x,'time interval [',f10.2,'-',f10.2, ']    nsteps=',i5,&
              ' nasteps=',i5,' nrsteps=',i5,/,2x,&
                'cputime=', f14.3,'seconds',/,2x,//)
 ! store results into a file
       open(unit=11,file='cbm4_ros_sol.txt',status='replace')
       do i=1,nvar
         write(11,2011) var(i)
  2011   format(e24.16)
       end do
       close(11)
 
  end subroutine output_statics

  subroutine output_sens(nvar,np,nadj,lambda,mu)
  implicit none
  integer :: nvar, nadj, np
  double precision :: lambda(nvar,nadj), mu(np,nadj)
  integer :: i,j

  open(unit=2012,file='cbm4_output_sens.txt',status='replace')

  do i=1,nadj
    write(2012,'(A,I3)') 'Lambda: column',i
    do j=1,nvar
      write(2012,'(E24.16)') lambda(j,i)
    end do
  end do

  do i=1,nadj
    write(2012,'(A,I3)') 'Mu: column',i
    do j=1,np
      write(2012,'(E24.16)') mu(j,i)
    end do
  end do
  close(2012)
  end subroutine output_sens

  program cbm4_ros_adj_dr
    use cbm4_parameters
    use ros_adj_f90_integrator
    implicit none
    external fun,jac,adjinit,jacp,hesstr_vec,hesstr_vec_f_py
    integer :: i,j,state,ntests
    integer, parameter :: np=29, nadj =nvar
    double precision :: lambda(nvar,nadj),mu(np,nadj)
    double precision :: atol(nvar),rtol(nvar), atol_adj(nvar,nadj),rtol_adj(nvar,nadj),var(nvar)
    double precision :: tstart, tend, t1, t2, tcpu
    double precision :: rstate(20),rcntrl(20)
    double precision :: a=30d0, cfactor=2.55d10
    integer ::istate(20),icntrl(20)
    
    do i=1,nvar
      rtol(i) = 1.0d-4
      atol(i) = 1.0d-6
    end do
 
    do i = 1,nvar
      rtol(i) = 0.1*rtol(i)
      atol(i) = 0.1*atol(i)
    end do

    do j=1,nadj
      do i=1,nvar
        rtol_adj(i,j) = rtol(i)
        atol_adj(i,j) = atol(i)
      end do
    end do

    ! initialization
    do i = 1,nfix
      fix(i) = 1.0d-8 * cfactor
    end do
    fix(1) = (1.25d+8)*cfactor
    ! constant rate coefficients
    rconst(4) = 9.3d-12
    rconst(11) = 2.2d-10
    rconst(18) = 1.3d-21
    rconst(21) = 4.39999d-40
    rconst(24) = 6.6d-12
    rconst(25) = 1d-20
    rconst(36) = 2.2d-13
    rconst(37) = 1d-11
    rconst(41) = 6.3d-16
    rconst(44) = 2.5d-15
    rconst(49) = 2d-12
    rconst(50) = 6.5d-12
    rconst(52) = 8.1d-13
    rconst(54) = 1600.d0
    rconst(55) = 1.5d-11
    rconst(59) = 7.7d-15
    rconst(64) = 8.1d-12
    rconst(65) = 4.2d0
    rconst(66) = 4.1d-11
    rconst(67) = 2.2d-11
    rconst(68) = 1.4d-11
    rconst(70) = 3d-11
    rconst(73) = 1.7d-11
    rconst(75) = 1.8d-11
    rconst(76) = 9.6d-11
    rconst(77) = 1.2d-17
    rconst(78) = 3.2d-13
    rconst(79) = 8.1d-12
    rconst(81) = 6.8d-13
    var(1)  = 3.652686375061191d-02
    var(2)  = 3.470772014720235d+11
    var(3)  = 2.558736854382598d+03
    var(4)  = 3.353179955426548d-23
    var(5)  = 5.287698343141944d-20
    var(6)  = 1.307511063007438d+07
    var(7)  = 0.0
    var(8)  = 1.985180650951722d+01
    var(9)  = 3.839565210778183d+08
    var(10) = 2.675425503948970d+08
    var(11) = 1.593461800730453d-24
    var(12) = 1.191924698309334d+12
    var(13) = 4.894101892069533d-05
    var(14) = 5.462490511678749d-21
    var(15) = 1.943997577140665d-23
    var(16) = 2.329863421033919d+12
    var(17) = 2.106517721289787d-30
    var(18) = 5.518957651363696d+08
    var(19) = 2.633932009302457d-21
    var(20) = 8.118206770655152d+03
    var(21) = 3.000632019118241d+10
    var(22) = 0.0d0
    var(23) = 0.0d0
    var(24) = 2.687919546821573d+02
    var(25) = 2.061185883675469d+12
    var(26) = 1.561906017211229d+10
    var(27) = 3.671987743784878d+07
    var(28) = 9.888485531199334d+08
    var(29) = 1.248793414327174d+04
    var(30) = 3.357869003391008d+07
    var(31) = 2.727638922137502d+09
    var(32) = 6.534644475169904d+00

    tstart = 12.0d0 * 3600.d0
    tend   = tstart + 72.d0 * 3600.d0
    temp   = 298.d0
    
    call update_sun(tstart)
    call update_rconst() 
    istate(:)=0
    icntrl(:)=0
    rcntrl(:)=0
!    icntrl(1)=1
    icntrl(3)=5
!    icntrl(9)=1
    !icntrl(3)=3
    ! change initial step size
    !      rcntrl(3)=1d-3
    call cpu_time(t1)
    call integrate_adj(tin=tstart, tout=tend, nvar=nvar, nnzero=nnz, np=np,  &
          nadj=nadj, y=var, lambda=lambda, mu=mu, atol_adj=atol_adj,         &
          rtol_adj=rtol_adj, jacp=jacp, rtol=rtol, atol=atol,fun=fun,jac=jac,&
          adjinit=adjinit, hesstr_vec=hesstr_vec,                            &
          hesstr_vec_f_py=hesstr_vec_f_py, rstatus_u=rstate,rcntrl_u=rcntrl, &
          istatus_u=istate, icntrl_u=icntrl)
    call cpu_time(t2)
    tcpu=t2-t1
    call output_statics(var,tstart,tend,tcpu,istate)
    call output_sens(nvar,np,nadj,lambda,mu)
  end program cbm4_ros_adj_dr
