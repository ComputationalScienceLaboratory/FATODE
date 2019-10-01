!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine jac(n,t,y,fjac)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use swe2dxy_parameters
    implicit none
    integer :: n
    double precision ::t
    double precision ::y(n),fjac(n,n)
    call vec2grid(y,u_par)
    call compute_f(u_par,f_par,feigvalues_par,geigvalues_par)
    call compute_jacf(u_par,feigvalues_par,geigvalues_par)
    fjac(:,:)=jf(:,:)
    end subroutine jac

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine fun(n,t, y, p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      use swe2dxy_parameters
      implicit none
      integer :: n
      double precision :: t
      double precision :: y(n), p(n)
      
      call vec2grid(y, temp_par)
      call compute_f(temp_par,fz_par,feigvalues_par,geigvalues_par)
      call grid2vec(fz_par, p)
      
      end subroutine fun

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
         lambda(k,k) = 1.0d0
       end do
      end subroutine   
   
      subroutine output_statics(var,tstart,tend, tcpu, istate)
      use swe2dxy_parameters
      implicit none
      double precision :: tstart, tend, tcpu
      double precision :: var(ndim)
      integer :: istate(20)
!      double precision :: tmp(ndim),delta(ndim)
!      double precision :: l2, l2tmp,error
!      integer :: i

!      open(unit=2010,file='swe_lsode_sol.txt',form='formatted')
!      do i=1,ndim
!        read(2010,"(e24.16)") tmp(i)
!      end do
!      close(2010)

!      do i=1,ndim
!         delta(i)=var(i)-tmp(i)
!      end do
!     *****errors are shown in the form ||delta(x)||/||x||*****
!      l2=0d0
!      l2tmp=0d0
!      do i=1,ndim
!         l2=l2+delta(i)*delta(i)
!         l2tmp =l2tmp+tmp(i)*tmp(i)
!      end do
!      error=sqrt(l2/l2tmp)
      write(6,2010) tstart,tend,istate(3:5),tcpu
 2010   format (2x,'time interval [',f5.2,'-',f5.2, ']',/,&
          2x,'nsteps=',i5,' nasteps=',i5,' nrsteps=',i5,/,&
               2x, 'cputime=', f10.3,//)

!      write(6,2010) tstart,tend,istate(3:5),tcpu, error
! 2010   format (2x,'time interval [',f5.2,'-',f5.2, ']',/,&
!          2x,'nsteps=',i5,' nasteps=',i5,' nrsteps=',i5,/,&
!               2x, 'cputime=', f10.3,/,2x,'rerror=',f8.5,//)

! store results into a file
!      open(unit=11,file='swe2d_sdirk_sol.txt',status='replace')
!      do i=1,ndim
!        write(11,2011) var(i)
! 2011   format(e24.16)
!      end do
!      close(11)

    end subroutine output_statics

    subroutine output_sen(ndim,nadj,y_adj,fno)
      implicit none
      double precision :: y_adj(ndim,nadj)
      integer       :: ndim,nadj,i,j,fno
      character(len=13) :: fname
      write(fname,'(a,i3.3,a)') 'ERK_ADJ_',fno,'.m'
      open(unit=20,file=fname,status='replace')
      do j=1,nadj
         do i=1,ndim
            write(20,'(E24.16)') y_adj(i,j)
         end do
      end do
    end subroutine output_sen
  

  program swe_sdirk_adj_dr
    use swe2dxy_parameters
    use erk_adj_f90_integrator
    implicit none
    external fun,jac,adjinit
    integer, parameter :: nadj = 1
    integer ::i,j,k,state,np = 0
    double precision :: atol(ndim),rtol(ndim), var(ndim)
    double precision :: y_adj(ndim,nadj) 
    double precision:: tstart, tend, t1, t2, tcpu
    double precision:: rstate(20),rcntrl(20)
    double precision :: a=30d0
    integer ::istate(20),icntrl(20)
    !~~> solution
    !   u(1,:,:) = h(t,x,y)
    !   u(2,:,:) = u(t,x,y)
    !   u(3,:,:) = v(t,x,y)
    allocate(f_par(3, mx+4, my+4), fz_par(3, mx+4, my+4), stat=state)
    if (state .ne. 0) then
      stop 'allocation error for f or fz'
    end if
   
    allocate(temp_par(3, mx+4, my+4), feigvalues_par(3, mx+4, my+4), geigvalues_par(3, mx+4,my+4), stat=state)
    if (state .ne. 0) then
      stop 'allocation error for temp or feigvalues or geigvalues'
    end if

    allocate(jf(ndim, ndim), jg(ndim, ndim), stat=state)
    if (state .ne. 0) then
      stop 'allocation error for jacf, jf or jg'
    end if

     do i=1,ndim
       rtol(i) = 0.1
       atol(i) = 0.1
     end do

     do j = 1,6
       do i=1,ndim
         rtol(i) = 0.1 * rtol(i)
         atol(i) = 0.1 * atol(i)
       end do

       call initializegaussiangrid(temp_par, a)
       temp_par(2,:,:) = 2d0
       temp_par(3,:,:) = 2d0
       call grid2vec(temp_par,var)

       tstart = 0.0d0
       tend = tstart + 5*2.1573758800627289e-003
       istate(:)=0
       icntrl(:)=0
       rcntrl(:)=0
! change initial step size
!       rcntrl(3)=1e-3
      
       call cpu_time(t1)
       call integrate_adj(tin=tstart, tout=tend, nvar=ndim, np=np,&
        lambda=y_adj, adjinit=adjinit,y=var, rtol=rtol, atol=atol, nadj=nadj,&
        fun=fun, jac=jac, rstatus_u=rstate,rcntrl_u=rcntrl, istatus_u=istate,&
         icntrl_u=icntrl)
       call cpu_time(t2)
       tcpu=t2-t1
       call output_statics(var,tstart,tend,tcpu,istate)
       call output_sen(ndim,nadj,y_adj,j)
       print *,'Number of rhs Evaluations:',istate(1)
     end do

     deallocate(f_par, fz_par, stat = state)
     if (state .ne. 0) then
       print *, 'deallocation error of f fz'
       stop 1
     end if
     deallocate(temp_par, feigvalues_par, geigvalues_par, stat = state)
     if (state .ne. 0) then
       print *, 'deallocation error of temp,feigvalues,geigvalues'
       stop 1
     end if

     deallocate(jf, jg, stat = state)
     if (state .ne. 0) then
       print *, 'deallocation error of jf jg'
       stop 1
     end if
 
end program swe_sdirk_adj_dr
