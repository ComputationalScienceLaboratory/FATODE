!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine hesstr_vec(n,t,y,u,k,tmp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      use swe2dxy_parameters
      integer :: n
      double precision :: t,y(n),u(n),k(n),tmp(n)
      double precision :: y_grid(3,mx+4,my+4),u_grid(3,mx+4,my+4),&
                          k_grid(3,mx+4,my+4),tmp_grid(3,mx+4,my+4)
      call vec2grid(y,y_grid)
      call vec2grid(u,u_grid)
      call vec2grid(k,k_grid)
      call vec2grid(tmp,tmp_grid)
      call g_g_compute_f(y_grid,u_grid,k_grid,tmp_grid)
      call grid2vec(tmp_grid,tmp)
      end subroutine hesstr_vec

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

      write(6,2010) tstart,tend,istate(3:5),tcpu
 2010   FORMAT (2x,'Time interval [',f5.2,'-',f5.2, ']',/,&
          2x,'Nsteps=',I5,' NAsteps=',I5,' NRsteps=',I5,/,&
               2x, 'CPUtime=', f14.3,' seconds',//)
! store results into a file
!      open(unit=11,file='swe2d_ros_sol.txt',status='replace')
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
      character(len=12) :: fname
      write(fname,'(a,i3.3,a)') 'RS_ADJ_',fno,'.m'
      open(unit=20,file=fname,status='replace')
      do j=1,nadj
         do i=1,ndim
            write(20,'(E24.16)') y_adj(i,j)
         end do
      end do
    end subroutine output_sen

  program swe_ros_adj_dr
    use swe2dxy_parameters
    use ros_adj_f90_integrator
    implicit none
    external fun,jac,adjinit,hesstr_vec
    integer, parameter :: nadj = 1
    integer ::i,j,k,state,np = 0
    double precision :: atol(ndim),rtol(ndim), var(ndim)
    double precision :: atol_adj(ndim,nadj),rtol_adj(ndim,nadj),y_adj(ndim,nadj) 
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

       do k=1,nadj
        do i=1,ndim
          rtol_adj(i,k) = rtol(i)
          atol_adj(i,k) = atol(i)
        end do
       end do


       call initializegaussiangrid(temp_par, a)
       temp_par(2,:,:) = 2d0
       temp_par(3,:,:) = 2d0
       call grid2vec(temp_par,var)

       tstart = 0.0d0
       tend = tstart + 5*2.1573758800627289e-003
       istate(:)=0
       icntrl(:)=0
!       icntrl(1)=1
       rcntrl(:)=0
! change initial step size
       rcntrl(3)=1e-4
      
       call cpu_time(t1)
       call integrate_adj(tin=tstart, tout=tend, nvar=ndim,np=np,nnzero=nnz,&
        lambda=y_adj, y=var, rtol=rtol, atol=atol, nadj=nadj,&
        atol_adj=atol_adj, rtol_adj=rtol_adj, fun=fun, jac=jac, adjinit=adjinit, &
        rstatus_u=rstate,rcntrl_u=rcntrl, istatus_u=istate,&
        hesstr_vec=hesstr_vec, icntrl_u=icntrl)
       call cpu_time(t2)
       tcpu=t2-t1
       call output_statics(var,tstart,tend,tcpu,istate)
       call output_sen(ndim,nadj,y_adj,j)
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
 
end program swe_ros_adj_dr
