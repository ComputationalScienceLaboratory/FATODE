!~~~> LAPACK implementation
#ifdef FULL_ALGEBRA
module lapack
      implicit none
      save
      integer :: nvar
      integer, allocatable :: ip(:)
      double precision, allocatable :: fjac(:,:),e(:,:)
      double precision, allocatable :: fjac1(:,:),djdt(:,:)
contains
    subroutine lapack_init(n)
      integer :: n,state
      nvar = n
      allocate(ip(nvar),fjac(nvar,nvar),fjac1(nvar,nvar),djdt(nvar,nvar),e(nvar,nvar),STAT=state)
      if(state .ne. 0) then
        stop 'ALlocation error in lapack_init'
      end if
    end subroutine lapack_init


    subroutine lapack_free
      integer :: state
      deallocate(fjac,fjac1,djdt,STAT=state)
      if(state .ne. 0) then
        stop 'Deallocation error for Jac'
      end if
      deallocate(e,STAT=state)
      if(state .ne. 0) then
        stop 'DeaLlocation error for e'
      end if
    end subroutine lapack_free


!~~~> decomposition of the left hand side (HGamma-Jac)
    subroutine lapack_decomp(hgamma,ising)
      integer :: ising, i, j
      double precision :: hgamma
! prepare left hand side
      do j=1,nvar
        do i=1,nvar
          e(i,j) = -fjac(i,j)
        end do
        e(j,j) = e(j,j) + hgamma
      end do
      call dgetrf( nvar, nvar, e, nvar, ip, ising )
    end subroutine lapack_decomp

!~~~>Solving the system (HGamma-Jac)*X = RHS
    subroutine lapack_solve(trans, rhs)
      double precision ::rhs(nvar)
      logical :: trans
      integer :: ising
      if(trans) then 
        call dgetrs( 't', nvar, 1, e, nvar, ip, rhs, nvar, ising)
      else
        call dgetrs( 'n', nvar, 1, e, nvar, ip, rhs, nvar, ising)
      end if
    end subroutine lapack_solve

end module lapack
#endif

#ifdef SPARSE_UMF
! UMFPACK implementation 
module umf
    implicit none 
    save
    integer :: nvar, nnz 
! ax(nnz),ai(nnz),ap(nvar+1)
    double precision, allocatable :: ax(:)
    double precision, allocatable :: fjac(:,:),fjac1(:,:),djdt(:,:)
    integer, allocatable :: ai(:),ap(:)
    double precision ::  control(20),info(90)
    integer*8 :: symbolic,numeric
contains
!~~~>Decomposition of left hand side (HGamma-Jac)
  subroutine umf_decomp(hgamma,ising)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,found,ising
    double precision :: hgamma
!   free previously used memory
    if(symbolic .ne. 0) then
      call umf4fsym(symbolic)
      call umf4fnum(numeric)
    end if

    idx = 1
    ! c index starts from 0, while fortran one starts from 1
    do j = 1,nvar
      found = 0
      do i = 1, nvar
        if ( fjac(i,j) .ne. 0d0) then
  !        ax(idx) = e(i,j)
          if( i .eq. j ) then
            ax(idx) = -fjac(i,j) + hgamma
          else 
            ax(idx) = -fjac(i,j)
          end if      
          ! write the row in the row indices vector ai
          ai(idx) = i - 1
          if(found .eq. 0) then
            found = idx
          end if
          idx = idx + 1
        end if
      end do
      !must also update the column pointers vector ap
      ap(j) = found - 1
    end do
    !last element in ap must be nnza+1
    ap(nvar+1) = nnz
    call umf4def(control)
    control(1) = 0
    call umf4sym(nvar, nvar, ap, ai, ax, symbolic, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sym: ', info(1)
      stop
    endif
    call umf4num(ap, ai, ax, symbolic, numeric, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4num: ', info(1)
      stop
    endif
    ising = 0
  end subroutine umf_decomp

!~~~>Solving the system (HGamma-Jac)*X=RHS
    subroutine umf_solve(trans, rhs)
    ! solve ax=b, without iterative refinement
    integer :: sys
    double precision :: x(nvar),rhs(nvar)
    logical :: trans
    if(trans) then
      sys = 1
    else
      sys = 0
    end if
    call umf4sol(sys, x, rhs, numeric, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sol: ', info(1)
      stop
    endif
    ! free the numeric factorization, to be continued
    rhs(:) = x(:)
    end subroutine umf_solve

      subroutine umf_init(n,nn)
      integer :: n,nn,state
      nvar = n
      nnz = nn
      symbolic = 0
      numeric = 0
      allocate(ax(nnz),ai(nnz),ap(nvar+1),fjac(nvar,nvar),fjac1(nvar,nvar),&
                                                djdt(nvar,nvar),STAT=state)
      if(state .ne. 0) then
        stop 'Allocation error in umf_init'
      end if
      end subroutine umf_init

!~~~>Clean objects created during the process
      subroutine umf_free
      integer :: state
      call umf4fsym(symbolic)
      call umf4fnum(numeric)
      deallocate(ax,ai,ap,fjac,fjac1,djdt,STAT=state)
      if(state .ne. 0) then
        stop 'Deallocation error in umf_free'
      end if
      end subroutine umf_free

end module umf
#endif

#ifdef SPARSE_LU
!~~~>SuperLU implementation
module superlu
    implicit none
    save
    integer :: nvar, nnz
    double precision,allocatable ::fjac(:,:),fjac1(:,:),djdt(:,:)
! ax(nnz),b(nvar),ai(nnz),ap(nvar+1)
    double precision,allocatable :: ax(:),b(:)
    integer,allocatable :: ai(:), ap(:)
    integer*8 factors
contains
!~~~>Decomposition of the left hand side (HGamma-Jac)
    subroutine superlu_decomp(hgamma,ising)
      !convert the matrix e from coordinate format to column compressed format
      integer :: i,j,idx,found,ising
      integer :: iopt, ldb, nrhs, tr
      double precision :: hgamma 
      nrhs = 1
      ldb = nvar
      tr = 0 
      if(factors .ne. 0) then
        iopt = 3
        call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, ising )
      end if

      idx = 1
      do j = 1,nvar
        found = 0
        do i = 1, nvar
          if ( fjac(i,j) .ne. 0d0) then
!            ax(idx) = e(i,j)
            if( i .eq. j ) then
              ax(idx) = -fjac(i,j) + hgamma
            else
            ax(idx) = -fjac(i,j)
            end if
            ! write the row in the row indices vector ai
            ai(idx) = i 
            if(found .eq. 0) then
              found = idx
            end if
            idx = idx + 1
          end if
        end do
        !must also update the column pointers vector ap
        ap(j) = found 
      end do
      !last element in ap must be nnza+1
      ap(nvar+1) = nnz + 1
  !   factorize the matrix. The factors are stored in *factors* handle.
      iopt = 1
      call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                          b, ldb, factors, ising)
  
      if(ising .ne. 0) then
        write(*,*) 'INFO from failed factorization = ', ising
      endif
    end subroutine superlu_decomp


!~~~>Solving the system (HGamma-Jac)*X=RHS
    subroutine superlu_solve(trans, rhs)
    double precision :: rhs(nvar)
    integer :: info, iopt, ldb, nrhs, tr
    logical :: trans
!   solve the system using the existing factors.
    iopt = 2
    nrhs = 1
    ldb = nvar
    if(trans) then
      tr = 1
    else
      tr = 0
    end if
    call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                          rhs, ldb, factors, info )
    if(info .ne. 0) then
      write(*,*) 'INFO from failed triangular solve = ', info
    end if
    end subroutine superlu_solve

      subroutine superlu_init(n,nn)
      integer :: n,nn,state
      nvar = n
      nnz = nn
      factors = 0
      allocate(ax(nnz),b(nvar),ai(nnz),ap(nvar+1),fjac(nvar,nvar),fjac1(nvar,nvar),djdt(nvar,nvar),STAT=state)
      if(state .ne. 0) then
        stop 'Allocation error for Jac'
      end if
      end subroutine superlu_init

      subroutine superlu_free
      integer :: info, iopt, ldb, nrhs, tr, state
      iopt = 3
      ldb = nvar
      nrhs = 1
      tr = 0
      call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, info )
      deallocate(ax,b,ai,ap,fjac,fjac1,djdt,STAT=state)
      if(state .ne. 0) then
        stop 'Deallocation error in superlu_free'
      end if
      end subroutine superlu_free

end module superlu
#endif

module ls_solver
#ifdef SPARSE_UMF
      use umf
#endif
#ifdef FULL_ALGEBRA
      use lapack
#endif
#ifdef SPARSE_LU
      use superlu
#endif
      implicit none
contains      
      subroutine lss_decomp(hgamma,ising)
      implicit none
      integer :: ising
      double precision :: hgamma
#ifdef SPARSE_UMF
      call umf_decomp(hgamma,ising)
#endif
#ifdef FULL_ALGEBRA
      call lapack_decomp(hgamma,ising)
#endif
#ifdef SPARSE_LU
      call superlu_decomp(hgamma,ising)
#endif
      end subroutine lss_decomp


      subroutine lss_solve(trans, rhs)
      double precision ::rhs(nvar)
      logical :: trans
#ifdef SPARSE_UMF  
      call umf_solve(trans, rhs)
#endif
#ifdef FULL_ALGEBRA
      call lapack_solve(trans, rhs)
#endif
#ifdef SPARSE_LU
      call superlu_solve(trans, rhs)
#endif
!      print *, rhs  
      end subroutine lss_solve    


      subroutine lss_init(n,nn)     
      integer :: n,nn
#ifdef SPARSE_UMF
      call umf_init(n,nn)
#endif
#ifdef FULL_ALGEBRA
      call lapack_init(n)
#endif
#ifdef SPARSE_LU
      call superlu_init(n,nn)
#endif
     end subroutine lss_init

      subroutine lss_free
#ifdef SPARSE_UMF
      call umf_free
#endif
#ifdef FULL_ALGEBRA
      call lapack_free
#endif
#ifdef SPARSE_LU
      call superlu_free
#endif
     end subroutine lss_free

     subroutine lss_jac_time_derivative(delta)
     double precision :: delta
     double precision :: a
     a = -1
     call daxpy(nvar*nvar,a,fjac,1,djdt,1)
     a = 1/delta
     call dscal(nvar*nvar,a,djdt,1)
     end subroutine lss_jac_time_derivative

!~~~> multiplication of the jacobian matrix and a vector
     subroutine lss_mul_jac(z,g)
      double precision :: z(nvar), g(nvar)
      z = matmul(fjac,g)
     end subroutine lss_mul_jac

!~~~> multiplication of the jacobian matrix and a vector
     subroutine lss_mul_jac1(z,g)
      double precision :: z(nvar), g(nvar)
      z = matmul(fjac1,g)
     end subroutine lss_mul_jac1

!~~~> multiplication of the jacobian matrix and a vector
     subroutine lss_mul_jac2(z,g)
      double precision :: z(nvar), g(nvar)
      z = matmul(djdt,g)
     end subroutine lss_mul_jac2
     
     subroutine lss_jac(t,y,jac)
     double precision :: t, y(nvar)
     external :: jac
     call jac(nvar,t,y,fjac)
     end subroutine lss_jac

     subroutine lss_jac1(t,y,jac)
     double precision :: t, y(nvar)
     external :: jac
     call jac(nvar,t,y,fjac1)
     end subroutine lss_jac1

     subroutine lss_jac2(t,y,jac)
     double precision :: t, y(nvar)
     external :: jac
     call jac(nvar,t,y,djdt)
     end subroutine lss_jac2

end module ls_solver
