#ifdef FULL_ALGEBRA
module lapack
      implicit none
      save
      integer :: nvar
      double precision,allocatable :: fjac(:,:),e1(:,:)
      integer, allocatable :: ip1(:), ip2(:)
      complex(kind=selected_real_kind(14,300)),allocatable :: e2(:,:)
contains
      subroutine lapack_decomp(ising,hgamma)
      integer :: ising,i,j
      double precision :: hgamma
      do j=1,nvar
        do i=1,nvar
          e1(i,j) = -fjac(i,j)
        end do
        e1(j,j) = e1(j,j) +hgamma
      end do
      call dgetrf( nvar, nvar, e1, nvar, ip1, ising )
      end subroutine lapack_decomp

      subroutine lapack_decomp_cmp(ising,alpha,beta)
      integer :: ising,i,j
      double precision :: alpha,beta
      do j=1,nvar
         do i=1,nvar
            e2(i,j) = dcmplx( -fjac(i,j), 0d0)
         end do
         e2(j,j) =e2(j,j) + cmplx( alpha, beta )
      end do
      call zgetrf( nvar, nvar, e2, nvar, ip2, ising )
      end subroutine lapack_decomp_cmp

      subroutine lapack_solve(rhs)
      double precision ::rhs(nvar)
      integer :: ising
      call dgetrs( 'n', nvar, 1, e1, nvar, ip1, rhs, nvar, ising)
      end subroutine lapack_solve

      subroutine lapack_solve_cmp(b,bz)
      double precision :: b(nvar), bz(nvar)
      complex(kind=selected_real_kind(14,300)) :: rhs(nvar)
      integer :: ising,i
      do i=1,nvar
         rhs(i) = dcmplx(b(i),bz(i))
      end do
      call zgetrs( 'n', nvar, 1, e2, nvar, ip2, rhs, nvar, ising)
      do i = 1,nvar
        b(i) = dble(rhs(i))
        bz(i) = aimag(rhs(i))
      end do
      end subroutine lapack_solve_cmp


      subroutine lapack_init(n)
      integer :: n,state
      nvar = n
      allocate(ip1(nvar),ip2(nvar),fjac(nvar,nvar),e1(nvar,nvar),&
                                       e2(nvar,nvar),STAT = state)
      if(state .ne. 0)stop 'Allocation error in lapack_init'
      end subroutine lapack_init

      subroutine lapack_free
      integer :: state
      deallocate(ip1,ip2,fjac,e1,e2,STAT=state)
      if(state .ne.0) stop 'Deallocation error in lapack_free'
      end subroutine lapack_free

end module lapack
#endif FULL_ALGEBRA

#ifdef SPARSE_UMF
module umf
    implicit none 
    save 
    integer :: nvar, nnz
    double precision, allocatable :: fjac(:,:)
    double precision, allocatable :: ax(:),axc(:),azc(:)
! ax(nnz),axc(nnz),azc(nnz),ai(nnz),ap(nvar+1)
    integer, allocatable :: ai(:), ap(:)
    double precision ::  control(20), info(90), zcontrol(20)
    integer*8 :: symbolic, numeric, symbolic_cmp, numeric_cmp
contains
    subroutine umf_decomp(ising,hgamma)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,found,ising
    double precision :: hgamma
!   free previously used memory
    if(symbolic .ne. 0) then
      call umf4fsym(symbolic)
      call umf4fnum(numeric)
    end if

    idx = 1
    do j = 1,nvar
      found = 0
      do i = 1, nvar
        if ( fjac(i,j) .ne. 0d0) then
!          ax(idx) = e(i,j)
          if( i .eq. j ) then
            ax(idx) = -fjac(i,j) + hgamma
          else 
            ax(idx) = -fjac(i,j)
          end if      
          ! write the row in the row indices vector ai
          ai(idx) = i -1
          if(found .eq. 0) then
            found = idx
          end if
          idx = idx + 1
        end if
      end do
      !must also update the column pointers vector ap
      ap(j) = found -1
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

    subroutine umf_solve(rhs)
    implicit none
    ! solve ax=b, without iterative refinement
    integer :: sys
    double precision :: x(nvar),rhs(nvar)
    sys = 0
    call umf4sol(sys, x, rhs, numeric, control, info)

    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sol: ', info(1)
      stop
    endif
    ! free the numeric factorization
    rhs(:) = x(:)
    end subroutine umf_solve
      
    subroutine umf_decomp_cmp(ising,alpha,beta)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,ising
    double precision :: alpha, beta
!   free previously used memory
    if(symbolic_cmp .ne. 0) then
      call umf4fsym(symbolic_cmp)
      call umf4fnum(numeric_cmp)
    end if

    idx = 1
    do j = 1,nvar
      do i = 1, nvar
        if ( fjac(i,j) .ne. 0d0) then
          if( i .eq. j ) then
            axc(idx) = -fjac(i,j) + alpha
            azc(idx) = beta
          else
            axc(idx) = -fjac(i,j)
            azc(idx) = 0d0
          end if
          idx = idx + 1
        end if
      end do
    end do

    call umf4zdef(zcontrol)
    zcontrol(1) = 0
    call umf4zsym(nvar, nvar, ap, ai, axc, azc, symbolic_cmp, zcontrol, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sym: ', info(1)
      stop
    endif

    call umf4znum(ap, ai, axc, azc, symbolic_cmp, numeric_cmp, zcontrol, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4num: ', info(1)
      stop
    endif

!    call umf4zfsym(symbolic_cmp)
    ising = 0
    end subroutine umf_decomp_cmp

    subroutine umf_solve_cmp(b,bz)
    ! solve ax=b, without iterative refinement complex version
    integer :: sys
    double precision :: x(nvar),xz(nvar),b(nvar),bz(nvar)
    
    sys = 0
    call umf4zsol(sys, x, xz, b, bz, numeric_cmp, zcontrol, info)

    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4zsol: ', info(1)
      stop
    endif

    b(:) = x(:)
    bz(:) = xz(:)
    end subroutine umf_solve_cmp

 
    subroutine umf_init(n,nn)
    integer :: n,nn,state
    nvar = n
    nnz = nn
    allocate(ax(nnz),axc(nnz),azc(nnz),ai(nnz),ap(nvar+1),&
                              fjac(nvar,nvar),STAT = state)
    if(state .ne. 0) stop 'Allocation error in umf_init'
    symbolic = 0
    numeric = 0
    symbolic_cmp = 0
    numeric_cmp = 0
    end subroutine umf_init

    subroutine umf_free
    implicit none
    integer :: state

    call umf4fsym(symbolic)
    call umf4fnum(numeric)
    call umf4zfsym(symbolic_cmp)
    call umf4zfnum(numeric_cmp)
    deallocate(ax,axc,azc,ai,ap,fjac,STAT = state)
    if(state .ne. 0) stop 'Deallocation error in umf_free'
    end subroutine

end module umf
#endif

#ifdef SPARSE_LU
module superlu
    implicit none
    save
    integer :: nvar, nnz
    double precision, allocatable :: fjac(:,:)
    double precision, allocatable :: ax(:),b(:)
    !ax(nnz),b(nvar),axc(nnz),bc(nvar),ai(nnz),ap(nvar+1)
    complex(kind=selected_real_kind(14,300)),allocatable :: axc(:), bc(:)
    integer, allocatable :: ai(:), ap(:)
    integer*8 :: factors,factors_cmp
contains
    subroutine superlu_decomp(ising,hgamma)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,found,ising
    integer :: tr, iopt, ldb, nrhs
    double precision :: hgamma
!   free previouslu used memory
    tr = 0
    nrhs = 1
    ldb = nvar
    if(factors .ne. 0) then
      iopt = 3
      call c_fortran_dgssv(tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, ising )
    end if

    idx = 1
    do j = 1,nvar
      found = 0
      do i = 1, nvar
        if ( fjac(i,j) .ne. 0d0) then
!          ax(idx) = e(i,j)
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
                          b, ldb, factors, ising )
    if(ising .ne. 0) then
      write(*,*) 'INFO from failed factorization = ', ising
    endif
    end subroutine superlu_decomp

    subroutine superlu_solve(rhs)
    double precision :: rhs(nvar)
    integer :: info, tr, iopt, ldb, nrhs
!   solve the system using the existing factors.
    iopt = 2
    tr = 0
    ldb = nvar
    nrhs = 1
    call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                          rhs, ldb, factors, info )
    if(info .ne. 0) then
      write(*,*) 'INFO from failed triangular solve = ', info
    endif
    end subroutine superlu_solve

    subroutine superlu_decomp_cmp(ising,alpha,beta)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,ising
    integer :: tr, iopt, ldb, nrhs
    double precision :: alpha, beta
!   free previouslu used memory
    nrhs = 1
    ldb = nvar
    tr = 0
    if(factors_cmp .ne. 0) then
      iopt = 3
      call c_fortran_zgssv( tr, iopt, nvar, nnz, nrhs, axc, ai, ap,&
                        bc, ldb, factors_cmp, ising )
    end if

    idx = 1
    do j = 1,nvar
      do i = 1, nvar
        if ( fjac(i,j) .ne. 0d0) then
          if( i .eq. j ) then
            axc(idx) = -fjac(i,j) + cmplx(alpha,beta)
          else
            axc(idx) = dcmplx(-fjac(i,j),0d0)
          end if
          idx = idx + 1
        end if
      end do
    end do

!   factorize the matrix. The factors are stored in *factors* handle.
    iopt = 1
    call c_fortran_zgssv( tr, iopt, nvar, nnz, nrhs, axc, ai, ap,&
                          bc, ldb, factors_cmp, ising )
    if(ising .ne. 0) then
      write(*,*) 'INFO from failed factorization = ', ising
    endif
    end subroutine superlu_decomp_cmp

    subroutine superlu_solve_cmp(rhs,zrhs)
    integer :: info, tr, iopt, ldb, nrhs
    integer :: i
    double precision :: rhs(nvar),zrhs(nvar)
!   solve the system using the existing factors.
    do i=1,nvar
      bc(i) = dcmplx(rhs(i),zrhs(i))
    end do
    iopt = 2
    tr = 0
    ldb = nvar
    nrhs = 1
    call c_fortran_zgssv( tr, iopt, nvar, nnz, nrhs, axc, ai, ap,&
                          bc, ldb, factors_cmp, info )
    do i = 1,nvar
      rhs(i) = dble(bc(i))
      zrhs(i) = aimag(bc(i))
    end do

    if(info .ne. 0) then
      write(*,*) 'INFO from failed triangular solve = ', info
    endif
    end subroutine superlu_solve_cmp


    subroutine superlu_init(n,nn)
    integer :: n,nn,state
    nvar = n
    nnz = nn
    factors = 0
    factors_cmp = 0
    allocate(ax(nnz),b(nvar),axc(nnz),bc(nvar),ai(nnz),ap(nvar+1),&
                           fjac(nvar,nvar),STAT=state)
    if(state .ne. 0) stop 'Allocation error in superlu_init'
    end subroutine superlu_init

    subroutine superlu_free
    integer :: state
    integer :: info, tr, iopt, ldb, nrhs

    tr = 0
    iopt = 3
    ldb = nvar
    nrhs = 1
    call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, info )
    iopt = 3
    call c_fortran_zgssv( tr, iopt, nvar, nnz, nrhs, axc, ai, ap,&
                       bc, ldb, factors_cmp, info )
    deallocate(ax,b,axc,bc,ai,ap,fjac,STAT=state)
    if(state .ne. 0) stop 'Deallocation error in superlu_free'
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
      public
contains      
      subroutine lss_jac(t,y,jac)
      double precision :: t,y(nvar)
      external :: jac
      call jac(nvar,t,y,fjac)
      end subroutine lss_jac

      subroutine lss_decomp(hgamma, ising)
      implicit none
      integer :: ising
      double precision :: hgamma
#ifdef SPARSE_UMF      
      call umf_decomp(ising, hgamma)
#endif
#ifdef FULL_ALGEBRA
      call lapack_decomp(ising,hgamma)
#endif
#ifdef SPARSE_LU
      call superlu_decomp(ising,hgamma)
#endif
      end subroutine lss_decomp

      subroutine lss_decomp_cmp(alpha, beta, ising)
      integer :: ising
      double precision :: alpha,beta

#ifdef SPARSE_UMF
      call umf_decomp_cmp(ising, alpha, beta)
#endif
#ifdef FULL_ALGEBRA
      call lapack_decomp_cmp(ising,alpha,beta)
#endif
#ifdef SPARSE_LU
      call superlu_decomp_cmp(ising,alpha,beta)
#endif
      end subroutine lss_decomp_cmp

      subroutine lss_solve(rhs)
      double precision ::rhs(nvar)
#ifdef SPARSE_UMF
      call umf_solve(rhs)
#endif 
#ifdef FULL_ALGEBRA
      call lapack_solve(rhs)
#endif
#ifdef SPARSE_LU
      call superlu_solve(rhs)
#endif
!      print *, rhs  
      end subroutine lss_solve

      subroutine lss_solve_cmp(b,bz)
      double precision  :: b(nvar), bz(nvar)
#ifdef SPARSE_UMF
      call umf_solve_cmp(b,bz)
#endif
#ifdef FULL_ALGEBRA
      call lapack_solve_cmp(b,bz)
#endif
#ifdef SPARSE_LU
      call superlu_solve_cmp(b,bz)
#endif
!      print *, rhs  
      end subroutine lss_solve_cmp

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
      end subroutine

end module ls_solver
