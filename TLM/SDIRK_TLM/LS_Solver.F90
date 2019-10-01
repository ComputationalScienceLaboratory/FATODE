#ifdef FULL_ALGEBRA
!~~~> LAPACK implementation
module lapack
      implicit none
      save
      integer :: nvar
! ip(nvar),ip_tlm(nvar)
      integer, allocatable :: ip(:),ip_tlm(:)
      double precision, allocatable :: fjac(:,:),e(:,:)
      double precision, allocatable :: fjac1(:,:),e_tlm(:,:)
contains
    subroutine lapack_init(n)
      integer :: n,state
      allocate(ip(nvar),ip_tlm(nvar),fjac(nvar,nvar),fjac1(nvar,nvar),&
                          e(nvar,nvar),e_tlm(nvar,nvar),STAT=state)
      if(state .ne. 0) then
        stop 'ALlocation error in lapack_init'
      end if
    end subroutine lapack_init

    subroutine lapack_free
      integer :: state
      deallocate(ip,ip_tlm,fjac,fjac1,e,e_tlm,STAT=state)
      if(state .ne. 0) then
        stop 'DeaLlocation error in lapack_free'
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

!~~~> decomposition of the left hand side (HGamma-Jac)
    subroutine lapack_decomp_tlm(hgamma,ising)
      integer :: ising, i, j
      double precision :: hgamma
! prepare left hand side
      do j=1,nvar
        do i=1,nvar
          e_tlm(i,j) = -fjac1(i,j)
        end do
        e_tlm(j,j) = e_tlm(j,j) + hgamma
      end do
      call dgetrf( nvar, nvar, e_tlm, nvar, ip_tlm, ising )
    end subroutine lapack_decomp_tlm

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

!~~~>Solving the system (HGamma-Jac)*X = RHS
    subroutine lapack_solve_tlm(trans, rhs)
      double precision ::rhs(nvar)
      logical :: trans
      integer :: ising
      if(trans) then
        call dgetrs( 't', nvar, 1, e_tlm, nvar, ip_tlm, rhs, nvar, ising)
      else
        call dgetrs( 'n', nvar, 1, e_tlm, nvar, ip_tlm, rhs, nvar, ising)
      end if
    end subroutine lapack_solve_tlm

end module lapack
#endif

#ifdef SPARSE_UMF
! UMFPACK implementation 
module umf
    implicit none 
    save 
    integer :: nvar, nnz
! ax(nnz),ax_tlm(nnz),ai(nnz),ap(nvar+1),ai_tlm(nnz),ap_tlm(nvar+1)
    double precision,allocatable :: ax(:),ax_tlm(:)
    double precision,allocatable :: fjac(:,:),fjac1(:,:)
    integer,allocatable :: ai(:),ap(:),ai_tlm(:),ap_tlm(:)
    double precision ::  control(20),info(90)
    integer*8 :: symbolic,numeric,symbolic_tlm,numeric_tlm
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

  subroutine umf_decomp_tlm(hgamma,ising)
    !convert the matrix e from coordinate format to column compressed format
    integer :: i,j,idx,found,ising
    double precision :: hgamma
!   free previously used memory
    if(symbolic_tlm .ne. 0) then
      call umf4fsym(symbolic_tlm)
      call umf4fnum(numeric_tlm)
    end if

    idx = 1
    ! c index starts from 0, while fortran one starts from 1
    do j = 1,nvar
      found = 0
      do i = 1, nvar
        if ( fjac1(i,j) .ne. 0d0) then
  !        ax(idx) = e(i,j)
          if( i .eq. j ) then
            ax_tlm(idx) = -fjac1(i,j) + hgamma
          else
            ax_tlm(idx) = -fjac1(i,j)
          end if
          ! write the row in the row indices vector ai
          ai_tlm(idx) = i - 1
          if(found .eq. 0) then
            found = idx
          end if
          idx = idx + 1
        end if
      end do
      !must also update the column pointers vector ap
      ap_tlm(j) = found - 1
    end do
    !last element in ap must be nnza+1
    ap_tlm(nvar+1) = nnz
    call umf4def(control)
    control(1) = 0
    call umf4sym(nvar, nvar, ap_tlm, ai_tlm, ax_tlm, symbolic_tlm, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sym: ', info(1)
      stop
    endif
    call umf4num(ap_tlm, ai_tlm, ax_tlm, symbolic_tlm, numeric_tlm, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4num: ', info(1)
      stop
    endif
    ising = 0
  end subroutine umf_decomp_tlm

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

!~~~>Solving the system (HGamma-Jac)*X=RHS
    subroutine umf_solve_tlm(trans, rhs)
    ! solve ax=b, without iterative refinement
    integer :: sys
    double precision :: x(nvar),rhs(nvar)
    logical :: trans
    if(trans) then
      sys = 1
    else
      sys = 0
    end if
    call umf4sol(sys, x, rhs, numeric_tlm, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sol: ', info(1)
      stop
    endif
    ! free the numeric factorization, to be continued
    rhs(:) = x(:)
    end subroutine umf_solve_tlm

      subroutine umf_init(n,nn)
      integer :: n,nn,state
      nvar = n
      nnz = nn
      symbolic = 0
      numeric = 0
      symbolic_tlm = 0
      numeric_tlm = 0
      allocate(ax(nnz),ax_tlm(nnz),ai(nnz),ap(nvar+1),ai_tlm(nnz),&
        ap_tlm(nvar+1),fjac(nvar,nvar),fjac1(nvar,nvar),STAT=state)
      if(state .ne. 0) then
        stop 'Allocation error in umf_init'
      end if
      end subroutine umf_init

!~~~>Clean objects created during the process
      subroutine umf_free
      integer :: state

      call umf4fsym(symbolic)
      call umf4fnum(numeric)
      if(symbolic_tlm .ne. 0) then
        call umf4fsym(symbolic_tlm)
        call umf4fnum(numeric_tlm)
      end if
      deallocate(ax,ax_tlm,ai,ap,ai_tlm,ap_tlm,fjac,fjac1,STAT=state)
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
    double precision,allocatable ::fjac(:,:),fjac1(:,:)
!ax(nnz),ax_tlm(nnz),b(nvar),ai(nnz), ai_tlm(nnz),ap(nvar+1),ap_tlm(nvar+1)
    double precision,allocatable :: ax(:),ax_tlm(:),b(:)
    integer,allocatable :: ai(:), ai_tlm(:),ap(:),ap_tlm(:)
    integer*8 factors,factors_tlm
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
                          b, ldb, factors, ising )
  
      if(ising .ne. 0) then
        write(*,*) 'INFO from failed factorization = ', ising
      endif
    end subroutine superlu_decomp

    subroutine superlu_decomp_tlm(hgamma,ising)
      implicit none
      !convert the matrix e from coordinate format to column compressed format
      integer :: i,j,idx,found,ising
      integer :: iopt, ldb, nrhs, tr
      double precision :: hgamma
      nrhs = 1
      ldb = nvar
      tr = 0
      if(factors_tlm .ne. 0) then
        iopt = 3
        call c_fortran_dgssv(tr,iopt,nvar,nnz,nrhs,ax_tlm,ai_tlm,ap_tlm,&
                       b, ldb, factors_tlm, ising )
      end if

      idx = 1
      do j = 1,nvar
        found = 0
        do i = 1, nvar
          if ( fjac1(i,j) .ne. 0d0) then
!            ax(idx) = e(i,j)
            if( i .eq. j ) then
              ax_tlm(idx) = -fjac1(i,j) + hgamma
            else
            ax_tlm(idx) = -fjac1(i,j)
            end if
            ! write the row in the row indices vector ai
            ai_tlm(idx) = i 
            if(found .eq. 0) then
              found = idx
            end if
            idx = idx + 1
          end if
        end do
        !must also update the column pointers vector ap
        ap_tlm(j) = found 
      end do
      !last element in ap must be nnza+1
      ap_tlm(nvar+1) = nnz + 1
  !   factorize the matrix. The factors are stored in *factors* handle.
      call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax_tlm, ai_tlm, ap_tlm,&
                          b, ldb, factors_tlm, ising )

      if(ising .ne. 0) then
        write(*,*) 'INFO from failed factorization = ', ising
      endif
    end subroutine superlu_decomp_tlm

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

!~~~>Solving the system (HGamma-Jac)*X=RHS
    subroutine superlu_solve_tlm(trans, rhs)
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
    call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax_tlm, ai_tlm, ap_tlm,&
                          rhs, ldb, factors_tlm, info )
    if(info .ne. 0) then
      write(*,*) 'INFO from failed triangular solve = ', info
    end if
    end subroutine superlu_solve_tlm

      subroutine superlu_init(n,nn)
      integer :: n,nn,state
      nvar = n
      nnz = nn
      factors = 0
      factors_tlm = 0
      allocate(ax(nnz),ax_tlm(nnz),b(nvar),ai(nnz),ai_tlm(nnz),ap(nvar+1),&
                ap_tlm(nvar+1),fjac(nvar,nvar),fjac1(nvar,nvar),STAT=state)
      if(state .ne. 0) then
        stop 'Allocation error in superlu_init'
      end if
      end subroutine superlu_init

      subroutine superlu_free
      integer :: state
      integer :: info, iopt, ldb, nrhs, tr     

      iopt = 3
      nrhs = 1
      tr = 0
      ldb = nvar
      call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b,ldb,factors,info )
      if (factors_tlm .ne. 0) then     
        call c_fortran_dgssv(tr,iopt,nvar,nnz,nrhs,ax_tlm,ai_tlm,ap_tlm,&
                       b,ldb,factors_tlm,info )
      end if
      deallocate(ax,ax_tlm,b,ai,ai_tlm,ap,ap_tlm,fjac,fjac1,STAT=state)
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

      subroutine lss_decomp_tlm(hgamma,ising)
      integer :: ising
      double precision :: hgamma
#ifdef SPARSE_UMF
      call umf_decomp_tlm(hgamma,ising)
#endif
#ifdef FULL_ALGEBRA
      call lapack_decomp_tlm(hgamma,ising)
#endif
#ifdef SPARSE_LU
      call superlu_decomp_tlm(hgamma,ising)
#endif
      end subroutine lss_decomp_tlm

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

      subroutine lss_solve_tlm(trans, rhs)
      double precision ::rhs(nvar)
      logical :: trans
#ifdef SPARSE_UMF  
      call umf_solve_tlm(trans, rhs)
#endif
#ifdef FULL_ALGEBRA
      call lapack_solve_tlm(trans, rhs)
#endif
#ifdef SPARSE_LU
      call superlu_solve_tlm(trans, rhs)
#endif
!      print *, rhs  
      end subroutine lss_solve_tlm    


      subroutine lss_init(n,nn)
      integer :: n ,nn     
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

     subroutine lss_mul_jac1(z,g)
      double precision :: z(nvar), g(nvar)
      z=matmul(fjac1,g)
     end subroutine lss_mul_jac1

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

end module ls_solver
