#ifdef FULL_ALGEBRA
!~~~> LAPACK implementation
module lapack
      implicit none
      save
      integer :: nvar
!ip1(nvar),ip2(nvar),ip_big(3*nvar)
      integer,allocatable :: ip1(:),ip2(:),ip_big(:)
      double precision ,allocatable:: fjac(:,:),e1(:,:)
      double precision,allocatable ::e_big(:,:), fjac1(:,:),fjac2(:,:),fjac3(:,:)
      complex(kind=selected_real_kind(14,300)),allocatable :: e2(:,:)
contains
    subroutine lapack_init(n)
      integer :: n,state
      nvar = n
      allocate(ip1(nvar),ip2(nvar),ip_big(3*nvar),fjac(nvar,nvar),        &
         fjac1(nvar,nvar),fjac2(nvar,nvar),fjac3(nvar,nvar),e1(nvar,nvar),&
                    e2(nvar,nvar),e_big(3*nvar,3*nvar),STAT=state)
      if(state .ne. 0) then
        stop 'ALlocation error in lapack_init'
      end if
    end subroutine lapack_init

    subroutine lapack_free
      integer :: state
      deallocate(ip1,ip2,ip_big,fjac,fjac1,fjac2,fjac3,e1,e2,e_big,STAT=state)
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
          e1(i,j) = -fjac(i,j)
        end do
        e1(j,j) = e1(j,j) + hgamma
      end do
      call dgetrf( nvar, nvar, e1, nvar, ip1, ising )
    end subroutine lapack_decomp

    subroutine lapack_decomp_cmp(alpha,beta,ising)
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


    subroutine lapack_decomp_big(H,rkA,ising)
      integer :: ising, i, j
      double precision :: H,rkA(3,3)
! prepare left hand side
      do j=1,nvar
        do i=1,nvar
          e_big(i,j)         = -H*rkA(1,1)*fjac1(i,j)
          e_big(nvar+i,j)       = -H*rkA(2,1)*fjac1(i,j)
          e_big(2*nvar+i,j)     = -H*rkA(3,1)*fjac1(i,j)
          e_big(i,nvar+j)       = -H*rkA(1,2)*fjac2(i,j)
          e_big(nvar+i,nvar+j)     = -H*rkA(2,2)*fjac2(i,j)
          e_big(2*nvar+i,nvar+j)   = -H*rkA(3,2)*fjac2(i,j)
          e_big(i,2*nvar+j)     = -H*rkA(1,3)*fjac3(i,j)
          e_big(nvar+i,2*nvar+j)   = -H*rkA(2,3)*fjac3(i,j)
          e_big(2*nvar+i,2*nvar+j) = -H*rkA(3,3)*fjac3(i,j)
        end do
      end do
      do i=1,3*nvar
         e_big(i,i) = e_big(i,i) + 1
      end do
      call dgetrf( 3*nvar, 3*nvar, e_big, 3*nvar, ip_big, ising )
    end subroutine lapack_decomp_big

!~~~>Solving the system (HGamma-Jac)*X = RHS
    subroutine lapack_solve(trans, rhs)
      double precision ::rhs(nvar)
      logical :: trans
      integer :: ising
      if(trans) then 
        call dgetrs( 't', nvar, 1, e1, nvar, ip1, rhs, nvar, ising)
      else
        call dgetrs( 'n', nvar, 1, e1, nvar, ip1, rhs, nvar, ising)
      end if
    end subroutine lapack_solve


    subroutine lapack_solve_cmp(trans,b,bz)
      double precision :: b(nvar), bz(nvar)
      complex(kind=selected_real_kind(14,300)) :: rhs(nvar)
      integer :: ising,i
      logical :: trans
      do i=1,nvar
         rhs(i) = dcmplx(b(i),bz(i))
      end do
      if(trans) then
        call zgetrs( 't', nvar, 1, e2, nvar, ip2, rhs, nvar, ising)
      else
        call zgetrs( 'n', nvar, 1, e2, nvar, ip2, rhs, nvar, ising)
      end if

      do i = 1,nvar
        b(i) = dble(rhs(i))
        bz(i) = aimag(rhs(i))
      end do
    end subroutine lapack_solve_cmp

    subroutine lapack_solve_big(trans, rhs)
      double precision ::rhs(3*nvar)
      logical :: trans
      integer :: ising
      if(trans) then
        call dgetrs( 't', 3*nvar, 1, e_big, 3*nvar, ip_big, rhs, 3*nvar, ising)
      else
        call dgetrs( 'n', 3*nvar, 1, e_big, 3*nvar, ip_big, rhs, 3*nvar, ising)
      end if
    end subroutine lapack_solve_big

end module lapack
#endif

#ifdef SPARSE_UMF
! UMFPACK implementation 
module umf
    implicit none 
    save 
    integer :: nvar, nnz
!ax(nnz),axc(nnz),azc(nnz),ai(nnz),ap(nvar+1),ax_big(9*nnz),ai_big(9*nnz),ap_big(3*nvar+1)
    double precision,allocatable :: ax(:),axc(:),azc(:),ax_big(:)
    double precision,allocatable :: fjac1(:,:),fjac2(:,:),fjac3(:,:),fjac(:,:)
    integer,allocatable :: ai(:),ap(:),ai_big(:),ap_big(:)
    double precision ::  control(20), zcontrol(20),info(90)
    integer*8 :: symbolic,numeric,symbolic_cmp, numeric_cmp,symbolic_big,numeric_big
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

    subroutine umf_decomp_cmp(alpha,beta,ising)
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
    !last element in ap must be nnza+1
    ap(nvar+1) = nnz
    call umf4zdef(zcontrol)
    zcontrol(1) = 0
    call umf4zsym(nvar, nvar, ap, ai, axc, azc, symbolic_cmp, zcontrol, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4zsym: ', info(1)
      stop
    endif
    call umf4znum(ap, ai, axc, azc, symbolic_cmp, numeric_cmp, zcontrol, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4znum: ', info(1)
      stop
    endif
!    call umf4zfsym(symbolic_cmp)
    ising = 0
    end subroutine umf_decomp_cmp

    subroutine umf_decomp_big(H,rkA,ising)
      !convert the matrix e from coordinate format to column compressed format
      integer :: i,j,idx,found,ising
      integer :: col,row
      double precision :: current
      double precision :: H,rkA(3,3)
  !   free previously used memory
      if(symbolic_big .ne. 0) then
        call umf4fsym(symbolic_big)
        call umf4fnum(numeric_big)
      end if
      
      idx = 1
    ! c index starts from 0, while fortran one starts from 1
      do j = 1,3*nvar
        found = 0
        do i = 1,3*nvar
          ! convert indices
          row = mod(i-1,nvar)+1
          col = mod(j-1,nvar)+1
          if(j .le. nvar) then
              current = -H*rkA((i-row)/nvar+1,1)*fjac1(row,col)
          else if (j .le. 2*nvar) then
              current = -H*rkA((i-row)/nvar+1,2)*fjac2(row,col)  
          else
              current = -H*rkA((i-row)/nvar+1,3)*fjac3(row,col) 
          end if
          if ( current .ne. 0d0) then
            if( i .eq. j ) then
              ax_big(idx) = current + 1
            else 
              ax_big(idx) = current
            end if      
          ! write the row in the row indices vector ai
            ai_big(idx) = i - 1
            if(found .eq. 0) then
              found = idx
            end if
            idx = idx + 1
          end if
        end do
      !must also update the column pointers vector ap
      ap_big(j) = found - 1
    end do
    !last element in ap must be nnza+1
    ap_big(3*nvar+1) = 9*nnz
    call umf4def(control)
    control(1) = 0
    call umf4sym(3*nvar, 3*nvar, ap_big, ai_big, ax_big, symbolic_big,&
                          control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sym: ', info(1)
      stop
    endif
    call umf4num(ap_big, ai_big, ax_big, symbolic_big, numeric_big, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4num: ', info(1)
      stop
    endif
    ising = 0
    end subroutine umf_decomp_big

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
    subroutine umf_solve_cmp(trans, b, bz)
      double precision :: b(nvar), bz(nvar)
      logical :: trans
      integer :: sys
      double precision :: x(nvar),xz(nvar)
      
    ! solve ax=b, without iterative refinement
      if(trans) then
        sys = 1
      else
        sys = 0
      end if
      call umf4zsol(sys, x, xz, b, bz, numeric_cmp, zcontrol, info)
      if (info(1) .lt. 0) then
        print *, 'error occurred in umf4sol: ', info(1)
        stop
      endif
      b(:) = x(:)
      bz(:) = xz(:)
    end subroutine umf_solve_cmp

    subroutine umf_solve_big(trans, rhs)
    ! solve ax=b, without iterative refinement
    integer :: sys
    double precision :: x(3*nvar),rhs(3*nvar)
    logical :: trans
    if(trans) then
      sys = 1
    else
      sys = 0
    end if
    call umf4sol(sys, x, rhs, numeric_big, control, info)
    if (info(1) .lt. 0) then
      print *, 'error occurred in umf4sol: ', info(1)
      stop
    endif
    ! free the numeric factorization, to be continued
    rhs(:) = x(:)
    end subroutine umf_solve_big

    subroutine umf_init(n,nn)
      integer :: n,nn,state
      nvar = n
      nnz = nn
      symbolic = 0
      numeric = 0
      symbolic_cmp = 0
      numeric_cmp = 0
      symbolic_big = 0
      numeric_big = 0
      allocate(ax(nnz),axc(nnz),azc(nnz),ai(nnz),ap(nvar+1),ax_big(9*nnz),&
          ai_big(9*nnz),ap_big(3*nvar+1),fjac(nvar,nvar),fjac1(nvar,nvar),&
           fjac2(nvar,nvar),fjac3(nvar,nvar),STAT=state)
      if(state .ne. 0) then
        stop 'Allocation error in umf_init'
      end if
    end subroutine umf_init


!~~~>Clean objects created during the process
    subroutine umf_free
    integer :: state

    call umf4fsym(symbolic)
    call umf4fnum(numeric)
    
    call umf4fsym(symbolic_cmp)
    call umf4fnum(numeric_cmp)
    if(symbolic_big .ne. 0) then
      call umf4fsym(symbolic_big)
      call umf4fnum(numeric_big)
    end if
    deallocate(ax,axc,azc,ai,ap,ax_big,ai_big,ap_big,fjac,fjac1,fjac2,fjac3,&
                                   STAT=state)
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
    double precision,allocatable ::fjac(:,:),fjac1(:,:),fjac2(:,:),fjac3(:,:)
!   ax(nnz),b(nvar),ai(nnz), ap(nvar+1),axc(nnz),bc(nvar),ax_big(9*nnz),b_big(3*nvar),ai_big(9*nnz), ap_big(3*nvar+1)
    double precision,allocatable :: ax(:),b(:),ax_big(:),b_big(:)
    integer,allocatable :: ai(:), ap(:),ai_big(:),ap_big(:)
    complex(kind=selected_real_kind(14,300)),allocatable :: axc(:), bc(:)
!    integer :: aic(nnz), apc(nvar+1)
    integer*8 factors,factors_cmp,factors_big
contains
!~~~>Decomposition of the left hand side (HGamma-Jac)
    subroutine superlu_decomp(hgamma,ising)
      !convert the matrix e from coordinate format to column compressed format
      integer :: i,j,idx,found,ising
      integer :: iopt, ldb, nrhs, tr
      double precision :: hgamma 

!   free previously used memory
      nrhs = 1
      ldb = nvar
      tr = 0
      if(factors .ne. 0) then
        iopt = 3
        call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, ising )
      end if
     
      iopt = 1
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
      call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                          b, ldb, factors, ising )
  
      if(ising .ne. 0) then
        write(*,*) 'INFO from failed factorization = ', ising
      endif
    end subroutine superlu_decomp

    subroutine superlu_decomp_cmp(alpha,beta,ising)
      !convert the matrix e from coordinate format to column compressed format
      integer :: i,j,idx,found,ising
      integer :: iopt, ldb, nrhs, tr
      double precision :: alpha,beta
!   free previouslu used memory
      nrhs = 1
      ldb = nvar
      tr = 0
      if(factors_cmp .ne. 0) then
        iopt = 3
        call c_fortran_zgssv( tr, iopt, nvar, nnz, nrhs, axc, ai, ap,&
                        bc, ldb, factors_cmp, ising )
      end if

      iopt = 1
      idx = 1
      do j = 1,nvar
        found = 0
        do i = 1, nvar
          if ( fjac(i,j) .ne. 0d0) then
!            ax(idx) = e(i,j)
            if( i .eq. j ) then
              axc(idx) = -fjac(i,j) + cmplx(alpha,beta)
            else
            axc(idx) = dcmplx(-fjac(i,j),0d0)
            end if
            ! write the row in the row indices vector ai
!            aic(idx) = i -1
            if(found .eq. 0) then
              found = idx
            end if
            idx = idx + 1
          end if
        end do
        !must also update the column pointers vector ap
!        apc(j) = found -1
      end do
      !last element in ap must be nnza+1
!      apc(nvar+1) = nnz
  !   factorize the matrix. The factors are stored in *factors* handle.
      call c_fortran_zgssv( tr, iopt, nvar, nnz, nrhs, axc, ai, ap,&
                          bc, ldb, factors_cmp, ising )

      if(ising .ne. 0) then
        write(*,*) 'INFO from failed factorization = ', ising
      endif
    end subroutine superlu_decomp_cmp

    subroutine superlu_decomp_big(H,rkA,ising)
      !convert the matrix e from coordinate format to column compressed format
      integer :: i,j,idx,found,ising
      double precision :: H,rkA(3,3)
      integer :: col, row
      integer :: iopt, ldb, nrhs, tr
      double precision :: current

!   free previouslu used memory
      nrhs = 1
      ldb = 3*nvar
      tr = 0
      if(factors_big .ne. 0) then
        iopt = 3
        call c_fortran_dgssv(tr,iopt,3*nvar,9*nnz,nrhs,ax_big,ai_big,ap_big,&
                        b_big, ldb, factors_big, ising )
      end if
     
      iopt = 1
     
      idx = 1
    ! c index starts from 0, while fortran one starts from 1
      do j = 1,3*nvar
        found = 0
        do i = 1,3*nvar
          ! convert indices
          row = mod(i-1,nvar)+1
          col = mod(j-1,nvar)+1
          if(j .le. nvar) then
              current = -H*rkA((i-row)/nvar+1,1)*fjac1(row,col)
          else if (j .le. 2*nvar) then
              current = -H*rkA((i-row)/nvar+1,2)*fjac2(row,col)
          else
              current = -H*rkA((i-row)/nvar+1,3)*fjac3(row,col)
          end if
          if ( current .ne. 0d0) then
            if( i .eq. j ) then
              ax_big(idx) = current + 1
            else
              ax_big(idx) = current
            end if
          ! write the row in the row indices vector ai
            ai_big(idx) = i 
            if(found .eq. 0) then
              found = idx
            end if
            idx = idx + 1
          end if
        end do
      !must also update the column pointers vector ap
        ap_big(j) = found 
      end do
    !last element in ap must be nnza+1
      ap_big(3*nvar+1) = 9*nnz + 1
      call c_fortran_dgssv(tr,iopt,3*nvar,9*nnz,nrhs,ax_big,ai_big,ap_big,&
                          b_big, ldb, factors_big, ising )
      if(ising .ne. 0) then
        write(*,*) 'INFO from failed factorization (big) = ', ising
      end if
    end subroutine superlu_decomp_big

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

    subroutine superlu_solve_cmp(trans, rhs, rhsz)
    double precision :: rhs(nvar),rhsz(nvar)
    integer :: info,iopt, ldb, nrhs, tr
    logical :: trans
    integer :: i
    iopt = 2
    nrhs = 1
    ldb = nvar
    if(trans) then
      tr = 1
    else
      tr = 0
    end if
    do i = 1,nvar
      bc(i) = cmplx(rhs(i),rhsz(i))
    end do
    call c_fortran_zgssv(tr,iopt,nvar,nnz,nrhs,axc, ai, ap,&
                          bc, ldb, factors_cmp, info )
    do i = 1,nvar
      rhs(i) = dble(bc(i))
      rhsz(i) = aimag(bc(i))
    end do
    end subroutine superlu_solve_cmp

    subroutine superlu_solve_big(trans, rhs)
    double precision :: rhs(3*nvar)
    integer :: info, iopt, ldb, nrhs, tr
    logical :: trans
!   solve the system using the existing factors.
    iopt = 2
    nrhs = 1
    ldb = 3*nvar
    if(trans) then
      tr = 1
    else
      tr = 0
    end if
    call c_fortran_dgssv(tr,iopt,3*nvar,9*nnz,nrhs,ax_big,ai_big,ap_big,&
                          rhs, ldb, factors_big, info )
    if(info .ne. 0) then
      write(*,*) 'INFO from failed triangular solve big = ', info
    endif
    end subroutine superlu_solve_big

    subroutine superlu_init(n,nn)
      integer :: n,nn,state
      nvar = n
      nnz = nn
      allocate( ax(nnz),b(nvar),ai(nnz),ap(nvar+1),axc(nnz),bc(nvar),     &
               ax_big(9*nnz),b_big(3*nvar),ai_big(9*nnz),ap_big(3*nvar+1),&
       fjac(nvar,nvar),fjac1(nvar,nvar),fjac2(nvar,nvar),fjac3(nvar,nvar),&
                                 STAT=state)
      if(state .ne. 0) then
        stop 'Allocation error in superlu_init'
      end if
      factors = 0
      factors_cmp = 0
      factors_big = 0
    end subroutine superlu_init

    subroutine superlu_free
      integer :: info, iopt, ldb, nrhs, tr
      integer :: state
      iopt = 3
      ldb = nvar
      tr = 0
      nrhs = 1
      call c_fortran_dgssv( tr, iopt, nvar, nnz, nrhs, ax, ai, ap,&
                       b, ldb, factors, info )
      if(factors_cmp .ne. 0) then
        call c_fortran_zgssv( tr, iopt, nvar, nnz, nrhs, axc, ai, ap,&
                       bc, ldb, factors_cmp, info )
      end if
      if(factors_big .ne. 0) then
        ldb = 3*nvar
        call c_fortran_dgssv( tr,iopt,3*nvar,9*nnz,nrhs,ax_big,ai_big,ap_big,&
                       b_big,ldb,factors_big, info )
      end if

      deallocate(ax,b,ai,ap,axc,bc,ax_big,b_big,ai_big,ap_big,fjac,fjac1,&
                                                 fjac2,fjac3,STAT=state)
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
      subroutine lss_jac(t,y,jac)
      double precision :: t, y(nvar)
      external :: jac
      call jac(nvar,t,y,fjac)
      end subroutine
       
	  subroutine lss_jac1(t,y,jac)
	  double precision :: t, y(nvar)
	  external :: jac
	  call jac(nvar,t,y,fjac1)
	  end subroutine lss_jac1

	  subroutine lss_jac2(t,y,jac)
	  double precision :: t, y(nvar)
	  external :: jac
	  call jac(nvar,t,y,fjac2)
	  end subroutine lss_jac2

	  subroutine lss_jac3(t,y,jac)
	  double precision :: t, y(nvar)
	  external :: jac
	  call jac(nvar,t,y,fjac3)
	  end subroutine lss_jac3

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

      subroutine lss_decomp_cmp(alpha, beta, ising)
      integer :: ising
      double precision :: alpha,beta
#ifdef SPARSE_UMF
      call umf_decomp_cmp(alpha, beta,ising)
#endif
#ifdef FULL_ALGEBRA
      call lapack_decomp_cmp(alpha,beta,ising)
#endif
#ifdef SPARSE_LU
      call superlu_decomp_cmp(alpha,beta,ising)
#endif
      end subroutine lss_decomp_cmp


      subroutine lss_decomp_big(h, rkA, ising)
      integer :: ising
      double precision :: h,rkA(3,3)
#ifdef SPARSE_UMF
      call umf_decomp_big(h,rkA,ising)
#endif
#ifdef FULL_ALGEBRA
      call lapack_decomp_big(h,rkA,ising)
#endif
#ifdef SPARSE_LU
      call superlu_decomp_big(h,rkA,ising)
#endif 
      end subroutine lss_decomp_big

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

      subroutine lss_solve_big(trans, rhs)
      double precision :: rhs(3*nvar)
      logical :: trans
#ifdef SPARSE_UMF  
      call umf_solve_big(trans, rhs)
#endif
#ifdef FULL_ALGEBRA
      call lapack_solve_big(trans, rhs)
#endif
#ifdef SPARSE_LU
      call superlu_solve_big(trans, rhs)
#endif
      end subroutine lss_solve_big

      subroutine lss_solve_cmp(trans,b,bz)
      double precision  :: b(nvar), bz(nvar)
      logical :: trans
#ifdef SPARSE_UMF
      call umf_solve_cmp(trans,b,bz)
#endif
#ifdef FULL_ALGEBRA
      call lapack_solve_cmp(trans,b,bz)
#endif
#ifdef SPARSE_LU
      call superlu_solve_cmp(trans,b,bz)
#endif
!      print *, rhs  
      end subroutine lss_solve_cmp

      subroutine lss_init(n,nn)     
      integer :: n , nn
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

!~~~> multiplication of the jacobian matrix transpose and a vector
     subroutine lss_mul_jactr1(z,g)
      double precision :: z(nvar), g(nvar)
          z = matmul(transpose(fjac1),g)
     end subroutine lss_mul_jactr1

     subroutine lss_mul_jactr2(z,g)
      double precision :: z(nvar), g(nvar)
          z = matmul(transpose(fjac2),g)
     end subroutine lss_mul_jactr2

     subroutine lss_mul_jactr3(z,g)
      double precision :: z(nvar), g(nvar)
          z = matmul(transpose(fjac3),g)
     end subroutine lss_mul_jactr3

end module ls_solver
