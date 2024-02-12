! Solve A x = b
! test sparse solver by using known b vector, e.g.
! multiply A with x vector of 2d0, which gives b. Use
! this b to solve for x again
subroutine solve_poisson

    use config, only: fstatus, fdebug
    use kernel, only: K_rho_phi_llp, K_rho_B_llp
    use sparse_mod, only: sp2fullComplex, sparse_solveComplex_b1, column_pointer2full, sparse_solve_suitesparseComplex_b1, &
                          sparse_solve_method
    use config, only: output_path
    use constants, only: sol, p_mass, e_charge, e_mass, pi
    use grid, only: npoib, rb
    use back_quants, only: vTi, vTe
    use plasma_parameter, only: Zi, Ai
    use setup, only: btor, cut_off_fac

    implicit none

    double complex, dimension(:), allocatable :: A_nz ! non-zero elements of A matrix
    double complex, dimension(:,:), allocatable :: A_mat ! A matrix
    double complex, dimension(:,:), allocatable :: A_sparse_check ! A matrix reconfigured from sparse matrix
    double complex, dimension(:), allocatable :: b_vec ! b vector and x vector
    integer, dimension(:), allocatable :: irow, pcol, icol
    integer :: nz_out, nrow, ncol
    integer :: i,j
    logical :: ex
    double precision :: rho_L
    integer, dimension(2) :: max_ind
    integer :: set_zero
    integer :: iopt

    if (fstatus == 1) write(*,*) 'Status: solve poisson equation'

    ! create A matrix (Laplacian + 4 pi K^phi) in l space
    allocate(A_mat(npoib, npoib))

    ! create Laplacian:
    A_mat = cmplx(0.0d0, 0.0d0)
    A_mat(1,1) = -1d0 / (rb(1)) - 1d0 / (rb(2) + rb(1))
    A_mat(npoib,npoib) = - 1d0 / (rb(npoib) - rb(npoib-1))
    A_mat(1,2) = 1d0 / (rb(2) - rb(1))
    A_mat(2,1) = 1d0 / (rb(2) - rb(1))

    do i = 2, npoib-1
        A_mat(i,i) = - 1d0 / (rb(i) - rb(i-1)) - 1d0 / (rb(i+1) - rb(i))
        A_mat(i, i+1) = 1d0 / (rb(i+1) - rb(i))
        A_mat(i, i-1) = 1d0 / (rb(i) - rb(i-1))
    end do

    if (fdebug == 1) then
        write(*,*) 'Debug: writing Laplacian A matrix before sparse'
        open(unit=77, file=trim(output_path)//'kernel/laplacian_re.dat')
        do i = 1,npoib
            do j = 1,npoib
                write(77,*) real(A_mat(i,j))
            end do
        end do
        close(77)
    end if

    ! make kernel matrix sparse with cut_off_fac * ion Larmor radius cut off
    ! max value of the ion Larmor radius
    !max_ind = maxloc(vTi)
    !rho_L = vTi(max_ind(1), max_ind(2)) * Ai(max_ind(1)) * p_mass * sol &
    !        / (e_charge * Zi(max_ind(1)) * abs(btor))

    !write(*,*) "max(rho_Li) = ", rho_L

    !do i = 1, npoib
    !    do j = 1, npoib
    !        if (abs(rb(i) - rb(j)) >= cut_off_fac * rho_L) then
    !            !write(*,*) 'set kernel zero'
    !            K_rho_phi_llp(i,j) = cmplx(0.0d0, 0.0d0)
    !            !K_rho_B_llp(i,j) = cmplx(0.0d0, 0.0d0)
    !            set_zero = set_zero + 1
    !        end if
    !        A_mat(i,j) = A_mat(i,j) + cmplx(4.0d0 * pi, 0.0d0) * K_rho_phi_llp(i,j)
    !        !write(*,*) (1.0d0, 0.0d0) * 4.0d0 * pi * K_rho_phi_llp(i,j)
    !    end do
    !end do
    !if (fdebug == 1) write(*,*) 'set_zero = ', set_zero

    A_mat = A_mat + 4d0 * pi * K_rho_phi_llp

    if (fdebug == 1) then
        write(*,*) 'Debug: writing A matrix before sparse'
        open(unit=77, file=trim(output_path)//'kernel/A_mat_before_re.dat')
        open(unit=78, file=trim(output_path)//'kernel/A_mat_before_im.dat')
        open(unit=79, file=trim(output_path)//'kernel/K_rho_phi_llp_sp_re.dat')
        open(unit=80, file=trim(output_path)//'kernel/K_rho_phi_llp_sp_im.dat')
        do i = 1,npoib
            do j = 1,npoib
                write(77,*) real(A_mat(i,j))
                !write(*,*) real(A_mat(i,j))
                write(78,*) dimag(A_mat(i,j))
                write(79,*) real(K_rho_phi_llp(i,j))
                write(80,*) dimag(K_rho_phi_llp(i,j))
            end do
        end do
        close(77)
        close(78)
        close(79)
        close(80)
    end if

    ! make A matrix sparse
    call dense_to_sparse(A_mat, irow, pcol, A_nz, nrow, ncol, nz_out)

    write(*,*) 'nz_out = ', nz_out, '; nrow = ', nrow, '; ncol = ', ncol

    !allocate(A_sparse_check(nrow,ncol))

    if (fdebug == 3) then
        write(*,*) 'Debug: writing A matrix after sparse'
        call sp2fullComplex(irow, pcol, A_nz, nrow, ncol, A_sparse_check)
        open(unit=77, file=trim(output_path)//'kernel/A_sparse_check_re.dat')
        open(unit=78, file=trim(output_path)//'kernel/A_sparse_check_im.dat')
        do i = 1,nrow
            do j = 1,ncol
                write(77,*) real(A_sparse_check(i,j))
                write(78,*) dimag(A_sparse_check(i,j))
            end do
        end do
        close(77)
        close(78)
        deallocate(A_sparse_check)
    end if

    ! multiply B kernel with B^r vector

    write(*,*) "create b vec"

    if (.not. allocated(b_vec)) allocate(b_vec(npoib))
    b_vec = cmplx(1.0d0, 0.0d0)
    b_vec = - 4d0 * pi * matmul(K_rho_B_llp, b_vec)
    !b_vec = matmul(A_mat, b_vec)

    write(*,*) "after b vec"
    ! solve matrix vector problem
    !call column_pointer2full(pcol, icol)
    !write(*,*) 'A_nz = ', A_nz
    iopt = 0
    sparse_solve_method = 1 ! this works, don't know why. Default value of 3 does not work. I.e. need 
    !to use superlu instead of suitesparse
    call sparse_solveComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, iopt)
    !call sparse_solve_suitesparseComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, 0)
    !write(*,*) b_vec

    inquire(file=trim(output_path)//'fields', exist=ex)
    if (.not. ex) then
        call system('mkdir -p '//trim(output_path)//'fields')
    end if
    open(unit = 77, file=trim(output_path)//'fields/phi_re.dat')
    open(unit = 78, file=trim(output_path)//'fields/phi_im.dat')
    do i = 1,npoib
        write(77,*) rb(i), real(b_vec(i))
        write(78,*) rb(i), dimag(b_vec(i))
    end do
    close(77)
    close(78)

    if (fdebug == 3) then
        write(*,*) 'Debug : write A matrix '
        open(unit = 80, file=trim(output_path)//'fields/A_mat_re.dat')
        open(unit = 81, file=trim(output_path)//'fields/A_mat_im.dat')
        do i = 1,npoib
            do j = 1,npoib
                write(80,*) real(A_mat(i,j))
                write(81,*) dimag(A_mat(i,j))
            end do
        end do
        close(80)
        close(81)
    end if

    contains

    ! transform kernel matrix in l space to sparse matrix
    ! ignores elements that are further apart than 5 times the
    ! (ion) Larmor radius
    subroutine dense_to_sparse(A, irow, pcol, A_nz, nrow, ncol, nz_out)

        use sparse_mod, only: column_full2pointer
        
        implicit none

        double complex, dimension(:,:), intent(in) :: A
        double complex, dimension(:), allocatable, intent(inout) :: A_nz
        integer, dimension(:), allocatable, intent(inout) :: irow, pcol
        integer, optional, intent(out) :: nz_out
        integer, intent(out) :: nrow, ncol

        integer, dimension(:), allocatable :: icol
        integer :: nz, nc, nr, n
        
        nrow = size(A,1)
        ncol = size(A,2)
        nz = 0

        do nc = 1, ncol
            do nr = 1, nrow
                if (.not.(A(nr,nc) == 0.0d0)) then
                    nz = nz + 1
                end if
            end do
        end do

        if (allocated(irow)) deallocate(irow)
        allocate(irow(nz))
        if (allocated(icol)) deallocate(icol)
        allocate(icol(nz))
        if (allocated(A_nz)) deallocate(A_nz)
        allocate(A_nz(nz))

        A_nz = 0.0d0

        n = 0
        do nc = 1, ncol
            do nr = 1, nrow
                if (.not.(A(nr, nc) == 0.0d0)) then
                    n = n + 1
                    irow(n) = nr
                    icol(n) = nc
                    A_nz(n) = A_nz(n) + A(nr, nc)
                end if
            end do
        end do

        call column_full2pointer(icol, pcol)

        if (present(nz_out)) nz_out = nz
        if (allocated(icol)) deallocate(icol)

    end subroutine

end subroutine


