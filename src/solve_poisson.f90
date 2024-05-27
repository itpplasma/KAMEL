module poisson_solver

    contains
    
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
        use setup, only: btor, cut_off_fac, type_br_field

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

        call prepare_Laplace_matrix

        call check_kernels_for_nans

        A_mat = (A_mat + 4d0 * pi * K_rho_phi_llp) !/ npoib

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

        call create_rhs_vector(type_br_field, b_vec)

        if (type_br_field == 1 .or. type_br_field == 3) then
            ! multiply with K_rho_B_llp if not point charge case
            write(*,*) "multiply with K_rho_B_llp"
            b_vec = - 4d0 * pi * matmul(K_rho_B_llp, b_vec)
        end if

        inquire(file=trim(output_path)//'fields', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'fields')
        end if
        open(unit = 79, file=trim(output_path)//'fields/Kbr_re.dat')
        open(unit = 80, file=trim(output_path)//'fields/Kbr_im.dat')
        do i = 1,npoib
            write(79,*) rb(i), real(b_vec(i))
            write(80,*) rb(i), dimag(b_vec(i))
        end do
        close(79)
        close(80)

        write(*,*) "after b vec"

        iopt = 0
        sparse_solve_method = 1 ! this works, don't know why. Default value of 3 does not work. I.e. need 
        !to use superlu instead of suitesparse
        call sparse_solveComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, iopt)
        !call sparse_solve_suitesparseComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, 0)
        !write(*,*) b_vec

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

        subroutine prepare_Laplace_matrix

            use grid, only: xl
            implicit none
            integer :: i,j
            double precision, allocatable :: dr(:)

            ! create A matrix (Laplacian + 4 pi K^phi) in l space
            allocate(A_mat(npoib, npoib))
            allocate(dr(npoib))

            call initialize_grid_spacing(dr, xl)
            ! create Laplacian:
            A_mat = cmplx(0.0d0, 0.0d0)

            do i = 2, npoib-1
                A_mat(i,i) = (- 1d0 / (rb(i) - rb(i-1)) - 1d0 / (rb(i+1) - rb(i)) ) !* 2 &
                !* (rb(i+1) - rb(i-1))
                A_mat(i, i+1) = 1d0 / (rb(i+1) - rb(i)) !* 2 * (rb(i+1) - rb(i-1))
                A_mat(i, i-1) = 1d0 / (rb(i) - rb(i-1)) !* 2 * (rb(i+1) - rb(i-1))
            end do

            ! boundary conditions:
            A_mat(1,1) = -1d0 / (rb(1)) !* 2 * (rb(2) - rb(1)) !- 1d0 / (rb(2) - rb(1))
            A_mat(npoib,npoib) = - 1d0 / (rb(npoib) - rb(npoib-1)) !* 2 * (rb(npoib) - rb(npoib-1))
            A_mat(1,2) = 0.0d0 !1d0 / (rb(2) - rb(1))
            A_mat(npoib, npoib-1) = 0.0d0


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

        end subroutine

        subroutine initialize_grid_spacing(dr, r_in)

            implicit none
            double precision, allocatable, intent(in) :: r_in(:)
            double precision, allocatable, intent(out) :: dr(:)

            integer :: i

            allocate(dr(size(r_in)))

            do i = 1, size(r_in)
                dr(i) = r_in(i+1) - r_in(i)
            end do
            dr(size(r_in)) = dr(size(r_in)-1)

        end subroutine

        ! transform kernel matrix in l space to sparse matrix
        ! ignores elements that are further apart than 5 times the
        ! (ion) Larmor radius
        subroutine check_kernels_for_nans

            ! check if there are nan's in the kernels
            if (any(isnan(real(K_rho_phi_llp)))) then
                print *, "K_rho_phi_llp contains NaN values."
            else
                print *, "K_rho_phi_llp does not contain NaN values."
            end if

            if (any(isnan(real(K_rho_B_llp)))) then
                print *, "K_rho_B_llp contains NaN values."
            else
                print *, "K_rho_B_llp does not contain NaN values."
            end if

            if (any(isnan(real(A_mat)))) then
                print *, "A_mat contains NaN values."
            else
                print *, "A_mat does not contain NaN values."
            end if

        end subroutine


        subroutine create_rhs_vector(type, rhs_vec)

            use resonances_mod, only: index_res

            implicit none

            integer, intent(in) :: type
            double complex, allocatable, intent(out) :: rhs_vec(:)
            integer :: i

            allocate(rhs_vec(npoib))
            rhs_vec = cmplx(0.0d0, 0.0d0)

            if (type ==1) then
                ! constant Br field
                rhs_vec = cmplx(1.0d0, 0.0d0)
            elseif(type == 2) then
                ! point charge like Br field
                rhs_vec(size(rhs_vec)/2) = cmplx(-4.0d0 * pi, 0.0d0)  * e_charge 
                rhs_vec(1) = cmplx(1.0d-13, 0.0d0)
                rhs_vec(npoib) = cmplx(1.0d-13, 0.0d0)
            elseif(type ==3) then
                do i = 5/6 *size(b_vec), size(b_vec)
                    rhs_vec(i) = (i - size(rhs_vec)/2) * 0.02d0 * cmplx(1.0d0, 0.0d0) - 0.2d0
                end do 
            elseif(type == 4) then
                ! put point charge at resonant surface
                rhs_vec(index_res) = cmplx(-4.0d0 * pi, 0.0d0) * e_charge 
            elseif(type == 5) then
                ! read from file
                print *, "Reading B vector from file not implemented yet"
                stop
            elseif(type == 6) then
                rhs_vec = cmplx(-4.0d0 * pi, 0.0d0) * e_charge * exp(- (rb - maxval(rb)/2)**2 / 0.1d0**2) &
                        * sqrt(pi / 0.1d0**2)
            end if

            rhs_vec = rhs_vec !/ npoib

            inquire(file=trim(output_path)//'fields', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'fields')
            end if
            open(unit = 79, file=trim(output_path)//'fields/br_re.dat')
            open(unit = 80, file=trim(output_path)//'fields/br_im.dat')
            do i = 1,npoib
                write(79,*) rb(i), real(b_vec(i))
                write(80,*) rb(i), dimag(b_vec(i))
            end do
            close(79)
            close(80)


        end subroutine

    end subroutine

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



end module