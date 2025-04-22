module poisson_solver

    contains
    
    ! Solve A x = b
    subroutine solve_poisson(K_rho_phi, K_rho_B, phi_sol)

        use config, only: fstatus, fdebug
        use sparse_mod, only: sp2fullComplex, sparse_solveComplex_b1, column_pointer2full, sparse_solve_suitesparseComplex_b1, &
                            sparse_solve_method
        use config, only: output_path
        use constants, only: sol, p_mass, e_charge, e_mass, pi
        use grid, only: xl_grid
        use back_quants, only: vTi, vTe
        use plasma_parameter, only: Zi, Ai
        use setup, only: btor, cut_off_fac, type_br_field
        use KIM_kinds, only: dp

        implicit none

        complex(dp), intent(in) :: K_rho_phi(:,:)
        complex(dp), intent(in) :: K_rho_B(:,:)
        complex(dp), dimension(:), allocatable :: A_nz ! non-zero elements of A matrix
        complex(dp), dimension(:,:), allocatable :: A_mat ! A matrix
        complex(dp), dimension(:), allocatable :: b_vec ! b vector and x vector
        complex(dp), dimension(:), allocatable :: phi_sol
        integer, dimension(:), allocatable :: irow, pcol, icol
        integer :: nz_out, nrow, ncol
        integer :: i,j
        logical :: exists
        real(dp) :: rho_L
        integer :: sparse_solver_option


        if (fstatus == 1) write(*,*) 'Status: solve poisson equation'

        allocate(A_mat(xl_grid%npts_b, xl_grid%npts_b))
        call prepare_Laplace_matrix(A_mat)

        call check_kernels_for_nans(K_rho_phi)
        call check_kernels_for_nans(K_rho_B)

        A_mat = (A_mat + 4.0d0 * pi * K_rho_phi) 

        if (fdebug == 3) then
            call write_A_matrix_sparse_check_to_file
        end if

        call dense_to_sparse(A_mat, irow, pcol, A_nz, nrow, ncol, nz_out)
        
        call create_rhs_vector(type_br_field, K_rho_B, b_vec)
        
        sparse_solver_option = 0
        sparse_solve_method = 1 ! this works, don't know why. Default value of 3 does not work. I.e. need 
        !to use superlu instead of suitesparse
        call sparse_solveComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, sparse_solver_option)
        !call sparse_solve_suitesparseComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, 0)
        phi_sol = b_vec
        call write_phi_to_file

        if (fdebug == 3) then
            call write_A_matrix_to_file
        end if

        contains

        subroutine write_A_matrix_to_file

            implicit none

            write(*,*) 'Debug : write A matrix '
            open(unit = 80, file=trim(output_path)//'fields/A_mat_re.dat')
            open(unit = 81, file=trim(output_path)//'fields/A_mat_im.dat')
            do i = 1, xl_grid%npts_b
                do j = 1, xl_grid%npts_b
                    write(80,*) real(A_mat(i,j))
                    write(81,*) dimag(A_mat(i,j))
                end do
            end do
            close(80)
            close(81)

        end subroutine

        subroutine write_phi_to_file

            implicit none

            open(unit = 77, file=trim(output_path)//'fields/phi_re.dat')
            open(unit = 78, file=trim(output_path)//'fields/phi_im.dat')
            do i = 1, xl_grid%npts_b
                write(77,*) xl_grid%xb(i), real(b_vec(i))
                write(78,*) xl_grid%xb(i), dimag(b_vec(i))
            end do
            close(77)
            close(78)

        end subroutine


        subroutine write_K_times_b_to_file

            implicit none

            inquire(file=trim(output_path)//'fields', exist=exists)
            if (.not. exists) then
                call system('mkdir -p '//trim(output_path)//'fields')
            end if
            open(unit = 79, file=trim(output_path)//'fields/Kbr_re.dat')
            open(unit = 80, file=trim(output_path)//'fields/Kbr_im.dat')
            do i = 1,xl_grid%npts_b
                write(79,*) xl_grid%xb(i), real(b_vec(i))
                write(80,*) xl_grid%xb(i), dimag(b_vec(i))
            end do
            close(79)
            close(80)

        end subroutine

        subroutine write_A_matrix_sparse_check_to_file

            implicit none
            double complex, dimension(:,:), allocatable :: A_sparse_check ! A matrix reconfigured from sparse matrix

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

        subroutine create_rhs_vector(type, K_rho_B, rhs_vec)

            use resonances_mod, only: index_rg_res
            use functions, only: varphi_l
            use grid, only: xl_grid
            use plotting, only: write_profile
            use KIM_kinds, only: dp

            implicit none

            integer, intent(in) :: type
            complex(dp), allocatable, intent(out) :: rhs_vec(:)
            complex(dp), intent(in) :: K_rho_B(:,:)

            integer :: i, idx
            real(dp) :: x0

            x0 = 35.0d0

            allocate(rhs_vec(xl_grid%npts_b))
            rhs_vec = cmplx(0.0d0, 0.0d0, dp)

            if (type ==1) then ! constant br
                rhs_vec = cmplx(1.0d0, 0.0d0, dp) * e_charge
            elseif(type == 2) then ! point charge like Br field

                idx = minloc(abs(xl_grid%xb - x0), dim=1)
                rhs_vec(idx) = cmplx(1.0d0, 0.0d0, dp) * e_charge ! * exp(- (rg_grid%xb(idx) - 35.0d0)**2 / 0.1d0**2)

            elseif(type ==3) then ! linear increase from the center of the plasma
                do i = int(5.0d0/6.0d0 *size(b_vec)), size(b_vec)
                    rhs_vec(i) = (i - size(rhs_vec)/2) * 0.02d0 * cmplx(1.0d0, 0.0d0, dp) - 0.2d0
                end do 
            elseif(type == 4) then ! point charge at resonant surface
                rhs_vec(index_rg_res) = cmplx(-4.0d0 * pi, 0.0d0, dp) * e_charge 
            elseif(type == 5) then ! read br from file (not implemented yet)
                print *, "Reading B vector from file not implemented yet"
                stop
            elseif(type == 6) then ! gaussian distribution
                rhs_vec = cmplx(1.0d0, 0.0d0, dp) * e_charge * exp(- (xl_grid%xb - x0)**2 / 0.1d0**2) &
                        * sqrt(pi / 0.1d0**2)
            elseif(type==7) then
                rhs_vec = cmplx(0.0d0, 0.0d0, dp)
                do i = 2, xl_grid%npts_b-1
                    rhs_vec(i) = e_charge * varphi_l(x0, xl_grid%xb(i-1), xl_grid%xb(i), xl_grid%xb(i+1))
                end do
            elseif(type==11) then
                rhs_vec = 1.0d0
                rhs_vec = matmul(K_rho_B, rhs_vec)
            elseif(type==12) then
                !rhs_vec = 1.0d0
                rhs_vec = cmplx(1.0d0, 0.0d0, dp) * exp(- (xl_grid%xb - x0)**2 / 1.0d0**2) &
                        * sqrt(pi / 0.1d0**2)
                rhs_vec = matmul(K_rho_B, rhs_vec)
            end if

            rhs_vec = - 4d0 * pi * rhs_vec

            call write_profile(xl_grid%xb, real(rhs_vec), xl_grid%npts_b, trim(output_path)//'fields/rhs_vec_re.dat')
            call write_profile(xl_grid%xb, dimag(rhs_vec), xl_grid%npts_b, trim(output_path)//'fields/rhs_vec_im.dat')

        end subroutine

    end subroutine

    subroutine check_kernels_for_nans(kernel)

            use KIM_kinds, only: dp

            implicit none

            complex(dp), intent(in) :: kernel(:,:)

            ! check if there are nan's in the kernels
            if (any(isnan(real(kernel)))) then
                print *, "kernel contains NaN values."
            else
                print *, "kernel does not contain NaN values."
            end if

        end subroutine

    subroutine prepare_Laplace_matrix(A_mat)

        use grid, only: xl_grid, rg_grid
        use KIM_kinds, only: dp
        use config, only: fdebug, output_path
        use constants, only: pi
        use plotting, only: write_matrix

        implicit none

        complex(dp), intent(inout) :: A_mat(:,:)
        real(dp) :: h
        integer :: i,j

        ! create Laplacian:
        A_mat = cmplx(0.0d0, 0.0d0, dp)

        !do i = 2, rg_grid%npts_b-1
            !A_mat(i,i) = (- 1d0 / (rg_grid%xb(i) - rg_grid%xb(i-1)) - 1d0 / (rg_grid%xb(i+1) - rg_grid%xb(i)) ) !* 2 &
            !!* (rg_grid%xb(i+1) - rg_grid%xb(i-1))
            !A_mat(i, i+1) = 1d0 / (rg_grid%xb(i+1) - rg_grid%xb(i)) !* 2 * (rg_grid%xb(i+1) - rg_grid%xb(i-1))
            !A_mat(i, i-1) = 1d0 / (rg_grid%xb(i) - rg_grid%xb(i-1)) !* 2 * (rg_grid%xb(i+1) - rg_grid%xb(i-1))
        !end do

        !! boundary conditions:
        !A_mat(1,1) = - 2d0 / (rg_grid%xb(2) - rg_grid%xb(1)) !* 2 * (rg_grid%xb(2) - rg_grid%xb(1)) !- 1d0 / (rg_grid%xb(2) - rg_grid%xb(1))
        !A_mat(rg_grid%npts_b,rg_grid%npts_b) = - 2d0 / (rg_grid%max_val - rg_grid%xb(rg_grid%npts_b-1)) !* 2 * (rg_grid%max_val - rg_grid%xl(rg_grid%npts_b-1))
        !A_mat(1,2) = 1d0 / (rg_grid%xb(2) - rg_grid%xb(1))
        !A_mat(rg_grid%npts_b, rg_grid%npts_b-1) = 1d0 / (rg_grid%max_val - rg_grid%xb(rg_grid%npts_b-1))

        A_mat = 0.0_dp

        do i = 1, xl_grid%npts_b - 1
            h = xl_grid%xb(i+1) - xl_grid%xb(i)
            A_mat(i,   i  ) = A_mat(i,   i  ) - 1.0_dp / h
            A_mat(i,   i+1) = A_mat(i,   i+1) + 1.0_dp / h
            A_mat(i+1, i  ) = A_mat(i+1, i  ) + 1.0_dp / h
            A_mat(i+1, i+1) = A_mat(i+1, i+1) - 1.0_dp / h
        end do

        !A_mat(1,:) = 0.0_dp
        !A_mat(:,1) = 0.0_dp
        !A_mat(1,1) = A_mat(1,1) - 1d0 / (xl_grid%xb(2) - xl_grid%xb(1))! 1.0_dp
        !A_mat(rg_grid%npts_b,rg_grid%npts_b) = A_mat(rg_grid%npts_b,rg_grid%npts_b) - 1d0 / (rg_grid%max_val - rg_grid%xb(rg_grid%npts_b-1)) !* 2 * (rg_grid%max_val - rg_grid%xl(rg_grid%npts_b-1))


        !if (fdebug == 1) then
            !write(*,*) 'Debug: writing Laplacian A matrix before sparse'
            !open(unit=77, file=trim(output_path)//'kernel/laplacian_re.dat')
            !do i = 1,rg_grid%npts_b
                !do j = 1,rg_grid%npts_b
                    !write(77,*) real(A_mat(i,j))
                !end do
            !end do
            !close(77)
        !end if

        call write_matrix(trim(output_path)//'kernel/laplacian_re.dat', real(A_mat), rg_grid%npts_b, rg_grid%npts_b)

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
