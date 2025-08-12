module poisson_solver

    contains
    
    ! Solve A x = b
    subroutine solve_poisson(K_rho_phi, K_rho_B, phi_sol)

        use config, only: fstatus, fdebug
        use sparse_mod, only: sp2fullComplex, sparse_solveComplex_b1, column_pointer2full, sparse_solve_suitesparseComplex_b1, &
                            sparse_solve_method
        use config, only: output_path
        use constants, only: pi, e_charge
        use grid, only: xl_grid
        use setup, only: type_br_field
        use KIM_kinds, only: dp

        implicit none

        complex(dp), intent(in) :: K_rho_phi(:,:)
        complex(dp), intent(in) :: K_rho_B(:,:)
        complex(dp), dimension(:), allocatable, intent(out) :: phi_sol
        complex(dp), dimension(:), allocatable :: A_nz ! non-zero elements of A matrix
        complex(dp), dimension(:,:), allocatable :: A_mat ! A matrix
        complex(dp), dimension(:), allocatable :: b_vec ! b vector and x vector
        integer, dimension(:), allocatable :: irow, pcol
        integer :: nz_out, nrow, ncol
        integer :: i,j
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
        sparse_solve_method = 1 
        call sparse_solveComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, sparse_solver_option)

        phi_sol = b_vec

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

        subroutine write_A_matrix_sparse_check_to_file

            implicit none

            complex(dp), dimension(:,:), allocatable :: A_sparse_check ! A matrix reconfigured from sparse matrix

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


        ! transform kernel matrix in l space to sparse matrix
        ! ignores elements that are further apart than 5 times the
        ! (ion) Larmor radius

        subroutine create_rhs_vector(type, K_rho_B, rhs_vec)

            use resonances_mod, only: index_rg_res, r_res
            use functions, only: varphi_l
            use grid, only: xl_grid
            use IO_collection, only: write_profile, write_complex_profile, plot_profile
            use KIM_kinds, only: dp
            use fields, only: EBdat, set_Br_field
            use config, only: output_path

            implicit none

            integer, intent(in) :: type
            complex(dp), allocatable, intent(out) :: rhs_vec(:)
            complex(dp), intent(in) :: K_rho_B(:,:)

            integer :: i, idx
            real(dp) :: x0

            x0 = 58.5d0

            allocate(rhs_vec(xl_grid%npts_b))
            rhs_vec = cmplx(0.0d0, 0.0d0, dp)

            if (type ==1) then ! constant br
                rhs_vec = cmplx(1.0d0, 0.0d0, dp) * e_charge
            elseif(type == 2) then ! point charge like Br field

                idx = minloc(abs(xl_grid%xb - r_res), dim=1)
                rhs_vec(idx) = cmplx(1.0d0, 0.0d0, dp) * e_charge ! * exp(- (rg_grid%xb(idx) - 35.0d0)**2 / 0.1d0**2)

            elseif(type ==3) then ! linear increase from the center of the plasma
                do i = int(5.0d0/6.0d0 *size(b_vec)), size(b_vec)
                    rhs_vec(i) = (i - size(rhs_vec)/2) * 0.02d0 * cmplx(1.0d0, 0.0d0, dp) - 0.2d0
                end do 
            elseif(type == 4) then ! point charge at resonant surface
                rhs_vec(index_rg_res) = cmplx(1.0, 0.0d0, dp) * e_charge 
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
            ! type > 10 uses Br kernel to determine right hand side of equation
            elseif(type==11) then 
                call set_Br_field(EBdat, 1) ! Br field from input
                rhs_vec = matmul(K_rho_B, EBdat%Br)
            elseif(type==12) then 
                call set_Br_field(EBdat, 0) ! constant Br field
                rhs_vec = matmul(K_rho_B, EBdat%Br)
            elseif(type==13) then ! linear increase including kernel
                rhs_vec = 0.0d0
                do i = 1, size(b_vec)
                    if (xl_grid%xb(i) > 50.0d0) then
                        rhs_vec(i) = (xl_grid%xb(i) - 50.0d0) * 0.1d0
                    end if
                end do 
                call write_complex_profile(xl_grid%xb, rhs_vec, xl_grid%npts_b, trim(output_path)//'fields/br_pert.dat')
                rhs_vec = matmul(K_rho_B, rhs_vec)
            end if

            rhs_vec = - 4d0 * pi * rhs_vec
            
            ! Enforce Dirichlet BC at right boundary: Phi_n = 0
            ! The RHS at the last point should be 0 for consistency
            rhs_vec(xl_grid%npts_b) = 0.0_dp

            call write_complex_profile(xl_grid%xb, rhs_vec, xl_grid%npts_b, trim(output_path)//'fields/rhs_vec.dat')

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

        use grid, only: xl_grid
        use KIM_kinds, only: dp
        use config, only: output_path
        use IO_collection, only: write_matrix

        implicit none

        complex(dp), intent(inout) :: A_mat(:,:)
        real(dp) :: h, hL, hR
        integer :: i
        integer :: n

        n = xl_grid%npts_b
        
        ! create Laplacian:
        A_mat = cmplx(0.0d0, 0.0d0, dp)
        A_mat = 0.0d0

        do i = 2, n-1
            hL = xl_grid%xb(i)   - xl_grid%xb(i-1)  ! left element size
            hR = xl_grid%xb(i+1) - xl_grid%xb(i)    ! right element size

            A_mat(i,i-1) = A_mat(i,i-1) - 1.0d0/hL
            A_mat(i,i)   = A_mat(i,i)   + 1.0d0/hL + 1.0d0/hR
            A_mat(i,i+1) = A_mat(i,i+1) - 1.0d0/hR
        end do
        ! Interior points (standard finite difference)
        ! do i = 2, n - 2
        !     h = xl_grid%xb(i+1) - xl_grid%xb(i)
        !     A_mat(i,   i  ) = A_mat(i,   i  ) - 1.0_dp / h
        !     A_mat(i,   i+1) = A_mat(i,   i+1) + 1.0_dp / h
        !     A_mat(i+1, i  ) = A_mat(i+1, i  ) + 1.0_dp / h
        !     A_mat(i+1, i+1) = A_mat(i+1, i+1) - 1.0_dp / h
        ! end do

        ! ---- Left boundary: Neumann BC (dPhi/dx = 0) ----
        ! This modifies only the first row/col to enforce derivative = 0
        hR = xl_grid%xb(2) - xl_grid%xb(1)
        A_mat(1,:) = cmplx(0.0d0, 0.0d0, dp)
        A_mat(1,1) = -1.0d0 / hR
        A_mat(1,2) =  1.0d0 / hR
        ! Keep symmetry in column 1
        A_mat(2,1) =  1.0d0 / hR

        ! ---- Right boundary: Dirichlet BC (Phi = 0) ----
        ! Overwrite last row to enforce Phi_n = 0
        A_mat(n,:) = cmplx(0.0d0, 0.0d0, dp)
        A_mat(n,n) = cmplx(1.0d0, 0.0d0, dp)
        
        call write_matrix(trim(output_path)//'kernel/laplacian_re.dat', real(A_mat), xl_grid%npts_b, xl_grid%npts_b)

    end subroutine


    subroutine dense_to_sparse(A, irow, pcol, A_nz, nrow, ncol, nz_out)

        use sparse_mod, only: column_full2pointer
        use KIM_kinds, only: dp
        
        implicit none

        complex(dp), dimension(:,:), intent(in) :: A
        complex(dp), dimension(:), allocatable, intent(inout) :: A_nz
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
