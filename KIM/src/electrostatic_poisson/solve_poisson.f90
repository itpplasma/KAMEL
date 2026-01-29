module poisson_solver_m

    implicit none

    contains

    ! Solve \Delta \Phi + 4\pi K_rho_phi \Phi = - 4 \pi K_rho_B B_r
    ! using sparse matrix solver
    subroutine solve_poisson(K_rho_phi, K_rho_B, phi_sol)

        use config_m, only: fstatus, fdebug
        use sparse_mod, only: sp2fullComplex, sparse_solveComplex_b1, sparse_solve_method, sparse_talk
        use config_m, only: output_path
        use constants_m, only: pi
        use grid_m, only: xl_grid, calc_mass_matrix, M_mat
        use setup_m, only: type_br_field
        use KIM_kinds_m, only: dp
        use IO_collection_m, only: write_matrix, write_complex_profile

        implicit none

        complex(dp), intent(in) :: K_rho_phi(:,:)
        complex(dp), intent(in) :: K_rho_B(:,:)
        complex(dp), dimension(:), allocatable, intent(out) :: phi_sol
        complex(dp), dimension(:), allocatable :: A_nz ! non-zero elements of A matrix
        complex(dp), dimension(:,:), allocatable :: A_mat ! A matrix (stiffness matrix in the beginning, then full right hand side matrix)

        complex(dp), dimension(:), allocatable :: b_vec ! b vector and x vector
        integer, dimension(:), allocatable :: irow, pcol
        integer :: nz_out, nrow, ncol
        integer :: i, j
        integer :: sparse_solver_option

        if (fstatus >= 1) write(*,*) 'Status: solve poisson equation'
        if (fdebug < 2) sparse_talk = .false. ! turn off sparse solver output for low fdebug

        allocate(A_mat(xl_grid%npts_b, xl_grid%npts_b))

        call prepare_Laplace_matrix(A_mat)

        if (.not. allocated(M_mat)) then
            allocate(M_mat(xl_grid%npts_b, xl_grid%npts_b))
            call calc_mass_matrix(M_mat)
        end if

        call write_matrix('kernel/mass_matrix.dat', M_mat, xl_grid%npts_b, xl_grid%npts_b, &
            'Mass matrix M', 'cm')

        A_mat = (A_mat + 4.0d0 * pi * K_rho_phi)

        if (fdebug == 3) then
            call write_A_matrix_sparse_check_to_file
        end if

        sparse_solve_method = 1
        sparse_solver_option = 0

        call create_rhs_vector(type_br_field, K_rho_B, b_vec)
        call impose_bc_on_matrix_and_rhs(A_mat, b_vec, K_rho_phi, K_rho_B)

        call write_matrix('kernel/A_mat', real(A_mat), xl_grid%npts_b, xl_grid%npts_b, &
            'Stiffness matrix A (left hand matrix of Ax=b)', '1/cm^2')
        ! call write_matrix('kernel/A_mat_im.dat', dimag(A_mat), xl_grid%npts_b, xl_grid%npts_b)
        call write_complex_profile(xl_grid%xb, b_vec, xl_grid%npts_b, 'kernel/b_vec', &
            'Right hand side vector b of Poisson equation Ax=b', 'statC/cm^3') !TODO: check units

        call dense_to_sparse(A_mat, irow, pcol, A_nz, nrow, ncol, nz_out)
        call sparse_solveComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, sparse_solver_option)

        phi_sol = b_vec

        contains

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

    end subroutine

    subroutine create_rhs_vector(type, K_rho_B, rhs_vec)

        use resonances_mod, only: index_rg_res, r_res
        use functions_m, only: varphi_l
        use grid_m, only: xl_grid
        use IO_collection_m, only: write_complex_profile, plot_profile
        use KIM_kinds_m, only: dp
        use fields_m, only: EBdat, set_Br_field
        use config_m, only: output_path
        use constants_m, only: pi, e_charge

        implicit none

        integer, intent(in) :: type
        complex(dp), allocatable, intent(out) :: rhs_vec(:)
        complex(dp), intent(in) :: K_rho_B(:,:)

        allocate(rhs_vec(xl_grid%npts_b))
        rhs_vec = cmplx(0.0d0, 0.0d0, dp)

        if (type ==1) then ! constant br
            rhs_vec = cmplx(1.0d0, 0.0d0, dp) * e_charge
        ! type > 10 uses Br kernel to determine right hand side of equation
        elseif(type==11) then
            call set_Br_field(EBdat, 1) ! Br field from input
            rhs_vec = matmul(K_rho_B, EBdat%Br)
        elseif(type==12) then
            call set_Br_field(EBdat, 0) ! constant Br field
            rhs_vec = matmul(K_rho_B, EBdat%Br)
        end if

        rhs_vec = - 4d0 * pi * rhs_vec

    end subroutine

    subroutine prepare_Laplace_matrix(A_mat)

        use grid_m, only: xl_grid
        use KIM_kinds_m, only: dp
        use config_m, only: output_path
        use IO_collection_m, only: write_matrix

        implicit none

        complex(dp), intent(inout) :: A_mat(:,:)
        real(dp) :: h, hL, hR
        integer :: i
        integer :: n

        n = xl_grid%npts_b

        A_mat = cmplx(0.0d0, 0.0d0, dp)

        do i = 2, n-1
            ! distinguish left and right element sizes for non-uniform grid
            hL = xl_grid%xb(i)   - xl_grid%xb(i-1)  ! left element size
            hR = xl_grid%xb(i+1) - xl_grid%xb(i)    ! right element size

            A_mat(i,i-1) = A_mat(i,i-1) + 1.0d0/hL
            A_mat(i,i)   = A_mat(i,i)   - 1.0d0/hL - 1.0d0/hR
            A_mat(i,i+1) = A_mat(i,i+1) + 1.0d0/hR
        end do

        call write_matrix('kernel/Laplace_in_FEM', real(A_mat), xl_grid%npts_b, xl_grid%npts_b, &
            'Laplace operator matrix in FEM discretization', '1/cm^2')

    end subroutine

    subroutine impose_bc_on_matrix_and_rhs(A_mat, b_vec, K_rho_phi, K_rho_B)

        use config_m, only: fdebug
        use grid_m, only: xl_grid
        use KIM_kinds_m, only: dp
        use setup_m, only: bc_type
        use fields_m, only: calculate_phi_aligned, EBdat
        use species_m, only: plasma

        implicit none

        complex(dp), intent(inout) :: A_mat(:,:), b_vec(:)
        complex(dp), intent(in) :: K_rho_phi(:,:), K_rho_B(:,:)
        complex(dp) :: phi_boundary_left, phi_boundary_right

        if(bc_type == 3) then ! zero misalignment field at boundaries

            call calculate_phi_aligned(plasma, EBdat)

            phi_boundary_left = - EBdat%phi_aligned(1)
            phi_boundary_right = - EBdat%phi_aligned(xl_grid%npts_b)

            if (fdebug == 1) print *, "Imposing BC of zero misalignment field: Phi_left = ", &
                phi_boundary_left, ", Phi_right = ", phi_boundary_right

            b_vec(2:xl_grid%npts_b-1) = b_vec(2:xl_grid%npts_b-1) &
                - A_mat(2:xl_grid%npts_b-1,1) * phi_boundary_left &
                - A_mat(2:xl_grid%npts_b-1,xl_grid%npts_b) * phi_boundary_right

            A_mat(:, 1) = 0.0d0
            A_mat(1, :) = 0.0d0
            A_mat(xl_grid%npts_b, :) = 0.0d0
            A_mat(:, xl_grid%npts_b) = 0.0d0
            A_mat(1, 1) = 1.0d0
            A_mat(xl_grid%npts_b, xl_grid%npts_b) = 1.0d0

            b_vec(1) = phi_boundary_left
            b_vec(xl_grid%npts_b) = phi_boundary_right
        end if

    end subroutine

    subroutine dense_to_sparse(A, irow, pcol, A_nz, nrow, ncol, nz_out)

        use sparse_mod, only: column_full2pointer
        use KIM_kinds_m, only: dp

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
                if (abs(A(nr,nc)) > 1.0e-14_dp) then
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
                if (abs(A(nr, nc)) > 1.0e-14_dp) then
                    n = n + 1
                    irow(n) = nr
                    icol(n) = nc
                    ! A_nz(n) = A_nz(n) + A(nr, nc)
                    A_nz(n) = A(nr, nc)
                end if
            end do
        end do

        call column_full2pointer(icol, pcol)

        if (present(nz_out)) nz_out = nz
        if (allocated(icol)) deallocate(icol)

    end subroutine

end module
