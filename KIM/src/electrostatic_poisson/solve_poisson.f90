module poisson_solver_m

    contains
    
    ! Solve \Delta \Phi + 4\pi K_rho_phi \Phi = - 4 \pi K_rho_B B_r
    ! using sparse matrix solver
    subroutine solve_poisson(K_rho_phi, K_rho_B, phi_sol)

        use config_m, only: fstatus, fdebug
        use sparse_mod, only: sp2fullComplex, sparse_solveComplex_b1, sparse_solve_method
        use config_m, only: output_path
        use constants_m, only: pi
        use grid_m, only: xl_grid, calc_mass_matrix
        use setup_m, only: type_br_field, bc_type
        use KIM_kinds_m, only: dp

        implicit none

        complex(dp), intent(in) :: K_rho_phi(:,:)
        complex(dp), intent(in) :: K_rho_B(:,:)
        complex(dp), dimension(:), allocatable, intent(out) :: phi_sol
        complex(dp), dimension(:), allocatable :: A_nz ! non-zero elements of A matrix
        complex(dp), dimension(:,:), allocatable :: A_mat ! A matrix (stiffness matrix in the beginning, then full right hand side matrix)
        real(dp), dimension(:,:), allocatable :: M_mat ! mass matrix
        complex(dp), dimension(:), allocatable :: b_vec ! b vector and x vector
        complex(dp), dimension(:,:), allocatable :: inv_K_rho_phi ! inverse of K_rho_phi
        integer, dimension(:), allocatable :: irow, pcol
        integer :: nz_out, nrow, ncol
        integer :: i,j
        integer :: sparse_solver_option

        if (fstatus == 1) write(*,*) 'Status: solve poisson equation'

        allocate(A_mat(xl_grid%npts_b, xl_grid%npts_b), M_mat(xl_grid%npts_b, xl_grid%npts_b))
        call prepare_Laplace_matrix(A_mat)
        call calc_mass_matrix(M_mat)

        A_mat = (A_mat + 4.0d0 * pi * K_rho_phi) 

        if (fdebug == 3) then
            call write_A_matrix_sparse_check_to_file
        end if

        call create_rhs_vector(type_br_field, K_rho_B, b_vec)
        call impose_bc_on_matrix_and_rhs(A_mat, b_vec, K_rho_phi, K_rho_B)

        call dense_to_sparse(A_mat, irow, pcol, A_nz, nrow, ncol, nz_out)

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

    end subroutine

    subroutine create_rhs_vector(type, K_rho_B, rhs_vec)

        use resonances_mod, only: index_rg_res, r_res
        use functions_m, only: varphi_l
        use grid_m, only: xl_grid
        use IO_collection_m, only: write_profile, write_complex_profile, plot_profile
        use KIM_kinds_m, only: dp
        use fields_m, only: EBdat, set_Br_field
        use config_m, only: output_path
        use setup_m, only: bc_type
        use constants_m, only: pi, e_charge

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
            do i = int(5.0d0/6.0d0 *xl_grid%npts_b), xl_grid%npts_b
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
            do i = 1, xl_grid%npts_b
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
        if (bc_type == 1) then
            rhs_vec(xl_grid%npts_b) = 0.0_dp
        end if

        call write_complex_profile(xl_grid%xb, rhs_vec, xl_grid%npts_b, trim(output_path)//'fields/rhs_vec.dat')

    end subroutine

    subroutine prepare_Laplace_matrix(A_mat)

        use grid_m, only: xl_grid
        use KIM_kinds_m, only: dp
        use config_m, only: output_path
        use IO_collection_m, only: write_matrix
        use setup_m, only: bc_type

        implicit none

        complex(dp), intent(inout) :: A_mat(:,:)
        real(dp) :: h, hL, hR
        integer :: i
        integer :: n

        n = xl_grid%npts_b
        
        ! create Laplacian:
        A_mat = cmplx(0.0d0, 0.0d0, dp)

        do i = 2, n-1
            hL = xl_grid%xb(i)   - xl_grid%xb(i-1)  ! left element size
            hR = xl_grid%xb(i+1) - xl_grid%xb(i)    ! right element size

            A_mat(i,i-1) = A_mat(i,i-1) + 1.0d0/hL
            A_mat(i,i)   = A_mat(i,i)   - 1.0d0/hL - 1.0d0/hR
            A_mat(i,i+1) = A_mat(i,i+1) + 1.0d0/hR
        end do

        if (bc_type == 1)then
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
        end if

        call write_matrix(trim(output_path)//'kernel/laplacian_re.dat', real(A_mat), xl_grid%npts_b, xl_grid%npts_b)

    end subroutine

    subroutine impose_bc_on_matrix_and_rhs(A_mat, b_vec, K_rho_phi, K_rho_B)

        use grid_m, only: xl_grid
        use KIM_kinds_m, only: dp
        use setup_m, only: bc_type
        use fields_m, only: calculate_phi_aligned, EBdat
        use species_m, only: plasma

        implicit none

        complex(dp), intent(inout) :: A_mat(:,:), b_vec(:)
        complex(dp), intent(in) :: K_rho_phi(:,:), K_rho_B(:,:)
        complex(dp) :: phi_boundary_left, phi_boundary_right

        if (bc_type == 2) then
            block
                complex(dp), allocatable :: inv_K_rho_phi(:,:)

                integer :: n
                complex(dp), allocatable :: phi_temp(:)

                n = xl_grid%npts_b

                allocate(phi_temp(n))

                ! Enforce Dirichlet BC at right boundary: Phi_n = (K_rho_phi^-1 * b_vec)_n

                call invert_complex_matrix(K_rho_phi, inv_K_rho_phi)

                phi_temp = matmul(inv_K_rho_phi, b_vec)
                phi_boundary_left = - phi_temp(1)
                phi_boundary_right = - phi_temp(n)
                deallocate(phi_temp)

                b_vec = b_vec - A_mat(:,1) * phi_boundary_left - A_mat(:,n) * phi_boundary_right

                A_mat(:,1) = 0.0d0
                A_mat(1,:) = 0.0d0
                A_mat(n,:) = 0.0d0
                A_mat(:,n) = 0.0d0
                A_mat(1,1) = 1.0d0
                A_mat(n,n) = 1.0d0

                b_vec(1) = phi_boundary_left
                b_vec(n) = phi_boundary_right
            end block
        else if(bc_type == 3) then ! zero misalignment field at boundaries
            block
                integer :: n

                n = xl_grid%npts_b

                call calculate_phi_aligned(plasma, EBdat)

                phi_boundary_left = - EBdat%phi_aligned(1)
                phi_boundary_right = - EBdat%phi_aligned(n)

                print *, "Imposing BC of zero misalignment field: Phi_left = ", phi_boundary_left, ", Phi_right = ", phi_boundary_right

                b_vec = b_vec - A_mat(:,1) * phi_boundary_left - A_mat(:,n) * phi_boundary_right

                A_mat(:,1) = 0.0d0
                A_mat(1,:) = 0.0d0
                A_mat(n,:) = 0.0d0
                A_mat(:,n) = 0.0d0
                A_mat(1,1) = 1.0d0
                A_mat(n,n) = 1.0d0

                b_vec(1) = phi_boundary_left
                b_vec(n) = phi_boundary_right
            end block
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
                if (.not.(A(nr, nc) == (0.0d0, 0.0d0))) then
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

    ! Invert a square complex(dp) matrix using LAPACK (zgetrf + zgetri)
    subroutine invert_complex_matrix(A_in, A_inv)

        use KIM_kinds_m, only: dp

        implicit none

        complex(dp), intent(in)  :: A_in(:,:)
        complex(dp), allocatable, intent(out) :: A_inv(:,:)

        complex(dp), allocatable :: A(:,:), work(:)
        integer, allocatable     :: ipiv(:)
        integer :: n, lda, info, lwork
        complex(dp) :: work_query(1)

        interface
            subroutine zgetrf(m, n, a, lda, ipiv, info)
                use KIM_kinds_m, only: dp
                integer, intent(in) :: m, n, lda
                integer, intent(out) :: ipiv(*)
                integer, intent(out) :: info
                complex(dp) :: a(lda,*)
            end subroutine zgetrf
            subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
                use KIM_kinds_m, only: dp
                integer, intent(in) :: n, lda, lwork
                integer, intent(in) :: ipiv(*)
                integer, intent(out) :: info
                complex(dp) :: a(lda,*), work(*)
            end subroutine zgetri
        end interface

        if (size(A_in,1) /= size(A_in,2)) then
            stop "invert_complex_matrix: input must be square"
        end if

        n = size(A_in,1)
        lda = max(1, n)

        allocate(A(n,n))
        A = A_in

        allocate(ipiv(n))

        call zgetrf(n, n, A, lda, ipiv, info)
        if (info /= 0) then
            stop "invert_complex_matrix: LU factorization failed"
        end if

        ! Workspace query
        lwork = -1
        call zgetri(n, A, lda, ipiv, work_query, lwork, info)
        lwork = max(1, int(real(work_query(1))))
        allocate(work(lwork))

        call zgetri(n, A, lda, ipiv, work, lwork, info)
        if (info /= 0) then
            stop "invert_complex_matrix: inversion failed"
        end if

        allocate(A_inv(n,n))
        A_inv = A

        deallocate(A, ipiv, work)

    end subroutine invert_complex_matrix

end module
