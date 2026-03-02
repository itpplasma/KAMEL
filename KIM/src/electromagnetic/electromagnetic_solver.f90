! Run type for electromagnetic model
! Solves coupled Poisson-Ampere block system for Phi and A_parallel self-consistently
module rt_electromagnetic_m

    use kim_base_m, only: kim_t

    implicit none

    type, extends(kim_t) :: electromagnetic_t
        contains
            procedure :: init => init_electromagnetic
            procedure :: run  => run_electromagnetic
    end type electromagnetic_t

    contains

    subroutine init_electromagnetic(this)

        use species_m, only: plasma, set_plasma_quantities
        use IO_collection_m, only: create_output_directories
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use grid_m, only: rg_grid

        implicit none

        class(electromagnetic_t), intent(inout) :: this

        this%run_type = "electromagnetic"

        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)

        print *, "..."//trim(this%run_type)//" model initialized."

    end subroutine

    subroutine run_electromagnetic(this)

        use kernel_m, only: FP_fill_kernels, kernel_spl_t, write_kernels
        use kernel_adaptive_m, only: FP_fill_kernels_adaptive
        use grid_m, only: xl_grid, calc_mass_matrix, M_mat, theta_integration
        use IO_collection_m, only: write_complex_profile_abs
        use poisson_solver_m, only: prepare_Laplace_matrix, dense_to_sparse
        use config_m, only: output_path, collision_model, fstatus, fdebug
        use fields_m, only: EBdat, postprocess_electric_field, &
                            calculate_charge_density, calculate_current_density
        use ampere_matrices_m, only: interpolate_equil_to_xl
        use setup_m, only: m_mode, n_mode, R0, Br_boundary_re, Br_boundary_im, bc_type
        use constants_m, only: pi, sol, com_unit
        use KIM_kinds_m, only: dp
        use species_m, only: plasma
        use sparse_mod, only: sparse_solveComplex_b1, sparse_solve_method, sparse_talk

        implicit none

        class(electromagnetic_t), intent(inout) :: this

        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp
        type(kernel_spl_t) :: kernel_j_phi_llp
        type(kernel_spl_t) :: kernel_j_B_llp

        complex(dp), allocatable :: A_block(:,:)
        complex(dp), allocatable :: b_block(:)
        complex(dp), allocatable :: A_nz(:)
        integer, allocatable :: irow(:), pcol(:)
        integer :: nz_out, nrow_sp, ncol_sp

        complex(dp), allocatable :: A_Phi(:,:)
        real(dp), allocatable :: laplace_perp(:,:)
        real(dp), allocatable :: hz_xl(:), hth_xl(:)
        real(dp), allocatable :: ks_xl(:)
        complex(dp), allocatable :: alpha(:)

        complex(dp), allocatable :: jpar(:), rho(:)
        complex(dp) :: Br_boundary

        integer :: N, i, j
        real(dp) :: kz, hL, hR
        character(8) :: date
        character(10) :: time
        character(5) :: zone
        integer, dimension(8) :: values

        ! Validate collision model
        if (trim(collision_model) /= "FokkerPlanck") then
            print *, "Error: electromagnetic run type requires collision_model = FokkerPlanck"
            print *, "Current collision_model: ", trim(collision_model)
            error stop
        end if

        N = xl_grid%npts_b
        kz = dble(n_mode) / R0
        Br_boundary = cmplx(Br_boundary_re, Br_boundary_im, dp)

        ! Initialize kernels
        call kernel_rho_phi_llp%init_kernel(N, N)
        call kernel_rho_B_llp%init_kernel(N, N)
        call kernel_j_phi_llp%init_kernel(N, N)
        call kernel_j_B_llp%init_kernel(N, N)

        ! Fill kernels
        call date_and_time(date, time, zone, values)
        write(*,*) "Start filling kernel at ", date, " ", time, " ..."

        select case (trim(theta_integration))
            case ("GaussLegendre")
                call FP_fill_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, &
                    kernel_j_phi_llp, kernel_j_B_llp)
            case ("RKF45", "QUADPACK")
                call FP_fill_kernels_adaptive(kernel_rho_phi_llp, kernel_rho_B_llp, &
                    kernel_j_phi_llp, kernel_j_B_llp)
            case default
                error stop "Error: theta integration method not recognized."
        end select

        call write_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, &
            kernel_j_phi_llp, kernel_j_B_llp)

        ! Compute mass matrix (needed by calculate_current_density)
        if (.not. allocated(M_mat)) then
            allocate(M_mat(N, N))
            call calc_mass_matrix(M_mat)
        end if

        ! Compute alpha(r) = i * (m/r * h_z - kz * h_theta) for A_par -> Br conversion
        call interpolate_equil_to_xl(hz_xl, hth_xl)
        allocate(alpha(N))
        do i = 1, N
            alpha(i) = com_unit * (dble(m_mode) / xl_grid%xb(i) * hz_xl(i) - kz * hth_xl(i))
        end do

        ! Compute ks(r) = (m * h_z - kz * h_theta) / r (perpendicular wavenumber)
        allocate(ks_xl(N))
        do i = 1, N
            ks_xl(i) = (dble(m_mode) * hz_xl(i) - kz * hth_xl(i)) / xl_grid%xb(i)
        end do

        ! Build laplace_perp: nabla^2_perp = d^2/dr^2 + (1/r)d/dr - ks^2
        allocate(laplace_perp(N, N))
        laplace_perp = 0.0d0
        do i = 2, N-1
            hL = xl_grid%xb(i) - xl_grid%xb(i-1)
            hR = xl_grid%xb(i+1) - xl_grid%xb(i)
            laplace_perp(i,i-1) = 1.0d0 / hL - 1.0d0 / (xl_grid%xb(i) * (hL + hR))
            laplace_perp(i,i)   = -(1.0d0 / hL + 1.0d0 / hR) - ks_xl(i)**2
            laplace_perp(i,i+1) = 1.0d0 / hR + 1.0d0 / (xl_grid%xb(i) * (hL + hR))
        end do

        ! Assemble Poisson LHS block: A_Phi = Delta + 4*pi * K_rho_phi
        allocate(A_Phi(N, N))
        call prepare_Laplace_matrix(A_Phi)
        A_Phi = A_Phi + 4.0d0 * pi * kernel_rho_phi_llp%Kllp

        ! Assemble 2N x 2N block system
        allocate(A_block(2*N, 2*N), b_block(2*N))
        A_block = cmplx(0.0d0, 0.0d0, dp)
        b_block = cmplx(0.0d0, 0.0d0, dp)

        ! Top-left: A_Phi
        A_block(1:N, 1:N) = A_Phi

        ! Top-right: 4*pi * K_rho_B * D_alpha  (since Br = alpha * A_par)
        do j = 1, N
            A_block(1:N, N+j) = 4.0d0 * pi * kernel_rho_B_llp%Kllp(1:N, j) * alpha(j)
        end do

        ! Bottom-left: (4*pi/c) * K_j_phi
        A_block(N+1:2*N, 1:N) = 4.0d0 * pi / sol * kernel_j_phi_llp%Kllp

        ! Bottom-right: laplace_perp + (4*pi/c) * K_j_B * D_alpha
        A_block(N+1:2*N, N+1:2*N) = cmplx(laplace_perp, 0.0d0, dp)
        do j = 1, N
            A_block(N+1:2*N, N+j) = A_block(N+1:2*N, N+j) &
                + 4.0d0 * pi / sol * kernel_j_B_llp%Kllp(1:N, j) * alpha(j)
        end do

        ! --- Boundary conditions ---

        ! Phi BCs: Dirichlet Phi = 0 at both boundaries
        A_block(1, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(1, 1) = cmplx(1.0d0, 0.0d0, dp)
        b_block(1) = cmplx(0.0d0, 0.0d0, dp)

        A_block(N, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(N, N) = cmplx(1.0d0, 0.0d0, dp)
        b_block(N) = cmplx(0.0d0, 0.0d0, dp)

        ! A_par BCs: A_par(left) = 0, A_par(right) = Br_boundary / alpha(N)
        A_block(N+1, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(N+1, N+1) = cmplx(1.0d0, 0.0d0, dp)
        b_block(N+1) = cmplx(0.0d0, 0.0d0, dp)

        A_block(2*N, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(2*N, 2*N) = cmplx(1.0d0, 0.0d0, dp)
        b_block(2*N) = Br_boundary / alpha(N)

        ! Zero-misalignment BC override for Phi
        if (bc_type == 3) then
            call apply_zero_misalignment_bc(A_block, b_block, N, Br_boundary)
        end if

        if (fstatus >= 1) write(*,*) 'Status: solving coupled Poisson-Ampere block system for (Phi, A_par) (2N =', 2*N, ')'

        ! Solve
        if (fdebug < 2) sparse_talk = .false.
        sparse_solve_method = 1

        call dense_to_sparse(A_block, irow, pcol, A_nz, nrow_sp, ncol_sp, nz_out)
        call sparse_solveComplex_b1(nrow_sp, ncol_sp, nz_out, irow, pcol, A_nz, b_block, 0)

        ! Extract solution: solve gives (Phi, A_par), recover Br = alpha * A_par
        allocate(EBdat%Phi(N), EBdat%Apar(N), EBdat%Br(N), EBdat%r_grid(N))
        EBdat%r_grid = xl_grid%xb
        EBdat%Phi  = b_block(1:N)
        EBdat%Apar = b_block(N+1:2*N)
        EBdat%Br   = alpha * EBdat%Apar

        ! Write fields
        call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, N, "/fields/Phi", &
            'Electrostatic potential perturbation Phi (self-consistent)', 'statV')
        call write_complex_profile_abs(xl_grid%xb, EBdat%Apar, N, "/fields/Apar", &
            'Parallel vector potential A_parallel (self-consistent)', 'G*cm')
        call write_complex_profile_abs(xl_grid%xb, EBdat%Br, N, "/fields/Br_selfconsistent", &
            'Self-consistent radial magnetic field perturbation Br', 'G')

        ! Current density
        call calculate_current_density(jpar, EBdat%Phi, EBdat%Br, &
            kernel_j_phi_llp%Kllp, kernel_j_B_llp%Kllp)
        call write_complex_profile_abs(xl_grid%xb, jpar, N, "/fields/jpar", &
            'Parallel current density from self-consistent solve', 'statA/cm^2')

        ! Electric field postprocessing
        call postprocess_electric_field(EBdat)

        ! Charge density
        call calculate_charge_density(rho, EBdat)
        call write_complex_profile_abs(xl_grid%xb, rho, N, "/fields/rho", &
            'Charge density from self-consistent solve', 'statC/cm^3')

        ! Cleanup
        deallocate(A_block, b_block, A_Phi, laplace_perp, alpha, hz_xl, hth_xl, ks_xl)

    end subroutine

    subroutine apply_zero_misalignment_bc(A_block, b_block, N, Br_boundary)

        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit
        use species_m, only: plasma
        use grid_m, only: xl_grid

        implicit none

        complex(dp), intent(inout) :: A_block(:,:), b_block(:)
        integer, intent(in) :: N
        complex(dp), intent(in) :: Br_boundary

        integer :: ir, ibeg, iend
        integer :: nlagr = 4, nder = 0
        real(dp), allocatable :: coef(:,:)
        real(dp) :: Er_right, kp_right, B0_right
        complex(dp) :: phi_right

        allocate(coef(0:nder, nlagr))

        ! Left boundary: Br = 0, so phi_aligned = 0, Phi_left = 0 (unchanged)

        ! Right boundary: Br = Br_boundary
        call binsrc(plasma%r_grid, 1, size(plasma%r_grid), xl_grid%xb(N), ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > size(plasma%r_grid)) then
            iend = size(plasma%r_grid)
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, xl_grid%xb(N), plasma%r_grid(ibeg:iend), coef)
        Er_right = sum(coef(0,:) * plasma%Er(ibeg:iend))
        kp_right = sum(coef(0,:) * plasma%kp(ibeg:iend))
        B0_right = sum(coef(0,:) * plasma%B0(ibeg:iend))

        phi_right = -com_unit * Er_right * Br_boundary / (B0_right * kp_right)

        b_block(N) = phi_right

        deallocate(coef)

    end subroutine

end module
