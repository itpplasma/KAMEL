! Run type for the enforced-periodicity electrostatic solve (issue #175):
! periodize the background around the resonant surface, assemble the
! discrete-Fourier kernel matrices, and solve the screened Poisson system.
! Supported perturbation: type_br_field = 12 (constant unit Br), matching
! the hat-basis kernel-weighted right-hand side.
module rt_fourier_periodic_m

    use kim_base_m, only: kim_t
    use KIM_kinds_m, only: dp

    implicit none

    type, extends(kim_t) :: fourier_periodic_t
        contains
            procedure :: init => init_fourier_periodic
            procedure :: run => run_fourier_periodic
    end type fourier_periodic_t

    integer, parameter :: n_output = 512

contains

    subroutine init_fourier_periodic(this)

        use species_m, only: plasma, set_plasma_quantities
        use IO_collection_m, only: create_output_directories
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use grid_m, only: rg_grid
        use setup_m, only: fp_r_res, fp_dr_layer, fp_dr_transition, &
                           fp_grid_points
        use config_m, only: rescale_density

        implicit none

        class(fourier_periodic_t), intent(inout) :: this

        this%run_type = "fourier_periodic"

        if (rescale_density) then
            error stop "fourier_periodic: rescale_density is not supported; &
                &the periodized collision frequencies would rescale twice"
        end if
        if (fp_r_res <= 0.0d0) then
            error stop "fourier_periodic: set fp_r_res > 0 in KIM_SETUP"
        end if
        if (fp_dr_layer <= 0.0d0 .or. fp_dr_transition <= 0.0d0) then
            error stop "fourier_periodic: set fp_dr_layer > 0 and &
                &fp_dr_transition > 0 in KIM_SETUP"
        end if
        if (fp_grid_points < 8) then
            error stop "fourier_periodic: fp_grid_points must be at least 8"
        end if

        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)

        call build_periodized_background(fp_r_res, fp_dr_layer, &
                                         fp_dr_transition, fp_grid_points)

        print *, "..."//trim(this%run_type)//" model initialized."

    end subroutine init_fourier_periodic

    ! Replace the global plasma state with its periodization around r_res:
    ! primitives (n, T, Er, B0, q, ks, kp) are periodized, om_E is
    ! recomputed algebraically, and every derived quantity is recomputed
    ! from the periodized primitives. Ordering invariant: gradients come
    ! from differentiating the periodized profiles, never from
    ! periodizing gradients. The hat-basis cell-center and susceptibility
    ! tables keep their original grids; this run type does not use them.
    subroutine build_periodized_background(r_res, dr_layer, dr_transition, &
                                           n_points)

        use species_m, only: plasma, calc_plasma_parameter_derivs, &
                             calculate_plasma_backs
        use periodize_m, only: periodized_profile_value
        use constants_m, only: sol, e_charge, ev

        implicit none

        real(dp), intent(in) :: r_res, dr_layer, dr_transition
        integer, intent(in) :: n_points

        real(dp), allocatable :: r_old(:), grid(:)
        real(dp) :: half_width
        integer :: i, sp

        half_width = dr_layer + dr_transition
        if (plasma%r_grid(1) > r_res - 3.0d0*half_width .or. &
            plasma%r_grid(plasma%grid_size) < r_res + 3.0d0*half_width) then
            print *, "fourier_periodic: plasma grid covers ", &
                plasma%r_grid(1), " to ", plasma%r_grid(plasma%grid_size)
            print *, "fourier_periodic: periodization needs ", &
                r_res - 3.0d0*half_width, " to ", r_res + 3.0d0*half_width
            error stop "fourier_periodic: background does not cover the &
                &periodization support"
        end if

        allocate (r_old(plasma%grid_size))
        r_old = plasma%r_grid
        allocate (grid(n_points))
        do i = 1, n_points
            grid(i) = r_res - 3.0d0*half_width &
                      + 6.0d0*half_width*(i - 1)/(n_points - 1)
        end do

        do sp = 0, plasma%n_species - 1
            call replace_periodized(plasma%spec(sp)%n)
            call replace_periodized(plasma%spec(sp)%T)
        end do
        call replace_periodized(plasma%Er)
        call replace_periodized(plasma%B0)
        call replace_periodized(plasma%q)
        call replace_periodized(plasma%ks)
        call replace_periodized(plasma%kp)

        deallocate (plasma%r_grid)
        allocate (plasma%r_grid(n_points))
        plasma%r_grid = grid
        plasma%grid_size = n_points

        deallocate (plasma%om_E)
        allocate (plasma%om_E(n_points))
        plasma%om_E = -sol*plasma%ks*plasma%Er/plasma%B0

        do sp = 0, plasma%n_species - 1
            associate (spec => plasma%spec(sp))
                if (allocated(spec%dndr)) deallocate (spec%dndr, spec%dTdr)
                if (allocated(spec%vT)) then
                    deallocate (spec%vT, spec%omega_c, spec%lambda_D, spec%nu)
                end if
                if (allocated(spec%z0)) deallocate (spec%z0)
                if (allocated(spec%rho_L)) deallocate (spec%rho_L)
                if (allocated(spec%A1)) deallocate (spec%A1, spec%A2)
            end associate
        end do
        if (allocated(plasma%dqdr)) deallocate (plasma%dqdr)

        call calc_plasma_parameter_derivs
        call calculate_plasma_backs(plasma)

        do sp = 0, plasma%n_species - 1
            associate (spec => plasma%spec(sp))
                allocate (spec%A1(n_points), spec%A2(n_points))
                do i = 1, n_points
                    spec%A1(i) = spec%dndr(i)/spec%n(i) &
                        - spec%Zspec*e_charge/(spec%T(i)*ev)*plasma%Er(i) &
                        - 3.0d0/(2.0d0*spec%T(i))*spec%dTdr(i)
                    spec%A2(i) = spec%dTdr(i)/spec%T(i)
                end do
            end associate
        end do

    contains

        subroutine replace_periodized(profile)

            implicit none

            real(dp), allocatable, intent(inout) :: profile(:)

            real(dp), allocatable :: old_values(:)
            integer :: j

            call move_alloc(profile, old_values)
            allocate (profile(n_points))
            do j = 1, n_points
                profile(j) = periodized_profile_value(r_old, old_values, &
                    r_res, dr_layer, dr_transition, grid(j))
            end do
            deallocate (old_values)

        end subroutine replace_periodized

    end subroutine build_periodized_background

    subroutine run_fourier_periodic(this)

        use fourier_periodic_m, only: periodic_k_grid, &
            assemble_periodic_kernel, solve_periodic_potential, &
            periodic_field_value
        use kernels_m, only: kernel_rho_phi_of_kr_krp_rg, &
                             kernel_rho_B_of_kr_krp_rg
        use setup_m, only: type_br_field, fp_r_res, fp_dr_layer, &
                           fp_dr_transition, fp_n_modes, fp_n_quad
        use constants_m, only: pi
        use IO_collection_m, only: write_complex_profile

        implicit none

        class(fourier_periodic_t), intent(inout) :: this

        real(dp) :: period, r_out(n_output)
        real(dp), allocatable :: k_modes(:)
        complex(dp), allocatable :: kernel_rho_phi_matrix(:, :)
        complex(dp), allocatable :: kernel_rho_B_matrix(:, :)
        complex(dp), allocatable :: br_modes(:), rhs(:), phi_modes(:)
        complex(dp) :: phi_out(n_output), br_out(n_output)
        integer :: i, n_k

        if (type_br_field /= 12) then
            error stop "fourier_periodic: only type_br_field = 12 &
                &(constant unit Br) is supported"
        end if

        period = 2.0d0*(fp_dr_layer + fp_dr_transition)
        n_k = 2*fp_n_modes + 1
        allocate (k_modes(n_k), br_modes(n_k))
        k_modes = periodic_k_grid(period, fp_n_modes)

        kernel_rho_phi_matrix = assemble_periodic_kernel( &
            kernel_rho_phi_of_kr_krp_rg, k_modes, fp_r_res, period, fp_n_quad)
        kernel_rho_B_matrix = assemble_periodic_kernel( &
            kernel_rho_B_of_kr_krp_rg, k_modes, fp_r_res, period, fp_n_quad)

        br_modes = (0.0d0, 0.0d0)
        br_modes(fp_n_modes + 1) = (1.0d0, 0.0d0)
        rhs = 4.0d0*pi*matmul(kernel_rho_B_matrix, br_modes)
        phi_modes = solve_periodic_potential(k_modes, &
                                             kernel_rho_phi_matrix, rhs)

        do i = 1, n_output
            r_out(i) = fp_r_res - 0.5d0*period &
                       + period*(i - 1)/(n_output - 1)
            phi_out(i) = periodic_field_value(k_modes, phi_modes, r_out(i))
            br_out(i) = periodic_field_value(k_modes, br_modes, r_out(i))
        end do

        call write_complex_profile(r_out, phi_out, n_output, &
            "/fields/Phi_periodic", &
            'Electrostatic potential perturbation Phi, periodic-Fourier &
            &solution of Poisson problem', 'statV')
        call write_complex_profile(r_out, br_out, n_output, &
            "/fields/Br_periodic", &
            'Radial magnetic field perturbation', 'G')

        print *, "...fourier_periodic solve finished with ", n_k, " modes."

    end subroutine run_fourier_periodic

end module rt_fourier_periodic_m
