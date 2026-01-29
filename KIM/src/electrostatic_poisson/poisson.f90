! run type for electrostatic model
! solves poisson's equation for electrostatic potential for given Br
module rt_electrostatic_m

    use kim_base_m, only: kim_t

    implicit none

    type, extends(kim_t) :: electrostatic_t
        contains
            procedure :: init => init_electrostatic
            procedure :: run => run_electrostatic
    end type electrostatic_t

    contains

    subroutine init_electrostatic(this)

        use species_m, only: plasma, set_plasma_quantities
        use IO_collection_m, only: create_output_directories
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use grid_m, only: rg_grid

        implicit none

        class(electrostatic_t), intent(inout) :: this

        this%run_type = "electrostatic"

        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)

        print *, "..."//trim(this%run_type)//" model initialized."

    end subroutine

    subroutine run_electrostatic(this)

        use kernel_m, only: Krook_fill_kernel_phi, FP_fill_kernels, fill_kernels_krook_fp, kernel_spl_t
        use kernel_adaptive_m, only: FP_fill_kernels_adaptive
        use grid_m, only: xl_grid
        use IO_collection_m, only: write_matrix, write_complex_profile, write_complex_profile_abs
        use poisson_solver_m, only: solve_poisson
        use config_m, only: output_path, collision_model, calculate_asymptotics
        use fields_m, only: EBdat, postprocess_electric_field,&
                            calculate_charge_density, calculate_current_density, calc_ideal_MA_phi
        use KIM_kinds_m, only: dp

        implicit none

        class(electrostatic_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp
        type(kernel_spl_t) :: kernel_j_phi_llp
        type(kernel_spl_t) :: kernel_j_B_llp

        complex(dp), allocatable :: rho(:)
        complex(dp), allocatable :: jpar(:)
        complex(dp), allocatable :: jpar_i(:)
        complex(dp), allocatable :: jpar_e(:)

        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        integer,dimension(8) :: values

        call date_and_time(date,time,zone,values)

        call kernel_rho_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
        call kernel_rho_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
        call kernel_j_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
        call kernel_j_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)

        write(*,*) "Start filling kernel at ", date, " ", time, " ..."

        if (collision_model == "Krook") then
            call run_Krook
        else if (collision_model == "FokkerPlanck") then
            call run_FP
        else if (collision_model == "Krook_FokkerPlanck") then
            call run_Krook_FP
            return
        else
            stop "Error: collision model not recognized. Options are: Krook, FokkerPlanck, Krook_FokkerPlanck."
        end if

    contains

        subroutine run_Krook
            call Krook_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)

            allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b),&
                    rho(xl_grid%npts_b))

            EBdat%r_grid = xl_grid%xb

            call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, "/fields/Phi", &
                'Electrostatic potential perturbation Phi, solution of Poisson problem', 'statV')

            call postprocess_electric_field(EBdat)

            call calculate_charge_density(rho, EBdat)
            call write_complex_profile_abs(xl_grid%xb, rho, xl_grid%npts_b, "/fields/rho", &
                'Charge density perturbation rho calculated from Poisson solution', 'statC/cm^3')
        end subroutine

        subroutine run_FP

            use grid_m, only: theta_integration
            use species_m, only: plasma
            use flr2_asymptotics_m, only: calc_flr2_asymptotic_Phi_MA, calc_hatK_Phi_in_Fourier
            use kernel_m, only: write_kernels
            use config_m, only: turn_off_electrons, turn_off_ions

            implicit none

            integer :: sp

            select case (trim(theta_integration))
                case ("GaussLegendre")
                    ! Fixed-order Gauss-Legendre path
                    call FP_fill_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, kernel_j_phi_llp, kernel_j_B_llp)
                case ("RKF45", "QUADPACK")
                    ! Adaptive path; actual integrator selected via theta_integration_method
                    call FP_fill_kernels_adaptive(kernel_rho_phi_llp, kernel_rho_B_llp, kernel_j_phi_llp, kernel_j_B_llp)
                case default
                    stop "Error: theta integration method not recognized."
            end select

            call write_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, kernel_j_phi_llp, kernel_j_B_llp)

            allocate(EBdat%Phi(xl_grid%npts_b), &
                    EBdat%Phi_i(xl_grid%npts_b), &
                    EBdat%Phi_e(xl_grid%npts_b), &
                    EBdat%Br(xl_grid%npts_b), &
                    EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), &
                    EBdat%E_perp(xl_grid%npts_b),&
                    rho(xl_grid%npts_b), &
                    jpar_i(xl_grid%npts_b), &
                    jpar_e(xl_grid%npts_b), &
                    jpar(xl_grid%npts_b))

            EBdat%r_grid = xl_grid%xb

            ! solve total Poisson problem and calculate total jpar
            call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, "/fields/Phi_m", &
                'Electrostatic potential perturbation Phi, solution of Poisson problem', 'statV')
            call calculate_current_density(jpar, EBdat%Phi, EBdat%Br, kernel_j_phi_llp%Kllp, kernel_j_B_llp%Kllp)
            call write_complex_profile_abs(xl_grid%xb, jpar, xl_grid%npts_b, "/fields/jpar", &
                'Parallel current density perturbation j_par calculated from Poisson solution', 'statA/cm^2')

            ! solve for electrons only, but use total phi in jpar calculation
            if (.not.turn_off_electrons) then
                call solve_poisson(kernel_rho_phi_llp%Kllp_e, kernel_rho_B_llp%Kllp_e, EBdat%Phi_e)
                call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_e, xl_grid%npts_b, "/fields/Phi_m_e", &
                    'Electrostatic potential perturbation Phi, solution of Poisson problem, only electrons', 'statV')
                call calculate_current_density(jpar_e, EBdat%Phi, EBdat%Br, kernel_j_phi_llp%Kllp_e, kernel_j_B_llp%Kllp_e)
                call write_complex_profile_abs(xl_grid%xb, jpar_e, xl_grid%npts_b, "/fields/jpar_e", &
                    'Parallel current density perturbation j_par calculated from Poisson solution', 'statA/cm^2')
            end if

            ! solve for ions only, but use total phi in jpar calculation
            do sp = 1, plasma%n_species - 1
                if (.not. turn_off_ions) then
                    call solve_poisson(kernel_rho_phi_llp%Kllp_i(:,:,sp), kernel_rho_B_llp%Kllp_i(:,:,sp), EBdat%Phi_i)
                    call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_i, xl_grid%npts_b, "/fields/Phi_m_"//trim(plasma%spec(sp)%name), &
                        'Electrostatic potential perturbation Phi, solution of Poisson problem for species '//trim(plasma%spec(sp)%name), 'statV')
                else
                    EBdat%Phi_i = (0.0d0, 0.0d0)
                    call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_i, xl_grid%npts_b, "/fields/Phi_m_"//trim(plasma%spec(sp)%name), &
                        'Electrostatic potential perturbation Phi, solution of Poisson problem for species '//trim(plasma%spec(sp)%name), 'statV')
                end if

                call calculate_current_density(jpar_i, EBdat%Phi, EBdat%Br, kernel_j_phi_llp%Kllp_i(:,:,sp), kernel_j_B_llp%Kllp_i(:,:,sp))
                call write_complex_profile_abs(xl_grid%xb, jpar_i, xl_grid%npts_b, "/fields/jpar_"//trim(plasma%spec(sp)%name), &
                    'Parallel current density perturbation j_par calculated from Poisson solution', 'statA/cm^2')
            end do


            call postprocess_electric_field(EBdat)

            call calculate_charge_density(rho, EBdat)
            call write_complex_profile_abs(xl_grid%xb, rho, xl_grid%npts_b, "/fields/rho", &
                'Charge density perturbation rho calculated from Poisson solution', 'statC/cm^3')

            if (calculate_asymptotics .eqv. .true.) then
                call calc_flr2_asymptotic_Phi_MA(plasma, EBdat)
                call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_MA_asymptotic, xl_grid%npts_b, &
                    "/fields/Phi_MA_asymptotic", 'Misalignment electrostatic potential perturbation in asymptotic limit', 'statV')

                call calc_ideal_MA_phi(EBdat, kernel_rho_phi_llp, kernel_rho_B_llp)
                call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_MA_ideal, xl_grid%npts_b, &
                    "/fields/Phi_MA_ideal", 'Electrostatic potential perturbation in ideal limit where E_perp_MA = 0', 'statV')

                call calc_hatK_Phi_in_Fourier(plasma)
            end if

        end subroutine

        subroutine run_Krook_FP
            call Krook_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)

            allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b),&
                    rho(xl_grid%npts_b))

            EBdat%r_grid = xl_grid%xb

            call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, "/fields/Phi_"//trim(collision_model), &
                'Electrostatic potential perturbation Phi, solution of Poisson problem', 'statV')

            call postprocess_electric_field(EBdat)

            call calculate_charge_density(rho, EBdat)
            call write_complex_profile_abs(xl_grid%xb, rho, xl_grid%npts_b, "/fields/rho_"//trim(collision_model), &
                'Charge density perturbation rho calculated from Poisson solution', 'statC/cm^3')
        end subroutine

    end subroutine

end module
