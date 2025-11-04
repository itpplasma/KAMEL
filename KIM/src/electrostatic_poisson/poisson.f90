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

        use species_m, only: init_plasma, plasma, set_plasma_quantities
        use IO_collection_m, only: create_output_directories
        use equilibrium_m, only: calculate_equil, interpolate_equil, write_equil
        use grid_m, only: rg_grid

        implicit none

        class(electrostatic_t), intent(inout) :: this

        this%run_type = "electrostatic"

        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)
        call write_equil

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

            implicit none

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

            call write_matrix("kernel/K_rho_phi", real(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b, &
                'Complex FLR2 benchmark kernel K_rho_phi', '1/cm^2')
            ! call write_matrix("kernel/K_rho_phi_im.dat", dimag(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)
            call write_matrix("kernel/K_rho_B.dat", real(kernel_rho_B_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b, &
                'Complex FLR2 benchmark kernel K_rho_B', '1/cm^2')
            ! call write_matrix("kernel/K_rho_B_im.dat", dimag(kernel_rho_B_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)

            allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b),&
                    rho(xl_grid%npts_b), jpar(xl_grid%npts_b))

            EBdat%r_grid = xl_grid%xb
            
            call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, "/fields/Phi", &
                'Electrostatic potential perturbation Phi, solution of Poisson problem', 'statV')

            call postprocess_electric_field(EBdat)

            call calculate_charge_density(rho, EBdat)
            call calculate_current_density(jpar, EBdat, kernel_j_phi_llp, kernel_j_B_llp)
            call write_complex_profile_abs(xl_grid%xb, rho, xl_grid%npts_b, "/fields/rho", &
                'Charge density perturbation rho calculated from Poisson solution', 'statC/cm^3')
            call write_complex_profile_abs(xl_grid%xb, jpar, xl_grid%npts_b, "/fields/jpar", &
                'Parallel current density perturbation j_par calculated from Poisson solution', 'statA/cm^2')

            if (calculate_asymptotics .eqv. .true.) then
                call calc_flr2_asymptotic_Phi_MA(plasma, EBdat)
                call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_MA_asymptotic, xl_grid%npts_b, &
                    "/fields/phi_MA_asymptotic", 'Misalignment electrostatic potential perturbation in asymptotic limit', 'statV')

                call calc_ideal_MA_phi(EBdat, kernel_rho_phi_llp, kernel_rho_B_llp)
                call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_MA_ideal, xl_grid%npts_b, &
                    "/fields/phi_MA_ideal", 'Electrostatic potential perturbation in ideal limit where E_perp_MA = 0', 'statV')

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
