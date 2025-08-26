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
        use equilibrium_m, only: calculate_equil

        implicit none

        class(electrostatic_t), intent(inout) :: this

        this%run_type = "electrostatic"

        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)

        print *, "... electrostatic model initialized."

    end subroutine

    subroutine run_electrostatic(this)

        use electrostatic_kernel_m, only: Krook_fill_kernel_phi, FP_fill_kernels, fill_kernels_krook_fp, kernel_spl_t
        use electrostatic_kernel_adaptive_mod, only: FP_fill_kernels_adaptive
        use grid_m, only: xl_grid
        use IO_collection_m, only: write_matrix, write_complex_profile, write_complex_profile_abs
        use poisson_solver_m, only: solve_poisson
        use config_m, only: output_path, collision_model
        use fields_m, only: EBdat, postprocess_electric_field, postprocess_electric_field_with_model,&
                            calculate_charge_density, calculate_current_density, calc_ideal_MA_phi
        use KIM_kinds_m, only: dp

        implicit none

        class(electrostatic_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp
        type(kernel_spl_t) :: kernel_j_phi_llp
        type(kernel_spl_t) :: kernel_j_B_llp
        type(kernel_spl_t) :: kernel_krook_rho_phi_llp
        type(kernel_spl_t) :: kernel_krook_rho_B_llp
        type(kernel_spl_t) :: kernel_fp_rho_phi_llp
        type(kernel_spl_t) :: kernel_fp_rho_B_llp

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
            stop "Error: collision model not recognized."
        end if

    contains

        subroutine run_Krook
            call Krook_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)

            allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b),&
                    rho(xl_grid%npts_b))

            EBdat%r_grid = xl_grid%xb
            
            call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_"//trim(collision_model)//".dat")

            call postprocess_electric_field(EBdat)

            call calculate_charge_density(rho, EBdat)
            call write_complex_profile_abs(xl_grid%xb, rho, xl_grid%npts_b, trim(output_path)//"/fields/rho_"//trim(collision_model)//".dat")
        end subroutine

        subroutine run_FP

            use grid_m, only: theta_integration
            use species_m, only: plasma
            use flr2_asymptotics_m, only: calc_flr2_asymptotic_Phi_MA

            implicit none

            if (trim(theta_integration) == "RKF45") then
                call FP_fill_kernels_adaptive(kernel_rho_phi_llp, kernel_rho_B_llp, kernel_j_phi_llp, kernel_j_B_llp)
            else if (trim(theta_integration) == "GaussLegendre") then
                call FP_fill_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, kernel_j_phi_llp, kernel_j_B_llp)
            else
                stop "Error: theta integration method not recognized."
            end if

            allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b),&
                    rho(xl_grid%npts_b), jpar(xl_grid%npts_b))

            EBdat%r_grid = xl_grid%xb
            
            call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_"//trim(collision_model)//".dat")

            call postprocess_electric_field(EBdat)

            call calculate_charge_density(rho, EBdat)
            call calculate_current_density(jpar, EBdat, kernel_j_phi_llp, kernel_j_B_llp)
            call write_complex_profile_abs(xl_grid%xb, rho, xl_grid%npts_b, trim(output_path)//"/fields/rho_"//trim(collision_model)//".dat")
            call write_complex_profile_abs(xl_grid%xb, jpar, xl_grid%npts_b, trim(output_path)//"/fields/jpar_"//trim(collision_model)//".dat")

            call calc_flr2_asymptotic_Phi_MA(plasma, EBdat)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_MA_asymptotic, xl_grid%npts_b, trim(output_path)//"/fields/phi_MA_asymptotic_"//trim(collision_model)//".dat")

            call calc_ideal_MA_phi(EBdat, kernel_rho_phi_llp, kernel_rho_B_llp)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_MA_ideal, xl_grid%npts_b, trim(output_path)//"/fields/phi_MA_ideal_"//trim(collision_model)//".dat")

        end subroutine
        
        subroutine run_Krook_FP
            call Krook_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)

            allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b),&
                    rho(xl_grid%npts_b))

            EBdat%r_grid = xl_grid%xb
            
            call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_"//trim(collision_model)//".dat")

            call postprocess_electric_field(EBdat)

            call calculate_charge_density(rho, EBdat)
            call write_complex_profile_abs(xl_grid%xb, rho, xl_grid%npts_b, trim(output_path)//"/fields/rho_"//trim(collision_model)//".dat")
        end subroutine

    end subroutine

end module
