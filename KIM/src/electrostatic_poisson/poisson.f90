module rt_electrostatic

    use kim_base, only: kim_t

    implicit none

    type, extends(kim_t) :: electrostatic_t
        contains
            procedure :: init => init_electrostatic
            procedure :: run => run_electrostatic
    end type electrostatic_t

    contains

    subroutine init_electrostatic(this)

        use species, only: init_plasma, plasma, set_plasma_quantities
        use IO_collection, only: create_output_directories
        use equilibrium, only: calculate_equil

        implicit none

        class(electrostatic_t), intent(inout) :: this

        this%run_type = "electrostatic"
        print *, " ____  __.___   _____  "
        print *, "|    |/ _|   | /     \  "
        print *, "|      < |   |/  \ /  \ "
        print *, "|    |  \|   /    Y    \"
        print *, "|____|__ \___\____|__  /"
        print *, "        \/           \/ "

        call create_output_directories
        call generate_grids
        !call init_plasma(plasma)
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)

        print *, "... electrostatic model initialized."

    end subroutine

    subroutine run_electrostatic(this)

        use electrostatic_kernel, only: Krook_fill_kernel_phi, FP_fill_kernel_phi, fill_kernels_krook_fp, kernel_spl_t
        use grid, only: xl_grid
        use IO_collection, only: write_matrix, write_complex_profile
        use poisson_solver, only: solve_poisson
        use config, only: output_path, collision_model
        use fields, only: EBdat, postprocess_electric_field, postprocess_electric_field_with_model

        implicit none

        class(electrostatic_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp
        type(kernel_spl_t) :: kernel_krook_rho_phi_llp
        type(kernel_spl_t) :: kernel_krook_rho_B_llp
        type(kernel_spl_t) :: kernel_fp_rho_phi_llp
        type(kernel_spl_t) :: kernel_fp_rho_B_llp

        call kernel_rho_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
        call kernel_rho_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)

        if (collision_model == "Krook") then
            call Krook_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)
        else if (collision_model == "FokkerPlanck") then
            call FP_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)
        else if (collision_model == "Krook_FokkerPlanck") then
            ! Initialize kernels for both models
            call kernel_krook_rho_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
            call kernel_krook_rho_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
            call kernel_fp_rho_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
            call kernel_fp_rho_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
            
            ! Fill both kernels using unified subroutine
            call fill_kernels_krook_fp(kernel_krook_rho_phi_llp, kernel_krook_rho_B_llp, &
                                      kernel_fp_rho_phi_llp, kernel_fp_rho_B_llp)
            
            ! Allocate EBdat fields
            allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                    EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b))
            EBdat%r_grid = xl_grid%xb
            
            ! Solve and write Krook solution
            call solve_poisson(kernel_krook_rho_phi_llp%Kllp, kernel_krook_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_Krook_sol.dat")
            call postprocess_electric_field_with_model(EBdat, "Krook")
            
            ! Solve and write Fokker-Planck solution
            call solve_poisson(kernel_fp_rho_phi_llp%Kllp, kernel_fp_rho_B_llp%Kllp, EBdat%Phi)
            call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_FokkerPlanck_sol.dat")
            call postprocess_electric_field_with_model(EBdat, "FokkerPlanck")
            
            return
        else
            stop "Error: collision model not recognized."
        end if

        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_re.dat", real(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)
        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_im.dat", dimag(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)

        allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b))

        EBdat%r_grid = xl_grid%xb
        
        call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
        
        ! Write phi solution with appropriate suffix
        if (collision_model == "Krook") then
            call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_krook_sol.dat")
        else if (collision_model == "FokkerPlanck") then
            call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_fp_sol.dat")
        else
            call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_sol.dat")
        end if

        call postprocess_electric_field(EBdat)
    
    end subroutine

end module