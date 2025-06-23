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

        use species, only: init_deuterium_plasma, set_deuterium_plasma, plasma, interpolate_plasma_backs
        use IO_collection, only: create_output_directories

        implicit none
        class(electrostatic_t), intent(inout) :: this

        this%run_type = "electrostatic"
        print *, " ____  __.___   _____  "
        print *, "|    |/ _|   | /     \  "
        print *, "|      < |   |/  \ /  \ "
        print *, "|    |  \|   /    Y    \"
        print *, "|____|__ \___\____|__  /"
        print *, "        \/           \/ "
        print *, "electrostatic model initialized."

        call create_output_directories
        call generate_grids
        call init_deuterium_plasma(plasma)
        call set_deuterium_plasma(plasma)

    end subroutine

    subroutine run_electrostatic(this)

        use electrostatic_kernel, only: Krook_fill_kernel_phi, FP_fill_kernel_phi, kernel_spl_t
        use grid, only: xl_grid
        use IO_collection, only: write_matrix, write_complex_profile
        use poisson_solver, only: solve_poisson
        use config, only: output_path, collision_model
        use fields, only: EBdat, postprocess_electric_field

        implicit none

        class(electrostatic_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp

        call kernel_rho_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
        call kernel_rho_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)

        if (collision_model == "Krook") then
            call Krook_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)
        else if (collision_model == "FokkerPlanck") then
            call FP_fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)
        else
            stop "Error: collision model not recognized."
        end if

        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_re.dat", real(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)
        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_im.dat", dimag(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)

        allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b-1))

        EBdat%r_grid = xl_grid%xb
        
        call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
        call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_sol.dat")

        call postprocess_electric_field(EBdat)

    
    end subroutine

end module