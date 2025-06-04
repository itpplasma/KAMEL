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
        use grid, only: rg_grid

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
        call interpolate_plasma_backs(plasma, rg_grid%xb)

    end subroutine

    subroutine run_electrostatic(this)

        use KIM_kinds, only: dp
        use electrostatic_kernel, only: fill_kernel_phi, kernel_spl_t
        use grid, only: xl_grid
        use IO_collection, only: write_matrix, write_complex_profile
        use poisson_solver, only: solve_poisson
        use config, only: output_path
        use fields, only: EBdat, calculate_E_perp_psi, set_Br_constant, calculate_E_perp, get_Br_from_txt, calculate_E_from_phi, calculate_E_in_rsp_from_cyl
        use species, only: plasma

        implicit none
        class(electrostatic_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp

        call kernel_rho_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
        call kernel_rho_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)

        call fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)
        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_re.dat", real(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)
        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_im.dat", dimag(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)

        allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b-1))
        EBdat%r_grid = xl_grid%xb
        
        call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)
        call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_sol.dat")

        call calculate_E_perp_psi(plasma, EBdat)
        call write_complex_profile(xl_grid%xb, EBdat%E_perp_psi, xl_grid%npts_b, trim(output_path)//"/fields/E_perp_psi.dat")
        !call calculate_E_perp(EBdat)
        !call write_complex_profile(xl_grid%xb(1:xl_grid%npts_b-1), EBdat%E_perp, xl_grid%npts_b-1, trim(output_path)//"/fields/E_perp.dat")

        !call calculate_E_from_phi(EBdat)
        !call calculate_E_in_rsp_from_cyl(EBdat)

        !call write_complex_profile(EBdat%r_grid, EBdat%Er, size(EBdat%r_grid), trim(output_path)//"/fields/Er.dat")
        !call write_complex_profile(EBdat%r_grid, EBdat%Etheta, size(EBdat%r_grid), trim(output_path)//"/fields/Etheta.dat")
        !call write_complex_profile(EBdat%r_grid, EBdat%Ez, size(EBdat%r_grid), trim(output_path)//"/fields/Ez.dat")

        !call write_complex_profile(EBdat%r_grid, EBdat%Es, size(EBdat%r_grid), trim(output_path)//"/fields/Es.dat")
        !call write_complex_profile(EBdat%r_grid, EBdat%Ep, size(EBdat%r_grid), trim(output_path)//"/fields/Ep.dat")

    
    end subroutine

end module