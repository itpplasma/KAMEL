module rt_reduced

    use kim_base, only: kim_t

    implicit none

    type, extends(kim_t) :: reduced_t
        contains
            procedure :: init => init_reduced
            procedure :: run => run_reduced
    end type reduced_t

    contains

    subroutine init_reduced(this)

        use species, only: init_deuterium_plasma, set_deuterium_plasma, plasma, interpolate_plasma_backs
        use IO_collection, only: create_output_directories
        use grid, only: rg_grid

        implicit none
        class(reduced_t), intent(inout) :: this

        this%run_type = "reduced"
        print *, " ____  __.___   _____  "
        print *, "|    |/ _|   | /     \  "
        print *, "|      < |   |/  \ /  \ "
        print *, "|    |  \|   /    Y    \"
        print *, "|____|__ \___\____|__  /"
        print *, "        \/           \/ "
        print *, "Reduced model initialized."

        call create_output_directories
        call generate_grids
        call init_deuterium_plasma(plasma)
        call set_deuterium_plasma(plasma)
        call interpolate_plasma_backs(plasma, rg_grid%xb)

    end subroutine

    subroutine run_reduced(this)

        use KIM_kinds, only: dp
        use reduced_kernel, only: fill_kernel_phi, kernel_spl_t
        use grid, only: xl_grid
        use IO_collection, only: write_matrix, write_complex_profile
        use poisson_solver, only: solve_poisson
        use config, only: output_path
        use fields, only: EBdat, calculate_E_perp_psi, set_Br_constant, calculate_E_perp, get_Br_from_txt
        use species, only: plasma

        implicit none
        class(reduced_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp
        !complex(dp) :: Br_const
        character(len=256) :: file_path

        call kernel_rho_phi_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
        call kernel_rho_B_llp%init_kernel(xl_grid%npts_b, xl_grid%npts_b)

        call fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)
        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_re.dat", real(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)
        call write_matrix(trim(output_path)//"kernel/kernel_phi_llp_im.dat", dimag(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b)

        allocate(EBdat%Phi(xl_grid%npts_b), EBdat%Br(xl_grid%npts_b), EBdat%E_perp_psi(xl_grid%npts_b), &
                EBdat%r_grid(xl_grid%npts_b), EBdat%E_perp(xl_grid%npts_b-1))
        EBdat%r_grid = xl_grid%xb
        !Br_const = 1.0d0
        !call set_Br_constant(EBdat, Br_const)
        file_path = './inp/Br_in.dat'
        call get_Br_from_txt(EBdat, file_path)
        call write_complex_profile(xl_grid%xb, EBdat%Br, xl_grid%npts_b, trim(output_path)//"/fields/Br.dat")

        call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, EBdat%Phi)

        call write_complex_profile(xl_grid%xb, EBdat%Phi, xl_grid%npts_b, trim(output_path)//"/fields/phi_sol.dat")

        call calculate_E_perp_psi(plasma, EBdat)
        call write_complex_profile(xl_grid%xb, EBdat%E_perp_psi, xl_grid%npts_b, trim(output_path)//"/fields/E_perp_psi.dat")
        call calculate_E_perp(EBdat)
        call write_complex_profile(xl_grid%xb, EBdat%E_perp, xl_grid%npts_b, trim(output_path)//"/fields/E_r.dat")
    
    end subroutine

end module