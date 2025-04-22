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
        use plotting, only: write_profile, plot_1D
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

        call generate_grids
        call init_deuterium_plasma(plasma)
        call set_deuterium_plasma(plasma)
        call interpolate_plasma_backs(plasma, rg_grid%xb)

    end subroutine

    subroutine run_reduced(this)

        use gauss_quad
        use KIM_kinds, only: dp
        use reduced_kernel, only: fill_kernel_phi, kernel_spl_t
        use grid, only: xl_grid
        use plotting, only: write_matrix, plot_matrix, write_complex_profile, plot_complex_1D
        use poisson_solver, only: solve_poisson
        use species, only: plasma

        implicit none
        class(reduced_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp
        integer :: npts_l, npts_lp
        complex(dp), allocatable :: phi_sol(:)
        
        allocate(phi_sol(xl_grid%npts_b))

        npts_l = xl_grid%npts_b
        npts_lp = npts_l

        print *, plasma%spec(0)%lambda_D(1)
        print *, plasma%spec(1)%lambda_D(1)

        call kernel_rho_phi_llp%init_kernel(npts_l, npts_lp)
        call kernel_rho_B_llp%init_kernel(npts_l, npts_lp)

        call fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)
        call write_matrix("kernel_phi_llp.txt", real(kernel_rho_phi_llp%Kllp), npts_l, npts_lp)
        !call plot_matrix("kernel_phi_llp.txt")

        call solve_poisson(kernel_rho_phi_llp%Kllp, kernel_rho_B_llp%Kllp, phi_sol)

        call write_complex_profile(xl_grid%xb, phi_sol, xl_grid%npts_b, "phi_sol.txt")
    
    end subroutine



end module