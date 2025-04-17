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

        use species, only: init_deuterium_plasma, set_deuterium_plasma, plasma
        use plotting, only: write_profile, plot_1D
        use plasma_parameter, only: r_prof, iprof_length

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

    end subroutine

    subroutine run_reduced(this)

        use gauss_quad
        use KIM_kinds, only: dp
        use reduced_kernel, only: fill_kernel_phi, kernel_spl_t
        use grid, only: xl_grid, rg_grid
        use plotting, only: write_matrix, plot_matrix, write_complex_profile, plot_complex_1D
        use poisson_solver, only: solve_poisson
        use species, only: plasma

        implicit none
        class(reduced_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_phi_llp
        integer :: npts_l, npts_lp
        complex(dp), allocatable :: phi_sol(:)
        
        allocate(phi_sol(rg_grid%npts_b))

        npts_l = xl_grid%npts_b
        npts_lp = npts_l

        print *, plasma%spec(0)%lambda_D(1)
        print *, plasma%spec(1)%lambda_D(1)

        call kernel_phi_llp%init_kernel(npts_l, npts_lp)

        call fill_kernel_phi(kernel_phi_llp)
        call write_matrix("kernel_phi_llp.txt", real(kernel_phi_llp%Kllp), npts_l, npts_lp)
        !call plot_matrix("kernel_phi_llp.txt")

        call solve_poisson(kernel_phi_llp%Kllp, phi_sol)

        call write_complex_profile(rg_grid%xb, phi_sol, rg_grid%npts_b, "phi_sol.txt")
    
    end subroutine



end module