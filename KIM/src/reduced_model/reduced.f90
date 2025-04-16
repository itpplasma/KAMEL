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

    end subroutine

    subroutine run_reduced(this)

        use gauss_quad
        use KIM_kinds, only: dp
        use reduced_kernel, only: fill_kernel_phi, kernel_spl_t
        use grid, only: xl_grid
        use plotting, only: write_profile, plot_1D
        use species, only: plasma_t, init_deuterium_plasma, set_deuterium_plasma
        use plasma_parameter, only: r_prof, iprof_length

        implicit none
        class(reduced_t), intent(inout) :: this
        type(kernel_spl_t) :: kernel_phi_llp
        integer :: npts_l, npts_lp
        type(plasma_t) :: plasma

        npts_l = xl_grid%npts_b
        npts_lp = npts_l

        call kernel_phi_llp%init_kernel(npts_l, npts_lp)

        call fill_kernel_phi(kernel_phi_llp)
        call write_matrix("kernel_phi_llp.txt", real(kernel_phi_llp%Kllp), npts_l, npts_lp)
        call plot_matrix("kernel_phi_llp.txt")
    
    end subroutine



end module