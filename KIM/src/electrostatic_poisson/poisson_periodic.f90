! run type for the electrostatic forced-periodicity model.
!
! On a given (m, n) case this run-type locates the resonant surface rm, sizes a
! periodic window from a representative Larmor radius rho_L(rm), builds the
! periodic plasma on that window (Phase-2.3), assembles the dense Fourier
! matrices K^{rhoPhi}/K^{rhoB} over one period (Phase-2.4), solves the periodic
! electrostatic system for the Fourier coefficients Phi_m and inverse-DFTs them
! into delta_Phi(r) on the window grid (Phase-2.5), then packs delta_Phi into
! fields_m::EBdat.
module rt_electrostatic_periodic_m

    use kim_base_m, only: kim_t

    implicit none

    type, extends(kim_t) :: electrostatic_periodic_t
        contains
            procedure :: init => init_electrostatic_periodic
            procedure :: run => run_electrostatic_periodic
    end type electrostatic_periodic_t

    contains

    !> Global setup, identical to the electrostatic run-type's init: build the
    !> grids and equilibrium and populate the GLOBAL plasma. run() reads the
    !> resonance rm and the representative rho_L(rm) from that global plasma
    !> before redirecting the grid onto the periodic window.
    subroutine init_electrostatic_periodic(this)

        use species_m, only: plasma, set_plasma_quantities
        use IO_collection_m, only: create_output_directories
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use grid_m, only: rg_grid
        use config_m, only: output_path, hdf5_output

        implicit none

        class(electrostatic_periodic_t), intent(inout) :: this
        logical :: ex

        this%run_type = "electrostatic_periodic"

        call create_output_directories

        ! Create setup/ directory for text output of setup parameters
        if (.not. hdf5_output) then
            inquire(file=trim(output_path)//'setup', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'setup')
            end if
        end if
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)

        print *, "..."//trim(this%run_type)//" model initialized."

    end subroutine init_electrostatic_periodic

    !> Locate the resonance, size and build the periodic window plasma, assemble
    !> the Fourier matrices, solve, reconstruct delta_Phi on the window grid, and
    !> pack it into fields_m::EBdat.
    subroutine run_electrostatic_periodic(this)

        use KIM_kinds_m, only: dp
        use constants_m, only: pi
        use config_m, only: periodic_dr_asis_scale, periodic_dr_tr_scale, &
                            periodic_kmax_scale, periodic_n_rg, hdf5_output
        use setup_m, only: Br_boundary_re, Br_boundary_im
        use species_m, only: plasma
        use grid_m, only: rg_grid
        use kim_resonances_m, only: r_res
        use periodic_background_m, only: build_periodic_plasma
        use periodic_assembly_m, only: assemble_periodic_matrices
        use periodic_solve_m, only: solve_periodic, reconstruct_delta_phi
        use fields_m, only: EBdat
        use IO_collection_m, only: write_complex_profile_abs

        implicit none

        class(electrostatic_periodic_t), intent(inout) :: this

        complex(dp), allocatable :: Kphi(:,:), KB(:,:), Phi_m(:), dPhi(:)
        complex(dp) :: Br_const
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr, L, k_max
        integer :: M, n_rg, info

        ! 1. Locate the resonant surface rm = r_res (q = |m/n|) on the global plasma.
        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, "Error (electrostatic_periodic): no resonance found, r_res = ", r_res
            error stop "electrostatic_periodic: resonance not found"
        end if
        rm = r_res

        ! 2. Representative Larmor radius rho_L(rm) (electrons) on the global
        ! plasma via 4-point Lagrange interpolation.
        rhoL_rm = interp_rho_L(rm)
        if (.not. (rhoL_rm > 0.0_dp)) then
            print *, "Error (electrostatic_periodic): rho_L(rm) not positive, = ", rhoL_rm
            error stop "electrostatic_periodic: invalid rho_L(rm)"
        end if

        ! 3. Window geometry and Fourier cutoff from rho_L(rm).
        dx_asis = periodic_dr_asis_scale * rhoL_rm
        dx_tr   = periodic_dr_tr_scale * rhoL_rm
        L       = 2.0_dp * (dx_asis + dx_tr)
        k_max   = periodic_kmax_scale / rhoL_rm
        M       = ceiling(k_max * L / (2.0_dp * pi))
        n_rg    = periodic_n_rg

        print *, "electrostatic_periodic: rm = ", rm, " rho_L(rm) = ", rhoL_rm
        print *, "electrostatic_periodic: L = ", L, " M = ", M, " n_rg = ", n_rg

        ! 4. Build the periodic window plasma (redirects rg_grid + plasma).
        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)

        ! 5. Assemble the dense periodic Fourier matrices over one period.
        call assemble_periodic_matrices(plasma, L, M, Kphi, KB)

        ! 6. Solve for the Fourier coefficients Phi_m under the constant Br drive.
        Br_const = cmplx(Br_boundary_re, Br_boundary_im, dp)
        call solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info)
        if (info /= 0) then
            print *, "Error (electrostatic_periodic): solve_periodic failed, info = ", info
            error stop "electrostatic_periodic: periodic solve failed"
        end if

        ! 7. Reconstruct delta_Phi on the window grid via inverse DFT.
        dPhi = reconstruct_delta_phi(Phi_m, L, M, rg_grid%xb)

        ! 8. Pack into EBdat (window grid + reconstructed potential).
        if (allocated(EBdat%r_grid)) deallocate(EBdat%r_grid)
        if (allocated(EBdat%Phi))    deallocate(EBdat%Phi)
        EBdat%r_grid = rg_grid%xb
        EBdat%Phi    = dPhi

        if (hdf5_output) then
            call write_complex_profile_abs(rg_grid%xb, EBdat%Phi, rg_grid%npts_b, &
                "/fields/Phi", &
                'Electrostatic potential perturbation Phi, forced-periodicity solution', &
                'statV')
        end if

    contains

        !> 4-point Lagrange interpolation of the global electron Larmor radius
        !> plasma%spec(0)%rho_L at radius x, using the same binsrc + plag_coeff
        !> stencil as the rest of KIM.
        real(dp) function interp_rho_L(x) result(rhoLx)
            real(dp), intent(in) :: x
            integer, parameter :: nlagr = 4, nder = 0
            real(dp) :: coef(0:nder, nlagr)
            integer :: gs, ir, ibeg, iend

            gs = plasma%grid_size
            call binsrc(plasma%r_grid, 1, gs, x, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend > gs) then
                iend = gs
                ibeg = iend - nlagr + 1
            end if
            call plag_coeff(nlagr, nder, x, plasma%r_grid(ibeg:iend), coef)
            rhoLx = sum(coef(0, :) * plasma%spec(0)%rho_L(ibeg:iend))
        end function interp_rho_L

    end subroutine run_electrostatic_periodic

end module rt_electrostatic_periodic_m
