! run type for the electrostatic forced-periodicity model.
!
! On a given (m, n) case this run-type locates the resonant surface rm, sizes a
! periodic window from a representative Larmor radius rho_L(rm), builds the
! periodic plasma on that window (Phase-2.3), assembles the dense Fourier
! matrices K^{rhoPhi}/K^{rhoB} over one period, solves for the periodic
! deviation from a physical aligned-potential lift, reconstructs delta_Phi(r)
! on the window grid, then packs it into fields_m::EBdat.
module rt_electrostatic_periodic_m

    use kim_base_m, only: kim_t
    use KIM_kinds_m, only: dp

    implicit none

    type, extends(kim_t) :: electrostatic_periodic_t
        contains
            procedure :: init => init_electrostatic_periodic
            procedure :: run => run_electrostatic_periodic
    end type electrostatic_periodic_t

    integer, parameter, public :: PERIODIC_SCALE_OK = 0
    integer, parameter, public :: PERIODIC_SCALE_NO_ACTIVE = 1
    integer, parameter, public :: PERIODIC_SCALE_INVALID_RHO = 2
    real(dp), parameter :: PERIODIC_LIFT_PROJECTION_TOLERANCE = 5.0e-2_dp

    public :: compute_periodic_delta_phi, select_periodic_reference_scale

    contains

    !> Reusable periodic core: build the window plasma, assemble the Fourier
    !> matrices, solve for deviation coefficients under a constant Br drive, and
    !> reconstruct physical delta_Phi on the EXPLICIT output grid r_out (not the window
    !> grid rg_grid%xb). All window sizing is passed in EXPLICITLY so both the
    !> run-type and the Phase-3 convergence harness can evaluate delta_Phi on a
    !> fixed diagnostic grid independent of the (n_rg, M, dx) discretization.
    !>
    !> Does NOT error stop on a singular solve: it returns info /= 0 so the
    !> caller can decide (the run-type error stops; the tests inspect info).
    subroutine compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, Br_const, &
                                          r_out, dPhi, info, lift_projection_residual)
        use KIM_kinds_m, only: dp
        use species_m, only: plasma
        use grid_m, only: rg_grid
        use periodic_background_m, only: build_periodic_plasma
        use periodic_assembly_m, only: assemble_periodic_matrices
        use periodic_solve_m, only: solve_periodic_deviation, &
            reconstruct_delta_phi_from_lift, build_physical_c2_lift, &
            PERIODIC_SOLVE_INVALID_INPUT

        real(dp),    intent(in)  :: rm, dx_asis, dx_tr
        integer,     intent(in)  :: M, n_rg
        complex(dp), intent(in)  :: Br_const
        real(dp),    intent(in)  :: r_out(:)
        complex(dp), allocatable, intent(out) :: dPhi(:)
        integer,     intent(out) :: info
        real(dp), intent(out), optional :: lift_projection_residual

        complex(dp), allocatable :: Kphi(:,:), KB(:,:), psi_m(:), Br_m(:)
        complex(dp), allocatable :: Phi_ref_grid(:), d2Phi_ref_grid(:)
        complex(dp), allocatable :: Phi_ref_out(:), d2Phi_ref_out(:)
        logical, allocatable :: projection_mask(:)
        complex(dp) :: phi_left, phi_right
        real(dp) :: L, projection_residual
        integer :: dim

        if (present(lift_projection_residual)) lift_projection_residual = huge(1.0_dp)
        if (dx_asis <= 0.0_dp .or. dx_tr <= 0.0_dp .or. M < 0 .or. n_rg <= 0) then
            info = PERIODIC_SOLVE_INVALID_INPUT
            return
        end if
        L = 2.0_dp * (dx_asis + dx_tr)

        ! Capture the two physical ideal-MHD edge values before periodization
        ! replaces the outer background by artificial periodic transitions.
        call aligned_lift_anchor(rm - 0.5_dp * L, Br_const, phi_left, info)
        if (info /= 0) return
        call aligned_lift_anchor(rm + 0.5_dp * L, Br_const, phi_right, info)
        if (info /= 0) return
        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)
        call assemble_periodic_matrices(plasma, L, M, Kphi, KB)

        allocate(Phi_ref_grid(size(rg_grid%xb)), d2Phi_ref_grid(size(rg_grid%xb)))
        call build_physical_c2_lift(rg_grid%xb, rm, dx_asis, dx_tr, phi_left, phi_right, &
                                    Phi_ref_grid, d2Phi_ref_grid)
        allocate(projection_mask(size(rg_grid%xb)))
        projection_mask = abs(rg_grid%xb - rm) <= dx_asis
        dim = 2 * M + 1
        allocate(Br_m(dim))
        Br_m = (0.0_dp, 0.0_dp)
        Br_m(M + 1) = Br_const
        call solve_periodic_deviation(Kphi, KB, L, M, Br_m, rg_grid%xb, &
            Phi_ref_grid, d2Phi_ref_grid, PERIODIC_LIFT_PROJECTION_TOLERANCE, &
            psi_m, projection_residual, info, projection_mask=projection_mask)
        if (present(lift_projection_residual)) lift_projection_residual = projection_residual
        if (info /= 0) return

        allocate(Phi_ref_out(size(r_out)), d2Phi_ref_out(size(r_out)))
        call build_physical_c2_lift(r_out, rm, dx_asis, dx_tr, phi_left, phi_right, &
                                    Phi_ref_out, d2Phi_ref_out)
        dPhi = reconstruct_delta_phi_from_lift(psi_m, L, M, r_out, Phi_ref_out)
    end subroutine compute_periodic_delta_phi

    subroutine aligned_lift_anchor(x, Br_const, phi_anchor, info)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit
        use periodic_background_m, only: kim_aperfuns, kim_geom_aperfuns
        use periodic_solve_m, only: PERIODIC_SOLVE_OK, PERIODIC_SOLVE_INVALID_INPUT

        real(dp), intent(in) :: x
        complex(dp), intent(in) :: Br_const
        complex(dp), intent(out) :: phi_anchor
        integer, intent(out) :: info

        real(dp) :: primitives(5), geometry(4), Er_x, kp_x, B0_x

        phi_anchor = (0.0_dp, 0.0_dp)
        info = PERIODIC_SOLVE_INVALID_INPUT
        ! Both callbacks sample the cached true background, so repeated
        ! convergence solves never derive anchors from an already-periodized
        ! live plasma.
        call kim_aperfuns(size(primitives), x, primitives)
        call kim_geom_aperfuns(size(geometry), x, geometry)
        Er_x = primitives(5)
        B0_x = geometry(1)
        kp_x = geometry(3)
        if (.not. ieee_is_finite(Er_x) .or. .not. ieee_is_finite(kp_x) .or. &
            .not. ieee_is_finite(B0_x) .or. abs(B0_x) <= tiny(1.0_dp) .or. &
            abs(kp_x) <= sqrt(epsilon(1.0_dp))) return

        ! Global bc_type=3 imposes Phi=-phi_aligned with
        ! phi_aligned=+i Er Br/(B0 k_parallel).
        phi_anchor = -com_unit * Er_x * Br_const / (B0_x * kp_x)
        info = PERIODIC_SOLVE_OK
    end subroutine aligned_lift_anchor

    !> Select the largest Larmor radius among species active in the FP kernel.
    subroutine select_periodic_reference_scale(plasma_in, x, electrons_active, ions_active, &
                                               rho_ref, reference_species, info)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use KIM_kinds_m, only: dp
        use species_m, only: plasma_t

        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: x
        logical, intent(in) :: electrons_active, ions_active
        real(dp), intent(out) :: rho_ref
        integer, intent(out) :: reference_species, info

        real(dp) :: rho_sp
        integer :: sp

        rho_ref = 0.0_dp
        reference_species = -1
        info = PERIODIC_SCALE_NO_ACTIVE
        if (plasma_in%n_species <= 0) return
        if (.not. allocated(plasma_in%spec) .or. .not. allocated(plasma_in%r_grid) &
            .or. plasma_in%grid_size < 4) then
            info = PERIODIC_SCALE_INVALID_RHO
            return
        end if

        do sp = 0, plasma_in%n_species - 1
            if (.not. electrons_active .and. sp == 0) cycle
            if (.not. ions_active .and. sp >= 1) cycle
            if (.not. allocated(plasma_in%spec(sp)%rho_L) &
                .or. size(plasma_in%spec(sp)%rho_L) < plasma_in%grid_size) then
                reference_species = sp
                info = PERIODIC_SCALE_INVALID_RHO
                return
            end if
            rho_sp = interp_species_rho_L(plasma_in, sp, x)
            if (.not. ieee_is_finite(rho_sp) .or. rho_sp <= 0.0_dp) then
                rho_ref = rho_sp
                reference_species = sp
                info = PERIODIC_SCALE_INVALID_RHO
                return
            end if
            if (rho_sp > rho_ref) then
                rho_ref = rho_sp
                reference_species = sp
            end if
        end do

        if (reference_species >= 0) info = PERIODIC_SCALE_OK
    end subroutine select_periodic_reference_scale

    real(dp) function interp_species_rho_L(plasma_in, sp, x) result(rhoLx)
        use KIM_kinds_m, only: dp
        use species_m, only: plasma_t

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp
        real(dp), intent(in) :: x
        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: gs, ir, ibeg, iend

        gs = plasma_in%grid_size
        call binsrc(plasma_in%r_grid, 1, gs, x, ir)
        ibeg = max(1, ir - nlagr / 2)
        iend = ibeg + nlagr - 1
        if (iend > gs) then
            iend = gs
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, x, plasma_in%r_grid(ibeg:iend), coef)
        rhoLx = sum(coef(0, :) * plasma_in%spec(sp)%rho_L(ibeg:iend))
    end function interp_species_rho_L

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
                            periodic_kmax_scale, periodic_n_rg, hdf5_output, &
                            turn_off_electrons, turn_off_ions
        use setup_m, only: Br_boundary_re, Br_boundary_im
        use species_m, only: plasma
        use grid_m, only: rg_grid
        use kim_resonances_m, only: r_res
        use fields_m, only: EBdat
        use IO_collection_m, only: write_complex_profile_abs, write_periodic_scale_metadata

        implicit none

        class(electrostatic_periodic_t), intent(inout) :: this

        complex(dp), allocatable :: dPhi(:)
        complex(dp) :: Br_const
        real(dp), allocatable :: r_win(:)
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr, L, k_max, lift_projection_residual
        integer :: M, n_rg, info, i, reference_species

        ! 1. Locate the resonant surface rm = r_res (q = |m/n|) on the global plasma.
        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, "Error (electrostatic_periodic): no resonance found, r_res = ", r_res
            error stop "electrostatic_periodic: resonance not found"
        end if
        rm = r_res

        ! 2. Select the largest Larmor radius among the species that the FP
        ! kernel will actually assemble. This keeps the full response of every
        ! active species resolved and is independent of species storage order.
        call select_periodic_reference_scale(plasma, rm, .not. turn_off_electrons, &
                                             .not. turn_off_ions, rhoL_rm, &
                                             reference_species, info)
        select case (info)
        case (PERIODIC_SCALE_NO_ACTIVE)
            error stop "electrostatic_periodic: no active kinetic species"
        case (PERIODIC_SCALE_INVALID_RHO)
            print *, "Error (electrostatic_periodic): invalid rho_L for species ", &
                     reference_species, " value = ", rhoL_rm
            error stop "electrostatic_periodic: invalid active-species rho_L(rm)"
        end select

        ! 3. Window geometry and Fourier cutoff from rho_L(rm).
        dx_asis = periodic_dr_asis_scale * rhoL_rm
        dx_tr   = periodic_dr_tr_scale * rhoL_rm
        L       = 2.0_dp * (dx_asis + dx_tr)
        k_max   = periodic_kmax_scale / rhoL_rm
        M       = ceiling(k_max * L / (2.0_dp * pi))
        n_rg    = periodic_n_rg

        print *, "electrostatic_periodic: reference species = ", reference_species, &
                 " Z = ", plasma%spec(reference_species)%Zspec, &
                 " mass [g] = ", plasma%spec(reference_species)%mass
        print *, "electrostatic_periodic: rm [cm] = ", rm, &
                 " rho_L(rm) [cm] = ", rhoL_rm
        print *, "electrostatic_periodic: dx_asis [cm] = ", dx_asis, &
                 " dx_tr [cm] = ", dx_tr, " L [cm] = ", L
        print *, "electrostatic_periodic: k_max [1/cm] = ", k_max, &
                 " M = ", M, " n_rg = ", n_rg
        call write_periodic_scale_metadata(reference_species, &
            plasma%spec(reference_species)%Zspec, &
            plasma%spec(reference_species)%mass, rhoL_rm, dx_asis, dx_tr, &
            k_max, M, n_rg)

        ! 4. Window output grid: n_rg equidistant points on [rm - L/2, rm + L/2],
        ! bit-identical to the grid build_periodic_plasma installs as rg_grid%xb
        ! via grid_generate_equidistant: spacing h = L/n_rg (NOT L/(n_rg-1)),
        ! endpoint-exclusive at the top (last node = rm + L/2 - L/n_rg).
        allocate(r_win(n_rg))
        do i = 1, n_rg
            r_win(i) = (rm - 0.5_dp * L) + real(i - 1, dp) * L / real(n_rg, dp)
        end do

        ! 5. Build -> assemble -> solve -> reconstruct on r_win via the reusable
        ! periodic core. It does NOT error stop on a singular solve; do it here.
        Br_const = cmplx(Br_boundary_re, Br_boundary_im, dp)
        call compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, Br_const, &
                                        r_win, dPhi, info, lift_projection_residual)
        if (info /= 0) then
            print *, "Error (electrostatic_periodic): solve_periodic failed, info = ", info
            error stop "electrostatic_periodic: periodic solve failed"
        end if
        print *, "electrostatic_periodic: aligned-lift projection residual = ", &
                 lift_projection_residual

        ! Lock: the core has installed rg_grid%xb via build_periodic_plasma; the
        ! run-type's output grid r_win MUST equal it so EBdat is on the window.
        if (maxval(abs(r_win - rg_grid%xb)) > 1.0e-9_dp) then
            error stop "run_electrostatic_periodic: r_win grid does not match rg_grid%xb"
        end if

        ! 6. Pack into EBdat (window grid + reconstructed potential).
        if (allocated(EBdat%r_grid)) deallocate(EBdat%r_grid)
        if (allocated(EBdat%Phi))    deallocate(EBdat%Phi)
        EBdat%r_grid = r_win
        EBdat%Phi    = dPhi

        if (hdf5_output) then
            call write_complex_profile_abs(EBdat%r_grid, EBdat%Phi, rg_grid%npts_b, &
                "/fields/Phi", &
                'Electrostatic potential perturbation Phi, forced-periodicity solution', &
                'statV')
        end if

    end subroutine run_electrostatic_periodic

end module rt_electrostatic_periodic_m
