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

    integer, parameter, public :: PERIODIC_SCALE_OK = 0
    integer, parameter, public :: PERIODIC_SCALE_NO_ACTIVE = 1
    integer, parameter, public :: PERIODIC_SCALE_INVALID_RHO = 2

    public :: compute_periodic_delta_phi, select_periodic_reference_scale

    contains

    !> Reusable periodic core: build the window plasma, assemble the Fourier
    !> matrices, solve for the coefficients Phi_m under a constant Br drive, and
    !> reconstruct delta_Phi on the EXPLICIT output grid r_out (NOT on the window
    !> grid rg_grid%xb). All window sizing is passed in EXPLICITLY so both the
    !> run-type and the Phase-3 convergence harness can evaluate delta_Phi on a
    !> fixed diagnostic grid independent of the (n_rg, M, dx) discretization.
    !>
    !> Does NOT error stop on a singular solve: it returns info /= 0 so the
    !> caller can decide (the run-type error stops; the tests inspect info).
    !>
    !> The optional jpar returns the total parallel current density perturbation
    !> on the same r_out, per thesis (11.7): j_par = K^{jPhi} Phi + K^{jB} Br.
    !> The optional jpar_species(:,sp) returns the contribution from each species,
    !> with electron index sp=0. Both are left unallocated when the solve fails.
    subroutine compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, Br_const, &
                                          r_out, dPhi, info, jpar, jpar_species, rho_B, rho_B_species, dPhi_dr)
        use KIM_kinds_m, only: dp
        use species_m, only: plasma
        use periodic_background_m, only: build_periodic_plasma
        use periodic_assembly_m, only: assemble_periodic_matrices
        use periodic_solve_m, only: solve_periodic, reconstruct_delta_phi, reconstruct_delta_phi_derivative, reconstruct_jpar
        use config_m, only: periodic_match_global_kernel_approximations
        use flr2_fourier_kernel_m, only: set_global_kernel_approximations

        real(dp),    intent(in)  :: rm, dx_asis, dx_tr
        integer,     intent(in)  :: M, n_rg
        complex(dp), intent(in)  :: Br_const
        real(dp),    intent(in)  :: r_out(:)
        complex(dp), allocatable, intent(out) :: dPhi(:)
        integer,     intent(out) :: info
        complex(dp), allocatable, intent(out), optional :: jpar(:)
        complex(dp), allocatable, intent(out), optional :: jpar_species(:,:)
        complex(dp), allocatable, intent(out), optional :: rho_B(:)
        complex(dp), allocatable, intent(out), optional :: rho_B_species(:,:)
        complex(dp), allocatable, intent(out), optional :: dPhi_dr(:)

        complex(dp), allocatable :: Kphi(:,:), KB(:,:), Kjphi(:,:), KjB(:,:), Phi_m(:)
        complex(dp), allocatable :: Kjphi_species(:,:,:), KjB_species(:,:,:)
        complex(dp), allocatable :: Kphi_species(:,:,:), KB_species(:,:,:)
        real(dp) :: L
        integer :: sp

        L = 2.0_dp * (dx_asis + dx_tr)

        call set_global_kernel_approximations(periodic_match_global_kernel_approximations)
        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)
        if (present(jpar_species) .and. present(rho_B_species)) then
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
                Kjphi_species, KjB_species, Kphi_species, KB_species)
        else if (present(jpar_species)) then
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
                Kjphi_species, KjB_species)
        else if (present(rho_B_species)) then
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
                Kphi_species=Kphi_species, KB_species=KB_species)
        else
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB)
        end if
        call solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info)
        if (info /= 0) return

        dPhi = reconstruct_delta_phi(Phi_m, L, M, r_out)
        if (present(dPhi_dr)) dPhi_dr = reconstruct_delta_phi_derivative(Phi_m, L, M, r_out)
        if (present(rho_B)) then
            rho_B = reconstruct_delta_phi(Br_const * KB(:, M + 1), L, M, r_out)
        end if

        if (present(rho_B_species)) then
            allocate(rho_B_species(size(r_out), 0:plasma%n_species - 1))
            do sp = 0, plasma%n_species - 1
                rho_B_species(:, sp) = reconstruct_delta_phi(&
                    Br_const * KB_species(:, M + 1, sp), L, M, r_out)
            end do
        end if

        if (present(jpar_species)) then
            allocate(jpar_species(size(r_out), 0:plasma%n_species - 1))
            do sp = 0, plasma%n_species - 1
                jpar_species(:, sp) = reconstruct_jpar(&
                    Kjphi_species(:, :, sp), KjB_species(:, :, sp), &
                    Phi_m, Br_const, L, M, r_out)
            end do
        end if

        if (present(jpar)) then
            jpar = reconstruct_jpar(Kjphi, KjB, Phi_m, Br_const, L, M, r_out)
        end if
    end subroutine compute_periodic_delta_phi

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
        use constants_m, only: pi, com_unit
        use config_m, only: periodic_dr_asis_scale, periodic_dr_tr_scale, &
                            periodic_kmax_scale, periodic_n_rg, hdf5_output, &
                            periodic_match_global_kernel_approximations, &
                            turn_off_electrons, turn_off_ions
        use setup_m, only: Br_boundary_re, Br_boundary_im, m_mode, n_mode, R0
        use species_m, only: plasma
        use grid_m, only: rg_grid
        use kim_resonances_m, only: r_res
        use fields_m, only: EBdat
        use fields_m, only: calculate_MA_field, calculate_E_in_rsp_from_cyl
        use IO_collection_m, only: itoa, write_complex_profile_abs, &
            write_periodic_scale_metadata

        implicit none

        class(electrostatic_periodic_t), intent(inout) :: this

        complex(dp), allocatable :: dPhi(:), dPhi_dr(:), jpar(:), jpar_species(:,:)
        complex(dp), allocatable :: rho_B(:), rho_B_species(:,:), rho_B_i(:)
        complex(dp) :: Br_const
        real(dp), allocatable :: r_win(:)
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr, L, k_max
        integer :: M, n_rg, info, i, reference_species, sp

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
        print *, "electrostatic_periodic: match global kernel approximations = ", &
                 periodic_match_global_kernel_approximations
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
                                        r_win, dPhi, info, jpar, jpar_species, &
                                        rho_B, rho_B_species, dPhi_dr)
        if (info /= 0) then
            print *, "Error (electrostatic_periodic): solve_periodic failed, info = ", info
            error stop "electrostatic_periodic: periodic solve failed"
        end if

        ! Lock: the core has installed rg_grid%xb via build_periodic_plasma; the
        ! run-type's output grid r_win MUST equal it so EBdat is on the window.
        if (maxval(abs(r_win - rg_grid%xb)) > 1.0e-9_dp) then
            error stop "run_electrostatic_periodic: r_win grid does not match rg_grid%xb"
        end if

        ! 6. Pack into EBdat (window grid + reconstructed potential + current).
        if (allocated(EBdat%r_grid)) deallocate(EBdat%r_grid)
        if (allocated(EBdat%Phi))    deallocate(EBdat%Phi)
        if (allocated(EBdat%jpar))   deallocate(EBdat%jpar)
        if (allocated(EBdat%jpar_e)) deallocate(EBdat%jpar_e)
        if (allocated(EBdat%jpar_i)) deallocate(EBdat%jpar_i)
        EBdat%r_grid = r_win
        EBdat%Phi    = dPhi
        allocate(EBdat%Br(size(dPhi)))
        EBdat%Br = Br_const
        allocate(EBdat%Bparallel(size(dPhi)))
        EBdat%Bparallel = (0.0_dp, 0.0_dp)
        allocate(EBdat%Er(size(dPhi)), EBdat%Etheta(size(dPhi)), EBdat%Ez(size(dPhi)))
        EBdat%Er = -dPhi_dr
        do i = 1, size(dPhi)
            EBdat%Etheta(i) = -com_unit * m_mode * dPhi(i) / r_win(i)
            EBdat%Ez(i) = -com_unit * n_mode * dPhi(i) / R0
        end do
        call calculate_MA_field(plasma, EBdat)
        call calculate_E_in_rsp_from_cyl(EBdat)
        EBdat%jpar   = jpar
        EBdat%jpar_e = jpar_species(:, 0)
        allocate(EBdat%jpar_i(size(jpar)))
        EBdat%jpar_i = (0.0_dp, 0.0_dp)
        if (plasma%n_species > 1) then
            EBdat%jpar_i = sum(jpar_species(:, 1:plasma%n_species - 1), dim=2)
        end if
        allocate(rho_B_i(size(rho_B)))
        rho_B_i = (0.0_dp, 0.0_dp)
        if (plasma%n_species > 1) then
            rho_B_i = sum(rho_B_species(:, 1:plasma%n_species - 1), dim=2)
        end if

        if (hdf5_output) then
            call write_complex_profile_abs(EBdat%r_grid, EBdat%Phi, rg_grid%npts_b, &
                "/fields/Phi", &
                'Electrostatic potential perturbation Phi, forced-periodicity solution', &
                'statV')
            call write_complex_profile_abs(EBdat%r_grid, EBdat%jpar, rg_grid%npts_b, &
                "/fields/jpar", &
                'Parallel current density perturbation j_par, forced-periodicity solution', &
                'statA/cm^2')
            call write_complex_profile_abs(EBdat%r_grid, EBdat%jpar_e, rg_grid%npts_b, &
                "/fields/jpar_e", &
                'Electron parallel current density, forced-periodicity solution', &
                'statA/cm^2')
            call write_complex_profile_abs(EBdat%r_grid, EBdat%jpar_i, rg_grid%npts_b, &
                "/fields/jpar_i", &
                'Summed ion parallel current density, forced-periodicity solution', &
                'statA/cm^2')
            call write_complex_profile_abs(EBdat%r_grid, rho_B, rg_grid%npts_b, &
                "/fields/rho_B", &
                'Charge density driven directly by imposed radial magnetic field', &
                'statC/cm^3')
            call write_complex_profile_abs(EBdat%r_grid, rho_B_species(:, 0), rg_grid%npts_b, &
                "/fields/rho_B_e", &
                'Electron charge density driven directly by imposed radial magnetic field', &
                'statC/cm^3')
            call write_complex_profile_abs(EBdat%r_grid, rho_B_i, rg_grid%npts_b, &
                "/fields/rho_B_i", &
                'Summed ion charge density driven directly by imposed radial magnetic field', &
                'statC/cm^3')
            do sp = 1, plasma%n_species - 1
                call write_complex_profile_abs(EBdat%r_grid, jpar_species(:, sp), &
                    rg_grid%npts_b, "/fields/jpar_i"//trim(itoa(sp)), &
                    'Ion-species parallel current density, forced-periodicity solution', &
                    'statA/cm^2')
                call write_complex_profile_abs(EBdat%r_grid, rho_B_species(:, sp), &
                    rg_grid%npts_b, "/fields/rho_B_i"//trim(itoa(sp)), &
                    'Ion-species charge density driven directly by imposed radial magnetic field', &
                    'statC/cm^3')
            end do
        end if

    end subroutine run_electrostatic_periodic

end module rt_electrostatic_periodic_m
