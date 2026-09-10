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
    integer, parameter, public :: PERIODIC_BPAR_OK = 0
    integer, parameter, public :: PERIODIC_BPAR_UNSUPPORTED_COLLISION = 1
    integer, parameter, public :: PERIODIC_BPAR_UNSUPPORTED_HARMONIC = 2
    integer, parameter, public :: PERIODIC_BPAR_UNSUPPORTED_APPROXIMATION = 3
    integer, parameter, public :: PERIODIC_BPAR_UNSUPPORTED_DEBYE_MODEL = 4
    integer, parameter, public :: PERIODIC_BPAR_INVALID_DRIVE = 5

    public :: compute_periodic_delta_phi, select_periodic_reference_scale
    public :: compute_periodic_ion_tensor, compute_periodic_ion_tensor_spectrum
    public :: periodic_bparallel_support_status

    contains

    pure integer function periodic_bparallel_support_status(Bparallel_drive, &
            electron_collision_model, ion_model, collisions_disabled, debye_case, &
            harmonic_max, match_global_approximations) result(status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use KIM_kinds_m, only: dp
        complex(dp), intent(in) :: Bparallel_drive
        character(*), intent(in) :: electron_collision_model, ion_model
        integer, intent(in) :: debye_case, harmonic_max
        logical, intent(in) :: collisions_disabled, match_global_approximations

        status = PERIODIC_BPAR_OK
        if (Bparallel_drive == (0.0_dp, 0.0_dp)) return
        if (.not. ieee_is_finite(real(Bparallel_drive, dp)) .or. &
                .not. ieee_is_finite(aimag(Bparallel_drive))) then
            status = PERIODIC_BPAR_INVALID_DRIVE
        else if (trim(electron_collision_model) /= 'FokkerPlanck' .or. &
                trim(ion_model) /= 'FokkerPlanck' .or. collisions_disabled) then
            status = PERIODIC_BPAR_UNSUPPORTED_COLLISION
        else if (debye_case /= 0 .and. debye_case /= 2) then
            status = PERIODIC_BPAR_UNSUPPORTED_DEBYE_MODEL
        else if (harmonic_max /= 0) then
            status = PERIODIC_BPAR_UNSUPPORTED_HARMONIC
        else if (match_global_approximations) then
            status = PERIODIC_BPAR_UNSUPPORTED_APPROXIMATION
        end if
    end function periodic_bparallel_support_status

    subroutine validate_periodic_bparallel(Bparallel_drive)
        use KIM_kinds_m, only: dp
        use config_m, only: artificial_debye_case, collision_model, &
            ion_collision_model, periodic_match_global_kernel_approximations
        use setup_m, only: collisions_off, mphi_max
        complex(dp), intent(in) :: Bparallel_drive
        integer :: status

        status = periodic_bparallel_support_status(Bparallel_drive, &
            collision_model, ion_collision_model, collisions_off, &
            artificial_debye_case, mphi_max, &
            periodic_match_global_kernel_approximations)
        select case (status)
        case (PERIODIC_BPAR_OK)
            return
        case (PERIODIC_BPAR_UNSUPPORTED_COLLISION)
            error stop 'Nonzero periodic Bparallel requires enabled FokkerPlanck collisions'
        case (PERIODIC_BPAR_UNSUPPORTED_DEBYE_MODEL)
            error stop 'Nonzero periodic Bparallel requires full or no-Debye FP response'
        case (PERIODIC_BPAR_UNSUPPORTED_HARMONIC)
            error stop 'Nonzero periodic Bparallel currently requires mphi_max=0'
        case (PERIODIC_BPAR_UNSUPPORTED_APPROXIMATION)
            error stop 'Nonzero periodic Bparallel cannot drop ks from gyrogeometry'
        case (PERIODIC_BPAR_INVALID_DRIVE)
            error stop 'Periodic Bparallel drive must be finite'
        case default
            error stop 'Unknown periodic Bparallel validation status'
        end select
    end subroutine validate_periodic_bparallel

    subroutine compute_periodic_ion_tensor(fields_s, ks_s, kr_s, kpar, vTi, nui, omega_ci, &
                                            omega_mode, om_E, B0, tensor)
        use KIM_kinds_m, only: dp
        use config_m, only: resolved_ion_ifunc_conservation_model
        use setup_m, only: mphi_max
        use species_m, only: evaluate_susceptibility
        use quasilinear_flr_m, only: calc_ion_flr_harmonic
        use constants_m, only: sol

        complex(dp), intent(in) :: fields_s(3)
        real(dp), intent(in) :: ks_s, kr_s, kpar, vTi, nui, omega_ci
        real(dp), intent(in) :: omega_mode, om_E, B0
        real(dp), intent(out) :: tensor(2,2)
        complex(dp) :: symbI(0:3,0:3), fields_o(3)
        real(dp) :: harmonic_tensor(2,2), x1, x2
        integer :: ell

        tensor = 0.0_dp
        if (nui <= tiny(1.0_dp) .or. abs(omega_ci) <= tiny(1.0_dp)) return
        fields_o = fields_s
        x1 = kpar*vTi/nui
        do ell = -mphi_max, mphi_max
            x2 = -(om_E + real(ell,dp)*omega_ci - omega_mode)/nui
            call evaluate_susceptibility(x1, x2, resolved_ion_ifunc_conservation_model, symbI)
            call calc_ion_flr_harmonic(ell, ks_s, kr_s, ks_s, kr_s, vTi, abs(omega_ci), &
                omega_ci, sol, B0, nui, fields_s, fields_o, symbI, harmonic_tensor)
            tensor = tensor + harmonic_tensor
        end do
    end subroutine compute_periodic_ion_tensor


    subroutine compute_periodic_ion_tensor_spectrum(phi_m, br, period, radius, ks, kpar, &
            vti, nui, omega_ci, omega_mode, om_e, b0, tensor, bparallel)
        use KIM_kinds_m, only: dp
        use config_m, only: resolved_ion_ifunc_conservation_model
        use setup_m, only: mphi_max
        use species_m, only: evaluate_susceptibility
        use quasilinear_flr_m, only: calc_ion_flr_harmonic
        use constants_m, only: pi, sol, com_unit

        complex(dp), intent(in) :: phi_m(:), br
        complex(dp), intent(in), optional :: bparallel
        real(dp), intent(in) :: period, radius, ks, kpar, vti, nui, omega_ci
        real(dp), intent(in) :: omega_mode, om_e, b0
        real(dp), intent(out) :: tensor(2,2)
        complex(dp) :: fields(3,size(phi_m)), moments(0:3,0:3)
        real(dp) :: kr(size(phi_m)), pair_tensor(2,2), x1, x2
        integer :: cutoff, i, source, observation, ell

        if (period <= 0.0_dp .or. size(phi_m) < 1 .or. mod(size(phi_m),2) /= 1) &
            error stop 'periodic ion tensor requires a positive period and odd Fourier spectrum'
        tensor = 0.0_dp
        if (nui <= tiny(1.0_dp) .or. abs(omega_ci) <= tiny(1.0_dp)) return
        cutoff = (size(phi_m)-1)/2
        fields = (0.0_dp, 0.0_dp)
        do i = 1, size(phi_m)
            kr(i) = 2.0_dp*pi*real(i-cutoff-1,dp)/period
            fields(1,i) = phi_m(i)*exp(com_unit*kr(i)*radius)
        end do
        fields(2,cutoff+1) = br
        if (present(bparallel)) fields(3,cutoff+1) = bparallel

        ! Phi_m are Fourier-series amplitudes: Phi(r)=sum Phi_m exp(i*k_m*r).
        ! Constant Br occupies only k_r=0. Put the spatial phase in each field
        ! BEFORE the kernel takes its real part, retaining complex interference.
        ! Each ordered pair is included once; no additional Fourier weight or
        ! factor of two belongs on the kernel's Hermitianized quadratic form.
        x1 = kpar*vti/nui
        do ell = -mphi_max, mphi_max
            x2 = -(om_e+real(ell,dp)*omega_ci-omega_mode)/nui
            call evaluate_susceptibility(x1, x2, resolved_ion_ifunc_conservation_model, moments)
            do observation = 1, size(phi_m)
                if (all(fields(:,observation) == (0.0_dp,0.0_dp))) cycle
                do source = 1, size(phi_m)
                    if (all(fields(:,source) == (0.0_dp,0.0_dp))) cycle
                    call calc_ion_flr_harmonic(ell, ks, kr(source), ks, kr(observation), &
                        vti, abs(omega_ci), omega_ci, sol, b0, nui, fields(:,source), &
                        fields(:,observation), moments, pair_tensor)
                    tensor = tensor + pair_tensor
                end do
            end do
        end do
    end subroutine compute_periodic_ion_tensor_spectrum

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
    !> on the same r_out. With a prescribed compression drive,
    !> j_par = K^{jPhi} Phi + K^{jB} Br + K^{jBparallel} Bparallel.
    !> The optional jpar_species(:,sp) returns the contribution from each species,
    !> with electron index sp=0. Both are left unallocated when the solve fails.
    subroutine compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, Br_const, &
            r_out, dPhi, info, jpar, jpar_species, rho_B, rho_B_species, &
            dPhi_dr, jrad, phi_spectrum, Bparallel_const, rho_Bparallel, &
            rho_Bparallel_species)
        use KIM_kinds_m, only: dp
        use species_m, only: plasma
        use periodic_background_m, only: build_periodic_plasma
        use periodic_assembly_m, only: assemble_periodic_matrices, &
            assemble_periodic_bparallel_matrices
        use periodic_solve_m, only: solve_periodic, reconstruct_delta_phi, &
            reconstruct_delta_phi_derivative, reconstruct_jpar, reconstruct_jrad
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
        complex(dp), allocatable, intent(out), optional :: dPhi_dr(:), phi_spectrum(:)
        complex(dp), allocatable, intent(out), optional :: jrad(:)
        complex(dp), intent(in), optional :: Bparallel_const
        complex(dp), allocatable, intent(out), optional :: rho_Bparallel(:)
        complex(dp), allocatable, intent(out), optional :: rho_Bparallel_species(:,:)

        complex(dp), allocatable :: Kphi(:,:), KB(:,:), Kjphi(:,:), KjB(:,:), Phi_m(:)
        complex(dp), allocatable :: Kjrphi(:,:), KjrB(:,:)
        complex(dp), allocatable :: Kjphi_species(:,:,:), KjB_species(:,:,:)
        complex(dp), allocatable :: Kphi_species(:,:,:), KB_species(:,:,:)
        complex(dp), allocatable :: KBparallel(:,:), KjBparallel(:,:), KjrBparallel(:,:)
        complex(dp), allocatable :: KBparallel_species(:,:,:)
        complex(dp), allocatable :: KjBparallel_species(:,:,:)
        complex(dp) :: Bparallel_drive
        real(dp) :: L
        integer :: sp
        logical :: bparallel_active, want_bparallel_species

        L = 2.0_dp * (dx_asis + dx_tr)
        Bparallel_drive = (0.0_dp, 0.0_dp)
        if (present(Bparallel_const)) Bparallel_drive = Bparallel_const
        call validate_periodic_bparallel(Bparallel_drive)
        bparallel_active = Bparallel_drive /= (0.0_dp, 0.0_dp)

        call set_global_kernel_approximations(periodic_match_global_kernel_approximations)
        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)
        if (present(jpar_species) .and. present(rho_B_species)) then
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
                Kjrphi, KjrB, Kjphi_species, KjB_species, Kphi_species, KB_species)
        else if (present(jpar_species)) then
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
                Kjrphi, KjrB, Kjphi_species, KjB_species)
        else if (present(rho_B_species)) then
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
                Kjrphi, KjrB, &
                Kphi_species=Kphi_species, KB_species=KB_species)
        else
            call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
                Kjrphi, KjrB)
        end if
        if (bparallel_active) then
            want_bparallel_species = present(jpar_species) .or. &
                present(rho_Bparallel_species)
            if (want_bparallel_species) then
                call assemble_periodic_bparallel_matrices(plasma, L, M, &
                    KBparallel, KjBparallel, KBparallel_species, &
                    KjBparallel_species, KjrBparallel)
            else
                call assemble_periodic_bparallel_matrices(plasma, L, M, &
                    KBparallel, KjBparallel, KjrBparallel=KjrBparallel)
            end if
            call solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info, &
                KBparallel=KBparallel, Bparallel_const=Bparallel_drive)
        else
            call solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info)
        end if
        if (info /= 0) return

        if (present(phi_spectrum)) phi_spectrum = Phi_m
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

        if (present(rho_Bparallel)) then
            allocate(rho_Bparallel(size(r_out)))
            rho_Bparallel = (0.0_dp, 0.0_dp)
            if (bparallel_active) then
                rho_Bparallel = reconstruct_delta_phi(&
                    Bparallel_drive * KBparallel(:, M + 1), L, M, r_out)
            end if
        end if

        if (present(rho_Bparallel_species)) then
            allocate(rho_Bparallel_species(size(r_out), 0:plasma%n_species - 1))
            rho_Bparallel_species = (0.0_dp, 0.0_dp)
            if (bparallel_active) then
                do sp = 0, plasma%n_species - 1
                    rho_Bparallel_species(:, sp) = reconstruct_delta_phi(&
                        Bparallel_drive * KBparallel_species(:, M + 1, sp), &
                        L, M, r_out)
                end do
            end if
        end if

        if (present(jpar_species)) then
            allocate(jpar_species(size(r_out), 0:plasma%n_species - 1))
            do sp = 0, plasma%n_species - 1
                if (bparallel_active) then
                    jpar_species(:, sp) = reconstruct_jpar(&
                        Kjphi_species(:, :, sp), KjB_species(:, :, sp), &
                        Phi_m, Br_const, L, M, r_out, &
                        KjBparallel=KjBparallel_species(:, :, sp), &
                        Bparallel_const=Bparallel_drive)
                else
                    jpar_species(:, sp) = reconstruct_jpar(&
                        Kjphi_species(:, :, sp), KjB_species(:, :, sp), &
                        Phi_m, Br_const, L, M, r_out)
                end if
            end do
        end if

        if (present(jpar)) then
            if (bparallel_active) then
                jpar = reconstruct_jpar(Kjphi, KjB, Phi_m, Br_const, L, M, &
                    r_out, KjBparallel=KjBparallel, &
                    Bparallel_const=Bparallel_drive)
            else
                jpar = reconstruct_jpar(Kjphi, KjB, Phi_m, Br_const, L, M, r_out)
            end if
        end if
        if (present(jrad)) then
            if (bparallel_active) then
                jrad = reconstruct_jrad(Kjrphi, KjrB, Phi_m, Br_const, L, M, r_out, &
                    KjrBparallel=KjrBparallel, Bparallel_const=Bparallel_drive)
            else
                jrad = reconstruct_jrad(Kjrphi, KjrB, Phi_m, Br_const, L, M, r_out)
            end if
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
                            periodic_Bparallel_ratio, turn_off_electrons, turn_off_ions
        use setup_m, only: Br_boundary_re, Br_boundary_im, m_mode, n_mode, R0, omega
        use species_m, only: plasma
        use grid_m, only: rg_grid
        use kim_resonances_m, only: r_res
        use fields_m, only: EBdat, EBdat_t
        use fields_m, only: calculate_MA_field, calculate_E_in_rsp_from_cyl
        use IO_collection_m, only: itoa, write_complex_profile_abs, &
            write_periodic_scale_metadata

        implicit none

        class(electrostatic_periodic_t), intent(inout) :: this

        complex(dp), allocatable :: dPhi(:), dPhi_dr(:), jpar(:), jpar_species(:,:)
        complex(dp), allocatable :: jrad(:)
        complex(dp), allocatable :: rho_B(:), rho_B_species(:,:), rho_B_i(:)
        complex(dp), allocatable :: rho_Bparallel(:), rho_Bparallel_species(:,:)
        complex(dp), allocatable :: rho_Bparallel_i(:)
        complex(dp), allocatable :: phi_spectrum(:)
        real(dp) :: tensor_local(2,2)
        complex(dp) :: Br_const, Bparallel_const
        real(dp), allocatable :: r_win(:)
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr, L, k_max
        integer :: M, n_rg, info, i, reference_species, sp

        ! 1. Locate the resonant surface rm = r_res (q = -m/n) on the global plasma.
        call kim_prepare_resonances
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
        Bparallel_const = periodic_Bparallel_ratio * Br_const
        call compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, Br_const, &
            r_win, dPhi, info, jpar, jpar_species, rho_B, rho_B_species, &
            dPhi_dr=dPhi_dr, jrad=jrad, phi_spectrum=phi_spectrum, &
            Bparallel_const=Bparallel_const, rho_Bparallel=rho_Bparallel, &
            rho_Bparallel_species=rho_Bparallel_species)
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
        ! A periodic solve may be repeated with a different window. Reset the
        ! complete allocatable record so newly added response fields cannot
        ! retain stale allocations or shapes between direct runs.
        EBdat = EBdat_t()
        EBdat%r_grid = r_win
        EBdat%r_resonance = rm
        EBdat%dx_asis = dx_asis
        EBdat%dx_transition = dx_tr
        EBdat%Phi    = dPhi
        allocate(EBdat%Br(size(dPhi)))
        EBdat%Br = Br_const
        allocate(EBdat%Bparallel(size(dPhi)))
        EBdat%Bparallel = Bparallel_const
        allocate(EBdat%Er(size(dPhi)), EBdat%Etheta(size(dPhi)), EBdat%Ez(size(dPhi)))
        EBdat%Er = -dPhi_dr
        do i = 1, size(dPhi)
            EBdat%Etheta(i) = -com_unit * m_mode * dPhi(i) / r_win(i)
            EBdat%Ez(i) = -com_unit * n_mode * dPhi(i) / R0
        end do
        call calculate_MA_field(plasma, EBdat)
        call calculate_E_in_rsp_from_cyl(EBdat)
        allocate(EBdat%D_ion(2, 2, size(dPhi)))
        EBdat%D_ion = 0.0_dp
        if (.not. turn_off_ions) then
            do i = 1, size(dPhi)
                do sp = 1, plasma%n_species - 1
                    call compute_periodic_ion_tensor_spectrum(phi_spectrum, Br_const, L, r_win(i), &
                        plasma%ks(i), plasma%kp(i), plasma%spec(sp)%vT(i), plasma%spec(sp)%nu(i), &
                        plasma%spec(sp)%omega_c(i), omega, plasma%om_E(i), &
                        plasma%B0(i), tensor_local, Bparallel_const)
                    EBdat%D_ion(:, :, i) = EBdat%D_ion(:, :, i) + tensor_local
                end do
            end do
        end if
        EBdat%jpar   = jpar
        EBdat%jpar_e = jpar_species(:, 0)
        allocate(EBdat%jpar_i(size(jpar)))
        EBdat%jpar_i = (0.0_dp, 0.0_dp)
        EBdat%jrad   = jrad
        if (plasma%n_species > 1) then
            EBdat%jpar_i = sum(jpar_species(:, 1:plasma%n_species - 1), dim=2)
        end if
        allocate(rho_B_i(size(rho_B)))
        rho_B_i = (0.0_dp, 0.0_dp)
        if (plasma%n_species > 1) then
            rho_B_i = sum(rho_B_species(:, 1:plasma%n_species - 1), dim=2)
        end if
        allocate(rho_Bparallel_i(size(rho_Bparallel)))
        rho_Bparallel_i = (0.0_dp, 0.0_dp)
        if (plasma%n_species > 1) then
            rho_Bparallel_i = sum(&
                rho_Bparallel_species(:, 1:plasma%n_species - 1), dim=2)
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
            call write_complex_profile_abs(EBdat%r_grid, EBdat%jrad, rg_grid%npts_b, &
                "/fields/jrad", &
                'Radial current density perturbation j_rad, forced-periodicity solution', &
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
            call write_complex_profile_abs(EBdat%r_grid, rho_Bparallel, rg_grid%npts_b, &
                "/fields/rho_Bparallel", &
                'Charge density driven directly by prescribed parallel magnetic field', &
                'statC/cm^3')
            call write_complex_profile_abs(EBdat%r_grid, rho_Bparallel_species(:, 0), &
                rg_grid%npts_b, "/fields/rho_Bparallel_e", &
                'Electron charge density driven by prescribed parallel magnetic field', &
                'statC/cm^3')
            call write_complex_profile_abs(EBdat%r_grid, rho_Bparallel_i, rg_grid%npts_b, &
                "/fields/rho_Bparallel_i", &
                'Summed ion charge density driven by prescribed parallel magnetic field', &
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
                call write_complex_profile_abs(EBdat%r_grid, rho_Bparallel_species(:, sp), &
                    rg_grid%npts_b, "/fields/rho_Bparallel_i"//trim(itoa(sp)), &
                    'Ion-species charge density driven by prescribed parallel magnetic field', &
                    'statC/cm^3')
            end do
        end if

    end subroutine run_electrostatic_periodic

end module rt_electrostatic_periodic_m
