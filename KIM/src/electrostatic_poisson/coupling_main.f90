program kim_coupling
    ! Thin exchange driver. Physics stays in the existing KIM assembly/solver.
    use KIM_kinds_m, only: dp
    use config_m, only: nml_config_path, periodic_electron_flr, collision_model, &
        resolved_electron_ifunc_conservation_model, artificial_debye_case, &
        periodic_match_global_kernel_approximations, turn_off_electrons, turn_off_ions
    use setup_m, only: m_mode, n_mode, omega, collisions_off, mphi_max
    use constants_m, only: sol
    use species_m, only: plasma
    use grid_m, only: rg_grid
    use kim_base_m, only: kim_t
    use kim_mod_m, only: from_kim_factory_get_kim
    use periodic_background_m, only: build_periodic_plasma
    use periodic_assembly_m, only: assemble_periodic_matrices
    use periodic_solve_m, only: solve_periodic, reconstruct_delta_phi, reconstruct_jpar
    use periodic_electron_response_m, only: electron_response_t
    use rt_electrostatic_periodic_m, only: compute_periodic_delta_phi
    use flr2_fourier_kernel_m, only: set_global_kernel_approximations, kern_zero_flr_electrons
    implicit none
    class(kim_t), allocatable :: instance
    type(electron_response_t) :: electrons
    character(1024) :: config_path, coupling_path = 'coupling.nml', action = 'export', header
    character(1024) :: prefix = 'exchange/', response_file = 'exchange/electron_response.dat'
    integer :: M = 3, n_rg = 96, unit, i, sp, info
    integer :: file_m, file_n, file_model, file_mmode, file_nmode
    real(dp) :: rm = -1.0_dp, dx_asis = 1.0_dp, dx_tr = 2.0_dp, L, re, im, file_L, file_rm, file_c
    complex(dp) :: Br = (1.0_dp, 0.0_dp)
    real(dp), allocatable :: background(:, :), file_background(:, :), r(:)
    complex(dp), allocatable :: kp(:, :), kb(:, :), jp(:, :), jb(:, :)
    complex(dp), allocatable :: kp_sp(:, :, :), kb_sp(:, :, :), jp_sp(:, :, :), jb_sp(:, :, :)
    complex(dp), allocatable :: phi_modes(:), phi(:), current(:), current_sp(:, :)
    namelist /COUPLING/ rm, dx_asis, dx_tr, M, n_rg, Br, prefix, response_file, action
    call get_command_argument(1, config_path)
    if (len_trim(config_path) == 0) &
        error stop 'usage: kim_coupling.x KIM_config.nml (with coupling.nml)'
    open (newunit=unit, file=trim(coupling_path), status='old', action='read')
    read (unit, nml=COUPLING)
    close (unit)
    if (trim(action) /= 'export' .and. trim(action) /= 'solve') &
        error stop 'expected export or solve'
    if (M < 0 .or. n_rg < 2*M + 1 .or. rm <= 0 .or. dx_asis <= 0 .or. dx_tr <= 0) &
        error stop 'invalid coupling radial window or resolution'
    nml_config_path = trim(config_path)
    call kim_init()
    if (periodic_match_global_kernel_approximations .or. turn_off_electrons .or. &
        turn_off_ions) error stop 'coupling requires full ion kernel and both species enabled'
    if (periodic_electron_flr) &
        error stop 'set periodic_electron_flr=.false. for drift-kinetic electrons'
    if (resolved_electron_ifunc_conservation_model /= 0) &
        error stop 'initial cylindrical provider requires number-conserving OU electrons (model 0)'
    if (omega /= 0.0_dp .or. collisions_off .or. artificial_debye_case /= 0) &
        error stop 'initial coupling requires static forcing, finite collisions and full charge'
    if (mphi_max /= 0) error stop 'initial coupling requires the zero cyclotron harmonic'
    if (trim(collision_model) /= 'FokkerPlanck') &
        error stop 'OU electron provider requires FokkerPlanck electron susceptibilities'
    call from_kim_factory_get_kim('electrostatic', instance)
    call instance%init()
    L = 2.0_dp*(dx_asis + dx_tr)
    call set_global_kernel_approximations(.false.)
    kern_zero_flr_electrons = .true.
    call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)
    allocate (background(n_rg, 13))
    background(:, 1) = rg_grid%xb
    background(:, 2) = plasma%ks
    background(:, 3) = plasma%kp
    background(:, 4) = plasma%om_E
    background(:, 5) = plasma%spec(0)%lambda_D
    background(:, 6) = plasma%spec(0)%nu
    background(:, 7) = plasma%spec(0)%vT
    background(:, 8) = plasma%spec(0)%omega_c
    background(:, 9) = plasma%spec(0)%A1
    background(:, 10) = plasma%spec(0)%A2
    background(:, 11) = plasma%spec(0)%n
    background(:, 12) = plasma%spec(0)%T
    background(:, 13) = plasma%B0
    r = rg_grid%xb
    if (trim(action) == 'export') then
        call assemble_periodic_matrices(plasma, L, M, kp, kb, jp, jb, &
                       Kphi_species=kp_sp, KB_species=kb_sp, Kjphi_species=jp_sp, KjB_species=jb_sp)
        open (newunit=unit, file=trim(prefix)//'background.dat', status='new', action='write')
        write (unit, '(a)') 'GK_BACKGROUND_V1'
        write (unit, *) M, n_rg, 0, m_mode, n_mode
        write (unit, '(3es26.17e3)') L, rm, sol
        do i = 1, n_rg
            write (unit, '(13es26.17e3)') background(i, :)
        end do
        close (unit)
        open (newunit=unit, file=trim(prefix)//'kim_blocks.dat', status='new', action='write')
        write (unit, '(a)') 'KIM_BLOCKS_V1'
        write (unit, *) 2*M + 1, plasma%n_species
        do sp = 0, plasma%n_species - 1
            call write_block(kp_sp(:, :, sp))
            call write_block(kb_sp(:, :, sp))
            call write_block(jp_sp(:, :, sp))
            call write_block(jb_sp(:, :, sp))
        end do
        close (unit)
        call solve_periodic(kp, kb, L, M, Br, phi_modes, info)
        if (info /= 0) error stop 'laminar reference solve failed'
        phi = reconstruct_delta_phi(phi_modes, L, M, r)
        current = reconstruct_jpar(jp, jb, phi_modes, Br, L, M, r)
        call write_fields('reference_fields.dat')
    else
        open (newunit=unit, file=trim(response_file), status='old', action='read')
        read (unit, '(a)') header
        if (trim(header) /= 'GK_RESPONSE_V1') error stop 'wrong electron response format'
        read (unit, *) file_m, file_n, file_model, file_mmode, file_nmode
        read (unit, *) file_L, file_rm, file_c
        if (file_m /= M .or. file_n /= n_rg .or. file_model /= 0 .or. &
            file_mmode /= m_mode .or. file_nmode /= n_mode) error stop 'response request mismatch'
        if (file_L /= L .or. file_rm /= rm .or. file_c /= sol) &
            error stop 'response geometry mismatch'
        allocate (file_background(n_rg, 13))
        do i = 1, n_rg
            read (unit, *) file_background(i, :)
        end do
        ! Decimal round trips preserve every background value. Reject stale blocks.
        if (any(file_background /= background)) error stop 'response background mismatch'
        call read_block(electrons%rho_phi)
        call read_block(electrons%rho_b)
        call read_block(electrons%j_phi)
        call read_block(electrons%j_b)
        close (unit)
        call compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, Br, r, phi, info, &
                                jpar=current, jpar_species=current_sp, external_electrons=electrons)
        if (info /= 0) error stop 'hybrid solve failed'
        call write_fields('hybrid_fields.dat')
    end if
contains
    subroutine write_block(block)
        complex(dp), intent(in) :: block(:, :)
        integer :: row, col
        do col = 1, size(block, 2)
            do row = 1, size(block, 1)
                write (unit, '(2es26.17e3)') real(block(row, col)), aimag(block(row, col))
            end do
        end do
    end subroutine
    subroutine read_block(block)
        complex(dp), allocatable, intent(out) :: block(:, :)
        integer :: row, col
        allocate (block(2*M + 1, 2*M + 1))
        do col = 1, 2*M + 1
            do row = 1, 2*M + 1
                read (unit, *) re, im
                block(row, col) = cmplx(re, im, dp)
            end do
        end do
    end subroutine
    subroutine write_fields(name)
        character(*), intent(in) :: name
        integer :: row
        open (newunit=unit, file=trim(prefix)//name, status='new', action='write')
        write (unit, '(a)') '# r[cm] RePhi[statV] ImPhi Rej[statA/cm2] Imj'
        do row = 1, size(r)
            write (unit, '(5es26.17e3)') r(row), real(phi(row)), aimag(phi(row)), &
                real(current(row)), aimag(current(row))
        end do
        close (unit)
    end subroutine
end program
