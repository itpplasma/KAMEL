module rt_flr2_m
    use kim_base_m, only: kim_t
    implicit none
    private

    type, extends(kim_t), public :: flr2_t
    contains
        procedure :: init => init_flr2
        procedure :: run => run_flr2
    end type flr2_t

contains

    subroutine init_flr2(this)
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use grid_m, only: rg_grid
        use IO_collection_m, only: create_output_directories
        use logger_m, only: log_info
        use species_m, only: plasma, set_plasma_quantities

        class(flr2_t), intent(inout) :: this

        this%run_type = 'flr2'
        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)
        call log_info('...flr2 model initialized with KIM background.')
    end subroutine init_flr2

    subroutine run_flr2(this)
        use config_m, only: collision_model, flr2_electron_current, &
                            flr2_electron_flr, flr2_electron_potential, &
                            flr2_include_potential_in_current, flr2_ion_current, &
                            flr2_ion_flr, flr2_ion_potential, turn_off_electrons, &
                            turn_off_ions, periodic_dr_asis_scale, output_path
        use equilibrium_m, only: B0z
        use fields_m, only: EBdat, get_Br_from_txt, postprocess_electric_field_no_write
        use flr2_response_m, only: flr2_options_t, solve_flr2_response, &
                                   FLR2_RESPONSE_OK
        use IO_collection_m, only: write_complex_profile_abs, &
                                   invalidate_kim_diagnostics, &
                                   KIM_DIAGNOSTICS_INVALIDATE_OK
        use kim_diagnostics_m, only: compute_and_write_diagnostics
        use rt_electrostatic_periodic_m, only: select_periodic_reference_scale, &
                                               PERIODIC_SCALE_OK
        use kim_resonances_m, only: r_res
        use KIM_kinds_m, only: dp
        use logger_m, only: log_error, log_info
        use setup_m, only: Br_boundary_im, Br_boundary_re, m_mode, n_mode, &
                           omega, R0, type_br_field
        use species_m, only: plasma

        class(flr2_t), intent(inout) :: this

        type(flr2_options_t) :: options
        complex(dp), allocatable :: bpsi_over_bphi(:), parcur_over_b0(:)
        complex(dp), allocatable :: parcur_e_over_b0(:), parcur_i_over_b0(:)
        complex(dp), allocatable :: principal_coefficient(:)
        complex(dp) :: br_constant
        character(len=160) :: error_message
        integer :: n, stat, reference_species, scale_info, i, iunit
        integer :: diagnostics_status, principal_status
        real(dp) :: rho_ref

        ! Invalidate before any validation or allocation.  An early FLR2
        ! failure must never leave a previous successful scan result behind.
        call invalidate_kim_diagnostics(trim(output_path)//'kim_diagnostics.dat', &
                                        diagnostics_status)
        call invalidate_kim_diagnostics(trim(output_path)//'principal_coefficient.dat', &
                                        principal_status)
        if (diagnostics_status /= KIM_DIAGNOSTICS_INVALIDATE_OK) then
            call log_error('FLR2 could not invalidate kim_diagnostics.dat before run.')
        end if
        if (principal_status /= KIM_DIAGNOSTICS_INVALIDATE_OK) then
            call log_error('FLR2 could not invalidate principal_coefficient.dat before run.')
        end if

        if (trim(collision_model) /= 'FokkerPlanck') then
            call log_error('FLR2 run type requires collision_model = FokkerPlanck.')
        end if
        if (abs(omega) > tiny(1.0_dp)) then
            call log_error('FLR2 run type currently requires omega = 0.')
        end if
        if (plasma%n_species /= 2) then
            call log_error('FLR2 run type currently requires exactly one ion species.')
        end if
        if (plasma%spec(1)%Zspec <= 0) then
            call log_error('FLR2 run type requires a positively charged ion species.')
        end if
        if (turn_off_electrons .and. turn_off_ions) then
            call log_error('FLR2 cannot disable both electrons and ions.')
        end if

        n = plasma%grid_size
        if (allocated(EBdat%r_grid)) deallocate(EBdat%r_grid)
        if (allocated(EBdat%Br)) deallocate(EBdat%Br)
        if (allocated(EBdat%Phi)) deallocate(EBdat%Phi)
        if (allocated(EBdat%jpar)) deallocate(EBdat%jpar)
        if (allocated(EBdat%jpar_e)) deallocate(EBdat%jpar_e)
        if (allocated(EBdat%jpar_i)) deallocate(EBdat%jpar_i)
        allocate(EBdat%r_grid(n), EBdat%Br(n), EBdat%Phi(n), EBdat%jpar(n), &
                 EBdat%jpar_e(n), EBdat%jpar_i(n))
        allocate(bpsi_over_bphi(n), parcur_over_b0(n), parcur_e_over_b0(n), &
                 parcur_i_over_b0(n), principal_coefficient(n))
        EBdat%r_grid = plasma%r_grid

        select case (type_br_field)
            case (11)
                call get_Br_from_txt(EBdat, './inp/Br_in.dat')
            case (12)
                br_constant = cmplx(Br_boundary_re, Br_boundary_im, dp)
                EBdat%Br = br_constant
            case default
                write(error_message, '(A,I0,A)') &
                    'FLR2 supports type_br_field = 11 (file) or 12 (constant), got ', &
                    type_br_field, '.'
                call log_error(trim(error_message))
        end select

        options%electron_flr = flr2_electron_flr
        options%ion_flr = flr2_ion_flr
        options%electron_potential = flr2_electron_potential .and. .not. turn_off_electrons
        options%ion_potential = flr2_ion_potential .and. .not. turn_off_ions
        options%electron_current = flr2_electron_current .and. .not. turn_off_electrons
        options%ion_current = flr2_ion_current .and. .not. turn_off_ions
        options%include_potential_in_current = flr2_include_potential_in_current
        if (.not. options%electron_potential .and. .not. options%ion_potential) then
            call log_error('FLR2 requires at least one species in the potential equation.')
        end if

        call kim_prepare_resonances
        EBdat%r_resonance = r_res
        if (r_res > 0.0_dp) then
            ! FLR2's shared window is defined by the ion Larmor radius,
            ! independently of which species are enabled in attribution.
            call select_periodic_reference_scale(plasma, r_res, .false., .true., &
                rho_ref, reference_species, scale_info)
            if (scale_info /= PERIODIC_SCALE_OK) then
                call log_error('FLR2 could not determine rho_L at its resonant surface.')
            end if
            EBdat%dx_asis = periodic_dr_asis_scale*rho_ref
        end if

        ! KiLCA-FLR2's radial-variable formulation uses B^r/B^phi, with
        ! B_phi = B0z*R0 in KIM's cylindrical equilibrium. Keep B0z signed.
        bpsi_over_bphi = EBdat%Br*R0/B0z
        call solve_flr2_response(plasma, m_mode, n_mode, R0, B0z, &
                                 bpsi_over_bphi, options, EBdat%Phi, &
                                 parcur_over_b0, stat, &
                                 jpar_e=parcur_e_over_b0, jpar_i=parcur_i_over_b0, &
                                 principal_coefficient=principal_coefficient)
        if (stat /= FLR2_RESPONSE_OK) then
            write(error_message, '(A,I0)') 'FLR2 response solve failed with status ', stat
            call log_error(trim(error_message))
        end if
        EBdat%jpar = parcur_over_b0*plasma%B0
        EBdat%jpar_e = parcur_e_over_b0*plasma%B0
        EBdat%jpar_i = parcur_i_over_b0*plasma%B0

        call postprocess_electric_field_no_write(EBdat)

        open(newunit=iunit, file=trim(output_path)//'principal_coefficient.dat', &
             status='replace', action='write')
        write(iunit, '(A)') '# r [cm]  Re(a2_out+a2_in)  Im(a2_out+a2_in)'
        do i = 1, n
            write(iunit, '(3(es23.15,1x))') EBdat%r_grid(i), &
                real(principal_coefficient(i), dp), aimag(principal_coefficient(i))
        end do
        close(iunit)

        call write_complex_profile_abs(EBdat%r_grid, EBdat%Br, n, &
            '/fields/Br', 'Radial magnetic-field perturbation', 'G')
        call write_complex_profile_abs(EBdat%r_grid, EBdat%Phi, n, &
            '/fields/Phi', 'FLR2 electrostatic potential perturbation', 'statV')
        call write_complex_profile_abs(EBdat%r_grid, EBdat%jpar, n, &
            '/fields/jpar', 'FLR2 parallel current density perturbation', 'statA/cm^2')
        call compute_and_write_diagnostics(require_window=.true.)
        call log_info('...'//trim(this%run_type)//' response solved.')
    end subroutine run_flr2

end module rt_flr2_m
