program test_periodic_multimode_feedback
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kim_wave_code_adapter_m, only: kim_initialize, kim_run_for_all_modes, &
        kim_update_profiles, kim_get_wave_fields, kim_get_current_densities, &
        kim_D_ion_modes, kim_transition_weights, kim_embedding_metadata, &
        kim_mode_m, kim_mode_n, kim_mode_status, kim_mode_resonance, &
        kim_periodic_scale_modes, kim_periodic_current_unit, kim_current_records
    use control_mod, only: wave_code, kim_config_path, kim_profiles_from_balance, &
        type_of_run, kim_run_type, kim_current_relaxation
    use wave_code_data, only: dim_mn, m_vals, n_vals, r, n, Te, Ti, q, &
        Vth, Vz, dPhi0, Es, Ep, Er, Et, Ez, Br, Bp, Jpe, Jpi, I_par_toroidal
    use plasma_parameters, only: params_b
    use grid_mod, only: Ercov
    use baseparam_mod, only: ev
    implicit none
    integer, parameter :: npts = 401
    type :: response_t
        complex(8), allocatable :: fields(:,:,:), currents(:,:,:), unit(:), scale(:)
        real(8), allocatable :: tensor(:,:,:,:), weights(:,:), window(:,:), resonance(:)
    end type
    type(response_t) :: reference, reordered, repeated, refreshed, fresh, single
    real(8) :: input_r(npts), saved_q(npts), saved_params(4,npts), saved_er(npts)
    character(len=32) :: argument
    integer :: point, scenario, profile
    external :: allocate_wave_code_data

    dim_mn = 2
    allocate(m_vals(2), n_vals(2))
    m_vals = [-6, -7]
    n_vals = 2
    do point = 1, npts
        input_r(point) = 2.0d0 + 58.0d0*real(point-1,8)/real(npts-1,8)
    end do
    call allocate_wave_code_data(npts, input_r)
    n = 2.0d13*(1.1d0-r/100.0d0)
    Te = 1.0d3*(1.2d0-r/100.0d0)
    Ti = Te
    q = 1.5d0+3.0d0*(r-2.0d0)/58.0d0
    dPhi0 = 0.5d0
    Vth = 0.0d0
    Vz = 0.0d0
    allocate(params_b(4,npts), Ercov(npts))
    params_b(1,:) = n
    params_b(2,:) = 0.0d0
    params_b(3,:) = Te*ev
    params_b(4,:) = Ti*ev
    Ercov = -dPhi0
    saved_q = q
    saved_params = params_b
    saved_er = Ercov
    wave_code = 'KIM'
    kim_run_type = 'electrostatic_periodic'
    kim_profiles_from_balance = .true.
    type_of_run = 'ParameterScan'
    call get_command_argument(1, argument)
    kim_config_path = 'KIM_config_multimode_feedback.nml'
    if (len_trim(argument)>0) kim_config_path = trim(argument)
    call write_config(trim(kim_config_path))
    call kim_initialize(npts, input_r)
    if (len_trim(argument)>0) then
        select case (trim(argument))
        case ('wrong_sign')
            m_vals = [6, 7]
        case ('outside')
            m_vals = [-20, -21]
        case default
            error stop 'unknown test argument'
        end select
        I_par_toroidal = 0.0d0
        call kim_run_for_all_modes()
        error stop 'invalid signed mode was accepted'
    end if
    I_par_toroidal = 0.25d0
    kim_current_relaxation = 1.0d0

    do scenario = 1, 2
        ! Same q profile: the second pair has closer rational surfaces.
        if (scenario == 1) then
            m_vals = [-6, -7]
            n_vals = 2
        else
            m_vals = [-12, -13]
            n_vals = 4
        end if
        call kim_update_profiles()
        call capture(reference)
        if (scenario == 1) then
            call require(.not. any(reference%weights(:,1)*reference%weights(:,2)>0), &
                'separated fixture has overlapping supports')
        else
            call require(any(reference%weights(:,1)*reference%weights(:,2)>0), &
                'overlap fixture has no jointly supported grid points')
        end if
        call capture(repeated)
        call compare(repeated, reference, [1,2], 'unchanged batch')
        m_vals = m_vals(2:1:-1)
        n_vals = n_vals(2:1:-1)
        call capture(reordered)
        call compare(reordered, reference, [2,1], 'mode permutation')
        m_vals = m_vals(2:1:-1)
        n_vals = n_vals(2:1:-1)

        ! Probe each profile separately; prescribed Br is never the feedback observable.
        do profile = 1, 5
            params_b = saved_params
            Ercov = saved_er
            q = saved_q
            select case (profile)
            case (1)
                params_b(1,:) = 1.2d0*saved_params(1,:)
            case (2)
                params_b(3,:) = 0.9d0*saved_params(3,:)
            case (3)
                params_b(4,:) = 1.1d0*saved_params(4,:)
            case (4)
                Ercov = 0.8d0*saved_er
            case (5)
                q = saved_q+0.1d0
            end select
            call kim_update_profiles()
            call capture(refreshed)
            call require(all(abs(kim_mode_resonance - &
                (2.0d0+58.0d0*(-real(m_vals,8)/n_vals-1.5d0 &
                -merge(0.1d0,0.0d0,profile==5))/3.0d0)) < 1.0d-9), &
                'profile refresh used a stale safety factor or wrong signed resonance')
            do point = 1, 2
                call require(maxval(abs(refreshed%fields(1,:,point) &
                    -reference%fields(1,:,point))) > &
                    1.0d-7*maxval(abs(reference%fields(1,:,point))), &
                    'updated profile did not change the solved Es response of each mode')
            end do
            ! A clean initialization from the updated QL profiles is an independent
            ! oracle for stale background/cache state in the repeated-refresh path.
            call kim_initialize(npts, input_r)
            call capture(fresh)
            call compare(refreshed, fresh, [1,2], 'refresh versus fresh initialization')
        end do
        params_b = saved_params
        Ercov = saved_er
        q = saved_q
        call kim_update_profiles()
        call capture(repeated)
        call compare(repeated, reference, [1,2], 'restored profiles')
    end do

    ! Shrink and grow the mode list: no response or metadata from a removed mode survives.
    dim_mn = 1
    m_vals = [-12]
    n_vals = [4]
    call capture(single)
    call compare(single, reference, [1], 'shrunk mode list')
    dim_mn = 2
    m_vals = [-12, -13]
    n_vals = [4, 4]
    call capture(repeated)
    call compare(repeated, reference, [1,2], 'grown mode list')
    print *, 'PASS: signed multimode ordering, supports, and independent profile feedback'

contains
    subroutine capture(output)
        type(response_t), intent(out) :: output
        integer :: mode
        call kim_run_for_all_modes()
        call require(all(kim_mode_m==m_vals) .and. all(kim_mode_n==n_vals), &
            'signed mode identity was rewritten')
        call require(all(kim_mode_status==0), 'mode solve failed')
        call require(size(kim_current_records)==dim_mn, 'stale current record count')
        allocate(output%fields(7,npts,dim_mn), output%currents(2,npts,dim_mn))
        do mode = 1, dim_mn
            call kim_get_wave_fields(mode)
            call kim_get_current_densities(mode)
            output%fields(:,:,mode) = transpose(reshape([Es,Ep,Er,Et,Ez,Br,Bp],[npts,7]))
            output%currents(1,:,mode) = Jpe
            output%currents(2,:,mode) = Jpi
            call require(kim_current_records(mode)%active, 'missing mode normalization record')
            call require(kim_current_records(mode)%status==0, 'normalization guard fired')
            call require(maxval(abs(kim_current_records(mode)%core &
                -kim_embedding_metadata(1:2,mode)))<1.0d-12, 'current belongs to wrong window')
        end do
        output%tensor = kim_D_ion_modes
        output%weights = kim_transition_weights
        output%window = kim_embedding_metadata
        output%resonance = kim_mode_resonance
        output%unit = kim_periodic_current_unit
        output%scale = kim_periodic_scale_modes
        call require(all(ieee_is_finite(real(output%fields))), 'nonfinite field')
        call require(all(ieee_is_finite(aimag(output%fields))), 'nonfinite field')
        call require(all(ieee_is_finite(output%tensor)), 'nonfinite tensor')
    end subroutine

    subroutine compare(actual, expected, mapping, context)
        type(response_t), intent(in) :: actual, expected
        integer, intent(in) :: mapping(:)
        character(*), intent(in) :: context
        integer :: mode, other, component, row, column
        real(8), parameter :: tol=2.0d-8
        do mode = 1, size(mapping)
            other = mapping(mode)
            do component = 1, 7
                call require(maxval(abs(actual%fields(component,:,mode) &
                    -expected%fields(component,:,other))) <= tol &
                    *max(1.0d-30,maxval(abs(expected%fields(component,:,other)))), &
                    context//': fields differ')
            end do
            do component = 1, 2
                call require(maxval(abs(actual%currents(component,:,mode) &
                    -expected%currents(component,:,other))) <= tol &
                    *max(1.0d-30,maxval(abs(expected%currents(component,:,other)))), &
                    context//': current differs')
            end do
            do row = 1, 2
                do column = 1, 2
                    call require(maxval(abs(actual%tensor(row,column,:,mode) &
                        -expected%tensor(row,column,:,other))) <= tol &
                        *max(1.0d-30,maxval(abs(expected%tensor(row,column,:,other)))), &
                        context//': tensor differs')
                end do
            end do
            call require(maxval(abs(actual%weights(:,mode)-expected%weights(:,other))) < tol, &
                context//': transition weights differ')
            call require(maxval(abs(actual%window(:,mode)-expected%window(:,other))) < tol, &
                context//': local window differs')
            call require(abs(actual%resonance(mode)-expected%resonance(other)) < tol, &
                context//': resonance differs')
            call require(abs(actual%unit(mode)-expected%unit(other)) <= &
                tol*abs(expected%unit(other)), context//': unit response differs')
            call require(abs(actual%scale(mode)-expected%scale(other)) <= &
                tol*abs(expected%scale(other)), context//': normalization amplitude differs')
        end do
    end subroutine

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine
    subroutine write_config(path)
        character(*), intent(in) :: path
        integer :: unit
        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(A)') &
            '&KIM_CONFIG', &
            ' number_of_ion_species = 1', &
            ' artificial_debye_case = 0', &
            " type_of_run = 'electrostatic_periodic'", &
            " collision_model = 'FokkerPlanck'", &
            " ion_collision_model = 'FokkerPlanck'", &
            ' read_species_from_namelist = .false.', &
            " plasma_type = 'D'", &
            ' turn_off_electrons = .false.', &
            ' turn_off_ions = .false.', &
            ' rescale_density = .false.', &
            ' ion_flr_scale_factor = 1.0', &
            '/', &
            '&WKB_DISPERSION', '/', &
            '&KIM_IO', &
            " profile_location = './'", &
            " output_path = './periodic_multimode_output/'", &
            ' hdf5_input = .false.', &
            ' hdf5_output = .false.', &
            ' log_level = 3', &
            ' data_verbosity = 0', &
            ' calculate_asymptotics = .false.', &
            ' write_diagnostics_dat = .false.', &
            '/', &
            '&KIM_SETUP', &
            ' btor = -18000.0', &
            ' R0 = 165.0', &
            ' m_mode = -6', &
            ' n_mode = 2', &
            ' omega = 0.0', &
            ' spline_base = 1', &
            ' type_br_field = 12', &
            ' collisions_off = .false.', &
            ' set_profiles_constant = 0', &
            ' bc_type = 3', &
            ' mphi_max = 0', &
            ' Br_boundary_re = 1.0', &
            ' Br_boundary_im = 0.0', &
            '/', &
            '&KIM_GRID', &
            " grid_spacing_rg = 'equidistant'", &
            " grid_spacing_xl = 'equidistant'", &
            ' l_space_dim = 32', &
            ' rg_space_dim = 32', &
            " theta_integration = 'GaussLegendre'", &
            ' Larmor_skip_factor = 5', &
            ' gauss_int_nodes_Ntheta = 17', &
            ' gauss_int_nodes_Nx = 31', &
            ' gauss_int_nodes_Nxp = 30', &
            ' r_plas = 50.0', &
            ' r_min = 10.0', &
            ' width_res = 0.5', &
            ' ampl_res = 15.0', &
            '/', &
            '&KIM_PROFILES', '/', &
            '&KIM_PERIODIC', &
            ' periodic_dr_asis_scale = 5.0', &
            ' periodic_dr_tr_scale = 10.0', &
            ' periodic_kmax_scale = 0.4', &
            ' periodic_n_rg = 64', &
            ' periodic_Bparallel_ratio = (0.0, 0.0)', &
            '/'
        close (unit)
    end subroutine write_config
end program test_periodic_multimode_feedback
