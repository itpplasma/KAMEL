
!> @brief subroutine get_dql. Calculates quasilinear diffusion coefficients.
subroutine get_dql

    use grid_mod, only: nbaleqs, npoib &
                        , deriv_coef &
                        , ipbeg, ipend, rb, reint_coef &
                        , sqrt_g_times_B_theta_over_c, Ercov &
                        , mwind &
                        , dqle11, dqle12, dqle21, dqle22 &
                        , dqli11, dqli12, dqli21, dqli22 &
                        , de11, de12, de21, de22, di11, di12, di21, di22 &
                        , rb_cut_in, rb_cut_out, re_cut_out, rb &
                        , r_resonant, d11_misalign, Es_pert_flux, Ipar
    use plasma_parameters
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, e_mass, ev, rtor, pi, rsepar
    use control_mod, only: irf, suppression_mode, misalign_diffusion, type_of_run, wave_code, &
                          jpar_method
    use time_evolution, only: save_prof_time_step, time_ind, br_formfactor, br_vac_res
    use h5mod
    use wave_code_data
    use kim_wave_code_adapter_m, only: kim_update_profiles, kim_run_for_all_modes, &
        kim_get_wave_fields, kim_get_wave_vectors, kim_vac_Br, kim_Br_modes, &
        kim_get_current_densities
    use QLBalance_diag, only: i_mn_loop
    use QLBalance_kinds, only: dp
    use PolyLagrangeInterpolation

    implicit none

    integer :: ipoi, ieq, i_mn, mwind_save
    real(dp), dimension(:), allocatable :: dummy
    real(dp), dimension(:), allocatable :: row_buffer

    real(dp), dimension(npoib) :: spec_weight
    real(dp) :: weight
    real(dp), dimension(npoib) :: vT_e, vT_i, nu_e, nu_i

    real(dp), DIMENSION(:), ALLOCATABLE :: dqle11_loc
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle12_loc
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle21_loc
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle22_loc
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli11_loc
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli12_loc
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli21_loc
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli22_loc
    complex(dp), DIMENSION(:), ALLOCATABLE :: Es_pert_flux_temp
    complex(dp), dimension(:), allocatable :: formfactor

    ! added variables for interpolation of Brvac
    integer :: ibrabsres, ind_begin_interp, ind_end_interp
    real(dp) :: MI_width

    ! curl(B) current density
    complex(dp), dimension(:), allocatable :: Jpar_curlB
    complex(dp) :: Ipar_conductivity, Ipar_curlB
    ! Temporary B field arrays for curl(B) computation (need fresh FLRE fields)
    complex(dp), dimension(:), allocatable :: Br_flre, Bt_flre, Bz_flre
    complex(dp), dimension(:), allocatable :: zeros_dim_r

    allocate (dqle11_loc(npoib))
    allocate (dqle12_loc(npoib))
    allocate (dqle21_loc(npoib))
    allocate (dqle22_loc(npoib))
    allocate (dqli11_loc(npoib))
    allocate (dqli12_loc(npoib))
    allocate (dqli21_loc(npoib))
    allocate (dqli22_loc(npoib))
    allocate (formfactor(npoib))
    if (.not. allocated(d11_misalign)) allocate (d11_misalign(npoib))
    if (.not. allocated(Es_pert_flux)) allocate (Es_pert_flux(npoib))
    if (.not. allocated(Es_pert_flux_temp)) allocate (Es_pert_flux_temp(npoib))
    allocate(Jpar_curlB(npoib))
    allocate(Br_flre(npoib), Bt_flre(npoib), Bz_flre(npoib))
    allocate(zeros_dim_r(npoib))
    zeros_dim_r = (0.0d0, 0.0d0)

    dqle11_loc = 0.0d0
    dqle12_loc = 0.0d0
    dqle21_loc = 0.0d0
    dqle22_loc = 0.0d0
    dqli11_loc = 0.0d0
    dqli12_loc = 0.0d0
    dqli21_loc = 0.0d0
    dqli22_loc = 0.0d0

    if (irf .eq. 0) then
        return
    elseif (irf .eq. 2) then
        dqle11 = 0.0d0
        dqle12 = 0.0d0
        dqle21 = 0.0d0
        dqle22 = 0.0d0
        dqli11 = 0.0d0
        dqli12 = 0.0d0
        dqli21 = 0.0d0
        dqli22 = 0.0d0
        return
    end if

    do ipoi = 1, npoib
        do ieq = 1, nbaleqs
            ! radial derivatives of equilibrium parameters at cell boundaries:
            ddr_params_nl(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
            ! equilibrium parameters at cell boundaries:
            params_b(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
        end do
    end do

    ! Smooth input for KILCA
    if (.true.) then
        allocate (dummy(npoib), row_buffer(npoib))
        do ieq = 1, nbaleqs
            row_buffer = ddr_params_nl(ieq, :)
            call smooth_array_gauss(npoib, mwind, row_buffer, dummy)
            ddr_params_nl(ieq, :) = dummy
            row_buffer = params_b(ieq, :)
            call smooth_array_gauss(npoib, mwind, row_buffer, dummy)
            params_b(ieq, :) = dummy
        end do
        deallocate (row_buffer)
        deallocate (dummy)
    end if

    ! Compute radial electric field:
    Ercov = sqrt_g_times_B_theta_over_c*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params_nl(1, :)/params_b(1, :) + ddr_params_nl(4, :)) &
            /(Z_i*e_charge)

    ! Compute diffusion coefficient matrices:

    select case (trim(wave_code))
    case ('KiLCA')
        if (irf .eq. 1) call update_background_files(path2profs)
        if (irf .eq. 1) call get_wave_code_data(1, dim_mn)
        if (irf .eq. 1) call get_background_magnetic_fields_from_wave_code(flre_cd_ptr(1), dim_r, r, B0t, B0z, B0)
        if (irf .eq. 1) call get_collision_frequences_from_wave_code(flre_cd_ptr(1), dim_r, r, nui, nue)
    case ('KIM')
        if (irf .eq. 1) call kim_update_profiles()
        if (irf .eq. 1) call kim_run_for_all_modes()
        ! B0 and collision freqs already set in kim_initialize
    case default
        error stop "Unknown wave_code: "//trim(wave_code)
    end select

    !  nu_e=15.4d-6*params_b(1,:)/sqrt(params_b(3,:)/ev)**3            &
    !      *(23.d0-0.5d0*log(params_b(1,:)/(params_b(3,:)/ev)**3))
    !  nu_i=1.d-7*params_b(1,:)/sqrt(params_b(4,:)/ev)**3              &
    !      *(23.d0-0.5d0*log(2.d0*params_b(1,:)/(params_b(4,:)/ev)**3))

    nu_e = nue
    nu_i = nui

    !initialization before summing up over modes:
    dqle11 = 0.0d0
    dqle12 = 0.0d0
    dqle21 = 0.0d0
    dqle22 = 0.0d0
    dqli11 = 0.0d0
    dqli12 = 0.0d0
    dqli21 = 0.0d0
    dqli22 = 0.0e0
    Es_pert_flux = 0.0d0

    !sum over modes:
    do i_mn = 1, dim_mn
        select case (trim(wave_code))
        case ('KiLCA')
            call get_wave_vectors_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                                    m_vals(i_mn), n_vals(i_mn), ks, kp)
            call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                                m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)
        case ('KIM')
            call kim_get_wave_vectors(i_mn)
            call kim_get_wave_fields(i_mn)
        case default
            error stop "Unknown wave_code: "//trim(wave_code)
        end select
        om_E = ks * c * dPhi0 / B0
        vT_e = sqrt(params_b(3, :)/e_mass)
        vT_i = sqrt(params_b(4, :)/p_mass/am)

        i_mn_loop = i_mn

        ! DIAGNOSTIC: compare field magnitudes between wave codes
        write(*,*) ''
        write(*,*) '========== FIELD DIAGNOSTICS (mode ', i_mn, ') =========='
        write(*,*) '  wave_code      = ', trim(wave_code)
        write(*,*) '  m/n            = ', m_vals(i_mn), '/', n_vals(i_mn)
        write(*,*) '  r_resonant     = ', r_resonant(i_mn)
        write(*,*) '  --- Peak field amplitudes ---'
        write(*,*) '  max|Es|        = ', maxval(abs(Es))
        write(*,*) '  max|Br|        = ', maxval(abs(Br))
        write(*,*) '  max|Er|        = ', maxval(abs(Er))
        write(*,*) '  max|Ep|        = ', maxval(abs(Ep))
        write(*,*) '  max|Et|        = ', maxval(abs(Et))
        write(*,*) '  max|Ez|        = ', maxval(abs(Ez))
        write(*,*) '  --- Background quantities ---'
        write(*,*) '  max B0         = ', maxval(B0)
        write(*,*) '  max |kp|       = ', maxval(abs(kp))
        write(*,*) '  max |ks|       = ', maxval(abs(ks))
        write(*,*) '  max |om_E|     = ', maxval(abs(om_E))
        write(*,*) '  max nue        = ', maxval(nue)
        write(*,*) '  max nui        = ', maxval(nui)
        write(*,*) '  --- D_22 input terms (at peak |Es| location) ---'
        write(*,*) '  c^2*max|Es|^2  = ', c**2 * maxval(abs(Es))**2
        write(*,*) '  vTe^2*max|Br|^2= ', maxval(vT_e)**2 * maxval(abs(Br))**2
        write(*,*) '  0.5/(nue*B0^2) = ', 0.5d0/(maxval(nue)*maxval(B0)**2)
        write(*,*) '================================================='
        write(*,*) ''

        ! TODO: add switch to choose calculation of collisionless transport coefficients.
        if (.false.) then
            call calc_transport_coeffs_collisionless(npoib, vT_e, de11, de12, de22)
            de21 = de12
            call calc_transport_coeffs_collisionless(npoib, vT_i, di11, di12, di22)
            di21 = di12
        else
            if (.true.) then
                call calc_transport_coeffs_ornuhl(npoib, vT_e, nu_e, de11, de12, de21, de22)
                call calc_transport_coeffs_ornuhl(npoib, vT_i, nu_i, di11, di12, di21, di22)
            else
                call calc_transport_coeffs_ornuhl_drift(1, npoib, de11, de12, de21, de22)
                call calc_transport_coeffs_ornuhl_drift(2, npoib, di11, di12, di21, di22)
            end if
        end if

        if (misalign_diffusion .eqv. .true.) then
            select case (trim(wave_code))
            case ('KiLCA')
                call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                                m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)
            case ('KIM')
                call kim_get_wave_fields(i_mn)  ! already loaded, but refresh
            case default
                error stop "Unknown wave_code: "//trim(wave_code)
            end select

            ! caluclate part of perpendicular electric field perturbation that comes from
            if (.not. allocated(coef)) allocate(coef(0:nder,nlagr))
            if (debug_mode) write(*,*) "at r_resonant(i_mn) = ", r_resonant(i_mn)
            call binsrc(rb, 1, npoib, r_resonant(i_mn), ibrabsres)
            if (debug_mode) write(*,*) "binary search found ibrabsres = ", ibrabsres

            call get_ind_Lagr_interp(ibrabsres, ind_begin_interp, ind_end_interp)
            call plag_coeff(nlagr, nder, r_resonant(i_mn), rb(ind_begin_interp:ind_end_interp), coef)
            call magnetic_island_width(coef, nder, nlagr, ind_begin_interp, ind_end_interp, m_vals(i_mn), MI_width)

            ! the perturbed flux surfaces
            !Es_pert_flux_temp = (-dPhi0) * Br * (m_vals(i_mn) * rtor**2d0 - n_vals(i_mn) * r**2d0 / qsaf) &
            !/ (B0 * r * rtor * (n_vals(i_mn) + (m_vals(i_mn)) / qsaf))
            Es_pert_flux_temp = (-dPhi0) * Br * ks / (B0 * kp)

            ! cut magnetic island from diffusion
            !do ipoi = 1, npoi
            !    if (r(ipoi) .gt. r_resonant(i_mn) - MI_width/2d0 .and. &
            !    r(ipoi) .lt. r_resonant(i_mn) + MI_width/2d0) then
            !        Es_pert_flux_temp(ipoi) = 0d0
            !    end if
            !end do
            Es_pert_flux = Es_pert_flux + Es_pert_flux_temp
        end if

        select case (trim(wave_code))
        case ('KiLCA')
            call get_wave_fields_from_wave_code(vac_cd_ptr(i_mn), dim_r, r, &
                                                m_vals(i_mn), n_vals(i_mn), Bz, Bz, Bz, Bz, Bz, Br, Bz, Bz, Bz, Bz)
        case ('KIM')
            Br = kim_vac_Br(:, i_mn)
        case default
            error stop "Unknown wave_code: "//trim(wave_code)
        end select
        formfactor = (1.d0, 0.d0)/Br

        ! DIAGNOSTIC: vacuum Br values
        write(*,*) '  --- Vacuum Br ---'
        write(*,*) '  |Br_vac| at r_min    = ', abs(Br(1)), ' at r=', r(1)
        write(*,*) '  |Br_vac| at r_max    = ', abs(Br(dim_r)), ' at r=', r(dim_r)
        write(*,*) '  max|Br_vac|          = ', maxval(abs(Br))

        ! In case that the tmhd code uses double sided Fourier series, it
        ! must be set to 4.0d0, since the factor 2.0d0 should occur in the fields.
        spec_weight = 1.0d0

        do ipoi = 1, npoib
            call localizer(1.d0, rb_cut_out, re_cut_out, r(ipoi), weight)
            spec_weight(ipoi) = spec_weight(ipoi)*weight
            call localizer(-1.d0, rb_cut_in, rb_cut_in, r(ipoi), weight)
            spec_weight(ipoi) = spec_weight(ipoi)*weight
        end do

        if (trim(type_of_run) .eq. "TimeEvolution") then !TODO: this is a very bad solution... Make it better
            if (time_ind>0) then
                call interp_rb_at_r0(Br, r_resonant(i_mn), br_vac_res(time_ind))
            end if
        end if

        select case (trim(wave_code))
        case ('KiLCA')
            call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                                m_vals(i_mn), n_vals(i_mn), Bz, Bz, Bz, Bz, Bz, Br, Bz, Bz, Bz, Bz)
        case ('KIM')
            Br = kim_Br_modes(:, i_mn)
        case default
            error stop "Unknown wave_code: "//trim(wave_code)
        end select
        formfactor = Br * formfactor

        ! DIAGNOSTIC: Br_total and formfactor = Br_total/Br_vac
        write(*,*) '  --- Br_total (self-consistent / full response) ---'
        write(*,*) '  |Br_total| at r_max  = ', abs(Br(dim_r)), ' at r=', r(dim_r)
        write(*,*) '  max|Br_total|        = ', maxval(abs(Br))
        write(*,*) '  |formfactor| at r_max= ', abs(formfactor(dim_r))
        write(*,*) '  max|formfactor|      = ', maxval(abs(formfactor))

        ! todo: interpolate formfactor at resonant surface and write out. This is Brtot/Brvac at the resonant surface
        if (trim(type_of_run) .eq. "TimeEvolution") then !TODO: this is a very bad solution... Make it better
            if (time_ind > 0) then
                call interp_rb_at_r0(formfactor, r_resonant(i_mn), br_formfactor(time_ind))
            end if
        end if

        dqle11_loc = dqle11_loc + de11*spec_weight
        dqle12_loc = dqle12_loc + de12*spec_weight
        dqle21_loc = dqle21_loc + de21*spec_weight
        dqle22_loc = dqle22_loc + de22*spec_weight
        dqli11_loc = dqli11_loc + di11*spec_weight
        dqli12_loc = dqli12_loc + di12*spec_weight
        dqli21_loc = dqli21_loc + di21*spec_weight
        dqli22_loc = dqli22_loc + di22*spec_weight

        ! --- Conductivity-based current (sigma*E) ---
        select case (trim(wave_code))
        case ('KiLCA')
            call get_current_densities_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                m_vals(i_mn), n_vals(i_mn), Jri, Jsi, Jpi, Jre, Jse, Jpe)
        case ('KIM')
            call kim_get_current_densities(i_mn)
        case default
            error stop "Unknown wave_code: "//trim(wave_code)
        end select

        write(*,*) '  --- Ipar (conductivity, sigma*E) ---'
        call integrate_parallel_current(dim_r, r, Jpe, Jpi, r_resonant(i_mn), Ipar_conductivity)

        ! --- curl(B)-based current (Ampere's law) ---
        ! Reload full FLRE B field components (Bt, Bz may have been overwritten above)
        select case (trim(wave_code))
        case ('KiLCA')
            call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br_flre, Bs, Bp, Bt_flre, Bz_flre)
        case ('KIM')
            Br_flre = kim_Br_modes(:, i_mn)
            ! KIM stores Et (=Btheta) and Ez via kim_get_wave_fields
            ! but Bt/Bz are not stored separately — for KIM, the EM solver
            ! does not output cylindrical B components other than Br.
            ! Use the conductivity-based current for KIM.
            Bt_flre = (0.0d0, 0.0d0)
            Bz_flre = (0.0d0, 0.0d0)
        end select

        call compute_jpar_curlB(dim_r, r, Br_flre, Bt_flre, Bz_flre, &
                                B0t, B0z, B0, m_vals(i_mn), n_vals(i_mn), rtor, &
                                Jpar_curlB)
        ! compute_jpar_curlB returns Jpar in c=1 Gaussian units.
        ! Convert to full CGS (multiply by c) so it goes through the
        ! same integration/antenna_factor pipeline as conductivity Jpar.
        Jpar_curlB = Jpar_curlB * c

        write(*,*) '  --- Ipar (curlB, Ampere law) ---'
        call integrate_parallel_current(dim_r, r, Jpar_curlB, zeros_dim_r, &
            r_resonant(i_mn), Ipar_curlB)

        ! Select which Ipar to use for antenna factor computation
        select case (trim(jpar_method))
        case ('curlB')
            Ipar = Ipar_curlB
            write(*,*) '  Using curlB Ipar for antenna factor'
        case ('conductivity')
            Ipar = Ipar_conductivity
            write(*,*) '  Using conductivity Ipar for antenna factor'
        case default
            error stop "Unknown jpar_method: "//trim(jpar_method)
        end select
    end do


    ! calculate diffusion due to misalignment of equipotentials and flux surfaces
    if (misalign_diffusion .eqv. .true.) then
        ! rsepar/rtor is the inverse aspect ratio
        d11_misalign = (16.0d0*sqrt(2.0d0) / (9.0d0 * pi**1.5d0)) * (c * abs(Es_pert_flux + Es) / B0)**2.0d0 &
         * (rsepar / rtor)**(1.5d0) / (0.5d0 * nu_e)
        !write(*,*) B0
        !write(*,*) "r         Es     abs(Es)     nu_e"
        !do ipoi = 1, npoib
        !    write(*,*) r(ipoi), Es(ipoi), abs(Es(ipoi)), nu_e(ipoi)
        !end do
        !dqle11_loc = dqle11_loc + d11_misalign
         !dqle12_loc = dqle12_loc + 3 * d11_misalign
         !dqle21_loc = dqle21_loc + 3 * d11_misalign
         !dqle22_loc = dqle22_loc + 12 * d11_misalign
    end if ! misalign_diffusion .eqv. .true.

    dqle11 = dqle11_loc
    dqle12 = dqle12_loc
    dqle21 = dqle21_loc
    dqle22 = dqle22_loc
    dqli11 = dqli11_loc
    dqli12 = dqli12_loc
    dqli21 = dqli21_loc
    dqli22 = dqli22_loc
    deallocate (dqle11_loc);
    deallocate (dqle12_loc);
    deallocate (dqle21_loc);
    deallocate (dqle22_loc);
    deallocate (dqli11_loc);
    deallocate (dqli12_loc);
    deallocate (dqli21_loc);
    deallocate (dqli22_loc);
    deallocate (formfactor)
    deallocate (Jpar_curlB)
    deallocate (Br_flre, Bt_flre, Bz_flre)
    deallocate (zeros_dim_r)

    call calc_parallel_current_directly
    call calc_ion_parallel_current_directly


    if (.true.) then
        mwind_save = mwind
        mwind = 30
        allocate (dummy(npoib))
        call smooth_array_gauss(npoib, mwind, dqle11, dummy)
        dqle11 = dummy
        call smooth_array_gauss(npoib, mwind, dqle12, dummy)
        dqle12 = dummy
        call smooth_array_gauss(npoib, mwind, dqle21, dummy)
        dqle21 = dummy
        call smooth_array_gauss(npoib, mwind, dqle22, dummy)
        dqle22 = dummy
        mwind = 30
        call smooth_array_gauss(npoib, mwind, dqli11, dummy)
        dqli11 = dummy
        call smooth_array_gauss(npoib, mwind, dqli12, dummy)
        dqli12 = dummy
        call smooth_array_gauss(npoib, mwind, dqli21, dummy)
        dqli21 = dummy
        call smooth_array_gauss(npoib, mwind, dqli22, dummy)
        dqli22 = dummy
        mwind = mwind_save
        deallocate (dummy)
    else
        ! set ion particle flux coefficients to zero
        mwind_save = mwind
        mwind = 30
        allocate (dummy(npoib))
        !call smooth_array_gauss(npoib, mwind, dqle11, dummy)
        dqle11 = 0.d0!dummy
        !call smooth_array_gauss(npoib, mwind, dqle12, dummy)
        dqle12 = 0.d0!dummy
        call smooth_array_gauss(npoib, mwind, dqle21, dummy)
        dqle21 = dummy
        call smooth_array_gauss(npoib, mwind, dqle22, dummy)
        dqle22 = dummy
        mwind = 30
        call smooth_array_gauss(npoib, mwind, dqli12, dummy)
        dqli12 = dummy
        call smooth_array_gauss(npoib, mwind, dqli21, dummy)
        dqli21 = dummy
        call smooth_array_gauss(npoib, mwind, dqli21, dummy)
        dqli21 = dummy
        call smooth_array_gauss(npoib, mwind, dqli22, dummy)
        dqli22 = dummy
        mwind = mwind_save
        deallocate (dummy)
    end if


    if (debug_mode) print *, "Debug: write_fields_currs_transp_coefs_to_h5"

    if (modulo(time_ind, save_prof_time_step) .eq. 0) then
        if (suppression_mode .eqv. .false.) then
            call write_fields_currs_transp_coefs_to_h5
            call write_D_one_over_nu_to_h5
        end if
    end if

    if (debug_mode) write(*,*) "Debug: going out of get_dql"

end subroutine get_dql

subroutine initialize_get_dql

    use control_mod, only: irf

    implicit none

    irf = 2
    call get_dql
    irf = 1

end subroutine


subroutine interp_rb_at_r0(func, r0, func_res)

    ! interpolate a function on the rb grid at r0 using Lagrange interpolation

    use PolyLagrangeInterpolation, only: nder, nlagr, binsrc, get_ind_Lagr_interp, plag_coeff
    use grid_mod, only: rb, npoib
    use QLBalance_kinds, only: dp

    implicit none

    complex(dp), intent(in) :: func(npoib)
    complex(dp), intent(out) :: func_res
    real(dp), intent(in) :: r0
    integer :: ibrabsres, ind_begin_interp, ind_end_interp
    real(dp), dimension(:, :), allocatable :: coef

    if (.not. allocated(coef)) allocate(coef(0:nder,nlagr))

    call binsrc(rb, 1, npoib, r0, ibrabsres)
    call get_ind_Lagr_interp(ibrabsres, ind_begin_interp, ind_end_interp)
    call plag_coeff(nlagr, nder, r0, rb(ind_begin_interp:ind_end_interp), coef)

    func_res = sum(coef(0,:) * func(ind_begin_interp:ind_end_interp))

end subroutine

subroutine get_Brvac(brvac_interp)

    use PolyLagrangeInterpolation
    use grid_mod, only: rb, npoib, r_resonant
    use QLBalance_kinds, only: dp
    use h5mod
    use wave_code_data, only: Br

    implicit none

    integer :: ibrabsres, ind_begin_interp, ind_end_interp
    complex(dp), intent(out) :: brvac_interp

    if (.not. allocated(coef)) allocate(coef(0:nder,nlagr))

    call binsrc(rb, 1, npoib, r_resonant(1), ibrabsres)
    call get_ind_Lagr_interp(ibrabsres, ind_begin_interp, ind_end_interp)
    call plag_coeff(nlagr, nder, r_resonant(1), rb(ind_begin_interp:ind_end_interp), coef)

    brvac_interp = sum(coef(0,:) * abs(Br(ind_begin_interp:ind_end_interp)))

end subroutine

subroutine write_Brvac(brvac_interp)

    use h5mod
    use QLBalance_kinds, only: dp
    use grid_mod, only: npoib, r_resonant
    use wave_code_data, only: Br, r
    use control_mod, only: ihdf5IO

    implicit none

    character(len=1024) :: tempch
    complex(dp), intent(in) :: brvac_interp
    integer :: ipoi

        if (ihdf5IO .eq. 1) then
            if (debug_mode) print *, "Debug: Writing Brvac to hdf5 file"
            if (debug_mode) write(*,*) "Debug: Brvac = ", brvac_interp, " at r_res = ", r_resonant(1)
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)
            !call create_group_if_not_existent("/"//trim(h5_mode_groupname))
            tempch = "/"//trim(h5_mode_groupname)//"/Brvac_res_real"
            if (debug_mode) write(*,*) "Debug: group name: ", trim(tempch)
            CALL h5_add_double_0(h5_id, trim(tempch), real(brvac_interp))

            tempch = "/"//trim(h5_mode_groupname)//"/Brvac_res_imag"
            if (debug_mode) write(*,*) "Debug: group name: ", trim(tempch)
            CALL h5_add_double_0(h5_id, trim(tempch), aimag(brvac_interp))

            if (diagnostics_output) then
                if (debug_mode) write (*, *) "Debug: writing Brvac.dat"
                tempch = "/"//trim(h5_mode_groupname)//"/Brvac.dat"
                CALL h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
                if (h5_exists_log) then
                CALL h5_delete(h5_id, trim(tempch))
                end if

                CALL h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                            H5T_NATIVE_DOUBLE, (/-1, 2/), dataset_id)
                CALL h5_append_double_1(dataset_id, r, 1)
                CALL h5_append_double_1(dataset_id, abs(Br), 2)
            end if

            CALL h5_close(h5_id)
            CALL h5_deinit()

        else
            open (7000, file='Brvac.dat')
            do ipoi = 1, npoib
                write (7000, *) r(ipoi), abs(Br(ipoi))
            end do
            close (7000)
        end if

end subroutine
