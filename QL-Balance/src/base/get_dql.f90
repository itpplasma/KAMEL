
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
    use control_mod, only: irf, suppression_mode, misalign_diffusion, type_of_run
    use time_evolution, only: save_prof_time_step, time_ind, br_formfactor, br_vac_res
    use h5mod
    use wave_code_data
    use parallelTools
    use QLBalance_diag, only: i_mn_loop
    use QLBalance_kinds, only: dp
    use PolyLagrangeInterpolation

    implicit none

    integer :: modpernode, imin, imax
    integer :: ipoi, ieq, i_mn, mwind_save
    real(dp), dimension(:), allocatable :: dummy

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

    if (debug_mode) write(*,*) "in det dql 1"

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
        allocate (dummy(npoib))
        do ieq = 1, nbaleqs
            call smooth_array_gauss(npoib, mwind, ddr_params_nl(ieq, :), dummy)
            ddr_params_nl(ieq, :) = dummy
            call smooth_array_gauss(npoib, mwind, params_b(ieq, :), dummy)
            params_b(ieq, :) = dummy
        end do
        deallocate (dummy)
    end if

    if (debug_mode) write(*,*) "in det dql 2"

    ! Compute radial electric field:
    Ercov = sqrt_g_times_B_theta_over_c*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params_nl(1, :)/params_b(1, :) + ddr_params_nl(4, :)) &
            /(Z_i*e_charge)

    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror)

    ! Compute diffusion coefficient matrices:

    call MPI_Comm_size(MPI_COMM_WORLD, np_num, ierror);
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror);
    !sum over modes:
    if (np_num .gt. dim_mn) then
        print *, ' '
        print *, 'Number of processes', np_num, 'is larger than number of modes', dim_mn
        call MPI_finalize(ierror)
        stop
    end if

    modpernode = ceiling(float(dim_mn)/float(np_num))
    imin = modpernode*irank + 1
    imax = min(dim_mn, modpernode*(irank + 1))

    if (debug_mode) write(*,*) "in det dql 3"
    if (irf .eq. 1) call update_background_files(path2profs)
    if (debug_mode) write(*,*) "in det dql 3.1"
    if (irf .eq. 1) call get_wave_code_data(imin, imax)
    if (debug_mode) write(*,*) "in det dql 3.2"
    if (irf .eq. 1) call get_background_magnetic_fields_from_wave_code(flre_cd_ptr(imin), dim_r, r, B0t, B0z, B0)
    if (debug_mode) write(*,*) "in det dql 3.3"
    if (irf .eq. 1) call get_collision_frequences_from_wave_code(flre_cd_ptr(imin), dim_r, r, nui, nue)
    if (debug_mode) write(*,*) "in det dql 4"

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
    do i_mn = imin, imax
        call get_wave_vectors_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                                m_vals(i_mn), n_vals(i_mn), ks, kp)
        call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)
        om_E = ks * c * dPhi0 / B0
        vT_e = sqrt(params_b(3, :)/e_mass)
        vT_i = sqrt(params_b(4, :)/p_mass/am)

        i_mn_loop = i_mn

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
            call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)

            ! caluclate part of perpendicular electric field perturbation that comes from
            if (.not. allocated(coef)) allocate(coef(0:nder,nlagr))
            if (debug_mode) write(*,*) "at r_resonant(i_mn) = ", r_resonant(i_mn)
            call binsrc(rb, 1, npoib, r_resonant(i_mn), ibrabsres)
            if (debug_mode) write(*,*) "binary search found ibrabsres = ", ibrabsres

            call get_ind_Lagr_interp(ibrabsres, ind_begin_interp, ind_end_interp)
            call plag_coeff(nlagr, nder, r_resonant(i_mn), rb(ind_begin_interp:ind_end_interp), coef)
            CALL magnetic_island_width(coef, nder, nlagr, ind_begin_interp, ind_end_interp, m_vals(i_mn), MI_width)

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

        call get_wave_fields_from_wave_code(vac_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Bz, Bz, Bz, Bz, Bz, Br, Bz, Bz, Bz, Bz)

        formfactor = (1.d0, 0.d0)/Br

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
            if (irank .eq. 0) then
                call interp_rb_at_r0(Br, r_resonant(i_mn), br_vac_res(time_ind))
            end if
        end if

        call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Bz, Bz, Bz, Bz, Bz, Br, Bz, Bz, Bz, Bz)
        formfactor = Br * formfactor

        ! todo: interpolate formfactor at resonant surface and write out. This is Brtot/Brvac at the resonant surface
        if (trim(type_of_run) .eq. "TimeEvolution") then !TODO: this is a very bad solution... Make it better
            if (irank .eq. 0) then
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

        call get_current_densities_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                            m_vals(i_mn), n_vals(i_mn), Jri, Jsi, Jpi, Jre, Jse, Jpe)

        call integrate_parallel_current(dim_r, r, Jpe, Jpi, Ipar)
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

    call MPI_Allreduce(dqle11_loc, dqle11, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqle12_loc, dqle12, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqle21_loc, dqle21, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqle22_loc, dqle22, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli11_loc, dqli11, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli12_loc, dqli12, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli21_loc, dqli21, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli22_loc, dqli22, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Barrier(MPI_COMM_WORLD, ierror);
    deallocate (dqle11_loc);
    deallocate (dqle12_loc);
    deallocate (dqle21_loc);
    deallocate (dqle22_loc);
    deallocate (dqli11_loc);
    deallocate (dqli12_loc);
    deallocate (dqli21_loc);
    deallocate (dqli22_loc);
    deallocate (formfactor)
    
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

    if (irank .eq. 0) then
        if (modulo(time_ind, save_prof_time_step) .eq. 0) then
            if (suppression_mode .eqv. .false.) then
                CALL write_fields_currs_transp_coefs_to_h5
                call write_D_one_over_nu_to_h5
            end if
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
