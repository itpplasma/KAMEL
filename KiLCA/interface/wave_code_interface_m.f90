!> Fortran interface to KiLCA library functions for the balance code
!> (QL-Balance/KIM), formerly interface/wave_code_interface.{h,cpp}. Callers
!> import native module procedures with compiler-checked argument types and ranks.
module kilca_wave_code_interface_m
    use kilca_legacy_interfaces_m, only: set_core_data_in_core_module
    use kilca_legacy_interfaces_m, only: set_wave_parameters_in_mode_data_module
    use kilca_legacy_interfaces_m, only: clear_all_data_in_mode_data_module
    use kilca_background_data_m, only: get_background_collision_freqs
    use kilca_background_data_m, only: get_background_magnetic_fields
    use kilca_background_settings_m, only: get_background_rtor
    use, intrinsic :: iso_c_binding, only: &
        c_int, c_intptr_t, c_double, c_double_complex, c_ptr, c_loc, c_f_pointer, &
        c_null_ptr
    use constants, only: dp
    use kilca_core_data_m, only: core_data_create_, core_data_destroy_, &
        core_data_calc_and_set_mode_independent_, &
        core_data_calc_and_set_mode_dependent_antenna_interface_, &
        core_data_calc_and_set_mode_dependent_antenna_interface_mn_, &
        core_data_get_dim_, core_data_get_mda_element_, core_data_get_bp_
    use kilca_mode_data_m, only: mode_data_get_wd_, mode_data_get_zone_handle_, &
        mode_data_eval_EB_fields_, mode_data_eval_diss_power_density_, &
        mode_data_eval_current_density_
    use kilca_wave_data_m, only: wave_data_get_m, wave_data_get_n, &
        wave_data_get_olab_re, wave_data_get_olab_im, &
        get_wave_data_obj_omov_re, get_wave_data_obj_omov_im
    use kilca_flre_zone_m, only: activate_fortran_modules_for_zone_, &
        deactivate_fortran_modules_for_zone_, flre_zone_get_flre_order_, flre_zone_get_cp_
    use kilca_cond_profiles_m, only: get_cond_dimx, get_cond_x_ptr, get_cond_k_ptr, get_cond_iks
    use kilca_transforms_m, only: transform_EB_from_cyl_to_rsp
    use kilca_background_data_m, only: background_interp_in_lab_frame, &
        eval_hthz_native => eval_hthz
    implicit none
    private

    public :: calc_wave_code_data_, calc_wave_code_data_for_mode_, clear_wave_code_data_
    public :: get_basic_background_profiles_from_wave_code_
    public :: get_wave_vectors_from_wave_code_
    public :: get_background_magnetic_fields_from_wave_code_
    public :: get_wave_fields_from_wave_code_
    public :: get_diss_power_density_from_wave_code_
    public :: get_antenna_spectrum_dim_, get_antenna_spectrum_numbers_
    public :: get_collision_frequences_from_wave_code_
    public :: get_current_densities_from_wave_code_
    public :: activate_kilca_modules_for_flre_zone_, deactivate_kilca_modules_for_flre_zone_
    public :: get_kilca_conductivity_array_, calc_conductivity_matrices_for_mode_
    public :: get_mode_parameters_, set_wave_parameters_, unset_wave_parameters_
contains

    !> --- entry points used by the antenna-driven (run-the-whole-code) path ---

    subroutine calc_wave_code_data_(cdptr, run_path, pathlength)
        integer(c_intptr_t), intent(out) :: cdptr
        character(len=*), intent(in) :: run_path
        integer(c_int), intent(in) :: pathlength

        cdptr = core_data_create_(build_path(run_path, pathlength))
        call set_core_data_in_core_module(cdptr)

        call core_data_calc_and_set_mode_independent_(cdptr)
        call core_data_calc_and_set_mode_dependent_antenna_interface_(cdptr)
    end subroutine calc_wave_code_data_

    subroutine calc_wave_code_data_for_mode_(cdptr, run_path, pathlength, m, n)
        integer(c_intptr_t), intent(out) :: cdptr
        character(len=*), intent(in) :: run_path
        integer(c_int), intent(in) :: pathlength, m, n

        cdptr = core_data_create_(build_path(run_path, pathlength))
        call set_core_data_in_core_module(cdptr)

        call core_data_calc_and_set_mode_independent_(cdptr)
        call core_data_calc_and_set_mode_dependent_antenna_interface_mn_(cdptr, m, n, 0_c_int)
    end subroutine calc_wave_code_data_for_mode_

    subroutine clear_wave_code_data_(cdptr)
        integer(c_intptr_t), intent(in) :: cdptr
        call core_data_destroy_(cdptr)
    end subroutine clear_wave_code_data_

    subroutine get_basic_background_profiles_from_wave_code_(cdptr, dim_r, r, q, n, &
                                                             Ti, Te, Vth, Vz, dPhi0)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_r
        real(c_double), intent(in) :: r(0:dim_r - 1)
        real(c_double), intent(out) :: q(0:dim_r - 1), n(0:dim_r - 1)
        real(c_double), intent(out) :: Ti(0:dim_r - 1), Te(0:dim_r - 1)
        real(c_double), intent(inout) :: Vth(0:dim_r - 1), Vz(0:dim_r - 1), dPhi0(0:dim_r - 1)

        call background_interp_in_lab_frame(dim_r, r, q, n, Ti, Te, Vth, Vz, dPhi0)
    end subroutine get_basic_background_profiles_from_wave_code_

    subroutine get_wave_vectors_from_wave_code_(cdptr, dim_r, r, m, n, ks, kp)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_r, m, n
        real(c_double), intent(in) :: r(0:dim_r - 1)
        real(c_double), intent(out) :: ks(0:dim_r - 1), kp(0:dim_r - 1)
        integer :: i
        real(dp) :: kth, kz, htz(0:1)

        do i = 0, dim_r - 1
            kth = real(m, dp)/r(i)
            kz = real(n, dp)/get_background_rtor()
            call eval_hthz_native(r(i), 0, 0, c_null_ptr, htz)
            kp(i) = htz(0)*kth + htz(1)*kz
            ks(i) = htz(1)*kth - htz(0)*kz
        end do
    end subroutine get_wave_vectors_from_wave_code_

    subroutine get_background_magnetic_fields_from_wave_code_(cdptr, dim_r, r, Bt, &
                                                              Bz, B0)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_r
        real(c_double), intent(in) :: r(0:dim_r - 1)
        real(c_double), intent(out) :: Bt(0:dim_r - 1), Bz(0:dim_r - 1), B0(0:dim_r - 1)
        integer :: i

        do i = 0, dim_r - 1
            call get_background_magnetic_fields(r(i), Bt(i), Bz(i), B0(i))
        end do
    end subroutine get_background_magnetic_fields_from_wave_code_

    subroutine get_collision_frequences_from_wave_code_(cdptr, dim_r, r, nui, nue)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_r
        real(c_double), intent(in) :: r(0:dim_r - 1)
        real(c_double), intent(out) :: nui(0:dim_r - 1), nue(0:dim_r - 1)
        integer :: i

        do i = 0, dim_r - 1
            call get_background_collision_freqs(r(i), nui(i), nue(i))
        end do
    end subroutine get_collision_frequences_from_wave_code_

    subroutine get_wave_fields_from_wave_code_(cdptr, dim_r, r, m, n, &
                                               Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_r, m, n
        real(c_double), intent(in) :: r(0:dim_r - 1)
        complex(c_double_complex), intent(out) :: Er(0:dim_r - 1), Es(0:dim_r - 1), &
                                                  Ep(0:dim_r - 1)
        complex(c_double_complex), intent(out) :: Et(0:dim_r - 1), Ez(0:dim_r - 1)
        complex(c_double_complex), intent(out) :: Br(0:dim_r - 1), Bs(0:dim_r - 1), &
                                                  Bp(0:dim_r - 1)
        complex(c_double_complex), intent(out) :: Bt(0:dim_r - 1), Bz(0:dim_r - 1)
        integer :: i, num, ind, comp
        real(c_double) :: EBflat(12)
        complex(c_double_complex) :: EBcyl(6), EBrsp(6)

        num = find_mode_index(cdptr, m, n)
        if (num == -1) then
            write (*, '(/,a,i0,a,i0,a)') &
                'warning: get_wave_fields_from_wave_code: (', m, ', ', n, &
                ') mode is not in the spectrum.'
            return
        end if

        do i = 0, dim_r - 1
            call mode_data_eval_EB_fields_(core_data_get_mda_element_(cdptr, num), r(i), EBflat)
            do comp = 1, 6
                EBcyl(comp) = cmplx(EBflat(2*comp - 1), EBflat(2*comp), c_double_complex)
            end do
            call transform_EB_from_cyl_to_rsp(core_data_get_bp_(cdptr), r(i), EBcyl, EBrsp)

            Er(i) = EBcyl(1)
            Es(i) = EBrsp(2)
            Ep(i) = EBrsp(3)
            Et(i) = EBcyl(2)
            Ez(i) = EBcyl(3)
            Br(i) = EBcyl(4)
            Bs(i) = EBrsp(5)
            Bp(i) = EBrsp(6)
            Bt(i) = EBcyl(5)
            Bz(i) = EBcyl(6)
        end do
    end subroutine get_wave_fields_from_wave_code_

    subroutine get_diss_power_density_from_wave_code_(cdptr, dim_r, r, m, n, &
                                                      ttype, spec, pdis)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_r, m, n, ttype, spec
        real(c_double), intent(in) :: r(0:dim_r - 1)
        real(c_double), intent(out) :: pdis(0:dim_r - 1)
        integer :: i, num
        integer(c_intptr_t) :: mdh

        num = find_mode_index(cdptr, m, n)
        if (num == -1) then
            !> Oracle copy-paste quirk preserved: this warning literally
            !> says "get_wave_fields_from_wave_code" in the C++ source too.
            write (*, '(/,a,i0,a,i0,a)') &
                'warning: get_wave_fields_from_wave_code: (', m, ', ', n, &
                ') mode is not in the spectrum.'
            return
        end if

        mdh = core_data_get_mda_element_(cdptr, num)
        do i = 0, dim_r - 1
            call mode_data_eval_diss_power_density_(mdh, r(i), ttype, spec, pdis(i))
        end do
    end subroutine get_diss_power_density_from_wave_code_

    subroutine get_antenna_spectrum_dim_(cdptr, dim_mn)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(out) :: dim_mn
        dim_mn = core_data_get_dim_(cdptr)
    end subroutine get_antenna_spectrum_dim_

    subroutine get_antenna_spectrum_numbers_(cdptr, dim_mn, m_vals, n_vals)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_mn
        integer(c_int), intent(out) :: m_vals(0:dim_mn - 1), n_vals(0:dim_mn - 1)
        integer :: i
        integer(c_intptr_t) :: wd

        if (dim_mn /= core_data_get_dim_(cdptr)) then
            write (*, '(a)') 'warning: get_antenna_spectrum_numbers: spectrum dimensions must match'
            return
        end if

        do i = 0, dim_mn - 1
            wd = mode_data_get_wd_(core_data_get_mda_element_(cdptr, i))
            m_vals(i) = wave_data_get_m(intptr_to_cptr(wd))
            n_vals(i) = wave_data_get_n(intptr_to_cptr(wd))
        end do
    end subroutine get_antenna_spectrum_numbers_

    subroutine get_current_densities_from_wave_code_(cdptr, dim_r, r, m, n, &
                                                     Jri, Jsi, Jpi, Jre, Jse, Jpe)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: dim_r, m, n
        real(c_double), intent(in) :: r(0:dim_r - 1)
        complex(c_double_complex), intent(out) :: Jri(0:dim_r - 1), &
                                                  Jsi(0:dim_r - 1), Jpi(0:dim_r - 1)
        complex(c_double_complex), intent(out) :: Jre(0:dim_r - 1), &
                                                  Jse(0:dim_r - 1), Jpe(0:dim_r - 1)
        integer :: i, num
        real(c_double) :: current(2)
        integer(c_intptr_t) :: mdh

        num = find_mode_index(cdptr, m, n)
        if (num == -1) then
            write (*, '(/,a,i0,a,i0,a)') &
                'warning: get_current_densities_from_wave_code: (', m, ', ', n, &
                ') mode is not in the spectrum.'
            return
        end if

        mdh = core_data_get_mda_element_(cdptr, num)
        do i = 0, dim_r - 1
            call mode_data_eval_current_density_(mdh, r(i), 0, 0, 0, current)
            Jri(i) = cmplx(current(1), current(2), c_double_complex)
            call mode_data_eval_current_density_(mdh, r(i), 0, 0, 1, current)
            Jsi(i) = cmplx(current(1), current(2), c_double_complex)
            call mode_data_eval_current_density_(mdh, r(i), 0, 0, 2, current)
            Jpi(i) = cmplx(current(1), current(2), c_double_complex)
            call mode_data_eval_current_density_(mdh, r(i), 0, 1, 0, current)
            Jre(i) = cmplx(current(1), current(2), c_double_complex)
            call mode_data_eval_current_density_(mdh, r(i), 0, 1, 1, current)
            Jse(i) = cmplx(current(1), current(2), c_double_complex)
            call mode_data_eval_current_density_(mdh, r(i), 0, 1, 2, current)
            Jpe(i) = cmplx(current(1), current(2), c_double_complex)
        end do
    end subroutine get_current_densities_from_wave_code_

    !> --- FLRE conductivity / QL-Balance coupling ---

    subroutine activate_kilca_modules_for_flre_zone_(cdptr)
        integer(c_intptr_t), intent(in) :: cdptr
        call activate_fortran_modules_for_zone_(mode_data_get_zone_handle_( &
            core_data_get_mda_element_(cdptr, 0), 0))
    end subroutine activate_kilca_modules_for_flre_zone_

    subroutine deactivate_kilca_modules_for_flre_zone_(cdptr)
        integer(c_intptr_t), intent(in) :: cdptr
        call deactivate_fortran_modules_for_zone_(mode_data_get_zone_handle_( &
            core_data_get_mda_element_(cdptr, 0), 0))
    end subroutine deactivate_kilca_modules_for_flre_zone_

    subroutine get_kilca_conductivity_array_(cdptr, m, n, zone_ind, spec, flreo, &
                                             dimv, &
                                             rptr, cptr)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: m, n, zone_ind, spec
        integer(c_int), intent(out) :: flreo, dimv
        type(c_ptr), intent(out) :: rptr, cptr
        integer :: num
        integer(c_intptr_t) :: zone, zone_cp
        integer(c_int) :: iks_off
        real(c_double), pointer :: kflat(:)
        type(c_ptr) :: kbase

        num = find_mode_index(cdptr, m, n)
        if (num == -1) then
            write (*, '(/,a,i0,a,i0,a)') &
                'warning: get_kilca_conductivity_array: (', m, ', ', n, &
                ') mode is not in the spectrum.'
            return
        end if

        zone = mode_data_get_zone_handle_(core_data_get_mda_element_(cdptr, num), zone_ind)
        zone_cp = flre_zone_get_cp_(zone)

        flreo = flre_zone_get_flre_order_(zone)
        dimv = get_cond_dimx(zone_cp)
        rptr = get_cond_x_ptr(zone_cp)
        iks_off = get_cond_iks(zone_cp, spec, 0, 0, 0, 0, 0, 0, 0)

        kbase = get_cond_k_ptr(zone_cp)
        call c_f_pointer(kbase, kflat, [iks_off + 1])
        cptr = c_loc(kflat(iks_off + 1))
    end subroutine get_kilca_conductivity_array_

    subroutine calc_conductivity_matrices_for_mode_(cdptr, run_path, pathlength, m, &
                                                    n)
        integer(c_intptr_t), intent(out) :: cdptr
        character(len=*), intent(in) :: run_path
        integer(c_int), intent(in) :: pathlength, m, n

        cdptr = core_data_create_(build_path(run_path, pathlength))
        call set_core_data_in_core_module(cdptr)

        call core_data_calc_and_set_mode_independent_(cdptr)
        call core_data_calc_and_set_mode_dependent_antenna_interface_mn_(cdptr, m, n, 1_c_int)
    end subroutine calc_conductivity_matrices_for_mode_

    subroutine get_mode_parameters_(cdptr, m, n, kz, omega_mov_re, omega_mov_im)
        use kilca_wave_data_m, only: get_wave_data_obj_omov_im
        use kilca_wave_data_m, only: get_wave_data_obj_omov_re
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: m, n
        real(c_double), intent(out) :: kz, omega_mov_re, omega_mov_im
        integer :: num
        integer(c_intptr_t) :: wd

        num = find_mode_index(cdptr, m, n)
        if (num == -1) then
            !> Oracle copy-paste quirk preserved: this warning literally
            !> says "get_kilca_conductivity_array" in the C++ source too.
            write (*, '(/,a,i0,a,i0,a)') &
                'warning: get_kilca_conductivity_array: (', m, ', ', n, &
                ') mode is not in the spectrum.'
            return
        end if

        wd = mode_data_get_wd_(core_data_get_mda_element_(cdptr, num))
        kz = real(wave_data_get_n(intptr_to_cptr(wd)), dp)/get_background_rtor()
        omega_mov_re = get_wave_data_obj_omov_re(intptr_to_cptr(wd))
        omega_mov_im = get_wave_data_obj_omov_im(intptr_to_cptr(wd))
    end subroutine get_mode_parameters_

    subroutine set_wave_parameters_(cdptr, m, n)
        use kilca_wave_data_m, only: get_wave_data_obj_omov_im
        use kilca_wave_data_m, only: get_wave_data_obj_omov_re
        use kilca_wave_data_m, only: wave_data_get_olab_im
        use kilca_wave_data_m, only: wave_data_get_olab_re
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: m, n
        integer :: num
        integer(c_intptr_t) :: wd
        integer(c_int) :: mm, nn
        real(c_double) :: olab_re, olab_im, omov_re, omov_im

        num = find_mode_index(cdptr, m, n)
        if (num == -1) then
            write (*, '(/,a,i0,a,i0,a)') &
                'warning: set_wave_parameters: (', m, ', ', n, &
                ') mode is not in the spectrum.'
            return
        end if

        wd = mode_data_get_wd_(core_data_get_mda_element_(cdptr, num))
        mm = wave_data_get_m(intptr_to_cptr(wd))
        nn = wave_data_get_n(intptr_to_cptr(wd))
        olab_re = wave_data_get_olab_re(intptr_to_cptr(wd))
        olab_im = wave_data_get_olab_im(intptr_to_cptr(wd))
        omov_re = get_wave_data_obj_omov_re(intptr_to_cptr(wd))
        omov_im = get_wave_data_obj_omov_im(intptr_to_cptr(wd))

        call set_wave_parameters_in_mode_data_module(mm, nn, olab_re, olab_im, omov_re, omov_im)
    end subroutine set_wave_parameters_

    subroutine unset_wave_parameters_()
        call clear_all_data_in_mode_data_module()
    end subroutine unset_wave_parameters_

    !> --- internal helpers ---

    !> Linear scan for the antenna-spectrum index of mode (m, n), mirroring
    !> the oracle's identical inline loop repeated in 6 of the entry points
    !> above. Returns a 0-based index, or -1 if (m, n) is not in cd's
    !> spectrum.
    function find_mode_index(cdptr, mval, nval) result(num)
        integer(c_intptr_t), intent(in) :: cdptr
        integer(c_int), intent(in) :: mval, nval
        integer :: num
        integer :: i, dimv
        integer(c_intptr_t) :: wd

        num = -1
        dimv = core_data_get_dim_(cdptr)
        do i = 0, dimv - 1
            wd = mode_data_get_wd_(core_data_get_mda_element_(cdptr, i))
            if (wave_data_get_m(intptr_to_cptr(wd)) == mval .and. &
                wave_data_get_n(intptr_to_cptr(wd)) == nval) then
                num = i
                return
            end if
        end do
    end function find_mode_index

    function intptr_to_cptr(h) result(cp)
        integer(c_intptr_t), intent(in) :: h
        type(c_ptr) :: cp
        cp = transfer(h, cp)
    end function intptr_to_cptr

    !> Reproduces calc_wave_code_data_'s C string handling exactly:
    !> `strcpy(path, run_path); path[*pathlength] = '\0';
    !> if (path[strlen(path)-1] != '/') strcat(path, "/");` - run_path is
    !> NOT null-terminated by the caller, pathlength is its true length.
    function build_path(run_path, pathlength) result(path)
        character(len=*), intent(in) :: run_path
        integer(c_int), intent(in) :: pathlength
        character(len=1024) :: path
        integer :: i

        path = ''
        do i = 1, pathlength
            path(i:i) = run_path(i:i)
        end do
        if (path(pathlength:pathlength) /= '/') then
            path(pathlength + 1:pathlength + 1) = '/'
        end if
    end function build_path

end module kilca_wave_code_interface_m
