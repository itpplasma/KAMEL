!> Per-mode wave description (harmonics, frequency, resonance location,
!> stitching determinant), formerly the C++ wave_data class. Per-instance
!> handle (same pattern as kilca_cond_profiles_m/kilca_background_data_m's
!> sibling modules): one wave_data per mode_data, shared by all zones of
!> that mode. mode_data (mode.cpp, still C++ until S6) owns the lifetime;
!> the legacy mode_data Fortran module (mode_m.f90) and the new zone_t
!> hierarchy both hold the same opaque handle, exactly as they held the raw
!> C++ pointer before.
module kilca_wave_data_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_ptr, &
        c_loc, c_f_pointer
    implicit none
    private

    public :: wave_data_t
    public :: wave_data_create, wave_data_destroy
    public :: wave_data_get_m, wave_data_get_n
    public :: wave_data_get_olab_re, wave_data_get_olab_im
    public :: wave_data_get_r_res, wave_data_set_r_res
    public :: wave_data_get_det_re, wave_data_get_det_im
    public :: get_wave_data_obj_omov_re, get_wave_data_obj_omov_im
    public :: set_det_in_wd_struct

    type :: wave_data_t
        integer(c_int) :: m, n
        complex(c_double) :: olab, omov, det
        real(c_double) :: r_res
    end type wave_data_t

contains

    function wave_data_create(m, n, olab_re, olab_im, omov_re, omov_im) &
        result(handle) bind(C, name="wave_data_create_")
        integer(c_int), value :: m, n
        real(c_double), value :: olab_re, olab_im, omov_re, omov_im
        integer(c_intptr_t) :: handle
        type(wave_data_t), pointer :: wd

        allocate (wd)
        wd%m = m
        wd%n = n
        wd%olab = cmplx(olab_re, olab_im, c_double)
        wd%omov = cmplx(omov_re, omov_im, c_double)
        wd%r_res = 0.0d0
        wd%det = (0.0d0, 0.0d0)

        handle = transfer(c_loc(wd), handle)
    end function wave_data_create

    subroutine wave_data_destroy(handle) bind(C, name="wave_data_destroy_")
        integer(c_intptr_t), value :: handle
        type(wave_data_t), pointer :: wd

        if (handle == 0_c_intptr_t) return
        call handle_to_wd(handle, wd)
        deallocate (wd)
    end subroutine wave_data_destroy

    function wave_data_get_m(handle) result(res) bind(C, name="wave_data_get_m_")
        type(c_ptr), value :: handle
        integer(c_int) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = wd%m
    end function wave_data_get_m

    function wave_data_get_n(handle) result(res) bind(C, name="wave_data_get_n_")
        type(c_ptr), value :: handle
        integer(c_int) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = wd%n
    end function wave_data_get_n

    function wave_data_get_olab_re(handle) result(res) bind(C, name="wave_data_get_olab_re_")
        type(c_ptr), value :: handle
        real(c_double) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = real(wd%olab, c_double)
    end function wave_data_get_olab_re

    function wave_data_get_olab_im(handle) result(res) bind(C, name="wave_data_get_olab_im_")
        type(c_ptr), value :: handle
        real(c_double) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = aimag(wd%olab)
    end function wave_data_get_olab_im

    function wave_data_get_r_res(handle) result(res) bind(C, name="wave_data_get_r_res_")
        type(c_ptr), value :: handle
        real(c_double) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = wd%r_res
    end function wave_data_get_r_res

    subroutine wave_data_set_r_res(handle, val) bind(C, name="wave_data_set_r_res_")
        type(c_ptr), value :: handle
        real(c_double), value :: val
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        wd%r_res = val
    end subroutine wave_data_set_r_res

    function wave_data_get_det_re(handle) result(res) bind(C, name="wave_data_get_det_re_")
        type(c_ptr), value :: handle
        real(c_double) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = real(wd%det, c_double)
    end function wave_data_get_det_re

    function wave_data_get_det_im(handle) result(res) bind(C, name="wave_data_get_det_im_")
        type(c_ptr), value :: handle
        real(c_double) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = aimag(wd%det)
    end function wave_data_get_det_im

    !> Pre-existing names/by-value c_ptr convention, called from Fortran
    !> (kilca_cond_profiles_m's per-point fallback path).
    function get_wave_data_obj_omov_re(handle) result(res) &
        bind(C, name="get_wave_data_obj_omov_re_")
        type(c_ptr), value :: handle
        real(c_double) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = real(wd%omov, c_double)
    end function get_wave_data_obj_omov_re

    function get_wave_data_obj_omov_im(handle) result(res) &
        bind(C, name="get_wave_data_obj_omov_im_")
        type(c_ptr), value :: handle
        real(c_double) :: res
        type(wave_data_t), pointer :: wd
        call c_f_pointer(handle, wd)
        res = aimag(wd%omov)
    end function get_wave_data_obj_omov_im

    !> Pre-existing name, BY-REFERENCE handle (Fortran caller convention):
    !> the legacy stitching-equations Fortran (calc_system_determinant_)
    !> calls this with the handle held in mode_m.f90's wd_ptr.
    subroutine set_det_in_wd_struct(wd_ptr, re_det, im_det) &
        bind(C, name="set_det_in_wd_struct_")
        integer(c_intptr_t), intent(in) :: wd_ptr
        real(c_double), intent(in) :: re_det, im_det
        type(wave_data_t), pointer :: wd
        call handle_to_wd(wd_ptr, wd)
        wd%det = cmplx(re_det, im_det, c_double)
    end subroutine set_det_in_wd_struct

    subroutine handle_to_wd(handle, wd)
        integer(c_intptr_t), value :: handle
        type(wave_data_t), pointer, intent(out) :: wd
        type(c_ptr) :: cwd
        cwd = transfer(handle, cwd)
        call c_f_pointer(cwd, wd)
    end subroutine handle_to_wd

end module kilca_wave_data_m
