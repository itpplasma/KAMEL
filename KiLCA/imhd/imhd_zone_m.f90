!> Ideal-MHD plasma zone, formerly the C++ imhd_zone class (imhd_zone.{h,cpp}).
!> imhd_zone adds no data members over the base zone in the C++ oracle - its
!> basis calculation, dispersion, and quantity hooks are split across
!> kilca_incompressible_m (version 0) and kilca_compressible_flow_m
!> (version 1), both of which operate on a bare class(zone_t) for the same
!> reason. The only extra state this Fortran type carries is own_handle,
!> bookkeeping the zone's own pool handle (the C++ oracle used `this`
!> directly; this port's handle-pool architecture has no equivalent
!> implicit self-reference, so read_settings needs somewhere to recover it
!> for set_imhd_data_module_).
module kilca_imhd_zone_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, c_ptr, &
        c_null_ptr
    use constants, only: dp, pi
    use kilca_wave_data_m, only: wave_data_t
    use kilca_zone_m, only: zone_t, zone_register, zone_read, zone_print, skip_line, &
        read_real_before_hash, read_int_before_hash, handle_to_zone
    use kilca_incompressible_m, only: calculate_basis_incompressible
    use kilca_compressible_flow_m, only: calculate_basis_flow
    implicit none
    private

    public :: imhd_zone_t
    public :: imhd_zone_create
    public :: calc_k_vals_

    complex(dp), parameter :: cmplx_one = (1.0_dp, 0.0_dp)

    type, extends(zone_t) :: imhd_zone_t
        integer(c_intptr_t) :: own_handle = 0_c_intptr_t
    contains
        procedure :: read_settings => imhd_read_settings
        procedure :: print_settings => imhd_print_settings
        procedure :: calc_basis_fields => imhd_calc_basis_fields
        procedure :: copy_E_and_B_fields => imhd_copy_E_and_B_fields
        procedure :: calc_final_fields => imhd_calc_final_fields
        procedure :: calc_dispersion => imhd_calc_dispersion
        procedure :: save_dispersion => imhd_save_dispersion
        procedure :: calc_all_quants => imhd_calc_all_quants
        procedure :: save_all_quants => imhd_save_all_quants
        procedure :: eval_diss_power_density => imhd_eval_diss_power_density
        procedure :: eval_current_density => imhd_eval_current_density
    end type imhd_zone_t

    interface
        function get_background_rtor() bind(C, name="get_background_rtor_") result(rtor)
            import :: c_double
            real(c_double) :: rtor
        end function get_background_rtor

        function get_background_mass(i) bind(C, name="get_background_mass_") result(m)
            import :: c_double, c_int
            integer(c_int), value :: i
            real(c_double) :: m
        end function get_background_mass
    end interface

    interface
        subroutine eval_B0_ht_hz_n0_Vz(rval, bp, R) bind(C, name="eval_B0_ht_hz_n0_Vz")
            import :: c_double, c_ptr
            real(c_double), value :: rval
            type(c_ptr), value :: bp
            real(c_double), intent(out) :: R(*)
        end subroutine eval_B0_ht_hz_n0_Vz
    end interface

    external :: set_imhd_data_module

contains

    function imhd_zone_create(sd_ptr, bp_ptr, wd_handle, path, index_p) &
        result(handle) bind(C, name="imhd_zone_create_")
        integer(c_intptr_t), value :: sd_ptr, bp_ptr, wd_handle
        character(kind=c_char), intent(in) :: path(*)
        integer(c_int), value :: index_p
        integer(c_intptr_t) :: handle

        type(imhd_zone_t), pointer :: iz
        class(zone_t), pointer :: zp
        type(c_ptr) :: wd_cptr

        allocate (iz)
        iz%bp = bp_ptr
        iz%index = int(index_p)
        iz%path = zone_c_string(path)
        wd_cptr = transfer(wd_handle, wd_cptr)
        call c_f_pointer_local(wd_cptr, iz%wd)

        zp => iz
        handle = zone_register(zp)
        iz%own_handle = handle
    end function imhd_zone_create

    subroutine c_f_pointer_local(cptr, wd)
        use, intrinsic :: iso_c_binding, only: c_f_pointer
        type(c_ptr), intent(in) :: cptr
        type(wave_data_t), pointer, intent(out) :: wd
        call c_f_pointer(cptr, wd)
    end subroutine c_f_pointer_local

    function zone_c_string(cstr) result(fstr)
        use, intrinsic :: iso_c_binding, only: c_null_char
        character(kind=c_char), intent(in) :: cstr(*)
        character(len=1024) :: fstr
        integer :: i
        fstr = ''
        i = 0
        do
            if (cstr(i + 1) == c_null_char .or. i >= 1024) exit
            fstr(i + 1:i + 1) = cstr(i + 1)
            i = i + 1
        end do
    end function zone_c_string

    !> Mirrors imhd_zone::read_settings: the base read() (8 lines), then a
    !> fresh open re-skipping those same 8 lines (matching hmedium_zone's
    !> double-open pattern), then Solution/Space out/Debugging settings -
    !> note imhd_zone has no "Medium settings" section. The failed-open
    !> error text below literally says "hmedium_zone" in the C++ oracle too
    !> (a copy-paste artifact carried over unchanged into imhd_zone.cpp).
    subroutine imhd_read_settings(self, file)
        class(imhd_zone_t), intent(inout) :: self
        character(len=*), intent(in) :: file
        integer :: unit, ios, k

        call zone_read(self, file)

        open (newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') &
                'error: hmedium_zone: read_settings: failed to open file ', trim(file)
            stop 1
        end if

        do k = 1, 8
            call skip_line(unit)
        end do

        call skip_line(unit)
        call read_int_before_hash(unit, self%max_dim)
        call read_real_before_hash(unit, self%eps_rel)
        call read_real_before_hash(unit, self%eps_abs)
        call skip_line(unit)

        call skip_line(unit)
        call read_int_before_hash(unit, self%deg)
        call read_real_before_hash(unit, self%reps)
        call read_real_before_hash(unit, self%aeps)
        call read_real_before_hash(unit, self%step)
        call skip_line(unit)

        call skip_line(unit)
        call read_int_before_hash(unit, self%flag_debug)
        call skip_line(unit)

        close (unit)

        self%Nwaves = 2
        self%Ncomps = 8

        call set_imhd_data_module(self%own_handle)

        if (self%flag_debug /= 0) call self%print_settings()
    end subroutine imhd_read_settings

    subroutine imhd_print_settings(self)
        class(imhd_zone_t), intent(in) :: self
        call zone_print(self)
    end subroutine imhd_print_settings

    subroutine imhd_calc_basis_fields(self, flag)
        class(imhd_zone_t), intent(inout) :: self
        integer, intent(in) :: flag

        select case (self%version)
        case (0)
            call calculate_basis_incompressible(self)
        case (1)
            call calculate_basis_flow(self)
        case default
            write (*, '(a)') &
                'error: imhd_zone::calc_basis_fields: unknown code version.'
            stop 1
        end select
    end subroutine imhd_calc_basis_fields

    !> Identical pattern to hmedium_zone_m's copy_E_and_B_fields: only the
    !> first 6 of imhd_zone's 8 state-vector components (Er, Et, Ez, Br, Bt,
    !> Bz) reach the lab-frame output; comps 6,7 ((r*zeta) and its
    !> derivative/pressure proxy) are dropped here, matching the oracle's
    !> `for (comp=0; comp<6; comp++)` loop bound.
    subroutine imhd_copy_E_and_B_fields(self, EB_out)
        class(imhd_zone_t), intent(in) :: self
        real(dp), intent(out) :: EB_out(*)
        integer :: node, comp

        do node = 0, self%dim - 1
            do comp = 0, 5
                EB_out(2*(comp + 6*node) + 1) = real(self%EB(comp + 1, node + 1), dp)
                EB_out(2*(comp + 6*node) + 2) = aimag(self%EB(comp + 1, node + 1))
            end do
        end do
    end subroutine imhd_copy_E_and_B_fields

    subroutine imhd_calc_final_fields(self)
        class(imhd_zone_t), intent(inout) :: self
    end subroutine imhd_calc_final_fields

    subroutine imhd_calc_dispersion(self)
        class(imhd_zone_t), intent(inout) :: self
    end subroutine imhd_calc_dispersion

    subroutine imhd_save_dispersion(self)
        class(imhd_zone_t), intent(inout) :: self
    end subroutine imhd_save_dispersion

    subroutine imhd_calc_all_quants(self)
        class(imhd_zone_t), intent(inout) :: self
    end subroutine imhd_calc_all_quants

    subroutine imhd_save_all_quants(self)
        class(imhd_zone_t), intent(inout) :: self
    end subroutine imhd_save_all_quants

    subroutine imhd_eval_diss_power_density(self, x, ttype, spec, dpd)
        class(imhd_zone_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(in) :: ttype, spec
        real(dp), intent(out) :: dpd(*)
        write (*, '(a)') &
            'error: eval_diss_power_density() is not implemented for the imhd_zone'
    end subroutine imhd_eval_diss_power_density

    subroutine imhd_eval_current_density(self, x, ttype, spec, comp, J)
        class(imhd_zone_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(in) :: ttype, spec, comp
        real(dp), intent(out) :: J(*)
        write (*, '(a)') &
            'error: eval_current_density() is not implemented for the imhd_zone'
    end subroutine imhd_eval_current_density

    !> bind(C, name="calc_k_vals_") so that imhd.f90's existing
    !> calc_k_vals_sub (an unchanged legacy Fortran external call) resolves
    !> to this routine. Inlines the oracle's free `calc_k_vals` (a helper
    !> with no other caller) directly, matching `calc_k_vals_`'s own
    !> C++ body. `E` in the oracle's `kfac = E - kA*kA/(kp*kp)` is
    !> constants.h's `const complex<double> E(1.0, 0.0)` (just the complex
    !> value 1+0i, unrelated to Euler's number or the elementary charge
    !> `e` also defined in that header) - inlined here as cmplx_one.
    subroutine calc_k_vals_(zone_ptr_handle, r, kt, kz, ks, kp, k2, kB, &
        re_kA, im_kA, re_kfac, im_kfac) bind(C, name="calc_k_vals_")
        integer(c_intptr_t), intent(in) :: zone_ptr_handle
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: kt, kz, ks, kp, k2, kB
        real(c_double), intent(out) :: re_kA, im_kA, re_kfac, im_kfac

        class(zone_t), pointer :: zone
        real(dp) :: Rb(5), B0, ht, hz, n0, Vz, VA
        complex(dp) :: kA, kfac

        call handle_to_zone(zone_ptr_handle, zone)

        kt = real(zone%wd%m, dp)/r
        kz = real(zone%wd%n, dp)/get_background_rtor()
        k2 = kt*kt + kz*kz

        call eval_B0_ht_hz_n0_Vz(r, c_null_ptr, Rb)
        B0 = Rb(1); ht = Rb(2); hz = Rb(3); n0 = Rb(4); Vz = Rb(5)

        ks = hz*kt - ht*kz
        kp = ht*kt + hz*kz

        kB = kp*B0

        VA = B0/sqrt(4.0_dp*pi*get_background_mass(0)*n0)

        kA = (zone%wd%olab - kz*Vz)/VA

        kfac = cmplx_one - kA*kA/(kp*kp)

        re_kA = real(kA, dp)
        im_kA = aimag(kA)
        re_kfac = real(kfac, dp)
        im_kfac = aimag(kfac)
    end subroutine calc_k_vals_

end module kilca_imhd_zone_m
