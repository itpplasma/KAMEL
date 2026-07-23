!> Vacuum/homogeneous-medium zone, formerly the C++ hmedium_zone class
!> (hmedium_zone.{h,cpp}). Simplest concrete zone_t extension: a closed-form
!> basis (eval_basis_in_hom_media, pre-existing legacy Fortran), and no
!> dispersion/quantities support (those methods are no-ops or error stubs
!> in the oracle, preserved as-is).
module kilca_hmedium_zone_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, c_ptr
    use constants, only: dp
    use kilca_wave_data_m, only: wave_data_t
    use kilca_zone_m, only: zone_t, zone_register, zone_read, zone_print, skip_line, &
        read_real_before_hash, read_int_before_hash, read_complex_before_hash
    implicit none
    private

    public :: hmedium_zone_t
    public :: hmedium_zone_create

    type, extends(zone_t) :: hmedium_zone_t
        complex(dp) :: sigma = (0, 0)
    contains
        procedure :: read_settings => hmedium_read_settings
        procedure :: print_settings => hmedium_print_settings
        procedure :: calc_basis_fields => hmedium_calc_basis_fields
        procedure :: copy_E_and_B_fields => hmedium_copy_E_and_B_fields
        procedure :: calc_final_fields => hmedium_calc_final_fields
        procedure :: calc_dispersion => hmedium_calc_dispersion
        procedure :: save_dispersion => hmedium_save_dispersion
        procedure :: calc_all_quants => hmedium_calc_all_quants
        procedure :: save_all_quants => hmedium_save_all_quants
        procedure :: eval_diss_power_density => hmedium_eval_diss_power_density
        procedure :: eval_current_density => hmedium_eval_current_density
    end type hmedium_zone_t

    interface
        subroutine eval_basis_in_hom_media(m, kz, omega, sigma, rval, EB)
            import :: dp
            integer, intent(in) :: m
            real(dp), intent(in) :: kz
            complex(dp), intent(in) :: omega, sigma
            real(dp), intent(in) :: rval
            complex(dp), intent(out) :: EB(6, 4)
        end subroutine eval_basis_in_hom_media

        function get_background_rtor() bind(C, name="get_background_rtor_") result(rtor)
            import :: c_double
            real(c_double) :: rtor
        end function get_background_rtor
    end interface

contains

    function hmedium_zone_create(sd_ptr, bp_ptr, wd_handle, path, index_p) &
        result(handle) bind(C, name="hmedium_zone_create_")
        integer(c_intptr_t), value :: sd_ptr, bp_ptr, wd_handle
        character(kind=c_char), intent(in) :: path(*)
        integer(c_int), value :: index_p
        integer(c_intptr_t) :: handle

        type(hmedium_zone_t), pointer :: hz
        class(zone_t), pointer :: zp
        type(c_ptr) :: wd_cptr

        allocate (hz)
        hz%bp = bp_ptr
        hz%index = int(index_p)
        hz%path = zone_c_string(path)
        wd_cptr = transfer(wd_handle, wd_cptr)
        call c_f_pointer_local(wd_cptr, hz%wd)

        zp => hz
        handle = zone_register(zp)
    end function hmedium_zone_create

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

    subroutine hmedium_read_settings(self, file)
        class(hmedium_zone_t), intent(inout) :: self
        character(len=*), intent(in) :: file
        integer :: unit, ios, k

        call zone_read(self, file)

        open (newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'error: hmedium_zone: read_settings: failed to open file ', trim(file)
            stop 1
        end if

        do k = 1, 8
            call skip_line(unit)
        end do

        call skip_line(unit)
        call read_complex_before_hash(unit, self%sigma)
        call skip_line(unit)

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

        self%Nwaves = 4
        self%Ncomps = 6

        if (self%flag_debug /= 0) call self%print_settings()
    end subroutine hmedium_read_settings

    subroutine hmedium_print_settings(self)
        class(hmedium_zone_t), intent(in) :: self
        call zone_print(self)
        write (*, '(a,es16.8e3,a,es16.8e3,a)') 'medium conductivity: (', &
            real(self%sigma, dp), ', ', aimag(self%sigma), ')'
    end subroutine hmedium_print_settings

    subroutine hmedium_calc_basis_fields(self, flag)
        class(hmedium_zone_t), intent(inout) :: self
        integer, intent(in) :: flag
        integer :: i, m
        real(dp) :: kz, dr

        self%dim = self%max_dim
        if (allocated(self%r)) deallocate (self%r)
        allocate (self%r(self%dim))
        if (allocated(self%basis)) deallocate (self%basis)
        allocate (self%basis(self%Ncomps, self%Nwaves, self%dim))

        m = self%wd%m
        kz = real(self%wd%n, dp)/get_background_rtor()

        dr = (self%r2 - self%r1)/(self%dim - 1)

        do i = 1, self%dim
            self%r(i) = self%r1 + dr*(i - 1)
            call eval_basis_in_hom_media(m, kz, self%wd%olab, self%sigma, self%r(i), self%basis(:, :, i))
        end do

        self%r(self%dim) = self%r2
        call eval_basis_in_hom_media(m, kz, self%wd%olab, self%sigma, self%r(self%dim), self%basis(:, :, self%dim))
    end subroutine hmedium_calc_basis_fields

    subroutine hmedium_copy_E_and_B_fields(self, EB_out)
        class(hmedium_zone_t), intent(in) :: self
        real(dp), intent(out) :: EB_out(*)
        integer :: node, comp

        do node = 0, self%dim - 1
            do comp = 0, 5
                EB_out(2*(comp + 6*node) + 1) = real(self%EB(comp + 1, node + 1), dp)
                EB_out(2*(comp + 6*node) + 2) = aimag(self%EB(comp + 1, node + 1))
            end do
        end do
    end subroutine hmedium_copy_E_and_B_fields

    subroutine hmedium_calc_final_fields(self)
        class(hmedium_zone_t), intent(inout) :: self
    end subroutine hmedium_calc_final_fields

    subroutine hmedium_calc_dispersion(self)
        class(hmedium_zone_t), intent(inout) :: self
    end subroutine hmedium_calc_dispersion

    subroutine hmedium_save_dispersion(self)
        class(hmedium_zone_t), intent(inout) :: self
    end subroutine hmedium_save_dispersion

    subroutine hmedium_calc_all_quants(self)
        class(hmedium_zone_t), intent(inout) :: self
    end subroutine hmedium_calc_all_quants

    subroutine hmedium_save_all_quants(self)
        class(hmedium_zone_t), intent(inout) :: self
    end subroutine hmedium_save_all_quants

    subroutine hmedium_eval_diss_power_density(self, x, ttype, spec, dpd)
        class(hmedium_zone_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(in) :: ttype, spec
        real(dp), intent(out) :: dpd(*)
        write (*, '(a)') 'error: eval_diss_power_density() is not implemented for the hmedium_zone'
    end subroutine hmedium_eval_diss_power_density

    subroutine hmedium_eval_current_density(self, x, ttype, spec, comp, J)
        class(hmedium_zone_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(in) :: ttype, spec, comp
        real(dp), intent(out) :: J(*)
        J(1) = 0.0d0
        J(2) = 0.0d0
    end subroutine hmedium_eval_current_density

end module kilca_hmedium_zone_m
