!> Abstract base for the zone hierarchy (zone.h/zone.cpp), a radial
!> interval over which one plasma model (vacuum/medium/imhd/flre) is
!> solved. F2003 type-extension polymorphism: concrete subtypes
!> hmedium_zone_t/imhd_zone_t/flre_zone_t (kilca_hmedium_zone_m/
!> kilca_imhd_zone_m/kilca_flre_zone_m) extend zone_t and override the
!> deferred bindings below.
!>
!> sd (settings*) is write-only at construction in the C++ oracle (zone.cpp
!> never dereferences it; `path` already holds a copy of sd->path2project),
!> so it is accepted by the create shims for ABI/call-site fidelity but not
!> stored. bp (background*) IS read by several subclasses (imhd_zone.cpp,
!> incompressible.cpp, compressible_flow.cpp, flre_zone.cpp all pass
!> `zone->bp` into eval_* calls), but background is now a Fortran singleton
!> (kilca_background_data_m) whose handle value is always ignored, so bp is
!> kept only as a pass-through sentinel for call-site shape fidelity.
!>
!> mode_data (mode.cpp/calc_mode.cpp, still C++ until S6) cannot hold or
!> dispatch through a class(zone_t) directly, so this module also exports a
!> bind(C) dispatch-shim layer: one zone_<method>_ function per
!> deferred/concrete binding, taking an opaque handle, recovering a
!> class(zone_t) pointer via c_f_pointer, and making a normal type-bound
!> call (dispatched by the compiler to the right override). Per-subtype
!> <name>_zone_create_ factories (in each subtype's own module) return that
!> handle. This shim is temporary: S6 removes it once mode_data itself
!> becomes Fortran and can hold class(zone_t) natively.
module kilca_zone_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_loc, c_f_pointer, c_null_char
    use constants, only: dp
    use kilca_wave_data_m, only: wave_data_t
    implicit none
    private

    public :: zone_t
    public :: PLASMA_MODEL_VACUUM, PLASMA_MODEL_MEDIUM, PLASMA_MODEL_IMHD, &
        PLASMA_MODEL_RMHD, PLASMA_MODEL_FLRE
    public :: BOUNDARY_CENTER, BOUNDARY_INFINITY, BOUNDARY_IDEALWALL, &
        BOUNDARY_INTERFACE, BOUNDARY_ANTENNA
    public :: bc_str, med_str, Nbc, Nmed
    public :: zone_read, zone_print
    public :: zone_register, handle_to_zone
    public :: skip_line, read_real_before_hash, read_int_before_hash
    public :: read_token_before_hash, read_complex_before_hash

    public :: zone_read_settings_c, zone_print_settings_c, zone_calc_basis_fields_c
    public :: zone_copy_E_and_B_fields_c, zone_calc_final_fields_c
    public :: zone_calc_dispersion_c, zone_save_dispersion_c
    public :: zone_calc_all_quants_c, zone_save_all_quants_c
    public :: zone_eval_diss_power_density_c, zone_eval_current_density_c
    public :: zone_get_r1_c, zone_get_r2_c
    public :: zone_get_dim_of_basis_c, zone_get_dim_of_basis_vector_c
    public :: zone_get_radial_grid_dimension_c, zone_get_code_version_c
    public :: zone_get_medium_c, zone_get_bc1_c, zone_get_bc2_c
    public :: zone_get_basis_at_left_boundary_c, zone_get_basis_at_right_boundary_c
    public :: zone_save_basis_fields_c, zone_calc_superposition_of_basis_fields_c
    public :: zone_save_final_fields_c, zone_copy_radial_grid_c
    public :: zone_destroy_c
    public :: get_right_boundary_of_zone, get_left_boundary_of_zone

    integer, parameter :: PLASMA_MODEL_VACUUM = 0
    integer, parameter :: PLASMA_MODEL_MEDIUM = 1
    integer, parameter :: PLASMA_MODEL_IMHD = 2
    integer, parameter :: PLASMA_MODEL_RMHD = 3
    integer, parameter :: PLASMA_MODEL_FLRE = 4

    integer, parameter :: BOUNDARY_CENTER = 0
    integer, parameter :: BOUNDARY_INFINITY = 1
    integer, parameter :: BOUNDARY_IDEALWALL = 2
    integer, parameter :: BOUNDARY_INTERFACE = 3
    integer, parameter :: BOUNDARY_ANTENNA = 4

    integer, parameter :: Nbc = 5
    character(len=16), parameter :: bc_str(0:Nbc - 1) = &
        [character(len=16) :: 'center', 'infinity', 'idealwall', 'interface', 'antenna']
    integer, parameter :: Nmed = 5
    character(len=16), parameter :: med_str(0:Nmed - 1) = &
        [character(len=16) :: 'vacuum', 'medium', 'imhd', 'rmhd', 'flre']

    type, abstract :: zone_t
        real(dp) :: r1 = 0, r2 = 0
        integer :: bc1 = -1, bc2 = -1
        integer :: medium = -1, version = 0
        integer :: index = 0
        integer(c_intptr_t) :: bp = 1_c_intptr_t
        type(wave_data_t), pointer :: wd => null()
        character(len=1024) :: path = ''
        integer :: dim = 0
        real(dp), allocatable :: r(:)
        complex(dp), allocatable :: basis(:, :, :)
        complex(dp), allocatable :: EB(:, :)
        complex(dp), allocatable :: S(:)
        integer :: Nwaves = 0, Ncomps = 0
        integer :: max_dim = 0
        real(dp) :: eps_rel = 0, eps_abs = 0
        integer :: deg = 0
        real(dp) :: reps = 0, aeps = 0, step = 0
        integer :: flag_debug = 0
    contains
        procedure :: get_r1 => zone_get_r1
        procedure :: get_r2 => zone_get_r2
        procedure :: get_dim_of_basis => zone_get_dim_of_basis
        procedure :: get_dim_of_basis_vector => zone_get_dim_of_basis_vector
        procedure :: get_code_version => zone_get_code_version
        procedure :: get_radial_grid_dimension => zone_get_radial_grid_dimension
        procedure :: save_basis_fields => zone_save_basis_fields
        procedure :: calc_superposition_of_basis_fields => zone_calc_superposition_of_basis_fields
        procedure :: save_final_fields => zone_save_final_fields
        procedure :: copy_radial_grid => zone_copy_radial_grid

        procedure(read_settings_if), deferred :: read_settings
        procedure(print_settings_if), deferred :: print_settings
        procedure(calc_basis_fields_if), deferred :: calc_basis_fields
        procedure(copy_E_and_B_fields_if), deferred :: copy_E_and_B_fields
        procedure(calc_final_fields_if), deferred :: calc_final_fields
        procedure(calc_dispersion_if), deferred :: calc_dispersion
        procedure(save_dispersion_if), deferred :: save_dispersion
        procedure(calc_all_quants_if), deferred :: calc_all_quants
        procedure(save_all_quants_if), deferred :: save_all_quants
        procedure(eval_diss_power_density_if), deferred :: eval_diss_power_density
        procedure(eval_current_density_if), deferred :: eval_current_density
    end type zone_t

    abstract interface
        subroutine read_settings_if(self, file)
            import :: zone_t
            class(zone_t), intent(inout) :: self
            character(len=*), intent(in) :: file
        end subroutine read_settings_if

        subroutine print_settings_if(self)
            import :: zone_t
            class(zone_t), intent(in) :: self
        end subroutine print_settings_if

        subroutine calc_basis_fields_if(self, flag)
            import :: zone_t
            class(zone_t), intent(inout) :: self
            integer, intent(in) :: flag
        end subroutine calc_basis_fields_if

        subroutine copy_E_and_B_fields_if(self, EB_out)
            import :: zone_t, dp
            class(zone_t), intent(in) :: self
            real(dp), intent(out) :: EB_out(*)
        end subroutine copy_E_and_B_fields_if

        subroutine calc_final_fields_if(self)
            import :: zone_t
            class(zone_t), intent(inout) :: self
        end subroutine calc_final_fields_if

        subroutine calc_dispersion_if(self)
            import :: zone_t
            class(zone_t), intent(inout) :: self
        end subroutine calc_dispersion_if

        subroutine save_dispersion_if(self)
            import :: zone_t
            class(zone_t), intent(inout) :: self
        end subroutine save_dispersion_if

        subroutine calc_all_quants_if(self)
            import :: zone_t
            class(zone_t), intent(inout) :: self
        end subroutine calc_all_quants_if

        subroutine save_all_quants_if(self)
            import :: zone_t
            class(zone_t), intent(inout) :: self
        end subroutine save_all_quants_if

        subroutine eval_diss_power_density_if(self, x, ttype, spec, dpd)
            import :: zone_t, dp
            class(zone_t), intent(in) :: self
            real(dp), intent(in) :: x
            integer, intent(in) :: ttype, spec
            real(dp), intent(out) :: dpd(*)
        end subroutine eval_diss_power_density_if

        subroutine eval_current_density_if(self, x, ttype, spec, comp, J)
            import :: zone_t, dp
            class(zone_t), intent(in) :: self
            real(dp), intent(in) :: x
            integer, intent(in) :: ttype, spec, comp
            real(dp), intent(out) :: J(*)
        end subroutine eval_current_density_if
    end interface

    interface
        subroutine eval_superposition_of_basis_functions(D, Nw, ndim, basis, S, EB)
            import :: dp
            integer, intent(in) :: D, Nw, ndim
            complex(dp), intent(in) :: basis(D, Nw, ndim)
            complex(dp), intent(in) :: S(Nw)
            complex(dp), intent(out) :: EB(D, ndim)
        end subroutine eval_superposition_of_basis_functions
    end interface

    !> c_loc/transfer cannot recover a polymorphic dynamic type from a bare
    !> address (the C pointer carries no type-descriptor information), so
    !> handles are plain 1-based indices into this fixed-size pool instead
    !> of memory addresses. Sized far above any plausible Nzones (single
    !> digits to a few dozen per mode, a handful of modes per run).
    integer, parameter :: max_zones = 4096
    type :: zone_box_t
        class(zone_t), pointer :: z => null()
    end type zone_box_t
    type(zone_box_t) :: zone_pool(max_zones)
    integer :: zone_pool_size = 0

contains

    !> Replicates zone::read: 8 lines (skip, r1, bc1-string, medium-string,
    !> version, bc2-string, r2, skip), each value taken as the text before
    !> '#' (matching read_line_2get_*'s strtok(buf, "#\n\0") convention).
    subroutine zone_read(self, file)
        class(zone_t), intent(inout) :: self
        character(len=*), intent(in) :: file
        integer :: unit, ios, k
        character(len=1024) :: line
        character(len=64) :: bstr1, bstr2, mstr

        open (newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'error: read_settings: failed to open file ', trim(file)
            stop 1
        end if

        read (unit, '(a)') line
        call read_real_before_hash(unit, self%r1)
        call read_token_before_hash(unit, bstr1)
        call read_token_before_hash(unit, mstr)
        call read_int_before_hash(unit, self%version)
        call read_token_before_hash(unit, bstr2)
        call read_real_before_hash(unit, self%r2)
        read (unit, '(a)') line

        close (unit)

        self%bc1 = -1
        self%bc2 = -1
        do k = 0, Nbc - 1
            if (trim(bstr1) == trim(bc_str(k))) self%bc1 = k
            if (trim(bstr2) == trim(bc_str(k))) self%bc2 = k
        end do

        if (self%bc1 == -1 .or. self%bc2 == -1) then
            write (*, '(a,a,a,a)') 'zone::read: BC type is not known: ', trim(bstr1), ' ', trim(bstr2)
            stop 1
        end if

        self%medium = -1
        do k = 0, Nmed - 1
            if (trim(mstr) == trim(med_str(k))) then
                self%medium = k
                exit
            end if
        end do

        if (self%medium == -1) then
            write (*, '(a,a)') 'zone::read: medium type is not known: ', trim(mstr)
            stop 1
        end if
    end subroutine zone_read

    subroutine read_real_before_hash(unit, val)
        integer, intent(in) :: unit
        real(dp), intent(out) :: val
        character(len=1024) :: line
        integer :: ipos
        read (unit, '(a)') line
        ipos = index(line, '#')
        if (ipos == 0) ipos = len_trim(line) + 1
        read (line(1:ipos - 1), *) val
    end subroutine read_real_before_hash

    subroutine read_int_before_hash(unit, val)
        integer, intent(in) :: unit
        integer, intent(out) :: val
        character(len=1024) :: line
        integer :: ipos
        read (unit, '(a)') line
        ipos = index(line, '#')
        if (ipos == 0) ipos = len_trim(line) + 1
        read (line(1:ipos - 1), *) val
    end subroutine read_int_before_hash

    subroutine read_token_before_hash(unit, val)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: val
        character(len=1024) :: line
        integer :: ipos
        read (unit, '(a)') line
        ipos = index(line, '#')
        if (ipos == 0) ipos = len_trim(line) + 1
        read (line(1:ipos - 1), *) val
    end subroutine read_token_before_hash

    !> Matches read_line_2get_complex: the text before '#' is "(re,im)".
    subroutine read_complex_before_hash(unit, val)
        integer, intent(in) :: unit
        complex(dp), intent(out) :: val
        character(len=1024) :: line
        integer :: ipos, p1, p2, p3
        real(dp) :: re, im
        read (unit, '(a)') line
        ipos = index(line, '#')
        if (ipos == 0) ipos = len_trim(line) + 1
        p1 = index(line(1:ipos - 1), '(')
        p2 = index(line(1:ipos - 1), ',')
        p3 = index(line(1:ipos - 1), ')')
        read (line(p1 + 1:p2 - 1), *) re
        read (line(p2 + 1:p3 - 1), *) im
        val = cmplx(re, im, dp)
    end subroutine read_complex_before_hash

    subroutine skip_line(unit)
        integer, intent(in) :: unit
        character(len=1024) :: line
        read (unit, '(a)') line
    end subroutine skip_line

    subroutine zone_print(self)
        class(zone_t), intent(in) :: self
        write (*, '(a,i0)') 'zone index = ', self%index
        write (*, '(a,es24.16e3)') 'r1 = ', self%r1
        write (*, '(a,es24.16e3)') 'r2 = ', self%r2
        write (*, '(a,i0)') 'bc1 = ', self%bc1
        write (*, '(a,i0)') 'bc2 = ', self%bc2
        write (*, '(a,i0)') 'medium = ', self%medium
        write (*, '(a,i0)') 'version = ', self%version
        write (*, '(a,i0)') 'number of waves: ', self%Nwaves
        write (*, '(a,i0)') 'number of field components: ', self%Ncomps
        write (*, '(a,i0)') 'max dimension of the radial grid for the solution: ', self%max_dim
        write (*, '(a,es16.8e3)') 'relative accuracy for the solver: ', self%eps_rel
        write (*, '(a,es16.8e3)') 'absolute accuracy for the solver: ', self%eps_abs
        write (*, '(a,i0)') 'degree of the polynomial used to space out the solution: ', self%deg
        write (*, '(a,es16.8e3)') 'relative accuracy of the sparse solution: ', self%reps
        write (*, '(a,es16.8e3)') 'absolute accuracy of the sparse solution: ', self%aeps
        write (*, '(a,es16.8e3)') 'max grid step in the solution: ', self%step
        write (*, '(a,i0)') 'flag for debugging mode: ', self%flag_debug
    end subroutine zone_print

    real(dp) function zone_get_r1(self) result(res)
        class(zone_t), intent(in) :: self
        res = self%r1
    end function zone_get_r1

    real(dp) function zone_get_r2(self) result(res)
        class(zone_t), intent(in) :: self
        res = self%r2
    end function zone_get_r2

    integer function zone_get_dim_of_basis(self) result(res)
        class(zone_t), intent(in) :: self
        res = self%Nwaves
    end function zone_get_dim_of_basis

    integer function zone_get_dim_of_basis_vector(self) result(res)
        class(zone_t), intent(in) :: self
        res = self%Ncomps
    end function zone_get_dim_of_basis_vector

    integer function zone_get_code_version(self) result(res)
        class(zone_t), intent(in) :: self
        res = self%version
    end function zone_get_code_version

    integer function zone_get_radial_grid_dimension(self) result(res)
        class(zone_t), intent(in) :: self
        res = self%dim
    end function zone_get_radial_grid_dimension

    subroutine zone_calc_superposition_of_basis_fields(self, S_p)
        class(zone_t), intent(inout) :: self
        real(dp), intent(in) :: S_p(0:2*self%Nwaves - 1)
        integer :: i

        if (allocated(self%EB)) deallocate (self%EB)
        allocate (self%EB(self%Ncomps, self%dim))

        if (allocated(self%S)) deallocate (self%S)
        allocate (self%S(self%Nwaves))

        do i = 1, self%Nwaves
            self%S(i) = cmplx(S_p(2*(i - 1)), S_p(2*(i - 1) + 1), dp)
        end do

        call eval_superposition_of_basis_functions(self%Ncomps, self%Nwaves, self%dim, self%basis, self%S, self%EB)
    end subroutine zone_calc_superposition_of_basis_fields

    subroutine zone_copy_radial_grid(self, r_p)
        class(zone_t), intent(in) :: self
        real(dp), intent(out) :: r_p(*)
        integer :: i
        do i = 1, self%dim
            r_p(i) = self%r(i)
        end do
    end subroutine zone_copy_radial_grid

    subroutine zone_save_basis_fields(self, path2linear)
        class(zone_t), intent(in) :: self
        character(len=*), intent(in) :: path2linear
        character(len=1024) :: fname
        integer :: unit, i, sol, comp

        write (fname, '(a,a,i0,a)') trim(path2linear), 'debug-data/zone_', self%index, '_basis.dat'
        open (newunit=unit, file=trim(fname), status='replace', action='write')

        do i = 1, self%dim
            write (unit, '(es24.16e3)', advance='no') self%r(i)
            do sol = 1, self%Nwaves
                do comp = 1, self%Ncomps
                    write (unit, '(a,es24.16e3,a,es24.16e3)', advance='no') &
                        char(9), real(self%basis(comp, sol, i), dp), char(9), aimag(self%basis(comp, sol, i))
                end do
            end do
            write (unit, *)
        end do

        close (unit)
    end subroutine zone_save_basis_fields

    subroutine zone_save_final_fields(self, path2linear)
        class(zone_t), intent(in) :: self
        character(len=*), intent(in) :: path2linear
        character(len=1024) :: fname
        integer :: unit, i, comp

        write (fname, '(a,a,i0,a)') trim(path2linear), 'zone_', self%index, '_EB.dat'
        open (newunit=unit, file=trim(fname), status='replace', action='write')

        do i = 1, self%dim
            write (unit, '(es24.16e3)', advance='no') self%r(i)
            do comp = 1, self%Ncomps
                write (unit, '(a,es24.16e3,a,es24.16e3)', advance='no') &
                    char(9), real(self%EB(comp, i), dp), char(9), aimag(self%EB(comp, i))
            end do
            write (unit, *)
        end do

        close (unit)
    end subroutine zone_save_final_fields

    !> Registers a newly-allocated concrete zone (called by each subtype's
    !> own *_zone_create_ factory) and returns its pool-index handle.
    function zone_register(z) result(handle)
        class(zone_t), pointer, intent(in) :: z
        integer(c_intptr_t) :: handle
        if (zone_pool_size >= max_zones) then
            write (*, '(a)') 'error: zone_register: zone pool exhausted'
            stop 1
        end if
        zone_pool_size = zone_pool_size + 1
        zone_pool(zone_pool_size)%z => z
        handle = int(zone_pool_size, c_intptr_t)
    end function zone_register

    subroutine handle_to_zone(handle, z)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer, intent(out) :: z
        z => zone_pool(int(handle))%z
    end subroutine handle_to_zone

    !> ---- bind(C) dispatch shims for still-C++ callers (mode.cpp,
    !> calc_mode.cpp, wave_code_interface.cpp) ----

    subroutine zone_read_settings_c(handle, file) bind(C, name="zone_read_settings_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char), intent(in) :: file(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%read_settings(c_string_to_fortran(file))
    end subroutine zone_read_settings_c

    subroutine zone_print_settings_c(handle) bind(C, name="zone_print_settings_")
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%print_settings()
    end subroutine zone_print_settings_c

    subroutine zone_calc_basis_fields_c(handle, flag) bind(C, name="zone_calc_basis_fields_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: flag
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%calc_basis_fields(int(flag))
    end subroutine zone_calc_basis_fields_c

    subroutine zone_copy_E_and_B_fields_c(handle, EB_out) bind(C, name="zone_copy_E_and_B_fields_")
        integer(c_intptr_t), value :: handle
        real(c_double), intent(out) :: EB_out(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%copy_E_and_B_fields(EB_out)
    end subroutine zone_copy_E_and_B_fields_c

    subroutine zone_calc_final_fields_c(handle) bind(C, name="zone_calc_final_fields_")
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%calc_final_fields()
    end subroutine zone_calc_final_fields_c

    subroutine zone_calc_dispersion_c(handle) bind(C, name="zone_calc_dispersion_")
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%calc_dispersion()
    end subroutine zone_calc_dispersion_c

    subroutine zone_save_dispersion_c(handle) bind(C, name="zone_save_dispersion_")
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%save_dispersion()
    end subroutine zone_save_dispersion_c

    subroutine zone_calc_all_quants_c(handle) bind(C, name="zone_calc_all_quants_")
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%calc_all_quants()
    end subroutine zone_calc_all_quants_c

    subroutine zone_save_all_quants_c(handle) bind(C, name="zone_save_all_quants_")
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%save_all_quants()
    end subroutine zone_save_all_quants_c

    subroutine zone_eval_diss_power_density_c(handle, x, ttype, spec, dpd) &
        bind(C, name="zone_eval_diss_power_density_")
        integer(c_intptr_t), value :: handle
        real(c_double), value :: x
        integer(c_int), value :: ttype, spec
        real(c_double), intent(out) :: dpd(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%eval_diss_power_density(x, int(ttype), int(spec), dpd)
    end subroutine zone_eval_diss_power_density_c

    subroutine zone_eval_current_density_c(handle, x, ttype, spec, comp, J) &
        bind(C, name="zone_eval_current_density_")
        integer(c_intptr_t), value :: handle
        real(c_double), value :: x
        integer(c_int), value :: ttype, spec, comp
        real(c_double), intent(out) :: J(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%eval_current_density(x, int(ttype), int(spec), int(comp), J)
    end subroutine zone_eval_current_density_c

    real(c_double) function zone_get_r1_c(handle) bind(C, name="zone_get_r1_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%get_r1()
    end function zone_get_r1_c

    real(c_double) function zone_get_r2_c(handle) bind(C, name="zone_get_r2_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%get_r2()
    end function zone_get_r2_c

    integer(c_int) function zone_get_dim_of_basis_c(handle) bind(C, name="zone_get_dim_of_basis_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%get_dim_of_basis()
    end function zone_get_dim_of_basis_c

    integer(c_int) function zone_get_dim_of_basis_vector_c(handle) &
        bind(C, name="zone_get_dim_of_basis_vector_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%get_dim_of_basis_vector()
    end function zone_get_dim_of_basis_vector_c

    integer(c_int) function zone_get_radial_grid_dimension_c(handle) &
        bind(C, name="zone_get_radial_grid_dimension_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%get_radial_grid_dimension()
    end function zone_get_radial_grid_dimension_c

    integer(c_int) function zone_get_code_version_c(handle) bind(C, name="zone_get_code_version_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%get_code_version()
    end function zone_get_code_version_c

    integer(c_int) function zone_get_medium_c(handle) bind(C, name="zone_get_medium_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%medium
    end function zone_get_medium_c

    integer(c_int) function zone_get_bc1_c(handle) bind(C, name="zone_get_bc1_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%bc1
    end function zone_get_bc1_c

    integer(c_int) function zone_get_bc2_c(handle) bind(C, name="zone_get_bc2_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = z%bc2
    end function zone_get_bc2_c

    type(c_ptr) function zone_get_basis_at_left_boundary_c(handle) &
        bind(C, name="zone_get_basis_at_left_boundary_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = c_loc(z%basis(1, 1, 1))
    end function zone_get_basis_at_left_boundary_c

    type(c_ptr) function zone_get_basis_at_right_boundary_c(handle) &
        bind(C, name="zone_get_basis_at_right_boundary_") result(res)
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        res = c_loc(z%basis(1, 1, z%dim))
    end function zone_get_basis_at_right_boundary_c

    subroutine zone_save_basis_fields_c(handle, path2linear) bind(C, name="zone_save_basis_fields_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char), intent(in) :: path2linear(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%save_basis_fields(c_string_to_fortran(path2linear))
    end subroutine zone_save_basis_fields_c

    subroutine zone_calc_superposition_of_basis_fields_c(handle, S_p) &
        bind(C, name="zone_calc_superposition_of_basis_fields_")
        integer(c_intptr_t), value :: handle
        real(c_double), intent(in) :: S_p(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%calc_superposition_of_basis_fields(S_p(1:2*z%Nwaves))
    end subroutine zone_calc_superposition_of_basis_fields_c

    subroutine zone_save_final_fields_c(handle, path2linear) bind(C, name="zone_save_final_fields_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char), intent(in) :: path2linear(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%save_final_fields(c_string_to_fortran(path2linear))
    end subroutine zone_save_final_fields_c

    subroutine zone_copy_radial_grid_c(handle, r_p) bind(C, name="zone_copy_radial_grid_")
        integer(c_intptr_t), value :: handle
        real(c_double), intent(out) :: r_p(*)
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        call z%copy_radial_grid(r_p)
    end subroutine zone_copy_radial_grid_c

    subroutine zone_destroy_c(handle) bind(C, name="zone_destroy_")
        integer(c_intptr_t), value :: handle
        class(zone_t), pointer :: z
        if (handle == 0_c_intptr_t) return
        call handle_to_zone(handle, z)
        deallocate (z)
        nullify (zone_pool(int(handle))%z)
    end subroutine zone_destroy_c

    !> Pre-existing names (zone.h), unchanged signatures: legacy stitching
    !> Fortran already treats zone**/subclass** purely as an opaque handle.
    subroutine get_right_boundary_of_zone(handle, rout) bind(C, name="get_right_boundary_of_zone_")
        integer(c_intptr_t), intent(in) :: handle
        real(c_double), intent(out) :: rout
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        rout = z%get_r2()
    end subroutine get_right_boundary_of_zone

    subroutine get_left_boundary_of_zone(handle, rout) bind(C, name="get_left_boundary_of_zone_")
        integer(c_intptr_t), intent(in) :: handle
        real(c_double), intent(out) :: rout
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        rout = z%get_r1()
    end subroutine get_left_boundary_of_zone

    function c_string_to_fortran(cstr) result(fstr)
        character(kind=c_char), intent(in) :: cstr(*)
        character(len=:), allocatable :: fstr
        integer :: i
        i = 0
        do
            if (cstr(i + 1) == c_null_char) exit
            i = i + 1
        end do
        allocate (character(len=i) :: fstr)
        do i = 1, len(fstr)
            fstr(i:i) = cstr(i)
        end do
    end function c_string_to_fortran

end module kilca_zone_m
