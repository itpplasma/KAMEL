!> FLR (finite Larmor radius) plasma-response zone, formerly the C++
!> flre_zone class (flre_zone.{h,cpp}). The most complex zone_t subtype:
!> orchestrates the (already-Fortran) per-zone Maxwell-equations/
!> conductivity/sysmat/dispersion/quants sibling modules, drives the ODE
!> integration of the basis solutions (calculate_field_profiles_orth), and
!> re-spaces/transforms the final solution to the lab frame
!> (calc_final_fields).
!>
!> Field renamed from the C++ class: `dp` -> `dp_handle` (collides with the
!> `dp` real-kind parameter from `constants`).
!>
!> Several free C functions declared in flre_zone.h are callbacks invoked BY
!> pre-existing legacy Fortran (flre_sett_m.f90's setup_flre_data_module,
!> KiLCA/flre/maxwell_eqs/flre.f90, KiLCA/interface/wave_code_interface.cpp)
!> with a zone handle passed BY REFERENCE (address-of-handle, matching the
!> oracle's `flre_zone **ptr` convention) -- NOT by VALUE like the
!> mode.cpp-facing zone_*_c dispatch shims in kilca_zone_m. These are defined
!> here as bind(C) subroutines taking `integer(c_intptr_t), intent(in)`
!> (no VALUE attribute).
!>
!> path2linear/path2dispersion: the oracle recomputes these via
!> eval_path_to_linear_data/eval_path_to_dispersion_data (mode.cpp, still
!> C++, NOT extern "C" so unreachable from Fortran). mode.cpp already
!> computes the identical strings once at mode construction and copies them
!> into the shared `mode_data` Fortran module (mode_m.f90) via
!> copy_mode_paths_to_mode_data_module_; this module reads them from there
!> instead of recomputing, which is both simpler and guarantees identical
!> strings to whatever mode.cpp itself uses for the same mode.
module kilca_flre_zone_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_funptr, c_funloc, c_loc, c_f_pointer, c_null_ptr, c_null_char
    use constants, only: dp, c, im
    use kilca_zone_m, only: zone_t, zone_register, handle_to_zone, zone_read, zone_print, &
        skip_line, read_real_before_hash, read_int_before_hash, read_token_before_hash, &
        read_complex_before_hash, &
        BOUNDARY_CENTER, BOUNDARY_INFINITY, BOUNDARY_IDEALWALL, BOUNDARY_INTERFACE, BOUNDARY_ANTENNA
    use kilca_maxwell_eqs_data_m, only: maxwell_eqs_data_create, maxwell_eqs_data_destroy, &
        get_me_num_vars, get_me_num_eqs, get_me_der_order, &
        get_me_dim_ersp_state, get_me_iersp_state, &
        get_me_dim_ersp_sys, get_me_iersp_sys, &
        get_me_dim_brsp_sys, get_me_ibrsp_sys, &
        get_me_sys_ind, get_me_nwaves
    use kilca_cond_profiles_m, only: cond_profiles_create, cond_profiles_destroy, &
        get_cond_flag_back, get_cond_path2linear, get_cond_nc
    use kilca_sysmat_profiles_m, only: sysmat_profiles_create, sysmat_profiles_destroy, &
        get_sysmat_dimx, get_sysmat_x_ptr
    use kilca_disp_profiles_m, only: disp_profiles_create, disp_profiles_destroy, &
        disp_profiles_calculate, disp_profiles_save
    use kilca_flre_quants_m, only: flre_quants_create, flre_quants_destroy, &
        flre_quants_calculate_jae, flre_quants_calculate_local_profiles, &
        flre_quants_calculate_integrated_profiles, flre_quants_transform_quants_to_lab_cyl_frame, &
        flre_quants_save_profiles, flre_quants_interp_diss_power_density, &
        flre_quants_interp_current_density
    use kilca_background_data_m, only: get_background_x0, get_background_xlast, eval_hthz
    use kilca_transforms_m, only: transform_EB_from_rsp_to_cyl
    use kilca_solver_m, only: rhs_func
    use kilca_interp_m, only: sparse_grid_polynom
    use mode_data, only: md_path2linear => path2linear, md_path2dispersion => path2dispersion
    implicit none
    private

    public :: flre_zone_t
    public :: flre_zone_create_
    public :: set_equations_settings_c_, set_conductivity_settings_c_, set_collisions_settings_c_
    public :: get_iersp_sys_array_, get_ibrsp_sys_array_, get_sys_ind_array_
    public :: activate_fortran_modules_for_zone_, deactivate_fortran_modules_for_zone_
    public :: calc_flre_basis_in_lab_cyl_frame_with_full_system_vectors_
    public :: flre_zone_get_flre_order_, flre_zone_get_cp_

    type, bind(C) :: solver_settings_local_t
        integer(c_int) :: Nort
        real(c_double) :: eps_rel
        real(c_double) :: eps_abs
        real(c_double) :: norm_fac
        integer(c_int) :: debug
    end type solver_settings_local_t

    !> Mirrors kilca_solver_m's private rhs_func_params_t field-for-field
    !> (cannot `use` it directly: not in that module's public list).
    type, bind(C) :: rhs_func_params_local_t
        integer(c_int) :: Nwaves
        integer(c_int) :: Nphys
        integer(c_int) :: Nfs
        type(c_ptr) :: Dmat
        integer(c_intptr_t) :: sp
    end type rhs_func_params_local_t

    type, extends(zone_t) :: flre_zone_t
        integer(c_intptr_t) :: me = 0, cp = 0, sp = 0, dp_handle = 0, qp = 0
        integer(c_intptr_t) :: self_handle = 0
        integer :: flre_order = 0, Nmax = 0, gal_corr = 0, N = 0, max_dim_c = 0
        real(dp) :: D = 0, eps_out = 0, eps_res = 0
        real(dp) :: dr_out = 0, dr_res = 0, del = 0
        integer :: hom_sys = 0
        integer :: Nort = 0
        real(dp) :: norm_fac = 0
        integer :: Nfs = 0
        complex(dp), allocatable :: EB_mov(:, :)
        integer :: collmod(2) = 0
        integer :: rsp = 0
    contains
        procedure :: read_settings => flre_read_settings
        procedure :: print_settings => flre_print_settings
        procedure :: calc_basis_fields => flre_calc_basis_fields
        procedure :: copy_E_and_B_fields => flre_copy_E_and_B_fields
        procedure :: calc_final_fields => flre_calc_final_fields
        procedure :: calc_dispersion => flre_calc_dispersion
        procedure :: save_dispersion => flre_save_dispersion
        procedure :: calc_all_quants => flre_calc_all_quants
        procedure :: save_all_quants => flre_save_all_quants
        procedure :: eval_diss_power_density => flre_eval_diss_power_density
        procedure :: eval_current_density => flre_eval_current_density
        final :: flre_zone_finalize
    end type flre_zone_t

    interface calc_deriv_product
        module procedure calc_deriv_product_rr
        module procedure calc_deriv_product_rc
    end interface calc_deriv_product

    !> Free-function / legacy-Fortran callbacks. Names without a trailing
    !> underscore are pre-existing plain Fortran subprograms (90k-LOC legacy
    !> physics code, predating this port); gfortran's default external-name
    !> mangling already matches the trailing-underscore C declarations in
    !> flre_zone.h, so they are called here by their bare Fortran name.
    !> Array dummies use assumed-size (F77 sequence association) to match
    !> the legacy routines' own explicit-shape, descriptor-free ABI.
    interface
        subroutine setup_flre_data_module(zone)
            import :: c_intptr_t
            integer(c_intptr_t), intent(in) :: zone
        end subroutine setup_flre_data_module

        subroutine clean_flre_data_module()
        end subroutine clean_flre_data_module

        subroutine calc_and_set_maxwell_system_parameters_module()
        end subroutine calc_and_set_maxwell_system_parameters_module

        subroutine clean_maxwell_system_parameters_module()
        end subroutine clean_maxwell_system_parameters_module

        subroutine allocate_and_set_conductivity_arrays()
        end subroutine allocate_and_set_conductivity_arrays

        subroutine deallocate_conductivity_arrays()
        end subroutine deallocate_conductivity_arrays

        subroutine set_cond_profiles_in_mode_data_module(cp)
            import :: c_intptr_t
            integer(c_intptr_t), intent(in) :: cp
        end subroutine set_cond_profiles_in_mode_data_module

        subroutine set_sysmat_profiles_in_mode_data_module(sp)
            import :: c_intptr_t
            integer(c_intptr_t), intent(in) :: sp
        end subroutine set_sysmat_profiles_in_mode_data_module

        subroutine calc_start_values_center_with_correct_asymptotic(rstart, zstart)
            import :: dp
            real(dp), intent(in) :: rstart
            complex(dp), intent(out) :: zstart(*)
        end subroutine calc_start_values_center_with_correct_asymptotic

        subroutine calc_start_values_anywhere_low_derivs(rstart, zstart)
            import :: dp
            real(dp), intent(in) :: rstart
            complex(dp), intent(out) :: zstart(*)
        end subroutine calc_start_values_anywhere_low_derivs

        subroutine state2sys(rpt, flagback, v_state, v_sys, rhs)
            import :: dp
            real(dp), intent(in) :: rpt
            character(len=*), intent(in) :: flagback
            complex(dp), intent(in) :: v_state(*)
            complex(dp), intent(out) :: v_sys(*)
            complex(dp), intent(in) :: rhs(*)
        end subroutine state2sys

        subroutine normalize_flre_basis(D, Nw, dim, basis, iErsp_sys, ind1, ind2, node)
            import :: dp
            integer, intent(in) :: D, Nw, dim, ind1, ind2, node
            complex(dp), intent(inout) :: basis(*)
            integer, intent(in) :: iErsp_sys(3)
        end subroutine normalize_flre_basis

        function get_background_rtor() bind(C, name="get_background_rtor_") result(rtor)
            import :: c_double
            real(c_double) :: rtor
        end function get_background_rtor

        function get_background_V_gal_sys() bind(C, name="get_background_V_gal_sys_") result(vgal)
            import :: c_double
            real(c_double) :: vgal
        end function get_background_V_gal_sys

        function get_background_flag_back() bind(C, name="get_background_flag_back_") result(ch)
            import :: c_char
            character(kind=c_char) :: ch
        end function get_background_flag_back

        integer(c_int) function get_output_flag_dispersion() &
            bind(C, name="get_output_flag_dispersion_")
            import :: c_int
        end function get_output_flag_dispersion

        !> Still-C++ shared.cpp free function (binomial_coefficients_,
        !> extern "C"); same interface duplicated per-module already by
        !> kilca_cond_profiles_m/kilca_flre_quants_m, since it is module-
        !> private there too.
        subroutine binomial_coefficients_local(N, BC) bind(C, name="binomial_coefficients_")
            import :: c_int, c_double
            integer(c_int), value :: N
            real(c_double), intent(out) :: BC(*)
        end subroutine binomial_coefficients_local

        !> Local redeclaration of kilca_solver_m's integrate_basis_vecs with
        !> a COMPLEX Smat dummy (bind(C) explicit-shape array dummies are
        !> passed as plain base addresses, so this is ABI-identical to the
        !> module's own REAL(c_double) declaration -- same established
        !> convention as LAPACK/legacy-Fortran complex-as-flat-real calls
        !> elsewhere in this port).
        function integrate_basis_vecs_local(f, Nfs, Nw, dim, rvec, Smat, ss_ptr, params) &
            result(ret) bind(C, name="integrate_basis_vecs")
            import :: c_funptr, c_int, c_double, c_ptr, dp
            type(c_funptr), value :: f
            integer(c_int), value :: Nfs, Nw, dim
            real(c_double), intent(in) :: rvec(0:dim - 1)
            complex(dp), intent(inout) :: Smat(0:Nfs*Nw*dim - 1)
            type(c_ptr), value :: ss_ptr, params
            integer(c_int) :: ret
        end function integrate_basis_vecs_local
    end interface

contains

    !> ---- construction / destruction ----

    function flre_zone_create_(sd_ptr, bp_ptr, wd_handle, path, index_p) &
        result(handle) bind(C, name="flre_zone_create_")
        integer(c_intptr_t), value :: sd_ptr, bp_ptr, wd_handle
        character(kind=c_char), intent(in) :: path(*)
        integer(c_int), value :: index_p
        integer(c_intptr_t) :: handle

        type(flre_zone_t), pointer :: fz
        class(zone_t), pointer :: zp
        type(c_ptr) :: wd_cptr
        integer :: i

        allocate (fz)
        fz%bp = bp_ptr
        fz%index = int(index_p)
        fz%path = ''
        i = 0
        do
            if (path(i + 1) == c_null_char .or. i >= len(fz%path)) exit
            fz%path(i + 1:i + 1) = path(i + 1)
            i = i + 1
        end do
        wd_cptr = transfer(wd_handle, wd_cptr)
        call c_f_pointer(wd_cptr, fz%wd)

        fz%me = 0
        fz%cp = 0
        fz%sp = 0
        fz%dp_handle = 0
        fz%qp = 0

        zp => fz
        handle = zone_register(zp)
        fz%self_handle = handle
    end function flre_zone_create_

    subroutine flre_zone_finalize(self)
        type(flre_zone_t), intent(inout) :: self
        if (self%me /= 0) call maxwell_eqs_data_destroy(self%me)
        if (self%cp /= 0) call cond_profiles_destroy(self%cp)
        if (self%sp /= 0) call sysmat_profiles_destroy(self%sp)
        if (self%dp_handle /= 0) call disp_profiles_destroy(self%dp_handle)
        if (self%qp /= 0) call flre_quants_destroy(self%qp)
    end subroutine flre_zone_finalize

    !> ---- settings ----

    subroutine flre_read_settings(self, file)
        class(flre_zone_t), intent(inout) :: self
        character(len=*), intent(in) :: file
        integer :: unit, ios, k

        call zone_read(self, file)

        open (newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'error: flre_zone: read_settings: failed to open file ', trim(file)
            stop 1
        end if

        do k = 1, 8
            call skip_line(unit)
        end do

        call skip_line(unit)
        call read_int_before_hash(unit, self%flre_order)
        call read_int_before_hash(unit, self%Nmax)
        call read_int_before_hash(unit, self%gal_corr)
        call read_int_before_hash(unit, self%N)
        call read_int_before_hash(unit, self%max_dim_c)
        call read_real_before_hash(unit, self%D)
        call read_real_before_hash(unit, self%eps_out)
        call read_real_before_hash(unit, self%eps_res)
        call read_int_before_hash(unit, self%hom_sys)
        call skip_line(unit)

        call skip_line(unit)
        call read_int_before_hash(unit, self%max_dim)
        call read_real_before_hash(unit, self%eps_rel)
        call read_real_before_hash(unit, self%eps_abs)
        call read_int_before_hash(unit, self%Nort)
        call read_real_before_hash(unit, self%norm_fac)
        call read_real_before_hash(unit, self%dr_out)
        call read_real_before_hash(unit, self%dr_res)
        call read_real_before_hash(unit, self%del)
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

        call skip_line(unit)
        call read_int_before_hash(unit, self%collmod(1))
        call read_int_before_hash(unit, self%collmod(2))
        call skip_line(unit)

        close (unit)

        if (self%hom_sys /= 0) then
            write (*, '(a)') 'warning: homogenious limit should not be normally used.'
            stop 1
        end if

        self%Nwaves = 6*self%flre_order - 2
        self%Ncomps = 6*self%flre_order + 7

        if (self%bc1 == BOUNDARY_CENTER .or. self%bc1 == BOUNDARY_IDEALWALL .or. &
            self%bc2 == BOUNDARY_INFINITY .or. self%bc2 == BOUNDARY_IDEALWALL) then
            self%Nfs = self%Nwaves/2
        else
            self%Nfs = self%Nwaves
        end if

        if (self%flag_debug /= 0) call self%print_settings()
    end subroutine flre_read_settings

    subroutine flre_print_settings(self)
        class(flre_zone_t), intent(in) :: self
        call zone_print(self)
        write (*, '(a)') ''
        write (*, '(a)') 'Check for flre specific settings below:'
        write (*, '(a,i0)') 'order of FLR expansion: ', self%flre_order
        write (*, '(a,i0)') 'highest cyclotron harmonic: ', self%Nmax
        write (*, '(a,i0)') 'flag if to use correction term in conductivity: ', self%gal_corr
        write (*, '(a,i0)') 'splines degree: ', self%N
        write (*, '(a,i0)') &
            'maximum dimension of the grid for conductivity martices: ', self%max_dim_c
        write (*, '(a,es16.8e3)') 'resonant layer width: ', self%D
        write (*, '(a,es16.8e3)') &
            'error parameter for adaptive radial grid outside the resonant layer: ', self%eps_out
        write (*, '(a,es16.8e3)') &
            'error parameter for the adaptive radial grid in the resonant layer: ', self%eps_res
        write (*, '(a,i0)') 'flag for homogenious system: ', self%hom_sys
        write (*, '(a,i0)') &
            'max number of orthonormalization steps (ONS) for the solver: ', self%Nort
        write (*, '(a,es16.8e3)') 'controlling factor for ONS by QR: ', self%norm_fac
        write (*, '(a,es16.8e3)') &
            'output grid step outside the resonance region for the ME solutions: ', self%dr_out
        write (*, '(a,es16.8e3)') &
            'output grid step inside the resonance region for the ME solutions: ', self%dr_res
        write (*, '(a,es16.8e3)') 'width of the resonance region: ', self%del
        write (*, '(a,i0,1x,i0)') 'collisions model flags: ', self%collmod(1), self%collmod(2)
    end subroutine flre_print_settings

    !> ---- basis-field calculation ----

    subroutine flre_calc_basis_fields(self, flag)
        class(flre_zone_t), intent(inout) :: self
        integer, intent(in) :: flag
        real(dp) :: a_, b_
        character(kind=c_char), target :: path2linear_buf(1025)
        character(kind=c_char), target :: flag_back_buf(2)
        character(kind=c_char), target :: path2linear_buf2(1025)

        if (flag /= 0) then
            self%rsp = 0
        else
            self%rsp = 1
        end if

        call setup_flre_data_module(self%self_handle)

        call calc_and_set_maxwell_system_parameters_module()

        self%me = maxwell_eqs_data_create(self%Nwaves)

        if (self%flag_debug /= 0) call print_maxwell_eqs_data_local(self%me)

        call allocate_and_set_conductivity_arrays()

        a_ = max(self%r1 - 1.0_dp, get_background_x0())
        b_ = min(self%r2 + 1.0_dp, get_background_xlast())

        call to_cstr(trim(md_path2linear), path2linear_buf)

        self%cp = cond_profiles_create(c_loc(path2linear_buf), self%flre_order, self%gal_corr, &
                                        self%N, self%max_dim_c, self%r1, self%r2, self%D, &
                                        self%eps_out, self%eps_res, a_, b_, self%wd%r_res, &
                                        real(self%wd%omov, dp), aimag(self%wd%omov), &
                                        self%flag_debug, flag)

        call set_cond_profiles_in_mode_data_module(self%cp)

        if (flag /= 0) then
            call deallocate_conductivity_arrays()
            call clean_maxwell_system_parameters_module()
            call clean_flre_data_module()
            return
        end if

        flag_back_buf(1) = get_cond_flag_back(self%cp)
        flag_back_buf(2) = c_null_char
        call get_cond_path2linear(self%cp, path2linear_buf2)

        self%sp = sysmat_profiles_create(self%Nwaves, c_loc(flag_back_buf), &
                                          c_loc(path2linear_buf2), get_cond_nc(self%cp), &
                                          self%max_dim_c, self%eps_out, self%flag_debug, &
                                          self%r1, self%r2, self%wd%r_res)

        call set_sysmat_profiles_in_mode_data_module(self%sp)

        if (get_output_flag_dispersion() > 1) then
            call self%calc_dispersion()
            call self%save_dispersion()
        end if

        call flre_calculate_field_profiles_orth(self)

        call deallocate_conductivity_arrays()
        call clean_maxwell_system_parameters_module()
        call clean_flre_data_module()
    end subroutine flre_calc_basis_fields

    !> Translates flre_zone::calculate_field_profiles_orth. Boundary-
    !> condition-dependent start values + integration direction, an
    !> adaptive grid thickened near wd%r_res via dr_out/dr_res/del, ODE
    !> integration of the fundamental solutions, and basis normalization.
    !> Error branches (BOUNDARY_INFINITY "not implemented", unknown BCs)
    !> stop exactly where the oracle exit(1)s, matching the oracle's own
    !> commented-out alternative start-value calls (left as comments, never
    !> activated).
    subroutine flre_calculate_field_profiles_orth(self)
        class(flre_zone_t), intent(inout) :: self
        complex(dp), allocatable :: y(:, :)
        integer :: ind1, ind2
        real(dp) :: ri, rf
        integer :: Neq
        real(dp), allocatable, target :: grid(:)
        real(dp) :: sgn, r_res, rcur, dr
        integer :: node, dim_local
        complex(dp), allocatable, target :: state(:, :, :)
        integer :: j, j1, j2, dj, k, ind, node_param
        real(dp), allocatable, target :: Dmat(:)
        type(solver_settings_local_t), target :: ss
        type(rhs_func_params_local_t), target :: params
        integer(c_int) :: ret
        integer :: iErsp_sys_local(0:2)

        allocate (y(self%Nwaves, self%Nwaves))
        y = (0.0_dp, 0.0_dp)

        if (self%bc1 == BOUNDARY_CENTER) then
            ri = self%r1
            rf = self%r2
            !calc_start_values_anywhere_low_derivs_ (&ri, y); -- oracle also has this commented out
            call calc_start_values_center_with_correct_asymptotic(ri, y(:, 1:self%Nwaves/2))
            self%Nfs = self%Nwaves/2
            ind1 = 0
            ind2 = self%Nfs - 1
        else if (self%bc1 == BOUNDARY_IDEALWALL) then
            ri = self%r1
            rf = self%r2
            call calc_start_values_anywhere_low_derivs(ri, y(:, 1:self%Nwaves/2))
            self%Nfs = self%Nwaves/2
            ind1 = 0
            ind2 = self%Nfs - 1
        else if (self%bc2 == BOUNDARY_INFINITY) then
            ri = self%r2
            rf = self%r1
            self%Nfs = self%Nwaves/2
            ind1 = self%Nfs
            ind2 = self%Nwaves - 1
            write (*, '(a)') 'flre_zone::calculate_field_profiles_orth: not implemented!'
            stop 1
        else if (self%bc2 == BOUNDARY_IDEALWALL) then
            ri = self%r2
            rf = self%r1
            call calc_start_values_anywhere_low_derivs(ri, y(:, 1:self%Nwaves/2))
            self%Nfs = self%Nwaves/2
            ind1 = 0
            ind2 = self%Nfs - 1
        else
            ri = self%r1
            rf = self%r2
            self%Nfs = self%Nwaves
            ind1 = 0
            ind2 = self%Nwaves - 1
            write (*, '(a)') 'flre_zone::calculate_field_profiles_orth: unknown BCs!'
            stop 1
        end if

        Neq = 2*self%Nwaves*self%Nfs

        sgn = real(signum_local(rf - ri), dp)

        if (self%r1 <= self%wd%r_res .and. self%r2 >= self%wd%r_res) then
            r_res = self%wd%r_res

            rcur = ri
            node = 0
            do while ((rcur - rf)*sgn < 0.0_dp)
                rcur = rcur + ((self%dr_out - self%dr_res)* &
                               (1.0_dp - exp(-(rcur - r_res)*(rcur - r_res)/self%del/self%del)) + &
                               self%dr_res)*sgn
                node = node + 1
            end do

            dim_local = node + 1
            allocate (grid(dim_local))

            rcur = ri
            node = 0
            do while ((rcur - rf)*sgn < 0.0_dp)
                node = node + 1
                grid(node) = rcur
                rcur = rcur + ((self%dr_out - self%dr_res)* &
                               (1.0_dp - exp(-(rcur - r_res)*(rcur - r_res)/self%del/self%del)) + &
                               self%dr_res)*sgn
            end do
            grid(node + 1) = rf
        else
            dim_local = int(abs(rf - ri)/self%dr_out)
            allocate (grid(dim_local))
            dr = (rf - ri)/real(dim_local - 1, dp)
            do node = 0, dim_local - 1
                grid(node + 1) = ri + dr*node
            end do
            grid(dim_local) = rf
        end if

        allocate (state(self%Nwaves, self%Nfs, dim_local))
        state = (0.0_dp, 0.0_dp)

        do j = 1, self%Nfs
            call galilean_transform_of_flre_state_vector(self, get_background_V_gal_sys(), &
                                                           self%wd%olab, ri, y(:, j), &
                                                           state(:, j, 1))
        end do

        allocate (Dmat(2*self%Nwaves*(self%Nwaves + 2*self%Nfs)))
        Dmat = 0.0_dp

        params%Nwaves = self%Nwaves
        params%Nphys = self%Nwaves
        params%Nfs = self%Nfs
        params%Dmat = c_loc(Dmat)
        params%sp = self%sp

        ss%Nort = self%Nort
        ss%eps_rel = self%eps_rel
        ss%eps_abs = self%eps_abs
        ss%norm_fac = self%norm_fac
        ss%debug = self%flag_debug

        ret = integrate_basis_vecs_local(c_funloc(rhs_func), self%Nfs, self%Nwaves, dim_local, &
                                          grid, state, c_loc(ss), c_loc(params))

        deallocate (Dmat)

        self%dim = dim_local
        if (allocated(self%r)) deallocate (self%r)
        allocate (self%r(self%dim))
        if (allocated(self%basis)) deallocate (self%basis)
        allocate (self%basis(self%Ncomps, self%Nwaves, self%dim))
        self%basis = (0.0_dp, 0.0_dp)

        if (grid(1) == self%r1 .and. grid(dim_local) == self%r2) then
            j1 = 1
            j2 = dim_local
            dj = 1
        else if (grid(1) == self%r2 .and. grid(dim_local) == self%r1) then
            j1 = dim_local
            j2 = 1
            dj = -1
        else
            write (*, '(a)') 'flre_zone::calculate_field_profiles_orth: error!'
            stop 1
        end if
        node_param = merge(dim_local - 1, 0, dj == 1)

        k = 0
        j = j1
        do while ((j - j2)*dj <= 0)
            k = k + 1
            self%r(k) = grid(j)
            do ind = ind1, ind2
                call state_to_system_copy(self, state(:, ind - ind1 + 1, j), &
                                           self%basis(:, ind + 1, k))
            end do
            j = j + dj
        end do

        deallocate (grid, state, y)

        iErsp_sys_local(0) = get_me_iersp_sys(self%me, 0)
        iErsp_sys_local(1) = get_me_iersp_sys(self%me, 1)
        iErsp_sys_local(2) = get_me_iersp_sys(self%me, 2)
        call normalize_flre_basis(self%Ncomps, self%Nwaves, self%dim, self%basis, &
                                   iErsp_sys_local, ind1, ind2, node_param)
    end subroutine flre_calculate_field_profiles_orth

    subroutine flre_copy_E_and_B_fields(self, EB_out)
        class(flre_zone_t), intent(in) :: self
        real(dp), intent(out) :: EB_out(*)
        complex(dp) :: EBrsp(6), EBcyl(6)
        integer :: node, comp, idx

        do node = 0, self%dim - 1
            do comp = 0, 5
                idx = iF_comp(self, comp)
                EBrsp(comp + 1) = self%EB(idx + 1, node + 1)
            end do
            call transform_EB_from_rsp_to_cyl(self%bp, self%r(node + 1), EBrsp, EBcyl)
            do comp = 0, 5
                EB_out(2*(comp + 6*node) + 1) = real(EBcyl(comp + 1), dp)
                EB_out(2*(comp + 6*node) + 2) = aimag(EBcyl(comp + 1))
            end do
        end do
    end subroutine flre_copy_E_and_B_fields

    !> Translates flre_zone::calc_final_fields: thins the grid via
    !> sparse_grid_polynom, evaluates full system vectors from the state
    !> vectors via the legacy state2sys, then Galilean-transforms back to
    !> the lab frame. self%r/self%EB are NOT reallocated (matching the
    !> oracle, which overwrites the existing buffers in place and only
    !> shrinks the logical `dim`); only self%EB_mov is freshly allocated at
    !> the new (smaller-or-equal) dimension.
    subroutine flre_calc_final_fields(self)
        class(flre_zone_t), intent(inout) :: self
        integer, allocatable :: ind(:)
        real(dp), allocatable :: rnew(:)
        real(dp), allocatable :: snew(:, :), sold(:, :)
        complex(dp) :: state_tmp(self%Nwaves)
        integer :: dimnew, num_vars, num_eqs, k, i, dim_old
        complex(dp), allocatable :: rhs(:), v_sys(:)
        character(kind=c_char) :: flag_back
        real(dp) :: aeps_local, reps_local

        dim_old = self%dim

        allocate (ind(dim_old), rnew(dim_old))
        allocate (sold(2*self%Nwaves, dim_old), snew(2*self%Nwaves, dim_old))

        do k = 1, dim_old
            call system_to_state_copy(self, self%EB(:, k), state_tmp)
            do i = 1, self%Nwaves
                sold(2*i - 1, k) = real(state_tmp(i), dp)
                sold(2*i, k) = aimag(state_tmp(i))
            end do
        end do

        aeps_local = self%aeps
        reps_local = self%reps
        call sparse_grid_polynom(self%r, dim_old, sold, 2*self%Nwaves, self%deg, &
                                  aeps_local, reps_local, self%step, dimnew, rnew, snew, ind)

        num_vars = get_me_num_vars(self%me)
        num_eqs = get_me_num_eqs(self%me)

        allocate (rhs(num_eqs))
        rhs = (0.0_dp, 0.0_dp)

        call flre_activate_fortran_modules(self)

        self%dim = dimnew
        if (allocated(self%EB_mov)) deallocate (self%EB_mov)
        allocate (self%EB_mov(self%Ncomps, self%dim))

        flag_back = get_background_flag_back()

        allocate (v_sys(num_vars))
        do k = 1, self%dim
            self%r(k) = rnew(k)
            do i = 1, self%Nwaves
                state_tmp(i) = cmplx(snew(2*i - 1, k), snew(2*i, k), dp)
            end do
            call state2sys(self%r(k), flag_back, state_tmp, v_sys, rhs)
            self%EB_mov(1:self%Ncomps, k) = v_sys(1:self%Ncomps)
        end do

        call flre_deactivate_fortran_modules(self)

        do k = 1, self%dim
            call galilean_transform_of_flre_system_vector(self, -get_background_V_gal_sys(), &
                                                            self%wd%omov, self%r(k), &
                                                            self%EB_mov(:, k), self%EB(:, k))
        end do

        !> if (DEBUG_FLAG) save_system_vector_in_mov_frame(); -- DEBUG_FLAG
        !> is hard-coded 0 in code_settings.h (matching the established
        !> flre_quants_m.f90 precedent for this exact macro), so this call
        !> is dead code in the oracle and is not invoked here either.

        deallocate (ind, rnew, sold, snew, rhs, v_sys)
    end subroutine flre_calc_final_fields

    !> Dead in the oracle (DEBUG_FLAG hard-coded 0 in code_settings.h), kept
    !> as a translated-but-uncalled private procedure for completeness, per
    !> the task's translation list.
    subroutine flre_save_system_vector_in_mov_frame(self)
        class(flre_zone_t), intent(in) :: self
        character(len=1024) :: fname
        integer :: unit, i, comp

        write (fname, '(a,a,i0,a)') trim(md_path2linear), 'zone_', self%index, '_EB_mov.dat'
        open (newunit=unit, file=trim(fname), status='replace', action='write')
        do i = 1, self%dim
            write (unit, '(es24.16e3)', advance='no') self%r(i)
            do comp = 1, self%Ncomps
                write (unit, '(a,es24.16e3,a,es24.16e3)', advance='no') &
                    char(9), real(self%EB_mov(comp, i), dp), char(9), aimag(self%EB_mov(comp, i))
            end do
            write (unit, *)
        end do
        close (unit)
    end subroutine flre_save_system_vector_in_mov_frame

    subroutine flre_calc_dispersion(self)
        class(flre_zone_t), intent(inout) :: self
        character(kind=c_char), target :: flag_back_buf(2)
        flag_back_buf(1) = get_background_flag_back()
        flag_back_buf(2) = c_null_char
        self%dp_handle = disp_profiles_create(self%Nwaves, get_sysmat_dimx(self%sp), &
                                               get_sysmat_x_ptr(self%sp), c_loc(flag_back_buf))
        call disp_profiles_calculate(self%dp_handle)
    end subroutine flre_calc_dispersion

    subroutine flre_save_dispersion(self)
        class(flre_zone_t), intent(inout) :: self
        character(kind=c_char), target :: filename(1025)
        character(len=1024) :: fname_str

        write (fname_str, '(a,a,i0,a)') trim(md_path2dispersion), 'zone_', self%index, '_kr.dat'
        call to_cstr(trim(fname_str), filename)
        call disp_profiles_save(self%dp_handle, c_loc(filename))
    end subroutine flre_save_dispersion

    subroutine flre_calc_all_quants(self)
        class(flre_zone_t), intent(inout) :: self
        class(zone_t), pointer :: zp
        type(c_ptr) :: x_cptr, ebmov_cptr
        character(kind=c_char), target :: path2linear_buf(1025)

        call handle_to_zone(self%self_handle, zp)
        select type (zp)
        type is (flre_zone_t)
            x_cptr = c_loc(zp%r)
            ebmov_cptr = c_loc(zp%EB_mov)
        end select

        call to_cstr(trim(md_path2linear), path2linear_buf)

        self%qp = flre_quants_create(self%cp, self%me, c_null_ptr, c_loc(path2linear_buf), &
                                      self%flre_order, self%dim, x_cptr, self%Ncomps, ebmov_cptr, &
                                      real(self%wd%omov, dp), aimag(self%wd%omov), &
                                      self%bc1, self%bc2, self%index)

        call flre_quants_calculate_jae(self%qp)
        call flre_quants_calculate_local_profiles(self%qp)
        call flre_quants_calculate_integrated_profiles(self%qp)
        call flre_quants_transform_quants_to_lab_cyl_frame(self%qp)
    end subroutine flre_calc_all_quants

    subroutine flre_save_all_quants(self)
        class(flre_zone_t), intent(inout) :: self
        call flre_quants_save_profiles(self%qp)
    end subroutine flre_save_all_quants

    subroutine flre_eval_diss_power_density(self, x, ttype, spec, dpd)
        class(flre_zone_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(in) :: ttype, spec
        real(dp), intent(out) :: dpd(*)
        call flre_quants_interp_diss_power_density(self%qp, x, ttype, spec, dpd)
    end subroutine flre_eval_diss_power_density

    subroutine flre_eval_current_density(self, x, ttype, spec, comp, J)
        class(flre_zone_t), intent(in) :: self
        real(dp), intent(in) :: x
        integer, intent(in) :: ttype, spec, comp
        real(dp), intent(out) :: J(*)
        call flre_quants_interp_current_density(self%qp, x, ttype, spec, comp, J)
    end subroutine flre_eval_current_density

    !> ---- activation helpers (shared by calc_final_fields and the
    !> bind(C) activate_/deactivate_fortran_modules_for_zone_ entry points
    !> called from wave_code_interface.cpp) ----

    subroutine flre_activate_fortran_modules(self)
        class(flre_zone_t), intent(inout) :: self
        call setup_flre_data_module(self%self_handle)
        call calc_and_set_maxwell_system_parameters_module()
        call allocate_and_set_conductivity_arrays()
        call set_cond_profiles_in_mode_data_module(self%cp)
        call set_sysmat_profiles_in_mode_data_module(self%sp)
    end subroutine flre_activate_fortran_modules

    subroutine flre_deactivate_fortran_modules(self)
        class(flre_zone_t), intent(inout) :: self
        call deallocate_conductivity_arrays()
        call clean_maxwell_system_parameters_module()
        call clean_flre_data_module()
    end subroutine flre_deactivate_fortran_modules

    !> ---- inline-helper translations (system_to_state_copy etc) ----

    !> Indexing function for Er,Es,Ep,Br,Bs,Bp (comp=0:5) -> the physical
    !> system-vector component index (0-based) that get_me_iersp_sys_/
    !> get_me_ibrsp_sys_ map it to. NOT an array-layout helper (basis/EB are
    !> native Fortran arrays, indexed directly elsewhere): this is purely
    !> the maxwell_eqs_data component lookup the oracle's flre_zone::iF(int)
    !> performed.
    integer function iF_comp(self, comp) result(idx)
        class(flre_zone_t), intent(in) :: self
        integer, intent(in) :: comp
        if (comp < 3) then
            idx = get_me_iersp_sys(self%me, comp)
        else
            idx = get_me_ibrsp_sys(self%me, comp - 3)
        end if
    end function iF_comp

    subroutine system_to_state_copy(self, system, state)
        class(flre_zone_t), intent(in) :: self
        complex(dp), intent(in) :: system(*)
        complex(dp), intent(out) :: state(*)
        integer :: k, i, dim_k, i_state, i_sys
        do k = 0, 2
            dim_k = get_me_dim_ersp_state(self%me, k)
            i_state = get_me_iersp_state(self%me, k)
            i_sys = get_me_iersp_sys(self%me, k)
            do i = 0, dim_k - 1
                state(i_state + i + 1) = system(i_sys + i + 1)
            end do
        end do
    end subroutine system_to_state_copy

    !> Only the Ersp slots of `system` are written (matching the oracle:
    !> "warning: basis array is filled (partly) by state vectors").
    subroutine state_to_system_copy(self, state, system)
        class(flre_zone_t), intent(in) :: self
        complex(dp), intent(in) :: state(*)
        complex(dp), intent(inout) :: system(*)
        integer :: k, i, dim_k, i_state, i_sys
        do k = 0, 2
            dim_k = get_me_dim_ersp_state(self%me, k)
            i_state = get_me_iersp_state(self%me, k)
            i_sys = get_me_iersp_sys(self%me, k)
            do i = 0, dim_k - 1
                system(i_sys + i + 1) = state(i_state + i + 1)
            end do
        end do
    end subroutine state_to_system_copy

    !> Galilean transform of a state vector (Er,Es,Ep with their radial
    !> derivative ladders). Assumes dim_Ersp_state(0)=dim_Ersp_state(1)=
    !> dim_Ersp_state(2) (= ho+1), matching the oracle's own implicit
    !> assumption (it sizes ht/hz/mor/ks/kp at ho+1 = dim_Ersp_state(2) and
    !> reads E1c[1]/E1c[2] up to index ho regardless of their own per-
    !> component allocation size).
    subroutine galilean_transform_of_flre_state_vector(self, V, omega, rpt, E1, E2)
        class(flre_zone_t), intent(in) :: self
        real(dp), intent(in) :: V, rpt
        complex(dp), intent(in) :: omega
        complex(dp), intent(in) :: E1(0:*)
        complex(dp), intent(out) :: E2(0:*)

        integer :: dim_Ersp_state(0:2), iErsp_state(0:2), mo(0:2)
        integer :: m, ho, o, k
        real(dp) :: kz
        real(dp), allocatable, target :: Cbin(:, :), htz(:), mor(:), ks(:), kp(:)
        real(dp), pointer :: ht(:), hz(:)
        real(dp) :: htmor, hzmor
        complex(dp), allocatable :: Br(:), Ez(:), htBr(:), hzBr(:)
        complex(dp) :: kpEs, ksEp, hzEp, htEs
        complex(dp), allocatable :: E1c(:, :), E2c(:, :)

        do k = 0, 2
            dim_Ersp_state(k) = get_me_dim_ersp_state(self%me, k)
            iErsp_state(k) = get_me_iersp_state(self%me, k)
        end do
        mo = dim_Ersp_state - 1

        m = self%wd%m
        kz = real(self%wd%n, dp)/get_background_rtor()

        ho = mo(2)

        allocate (Cbin(0:ho, 0:ho))
        call binomial_coefficients_local(ho, Cbin)

        allocate (htz(0:2*ho + 1))
        call eval_hthz(rpt, 0, ho, c_null_ptr, htz)
        ht => htz(0:ho)
        hz => htz(ho + 1:2*ho + 1)

        allocate (mor(0:ho))
        mor(0) = real(m, dp)/rpt
        do o = 1, ho
            mor(o) = (-real(o, dp)/rpt)*mor(o - 1)
        end do

        allocate (ks(0:ho), kp(0:ho))
        do o = 0, ho
            htmor = calc_deriv_product(ho, Cbin, o, ht, mor)
            hzmor = calc_deriv_product(ho, Cbin, o, hz, mor)
            ks(o) = hzmor - ht(o)*kz
            kp(o) = htmor + hz(o)*kz
        end do

        allocate (E1c(0:2, 0:ho), E2c(0:2, 0:ho))
        E1c = (0.0_dp, 0.0_dp)
        E2c = (0.0_dp, 0.0_dp)

        do k = 0, 2
            do o = 0, mo(k)
                E1c(k, o) = E1(iErsp_state(k) + o)
            end do
        end do

        allocate (Br(0:ho), Ez(0:ho), htBr(0:ho), hzBr(0:ho))

        do o = 0, ho
            kpEs = calc_deriv_product(ho, Cbin, o, kp, E1c(1, :))
            ksEp = calc_deriv_product(ho, Cbin, o, ks, E1c(2, :))

            Br(o) = (c/omega)*(ksEp - kpEs)

            htBr(o) = calc_deriv_product(ho, Cbin, o, ht, Br)
            hzBr(o) = calc_deriv_product(ho, Cbin, o, hz, Br)

            htEs = calc_deriv_product(ho, Cbin, o, ht, E1c(1, :))
            hzEp = calc_deriv_product(ho, Cbin, o, hz, E1c(2, :))

            Ez(o) = hzEp - htEs
        end do

        do o = 0, mo(0)
            E2c(0, o) = E1c(0, o) - V/omega*(kz*E1c(0, o) + im*Ez(o + 1))
        end do
        do o = 0, mo(1)
            E2c(1, o) = E1c(1, o) + (V/c)*hzBr(o)
        end do
        do o = 0, mo(2)
            E2c(2, o) = E1c(2, o) + (V/c)*htBr(o)
        end do

        do k = 0, 2
            do o = 0, mo(k)
                E2(iErsp_state(k) + o) = E2c(k, o)
            end do
        end do
    end subroutine galilean_transform_of_flre_state_vector

    !> Galilean transform of a full (E,B) system vector. Same uniform-ho
    !> sizing assumption as galilean_transform_of_flre_state_vector, also
    !> covering dim_Brsp_sys <= dim_Ersp_sys(2). B2c is set equal to B1c
    !> (the oracle's V/c correction terms for B are commented out in the
    !> source and never activated); Et/htEr/hzEr are computed but unused,
    !> matching the oracle's own dead-but-computed values exactly.
    subroutine galilean_transform_of_flre_system_vector(self, V, omega, rpt, EB1, EB2)
        class(flre_zone_t), intent(in) :: self
        real(dp), intent(in) :: V, rpt
        complex(dp), intent(in) :: omega
        complex(dp), intent(in) :: EB1(0:*)
        complex(dp), intent(out) :: EB2(0:*)

        integer :: dim_Ersp_sys(0:2), iErsp_sys(0:2)
        integer :: dim_Brsp_sys(0:2), iBrsp_sys(0:2)
        integer :: moe(0:2), mob(0:2)
        integer :: m, ho, o, k
        real(dp) :: kz
        real(dp), allocatable, target :: Cbin(:, :), htz(:), mor(:), ks(:), kp(:)
        real(dp), pointer :: ht(:), hz(:)
        real(dp) :: htmor, hzmor
        complex(dp), allocatable :: E1c(:, :), E2c(:, :), B1c(:, :), B2c(:, :)
        complex(dp), allocatable :: Br(:), Ez(:), Et(:), htBr(:), hzBr(:), htEr(:), hzEr(:)
        complex(dp) :: kpEs, ksEp, hzEp, htEs, hzEs, htEp

        do k = 0, 2
            dim_Ersp_sys(k) = get_me_dim_ersp_sys(self%me, k)
            iErsp_sys(k) = get_me_iersp_sys(self%me, k)
            dim_Brsp_sys(k) = get_me_dim_brsp_sys(self%me, k)
            iBrsp_sys(k) = get_me_ibrsp_sys(self%me, k)
        end do
        moe = dim_Ersp_sys - 1
        mob = dim_Brsp_sys - 1

        m = self%wd%m
        kz = real(self%wd%n, dp)/get_background_rtor()

        ho = moe(2)

        allocate (Cbin(0:ho, 0:ho))
        call binomial_coefficients_local(ho, Cbin)

        allocate (htz(0:2*ho + 1))
        call eval_hthz(rpt, 0, ho, c_null_ptr, htz)
        ht => htz(0:ho)
        hz => htz(ho + 1:2*ho + 1)

        allocate (mor(0:ho))
        mor(0) = real(m, dp)/rpt
        do o = 1, ho
            mor(o) = (-real(o, dp)/rpt)*mor(o - 1)
        end do

        allocate (ks(0:ho), kp(0:ho))
        do o = 0, ho
            htmor = calc_deriv_product(ho, Cbin, o, ht, mor)
            hzmor = calc_deriv_product(ho, Cbin, o, hz, mor)
            ks(o) = hzmor - ht(o)*kz
            kp(o) = htmor + hz(o)*kz
        end do

        allocate (E1c(0:2, 0:ho), E2c(0:2, 0:ho), B1c(0:2, 0:ho), B2c(0:2, 0:ho))
        E1c = (0.0_dp, 0.0_dp)
        E2c = (0.0_dp, 0.0_dp)
        B1c = (0.0_dp, 0.0_dp)
        B2c = (0.0_dp, 0.0_dp)

        do k = 0, 2
            do o = 0, moe(k)
                E1c(k, o) = EB1(iErsp_sys(k) + o)
            end do
            do o = 0, mob(k)
                B1c(k, o) = EB1(iBrsp_sys(k) + o)
            end do
        end do

        allocate (Br(0:ho), Ez(0:ho), Et(0:ho), htBr(0:ho), hzBr(0:ho), htEr(0:ho), hzEr(0:ho))

        do o = 0, ho
            kpEs = calc_deriv_product(ho, Cbin, o, kp, E1c(1, :))
            ksEp = calc_deriv_product(ho, Cbin, o, ks, E1c(2, :))

            Br(o) = (c/omega)*(ksEp - kpEs)

            htBr(o) = calc_deriv_product(ho, Cbin, o, ht, Br)
            hzBr(o) = calc_deriv_product(ho, Cbin, o, hz, Br)

            htEs = calc_deriv_product(ho, Cbin, o, ht, E1c(1, :))
            hzEp = calc_deriv_product(ho, Cbin, o, hz, E1c(2, :))

            htEp = calc_deriv_product(ho, Cbin, o, ht, E1c(2, :))
            hzEs = calc_deriv_product(ho, Cbin, o, hz, E1c(1, :))

            Et(o) = hzEs + htEp
            Ez(o) = hzEp - htEs
        end do

        do o = 0, mob(2)
            htEr(o) = calc_deriv_product(ho, Cbin, o, ht, E1c(0, :))
            hzEr(o) = calc_deriv_product(ho, Cbin, o, hz, E1c(0, :))
        end do

        do o = 0, moe(0)
            E2c(0, o) = E1c(0, o) - V/omega*(kz*E1c(0, o) + im*Ez(o + 1))
        end do
        do o = 0, moe(1)
            E2c(1, o) = E1c(1, o) + (V/c)*hzBr(o)
        end do
        do o = 0, moe(2)
            E2c(2, o) = E1c(2, o) + (V/c)*htBr(o)
        end do

        do o = 0, mob(0)
            B2c(0, o) = B1c(0, o)
        end do
        do o = 0, mob(1)
            B2c(1, o) = B1c(1, o)
        end do
        do o = 0, mob(2)
            B2c(2, o) = B1c(2, o)
        end do

        do k = 0, 2
            do o = 0, moe(k)
                EB2(iErsp_sys(k) + o) = E2c(k, o)
            end do
            do o = 0, mob(k)
                EB2(iBrsp_sys(k) + o) = B2c(k, o)
            end do
        end do
    end subroutine galilean_transform_of_flre_system_vector

    !> Transform of a full system vector from (Er,Es,Ep,Br,Bs,Bp) to
    !> (Er,Et,Ez,Br,Bt,Bz). Same uniform-ho sizing assumption as the two
    !> galilean transforms above.
    subroutine transform_of_flre_system_vector_to_cyl_coordinates(self, rpt, EB1, EB2)
        class(flre_zone_t), intent(in) :: self
        real(dp), intent(in) :: rpt
        complex(dp), intent(in) :: EB1(0:*)
        complex(dp), intent(out) :: EB2(0:*)

        integer :: dim_Ersp_sys(0:2), iErsp_sys(0:2)
        integer :: dim_Brsp_sys(0:2), iBrsp_sys(0:2)
        integer :: moe(0:2), mob(0:2)
        integer :: ho, o, k
        real(dp), allocatable, target :: Cbin(:, :), htz(:)
        real(dp), pointer :: ht(:), hz(:)
        complex(dp), allocatable :: E1c(:, :), E2c(:, :), B1c(:, :), B2c(:, :)
        complex(dp) :: hzAp, htAs, hzAs, htAp

        do k = 0, 2
            dim_Ersp_sys(k) = get_me_dim_ersp_sys(self%me, k)
            iErsp_sys(k) = get_me_iersp_sys(self%me, k)
            dim_Brsp_sys(k) = get_me_dim_brsp_sys(self%me, k)
            iBrsp_sys(k) = get_me_ibrsp_sys(self%me, k)
        end do
        moe = dim_Ersp_sys - 1
        mob = dim_Brsp_sys - 1

        ho = moe(2)

        allocate (Cbin(0:ho, 0:ho))
        call binomial_coefficients_local(ho, Cbin)

        allocate (htz(0:2*ho + 1))
        call eval_hthz(rpt, 0, ho, c_null_ptr, htz)
        ht => htz(0:ho)
        hz => htz(ho + 1:2*ho + 1)

        allocate (E1c(0:2, 0:ho), E2c(0:2, 0:ho), B1c(0:2, 0:ho), B2c(0:2, 0:ho))
        E1c = (0.0_dp, 0.0_dp)
        E2c = (0.0_dp, 0.0_dp)
        B1c = (0.0_dp, 0.0_dp)
        B2c = (0.0_dp, 0.0_dp)

        do k = 0, 2
            do o = 0, moe(k)
                E1c(k, o) = EB1(iErsp_sys(k) + o)
            end do
            do o = 0, mob(k)
                B1c(k, o) = EB1(iBrsp_sys(k) + o)
            end do
        end do

        do o = 0, moe(0)
            E2c(0, o) = E1c(0, o)
        end do
        do o = 0, moe(1)
            hzAs = calc_deriv_product(ho, Cbin, o, hz, E1c(1, :))
            htAp = calc_deriv_product(ho, Cbin, o, ht, E1c(2, :))
            E2c(1, o) = hzAs + htAp
        end do
        do o = 0, moe(2)
            hzAp = calc_deriv_product(ho, Cbin, o, hz, E1c(2, :))
            htAs = calc_deriv_product(ho, Cbin, o, ht, E1c(1, :))
            E2c(2, o) = hzAp - htAs
        end do

        do o = 0, mob(0)
            B2c(0, o) = B1c(0, o)
        end do
        do o = 0, mob(1)
            hzAs = calc_deriv_product(ho, Cbin, o, hz, B1c(1, :))
            htAp = calc_deriv_product(ho, Cbin, o, ht, B1c(2, :))
            B2c(1, o) = hzAs + htAp
        end do
        do o = 0, mob(2)
            hzAp = calc_deriv_product(ho, Cbin, o, hz, B1c(2, :))
            htAs = calc_deriv_product(ho, Cbin, o, ht, B1c(1, :))
            B2c(2, o) = hzAp - htAs
        end do

        do k = 0, 2
            do o = 0, moe(k)
                EB2(iErsp_sys(k) + o) = E2c(k, o)
            end do
            do o = 0, mob(k)
                EB2(iBrsp_sys(k) + o) = B2c(k, o)
            end do
        end do
    end subroutine transform_of_flre_system_vector_to_cyl_coordinates

    !> calc_derive_of_func_product<T1,T2,T3>: fg[n] = sum_k C(n,k)*f(k)*g(n-k).
    !> The oracle's only live instantiations are (real,real->real) and
    !> (real,complex->complex).
    real(dp) function calc_deriv_product_rr(N, Cbin, ordn, f, g) result(fg)
        integer, intent(in) :: N, ordn
        real(dp), intent(in) :: Cbin(0:N, 0:N), f(0:N), g(0:N)
        integer :: k
        fg = 0.0_dp
        do k = 0, ordn
            fg = fg + Cbin(ordn, k)*f(k)*g(ordn - k)
        end do
    end function calc_deriv_product_rr

    function calc_deriv_product_rc(N, Cbin, ordn, f, g) result(fg)
        integer, intent(in) :: N, ordn
        real(dp), intent(in) :: Cbin(0:N, 0:N), f(0:N)
        complex(dp), intent(in) :: g(0:N)
        complex(dp) :: fg
        integer :: k
        fg = (0.0_dp, 0.0_dp)
        do k = 0, ordn
            fg = fg + Cbin(ordn, k)*f(k)*g(ordn - k)
        end do
    end function calc_deriv_product_rc

    integer function signum_local(x) result(s)
        real(dp), intent(in) :: x
        if (x < 0.0_dp) then
            s = -1
        else if (x == 0.0_dp) then
            s = 0
        else
            s = 1
        end if
    end function signum_local

    !> Mirrors maxwell_eqs_data.h's inline print_maxwell_eqs_data (header-
    !> only C++, not extern "C", so unreachable from Fortran): reads the
    !> same fields via this module's own get_me_* getters.
    subroutine print_maxwell_eqs_data_local(me)
        integer(c_intptr_t), intent(in) :: me
        integer :: k, Nw

        write (*, '(a)') ''
        write (*, '(a)') 'check Maxwell system parameters:'
        write (*, '(a,i0)') 'num_vars=', get_me_num_vars(me)
        write (*, '(a,i0)') 'num_eqs=', get_me_num_eqs(me)
        write (*, '(a,i0,1x,i0,1x,i0)') 'dim_Ersp_state: ', get_me_dim_ersp_state(me, 0), &
            get_me_dim_ersp_state(me, 1), get_me_dim_ersp_state(me, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'iErsp_state: ', get_me_iersp_state(me, 0), &
            get_me_iersp_state(me, 1), get_me_iersp_state(me, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'dim_Ersp_sys: ', get_me_dim_ersp_sys(me, 0), &
            get_me_dim_ersp_sys(me, 1), get_me_dim_ersp_sys(me, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'iErsp_sys: ', get_me_iersp_sys(me, 0), &
            get_me_iersp_sys(me, 1), get_me_iersp_sys(me, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'dim_Brsp_sys: ', get_me_dim_brsp_sys(me, 0), &
            get_me_dim_brsp_sys(me, 1), get_me_dim_brsp_sys(me, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'iBrsp_sys: ', get_me_ibrsp_sys(me, 0), &
            get_me_ibrsp_sys(me, 1), get_me_ibrsp_sys(me, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'jr der_orders: ', get_me_der_order(me, 0, 0), &
            get_me_der_order(me, 0, 1), get_me_der_order(me, 0, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'js der_orders: ', get_me_der_order(me, 1, 0), &
            get_me_der_order(me, 1, 1), get_me_der_order(me, 1, 2)
        write (*, '(a,i0,1x,i0,1x,i0)') 'jp der_orders: ', get_me_der_order(me, 2, 0), &
            get_me_der_order(me, 2, 1), get_me_der_order(me, 2, 2)

        Nw = get_me_nwaves(me)
        do k = 0, Nw - 1
            write (*, '(a,i0,a,i0)') 'sys_ind[', k, '] = ', get_me_sys_ind(me, k)
        end do
    end subroutine print_maxwell_eqs_data_local

    subroutine to_cstr(str, buf)
        character(len=*), intent(in) :: str
        character(kind=c_char), intent(out) :: buf(*)
        integer :: i, n
        n = len_trim(str)
        do i = 1, n
            buf(i) = str(i:i)
        end do
        buf(n + 1) = c_null_char
    end subroutine to_cstr

    !> ---- bind(C) free-function callbacks (called BY legacy Fortran /
    !> still-C++ callers with a zone handle BY REFERENCE) ----

    subroutine set_equations_settings_c_(ptr, hom_sys, Nwaves, Nfs, Nphys, flag_debug) &
        bind(C, name="set_equations_settings_c_")
        integer(c_intptr_t), intent(in) :: ptr
        integer(c_int), intent(out) :: hom_sys, Nwaves, Nfs, Nphys, flag_debug
        class(zone_t), pointer :: z
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            hom_sys = z%hom_sys
            Nwaves = z%Nwaves
            Nfs = z%Nfs
            Nphys = 0
            flag_debug = z%flag_debug
        end select
    end subroutine set_equations_settings_c_

    subroutine set_conductivity_settings_c_(ptr, flre_order, Nmax, gal_corr, rsp) &
        bind(C, name="set_conductivity_settings_c_")
        integer(c_intptr_t), intent(in) :: ptr
        integer(c_int), intent(out) :: flre_order, Nmax, gal_corr, rsp
        class(zone_t), pointer :: z
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            flre_order = z%flre_order
            Nmax = z%Nmax
            gal_corr = z%gal_corr
            rsp = z%rsp
        end select
    end subroutine set_conductivity_settings_c_

    subroutine set_collisions_settings_c_(ptr, collmod) bind(C, name="set_collisions_settings_c_")
        integer(c_intptr_t), intent(in) :: ptr
        integer(c_int), intent(out) :: collmod(2)
        class(zone_t), pointer :: z
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            collmod = z%collmod
        end select
    end subroutine set_collisions_settings_c_

    subroutine get_iersp_sys_array_(ptr, iErsp_sys) bind(C, name="get_iersp_sys_array_")
        integer(c_intptr_t), intent(in) :: ptr
        integer(c_int), intent(out) :: iErsp_sys(*)
        class(zone_t), pointer :: z
        integer :: k
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            do k = 0, 2
                iErsp_sys(k + 1) = get_me_iersp_sys(z%me, k) + 1
            end do
        end select
    end subroutine get_iersp_sys_array_

    subroutine get_ibrsp_sys_array_(ptr, iBrsp_sys) bind(C, name="get_ibrsp_sys_array_")
        integer(c_intptr_t), intent(in) :: ptr
        integer(c_int), intent(out) :: iBrsp_sys(*)
        class(zone_t), pointer :: z
        integer :: k
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            do k = 0, 2
                iBrsp_sys(k + 1) = get_me_ibrsp_sys(z%me, k) + 1
            end do
        end select
    end subroutine get_ibrsp_sys_array_

    subroutine get_sys_ind_array_(ptr, sys_ind) bind(C, name="get_sys_ind_array_")
        integer(c_intptr_t), intent(in) :: ptr
        integer(c_int), intent(out) :: sys_ind(*)
        class(zone_t), pointer :: z
        integer :: k
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            do k = 0, z%Nwaves - 1
                sys_ind(k + 1) = get_me_sys_ind(z%me, k) + 1
            end do
        end select
    end subroutine get_sys_ind_array_

    subroutine activate_fortran_modules_for_zone_(ptr) &
        bind(C, name="activate_fortran_modules_for_zone_")
        integer(c_intptr_t), intent(in) :: ptr
        class(zone_t), pointer :: z
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            call flre_activate_fortran_modules(z)
        end select
    end subroutine activate_fortran_modules_for_zone_

    subroutine deactivate_fortran_modules_for_zone_(ptr) &
        bind(C, name="deactivate_fortran_modules_for_zone_")
        integer(c_intptr_t), intent(in) :: ptr
        class(zone_t), pointer :: z
        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            call flre_deactivate_fortran_modules(z)
        end select
    end subroutine deactivate_fortran_modules_for_zone_

    !> rpt/EB1/EB2: a single radial point and the FULL set of Nwaves basis
    !> system vectors there (each Ncomps components), transformed from the
    !> moving frame to the lab cylindrical frame -- matches the legacy
    !> Fortran caller's `complex(8), dimension(D1,Nw1)` convention
    !> (flre.f90's stitching_equations_flre_N1_hommed), so EB1/EB2 here are
    !> COMPLEX, not flat real, despite the oracle's `double *` signature.
    subroutine calc_flre_basis_in_lab_cyl_frame_with_full_system_vectors_(ptr, rpt, EB1, EB2) &
        bind(C, name="calc_flre_basis_in_lab_cyl_frame_with_full_system_vectors_")
        integer(c_intptr_t), intent(in) :: ptr
        real(dp), intent(in) :: rpt
        complex(dp), intent(in) :: EB1(*)
        complex(dp), intent(out) :: EB2(*)
        class(zone_t), pointer :: z
        complex(dp), allocatable :: rhs(:), state(:), sys_mov(:), sys_lab(:)
        integer :: num_eqs, k, base
        character(kind=c_char) :: flag_back

        call handle_to_zone(ptr, z)
        select type (z)
        type is (flre_zone_t)
            num_eqs = get_me_num_eqs(z%me)
            allocate (rhs(num_eqs))
            rhs = (0.0_dp, 0.0_dp)
            allocate (state(z%Nwaves), sys_mov(z%Ncomps), sys_lab(z%Ncomps))

            call flre_activate_fortran_modules(z)

            flag_back = get_background_flag_back()

            do k = 0, z%Nwaves - 1
                base = k*z%Ncomps
                call system_to_state_copy(z, EB1(base + 1:base + z%Ncomps), state)
                call state2sys(rpt, flag_back, state, sys_mov, rhs)
                call galilean_transform_of_flre_system_vector(z, -get_background_V_gal_sys(), &
                                                                z%wd%omov, rpt, sys_mov, sys_lab)
                call transform_of_flre_system_vector_to_cyl_coordinates( &
                    z, rpt, sys_lab, EB2(base + 1:base + z%Ncomps))
            end do

            call flre_deactivate_fortran_modules(z)

            deallocate (rhs, state, sys_mov, sys_lab)
        end select
    end subroutine calc_flre_basis_in_lab_cyl_frame_with_full_system_vectors_

    !> wave_code_interface.cpp's get_kilca_conductivity_array_ used to
    !> `static_cast<flre_zone*>` a zone* and read ->flre_order/->cp
    !> directly; with the hierarchy now Fortran, it gets the same data
    !> through these two getters instead.
    function flre_zone_get_flre_order_(handle) bind(C, name="flre_zone_get_flre_order_") result(res)
        integer(c_intptr_t), value :: handle
        integer(c_int) :: res
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        select type (z)
        type is (flre_zone_t)
            res = z%flre_order
        end select
    end function flre_zone_get_flre_order_

    function flre_zone_get_cp_(handle) bind(C, name="flre_zone_get_cp_") result(res)
        integer(c_intptr_t), value :: handle
        integer(c_intptr_t) :: res
        class(zone_t), pointer :: z
        call handle_to_zone(handle, z)
        select type (z)
        type is (flre_zone_t)
            res = z%cp
        end select
    end function flre_zone_get_cp_

end module kilca_flre_zone_m
