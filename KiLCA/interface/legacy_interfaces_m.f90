!> Native interfaces for legacy external Fortran procedures used by the ported modules.
!> Keep these declarations aligned with their defining routines. Strings use the
!> compiler's Fortran convention, and complex arrays retain their actual element types.
module kilca_legacy_interfaces_m
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use constants, only: dp, dpc, pp
    implicit none
    private

    public :: binomial_coefficients, get_pointer_precision
    public :: clear_all_data_in_mode_data_module, set_settings_in_core_module
    public :: eval_and_set_background_parameters_spec_independent
    public :: eval_and_set_background_parameters_spec_dependent
    public :: eval_and_set_wave_parameters, eval_and_set_f0_parameters_nu_and_derivs
    public :: eval_a_matrix, calc_dem_djmi_arrays, eval_electric_drift_velocities
    public :: eval_fgi_arrays, calc_w2_array, calc_d_array, calc_k_matrices, &
              calc_k1_matrices
    public :: calc_and_add_galilelian_correction, get_flre_order, get_gal_corr
    public :: current_density, cyl2rsp, calc_current_density_r_s_p
    public :: get_wave_parameters, get_magnetic_field_parameters, &
              calc_diff_sys_matrix
    public :: copy_module_data_to_maxwell_eqs_data_struct_f
    public :: get_ersp_state_indices_and_dims_f, get_sys_ind_array_f
    public :: get_background_dimension_from_balance, &
              get_background_profiles_from_balance

    public :: set_core_data_in_core_module, set_wave_parameters_in_mode_data_module
    interface
        subroutine set_core_data_in_core_module(cd)
            import :: pp
            integer(pp), intent(in) :: cd
        end subroutine
        subroutine set_wave_parameters_in_mode_data_module(m, n, olab_re, olab_im, &
                                                           omov_re, omov_im)
            import :: dp
            integer, intent(in) :: m, n
            real(dp), intent(in) :: olab_re, olab_im, omov_re, omov_im
        end subroutine
        subroutine binomial_coefficients(n, bc)
            import :: c_int, c_double
            integer(c_int), value :: n
            real(c_double), intent(out) :: bc(0:n, 0:n)
        end subroutine binomial_coefficients

        subroutine get_pointer_precision(ppp)
            integer, intent(out) :: ppp
        end subroutine get_pointer_precision

        subroutine clear_all_data_in_mode_data_module()
        end subroutine clear_all_data_in_mode_data_module

        subroutine set_settings_in_core_module(sd)
            import :: pp
            integer(pp), intent(in) :: sd
        end subroutine set_settings_in_core_module

        subroutine eval_and_set_background_parameters_spec_independent(r, flag_back)
            import :: dp
            real(dp), intent(in) :: r
            character(len=*), intent(in) :: flag_back
        end subroutine eval_and_set_background_parameters_spec_independent

        subroutine eval_and_set_background_parameters_spec_dependent(r, spec, &
                                                                     flag_back)
            import :: dp
            real(dp), intent(in) :: r
            integer, intent(in) :: spec
            character(len=*), intent(in) :: flag_back
        end subroutine eval_and_set_background_parameters_spec_dependent

        subroutine eval_and_set_wave_parameters(r, flag_back)
            import :: dp
            real(dp), intent(in) :: r
            character(len=*), intent(in) :: flag_back
        end subroutine eval_and_set_wave_parameters

        subroutine eval_and_set_f0_parameters_nu_and_derivs(r, spec, flag_back)
            import :: dp
            real(dp), intent(in) :: r
            integer, intent(in) :: spec
            character(len=*), intent(in) :: flag_back
        end subroutine eval_and_set_f0_parameters_nu_and_derivs

        subroutine eval_a_matrix()
        end subroutine eval_a_matrix

        subroutine calc_dem_djmi_arrays(r)
            import :: dp
            real(dp), intent(in) :: r
        end subroutine calc_dem_djmi_arrays

        subroutine eval_electric_drift_velocities()
        end subroutine eval_electric_drift_velocities

        subroutine eval_fgi_arrays()
        end subroutine eval_fgi_arrays

        subroutine calc_w2_array(spec)
            integer, intent(in) :: spec
        end subroutine calc_w2_array

        subroutine calc_d_array()
        end subroutine calc_d_array

        subroutine calc_k_matrices(kmat)
            use flre_sett, only: flre_order
            import :: dpc
            complex(dpc), intent(out) :: kmat(1:3, 1:3, 0:flre_order, 0:flre_order)
        end subroutine calc_k_matrices

        subroutine calc_k1_matrices(kmat)
            use flre_sett, only: flre_order
            import :: dpc
            complex(dpc), intent(out) :: kmat(1:3, 1:3, 0:flre_order, 0:flre_order)
        end subroutine calc_k1_matrices

        subroutine calc_and_add_galilelian_correction(r, spec, flag_back, ct)
            use flre_sett, only: flre_order
            import :: dp, dpc
            real(dp), intent(in) :: r
            integer, intent(in) :: spec
            character(len=1), intent(in) :: flag_back
            complex(dpc), intent(inout) :: ct(1:3, 1:3, 0:2 * flre_order)
        end subroutine calc_and_add_galilelian_correction

        subroutine get_flre_order(flre_order_p)
            integer, intent(out) :: flre_order_p
        end subroutine get_flre_order

        subroutine get_gal_corr(gal_corr_p)
            integer, intent(out) :: gal_corr_p
        end subroutine get_gal_corr

        subroutine current_density(j_surf)
            import :: dpc
            complex(dpc), intent(out) :: j_surf(2)
        end subroutine current_density

        subroutine cyl2rsp(r, th, z, per, par)
            import :: dpc
            real(dpc), intent(in) :: r
            complex(dpc), intent(in) :: th, z
            complex(dpc), intent(out) :: per, par
        end subroutine cyl2rsp

        subroutine calc_current_density_r_s_p(r, ja_rsp)
            import :: dp, dpc
            real(dp), intent(in) :: r
            complex(dpc), intent(out) :: ja_rsp(2)
        end subroutine calc_current_density_r_s_p

        subroutine get_wave_parameters(kvals)
            import :: dp
            real(dp), intent(out) :: kvals(3)
        end subroutine get_wave_parameters

        subroutine get_magnetic_field_parameters(hvals)
            import :: dp
            real(dp), intent(out) :: hvals(3)
        end subroutine get_magnetic_field_parameters

        subroutine calc_diff_sys_matrix(r, flagback, dmat)
            use flre_sett, only: nwaves
            import :: dp, dpc
            real(dp), intent(in) :: r
            character(len=*), intent(in) :: flagback
            complex(dpc), intent(out) :: dmat(nwaves, nwaves)
        end subroutine calc_diff_sys_matrix

        subroutine copy_module_data_to_maxwell_eqs_data_struct_f(num_vars_p, &
                                                                 num_eqs_p, &
                              dim_ersp_sys_p, iersp_sys_p, dim_brsp_sys_p, ibrsp_sys_p, der_order_p)
            integer, intent(out) :: num_vars_p, num_eqs_p
            integer, intent(out) :: dim_ersp_sys_p(3), iersp_sys_p(3)
            integer, intent(out) :: dim_brsp_sys_p(3), ibrsp_sys_p(3)
            integer, intent(out) :: der_order_p(3, 3)
        end subroutine copy_module_data_to_maxwell_eqs_data_struct_f

        subroutine get_ersp_state_indices_and_dims_f(dim_ersp_state_p, iersp_state_p)
            integer, intent(out) :: dim_ersp_state_p(3), iersp_state_p(3)
        end subroutine get_ersp_state_indices_and_dims_f

        subroutine get_sys_ind_array_f(sys_ind_p)
            use flre_sett, only: nwaves
            integer, intent(out) :: sys_ind_p(nwaves)
        end subroutine get_sys_ind_array_f

        subroutine get_background_dimension_from_balance(dim_p)
            integer, intent(out) :: dim_p
        end subroutine get_background_dimension_from_balance

        !> Callers obtain dim_p from get_background_dimension_from_balance first.
        !> The implementation supplies the extent through that preceding call.
        !> Assumed-size vectors match the explicit-shape wave_code_data API without
        !> requiring a dependency on its global profile dimension.
        subroutine get_background_profiles_from_balance(dim_p, r_p, q_p, n_p, ti_p, &
                                                        te_p, &
                                                        vth_p, vz_p, er_p)
            integer, intent(out) :: dim_p
            real(8), intent(out) :: r_p(*), q_p(*), n_p(*), ti_p(*), te_p(*)
            real(8), intent(out) :: vth_p(*), vz_p(*), er_p(*)
        end subroutine get_background_profiles_from_balance
    end interface
end module kilca_legacy_interfaces_m
