!> Fake implementations of the C-ABI data-fill entry points that
!> kilca_maxwell_eqs_data_m's constructor calls, used only by
!> test_maxwell_eqs_data so the getter index mapping can be checked with known
!> values without needing a real Maxwell-equations zone set up. Linked instead
!> of (not alongside) the real definitions in maxwell_eqs_m.f90.
subroutine copy_module_data_to_maxwell_eqs_data_struct_f(num_vars_p, num_eqs_p, &
    dim_Ersp_sys_p, iErsp_sys_p, dim_Brsp_sys_p, iBrsp_sys_p, der_order_p) &
    bind(C, name="copy_module_data_to_maxwell_eqs_data_struct_f_")
    use, intrinsic :: iso_c_binding, only: c_int
    integer(c_int), intent(out) :: num_vars_p, num_eqs_p
    integer(c_int), intent(out) :: dim_Ersp_sys_p(3), iErsp_sys_p(3)
    integer(c_int), intent(out) :: dim_Brsp_sys_p(3), iBrsp_sys_p(3)
    integer(c_int), intent(out) :: der_order_p(9)
    integer(c_int) :: k

    num_vars_p = 11
    num_eqs_p = 22
    do k = 1, 3
        dim_Ersp_sys_p(k) = 30 + (k - 1)
        iErsp_sys_p(k) = 40 + (k - 1)
        dim_Brsp_sys_p(k) = 50 + (k - 1)
        iBrsp_sys_p(k) = 60 + (k - 1)
    end do
    ! C row-major [3][3]: flat offset i*3+j holds 100+i*3+j (i,j 0-based).
    do k = 1, 9
        der_order_p(k) = 100 + (k - 1)
    end do
end subroutine

subroutine get_ersp_state_indices_and_dims_f(dim_Ersp_state_p, iErsp_state_p) &
    bind(C, name="get_ersp_state_indices_and_dims_f_")
    use, intrinsic :: iso_c_binding, only: c_int
    integer(c_int), intent(out) :: dim_Ersp_state_p(3), iErsp_state_p(3)
    integer(c_int) :: k
    do k = 1, 3
        dim_Ersp_state_p(k) = 70 + (k - 1)
        iErsp_state_p(k) = 80 + (k - 1)
    end do
end subroutine

subroutine get_sys_ind_array_f(sys_ind_p, n) bind(C, name="get_sys_ind_array_f_")
    use, intrinsic :: iso_c_binding, only: c_int
    integer(c_int), value :: n
    integer(c_int), intent(out) :: sys_ind_p(n)
    integer(c_int) :: k
    do k = 1, n
        sys_ind_p(k) = 90 + (k - 1)
    end do
end subroutine
