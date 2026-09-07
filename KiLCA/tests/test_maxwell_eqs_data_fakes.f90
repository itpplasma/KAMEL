!> Minimal domain parameters for this isolated native-interface test.
module constants
    use iso_c_binding, only: c_double, c_double_complex, c_intptr_t
    implicit none
    integer, parameter :: dp = c_double, dpc = c_double_complex, pp = c_intptr_t
end module constants

module flre_sett
    implicit none
    integer :: nwaves = 3, flre_order = 1
end module flre_sett

subroutine copy_module_data_to_maxwell_eqs_data_struct_f(num_vars_p, num_eqs_p, &
                                                         dim_Ersp_sys_p, &
                                              iErsp_sys_p, dim_Brsp_sys_p, iBrsp_sys_p, der_order_p)
    use, intrinsic :: iso_c_binding, only: c_int
    integer(c_int), intent(out) :: num_vars_p, num_eqs_p
    integer(c_int), intent(out) :: dim_Ersp_sys_p(3), iErsp_sys_p(3)
    integer(c_int), intent(out) :: dim_Brsp_sys_p(3), iBrsp_sys_p(3)
    integer(c_int), intent(out) :: der_order_p(3, 3)
    integer(c_int) :: k, i, j

    num_vars_p = 11
    num_eqs_p = 22
    do k = 1, 3
        dim_Ersp_sys_p(k) = 30 + (k - 1)
        iErsp_sys_p(k) = 40 + (k - 1)
        dim_Brsp_sys_p(k) = 50 + (k - 1)
        iBrsp_sys_p(k) = 60 + (k - 1)
    end do
    ! C row-major [3][3]: flat offset i*3+j holds 100+i*3+j (i,j 0-based).
    do j = 1, 3
        do i = 1, 3
            der_order_p(i, j) = 100 + (j - 1) * 3 + i - 1
        end do
    end do
end subroutine

subroutine get_ersp_state_indices_and_dims_f(dim_Ersp_state_p, iErsp_state_p)
    use, intrinsic :: iso_c_binding, only: c_int
    integer(c_int), intent(out) :: dim_Ersp_state_p(3), iErsp_state_p(3)
    integer(c_int) :: k
    do k = 1, 3
        dim_Ersp_state_p(k) = 70 + (k - 1)
        iErsp_state_p(k) = 80 + (k - 1)
    end do
end subroutine

subroutine get_sys_ind_array_f(sys_ind_p)
    use, intrinsic :: iso_c_binding, only: c_int
    use flre_sett, only: nwaves
    integer(c_int), intent(out) :: sys_ind_p(nwaves)
    integer(c_int) :: k
    do k = 1, nwaves
        sys_ind_p(k) = 90 + (k - 1)
    end do
end subroutine
