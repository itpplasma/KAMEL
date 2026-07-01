!> Unit test for kilca_maxwell_eqs_data_m, isolating the trickiest part of the
!> translation: der_order's flat C row-major [3][3] layout vs. Fortran 2D
!> indexing. Provides fake implementations of the three C-ABI data-fill
!> entry points (copy_module_data_to_maxwell_eqs_data_struct_f_,
!> get_ersp_state_indices_and_dims_f_, get_sys_ind_array_f_) with known,
!> all-distinct values, so each getter's index mapping can be checked exactly
!> -- independent of having a real Maxwell-equations zone set up.
program test_maxwell_eqs_data
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t
    implicit none

    interface
        function maxwell_eqs_data_create(Nwaves) result(handle) &
            bind(C, name="maxwell_eqs_data_create_")
            import :: c_int, c_intptr_t
            integer(c_int), value :: Nwaves
            integer(c_intptr_t) :: handle
        end function
        subroutine maxwell_eqs_data_destroy(handle) bind(C, name="maxwell_eqs_data_destroy_")
            import :: c_intptr_t
            integer(c_intptr_t), value :: handle
        end subroutine
        integer(c_int) function get_me_num_vars(handle) bind(C, name="get_me_num_vars_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
        end function
        integer(c_int) function get_me_num_eqs(handle) bind(C, name="get_me_num_eqs_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
        end function
        integer(c_int) function get_me_der_order(handle, i, j) bind(C, name="get_me_der_order_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: i, j
        end function
        integer(c_int) function get_me_dim_ersp_state(handle, k) bind(C, name="get_me_dim_ersp_state_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: k
        end function
        integer(c_int) function get_me_iersp_state(handle, k) bind(C, name="get_me_iersp_state_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: k
        end function
        integer(c_int) function get_me_dim_ersp_sys(handle, k) bind(C, name="get_me_dim_ersp_sys_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: k
        end function
        integer(c_int) function get_me_iersp_sys(handle, k) bind(C, name="get_me_iersp_sys_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: k
        end function
        integer(c_int) function get_me_dim_brsp_sys(handle, k) bind(C, name="get_me_dim_brsp_sys_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: k
        end function
        integer(c_int) function get_me_ibrsp_sys(handle, k) bind(C, name="get_me_ibrsp_sys_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: k
        end function
        integer(c_int) function get_me_sys_ind(handle, k) bind(C, name="get_me_sys_ind_")
            import :: c_int, c_intptr_t
            integer(c_intptr_t), value :: handle
            integer(c_int), value :: k
        end function
    end interface

    integer(c_intptr_t) :: handle
    integer(c_int) :: i, j, failures

    failures = 0

    handle = maxwell_eqs_data_create(3_c_int)

    call check_i("num_vars", get_me_num_vars(handle), 11)
    call check_i("num_eqs", get_me_num_eqs(handle), 22)

    ! der_order[i][j] (C 0-based) must equal the value the fake C-ABI writer
    ! placed at flat row-major offset i*3+j: 100 + i*3+j.
    do i = 0, 2
        do j = 0, 2
            call check_i("der_order", get_me_der_order(handle, i, j), 100 + i*3 + j)
        end do
    end do

    do i = 0, 2
        call check_i("dim_Ersp_sys", get_me_dim_ersp_sys(handle, i), 30 + i)
        call check_i("iErsp_sys", get_me_iersp_sys(handle, i), 40 + i - 1)  ! -1: Fortran->C index shift
        call check_i("dim_Brsp_sys", get_me_dim_brsp_sys(handle, i), 50 + i)
        call check_i("iBrsp_sys", get_me_ibrsp_sys(handle, i), 60 + i - 1)
        call check_i("dim_Ersp_state", get_me_dim_ersp_state(handle, i), 70 + i)
        call check_i("iErsp_state", get_me_iersp_state(handle, i), 80 + i - 1)
    end do

    do i = 0, 2
        call check_i("sys_ind", get_me_sys_ind(handle, i), 90 + i - 1)
    end do

    call maxwell_eqs_data_destroy(handle)

    if (failures == 0) then
        write (*, '(a)') "PASS: kilca_maxwell_eqs_data_m getters match expected layout"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if

contains

    subroutine check_i(label, got, want)
        character(*), intent(in) :: label
        integer(c_int), intent(in) :: got, want
        if (got /= want) then
            write (*, '(a,a,2(1x,i0))') "FAIL ", label, got, want
            failures = failures + 1
        end if
    end subroutine

end program test_maxwell_eqs_data
