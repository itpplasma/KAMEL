!> Unit test for kilca_maxwell_eqs_data_m, isolating the trickiest part of the
!> translation: der_order's flat C row-major [3][3] layout vs. Fortran 2D
!> indexing. Provides fake implementations of the three native data-fill
!> entry points (copy_module_data_to_maxwell_eqs_data_struct_f,
!> get_ersp_state_indices_and_dims_f, get_sys_ind_array_f) with known,
!> all-distinct values, so each getter's index mapping can be checked exactly
!> -- independent of having a real Maxwell-equations zone set up.
program test_maxwell_eqs_data
    use kilca_maxwell_eqs_data_m, only: get_me_der_order
    use kilca_maxwell_eqs_data_m, only: get_me_dim_brsp_sys
    use kilca_maxwell_eqs_data_m, only: get_me_dim_ersp_state
    use kilca_maxwell_eqs_data_m, only: get_me_dim_ersp_sys
    use kilca_maxwell_eqs_data_m, only: get_me_ibrsp_sys
    use kilca_maxwell_eqs_data_m, only: get_me_iersp_state
    use kilca_maxwell_eqs_data_m, only: get_me_iersp_sys
    use kilca_maxwell_eqs_data_m, only: get_me_num_eqs
    use kilca_maxwell_eqs_data_m, only: get_me_num_vars
    use kilca_maxwell_eqs_data_m, only: get_me_sys_ind
    use kilca_maxwell_eqs_data_m, only: maxwell_eqs_data_create
    use kilca_maxwell_eqs_data_m, only: maxwell_eqs_data_destroy
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t
    implicit none
    integer(c_intptr_t) :: handle
    integer(c_int) :: i, j, failures

    failures = 0

    handle = maxwell_eqs_data_create(3_c_int)

    call check_i("num_vars", get_me_num_vars(handle), 11)
    call check_i("num_eqs", get_me_num_eqs(handle), 22)

    ! der_order[i][j] (C 0-based) must equal the value the fake native writer
    ! placed at flat row-major offset i*3+j: 100 + i*3+j.
    do i = 0, 2
        do j = 0, 2
            call check_i("der_order", get_me_der_order(handle, i, j), 100 + i*3 + j)
        end do
    end do

    do i = 0, 2
        call check_i("dim_Ersp_sys", get_me_dim_ersp_sys(handle, i), 30 + i)
        call check_i("iErsp_sys", get_me_iersp_sys(handle, i), &
                     40 + i - 1) ! -1: Fortran->C index shift
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
