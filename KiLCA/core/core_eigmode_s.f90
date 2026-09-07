!> Eigenmode orchestration is implemented separately to keep core/search imports acyclic.
submodule(kilca_core_data_m) kilca_core_eigmode_s
    use kilca_eigmode_solve_m, only: loop_over_frequences, find_det_zeros, &
                                     find_eigmodes
    use, intrinsic :: iso_c_binding, only: c_null_ptr
    implicit none
contains
    module procedure core_data_calc_and_set_mode_dependent_eigmode_
    type(core_data_t), pointer :: cd
    integer :: ind, m, n, stat_unused

    call c_f_pointer(transfer(handle, c_null_ptr), cd)

    cd%dim = get_antenna_dma()
    if (allocated(cd%mda)) deallocate (cd%mda)
    allocate (cd%mda(cd%dim))
    cd%mda = 0

    do ind = 1, cd%dim
        call get_antenna_mode(int(ind - 1, c_int), m, n)

        select case (get_eigmode_search_flag())
        case (1)
            stat_unused = loop_over_frequences(int(ind - 1, c_int), m, n, handle)
        case (0)
            stat_unused = find_det_zeros(int(ind - 1, c_int), m, n, handle)
        case (-1)
            stat_unused = find_eigmodes(int(ind - 1, c_int), m, n, handle)
        case default
            write (*, '(a,i0,a)') 'Error: unknown search_flag in eigmode options file: ', &
                get_eigmode_search_flag(), '.'
            stop 1
        end select
    end do
    end procedure core_data_calc_and_set_mode_dependent_eigmode_
end submodule kilca_core_eigmode_s
