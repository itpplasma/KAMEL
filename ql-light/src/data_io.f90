module data_io

    use h5mod

    implicit none

    double precision, dimension(:), allocatable :: ne, Te, Ti, Vz, Er, q, r
    double precision, dimension(:), allocatable :: Dqle22

    contains

    subroutine read_input_data_h5

        integer :: lb, ub

        CALL h5_init()
        CALL h5_open_ro(path2inp, h5_id)
        CALL h5_get_bounds_1(h5_id, "/input_data/Te", lb, ub)
        write (*, *) "lower bound Te", lb, " upper bound ", ub
        allocate(Te(ub))
        CALL h5_get_double_1(h5_id, "/input_data/Te", Te)
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

    subroutine read_input_data_txt

        implicit none

    end subroutine

end module data_io