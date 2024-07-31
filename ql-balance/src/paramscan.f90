
subroutine getfactors

    use paramscan_mod
    use wave_code_data
    use h5mod

    integer :: lb, ub
    CALL h5_init()
    CALL h5_open_rw(path2out, h5_id)
    CALL h5_get_bounds_1(h5_id, "/factors/fac_n", lb, ub)
    write (*, *) "lower bound fac_n", lb, " upper bound ", ub
    allocate(fac_n(ub))
    CALL h5_get_bounds_1(h5_id, "/factors/fac_Te", lb, ub)
    write (*, *) "lower bound fac_Te", lb, " upper bound ", ub
    allocate(fac_Te(ub))
    CALL h5_get_bounds_1(h5_id, "/factors/fac_Ti", lb, ub)
    write (*, *) "lower bound fac_Ti", lb, " upper bound ", ub
    allocate(fac_Ti(ub))
    CALL h5_get_bounds_1(h5_id, "/factors/fac_vz", lb, ub)
    write (*, *) "lower bound fac_vz", lb, " upper bound ", ub
    allocate(fac_vz(ub))

    CALL h5_get_double_1(h5_id, "/factors/fac_n", fac_n)
    CALL h5_get_double_1(h5_id, "/factors/fac_Te", fac_Te)
    CALL h5_get_double_1(h5_id, "/factors/fac_Ti", fac_Ti)
    CALL h5_get_double_1(h5_id, "/factors/fac_vz", fac_vz)

    CALL h5_close(h5_id)
    CALL h5_deinit()

end subroutine
