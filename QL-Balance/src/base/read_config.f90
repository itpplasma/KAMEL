subroutine read_config
    use baseparam_mod, only: btor, rtor, rsepar, dperp, Z_i, am, urelax
    use control_mod, only: eps, paramscan, data_verbosity, suppression_mode, log_level, &
                           readfromtimestep, temperature_limit, gyro_current_study, &
                           misalign_diffusion, equil_path, ihdf5IO, wave_code, &
                           kim_config_path, kim_profiles_from_balance, &
                           kim_n_modes, kim_m_list, kim_n_list, &
                           jpar_method
    use grid_mod, only: rmin, rmax, npoimin, gg_factor, gg_width, gg_r_res, iboutype, rb_cut_in, &
                        re_cut_in, rb_cut_out, re_cut_out
    use h5mod, only: path2inp, path2out, path2time
    use paramscan_mod, only: viscosity_factor
    use time_evolution
    use wave_code_data, only: flre_path, vac_path, antenna_factor, I_par_toroidal
    use logger_m, only: set_log_level, log_info, fmt_val

    implicit none

    integer :: u, ios

    character(len=*), parameter :: config_file = "balance_conf.nml"

    !> namelist for the balance configuration input
    namelist /BALANCENML/ flre_path, vac_path, btor, rtor, rmin, rmax, rsepar, npoimin, gg_factor, &
        gg_width, gg_r_res, Nstorage, tmax_factor, antenna_factor, iboutype, eps, dperp, Z_i, am, &
        rb_cut_in, re_cut_in, rb_cut_out, re_cut_out, stop_time_step, path2inp, path2out, &
        timstep_min, paramscan, save_prof_time_step, data_verbosity, br_stopping, &
        suppression_mode, log_level, readfromtimestep, path2time, ramp_up_mode, t_max_ramp_up, &
        temperature_limit, antenna_max_stopping, gyro_current_study, viscosity_factor, &
        misalign_diffusion, equil_path, ihdf5IO, type_of_run, wave_code, &
        set_constant_time_step, constant_time_step, urelax, kim_config_path, &
        kim_profiles_from_balance, kim_n_modes, kim_m_list, kim_n_list, &
        I_par_toroidal, jpar_method

    ! read the parameters from namelist file
    open (newunit=u, file=config_file, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Failed to open config file"
    read (u, nml=BALANCENML, iostat=ios)
    if (ios /= 0) error stop "Failed to read namelist"
    close (u)

    call set_log_level(log_level)

    call log_info("=================================================================================")
    call log_info(fmt_val("    Type of Run", type_of_run))
    call log_info(fmt_val("    Wave Code", trim(adjustl(wave_code))))
    call log_info("")
    call log_info("    Parameters from " // config_file // ":")
    call log_info("   ------------------------------------------------------------------------------")
    call log_info(fmt_val("    flre path", trim(adjustl(flre_path))))
    call log_info(fmt_val("    vac path", trim(adjustl(vac_path))))
    call log_info(fmt_val("    B_tor", btor, "G"))
    call log_info(fmt_val("    R_tor", rtor, "cm"))
    call log_info(fmt_val("    r_min", rmin, "cm"))
    call log_info(fmt_val("    r_max", rmax, "cm"))
    call log_info(fmt_val("    npoimin", npoimin))
    call log_info(fmt_val("    gg_factor", gg_factor))
    call log_info(fmt_val("    gg_width", gg_width))
    call log_info(fmt_val("    gg_r_res", gg_r_res))
    call log_info(fmt_val("    Nstorage", Nstorage))
    call log_info(fmt_val("    tmax_factor", tmax_factor))
    call log_info(fmt_val("    timstep_min", timstep_min))
    call log_info(fmt_val("    antenna_factor", antenna_factor))
    call log_info(fmt_val("    I_par_toroidal", I_par_toroidal))
    call log_info(fmt_val("    iboutype", iboutype))
    call log_info(fmt_val("    eps", eps))
    call log_info(fmt_val("    dperp", dperp))
    call log_info(fmt_val("    Z_i", Z_i))
    call log_info(fmt_val("    am", am))
    call log_info(fmt_val("    stop_time_step", stop_time_step))
    call log_info(fmt_val("    path2inp", trim(adjustl(path2inp))))
    call log_info(fmt_val("    path2out", trim(adjustl(path2out))))
    call log_info(fmt_val("    paramscan", paramscan))
    call log_info(fmt_val("    data_verbosity", data_verbosity))
    call log_info(fmt_val("    br_stopping", br_stopping))
    call log_info(fmt_val("    log_level", log_level))
    call log_info(fmt_val("    readfromtimestep", readfromtimestep))
    call log_info(fmt_val("    suppression_mode", suppression_mode))
    call log_info(fmt_val("    ramp_up_mode", ramp_up_mode))
    call log_info(fmt_val("    t_max_ramp_up", t_max_ramp_up, "s"))
    call log_info(fmt_val("    temperature_limit", temperature_limit, "eV"))
    call log_info(fmt_val("    antenna_max_stopping", antenna_max_stopping))
    call log_info(fmt_val("    gyro_current_study", gyro_current_study))
    call log_info(fmt_val("    viscosity_factor", viscosity_factor))
    call log_info(fmt_val("    misalign_diffusion", misalign_diffusion))
    call log_info(fmt_val("    equil_path", trim(adjustl(equil_path))))
    call log_info(fmt_val("    ihdf5IO", ihdf5IO))
    call log_info(fmt_val("    set_constant_time_step", set_constant_time_step))
    call log_info(fmt_val("    constant_time_step", constant_time_step, "s"))
    call log_info(fmt_val("    urelax", urelax))
    call log_info(fmt_val("    jpar_method", trim(adjustl(jpar_method))))
    call log_info(fmt_val("    kim_config_path", trim(adjustl(kim_config_path))))
    call log_info(fmt_val("    kim_profiles_from_balance", kim_profiles_from_balance))
    call log_info(fmt_val("    kim_n_modes", kim_n_modes))
    if (kim_n_modes > 0) then
        write (*, "(A,100I5)") "    kim_m_list = ", kim_m_list(1:kim_n_modes)
        write (*, "(A,100I5)") "    kim_n_list = ", kim_n_list(1:kim_n_modes)
    end if
    call log_info("")
    call log_info("=================================================================================")
end subroutine read_config
