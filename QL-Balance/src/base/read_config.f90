subroutine read_config
    use baseparam_mod, only: btor, rtor, rsepar, dperp, Z_i, am, urelax
    use control_mod, only: eps, paramscan, diagnostics_output, suppression_mode, debug_mode, &
                           readfromtimestep, temperature_limit, gyro_current_study, &
                           misalign_diffusion, equil_path, ihdf5IO, wave_code
    use grid_mod, only: rmin, rmax, npoimin, gg_factor, gg_width, gg_r_res, iboutype, rb_cut_in, &
                        re_cut_in, rb_cut_out, re_cut_out
    use h5mod, only: path2inp, path2out, path2time
    use paramscan_mod, only: viscosity_factor
    use time_evolution
    use wave_code_data, only: flre_path, vac_path, antenna_factor

    implicit none

    integer :: u, ios

    character(len=*), parameter :: config_file = "balance_conf.nml"

    !> namelist for the balance configuration input
    namelist /BALANCENML/ flre_path, vac_path, btor, rtor, rmin, rmax, rsepar, npoimin, gg_factor, &
        gg_width, gg_r_res, Nstorage, tmax_factor, antenna_factor, iboutype, eps, dperp, Z_i, am, &
        rb_cut_in, re_cut_in, rb_cut_out, re_cut_out, stop_time_step, path2inp, path2out, &
        timstep_min, paramscan, save_prof_time_step, diagnostics_output, br_stopping, &
        suppression_mode, debug_mode, readfromtimestep, path2time, ramp_up_mode, t_max_ramp_up, &
        temperature_limit, antenna_max_stopping, gyro_current_study, viscosity_factor, &
        misalign_diffusion, equil_path, ihdf5IO, type_of_run, wave_code, &
        set_constant_time_step, constant_time_step, urelax

    ! read the parameters from namelist file
    open (newunit=u, file=config_file, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "Failed to open config file"
    read (u, nml=BALANCENML, iostat=ios)
    if (ios /= 0) error stop "Failed to read namelist"
    close (u)

    write (*, *) ""
    write (*, *) "================================================================================="
    write (*, "(A,A)") "    Type of Run: ", type_of_run
    write (*, "(A,A)") "    Wave Code:   ", trim(adjustl(wave_code))
    write (*, *) ""
    write (*, "(A,A,A)") "    Parameters from ", config_file, ":"
    write (*, *) "   ------------------------------------------------------------------------------"
    write (*, "(A,A)") "    flre path: ", trim(adjustl(flre_path))
    write (*, "(A,A)") "    vac path: ", trim(adjustl(vac_path))
    write (*, "(A,ES15.8,A)") "    B_tor = ", btor, " G"
    write (*, "(A,ES15.8,A)") "    R_tor = ", rtor, " cm"
    write (*, "(A,ES15.8,A)") "    r_min = ", rmin, " cm"
    write (*, "(A,ES15.8,A)") "    r_max = ", rmax, " cm"
    write (*, "(A,I0)") "    npoimin = ", npoimin
    write (*, "(A,ES15.8)") "    gg_factor = ", gg_factor
    write (*, "(A,ES15.8)") "    gg_width = ", gg_width
    write (*, "(A,ES15.8)") "    gg_r_res = ", gg_r_res
    write (*, "(A,I0)") "    Nstorage = ", Nstorage
    write (*, "(A,ES15.8)") "    tmax_factor = ", tmax_factor
    write (*, "(A,ES15.8)") "    timstep_min = ", timstep_min
    write (*, "(A,ES15.8)") "    antenna_factor = ", antenna_factor
    write (*, "(A,I0)") "    iboutype = ", iboutype
    write (*, "(A,ES15.8)") "    eps = ", eps
    write (*, "(A,ES15.8)") "    dperp = ", dperp
    write (*, "(A,ES15.8)") "    Z_i = ", Z_i
    write (*, "(A,ES15.8)") "    am = ", am
    write (*, "(A,ES12.4)") "    stop_time_step = ", stop_time_step
    write (*, "(A,A)") "    path2inp = ", trim(adjustl(path2inp))
    write (*, "(A,A)") "    path2out = ", trim(adjustl(path2out))
    write (*, "(A,L0)") "    paramscan = ", paramscan
    write (*, "(A,L0)") "    diagnostics_output = ", diagnostics_output
    write (*, "(A,L0)") "    br_stopping = ", br_stopping
    write (*, "(A,L0)") "    debug_mode = ", debug_mode
    write (*, "(A,I0)") "    readfromtimestep = ", readfromtimestep
    write (*, "(A,L0)") "    suppression_mode = ", suppression_mode
    write (*, "(A,I0)") "    ramp_up_mode = ", ramp_up_mode
    write (*, "(A,ES15.8,A)") "    t_max_ramp_up = ", t_max_ramp_up, " s"
    write (*, "(A,ES15.8,A)") "    temperature_limit = ", temperature_limit, " eV"
    write (*, "(A,ES15.8)") "    antenna_max_stopping = ", antenna_max_stopping
    write (*, "(A,I0)") "    gyro_current_study = ", gyro_current_study
    write (*, "(A,ES15.8)") "    viscosity_factor = ", viscosity_factor
    write (*, "(A,L0)") "    misalign_diffusion = ", misalign_diffusion
    write (*, "(A,A)") "    equil_path = ", trim(adjustl(equil_path))
    write (*, "(A,I0)") "    ihdf5IO = ", ihdf5IO
    write (*, "(A,L0)") "    set_constant_time_step = ", set_constant_time_step
    write (*, "(A,ES15.8,A)") "    constant_time_step = ", constant_time_step, " s"
    write (*, "(A,ES15.8)") "    urelax = ", urelax
    write (*, *) ""
    write (*, *) "================================================================================="
end subroutine read_config
