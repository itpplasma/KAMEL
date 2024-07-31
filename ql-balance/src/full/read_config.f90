subroutine read_config

    use wave_code_data, only: flre_path, vac_path, antenna_factor
    use baseparam_mod, only: btor, rtor, rsepar, dperp, Z_i, am
    use grid_mod, only: rmin, rmax, npoimin, gg_factor, gg_width, &
        gg_r_res, iboutype, rb_cut_in, re_cut_in, rb_cut_out, re_cut_out
    use control_mod, only: Nstorage, tmax_factor, iwrite, eps, icoll, &
        write_formfactors, flag_run_time_evolution, stop_time_step, &
        timstep_min, paramscan, save_prof_time_step, diagnostics_output, br_stopping, &
        suppression_mode, debug_mode, readfromtimestep, ramp_up_mode, &
        t_max_ramp_up, temperature_limit, antenna_max_stopping, gyro_current_study, &
        misalign_diffusion, equil_path
    use h5mod, only: path2inp, path2out, path2time
    use paramscan_mod, only: viscosity_factor

    implicit none

    !> namelist for the balance configuration input
    NAMELIST /BALANCENML/ flre_path, vac_path, btor, rtor, rmin, rmax, &
        rsepar, npoimin, gg_factor, gg_width, gg_r_res, Nstorage, &
        tmax_factor, antenna_factor, iboutype, iwrite, eps, dperp, &
        icoll, Z_i, am, rb_cut_in, re_cut_in, rb_cut_out, re_cut_out, &
        write_formfactors, flag_run_time_evolution, stop_time_step, &
        path2inp, path2out, timstep_min, paramscan, save_prof_time_step, &
        diagnostics_output, br_stopping, suppression_mode, debug_mode, &
        readfromtimestep, path2time, ramp_up_mode, t_max_ramp_up, temperature_limit, &
        antenna_max_stopping, gyro_current_study, viscosity_factor, misalign_diffusion, &
        equil_path

    ! read the parameters from namelist file
    open (22, file='balance_conf.nml');
    read (22, NML=BALANCENML)
    close (22);


    write(*,*) ''
    write(*,*) 'Balance code'
    write(*,*) '====================================================================='
    write(*,*) 'Run time evolution: ', flag_run_time_evolution
    write(*,*) ''
    write(*,*) 'Parameters from balance_conf.nml:'
    write(*,*) '---------------------------------------------------------------------'
    write(*,*) 'flre path: ', trim(flre_path)
    write(*,*) 'vac path: ', trim(vac_path)
    write(*,*) 'B_tor = ', btor
    write(*,*) 'R_tor = ', rtor
    write(*,*) 'r_min = ', rmin
    write(*,*) 'r_max = ', rmax
    write(*,*) 'npoimin = ', npoimin
    write(*,*) 'gg_factor = ', gg_factor
    write(*,*) 'gg_width = ', gg_width
    write(*,*) 'gg_r_res = ', gg_r_res
    write(*,*) 'Nstorage = ', Nstorage
    write(*,*) 'tmax_factor = ', tmax_factor
    write(*,*) "timstep_min = ", timstep_min
    write(*,*) 'antenna_factor = ', antenna_factor
    write(*,*) 'iboutype = ', iboutype
    write(*,*) 'iwrite = ', iwrite
    write(*,*) 'eps = ', eps
    write(*,*) 'dperp = ', dperp
    write(*,*) 'icoll = ', icoll
    write(*,*) 'Z_i = ', Z_i
    write(*,*) 'am = ', am
    write(*,*) 'stop_time_step = ', stop_time_step
    write(*,*) 'path2inp = ', trim(path2inp)
    write(*,*) 'path2out = ', trim(path2out)
    write(*,*) 'paramscan = ', paramscan
    write(*,*) 'diagnostics_output = ', diagnostics_output
    write(*,*) 'br_stopping = ', br_stopping
    write(*,*) 'debug_mode = ', debug_mode
    write(*,*) 'readfromtimestep = ', readfromtimestep
    write(*,*) 'suppression_mode = ', suppression_mode
    write(*,*) "ramp_up_mode = ", ramp_up_mode
    write(*,*) "t_max_ramp_up = ", t_max_ramp_up
    write(*,*) "temperature_limit = ", temperature_limit
    write(*,*) "antenna_max_stopping = ", antenna_max_stopping
    write(*,*) "gyro_current_study = ", gyro_current_study
    write(*,*) "viscosity_factor = ", viscosity_factor
    write(*,*) "misalign_diffusion = ", misalign_diffusion
    write(*,*) "equil_path = ", trim(equil_path)
    write(*,*) ''
    write(*,*) '====================================================================='


end subroutine