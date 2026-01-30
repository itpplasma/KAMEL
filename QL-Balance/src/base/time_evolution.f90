module time_evolution

    use control_mod
    use h5mod
    use balance_base, only: balance_t
    use QLBalance_kinds, only: dp

    implicit none

    logical :: br_stopping ! trigger Br stopping criterion
    logical :: discr_reached = .false. ! variable to say if discrepancy to linear regression
    logical :: scratch

    integer :: Nstorage
    integer :: ramp_up_mode !> control ramp up mode of the RMP coil current amplitude
    integer :: save_prof_time_step
    integer :: iexit ! used for ramp-up skipping of saving
    integer :: ramp_up_down = 0 !> used in hysteresis mode, tells if ramp-up (0) or ramp-down (1)
    integer :: time_ind

    real(dp) :: tmax, timescale
    real(dp) :: br_beta = 0
    real(dp) :: br_predicted

    real(dp) :: tmax_factor!, antenna_factor
    real(dp) :: stop_time_step
    real(dp) :: timstep_min
    real(dp) :: t_max_ramp_up = 1e-2 !> 10ms ramp up until antenna_factor_max is reached
    real(dp) :: timstep
    real(dp) :: time
    real(dp) :: t_hysteresis_turn = 0
    real(dp) :: constant_time_step
    logical :: set_constant_time_step

    real(dp), dimension(:), allocatable :: timscal

    real(dp), DIMENSION(:), ALLOCATABLE :: yprev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle11_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle12_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle21_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle22_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli11_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli12_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli21_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli22_prev

    real(dp) :: antenna_factor_max
    real(dp) :: antenna_max_stopping

    !needed for interpolation of br abs and stopping criterion
    real(dp), DIMENSION(:), ALLOCATABLE :: br_abs
    complex(dp), DIMENSION(:), ALLOCATABLE :: br_formfactor
    complex(dp), DIMENSION(:), ALLOCATABLE :: br_vac_res
    real(dp), DIMENSION(:), ALLOCATABLE :: br_abs_time
    real(dp), DIMENSION(:), ALLOCATABLE :: br_abs_antenna_factor
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle22_res_time
    real(dp), DIMENSION(:), ALLOCATABLE :: dae22_res_time
    real(dp), DIMENSION(:), ALLOCATABLE :: bif_criterion
    complex(dp), dimension(:), allocatable :: Ipar_time

    logical :: firstiterationdone = .false. !Some steps in saving the data to hdf5 file
    !need to be done only the first time iteration

    integer(HID_T) :: time_dataset_id !> variable to save the time dataset id

    !integer :: time_ind

    type, extends(balance_t) :: TimeEvolution_t
        contains
            procedure :: init_balance => initTimeEvolution
            procedure :: run_balance => runTimeEvolution
    end type

    private :: initTimeEvolution
    private :: runTimeEvolution

    contains

    subroutine initTimeEvolution(this)

        use recstep_mod, only: tol
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use grid_mod, only: mwind, rmax, rmin, set_boundary_condition, npoib, rb
        use baseparam_mod, only: dperp, tol_max
        use QLbalance_diag, only: write_diag, write_diag_b
        use KAMEL_hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, debug_mode, &
                        ihdf5IO
        use wave_code_data, only: m_vals, n_vals
        use plasma_parameters, only: write_initial_parameters, alloc_hold_parameters, &
                                params, params_begbeg, init_background_profiles
        use resonances_mod, only: write_resonant_radii_to_hdf5

        implicit none

        class(TimeEvolution_t), intent(inout) :: this
        this%runType = "TimeEvolution"

        iexit = 0 ! 0 - dont skip, 1 - skip, 2 - stop
        mwind = 10
        write_diag = .false.
        write_diag_b = .false.

        if (gyro_current_study .ne. 0) then
            write_gyro_current = .true.
        else
            write_gyro_current = .false.
        end if

        !call read_config

        timescale = (rmax - rmin)**2 / dperp
        tmax = timescale * tmax_factor
        timstep = tmax / Nstorage
        time = 0.0d0
        tol = tol_max

        call gengrid
        call set_boundary_condition

        call initialize_wave_code_interface(npoib, rb)

        mode_m = m_vals(1)
        mode_n = n_vals(1)
        if (ihdf5IO .eq. 1) then
            call create_group_structure_timeevol
        end if
        if (debug_mode) write(*,*) 'Debug: mode_m = ', mode_m, 'mode_n = ', mode_n

        call write_resonant_radii_to_hdf5

        call allocate_prev_variables
        call init_background_profiles
        call write_initial_parameters
        !call alloc_hold_parameters

        call calc_geometric_parameter_profiles
        call initialize_get_dql
        call initialize_antenna_factor
        call det_balance_eqs_source_terms

        if (.not. suppression_mode) call write_kin_prof_data_to_disk

        call allocate_timscal_and_params
        timstep = timstep*tol
        scratch = .true.

        call reset_timstep_arr_w_timstep
        !timstep_arr = timstep
        !tim_stack = timstep_arr

        call alloc_Br_Dqle_for_timeevol
        call get_dql
        call rescale_transp_coeffs_by_ant_fac
        call hold_prev_transp_coeffs

        params_begbeg = params

    end subroutine

    subroutine runTimeEvolution(this)
        class(TimeEvolution_t), intent(inout) :: this

        do time_ind = 1, Nstorage
            call doStep(this)
        end do
    end subroutine runTimeEvolution

    subroutine doStep(this)
        use baseparam_mod, only: factolmax, factolred
        use control_mod, only: debug_mode
        use plasma_parameters, only: params, params_beg, params_begbeg, limit_temps_from_below
        use recstep_mod, only: timstep_arr
        use recstep_mod, only: tol
        use restart_mod, only: redostep
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac

        implicit none

        class(TimeEvolution_t), intent(inout) :: this

        integer :: iredo

        call copy_kin_profs_to_yprev
        redostep = .false.

        call get_dql
        call rescale_transp_coeffs_by_ant_fac
        call interp_Br_Dql_at_resonance_timeevol
        call determine_Dql_diagnostic

        call write_br_dqle22_time_data
        call message_Br_Dqle_values

        if (diagnostics_output)then
            call writefort9999
        end if

        if (.true.) then
            call hold_prev_transp_coeffs
            params_begbeg = params
        else
            call redoTimeStep
        end if

        iredo = 0
        do ! redo step loop
            iredo = iredo + 1
            params_beg = params

            print *, ""
            if (debug_mode) write(*, "(a, i0, a, f12.6)") &
                            "Debug: Timstep before evolvestep is ", timstep, " eps = " , eps

            call evolvestep(timstep, eps)

            call limit_temps_from_below

            call calc_params_num_and_denom
            call smooth_params_num_and_denom
            call determine_timscal

            if (maxval(timscal) .lt. tol * factolmax) then
                exit
            end if

            timstep_arr = timstep_arr * factolred
            params = params_beg

            if (debug_mode) then
                print *, "Redoing step: Maxval(timscal) is not lesser than tol * factolmax"
                print *, "Maxval(timscal) = ", maxval(timscal)
                print *, "tol = ", tol
                print *, "factolmax = ", factolmax
                print *, "tol * factolmax = ", tol * factolmax
                print *, "timstep_arr(1) = ", timstep_arr(1)
                print *, "timstep_arr(100) = ", timstep_arr(100)
                print *, ""
            end if
            if (iredo > 100) then
                stop "Redoing step: Maxval(timscal) is not lesser than tol * factolmax " // &
                        "after 100 redos"
            end if
        end do

        call rescale_time_step_array
        call set_time_step
        call stop_if_time_step_too_small
        call reset_timstep_arr_w_timstep
        call write_time_info
        call relax_plasma_parameters

        timstep_arr = 0.0d0
        call evolvestep(timstep, eps)
        timstep_arr = timstep
        time = time + timstep

        if (debug_mode) call msg_time_info
        if (.not. suppression_mode) call write_kin_profile_at_time_index
        call set_first_iteration_true
        call calculate_total_toroidal_torque(time_ind)
        call write_total_toroidal_torque_to_file(time_ind)
        call check_linear_discr_pen_ratio
        call stop_if_antenna_fac_max_reached

        call ramp_coil
    end subroutine doStep

    subroutine allocate_prev_variables

        use recstep_mod, only: timstep_arr, tim_stack
        use grid_mod, only: neqset, npoib

        implicit none

        allocate (yprev(neqset))
        allocate (dqle11_prev(npoib))
        allocate (dqle12_prev(npoib))
        allocate (dqle21_prev(npoib))
        allocate (dqle22_prev(npoib))
        allocate (dqli11_prev(npoib))
        allocate (dqli12_prev(npoib))
        allocate (dqli21_prev(npoib))
        allocate (dqli22_prev(npoib))
        allocate (timstep_arr(neqset), tim_stack(neqset))

    end subroutine

    subroutine copy_kin_profs_to_yprev

        use grid_mod, only: npoi, nbaleqs
        use plasma_parameters, only: params

        implicit none

        integer :: ipoi, ieq, k

        if (debug_mode) write(*,*) 'Debug: yprev loop'
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                yprev(k) = params(ieq, ipoi)
            end do
        end do

    end subroutine

    subroutine initialize_antenna_factor

        use wave_code_data, only: antenna_factor
        implicit none

        antenna_factor_max = antenna_factor
        if (ramp_up_mode .eq. 4) then
            antenna_factor = 0d0
        else
            antenna_factor = 1.d-4
        end if

    end subroutine

    subroutine alloc_Br_Dqle_for_timeevol

        use grid_mod, only: T_tot_phi_e, T_tot_phi_i

        implicit none

        allocate(br_abs(Nstorage))
        allocate(br_formfactor(Nstorage))
        allocate(br_vac_res(Nstorage))
        allocate(br_abs_antenna_factor(Nstorage))
        allocate(br_abs_time(Nstorage))
        allocate(dqle22_res_time(Nstorage))
        allocate(dae22_res_time(Nstorage))
        allocate(bif_criterion(Nstorage))
        allocate(Ipar_time(Nstorage))
        allocate(T_tot_phi_e(Nstorage))
        allocate(T_tot_phi_i(Nstorage))

    end subroutine

    subroutine hold_prev_transp_coeffs

        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
                            dqli11, dqli12, dqli21, dqli22
        implicit none

        dqle11_prev = dqle11
        dqle12_prev = dqle12
        dqle21_prev = dqle21
        dqle22_prev = dqle22
        dqli11_prev = dqli11
        dqli12_prev = dqli12
        dqli21_prev = dqli21
        dqli22_prev = dqli22

    end subroutine

    subroutine allocate_timscal_and_params

        use grid_mod, only: npoi, npoic, nbaleqs, dummy
        use plasma_parameters, only: params_beg, params_num, params_denom, params_begbeg

        implicit none

        allocate(timscal(npoi), dummy(npoic))
        allocate(params_beg(nbaleqs, npoic), params_num(nbaleqs, npoic))
        allocate(params_denom(nbaleqs, npoic))
        allocate(params_begbeg(nbaleqs, npoic))

    end subroutine


    !> @brief subroutine write_br_dqle22_time_data. Writes radial magnetic field perturbation evaluated at the resonant
    !> surface, the antenna factor, the time and Dqle22 evaluated at the resonant surface for a given
    !> time step to the hdf5 file
    !> @param[in] i Integer of time step to which the data will be saved. Goes from 1:i.
    !> @param[in] br_abs_time Time value of the time evolution.
    !> @param[in] br_abs_antenna_factor Value of the antenna factor, i.e. the RMP coil current.
    !> @param[in] br_abs Absolute value of the radial magnetic field evaluated at the resonant surface in question.
    !> @param[in] dqle22_res_time Value of Dqle22 evaluated at the resonant surface during the time evolution.
    subroutine write_br_dqle22_time_data

        use control_mod
        use baseparam_mod
        use h5mod
        use KAMEL_hdf5_tools
        use wave_code_data, only: antenna_factor
        use resonances_mod, only: r_res

        implicit none

        if (debug_mode) write(*,*) "Debug: writing out br time evolution data"

        if (ihdf5IO .eq. 1) then
        !if (.false.) then
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

            h5overwrite = .true.

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_time"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_time(1:time_ind), &
                lbound(br_abs_time(1:time_ind)), ubound(br_abs_time(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_vac_res"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), abs(br_vac_res(1:time_ind)), &
                lbound(br_vac_res(1:time_ind)), ubound(br_vac_res(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_antenna_factor"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_antenna_factor(1:time_ind), &
                lbound(br_abs_antenna_factor(1:time_ind)), ubound(br_abs_antenna_factor(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_res"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs(1:time_ind), &
                lbound(br_abs(1:time_ind)), ubound(br_abs(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/dqle22_res_time"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), dqle22_res_time(1:time_ind), &
                lbound(dqle22_res_time(1:time_ind)), ubound(dqle22_res_time(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/dae22_res_time"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), dae22_res_time(1:time_ind), &
                lbound(dae22_res_time(1:time_ind)), ubound(dae22_res_time(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/bifurcation_criterion"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), bif_criterion(1:time_ind), &
                lbound(bif_criterion(1:time_ind)), ubound(bif_criterion(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_formfactor"
            CALL h5_add_complex_1(h5_id, trim(h5_currentgrp), br_formfactor(1:time_ind), &
                lbound(real(br_formfactor(1:time_ind))), ubound(real(br_formfactor(1:time_ind))))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/Ipar"
            CALL h5_add_complex_1(h5_id, trim(h5_currentgrp), Ipar_time(1:time_ind), &
                lbound(real(Ipar_time(1:time_ind))), ubound(real(Ipar_time(1:time_ind))))

            h5overwrite = .false.

            CALL h5_close(h5_id)
            CALL h5_deinit()

        else
            open (777, file='br_abs_res.dat', position='append')
            write (777, *) time_ind, time, antenna_factor, br_abs(time_ind)
            close (777)
        end if
    end subroutine ! write_br_dqle22_time_data


    subroutine check_linear_discr_pen_ratio

        implicit none

        if (time_ind .gt. 50 .and. .not. discr_reached) then
            ! calculate beta only once
            if (br_beta .eq. 0) then
                ! Calculate slope from the data until this time step. (Simple linear regression algorithm)
                br_beta = sum(br_abs(3:time_ind)*br_abs_time(3:time_ind))/sum(br_abs_time(1:time_ind)**2)
                write(*,*) 'br_beta = ', br_beta
            end if
            ! calculate the from the linear regression predicted value of Br_abs
            br_predicted = br_beta*br_abs_time(time_ind)
            write(*,*) "Delta = ", abs(br_abs(time_ind) - br_predicted)
            if (abs(br_abs(time_ind) - br_predicted) .gt. 0.1) then
                write(*,*) 'discrepancy to linearly predicted value of Br_abs_res > delta'
                if (modulo(time_ind, save_prof_time_step) .ne. 0) then
                    if (suppression_mode .eqv. .false.) then
                        CALL write_fields_currs_transp_coefs_to_h5
                    end if
                end if
                if (suppression_mode .eqv. .false.) then
                    CALL write_kin_prof_data_to_disk
                end if
                if (br_stopping) then

                    call write_reason_for_stop_to_h5("discrepancy to " //&
                        "linearly predicted value of Br_abs_res > delta")
                    CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                    stop "Finished time evolution: br_stopping"

                else
                    call write_br_discrepancy_reached_info
                    discr_reached = .true.
                end if
            end if
        end if

    end subroutine

    subroutine write_br_discrepancy_reached_info

        use h5mod

        implicit none

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
            '/info', 'discrepancy to linearly predicted value of Br_abs_res > delta')
        CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)// &
            '/discrep_time', (/time_ind*1.d0, time/), (/1/), (/2/))
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine


    !> @brief subroutine write_kin_prof_data_to_disk(time_ind). Writes the profile data to hdf5 files.
    !> Formerly, this data was written to fort.1xxx ascii files.
    !> This routine was added because of the change that only every
    !>  "save_prof_time_step"th timestep is written. If the program is to be stopped
    !> because a stopping criterion was met, the profiles should be written for that
    !> last time step. Because this occurs more than once, it is more convenient to
    !> summarize this in a subroutine.
    !> @author Markus Markl
    !> @date 12.03.2021
    !> @param[in] time_ind Current step of the time evolution. Used to name the fort.1000 group in which
    !> the data is written.
    subroutine write_kin_prof_data_to_disk

        use grid_mod
        use plasma_parameters
        use control_mod
        use baseparam_mod
        use h5mod
        use KAMEL_hdf5_tools
        use wave_code_data, only: Vth

        implicit none
        integer :: ipoi

        if (ihdf5IO .eq. 1) then
            if (debug_mode) print *, "Debug: Write kinetic profiles"
            do ipoi = 1, npoic
                sqg_bthet_overcavg(ipoi) = 0.5d0*(sqrt_g_times_B_theta_over_c(ipoi) &
                                                + sqrt_g_times_B_theta_over_c(ipoi + 1))
                Ercovavg(ipoi) = 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            ! h5_mode_groupname
            !h5_currentgrp = "/"//trim(h5_mode_groupname) &
            h5_currentgrp = trim(h5_mode_groupname) &
                            //"/KinProfiles"

            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

            write (h5_currentgrp, "(A,A,I4,A)") trim(h5_currentgrp), &
                "/", 1000 + time_ind, "/"

            if (debug_mode) write (*, *) "Debug: h5_currentgrp ", trim(h5_currentgrp)
            if (debug_mode) write (*, *) "Debug: defining KinProfiles/1000 group ", 1000 + time_ind

            CALL h5_obj_exists(h5_id, trim(h5_currentgrp), h5_exists_log)
            if (.not. h5_exists_log) then
                CALL h5_create_parent_groups(h5_id, trim(h5_currentgrp))
            end if

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"rc", &
                                real(rc), lbound(rc), ubound(rc))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"n", &
                                real(params(1, :)), lbound(params(1, :)), ubound(params(1, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Vz", &
                                real(params(2, :)), lbound(params(2, :)), ubound(params(2, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Te", &
                                real(params(3, :)/ev), lbound(params(3, :)), ubound(params(3, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Ti", &
                                real(params(4, :)/ev), lbound(params(4, :)), ubound(params(4, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Er", &
                                real(Ercovavg), lbound(Ercovavg), ubound(Ercovavg))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"sqg_btheta_overc", &
                                real(sqg_bthet_overcavg), lbound(sqg_bthet_overcavg), &
                                ubound(sqg_bthet_overcavg))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Vth", &
                                real(Vth), lbound(Vth), ubound(Vth))

            CALL h5_close(h5_id)
            CALL h5_deinit()
            if (debug_mode) write (*, *) "Debug: finished writing KinProfiles"

        else
            do ipoi = 1, npoic
                write (1000 + time_ind, *) rc(ipoi), params(1:2, ipoi) &
                    , params(3, ipoi)/ev &
                    , params(4, ipoi)/ev &
                    , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1)) &
                    , 0.5d0*(sqrt_g_times_B_theta_over_c(ipoi) + &
                            sqrt_g_times_B_theta_over_c(ipoi + 1))
            end do
            close (1000 + time_ind)
        end if

    end subroutine write_kin_prof_data_to_disk


    subroutine write_kin_profile_at_time_index

        implicit none

        if (debug_mode) write(*,*) "Debug: Write kinetic profiles at time index: ", time_ind
        if (modulo(time_ind, save_prof_time_step) .eq. 0) then
            call write_kin_prof_data_to_disk
        end if

    end subroutine


    subroutine interp_Br_Dql_at_resonance_timeevol

        use PolyLagrangeInterpolation
        use grid_mod, only: npoib, r_resonant, rb, dqle22, dae22
        use wave_code_data, only: antenna_factor, Br

        implicit none

        integer :: indResRadius, ind_begin_interp, ind_end_interp

        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call get_ind_Lagr_interp(indResRadius, ind_begin_interp, ind_end_interp)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(ind_begin_interp:ind_end_interp), coef)

        br_abs(time_ind) = sum(coef(0, :)*abs(Br(ind_begin_interp:ind_end_interp)))!*sqrt(antenna_factor)
        dqle22_res_time(time_ind) = sum(coef(0, :)*dqle22(ind_begin_interp:ind_end_interp))
        dae22_res_time(time_ind) = sum(coef(0, :)*dae22(ind_begin_interp:ind_end_interp))
        bif_criterion(time_ind) = dqle22_res_time(time_ind)/dae22_res_time(time_ind)

        if (bif_criterion(time_ind) .ge. 1.0d0) then
            write(*,*) "!!! Bifurcation criterion met !!!"
            write(*,*) "bif_criterion = ", bif_criterion(time_ind)
        end if

        ! save the time for the improved stopping criterion
        br_abs_time(time_ind) = time
        br_abs_antenna_factor(time_ind) = antenna_factor

    end subroutine


    subroutine message_Br_Dqle_values

        use wave_code_data, only: antenna_factor

        implicit none

        write(*,*) " "
        write(*,*) "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +"
        WRITE(*,'(A9,F10.4,A16,I4)') '  time = ', br_abs_time(time_ind), " s, time_ind = ", time_ind
        WRITE(*,'(A23,E10.5,A11,F6.2,A12)') '    Antenna_factor   = ', antenna_factor, " which are ", &
            antenna_factor/antenna_factor_max*100, "% of the max"
        WRITE(*,'(A23,F10.5,A2)') '    Br abs res * C_mn= ', br_abs(time_ind) * SQRT(antenna_factor), " G"
        WRITE(*,'(A23,F10.5,A2)') '    Br abs res       = ', br_abs(time_ind), " G"
        WRITE(*,'(A23,F10.3,A7)') '    Dqle22 res       = ', dqle22_res_time(time_ind), " cm^2/s"
        WRITE(*,'(A23,F10.5)')    '    bif crit         = ', bif_criterion(time_ind)
        write(*,*) '   Form factor      = ', abs(br_formfactor(time_ind))
        write(*,*) '   Ipar             = ', abs(Ipar_time(time_ind))
        write(*,*) "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +"
        write(*,*) " "

    end subroutine

    subroutine stop_if_time_step_too_small

        use h5mod

        implicit none

        character(*), parameter :: reason = 'timestep < stop time step'

        if (timstep .lt. stop_time_step .and. time .gt. 1.0d-3) then
            write(*,*) 'stop: timestep smaller than stop limit'
            if (suppression_mode .eqv. .false.) then
                call write_kin_prof_data_to_disk
            end if

            call write_reason_for_stop_to_h5(reason)
            stop
        end if

    end subroutine


    subroutine calc_params_num_and_denom

        use plasma_parameters, only: params_num, params_denom, params, params_beg

        implicit none

        params_num = (params - params_beg)**2
        params_denom = params**2 + params_beg**2

    end subroutine

    subroutine smooth_params_num_and_denom

        use grid_mod, only: npoi, nbaleqs, mwind, dummy
        use plasma_parameters, only: params_num, params_denom
        use QLBalance_kinds, only: dp

        implicit none

        integer :: ieq
        real(dp), allocatable :: row_buffer(:)

        allocate(row_buffer(npoi))
        do ieq = 1, nbaleqs
            row_buffer = params_num(ieq, :)
            call smooth_array_gauss(npoi, mwind, row_buffer, dummy)
            params_num(ieq, :) = dummy
            row_buffer = params_denom(ieq, :)
            call smooth_array_gauss(npoi, mwind, row_buffer, dummy)
            params_denom(ieq, :) = dummy
        end do
        deallocate(row_buffer)

    end subroutine

    subroutine determine_timscal

        use grid_mod, only: npoi, rc
        use plasma_parameters, only: params_num, params_denom
        use baseparam_mod, only: factolmax
        use recstep_mod, only: tol

        implicit none

        integer :: ipoi

        do ipoi = 1, npoi
            if (rc(ipoi) .lt. 0.95d0*rc(npoi)) then
                timscal(ipoi) = sum(sqrt(params_num(3:4, ipoi)/params_denom(3:4, ipoi)))
            else
                timscal(ipoi) = 1d-30 !0.d0
            end if
        end do

        if (debug_mode) then
            write(*,*) "maxval(timscal) = ", maxval(timscal)
            write(*,*) "tol*factolmax = ", tol*factolmax
        end if

    end subroutine

    subroutine rescale_time_step_array

        use grid_mod, only: npoi, nbaleqs
        use recstep_mod, only: tim_stack, timstep_arr, tol
        use QLbalance_diag, only: timscal_dql

        implicit none

        integer :: ipoi, ieq, k

        timscal = timscal + timscal_dql

        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                !timstep_arr(k)=timstep_arr(k)/timscal(ipoi)*tol
                timstep_arr(k) = timstep_arr(k)/max(timscal(ipoi), epsilon(1.d0))*tol
                ! steady state solution:
                !if (ieq .gt. 1 .and. r(ipoi) .gt. rsepar-0.5d0) then
                !    timstep_arr(k) = 0d0
                !end if
            end do
        end do

        timstep_arr = timstep_arr*timescale/(timstep_arr + timescale)
        if (scratch) then
            scratch = .false.
            tim_stack = timstep_arr
        end if
        timstep_arr = 2.d0*timstep_arr*tim_stack/(timstep_arr + tim_stack)

    end subroutine

    subroutine set_time_step

        use recstep_mod, only: timstep_arr, tol

        implicit none

        timstep = minval(timstep_arr)

        if (.not. set_constant_time_step) then
            ! limit time step from below:
            timstep = max(timstep, timstep_min)
            ! limit timestep from above:
            !if (ramp_up_mode .ne. 0) timstep = min(timstep,0.1)
            !timstep = min(timstep,0.005)
        else
            ! use for constant time step:
            timstep = constant_time_step
            write(*,*) "constant time step = ", timstep
        end if

        if (debug_mode) write(*,*) 'Debug: timstep', real(timstep), '   timescale', real(timescale), &
            'tolerance', real(tol)

    end subroutine

    subroutine reset_timstep_arr_w_timstep

        use recstep_mod, only: tim_stack, timstep_arr

        implicit none

        timstep_arr = 0.0d0
        timstep_arr = timstep
        tim_stack = timstep_arr

    end subroutine

    subroutine write_time_info

        implicit none

        if (ihdf5IO .eq. 1) then
            call write_time_info_to_h5
        else
            call write_time_info_to_txt
        end if

    end subroutine

    subroutine write_time_info_to_txt

        use QLbalance_diag, only: rate_dql, timscal_dql

        implicit none

        open (4321, file='timstep_evol.dat', position='append')
        write (4321, *) time_ind, timstep, timscal_dql, timscal(1), rate_dql, time
        close (4321)

    end subroutine

    subroutine write_time_info_to_h5

        use QLbalance_diag, only: rate_dql, timscal_dql

        implicit none

        h5_currentgrp = trim("/"//trim(h5_mode_groupname)//"/timstep_evol.dat")

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)

        if (.not. firstiterationdone) then
            if (diagnostics_output) then
                CALL h5_define_unlimited_matrix(h5_id, trim(h5_currentgrp), &
                                                H5T_NATIVE_DOUBLE, (/6, -1/), time_dataset_id)
            else
                CALL h5_define_unlimited_matrix(h5_id, trim(h5_currentgrp), &
                                                H5T_NATIVE_DOUBLE, (/3, -1/), time_dataset_id)
            end if
        end if
        ! this output works differently, because it is done for each
        ! time step and not all at once. For every time step a row
        ! is appended to the data. This is different to the case
        ! where the whole data is written at once. In the latter
        ! the data for one variable is saved in one row, i.e. there
        ! are length(variable) number of columns.
        if (diagnostics_output) then
            CALL h5_append_double_1(time_dataset_id, (/time_ind*1.d0, timstep, &
                                                        timscal_dql, timscal(1), rate_dql, time/), time_ind)
        else
            CALL h5_append_double_1(time_dataset_id, (/time_ind*1.d0, timstep, &
                                                        time/), time_ind)
        end if
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

    subroutine relax_plasma_parameters

        use grid_mod, only: npoi, nbaleqs
        use plasma_parameters, only: params
        use baseparam_mod, only: urelax

        implicit none

        integer :: ipoi, ieq, k

        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                params(ieq, ipoi) = yprev(k)*urelax + params(ieq, ipoi)*(1.d0 - urelax)
            end do
        end do

    end subroutine

    subroutine msg_time_info

        implicit none

        write(*,*) ' '
        write(*,*) 'Debug: i = ', int2(time_ind), 'time = ', real(time)
        write(*,*) ' '

    end subroutine

    subroutine set_first_iteration_true

        implicit none

        if (firstiterationdone .eqv. .false.) firstiterationdone = .true.

    end subroutine


    subroutine determine_Dql_diagnostic

        use grid_mod, only: dqle11, dqli11
        use QLbalance_diag

        implicit none

        timscal_dql = maxval(abs(dqle11_prev - dqle11))/maxval(dqle11_prev + dqle11)
        ind_dqle = maxloc(abs(dqle11_prev - dqle11))
        timscal_dqli = maxval(abs(dqli11_prev - dqli11))/maxval(dqli11_prev + dqli11)
        ind_dqli = maxloc(abs(dqli11_prev - dqli11))
        rate_dql = timscal_dql/timstep

    end subroutine

    subroutine create_group_structure_timeevol

        use control_mod
        use wave_code_data, only: m_vals, n_vals
        use resonances_mod, only: numres
        use h5mod

        implicit none

        if (debug_mode) write (*, *) "Debug: Creating group structure for TimeEvol"

        if (numres .eq. 1) then
            write (h5_mode_groupname, "(A,I0,A,I0)") "f_", m_vals(1), "_", n_vals(1)
        else
            write (h5_mode_groupname, "(A)") "multi_mode"
        end if

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)

        if (.not. suppression_mode) then
            if (debug_mode) write(*,*) "Debug: h5_mode_groupname ", trim(h5_mode_groupname)
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname) //'/')
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname)//"/KinProfiles/")

            call create_group_if_not_existent(trim(h5_mode_groupname)//"/LinearProfiles/")
            call create_group_if_not_existent("/init_params")
        else
            if (debug_mode) write (*,*) "Debug: h5_mode_groupname: ", trim(h5_mode_groupname)
            call create_group_if_not_existent(trim(h5_mode_groupname))
            call create_group_if_not_existent("/init_params")
        end if

        CALL h5_close(h5_id)
        CALL h5_deinit()

        if (debug_mode) write (*, *) "Debug: finished creating group structure for TimeEvol"
    end subroutine



    subroutine redoTimeStep

        use recstep_mod, only: timstep_arr
        use grid_mod, only: npoic, rc, Ercov
        use plasma_parameters, only: params, params_begbeg
        use baseparam_mod, only: eV, factolmax
        use restart_mod

        implicit none

        integer :: ipoi

        print *, 'redo step with old DQL'

        call hold_prev_transp_coeffs
        iunit_redo = 137

        open (iunit_redo, file='params_redostep.after')
        do ipoi = 1, npoic
            write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                , params(3, ipoi)/ev &
                , params(4, ipoi)/ev &
                , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
        end do
        close (iunit_redo)
        params = params_begbeg
        open (iunit_redo, file='params_redostep.before')
        do ipoi = 1, npoic
            write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                , params(3, ipoi)/ev &
                , params(4, ipoi)/ev &
                , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
        end do
        close (iunit_redo)

        timstep = timstep/factolmax
        timstep_arr = timstep

    end subroutine

end module
