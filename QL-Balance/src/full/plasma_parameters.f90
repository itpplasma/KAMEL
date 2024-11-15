module plasma_parameters

    use control_mod

    implicit none

    double precision, dimension(:,:), allocatable :: params,dot_params
    double precision, dimension(:,:), allocatable :: ddr_params, ddr_params_nl
    double precision, dimension(:, :), allocatable :: params_beg, params_begbeg
    double precision, dimension(:, :), allocatable :: params_num, params_denom
    double precision, dimension(:,:), allocatable :: params_b
    double precision, dimension(:,:), allocatable :: init_params
    double precision, dimension(:,:), allocatable :: params_lin,params_b_lin
    double precision, dimension(:),   allocatable :: qsafb,qsaf

    double precision, dimension(:), allocatable :: hold_n, hold_Te, hold_Ti, hold_Vz, hold_dphi0! variables to hold the initial bg profiles

    contains


    subroutine alloc_hold_parameters
        
        use grid_mod, only: npoib
        use wave_code_data, only: idPhi0

        implicit none

        allocate(hold_n(npoib))
        allocate(hold_Vz(npoib))
        allocate(hold_Te(npoib))
        allocate(hold_Ti(npoib))
        allocate(hold_dphi0(npoib))
        hold_n = params(1, :)
        hold_Vz = params(2, :)
        hold_Te = params(3, :)
        hold_Ti = params(4, :)
        hold_dphi0 = idPhi0

    end subroutine

    subroutine limit_temps_from_below

        use grid_mod, only: npoic
        use wave_code_data, only: r
        use baseparam_mod, only: ev, rsepar

        implicit none

        integer :: ipoi

        do ipoi = 1, npoic
            if (.true.) then
                params(3, ipoi) = max(params(3, ipoi), temperature_limit*ev)
                params(4, ipoi) = max(params(4, ipoi), temperature_limit*ev)
            else 
                !> Quick fix of steady state solution. Keep boundary inside the separatrix.
                if (r(ipoi) > rsepar-0.5d0) then
                    params(3, ipoi) = hold_Te(ipoi)
                    params(4, ipoi) = hold_Ti(ipoi)
                end if
            end if
        end do


    end subroutine limit_temps_from_below



    subroutine init_background_profiles

        use grid_mod, only: npoic
        use wave_code_data, only: q, n, Vz, Te, Ti
        use baseparam_mod, only: ev, rtor

        implicit none    
        integer :: ipoi                    !initial background profiles:

        do ipoi = 1, npoic
            !safety factor:
            qsaf(ipoi) = 0.5*(q(ipoi) + q(ipoi + 1))
            !electron density :
            params(1, ipoi) = 0.5*(n(ipoi) + n(ipoi + 1))
            !toroidal rotation frequency :
            params(2, ipoi) = 0.5*(Vz(ipoi) + Vz(ipoi + 1))/rtor
            !electron temeperature :
            params(3, ipoi) = 0.5*(Te(ipoi) + Te(ipoi + 1))*ev
            !ion temeperature :
            params(4, ipoi) = 0.5*(Ti(ipoi) + Ti(ipoi + 1))*ev
        end do

    end subroutine

    !> @brief subroutine write_init_profiles. Write initial profiles to hdf5 or ascii.
    !> @author Markus Markl
    !> @date 05.10.2022
    subroutine write_initial_parameters

        use control_mod, only: debug_mode, ihdf5IO
        use wave_code_data, only: r
        use baseparam_mod, only: ev
        use h5mod

        implicit none

        if (debug_mode) write(*,*) "Debug: writing initial background profiles"
        if (ihdf5IO .eq. 1) then
            CALL h5_init()
            ! open hdf5 file
            CALL h5_open_rw(path2out, h5_id)
            CALL h5_obj_exists(h5_id, "/init_params/n", h5_exists_log)
            if (.not. h5_exists_log) then
                CALL h5_add_double_1(h5_id, "/init_params/n", &
                    params(1, :), lbound(params(1, :)), ubound(params(1, :)))
                CALL h5_add_double_1(h5_id, "/init_params/Vz", &
                    params(2, :), lbound(params(2, :)), ubound(params(2, :)))
                CALL h5_add_double_1(h5_id, "/init_params/Te", &
                    params(3, :)/ev, lbound(params(3, :)), ubound(params(3, :)))
                CALL h5_add_double_1(h5_id, "/init_params/Ti", &
                    params(4, :)/ev, lbound(params(4, :)), ubound(params(4, :)))
                CALL h5_add_double_1(h5_id, "/init_params/qsaf", &
                    qsaf(:), lbound(qsaf(:)), ubound(qsaf(:)))
                CALL h5_add_double_1(h5_id, "/init_params/r", &
                    r, lbound(r), ubound(r))
            else
                if (debug_mode) write(*,*) "Debug: they are already there -> skiping"
            end if

            CALL h5_close(h5_id)
            CALL h5_deinit()
            if (debug_mode) write(*,*) "Debug: finished writing initial background profiles"
                                    !stop ! for test purposes

        else
            open (123, form='unformatted', file='init_params.dat')
            write (123) params
            close (123)
        end if

    end subroutine



end module