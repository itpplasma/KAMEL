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

    contains

    subroutine limitTemperaturesFromBelow

        use grid_mod, only: npoic
        use wave_code_data, only: r
        use baseparam_mod, only: ev, rsepar
        use paramscan_mod, only: hold_Te, hold_Ti

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


    end subroutine limitTemperaturesFromBelow



end module