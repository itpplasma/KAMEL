  module time_evolution

    implicit none

    integer :: Nstorage
    double precision :: tmax_factor!, antenna_factor
    logical :: flag_run_time_evolution !Added by Philipp Ulbl 12.05.2020
    double precision :: stop_time_step !Added by Philipp Ulbl 13.05.2020
    double precision :: timstep_min
    logical :: br_stopping ! trigger Br stopping criterion
    integer :: ramp_up_mode !> control ramp up mode of the RMP coil current amplitude
    DOUBLE PRECISION :: t_max_ramp_up = 1e-2 !> 10ms ramp up until antenna_factor_max is reached
    integer :: save_prof_time_step ! added by Markus Markl 11.03.2021
    double precision :: timstep
    
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: yprev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle11_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle12_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle21_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle22_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli11_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli12_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli21_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli22_prev

    double precision :: antenna_factor_max
    DOUBLE PRECISION :: antenna_max_stopping

    !needed for interpolation of br abs and stopping criterion
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: br_abs
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: br_abs_time
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: br_abs_antenna_factor
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle22_res_time

    contains

    subroutine initAntennaFactor
        !For time evolution mode use antenna_factor as maximum
        !and start with a very small value and ramp this up

        use wave_code_data, only: antenna_factor
        implicit none

        antenna_factor_max = antenna_factor
        if (ramp_up_mode .eq. 4) then 
            antenna_factor = 0d0
        else
            antenna_factor = 1.d-4
        end if

    end subroutine

    subroutine allocateBrAndDqleForTimeEvolution

        implicit none

        allocate(br_abs(Nstorage))
        allocate(br_abs_antenna_factor(Nstorage))
        allocate(br_abs_time(Nstorage))
        allocate(dqle22_res_time(Nstorage))

    end subroutine

    subroutine savePrevTranspCoefficients

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

    subroutine rescaleTranspCoefficientsByAntennaFac

        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
                            dqli11, dqli12, dqli21, dqli22
        use wave_code_data, only: antenna_factor
        
        implicit none

        dqle11 = dqle11*antenna_factor
        dqle12 = dqle12*antenna_factor
        dqle21 = dqle21*antenna_factor
        dqle22 = dqle22*antenna_factor
        dqli11 = dqli11*antenna_factor
        dqli12 = dqli12*antenna_factor
        dqli21 = dqli21*antenna_factor
        dqli22 = dqli22*antenna_factor

    end subroutine

  end module