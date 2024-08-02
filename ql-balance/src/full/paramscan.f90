
module paramscan_mod
    integer :: ifac_n, ifac_Te, ifac_Ti, ifac_vz ! counter for the do loops
    integer :: numoffac                          ! total number of factors
    double precision, dimension(:), allocatable :: fac_n, fac_Te, &
        fac_Ti, fac_vz
    character(len=1024) :: parscan_str
    double precision :: viscosity_factor
    DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Er_res
    DOUBLE PRECISION, DIMENSION(:, :, :, :), ALLOCATABLE :: br_abs_res_parscan
    DOUBLE PRECISION, DIMENSION(:, :, :, :), ALLOCATABLE :: dqle22_res


    double precision, dimension(:), allocatable :: hold_n, hold_Te, hold_Ti, hold_Vz, hold_dphi0! variables to hold the initial bg profiles

    contains

    !> @brief subroutine initialize_parameter_scan_vars. Read factors for parameter scan and allocate variables. Still needed if no parameter scan is done.
    !> @author Markus Markl
    !> @date 05.10.2022
    subroutine initialize_parameter_scan_vars

        use control_mod, only: paramscan

        implicit none

        if (paramscan) then
            write(*,*) "Parameter scan: fetch factors for parameter scan"
            CALL getfactors
            write(*,*) "Got out of getfactors"
            if (size(fac_vz) .ne. 1) allocate(Er_res(size(fac_n), size(fac_Te), size(fac_Ti), size(fac_vz)))
            write(*,*) "allocated Er_res"
        else
            allocate(fac_n(1))
            allocate(fac_Ti(1))
            allocate(fac_Te(1))
            allocate(fac_vz(1))
            fac_n = (/1.d0/)
            fac_Ti = (/1.d0/)
            fac_Te = (/1.d0/)
            fac_vz = (/1.d0/)
        end if

        allocate(dqle22_res(size(fac_n), size(fac_Te), size(fac_Ti), size(fac_vz)))
        allocate(br_abs_res_parscan(size(fac_n), size(fac_Te), size(fac_Ti), size(fac_vz)))
        write(*,*) "Finished initialize parameter scan vars"

    end subroutine ! initialize_parameter_scan_vars


    subroutine getfactors

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

    subroutine alloc_hold_parameters
        
        use grid_mod, only: npoib, params
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

    !> @brief subroutine rescale_profiles. Rescales kinetic profiles (n,Vz,Te,Ti).
    !> @author Markus Markl
    !> @date 05.10.2022
    subroutine rescale_profiles

        use grid_mod, only: params
        use control_mod, only: debug_mode
        use h5mod
        use wave_code_data, only: idPhi0

        implicit none

        double precision, dimension(:), allocatable :: ErVzfac ! factor to rescale Er
        integer :: lowerBound, upperBound

        if (debug_mode) write(*,*) "Debug: coming into rescaling profiles"

        params(1, :) = hold_n * fac_n(ifac_n)
        params(2, :) = hold_vz * fac_vz(ifac_vz)
        params(3, :) = hold_Te * fac_Te(ifac_Te)
        params(4, :) = hold_Ti * fac_Ti(ifac_Ti)

        if (fac_vz(ifac_vz) .ne. 1.d0) then
            if (debug_mode) write(*,*) "Debug: fac_vz not equal 1. need to rescale Er as well"
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)
            CALL h5_get_bounds_1(h5_id, '/factors/ErVzfac', lowerBound, upperBound)
            allocate(ErVzfac(upperBound))
            CALL h5_get_double_1(h5_id, '/factors/ErVzfac', ErVzfac)
            CALL h5_close(h5_id)
            CALL h5_deinit()
            write (*, *) "rescale Er"
            idPhi0 = hold_dphi0 + ErVzfac*params(2, :)*(fac_vz(ifac_vz) - 1.d0)
            deallocate (ErVzfac)
        end if

        write(*,*) "Parameter scan, current factors: "
        write(*,*) "fac_n = ", fac_n(ifac_n), "   ", ifac_n, " of ", size(fac_n)
        write(*,*) "fac_Ti = ", fac_Ti(ifac_Ti), "   ", ifac_Ti, " of ", size(fac_Ti)
        write(*,*) "fac_Te = ", fac_Te(ifac_Te), "   ", ifac_Te, " of ", size(fac_Te)
        write(*,*) "fac_vz = ", fac_vz(ifac_vz), "   ", ifac_vz, " of ", size(fac_vz)

        if (debug_mode) write(*,*) "Debug: going out of rescaling profiles"

    end subroutine !rescale_profiles


    subroutine interpBrAndDqlAtResonanceParamScan

    end subroutine interpBrAndDqlAtResonanceParamScan


    subroutine run_parameter_scan


    end subroutine

end module