module kim_wave_code_adapter_m
    !! Adapter module that wraps KIM library calls to populate
    !! wave_code_data module arrays, matching the contract expected
    !! by get_dql.f90 and diff_coeffs.f90.

    implicit none
    private

    public :: kim_initialize
    public :: kim_run_for_all_modes
    public :: kim_update_profiles
    public :: kim_get_wave_fields
    public :: kim_get_wave_vectors
    public :: kim_get_background_magnetic_fields
    public :: kim_get_collision_frequencies

contains

    subroutine kim_initialize(nrad, r_grid)
        !! Initialize KIM backend: read config, profiles, grids.
        !! Populate wave_code_data module arrays with background quantities.
        integer, intent(in) :: nrad
        real(8), intent(in) :: r_grid(nrad)
        ! TODO: implement in Task 4
    end subroutine

    subroutine kim_run_for_all_modes()
        !! Run KIM solver for all (m,n) modes.
        ! TODO: implement in Task 5
    end subroutine

    subroutine kim_update_profiles()
        !! Update KIM's internal profiles from QL-Balance's evolved parameters.
        ! TODO: implement in Task 5
    end subroutine

    subroutine kim_get_wave_fields(i_mn)
        !! Extract wave fields from KIM's EBdat into wave_code_data arrays.
        integer, intent(in) :: i_mn
        ! TODO: implement in Task 6
    end subroutine

    subroutine kim_get_wave_vectors()
        !! Extract wave vectors from KIM's plasma into wave_code_data arrays.
        ! TODO: implement in Task 6
    end subroutine

    subroutine kim_get_background_magnetic_fields()
        !! Extract background B from KIM's equilibrium into wave_code_data arrays.
        ! TODO: implement in Task 6
    end subroutine

    subroutine kim_get_collision_frequencies()
        !! Extract collision frequencies from KIM's plasma into wave_code_data arrays.
        ! TODO: implement in Task 6
    end subroutine

end module kim_wave_code_adapter_m
