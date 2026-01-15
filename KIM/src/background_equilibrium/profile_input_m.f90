module profile_input_m
    !> Profile input preprocessing module for KIM
    !> Handles coordinate detection, transformation, and Er calculation

    use KIM_kinds_m, only: dp
    use config_m, only: coord_type, input_profile_dir, equil_file, profile_location
    use setup_m, only: btor, R0

    implicit none
    private

    public :: prepare_profiles

    ! Coordinate detection threshold
    real(dp), parameter :: COORD_THRESHOLD = 2.0_dp

contains

    subroutine prepare_profiles()
        !> Main entry point - called before read_profiles() in KIM_init
        implicit none

        character(20) :: detected_coord_type
        logical :: er_exists

        ! 1. Detect coordinate type if auto
        if (trim(coord_type) == 'auto') then
            call detect_coordinate_type(detected_coord_type)
        else
            detected_coord_type = coord_type
        end if

        ! 2. Process based on coordinate type
        if (trim(detected_coord_type) == 'sqrt_psiN') then
            call run_preprocessing()
        end if

        ! 3. Check for Er.dat and calculate if missing
        call check_and_calculate_er(er_exists)

        ! 4. Validate btor/R0 against btor_rbig.dat
        call validate_btor_rbig()

    end subroutine prepare_profiles

    subroutine detect_coordinate_type(detected_type)
        !> Auto-detect coordinate type from input profile
        !> Reads first column of n profile and checks max value
        character(20), intent(out) :: detected_type

        character(256) :: filename
        real(dp) :: r_val, dummy
        real(dp) :: max_r
        integer :: ios, iunit
        logical :: file_exists

        ! Try n_of_psiN.dat first, then n.dat
        filename = trim(input_profile_dir) // '/n_of_psiN.dat'
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            filename = trim(input_profile_dir) // '/n.dat'
            inquire(file=trim(filename), exist=file_exists)
        end if

        if (.not. file_exists) then
            write(*,*) 'ERROR: No density profile found for coordinate detection'
            write(*,*) '  Searched: ', trim(input_profile_dir), '/n_of_psiN.dat'
            write(*,*) '       and: ', trim(input_profile_dir), '/n.dat'
            stop 1
        end if

        ! Read first column and find max
        max_r = 0.0_dp
        open(newunit=iunit, file=trim(filename), status='old', action='read')
        do
            read(iunit, *, iostat=ios) r_val, dummy
            if (ios /= 0) exit
            if (r_val > max_r) max_r = r_val
        end do
        close(iunit)

        ! Determine coordinate type
        if (max_r > COORD_THRESHOLD) then
            detected_type = 'r_eff'
            write(*,*) 'Auto-detected coordinate type: r_eff (max_r = ', max_r, ' cm)'
        else
            detected_type = 'sqrt_psiN'
            write(*,*) 'Auto-detected coordinate type: sqrt_psiN (max_r = ', max_r, ')'
        end if

    end subroutine detect_coordinate_type

    subroutine run_preprocessing()
        !> Run profile_preprocessor for sqrt_psiN -> r_eff transformation
        ! Placeholder - will implement
    end subroutine run_preprocessing

    subroutine check_and_calculate_er(er_exists)
        !> Check for Er.dat, calculate from force balance if missing
        logical, intent(out) :: er_exists

        er_exists = .true.  ! Placeholder - will implement
    end subroutine check_and_calculate_er

    subroutine validate_btor_rbig()
        !> Compare namelist btor/R0 with btor_rbig.dat, warn on mismatch
        ! Placeholder - will implement
    end subroutine validate_btor_rbig

end module profile_input_m
