module profile_input_m
    !> Profile input preprocessing module for KIM
    !> Handles coordinate detection, transformation, and Er calculation

    use KIM_kinds_m, only: dp
    use config_m, only: coord_type, input_profile_dir, equil_file, profile_location
    use setup_m, only: btor, R0
    use constants_m, only: e_charge, sol

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
        !> Er = (Ti/ei*ni)*dni/dr + (1/ei)*dTi/dr + (r*B0*Vz)/(c*q*R0)
        !> Uses k=0 (no poloidal rotation)
        implicit none

        logical, intent(out) :: er_exists
        character(256) :: er_filename

        er_filename = trim(profile_location) // '/Er.dat'
        inquire(file=trim(er_filename), exist=er_exists)

        if (er_exists) then
            write(*,*) 'Found Er.dat, using provided radial electric field'
            return
        end if

        write(*,*) 'WARNING: Er.dat not found'
        write(*,*) 'Calculating Er from force balance without V_pol (k=0)'
        call calculate_er_from_force_balance()

    end subroutine check_and_calculate_er

    subroutine validate_btor_rbig()
        !> Compare namelist btor/R0 with btor_rbig.dat, warn on mismatch >1%
        implicit none

        character(256) :: filename
        real(dp) :: file_btor, file_rbig
        real(dp) :: btor_diff, rbig_diff
        integer :: iunit, ios
        logical :: file_exists
        real(dp), parameter :: MISMATCH_THRESHOLD = 0.01_dp  ! 1%

        filename = trim(profile_location) // '/btor_rbig.dat'
        inquire(file=trim(filename), exist=file_exists)

        if (.not. file_exists) then
            write(*,*) 'Note: btor_rbig.dat not found, using namelist values only'
            return
        end if

        ! Read btor_rbig.dat
        open(newunit=iunit, file=trim(filename), status='old', action='read')
        read(iunit, *, iostat=ios) file_btor, file_rbig
        close(iunit)

        if (ios /= 0) then
            write(*,*) 'WARNING: Could not read btor_rbig.dat'
            return
        end if

        ! Compare with namelist values
        if (abs(btor) > 1.0e-10_dp) then
            btor_diff = abs(file_btor - btor) / abs(btor)
            if (btor_diff > MISMATCH_THRESHOLD) then
                write(*,*) 'WARNING: btor mismatch > 1%'
                write(*,*) '  Namelist btor = ', btor
                write(*,*) '  File btor     = ', file_btor
                write(*,*) '  Difference    = ', btor_diff * 100.0_dp, '%'
            end if
        end if

        if (abs(R0) > 1.0e-10_dp) then
            rbig_diff = abs(file_rbig - R0) / abs(R0)
            if (rbig_diff > MISMATCH_THRESHOLD) then
                write(*,*) 'WARNING: R0 mismatch > 1%'
                write(*,*) '  Namelist R0 = ', R0
                write(*,*) '  File rbig   = ', file_rbig
                write(*,*) '  Difference  = ', rbig_diff * 100.0_dp, '%'
            end if
        end if

    end subroutine validate_btor_rbig

    subroutine calculate_er_from_force_balance()
        !> Calculate Er from radial force balance (k=0, no V_pol)
        !> Er = (Ti/ei*ni)*dni/dr + (1/ei)*dTi/dr + (r*B0*Vz)/(c*q*R0)
        !> Units: CGS (statV/cm)
        implicit none

        real(dp), allocatable :: r(:), n(:), Ti(:), Vz(:), q(:), Er(:)
        real(dp), allocatable :: dn_dr(:), dTi_dr(:)
        integer :: npts, i, iunit
        character(256) :: filename
        logical :: vz_exists

        ! Read required profiles
        call read_profile_data(trim(profile_location)//'/n.dat', r, n, npts)
        call read_profile_data(trim(profile_location)//'/Ti.dat', r, Ti, npts)
        call read_profile_data(trim(profile_location)//'/q.dat', r, q, npts)

        ! Check for Vz profile
        filename = trim(profile_location) // '/Vz.dat'
        inquire(file=trim(filename), exist=vz_exists)
        if (vz_exists) then
            call read_profile_data(filename, r, Vz, npts)
        else
            write(*,*) 'Note: Vz.dat not found, assuming Vz=0'
            allocate(Vz(npts))
            Vz = 0.0_dp
        end if

        ! Calculate derivatives
        allocate(dn_dr(npts), dTi_dr(npts), Er(npts))
        call calculate_derivative(r, n, dn_dr, npts)
        call calculate_derivative(r, Ti, dTi_dr, npts)

        ! Calculate Er from force balance (CGS units)
        ! Ti is in eV, convert to erg: Ti_erg = Ti_eV * e_charge
        ! Er in statV/cm
        do i = 1, npts
            if (abs(n(i)) > 1.0e-20_dp) then
                ! Term 1: (Ti/ei*ni)*dni/dr
                Er(i) = (Ti(i) * e_charge / (e_charge * n(i))) * dn_dr(i)
                ! Term 2: (1/ei)*dTi/dr
                Er(i) = Er(i) + (1.0_dp / e_charge) * dTi_dr(i) * e_charge
                ! Term 3: (r*B0*Vz)/(c*q*R0)
                if (abs(q(i)) > 1.0e-10_dp .and. abs(R0) > 1.0e-10_dp) then
                    Er(i) = Er(i) + (r(i) * btor * Vz(i)) / (sol * q(i) * R0)
                end if
            else
                Er(i) = 0.0_dp
            end if
        end do

        ! Write to Er_no_Vpol.dat
        filename = trim(profile_location) // '/Er_no_Vpol.dat'
        open(newunit=iunit, file=trim(filename), status='replace', action='write')
        do i = 1, npts
            write(iunit, '(2E20.12)') r(i), Er(i)
        end do
        close(iunit)

        write(*,*) 'Er calculated from force balance, written to Er_no_Vpol.dat'

        ! Also create Er.dat for read_profiles to find
        filename = trim(profile_location) // '/Er.dat'
        open(newunit=iunit, file=trim(filename), status='replace', action='write')
        do i = 1, npts
            write(iunit, '(2E20.12)') r(i), Er(i)
        end do
        close(iunit)

        deallocate(r, n, Ti, Vz, q, Er, dn_dr, dTi_dr)

    end subroutine calculate_er_from_force_balance

    subroutine read_profile_data(filename, r, data, npts)
        !> Read two-column profile data file
        character(*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: r(:), data(:)
        integer, intent(out) :: npts

        integer :: iunit, ios
        real(dp) :: r_val, data_val
        real(dp), allocatable :: r_temp(:), data_temp(:)
        integer, parameter :: MAX_PTS = 10000

        allocate(r_temp(MAX_PTS), data_temp(MAX_PTS))
        npts = 0

        open(newunit=iunit, file=trim(filename), status='old', action='read')
        do
            read(iunit, *, iostat=ios) r_val, data_val
            if (ios /= 0) exit
            npts = npts + 1
            if (npts > MAX_PTS) then
                write(*,*) 'ERROR: Too many points in ', trim(filename)
                stop 1
            end if
            r_temp(npts) = r_val
            data_temp(npts) = data_val
        end do
        close(iunit)

        allocate(r(npts), data(npts))
        r = r_temp(1:npts)
        data = data_temp(1:npts)
        deallocate(r_temp, data_temp)

    end subroutine read_profile_data

    subroutine calculate_derivative(x, y, dydx, n)
        !> Calculate derivative using central differences
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: dydx(:)
        integer, intent(in) :: n
        integer :: i

        if (n < 2) then
            dydx = 0.0_dp
            return
        end if

        ! Forward difference at first point
        dydx(1) = (y(2) - y(1)) / (x(2) - x(1))

        ! Central differences for interior
        do i = 2, n-1
            dydx(i) = (y(i+1) - y(i-1)) / (x(i+1) - x(i-1))
        end do

        ! Backward difference at last point
        dydx(n) = (y(n) - y(n-1)) / (x(n) - x(n-1))

    end subroutine calculate_derivative

end module profile_input_m
