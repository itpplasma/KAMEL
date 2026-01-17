module profile_input_m
    !> Profile input preprocessing module for KIM
    !> Handles coordinate detection, transformation, and Er calculation

    use KIM_kinds_m, only: dp
    use config_m, only: coord_type, input_profile_dir, equil_file, geqdsk_file, profile_location, &
                        n_input_file, Te_input_file, Ti_input_file, Vz_input_file, &
                        n_file, Te_file, Ti_file, Vz_file, Er_file, q_file
    use setup_m, only: btor, R0
    use constants_m, only: e_charge, sol, ev
    use grid_m, only: r_min, r_plas
    use profile_preprocessor_m, only: profile_preprocessor_t

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
        else if (trim(detected_coord_type) == 'r_eff') then
            ! Profiles already in r_eff coordinates - copy to profile_location if needed
            call copy_profiles_if_needed()
        end if

        ! 3. Check for Er.dat and calculate if missing
        call check_and_calculate_er(er_exists)

        ! 4. Validate btor/R0 against btor_rbig.dat
        call validate_btor_rbig()

        ! 5. Validate r_min/r_plas are within profile range
        call validate_radial_range()

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

        ! Try input file first, then output file
        filename = trim(input_profile_dir) // '/' // trim(n_input_file)
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            filename = trim(input_profile_dir) // '/' // trim(n_file)
            inquire(file=trim(filename), exist=file_exists)
        end if

        if (.not. file_exists) then
            write(*,*) 'ERROR: No density profile found for coordinate detection'
            write(*,*) '  Searched: ', trim(input_profile_dir), '/', trim(n_input_file)
            write(*,*) '       and: ', trim(input_profile_dir), '/', trim(n_file)
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
        implicit none

        type(profile_preprocessor_t) :: preprocessor
        character(256) :: equil_path
        logical :: equil_exists, geqdsk_exists
        character(256) :: temp_nml_file
        integer :: iunit

        ! Determine equilibrium file path
        if (len_trim(equil_file) > 0) then
            equil_path = equil_file
        else
            equil_path = trim(input_profile_dir) // '/equil_r_q_psi.dat'
        end if

        inquire(file=trim(equil_path), exist=equil_exists)
        if (.not. equil_exists) then
            ! Check if geqdsk_file is provided for equilibrium computation
            if (len_trim(geqdsk_file) > 0) then
                inquire(file=trim(geqdsk_file), exist=geqdsk_exists)
                if (.not. geqdsk_exists) then
                    write(*,*) 'ERROR: GEQDSK file not found: ', trim(geqdsk_file)
                    stop 1
                end if
                write(*,*) 'Equilibrium file not found, will compute from GEQDSK'
                write(*,*) '  GEQDSK file: ', trim(geqdsk_file)
                ! Write field_divB0.inp for libneo to read geqdsk
                call write_field_divB0_inp(geqdsk_file)
                ! Pass empty equil_path so preprocessor computes equilibrium
                equil_path = ''
            else
                write(*,*) 'ERROR: Equilibrium file not found: ', trim(equil_path)
                write(*,*) '  and no geqdsk_file specified for computation'
                write(*,*) 'Cannot perform sqrt_psiN -> r_eff transformation'
                stop 1
            end if
        end if

        write(*,*) 'Running profile preprocessing (sqrt_psiN -> r_eff)'
        write(*,*) '  Input directory: ', trim(input_profile_dir)
        write(*,*) '  Equilibrium: ', trim(equil_path)
        write(*,*) '  Output directory: ', trim(profile_location)

        ! Create output directory if it doesn't exist
        call execute_command_line('mkdir -p ' // trim(profile_location))

        ! Create temporary namelist file for the preprocessor
        temp_nml_file = trim(profile_location) // '/.profile_preproc_temp.nml'
        call write_preprocessor_namelist(temp_nml_file, equil_path)

        ! Initialize and run preprocessor
        call preprocessor%init(trim(temp_nml_file))
        call preprocessor%process()
        call preprocessor%write_output()
        call preprocessor%cleanup()

        ! Remove temporary namelist file
        open(newunit=iunit, file=trim(temp_nml_file), status='old')
        close(iunit, status='delete')

        write(*,*) 'Profile preprocessing complete'

    end subroutine run_preprocessing

    subroutine write_preprocessor_namelist(filename, equil_path)
        !> Write temporary namelist file for profile_preprocessor
        character(*), intent(in) :: filename
        character(*), intent(in) :: equil_path

        integer :: iunit

        open(newunit=iunit, file=trim(filename), status='replace', action='write')
        write(iunit, '(A)') '&profile_preprocessor'
        write(iunit, '(A)') '  equil_file = "' // trim(equil_path) // '"'
        write(iunit, '(A)') '  input_dir = "' // trim(input_profile_dir) // '"'
        write(iunit, '(A)') '  output_dir = "' // trim(profile_location) // '"'
        write(iunit, '(A)') '  coord_type = 1'  ! COORD_SQRT_PSIN
        write(iunit, '(A)') '  n_input_file = "' // trim(n_input_file) // '"'
        write(iunit, '(A)') '  Te_input_file = "' // trim(Te_input_file) // '"'
        write(iunit, '(A)') '  Ti_input_file = "' // trim(Ti_input_file) // '"'
        write(iunit, '(A)') '  Vz_input_file = "' // trim(Vz_input_file) // '"'
        write(iunit, '(A)') '/'
        close(iunit)

    end subroutine write_preprocessor_namelist

    subroutine write_field_divB0_inp(gfile_path)
        !> Write field_divB0.inp for libneo equilibrium computation
        !> Uses equilibrium-only mode (ipert=0, iequil=1)
        character(*), intent(in) :: gfile_path

        integer :: iunit
        character(512) :: code_path, convexwall_path

        ! Get KAMEL path from $CODE environment variable
        call get_environment_variable('CODE', code_path)
        if (len_trim(code_path) == 0) then
            write(*,*) 'ERROR: Environment variable $CODE is not set'
            write(*,*) '  Required to locate convexwall file for equilibrium computation'
            stop 1
        end if
        ! Remove trailing slash if present
        if (code_path(len_trim(code_path):len_trim(code_path)) == '/') then
            code_path = code_path(1:len_trim(code_path)-1)
        end if
        convexwall_path = trim(code_path) // '/KAMEL/common/equil/convexwall/convexwall.asdex'

        ! Warn user about ASDEX Upgrade specific convexwall
        write(*,*) 'WARNING: Using ASDEX Upgrade convexwall file for equilibrium computation'
        write(*,*) '  ', trim(convexwall_path)

        open(newunit=iunit, file='field_divB0.inp', status='replace', action='write')
        write(iunit, '(A)') '0                                 ipert        ! 0=eq only'
        write(iunit, '(A)') '1                                 iequil       ! 1=with equil.'
        write(iunit, '(A)') '1.00                              ampl         ! amplitude'
        write(iunit, '(A)') '72                                ntor         ! toroidal harmonics'
        write(iunit, '(A)') '0.99                              cutoff       ! inner cutoff'
        write(iunit, '(A)') '4                                 icftype      ! coil file type'
        write(iunit, '(A)') "'" // trim(gfile_path) // "'"
        write(iunit, '(A)') "''"
        write(iunit, '(A)') "'" // trim(convexwall_path) // "'"
        write(iunit, '(A)') "''"
        write(iunit, '(A)') '0                                 nwindow_r'
        write(iunit, '(A)') '0                                 nwindow_z'
        write(iunit, '(A)') '1                                 ieqfile      ! 1=EFIT format'
        close(iunit)

        write(*,*) 'Created field_divB0.inp for equilibrium computation'

    end subroutine write_field_divB0_inp

    subroutine check_and_calculate_er(er_exists)
        !> Check for Er.dat, calculate from force balance if missing
        !> Er = (Ti/ei*ni)*dni/dr + (1/ei)*dTi/dr + (r*B0*Vz)/(c*q*R0)
        !> Uses k=0 (no poloidal rotation)
        implicit none

        logical, intent(out) :: er_exists
        character(256) :: er_filename

        er_filename = trim(profile_location) // '/' // trim(Er_file)
        inquire(file=trim(er_filename), exist=er_exists)

        if (er_exists) then
            write(*,*) 'Found ', trim(Er_file), ', using provided radial electric field'
            return
        end if

        write(*,*) 'WARNING: ', trim(Er_file), ' not found'
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
        call read_profile_data(trim(profile_location)//'/'//trim(n_file), r, n, npts)
        call read_profile_data(trim(profile_location)//'/'//trim(Ti_file), r, Ti, npts)
        call read_profile_data(trim(profile_location)//'/'//trim(q_file), r, q, npts)

        ! Check for Vz profile
        filename = trim(profile_location) // '/' // trim(Vz_file)
        inquire(file=trim(filename), exist=vz_exists)
        if (vz_exists) then
            call read_profile_data(filename, r, Vz, npts)
        else
            write(*,*) 'Note: ', trim(Vz_file), ' not found, assuming Vz=0'
            allocate(Vz(npts))
            Vz = 0.0_dp
        end if

        ! Calculate derivatives
        allocate(dn_dr(npts), dTi_dr(npts), Er(npts))
        call calculate_derivative(r, n, dn_dr, npts)
        call calculate_derivative(r, Ti, dTi_dr, npts)

        ! Calculate Er from radial force balance (CGS units)
        ! From momentum balance: E_r = (T_i/e_i*n_i)*dn_i/dr + ((1-k)/e_i)*dT_i/dr + r*B0*Vtor/(c*q*R0)
        ! With k=0 (no poloidal rotation):
        ! E_r = (T_i*ev)/(e*n)*dn/dr + (ev/e)*dT/dr + r*B0*Vz/(c*q*R0)
        ! Ti is in eV, convert to erg using ev = 1.6022e-12 erg/eV
        ! e_charge = 4.803e-10 statcoulomb (CGS)
        ! Er in statV/cm
        do i = 1, npts
            if (abs(n(i)) > 1.0e-20_dp) then
                ! Term 1: (Ti*ev)/(e*n) * dn/dr (density gradient contribution)
                Er(i) = (Ti(i) * ev / (e_charge * n(i))) * dn_dr(i)
                ! Term 2: (ev/e) * dTi/dr (temperature gradient contribution, k=0)
                Er(i) = Er(i) + (ev / e_charge) * dTi_dr(i)
                ! Term 3: (r*B0*Vz)/(c*q*R0) (toroidal rotation contribution)
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
        filename = trim(profile_location) // '/' // trim(Er_file)
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
        integer, parameter :: MAX_PTS = 20000

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

    subroutine copy_profiles_if_needed()
        !> Copy Er.dat from input_profile_dir to profile_location when
        !> profiles are already in r_eff coordinates
        !> This ensures check_and_calculate_er() finds the provided Er file
        implicit none

        character(256) :: src_file, dst_file
        integer :: src_unit, dst_unit, ios
        logical :: src_exists, dst_exists, dirs_same
        character(4096) :: line

        ! Check if directories are the same (no copy needed)
        dirs_same = (trim(input_profile_dir) == trim(profile_location))
        if (dirs_same) then
            write(*,*) 'Profiles already in r_eff coordinates at: ', trim(profile_location)
            return
        end if

        ! Only Er.dat needs to be copied - other profiles are read from input_profile_dir
        src_file = trim(input_profile_dir) // '/' // trim(Er_file)
        dst_file = trim(profile_location) // '/' // trim(Er_file)

        inquire(file=trim(src_file), exist=src_exists)
        if (.not. src_exists) then
            ! No Er.dat provided - will be calculated later
            return
        end if

        ! Check if destination already exists
        inquire(file=trim(dst_file), exist=dst_exists)
        if (dst_exists) then
            write(*,*) 'Er.dat already exists at ', trim(profile_location)
            return
        end if

        ! Create output directory if it doesn't exist
        call execute_command_line('mkdir -p ' // trim(profile_location))

        ! Copy Er.dat
        open(newunit=src_unit, file=trim(src_file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'WARNING: Failed to open ', trim(src_file)
            return
        end if

        open(newunit=dst_unit, file=trim(dst_file), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            close(src_unit)
            write(*,*) 'WARNING: Failed to create ', trim(dst_file)
            return
        end if

        do
            read(src_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            write(dst_unit, '(A)') trim(line)
        end do

        close(src_unit)
        close(dst_unit)
        write(*,*) 'Copied Er.dat from ', trim(input_profile_dir), ' to ', trim(profile_location)

    end subroutine copy_profiles_if_needed

    subroutine validate_radial_range()
        !> Validate that r_min and r_plas are within the plasma profile region
        !> Reads n.dat (density profile) to determine the valid r_eff range
        implicit none

        real(dp), allocatable :: r(:), n(:)
        real(dp) :: r_min_profile, r_max_profile
        integer :: npts
        character(256) :: filename
        logical :: file_exists

        ! Read density profile to get radial range
        filename = trim(profile_location) // '/' // trim(n_file)
        inquire(file=trim(filename), exist=file_exists)

        if (.not. file_exists) then
            write(*,*) 'WARNING: Cannot validate radial range - density profile not found'
            write(*,*) '  Expected: ', trim(filename)
            return
        end if

        call read_profile_data(filename, r, n, npts)

        if (npts < 2) then
            write(*,*) 'WARNING: Cannot validate radial range - profile has < 2 points'
            deallocate(r, n)
            return
        end if

        r_min_profile = r(1)
        r_max_profile = r(npts)

        write(*,*) 'Validating radial range against plasma profiles:'
        write(*,*) '  Profile r_eff range: [', r_min_profile, ', ', r_max_profile, '] cm'
        write(*,*) '  Requested r_min    : ', r_min, ' cm'
        write(*,*) '  Requested r_plas   : ', r_plas, ' cm'

        ! Check r_min
        if (r_min < r_min_profile) then
            write(*,*) 'ERROR: r_min is below the profile region'
            write(*,*) '  r_min = ', r_min, ' cm'
            write(*,*) '  Profile starts at r_eff = ', r_min_profile, ' cm'
            write(*,*) '  Please increase r_min in KIM_GRID namelist'
            stop 1
        end if

        if (r_min > r_max_profile) then
            write(*,*) 'ERROR: r_min is beyond the profile region'
            write(*,*) '  r_min = ', r_min, ' cm'
            write(*,*) '  Profile ends at r_eff = ', r_max_profile, ' cm'
            write(*,*) '  Please decrease r_min in KIM_GRID namelist'
            stop 1
        end if

        ! Check r_plas
        if (r_plas < r_min_profile) then
            write(*,*) 'ERROR: r_plas is below the profile region'
            write(*,*) '  r_plas = ', r_plas, ' cm'
            write(*,*) '  Profile starts at r_eff = ', r_min_profile, ' cm'
            write(*,*) '  Please increase r_plas in KIM_GRID namelist'
            stop 1
        end if

        if (r_plas > r_max_profile) then
            write(*,*) 'ERROR: r_plas is beyond the profile region'
            write(*,*) '  r_plas = ', r_plas, ' cm'
            write(*,*) '  Profile ends at r_eff = ', r_max_profile, ' cm'
            write(*,*) '  Please decrease r_plas in KIM_GRID namelist'
            stop 1
        end if

        ! Additional check: r_min should be less than r_plas
        if (r_min >= r_plas) then
            write(*,*) 'ERROR: r_min must be less than r_plas'
            write(*,*) '  r_min = ', r_min, ' cm'
            write(*,*) '  r_plas = ', r_plas, ' cm'
            stop 1
        end if

        write(*,*) '  Radial range validation: PASSED'

        deallocate(r, n)

    end subroutine validate_radial_range

end module profile_input_m
