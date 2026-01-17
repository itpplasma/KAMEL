!> @file profile_preprocessor.f90
!> @brief Profile preprocessor for mapping plasma profiles from psi_N to r_eff
!>
!> This module maps plasma profiles (density, temperatures, rotation) from
!> normalized poloidal flux coordinates to effective radius coordinates.
!> It reads input profiles as functions of psi_N or sqrt(psi_N), uses
!> equilibrium data for the mapping, and outputs profiles in the format
!> expected by KiLCA/KIM/QL-Balance.

module profile_preprocessor_m

    use equil_profiles_m, only: equil_profiles_t, init_equilibrium_field

    implicit none
    private

    double precision, parameter :: pi = 3.14159265358979d0

    !> Coordinate type constants (public for namelist use)
    integer, parameter, public :: COORD_SQRT_PSIN = 1  !< sqrt(psi/psi_max) - default
    integer, parameter, public :: COORD_PSIN = 2       !< psi/psi_max

    !> Profile preprocessor derived type
    type, public :: profile_preprocessor_t
        !> Configuration from namelist
        character(len=512) :: equil_file = ''       !< Path to equil_r_q_psi.dat
        character(len=512) :: output_dir = '.'      !< Where to write output
        character(len=512) :: input_dir = '.'       !< Where to read input profiles
        integer :: coord_type = COORD_SQRT_PSIN     !< Input coordinate convention

        !> Input filenames (configurable)
        character(len=128) :: n_input_file = 'n_of_psiN.dat'
        character(len=128) :: Te_input_file = 'Te_of_psiN.dat'
        character(len=128) :: Ti_input_file = 'Ti_of_psiN.dat'
        character(len=128) :: Vz_input_file = 'Vz_of_psiN.dat'

        !> Equilibrium recomputation parameters
        integer :: nsqpsi = 100     !< Grid points for equilibrium
        integer :: nstep = 1000     !< Integration steps per circuit
        integer :: nsurfmax = 200   !< Surface search points
        integer :: niter = 10       !< Newton iterations

        !> Equilibrium mapping arrays
        integer :: nrad = 0
        double precision, allocatable :: r_eff(:)    !< Effective radius grid
        double precision, allocatable :: psi_n(:)    !< Normalized poloidal flux
        double precision :: psi_max = 0.d0           !< Normalization value (psi at separatrix)

        !> Output profiles on r_eff grid
        double precision, allocatable :: n(:)        !< Density profile
        double precision, allocatable :: Te(:)       !< Electron temperature
        double precision, allocatable :: Ti(:)       !< Ion temperature
        double precision, allocatable :: Vz(:)       !< Toroidal rotation
        double precision, allocatable :: q(:)        !< Safety factor profile

        !> Flags for which profiles are loaded
        logical :: has_n = .false.
        logical :: has_Te = .false.
        logical :: has_Ti = .false.
        logical :: has_Vz = .false.
        logical :: has_q = .false.

        logical :: initialized = .false.

    contains
        procedure :: init => profile_preprocessor_init
        procedure :: process => profile_preprocessor_process
        procedure :: write_output => profile_preprocessor_write_output
        procedure :: cleanup => profile_preprocessor_cleanup
        procedure :: get_n => profile_preprocessor_get_n
        procedure :: get_Te => profile_preprocessor_get_Te
        procedure :: get_Ti => profile_preprocessor_get_Ti
        procedure :: get_Vz => profile_preprocessor_get_Vz
        ! Private helper procedures
        procedure, private :: read_equil_file
        procedure, private :: load_from_equil_profiles
        procedure, private :: compute_equilibrium
        procedure, private :: process_single_profile
    end type profile_preprocessor_t

contains

    !---------------------------------------------------------------------------
    !> Initialize the profile preprocessor
    !>
    !> @param[in] namelist_file Path to namelist configuration file
    !> @param[in] equil_in Optional pre-computed equilibrium profiles
    !---------------------------------------------------------------------------
    subroutine profile_preprocessor_init(self, namelist_file, equil_in)
        implicit none
        class(profile_preprocessor_t), intent(inout) :: self
        character(len=*), intent(in) :: namelist_file
        type(equil_profiles_t), intent(in), optional :: equil_in

        type(equil_profiles_t) :: equil_local
        integer :: iunit, ios
        logical :: file_exists
        character(len=512) :: equil_file, output_dir, input_dir
        character(len=128) :: n_input_file, Te_input_file, Ti_input_file, Vz_input_file
        integer :: coord_type, nsqpsi, nstep, nsurfmax, niter

        namelist /profile_preprocessor/ equil_file, output_dir, input_dir, &
            coord_type, n_input_file, Te_input_file, Ti_input_file, Vz_input_file, &
            nsqpsi, nstep, nsurfmax, niter

        ! Clean up any previous state
        call self%cleanup()

        ! Set defaults for namelist variables
        equil_file = self%equil_file
        output_dir = self%output_dir
        input_dir = self%input_dir
        coord_type = self%coord_type
        n_input_file = self%n_input_file
        Te_input_file = self%Te_input_file
        Ti_input_file = self%Ti_input_file
        Vz_input_file = self%Vz_input_file
        nsqpsi = self%nsqpsi
        nstep = self%nstep
        nsurfmax = self%nsurfmax
        niter = self%niter

        ! Read namelist
        inquire(file=namelist_file, exist=file_exists)
        if (.not. file_exists) then
            write(*,*) '[profile_preprocessor_m:init] ERROR: Namelist file not found: ', &
                       trim(namelist_file)
            stop 1
        end if

        open(newunit=iunit, file=namelist_file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[profile_preprocessor_m:init] ERROR: Cannot open namelist file: ', &
                       trim(namelist_file)
            stop 1
        end if

        read(iunit, nml=profile_preprocessor, iostat=ios)
        if (ios /= 0 .and. ios /= -1) then  ! -1 is EOF, which is OK
            write(*,*) '[profile_preprocessor_m:init] ERROR: Error reading namelist, iostat=', ios
            close(iunit)
            stop 1
        end if
        close(iunit)

        ! Store configuration
        self%equil_file = equil_file
        self%output_dir = output_dir
        self%input_dir = input_dir
        self%coord_type = coord_type
        self%n_input_file = n_input_file
        self%Te_input_file = Te_input_file
        self%Ti_input_file = Ti_input_file
        self%Vz_input_file = Vz_input_file
        self%nsqpsi = nsqpsi
        self%nstep = nstep
        self%nsurfmax = nsurfmax
        self%niter = niter

        ! Validate coord_type
        if (coord_type /= COORD_SQRT_PSIN .and. coord_type /= COORD_PSIN) then
            write(*,*) '[profile_preprocessor_m:init] ERROR: Invalid coord_type:', coord_type
            write(*,*) '  Use 1 for sqrt(psi_N) or 2 for psi_N'
            stop 1
        end if

        ! Obtain equilibrium mapping
        if (present(equil_in)) then
            ! Use provided equilibrium
            call self%load_from_equil_profiles(equil_in)
        else if (len_trim(equil_file) > 0) then
            ! Try to read from file
            inquire(file=equil_file, exist=file_exists)
            if (file_exists) then
                call self%read_equil_file(equil_file)
            else
                write(*,*) '[profile_preprocessor_m:init] WARNING: equil_file not found: ', &
                           trim(equil_file)
                write(*,*) '  Falling back to equilibrium recomputation'
                call self%compute_equilibrium(equil_local)
                call self%load_from_equil_profiles(equil_local)
                call equil_local%cleanup()
            end if
        else
            ! Compute equilibrium from scratch
            call self%compute_equilibrium(equil_local)
            call self%load_from_equil_profiles(equil_local)
            call equil_local%cleanup()
        end if

        ! Allocate output profile arrays
        allocate(self%n(self%nrad))
        allocate(self%Te(self%nrad))
        allocate(self%Ti(self%nrad))
        allocate(self%Vz(self%nrad))
        self%n = 0.d0
        self%Te = 0.d0
        self%Ti = 0.d0
        self%Vz = 0.d0

        self%initialized = .true.

    end subroutine profile_preprocessor_init

    !---------------------------------------------------------------------------
    !> Read equilibrium mapping from equil_r_q_psi.dat file
    !---------------------------------------------------------------------------
    subroutine read_equil_file(self, filename)
        implicit none
        class(profile_preprocessor_t), intent(inout) :: self
        character(len=*), intent(in) :: filename

        integer :: iunit, ios, i, nlines
        character(len=1024) :: line
        double precision :: r_eff_val, q_val, psi_val
        double precision, allocatable :: r_eff_tmp(:), psi_tmp(:), q_tmp(:)

        ! First pass: count data lines
        open(newunit=iunit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[profile_preprocessor_m:read_equil_file] ERROR: Cannot open file: ', &
                       trim(filename)
            stop 1
        end if

        nlines = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            nlines = nlines + 1
        end do
        close(iunit)

        if (nlines < 2) then
            write(*,*) '[profile_preprocessor_m:read_equil_file] ERROR: File has < 2 data points'
            stop 1
        end if

        ! Allocate arrays
        self%nrad = nlines
        allocate(self%r_eff(nlines))
        allocate(self%psi_n(nlines))
        allocate(self%q(nlines))
        allocate(r_eff_tmp(nlines), psi_tmp(nlines), q_tmp(nlines))

        ! Second pass: read data
        open(newunit=iunit, file=filename, status='old', action='read', iostat=ios)
        i = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle

            i = i + 1
            ! Format: r_eff, q, psi, phi, dphidpsi, rsmall, volume [, R_beg, Z_beg, R_min, R_max]
            read(line, *, iostat=ios) r_eff_val, q_val, psi_val
            if (ios /= 0) then
                write(*,*) '[profile_preprocessor_m:read_equil_file] ERROR: Parse error at line ', i
                stop 1
            end if

            r_eff_tmp(i) = r_eff_val
            q_tmp(i) = q_val
            psi_tmp(i) = psi_val
        end do
        close(iunit)

        ! Normalize psi
        self%psi_max = psi_tmp(nlines)
        if (abs(self%psi_max) < 1.d-30) then
            write(*,*) '[profile_preprocessor_m:read_equil_file] ERROR: psi_max is zero'
            stop 1
        end if

        do i = 1, nlines
            self%r_eff(i) = r_eff_tmp(i)
            self%q(i) = q_tmp(i)
            self%psi_n(i) = psi_tmp(i) / self%psi_max
        end do

        self%has_q = .true.

        deallocate(r_eff_tmp, psi_tmp, q_tmp)

    end subroutine read_equil_file

    !---------------------------------------------------------------------------
    !> Load equilibrium mapping from equil_profiles_t
    !---------------------------------------------------------------------------
    subroutine load_from_equil_profiles(self, equil)
        implicit none
        class(profile_preprocessor_t), intent(inout) :: self
        type(equil_profiles_t), intent(in) :: equil

        integer :: i

        self%nrad = equil%nsqpsi
        allocate(self%r_eff(self%nrad))
        allocate(self%psi_n(self%nrad))
        allocate(self%q(self%nrad))

        ! Validate equilibrium data BEFORE using it
        if (abs(equil%btor) < 1.d-10) then
            write(*,*) '[profile_preprocessor_m:load_from_equil] ERROR: btor is zero or too small'
            write(*,*) '  btor = ', equil%btor
            write(*,*) '  The equilibrium computation may have failed.'
            stop 1
        end if

        ! r_eff from equil is sqrt(2*phi_tor/B_tor)
        ! psi is the poloidal flux relative to axis
        self%psi_max = equil%psisurf(equil%nsqpsi)

        if (abs(self%psi_max) < 1.d-20) then
            write(*,*) '[profile_preprocessor_m:load_from_equil] ERROR: psi_max is zero or too small'
            write(*,*) '  psi_max = ', self%psi_max
            write(*,*) '  psisurf(nsqpsi) = ', equil%psisurf(equil%nsqpsi)
            write(*,*) '  The equilibrium computation may have failed.'
            write(*,*) '  Check that the GEQDSK file is valid and the convexwall file is correct.'
            stop 1
        end if

        do i = 1, self%nrad
            self%r_eff(i) = sqrt(2.d0 * abs(equil%phitor(i) / equil%btor))
            self%psi_n(i) = equil%psisurf(i) / self%psi_max
            self%q(i) = equil%qsaf(i)
        end do

        self%has_q = .true.

        ! Check for NaN in computed values
        do i = 1, self%nrad
            if (self%r_eff(i) /= self%r_eff(i) .or. self%psi_n(i) /= self%psi_n(i)) then
                write(*,*) '[profile_preprocessor_m:load_from_equil] ERROR: NaN at index ', i
                write(*,*) '  r_eff = ', self%r_eff(i), ', psi_n = ', self%psi_n(i)
                write(*,*) '  phitor = ', equil%phitor(i), ', psisurf = ', equil%psisurf(i)
                stop 1
            end if
        end do

    end subroutine load_from_equil_profiles

    !---------------------------------------------------------------------------
    !> Compute equilibrium profiles from scratch
    !---------------------------------------------------------------------------
    subroutine compute_equilibrium(self, equil)
        implicit none
        class(profile_preprocessor_t), intent(in) :: self
        type(equil_profiles_t), intent(out) :: equil

        character(len=512) :: output_file

        write(*,*) '[profile_preprocessor_m:compute_equilibrium] Computing equilibrium profiles...'

        ! Initialize the magnetic field from geqdsk
        call init_equilibrium_field()

        ! Initialize and compute equilibrium
        call equil%init(self%nsqpsi)
        call equil%find_axis(self%nstep, 10, self%niter)
        call equil%compute_profiles(self%nstep, self%nsurfmax, self%niter)

        ! Write output file
        output_file = trim(self%output_dir) // '/equil_r_q_psi.dat'
        call equil%write_output(self%output_dir)

        write(*,*) '[profile_preprocessor_m:compute_equilibrium] Wrote: ', trim(output_file)

    end subroutine compute_equilibrium

    !---------------------------------------------------------------------------
    !> Process input profiles: read, interpolate, and store
    !---------------------------------------------------------------------------
    subroutine profile_preprocessor_process(self)
        implicit none
        class(profile_preprocessor_t), intent(inout) :: self

        if (.not. self%initialized) then
            write(*,*) '[profile_preprocessor_m:process] ERROR: Not initialized'
            stop 1
        end if

        ! Process each profile
        call self%process_single_profile(self%n_input_file, self%n, self%has_n, 'n')
        call self%process_single_profile(self%Te_input_file, self%Te, self%has_Te, 'Te')
        call self%process_single_profile(self%Ti_input_file, self%Ti, self%has_Ti, 'Ti')
        call self%process_single_profile(self%Vz_input_file, self%Vz, self%has_Vz, 'Vz')

    end subroutine profile_preprocessor_process

    !---------------------------------------------------------------------------
    !> Process a single profile file
    !---------------------------------------------------------------------------
    subroutine process_single_profile(self, filename, profile_out, has_profile, profile_name)
        implicit none
        class(profile_preprocessor_t), intent(inout) :: self
        character(len=*), intent(in) :: filename
        double precision, intent(out) :: profile_out(:)
        logical, intent(out) :: has_profile
        character(len=*), intent(in) :: profile_name

        character(len=1024) :: filepath, line
        integer :: iunit, ios, nlines, i, n_in
        double precision, allocatable :: coord_in(:), value_in(:), psi_n_in(:)
        double precision :: coord_val, value_val
        double precision :: psi_n_min_in, psi_n_max_in
        double precision :: psi_n_min_eq, psi_n_max_eq
        logical :: file_exists, truncated, has_negative

        has_profile = .false.
        profile_out = 0.d0

        ! Construct full path
        filepath = trim(self%input_dir) // '/' // trim(filename)

        inquire(file=filepath, exist=file_exists)
        if (.not. file_exists) then
            write(*,*) '[profile_preprocessor_m:process] WARNING: File not found: ', trim(filepath)
            write(*,*) '  Skipping profile: ', trim(profile_name)
            return
        end if

        ! First pass: count data lines
        open(newunit=iunit, file=filepath, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[profile_preprocessor_m:process] WARNING: Cannot open: ', trim(filepath)
            return
        end if

        nlines = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (is_data_line(line)) nlines = nlines + 1
        end do
        close(iunit)

        if (nlines < 2) then
            write(*,*) '[profile_preprocessor_m:process] ERROR: Profile ', trim(profile_name), &
                       ' has < 2 data points'
            stop 1
        end if

        ! Allocate input arrays
        n_in = nlines
        allocate(coord_in(n_in), value_in(n_in), psi_n_in(n_in))

        ! Second pass: read data
        open(newunit=iunit, file=filepath, status='old', action='read', iostat=ios)
        i = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (.not. is_data_line(line)) cycle

            i = i + 1
            read(line, *, iostat=ios) coord_val, value_val
            if (ios /= 0) then
                write(*,*) '[profile_preprocessor_m:process] ERROR: Parse error in ', &
                           trim(profile_name), ' at data line ', i
                stop 1
            end if

            ! Check for NaN/Inf
            if (value_val /= value_val .or. abs(value_val) > huge(1.d0)/2) then
                write(*,*) '[profile_preprocessor_m:process] ERROR: NaN/Inf in ', &
                           trim(profile_name), ' at line ', i
                stop 1
            end if

            coord_in(i) = coord_val
            value_in(i) = value_val
        end do
        close(iunit)

        ! Check monotonicity
        do i = 2, n_in
            if (coord_in(i) <= coord_in(i-1)) then
                write(*,*) '[profile_preprocessor_m:process] ERROR: Non-monotonic coordinate in ', &
                           trim(profile_name), ' at line ', i
                stop 1
            end if
        end do

        ! Convert to psi_n based on coord_type
        if (self%coord_type == COORD_SQRT_PSIN) then
            psi_n_in = coord_in**2
        else
            psi_n_in = coord_in
        end if

        ! Check for negative values
        has_negative = .false.
        do i = 1, n_in
            if (value_in(i) < 0.d0) has_negative = .true.
        end do

        if (has_negative) then
            if (profile_name == 'n') then
                write(*,*) '[profile_preprocessor_m:process] WARNING: Negative values in density profile'
            else if (profile_name == 'Te' .or. profile_name == 'Ti') then
                write(*,*) '[profile_preprocessor_m:process] WARNING: Negative values in ', &
                           trim(profile_name), ' profile'
            end if
        end if

        ! Check bounds
        psi_n_min_in = psi_n_in(1)
        psi_n_max_in = psi_n_in(n_in)
        psi_n_min_eq = self%psi_n(1)
        psi_n_max_eq = self%psi_n(self%nrad)

        truncated = .false.
        if (psi_n_max_in > psi_n_max_eq * 1.001d0) then
            write(*,*) '[profile_preprocessor_m:process] WARNING: Profile ', trim(profile_name), &
                       ' extends beyond equilibrium'
            write(*,*) '  Input psi_n_max = ', psi_n_max_in, ', equilibrium psi_n_max = ', psi_n_max_eq
            write(*,*) '  Truncating to equilibrium boundary'
            truncated = .true.
        end if

        if (psi_n_min_in > psi_n_min_eq + 0.1d0 * (psi_n_max_eq - psi_n_min_eq)) then
            write(*,*) '[profile_preprocessor_m:process] ERROR: Profile ', trim(profile_name), &
                       ' does not cover equilibrium range'
            write(*,*) '  Gap at start exceeds 10%'
            stop 1
        end if

        ! Debug: print psi_n ranges
        write(*,*) '[profile_preprocessor_m:process] Input psi_n range: ', &
                   psi_n_in(1), ' to ', psi_n_in(n_in)
        write(*,*) '[profile_preprocessor_m:process] Equil psi_n range: ', &
                   self%psi_n(1), ' to ', self%psi_n(self%nrad)

        ! Interpolate using cubic spline
        call cubic_spline_interpolate(psi_n_in, value_in, n_in, &
                                      self%psi_n, profile_out, self%nrad)

        ! Check for NaN in interpolated output
        do i = 1, self%nrad
            if (profile_out(i) /= profile_out(i)) then
                write(*,*) '[profile_preprocessor_m:process] ERROR: NaN in interpolated ', &
                           trim(profile_name), ' at index ', i
                write(*,*) '  psi_n(i) = ', self%psi_n(i), ', r_eff(i) = ', self%r_eff(i)
                write(*,*) '  This may indicate a coordinate system mismatch between'
                write(*,*) '  the input profiles and the computed equilibrium.'
                stop 1
            end if
        end do

        has_profile = .true.
        deallocate(coord_in, value_in, psi_n_in)

        write(*,*) '[profile_preprocessor_m:process] Loaded profile: ', trim(profile_name), &
                   ' (', n_in, ' points)'

    end subroutine process_single_profile

    !---------------------------------------------------------------------------
    !> Check if a line contains data (not a header/comment)
    !---------------------------------------------------------------------------
    function is_data_line(line) result(is_data)
        implicit none
        character(len=*), intent(in) :: line
        logical :: is_data
        character(len=1024) :: trimmed
        character(len=1) :: first_char

        is_data = .false.
        trimmed = adjustl(line)

        if (len_trim(trimmed) == 0) return
        first_char = trimmed(1:1)

        ! Skip comments and headers
        if (first_char == '#') return
        if (first_char == '!') return
        if (first_char == '%') return

        ! Check if first character could start a number
        if (first_char == '-' .or. first_char == '+' .or. first_char == '.' .or. &
            (ichar(first_char) >= ichar('0') .and. ichar(first_char) <= ichar('9'))) then
            is_data = .true.
        end if

    end function is_data_line

    !---------------------------------------------------------------------------
    !> Cubic spline interpolation
    !>
    !> Natural cubic spline with zero second derivatives at boundaries
    !---------------------------------------------------------------------------
    subroutine cubic_spline_interpolate(x_in, y_in, n_in, x_out, y_out, n_out)
        implicit none
        integer, intent(in) :: n_in, n_out
        double precision, intent(in) :: x_in(n_in), y_in(n_in)
        double precision, intent(in) :: x_out(n_out)
        double precision, intent(out) :: y_out(n_out)

        double precision, allocatable :: h(:), alpha(:), l(:), mu(:), z(:), c(:), b(:), d(:)
        double precision :: x, dx
        integer :: i, j, k

        ! Compute spline coefficients
        allocate(h(n_in-1), alpha(n_in), l(n_in), mu(n_in), z(n_in))
        allocate(c(n_in), b(n_in-1), d(n_in-1))

        ! Step 1: Compute h_i = x_{i+1} - x_i
        do i = 1, n_in - 1
            h(i) = x_in(i+1) - x_in(i)
        end do

        ! Step 2: Compute alpha
        alpha(1) = 0.d0
        alpha(n_in) = 0.d0
        do i = 2, n_in - 1
            alpha(i) = 3.d0/h(i) * (y_in(i+1) - y_in(i)) - &
                       3.d0/h(i-1) * (y_in(i) - y_in(i-1))
        end do

        ! Step 3: Solve tridiagonal system for c
        l(1) = 1.d0
        mu(1) = 0.d0
        z(1) = 0.d0

        do i = 2, n_in - 1
            l(i) = 2.d0 * (x_in(i+1) - x_in(i-1)) - h(i-1) * mu(i-1)
            mu(i) = h(i) / l(i)
            z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i)
        end do

        l(n_in) = 1.d0
        z(n_in) = 0.d0
        c(n_in) = 0.d0

        do j = n_in - 1, 1, -1
            c(j) = z(j) - mu(j) * c(j+1)
            b(j) = (y_in(j+1) - y_in(j)) / h(j) - h(j) * (c(j+1) + 2.d0*c(j)) / 3.d0
            d(j) = (c(j+1) - c(j)) / (3.d0 * h(j))
        end do

        ! Evaluate spline at output points
        do i = 1, n_out
            x = x_out(i)

            ! Clamp to input range
            if (x <= x_in(1)) then
                y_out(i) = y_in(1)
                cycle
            else if (x >= x_in(n_in)) then
                y_out(i) = y_in(n_in)
                cycle
            end if

            ! Find interval (binary search)
            k = 1
            do j = 1, n_in - 1
                if (x >= x_in(j) .and. x < x_in(j+1)) then
                    k = j
                    exit
                end if
            end do

            ! Evaluate cubic polynomial
            dx = x - x_in(k)
            y_out(i) = y_in(k) + b(k)*dx + c(k)*dx**2 + d(k)*dx**3
        end do

        deallocate(h, alpha, l, mu, z, c, b, d)

    end subroutine cubic_spline_interpolate

    !---------------------------------------------------------------------------
    !> Write output profiles to files
    !---------------------------------------------------------------------------
    subroutine profile_preprocessor_write_output(self)
        implicit none
        class(profile_preprocessor_t), intent(in) :: self

        character(len=512) :: filepath

        if (.not. self%initialized) then
            write(*,*) '[profile_preprocessor_m:write_output] ERROR: Not initialized'
            stop 1
        end if

        ! Check if any profiles are loaded
        if (.not. (self%has_n .or. self%has_Te .or. self%has_Ti .or. self%has_Vz .or. self%has_q)) then
            write(*,*) '[profile_preprocessor_m:write_output] WARNING: No profiles to write'
            return
        end if

        ! Write each profile that was loaded
        if (self%has_n) then
            filepath = trim(self%output_dir) // '/n.dat'
            call write_profile_file(filepath, self%r_eff, self%n, self%nrad)
            write(*,*) '[profile_preprocessor_m:write_output] Wrote: ', trim(filepath)
        end if

        if (self%has_Te) then
            filepath = trim(self%output_dir) // '/Te.dat'
            call write_profile_file(filepath, self%r_eff, self%Te, self%nrad)
            write(*,*) '[profile_preprocessor_m:write_output] Wrote: ', trim(filepath)
        end if

        if (self%has_Ti) then
            filepath = trim(self%output_dir) // '/Ti.dat'
            call write_profile_file(filepath, self%r_eff, self%Ti, self%nrad)
            write(*,*) '[profile_preprocessor_m:write_output] Wrote: ', trim(filepath)
        end if

        if (self%has_Vz) then
            filepath = trim(self%output_dir) // '/Vz.dat'
            call write_profile_file(filepath, self%r_eff, self%Vz, self%nrad)
            write(*,*) '[profile_preprocessor_m:write_output] Wrote: ', trim(filepath)
        end if

        ! Write q.dat (safety factor profile)
        if (self%has_q) then
            filepath = trim(self%output_dir) // '/q.dat'
            call write_profile_file(filepath, self%r_eff, self%q, self%nrad)
            write(*,*) '[profile_preprocessor_m:write_output] Wrote: ', trim(filepath)
        end if

    end subroutine profile_preprocessor_write_output

    !---------------------------------------------------------------------------
    !> Write a single profile to file
    !---------------------------------------------------------------------------
    subroutine write_profile_file(filepath, r_eff, profile, n)
        implicit none
        character(len=*), intent(in) :: filepath
        double precision, intent(in) :: r_eff(n), profile(n)
        integer, intent(in) :: n

        integer :: iunit, ios, i

        open(newunit=iunit, file=filepath, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write(*,*) '[profile_preprocessor_m:write_profile_file] ERROR: Cannot write to: ', &
                       trim(filepath)
            stop 1
        end if

        do i = 1, n
            write(iunit, '(ES23.15, 1X, ES23.15)') r_eff(i), profile(i)
        end do

        close(iunit)

    end subroutine write_profile_file

    !---------------------------------------------------------------------------
    !> Get density at arbitrary r_eff (interpolation)
    !---------------------------------------------------------------------------
    function profile_preprocessor_get_n(self, r) result(val)
        implicit none
        class(profile_preprocessor_t), intent(in) :: self
        double precision, intent(in) :: r
        double precision :: val

        if (.not. self%initialized) then
            write(*,*) '[profile_preprocessor_m:get_n] ERROR: Not initialized'
            stop 1
        end if

        if (.not. self%has_n) then
            write(*,*) '[profile_preprocessor_m:get_n] ERROR: Density profile not loaded'
            stop 1
        end if

        call interpolate_at_r(self%r_eff, self%n, self%nrad, r, val)

    end function profile_preprocessor_get_n

    !---------------------------------------------------------------------------
    !> Get electron temperature at arbitrary r_eff
    !---------------------------------------------------------------------------
    function profile_preprocessor_get_Te(self, r) result(val)
        implicit none
        class(profile_preprocessor_t), intent(in) :: self
        double precision, intent(in) :: r
        double precision :: val

        if (.not. self%initialized) then
            write(*,*) '[profile_preprocessor_m:get_Te] ERROR: Not initialized'
            stop 1
        end if

        if (.not. self%has_Te) then
            write(*,*) '[profile_preprocessor_m:get_Te] ERROR: Te profile not loaded'
            stop 1
        end if

        call interpolate_at_r(self%r_eff, self%Te, self%nrad, r, val)

    end function profile_preprocessor_get_Te

    !---------------------------------------------------------------------------
    !> Get ion temperature at arbitrary r_eff
    !---------------------------------------------------------------------------
    function profile_preprocessor_get_Ti(self, r) result(val)
        implicit none
        class(profile_preprocessor_t), intent(in) :: self
        double precision, intent(in) :: r
        double precision :: val

        if (.not. self%initialized) then
            write(*,*) '[profile_preprocessor_m:get_Ti] ERROR: Not initialized'
            stop 1
        end if

        if (.not. self%has_Ti) then
            write(*,*) '[profile_preprocessor_m:get_Ti] ERROR: Ti profile not loaded'
            stop 1
        end if

        call interpolate_at_r(self%r_eff, self%Ti, self%nrad, r, val)

    end function profile_preprocessor_get_Ti

    !---------------------------------------------------------------------------
    !> Get toroidal rotation at arbitrary r_eff
    !---------------------------------------------------------------------------
    function profile_preprocessor_get_Vz(self, r) result(val)
        implicit none
        class(profile_preprocessor_t), intent(in) :: self
        double precision, intent(in) :: r
        double precision :: val

        if (.not. self%initialized) then
            write(*,*) '[profile_preprocessor_m:get_Vz] ERROR: Not initialized'
            stop 1
        end if

        if (.not. self%has_Vz) then
            write(*,*) '[profile_preprocessor_m:get_Vz] ERROR: Vz profile not loaded'
            stop 1
        end if

        call interpolate_at_r(self%r_eff, self%Vz, self%nrad, r, val)

    end function profile_preprocessor_get_Vz

    !---------------------------------------------------------------------------
    !> Interpolate profile at given r using cubic spline
    !---------------------------------------------------------------------------
    subroutine interpolate_at_r(r_grid, profile, n, r, val)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: r_grid(n), profile(n)
        double precision, intent(in) :: r
        double precision, intent(out) :: val

        double precision :: r_arr(1), val_arr(1)

        r_arr(1) = r
        call cubic_spline_interpolate(r_grid, profile, n, r_arr, val_arr, 1)
        val = val_arr(1)

        ! Warning for out-of-range
        if (r < r_grid(1) .or. r > r_grid(n)) then
            write(*,*) '[profile_preprocessor_m:interpolate_at_r] WARNING: r=', r, &
                       ' outside grid range [', r_grid(1), ',', r_grid(n), ']'
        end if

    end subroutine interpolate_at_r

    !---------------------------------------------------------------------------
    !> Clean up allocated memory
    !---------------------------------------------------------------------------
    subroutine profile_preprocessor_cleanup(self)
        implicit none
        class(profile_preprocessor_t), intent(inout) :: self

        if (allocated(self%r_eff)) deallocate(self%r_eff)
        if (allocated(self%psi_n)) deallocate(self%psi_n)
        if (allocated(self%n)) deallocate(self%n)
        if (allocated(self%Te)) deallocate(self%Te)
        if (allocated(self%Ti)) deallocate(self%Ti)
        if (allocated(self%Vz)) deallocate(self%Vz)
        if (allocated(self%q)) deallocate(self%q)

        self%has_n = .false.
        self%has_Te = .false.
        self%has_Ti = .false.
        self%has_Vz = .false.
        self%has_q = .false.
        self%initialized = .false.
        self%nrad = 0

    end subroutine profile_preprocessor_cleanup

end module profile_preprocessor_m
