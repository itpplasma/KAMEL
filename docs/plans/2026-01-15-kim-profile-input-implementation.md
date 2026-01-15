# KIM Profile Input Module Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Integrate profile_preprocessor into KIM to support input profiles in sqrt(psi_N) coordinates with automatic coordinate detection and Er force balance calculation.

**Architecture:** New `profile_input_m` module called early in KIM initialization (before `read_profiles()`). Detects coordinate type, runs preprocessing if needed, calculates Er from force balance when missing.

**Tech Stack:** Fortran 90, kamel_equil library (profile_preprocessor_m), CGS units

---

## Task 1: Add KIM_profiles Namelist Variables

**Files:**
- Modify: `KIM/src/setup/config_mod.f90:27-40`

**Step 1: Write the test (manual verification)**

No unit test needed - namelist variables are data declarations. Verification via compilation.

**Step 2: Add namelist variables to config_mod.f90**

Add after line 40 (after `ion_flr_scale_factor`):

```fortran
    ! KIM_PROFILES namelist variables
    character(20) :: coord_type = 'auto'           ! 'auto', 'sqrt_psiN', or 'r_eff'
    character(256) :: input_profile_dir = './'     ! Directory for raw input profiles
    character(256) :: equil_file = ''              ! Path to equilibrium file (empty = compute)
```

**Step 3: Build to verify compilation**

Run: `make KIM`
Expected: BUILD SUCCESS

**Step 4: Commit**

```bash
git add KIM/src/setup/config_mod.f90
git commit -m "feat(KIM): add KIM_profiles namelist variables to config_mod"
```

---

## Task 2: Add KIM_profiles Namelist Reading

**Files:**
- Modify: `KIM/src/general/read_config.f90:16-31` (namelist declarations)
- Modify: `KIM/src/general/read_config.f90:56-62` (namelist reading)

**Step 1: Add namelist declaration**

After the KIM_GRID namelist declaration (around line 40), add:

```fortran
    namelist /KIM_PROFILES/ coord_type, input_profile_dir, equil_file
```

**Step 2: Add namelist reading**

After line 61 (`read(unit = 77, nml = KIM_GRID)`), add:

```fortran
    read(unit = 77, nml = KIM_PROFILES)
```

**Step 3: Build to verify compilation**

Run: `make KIM`
Expected: BUILD SUCCESS

**Step 4: Commit**

```bash
git add KIM/src/general/read_config.f90
git commit -m "feat(KIM): read KIM_profiles namelist in read_config"
```

---

## Task 3: Add KIM_profiles to Default Namelist File

**Files:**
- Modify: `KIM/nmls/KIM_config.nml`

**Step 1: Add namelist group**

Add at the end of the file:

```fortran
&KIM_PROFILES
    ! Coordinate type for input profiles
    ! 'auto'      - detect automatically (default)
    ! 'sqrt_psiN' - profiles are in sqrt(psi/psi_max) coordinates
    ! 'r_eff'     - profiles are already in effective radius [cm]
    coord_type = 'auto'

    ! Input directory for raw profiles (n_of_psiN.dat, etc.)
    input_profile_dir = './'

    ! Equilibrium file path (empty = compute from geqdsk)
    equil_file = ''
/
```

**Step 2: Verify namelist file syntax**

Run: `make KIM`
Expected: BUILD SUCCESS (namelist file not read at compile time, but confirms no syntax errors in Fortran)

**Step 3: Commit**

```bash
git add KIM/nmls/KIM_config.nml
git commit -m "feat(KIM): add KIM_profiles namelist group to default config"
```

---

## Task 4: Create profile_input_m Module - Core Structure

**Files:**
- Create: `KIM/src/background_equilibrium/profile_input_m.f90`

**Step 1: Create module skeleton**

```fortran
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
        character(20), intent(out) :: detected_type

        detected_type = 'r_eff'  ! Placeholder - will implement
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
```

**Step 2: Build to verify compilation**

Run: `make KIM`
Expected: BUILD SUCCESS (module compiles but not yet linked)

**Step 3: Commit**

```bash
git add KIM/src/background_equilibrium/profile_input_m.f90
git commit -m "feat(KIM): create profile_input_m module skeleton"
```

---

## Task 5: Update CMakeLists.txt - Add Source and Link kamel_equil

**Files:**
- Modify: `KIM/src/CMakeLists.txt`

**Step 1: Read current CMakeLists.txt to find source list**

Read file to locate where to add new source.

**Step 2: Add source file to KIM_lib**

Find the `add_library(KIM_lib ...)` or source file list and add:
`background_equilibrium/profile_input_m.f90`

**Step 3: Link kamel_equil library**

Add after `target_link_libraries(KIM_lib ...`:
```cmake
target_link_libraries(KIM_lib PUBLIC kamel_equil)
```

**Step 4: Build to verify**

Run: `make KIM`
Expected: BUILD SUCCESS

**Step 5: Commit**

```bash
git add KIM/src/CMakeLists.txt
git commit -m "feat(KIM): add profile_input_m to build and link kamel_equil"
```

---

## Task 6: Implement detect_coordinate_type

**Files:**
- Modify: `KIM/src/background_equilibrium/profile_input_m.f90`

**Step 1: Write the failing test**

Create `KIM/tests/test_profile_input.f90`:

```fortran
program test_profile_input
    use KIM_kinds_m, only: dp
    implicit none

    logical :: all_passed
    all_passed = .true.

    call test_coord_detection_reff(all_passed)
    call test_coord_detection_psiN(all_passed)

    if (all_passed) then
        print *, 'All tests PASSED'
        stop 0
    else
        print *, 'Some tests FAILED'
        stop 1
    end if

contains

    subroutine test_coord_detection_reff(passed)
        logical, intent(inout) :: passed
        real(dp) :: test_data(5)
        character(20) :: result

        ! r_eff data: max > 2.0
        test_data = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 45.0_dp]
        call detect_from_data(test_data, 5, result)

        if (trim(result) /= 'r_eff') then
            print *, 'FAIL: test_coord_detection_reff'
            print *, '  Expected: r_eff'
            print *, '  Got: ', trim(result)
            passed = .false.
        else
            print *, 'PASS: test_coord_detection_reff'
        end if
    end subroutine

    subroutine test_coord_detection_psiN(passed)
        logical, intent(inout) :: passed
        real(dp) :: test_data(5)
        character(20) :: result

        ! sqrt_psiN data: max <= 2.0
        test_data = [0.1_dp, 0.3_dp, 0.5_dp, 0.7_dp, 0.95_dp]
        call detect_from_data(test_data, 5, result)

        if (trim(result) /= 'sqrt_psiN') then
            print *, 'FAIL: test_coord_detection_psiN'
            print *, '  Expected: sqrt_psiN'
            print *, '  Got: ', trim(result)
            passed = .false.
        else
            print *, 'PASS: test_coord_detection_psiN'
        end if
    end subroutine

    subroutine detect_from_data(data, n, coord_type)
        real(dp), intent(in) :: data(:)
        integer, intent(in) :: n
        character(20), intent(out) :: coord_type
        real(dp) :: max_val
        real(dp), parameter :: COORD_THRESHOLD = 2.0_dp

        max_val = maxval(data(1:n))
        if (max_val > COORD_THRESHOLD) then
            coord_type = 'r_eff'
        else
            coord_type = 'sqrt_psiN'
        end if
    end subroutine detect_from_data

end program test_profile_input
```

**Step 2: Add test to CMakeLists.txt**

Add to `KIM/tests/CMakeLists.txt`:
```cmake
add_executable(test_profile_input test_profile_input.f90)
target_link_libraries(test_profile_input KIM_lib)
add_test(NAME test_profile_input COMMAND test_profile_input)
```

**Step 3: Run test to verify it compiles and passes**

Run: `make test`
Expected: test_profile_input PASSED

**Step 4: Implement detect_coordinate_type in module**

Update `profile_input_m.f90`:

```fortran
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
```

**Step 5: Build and run tests**

Run: `make KIM && make test`
Expected: BUILD SUCCESS, test_profile_input PASSED

**Step 6: Commit**

```bash
git add KIM/src/background_equilibrium/profile_input_m.f90 KIM/tests/test_profile_input.f90 KIM/tests/CMakeLists.txt
git commit -m "feat(KIM): implement coordinate type detection"
```

---

## Task 7: Implement validate_btor_rbig

**Files:**
- Modify: `KIM/src/background_equilibrium/profile_input_m.f90`

**Step 1: Implement btor/R0 validation**

```fortran
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
```

**Step 2: Build to verify**

Run: `make KIM`
Expected: BUILD SUCCESS

**Step 3: Commit**

```bash
git add KIM/src/background_equilibrium/profile_input_m.f90
git commit -m "feat(KIM): implement btor/R0 validation against btor_rbig.dat"
```

---

## Task 8: Implement Er Force Balance Calculation

**Files:**
- Modify: `KIM/src/background_equilibrium/profile_input_m.f90`

**Step 1: Add constants module use**

Add to module header:
```fortran
    use constants_m, only: e_charge, c_light
```

**Step 2: Implement check_and_calculate_er**

```fortran
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
                    Er(i) = Er(i) + (r(i) * btor * Vz(i)) / (c_light * q(i) * R0)
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

        ! Also create Er.dat symlink or copy for read_profiles to find
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

        integer :: iunit, ios, i
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
```

**Step 3: Build to verify**

Run: `make KIM`
Expected: BUILD SUCCESS

**Step 4: Commit**

```bash
git add KIM/src/background_equilibrium/profile_input_m.f90
git commit -m "feat(KIM): implement Er force balance calculation"
```

---

## Task 9: Implement run_preprocessing (Profile Transformation)

**Files:**
- Modify: `KIM/src/background_equilibrium/profile_input_m.f90`

**Step 1: Add profile_preprocessor module use**

Add to module header:
```fortran
    use profile_preprocessor_m, only: profile_preprocessor_t
```

**Step 2: Implement run_preprocessing**

```fortran
    subroutine run_preprocessing()
        !> Run profile_preprocessor for sqrt_psiN -> r_eff transformation
        implicit none

        type(profile_preprocessor_t) :: preprocessor
        character(256) :: equil_path
        logical :: equil_exists

        ! Determine equilibrium file path
        if (len_trim(equil_file) > 0) then
            equil_path = equil_file
        else
            equil_path = trim(input_profile_dir) // '/equil_r_q_psi.dat'
        end if

        inquire(file=trim(equil_path), exist=equil_exists)
        if (.not. equil_exists) then
            write(*,*) 'ERROR: Equilibrium file not found: ', trim(equil_path)
            write(*,*) 'Cannot perform sqrt_psiN -> r_eff transformation'
            stop 1
        end if

        write(*,*) 'Running profile preprocessing (sqrt_psiN -> r_eff)'
        write(*,*) '  Input directory: ', trim(input_profile_dir)
        write(*,*) '  Equilibrium: ', trim(equil_path)
        write(*,*) '  Output directory: ', trim(profile_location)

        ! Initialize preprocessor
        call preprocessor%init(trim(equil_path))

        ! Process each profile
        call process_profile(preprocessor, 'n')
        call process_profile(preprocessor, 'Te')
        call process_profile(preprocessor, 'Ti')
        call process_profile(preprocessor, 'Vz')

        write(*,*) 'Profile preprocessing complete'

    end subroutine run_preprocessing

    subroutine process_profile(preprocessor, profile_name)
        !> Process a single profile through the preprocessor
        type(profile_preprocessor_t), intent(inout) :: preprocessor
        character(*), intent(in) :: profile_name

        character(256) :: input_file, output_file
        logical :: file_exists
        real(dp), allocatable :: psi_in(:), data_in(:), r_out(:), data_out(:)
        integer :: npts, i, iunit

        ! Check for input file
        input_file = trim(input_profile_dir) // '/' // trim(profile_name) // '_of_psiN.dat'
        inquire(file=trim(input_file), exist=file_exists)

        if (.not. file_exists) then
            if (trim(profile_name) == 'Vz') then
                write(*,*) '  Note: ', trim(profile_name), '_of_psiN.dat not found, skipping'
                return
            else
                write(*,*) 'ERROR: Required profile not found: ', trim(input_file)
                stop 1
            end if
        end if

        ! Read input profile
        call read_profile_data(input_file, psi_in, data_in, npts)

        ! Transform coordinates
        allocate(r_out(npts), data_out(npts))
        do i = 1, npts
            r_out(i) = preprocessor%psi_to_r(psi_in(i))
            data_out(i) = data_in(i)
        end do

        ! Write output
        output_file = trim(profile_location) // '/' // trim(profile_name) // '.dat'
        open(newunit=iunit, file=trim(output_file), status='replace', action='write')
        do i = 1, npts
            write(iunit, '(2E20.12)') r_out(i), data_out(i)
        end do
        close(iunit)

        write(*,*) '  Processed: ', trim(profile_name), ' (', npts, ' points)'

        deallocate(psi_in, data_in, r_out, data_out)

    end subroutine process_profile
```

**Step 3: Build to verify**

Run: `make KIM`
Expected: BUILD SUCCESS

**Step 4: Commit**

```bash
git add KIM/src/background_equilibrium/profile_input_m.f90
git commit -m "feat(KIM): implement profile preprocessing transformation"
```

---

## Task 10: Integrate into KIM_init.f90

**Files:**
- Modify: `KIM/src/general/KIM_init.f90`

**Step 1: Add module use**

Add to use statements:
```fortran
    use profile_input_m, only: prepare_profiles
```

**Step 2: Call prepare_profiles before read_profiles**

After the HDF5 initialization block and before `call allocate_plasma`, add:
```fortran
    call prepare_profiles()
```

The modified KIM_init should look like:
```fortran
subroutine kim_init

    use species_m, only: read_profiles, allocate_plasma, init_plasma, plasma
    use IO_collection_m, only: initialize_hdf5_output, write_KIM_namelist_to_hdf5
    use config_m, only: hdf5_output
    use profile_input_m, only: prepare_profiles

    implicit none

    call read_config

    if (hdf5_output) then
        call initialize_hdf5_output()
        call write_KIM_namelist_to_hdf5()
    end if

    call prepare_profiles()

    call allocate_plasma
    call init_plasma(plasma)
    call read_profiles()

end subroutine
```

**Step 3: Build to verify**

Run: `make KIM`
Expected: BUILD SUCCESS

**Step 4: Run tests**

Run: `make test`
Expected: All tests PASS

**Step 5: Commit**

```bash
git add KIM/src/general/KIM_init.f90
git commit -m "feat(KIM): integrate prepare_profiles into KIM initialization"
```

---

## Task 11: Add Integration Test

**Files:**
- Create: `KIM/tests/test_profile_input_integration.f90`
- Create: `KIM/tests/test_data/` directory with test profiles

**Step 1: Create test data directory and files**

Create minimal test profile files for integration testing.

**Step 2: Create integration test**

Test that:
1. Auto-detection works correctly
2. Er calculation produces valid output
3. btor/R0 validation runs without error

**Step 3: Add to CMakeLists.txt**

**Step 4: Run test**

Run: `make test`
Expected: All tests PASS

**Step 5: Commit**

```bash
git add KIM/tests/
git commit -m "test(KIM): add profile_input integration tests"
```

---

## Task 12: Final Verification and Cleanup

**Step 1: Run full test suite**

Run: `make clean && make all && make test`
Expected: BUILD SUCCESS, All tests PASS

**Step 2: Verify no compiler warnings**

Check build output for warnings related to profile_input_m.

**Step 3: Update CLAUDE.md if needed**

Document new namelist group KIM_profiles.

**Step 4: Final commit**

```bash
git add -A
git commit -m "docs(KIM): document KIM_profiles namelist configuration"
```

---

## Summary

| Task | Description | Estimated Complexity |
|------|-------------|---------------------|
| 1 | Add namelist variables | Simple |
| 2 | Add namelist reading | Simple |
| 3 | Update default namelist | Simple |
| 4 | Create module skeleton | Medium |
| 5 | Update CMakeLists.txt | Simple |
| 6 | Implement coord detection | Medium |
| 7 | Implement btor/R0 validation | Medium |
| 8 | Implement Er calculation | Complex |
| 9 | Implement preprocessing | Complex |
| 10 | Integrate into KIM_init | Simple |
| 11 | Add integration tests | Medium |
| 12 | Final verification | Simple |
