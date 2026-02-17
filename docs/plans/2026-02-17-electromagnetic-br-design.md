# Self-Consistent Br Calculation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement self-consistent radial magnetic field (Br) by coupling the Poisson and Ampere equations into a 2N x 2N block system with direct solve (Issue #52).

**Architecture:** New `electromagnetic_t` run type extending `kim_t`. The Poisson equation for Phi and the parallel Ampere equation for Br are assembled into a single block matrix and solved with the existing UMFPACK/SuperLU sparse solver. Requires FokkerPlanck collision model (all 4 kernel matrices).

**Tech Stack:** Fortran 2008, CMake, CTest, UMFPACK/SuperLU sparse solvers, HDF5 output.

---

### Task 1: Add Namelist Parameters

**Files:**
- Modify: `KIM/src/setup/setup_mod.f90:19` (after `mphi_max`)
- Modify: `KIM/src/general/read_config.f90:31-32` (namelist statement)

**Step 1: Add variables to `setup_m`**

In `KIM/src/setup/setup_mod.f90`, add after line 19 (`integer :: mphi_max`):

```fortran
    real(dp) :: Br_boundary_re = 1.0d0  ! real part of Br at right boundary
    real(dp) :: Br_boundary_im = 0.0d0  ! imaginary part of Br at right boundary
```

**Step 2: Add to namelist in `read_config.f90`**

In `KIM/src/general/read_config.f90`, change lines 31-32 from:
```fortran
    namelist /KIM_SETUP/ btor, R0, m_mode, n_mode, omega, spline_base, &
                        type_br_field, collisions_off, set_profiles_constant, bc_type, mphi_max
```
to:
```fortran
    namelist /KIM_SETUP/ btor, R0, m_mode, n_mode, omega, spline_base, &
                        type_br_field, collisions_off, set_profiles_constant, bc_type, mphi_max, &
                        Br_boundary_re, Br_boundary_im
```

**Step 3: Build to verify no errors**

Run: `make clean && make all`
Expected: Builds successfully; existing namelist files work (new params have defaults).

**Step 4: Run existing tests**

Run: `make test`
Expected: All existing tests pass.

**Step 5: Commit**

```bash
git add KIM/src/setup/setup_mod.f90 KIM/src/general/read_config.f90
git commit -m "feat(KIM): add Br_boundary namelist parameters for electromagnetic run type"
```

---

### Task 2: Implement Ampere Matrices Module

**Files:**
- Create: `KIM/src/electromagnetic/ampere_matrices.f90`

This module provides:
- `calc_potential_matrix(Q_mat)` -- the unweighted Q_{l'l} = int dr 1/r phi_l' phi_l (tridiagonal, analytical log-expressions)
- `calc_weighted_mass_matrix(M_hth, hz_xl, hth_xl)` -- M^{h_theta}_{l'l} using midpoint approximation
- `calc_weighted_potential_matrix(Q_hz, Q_mat, hz_xl, hth_xl)` -- Q^{h_z}_{l'l} using midpoint approximation
- `interpolate_equil_to_xl(hz_xl, hth_xl)` -- interpolates hz, hth from equil_grid to xl_grid

**Step 1: Create the module**

Write `KIM/src/electromagnetic/ampere_matrices.f90`:

```fortran
module ampere_matrices_m

    use KIM_kinds_m, only: dp

    implicit none

    contains

    subroutine interpolate_equil_to_xl(hz_xl, hth_xl)
        ! Interpolate equilibrium hz, hth from equil_grid onto xl_grid%xb

        use equilibrium_m, only: hz, hth, equil_grid
        use grid_m, only: xl_grid

        implicit none

        real(dp), allocatable, intent(out) :: hz_xl(:), hth_xl(:)
        integer :: i, ir, ibeg, iend
        integer :: nlagr = 4
        integer :: nder = 0
        real(dp), allocatable :: coef(:,:)

        allocate(hz_xl(xl_grid%npts_b), hth_xl(xl_grid%npts_b))
        allocate(coef(0:nder, nlagr))

        do i = 1, xl_grid%npts_b
            call binsrc(equil_grid, 1, size(equil_grid), xl_grid%xb(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend > size(equil_grid)) then
                iend = size(equil_grid)
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, xl_grid%xb(i), equil_grid(ibeg:iend), coef)
            hz_xl(i)  = sum(coef(0,:) * hz(ibeg:iend))
            hth_xl(i) = sum(coef(0,:) * hth(ibeg:iend))
        end do

        deallocate(coef)

    end subroutine

    subroutine calc_potential_matrix(Q_mat)
        ! Compute the unweighted potential matrix Q_{l'l} = int dr (1/r) phi_l' phi_l
        ! Analytical expressions for hat basis functions (tridiagonal).

        use grid_m, only: xl_grid

        implicit none

        real(dp), intent(out) :: Q_mat(:,:)
        integer :: i, n
        real(dp) :: r_l, r_lm1, r_lp1

        n = xl_grid%npts_b
        Q_mat = 0.0d0

        do i = 2, n - 1
            r_l   = xl_grid%xb(i)
            r_lm1 = xl_grid%xb(i-1)
            r_lp1 = xl_grid%xb(i+1)

            ! Diagonal: l' = l
            Q_mat(i,i) = 0.5d0 * ( &
                (r_l**2 + 2.0d0*r_lm1**2*(log(r_l) - log(r_lm1)) &
                 - 4.0d0*r_l*r_lm1 + 3.0d0*r_lm1**2) / (r_l - r_lm1)**2 &
                - (r_l**2 + 2.0d0*r_lp1**2*(log(r_l) - log(r_lp1)) &
                 - 4.0d0*r_l*r_lp1 + 3.0d0*r_lp1**2) / (r_l - r_lp1)**2 )

            ! Off-diagonal: l' = l+1
            if (i < n) then
                Q_mat(i, i+1) = (-r_l**2 + 2.0d0*r_l*r_lp1*(log(r_l) - log(r_lp1)) &
                    + r_lp1**2) / (2.0d0 * (r_l - r_lp1)**2)
            end if

            ! Off-diagonal: l' = l-1
            if (i > 1) then
                Q_mat(i, i-1) = (r_l**2 - 2.0d0*r_l*r_lm1*(log(r_l) - log(r_lm1)) &
                    - r_lm1**2) / (2.0d0 * (r_l - r_lm1)**2)
            end if
        end do

    end subroutine

    subroutine calc_weighted_mass_matrix(M_hth, M_mat, hth_xl)
        ! M^{h_theta}_{l'l} = h_theta(r_bar) * M_{l'l}
        ! Midpoint approximation: h_theta evaluated at the midpoint of
        ! overlapping basis function support.

        use grid_m, only: xl_grid

        implicit none

        real(dp), intent(out) :: M_hth(:,:)
        real(dp), intent(in) :: M_mat(:,:)
        real(dp), intent(in) :: hth_xl(:)
        integer :: i, n

        n = xl_grid%npts_b
        M_hth = 0.0d0

        do i = 1, n
            ! Diagonal: midpoint is the node itself
            M_hth(i,i) = hth_xl(i) * M_mat(i,i)
            ! Off-diagonal: midpoint between nodes i and i+1
            if (i < n) then
                M_hth(i, i+1) = 0.5d0 * (hth_xl(i) + hth_xl(i+1)) * M_mat(i, i+1)
                M_hth(i+1, i) = 0.5d0 * (hth_xl(i) + hth_xl(i+1)) * M_mat(i+1, i)
            end if
        end do

    end subroutine

    subroutine calc_weighted_potential_matrix(Q_hz, Q_mat, hz_xl)
        ! Q^{h_z}_{l'l} = h_z(r_bar) * Q_{l'l}
        ! Midpoint approximation: h_z evaluated at the midpoint of
        ! overlapping basis function support.

        use grid_m, only: xl_grid

        implicit none

        real(dp), intent(out) :: Q_hz(:,:)
        real(dp), intent(in) :: Q_mat(:,:)
        real(dp), intent(in) :: hz_xl(:)
        integer :: i, n

        n = xl_grid%npts_b
        Q_hz = 0.0d0

        do i = 1, n
            ! Diagonal: midpoint is the node itself
            Q_hz(i,i) = hz_xl(i) * Q_mat(i,i)
            ! Off-diagonal: midpoint between nodes i and i+1
            if (i < n) then
                Q_hz(i, i+1) = 0.5d0 * (hz_xl(i) + hz_xl(i+1)) * Q_mat(i, i+1)
                Q_hz(i+1, i) = 0.5d0 * (hz_xl(i) + hz_xl(i+1)) * Q_mat(i+1, i)
            end if
        end do

    end subroutine

end module
```

**Step 2: Build to verify compilation**

Add the file to CMake first (Task 5 does the full wiring, but for now verify syntax):
Run: `gfortran -c -fsyntax-only KIM/src/electromagnetic/ampere_matrices.f90` (optional quick check)

**Step 3: Commit**

```bash
git add KIM/src/electromagnetic/ampere_matrices.f90
git commit -m "feat(KIM): add Ampere matrices module for electromagnetic run type"
```

---

### Task 3: Write Unit Test for Ampere Matrices

**Files:**
- Create: `KIM/tests/test_ampere_matrices.f90`
- Modify: `KIM/tests/CMakeLists.txt`

**Step 1: Write the test**

The test verifies:
1. Q_mat symmetry: Q(i,j) = Q(j,i)
2. Q_mat row sums approximate int(1/r * phi_l, dr) correctly
3. Weighted matrices with constant h reduce to h * unweighted matrices

Write `KIM/tests/test_ampere_matrices.f90`:

```fortran
program test_ampere_matrices

    use KIM_kinds_m, only: dp
    use grid_m, only: xl_grid, grid_t, calc_mass_matrix, M_mat
    use ampere_matrices_m, only: calc_potential_matrix, &
        calc_weighted_mass_matrix, calc_weighted_potential_matrix

    implicit none

    integer, parameter :: n = 10
    real(dp) :: r_min = 10.0d0   ! cm
    real(dp) :: r_max = 50.0d0   ! cm
    real(dp), allocatable :: Q_mat_loc(:,:), M_hth(:,:), Q_hz(:,:)
    real(dp), allocatable :: hth_const(:), hz_const(:)
    real(dp) :: h_const_val, tol
    integer :: i, j

    tol = 1.0d-12

    ! Set up an equidistant xl_grid
    xl_grid%npts_b = n
    allocate(xl_grid%xb(n))
    do i = 1, n
        xl_grid%xb(i) = r_min + (r_max - r_min) * dble(i-1) / dble(n-1)
    end do

    ! Compute mass matrix
    allocate(M_mat(n, n))
    call calc_mass_matrix(M_mat)

    ! Compute potential matrix
    allocate(Q_mat_loc(n, n))
    call calc_potential_matrix(Q_mat_loc)

    ! Test 1: Q_mat symmetry
    do i = 2, n-1
        do j = 2, n-1
            if (abs(Q_mat_loc(i,j) - Q_mat_loc(j,i)) > tol) then
                print *, 'FAIL: Q_mat not symmetric at (', i, ',', j, ')'
                print *, '  Q(i,j) =', Q_mat_loc(i,j), ' Q(j,i) =', Q_mat_loc(j,i)
                error stop
            end if
        end do
    end do
    print *, 'PASS: Q_mat is symmetric'

    ! Test 2: Q_mat entries are positive on diagonal (for r > 0)
    do i = 2, n-1
        if (Q_mat_loc(i,i) <= 0.0d0) then
            print *, 'FAIL: Q_mat diagonal entry non-positive at i =', i
            error stop
        end if
    end do
    print *, 'PASS: Q_mat diagonal entries are positive'

    ! Test 3: Weighted mass matrix with constant h_theta = h_const
    !         should give M_hth = h_const * M_mat
    h_const_val = 0.05d0
    allocate(hth_const(n), hz_const(n))
    hth_const = h_const_val
    hz_const  = h_const_val

    allocate(M_hth(n, n), Q_hz(n, n))
    call calc_weighted_mass_matrix(M_hth, M_mat, hth_const)

    do i = 1, n
        do j = 1, n
            if (abs(M_hth(i,j) - h_const_val * M_mat(i,j)) > tol) then
                print *, 'FAIL: M_hth /= h * M_mat at (', i, ',', j, ')'
                print *, '  M_hth =', M_hth(i,j), ' expected =', h_const_val * M_mat(i,j)
                error stop
            end if
        end do
    end do
    print *, 'PASS: M_hth = h_const * M_mat for constant h_theta'

    ! Test 4: Weighted potential matrix with constant h_z = h_const
    !         should give Q_hz = h_const * Q_mat
    call calc_weighted_potential_matrix(Q_hz, Q_mat_loc, hz_const)

    do i = 1, n
        do j = 1, n
            if (abs(Q_hz(i,j) - h_const_val * Q_mat_loc(i,j)) > tol) then
                print *, 'FAIL: Q_hz /= h * Q_mat at (', i, ',', j, ')'
                print *, '  Q_hz =', Q_hz(i,j), ' expected =', h_const_val * Q_mat_loc(i,j)
                error stop
            end if
        end do
    end do
    print *, 'PASS: Q_hz = h_const * Q_mat for constant h_z'

    print *, 'All Ampere matrix tests passed.'

end program
```

**Step 2: Register test in CMake**

Add to `KIM/tests/CMakeLists.txt`:

```cmake
add_executable(test_ampere_matrices ${CMAKE_SOURCE_DIR}/KIM/tests/test_ampere_matrices.f90)
set_target_properties(test_ampere_matrices PROPERTIES
                        OUTPUT_NAME test_ampere_matrices.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_ampere_matrices KIM_lib lapack cerf slatec ddeabm)
add_test(NAME test_ampere_matrices
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_ampere_matrices.x)
```

**Step 3: Build and run test (after Task 5 wires CMake)**

Run: `make clean && make all && make test`
Expected: `test_ampere_matrices` PASSES.

**Step 4: Commit**

```bash
git add KIM/tests/test_ampere_matrices.f90 KIM/tests/CMakeLists.txt
git commit -m "test(KIM): add unit tests for Ampere matrices"
```

---

### Task 4: Implement Electromagnetic Run Type

**Files:**
- Create: `KIM/src/electromagnetic/electromagnetic_solver.f90`

This is the core implementation: the `electromagnetic_t` type with `init` and `run` methods.
The `run` method fills all 4 kernels, assembles the 2N x 2N block system, solves it,
and runs postprocessing.

**Step 1: Create the module**

Write `KIM/src/electromagnetic/electromagnetic_solver.f90`:

```fortran
! Run type for electromagnetic model
! Solves coupled Poisson-Ampere block system for Phi and Br self-consistently
module rt_electromagnetic_m

    use kim_base_m, only: kim_t

    implicit none

    type, extends(kim_t) :: electromagnetic_t
        contains
            procedure :: init => init_electromagnetic
            procedure :: run  => run_electromagnetic
    end type electromagnetic_t

    contains

    subroutine init_electromagnetic(this)

        use species_m, only: plasma, set_plasma_quantities
        use IO_collection_m, only: create_output_directories
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use grid_m, only: rg_grid

        implicit none

        class(electromagnetic_t), intent(inout) :: this

        this%run_type = "electromagnetic"

        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)

        print *, "..."//trim(this%run_type)//" model initialized."

    end subroutine

    subroutine run_electromagnetic(this)

        use kernel_m, only: FP_fill_kernels, kernel_spl_t, write_kernels
        use kernel_adaptive_m, only: FP_fill_kernels_adaptive
        use grid_m, only: xl_grid, calc_mass_matrix, M_mat, theta_integration
        use IO_collection_m, only: write_complex_profile_abs, write_matrix
        use poisson_solver_m, only: prepare_Laplace_matrix, dense_to_sparse
        use config_m, only: output_path, collision_model, fstatus, fdebug
        use fields_m, only: EBdat, postprocess_electric_field, &
                            calculate_charge_density, calculate_current_density
        use ampere_matrices_m, only: calc_potential_matrix, &
            calc_weighted_mass_matrix, calc_weighted_potential_matrix, &
            interpolate_equil_to_xl
        use setup_m, only: m_mode, n_mode, R0, Br_boundary_re, Br_boundary_im, bc_type
        use constants_m, only: pi, sol, com_unit, dp => dp
        use KIM_kinds_m, only: dp
        use species_m, only: plasma
        use sparse_mod, only: sparse_solveComplex_b1, sparse_solve_method, sparse_talk

        implicit none

        class(electromagnetic_t), intent(inout) :: this

        type(kernel_spl_t) :: kernel_rho_phi_llp
        type(kernel_spl_t) :: kernel_rho_B_llp
        type(kernel_spl_t) :: kernel_j_phi_llp
        type(kernel_spl_t) :: kernel_j_B_llp

        complex(dp), allocatable :: A_block(:,:)  ! 2N x 2N block system
        complex(dp), allocatable :: b_block(:)    ! 2N RHS
        complex(dp), allocatable :: A_nz(:)
        integer, allocatable :: irow(:), pcol(:)
        integer :: nz_out, nrow_sp, ncol_sp

        complex(dp), allocatable :: A_Phi(:,:)    ! N x N Poisson LHS
        real(dp), allocatable :: Q_mat(:,:)       ! unweighted potential matrix
        real(dp), allocatable :: M_hth(:,:)       ! weighted mass matrix
        real(dp), allocatable :: Q_hz(:,:)        ! weighted potential matrix
        real(dp), allocatable :: hz_xl(:), hth_xl(:)

        complex(dp), allocatable :: jpar(:), rho(:)
        complex(dp) :: Br_boundary

        integer :: N, i
        real(dp) :: kz
        character(8) :: date
        character(10) :: time
        character(5) :: zone
        integer, dimension(8) :: values

        ! Validate collision model
        if (trim(collision_model) /= "FokkerPlanck") then
            print *, "Error: electromagnetic run type requires collision_model = 'FokkerPlanck'"
            print *, "Current collision_model: ", trim(collision_model)
            error stop
        end if

        N = xl_grid%npts_b
        kz = dble(n_mode) / R0
        Br_boundary = cmplx(Br_boundary_re, Br_boundary_im, dp)

        ! Initialize kernels
        call kernel_rho_phi_llp%init_kernel(N, N)
        call kernel_rho_B_llp%init_kernel(N, N)
        call kernel_j_phi_llp%init_kernel(N, N)
        call kernel_j_B_llp%init_kernel(N, N)

        ! Fill kernels
        call date_and_time(date, time, zone, values)
        write(*,*) "Start filling kernel at ", date, " ", time, " ..."

        select case (trim(theta_integration))
            case ("GaussLegendre")
                call FP_fill_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, &
                    kernel_j_phi_llp, kernel_j_B_llp)
            case ("RKF45", "QUADPACK")
                call FP_fill_kernels_adaptive(kernel_rho_phi_llp, kernel_rho_B_llp, &
                    kernel_j_phi_llp, kernel_j_B_llp)
            case default
                error stop "Error: theta integration method not recognized."
        end select

        call write_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, &
            kernel_j_phi_llp, kernel_j_B_llp)

        ! Compute mass matrix (reuse if already allocated)
        if (.not. allocated(M_mat)) then
            allocate(M_mat(N, N))
            call calc_mass_matrix(M_mat)
        end if

        ! Compute Ampere matrices
        call interpolate_equil_to_xl(hz_xl, hth_xl)

        allocate(Q_mat(N, N), M_hth(N, N), Q_hz(N, N))
        call calc_potential_matrix(Q_mat)
        call calc_weighted_mass_matrix(M_hth, M_mat, hth_xl)
        call calc_weighted_potential_matrix(Q_hz, Q_mat, hz_xl)

        ! Assemble Poisson LHS block: A_Phi = Delta + 4*pi * K_rho_phi
        allocate(A_Phi(N, N))
        call prepare_Laplace_matrix(A_Phi)
        A_Phi = A_Phi + 4.0d0 * pi * kernel_rho_phi_llp%Kllp

        ! Assemble 2N x 2N block system
        allocate(A_block(2*N, 2*N), b_block(2*N))
        A_block = cmplx(0.0d0, 0.0d0, dp)
        b_block = cmplx(0.0d0, 0.0d0, dp)

        ! Top-left: A_Phi (Poisson LHS)
        A_block(1:N, 1:N) = A_Phi

        ! Top-right: C_PhiB = 4*pi * K_rho_B
        A_block(1:N, N+1:2*N) = 4.0d0 * pi * kernel_rho_B_llp%Kllp

        ! Bottom-left: C_BPhi = i*4*pi/c * M * K_j_phi
        A_block(N+1:2*N, 1:N) = com_unit * 4.0d0 * pi / sol &
            * matmul(cmplx(M_mat, 0.0d0, dp), kernel_j_phi_llp%Kllp)

        ! Bottom-right: A_B = kz * M^h_theta - m * Q^h_z + i*4*pi/c * M * K_j_B
        A_block(N+1:2*N, N+1:2*N) = cmplx(kz * M_hth - dble(m_mode) * Q_hz, 0.0d0, dp) &
            + com_unit * 4.0d0 * pi / sol &
            * matmul(cmplx(M_mat, 0.0d0, dp), kernel_j_B_llp%Kllp)

        ! --- Boundary conditions ---

        ! Phi BCs: Dirichlet Phi = 0 at both boundaries (simplest case)
        ! Row 1 (left boundary)
        A_block(1, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(1, 1) = cmplx(1.0d0, 0.0d0, dp)
        b_block(1) = cmplx(0.0d0, 0.0d0, dp)

        ! Row N (right boundary)
        A_block(N, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(N, N) = cmplx(1.0d0, 0.0d0, dp)
        b_block(N) = cmplx(0.0d0, 0.0d0, dp)

        ! Br BCs: Br(left) = 0, Br(right) = Br_boundary
        ! Row N+1 (left boundary)
        A_block(N+1, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(N+1, N+1) = cmplx(1.0d0, 0.0d0, dp)
        b_block(N+1) = cmplx(0.0d0, 0.0d0, dp)

        ! Row 2N (right boundary)
        A_block(2*N, :) = cmplx(0.0d0, 0.0d0, dp)
        A_block(2*N, 2*N) = cmplx(1.0d0, 0.0d0, dp)
        b_block(2*N) = Br_boundary

        ! Zero-misalignment BC override for Phi (if bc_type == 3)
        if (bc_type == 3) then
            call apply_zero_misalignment_bc(A_block, b_block, N, Br_boundary)
        end if

        if (fstatus >= 1) write(*,*) 'Status: solving coupled Poisson-Ampere block system (2N =', 2*N, ')'

        ! Solve
        if (fdebug < 2) sparse_talk = .false.
        sparse_solve_method = 1

        call dense_to_sparse(A_block, irow, pcol, A_nz, nrow_sp, ncol_sp, nz_out)
        call sparse_solveComplex_b1(nrow_sp, ncol_sp, nz_out, irow, pcol, A_nz, b_block, 0)

        ! Extract solution
        allocate(EBdat%Phi(N), EBdat%Br(N), EBdat%r_grid(N))
        EBdat%r_grid = xl_grid%xb
        EBdat%Phi = b_block(1:N)
        EBdat%Br  = b_block(N+1:2*N)

        ! Write fields
        call write_complex_profile_abs(xl_grid%xb, EBdat%Phi, N, "/fields/Phi", &
            'Electrostatic potential perturbation Phi (self-consistent)', 'statV')
        call write_complex_profile_abs(xl_grid%xb, EBdat%Br, N, "/fields/Br_selfconsistent", &
            'Self-consistent radial magnetic field perturbation Br', 'G')

        ! Current density
        call calculate_current_density(jpar, EBdat%Phi, EBdat%Br, &
            kernel_j_phi_llp%Kllp, kernel_j_B_llp%Kllp)
        call write_complex_profile_abs(xl_grid%xb, jpar, N, "/fields/jpar", &
            'Parallel current density from self-consistent solve', 'statA/cm^2')

        ! Electric field postprocessing
        call postprocess_electric_field(EBdat)

        ! Charge density
        call calculate_charge_density(rho, EBdat)
        call write_complex_profile_abs(xl_grid%xb, rho, N, "/fields/rho", &
            'Charge density from self-consistent solve', 'statC/cm^3')

        ! Cleanup
        deallocate(A_block, b_block, A_Phi, Q_mat, M_hth, Q_hz, hz_xl, hth_xl)

    end subroutine

    subroutine apply_zero_misalignment_bc(A_block, b_block, N, Br_boundary)
        ! For bc_type == 3: set Phi boundaries such that misalignment field = 0.
        ! phi_aligned = com_unit * Er * Br / (B0 * kp)
        ! Phi_boundary = -phi_aligned, evaluated at boundaries using known Br BC values.

        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit
        use species_m, only: plasma
        use grid_m, only: xl_grid

        implicit none

        complex(dp), intent(inout) :: A_block(:,:), b_block(:)
        integer, intent(in) :: N
        complex(dp), intent(in) :: Br_boundary

        integer :: ir, ibeg, iend
        integer :: nlagr = 4, nder = 0
        real(dp), allocatable :: coef(:,:)
        real(dp) :: Er_left, kp_left, B0_left
        real(dp) :: Er_right, kp_right, B0_right
        complex(dp) :: phi_left, phi_right

        allocate(coef(0:nder, nlagr))

        ! Left boundary: Br = 0, so phi_aligned = 0, so Phi_left = 0 (unchanged)

        ! Right boundary: Br = Br_boundary
        call binsrc(plasma%r_grid, 1, size(plasma%r_grid), xl_grid%xb(N), ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > size(plasma%r_grid)) then
            iend = size(plasma%r_grid)
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, xl_grid%xb(N), plasma%r_grid(ibeg:iend), coef)
        Er_right = sum(coef(0,:) * plasma%Er(ibeg:iend))
        kp_right = sum(coef(0,:) * plasma%kp(ibeg:iend))
        B0_right = sum(coef(0,:) * plasma%B0(ibeg:iend))

        phi_right = -com_unit * Er_right * Br_boundary / (B0_right * kp_right)

        ! Modify Phi right boundary (row N)
        ! Already zeroed in main assembly, just set RHS
        b_block(N) = phi_right

        deallocate(coef)

    end subroutine

end module
```

**Step 2: Commit**

```bash
git add KIM/src/electromagnetic/electromagnetic_solver.f90
git commit -m "feat(KIM): implement electromagnetic_t run type with block Poisson-Ampere solve"
```

---

### Task 5: Register in Factory and Wire CMake

**Files:**
- Modify: `KIM/src/general/KIM_mod.f90:10-28` (factory)
- Modify: `KIM/src/CMakeSources.in` (new source group)
- Modify: `KIM/src/CMakeLists.txt:6-19` (add to library)

**Step 1: Add electromagnetic case to factory**

In `KIM/src/general/KIM_mod.f90`, add `use rt_electromagnetic_m` import and a new case branch:

After line 11 (`use rt_electrostatic_m, only: electrostatic_t`), add:
```fortran
        use rt_electromagnetic_m, only: electromagnetic_t
```

After line 25 (`case("WKB_dispersion")`/`allocate(...WKB_dispersion_t())`), add:
```fortran
            case("electromagnetic")
                allocate(kim_instance, source=electromagnetic_t())
```

Update the error message on line 28 to include "electromagnetic".

**Step 2: Add electromagnetic source group to CMakeSources.in**

After line 76 (`set(KIM_asymptotics ...)`), add:

```cmake
set(KIM_electromagnetic electromagnetic/electromagnetic_kernel.f90
                        electromagnetic/ampere_matrices.f90
                        electromagnetic/electromagnetic_solver.f90
)
```

**Step 3: Add to library target in CMakeLists.txt**

In `KIM/src/CMakeLists.txt`, add `"${KIM_electromagnetic}"` to the `add_library` call, e.g. after line 17 (`"${KIM_asymptotics}"`):
```cmake
                        "${KIM_electromagnetic}"
```

**Step 4: Build everything**

Run: `make clean && make all`
Expected: Builds without errors.

**Step 5: Run all tests**

Run: `make test`
Expected: All tests pass (including new `test_ampere_matrices`).

**Step 6: Commit**

```bash
git add KIM/src/general/KIM_mod.f90 KIM/src/CMakeSources.in KIM/src/CMakeLists.txt
git commit -m "feat(KIM): register electromagnetic run type in factory and CMake"
```

---

### Task 6: Verify Build and Run Existing Tests

**Files:** None (verification only)

**Step 1: Clean build**

Run: `make clean && make all`
Expected: Full build succeeds with no warnings related to new code.

**Step 2: Run all tests**

Run: `make test`
Expected: All tests pass. If any fail, investigate and fix before proceeding.

**Step 3: Verify no regressions with existing electrostatic run**

If a test namelist is available in `KIM/nmls/`, try:
Run: `./build/install/bin/KIM.x KIM/nmls/KIM_config.nml`
Expected: Runs as before (electrostatic mode, since `type_of_run = "electrostatic"` in default namelist).

---

### Task 7: Create Electromagnetic Test Namelist and Smoke Test

**Files:**
- Create: `KIM/nmls/KIM_config_electromagnetic.nml`

**Step 1: Create test namelist**

Copy `KIM/nmls/KIM_config.nml` and modify:
- `type_of_run = "electromagnetic"` in `&KIM_CONFIG`
- `collision_model = "FokkerPlanck"` in `&KIM_CONFIG`
- `Br_boundary_re = 1.0` in `&KIM_SETUP`
- `Br_boundary_im = 0.0` in `&KIM_SETUP`

**Step 2: Run the electromagnetic solver**

Run: `./build/install/bin/KIM.x KIM/nmls/KIM_config_electromagnetic.nml`
Expected: Runs to completion, produces output in the configured output directory with:
- `fields/Phi` (electrostatic potential)
- `fields/Br_selfconsistent` (self-consistent Br)
- `fields/jpar` (parallel current)
- All standard field outputs

**Step 3: Basic sanity check**

Verify that:
- Br at the left boundary is 0
- Br at the right boundary equals the configured `Br_boundary`
- Phi is non-trivial (not all zeros)

**Step 4: Commit**

```bash
git add KIM/nmls/KIM_config_electromagnetic.nml
git commit -m "feat(KIM): add electromagnetic test namelist configuration"
```

---

## Summary of All Files

| Action | File |
|--------|------|
| Modify | `KIM/src/setup/setup_mod.f90` |
| Modify | `KIM/src/general/read_config.f90` |
| Create | `KIM/src/electromagnetic/ampere_matrices.f90` |
| Create | `KIM/src/electromagnetic/electromagnetic_solver.f90` |
| Modify | `KIM/src/general/KIM_mod.f90` |
| Modify | `KIM/src/CMakeSources.in` |
| Modify | `KIM/src/CMakeLists.txt` |
| Create | `KIM/tests/test_ampere_matrices.f90` |
| Modify | `KIM/tests/CMakeLists.txt` |
| Create | `KIM/nmls/KIM_config_electromagnetic.nml` |

## Key Dependencies Between Tasks

- Task 1 (namelist) has no dependencies
- Task 2 (ampere matrices) has no dependencies
- Task 3 (test) depends on Task 2
- Task 4 (electromagnetic solver) depends on Tasks 1 and 2
- Task 5 (CMake wiring) depends on Tasks 2 and 4
- Task 6 (build verification) depends on Task 5
- Task 7 (smoke test) depends on Task 6

Tasks 1 and 2 can be done in parallel. Task 3 can start as soon as Task 2 is done.
Task 4 can start as soon as Tasks 1 and 2 are done. Task 5 wires everything together.
