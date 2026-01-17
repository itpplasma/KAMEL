module rt_WKB_dispersion_m

    !==========================================================================================
    ! module for running the WKB dispersion relation solver
    ! Algorithms for root finding are Muller and ZEAL, both with per-branch tracking
    ! Dispersion functions implemented: KIM and FLRE to second order
    !
    ! Configuration namelist is in KIM_config.nml named &wkb_dispersion
    !==========================================================================================

    use kim_base_m, only: kim_t

    implicit none

    type, extends(kim_t) :: WKB_dispersion_t
        contains
            procedure :: init => init_wkb_dispersion
            procedure :: run => run_wkb_dispersion
    end type WKB_dispersion_t

    contains

    subroutine init_wkb_dispersion(this)

        use species_m, only: plasma, set_plasma_quantities
        use IO_collection_m, only: create_output_directories
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use grid_m, only: rg_grid

        implicit none

        class(WKB_dispersion_t), intent(inout) :: this

        this%run_type = "WKB_dispersion"

        call create_output_directories
        call generate_grids
        call calculate_equil(.true.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)

        print *, "..."//trim(this%run_type)//" model initialized."

    end subroutine

    subroutine run_wkb_dispersion(this)

        use KIM_kinds_m, only: dp
        use muller_root_finding
        use grid_m, only: rg_grid
        use IO_collection_m, only: write_complex_profile_abs, ensure_dispersion_dir_exists
        use config_m, only: WKB_dispersion_solver, dispersion_output_path

        implicit none

        class(WKB_dispersion_t), intent(inout) :: this
        character(len=20) :: out_str
        complex(dp), allocatable :: found_roots(:,:)
        integer :: km, i

        print *, "Running "//trim(this%run_type)//" model ..."

        ! Ensure dispersion output directory exists
        call ensure_dispersion_dir_exists()

        km = 1 ! max number of roots to find
        allocate(found_roots(km, rg_grid%npts_b))

        select case (WKB_dispersion_solver)
            case ('Muller')
                call run_Muller_dispersion(km, found_roots)
                do i = 1, km
                    write(out_str, '(A,I0)') "kr", i
                    call write_complex_profile_abs(rg_grid%xb, found_roots(:,i), rg_grid%npts_b, &
                        trim(out_str), "WKB dispersion relation", "cm^{-2}", dispersion_output_path)
                end do
            case ('ZEAL')
                call run_ZEAL_dispersion()
            case default
                stop "Unknown WKB dispersion solver"
        end select

    end subroutine

    subroutine run_Muller_dispersion(km, found_roots)

        use KIM_kinds_m, only: dp
        use muller_root_finding
        use grid_m, only: rg_grid
        use IO_collection_m, only: write_complex_profile_abs, &
            track_root_branches, write_tracked_roots
        use config_m, only: WKB_dispersion_mode, dispersion_output_path
        use Function_Input_Module, only: f_ptr

        implicit none

        complex(dp), intent(out) :: found_roots(:,:)
        integer, intent(inout) :: km
        complex(dp), allocatable :: fnv(:,:)
        integer :: known, nmore, guess, maxits, nfound, ifail
        logical :: realrt
        real(dp) :: ep1, ep2
        integer :: i, j
        character(len=20) :: out_str

        ! Branch tracking arrays
        integer, allocatable :: n_roots_per_point(:)
        integer, allocatable :: all_multiplicities(:,:)
        integer, allocatable :: branch_id(:,:)
        integer :: n_branches

        known = 0 ! number of known roots
        nmore = 1 ! number of additional roots to find
        guess = 0 ! number of initial guesses provided
        maxits = 10000
        ep1 = 0.5d-6
        ep2 = 1.0d-7
        realrt = .false.

        allocate(fnv(km, rg_grid%npts_b))
        allocate(n_roots_per_point(rg_grid%npts_b))
        allocate(all_multiplicities(km, rg_grid%npts_b))
        allocate(branch_id(km, rg_grid%npts_b))

        found_roots = (0.0_dp, 0.0_dp)
        fnv = (0.0_dp, 0.0_dp)
        all_multiplicities = 1  ! Muller doesn't compute multiplicities, assume 1
        branch_id = 0

        if (nmore > km) then
            print *, "Error: nmore exceeds km in WKB dispersion root finding."
            stop
        end if

        select case (WKB_dispersion_mode)
            case ('KIM')
                f_ptr => dispersion_KIM
            case ('FLRE')
                f_ptr => dispersion_FLRE
            case default
                stop "Unknown WKB dispersion mode"
        end select

        do j = 1, rg_grid%npts_b

            ! initial guesses (optional)
            if (j > 1) then
                do i=1, min(nmore, km)
                    found_roots(i, j) = found_roots(i, j-1)
                end do
                guess = min(nmore, km)
            end if

            call roots(f_ptr, known, nmore, km, realrt, ep1, ep2, &
                    guess, maxits, found_roots(:, j), nfound, fnv(:, j), ifail)

            if (ifail == 1) then
                print *, 'Known roots are not actual roots at grid point ', j
            else if (ifail == 2) then
                print *, 'Bad input parameters at grid point ', j
            else if (ifail == 3) then
                print *, 'Maximum iterations exceeded without convergence at grid point ', j
            end if

            ! Store number of roots found at this grid point
            n_roots_per_point(j) = nfound

        end do

        ! Track root branches across grid points
        print *
        print *, '=== Muller Solver Summary ==='
        print *, 'Total grid points: ', rg_grid%npts_b
        print *, 'Total roots stored: ', sum(n_roots_per_point)
        print *
        print *, 'Tracking root branches...'

        call track_root_branches(found_roots, n_roots_per_point, rg_grid%npts_b, &
            km, branch_id, n_branches)

        print *, 'Found ', n_branches, ' distinct branches.'

        ! Write tracked branches to file
        call write_tracked_roots(rg_grid%xb, rg_grid%npts_b, &
            found_roots, fnv, all_multiplicities, &
            n_roots_per_point, branch_id, n_branches, &
            'muller_branches', 'Muller tracked dispersion branches', dispersion_output_path)

        print *, 'Tracked branches written to output files.'

        ! Clean up branch tracking arrays
        deallocate(n_roots_per_point, all_multiplicities, branch_id)

        contains

        subroutine dispersion_KIM(kr, f)

            use species_m, only: plasma
            use constants_m, only: com_unit
            use KIM_kinds_m, only: dp
            use grid_m, only: rg_grid

            implicit none

            complex(dp), intent(in)  :: kr
            complex(dp), intent(out) :: f
            integer :: sp
            integer :: ifail_bess, nz
            complex(dp) :: bess0, bessm1
            real(dp) :: bess0_re(2), bess0_im(2)
            complex(dp) :: kr_rho_squared

            bess0_re = 0.0d0
            bess0_im = 0.0d0
            f = (0.0d0, 0.0d0)

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr**2.0d0 * plasma%spec(sp)%rho_L(j)**2.0d0

                ! calculate modified Bessel functions I0 and I1 with complex argument
                call zbesi(real(kr_rho_squared, kind=dp), &
                           dimag(kr_rho_squared), &
                           0.0d0, &! initial order
                           2, &! KODE, 1= no exponential scaling, 2= with scaling
                           2, &! n  number of terms, with 2 goes to order +1 (BesselI is symmetric in order)
                           bess0_re, &
                           bess0_im, &
                           nz, & ! number of underflows set to zero
                           ifail_bess)

                if (ifail_bess == 1) then
                    print *, 'Warning: Bessel function calculation, I0 and I1, bad input - no computation'
                else if (ifail_bess == 2) then
                    print *, 'Warning: Bessel function calculation, I0 and I1, overflow occurred - no computation'
                else if (ifail_bess == 3) then
                    print *, 'Warning: Bessel function calculation, I0 and I1, precision warning - computation completed'
                else if (ifail_bess == 4) then
                    print *, 'Error: Bessel function calculation, I0 and I1, precision error - no computation'
                else if (ifail_bess == 5) then
                    print *, 'Error: Bessel function calculation, I0 and I1, algorithmic error - no computation'
                end if

                bess0 = cmplx(bess0_re(1), bess0_im(1), kind=dp) ! zero order Bessel function
                bessm1 = cmplx(bess0_re(2), bess0_im(2), kind=dp) ! first order Bessel function, is symmetric for integer order

                f = f + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 * ( -1.0d0 &
                    + com_unit * plasma%spec(sp)%vT(j)**2.0d0 * plasma%ks(j) &
                    / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * &
                    exp(abs(real(kr_rho_squared, kind=dp)) - kr_rho_squared) * & ! correct the scaling of the bessel functions
                    (&
                        plasma%spec(sp)%I00(j, 0) * (&
                            bess0 * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1.0d0 - kr_rho_squared)) &
                            + plasma%spec(sp)%A2(j) * kr_rho_squared * bessm1 &
                        )&
                         + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) * bess0 &
                    )&
                )

            end do

            ! more precise, includes first order derivative:
            ! f = kr**2.0d0 - com_unit / rg_grid%xb(j) * kr + plasma%kp(j)**2.0d0 - f
            
            ! neglect linear term, accurate enough
            f = kr**2.0d0 + plasma%kp(j)**2.0d0 - f

        end subroutine dispersion_KIM

        subroutine dispersion_FLRE(kr, f)
            ! dispersion equation for second order finite Larmor radius expansion

            use species_m, only: plasma
            use constants_m, only: com_unit
            use KIM_kinds_m, only: dp
            use grid_m, only: rg_grid

            implicit none

            complex(dp), intent(in)  :: kr
            complex(dp), intent(out) :: f
            integer :: sp
            complex(dp) :: kr_rho_squared

            f = (0.0d0, 0.0d0)

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr**2.0d0 * plasma%spec(sp)%rho_L(j)**2.0d0

                f = f + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 * ( -1.0d0 &
                    + com_unit * plasma%spec(sp)%vT(j)**2.0d0 * plasma%ks(j) &
                    / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * &
                    (&
                        plasma%spec(sp)%I00(j, 0) * ((1.0d0 - kr_rho_squared) * &
                            (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1.0d0 - kr_rho_squared)) &
                            !+ plasma%spec(sp)%A2(j) * kr_rho_squared**2.0d0 / 2.0d0 & ! 4th order in the FLRE
                        )&
                        + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) * (1.0d0 - kr_rho_squared) &
                    )&
                )
            end do

            ! more accureate, includes first order derivative:
            ! f = kr**2.0d0 - com_unit / rg_grid%xb(j) * kr + plasma%kp(j)**2.0d0 - f

            ! use the version that neglects the first order derivative term (is accurate enough)
            f = kr**2.0d0 + plasma%kp(j)**2.0d0 - f

        end subroutine dispersion_FLRE


    end subroutine


    subroutine test_dispersion_KIM(kr, f)
        ! Standalone version of dispersion_KIM for testing/comparison
        use species_m, only: plasma
        use constants_m, only: com_unit
        use KIM_kinds_m, only: dp
        use Function_Input_Module, only: rg_index

        implicit none

        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: f
        integer :: sp, j
        integer :: ifail_bess, nz
        complex(dp) :: bess0, bessm1
        real(dp) :: bess0_re(2), bess0_im(2)
        complex(dp) :: kr_rho_squared

        j = rg_index

        bess0_re = 0.0d0
        bess0_im = 0.0d0
        f = (0.0d0, 0.0d0)

        do sp = 0, plasma%n_species-1

            kr_rho_squared = kr**2.0d0 * plasma%spec(sp)%rho_L(j)**2.0d0

            call zbesi(real(kr_rho_squared, kind=dp), &
                       dimag(kr_rho_squared), &
                       0.0d0, 2, 2, bess0_re, bess0_im, nz, ifail_bess)

            bess0 = cmplx(bess0_re(1), bess0_im(1), kind=dp)
            bessm1 = cmplx(bess0_re(2), bess0_im(2), kind=dp)

            f = f + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 * ( -1.0d0 &
                + com_unit * plasma%spec(sp)%vT(j)**2.0d0 * plasma%ks(j) &
                / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * &
                exp(abs(real(kr_rho_squared, kind=dp)) - kr_rho_squared) * &
                (&
                    plasma%spec(sp)%I00(j, 0) * (&
                        bess0 * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1.0d0 - kr_rho_squared)) &
                        + plasma%spec(sp)%A2(j) * kr_rho_squared * bessm1 &
                    )&
                     + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) * bess0 &
                )&
            )

        end do

        f = kr**2.0d0 + plasma%kp(j)**2.0d0 - f

    end subroutine test_dispersion_KIM


    subroutine run_ZEAL_dispersion()
        !-----------------------------------------------------------------------
        ! Per-branch tracking ZEAL dispersion solver.
        ! Tracks up to WKB_max_tracked_branches branches, each with its own search window.
        ! 
        ! The search window of the grid point j>1 is based on the roots found at grid point j-1.
        ! It is centered around the roots with a specified half width.
        ! Additional branches within this window might be found. Set 'do_broad_search = true' in 
        ! KIM_config.nml to perform a broad search at certain grid points to find new branches.
        !-----------------------------------------------------------------------
        use Zeal_Module
        use KIM_kinds_m, only: dp
        use Function_Input_Module, only: rg_index, test_FDF_derivative, init_dispersion_mode
        use grid_m, only: rg_grid
        use IO_collection_m, only: write_tracked_roots
        use Zeal_Input_Module, only: set_zeal_search_region
        use config_m, only: WKB_max_tracked_branches, WKB_branch_search_halfwidth, &
            WKB_broad_search_halfwidth, WKB_broad_search_interval, WKB_root_tolerance, &
            WKB_dispersion_mode, dispersion_output_path

        implicit none

        integer :: j, i, b, k
        integer :: totalnumber, distinctnumber, refinednumber
        integer, dimension(:), pointer            :: multiplicities => null()
        logical, dimension(:), pointer            :: refinement_ok => null()
        complex(kind=dp), dimension(:), pointer   :: zeros => null()
        complex(kind=dp), dimension(:), pointer   :: fzeros => null()

        ! Per-branch tracking data
        complex(dp), allocatable :: branch_center(:)      ! Last known root position for each branch
        complex(dp), allocatable :: branch_fvalue(:)      ! Function value at branch root
        integer, allocatable :: branch_mult(:)            ! Multiplicity of each branch
        logical, allocatable :: branch_active(:)          ! Is branch currently being tracked?
        integer, allocatable :: branch_miss_count(:)      ! Consecutive misses (no root found)
        integer :: n_active_branches

        ! Storage arrays for all grid points (one root per branch per grid point)
        complex(dp), allocatable :: all_zeros(:,:)        ! (max_branches, n_grid)
        complex(dp), allocatable :: all_fzeros(:,:)
        integer, allocatable :: all_multiplicities(:,:)
        integer, allocatable :: all_branch_ids(:,:)       ! Branch ID for each stored root
        integer, allocatable :: n_roots_per_point(:)

        ! Temporary storage for ZEAL results
        complex(dp), allocatable :: temp_zeros(:)
        complex(dp), allocatable :: temp_fzeros(:)
        integer, allocatable :: temp_mults(:)
        logical, allocatable :: temp_valid(:)
        integer :: n_temp_valid

        ! Working variables
        real(dp) :: search_center_re, search_center_im
        real(dp) :: dist
        integer :: best_idx
        logical :: do_broad_search
        integer, parameter :: MAX_MISS_COUNT = 3  ! Deactivate branch after this many consecutive misses

        ! Allocate branch tracking arrays
        allocate(branch_center(WKB_max_tracked_branches))
        allocate(branch_fvalue(WKB_max_tracked_branches))
        allocate(branch_mult(WKB_max_tracked_branches))
        allocate(branch_active(WKB_max_tracked_branches))
        allocate(branch_miss_count(WKB_max_tracked_branches))

        ! Allocate storage for results
        allocate(all_zeros(WKB_max_tracked_branches, rg_grid%npts_b))
        allocate(all_fzeros(WKB_max_tracked_branches, rg_grid%npts_b))
        allocate(all_multiplicities(WKB_max_tracked_branches, rg_grid%npts_b))
        allocate(all_branch_ids(WKB_max_tracked_branches, rg_grid%npts_b))
        allocate(n_roots_per_point(rg_grid%npts_b))

        ! Temporary arrays for broad search
        allocate(temp_zeros(20))
        allocate(temp_fzeros(20))
        allocate(temp_mults(20))
        allocate(temp_valid(20))

        ! Initialize
        branch_center = (0.0_dp, 0.0_dp)
        branch_fvalue = (0.0_dp, 0.0_dp)
        branch_mult = 1
        branch_active = .false.
        branch_miss_count = 0
        n_active_branches = 0

        all_zeros = (0.0_dp, 0.0_dp)
        all_fzeros = (0.0_dp, 0.0_dp)
        all_multiplicities = 0
        all_branch_ids = 0
        n_roots_per_point = 0

        ! Initialize dispersion function pointer (KIM or FLRE)
        call init_dispersion_mode()

        ! Test derivative at first grid point
        rg_index = 1
        print *, 'Testing FDF derivative at grid point 1...'
        call test_FDF_derivative(cmplx(1.0d0, 0.5d0, dp))

        print *
        print *, '=== Per-Branch ZEAL Tracking ==='
        print *, 'Dispersion mode: ', trim(WKB_dispersion_mode)
        print *, 'Max tracked branches: ', WKB_max_tracked_branches
        print *, 'Branch search half-width: ', WKB_branch_search_halfwidth
        print *, 'Broad search half-width: ', WKB_broad_search_halfwidth
        print *

        ! Main loop over grid points
        do j = 1, rg_grid%npts_b

            print *
            print *, '--- Grid point ', j, ' at x = ', rg_grid%xb(j), ' ---'

            rg_index = j

            ! Decide whether to do broad search
            do_broad_search = (j == 1)  ! Always at first point
            if (WKB_broad_search_interval > 0 .and. j > 1) then
                if (mod(j-1, WKB_broad_search_interval) == 0) do_broad_search = .true.
            end if
            ! Also do broad search if we have no active branches
            if (n_active_branches == 0) do_broad_search = .true.

            !-------------------------------------------------------------------
            ! Step 1: Track existing branches with focused searches
            !-------------------------------------------------------------------
            do b = 1, WKB_max_tracked_branches
                if (.not. branch_active(b)) cycle

                ! Set small search window centered on this branch's last position
                search_center_re = real(branch_center(b), dp)
                search_center_im = aimag(branch_center(b))
                call set_zeal_search_region(search_center_re, search_center_im, WKB_branch_search_halfwidth)

                print *, '  Branch ', b, ': searching near (', search_center_re, ',', search_center_im, ')'

                ! Run ZEAL
                nullify(zeros, fzeros, multiplicities, refinement_ok)
                call zeal(totalnumber, distinctnumber, zeros, fzeros, multiplicities, refinednumber, refinement_ok)

                ! Find best valid root (smallest |f(z)|)
                best_idx = 0
                if (distinctnumber > 0 .and. associated(zeros)) then
                    do i = 1, distinctnumber
                        if (abs(fzeros(i)) < WKB_root_tolerance .or. refinement_ok(i)) then
                            if (best_idx == 0) then
                                best_idx = i
                            else if (abs(fzeros(i)) < abs(fzeros(best_idx))) then
                                best_idx = i
                            end if
                        end if
                    end do
                end if

                if (best_idx > 0) then
                    ! Found valid root - update branch
                    branch_center(b) = zeros(best_idx)
                    branch_fvalue(b) = fzeros(best_idx)
                    branch_mult(b) = multiplicities(best_idx)
                    branch_miss_count(b) = 0

                    ! Store in results
                    n_roots_per_point(j) = n_roots_per_point(j) + 1
                    all_zeros(n_roots_per_point(j), j) = zeros(best_idx)
                    all_fzeros(n_roots_per_point(j), j) = fzeros(best_idx)
                    all_multiplicities(n_roots_per_point(j), j) = multiplicities(best_idx)
                    all_branch_ids(n_roots_per_point(j), j) = b

                    print *, '    Found: z = (', real(zeros(best_idx)), ',', aimag(zeros(best_idx)), &
                        '), |f| = ', abs(fzeros(best_idx))
                else
                    ! No valid root found
                    branch_miss_count(b) = branch_miss_count(b) + 1
                    print *, '    No valid root found (miss count: ', branch_miss_count(b), ')'

                    if (branch_miss_count(b) >= MAX_MISS_COUNT) then
                        branch_active(b) = .false.
                        n_active_branches = n_active_branches - 1
                        print *, '    Branch ', b, ' deactivated after ', MAX_MISS_COUNT, ' consecutive misses'
                    end if
                end if
            end do

            !-------------------------------------------------------------------
            ! Step 2: Broad search to find new branches (if needed)
            !-------------------------------------------------------------------
            if (do_broad_search) then
                print *, '  Running broad search...'

                ! Use broad search window (centered at origin or average of active branches)
                if (n_active_branches > 0) then
                    search_center_re = 0.0_dp
                    search_center_im = 0.0_dp
                    k = 0
                    do b = 1, WKB_max_tracked_branches
                        if (branch_active(b)) then
                            search_center_re = search_center_re + real(branch_center(b), dp)
                            search_center_im = search_center_im + aimag(branch_center(b))
                            k = k + 1
                        end if
                    end do
                    if (k > 0) then
                        search_center_re = search_center_re / real(k, dp)
                        search_center_im = search_center_im / real(k, dp)
                    end if
                else
                    search_center_re = 0.0_dp
                    search_center_im = -1.0_dp  ! Default center
                end if

                call set_zeal_search_region(search_center_re, search_center_im, WKB_broad_search_halfwidth)
                print *, '  Broad search centered at (', search_center_re, ',', search_center_im, &
                    '), half-width = ', WKB_broad_search_halfwidth

                ! Run ZEAL
                nullify(zeros, fzeros, multiplicities, refinement_ok)
                call zeal(totalnumber, distinctnumber, zeros, fzeros, multiplicities, refinednumber, refinement_ok)

                print *, '  Broad search found ', distinctnumber, ' zeros'

                ! Collect valid roots
                n_temp_valid = 0
                if (distinctnumber > 0 .and. associated(zeros)) then
                    do i = 1, min(distinctnumber, 20)
                        if (abs(fzeros(i)) < WKB_root_tolerance .or. refinement_ok(i)) then
                            n_temp_valid = n_temp_valid + 1
                            temp_zeros(n_temp_valid) = zeros(i)
                            temp_fzeros(n_temp_valid) = fzeros(i)
                            temp_mults(n_temp_valid) = multiplicities(i)
                            temp_valid(n_temp_valid) = .true.
                        end if
                    end do
                end if

                print *, '  Valid roots from broad search: ', n_temp_valid

                ! Try to assign new roots to inactive branch slots
                do i = 1, n_temp_valid
                    if (.not. temp_valid(i)) cycle

                    ! Check if this root is already tracked by an active branch
                    do b = 1, WKB_max_tracked_branches
                        if (branch_active(b)) then
                            dist = abs(temp_zeros(i) - branch_center(b))
                            if (dist < WKB_branch_search_halfwidth) then
                                temp_valid(i) = .false.  ! Already tracked
                                exit
                            end if
                        end if
                    end do

                    if (.not. temp_valid(i)) cycle

                    ! Find an inactive branch slot
                    do b = 1, WKB_max_tracked_branches
                        if (.not. branch_active(b)) then
                            ! Activate new branch
                            branch_active(b) = .true.
                            branch_center(b) = temp_zeros(i)
                            branch_fvalue(b) = temp_fzeros(i)
                            branch_mult(b) = temp_mults(i)
                            branch_miss_count(b) = 0
                            n_active_branches = n_active_branches + 1

                            ! Store in results (if not already stored by focused search)
                            n_roots_per_point(j) = n_roots_per_point(j) + 1
                            all_zeros(n_roots_per_point(j), j) = temp_zeros(i)
                            all_fzeros(n_roots_per_point(j), j) = temp_fzeros(i)
                            all_multiplicities(n_roots_per_point(j), j) = temp_mults(i)
                            all_branch_ids(n_roots_per_point(j), j) = b

                            print *, '  New branch ', b, ' at (', real(temp_zeros(i)), ',', &
                                aimag(temp_zeros(i)), ')'
                            temp_valid(i) = .false.
                            exit
                        end if
                    end do
                end do
            end if

            print *, '  Active branches: ', n_active_branches, ', Roots stored: ', n_roots_per_point(j)

        end do  ! End grid loop

        !-----------------------------------------------------------------------
        ! Summary and output
        !-----------------------------------------------------------------------
        print *
        print *, '=== ZEAL Per-Branch Tracking Summary ==='
        print *, 'Total grid points: ', rg_grid%npts_b
        print *, 'Final active branches: ', n_active_branches
        print *, 'Total roots stored: ', sum(n_roots_per_point)

        do b = 1, WKB_max_tracked_branches
            if (branch_active(b)) then
                print *, '  Branch ', b, ': final position (', real(branch_center(b)), ',', &
                    aimag(branch_center(b)), ')'
            end if
        end do

        ! Write tracked branches to file
        call write_tracked_roots(rg_grid%xb, rg_grid%npts_b, &
            all_zeros, all_fzeros, all_multiplicities, &
            n_roots_per_point, all_branch_ids, WKB_max_tracked_branches, &
            'zeal_branches', 'ZEAL per-branch tracked dispersion', dispersion_output_path)

        print *
        print *, 'Tracked branches written to output files.'

        ! Clean up
        deallocate(branch_center, branch_fvalue, branch_mult, branch_active, branch_miss_count)
        deallocate(all_zeros, all_fzeros, all_multiplicities, all_branch_ids, n_roots_per_point)
        deallocate(temp_zeros, temp_fzeros, temp_mults, temp_valid)

    end subroutine

end module 

