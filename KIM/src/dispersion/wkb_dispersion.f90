module rt_WKB_dispersion_m

    !==========================================================================================
    ! module for running the WKB dispersion relation solver
    ! Algorithms for root finding are Muller and ZEAL, both with per-branch tracking
    ! Dispersion functions implemented: KIM and FLRE to second order
    !
    ! Configuration namelist is in KIM_config.nml named &wkb_dispersion
    !==========================================================================================

    use kim_base_m, only: kim_t
    use KIM_kinds_m, only: dp

    implicit none

    integer, parameter :: ZEAL_M_MAX = 5
    integer, parameter :: ZEAL_MAX_MISS_COUNT = 3

    type zeal_branch_state_t
        complex(dp), allocatable :: center(:)
        complex(dp), allocatable :: fvalue(:)
        integer, allocatable :: mult(:)
        logical, allocatable :: active(:)
        integer, allocatable :: miss_count(:)
        integer :: n_active = 0
    end type zeal_branch_state_t

    type zeal_root_store_t
        complex(dp), allocatable :: zeros(:,:)
        complex(dp), allocatable :: fzeros(:,:)
        integer, allocatable :: multiplicities(:,:)
        integer, allocatable :: branch_id(:,:)
        integer, allocatable :: n_per_point(:)
    end type zeal_root_store_t

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
        use config_m, only: WKB_dispersion_solver, dispersion_output_path, WKB_solve_for_kr_squared

        implicit none

        class(WKB_dispersion_t), intent(inout) :: this
        character(len=20) :: out_str
        complex(dp), allocatable :: found_roots(:,:)
        integer :: km, i

        print *, "Running "//trim(this%run_type)//" model ..."
        if (WKB_solve_for_kr_squared) then
            print *, "  Solver mode: solving for kr^2 (D(kr^2)=0)"
        else
            print *, "  Solver mode: solving for kr (D(kr)=0)"
        end if

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
        use config_m, only: WKB_dispersion_mode, dispersion_output_path, WKB_solve_for_kr_squared
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
        print *, 'Dispersion mode: ', trim(WKB_dispersion_mode)
        if (WKB_solve_for_kr_squared) then
            print *, 'Solver variable: kr^2 (solving D(kr^2)=0)'
        else
            print *, 'Solver variable: kr (solving D(kr)=0)'
        end if
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
            ! When WKB_solve_for_kr_squared=.false., kr is the unknown (solve D(kr)=0)
            ! When WKB_solve_for_kr_squared=.true., kr represents kr^2 (solve D(kr^2)=0)

            use species_m, only: plasma
            use constants_m, only: com_unit
            use KIM_kinds_m, only: dp
            use grid_m, only: rg_grid
            use config_m, only: WKB_solve_for_kr_squared

            implicit none

            complex(dp), intent(in)  :: kr
            complex(dp), intent(out) :: f
            integer :: sp
            integer :: ifail_bess, nz
            complex(dp) :: bess0, bessm1
            real(dp) :: bess0_re(2), bess0_im(2)
            complex(dp) :: kr_rho_squared
            complex(dp) :: kr_squared  ! Either kr^2 or kr depending on mode

            bess0_re = 0.0d0
            bess0_im = 0.0d0
            f = (0.0d0, 0.0d0)

            ! Determine kr^2 based on solver mode
            if (WKB_solve_for_kr_squared) then
                kr_squared = kr  ! Input is already kr^2
            else
                kr_squared = kr**2.0d0  ! Input is kr, compute kr^2
            end if

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr_squared * plasma%spec(sp)%rho_L(j)**2.0d0

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

            ! neglect linear term, accurate enough
            f = kr_squared + plasma%kp(j)**2.0d0 - f

        end subroutine dispersion_KIM

        subroutine dispersion_FLRE(kr, f)
            ! dispersion equation for second order finite Larmor radius expansion
            ! When WKB_solve_for_kr_squared=.false., kr is the unknown (solve D(kr)=0)
            ! When WKB_solve_for_kr_squared=.true., kr represents kr^2 (solve D(kr^2)=0)

            use species_m, only: plasma
            use constants_m, only: com_unit
            use KIM_kinds_m, only: dp
            use grid_m, only: rg_grid
            use config_m, only: WKB_solve_for_kr_squared

            implicit none

            complex(dp), intent(in)  :: kr
            complex(dp), intent(out) :: f
            integer :: sp
            complex(dp) :: kr_rho_squared
            complex(dp) :: kr_squared  ! Either kr^2 or kr depending on mode

            f = (0.0d0, 0.0d0)

            ! Determine kr^2 based on solver mode
            if (WKB_solve_for_kr_squared) then
                kr_squared = kr  ! Input is already kr^2
            else
                kr_squared = kr**2.0d0  ! Input is kr, compute kr^2
            end if

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr_squared * plasma%spec(sp)%rho_L(j)**2.0d0

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

            ! use the version that neglects the first order derivative term (is accurate enough)
            f = kr_squared + plasma%kp(j)**2.0d0 - f

        end subroutine dispersion_FLRE


    end subroutine


    subroutine test_dispersion_KIM(kr, f)
        ! Standalone version of dispersion_KIM for testing/comparison
        ! When WKB_solve_for_kr_squared=.false., kr is the unknown (solve D(kr)=0)
        ! When WKB_solve_for_kr_squared=.true., kr represents kr^2 (solve D(kr^2)=0)
        use species_m, only: plasma
        use constants_m, only: com_unit
        use KIM_kinds_m, only: dp
        use Function_Input_Module, only: rg_index
        use config_m, only: WKB_solve_for_kr_squared

        implicit none

        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: f
        integer :: sp, j
        integer :: ifail_bess, nz
        complex(dp) :: bess0, bessm1
        real(dp) :: bess0_re(2), bess0_im(2)
        complex(dp) :: kr_rho_squared
        complex(dp) :: kr_squared  ! Either kr^2 or kr depending on mode

        j = rg_index

        bess0_re = 0.0d0
        bess0_im = 0.0d0
        f = (0.0d0, 0.0d0)

        ! Determine kr^2 based on solver mode
        if (WKB_solve_for_kr_squared) then
            kr_squared = kr  ! Input is already kr^2
        else
            kr_squared = kr**2.0d0  ! Input is kr, compute kr^2
        end if

        do sp = 0, plasma%n_species-1

            kr_rho_squared = kr_squared * plasma%spec(sp)%rho_L(j)**2.0d0

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

        f = kr_squared + plasma%kp(j)**2.0d0 - f

    end subroutine test_dispersion_KIM


    subroutine run_ZEAL_dispersion()
        use Function_Input_Module, only: rg_index, init_dispersion_mode
        use grid_m, only: rg_grid
        use IO_collection_m, only: write_tracked_roots
        use config_m, only: dispersion_output_path, WKB_verbose

        implicit none

        type(zeal_branch_state_t) :: branches
        type(zeal_root_store_t) :: roots
        integer :: j, n_output_branches

        call allocate_zeal_tracking(branches, roots, rg_grid%npts_b)
        call init_dispersion_mode()
        call print_zeal_start()

        do j = 1, rg_grid%npts_b
            rg_index = j
            if (WKB_verbose) print *, '--- Grid point ', j, ' at x = ', rg_grid%xb(j), ' ---'

            call track_zeal_branches(j, branches, roots)
            if (zeal_needs_broad_search(j, branches)) call run_zeal_broad_search(j, branches, roots)

            if (WKB_verbose) then
                print *, '  Active branches: ', branches%n_active, &
                    ', Roots stored: ', roots%n_per_point(j)
                print *
            end if
        end do

        call compact_zeal_branch_ids(roots, n_output_branches)
        call print_zeal_summary(branches, roots, n_output_branches)
        call write_tracked_roots(rg_grid%xb, rg_grid%npts_b, &
            roots%zeros, roots%fzeros, roots%multiplicities, &
            roots%n_per_point, roots%branch_id, n_output_branches, &
            'zeal_branches', 'ZEAL per-branch tracked dispersion', dispersion_output_path)
        call deallocate_zeal_tracking(branches, roots)

        print *, 'Tracked branches written to output files.'

    end subroutine run_ZEAL_dispersion

    subroutine allocate_zeal_tracking(branches, roots, n_grid)
        use config_m, only: WKB_max_tracked_branches

        implicit none

        type(zeal_branch_state_t), intent(out) :: branches
        type(zeal_root_store_t), intent(out) :: roots
        integer, intent(in) :: n_grid

        allocate(branches%center(WKB_max_tracked_branches))
        allocate(branches%fvalue(WKB_max_tracked_branches))
        allocate(branches%mult(WKB_max_tracked_branches))
        allocate(branches%active(WKB_max_tracked_branches))
        allocate(branches%miss_count(WKB_max_tracked_branches))
        allocate(roots%zeros(WKB_max_tracked_branches, n_grid))
        allocate(roots%fzeros(WKB_max_tracked_branches, n_grid))
        allocate(roots%multiplicities(WKB_max_tracked_branches, n_grid))
        allocate(roots%branch_id(WKB_max_tracked_branches, n_grid))
        allocate(roots%n_per_point(n_grid))

        branches%center = (0.0_dp, 0.0_dp)
        branches%fvalue = (0.0_dp, 0.0_dp)
        branches%mult = 1
        branches%active = .false.
        branches%miss_count = 0
        branches%n_active = 0
        roots%zeros = (0.0_dp, 0.0_dp)
        roots%fzeros = (0.0_dp, 0.0_dp)
        roots%multiplicities = 0
        roots%branch_id = 0
        roots%n_per_point = 0

    end subroutine allocate_zeal_tracking

    subroutine deallocate_zeal_tracking(branches, roots)
        implicit none

        type(zeal_branch_state_t), intent(inout) :: branches
        type(zeal_root_store_t), intent(inout) :: roots

        deallocate(branches%center, branches%fvalue, branches%mult)
        deallocate(branches%active, branches%miss_count)
        deallocate(roots%zeros, roots%fzeros, roots%multiplicities)
        deallocate(roots%branch_id, roots%n_per_point)

    end subroutine deallocate_zeal_tracking

    subroutine print_zeal_start()
        use config_m, only: WKB_max_tracked_branches, WKB_branch_search_halfwidth, &
            WKB_broad_search_halfwidth, WKB_broad_search_interval, WKB_dispersion_mode, &
            WKB_solve_for_kr_squared, WKB_verbose

        implicit none

        print *
        print *, '=== Per-Branch ZEAL Tracking ==='
        print *, 'Dispersion mode: ', trim(WKB_dispersion_mode)
        if (WKB_solve_for_kr_squared) then
            print *, 'Solver variable: kr^2 (solving D(kr^2)=0)'
        else
            print *, 'Solver variable: kr (solving D(kr)=0)'
        end if
        print *, 'Max tracked branches: ', WKB_max_tracked_branches
        print *, 'Branch search half-width: ', WKB_branch_search_halfwidth
        print *, 'Broad search half-width: ', WKB_broad_search_halfwidth
        print *, 'Broad search interval: ', WKB_broad_search_interval
        print *, 'Verbose branch tracing: ', WKB_verbose
        print *

    end subroutine print_zeal_start

    subroutine track_zeal_branches(j, branches, roots)
        use config_m, only: WKB_max_tracked_branches, WKB_branch_search_halfwidth, WKB_verbose

        implicit none

        integer, intent(in) :: j
        type(zeal_branch_state_t), intent(inout) :: branches
        type(zeal_root_store_t), intent(inout) :: roots
        integer :: b, best_idx, distinctnumber
        complex(dp), allocatable :: zeros(:), fzeros(:)
        integer, allocatable :: multiplicities(:)
        complex(dp) :: ll, ur

        do b = 1, WKB_max_tracked_branches
            if (.not. branches%active(b)) cycle

            ll = branches%center(b) - cmplx(WKB_branch_search_halfwidth, &
                WKB_branch_search_halfwidth, dp)
            ur = branches%center(b) + cmplx(WKB_branch_search_halfwidth, &
                WKB_branch_search_halfwidth, dp)
            if (WKB_verbose) print *, '  Branch ', b, ': searching near ', branches%center(b)

            call zeal_region_search('branch', j, ll, ur, zeros, fzeros, &
                multiplicities, distinctnumber)
            call best_zeal_root(fzeros, distinctnumber, best_idx)
            if (best_idx > 0) then
                call update_zeal_branch(branches, roots, j, b, zeros(best_idx), &
                    fzeros(best_idx), multiplicities(best_idx))
            else
                call record_zeal_miss(branches, b)
            end if
        end do

    end subroutine track_zeal_branches

    subroutine run_zeal_broad_search(j, branches, roots)
        use config_m, only: WKB_broad_search_halfwidth, WKB_root_tolerance, WKB_verbose

        implicit none

        integer, intent(in) :: j
        type(zeal_branch_state_t), intent(inout) :: branches
        type(zeal_root_store_t), intent(inout) :: roots
        integer :: i, distinctnumber
        complex(dp), allocatable :: zeros(:), fzeros(:)
        integer, allocatable :: multiplicities(:)
        complex(dp) :: center, ll, ur

        center = zeal_broad_search_center(branches)
        ll = center - cmplx(WKB_broad_search_halfwidth, WKB_broad_search_halfwidth, dp)
        ur = center + cmplx(WKB_broad_search_halfwidth, WKB_broad_search_halfwidth, dp)

        if (WKB_verbose) then
            print *, '  Running broad search centered at ', center, &
                ' with half-width ', WKB_broad_search_halfwidth
        end if

        call zeal_region_search('broad', j, ll, ur, zeros, fzeros, multiplicities, distinctnumber)
        if (WKB_verbose) print *, '  Broad search found ', distinctnumber, ' zeros'
        if (distinctnumber <= 0) return
        if (.not. allocated(zeros)) return

        do i = 1, min(distinctnumber, size(zeros))
            if (abs(fzeros(i)) >= WKB_root_tolerance) cycle
            if (zeal_root_is_tracked(zeros(i), branches)) cycle
            if (.not. activate_zeal_branch(branches, roots, j, zeros(i), &
                    fzeros(i), multiplicities(i))) then
                print *, 'Warning: ZEAL broad search found more valid roots than tracked ', &
                    'branch slots.'
                exit
            end if
        end do

    end subroutine run_zeal_broad_search

    subroutine zeal_region_search(label, j, ll, ur, zeros, fzeros, multiplicities, distinctnumber)
        use Function_Input_Module, only: dispersion_region_fn
        use fortnum_roots_complex, only: complex_region_roots
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK

        implicit none

        character(len=*), intent(in) :: label
        integer, intent(in) :: j
        complex(dp), intent(in) :: ll, ur
        complex(dp), allocatable, intent(out) :: zeros(:), fzeros(:)
        integer, allocatable, intent(out) :: multiplicities(:)
        integer, intent(out) :: distinctnumber
        type(fortnum_status_t) :: rstatus

        call complex_region_roots(dispersion_region_fn, ll, ur, &
            zeros, fzeros, multiplicities, distinctnumber, rstatus, m_max=ZEAL_M_MAX)
        if (rstatus%code == FORTNUM_OK) return

        print *, 'Warning: ZEAL ', trim(label), ' search failed at grid point ', j, ': ', &
            trim(rstatus%msg)
        distinctnumber = 0

    end subroutine zeal_region_search

    subroutine best_zeal_root(fzeros, distinctnumber, best_idx)
        use config_m, only: WKB_root_tolerance

        implicit none

        complex(dp), allocatable, intent(in) :: fzeros(:)
        integer, intent(in) :: distinctnumber
        integer, intent(out) :: best_idx
        integer :: i, n_search

        best_idx = 0
        if (distinctnumber <= 0) return
        if (.not. allocated(fzeros)) return

        n_search = min(distinctnumber, size(fzeros))
        do i = 1, n_search
            if (abs(fzeros(i)) < WKB_root_tolerance) then
                if (best_idx == 0) then
                    best_idx = i
                else if (abs(fzeros(i)) < abs(fzeros(best_idx))) then
                    best_idx = i
                end if
            end if
        end do

    end subroutine best_zeal_root

    subroutine update_zeal_branch(branches, roots, j, branch_id, zero, fzero, multiplicity)
        use config_m, only: WKB_verbose

        implicit none

        type(zeal_branch_state_t), intent(inout) :: branches
        type(zeal_root_store_t), intent(inout) :: roots
        integer, intent(in) :: j, branch_id, multiplicity
        complex(dp), intent(in) :: zero, fzero

        branches%center(branch_id) = zero
        branches%fvalue(branch_id) = fzero
        branches%mult(branch_id) = multiplicity
        branches%miss_count(branch_id) = 0
        call store_zeal_root(roots, j, branch_id, zero, fzero, multiplicity)

        if (WKB_verbose) then
            print *, '    Found branch ', branch_id, ' root ', zero, ' |f| = ', abs(fzero)
        end if

    end subroutine update_zeal_branch

    subroutine record_zeal_miss(branches, branch_id)
        use config_m, only: WKB_verbose

        implicit none

        type(zeal_branch_state_t), intent(inout) :: branches
        integer, intent(in) :: branch_id

        branches%miss_count(branch_id) = branches%miss_count(branch_id) + 1
        if (WKB_verbose) then
            print *, '    No valid root found for branch ', branch_id, &
                ' (miss count: ', branches%miss_count(branch_id), ')'
        end if
        if (branches%miss_count(branch_id) < ZEAL_MAX_MISS_COUNT) return

        branches%active(branch_id) = .false.
        branches%n_active = branches%n_active - 1
        if (WKB_verbose) then
            print *, '    Branch ', branch_id, ' deactivated after ', ZEAL_MAX_MISS_COUNT, &
                ' consecutive misses'
        end if

    end subroutine record_zeal_miss

    logical function zeal_needs_broad_search(j, branches)
        use config_m, only: WKB_broad_search_interval

        implicit none

        integer, intent(in) :: j
        type(zeal_branch_state_t), intent(in) :: branches

        zeal_needs_broad_search = (j == 1)
        if (.not. zeal_needs_broad_search) then
            if (WKB_broad_search_interval > 0) then
                if (mod(j - 1, WKB_broad_search_interval) == 0) zeal_needs_broad_search = .true.
            end if
        end if
        if (branches%n_active == 0) zeal_needs_broad_search = .true.

    end function zeal_needs_broad_search

    complex(dp) function zeal_broad_search_center(branches)
        use config_m, only: WKB_max_tracked_branches

        implicit none

        type(zeal_branch_state_t), intent(in) :: branches
        integer :: b, n_centers

        zeal_broad_search_center = (0.0_dp, -1.0_dp)
        if (branches%n_active <= 0) return

        zeal_broad_search_center = (0.0_dp, 0.0_dp)
        n_centers = 0
        do b = 1, WKB_max_tracked_branches
            if (branches%active(b)) then
                zeal_broad_search_center = zeal_broad_search_center + branches%center(b)
                n_centers = n_centers + 1
            end if
        end do
        if (n_centers > 0) zeal_broad_search_center = zeal_broad_search_center / real(n_centers, dp)

    end function zeal_broad_search_center

    logical function zeal_root_is_tracked(zero, branches)
        use config_m, only: WKB_max_tracked_branches, WKB_branch_search_halfwidth

        implicit none

        complex(dp), intent(in) :: zero
        type(zeal_branch_state_t), intent(in) :: branches
        integer :: b

        zeal_root_is_tracked = .false.
        do b = 1, WKB_max_tracked_branches
            if (branches%active(b)) then
                if (abs(zero - branches%center(b)) < WKB_branch_search_halfwidth) then
                    zeal_root_is_tracked = .true.
                    return
                end if
            end if
        end do

    end function zeal_root_is_tracked

    logical function activate_zeal_branch(branches, roots, j, zero, fzero, multiplicity)
        use config_m, only: WKB_max_tracked_branches, WKB_verbose

        implicit none

        type(zeal_branch_state_t), intent(inout) :: branches
        type(zeal_root_store_t), intent(inout) :: roots
        integer, intent(in) :: j, multiplicity
        complex(dp), intent(in) :: zero, fzero
        integer :: b

        activate_zeal_branch = .false.
        do b = 1, WKB_max_tracked_branches
            if (branches%active(b)) cycle

            branches%active(b) = .true.
            branches%center(b) = zero
            branches%fvalue(b) = fzero
            branches%mult(b) = multiplicity
            branches%miss_count(b) = 0
            branches%n_active = branches%n_active + 1
            call store_zeal_root(roots, j, b, zero, fzero, multiplicity)
            activate_zeal_branch = .true.
            if (WKB_verbose) print *, '  New branch ', b, ' at ', zero
            return
        end do

    end function activate_zeal_branch

    subroutine store_zeal_root(roots, j, branch_id, zero, fzero, multiplicity)
        implicit none

        type(zeal_root_store_t), intent(inout) :: roots
        integer, intent(in) :: j, branch_id, multiplicity
        complex(dp), intent(in) :: zero, fzero
        integer :: slot

        if (roots%n_per_point(j) >= size(roots%zeros, 1)) then
            print *, 'Warning: ZEAL root storage full at grid point ', j, &
                '; dropping branch ', branch_id
            return
        end if

        slot = roots%n_per_point(j) + 1
        roots%n_per_point(j) = slot
        roots%zeros(slot, j) = zero
        roots%fzeros(slot, j) = fzero
        roots%multiplicities(slot, j) = multiplicity
        roots%branch_id(slot, j) = branch_id

    end subroutine store_zeal_root

    subroutine compact_zeal_branch_ids(roots, n_branches)
        implicit none

        type(zeal_root_store_t), intent(inout) :: roots
        integer, intent(out) :: n_branches
        integer, allocatable :: branch_map(:)
        integer :: i, j, branch_id, max_branches

        max_branches = size(roots%zeros, 1)
        allocate(branch_map(max_branches))
        branch_map = 0
        n_branches = 0

        do j = 1, size(roots%zeros, 2)
            do i = 1, roots%n_per_point(j)
                branch_id = roots%branch_id(i, j)
                if (branch_id < 1) then
                    roots%branch_id(i, j) = 0
                else if (branch_id > max_branches) then
                    roots%branch_id(i, j) = 0
                else
                    if (branch_map(branch_id) == 0) then
                        n_branches = n_branches + 1
                        branch_map(branch_id) = n_branches
                    end if
                end if
            end do
        end do

        do j = 1, size(roots%zeros, 2)
            do i = 1, roots%n_per_point(j)
                branch_id = roots%branch_id(i, j)
                if (branch_id > 0) roots%branch_id(i, j) = branch_map(branch_id)
            end do
        end do

        deallocate(branch_map)

    end subroutine compact_zeal_branch_ids

    subroutine print_zeal_summary(branches, roots, n_output_branches)
        use config_m, only: WKB_max_tracked_branches, WKB_verbose

        implicit none

        type(zeal_branch_state_t), intent(in) :: branches
        type(zeal_root_store_t), intent(in) :: roots
        integer, intent(in) :: n_output_branches
        integer :: b

        print *
        print *, '=== ZEAL Per-Branch Tracking Summary ==='
        print *, 'Total grid points: ', size(roots%n_per_point)
        print *, 'Final active branches: ', branches%n_active
        print *, 'Output branches with stored roots: ', n_output_branches
        print *, 'Total roots stored: ', sum(roots%n_per_point)

        if (WKB_verbose) then
            do b = 1, WKB_max_tracked_branches
                if (branches%active(b)) then
                    print *, '  Branch ', b, ': final position ', branches%center(b)
                end if
            end do
        end if
        print *

    end subroutine print_zeal_summary

end module
