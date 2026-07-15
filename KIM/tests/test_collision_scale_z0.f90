program test_collision_scale_z0
    ! Two background-quantity invariants, extracted from the retired wp5 stack
    ! (PRs #182/#183) onto the forced-periodicity lineage:
    !
    !   (a) interpolate_plasma_backs recomputes the derived Krook argument z0
    !       from interpolated primitives instead of interpolating z0 itself.
    !       z0 carries a 1/|k_parallel| singularity at the resonance, so a
    !       four-point Lagrange stencil spanning k_parallel=0 can return values
    !       many orders of magnitude away from the defining formula
    !         z0 = -(omega_E - omega - i nu)/(|k_parallel| sqrt(2) v_T).
    !       Asserted both on a stencil that straddles the resonance (the
    !       regression proper) and on a smooth off-resonance point.
    !
    !   (b) collision_frequency_scale multiplies the calculated collision
    !       frequencies, and the default 1.0 leaves them untouched. The scaled
    !       nu must reach z0, since z0 is assembled after the scaling.

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use constants_m, only: com_unit
    use species_m, only: plasma, plasma_t, init_electron_species, &
                         init_deuterium_species, calculate_plasma_backs, &
                         interpolate_plasma_backs, compute_rg_cell_centers, &
                         calculate_thermodynamic_forces_and_susc
    use config_m, only: number_of_ion_species, collision_frequency_scale, &
                        rescale_density, ion_flr_scale_factor
    use grid_m, only: rg_grid
    use setup_m, only: omega, collisions_off, mphi_max

    implicit none

    call test_z0_recomputed_across_resonance()
    call test_z0_matches_formula_off_resonance()
    call test_collision_scale_reaches_production_paths()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine build_coarse_plasma(coarse, nu_scale)
        ! Minimal single-species background whose k_parallel changes sign
        ! between grid points 2 and 3, i.e. a resonance inside the stencil.
        type(plasma_t), intent(out) :: coarse
        real(dp), intent(in) :: nu_scale
        integer :: n

        n = 4
        coarse%n_species = 1
        coarse%grid_size = n
        allocate(coarse%spec(0:0))
        call init_electron_species(coarse%spec(0))
        allocate(coarse%r_grid(n), coarse%ks(n), coarse%kp(n), coarse%om_E(n), &
                 coarse%B0(n), coarse%q(n), coarse%dqdr(n), coarse%Er(n))
        allocate(coarse%spec(0)%n(n), coarse%spec(0)%dndr(n), coarse%spec(0)%T(n), &
                 coarse%spec(0)%dTdr(n), coarse%spec(0)%nu(n), coarse%spec(0)%vT(n), &
                 coarse%spec(0)%omega_c(n), coarse%spec(0)%lambda_D(n), &
                 coarse%spec(0)%rho_L(n), coarse%spec(0)%z0(n))

        coarse%r_grid = [10.0d0, 20.0d0, 30.0d0, 40.0d0]
        coarse%kp     = [3.0d-3, 1.0d-3, -1.0d-3, -3.0d-3]
        coarse%om_E   = [9.3d4, 8.0d4, 7.0d4, 6.0d4]
        coarse%ks     = 1.0d-2
        coarse%B0     = 2.0d4
        coarse%q      = 1.5d0
        coarse%dqdr   = 1.0d-2
        coarse%Er     = 1.0d0

        coarse%spec(0)%n       = 1.0d13
        coarse%spec(0)%dndr    = -1.0d11
        coarse%spec(0)%T       = 1.0d3
        coarse%spec(0)%dTdr    = -1.0d1
        coarse%spec(0)%vT      = 5.9d8
        coarse%spec(0)%omega_c = -3.5d11
        coarse%spec(0)%lambda_D = 7.4d-3
        coarse%spec(0)%rho_L   = 1.7d-3
        coarse%spec(0)%nu      = 1.0d4 * nu_scale

        ! z0 consistent with the coarse primitives; the point of the fix is
        ! that this array is NOT carried through the interpolation.
        coarse%spec(0)%z0 = -(coarse%om_E - omega - com_unit*coarse%spec(0)%nu) &
                            /(abs(coarse%kp)*sqrt(2d0)*coarse%spec(0)%vT)
    end subroutine build_coarse_plasma

    complex(dp) function z0_formula(p, i) result(z0)
        type(plasma_t), intent(in) :: p
        integer, intent(in) :: i

        z0 = -(p%om_E(i) - omega - com_unit*p%spec(0)%nu(i)) &
             /(abs(p%kp(i))*sqrt(2d0)*p%spec(0)%vT(i))
    end function z0_formula

    subroutine require_z0_matches_formula(p, i, label)
        type(plasma_t), intent(in) :: p
        integer, intent(in) :: i
        character(*), intent(in) :: label
        complex(dp) :: want
        real(dp) :: err

        want = z0_formula(p, i)
        err = abs(p%spec(0)%z0(i) - want)/max(abs(want), 1.0d0)
        if (.not. ieee_is_finite(real(p%spec(0)%z0(i), dp)) .or. &
                .not. ieee_is_finite(aimag(p%spec(0)%z0(i))) .or. &
                .not. ieee_is_finite(real(want, dp)) .or. &
                .not. ieee_is_finite(aimag(want)) .or. &
                .not. ieee_is_finite(err)) then
            print *, 'FAIL: ', label
            print *, '  k_par   = ', p%kp(i)
            print *, '  z0 got  = ', p%spec(0)%z0(i)
            print *, '  z0 want = ', want
            error stop 'z0 comparison contains a non-finite value'
        end if
        if (err > 1.0d-12) then
            print *, 'FAIL: ', label
            print *, '  r        = ', p%r_grid(i)
            print *, '  k_par    = ', p%kp(i)
            print *, '  z0 got   = ', p%spec(0)%z0(i)
            print *, '  z0 want  = ', want
            print *, '  rel err  = ', err
            error stop 'z0 violates its defining formula after interpolation'
        end if
    end subroutine require_z0_matches_formula

    subroutine test_z0_recomputed_across_resonance()
        ! Target grid point sits where k_parallel is near zero, so the old
        ! interpolated-z0 path blew up while the recomputed one stays exact.
        type(plasma_t) :: p
        real(dp) :: fine_grid(4)
        integer :: i

        call build_coarse_plasma(p, 1.0d0)
        ! Sample finite points on both sides of the singular surface. The exact
        ! k_parallel=0 point is mathematically singular and is not a valid
        ! finite-value invariant check.
        fine_grid = [24.0d0, 24.9d0, 25.1d0, 26.0d0]
        call interpolate_plasma_backs(p, fine_grid)

        do i = 1, size(fine_grid)
            call require_z0_matches_formula(p, i, 'z0 across the resonance')
        end do
        print *, 'PASS: z0 recomputed from primitives across k_parallel = 0'
    end subroutine test_z0_recomputed_across_resonance

    subroutine test_z0_matches_formula_off_resonance()
        type(plasma_t) :: p
        real(dp) :: fine_grid(2)
        integer :: i

        call build_coarse_plasma(p, 1.0d0)
        fine_grid = [12.0d0, 14.0d0]
        call interpolate_plasma_backs(p, fine_grid)

        do i = 1, size(fine_grid)
            call require_z0_matches_formula(p, i, 'z0 away from the resonance')
        end do
        print *, 'PASS: z0 obeys its formula away from the resonance'
    end subroutine test_z0_matches_formula_off_resonance

    subroutine require_close(got, want, label)
        real(dp), intent(in) :: got, want
        character(*), intent(in) :: label

        if (abs(got - want) > 1.0d-12*max(abs(want), 1.0d0)) then
            print *, 'FAIL: ', label
            print *, '  got  = ', got
            print *, '  want = ', want
            error stop 'scaled production quantity has the wrong value'
        end if
    end subroutine require_close

    subroutine build_collision_plasma(template)
        type(plasma_t), intent(out) :: template
        integer, parameter :: n = 4
        integer :: sp

        template%n_species = 2
        template%grid_size = n
        allocate(template%spec(0:1))
        call init_electron_species(template%spec(0))
        call init_deuterium_species(template%spec(1))
        allocate(template%r_grid(n), template%ks(n), template%kp(n), &
                 template%om_E(n), template%B0(n), template%q(n), &
                 template%dqdr(n), template%Er(n))
        do sp = 0, 1
            allocate(template%spec(sp)%n(n), template%spec(sp)%dndr(n), &
                     template%spec(sp)%T(n), template%spec(sp)%dTdr(n))
        end do

        template%r_grid = [10.0d0, 20.0d0, 30.0d0, 40.0d0]
        template%kp = [3.0d-3, 2.0d-3, 1.0d-3, 5.0d-4]
        template%om_E = [9.3d4, 8.0d4, 7.0d4, 6.0d4]
        template%ks = 1.0d-2
        template%B0 = 2.0d4
        template%q = 1.5d0
        template%dqdr = 1.0d-2
        template%Er = 0.0d0
        template%spec(0)%n = 1.0d13
        template%spec(1)%n = 1.0d13
        template%spec(0)%T = 1.0d3
        template%spec(1)%T = 8.0d2
        do sp = 0, 1
            template%spec(sp)%dndr = -1.0d11
            template%spec(sp)%dTdr = -1.0d1
        end do

        rg_grid%npts_b = n
        rg_grid%npts_c = n - 1
        if (allocated(rg_grid%xb)) deallocate(rg_grid%xb)
        if (allocated(rg_grid%xc)) deallocate(rg_grid%xc)
        allocate(rg_grid%xb(n), rg_grid%xc(n - 1))
        rg_grid%xb = template%r_grid
        rg_grid%xc = 0.5d0*(rg_grid%xb(:n - 1) + rg_grid%xb(2:))
    end subroutine build_collision_plasma

    subroutine run_collision_scale_case(template, scale, nu, z0, x1, x2)
        type(plasma_t), intent(in) :: template
        real(dp), intent(in) :: scale
        real(dp), intent(out) :: nu(0:1), x1(0:1), x2(0:1)
        complex(dp), intent(out) :: z0(0:1)
        integer :: sp

        plasma = template
        collision_frequency_scale = scale
        call calculate_plasma_backs(plasma)
        call compute_rg_cell_centers(plasma)
        call calculate_thermodynamic_forces_and_susc(plasma)

        do sp = 0, 1
            nu(sp) = plasma%spec(sp)%nu(2)
            z0(sp) = plasma%spec(sp)%z0(2)
            x1(sp) = plasma%spec(sp)%x1(2)
            x2(sp) = plasma%spec(sp)%x2(2, 0)
        end do
    end subroutine run_collision_scale_case

    subroutine test_collision_scale_reaches_production_paths()
        type(plasma_t) :: template
        real(dp), parameter :: scale = 3.0d0
        real(dp) :: nu_base(0:1), nu_scaled(0:1)
        real(dp) :: x1_base(0:1), x1_scaled(0:1)
        real(dp) :: x2_base(0:1), x2_scaled(0:1)
        real(dp) :: lee, lei, expected_nue
        complex(dp) :: z0_base(0:1), z0_scaled(0:1)
        integer :: sp

        number_of_ion_species = 1
        rescale_density = .false.
        ion_flr_scale_factor = 1.0d0
        collisions_off = .false.
        mphi_max = 0
        omega = 1.0d5
        call build_collision_plasma(template)
        call run_collision_scale_case(template, 1.0d0, nu_base, z0_base, &
                                      x1_base, x2_base)
        call run_collision_scale_case(template, scale, nu_scaled, z0_scaled, &
                                      x1_scaled, x2_scaled)

        lee = 23.5d0 - log(sqrt(1.0d13)/1.0d3**1.25d0) &
            - sqrt(1.0d-5 + (log(1.0d3) - 2.0d0)**2/16.0d0)
        lei = 24.0d0 - log(sqrt(1.0d13)/1.0d3)
        expected_nue = real(5.8e-6, dp)*1.0d13*lee/1.0d3**1.5d0 &
            + 7.7d-6*1.0d13*lei/1.0d3**1.5d0
        call require_close(nu_base(0), expected_nue, &
                           'default scale preserves the electron collision formula')

        do sp = 0, 1
            call require_close(nu_scaled(sp), scale*nu_base(sp), &
                               'collision scale reaches electron/ion nu')
            call require_close(real(z0_scaled(sp), dp), real(z0_base(sp), dp), &
                               'collision scale leaves real(z0) unchanged')
            call require_close(aimag(z0_scaled(sp)), scale*aimag(z0_base(sp)), &
                               'collision scale reaches imag(z0)')
            call require_close(x1_scaled(sp), x1_base(sp)/scale, &
                               'scaled nu reaches FP x1')
            call require_close(x2_scaled(sp), x2_base(sp)/scale, &
                               'scaled nu reaches FP x2')
        end do
        print *, 'PASS: collision scale reaches electron/ion nu, z0, and FP inputs'
    end subroutine test_collision_scale_reaches_production_paths

end program test_collision_scale_z0
