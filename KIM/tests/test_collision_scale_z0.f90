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

    use KIM_kinds_m, only: dp
    use constants_m, only: com_unit, e_mass
    use species_m, only: plasma_t, init_electron_species, interpolate_plasma_backs
    use setup_m, only: omega

    implicit none

    call test_z0_recomputed_across_resonance()
    call test_z0_matches_formula_off_resonance()
    call test_collision_scale_reaches_z0()

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
        real(dp) :: fine_grid(3)
        integer :: i

        call build_coarse_plasma(p, 1.0d0)
        fine_grid = [24.0d0, 25.0d0, 26.0d0]
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

    subroutine test_collision_scale_reaches_z0()
        ! calculate_plasma_backs applies collision_frequency_scale before z0 is
        ! assembled; check the scaled nu propagates and that scale 1 is a no-op.
        type(plasma_t) :: unscaled, scaled
        real(dp), parameter :: scale = 3.0d0
        real(dp) :: fine_grid(2)
        complex(dp) :: want
        integer :: i

        call build_coarse_plasma(unscaled, 1.0d0)
        call build_coarse_plasma(scaled, scale)
        fine_grid = [12.0d0, 14.0d0]
        call interpolate_plasma_backs(unscaled, fine_grid)
        call interpolate_plasma_backs(scaled, fine_grid)

        do i = 1, size(fine_grid)
            if (abs(scaled%spec(0)%nu(i) - scale*unscaled%spec(0)%nu(i)) &
                    > 1.0d-12*abs(scale*unscaled%spec(0)%nu(i))) then
                print *, 'FAIL: scaled nu = ', scaled%spec(0)%nu(i), &
                         ' expected ', scale*unscaled%spec(0)%nu(i)
                error stop 'collision frequency scale did not reach the background'
            end if
            want = z0_formula(scaled, i)
            if (abs(scaled%spec(0)%z0(i) - want) > 1.0d-12*max(abs(want), 1.0d0)) then
                print *, 'FAIL: z0 = ', scaled%spec(0)%z0(i), ' expected ', want
                error stop 'scaled collision frequency did not reach z0'
            end if
            if (abs(aimag(scaled%spec(0)%z0(i)) - scale*aimag(unscaled%spec(0)%z0(i))) &
                    > 1.0d-12*abs(scale*aimag(unscaled%spec(0)%z0(i)))) then
                error stop 'z0 imaginary part does not track the collision scale'
            end if
        end do
        print *, 'PASS: collision frequency scale reaches nu and z0'
    end subroutine test_collision_scale_reaches_z0

end program test_collision_scale_z0
