module kernel_test_background_m
    ! Uniform two-species background written directly into species_m state,
    ! so kernel-test expectations derive from the same first-principles
    ! inputs without running the profile pipeline.

    use KIM_kinds_m, only: dp
    use constants_m, only: pi, ev, e_charge, sol, com_unit
    use species_m, only: plasma, init_electron_species, init_hydrogen_species

    implicit none

    integer, parameter :: npts = 100
    real(dp), parameter :: r_inner = 3.0d0
    real(dp), parameter :: r_outer = 70.0d0
    real(dp), parameter :: n0 = 2.0d13
    real(dp), parameter :: T0 = 1.0d3
    real(dp), parameter :: b0 = 2.0d4
    real(dp), parameter :: ks0 = 5.0d-2
    real(dp), parameter :: kp0 = 1.0d-3
    real(dp), parameter :: nu0 = 1.0d4

contains

    subroutine setup_uniform_background()

        implicit none

        integer :: sp, i

        if (.not. allocated(plasma%spec)) then
            plasma%n_species = 2
            allocate (plasma%spec(0:1))
            call init_electron_species(plasma%spec(0))
            call init_hydrogen_species(plasma%spec(1))
        end if

        plasma%grid_size = npts
        if (.not. allocated(plasma%r_grid)) then
            allocate (plasma%r_grid(npts), plasma%ks(npts), plasma%kp(npts), &
                      plasma%om_E(npts))
        end if
        do i = 1, npts
            plasma%r_grid(i) = r_inner + (r_outer - r_inner)*(i - 1)/(npts - 1)
        end do
        plasma%ks = ks0
        plasma%kp = kp0
        plasma%om_E = 0.0d0

        do sp = 0, plasma%n_species - 1
            associate (spec => plasma%spec(sp))
                if (.not. allocated(spec%vT)) then
                    allocate (spec%vT(npts), spec%omega_c(npts), &
                              spec%lambda_D(npts), spec%nu(npts), &
                              spec%z0(npts), spec%A1(npts), spec%A2(npts))
                end if
                spec%vT = sqrt(T0*ev/spec%mass)
                spec%omega_c = spec%Zspec*e_charge*b0/(spec%mass*sol)
                spec%lambda_D = sqrt(T0*ev/(4.0d0*pi*n0 &
                                             *(spec%Zspec*e_charge)**2))
                spec%nu = nu0
                spec%z0 = -(plasma%om_E - com_unit*nu0) &
                          /(abs(kp0)*sqrt(2.0d0)*spec%vT)
                spec%A1 = 0.0d0
                spec%A2 = 0.0d0
            end associate
        end do

    end subroutine setup_uniform_background

    subroutine set_forces(sp, a1_value, a2_value)

        implicit none

        integer, intent(in) :: sp
        real(dp), intent(in) :: a1_value, a2_value

        plasma%spec(sp)%A1 = a1_value
        plasma%spec(sp)%A2 = a2_value

    end subroutine set_forces

    subroutine require_close(name, got, want, rel_tol)

        implicit none

        character(*), intent(in) :: name
        complex(dp), intent(in) :: got, want
        real(dp), intent(in) :: rel_tol

        real(dp) :: scale

        scale = max(abs(want), 1.0d-300)
        if (abs(got - want) <= rel_tol*scale) then
            print *, "PASS ", name
        else
            print *, "FAIL ", name
            print *, "  got  = ", got
            print *, "  want = ", want
            print *, "  rel  = ", abs(got - want)/scale
            error stop
        end if

    end subroutine require_close

end module kernel_test_background_m
