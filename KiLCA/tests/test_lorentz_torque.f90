! Preserve the Lorentz-force identities introduced on main before the migration.
program test_lorentz_torque
    use, intrinsic :: iso_c_binding, only: c_double
    use kilca_flre_quants_m, only: calc_time_averaged_lorentz_force, &
                                   calc_cylindrical_torque_density, integrate_over_cylinder
    implicit none

    real(c_double), parameter :: cspeed = 29979245800.0d0
    complex(c_double) :: density, electric(3), current(3), magnetic(3)
    real(c_double) :: force(3), torque(3), integrated(3)
    integer :: failures

    failures = 0
    density = (1.0d0, 2.0d0)
    electric = [(3.0d0, 4.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0)]
    current = [(1.0d0, 1.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0)]
    magnetic = [(0.0d0, 0.0d0), (0.0d0, 0.0d0), (2.0d0, 3.0d0)]
    call calc_time_averaged_lorentz_force(2.0d0, density, electric, current, magnetic, force)
    call check('electrostatic force', force(1), 11.0d0)
    call check('magnetic cross-product sign', force(2) * cspeed, -2.5d0)
    call check('zero force component', force(3), 0.0d0)
    call calc_cylindrical_torque_density(5.0d0, 7.0d0, force, torque)
    call check('radial component', torque(1), 11.0d0)
    call check('poloidal torque moment arm', torque(2) * cspeed, -12.5d0)
    call check('zero toroidal torque', torque(3), 0.0d0)
    call calc_cylindrical_torque_density(5.0d0, 7.0d0, [1.0d0, 2.0d0, -3.0d0], torque)
    call check('toroidal torque moment arm', torque(3), -21.0d0)
    ! A common phase rotation cannot change a time-averaged force.
    call calc_time_averaged_lorentz_force(2.0d0, density * (0.0d0, 1.0d0), &
             electric * (0.0d0, 1.0d0), current * (0.0d0, 1.0d0), magnetic * (0.0d0, 1.0d0), torque)
    call check('phase-invariant electric force', torque(1), force(1))
    call check('phase-invariant magnetic force', torque(2) * cspeed, force(2) * cspeed)
    call integrate_over_cylinder(3, [1.0d0, 2.0d0, 4.0d0], [1.0d0, -2.0d0, 3.0d0], &
                                 10.0d0, integrated)
    call check('integral origin', integrated(1), 0.0d0)
    call check('integral first interval', integrated(2), -15.0d0)
    call check('integral signed total', integrated(3), 65.0d0)
    if (failures /= 0) stop 1
    print *, 'PASS KiLCA Lorentz torque identities'

contains

    subroutine check(label, actual, expected)
        character(len=*), intent(in) :: label
        real(c_double), intent(in) :: actual, expected
        if (abs(actual - expected) <= 1.0d-14) return
        print *, 'FAIL ', label, ': expected', expected, 'got', actual
        failures = failures + 1
    end subroutine check
end program test_lorentz_torque
