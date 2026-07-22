program probe_lorentz_torque_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use kilca_flre_quants_m, only: integrate_over_cylinder, vec_product_3d

    implicit none

    real(c_double), parameter :: pi = 3.141592653589793238462643383279502884197d0
    real(c_double), parameter :: cspeed = 2.9979245800d10
    real(c_double), parameter :: echarge = 4.8032d-10
    real(c_double), parameter :: rtor = 3.0d0
    real(c_double) :: force(0:2), torque_z, expected_total
    real(c_double) :: radius(3), torque_density(3), integrated(3)
    complex(c_double) :: current(0:2), magnetic_conjugate(0:2), cross(0:2)
    integer(c_int), parameter :: dim = 3

    current = [cmplx(1.0d0, 2.0d0, c_double), &
        cmplx(-0.5d0, 0.25d0, c_double), &
        cmplx(0.75d0, -1.25d0, c_double)]
    magnetic_conjugate = conjg([cmplx(0.3d0, -0.4d0, c_double), &
        cmplx(-1.1d0, 0.2d0, c_double), &
        cmplx(0.6d0, 0.9d0, c_double)])
    call vec_product_3d(current, magnetic_conjugate, cross)
    force = 0.5d0*echarge/cspeed*real(cross, c_double)

    torque_z = rtor*force(2)
    radius = [0.0d0, 0.5d0, 1.0d0]
    torque_density = torque_z
    call integrate_over_cylinder(dim, radius, torque_density, &
        4.0d0*pi*pi*rtor, integrated)
    expected_total = 4.0d0*pi*pi*rtor*0.5d0*torque_z

    write (*, '(a,3(1x,es24.16))') 'KILCA_FORCE', force
    write (*, '(a,1x,es24.16)') 'KILCA_TORQUE_Z', torque_z
    write (*, '(a,1x,es24.16)') 'KILCA_TOTAL', integrated(dim)
    write (*, '(a,1x,es24.16)') 'KILCA_EXPECTED', expected_total
end program probe_lorentz_torque_kernel
