!> Regression test for the Fortran vel_integral port (vel_integral_m).
!>
!> Reference values are the original vel_integral.cpp output (same fortnum
!> semi-infinite quadrature) for the three index branches; the port reproduces
!> them exactly.
program test_vel_integral
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use vel_integral_m, only: calc_velocity_integral
    implicit none

    real(c_double), parameter :: tol = 1.0d-10
    real(c_double) :: Es(2), Ep(2), res(1)
    real(c_double) :: ref(3)
    integer :: ind, failures

    Es = [0.7d0, -0.3d0]
    Ep = [0.2d0, 0.5d0]
    ref = [4.9943549094035191d0, 6.6097086357211552d0, 1.5089766708962472d1]

    failures = 0
    do ind = 1, 3
        res = 0.0d0
        call calc_velocity_integral(ind, 1.3d0, 0.8d0, 0.5d0, 0.25d0, 0.1d0, Es, Ep, res)
        if (abs(res(1) - ref(ind)) > tol) then
            write (*, '(a,i0,2(1x,es23.16))') "FAIL ind=", ind, res(1), ref(ind)
            failures = failures + 1
        end if
    end do

    if (failures == 0) then
        write (*, '(a)') "PASS: vel_integral_m matches reference values"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if
end program test_vel_integral
