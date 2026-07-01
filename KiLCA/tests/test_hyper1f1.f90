!> Regression test for the Fortran 1F1 port (kilca_hyper1f1_m).
!>
!> Reference values are the original hyper1F1.cpp output for the dispatcher
!> hypergeometric1f1_cont_fract_1_modified_0_ada, captured before the port, at
!> inputs that exercise both branches: |z/b| < 0.1 (Kummer series) and |z/b|
!> >= 0.1 (continued fraction). The port reproduces them to ~1e-15.
program test_hyper1f1
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none

    interface
        subroutine disp(br, bi, zr, zi, fr, fi) &
            bind(C, name="hypergeometric1f1_cont_fract_1_modified_0_ada_")
            import :: c_double
            real(c_double) :: br, bi, zr, zi, fr, fi
        end subroutine disp
    end interface

    real(c_double), parameter :: tol = 1.0d-12
    real(c_double) :: bs(2, 4), zs(2, 4), ref(2, 4)
    real(c_double) :: br, bi, zr, zi, fr, fi
    integer :: i, failures

    bs = reshape([2.0d0, 0.5d0, 1.5d0, -0.3d0, 3.0d0, 1.0d0, 0.7d0, 0.2d0], [2, 4])
    zs = reshape([0.05d0, 0.02d0, 0.8d0, 0.4d0, 2.0d0, -1.0d0, 0.1d0, 0.05d0], [2, 4])
    ref = reshape([ &
        1.3046957063426200d-02, 3.4588099438105468d-03, &
        2.4008192680132767d-01, 1.8665462901964203d-01, &
        3.5916458035555343d-01, -4.8084640665540940d-01, &
        3.9077977106348083d-02, 1.6610626441553922d-02], [2, 4])

    failures = 0
    do i = 1, 4
        br = bs(1, i); bi = bs(2, i); zr = zs(1, i); zi = zs(2, i)
        fr = 0.0d0; fi = 0.0d0
        call disp(br, bi, zr, zi, fr, fi)
        if (abs(fr - ref(1, i)) > tol .or. abs(fi - ref(2, i)) > tol) then
            write (*, '(a,i0,4(1x,es23.16))') "FAIL case ", i, fr, ref(1, i), fi, ref(2, i)
            failures = failures + 1
        end if
    end do

    if (failures == 0) then
        write (*, '(a)') "PASS: kilca_hyper1f1_m matches reference 1F1 values"
    else
        write (*, '(a,i0)') "FAILED cases: ", failures
        stop 1
    end if
end program test_hyper1f1
