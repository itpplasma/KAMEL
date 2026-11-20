!> The fortnum-backed 1F1 evaluations must not return a value fortnum flagged as invalid.
!> At a pole in b (b+2 = -1 for the modified form, b = -1 for M itself) fortnum reports a
!> domain error; the KiLCA wrapper has to stop with the routine name, b and z instead of
!> passing the result on. The routine under test is selected by the first argument.
program test_hyper1f1_status
    use kilca_hyper1f1_m, only: h_cont_fract_1_modified_0_ada, h_quad
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none

    character(len=32) :: which
    real(c_double) :: f_re, f_im

    call get_command_argument(1, which)
    select case (trim(which))
    case ('modified')
        call h_cont_fract_1_modified_0_ada(-3.0_c_double, 0.0_c_double, 1.0_c_double, &
                                           0.0_c_double, f_re, f_im)
    case ('quad')
        call h_quad(-1.0_c_double, 0.0_c_double, 1.0_c_double, 0.0_c_double, f_re, f_im)
    case default
        error stop 'usage: test_hyper1f1_status.x modified|quad'
    end select

    print *, 'NO_ERROR_STOP: returned', f_re, f_im
end program test_hyper1f1_status
