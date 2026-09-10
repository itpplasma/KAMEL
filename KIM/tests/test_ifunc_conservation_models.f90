program test_ifunc_conservation_models

    use getIfunc_config_m, only: boole_energy_conservation
    use KIM_kinds_m, only: dp

    implicit none

    interface
        subroutine getIfunc(x1, x2, symbI)
            import dp
            real(dp), intent(in) :: x1, x2
            complex(dp), intent(out) :: symbI(0:3, 0:3)
        end subroutine getIfunc

        subroutine getIfunc_model(x1, x2, conservation_model, symbI)
            import dp
            real(dp), intent(in) :: x1, x2
            integer, intent(in) :: conservation_model
            complex(dp), intent(out) :: symbI(0:3, 0:3)
        end subroutine getIfunc_model

        subroutine W2_arr(x1, x2, Imn)
            import dp
            real(dp), intent(in) :: x1, x2
            complex(dp), intent(out) :: Imn(0:3, 0:3)
        end subroutine W2_arr
    end interface

    real(dp), parameter :: tolerance = 5.0e-13_dp
    real(dp), parameter :: x1 = 0.73_dp
    real(dp), parameter :: x2 = -1.21_dp
    complex(dp) :: actual(0:3, 0:3), expected(0:3, 0:3)
    complex(dp) :: imn(0:3, 0:3), legacy(0:3, 0:3)
    integer :: conservation_model

    call W2_arr(x1, x2, imn)

    do conservation_model = 0, 3
        call expected_model(imn, conservation_model, expected)
        call getIfunc_model(x1, x2, conservation_model, actual)
        call assert_matrix_close(actual, expected, tolerance, conservation_model)
    end do

    boole_energy_conservation = .false.
    call getIfunc(x1, x2, legacy)
    call getIfunc_model(x1, x2, 0, actual)
    call assert_matrix_close(actual, legacy, tolerance, 0)

    boole_energy_conservation = .true.
    call getIfunc(x1, x2, legacy)
    call getIfunc_model(x1, x2, 1, actual)
    call assert_matrix_close(actual, legacy, tolerance, 1)

    print *, 'PASS: I-function conservation models match the KiLCA FPGEN formulas'

contains

    subroutine expected_model(bare, model, result)

        complex(dp), intent(in) :: bare(0:3, 0:3)
        integer, intent(in) :: model
        complex(dp), intent(out) :: result(0:3, 0:3)

        complex(dp) :: energy_denom, momentum_denom, coupling, determinant
        integer :: m, n

        energy_denom = 1.0_dp - bare(0, 0) + 2.0_dp * bare(2, 0) - bare(2, 2)
        momentum_denom = 1.0_dp - bare(1, 1)
        coupling = bare(1, 2) - bare(1, 0)
        determinant = energy_denom * momentum_denom - coupling**2

        select case (model)
        case (0)
            result = bare
        case (1)
            do m = 0, 3
                do n = 0, 3
                    result(m, n) = bare(m, n) &
                        + (bare(m, 0) - bare(m, 2)) &
                        * (bare(n, 0) - bare(n, 2)) / energy_denom
                end do
            end do
        case (2)
            do m = 0, 3
                do n = 0, 3
                    result(m, n) = bare(m, n) &
                        + bare(m, 1) * bare(n, 1) / momentum_denom
                end do
            end do
        case (3)
            do m = 0, 3
                do n = 0, 3
                    result(m, n) = bare(m, n) + ( &
                        coupling * bare(m, 1) * (bare(n, 2) - bare(n, 0)) &
                        + energy_denom * bare(m, 1) * bare(n, 1) &
                        + momentum_denom * (bare(m, 2) - bare(m, 0)) &
                        * (bare(n, 2) - bare(n, 0)) &
                        + coupling * (bare(m, 2) - bare(m, 0)) * bare(n, 1) &
                        ) / determinant
                end do
            end do
        case default
            error stop 'invalid model in expected_model'
        end select

    end subroutine expected_model

    subroutine assert_matrix_close(actual, expected, tol, model)

        complex(dp), intent(in) :: actual(0:3, 0:3), expected(0:3, 0:3)
        real(dp), intent(in) :: tol
        integer, intent(in) :: model

        real(dp) :: error

        error = maxval(abs(actual - expected) &
            / max(abs(expected), tiny(1.0_dp)))
        if (error > tol) then
            print *, 'FAIL: conservation model ', model, ' relative error = ', error
            error stop
        end if

    end subroutine assert_matrix_close

end program test_ifunc_conservation_models
