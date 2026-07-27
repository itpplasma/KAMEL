program test_ion_temperature_gradient_model

    use KIM_kinds_m, only: dp
    use config_m, only: ion_temperature_gradient_model_is_valid, &
                        temperature_gradient_force_terms

    implicit none

    real(dp), parameter :: gradient = 0.4_dp
    real(dp) :: a1_temperature, a2

    call assert_valid('full', .true.)
    call assert_valid('zero_A2', .true.)
    call assert_valid('zero_Tprime', .true.)
    call assert_valid('unknown', .false.)

    call temperature_gradient_force_terms(0, 'zero_Tprime', gradient, &
        a1_temperature, a2)
    call assert_close(a1_temperature, -0.6_dp, 'electron A1 is unchanged')
    call assert_close(a2, 0.4_dp, 'electron A2 is unchanged')

    call temperature_gradient_force_terms(1, 'full', gradient, &
        a1_temperature, a2)
    call assert_close(a1_temperature, -0.6_dp, 'full ion A1')
    call assert_close(a2, 0.4_dp, 'full ion A2')

    call temperature_gradient_force_terms(1, 'zero_A2', gradient, &
        a1_temperature, a2)
    call assert_close(a1_temperature, -0.6_dp, 'zero_A2 ion A1')
    call assert_close(a2, 0.0_dp, 'zero_A2 ion A2')

    call temperature_gradient_force_terms(1, 'zero_Tprime', gradient, &
        a1_temperature, a2)
    call assert_close(a1_temperature, 0.0_dp, 'zero_Tprime ion A1')
    call assert_close(a2, 0.0_dp, 'zero_Tprime ion A2')

    print *, 'PASS: ion temperature-gradient force models'

contains

    subroutine assert_valid(model, expected)

        character(*), intent(in) :: model
        logical, intent(in) :: expected

        if (ion_temperature_gradient_model_is_valid(model) .neqv. expected) then
            error stop 'ion temperature-gradient model validation mismatch'
        end if

    end subroutine assert_valid

    subroutine assert_close(actual, expected, label)

        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label

        if (abs(actual - expected) > 10.0_dp * epsilon(1.0_dp) &
                * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label), ': expected ', expected, ', got ', actual
            error stop
        end if

    end subroutine assert_close

end program test_ion_temperature_gradient_model
