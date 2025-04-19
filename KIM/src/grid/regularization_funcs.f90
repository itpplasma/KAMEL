module regularization_funcs

    use KIM_kinds, only: dp
    use resonances_mod, only: r_res

    implicit none

    real(dp) :: gauss_width = 0.1d0
    real(dp) :: smooth_step_width = 5.0d0
    real(dp) :: delta_rs = 3.0d0

    contains

    real(dp) function theta_middle(r)

        implicit none
        real(dp), intent(in) :: r

        !theta_middle = exp(-gauss_width *(r - r_res)**2.0d0)
        theta_middle = theta_right(r) * theta_left(r)

    end function

    real(dp) function theta_right(r)

        use constants, only: pi

        implicit none

        real(dp), intent(in) :: r

        theta_right = atan((r - (r_res - delta_rs))/0.1d0) / pi + 0.5d0
        if (r < (r_res-delta_rs)) then
            theta_right = theta_right * exp(-0.5d0 * (r - (r_res - delta_rs))**2.0d0)
        end if

    end function


    real(dp) function theta_left(r)

        use constants, only: pi

        implicit none

        real(dp), intent(in) :: r

        theta_left = -atan((r - (r_res + delta_rs))/0.1d0) / pi + 0.5d0
        if (r > (r_res+delta_rs)) then
            theta_left = theta_left * exp(-0.5d0 * (r - (r_res + delta_rs))**2.0d0)
        end if

    end function

end module