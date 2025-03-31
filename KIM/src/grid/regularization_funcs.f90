module regularization_funcs

    use KIM_kinds, only: dp
    use resonances_mod, only: r_res

    implicit none

    real(dp) :: gauss_width = 0.1d0

    contains

    real(dp) function theta_middle(r)

        implicit none
        real(dp), intent(in) :: r

        theta_middle = exp(-gauss_width *(r - r_res)**2.0d0)

    end function


end module