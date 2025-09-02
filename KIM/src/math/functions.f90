module functions_m

    contains

    function varphi_l(x, x_lm1, x_l, x_lp1) result(phi)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), intent(in) :: x        ! Evaluation point
        real(dp), intent(in) :: x_lm1    ! x_{l-1}
        real(dp), intent(in) :: x_l      ! x_l
        real(dp), intent(in) :: x_lp1    ! x_{l+1}
        real(dp) :: h_lm1, h_l
        real(dp) :: phi

        h_lm1 = x_l - x_lm1
        h_l   = x_lp1 - x_l

        if (x >= x_lm1 .and. x < x_l) then
            phi = (x - x_lm1) / h_lm1
        else if (x >= x_l .and. x <= x_lp1) then  ! Changed < to <= to include boundary
            phi = (x_lp1 - x) / h_l
        else
            phi = 0.0d0
        end if

    end function varphi_l

    function dvarphi_l_dx(x, x_lm1, x_l, x_lp1) result(dphi)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), intent(in) :: x        ! Evaluation point
        real(dp), intent(in) :: x_lm1    ! r_{l-1}
        real(dp), intent(in) :: x_l      ! r_l
        real(dp), intent(in) :: x_lp1    ! r_{l+1}
        real(dp) :: h_lm1, h_l
        real(dp) :: dphi

        h_lm1 = x_l - x_lm1
        h_l   = x_lp1 - x_l

        if (x >= x_lm1 .and. x < x_l) then
            dphi = 1.0d0 / h_lm1
        else if (x >= x_l .and. x <= x_lp1) then  ! Changed < to <= to include boundary
            dphi = -1.0d0 / h_l
        else
            dphi = 0.0d0
        end if

    end function dvarphi_l_dx

end module