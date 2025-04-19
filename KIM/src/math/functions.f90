module functions

    contains

    function varphi_l(r, r_lm1, r_l, r_lp1) result(phi)

        use KIM_kinds, only: dp

        implicit none

        real(dp), intent(in) :: r        ! Evaluation point
        real(dp), intent(in) :: r_lm1    ! r_{l-1}
        real(dp), intent(in) :: r_l      ! r_l
        real(dp), intent(in) :: r_lp1    ! r_{l+1}
        real(dp) :: h_lm1, h_l
        real(dp) :: phi

        h_lm1 = r_l - r_lm1
        h_l   = r_lp1 - r_l

        if (r >= r_lm1 .and. r < r_l) then
            phi = (r - r_lm1) / h_lm1
        else if (r >= r_l .and. r < r_lp1) then
            phi = (r_lp1 - r) / h_l
        else
            phi = 0.0d0
        end if
    end function varphi_l
end module