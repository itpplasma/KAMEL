module periodic_current_normalization_m
    !! Target-current normalization for one linear periodic response.
    !! integrate_trusted_current returns int(J_parallel*r dr); the scale
    !! routine applies the documented 2*pi cylindrical-current convention.
    use QLBalance_kinds, only: dp
    implicit none
    private
    public :: integrate_trusted_current, periodic_drive_scale

contains

    complex(dp) function integrate_trusted_current(r, jpar, core_lo, core_hi) result(current)
        real(dp), intent(in) :: r(:), core_lo, core_hi
        complex(dp), intent(in) :: jpar(:)
        integer :: i
        real(dp) :: left, right, midpoint
        if (size(r) /= size(jpar) .or. size(r) < 2) error stop 'current integration grid mismatch'
        current = (0.0_dp, 0.0_dp)
        do i = 1, size(r)-1
            midpoint = 0.5_dp*(r(i)+r(i+1))
            if (midpoint < core_lo .or. midpoint > core_hi) cycle
            left = max(r(i), core_lo)
            right = min(r(i+1), core_hi)
            if (right > left) current = current + 0.5_dp*(r(i+1)*jpar(i+1)+r(i)*jpar(i))*(right-left)
        end do
    end function integrate_trusted_current

    subroutine periodic_drive_scale(target_current, unit_current, c_light, current_floor, &
                                    max_scale_ratio, relaxation, scale, status)
        real(dp), intent(in) :: target_current, c_light, current_floor, max_scale_ratio, relaxation
        complex(dp), intent(in) :: unit_current
        complex(dp), intent(out) :: scale
        integer, intent(out) :: status
        real(dp) :: ratio
        status = 0
        scale = (1.0_dp, 0.0_dp)
        if (target_current <= 0.0_dp .or. c_light <= 0.0_dp .or. current_floor <= 0.0_dp &
            .or. max_scale_ratio <= 0.0_dp .or. relaxation <= 0.0_dp .or. relaxation > 1.0_dp) then
            status = 1
            return
        end if
        if (abs(unit_current) <= current_floor) then
            status = 2
            return
        end if
        scale = relaxation * (target_current*c_light/(unit_current)) / (2.0_dp*acos(-1.0_dp)) &
            + (1.0_dp-relaxation)*(1.0_dp,0.0_dp)
        ratio = abs(scale)
        if (.not. (ratio == ratio) .or. ratio > max_scale_ratio) then
            status = 3
            scale = (1.0_dp, 0.0_dp)
        end if
    end subroutine periodic_drive_scale

end module periodic_current_normalization_m
