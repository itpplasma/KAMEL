module periodic_current_normalization_m
    !! Target-current normalization for one linear periodic response.
    !! integrate_trusted_current returns int(J_parallel*r dr); the scale
    !! routine applies the documented 2*pi cylindrical-current convention.
    use QLBalance_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: integrate_trusted_current, periodic_drive_scale

contains

    complex(dp) function integrate_trusted_current(r, jpar, core_lo, core_hi) result(current)
        real(dp), intent(in) :: r(:), core_lo, core_hi
        complex(dp), intent(in) :: jpar(:)
        integer :: i
        real(dp) :: left, right, cell_width, left_fraction, right_fraction
        complex(dp) :: integrand_left, integrand_right, clipped_left, clipped_right
        if (size(r) /= size(jpar) .or. size(r) < 2) error stop 'current integration grid mismatch'
        if (.not. all(ieee_is_finite(r)) .or. .not. ieee_is_finite(core_lo) .or. &
            .not. ieee_is_finite(core_hi)) error stop 'current integration nonfinite geometry'
        if (any(r(2:) <= r(:size(r) - 1))) error stop 'current integration grid must increase'
        if (core_lo >= core_hi .or. core_lo < r(1) .or. core_hi > r(size(r))) &
            error stop 'current integration core must lie within sampled grid'
        ! A failed physical response is handled by the caller's normalization
        ! guard. Geometry errors remain hard errors instead of silent clipping.
        if (.not. all(ieee_is_finite(real(jpar, dp))) .or. &
            .not. all(ieee_is_finite(aimag(jpar)))) then
            current = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
            return
        end if
        current = (0.0_dp, 0.0_dp)
        do i = 1, size(r) - 1
            left = max(r(i), core_lo)
            right = min(r(i + 1), core_hi)
            if (right <= left) cycle
            ! Interpolate the cylindrical integrand r*J, consistently with
            ! the full-cell trapezoid, before integrating a clipped interval.
            cell_width = r(i + 1) - r(i)
            left_fraction = (left - r(i)) / cell_width
            right_fraction = (right - r(i)) / cell_width
            integrand_left = r(i) * jpar(i)
            integrand_right = r(i + 1) * jpar(i + 1)
            clipped_left = (1.0_dp - left_fraction) * integrand_left &
                + left_fraction * integrand_right
            clipped_right = (1.0_dp - right_fraction) * integrand_left &
                + right_fraction * integrand_right
            current = current + 0.5_dp * (clipped_left + clipped_right) * (right - left)
        end do
    end function integrate_trusted_current

    subroutine periodic_drive_scale(target_current, unit_current, c_light, current_floor, &
                                    max_scale_ratio, relaxation, scale, status)
        !! status = 0: normalized response
        !! status = 1: invalid finite configuration (response suppressed)
        !! status = 2: unit current below the trusted floor (response suppressed)
        !! status = 3: non-finite or excessive scale (response suppressed)
        !! Relaxation mixes the desired scale with the unit scale for this
        !! one response; it is not an update from a previous time step.
        real(dp), intent(in) :: target_current, c_light, current_floor, max_scale_ratio, relaxation
        complex(dp), intent(in) :: unit_current
        complex(dp), intent(out) :: scale
        integer, intent(out) :: status
        real(dp) :: ratio
        status = 3
        scale = (0.0_dp, 0.0_dp)
        if (.not. all(ieee_is_finite([target_current, c_light, current_floor, &
                 max_scale_ratio, relaxation, real(unit_current, dp), aimag(unit_current)]))) return
        if (target_current <= 0.0_dp .or. c_light <= 0.0_dp .or. current_floor <= 0.0_dp &
            .or. max_scale_ratio <= 0.0_dp .or. relaxation <= 0.0_dp .or. relaxation > 1.0_dp) then
            status = 1
            return
        end if
        if (abs(unit_current) <= current_floor) then
            status = 2
            return
        end if
        scale = relaxation * (target_current * c_light / unit_current) &
                / (2.0_dp * acos(-1.0_dp)) &
                + (1.0_dp - relaxation) * (1.0_dp, 0.0_dp)
        if (.not. ieee_is_finite(real(scale, dp)) .or. .not. ieee_is_finite(aimag(scale))) then
            scale = (0.0_dp, 0.0_dp)
            return
        end if
        ratio = abs(scale)
        if (ratio > max_scale_ratio) then
            scale = (0.0_dp, 0.0_dp)
            return
        end if
        status = 0
    end subroutine periodic_drive_scale

end module periodic_current_normalization_m
