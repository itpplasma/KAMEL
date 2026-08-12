module periodic_embedding_m
    !! Compact local-to-global transition for periodic KIM responses.
    !!
    !! The returned weight is one on the trusted core, decreases smoothly
    !! across each transition band, and is exactly zero outside bounded
    !! support.  Physical fields must be evaluated before multiplying by
    !! this weight; no derivative of the window is part of the contract.
    use QLBalance_kinds, only: dp
    implicit none
    private
    public :: compact_transition, embed_complex_profile, embed_tensor_profile

contains

    real(dp) function compact_transition(r, core_lo, core_hi, width) result(w)
        real(dp), intent(in) :: r, core_lo, core_hi, width
        real(dp) :: t

        if (width <= 0.0_dp .or. core_hi < core_lo) then
            error stop 'compact_transition: invalid core or transition width'
        end if

        if (r < core_lo - width .or. r > core_hi + width) then
            w = 0.0_dp
        elseif (r >= core_lo .and. r <= core_hi) then
            w = 1.0_dp
        elseif (r < core_lo) then
            t = (r - (core_lo - width)) / width
            w = smooth_ramp(t)
        else
            t = ((core_hi + width) - r) / width
            w = smooth_ramp(t)
        end if
    end function compact_transition

    subroutine embed_complex_profile(local_r, local_z, global_r, core_lo, core_hi, width, &
                                     global_z, weights)
        !! Embed a locally computed physical complex field on a global grid.
        !! The field is interpolated before the common compact transition is
        !! applied; no derivative of the transition is ever introduced.
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(dp), intent(in) :: local_r(:), global_r(:)
        complex(dp), intent(in) :: local_z(:)
        real(dp), intent(in) :: core_lo, core_hi, width
        complex(dp), intent(out) :: global_z(:)
        real(dp), intent(out) :: weights(:)
        integer :: i
        real(dp) :: re, im

        call validate_embedding_grids(local_r, global_r, core_lo, core_hi, width)
        if (size(local_z) /= size(local_r) .or. size(global_z) /= size(global_r) &
            .or. size(weights) /= size(global_r)) then
            error stop 'embed_complex_profile: inconsistent array sizes'
        end if
        do i = 1, size(global_r)
            weights(i) = compact_transition(global_r(i), core_lo, core_hi, width)
            if (weights(i) == 0.0_dp) then
                global_z(i) = (0.0_dp, 0.0_dp)
            else
                re = linear_profile(local_r, real(local_z), global_r(i))
                im = linear_profile(local_r, aimag(local_z), global_r(i))
                if (.not. ieee_is_finite(re) .or. .not. ieee_is_finite(im)) &
                    error stop 'embed_complex_profile: non-finite local field'
                global_z(i) = weights(i) * cmplx(re, im, dp)
            end if
        end do
    end subroutine embed_complex_profile

    subroutine embed_tensor_profile(local_r, local_d, global_r, core_lo, core_hi, width, global_d)
        !! Embed a local real 2x2 tensor.  Since the tensor is quadratic in
        !! perturbation amplitude, the physical field transition enters as w^2.
        real(dp), intent(in) :: local_r(:), local_d(:,:,:), global_r(:)
        real(dp), intent(in) :: core_lo, core_hi, width
        real(dp), intent(out) :: global_d(:,:,:)
        real(dp) :: w
        integer :: i, j, k

        call validate_embedding_grids(local_r, global_r, core_lo, core_hi, width)
        if (size(local_d, 1) /= 2 .or. size(local_d, 2) /= 2 .or. &
            size(local_d, 3) /= size(local_r) .or. size(global_d, 1) /= 2 .or. &
            size(global_d, 2) /= 2 .or. size(global_d, 3) /= size(global_r)) &
            error stop 'embed_tensor_profile: inconsistent tensor dimensions'
        do i = 1, size(global_r)
            w = compact_transition(global_r(i), core_lo, core_hi, width)
            do j = 1, 2
                do k = 1, 2
                    if (w == 0.0_dp) then
                        global_d(j, k, i) = 0.0_dp
                    else
                        global_d(j, k, i) = w*w*linear_profile(local_r, local_d(j, k, :), global_r(i))
                    end if
                end do
            end do
        end do
    end subroutine embed_tensor_profile

    subroutine validate_embedding_grids(local_r, global_r, core_lo, core_hi, width)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(dp), intent(in) :: local_r(:), global_r(:), core_lo, core_hi, width
        integer :: i

        if (size(local_r) < 2 .or. size(global_r) < 1) &
            error stop 'periodic embedding requires non-empty grids'
        if (width <= 0.0_dp .or. core_hi < core_lo) &
            error stop 'periodic embedding has invalid core or width'
        if (local_r(1) > core_lo - width .or. local_r(size(local_r)) < core_hi + width) &
            error stop 'periodic embedding local grid does not cover compact support'
        do i = 2, size(local_r)
            if (.not. ieee_is_finite(local_r(i)) .or. local_r(i) <= local_r(i-1)) &
                error stop 'periodic embedding local grid must be finite and increasing'
        end do
    end subroutine validate_embedding_grids

    real(dp) function linear_profile(x, y, xq) result(value)
        real(dp), intent(in) :: x(:), y(:), xq
        integer :: i
        real(dp) :: t

        if (size(x) /= size(y)) error stop 'periodic embedding interpolation size mismatch'
        if (xq <= x(1)) then
            value = y(1)
            return
        end if
        if (xq >= x(size(x))) then
            value = y(size(y))
            return
        end if
        do i = 1, size(x)-1
            if (xq <= x(i+1)) then
                t = (xq-x(i))/(x(i+1)-x(i))
                value = (1.0_dp-t)*y(i) + t*y(i+1)
                return
            end if
        end do
        error stop 'periodic embedding interpolation failed'
    end function linear_profile

    pure real(dp) function smooth_ramp(t) result(w)
        real(dp), intent(in) :: t
        real(dp) :: x

        x = min(1.0_dp, max(0.0_dp, t))
        if (x <= 0.0_dp) then
            w = 0.0_dp
        elseif (x >= 1.0_dp) then
            w = 1.0_dp
        else
            ! Standard C-infinity bump transition, normalized to endpoints.
            w = exp(-1.0_dp/x)
            w = w / (w + exp(-1.0_dp/(1.0_dp-x)))
        end if
    end function smooth_ramp

end module periodic_embedding_m
