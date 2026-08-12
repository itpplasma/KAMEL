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
    public :: compact_transition

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
