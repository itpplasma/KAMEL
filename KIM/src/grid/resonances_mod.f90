module kim_resonances_m

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp

    implicit none

    integer :: iunit_res
    real(dp) :: r_res
    integer :: index_rg_res
    logical :: prop = .true.
    real(dp) :: prescribed_r_res = 0.0_dp
    logical :: prescribed_r_res_active = .false.

    integer, parameter :: KIM_RESONANCE_OK = 0
    integer, parameter :: KIM_RESONANCE_INVALID_INPUT = 1
    integer, parameter :: KIM_RESONANCE_NOT_FOUND = 2
    integer, parameter :: KIM_RESONANCE_AMBIGUOUS = 3

contains

    pure subroutine locate_periodic_resonance(radius, q_profile, m_mode, &
            n_mode, resonance, status)
        real(dp), intent(in) :: radius(:), q_profile(:)
        integer, intent(in) :: m_mode, n_mode
        real(dp), intent(out) :: resonance
        integer, intent(out) :: status
        real(dp) :: q_left, q_right, q_target, weight
        integer :: i, root_count

        resonance = 0.0_dp
        if (size(radius) < 2 .or. size(q_profile) /= size(radius) .or. &
                n_mode == 0) then
            status = KIM_RESONANCE_INVALID_INPUT
            return
        end if
        if (.not. all(ieee_is_finite(radius)) .or. &
                .not. all(ieee_is_finite(q_profile)) .or. &
                any(radius(2:) <= radius(:size(radius) - 1))) then
            status = KIM_RESONANCE_INVALID_INPUT
            return
        end if

        ! KiLCA and the manuscript use the signed resonance condition
        ! q(r_res) = -m/n. Absolute values would invent a surface for a mode
        ! whose helicity is incompatible with the equilibrium field.
        q_target = -real(m_mode, dp) / real(n_mode, dp)
        root_count = 0
        do i = 1, size(radius) - 1
            q_left = q_profile(i) - q_target
            q_right = q_profile(i + 1) - q_target
            if (q_left == 0.0_dp .and. q_right == 0.0_dp) then
                resonance = 0.0_dp
                status = KIM_RESONANCE_AMBIGUOUS
                return
            end if
            if (q_left == 0.0_dp) then
                root_count = root_count + 1
                resonance = radius(i)
            else if ((q_left < 0.0_dp .and. q_right > 0.0_dp) .or. &
                    (q_left > 0.0_dp .and. q_right < 0.0_dp)) then
                root_count = root_count + 1
                weight = -q_left / (q_right - q_left)
                resonance = (1.0_dp - weight) * radius(i) + &
                    weight * radius(i + 1)
            end if
            if (root_count > 1) then
                resonance = 0.0_dp
                status = KIM_RESONANCE_AMBIGUOUS
                return
            end if
        end do
        if (q_profile(size(q_profile)) == q_target) then
            root_count = root_count + 1
            resonance = radius(size(radius))
        end if

        select case (root_count)
        case (0)
            status = KIM_RESONANCE_NOT_FOUND
        case (1)
            status = KIM_RESONANCE_OK
        case default
            resonance = 0.0_dp
            status = KIM_RESONANCE_AMBIGUOUS
        end select
    end subroutine locate_periodic_resonance

end module kim_resonances_m
