module wavenumber_geometry_m

    use KIM_kinds_m, only: dp
    use constants_m, only: sol

    implicit none
    private

    public :: parallel_wavenumber
    public :: perpendicular_wavenumber
    public :: exb_rotation_frequency

contains

    pure elemental real(dp) function parallel_wavenumber(m_mode, n_mode, r, R0, &
                                                         hth, hz) result(kp)
        integer, intent(in) :: m_mode, n_mode
        real(dp), intent(in) :: r, R0, hth, hz

        kp = hth * real(m_mode, dp) / r + hz * real(n_mode, dp) / R0
    end function parallel_wavenumber

    pure elemental real(dp) function perpendicular_wavenumber(m_mode, n_mode, r, &
                                                              R0, hth, hz) result(ks)
        integer, intent(in) :: m_mode, n_mode
        real(dp), intent(in) :: r, R0, hth, hz

        ks = hz * real(m_mode, dp) / r - hth * real(n_mode, dp) / R0
    end function perpendicular_wavenumber

    pure elemental real(dp) function exb_rotation_frequency(ks, Er, B0) result(omega_E)
        real(dp), intent(in) :: ks, Er, B0

        omega_E = -sol * Er * ks / B0
    end function exb_rotation_frequency

end module wavenumber_geometry_m
