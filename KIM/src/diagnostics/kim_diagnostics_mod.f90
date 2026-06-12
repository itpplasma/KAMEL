module kim_diagnostics_m

    use KIM_kinds_m, only: dp

    implicit none

    private
    public :: integrate_Ipar

contains

    pure function integrate_Ipar(npts, r, jpar) result(Ipar)
        ! Integrated parallel current I_par = 2*pi * int(r' * jpar(r') dr')
        ! over the full passed grid (complex trapezoid rule). The caller is
        ! expected to pass the restricted layer domain.
        integer, intent(in) :: npts
        real(dp), intent(in) :: r(npts)
        complex(dp), intent(in) :: jpar(npts)
        complex(dp) :: Ipar

        real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
        integer :: i

        Ipar = (0.0_dp, 0.0_dp)
        do i = 1, npts - 1
            Ipar = Ipar + (r(i) * jpar(i) + r(i+1) * jpar(i+1)) &
                          * (r(i+1) - r(i)) / 2.0_dp
        end do
        Ipar = 2.0_dp * pi * Ipar

    end function

end module
