subroutine integrate_parallel_current(dim_r, r, Jpe, Jpi, Ipar)

    use QLBalance_kinds, only: dp
    use grid_mod, only: r_resonant
    use resonances_mod, only: numres

    implicit none

    integer, intent(in) :: dim_r
    real(dp), dimension(dim_r), intent(in) :: r
    complex(dp), dimension(dim_r), intent(in) :: Jpe, Jpi
    complex(dp), intent(out) :: Ipar

    real(dp) :: int_range_1, int_range_2

    integer :: imn, ipoi, ind_r1, ind_r2

    Ipar = 0.0d0
    ind_r1 = 1
    ind_r2 = dim_r
    int_range_1 = 2.5d0
    int_range_2 = 2.5d0

    do imn = 1, numres
        do ipoi = 1, dim_r
            if (r(ipoi) < r_resonant(imn) - int_range_1) then
                ind_r1 = ipoi
            else if (r(ipoi) > r_resonant(imn) + int_range_2) then
                ind_r2 = ipoi
                exit
            end if
        end do

        do ipoi=1, dim_r-1
            Ipar = Ipar + (Jpi(ipoi)-Jpe(ipoi) + Jpi(ipoi+1) - Jpe(ipoi+1)) * r(ipoi) * (r(ipoi+1) - r(ipoi)) / 2.0d0
        end do
    end do

end subroutine
