subroutine integrate_parallel_current(dim_r, r, Jpe, Jpi, r_res_in, Ipar)
!
!  Integrate (Jpi + Jpe) * r * dr over a Gaussian-fitted layer around
!  the resonant surface to compute the cylindrical shielding current Ipar.
!
!  KiLCA convention: total current = ions + electrons (charge sign is
!  already built into the conductivity tensor, so species add directly).
!  This matches KiLCA's internal spec=2 = spec=0 + spec=1.
!
!  Method (matches KiLCA postprocessor):
!    1. Compute Jpar_net = Jpi + Jpe
!    2. Find peak |Jpar_net| within ±5 cm of r_res
!    3. Fit Gaussian to |Jpar_net| via log-linearization
!    4. Integration range = r_res ± d/2 where d = 5*sigma
!    5. Trapezoidal integration within that range
!
    use QLBalance_kinds, only: dp

    implicit none

    integer, intent(in) :: dim_r
    real(dp), dimension(dim_r), intent(in) :: r
    complex(dp), dimension(dim_r), intent(in) :: Jpe, Jpi
    real(dp), intent(in) :: r_res_in
    complex(dp), intent(out) :: Ipar

    ! Local variables
    complex(dp), dimension(dim_r) :: Jpar_net
    real(dp), dimension(dim_r) :: absJ
    real(dp) :: search_half_width, threshold, peak_val
    real(dp) :: sigma, d_layer, r_lo, r_hi
    integer :: ipeak, ipoi
    integer :: ind_lo, ind_hi

    ! Gaussian fit variables
    integer :: nfit
    real(dp) :: xi, yi, sum_x, sum_y, sum_xx, sum_xy, denom, a1
    real(dp), parameter :: sigma_min = 0.01d0   ! cm
    real(dp), parameter :: sigma_max = 5.0d0    ! cm
    real(dp), parameter :: d_min = 0.1d0        ! cm
    real(dp), parameter :: d_max = 10.0d0       ! cm
    real(dp), parameter :: fallback_half = 2.5d0 ! cm fallback half-width
    logical :: fit_ok

    ! Step 1: Compute total parallel current density (KiLCA: J_total = J_i + J_e)
    Jpar_net = Jpi + Jpe
    absJ = abs(Jpar_net)

    ! Step 2: Find peak |Jpar_net| within ±5 cm of r_res
    search_half_width = 5.0d0
    ipeak = 0
    peak_val = 0.0d0
    do ipoi = 1, dim_r
        if (abs(r(ipoi) - r_res_in) <= search_half_width) then
            if (absJ(ipoi) > peak_val) then
                peak_val = absJ(ipoi)
                ipeak = ipoi
            end if
        end if
    end do

    ! If no peak found or peak is zero, return zero
    if (ipeak == 0 .or. peak_val < 1.0d-30) then
        Ipar = (0.0d0, 0.0d0)
        write(*,*) '[integrate_parallel_current] WARNING: no significant Jpar peak found near r_res=', r_res_in
        return
    end if

    ! Step 3: Gaussian fit via log-linearization
    !   log(|Jpar(r)|) = a0 + a1*(r - r_res)^2
    !   where a1 = -1/(2*sigma^2)
    threshold = 0.1d0 * peak_val

    ! Count fitting points
    nfit = 0
    do ipoi = 1, dim_r
        if (abs(r(ipoi) - r_res_in) <= search_half_width .and. absJ(ipoi) > threshold) then
            nfit = nfit + 1
        end if
    end do

    fit_ok = .false.
    if (nfit >= 3) then
        ! Linear least squares: y = a0 + a1*x where x=(r-r_res)^2, y=log(|Jpar|)
        sum_x  = 0.0d0
        sum_y  = 0.0d0
        sum_xx = 0.0d0
        sum_xy = 0.0d0

        do ipoi = 1, dim_r
            if (abs(r(ipoi) - r_res_in) <= search_half_width .and. absJ(ipoi) > threshold) then
                xi = (r(ipoi) - r_res_in)**2
                yi = log(absJ(ipoi))
                sum_x  = sum_x  + xi
                sum_y  = sum_y  + yi
                sum_xx = sum_xx + xi * xi
                sum_xy = sum_xy + xi * yi
            end if
        end do

        denom = dble(nfit) * sum_xx - sum_x * sum_x
        if (abs(denom) > 1.0d-30) then
            a1 = (dble(nfit) * sum_xy - sum_x * sum_y) / denom
            if (a1 < 0.0d0) then
                sigma = sqrt(-1.0d0 / (2.0d0 * a1))
                if (sigma >= sigma_min .and. sigma <= sigma_max) then
                    fit_ok = .true.
                end if
            end if
        end if
    end if

    ! Step 4: Set integration range
    if (fit_ok) then
        d_layer = 5.0d0 * sigma
        ! Clamp to safety bounds
        d_layer = max(d_min, min(d_max, d_layer))
        r_lo = r_res_in - d_layer / 2.0d0
        r_hi = r_res_in + d_layer / 2.0d0

        write(*,*) '[integrate_parallel_current] Gaussian fit succeeded:'
        write(*,*) '  r_res   = ', r_res_in, ' cm'
        write(*,*) '  sigma   = ', sigma, ' cm'
        write(*,*) '  d=5sig  = ', d_layer, ' cm'
        write(*,*) '  range   = [', r_lo, ',', r_hi, '] cm'
        write(*,*) '  nfit    = ', nfit
    else
        ! Fallback: fixed ±2.5 cm window
        r_lo = r_res_in - fallback_half
        r_hi = r_res_in + fallback_half

        write(*,*) '[integrate_parallel_current] WARNING: Gaussian fit failed, using fallback +-2.5 cm'
        write(*,*) '  r_res   = ', r_res_in, ' cm'
        write(*,*) '  nfit    = ', nfit
        if (nfit >= 3) write(*,*) '  a1      = ', a1
        write(*,*) '  range   = [', r_lo, ',', r_hi, '] cm'
    end if

    ! Clamp to grid bounds
    r_lo = max(r_lo, r(1))
    r_hi = min(r_hi, r(dim_r))

    ! Find index bounds
    ind_lo = 1
    ind_hi = dim_r
    do ipoi = 1, dim_r
        if (r(ipoi) < r_lo) then
            ind_lo = ipoi
        end if
        if (r(ipoi) > r_hi) then
            ind_hi = ipoi
            exit
        end if
    end do

    ! Step 5: Trapezoidal integration of (Jpi + Jpe) * r * dr within bounds
    Ipar = (0.0d0, 0.0d0)
    do ipoi = ind_lo, ind_hi - 1
        Ipar = Ipar + (Jpar_net(ipoi) * r(ipoi) + Jpar_net(ipoi+1) * r(ipoi+1)) &
                     * (r(ipoi+1) - r(ipoi)) / 2.0d0
    end do

    ! Step 6: Diagnostic output
    write(*,*) '  |Ipar|  = ', abs(Ipar)
    write(*,*) '  Re(Ipar)= ', real(Ipar), ' Im(Ipar)= ', aimag(Ipar)

end subroutine
