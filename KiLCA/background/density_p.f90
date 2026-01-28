!evals density f0 parameter for given density moment n0 and other f0 parameters

subroutine dens_par(n, res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)

real(dp), intent(in) :: n
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t7 = 1/r_
      t13 = hz_**2
      t16 = omc_**2
      val(1) = 1/omc_*t7*(r_*omc_+Vp_*dpsi_*r_+Vp_*hz_*ht_)+q_/m_ &
    *(dPhi0_*t7/t16/omc_*(-domc_*r_+omc_*t13)+1/t16*ddPhi0_)
    res = n/val(1)

end subroutine dens_par
