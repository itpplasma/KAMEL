!evals density f0 parameter for given density moment n0 and other f0 parameters

subroutine dens_par(n, res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)

real(dp), intent(in) :: n
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

real(dp) :: vE_

vE_ = q_ * dPhi0_ /m_ / omc_

res = n * exp(- 0.5 * (vE_ / vT_)**2)

end subroutine dens_par
