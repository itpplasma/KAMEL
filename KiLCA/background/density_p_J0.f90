!evals density f0 parameter for given density moment n0 and other f0 parameters

subroutine dens_par(n, res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)

real(dp), intent(in) :: n
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

    res = n

end subroutine dens_par
