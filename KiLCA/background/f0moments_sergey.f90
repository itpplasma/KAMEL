
!Here is an automatically generated fortran subroutine:

subroutine dens_mom(res) !<f_0 moment: number density

use constants, only: dp, vp
use conduct_parameters

implicit real(vp) (s-t)
real(dp), intent(out) :: res
real(vp), dimension(1) :: val

real(dp) :: vE_

vE_ = q_ * dPhi0_ /m_ / omc_

res = n_ * exp(0.5 * (vE_ / vT_)**2)

end subroutine dens_mom

!Here is an automatically generated fortran subroutine:

subroutine Vth_mom(res) !<f_0 moment: V_theta velocity

use constants, only: dp, vp
use conduct_parameters

implicit real(vp) (s-t)
real(dp), intent(out) :: res
real(vp), dimension(1) :: val

      t5 = 1/r_
      t7 = omc_**2
      t9 = vT_**2
      t12 = 1/omc_
      val(1) = n_*(-t9/t7*t5*(hz_*domc_+ht_*dpsi_*omc_)+2*vT_*t5*dvT_*hz_ &
    *t12+Vp_*t5*ht_)+t5*dn_*hz_*t12*t9+t5*t12*hz_/m_*n_*q_*dPhi0_

res = val(1)

end subroutine Vth_mom

!Here is an automatically generated fortran subroutine:

subroutine Vz_mom(res) !<f_0 moment: V_z velocity

use constants, only: dp, vp
use conduct_parameters

implicit real(vp) (s-t)
real(dp), intent(out) :: res
real(vp), dimension(1) :: val

      t8 = omc_**2
      t12 = vT_**2
      t15 = 1/omc_
      val(1) = n_*(t12/r_/t8*(ht_*domc_*r_-omc_*hz_*dpsi_*r_-omc_*ht_)-2 &
    *dvT_*vT_*ht_*t15+hz_*Vp_)-ht_*dn_*t15*t12-t15/m_*n_*q_*dPhi0_*ht_

res = val(1)

end subroutine Vz_mom

!Here is an automatically generated fortran subroutine:

subroutine Vs_mom(res) !<f_0 moment: V_s velocity

use constants, only: dp, vp
use conduct_parameters

implicit real(vp) (s-t)
real(dp), intent(out) :: res
real(vp), dimension(1) :: val

      t1 = hz_**2
      t7 = omc_**2
      t9 = vT_**2
      t12 = 1/omc_
      val(1) = n_*(-t9/t7/r_*(-omc_+t1*omc_+domc_*r_)+2*vT_*t12*dvT_)+t12 &
    *t9*dn_+1/m_*t12*dPhi0_*n_*q_

res = val(1)

end subroutine Vs_mom

!Here is an automatically generated fortran subroutine:

subroutine Vp_mom(res) !<f_0 moment: V_p velocity

use constants, only: dp, vp
use conduct_parameters

implicit real(vp) (s-t)
real(dp), intent(out) :: res
real(vp), dimension(1) :: val

      t7 = vT_**2
      val(1) = n_*(-t7/omc_*(dpsi_*r_+hz_*ht_)/r_+Vp_)

res = val(1)

end subroutine Vp_mom

!Here is an automatically generated fortran subroutine:

subroutine Etot_mom(res) !<f_0 moment: total energy

use constants, only: dp, vp
use conduct_parameters

implicit real(vp) (s-t)
real(dp), intent(out) :: res
real(vp), dimension(1) :: val

      t1 = vT_**2
      t3 = Vp_**2
      val(1) = m_*n_*(3.D0/2.D0*t1+t3/2)

res = val(1)

end subroutine Etot_mom

!Here is an automatically generated fortran subroutine:

subroutine Eterm_mom(res) !<f_0 moment: termal energy

use constants, only: dp, vp
use conduct_parameters

implicit real(vp) (s-t)
real(dp), intent(out) :: res
real(vp), dimension(1) :: val

      t2 = vT_**2
      val(1) = 3.D0/2.D0*t2*m_*n_

res = val(1)

end subroutine Eterm_mom
