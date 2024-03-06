!<Subroutines used to evaluate galilean correction term from Sergey.

!Here is an automatically generated fortran subroutine:

subroutine delC0_11(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC0_11

!Here is an automatically generated fortran subroutine:

subroutine delC1_11(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC1_11

!Here is an automatically generated fortran subroutine:

subroutine delC2_11(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_11

!Here is an automatically generated fortran subroutine:

subroutine delC0_12(res)

use constants, only: dpc, dp, vpc
use conduct_parameters
use gal_corr, only: dVeq

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      t2 = hz_**2
      t3 = 1.0d0/omega_
      t9 = vT_**2
      t10 = t9*t3
      t14 = hz_*ht_
      t15 = r_*dpsi_
      t22 = Vp_*ks_
      t25 = omega_*t2+t14*t22-t15*t22
      val(1) = m_*(n_*(kp_*(-2*vT_*t3*t2*dvT_*Vp_-t10*t2*dVeq)-t3*t9*ks_ &
	*(-t14+t15)*dVeq+2*vT_*t3*t25*dvT_)-Vp_*t9*kp_*t3*t2*dn_+t10*t25*dn_)

res = val(1)

end subroutine delC0_12

!Here is an automatically generated fortran subroutine:

subroutine delC1_12(res)

use constants, only: dpc, dp, vpc
use conduct_parameters
use gal_corr, only: dVeq

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      t2 = 1.0d0/omega_
      t3 = vT_**2
      t6 = r_*dvT_
      t17 = dn_*r_
      val(1) = m_*(n_*(kp_*(-t3*t2*r_*dVeq-2*vT_*t2*Vp_*t6)+2*vT_*t6)-Vp_ &
	*t3*kp_*t2*t17+t3*t17)

res = val(1)

end subroutine delC1_12

!Here is an automatically generated fortran subroutine:

subroutine delC2_12(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_12

!Here is an automatically generated fortran subroutine:

subroutine delC0_13(res)

use constants, only: dpc, dp, vpc
use conduct_parameters
use gal_corr, only: dVeq

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      t3 = r_*dpsi_+hz_*ht_
      t6 = 1.0d0/omega_
      t11 = vT_**2
      t12 = t6*t11
      t16 = hz_**2
      t26 = ks_*Vp_
      t28 = omega_*r_*dpsi_+omega_*hz_*ht_-t16*t26+t26
      val(1) = m_*(n_*(kp_*(-2*t6*vT_*Vp_*t3*dvT_-t12*t3*dVeq)-t6*t11*ks_ &
	*(t16-1)*dVeq+2*vT_*t6*t28*dvT_)-t6*Vp_*t11*kp_*t3*dn_+t12*t28*dn_)

res = val(1)

end subroutine delC0_13

!Here is an automatically generated fortran subroutine:

subroutine delC1_13(res)

use constants, only: dpc, dp, vpc
use conduct_parameters
use gal_corr, only: dVeq

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      t2 = 1.0d0/omega_
      t3 = ks_*t2
      t4 = vT_**2
      val(1) = m_*(n_*(t4*t3*r_*dVeq+2*vT_*t3*r_*dvT_*Vp_)+Vp_*t4*ks_*t2 &
	*dn_*r_)

res = val(1)

end subroutine delC1_13

!Here is an automatically generated fortran subroutine:

subroutine delC2_13(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_13

!Here is an automatically generated fortran subroutine:

subroutine delC0_21(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC0_21

!Here is an automatically generated fortran subroutine:

subroutine delC1_21(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC1_21

!Here is an automatically generated fortran subroutine:

subroutine delC2_21(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_21

!Here is an automatically generated fortran subroutine:

subroutine delC0_22(res)

use constants, only: dpc, dp, vpc
use conduct_parameters
use gal_corr, only: dVeq

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      t1 = cmplx(0.D0,-1.D0,DP)
      t2 = r_*t1
      t4 = 1.0d0/omega_
      t6 = vT_**2
      t9 = cmplx(0.D0,-2.D0,DP)
      t18 = cmplx(0.D0,2.D0,DP)
      t31 = cmplx(0.D0,1.D0,DP)
      val(1) = m_*(n_*(kp_*(t6*ks_*t4*dVeq*t2+ks_*vT_*t4*Vp_*dvT_*r_*t9) &
	+dvT_*ks_*vT_*r_*t18)+Vp_*t6*kp_*ks_*t4*dn_*t2+t6*ks_*dn_*r_*t31)

res = val(1)

end subroutine delC0_22

!Here is an automatically generated fortran subroutine:

subroutine delC1_22(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC1_22

!Here is an automatically generated fortran subroutine:

subroutine delC2_22(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_22

!Here is an automatically generated fortran subroutine:

subroutine delC0_23(res)

use constants, only: dpc, dp, vpc
use conduct_parameters
use gal_corr, only: dVeq

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      t1 = cmplx(0.D0,1.D0,DP)
      t2 = r_*t1
      t4 = 1.0d0/omega_
      t5 = ks_**2
      t6 = t5*t4
      t7 = vT_**2
      t10 = cmplx(0.D0,2.D0,DP)
      val(1) = m_*(n_*(t7*t6*dVeq*t2+vT_*t5*t4*dvT_*Vp_*r_*t10)+Vp_*t7*t6 &
	*dn_*t2)

res = val(1)

end subroutine delC0_23

!Here is an automatically generated fortran subroutine:

subroutine delC1_23(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC1_23

!Here is an automatically generated fortran subroutine:

subroutine delC2_23(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_23

!Here is an automatically generated fortran subroutine:

subroutine delC0_31(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC0_31

!Here is an automatically generated fortran subroutine:

subroutine delC1_31(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC1_31

!Here is an automatically generated fortran subroutine:

subroutine delC2_31(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_31

!Here is an automatically generated fortran subroutine:

subroutine delC0_32(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC0_32

!Here is an automatically generated fortran subroutine:

subroutine delC1_32(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC1_32

!Here is an automatically generated fortran subroutine:

subroutine delC2_32(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_32

!Here is an automatically generated fortran subroutine:

subroutine delC0_33(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC0_33

!Here is an automatically generated fortran subroutine:

subroutine delC1_33(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC1_33

!Here is an automatically generated fortran subroutine:

subroutine delC2_33(res)

use constants, only: dpc, dp, vpc
use conduct_parameters

implicit complex(vpc) (s-t)
complex(dpc), intent(out) :: res
complex(vpc), dimension(1) :: val

      val(1) = 0

res = val(1)

end subroutine delC2_33
