
!Here is an automatically generated fortran subroutine:
!see also f0moments_J0.f90

subroutine dens_mom(res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t7 = 1/r_
      t13 = hz_**2
      t16 = omc_**2
      val(1) = 1/omc_*n_*t7*(r_*omc_+Vp_*dpsi_*r_+Vp_*hz_*ht_)+q_/m_*n_ &
    *(dPhi0_*t7/t16/omc_*(-domc_*r_+omc_*t13)+1/t16*ddPhi0_)

res = val(1)

end subroutine dens_mom

!Here is an automatically generated fortran subroutine:

subroutine Vth_mom(res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t1 = hz_**2
      t2 = t1*hz_
      t3 = r_*t2
      t4 = omc_**2
      t6 = r_*hz_
      t7 = t4*t6
      t8 = r_**2
      t9 = t8*hz_
      t13 = ht_*omc_
      t19 = Vp_*omc_
      t20 = dpsi_**2
      t28 = dpsi_*t6*t19
      t34 = ht_*t1*t19
      t35 = hz_*domc_
      t41 = r_*Vp_
      t45 = t4*t3-t7+omc_*domc_*t9-dVp_*t13*r_*t1-omc_*dVp_*dpsi_*t9+t8 &
    *t20*ht_*t19-3*dpsi_*t3*t19+2*t28-t8*ddpsi_*hz_*t19+t34+2*dpsi_*t8*Vp_*t35 &
    +2*ht_*t41*domc_*t1
      t47 = 1/t8/r_
      t50 = 1/t4/omc_
      t51 = vT_**2
      t54 = dvT_*hz_
      t55 = r_*omc_
      t58 = Vp_*hz_
      t60 = t55+Vp_*dpsi_*r_+ht_*t58
      t62 = 1/t8
      t63 = 1/t4
      t71 = dpsi_*r_
      t81 = dn_*hz_
      t86 = omc_*t8
      t87 = ht_*domc_
      t101 = domc_**2
      t107 = t4**2
      t118 = 1/t107
      t122 = 1/r_
      t123 = t50*t122
      t129 = -domc_*r_+omc_*t1
      t131 = t118*t62
      t132 = dPhi0_*t131
      t135 = t50*t122*ddPhi0_
      t161 = dPhi0_**2
      t169 = q_**2
      t170 = m_**2
      val(1) = n_*(-t51*t50*t47*t45+2*vT_*t63*t62*t60*t54+1/omc_*t62*(ht_ &
    *r_*omc_+t71*Vp_*ht_-t2*Vp_+t58)*Vp_)+t51*t63*t60*t62*t81+1/m_*q_*(n_*(t51 &
    *(-dPhi0_/t107/omc_*t47*(-dpsi_*t87*t86+3*t2*domc_*t55+hz_*ddomc_*t86+3 &
    *ht_*dpsi_*r_*t4*t1+t4*t2-4*t8*t101*hz_)-ddPhi0_*t118*t62*(t71*t13-t2*omc_ &
    +4*r_*t35)+t123*hz_*dddPhi0_)+2*vT_*(t132*t129*t54+t135*t54)+dPhi0_*t50 &
    *t62*(t28+2*t34-t41*t87+t7)+Vp_*t63*t122*ht_*ddPhi0_)+t51*(t132*t129*t81 &
    +t135*t81))+1/t170*t169*n_*(t161*t131*t129*hz_+dPhi0_*t123*hz_*ddPhi0_)

res = val(1)

end subroutine Vth_mom

!Here is an automatically generated fortran subroutine:

subroutine Vz_mom(res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t1 = omc_**2
      t2 = t1*ht_
      t3 = hz_**2
      t5 = omc_*ht_
      t6 = domc_*r_
      t8 = hz_*omc_
      t10 = t3*hz_
      t16 = Vp_*omc_
      t17 = dpsi_**2
      t21 = ht_*t3
      t28 = domc_*Vp_
      t30 = ht_*dpsi_*r_
      t37 = -t2+t3*t2+t6*t5-dVp_*t8+dVp_*t10*omc_-dVp_*dpsi_*r_*t5-r_*t17 &
    *hz_*t16-3*dpsi_*t21*t16-ht_*ddpsi_*r_*t16+2*t30*t28+2*hz_*t28-2*t10*t28
      t38 = 1/r_
      t41 = 1/t1/omc_
      t42 = vT_**2
      t45 = r_*ht_
      t48 = dpsi_*r_
      t51 = Vp_*hz_
      t52 = omc_*t45+t48*Vp_*ht_-t10*Vp_+t51
      t54 = 1/t1
      t55 = t54*t38
      t73 = omc_*domc_
      t74 = dpsi_*hz_
      t87 = domc_**2
      t93 = t1**2
      t105 = 1/t93
      t114 = -t6+omc_*t3
      t116 = t105*t38
      t117 = dPhi0_*t116
      t133 = dPhi0_*t41
      t140 = dn_*ht_
      t152 = dPhi0_**2
      t159 = q_**2
      t160 = m_**2
      val(1) = n_*(t42*t41*t38*t37-2*vT_*t55*t52*dvT_+1/omc_*t38*(r_*omc_ &
    +Vp_*dpsi_*r_+ht_*t51)*t51)-t42*t55*t52*dn_+1/m_*q_*(n_*(t42*(-dPhi0_/t93 &
    /omc_*t38*(-r_*t74*t73-3*t21*t73-domc_*t5-ddomc_*r_*t5+3*t1*dpsi_*t10-2 &
    *t1*t74+4*r_*t87*ht_)-ddPhi0_*t105*t38*(t48*t8+t5+t3*t5-4*ht_*domc_*r_) &
    -t41*dddPhi0_*ht_)+2*vT_*(-t117*t114*dvT_*ht_-t41*ht_*dvT_*ddPhi0_)-t133 &
    *t38*(t30*t16+hz_*t16-2*t10*t16+hz_*r_*t28+t1*t45)+Vp_*t54*hz_*ddPhi0_) &
    +t42*(-t117*t114*t140-t41*ddPhi0_*t140))+1/t160*t159*n_*(-t152*t116*t114 &
    *ht_-t133*ht_*ddPhi0_)

res = val(1)

end subroutine Vz_mom

!Here is an automatically generated fortran subroutine:

subroutine Vs_mom(res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t1 = Vp_*omc_
      t2 = hz_**2
      t6 = r_**2
      t12 = domc_*Vp_
      t20 = omc_*t6
      t24 = r_*omc_
      t28 = omc_**2
      t29 = t28*r_
      t32 = 1/t6
      t35 = 1/t28/omc_
      t36 = vT_**2
      t43 = t24+Vp_*dpsi_*r_+Vp_*hz_*ht_
      t45 = 1/r_
      t46 = 1/t28
      t62 = domc_**2
      t65 = t2**2
      t74 = t28**2
      t79 = domc_*r_
      t83 = 1/t74
      t90 = -t79+omc_*t2
      t93 = dPhi0_*t83*t45
      t115 = dPhi0_**2
      t122 = q_**2
      t123 = m_**2
      val(1) = n_*(-t36*t35*t32*(-r_*dpsi_*t2*t1-t6*ddpsi_*t1+ht_*t2*hz_ &
    *t1+2*dpsi_*t6*t12+2*r_*hz_*ht_*t12-dVp_*dpsi_*t20+domc_*t20-hz_*ht_*dVp_ &
    *t24+t2*t29-t29)+2*vT_*t46*t45*t43*dvT_)+t36*t46*t43*t45*dn_+1/m_*q_*(n_ &
    *(t36*(-dPhi0_/t74/omc_*t32*(ddomc_*t20+domc_*t24+2*t2*domc_*t24-4*t6*t62 &
    +t28*t65+2*ht_*dpsi_*r_*t28*hz_)+ddPhi0_*t83*t45*(omc_-4*t79)+t35*dddPhi0_) &
    +2*vT_*(t93*t90*dvT_+t35*dvT_*ddPhi0_)+dPhi0_*t46*t45*t43)+t36*(t93*t90 &
    *dn_+t35*ddPhi0_*dn_))+1/t123*t122*n_*(t115*t83*t45*t90+dPhi0_*t35*ddPhi0_)

res = val(1)

end subroutine Vs_mom

!Here is an automatically generated fortran subroutine:

subroutine Vp_mom(res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t1 = dpsi_**2
      t2 = r_**2
      t4 = ht_*hz_
      t5 = dpsi_*r_
      t8 = hz_**2
      t9 = t8**2
      t11 = vT_**2
      t13 = 1/t2
      t15 = omc_**2
      t16 = 1/t15
      t26 = 1/r_
      t32 = t5+t4
      t35 = -domc_*r_+omc_*t8
      t37 = t15**2
      t44 = 1/t15/omc_
      val(1) = n_*(-t16*t13*Vp_*t11*(t2*t1+2*t5*t4+t8-t9)+1/omc_*t26*Vp_ &
    *(r_*omc_+Vp_*dpsi_*r_+Vp_*hz_*ht_))+q_/m_*n_*(t11*(-dPhi0_/t37*t13*t35 &
    *t32-t44*ddPhi0_*t32*t26)+dPhi0_*t44*t26*t35*Vp_+t16*Vp_*ddPhi0_)

res = val(1)

end subroutine Vp_mom

!Here is an automatically generated fortran subroutine:

subroutine Etot_mom(res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t1 = r_*omc_
      t4 = Vp_*dpsi_*r_
      t7 = Vp_*hz_*ht_
      t10 = 1/r_
      t12 = 1/omc_
      t13 = vT_**2
      t16 = Vp_**2
      t25 = hz_**2
      t27 = -domc_*r_+omc_*t25
      t28 = omc_**2
      t30 = 1/t28/omc_
      t34 = 1/t28
      val(1) = m_*n_*(t13*t12*t10*(3*t1+5*t4+5*t7)+t12*t10*(t1+t4+t7)*t16) &
    /2+q_*n_*(3.D0/2.D0*t13*(dPhi0_*t10*t30*t27+t34*ddPhi0_)+dPhi0_*t30*t10 &
    *t27*t16/2+t34*t16*ddPhi0_/2)

res = val(1)

end subroutine Etot_mom

!Here is an automatically generated fortran subroutine:

subroutine Eterm_mom(res)

use constants, only: dpc, dp
use conduct_parameters

implicit real(dp) (s-t)
real(dp), intent(out) :: res
real(dp), dimension(1) :: val

      t7 = 1/r_
      t11 = vT_**2
      t17 = hz_**2
      t20 = omc_**2
      val(1) = 3.D0/2.D0*m_*n_*t11/omc_*t7*(r_*omc_+Vp_*dpsi_*r_+Vp_*hz_ &
    *ht_)+3.D0/2.D0*n_*q_*t11*(dPhi0_*t7/t20/omc_*(-domc_*r_+omc_*t17)+1/t20 &
    *ddPhi0_)

res = val(1)

end subroutine Eterm_mom
