!<The module contains all parameters needed for evaluation of conductivity matrices

module conduct_parameters

use constants, only: dpc, dp;

implicit none;

!array contains quants evaluated by splines: target for pointers.
real(dp), dimension(0:30), target :: back_data

!background related variables:
!background quants which do not depend on a sort of particles
real(dp), pointer :: ht_, dht_, ddht_, hz_, dhz_, ddhz_
real(dp), pointer :: B, dB, ddB

real(dp) :: r_, dpsi_, ddpsi_

!background quants which depend on a sort of particles
real(dp), pointer :: nu_, dnu_, ddnu_

!parameters of dstr. function:
real(dp), pointer :: n_, dn_, ddn_, vT_, dvT_, ddvT_, Vp_, dVp_, ddVp_

!electric field;
real(dp), pointer :: dPhi0_, ddPhi0_, dddPhi0_

real(dp) :: h_t, h_z

!background quants which depend on a sort of particles
real(dp) :: q_, m_, omc_, domc_, ddomc_

!quants which are related to wave and particles
real(dp) :: kp_, dkp_, ddkp_, ks_, dks_, ddks_
complex(dpc) :: omega_

complex(dpc) :: oms_, doms_, ddoms_
complex(dpc) :: oms_c_, doms_c_, ddoms_c_

!for plasma function: allocated in set_back_aliases_in_conductivity_parameters()
complex(dpc), dimension(-100:100) :: zv, Wv, Wmv
complex(dpc), dimension(1:2,-100:100) :: dzv, dWv, dWmv

end module

!------------------------------------------------------------------------------

subroutine set_back_aliases_in_conductivity_parameters ()

use conduct_parameters;

implicit none;

integer :: ind;

ind = 0;

B => back_data(ind); ind = ind+1;

dB => back_data(ind); ind = ind+1;

ddB => back_data(ind); ind = ind+1;

ht_ => back_data(ind); ind = ind+1;

dht_ => back_data(ind); ind = ind+1;

ddht_ => back_data(ind); ind = ind+1;

hz_ => back_data(ind); ind = ind+1;

dhz_ => back_data(ind); ind = ind+1;

ddhz_ => back_data(ind); ind = ind+1;

dPhi0_ => back_data(ind); ind = ind+1;

ddPhi0_ => back_data(ind); ind = ind+1;

dddPhi0_ => back_data(ind); ind = ind+1;

n_ => back_data(ind); ind = ind+1;

dn_ => back_data(ind); ind = ind+1;

ddn_ => back_data(ind); ind = ind+1;

Vp_ => back_data(ind); ind = ind+1;

dVp_ => back_data(ind); ind = ind+1;

ddVp_ => back_data(ind); ind = ind+1;

vT_ => back_data(ind); ind = ind+1;

dvT_ => back_data(ind); ind = ind+1;

ddvT_ => back_data(ind); ind = ind+1;

nu_ => back_data(ind); ind = ind+1;

dnu_ => back_data(ind); ind = ind+1;

ddnu_ => back_data(ind); ind = ind+1;

!Nmax must match to maximum cyclotron harmonic used in conductivity sources <=100
!allocate (zv(-100:100), Wv(-100:100), Wmv(-100:100));
!allocate (dzv(1:2,-100:100), dWv(1:2,-100:100), dWmv(1:2,-100:100));

end subroutine

!------------------------------------------------------------------------------

subroutine eval_and_set_conduct_parameters_spec_independent(r, flag_back)

use constants, only: dp;

implicit none;

real(dp), intent(in) :: r
!character(1), intent(in) :: flag_back
character(*), intent(in) :: flag_back

call eval_and_set_background_parameters_spec_independent (r, flag_back);
call eval_and_set_wave_parameters (r, flag_back);

end subroutine

!------------------------------------------------------------------------------

subroutine eval_and_set_conduct_parameters_spec_dependent(r, spec, flag_back)

use constants, only: dp;

implicit none;

real(dp), intent(in) :: r
integer, intent(in) :: spec
!character(1), intent(in) :: flag_back
character(*), intent(in) :: flag_back

call eval_and_set_background_parameters_spec_dependent (r, spec, flag_back);
call eval_and_set_f0_parameters_nu_and_derivs (r, spec, flag_back);
call eval_z_and_oms ();
call eval_w_array ();

end subroutine

!------------------------------------------------------------------------------

subroutine eval_and_set_all_conduct_parameters(r, spec, flag_back)

use constants, only: dp;

implicit none;

real(dp), intent(in) :: r
integer, intent(in) :: spec
!character(1), intent(in) :: flag_back
character(*), intent(in) :: flag_back

call eval_and_set_conduct_parameters_spec_independent(r, flag_back);
call eval_and_set_conduct_parameters_spec_dependent(r, spec, flag_back);

end subroutine

!------------------------------------------------------------------------------

subroutine eval_and_set_background_parameters_spec_independent (r, flag_back)

!evaluates and sets background quants (B, hth, hz, dPhi0 + derivs) which do not depend on a sort of particles
use constants, only: dp;
use core, only: bp_ptr;
use conduct_parameters;
use background, only: B0, huge_factor;

implicit none;

external :: eval_background_spec_independent

real(dp), intent(in) :: r
character(*), intent(in) :: flag_back

r_ = r

call eval_background_spec_independent (r, bp_ptr, B); !fills a part of the back_data array

h_t = ht_;
h_z = hz_;

dpsi_ = dht_/hz_;
ddpsi_ = (ddht_-dhz_*dpsi_)/hz_;

if (flag_back == 'w') then
    !keep rotational transform and set only derivs to zero
    dhz_ = 0.0d0;
    ddhz_ = 0.0d0;
    dht_ = 0.0d0;
    ddht_ = 0.0d0;
    dpsi_ = 0.0d0;
    ddpsi_ = 0.0d0;
    dB = 0.0d0;
    ddB = 0.0d0;
    ddPhi0_ = 0.0d0;
    dddPhi0_ = 0.0d0;
    r_ = r_*huge_factor
else if (flag_back == 'h') then
    !set derivs and hth to zero, hz=1
    ht_ = 0.0d0;
    hz_ = 1.0d0;
    dht_ = 0.0d0;
    dhz_ = 0.0d0;
    ddht_ = 0.0d0;
    ddhz_ = 0.0d0;
    dpsi_ = 0.0d0;
    ddpsi_ = 0.0d0;
    B = B0; !use B at the center!
    dB = 0.0d0;
    ddB = 0.0d0;
    dPhi0_ = 0.0d0;
    ddPhi0_ = 0.0d0;
    dddPhi0_ = 0.0d0;
    r_ = r_*huge_factor
else if (flag_back == 's') then
    print *, 'warning: this flag_back is not implemented!..'
end if

end subroutine

!------------------------------------------------------------------------------

subroutine get_magnetic_field_parameters (hvals)

use conduct_parameters, only: ht_, hz_;
use constants, only: dp;

implicit none;

real(dp), dimension(3), intent(out) :: hvals

hvals(1) = 0.0d0;
hvals(2) = ht_;
hvals(3) = hz_;

end subroutine

!------------------------------------------------------------------------------

subroutine eval_and_set_background_parameters_spec_dependent (r, spec, flag_back)

!evaluates and sets back quants which depend on a sort of particles

use constants, only: dp, c;
use conduct_parameters;
use background, only: mass, charge;

implicit none;

real(dp), intent(in) :: r
integer, intent(in) :: spec
!character(1), intent(in) :: flag_back
character(*), intent(in) :: flag_back

m_ = mass(spec)
q_ = charge(spec)

omc_ = q_*B/m_/c
domc_ = q_*dB/m_/c
ddomc_ = q_*ddB/m_/c

end subroutine

!------------------------------------------------------------------------------

subroutine eval_and_set_wave_parameters (r, flag_back)

!evaluates and sets quants which are related to wave
!must be called after eval_back_spec_independent only!!!

use conduct_parameters;
use constants, only: dp;
use background, only: rtor;
use mode_data, only: m, n, omov;

implicit none;

real(dp), intent(in) :: r
!character(1), intent(in) :: flag_back
character(*), intent(in) :: flag_back

real(dp) :: kth, dkth, ddkth, kz, dkz, ddkz

kth = m/r
dkth = -m/r**2
ddkth = 2.d0*m/r**3

kz = n/rtor
dkz = 0.d0
ddkz = 0.d0

omega_ = omov;

kp_ = ht_*kth + hz_*kz
ks_ = hz_*kth - ht_*kz

dkp_ = dht_*kth + dhz_*kz + ht_*dkth + hz_*dkz
dks_ = dhz_*kth - dht_*kz + hz_*dkth - ht_*dkz

ddkp_ = ddht_*kth + dht_*dkth + ddhz_*kz + dhz_*dkz + &
         dht_*dkth + ht_*ddkth + dhz_*dkz + hz_*ddkz

ddks_ = ddhz_*kth + dhz_*dkth - ddht_*kz - dht_*dkz + &
         dhz_*dkth + hz_*ddkth - dht_*dkz - ht_*ddkz

if (flag_back == 'w') then
    !keep rotational transform and set only derivs to zero
    !ks_ = 0.0d0;
    dks_ = 0.0d0;
    ddks_ = 0.0d0;
    dkp_ = 0.0d0;
    ddkp_ = 0.0d0;
else if (flag_back == 'h') then
    !ks_ = 1.0d-6;
    !ks_ = 1.0d-2;
    !ks_ = h_z*kth - h_t*kz;
    dks_ = 0.0d0;
    ddks_ = 0.0d0;
    kp_ = 1.0d-3; !kz; !1.0d-4;
    !kp_ = h_t*kth + h_z*kz;
    dkp_ = 0.0d0;
    ddkp_ = 0.0d0;
else if (flag_back == 's') then
    print *, 'warning: the flag_back is not implemented!..'
end if

end subroutine

!------------------------------------------------------------------------------

subroutine get_wave_parameters (kvals)

use conduct_parameters, only: ks_, kp_;
use constants, only: dp;

implicit none;

real(dp), dimension(3), intent(out) :: kvals

kvals(1) = 0.0d0;
kvals(2) = ks_;
kvals(3) = kp_;

end subroutine

!------------------------------------------------------------------------------

subroutine eval_and_set_f0_parameters_nu_and_derivs (r, spec, flag_back)

!evaluates and sets (n_p, Vp_p, Vt_p, nu + derivs) which depend on a sort of particles

use core, only: bp_ptr;
use conduct_parameters;
use constants, only: dp;

implicit none;

external :: eval_f0_parameters_nu_and_derivs

real(dp), intent(in) :: r
integer, intent(in) :: spec
!character(1), intent(in) :: flag_back
character(*), intent(in) :: flag_back

call eval_f0_parameters_nu_and_derivs (r, spec, bp_ptr, n_);

if (flag_back == 'w') then !keep rotational transform and set only derivs to zero
    dn_ = 0.0d0;
    ddn_ = 0.0d0;
    dVp_ = 0.0d0;
    ddVp_ = 0.0d0;
    dvT_ = 0.0d0;
    ddvT_ = 0.0d0;
    dnu_ = 0.0d0;
    ddnu_ = 0.0d0;
else if (flag_back == 'h') then
    dn_ = 0.0d0; !comment if you need n gradient
    ddn_ = 0.0d0;
    Vp_ = 0.0d0;
    dVp_ = 0.0d0;
    ddVp_ = 0.0d0;
    dvT_ = 0.0d0;
    ddvT_ = 0.0d0;
    dnu_ = 0.0d0;
    ddnu_ = 0.0d0;
else if (flag_back == 's') then
    print *, 'warning: the flag_back is not implemented!..'
end if

end subroutine

!------------------------------------------------------------------------------

subroutine set_n_vp_dphi0 (n, Vp, dPhi0)

use conduct_parameters;
use constants, only: dp;

implicit none;

real (dp), intent(in) :: n, Vp, dPhi0

n_ = n
Vp_ = Vp
dPhi0_ = dPhi0

end subroutine

!------------------------------------------------------------------------------

subroutine eval_z_and_oms ()

use conduct_parameters;
use constants, only: dpc, dp, im;
use flre_sett, only: Nmax;

implicit complex(dpc) (t)

integer :: k
complex(dpc) :: dt_num, ddt_num, dt_den, ddt_den

oms_ = omega_-kp_*Vp_-q_*dPhi0_*ks_/m_/omc_
oms_c_ = oms_ + IM*nu_

t7 = omc_**2
doms_ = -dkp_*Vp_-kp_*dVp_+1.0d0/m_*q_*((-ddPhi0_*ks_-dPhi0_*dks_)/omc_+domc_/t7*dPhi0_*ks_)
doms_c_ = doms_ + IM*dnu_

t14 = dPhi0_*ks_
t17 = omc_**2
t22 = domc_**2
ddoms_ = -ddkp_*Vp_-kp_*ddVp_-2.d0*dkp_*dVp_+ &
    1.0d0/m_*q_*((-dddPhi0_*ks_-dPhi0_*ddks_-2.0d0*ddPhi0_*dks_)/omc_ + &
    1.0d0/t17*(2.0d0*dPhi0_*dks_*domc_+2.0d0*ks_*ddPhi0_*domc_+ddomc_*t14)-2.0d0*t22/t17/omc_*t14)
ddoms_c_ = ddoms_ + IM*ddnu_

!z-array:
t1 = sqrt (2.d0)
t_den = kp_*vT_
dt_den = dkp_*vT_ + kp_*dvT_
ddt_den = ddkp_*vT_ + 2.0d0*dkp_*dvT_ + kp_*ddvT_

do k = -Nmax,Nmax
    t_num = oms_c_ - k*omc_
    zv(k) = t_num/t_den/t1

    dt_num = doms_c_ - k*domc_
    dzv(1,k) = (dt_num*t_den - t_num*dt_den)/t_den/t_den/t1

    ddt_num = ddoms_c_ - k*ddomc_
    dzv(2,k) = (-2.0d0*dt_den*(dt_num*t_den - t_num*dt_den) + &
           t_den*(ddt_num*t_den-t_num*ddt_den))/t_den/t_den/t_den/t1
end do

end subroutine

!------------------------------------------------------------------------------

subroutine eval_w_array ()

use constants, only: dpc, dp, pi, im
use flre_sett, only: Nmax
use conduct_parameters

implicit none;

external :: wm_asympt, wfunc;

integer :: k
complex(dpc) :: isqpi, z2, t1

!I do not see high derivatives in the tensor components:
!if appears it must be segmentaion fault, because arrays are 1:2
!if(flre_order > 2) then
!   print *, 'W_array: warning: W deriv > 2 is not implemented!..'
!end if

isqpi = im/sqrt(pi);

do k = -Nmax,Nmax
    call wfunc (kp_, zv(k), Wv(k));
    dWv(1,k) = cmplx(2.0d0,0.0d0)*(isqpi - zv(k)*Wv(k));
    dWv(2,k) = cmplx(-2.0d0,0.0d0)*(Wv(k) + zv(k)*dWv(1,k));

    if (abs(zv(k)) < 10.0d0) then
        t1 = zv(k)**2;
        Wmv(k) = t1*(t1*zv(k)*Wv(k)-im/sqrt(pi)*(t1+0.5d0));
    else
        call wm_asympt (zv(k), Wmv(k));
    end if;

    !Also asymptotics can be used: not really needed in the code!
    z2 = zv(k)**2;
    dWmv(1,k) = zv(k)*(z2*zv(k)*(5.0d0*Wv(k)+zv(k)*dWv(1,k)) - &
        isqpi*(4.0d0*z2+1.0d0));

    dWmv(2,k) = z2*zv(k)*(20.0d0*Wv(k)+zv(k)* &
    (10.0d0*dWv(1,k)+zv(k)*dWv(2,k))) - isqpi*(12.0d0*z2+1.0d0);
end do

!call print_conduct_params ();

end subroutine

!------------------------------------------------------------------------------

subroutine eval_w_array_new ()
!more correct for complex omega then wm_asympt!
!warning Wv is undefined!
!definition is different with one in the old code!
use constants, only: dpc, dp, pi, im
use flre_sett, only: Nmax
use conduct_parameters, only: kp_, zv, Wmv

implicit none;

external :: wmfunc;

integer :: k

do k = -Nmax,Nmax
        call wmfunc (kp_, zv(k), Wmv(k));
end do

!call print_conduct_params ();

end subroutine

!------------------------------------------------------------------------------

subroutine print_conduct_params ()

use conduct_parameters;
use flre_sett, only: Nmax

implicit none;

integer :: k

print *
print *, 'conductivity parameters and derivatives:'
print *, 'r_:', r_
print *, 'ht_:', ht_, dht_, ddht_
print *, 'hz_:', hz_, dhz_, ddhz_
print *, 'dpsi_:', dpsi_, ddpsi_
print *, 'B:', B, dB, ddB

print *, 'q_:', q_
print *, 'm_:', m_
print *, 'omc_:', omc_, domc_, ddomc_
print *, 'nu_:', nu_, dnu_, ddnu_

print *, 'omega_:', omega_
print *, 'kp_:', kp_, dkp_, ddkp_
print *, 'ks_:', ks_, dks_, ddks_
print *, 'oms_:',oms_, doms_, ddoms_
print *, 'oms_c_:', oms_c_, doms_c_, ddoms_c_

!plasma function and z:
do k=-Nmax,Nmax
    print *,'zv(',k,'):', zv(k),dzv(1,k),dzv(2,k)
    print *,'Wv(',k,'):', Wv(k),dWv(1,k),dWv(2,k)
    print *,'Wmv(',k,'):', Wmv(k),dWmv(1,k),dWmv(2,k)
end do

!parameters of dstr. function:
print *,'n_:', n_, dn_, ddn_
print *,'Vp_:', Vp_, dVp_, ddVp_
print *,'vT_:', vT_, dvT_, ddvT_

!electric field;
print *,'dPhi0_:', dPhi0_, ddPhi0_, dddPhi0_

call flush(101)

end subroutine

!------------------------------------------------------------------------------
