!--------------------------------------------------------------------

subroutine calc_dispersion (r, flagback, flagprint, kval, polvec)

use constants, only: im
use flre_sett, only: Nwaves
use mode_data, only: sp_ptr
use background, only: flag_back

implicit none;

integer, parameter :: dp = 8, dpc = 8

real(dp), intent(in) :: r
character(*), intent(in) :: flagback
integer, intent(in) :: flagprint
complex(dpc), dimension(Nwaves), intent(out) :: kval
complex(dpc), dimension(Nwaves,Nwaves), intent(out) :: polvec !(components, waves)

integer :: i
complex(dpc), dimension(Nwaves,Nwaves) :: cc, evc
complex(dpc), dimension(Nwaves) :: evl
complex(dpc), dimension(Nwaves) :: tmp ! temp for sorting
integer, dimension(Nwaves) :: iperm
integer :: ier

if (flagback == flag_back) then
    !call eval_diff_sys_matrix_f (sp_ptr, r, cc);
    call calc_diff_sys_matrix (r, flagback, cc)
else
    call calc_diff_sys_matrix (r, flagback, cc)
end if

call eigsys (Nwaves, cc, evl, evc)

kval = -im*evl
polvec = evc

!sorting: might be important for start_vals!
call dpsort (abs(kval), Nwaves, iperm, 1, ier)
if (ier /=0) then
    print *, 'dispersion: error: sorting failed!'
end if

!rearranging:
tmp = kval
do i=1,Nwaves
    kval(i) = tmp(iperm(i))
end do

cc=polvec
do i=1, Nwaves
    polvec(:,i) = cc(:,iperm(i))
end do

if (flagprint /= 0) then
    print *, 'r=',r
    do i=1,Nwaves
        print *, kval(i)
    end do
end if

end subroutine

!--------------------------------------------------------------------

subroutine calc_dispersion_grid (dimr, rvec, flagback) !old code - not used now

use constants, only: im
use conduct_parameters, only: kp_, ks_
use flre_sett, only: Nwaves
use mode_data, only: path2dispersion, m

implicit none;

integer, parameter :: dp = 8, dpc = 8

integer, intent(in) :: dimr
real(dp), dimension(dimr), intent(in) :: rvec
character(*), intent(in) :: flagback

real(dp), dimension(dimr,Nwaves+1) :: aknr, akrr
real(dp), dimension(dimr) :: kp_v, ks_v

complex(dp), dimension(dimr,Nwaves+1) :: knr, krr
complex(dpc), dimension(Nwaves,Nwaves) :: polvec
complex(dpc), dimension(Nwaves) :: kval

integer :: i
real(dp) :: r

do i=1,dimr
    r = rvec(i)

    call calc_dispersion (r, flagback, 0, kval, polvec)

    if (flagback /= 'f') then
        knr(i,1) = r + im*r
        knr(i,2:Nwaves+1) = kval

        aknr(i,1) = r
        aknr(i,2:Nwaves+1) = abs(kval)

        krr(i,1) = r + im*r
        krr(i,2:Nwaves+1) = sqrt(kval**2 - m**2/r**2)

        akrr(i,1) = r
        akrr(i,2:Nwaves+1) = abs(sqrt(kval**2 - m**2/r**2))

    else
        krr(i,1) = r + im*r
        krr(i,2:Nwaves+1) = kval

        akrr(i,1) = r
        akrr(i,2:Nwaves+1) = abs(kval)

        knr(i,1) = r + im*r
        knr(i,2:Nwaves+1) = sqrt(kval**2 + m**2/r**2)

        aknr(i,1) = r
        aknr(i,2:Nwaves+1) = abs(sqrt(kval**2 + m**2/r**2))
       end if

       kp_v(i) = kp_
       ks_v(i) = ks_

end do

end subroutine

!--------------------------------------------------------------------
