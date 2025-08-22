!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ctensor (spec, tip, r, backflag, ct)

use constants, only: dp, dpc, pp;
use flre_sett, only: flre_order;
use background, only: flag_back;
use mode_data, only: sd_ptr, bp_ptr, cp_ptr, wd_ptr;

implicit none;

integer, intent(in) :: spec, tip;
real(dp), intent(in) :: r;
character(*), intent(in) :: backflag;
complex(dpc), dimension(1:3, 1:3, 0:2*flre_order, 0:1), intent(out) :: ct;

integer :: Dmin = 0, Dmax = 1;

integer(pp) :: cp_point;

if (backflag == flag_back) then
    call eval_c_matrices_f (cp_ptr, spec, tip, Dmin, Dmax, r, ct);
    return;
end if

!print *, 'warning: ctensor: exact sub is called:', ' r =', r, ' flagback = ', backflag

call calc_and_spline_conductivity_for_point (sd_ptr, bp_ptr, wd_ptr, backflag, r, cp_point);

call eval_c_matrices_f (cp_point, spec, tip, Dmin, Dmax, r, ct);

call delete_conductivity_profiles_f (cp_point);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kmatrices (spec, tip, Dmin, Dmax, r, backflag, K)

use constants, only: dp, dpc, pp;
use flre_sett, only: flre_order;
use background, only: flag_back;
use mode_data, only: sd_ptr, bp_ptr, cp_ptr, wd_ptr;

implicit none;

integer, intent(in) :: spec, tip, Dmin, Dmax;
real(dp), intent(in) :: r;
character(*), intent(in) :: backflag;
complex(dpc), dimension(1:3, 1:3, 0:flre_order, 0:flre_order, Dmin:Dmax), intent(out) :: K; !j,i,q,p,order

integer(pp) :: cp_point;

if (backflag == flag_back) then
    call eval_k_matrices_f (cp_ptr, spec, tip, Dmin, Dmax, r, K);
    return;
end if

!print *, 'warning: kmatrices: exact sub is called:', ' r =', r, ' flagback = ', backflag

call calc_and_spline_conductivity_for_point (sd_ptr, bp_ptr, wd_ptr, backflag, r, cp_point);

call eval_k_matrices_f (cp_point, spec, tip, Dmin, Dmax, r, K);

call delete_conductivity_profiles_f (cp_point);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
