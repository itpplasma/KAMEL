subroutine kilca_matrix (r, amat)

use constants, only: dp, dpc;
use eqs_sett, only: Nwaves;

implicit none;

real(dp), intent(in) :: r;
complex(dpc), dimension(Nwaves, Nwaves), intent(out) :: amat;

integer, save :: flg = 0;

if (flg == 0) then
    call init_kilca_data ();
    print *, 'Nwaves=', Nwaves
    flg = 1;
end if

call calc_diff_sys_matrix (r, 'f', amat);

end subroutine
