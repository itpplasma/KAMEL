!<Calculates array of D special functions.

!------------------------------------------------------------------------------

subroutine calc_D_array ()

use flre_sett, only: flre_order;

implicit none;

if (flre_order == 1) then

    call calc_D_array_1 ();

else if (flre_order == 3) then

    call calc_D_array_3 ();

else if (flre_order == 5) then

    call calc_D_array_5 ();

else

    print *, 'calc_D_array: error: function for this order of expansion is not implemented:', flre_order;
    stop;

end if

end subroutine

!------------------------------------------------------------------------------
