!------------------------------------------------------------------------------

subroutine calc_W1_array ()

!only to be called after vT_ is updated in cond_parameters module for given r_!

use constants, only: sqrt2p;
use conduct_parameters, only: vT_;
use conduct_arrays, only: W1;

implicit none;

W1(0) = sqrt2p * vT_;
W1(1) = 0.0d0;
W1(2) = W1(0) * vT_ * vT_;
W1(3) = 0.0d0;

end subroutine

!------------------------------------------------------------------------------
