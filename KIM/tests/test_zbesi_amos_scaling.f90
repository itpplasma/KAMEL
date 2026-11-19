program test_zbesi_amos_scaling
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    real(dp) :: pos_re(1), pos_im(1), neg_re(1), neg_im(1)
    integer :: nz, ierr

    call zbesi(1.0_dp, 0.0_dp, 0.0_dp, 2, 1, pos_re, pos_im, nz, ierr)
    if (ierr /= 0) error stop 'zbesi failed for positive real argument'

    call zbesi(-1.0_dp, 0.0_dp, 0.0_dp, 2, 1, neg_re, neg_im, nz, ierr)
    if (ierr /= 0) error stop 'zbesi failed for negative real argument'

    if (abs(pos_re(1) - exp(-1.0_dp) * 1.2660658777520082_dp) > 1.0e-12_dp) then
        error stop 'KODE=2 must return exp(-1) I0(1) for z = 1'
    end if
    if (abs(pos_re(1) - neg_re(1)) > 1.0e-12_dp) then
        error stop 'KODE=2 must scale I0 by exp(-abs(real(z)))'
    end if
    if (abs(pos_im(1)) > 1.0e-12_dp .or. abs(neg_im(1)) > 1.0e-12_dp) then
        error stop 'I0 of a real argument must be real'
    end if
end program test_zbesi_amos_scaling
