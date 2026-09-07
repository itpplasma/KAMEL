!<The module is designed to store some FLRE zone properties used in Fortran part of the code.

module flre_sett

    use constants, only: dp, dpc

    implicit none;
!Conductivity settings:
    integer :: flre_order; !order of FLR expansion
    integer :: Nmax; !highest cyclotron harmonic
    integer :: gal_corr; !flag if use correction term in conductivity
    integer :: Nbmax; !max number of terms in Bessel expansions
    integer :: rsp; !use cylindrical components if rsp = 0 and rsp components otherwise

!Equations settings:
    integer :: hom_sys; !flag if system of equations should be used in homogenious limit
    integer :: Nwaves; !number of waves
    integer :: Nfs; !number of fundamental solutions
    integer :: Nphys; !number of physical modes

    integer :: flag_debug; !flag for debugging mode

    integer, dimension(0:1) :: collmod; ! collisions model flags for ions and electrons

end module
