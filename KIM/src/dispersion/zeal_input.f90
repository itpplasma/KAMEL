!-----------------------------------------------------------------------
! File copied from ZEAL package. Adapted for KIM.
! Comments kept from the original.
!-----------------------------------------------------------------------

module Zeal_Input_Module

     use Precision_Module
     use config_m, only: VERBOSE => WKB_verbose

     implicit none

     public

!-----------------------------------------------------------------------
!  specify the geometry of the rectangular region.
!  the lower left vertex has coordinates lv(1) and lv(2). the stepsizes
!  along the x- and y-direction are h(1) and h(2), respectively.
!
!  note: these are now variables (not parameters) to allow adaptive
!  search regions centered on previous roots.
!
     real(kind=dp), dimension(2) :: LV = (/-5.0_dp, -5.0_dp/)
     real(kind=dp), dimension(2) :: H  = (/ 5.0_dp, 5.0_dp/)

!-----------------------------------------------------------------------
!  The value of M determines the maximum number of zeros
!  (counting multiplicities) that are considered within a subregion.
!  M has to be larger than the maximum of the multiplicities of the
!  zeros. A recommended value is 5.
!
     integer, parameter :: M = 5

!-----------------------------------------------------------------------
!  Specify the values of ICON and NR.
!  The parameter ICON may take the values 1, 2, 3 and 4.
!
!    1  =>  calculation of the total number of zeros, only.
!
!    2  =>  calculation of the total number of zeros and isolation
!           of subregions that contain at most M zeros.
!
!    3  =>  calculation of the total number of zeros, isolation of
!           subregions that contain at most M zeros, and computation
!           of all the zeros and their multiplicity.
!
!    4  =>  calculation of the total number of zeros, isolation of
!           subregions that contain at most M zeros, and computation
!           of NR mutually distinct zeros and their multiplicity.
!
     integer, parameter :: ICON = 3  ! Find ALL zeros (not just NR)
     integer            :: NR   = 2  ! Only used if ICON=4 (ZEAL modifies this internally)

!-----------------------------------------------------------------------
!  VERBOSE is imported from config_m (WKB_verbose) and controlled via
!  the WKB_DISPERSION namelist in KIM_config.nml.

!-----------------------------------------------------------------------
!  Specify the value of FILES.
!  If FILES is set equal to .TRUE. then ZEAL generates the files
!  "zeros.dat" and "mult.dat". They contain the computed approximations
!  for the zeros as well as their respective multiplicities. ZEAL also
!  writes the file "fzeros.dat", which contains the values that the
!  function takes at the computed approximations for the zeros.
!
     logical, parameter :: FILES = .FALSE.

!-----------------------------------------------------------------------
!  Specify the value of IFAIL.
!  This variable determines how errors are to be handled.
!  We follow the NAG convention:
!
!    1  =>  "soft silent error"
!           Control returned to calling program.
!
!   -1  =>  "soft noisy error"
!           Error message is printed.
!           Control returned to calling program.
!
!    0  =>  "hard noisy error"
!           Error message is printed and program is stopped.
!
     integer, parameter :: IFAIL = -1


contains

     subroutine set_zeal_search_region(center_re, center_im, halfwidth)
     !> Set the ZEAL search region centered on a given point.
         real(kind=dp), intent(in) :: center_re, center_im
         real(kind=dp), intent(in) :: halfwidth

         ! Set lower-left vertex and box size
         LV(1) = center_re - halfwidth
         LV(2) = center_im - halfwidth
         H(1) = 2.0_DP * halfwidth
         H(2) = 2.0_DP * halfwidth
     end subroutine set_zeal_search_region

end module Zeal_Input_Module
