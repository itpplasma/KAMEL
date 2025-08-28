MODULE Function_Input_Module

     USE Precision_Module

     IMPLICIT NONE

!-----------------------------------------------------------------------
!  If ICON = 3 or 4 (as specified in Zeal_Input_Module), then specify 
!  the values of NEWTONZ and NEWTONF.
!  These variables are used as follows. The modified Newton's method,
!  which takes into account the multiplicity of a zero and converges
!  quadratically, is used to refine the calculated approximations for
!  the zeros. The iteration stops if
!    - the relative distance between two successive approximations is
!      at most NEWTONZ
!  or
!    - the absolute value of the function at the last approximation is
!      at most NEWTONF
!  or if a maximum number of iterations is exceeded.
!
     double precision :: NEWTONZ = 5.0E-08
     double precision :: NEWTONF = 1.0E-14

     CONTAINS


     subroutine FDF(Z, F, DF)

          use config_m, only: type_of_run

          implicit none
          double complex, intent(in) :: Z
          double complex, intent(out) :: F, DF

          select case (trim(type_of_run))
               case ("WKB_dispersion")
                    call FDF_WKB(Z, F, DF)
               case default 
                    stop "Invalid run type"
          end select

     end subroutine


!     SUBROUTINE FDF(Z,F,DF)
!!-----------------------------------------------------------------------
!!**PURPOSE
!!  Define the function F and its derivative DF.
!!-----------------------------------------------------------------------
     !COMPLEX(KIND=DP), INTENT(IN)   :: Z
     !COMPLEX(KIND=DP), INTENT(OUT)  :: F, DF

!!     F = EXP(3.0_DP*Z) + TWO*Z*COS(Z) - ONE
!!     DF = 3.0_DP*EXP(3.0_DP*Z) + TWO*COS(Z) - TWO*Z*SIN(Z)

!!
!!  Other examples
!!
    !F = (Z**2)*(Z-ONE)*(Z-TWO)*(Z-3.0_DP)*(Z-4.0_DP) + Z*SIN(Z)
    !DF = 6.0_DP*Z**5 - 50.0_DP*Z**4 + 140.0_DP*Z**3 - 150.0_DP*Z**2 &
         !+ 48.0_DP*Z + Z*COS(Z) + SIN(Z)

!!    F = (Z**2)*((Z-TWO)**2)*(COS(Z)*EXP(TWO*Z)+Z**3-ONE-SIN(Z))
!!    DF = TWO*Z*((Z-TWO)**2)*(COS(Z)*EXP(TWO*Z)+Z**3-ONE-SIN(Z))    +  &
!!         (Z**2)*TWO*(Z-TWO)*(COS(Z)*EXP(TWO*Z)+Z**3-ONE-SIN(Z))  +  &
!!         (Z**2)*((Z-TWO)**2)*(-SIN(Z)*EXP(TWO*Z)                 +  &
!!         TWO*COS(Z)*EXP(TWO*Z)+3.0_DP*Z**2-COS(Z))

     !END SUBROUTINE FDF


     FUNCTION VALREG(LV,H)
!-----------------------------------------------------------------------
!**PURPOSE
!  Given a rectangular region specified by LV and H (cf. the module
!  Zeal_Input_Module), decide whether the function is analytic inside 
!  this region or not.
!-----------------------------------------------------------------------
     LOGICAL VALREG
     REAL(KIND=DP), INTENT(IN) :: LV(2), H(2)

     VALREG = .TRUE.

!  The following statement can be used for functions that have a
!  branch cut along the non-positive real axis.
!
!    VALREG = .NOT. ( LV(2)*(LV(2)+H(2)) <= ZERO .AND. LV(1) <= ZERO )

     END FUNCTION VALREG


END MODULE Function_Input_Module


