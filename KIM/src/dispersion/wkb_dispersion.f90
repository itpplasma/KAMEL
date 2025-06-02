module rt_WKB_dispersion

    use kim_base, only: kim_t

    implicit none

    type, extends(kim_t) :: WKB_dispersion_t
        contains
            procedure :: init => init_wkb_dispersion
            procedure :: run => run_wkb_dispersion
    end type WKB_dispersion_t

    contains

    subroutine init_wkb_dispersion(this)

        implicit none

        class(WKB_dispersion_t), intent(inout) :: this

        this%run_type = "WKB_dispersion"
        print *, " ____  __.___   _____     __      __ ____  __.__________ "
        print *, "|    |/ _|   | /     \   /  \    /  \    |/ _|\______   \\"
        print *, "|      < |   |/  \ /  \  \   \/\/   /      <   |    |  _/"
        print *, "|    |  \|   /    Y    \  \        /|    |  \  |    |   \\"
        print *, "|____|__ \___\____|__  /   \__/\  / |____|__ \ |______  /"
        print *, "        \/           \/         \/          \/        \/ "

        call generate_grids

    end subroutine

    subroutine run_wkb_dispersion(this)

        use Precision_Module
        use Zeal_Module

        implicit none

        integer :: totalnumber, distinctnumber, refinednumber, j
        integer, dimension(:), pointer :: multiplicities
        logical, dimension(:), pointer :: refinement_ok
        double complex, dimension(:), pointer :: zeros, fzeros
        logical :: output

        class(WKB_dispersion_t), intent(inout) :: this

        print *, "Running WKB_dispersion"

        call zeal(totalnumber, distinctnumber, zeros, fzeros, &
                    multiplicities, refinednumber, refinement_ok)
        
        output = .true.
        if (output) then
            print *
            print *, 'Total number of zeros = ', totalnumber
            print *, 'Number of mutually distinct zeros = ', distinctnumber
            print *, 'Approximations for the zeros and verification:'
            print *, zeros
            print *, size(zeros)
            print *, refinednumber
            print *, size(refinement_ok)
            do j = 1, distinctnumber
                print *, "zero No. ", j
                print *, 'z = ',real(zeros(j)), aimag(zeros(j)), &
                        'f(z) = ',real(fzeros(j)), aimag(fzeros(j)), &
                        multiplicities(j)
            end do
        end if

    end subroutine



end module


SUBROUTINE FDF_WKB(Z,F,DF)
!-----------------------------------------------------------------------
!**PURPOSE
!  Define the function F and its derivative DF.
!-----------------------------------------------------------------------
    implicit none
    double complex, intent(in)   :: Z
    double complex, intent(out)  :: F, DF

    F = Z**3.0d0 - 5.0d0 * Z + 4.0d0!EXP(3.0*Z) + 2.0d0*Z*COS(Z) - 1.0d0
    DF = 3.0d0 * Z**2 - 5.0d0!3.0*EXP(3.0*Z) + 2.0d0*COS(Z) - 2.0d0*Z*SIN(Z)
!
!  Other examples
!
    !F = (Z**2)*(Z-1.0)*(Z-2.0)*(Z-3.0d0)*(Z-4.0d0) + Z*SIN(Z)
    !DF = 6.0d0*Z**5 - 50.0d0*Z**4 + 140.0d0*Z**3 - 150.0d0*Z**2 &
        !+ 48.0d0*Z + Z*COS(Z) + SIN(Z)

    !F = (Z**2)*((Z-2.0d0)**2)*(COS(Z)*EXP(2.0d0*Z)+Z**3-1.0d0-SIN(Z))
    !DF = 2.0d0*Z*((Z-2.0d0)**2)*(COS(Z)*EXP(2.0d0*Z)+Z**3-1.0d0-SIN(Z))    +  &
    !    (Z**2)*2.0d0*(Z-2.0d0)*(COS(Z)*EXP(2.0d0*Z)+Z**3-1.0d0-SIN(Z))  +  &
    !    (Z**2)*((Z-2.0d0)**2)*(-SIN(Z)*EXP(2.0d0*Z)                 +  &
    !    2.0d0*COS(Z)*EXP(2.0d0*Z)+3.0*Z**2-COS(Z))

END SUBROUTINE FDF_WKB