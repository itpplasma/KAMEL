program test_zeal

    USE Precision_Module
    USE Zeal_Module

    implicit none

    integer :: J 
    logical :: OUTPUT

    integer :: TOTALNUMBER, DISTINCTNUMBER, REFINEDNUMBER 
    integer, dimension(:), POINTER            :: MULTIPLICITIES
    logical, dimension(:), POINTER            :: REFINEMENT_OK
    COMPLEX(KIND=DP), DIMENSION(:), POINTER   :: ZEROS, FZEROS

    CALL ZEAL(TOTALNUMBER,DISTINCTNUMBER,ZEROS,FZEROS,   &
                MULTIPLICITIES,REFINEDNUMBER,REFINEMENT_OK)

    OUTPUT = .true.
    IF ( OUTPUT ) THEN
        PRINT *
        PRINT *
        PRINT *, 'Total number of zeros = ', TOTALNUMBER
        PRINT *
        PRINT *, 'Number of mutually distinct zeros = ', DISTINCTNUMBER
        PRINT *
        PRINT *, 'Approximations for the zeros and verification:'
        DO J = 1, DISTINCTNUMBER
        PRINT 1, REAL(ZEROS(J),DP),  AIMAG(ZEROS(J)),  &
                REAL(FZEROS(J),DP), AIMAG(FZEROS(J)), &
                MULTIPLICITIES(J)
        IF ( .NOT. REFINEMENT_OK(J) ) PRINT 2
        END DO 
    END IF

    1 FORMAT (/3X, 'z    = (', G22.15, ',', G22.15, ' )', &
            /3X, 'f(z) = (', G22.15, ',', G22.15, ' )', &
            /3X, 'multiplicity = ', I5)
    2 FORMAT (3X,'refinement not succesful')


end program