module initialize

    implicit none


    contains 

    subroutine init

        implicit none

        call gengrid(npoimin)

    end subroutine


end module