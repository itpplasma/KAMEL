module findIndex_m
    contains

    subroutine findClosestIndex(array, target, closest_index)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), dimension(:), intent(in) :: array
        real(dp), intent(in) :: target
        integer, intent(out) :: closest_index
        real(dp) :: min_difference
        integer :: i

        ! Initialize with a large value
        min_difference = huge(1.0)

        ! Loop through the array to find the closest element
        do i = 1, size(array)
            if (abs(array(i) - target) < min_difference) then
                min_difference = abs(array(i) - target)
                closest_index = i
            end if
        end do

    end subroutine

end module