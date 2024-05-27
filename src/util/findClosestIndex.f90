integer function findClosestIndex(array, target)

    implicit none
    real(kind=8), intent(in) :: array(:)
    real(kind=8), intent(in) :: target
    integer :: closest_index
    real(kind=8) :: min_difference
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

    ! Return the index of the closest element
    findClosestIndex = closest_index
end function findClosestIndex

