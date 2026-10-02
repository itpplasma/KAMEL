module ampere_operator_m
    use KIM_kinds_m, only: dp
    implicit none
    private
    public :: assemble_ampere_weak

contains

    subroutine assemble_ampere_weak(radius, ks, matrix)
        ! Unweighted hat-basis moments of A'' + A'/r - ks**2 A.
        ! Current kernels already contain the same dr test-function integral.
        ! Boundary rows are imposed by the coupled solver after assembly.
        real(dp), intent(in) :: radius(:), ks(:)
        real(dp), intent(out) :: matrix(:, :)
        real(dp), parameter :: nodes(3) = &
            [0.11270166537925831148_dp, 0.5_dp, 0.88729833462074168852_dp]
        real(dp), parameter :: weights(3) = &
            [5.0_dp/18.0_dp, 4.0_dp/9.0_dp, 5.0_dp/18.0_dp]
        real(dp) :: h, r, k, phi(2), deriv(2), block(2, 2), weight
        integer :: n, element, point, i, j

        n = size(radius)
        if (n < 2) error stop 'Ampere grid needs at least two points'
        if (size(ks) /= n) error stop 'Ampere wavenumber shape mismatch'
        if (size(matrix, 1) /= n) error stop 'Ampere matrix row mismatch'
        if (size(matrix, 2) /= n) error stop 'Ampere matrix column mismatch'
        if (any(radius < 0.0_dp)) error stop 'Ampere radius must be nonnegative'
        if (any(radius(2:n) <= radius(1:n-1))) error stop 'Ampere grid must increase'

        matrix = 0.0_dp
        do element = 1, n-1
            h = radius(element+1) - radius(element)
            deriv = [-1.0_dp/h, 1.0_dp/h]
            do j = 1, 2
                do i = 1, 2
                    block(i, j) = -h*deriv(i)*deriv(j)
                end do
            end do
            do point = 1, 3
                phi = [1.0_dp-nodes(point), nodes(point)]
                r = sum(phi*radius(element:element+1))
                k = sum(phi*ks(element:element+1))
                weight = h*weights(point)
                do j = 1, 2
                    do i = 1, 2
                        block(i, j) = block(i, j) + weight* &
                            (phi(i)*deriv(j)/r - k*k*phi(i)*phi(j))
                    end do
                end do
            end do
            matrix(element:element+1, element:element+1) = &
                matrix(element:element+1, element:element+1) + block
        end do
    end subroutine assemble_ampere_weak
end module ampere_operator_m
