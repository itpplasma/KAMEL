subroutine getIfunc(x1, x2, symbI)
    integer, parameter :: mnmax = 3
    integer :: m, n
    double precision, intent(in) :: x1, x2
    double precision :: z
    double complex :: denom
    double complex, dimension(0:mnmax, 0:mnmax), intent(out) :: symbI
    double complex, dimension(0:mnmax, 0:mnmax) :: Imn

    if (.false.) then
        ! collisionless case:
        symbI = (0.d0, 0.d0)
        z = x2/(sqrt(2.d0)*x1)

        symbI(0, 0) = sqrt(2.d0)*exp(-z**2)/abs(x1)
        symbI(1, 0) = symbI(0, 0)*x2/x1
        symbI(1, 1) = symbI(1, 0)*x2/x1

        symbI(2, 0) = symbI(0, 0)*(x2/x1)**2
        symbI(2, 1) = symbI(1, 0)*(x2/x1)**2
        symbI(3, 0) = symbI(2, 1)
        symbI(3, 1) = symbI(1, 1)*(x2/x1)**2

        symbI(2, 2) = symbI(2, 0)*(x2/x1)**2
        symbI(3, 2) = symbI(2, 1)*(x2/x1)**2
        symbI(3, 3) = symbI(3, 1)*(x2/x1)**2
    else
        ! collisional case:
        call W2_arr(x1, x2, Imn)

        if (.true.) then
            ! energy conservation:
            denom = (1.d0, 0.d0) - Imn(0, 0) + (2.d0, 0.d0)*Imn(2, 0) - Imn(2, 2)
            do m = 0, 3
                do n = 0, 3
                    symbI(m, n) = Imn(m, n) + (Imn(m, 0) - Imn(m, 2))*(Imn(n, 0) - Imn(n, 2))/denom
                end do
            end do
        else
            symbI = Imn
        end if
    end if

end subroutine getIfunc
