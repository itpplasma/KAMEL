subroutine smooth_array_gauss(dimx, ngauss, y, ys)

    implicit none

    integer :: dimx, ngauss, i
    real(8) :: sgauss, dummy
    real(8), dimension(dimx) :: y, ys
    double precision, dimension(:), allocatable :: wgauss
    real(8) :: sum

    sgauss = dfloat(ngauss)/5.d0
    allocate (wgauss(-ngauss:ngauss))
    do i = -ngauss, ngauss
        wgauss(i) = exp(-(dfloat(i)/sgauss)**2)
    end do
    dummy = sum(wgauss)
    wgauss = wgauss/dummy

    do i = ngauss + 1, dimx - ngauss
        ys(i) = sum(wgauss*y(i - ngauss:i + ngauss))
    end do
    do i = 1, ngauss
        ys(i) = sum(wgauss(1 - i:ngauss)*y(1:i + ngauss))/sum(wgauss(1 - i:ngauss))
    end do
    do i = dimx - ngauss + 1, dimx
        ys(i) = sum(wgauss(-ngauss:dimx - i)*y(i - ngauss:dimx))/sum(wgauss(-ngauss:dimx - i))
    end do

end subroutine
