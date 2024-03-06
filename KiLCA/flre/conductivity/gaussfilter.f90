subroutine smooth_array_gauss(dimx, ngauss, y)

implicit none

integer :: dimx, ngauss
real(8), dimension(dimx) :: y
real(8), dimension(:), allocatable :: wgauss, ys
real(8) :: sgauss, dummy
integer :: i, j, k

sgauss=dfloat(ngauss)/5.d0

allocate(wgauss(-ngauss:ngauss), ys(dimx))

do i=-ngauss,ngauss
    wgauss(i)=exp(-(dfloat(i)/sgauss)**2)
enddo

dummy=sum(wgauss)
wgauss=wgauss/dummy

do i=ngauss+1,dimx-ngauss
    ys(i)=sum(wgauss*y(i-ngauss:i+ngauss))
enddo

do i=1,ngauss
    ys(i)=sum(wgauss(1-i:ngauss)*y(1:i+ngauss))/sum(wgauss(1-i:ngauss))
enddo

do i=dimx-ngauss+1,dimx
    ys(i)=sum(wgauss(-ngauss:dimx-i)*y(i-ngauss:dimx))/sum(wgauss(-ngauss:dimx-i))
enddo

y=ys

deallocate(wgauss, ys)

end subroutine
