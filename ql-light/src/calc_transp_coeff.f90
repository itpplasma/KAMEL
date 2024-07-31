subroutine calc_transport_coeffs_ornuhl(dim, vT, nu, D_11, D_12, D_21, D_22)

    use baseparam_mod
    use wave_code_data, only: om_E, kp, Es, Br, B0, r, ks, Ep
    use diag_mod, only: i_mn_loop
    use grid_mod, only: r_resonant, gg_width, rb

    integer, parameter :: mnmax = 3
    integer, intent(in) :: dim
    real(8), dimension(dim), intent(in) :: vT, nu
    real(8), dimension(dim), intent(out) :: D_11, D_12, D_21, D_22
    double precision, dimension(:), allocatable :: x1, x2, comfac, d_12a
    double precision, dimension(:), allocatable :: epm2, brm2, epbr_re, epbr_im
    double complex, dimension(:, :, :), allocatable :: symbI

    allocate (comfac(dim), d_12a(dim), epm2(dim), brm2(dim), epbr_re(dim), epbr_im(dim))
    allocate (x1(dim), x2(dim), symbI(0:mnmax, 0:mnmax, dim))
!
!    if  Br=c*kp*Es/om_E diffusion tensor iz zero

    comfac = 0.5d0/(nu*B0**2)
    epm2 = c**2*abs(Es)**2
    brm2 = vT**2*abs(Br)**2
    epbr_re = 2.d0*c*vT*real(conjg(Es)*Br)
    epbr_im = 2.d0*c*vT*dimag(conjg(Es)*Br)
!epm2=0.0d0 !c**2*abs(Es)**2
!brm2=1.0d0 !vT**2*abs(Br)**2
!epbr_re=0.0d0 !2.d0*c*vT*real(conjg(Es)*Br)
!epbr_im=0.0d0 !2.d0*c*vT*dimag(conjg(Es)*Br)

    x1 = kp*vT/nu
    x2 = -om_E/nu
    symbI = (0.d0, 0.d0)
    do i = 1, dim
        if (rb(i) .lt. r_resonant(i_mn_loop) - 2.d0*gg_width) cycle
        if (rb(i) .gt. r_resonant(i_mn_loop) + 2.d0*gg_width) cycle
        call getIfunc(x1(i), x2(i), symbI(:, :, i))
    end do

    D_11 = comfac*(epm2*real(symbI(0, 0, :)) &
                   + epbr_re*real(symbI(1, 0, :)) &
                   + brm2*real(symbI(1, 1, :)))
    D_12 = comfac*(epm2*real(symbI(0, 0, :) + 0.5d0*symbI(2, 0, :)) &
                   + epbr_re*real(symbI(1, 0, :) + 0.25d0*(symbI(3, 0, :) + symbI(2, 1, :))) &
                   + brm2*real(symbI(1, 1, :) + 0.5d0*symbI(3, 1, :)))
    D_21 = D_12
    D_22 = comfac*(epm2*real(2.d0*symbI(0, 0, :) + symbI(2, 0, :) &
                             + 0.25d0*symbI(2, 2, :)) &
                   + epbr_re*real(2.d0*symbI(1, 0, :) &
                                  + 0.5d0*(symbI(3, 0, :) + symbI(2, 1, :)) &
                                  + 0.25d0*symbI(3, 2, :)) &
                   + brm2*real(2.d0*symbI(1, 1, :) + symbI(3, 1, :) &
                               + 0.25d0*symbI(3, 3, :)))

    D_12a = comfac*epbr_im*0.25d0*dimag(symbI(2, 1, :) - symbI(3, 0, :))

    D_12 = D_12 + D_12a
    D_21 = D_21 - D_12a

!D_11=comfac*epm2*real(symbI(0,0,:))
!D_12=comfac*epbr_re*real(symbI(1,0,:))
!D_22=comfac*brm2*real(symbI(1,1,:))
!do i=1,dim
!write(7000,*) r(i),D_11(i),D_12(i),D_22(i),D_11(i)+D_12(i)+D_22(i)
!enddo
!D_11=comfac*real(symbI(0,0,:))
!D_12=comfac*real(symbI(1,0,:))
!D_22=comfac*real(symbI(1,1,:))
!do i=1,dim
!write(7001,*) r(i),D_11(i),D_12(i),D_22(i),x1(i),x2(i)
!enddo
!stop

    deallocate (x1, x2, symbI)
    deallocate (comfac, d_12a, epm2, brm2, epbr_re, epbr_im)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
