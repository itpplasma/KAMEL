subroutine calc_equil_diffusion_coeffs
    
    use grid_mod, only: npoib, dae11, dae12, dae22, dai11, dai12, dai22 &
                        , dni22, visca, rb, cneo
    use plasma_parameters, only: params_b
    use baseparam_mod, only: rsepar
    use h5mod
    use control_mod, only: ihdf5IO, debug_mode
    use paramscan_mod, only: viscosity_factor
    use PolyLagrangeInterpolation
    use QLBalance_kinds, only: dp
    
    implicit none
    
    integer :: ipoi
    real(dp) :: weight

    character(1024) :: fname;
    integer :: nr, i, iunit_res
    real(dp) :: r
    real(dp), dimension(:), allocatable :: r_raw
    real(dp), dimension(:), allocatable :: Da_raw
    integer :: lb, ub, ind_begin_interp, ind_end_interp

    !This subroutine was changed to include estimated Da from outside of this code
    !file Da.dat which is located in profiles is read and the other Da.. calculated based on this
    !Changed by Philipp Ulbl 18.05.2020
    ! Added the option to read from hdf5 file by Markus Markl 12.04.2021

    if (ihdf5IO .eq. 1) then
        ! read Da data from hdf5 input file
        fname = "/da_estimation/" ! fname is used for the group name in the hdf5 version
        CALL h5_init()
        call h5_open(path2inp, h5_id)
        CALL h5_get_bounds_1(h5_id, trim(fname)//'Da', lb, ub)
        allocate (r_raw(ub), Da_raw(ub))
        CALL h5_get_double_1(h5_id, trim(fname)//'Da', Da_raw)
        CALL h5_get_double_1(h5_id, trim(fname)//'r', r_raw)
        CALL h5_close(h5_id)
        CALL h5_deinit()
        nr = ub
    else
        !code from gengrid to read Da file
        fname = 'profiles/Da.dat';
        iunit_res = 157
        nr = 0
        open (iunit_res, file=fname)
        do
            read (iunit_res, *, end=1)
            nr = nr + 1
        end do
1       continue
        close (iunit_res)
        allocate (r_raw(nr), Da_raw(nr))
        !
        open (iunit_res, file=fname)
        do i = 1, nr
            read (iunit_res, *) r_raw(i), Da_raw(i)
        end do
        close (iunit_res)
    end if
    !code from amn_of_r to interpolate (+do loop)

    !interpolate on balance grid
    if (.not. allocated(coef)) allocate (coef(0:nder, nlagr))

    do ipoi = 1, npoib
        r = rb(ipoi);
        if (r .gt. r_raw(nr)) then
            r = r_raw(nr)
        end if

        !binsearch
        call binsrc(r_raw, 1, nr, r, i)
        call get_ind_Lagr_interp(i, ind_begin_interp, ind_end_interp)
        !lagrange interpolation with order 4 only for function (0)

        call plag_coeff(nlagr, nder, r, r_raw(ind_begin_interp:ind_end_interp), coef)
        !
        !Da estimated is dae12 -> see notes on conversion
        dae12(ipoi) = sum(coef(0, :)*Da_raw(ind_begin_interp:ind_end_interp))
        !call localizer(-1.d0, rsepar, rsepar + 0.5d0, rb(ipoi), weight)
        !dae12(ipoi) = dae12(ipoi)*(1.d0 - weight) + weight*1d6
        !write(77, *) r, dae12(ipoi)
    end do
    !call localizer(-1.d0, rsepar, rsepar + 0.5d0, rb(ipoi), weight)
    !dae12(ipoi) = dae12(ipoi)*(1.d0 - weight) + weight*1d6

    !get other da
    ! [Heyn2014 below (68)]
    ! D11 = D_perb
    dae11 = dae12/1.499999d0 !previously used instead of 1.5d0, no idea why
    dae22 = 3.75d0*dae11
    dai11 = dae11
    dai12 = dae12
    dai22 = dae22
    visca = dae11 * viscosity_factor
    
    dni22 = cneo * params_b(1, :)/sqrt(abs(params_b(4, :)))

end subroutine calc_equil_diffusion_coeffs


subroutine calc_transport_coeffs_collisionless(dim, vT, D_11, D_12, D_22)

    use baseparam_mod
    use wave_code_data, only: om_E, kp, ks, Ep, Br
    use QLBalance_kinds, only: dp
    implicit none

    integer, intent(in) :: dim
    real(dp), dimension(dim), intent(in) :: vT
    real(dp), dimension(dim), intent(out) :: D_11, D_12, D_22
    real(dp), dimension(dim) :: Z
    real(dp), dimension(dim) :: field_fac

    Z = om_E/sqrt(2.0d0)/kp/vT

    field_fac = abs(c*ks*Ep - om_E*Br)**2

    D_11 = sqrt(pi)*vT**2.0d0/btor**2.0d0*(abs(Z/om_E))**3.0d0*exp(-Z**2.0d0)*field_fac
    D_12 = (1.0d0 + Z**2.0d0)*D_11
    D_22 = (1.0d0 + (1.0d0 + Z**2.0d0)**2.0d0)*D_11

end subroutine


subroutine calc_transport_coeffs_ornuhl(dim, vT, nu, D_11, D_12, D_21, D_22)

    use baseparam_mod
    use wave_code_data, only: om_E, kp, Es, Br, B0
    use QLbalance_diag, only: i_mn_loop
    use grid_mod, only: r_resonant, gg_width, rb
    use QLBalance_kinds, only: dp

    implicit none

    integer :: i
    integer, parameter :: mnmax = 3
    integer, intent(in) :: dim
    real(dp), dimension(dim), intent(in) :: vT, nu
    real(dp), dimension(dim), intent(out) :: D_11, D_12, D_21, D_22
    real(dp), dimension(:), allocatable :: x1, x2, comfac, d_12a
    real(dp), dimension(:), allocatable :: epm2, brm2, epbr_re, epbr_im
    complex(dp), dimension(:, :, :), allocatable :: symbI

    allocate (comfac(dim), d_12a(dim), epm2(dim), brm2(dim), epbr_re(dim), epbr_im(dim))
    allocate (x1(dim), x2(dim), symbI(0:mnmax, 0:mnmax, dim))

    symbI = 0.d0
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

    deallocate (x1, x2, symbI)
    deallocate (comfac, d_12a, epm2, brm2, epbr_re, epbr_im)

end subroutine calc_transport_coeffs_ornuhl
