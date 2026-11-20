!> The anomalous diffusion profile Da(r) from the input. It does not change during a run, but
!> calc_equil_diffusion_coeffs runs at every right-hand-side evaluation, so the profile is read
!> once and kept.
module da_profile_m

    use QLBalance_kinds, only: dp

    implicit none

    private
    public :: r_raw, Da_raw, load_da_profile, forget_da_profile

    real(dp), allocatable :: r_raw(:), Da_raw(:)

contains

    !> Read Da from path2inp (HDF5 input) or profiles/Da.dat, unless it is already loaded.
    !> Profile from outside this code by Philipp Ulbl (2020), HDF5 input by Markus Markl (2021).
    subroutine load_da_profile

        use h5mod
        use control_mod, only: ihdf5IO

        character(1024) :: fname
        integer :: nr, i, iunit_res
        integer :: lb, ub

        if (allocated(r_raw)) return

        if (ihdf5IO .eq. 1) then
            fname = "/da_estimation/"
            CALL h5_init()
            call h5_open(path2inp, h5_id)
            CALL h5_get_bounds_1(h5_id, trim(fname)//'Da', lb, ub)
            allocate (r_raw(ub), Da_raw(ub))
            CALL h5_get_double_1(h5_id, trim(fname)//'Da', Da_raw)
            CALL h5_get_double_1(h5_id, trim(fname)//'r', r_raw)
            CALL h5_close(h5_id)
            CALL h5_deinit()
        else
            fname = 'profiles/Da.dat'
            nr = 0
            open (newunit=iunit_res, file=fname)
            do
                read (iunit_res, *, end=1)
                nr = nr + 1
            end do
1           continue
            close (iunit_res)
            allocate (r_raw(nr), Da_raw(nr))
            open (newunit=iunit_res, file=fname)
            do i = 1, nr
                read (iunit_res, *) r_raw(i), Da_raw(i)
            end do
            close (iunit_res)
        end if

    end subroutine load_da_profile

    !> Drop the loaded profile, so the next load_da_profile reads the input again.
    subroutine forget_da_profile
        if (allocated(r_raw)) deallocate (r_raw)
        if (allocated(Da_raw)) deallocate (Da_raw)
    end subroutine forget_da_profile

end module da_profile_m

subroutine calc_equil_diffusion_coeffs

    use grid_mod, only: npoib, dae11, dae12, dae22, dai11, dai12, dai22 &
                        , dni22, visca, rb, cneo
    use plasma_parameters, only: params_b
    use baseparam_mod, only: rsepar
    use paramscan_mod, only: viscosity_factor
    use PolyLagrangeInterpolation
    use QLBalance_kinds, only: dp
    use da_profile_m, only: r_raw, Da_raw, load_da_profile

    implicit none

    integer :: ipoi
    real(dp) :: weight

    integer :: nr, i
    real(dp) :: r
    integer :: ind_begin_interp, ind_end_interp

    call load_da_profile
    nr = size(r_raw)
    !code from amn_of_r to interpolate (+do loop)

    !interpolate on balance grid
    if (nr < nlagr) error stop "Da grid has fewer points than the interpolation stencil"
    if (.not. allocated(coef)) allocate (coef(0:nder, nlagr))

    do ipoi = 1, npoib
        r = rb(ipoi);
        if (r .gt. r_raw(nr)) then
            r = r_raw(nr)
        end if

        !binsearch
        call binsrc(r_raw, 1, nr, r, i)
        ! A stencil on the Da grid r_raw, which may be shorter (e.g. ending at the
        ! separatrix) or longer than the balance grid rb.
        call get_lagr_stencil(i, nr, ind_begin_interp, ind_end_interp)
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
    ! D11 = D_perp
    ! D12 = 1.499999 D11 on purpose, not exactly 3/2: with 3/2 the thermodiffusive
    ! terms of the anomalous particle flux cancel and their Jacobian entries become
    ! exactly zero. rhs_balance takes the sparsity pattern from the values
    ! (abs(dy) > tiny), so its count pass and fill pass then disagree and the run
    ! stops. The deviation from 3/2 is 7e-7 relative and leaves a negligible
    ! non-ambipolar remainder in the anomalous flux.
    dae11 = dae12/1.499999d0
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

    integer :: i, i_res
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

    ! DIAGNOSTIC: print D_22 terms at resonance
    ! Find grid index closest to resonance
    i_res = 1
    do i = 2, dim
        if (abs(rb(i) - r_resonant(i_mn_loop)) .lt. abs(rb(i_res) - r_resonant(i_mn_loop))) i_res = i
    end do
    write(*,*) '--- D_22 terms at resonance (i=', i_res, ', r=', rb(i_res), ') ---'
    write(*,*) '  |Es| at res    = ', abs(Es(i_res))
    write(*,*) '  |Br| at res    = ', abs(Br(i_res))
    write(*,*) '  B0  at res     = ', B0(i_res)
    write(*,*) '  nu  at res     = ', nu(i_res)
    write(*,*) '  vT  at res     = ', vT(i_res)
    write(*,*) '  kp  at res     = ', kp(i_res)
    write(*,*) '  om_E at res    = ', om_E(i_res)
    write(*,*) '  comfac at res  = ', comfac(i_res)
    write(*,*) '  epm2 at res    = ', epm2(i_res), '  (c^2*|Es|^2)'
    write(*,*) '  brm2 at res    = ', brm2(i_res), '  (vT^2*|Br|^2)'
    write(*,*) '  epbr_re at res = ', epbr_re(i_res)
    write(*,*) '  --- Approximate D_22 ~ comfac*(epm2 + brm2) ---'
    write(*,*) '  epm2 contrib   = ', comfac(i_res)*epm2(i_res)
    write(*,*) '  brm2 contrib   = ', comfac(i_res)*brm2(i_res)
    write(*,*) ''

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
