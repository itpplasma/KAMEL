!> calc_equil_diffusion_coeffs interpolates Da from profiles/Da.dat onto rb. The
!> interpolation stencil must stay inside the Da grid, whether that grid ends
!> before rb (Da is then held at its last value) or has more points than rb.
!> The Da input is read once and reused until forget_da_profile is called.
program test_diff_coeffs
    use QLBalance_kinds, only: dp
    use control_mod, only: ihdf5IO
    use paramscan_mod, only: viscosity_factor
    use PolyLagrangeInterpolation, only: get_lagr_stencil
    use da_profile_m, only: forget_da_profile
    implicit none

    ihdf5IO = 0
    viscosity_factor = 1.0_dp

    call check_stencil(5, 10, 3, 6, "interior")
    call check_stencil(1, 10, 1, 4, "first point")
    call check_stencil(10, 10, 7, 10, "last point")

    ! Da grid ends at r = 40 inside rb, which reaches 60.
    call check_interpolation(nr=40, r_end=40.0_dp, npoib=50, label="Da grid shorter than rb")
    ! Da grid with more points than rb, covering all of it.
    call check_interpolation(nr=400, r_end=70.0_dp, npoib=50, label="Da grid longer than rb")

    call check_da_read_once

contains

    !> A smooth Da that a 4-point Lagrange stencil reproduces only near its nodes.
    elemental function da_model(r) result(da)
        real(dp), intent(in) :: r
        real(dp) :: da
        da = 1.0e4_dp * (2.0_dp + sin(r / 5.0_dp))
    end function da_model

    subroutine check_stencil(i, n, expected_begin, expected_end, label)
        integer, intent(in) :: i, n, expected_begin, expected_end
        character(len=*), intent(in) :: label
        integer :: ind_begin, ind_end

        call get_lagr_stencil(i, n, ind_begin, ind_end)
        if (ind_begin /= expected_begin .or. ind_end /= expected_end) then
            print '(a,a,4(a,i0))', "FAIL: stencil at ", label, ": ", ind_begin, ":", ind_end, &
                ", expected ", expected_begin, ":", expected_end
            error stop 1
        end if
    end subroutine check_stencil

    subroutine check_interpolation(nr, r_end, npoib, label)
        use grid_mod, only: npoib_ => npoib, rb, dae11, dae12, dae22, dai11, dai12, dai22, &
                            dni22, visca, cneo
        use plasma_parameters, only: params_b
        integer, intent(in) :: nr, npoib
        real(dp), intent(in) :: r_end
        character(len=*), intent(in) :: label
        real(dp), parameter :: r_start = 1.0_dp, rb_end = 60.0_dp, tolerance = 1.0e-4_dp
        real(dp) :: r_raw(nr), expected(npoib)
        integer :: i

        r_raw = [(r_start + (r_end - r_start) * (i - 1) / (nr - 1), i = 1, nr)]
        call write_da_file(r_raw, da_model(r_raw))
        call forget_da_profile

        npoib_ = npoib
        if (allocated(rb)) deallocate (rb, dae11, dae12, dae22, dai11, dai12, dai22, dni22, &
                                       visca, cneo, params_b)
        allocate (rb(npoib), dae11(npoib), dae12(npoib), dae22(npoib), dai11(npoib), &
                  dai12(npoib), dai22(npoib), dni22(npoib), visca(npoib), cneo(npoib), &
                  params_b(4, npoib))
        rb = [(r_start + (rb_end - r_start) * (i - 1) / (npoib - 1), i = 1, npoib)]
        cneo = 1.0_dp
        params_b = 1.0_dp

        call calc_equil_diffusion_coeffs

        ! Beyond the Da grid, Da keeps its value at the last point.
        expected = da_model(min(rb, r_raw(nr)))
        if (any(abs(dae12 - expected) > tolerance * abs(expected))) then
            i = maxloc(abs(dae12 - expected), 1)
            print '(a,a)', "FAIL: Da interpolation, ", label
            print '(a,f8.3,2(a,es14.6))', "  at r = ", rb(i), ": ", dae12(i), &
                ", expected ", expected(i)
            error stop 1
        end if
    end subroutine check_interpolation

    subroutine write_da_file(r, da)
        real(dp), intent(in) :: r(:), da(:)
        integer :: i, unit

        open (newunit=unit, file="profiles/Da.dat", status="replace", action="write")
        do i = 1, size(r)
            write (unit, *) r(i), da(i)
        end do
        close (unit)
    end subroutine write_da_file

    !> calc_equil_diffusion_coeffs runs at every right-hand-side evaluation, so
    !> it must not re-read the Da input each time.
    subroutine check_da_read_once
        use grid_mod, only: dae12
        integer, parameter :: nr = 400, npoib = 50
        real(dp), parameter :: r_end = 70.0_dp
        real(dp) :: r_raw(nr), first(npoib)
        integer :: i

        call check_interpolation(nr=nr, r_end=r_end, npoib=npoib, label="before rewrite")
        first = dae12

        r_raw = [(1.0_dp + (r_end - 1.0_dp) * (i - 1) / (nr - 1), i = 1, nr)]
        call write_da_file(r_raw, 2.0_dp * da_model(r_raw))

        call calc_equil_diffusion_coeffs
        if (any(dae12 /= first)) then
            print '(a)', "FAIL: Da was re-read although the profile was already loaded"
            error stop 1
        end if

        call forget_da_profile
        call calc_equil_diffusion_coeffs
        if (any(abs(dae12 - 2.0_dp * first) > 1.0e-12_dp * abs(first))) then
            print '(a)', "FAIL: forget_da_profile did not make the next call reload Da"
            error stop 1
        end if
    end subroutine check_da_read_once

end program test_diff_coeffs
