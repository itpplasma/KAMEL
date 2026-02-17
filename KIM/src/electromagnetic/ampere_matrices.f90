module ampere_matrices_m

    use KIM_kinds_m, only: dp

    implicit none

    contains

    subroutine interpolate_equil_to_xl(hz_xl, hth_xl)
        ! Interpolate equilibrium hz, hth from equil_grid onto xl_grid%xb

        use equilibrium_m, only: hz, hth, equil_grid
        use grid_m, only: xl_grid

        implicit none

        real(dp), allocatable, intent(out) :: hz_xl(:), hth_xl(:)
        integer :: i, ir, ibeg, iend
        integer :: nlagr = 4
        integer :: nder = 0
        real(dp), allocatable :: coef(:,:)

        allocate(hz_xl(xl_grid%npts_b), hth_xl(xl_grid%npts_b))
        allocate(coef(0:nder, nlagr))

        do i = 1, xl_grid%npts_b
            call binsrc(equil_grid, 1, size(equil_grid), xl_grid%xb(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend > size(equil_grid)) then
                iend = size(equil_grid)
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, xl_grid%xb(i), equil_grid(ibeg:iend), coef)
            hz_xl(i)  = sum(coef(0,:) * hz(ibeg:iend))
            hth_xl(i) = sum(coef(0,:) * hth(ibeg:iend))
        end do

        deallocate(coef)

    end subroutine

    subroutine calc_potential_matrix(Q_mat)
        ! Compute the unweighted potential matrix Q_{l'l} = int dr (1/r) phi_l' phi_l
        ! Analytical expressions for hat basis functions (tridiagonal).

        use grid_m, only: xl_grid

        implicit none

        real(dp), intent(out) :: Q_mat(:,:)
        integer :: i, n
        real(dp) :: r_l, r_lm1, r_lp1

        n = xl_grid%npts_b
        Q_mat = 0.0d0

        do i = 2, n - 1
            r_l   = xl_grid%xb(i)
            r_lm1 = xl_grid%xb(i-1)
            r_lp1 = xl_grid%xb(i+1)

            ! Diagonal: l' = l
            Q_mat(i,i) = 0.5d0 * ( &
                (r_l**2 + 2.0d0*r_lm1**2*(log(r_l) - log(r_lm1)) &
                 - 4.0d0*r_l*r_lm1 + 3.0d0*r_lm1**2) / (r_l - r_lm1)**2 &
                - (r_l**2 + 2.0d0*r_lp1**2*(log(r_l) - log(r_lp1)) &
                 - 4.0d0*r_l*r_lp1 + 3.0d0*r_lp1**2) / (r_l - r_lp1)**2 )

            ! Off-diagonal: l' = l+1
            if (i < n) then
                Q_mat(i, i+1) = (-r_l**2 + 2.0d0*r_l*r_lp1*(log(r_l) - log(r_lp1)) &
                    + r_lp1**2) / (2.0d0 * (r_l - r_lp1)**2)
            end if

            ! Off-diagonal: l' = l-1
            if (i > 1) then
                Q_mat(i, i-1) = (r_l**2 - 2.0d0*r_l*r_lm1*(log(r_l) - log(r_lm1)) &
                    - r_lm1**2) / (2.0d0 * (r_l - r_lm1)**2)
            end if
        end do

    end subroutine

    subroutine calc_weighted_mass_matrix(M_hth, M_mat, hth_xl)
        ! M^{h_theta}_{l'l} = h_theta(r_bar) * M_{l'l}
        ! Midpoint approximation: h_theta evaluated at the midpoint of
        ! overlapping basis function support.

        use grid_m, only: xl_grid

        implicit none

        real(dp), intent(out) :: M_hth(:,:)
        real(dp), intent(in) :: M_mat(:,:)
        real(dp), intent(in) :: hth_xl(:)
        integer :: i, n

        n = xl_grid%npts_b
        M_hth = 0.0d0

        do i = 1, n
            ! Diagonal: midpoint is the node itself
            M_hth(i,i) = hth_xl(i) * M_mat(i,i)
            ! Off-diagonal: midpoint between nodes i and i+1
            if (i < n) then
                M_hth(i, i+1) = 0.5d0 * (hth_xl(i) + hth_xl(i+1)) * M_mat(i, i+1)
                M_hth(i+1, i) = 0.5d0 * (hth_xl(i) + hth_xl(i+1)) * M_mat(i+1, i)
            end if
        end do

    end subroutine

    subroutine calc_weighted_potential_matrix(Q_hz, Q_mat, hz_xl)
        ! Q^{h_z}_{l'l} = h_z(r_bar) * Q_{l'l}
        ! Midpoint approximation: h_z evaluated at the midpoint of
        ! overlapping basis function support.

        use grid_m, only: xl_grid

        implicit none

        real(dp), intent(out) :: Q_hz(:,:)
        real(dp), intent(in) :: Q_mat(:,:)
        real(dp), intent(in) :: hz_xl(:)
        integer :: i, n

        n = xl_grid%npts_b
        Q_hz = 0.0d0

        do i = 1, n
            ! Diagonal: midpoint is the node itself
            Q_hz(i,i) = hz_xl(i) * Q_mat(i,i)
            ! Off-diagonal: midpoint between nodes i and i+1
            if (i < n) then
                Q_hz(i, i+1) = 0.5d0 * (hz_xl(i) + hz_xl(i+1)) * Q_mat(i, i+1)
                Q_hz(i+1, i) = 0.5d0 * (hz_xl(i) + hz_xl(i+1)) * Q_mat(i+1, i)
            end if
        end do

    end subroutine

end module
