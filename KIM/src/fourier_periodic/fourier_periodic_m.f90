module fourier_periodic_m
    ! Discrete-Fourier solve on the periodized background (issue #175):
    ! modes k_m = 2 pi m / L with L = 2 (dr_layer + dr_transition), matrix
    ! elements K_{m,m'} = (2 pi / L) * one-period r_g quadrature of the
    ! restored kernels, screened Poisson system
    !     (k_m^2 delta_{m,m'} - 4 pi K^{rho Phi}_{m,m'}) Phi_{m'}
    !         = 4 pi (K^{rho B} B_r)_m + s_m
    ! solved densely with LAPACK zgesv, and the inverse expansion
    ! f(r) = sum_m f_m exp(i k_m r). Absolute-coordinate phases match the
    ! exp(i (k_r - k_r') r_g) factor carried by the kernels, so a uniform
    ! background gives an exactly diagonal matrix through the equidistant
    ! discrete orthogonality of the quadrature.

    use KIM_kinds_m, only: dp
    use constants_m, only: pi, com_unit

    implicit none

contains

    ! Modes k_m = 2 pi m / L for m = -n_modes .. n_modes.
    function periodic_k_grid(period, n_modes) result(k_modes)

        implicit none

        real(dp), intent(in) :: period
        integer, intent(in) :: n_modes
        real(dp) :: k_modes(2*n_modes + 1)

        integer :: m

        do m = -n_modes, n_modes
            k_modes(m + n_modes + 1) = 2.0d0*pi*m/period
        end do

    end function periodic_k_grid

    ! One-period equidistant quadrature of a kernel function
    ! G(k_m, k_m', r_g): K_{m,m'} = (2 pi / L) sum_g G dr_g with the
    ! midpoint rule, which is spectrally accurate for the periodic
    ! integrand of the periodized background.
    function assemble_periodic_kernel(kernel, k_modes, r_mid, period, &
                                      n_quad) result(kernel_matrix)

        implicit none

        interface
            complex(dp) function kernel(val_kr, val_krp, val_rg)
                import :: dp
                real(dp), intent(in) :: val_kr, val_krp, val_rg
            end function kernel
        end interface
        real(dp), intent(in) :: k_modes(:), r_mid, period
        integer, intent(in) :: n_quad
        complex(dp) :: kernel_matrix(size(k_modes), size(k_modes))

        real(dp) :: r_g, dr_g
        integer :: m, mp, ig

        dr_g = period/n_quad
        kernel_matrix = (0.0d0, 0.0d0)
        do ig = 1, n_quad
            r_g = r_mid - 0.5d0*period + (ig - 0.5d0)*dr_g
            do mp = 1, size(k_modes)
                do m = 1, size(k_modes)
                    kernel_matrix(m, mp) = kernel_matrix(m, mp) &
                        + kernel(k_modes(m), k_modes(mp), r_g)
                end do
            end do
        end do
        kernel_matrix = kernel_matrix*2.0d0*pi*dr_g/period

    end function assemble_periodic_kernel

    ! Dense solve of (k_m^2 delta - 4 pi K^{rho Phi}) Phi = rhs.
    function solve_periodic_potential(k_modes, kernel_rho_phi_matrix, &
                                      rhs) result(phi_modes)

        implicit none

        real(dp), intent(in) :: k_modes(:)
        complex(dp), intent(in) :: kernel_rho_phi_matrix(:, :)
        complex(dp), intent(in) :: rhs(:)
        complex(dp) :: phi_modes(size(rhs))

        complex(dp) :: system(size(rhs), size(rhs))
        integer :: ipiv(size(rhs)), info, m

        system = -4.0d0*pi*kernel_rho_phi_matrix
        do m = 1, size(k_modes)
            system(m, m) = system(m, m) + k_modes(m)**2
        end do
        phi_modes = rhs
        call zgesv(size(rhs), 1, system, size(rhs), ipiv, phi_modes, &
                   size(rhs), info)
        if (info /= 0) then
            print *, "fourier_periodic_m: zgesv failed with info = ", info
            error stop
        end if

    end function solve_periodic_potential

    ! f(r) = sum_m f_m exp(i k_m r).
    complex(dp) function periodic_field_value(k_modes, field_modes, r)

        implicit none

        real(dp), intent(in) :: k_modes(:), r
        complex(dp), intent(in) :: field_modes(:)

        periodic_field_value = sum(field_modes &
                                   *exp(com_unit*k_modes*r))

    end function periodic_field_value

end module fourier_periodic_m
