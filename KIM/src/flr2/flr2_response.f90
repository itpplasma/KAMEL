module flr2_response_m
    ! Second-order local response adapted from KiLCA-FLR2 response_current.f90.
    ! Profile preparation is intentionally absent: KIM's plasma_t is the single
    ! source of equilibrium, collision, force, FLR, and susceptibility data.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use constants_m, only: com_unit, e_charge, pi
    use species_m, only: plasma_t
    use flr2_tridiagonal_m, only: solve_flr2_tridiagonal, FLR2_SOLVE_OK
    implicit none
    private

    integer, parameter, public :: FLR2_RESPONSE_OK = 0
    integer, parameter, public :: FLR2_RESPONSE_BAD_INPUT = 10

    type, public :: flr2_options_t
        logical :: electron_flr = .true.
        logical :: ion_flr = .true.
        logical :: electron_potential = .true.
        logical :: ion_potential = .true.
        logical :: electron_current = .true.
        logical :: ion_current = .true.
        logical :: include_potential_in_current = .true.
    end type flr2_options_t

    public :: solve_flr2_response

contains

    subroutine solve_flr2_response(background, mpol, ntor, major_radius, &
                                   toroidal_field, bpsi_over_bphi, options, &
                                   phi, parcur_over_b0, stat)
        type(plasma_t), intent(in) :: background
        integer, intent(in) :: mpol, ntor
        real(dp), intent(in) :: major_radius
        real(dp), intent(in) :: toroidal_field(:)
        complex(dp), intent(in) :: bpsi_over_bphi(:)
        type(flr2_options_t), intent(in) :: options
        complex(dp), intent(out) :: phi(:), parcur_over_b0(:)
        integer, intent(out) :: stat

        integer :: i, include_phi, n, solver_stat, sp
        real(dp) :: charge, current_switch, flr_switch, potential_switch
        real(dp) :: dphi_dr, mode_factor, rho2_factor
        complex(dp) :: factor_of_phi
        complex(dp) :: f_species, g_species, h_species
        complex(dp) :: f_current, f_potential
        complex(dp) :: g_current, g_potential, h_potential
        complex(dp), allocatable :: a2_in(:), a2_out(:), a0(:)
        complex(dp), allocatable :: b2_in(:), b2_out(:), b0(:)
        complex(dp), allocatable :: c2_in(:), c2_out(:), c0(:)
        complex(dp), allocatable :: d2_in(:), d2_out(:), d0(:)

        phi = (0.0_dp, 0.0_dp)
        parcur_over_b0 = (0.0_dp, 0.0_dp)
        stat = FLR2_RESPONSE_BAD_INPUT
        n = background%grid_size

        if (.not. valid_background(background, n)) return
        if (.not. ieee_is_finite(major_radius) .or. major_radius <= 0.0_dp) return
        if (size(toroidal_field) /= n .or. size(bpsi_over_bphi) /= n &
            .or. size(phi) /= n .or. size(parcur_over_b0) /= n) return
        if (.not. all(ieee_is_finite(toroidal_field))) return
        if (.not. finite_complex_1d(bpsi_over_bphi)) return

        allocate(a2_in(n), a2_out(n), a0(n))
        allocate(b2_in(n), b2_out(n), b0(n))
        allocate(c2_in(n), c2_out(n), c0(n))
        allocate(d2_in(n), d2_out(n), d0(n))

        do i = 1, n
            mode_factor = real(mpol, dp) + background%q(i)*real(ntor, dp)
            dphi_dr = -background%Er(i)
            if (abs(dphi_dr) <= tiny(1.0_dp) &
                .or. abs(background%q(i)) <= tiny(1.0_dp) &
                .or. abs(background%om_E(i)) <= tiny(1.0_dp) &
                .or. abs(toroidal_field(i)) <= tiny(1.0_dp)) return

            factor_of_phi = com_unit*mode_factor/(background%q(i)*dphi_dr)
            f_potential = (0.0_dp, 0.0_dp)
            g_potential = (0.0_dp, 0.0_dp)
            h_potential = (0.0_dp, 0.0_dp)
            f_current = (0.0_dp, 0.0_dp)
            g_current = (0.0_dp, 0.0_dp)

            do sp = 0, background%n_species - 1
                if (background%spec(sp)%nu(i) <= 0.0_dp) return

                if (sp == 0) then
                    flr_switch = logical_factor(options%electron_flr)
                    potential_switch = logical_factor(options%electron_potential)
                    current_switch = logical_factor(options%electron_current)
                else
                    flr_switch = logical_factor(options%ion_flr)
                    potential_switch = logical_factor(options%ion_potential)
                    current_switch = logical_factor(options%ion_current)
                end if

                charge = real(background%spec(sp)%Zspec, dp)*e_charge
                rho2_factor = (background%spec(sp)%rho_L(i)/major_radius)**2*flr_switch

                f_species = charge*background%spec(sp)%n(i)*background%spec(sp)%vT(i)**2 &
                            /background%spec(sp)%nu(i) &
                            *(background%spec(sp)%I11(i, 0) &
                              *(background%spec(sp)%A1(i) + background%spec(sp)%A2(i)) &
                              + 0.5_dp*background%spec(sp)%I13(i, 0)*background%spec(sp)%A2(i))
                g_species = 0.5_dp*rho2_factor*charge*background%spec(sp)%n(i) &
                            *background%spec(sp)%vT(i)**2/background%spec(sp)%nu(i) &
                            *(background%spec(sp)%I11(i, 0) &
                              *(background%spec(sp)%A1(i) + 2.0_dp*background%spec(sp)%A2(i)) &
                              + 0.5_dp*background%spec(sp)%I13(i, 0)*background%spec(sp)%A2(i))
                h_species = 0.5_dp*rho2_factor*charge*background%spec(sp)%n(i) &
                            *major_radius**2/dphi_dr &
                            *(background%spec(sp)%A1(i) + 2.5_dp*background%spec(sp)%A2(i))

                f_potential = f_potential + potential_switch*f_species
                g_potential = g_potential + potential_switch*g_species
                h_potential = h_potential + potential_switch*h_species
                f_current = f_current + current_switch*f_species
                g_current = g_current + current_switch*g_species
            end do

            b0(i) = -4.0_dp*pi*mode_factor*f_potential &
                    /(background%q(i)*background%om_E(i)*major_radius**2)
            b2_in(i) = -4.0_dp*pi*mode_factor*g_potential &
                       /(background%q(i)*background%om_E(i))
            b2_out(i) = b2_in(i)

            a0(i) = -b0(i)*factor_of_phi
            a2_in(i) = -b2_in(i)*factor_of_phi + 4.0_dp*pi*h_potential
            a2_out(i) = cmplx(1.0_dp, 0.0_dp, dp) + a2_in(i)

            d0(i) = -f_current/(toroidal_field(i)*major_radius)
            d2_in(i) = -g_current*major_radius/toroidal_field(i)
            d2_out(i) = d2_in(i)
            c0(i) = d0(i)*factor_of_phi
            c2_in(i) = d2_in(i)*factor_of_phi
            c2_out(i) = c2_in(i)
        end do

        include_phi = merge(1, 0, options%include_potential_in_current)
        call solve_flr2_tridiagonal(include_phi, background%r_grid, bpsi_over_bphi, &
                                    a2_in, a2_out, a0, b2_in, b2_out, b0, &
                                    c2_in, c2_out, c0, d2_in, d2_out, d0, &
                                    phi, parcur_over_b0, solver_stat)
        if (solver_stat /= FLR2_SOLVE_OK) then
            stat = solver_stat
            return
        end if
        stat = FLR2_RESPONSE_OK
    end subroutine solve_flr2_response

    logical function valid_background(background, n)
        type(plasma_t), intent(in) :: background
        integer, intent(in) :: n
        integer :: sp

        valid_background = .false.
        if (n < 3 .or. background%n_species < 1) return
        if (.not. allocated(background%spec)) return
        if (.not. allocated(background%r_grid) .or. size(background%r_grid) /= n) return
        if (.not. allocated(background%q) .or. size(background%q) /= n) return
        if (.not. allocated(background%B0) .or. size(background%B0) /= n) return
        if (.not. allocated(background%Er) .or. size(background%Er) /= n) return
        if (.not. allocated(background%om_E) .or. size(background%om_E) /= n) return
        if (.not. all(ieee_is_finite(background%r_grid))) return
        if (.not. all(ieee_is_finite(background%q))) return
        if (.not. all(ieee_is_finite(background%B0)) &
            .or. any(background%B0 <= 0.0_dp)) return
        if (.not. all(ieee_is_finite(background%Er))) return
        if (.not. all(ieee_is_finite(background%om_E))) return

        do sp = 0, background%n_species - 1
            if (.not. allocated(background%spec(sp)%n)) return
            if (.not. allocated(background%spec(sp)%vT)) return
            if (.not. allocated(background%spec(sp)%nu)) return
            if (.not. allocated(background%spec(sp)%rho_L)) return
            if (.not. allocated(background%spec(sp)%A1)) return
            if (.not. allocated(background%spec(sp)%A2)) return
            if (.not. allocated(background%spec(sp)%I11)) return
            if (.not. allocated(background%spec(sp)%I13)) return
            if (size(background%spec(sp)%n) /= n) return
            if (size(background%spec(sp)%vT) /= n) return
            if (size(background%spec(sp)%nu) /= n) return
            if (size(background%spec(sp)%rho_L) /= n) return
            if (size(background%spec(sp)%A1) /= n) return
            if (size(background%spec(sp)%A2) /= n) return
            if (size(background%spec(sp)%I11, 1) /= n) return
            if (size(background%spec(sp)%I13, 1) /= n) return
            if (lbound(background%spec(sp)%I11, 2) > 0 &
                .or. ubound(background%spec(sp)%I11, 2) < 0) return
            if (lbound(background%spec(sp)%I13, 2) > 0 &
                .or. ubound(background%spec(sp)%I13, 2) < 0) return
            if (.not. all(ieee_is_finite(background%spec(sp)%n)) &
                .or. any(background%spec(sp)%n <= 0.0_dp)) return
            if (.not. all(ieee_is_finite(background%spec(sp)%vT)) &
                .or. any(background%spec(sp)%vT <= 0.0_dp)) return
            if (.not. all(ieee_is_finite(background%spec(sp)%nu)) &
                .or. any(background%spec(sp)%nu <= 0.0_dp)) return
            if (.not. all(ieee_is_finite(background%spec(sp)%rho_L)) &
                .or. any(background%spec(sp)%rho_L < 0.0_dp)) return
            if (.not. all(ieee_is_finite(background%spec(sp)%A1))) return
            if (.not. all(ieee_is_finite(background%spec(sp)%A2))) return
            if (.not. finite_complex_2d(background%spec(sp)%I11)) return
            if (.not. finite_complex_2d(background%spec(sp)%I13)) return
        end do
        valid_background = .true.
    end function valid_background

    logical function finite_complex_1d(values)
        complex(dp), intent(in) :: values(:)

        finite_complex_1d = all(ieee_is_finite(real(values, dp))) &
                            .and. all(ieee_is_finite(aimag(values)))
    end function finite_complex_1d

    logical function finite_complex_2d(values)
        complex(dp), intent(in) :: values(:, :)

        finite_complex_2d = all(ieee_is_finite(real(values, dp))) &
                            .and. all(ieee_is_finite(aimag(values)))
    end function finite_complex_2d

    pure real(dp) function logical_factor(enabled)
        logical, intent(in) :: enabled

        logical_factor = merge(1.0_dp, 0.0_dp, enabled)
    end function logical_factor

end module flr2_response_m
