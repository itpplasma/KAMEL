module kim_qldiff_m

    use KIM_kinds_m, only: dp

    implicit none

    private
    public :: calc_dqle22

contains

    function calc_dqle22(vTe, nue, om_E, B0, kpar, Es, Br) result(dqle22)
        ! Quasilinear electron heat diffusion coefficient D_ql,e22, the
        ! (2,2) component of the electron Onsager matrix, evaluated from
        ! the local wave fields and plasma parameters.
        !
        ! Faithful port of the electron (2,2) component of
        ! calc_transport_coeffs_ornuhl (QL-Balance/src/base/diff_coeffs.f90),
        ! using the same generalized plasma dispersion functions I^kl via
        ! getIfunc (QL-Balance/src/base/getIfunc.f90, compiled into KIM_lib
        ! through the KIM_fokkerplanck source group).
        !
        ! All quantities in CGS units:
        !   vTe  - electron thermal velocity [cm/s]
        !   nue  - electron collision frequency [1/s]
        !   om_E - ExB rotation frequency [rad/s]
        !   B0   - equilibrium magnetic field [G]
        !   kpar - parallel wave number [1/cm] (-> 0 at the resonant surface)
        !   Es   - complex perpendicular electric field amplitude [statV/cm]
        !   Br   - complex radial magnetic perturbation amplitude [G]
        use constants_m, only: sol

        real(dp), intent(in) :: vTe, nue, om_E, B0, kpar
        complex(dp), intent(in) :: Es, Br
        real(dp) :: dqle22

        real(dp) :: x1, x2, comfac, epm2, brm2, epbr_re
        complex(dp) :: symbI(0:3, 0:3)

        interface
            subroutine getIfunc(x1, x2, symbI)
                double precision, intent(in) :: x1, x2
                double complex, dimension(0:3, 0:3), intent(out) :: symbI
            end subroutine
        end interface

        ! Normalized distance to resonance and inverse normalized
        ! collisionality (static perturbation, om = 0)
        x1 = kpar * vTe / nue
        x2 = -om_E / nue

        call getIfunc(x1, x2, symbI)

        comfac = 0.5_dp / (nue * B0**2)
        epm2 = sol**2 * abs(Es)**2
        brm2 = vTe**2 * abs(Br)**2
        epbr_re = 2.0_dp * sol * vTe * real(conjg(Es) * Br, dp)

        dqle22 = comfac * (epm2 * real(2.0_dp * symbI(0, 0) + symbI(2, 0) &
                                       + 0.25_dp * symbI(2, 2), dp) &
                           + epbr_re * real(2.0_dp * symbI(1, 0) &
                                            + 0.5_dp * (symbI(3, 0) + symbI(2, 1)) &
                                            + 0.25_dp * symbI(3, 2), dp) &
                           + brm2 * real(2.0_dp * symbI(1, 1) + symbI(3, 1) &
                                         + 0.25_dp * symbI(3, 3), dp))

    end function

end module
