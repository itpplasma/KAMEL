module kim_qldiff_m

    use KIM_kinds_m, only: dp

    implicit none

    private
    public :: calc_dqle22, calc_dqli11_phi, calc_dqli_tensor
    public :: calc_dqli_integral_harmonic

contains

    subroutine calc_dqli_integral_harmonic(ell, ks_s, kr_s, ks_o, kr_o, &
            vTi, nui, omega_ci, omega_mode, om_E, B0, kpar, fields_s, &
            fields_o, tensor)
        !! Production entry point for one ion cyclotron harmonic and one
        !! ordered radial-wave pair.  Periodic KIM's spectral assembly sums
        !! this building block over ell, kr_s and kr_o; #290 consumes that
        !! completed local tensor without duplicating KIM's conventions.
        use constants_m, only: sol
        use config_m, only: resolved_ion_ifunc_conservation_model
        use quasilinear_integral_m, only: calc_ion_integral_harmonic
        use species_m, only: evaluate_susceptibility
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks_s, kr_s, ks_o, kr_o
        real(dp), intent(in) :: vTi, nui, omega_ci, omega_mode, om_E, B0, kpar
        complex(dp), intent(in) :: fields_s(3), fields_o(3)
        real(dp), intent(out) :: tensor(2,2)
        real(dp) :: x1, x2
        complex(dp) :: symbI(0:3,0:3)
        if (abs(omega_ci) <= tiny(1.0_dp)) &
            error stop 'calc_dqli_integral_harmonic: zero signed cyclotron frequency'
        x1 = kpar*vTi/nui
        x2 = -(om_E+real(ell,dp)*omega_ci-omega_mode)/nui
        call evaluate_susceptibility(x1, x2, &
            resolved_ion_ifunc_conservation_model, symbI)
        call calc_ion_integral_harmonic(ell, ks_s, kr_s, ks_o, kr_o, &
            vTi, abs(omega_ci), omega_ci, sol, B0, nui, fields_s, &
            fields_o, symbI, tensor)
    end subroutine calc_dqli_integral_harmonic

    function calc_dqli11_phi(vTi, nui, om_E, B0, kpar, Es) result(dqli11)
        ! Ion Phi-only integral coefficient (the D11 tracer bullet).  This is
        ! the direct integral-formalism analogue of the electrostatic term in
        ! the Onsager tensor; magnetic and cross terms are deliberately absent.
        use constants_m, only: sol
        use config_m, only: resolved_ion_ifunc_conservation_model
        real(dp), intent(in) :: vTi, nui, om_E, B0, kpar
        complex(dp), intent(in) :: Es
        real(dp) :: dqli11, x1, x2, comfac
        complex(dp) :: symbI(0:3,0:3)
        interface
            subroutine getIfunc_model(x1, x2, conservation_model, symbI)
                double precision, intent(in) :: x1, x2
                integer, intent(in) :: conservation_model
                double complex, dimension(0:3,0:3), intent(out) :: symbI
            end subroutine getIfunc_model
        end interface
        x1 = kpar*vTi/nui
        x2 = -om_E/nui
        call getIfunc_model(x1, x2, resolved_ion_ifunc_conservation_model, symbI)
        comfac = 0.5_dp/(nui*B0**2)
        dqli11 = comfac*sol**2*abs(Es)**2*real(symbI(0,0),dp)
    end function calc_dqli11_phi

    subroutine calc_dqli_tensor(vTi, nui, om_E, B0, kpar, Es, Br, D11, D12, D21, D22)
        ! Complete ion Phi/Br Onsager tensor in the integral formalism.
        ! This is the scalar, local counterpart of QL-Balance's
        ! calc_transport_coeffs_ornuhl; electrons intentionally retain their
        ! existing drift-kinetic path.
        use constants_m, only: sol
        use config_m, only: resolved_ion_ifunc_conservation_model
        real(dp), intent(in) :: vTi, nui, om_E, B0, kpar
        complex(dp), intent(in) :: Es, Br
        real(dp), intent(out) :: D11, D12, D21, D22
        real(dp) :: x1, x2, comfac, epm2, brm2, epbr_re, epbr_im, d12a
        complex(dp) :: symbI(0:3,0:3)
        interface
            subroutine getIfunc_model(x1, x2, conservation_model, symbI)
                double precision, intent(in) :: x1, x2
                integer, intent(in) :: conservation_model
                double complex, dimension(0:3,0:3), intent(out) :: symbI
            end subroutine getIfunc_model
        end interface

        x1 = kpar*vTi/nui
        x2 = -om_E/nui
        call getIfunc_model(x1, x2, resolved_ion_ifunc_conservation_model, symbI)

        comfac = 0.5_dp/(nui*B0**2)
        epm2 = sol**2*abs(Es)**2
        brm2 = vTi**2*abs(Br)**2
        epbr_re = 2.0_dp*sol*vTi*real(conjg(Es)*Br,dp)
        epbr_im = 2.0_dp*sol*vTi*aimag(conjg(Es)*Br)

        D11 = comfac*(epm2*real(symbI(0,0),dp) + epbr_re*real(symbI(1,0),dp) &
             + brm2*real(symbI(1,1),dp))
        D12 = comfac*(epm2*real(symbI(0,0)+0.5_dp*symbI(2,0),dp) &
             + epbr_re*real(symbI(1,0)+0.25_dp*(symbI(3,0)+symbI(2,1)),dp) &
             + brm2*real(symbI(1,1)+0.5_dp*symbI(3,1),dp))
        D21 = D12
        D22 = comfac*(epm2*real(2.0_dp*symbI(0,0)+symbI(2,0)+0.25_dp*symbI(2,2),dp) &
             + epbr_re*real(2.0_dp*symbI(1,0)+0.5_dp*(symbI(3,0)+symbI(2,1)) &
                              +0.25_dp*symbI(3,2),dp) &
             + brm2*real(2.0_dp*symbI(1,1)+symbI(3,1)+0.25_dp*symbI(3,3),dp))

        d12a = comfac*epbr_im*0.25_dp*aimag(symbI(2,1)-symbI(3,0))
        D12 = D12 + d12a
        D21 = D21 - d12a
    end subroutine calc_dqli_tensor

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
        use config_m, only: resolved_electron_ifunc_conservation_model

        real(dp), intent(in) :: vTe, nue, om_E, B0, kpar
        complex(dp), intent(in) :: Es, Br
        real(dp) :: dqle22

        real(dp) :: x1, x2, comfac, epm2, brm2, epbr_re
        complex(dp) :: symbI(0:3, 0:3)

        interface
            subroutine getIfunc_model(x1, x2, conservation_model, symbI)
                double precision, intent(in) :: x1, x2
                integer, intent(in) :: conservation_model
                double complex, dimension(0:3, 0:3), intent(out) :: symbI
            end subroutine getIfunc_model
        end interface

        ! Normalized distance to resonance and inverse normalized
        ! collisionality (static perturbation, om = 0)
        x1 = kpar * vTe / nue
        x2 = -om_E / nue

        call getIfunc_model(x1, x2, &
            resolved_electron_ifunc_conservation_model, symbI)

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
