program test_dqli_limit_benchmark
    !! Algebraic limiting-case benchmark for the integral ion tensor.
    !! At ell=0, B_parallel=0, and k_perp rho_i -> 0 the gyrokinetic
    !! integral must reduce to the established Heyn/Markl tensor.
    use KIM_kinds_m, only: dp
    use quasilinear_integral_m, only: calc_ion_integral_harmonic
    implicit none

    complex(dp) :: fields(3), ifunc(0:3,0:3), es_field
    real(dp) :: integral_tensor(2,2), drift_tensor(2,2)
    real(dp) :: vT, nu, omega_c, c_light, B0, epm2, brm2, epbr_re, epbr_im, comfac
    real(dp) :: ks, limit_error, scale, d12a, d12base
    real(dp), parameter :: ks_values(2) = [0.0_dp, 1.0e-4_dp]
    integer :: ik

    vT = 2.3e7_dp
    nu = 3.7e5_dp
    omega_c = 4.1e7_dp
    c_light = 2.99792458e10_dp
    B0 = 1.8e4_dp
    fields = (0.0_dp, 0.0_dp)
    fields(2) = cmplx(0.8_dp, 0.3_dp, dp)
    ! The limiting benchmark is deliberately B_parallel-free.
    fields(3) = (0.0_dp, 0.0_dp)

    ifunc = (0.0_dp, 0.0_dp)
    ! Local diagonal radial-wave-number reduction retains the three moments
    ! entering the zero-FLR drift tensor; higher energy moments are O(k_perp^2)
    ! and vanish in this benchmark limit.
    ifunc(0,0) = cmplx(0.17_dp, 0.02_dp, dp)
    ifunc(1,0) = cmplx(0.23_dp, -0.01_dp, dp)
    ifunc(1,1) = cmplx(0.31_dp, 0.03_dp, dp)

    do ik = 1, size(ks_values)
        ks = ks_values(ik)
        if (ik == 1) then
            ! At exact k_perp=0 the electrostatic potential channel has
            ! vanishing E_perp; retain the magnetic channel for the regular
            ! zero-FLR branch.
            es_field = (0.0_dp, 0.0_dp)
            fields(1) = (0.0_dp, 0.0_dp)
        else
            es_field = cmplx(1.2e-3_dp, -0.7e-3_dp, dp)
            ! Hold -i c k_s Phi fixed while k_s approaches zero.
            fields(1) = cmplx(0.0_dp, 1.0_dp, dp)*es_field/(c_light*ks)
        end if
        comfac = 0.5_dp/(nu*B0**2)
        epm2 = abs(es_field)**2
        brm2 = vT**2*abs(fields(2))**2
        epbr_re = 2.0_dp*vT*real(conjg(es_field)*fields(2), dp)
        epbr_im = 2.0_dp*vT*aimag(conjg(es_field)*fields(2))
        drift_tensor(1,1) = comfac*(epm2*real(ifunc(0,0),dp) + epbr_re*real(ifunc(1,0),dp) &
            + brm2*real(ifunc(1,1),dp))
        d12base = comfac*(epm2*real(ifunc(0,0)+0.5_dp*ifunc(2,0),dp) &
            + epbr_re*real(ifunc(1,0)+0.25_dp*(ifunc(3,0)+ifunc(2,1)),dp) &
            + brm2*real(ifunc(1,1)+0.5_dp*ifunc(3,1),dp))
        d12a = comfac*epbr_im*0.25_dp*aimag(ifunc(2,1)-ifunc(3,0))
        drift_tensor(1,2) = d12base + d12a
        drift_tensor(2,1) = d12base - d12a
        drift_tensor(2,2) = comfac*(epm2*real(2.0_dp*ifunc(0,0)+ifunc(2,0)+0.25_dp*ifunc(2,2),dp) &
            + epbr_re*real(2.0_dp*ifunc(1,0)+0.5_dp*(ifunc(3,0)+ifunc(2,1)) &
                + 0.25_dp*ifunc(3,2),dp) &
            + brm2*real(2.0_dp*ifunc(1,1)+ifunc(3,1)+0.25_dp*ifunc(3,3),dp))
        call calc_ion_integral_harmonic(0, ks, 0.0_dp, ks, 0.0_dp, vT, omega_c, -omega_c, &
            c_light, B0, nu, fields, fields, ifunc, integral_tensor)
        scale = max(1.0_dp, maxval(abs(drift_tensor)))
        limit_error = maxval(abs(integral_tensor-drift_tensor))/scale
        if (ik == 1 .and. limit_error > 2.0e-12_dp) &
            then
            print *, 'integral tensor = ', integral_tensor
            print *, 'drift tensor    = ', drift_tensor
            print *, 'relative error  = ', limit_error
            error stop 'exact zero-FLR integral/drift benchmark mismatch'
        end if
        if (ik == 2 .and. limit_error > 2.0e-7_dp) &
            error stop 'finite small-FLR integral/drift convergence mismatch'
    end do

    print *, 'drift-kinetic limiting benchmark passed; max relative error = ', limit_error
end program test_dqli_limit_benchmark
