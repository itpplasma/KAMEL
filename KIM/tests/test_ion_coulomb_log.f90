program test_ion_coulomb_log
    ! Ion-ion Coulomb logarithm against the NRL Plasma Formulary (2019, p. 34):
    !
    !   lambda_ii' = 23 - ln[ Z Z' (mu + mu') / (mu T_i' + mu' T_i)
    !                         * (n_i Z^2 / T_i + n_i' Z'^2 / T_i')^(1/2) ]
    !
    ! with n in cm^-3, T in eV and mu the mass number. For like ions it reduces
    ! to 23 - ln[ Z^3 T^(-3/2) (2 n)^(1/2) ], the form KiLCA uses in
    ! background_data_m.f90, which serves as an independent check.

    use KIM_kinds_m, only: dp
    use species_m, only: ion_ion_coulomb_log

    implicit none

    call test_like_ions_match_kilca_form()
    call test_like_ions_reference_value()
    call test_mixed_ions_reference_value()
    call test_symmetric_in_species()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine require_close(got, want, tol, label)
        real(dp), intent(in) :: got, want, tol
        character(*), intent(in) :: label

        if (abs(got - want) > tol * abs(want)) then
            print *, 'FAIL: ', label
            print *, '  got  = ', got
            print *, '  want = ', want
            error stop 'ion-ion Coulomb logarithm has the wrong value'
        end if
        print *, 'PASS: ', label
    end subroutine require_close

    subroutine test_like_ions_match_kilca_form()
        integer, parameter :: Z = 1, A = 2
        real(dp), parameter :: n_values(3) = [1.0d12, 1.0d13, 5.0d14]
        real(dp), parameter :: T_values(3) = [5.0d1, 1.0d3, 1.0d4]
        real(dp) :: kilca_form
        integer :: j, k

        do j = 1, size(n_values)
            do k = 1, size(T_values)
                kilca_form = 23.0d0 - log(real(Z, dp)**3 / T_values(k)**1.5d0 &
                                          * sqrt(2.0d0 * n_values(j)))
                call require_close(ion_ion_coulomb_log(Z, A, n_values(j), T_values(k), &
                                                       Z, A, n_values(j), T_values(k)), &
                                   kilca_form, 1.0d-13, 'like ions match the KiLCA form')
            end do
        end do
    end subroutine test_like_ions_match_kilca_form

    subroutine test_like_ions_reference_value()
        ! Deuterium, n = 1e13 cm^-3, T = 1 keV: 23 - ln(sqrt(2e13) / 1000^1.5).
        call require_close(ion_ion_coulomb_log(1, 2, 1.0d13, 1.0d3, 1, 2, 1.0d13, 1.0d3), &
                           18.048256d0, 1.0d-7, 'D-D at 1e13 cm^-3 and 1 keV')
    end subroutine test_like_ions_reference_value

    subroutine test_mixed_ions_reference_value()
        ! D (1e13 cm^-3, 1 keV) with C6+ (1e11 cm^-3, 800 eV):
        ! 23 - ln[ 6 * 14 / (2 * 800 + 12 * 1000) * sqrt(1e10 + 4.5e9) ].
        call require_close(ion_ion_coulomb_log(1, 2, 1.0d13, 1.0d3, 6, 12, 1.0d11, 8.0d2), &
                           16.388301d0, 1.0d-7, 'D-C6+ mixed ions')
    end subroutine test_mixed_ions_reference_value

    subroutine test_symmetric_in_species()
        real(dp) :: forward, backward

        forward = ion_ion_coulomb_log(1, 2, 1.0d13, 1.0d3, 6, 12, 1.0d11, 8.0d2)
        backward = ion_ion_coulomb_log(6, 12, 1.0d11, 8.0d2, 1, 2, 1.0d13, 1.0d3)
        call require_close(backward, forward, 1.0d-14, 'lambda_ii'' = lambda_i''i')
    end subroutine test_symmetric_in_species

end program test_ion_coulomb_log
