module transp_coeffs_mod

    implicit none

    contains

    subroutine rescale_transp_coeffs_by_ant_fac

        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
                            dqli11, dqli12, dqli21, dqli22
        use wave_code_data, only: antenna_factor

        implicit none

        dqle11 = dqle11*antenna_factor
        dqle12 = dqle12*antenna_factor
        dqle21 = dqle21*antenna_factor
        dqle22 = dqle22*antenna_factor
        dqli11 = dqli11*antenna_factor
        dqli12 = dqli12*antenna_factor
        dqli21 = dqli21*antenna_factor
        dqli22 = dqli22*antenna_factor

    end subroutine

end module
