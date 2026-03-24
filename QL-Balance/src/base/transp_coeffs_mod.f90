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

    subroutine compute_antenna_factor_from_Ipar
        !
        ! Compute antenna_factor = (I_par_toroidal / I_par_cylindrical)^2
        ! from the GPEC toroidal shielding current and the wave-code
        ! cylindrical shielding current integrated at runtime.
        !
        ! antenna_factor is the squared ratio of currents (= squared field ratio).
        ! It is applied once to D_ql in rescale_transp_coeffs_by_ant_fac,
        ! consistent with the legacy namelist convention.
        !
        ! Unit conversion:
        !   I_par_toroidal is in CGS with c=1 (Gaussian-like, from GPEC)
        !   Ipar from integrate_parallel_current is int(J*r*dr) in full CGS
        !     without the 2*pi factor, i.e. physical current = 2*pi*Ipar
        !   To convert I_par_toroidal to full CGS: multiply by c
        !   So: antenna_factor = (I_par_toroidal * c / (2*pi * |Ipar|))^2
        !
        use wave_code_data, only: antenna_factor, I_par_toroidal
        use grid_mod, only: Ipar

        implicit none

        double precision, parameter :: c_light = 2.99792458d10   ! speed of light [cm/s]
        double precision, parameter :: twopi = 6.283185307179586d0
        double precision :: I_tor_cgs, I_cyl_cgs

        if (I_par_toroidal <= 0.0d0) return  ! legacy mode: use namelist antenna_factor

        if (abs(Ipar) < 1.0d-30) then
            write(*,*) 'WARNING: I_par from wave code is ~0; keeping namelist antenna_factor'
            return
        end if

        ! Convert GPEC current (c=1 CGS) to full CGS
        I_tor_cgs = I_par_toroidal * c_light
        ! Physical cylindrical current = 2*pi * int(J*r*dr)
        I_cyl_cgs = twopi * abs(Ipar)

        antenna_factor = (I_tor_cgs / I_cyl_cgs)**2

        write(*,*) 'Antenna factor computed from I_par_toroidal:'
        write(*,*) '  I_par_toroidal (c=1 CGS) = ', I_par_toroidal
        write(*,*) '  I_par_toroidal (full CGS) = ', I_tor_cgs
        write(*,*) '  I_par_cylindrical (full CGS) = ', I_cyl_cgs
        write(*,*) '  I_par_cylindrical (c=1 CGS) = ', I_cyl_cgs / c_light
        write(*,*) '  antenna_factor  = ', antenna_factor

    end subroutine

end module
