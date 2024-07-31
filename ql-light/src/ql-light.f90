program ql_light

    use h5mod

    implicit none

    write(*,*) 'Hello, World!'

    ! read in profiles and other data (txt and/or hdf5)
    
    ! determine gradients
 
    ! call wave code to calculate EM field perturbations
    ! from wave_code_data_64bit (uses kilca to calculate wave fields, needs maybe linking of kilca)
        ! something like thsi
        !if (irf .eq. 1) call update_background_files(path2profs);
        !if (irf .eq. 1) call get_wave_code_data(imin, imax);
        !if (irf .eq. 1) call get_background_magnetic_fields_from_wave_code(flre_cd_ptr(imin), dim_r, r, B0t, B0z, B0);
        !if (irf .eq. 1) call get_collision_frequences_from_wave_code(flre_cd_ptr(imin), dim_r, r, nui, nue);


        !call get_wave_vectors_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
        !                                     m_vals(i_mn), n_vals(i_mn), ks, kp)
        !call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
        !                                    m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)

    ! calculate quasilinear diffusion coefficient Dqle22
    ! From get_dql routine in rhs_balance.f90

        !if (.true.) then
                !call calc_transport_coeffs_ornuhl(npoib, vT_e, nu_e, de11, de12, de21, de22)
                !call calc_transport_coeffs_ornuhl(npoib, vT_i, nu_i, di11, di12, di21, di22)
            !else
                !call calc_transport_coeffs_ornuhl_drift(1, npoib, de11, de12, de21, de22)
                !call calc_transport_coeffs_ornuhl_drift(2, npoib, di11, di12, di21, di22)
            !end if
        !end if

        ! localizer ?!?!??!
    
    ! interpolate Dqle22

    ! write out Dqle22 to file (txt and/or hdf5)



end program