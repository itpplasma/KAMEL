subroutine kim_prepare_resonances

    use kim_resonances_m
    use setup_m, only: m_mode, n_mode, type_br_field
    use species_m, only: plasma

    implicit none

    integer :: resonance_status

    iunit_res=157

    if (prescribed_r_res_active) then
        if (prescribed_r_res < minval(plasma%r_grid) .or. &
                prescribed_r_res > maxval(plasma%r_grid)) then
            error stop 'Prescribed resonance is outside the KIM profile grid'
        end if
        r_res = prescribed_r_res
    else
        call locate_periodic_resonance(plasma%r_grid, plasma%q, m_mode, n_mode, &
            r_res, resonance_status)
        if (resonance_status /= KIM_RESONANCE_OK) then
            write(*,*) "Resonance location not found in q"
            r_res = 0.0d0
            return
        end if
    end if

    if (.not. prescribed_r_res_active .and. type_br_field == 2) then
        r_res = plasma%r_grid(plasma%grid_size)/2
    end if

    write(*,*) 'resonant radius: ',r_res

end subroutine kim_prepare_resonances
