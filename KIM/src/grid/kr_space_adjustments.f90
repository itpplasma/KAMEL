
subroutine kr_space_adjustments

    use kr_grid, only: kr, k_space_dim
    use setup_m, only: kr_cut_off_fac
    use config_m, only: fstatus
    use findIndex_m, only: findClosestIndex
    use KIM_kinds_m, only: dp

    implicit none

    real(dp) :: kr_cutoff

    ! determine cut-off in kr and corresponding indices
    kr_cutoff = kr_cut_off_fac !* rho_L

    call generate_k_space_grid(k_space_dim, .true., kr_cutoff)

    !call findClosestIndex(kr, kr_cutoff, closest_kr_ind_upper)
    !call findClosestIndex(-kr, kr_cutoff, closest_kr_ind_lower)

    if (fstatus == 1) then
        write(*,*) '- - - kr space adjustment - - -'
        write(*,*) ' rho_L      = ', rho_L, ' cm'
        write(*,*) ' kr cut-off = ', kr_cutoff
        write(*,*) ' min kr = ', kr(1), ', max kr = ', kr(k_space_dim)
        write(*,*) ' minval kr = ', minval(kr), ', maxval kr = ', maxval(kr)
        write(*,*) '- - - - - - - - - - - - - - - -'
    end if

end subroutine