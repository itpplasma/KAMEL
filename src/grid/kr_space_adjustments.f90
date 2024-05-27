
subroutine kr_space_adjustments

    use kr_grid, only: kr, k_space_dim
    use setup, only: kr_cut_off_fac
    use plasma_parameter, only: rho_L
    use kr_grid, only: closest_kr_ind_upper, closest_kr_ind_lower
    use config, only: fstatus

    implicit none

    double precision :: kr_cutoff
    integer :: findClosestIndex

    ! determine cut-off in kr and corresponding indices
    kr_cutoff = kr_cut_off_fac / rho_L

    call generate_k_space_grid(k_space_dim, .true., kr_cutoff)

    closest_kr_ind_upper = findClosestIndex(kr, kr_cutoff)
    closest_kr_ind_lower = findClosestIndex(-kr, kr_cutoff)

    if (fstatus == 1) then
        write(*,*) ' rho_L      = ', rho_L
        write(*,*) ' kr cut-off = ', kr_cutoff
        write(*,*) ' closest index lower = ', closest_kr_ind_lower, ', closest index upper = ', closest_kr_ind_upper
        write(*,*) ' closest lower = ', kr(closest_kr_ind_lower), ', closest upper = ', kr(closest_kr_ind_upper)
    end if
end subroutine