subroutine generate_grids

    use grid_m, only: rg_grid, xl_grid, reduced_rg_dim, l_space_dim, grid_spacing_rg, grid_spacing_xl, r_min, r_plas
    use species_m, only: plasma
    use config_m, only: fdebug

    implicit none

    call rg_grid%grid_init(reduced_rg_dim, 3.0d0, plasma%r_grid(size(plasma%r_grid)), 'rg')
    call xl_grid%grid_init(l_space_dim, 3.0d0, plasma%r_grid(size(plasma%r_grid)), 'xl')

    ! Generate rg_grid according to requested spacing
    select case (trim(adjustl(grid_spacing_rg)))
    case ("equidistant")
        call rg_grid%grid_generate_linear()
    case ("non-equidistant", "adaptive")
        call rg_grid%grid_generate()
    case default
        call rg_grid%grid_generate()
    end select

    ! Generate xl_grid according to requested spacing
    select case (trim(adjustl(grid_spacing_xl)))
    case ("equidistant")
        call xl_grid%grid_generate_linear()
    case ("non-equidistant", "adaptive")
        call xl_grid%grid_generate()
    case default
        call xl_grid%grid_generate()
    end select

    if (fdebug == 1) then
        write(*,*) " Generated Grid number of points:"
        write(*,*) ' Nrg = ', rg_grid%npts_b, ", Nl = ", xl_grid%npts_b
        write(*,*) ''
    end if

end subroutine
