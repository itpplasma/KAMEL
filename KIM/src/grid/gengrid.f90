subroutine generate_grids

    use grid_m, only: rg_grid, xl_grid, reduced_rg_dim, l_space_dim, grid_spacing, r_min, r_plas
    use species_m, only: plasma
    use config_m, only: fdebug

    implicit none

    !call rg_grid%grid_init(reduced_rg_dim, plasma%r_grid(1), plasma%r_grid(size(plasma%r_grid)), 'rg')
    !call xl_grid%grid_init(l_space_dim, plasma%r_grid(1), plasma%r_grid(size(plasma%r_grid)), 'xl')
    call rg_grid%grid_init(reduced_rg_dim, 3.0d0, plasma%r_grid(size(plasma%r_grid)), 'rg')
    call xl_grid%grid_init(l_space_dim, 3.0d0, plasma%r_grid(size(plasma%r_grid)), 'xl')


    if (grid_spacing == 2) then
        write(*,*) "Generating linear grids"
        call rg_grid%grid_generate_linear()
        call xl_grid%grid_generate_linear()
    else if (grid_spacing == 3) then
        call rg_grid%grid_generate()
        call xl_grid%grid_generate()
    else
        call rg_grid%grid_generate()
        call xl_grid%grid_generate()
    end if

    if (fdebug == 1) then
        write(*,*) " Generated Grid number of points:"
        write(*,*) ' Nrg = ', rg_grid%npts_b, ", Nl = ", xl_grid%npts_b
        write(*,*) ''
    end if

end subroutine
