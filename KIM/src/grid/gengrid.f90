subroutine generate_grids

    use grid, only: rg_grid, xl_grid, reduced_rg_dim, l_space_dim, grid_spacing
    use species, only: plasma

    implicit none

    call rg_grid%grid_init(reduced_rg_dim, plasma%r_grid(1), plasma%r_grid(size(plasma%r_grid)), 'rg')
    call xl_grid%grid_init(l_space_dim, plasma%r_grid(1), plasma%r_grid(size(plasma%r_grid)), 'xl')

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

end subroutine
