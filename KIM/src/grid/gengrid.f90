subroutine generate_grids

    use grid
    use config, only: fdebug, output_path
    use plasma_parameter, only: r_prof, iprof_length
    use setup, only: kr_cut_off_fac

    implicit none

    call rg_grid%grid_init(reduced_rg_dim, r_prof(1), r_prof(size(r_prof)), 'rg')
    call xl_grid%grid_init(l_space_dim, r_prof(1), r_prof(size(r_prof)), 'xl')

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
