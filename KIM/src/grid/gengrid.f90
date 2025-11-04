subroutine generate_grids

    use grid_m, only: rg_grid, xl_grid, rg_space_dim, l_space_dim, grid_spacing_rg, grid_spacing_xl, &
        r_min, r_plas
    use species_m, only: plasma
    use config_m, only: fdebug
    use IO_collection_m, only: write_profile

    implicit none

    ! Generate rg_grid according to requested spacing
    select case (trim(adjustl(grid_spacing_rg)))
    case ("equidistant")
        call rg_grid%grid_init_equidistant(rg_space_dim, r_min, r_plas, 'rg')
        call rg_grid%grid_generate_equidistant()
    case ("non-equidistant", "adaptive")
        call rg_grid%grid_init(rg_space_dim, r_min, r_plas, 'rg')
        call rg_grid%grid_generate()
    case default
        call rg_grid%grid_init(rg_space_dim, r_min, r_plas, 'rg')
        call rg_grid%grid_generate()
    end select

    ! Generate xl_grid according to requested spacing
    select case (trim(adjustl(grid_spacing_xl)))
    case ("equidistant")
        call xl_grid%grid_init_equidistant(l_space_dim, r_min, r_plas, 'xl')
        call xl_grid%grid_generate_equidistant()
    case ("non-equidistant", "adaptive")
        call xl_grid%grid_init(l_space_dim, r_min, r_plas, 'xl')
        call xl_grid%grid_generate()
    case default
        call xl_grid%grid_init(l_space_dim, r_min, r_plas, 'xl')
        call xl_grid%grid_generate()
    end select

    if (fdebug == 1) then
        write(*,*) " Generated Grid number of points:"
        write(*,*) ' Nrg = ', rg_grid%npts_b, ", Nl = ", xl_grid%npts_b
        write(*,*) " rg grid h = ", rg_grid%xb(2) - rg_grid%xb(1)
        write(*,*) " xl grid h = ", xl_grid%xb(2) - xl_grid%xb(1)
        write(*,*) " xl / rg grid h ratio = ", (xl_grid%xb(2) - xl_grid%xb(1)) / (rg_grid%xb(2) - rg_grid%xb(1))
        write(*,*) ''
    end if

    call write_profile(xl_grid%xb, xl_grid%xb, xl_grid%npts_b, 'grid/'//trim(xl_grid%name)//'_xb', &
        'Cell boundary points of grid', 'cm')
    call write_profile(xl_grid%xc, xl_grid%xc, xl_grid%npts_c, 'grid/'//trim(xl_grid%name)//'_xc', &
        'Cell center points of grid', 'cm')

    call write_profile(rg_grid%xb, rg_grid%xb, rg_grid%npts_b, 'grid/'//trim(rg_grid%name)//'_xb', &
        'Cell boundary points of grid', 'cm')
    call write_profile(rg_grid%xc, rg_grid%xc, rg_grid%npts_c, 'grid/'//trim(rg_grid%name)//'_xc', &
        'Cell center points of grid', 'cm')

end subroutine
