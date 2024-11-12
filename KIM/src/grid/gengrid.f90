subroutine generate_grids

    use grid
    use config, only: fdebug, output_path
    use plasma_parameter, only: r_prof, iprof_length
    use setup, only: kr_cut_off_fac

    implicit none

    call rg_grid%grid_init(reduced_rg_dim, r_prof(1), r_prof(iprof_length), 'rg')
    call xl_grid%grid_init(l_space_dim, r_prof(1), r_prof(iprof_length), 'xl')
    call kr_grid%grid_init(k_space_dim, -kr_cut_off_fac, kr_cut_off_fac, 'kr')
    call krp_grid%grid_init(k_space_dim, -kr_cut_off_fac, kr_cut_off_fac, 'krp')

    if (grid_spacing == 2) then
        write(*,*) "Generating linear grids"
        call rg_grid%grid_generate_linear()
        call xl_grid%grid_generate_linear()
        call kr_grid%grid_generate_linear()
        call krp_grid%grid_generate_linear()
    else if (grid_spacing == 3) then
        call rg_grid%grid_generate()
        call xl_grid%grid_generate()
        call kr_grid%grid_generate_integer()
        call krp_grid%grid_generate_integer()
    else
        call rg_grid%grid_generate()
        call xl_grid%grid_generate()
        call kr_grid%grid_generate()
        call krp_grid%grid_generate()
    end if

end subroutine

!!
  !subroutine gengrid(npoimin, write_out)
    !! Generates the grid for two types of boundary conditions at the outer
    !! boundary: iboutype=1 - fixed parameters, iboutype=2 - fixed fluxes
    !! At the inner boundary fixed fluxes = 0 are always assumed

    !use grid
    !use config, only: fdebug, output_path
    !use plasma_parameter, only: r_prof
    !use resonances_mod, only: index_rg_res, r_res, width_res, ampl_res

    !implicit none

    !integer, intent(in) :: npoimin
    !logical, intent(in) :: write_out
    !integer :: ipoib, nder, ipb, ipe
    !double precision :: hrmax, r, rnext, recnsp, rscale
    !!double precision, dimension(:),   allocatable :: x
    !double precision, dimension(:,:), allocatable :: coef

    !if (fdebug == 1) write(*,*) "Debug: entering gengrid"
    !!nbaleqs=4

    !nder=1
    !npoi_der=4
    !allocate(coef(0:nder,npoi_der))
    !!allocate(x(npoi_der))

    !! set parameters for grid spacing. grid_spacing=1: quidistant grid, grid_spacing=2: non-equidistant grid
    !if (grid_spacing == 1) then
        !width_res = 1.0
        !ampl_res = 0.0
    !elseif (grid_spacing == 2) then
        !width_res = 3.0
        !ampl_res = 0.3
    !end if

    !rmin = minval(r_prof)
    !rmax = maxval(r_prof)

    !hrmax = (rmax - rmin)/(npoimin+1)

    !npoib = 1
    !r = rmin

    !do while(r .lt. rmax)
        !call recnsplit(r,recnsp)
        !rnext = r + hrmax / recnsp
        !call recnsplit(rnext,recnsp)
        !r = 0.5d0 * (rnext + r + hrmax / recnsp)
        !npoib= npoib + 1
    !enddo

    !npoic = npoib - 1

    !allocate(rb(npoib), rc(npoic))
    !allocate(Sb(npoib), Sc(npoic))

    !allocate(xl(npoib)) 

    !r = rmin
    !rb(1) = r

    !do ipoib=2,npoib
        !call recnsplit(r,recnsp)
        !rnext = r + hrmax / recnsp
        !call recnsplit(rnext,recnsp)
        !r = 0.5d0 * (rnext + r + hrmax / recnsp)
        !rb(ipoib) = r
        !rc(ipoib-1) = 0.5 * (rb(ipoib-1) + rb(ipoib))
    !enddo

    !write(*,*) " - - - r grid: - - - "
    !write(*,*) "    h = ", rb(2) - rb(1)
    !write(*,*) '    Number points r (l) grid: ', npoib
    !write(*,*) " - - - - - - - - - - "

    !! modifie this if xl and rg should have different grids
    !xl = rb + 0.3
    !l_space_dim = npoib

    !! get index for resonant radius
    !call binsrc(abs(xl), 1, npoib, abs(r_res), index_rg_res)

    !if(iboutype .eq. 1) then
        !rscale = (rmax - rmin) / (rc(npoic) - rmin)
    !else
        !rscale = (rmax - rmin) / (rb(npoib) - rmin)
    !endif
  !!rb=rmin+rscale*(rb-rmin)
  !!rc=rmin+rscale*(rc-rmin)

    !if(npoi_der .gt. npoic) then
        !write(*,*) '! Error : not enough grid points for derivatives'
        !stop
    !endif

  !allocate(deriv_coef(npoi_der, npoib), ipbeg(npoib), ipend(npoib))
  !allocate(reint_coef(npoi_der, npoib))

    !do ipoib = 1, npoib
        !ipb = ipoib - npoi_der / 2
        !ipe = ipb + npoi_der - 1
        !if(ipb .lt. 1) then
            !ipb = 1
            !ipe = ipb + npoi_der - 1
        !elseif(ipe .gt. npoic) then
          !ipe = npoic
          !ipb = ipe - npoi_der + 1
        !endif
        !ipbeg(ipoib) = ipb
        !ipend(ipoib) = ipe
        !call plag_coeff(npoi_der, nder, rb(ipoib), rc(ipb:ipe), coef)
        !deriv_coef(:, ipoib) = coef(1,:)
        !reint_coef(:, ipoib) = coef(0,:)
    !enddo

    !deallocate(coef)

    !if (write_out) call write_new_grid

    !if (fdebug == 1) write(*,*) "Debug: exiting gengrid"


    !contains

    !subroutine write_new_grid

        !implicit none
        !integer :: i
        !logical :: ex

        !inquire(file = trim(output_path)//'grid', exist = ex)
        !if (.not. ex) then
            !call system('mkdir -p '//trim(output_path)//'grid')
        !end if

        !open(unit = 77, file=trim(output_path)//'grid/rb.dat')
        !open(unit = 78, file=trim(output_path)//'grid/xl.dat')
        !do i = 1, npoib
            !write(77,*) i, rb(i)
            !write(78,*) i, xl(i)
        !end do
        !close(77)
        !close(78)

        !open(unit = 78, file=trim(output_path)//'grid/rc.dat')
        !do i = 1, npoic
            !write(78,*) i, rc(i)
        !end do
        !close(78)


    !end subroutine

!end subroutine gengrid

!subroutine generate_rg_grid(write_out)
    !! Generates the grid for two types of boundary conditions at the outer
    !! boundary: iboutype=1 - fixed parameters, iboutype=2 - fixed fluxes
    !! At the inner boundary fixed fluxes = 0 are always assumed

    !use grid, rg_b, rg_c, number_points_rg_b, number_points_rg_c, reduced_rg_dim
    !use config, only: fdebug, output_path
    !use plasma_parameter, only: r_prof
    !use resonances_mod, only: index_rg_res, r_res, width_res, ampl_res

    !implicit none

    !logical, intent(in) :: write_out
    !integer :: ipoib, nder, ipb, ipe
    !double precision :: hrmax, r, rnext, recnsp, rscale
    !double precision, dimension(:,:), allocatable :: coef

    !if (fdebug == 1) write(*,*) "Debug: entering generate_rg_grid"

    !nder=1
    !npoi_der=4
    !allocate(coef(0:nder,npoi_der))

    !! set parameters for grid spacing. grid_spacing=1: quidistant grid, grid_spacing=2: non-equidistant grid
    !if (grid_spacing == 1) then
        !width_res = 1.0
        !ampl_res = 0.0
    !elseif (grid_spacing == 2) then
        !width_res = 3.0
        !ampl_res = 0.3
    !end if

    !rmin = minval(r_prof)
    !rmax = maxval(r_prof)

    !hrmax = (rmax - rmin) / (reduced_rg_dim + 1)

    !number_points_rg_b = 1
    !r = rmin

    !do while(r .lt. rmax)
        !call recnsplit(r,recnsp)
        !rnext = r + hrmax / recnsp
        !call recnsplit(rnext,recnsp)
        !r = 0.5d0 * (rnext + r + hrmax / recnsp)
        !number_points_rg_b = number_points_rg_b + 1
    !enddo

    !number_points_rg_c = number_points_rg_b - 1

    !allocate(rg_b(number_points_rg_b), rg_c(number_points_rg_c))
    !!allocate(S_b(number_points_rg_b), S_c(npoic))

    !r = rmin
    !rg_b(1) = r

    !do ipoib=2, number_points_rg_b
        !call recnsplit(r,recnsp)
        !rnext = r + hrmax / recnsp
        !call recnsplit(rnext,recnsp)
        !r = 0.5d0 * (rnext + r + hrmax / recnsp)
        !rg_b(ipoib) = r
        !rg_c(ipoib-1) = 0.5 * (rg_b(ipoib-1) + rg_b(ipoib))
    !enddo

    !write(*,*) " - - - r grid: - - - "
    !write(*,*) "    h = ", rb(2) - rb(1)
    !write(*,*) '    Number points r grid (boundary): ', number_points_rg_b
    !write(*,*) " - - - - - - - - - - "

    !! get index for resonant radius
    !call binsrc(abs(rg_b), 1, number_points_rg_b, abs(r_res), index_rg_res)

    !if(npoi_der .gt. number_points_rg_c) then
        !write(*,*) '! Error : not enough grid points for derivatives'
        !stop
    !endif


    !! determine coefficients for derivatives
    !allocate(deriv_coef(npoi_der, number_points_rg_b), ipbeg(number_points_rg_b), ipend(number_points_rg_b))
    !allocate(reint_coef(npoi_der, number_points_rg_b))

    !do ipoib = 1, number_points_rg_b

        !ipb = ipoib - npoi_der / 2
        !ipe = ipb + npoi_der - 1

        !if(ipb .lt. 1) then
            !ipb = 1
            !ipe = ipb + npoi_der - 1
        !elseif(ipe .gt. npoic) then
            !ipe = number_points_rg_c
            !ipb = ipe - npoi_der + 1
        !endif

        !ipbeg(ipoib) = ipb
        !ipend(ipoib) = ipe

        !call plag_coeff(npoi_der, nder, rg_b(ipoib), rg_c(ipb:ipe), coef)
        !deriv_coef(:, ipoib) = coef(1,:)
        !reint_coef(:, ipoib) = coef(0,:)

    !enddo

    !deallocate(coef)

    !if (write_out) call write_new_grid

    !if (fdebug == 1) write(*,*) "Debug: exiting generate_rg_grid"

    !contains

    !subroutine write_new_grid

        !implicit none
        !integer :: i
        !logical :: ex

        !inquire(file = trim(output_path)//'grid', exist = ex)
        !if (.not. ex) then
            !call system('mkdir -p '//trim(output_path)//'grid')
        !end if

        !open(unit = 77, file=trim(output_path)//'grid/rg_b.dat')
        !open(unit = 78, file=trim(output_path)//'grid/rg_c.dat')
        !do i = 1, npoib
            !write(77,*) i, rg_b(i)
            !write(78,*) i, rg_c(i)
        !end do
        !close(77)
        !close(78)

    !end subroutine

!end subroutine generate_rg_grid

