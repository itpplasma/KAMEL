module grid

    use KIM_kinds, only: dp

    implicit none

    integer :: l_space_dim ! dimension of spline grid
    integer :: r_space_dim ! dimension of r grid
    integer :: k_space_dim ! dimension of kr space grid
    logical :: reduce_r
    integer :: reduced_rg_dim
    integer :: spline_base
    integer :: grid_spacing
    integer :: num_gengrid_points

    integer :: nder=1
    integer :: npoi_der=4

    real(dp), dimension(:), allocatable :: xl  ! xl grid (real space)

    complex(dp), dimension(:,:), allocatable :: varphi_lkr

!    integer :: number_points_rg_b, number_points_rg_c
    real(dp) :: gg_factor = 1.0
    real(dp) :: gg_width = 0.0
    real(dp) :: gg_r_res = 0.0!95.34
    integer, dimension(:),   allocatable :: ipbeg, ipend

    ! k-space specific
    real(dp) :: kr_grid_ampl_res
    real(dp) :: kr_grid_width_res
    real(dp) :: kr_res = 0.0d0
    
    type grid_type
        integer :: npts_b, npts_c, npts
        real(dp) :: min_val
        real(dp) :: max_val
        real(dp), dimension(:), allocatable :: xb
        real(dp), dimension(:), allocatable :: xc
        real(dp), dimension(:,:), allocatable :: deriv_coef
        real(dp), dimension(:,:), allocatable :: reint_coef
        character(len=:), allocatable :: name
        contains
            procedure :: grid_init
            procedure :: grid_generate
            procedure :: grid_generate_linear
            procedure :: grid_generate_integer
    end type grid_type

    type(grid_type) :: rg_grid, xl_grid, kr_grid, krp_grid

    contains

    subroutine grid_init(this, npts, min_val, max_val, name)

        use resonances_mod, only: width_res, ampl_res
        implicit none

        class(grid_type), intent(inout) :: this

        integer, intent(inout) :: npts
        real(dp), intent(in) :: min_val, max_val
        character(len=*), intent(in) :: name

        real(dp) :: hrmax
        real(dp) :: x_current, x_next
        real(dp) :: recnsp

        this%npts = npts
        this%npts_b = npts
        this%min_val = min_val
        this%max_val = max_val
        allocate(character(len=len(name)) :: this%name)
        this%name = name

        ! set parameters for grid spacing. grid_spacing=1: quidistant grid, grid_spacing=2: non-equidistant grid
        if (grid_spacing == 1) then
            width_res = 1.0
            ampl_res = 0.0
        elseif (grid_spacing == 2) then
            width_res = 3.0
            ampl_res = 0.3
        end if

        hrmax = (this%max_val - this%min_val) / (this%npts_b)

        this%npts_b = 1
        x_current = this%min_val

        do while(x_current .lt. this%max_val)
            call recnsplit(x_current, recnsp)
            x_next = x_current + hrmax / recnsp
            call recnsplit(x_next, recnsp)
            x_current = 0.5d0 * (x_next + x_current + hrmax / recnsp)
            this%npts_b = this%npts_b + 1
        enddo

        this%npts_c = this%npts_b - 1

    end subroutine

    subroutine grid_generate(this)

        use resonances_mod, only: r_res, index_rg_res
        use config, only: fdebug, output_path

        implicit none

        class(grid_type), intent(inout) :: this

        real(dp) :: x_current, x_next
        real(dp) :: hrmax
        integer :: ipoib, ipb, ipe
        real(dp), dimension(:,:), allocatable :: coef
        real(dp) :: recnsp
        
        allocate(this%xb(this%npts_b), this%xc(this%npts_c))
        allocate(coef(0:nder,npoi_der))

        hrmax = (this%max_val - this%min_val) / (this%npts_b + 1)

        x_current = this%min_val
        this%xb(1) = x_current

        do ipoib=2, this%npts_b
            call recnsplit(x_current, recnsp)
            x_next = x_current + hrmax / recnsp
            call recnsplit(x_next, recnsp)
            x_current = 0.5d0 * (x_next + x_current + hrmax / recnsp)
            this%xb(ipoib) = x_current
            this%xc(ipoib-1) = 0.5 * (this%xb(ipoib-1) + this%xb(ipoib))
        enddo

        write(*,*) " - - - grid ", this%name, ": - - - "
        write(*,*) "    h = ", this%xb(2) - this%xb(1)
        write(*,*) '    Number points r (l) grid: ', this%npts_b
        write(*,*) " - - - - - - - - - - "

        ! get index for resonant radius
        call binsrc(abs(this%xb), 1, this%npts_b, abs(r_res), index_rg_res)

        if(npoi_der .gt. this%npts_c) then
            write(*,*) '! Error : not enough grid points for derivatives'
            stop
        endif

        if (.not. allocated(ipbeg)) allocate(ipbeg(this%npts_b), ipend(this%npts_b))
        allocate(this%deriv_coef(npoi_der, this%npts_b))
        allocate(this%reint_coef(npoi_der, this%npts_b))

        do ipoib = 1, this%npts_c
            ipb = ipoib - npoi_der / 2
            ipe = ipb + npoi_der - 1
            if(ipb .lt. 1) then
                ipb = 1
                ipe = ipb + npoi_der - 1
            elseif(ipe .gt. this%npts_c) then
                ipe = this%npts_c
                ipb = ipe - npoi_der + 1
            endif
            ipbeg(ipoib) = ipb
            ipend(ipoib) = ipe
            call plag_coeff(npoi_der, nder, this%xb(ipoib), this%xc(ipb:ipe), coef)

            this%deriv_coef(:, ipoib) = coef(1,:)
            this%reint_coef(:, ipoib) = coef(0,:)

        enddo

        deallocate(coef)

        call write_new_grid

        if (fdebug == 1) write(*,*) "Debug: exiting gengrid"

        contains

            subroutine write_new_grid

                implicit none
                integer :: i
                logical :: ex

                inquire(file = trim(output_path)//'grid', exist = ex)
                if (.not. ex) then
                    call system('mkdir -p '//trim(output_path)//'grid')
                end if
                
                open(unit = 77, file=trim(output_path)//'grid/'//trim(this%name)//'_xb.dat')
                open(unit = 78, file=trim(output_path)//'grid/'//trim(this%name)//'_xc.dat')
                do i = 1, this%npts_b
                    write(77,*) i, this%xb(i)
                    write(78,*) i, this%xc(i)
                end do
                close(77)
                close(78)

            end subroutine

    end subroutine grid_generate


    subroutine grid_generate_linear(this)

        use resonances_mod, only: r_res, index_rg_res
        use config, only: fdebug, output_path

        implicit none

        class(grid_type), intent(inout) :: this

        real(dp) :: h
        integer :: ipoib, ipb, ipe
        real(dp), dimension(:,:), allocatable :: coef

        allocate(this%xb(this%npts_b), this%xc(this%npts_c))

        h = (this%max_val - this%min_val) / (this%npts_b-1)

        this%xb(1) = this%min_val
        do ipoib=2, this%npts_b
            this%xb(ipoib) = this%min_val + (ipoib - 1) * h
            this%xc(ipoib-1) = 0.5 * (this%xb(ipoib-1) + this%xb(ipoib))
        end do
        
        allocate(coef(0:nder,npoi_der))

        write(*,*) " - - - grid ", this%name, ": - - - "
        write(*,*) "    h = ", this%xb(2) - this%xb(1)
        write(*,*) "    max = ", this%max_val
        write(*,*) '    Number points r (l) grid: ', this%npts_b
        write(*,*) "    generating linear grid..."

        ! get index for resonant radius
        call binsrc(abs(this%xb), 1, this%npts_b, abs(r_res), index_rg_res)

        if(npoi_der .gt. this%npts_c) then
            write(*,*) '! Error : not enough grid points for derivatives'
            stop
        endif

        if (.not. allocated(ipbeg)) allocate(ipbeg(this%npts_c), ipend(this%npts_c))
        allocate(this%deriv_coef(npoi_der, this%npts_c))
        allocate(this%reint_coef(npoi_der, this%npts_c))

        do ipoib = 1, this%npts_c
            ipb = ipoib - npoi_der / 2
            ipe = ipb + npoi_der - 1
            if(ipb .lt. 1) then
                ipb = 1
                ipe = ipb + npoi_der - 1
            elseif(ipe .gt. this%npts_c) then
                ipe = this%npts_c
                ipb = ipe - npoi_der + 1
            endif
            ipbeg(ipoib) = ipb
            ipend(ipoib) = ipe
            
            call plag_coeff(npoi_der, nder, this%xb(ipoib), this%xc(ipb:ipe), coef)

            this%deriv_coef(:, ipoib) = coef(1,:)
            this%reint_coef(:, ipoib) = coef(0,:)

        enddo

        deallocate(coef, ipbeg, ipend)

        call write_new_grid

        if (fdebug == 1) write(*,*) "Debug: exiting gengrid"

        contains

            subroutine write_new_grid

                implicit none
                integer :: i
                logical :: ex

                inquire(file = trim(output_path)//'grid', exist = ex)
                if (.not. ex) then
                    call system('mkdir -p '//trim(output_path)//'grid')
                end if
                
                open(unit = 77, file=trim(output_path)//'grid/'//trim(this%name)//'_xb.dat')
                open(unit = 78, file=trim(output_path)//'grid/'//trim(this%name)//'_xc.dat')
                do i = 1, this%npts_c
                    write(77,*) i, this%xb(i)
                    write(78,*) i, this%xc(i)
                end do
                close(77)
                close(78)

            end subroutine

    end subroutine grid_generate_linear


    subroutine grid_generate_integer(this)

        use resonances_mod, only: r_res, index_rg_res
        use config, only: fdebug, output_path

        implicit none

        class(grid_type), intent(inout) :: this

        integer :: ipoib, ipb, ipe
        real(dp), dimension(:,:), allocatable :: coef

        this%npts = this%npts +1 
        this%npts_b = this%npts
        this%npts_c = this%npts

        allocate(this%xb(this%npts), this%xc(this%npts -1))

        this%xb(1) = -(this%npts-1)/2
        do ipoib=2, this%npts
            this%xb(ipoib) = this%xb(1) + ipoib - 1 
            this%xc(ipoib-1) = 0.5 * (this%xb(ipoib-1) + this%xb(ipoib))
        end do
        
        allocate(coef(0:nder,npoi_der))

        this%npts_b = this%npts

        write(*,*) " - - - grid ", this%name, ": - - - "
        write(*,*) "    h = ", this%xb(2) - this%xb(1)
        write(*,*) '    Number points r (l) grid: ', this%npts_b
        write(*,*) " - - - - - - - - - - "

        ! get index for resonant radius
        call binsrc(abs(this%xb), 1, this%npts_b, abs(r_res), index_rg_res)

        if(npoi_der .gt. this%npts_c) then
            write(*,*) '! Error : not enough grid points for derivatives'
            stop
        endif

        if (.not. allocated(ipbeg)) allocate(ipbeg(this%npts_b), ipend(this%npts_b))
        allocate(this%deriv_coef(npoi_der, this%npts_b))
        allocate(this%reint_coef(npoi_der, this%npts_b))

        do ipoib = 1, this%npts_b
            ipb = ipoib - npoi_der / 2
            ipe = ipb + npoi_der - 1
            if(ipb .lt. 1) then
                ipb = 1
                ipe = ipb + npoi_der - 1
            elseif(ipe .gt. this%npts_c) then
                ipe = this%npts_c
                ipb = ipe - npoi_der + 1
            endif
            ipbeg(ipoib) = ipb
            ipend(ipoib) = ipe
            call plag_coeff(npoi_der, nder, this%xb(ipoib), this%xc(ipb:ipe), coef)

            this%deriv_coef(:, ipoib) = coef(1,:)
            this%reint_coef(:, ipoib) = coef(0,:)

        enddo

        deallocate(coef)

        call write_new_grid

        if (fdebug == 1) write(*,*) "Debug: exiting gengrid"

        contains

            subroutine write_new_grid

                implicit none
                integer :: i
                logical :: ex

                inquire(file = trim(output_path)//'grid', exist = ex)
                if (.not. ex) then
                    call system('mkdir -p '//trim(output_path)//'grid')
                end if
                
                open(unit = 77, file=trim(output_path)//'grid/'//trim(this%name)//'_xb.dat')
                open(unit = 78, file=trim(output_path)//'grid/'//trim(this%name)//'_xc.dat')
                do i = 1, this%npts_b
                    write(77,*) i, this%xb(i)
                    write(78,*) i, this%xc(i)
                end do
                close(77)
                close(78)

            end subroutine

    end subroutine grid_generate_integer

end module