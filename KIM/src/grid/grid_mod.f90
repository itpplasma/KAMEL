module grid_m

    use KIM_kinds_m, only: dp

    implicit none

    real(dp) :: r_min
    real(dp) :: r_plas
    integer :: l_space_dim ! dimension of spline grid
    integer :: r_space_dim ! dimension of r grid
    integer :: k_space_dim ! dimension of kr space grid
    integer :: rg_space_dim
    integer :: spline_base
    ! Grid spacing modes (strings): "equidistant", "non-equidistant", "adaptive"
    character(len=32) :: grid_spacing_rg = "adaptive"
    character(len=32) :: grid_spacing_xl = "adaptive"
    integer :: gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp
    real(dp):: Larmor_skip_factor
    real(dp):: width_res, ampl_res, hrmax_scaling
    character(len=64) :: theta_integration ! RKF45, GaussLegendre, or QUADPACK
    character(len=64) :: theta_integration_method = "RKF45" ! Default to RKF45 for backward compatibility
    real(dp) :: rkf45_atol = 1.0d-9  ! Absolute tolerance for RKF45 adaptive integration
    real(dp) :: rkf45_rtol = 1.0d-6  ! Relative tolerance for RKF45 adaptive integration
    real(dp) :: kernel_taper_skip_threshold = 1.0d-6  ! Skip element calc when taper weight below this
    
    ! QUADPACK integration parameters
    character(len=32) :: quadpack_algorithm = "QAG"  ! QAG or QAGS
    integer :: quadpack_key = 6  ! Gauss-Kronrod rule: 1-6 for 15-61 points
    integer :: quadpack_limit = 500  ! Maximum number of subdivisions
    real(dp) :: quadpack_epsabs = 1.0d-10  ! Absolute tolerance for QUADPACK
    real(dp) :: quadpack_epsrel = 1.0d-10  ! Relative tolerance for QUADPACK
    logical :: quadpack_use_u_substitution = .true.  ! Use u=sin(theta/2) transformation

    integer :: nder=2
    integer :: npoi_der=4

    real(dp), dimension(:), allocatable :: xl  ! xl grid (real space)

    complex(dp), dimension(:,:), allocatable :: varphi_lkr

    real(dp) :: gg_factor = 1.0
    real(dp) :: gg_width = 0.0
    real(dp) :: gg_r_res = 0.0!95.34

    ! k-space specific
    real(dp) :: kr_grid_ampl_res
    real(dp) :: kr_grid_width_res
    real(dp) :: kr_res = 0.0d0
    
    type grid_type
        integer :: npts_b, npts_c, npts
        integer, dimension(:), allocatable :: ipbeg, ipend
        real(dp) :: min_val
        real(dp) :: max_val
        real(dp) :: hrmax
        real(dp), dimension(:), allocatable :: xb
        real(dp), dimension(:), allocatable :: xc
        real(dp), dimension(:,:), allocatable :: deriv_coef
        real(dp), dimension(:,:), allocatable :: deriv2_coef
        real(dp), dimension(:,:), allocatable :: reint_coef
        character(len=:), allocatable :: name
        contains
            procedure :: grid_init
            procedure :: grid_generate
            procedure :: grid_generate_linear
    end type grid_type

    type(grid_type) :: rg_grid, xl_grid, kr_grid, krp_grid

    contains

    subroutine ensure_node_at_r_res(this)
        use resonances_mod, only: r_res
        use KIM_kinds_m, only: dp
        implicit none
        class(grid_type), intent(inout) :: this

        integer :: i, insert_pos
        real(dp), allocatable :: new_xb(:), new_xc(:)
        real(dp), parameter :: tol = 1.0d-12

        if (.not. allocated(this%xb)) return
        if (.not. allocated(this%xc)) return

        if (r_res <= this%min_val + tol) return
        if (r_res >= this%max_val - tol) return

        do i = 1, this%npts_b
            if (abs(this%xb(i) - r_res) <= tol) return
        end do

        insert_pos = -1
        do i = 1, this%npts_b - 1
            if (this%xb(i) < r_res .and. r_res < this%xb(i+1)) then
                insert_pos = i
                exit
            end if
        end do
        if (insert_pos < 0) return

        allocate(new_xb(this%npts_b + 1))
        new_xb(1:insert_pos) = this%xb(1:insert_pos)
        new_xb(insert_pos+1) = r_res
        new_xb(insert_pos+2:this%npts_b+1) = this%xb(insert_pos+1:this%npts_b)

        allocate(new_xc(this%npts_b))
        do i = 1, size(new_xc)
            new_xc(i) = 0.5d0 * (new_xb(i) + new_xb(i+1))
        end do

        deallocate(this%xb)
        deallocate(this%xc)
        this%npts_b = this%npts_b + 1
        this%npts_c = this%npts_b - 1
        allocate(this%xb(this%npts_b), this%xc(this%npts_c))
        this%xb = new_xb
        this%xc = new_xc

        deallocate(new_xb)
        deallocate(new_xc)
    end subroutine ensure_node_at_r_res

    subroutine grid_init(this, npts, min_val, max_val, name)

        implicit none

        class(grid_type), intent(inout) :: this

        integer, intent(inout) :: npts
        real(dp), intent(in) :: min_val, max_val
        character(len=*), intent(in) :: name

        real(dp) :: x_current, x_next
        real(dp) :: recnsp

        this%npts = npts
        this%npts_b = npts
        this%min_val = min_val
        this%max_val = max_val
        allocate(character(len=len(name)) :: this%name)
        this%name = name

        ! set parameters for grid spacing. grid_spacing=1: quidistant grid, grid_spacing=2: non-equidistant grid
        !if (grid_spacing == 1) then
            !width_res = 1.0
            !ampl_res = 0.0
        !elseif (grid_spacing == 2) then
            !width_res = 3.0
            !ampl_res = 0.3
        !else 
            !width_res = 0.2
            !ampl_res = 15.0
        !end if

        this%hrmax = hrmax_scaling * (this%max_val - this%min_val) / (this%npts_b)

        this%npts_b = 1
        x_current = this%min_val

        do while(x_current .lt. this%max_val)
            call recnsplit(x_current, recnsp)
            x_next = x_current + this%hrmax / recnsp
            call recnsplit(x_next, recnsp)
            x_current = 0.5d0 * (x_next + x_current + this%hrmax / recnsp)
            this%npts_b = this%npts_b + 1
        enddo

        this%npts_c = this%npts_b - 1

    end subroutine

    subroutine grid_generate(this)

        use resonances_mod, only: r_res, index_rg_res
        use config_m, only: fdebug, output_path

        implicit none

        class(grid_type), intent(inout) :: this

        real(dp) :: x_current, x_next
        integer :: ipoib, ipb, ipe
        real(dp), dimension(:,:), allocatable :: coef
        real(dp) :: recnsp
        
        allocate(this%xb(this%npts_b), this%xc(this%npts_c))
        allocate(coef(0:nder,npoi_der))


        x_current = this%min_val
        this%xb(1) = x_current

        do ipoib=2, this%npts_b
            call recnsplit(x_current, recnsp)
            x_next = x_current + this%hrmax / recnsp
            call recnsplit(x_next, recnsp)
            x_current = 0.5d0 * (x_next + x_current + this%hrmax / recnsp)
            this%xb(ipoib) = x_current
            this%xc(ipoib-1) = 0.5 * (this%xb(ipoib-1) + this%xb(ipoib))
        enddo

        ! call ensure_node_at_r_res(this)

        ! get index for resonant radius
        call binsrc(abs(this%xb), 1, this%npts_b, abs(r_res), index_rg_res)

        if(npoi_der .gt. this%npts_c) then
            write(*,*) '! Error : not enough grid points for derivatives'
            stop
        endif

        allocate(this%deriv_coef(npoi_der, this%npts_b))
        allocate(this%deriv2_coef(npoi_der, this%npts_b))
        allocate(this%reint_coef(npoi_der, this%npts_b))
        allocate(this%ipbeg(this%npts_b))
        allocate(this%ipend(this%npts_b))

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
            this%ipbeg(ipoib) = ipb
            this%ipend(ipoib) = ipe
            call plag_coeff(npoi_der, nder, this%xb(ipoib), this%xc(ipb:ipe), coef)

            this%reint_coef(:, ipoib) = coef(0,:)
            this%deriv_coef(:, ipoib) = coef(1,:)
            this%deriv2_coef(:, ipoib) = coef(2,:)

        enddo

        deallocate(coef)

        call write_new_grid

        contains

            subroutine write_new_grid

                implicit none
                integer :: i, ierr
                logical :: dir_exists

                inquire(file=trim(output_path)//'grid/', exist=dir_exists)

                if (.not. dir_exists) then
                    ! Try to create the directory using a system call (POSIX mkdir)
                    call execute_command_line("mkdir -p " // trim(output_path)//'grid', exitstat=ierr)

                    if (ierr /= 0) then
                    print *, "Error: Could not create directory. Exit status = ", ierr
                    else
                    print *, "Directory created successfully."
                    end if
                end if

                open(unit = 77, file=trim(output_path)//'grid/'//trim(this%name)//'_xb.dat')
                open(unit = 78, file=trim(output_path)//'grid/'//trim(this%name)//'_xc.dat')
                do i = 1, size(this%xb)
                    write(77,*) i, this%xb(i)
                    if (i > this%npts_c) cycle
                    write(78,*) i, this%xc(i)
                end do
                close(77)
                close(78)

            end subroutine

    end subroutine grid_generate

    subroutine grid_generate_linear(this)

        use resonances_mod, only: r_res, index_rg_res
        use config_m, only: fdebug, output_path

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

        ! call ensure_node_at_r_res(this) ! could be used for adding r_res point in grid, but introduces some small oscillations

        ! get index for resonant radius
        call binsrc(abs(this%xb), 1, this%npts_b, abs(r_res), index_rg_res)

        if(npoi_der .gt. this%npts_c) then
            write(*,*) '! Error : not enough grid points for derivatives'
            stop
        endif

        allocate(this%deriv_coef(npoi_der, this%npts_c))
        allocate(this%deriv2_coef(npoi_der, this%npts_c))
        allocate(this%reint_coef(npoi_der, this%npts_c))
        allocate(this%ipbeg(this%npts_b))
        allocate(this%ipend(this%npts_b))

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
            this%ipbeg(ipoib) = ipb
            this%ipend(ipoib) = ipe
            
            call plag_coeff(npoi_der, nder, this%xb(ipoib), this%xc(ipb:ipe), coef)

            this%reint_coef(:, ipoib) = coef(0,:)
            this%deriv_coef(:, ipoib) = coef(1,:)
            this%deriv2_coef(:, ipoib) = coef(2,:)

        enddo

        deallocate(coef)

        call write_new_grid

        if (fdebug == 1) write(*,*) "Debug: exiting gengrid"

        contains

            subroutine write_new_grid

                implicit none
                integer :: i
                
                open(unit = 77, file=trim(output_path)//'grid/'//trim(this%name)//'_xb.dat')
                open(unit = 78, file=trim(output_path)//'grid/'//trim(this%name)//'_xc.dat')
                do i = 1, this%npts_b
                    write(77,*) i, this%xb(i)
                    if (i > this%npts_c) cycle
                    write(78,*) i, this%xc(i)
                end do
                close(77)
                close(78)

            end subroutine

    end subroutine grid_generate_linear

    subroutine calc_mass_matrix(M_mat)

        use KIM_kinds_m, only: dp
        use config_m, only: output_path
        use IO_collection_m, only: write_matrix

        implicit none

        real(dp), intent(inout) :: M_mat(:,:)
        real(dp) :: h

        integer :: i, n

        M_mat = 0.0d0
        n = xl_grid%npts_b

        do i = 1, xl_grid%npts_b-1
            h = xl_grid%xb(i+1) - xl_grid%xb(i)

            M_mat(i,  i  ) = M_mat(i,  i  ) + 2.0d0*h/6.0d0
            M_mat(i,  i+1) = M_mat(i,  i+1) + 1.0d0*h/6.0d0
            M_mat(i+1,i  ) = M_mat(i+1,i  ) + 1.0d0*h/6.0d0
            M_mat(i+1,i+1) = M_mat(i+1,i+1) + 2.0d0*h/6.0d0
        end do

        ! Enforce Dirichlet BC at right boundary: Phi_n = 0
        ! This ensures consistency with A_mat boundary conditions
        M_mat(n,:) = 0.0d0
        M_mat(:,n) = 0.0d0
        M_mat(n,n) = 1.0d0

        call write_matrix(trim(output_path)//'kernel/mass_matrix.dat', M_mat, xl_grid%npts_b, xl_grid%npts_b)
    end subroutine


end module
