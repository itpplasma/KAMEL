module grid_m

    use KIM_kinds_m, only: dp

    implicit none
    private :: adaptive_next_node

    real(dp) :: r_min
    real(dp) :: r_plas
    integer :: l_space_dim ! dimension of spline grid
    integer :: r_space_dim ! dimension of r grid
    integer :: rg_space_dim
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
    real(dp), dimension(:,:), allocatable :: M_mat ! mass matrix

    complex(dp), dimension(:,:), allocatable :: varphi_lkr

    real(dp) :: gg_factor = 1.0
    real(dp) :: gg_width = 0.0
    real(dp) :: gg_r_res = 0.0!95.34

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
            procedure :: grid_init_equidistant
            procedure :: grid_generate
            procedure :: grid_generate_equidistant
    end type grid_type

    type(grid_type) :: rg_grid, xl_grid, kr_grid, krp_grid

    contains

    subroutine ensure_node_at_r_res(this)
        use kim_resonances_m, only: r_res
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

    subroutine grid_init_equidistant(this, npts, min_val, max_val, name)

        implicit none

        class(grid_type), intent(inout) :: this

        integer, intent(inout) :: npts
        real(dp), intent(in) :: min_val, max_val
        character(len=*), intent(in) :: name

        real(dp) :: x_current, x_next
        real(dp) :: recnsp

        this%npts = npts
        this%npts_b = npts
        this%npts_c = this%npts_b - 1
        this%min_val = min_val
        this%max_val = max_val
        if (allocated(this%name)) deallocate(this%name)
        allocate(character(len=len(name)) :: this%name)
        this%name = name

    end subroutine



    subroutine grid_init(this, npts, min_val, max_val, name)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        implicit none

        class(grid_type), intent(inout) :: this

        integer, intent(inout) :: npts
        real(dp), intent(in) :: min_val, max_val
        character(len=*), intent(in) :: name

        real(dp) :: x_current, path_sum

        if (npts < 2) error stop 'Adaptive grid needs two nominal nodes'
        if (.not. ieee_is_finite(min_val)) error stop 'Adaptive lower bound must be finite'
        if (.not. ieee_is_finite(max_val)) error stop 'Adaptive upper bound must be finite'
        if (max_val <= min_val) error stop 'Adaptive grid bounds must increase'
        if (.not. ieee_is_finite(width_res)) error stop 'Adaptive width must be finite'
        if (width_res <= 0.0_dp) error stop 'Adaptive width must be positive'
        if (.not. ieee_is_finite(ampl_res)) error stop 'Adaptive amplitude must be finite'
        if (ampl_res <= -1.0_dp) error stop 'Adaptive spacing factor must be positive'
        if (.not. ieee_is_finite(hrmax_scaling)) &
            error stop 'Adaptive spacing scale must be finite'
        if (hrmax_scaling <= 0.0_dp) error stop 'Adaptive spacing scale must be positive'

        this%npts = npts
        this%npts_b = npts
        this%min_val = min_val
        this%max_val = max_val
        if (allocated(this%name)) deallocate(this%name)
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
        if (.not. ieee_is_finite(this%hrmax)) error stop 'Adaptive spacing must be finite'
        if (this%hrmax <= 0.0_dp) error stop 'Adaptive spacing must be positive'

        this%npts_b = 1
        x_current = this%min_val
        path_sum = 0.0_dp

        do while(x_current .lt. this%max_val)
            x_current = adaptive_next_node(this, x_current, this%npts_b, path_sum)
            this%npts_b = this%npts_b + 1
        enddo

        this%npts_c = this%npts_b - 1

    end subroutine

    function adaptive_next_node(this, current, nsteps, path_sum) result(next)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none
        class(grid_type), intent(in) :: this
        real(dp), intent(in) :: current
        integer, intent(in) :: nsteps
        real(dp), intent(inout) :: path_sum
        real(dp) :: next, factor, step, increment, gamma, unit_roundoff, allowance

        call kim_recnsplit(current, factor)
        if (.not. ieee_is_finite(factor)) error stop 'Adaptive factor must be finite'
        if (factor <= 0.0_dp) error stop 'Adaptive factor must be positive'
        step = this%hrmax/factor
        next = current + step
        if (.not. ieee_is_finite(next)) error stop 'Adaptive predictor must be finite'
        call kim_recnsplit(next, factor)
        if (.not. ieee_is_finite(factor)) error stop 'Adaptive factor must be finite'
        if (factor <= 0.0_dp) error stop 'Adaptive factor must be positive'
        increment = 0.5_dp*(step + this%hrmax/factor)
        next = current + increment
        if (.not. ieee_is_finite(next)) error stop 'Adaptive node must be finite'
        path_sum = path_sum + increment
        if (.not. ieee_is_finite(path_sum)) error stop 'Adaptive path must be finite'
        if (next >= this%max_val) then
            next = this%max_val
        else
            ! Bound accumulated coordinate and positive-path summation roundoff.
            ! This is not a forward-error bound for the nonlinear spacing monitor.
            unit_roundoff = 0.5_dp*epsilon(1.0_dp)
            gamma = (2.0_dp*real(nsteps, dp) + 6.0_dp)*unit_roundoff
            if (gamma >= 1.0_dp) error stop 'Adaptive roundoff bound is unresolved'
            gamma = gamma/(1.0_dp - gamma)
            allowance = gamma*abs(this%min_val) + gamma*path_sum &
                + 2.0_dp*spacing(max(abs(this%min_val), abs(this%max_val)))
            if (this%max_val - next <= allowance) then
                if (allowance >= sqrt(epsilon(1.0_dp))*increment) &
                    error stop 'Adaptive terminal cell is numerically unresolved'
                next = this%max_val
            end if
        end if
        if (next <= current) error stop 'Adaptive grid cannot make positive progress'
    end function adaptive_next_node

    subroutine grid_generate(this)

        use kim_resonances_m, only: r_res, index_rg_res
        use config_m, only: output_path
        use logger_m, only: log_error

        implicit none

        class(grid_type), intent(inout) :: this

        real(dp) :: x_current, path_sum
        integer :: ipoib, ipb, ipe
        real(dp), dimension(:,:), allocatable :: coef

        if (allocated(this%xb)) deallocate(this%xb)
        if (allocated(this%xc)) deallocate(this%xc)
        allocate(this%xb(this%npts_b), this%xc(this%npts_c))
        allocate(coef(0:nder,npoi_der))


        x_current = this%min_val
        path_sum = 0.0_dp
        this%xb(1) = x_current

        do ipoib=2, this%npts_b
            x_current = adaptive_next_node(this, x_current, ipoib - 1, path_sum)
            this%xb(ipoib) = x_current
            this%xc(ipoib-1) = this%xb(ipoib-1) &
                + 0.5_dp*(this%xb(ipoib) - this%xb(ipoib-1))
            if (this%xc(ipoib-1) <= this%xb(ipoib-1)) &
                error stop 'Adaptive cell center must be inside its cell'
            if (this%xc(ipoib-1) >= this%xb(ipoib)) &
                error stop 'Adaptive cell center must be inside its cell'
        enddo
        if (this%xb(this%npts_b) /= this%max_val) &
            error stop 'Adaptive grid changed after node counting'

        ! call ensure_node_at_r_res(this)

        ! get index for resonant radius
        call binsrc(abs(this%xb), 1, this%npts_b, abs(r_res), index_rg_res)

        if(npoi_der .gt. this%npts_c) then
            call log_error('Not enough grid points for derivatives')
        endif

        if (allocated(this%deriv_coef)) deallocate(this%deriv_coef)
        if (allocated(this%deriv2_coef)) deallocate(this%deriv2_coef)
        if (allocated(this%reint_coef)) deallocate(this%reint_coef)
        if (allocated(this%ipbeg)) deallocate(this%ipbeg)
        if (allocated(this%ipend)) deallocate(this%ipend)
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

    end subroutine grid_generate

    subroutine grid_generate_equidistant(this)

        use kim_resonances_m, only: r_res, index_rg_res
        use config_m, only: output_path
        use logger_m, only: log_debug, log_error, fmt_val

        implicit none

        class(grid_type), intent(inout) :: this

        real(dp) :: h
        integer :: ipoib, ipb, ipe
        real(dp), dimension(:,:), allocatable :: coef

        if (allocated(this%xb)) deallocate(this%xb)
        if (allocated(this%xc)) deallocate(this%xc)
        allocate(this%xb(this%npts_b), this%xc(this%npts_c))

        h = (this%max_val - this%min_val) / this%npts_b
        call log_debug(trim(fmt_val('Equidistant grid spacing h', h, 'cm')))

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
            call log_error('Not enough grid points for derivatives')
        endif

        if (allocated(this%deriv_coef)) deallocate(this%deriv_coef)
        if (allocated(this%deriv2_coef)) deallocate(this%deriv2_coef)
        if (allocated(this%reint_coef)) deallocate(this%reint_coef)
        if (allocated(this%ipbeg)) deallocate(this%ipbeg)
        if (allocated(this%ipend)) deallocate(this%ipend)
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

    end subroutine grid_generate_equidistant

    subroutine calc_mass_matrix(M_mat)

        use KIM_kinds_m, only: dp
        use config_m, only: output_path

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

    end subroutine


end module
