!> Zero-finding module for locating all zeros of complex analytic functions
!! within a rectangular region using the argument principle and Newton's method
module kilca_zerofind_m
    use kilca_types_m
    implicit none
    
    private
    
    ! Public types
    public :: zerofind_settings_t
    public :: newton_settings_t
    public :: rectangle_t
    public :: rectangle_tree_t
    
    ! Public procedures
    public :: zerofind_settings_create, zerofind_settings_destroy, zerofind_settings_validate
    public :: zerofind_find_roots
    public :: zerofind_calc_winding_number
    public :: zerofind_calc_arg_change
    public :: zerofind_check_convergence
    public :: newton_settings_create, newton_settings_destroy
    public :: newton_solve
    public :: rectangle_create, rectangle_destroy, rectangle_contains, rectangle_subdivide
    public :: rect_tree_create, rect_tree_destroy, rect_tree_adaptive_subdivide, rect_tree_count_leaves
    
    !> Settings for the zero-finding algorithm
    type :: zerofind_settings_t
        logical :: use_winding = .true.              !< Use winding number evaluation
        integer :: max_partition_level = 128         !< Maximum recursion level
        integer :: min_recursion_level = 4           !< Minimum recursion level
        integer :: target_n_zeros = 20               !< Expected number of zeros
        integer :: n_split_x = 4                     !< Grid points in x for Newton starts
        integer :: n_split_y = 4                     !< Grid points in y for Newton starts
        real(dp) :: eps_abs = 1.0e-12_dp    !< Absolute tolerance
        real(dp) :: eps_rel = 1.0e-12_dp    !< Relative tolerance
        real(dp) :: eps_residual = 1.0e-10_dp !< Residual tolerance
        real(dp) :: interpolation_error = 1.0e-2_dp !< Interpolation tolerance
        real(dp) :: jump_error = 1.0e-2_dp   !< Discontinuity tolerance
        logical :: handle_multiplicities = .false.   !< Handle multiple roots
        integer :: debug_level = 0                   !< Debug output level
        integer :: print_level = 1                   !< Print output level
    end type zerofind_settings_t
    
    !> Settings for Newton's method
    type :: newton_settings_t
        integer :: max_iter = 50                     !< Maximum iterations
        real(dp) :: eps_abs = 1.0e-12_dp    !< Absolute tolerance
        real(dp) :: eps_rel = 1.0e-12_dp    !< Relative tolerance
        real(dp) :: eps_func = 1.0e-12_dp   !< Function value tolerance
        real(dp) :: delta = 1.0e-8_dp       !< Finite difference step
        logical :: use_numerical_deriv = .false.     !< Use numerical derivative
    end type newton_settings_t
    
    !> Rectangular region in complex plane
    type :: rectangle_t
        real(dp) :: xmin, xmax, ymin, ymax
        integer :: level = 0                         !< Subdivision level
        logical :: is_subdivided = .false.
        type(rectangle_t), pointer :: sub_rects(:) => null()
    end type rectangle_t
    
    !> Tree structure for adaptive subdivision
    type :: rectangle_tree_t
        type(rectangle_t) :: root
        integer :: max_level = 64
        integer :: n_rectangles = 0
    end type rectangle_tree_t
    
    ! Interface for user-defined complex functions
    abstract interface
        function complex_func_interface(z) result(f)
            import :: dp
            complex(dp), intent(in) :: z
            complex(dp) :: f
        end function complex_func_interface
    end interface
    
contains
    
    !---------------------------------------------------------------------------
    ! Zero-finding settings management
    !---------------------------------------------------------------------------
    
    !> Create and initialize zero-finding settings
    subroutine zerofind_settings_create(settings, ierr)
        type(zerofind_settings_t), intent(out) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Initialize with default values (already set in type definition)
        settings%use_winding = .true.
        settings%max_partition_level = 128
        settings%min_recursion_level = 4
        settings%target_n_zeros = 20
        settings%n_split_x = 4
        settings%n_split_y = 4
        settings%eps_abs = 1.0e-12_dp
        settings%eps_rel = 1.0e-12_dp
        settings%eps_residual = 1.0e-10_dp
        settings%interpolation_error = 1.0e-2_dp
        settings%jump_error = 1.0e-2_dp
        settings%handle_multiplicities = .false.
        settings%debug_level = 0
        settings%print_level = 1
        
    end subroutine zerofind_settings_create
    
    !> Destroy zero-finding settings
    subroutine zerofind_settings_destroy(settings, ierr)
        type(zerofind_settings_t), intent(inout) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        ! Nothing to deallocate for this simple type
        
    end subroutine zerofind_settings_destroy
    
    !> Validate zero-finding settings
    subroutine zerofind_settings_validate(settings, ierr)
        type(zerofind_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (settings%max_partition_level < 1) then
            ierr = -1
            return
        end if
        
        if (settings%min_recursion_level < 0 .or. &
            settings%min_recursion_level > settings%max_partition_level) then
            ierr = -2
            return
        end if
        
        if (settings%eps_abs < 0.0_dp .or. settings%eps_rel < 0.0_dp) then
            ierr = -3
            return
        end if
        
        if (settings%n_split_x < 1 .or. settings%n_split_y < 1) then
            ierr = -4
            return
        end if
        
    end subroutine zerofind_settings_validate
    
    !---------------------------------------------------------------------------
    ! Newton settings management
    !---------------------------------------------------------------------------
    
    !> Create Newton method settings
    subroutine newton_settings_create(settings, ierr)
        type(newton_settings_t), intent(out) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        settings%max_iter = 50
        settings%eps_abs = 1.0e-12_dp
        settings%eps_rel = 1.0e-12_dp
        settings%eps_func = 1.0e-12_dp
        settings%delta = 1.0e-8_dp
        settings%use_numerical_deriv = .false.
        
    end subroutine newton_settings_create
    
    !> Destroy Newton method settings
    subroutine newton_settings_destroy(settings, ierr)
        type(newton_settings_t), intent(inout) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        ! Nothing to deallocate
        
    end subroutine newton_settings_destroy
    
    !---------------------------------------------------------------------------
    ! Rectangle management
    !---------------------------------------------------------------------------
    
    !> Create a rectangle
    subroutine rectangle_create(rect, xmin, xmax, ymin, ymax, ierr)
        type(rectangle_t), intent(out) :: rect
        real(dp), intent(in) :: xmin, xmax, ymin, ymax
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (xmin >= xmax .or. ymin >= ymax) then
            ierr = -1
            return
        end if
        
        rect%xmin = xmin
        rect%xmax = xmax
        rect%ymin = ymin
        rect%ymax = ymax
        rect%level = 0
        rect%is_subdivided = .false.
        nullify(rect%sub_rects)
        
    end subroutine rectangle_create
    
    !> Destroy a rectangle and its subdivisions
    recursive subroutine rectangle_destroy(rect, ierr)
        type(rectangle_t), intent(inout) :: rect
        integer, intent(out) :: ierr
        integer :: i
        
        ierr = 0
        
        if (rect%is_subdivided .and. associated(rect%sub_rects)) then
            do i = 1, 4
                call rectangle_destroy(rect%sub_rects(i), ierr)
            end do
            deallocate(rect%sub_rects)
            nullify(rect%sub_rects)
        end if
        
        rect%is_subdivided = .false.
        
    end subroutine rectangle_destroy
    
    !> Check if rectangle contains a point
    subroutine rectangle_contains(rect, z, contains_point)
        type(rectangle_t), intent(in) :: rect
        complex(dp), intent(in) :: z
        logical, intent(out) :: contains_point
        real(dp) :: x, y
        
        x = real(z, dp)
        y = aimag(z)
        
        contains_point = (x >= rect%xmin .and. x <= rect%xmax .and. &
                         y >= rect%ymin .and. y <= rect%ymax)
        
    end subroutine rectangle_contains
    
    !> Subdivide rectangle into 4 sub-rectangles
    subroutine rectangle_subdivide(rect, ierr)
        type(rectangle_t), intent(inout) :: rect
        integer, intent(out) :: ierr
        real(dp) :: xmid, ymid
        
        ierr = 0
        
        if (rect%is_subdivided) then
            ierr = -1
            return
        end if
        
        allocate(rect%sub_rects(4), stat=ierr)
        if (ierr /= 0) return
        
        xmid = 0.5_dp * (rect%xmin + rect%xmax)
        ymid = 0.5_dp * (rect%ymin + rect%ymax)
        
        ! Bottom-left
        call rectangle_create(rect%sub_rects(1), rect%xmin, xmid, rect%ymin, ymid, ierr)
        rect%sub_rects(1)%level = rect%level + 1
        
        ! Bottom-right
        call rectangle_create(rect%sub_rects(2), xmid, rect%xmax, rect%ymin, ymid, ierr)
        rect%sub_rects(2)%level = rect%level + 1
        
        ! Top-left
        call rectangle_create(rect%sub_rects(3), rect%xmin, xmid, ymid, rect%ymax, ierr)
        rect%sub_rects(3)%level = rect%level + 1
        
        ! Top-right
        call rectangle_create(rect%sub_rects(4), xmid, rect%xmax, ymid, rect%ymax, ierr)
        rect%sub_rects(4)%level = rect%level + 1
        
        rect%is_subdivided = .true.
        
    end subroutine rectangle_subdivide
    
    !---------------------------------------------------------------------------
    ! Rectangle tree management
    !---------------------------------------------------------------------------
    
    !> Create rectangle tree
    subroutine rect_tree_create(tree, xmin, xmax, ymin, ymax, ierr)
        type(rectangle_tree_t), intent(out) :: tree
        real(dp), intent(in) :: xmin, xmax, ymin, ymax
        integer, intent(out) :: ierr
        
        call rectangle_create(tree%root, xmin, xmax, ymin, ymax, ierr)
        if (ierr /= 0) return
        
        tree%max_level = 64
        tree%n_rectangles = 1
        
    end subroutine rect_tree_create
    
    !> Destroy rectangle tree
    subroutine rect_tree_destroy(tree, ierr)
        type(rectangle_tree_t), intent(inout) :: tree
        integer, intent(out) :: ierr
        
        call rectangle_destroy(tree%root, ierr)
        tree%n_rectangles = 0
        
    end subroutine rect_tree_destroy
    
    !> Adaptive subdivision based on function properties
    subroutine rect_tree_adaptive_subdivide(tree, func, max_level, ierr)
        type(rectangle_tree_t), intent(inout) :: tree
        procedure(complex_func_interface) :: func
        integer, intent(in) :: max_level
        integer, intent(out) :: ierr
        
        ierr = 0
        
        call subdivide_recursive(tree%root, func, max_level, tree%n_rectangles, ierr)
        
    end subroutine rect_tree_adaptive_subdivide
    
    !> Recursive subdivision helper
    recursive subroutine subdivide_recursive(rect, func, max_level, n_rects, ierr)
        type(rectangle_t), intent(inout) :: rect
        procedure(complex_func_interface) :: func
        integer, intent(in) :: max_level
        integer, intent(inout) :: n_rects
        integer, intent(out) :: ierr
        integer :: i
        logical :: should_subdivide
        
        ierr = 0
        
        if (rect%level >= max_level) return
        
        ! Check if subdivision is needed based on function behavior
        call check_subdivision_needed(rect, func, should_subdivide)
        
        if (should_subdivide) then
            call rectangle_subdivide(rect, ierr)
            if (ierr /= 0) return
            
            n_rects = n_rects + 3  ! Added 3 more rectangles (4 total minus original)
            
            do i = 1, 4
                call subdivide_recursive(rect%sub_rects(i), func, max_level, n_rects, ierr)
                if (ierr /= 0) return
            end do
        end if
        
    end subroutine subdivide_recursive
    
    !> Check if subdivision is needed based on function variation
    subroutine check_subdivision_needed(rect, func, should_subdivide)
        type(rectangle_t), intent(in) :: rect
        procedure(complex_func_interface) :: func
        logical, intent(out) :: should_subdivide
        complex(dp) :: z1, z2, z3, z4, f1, f2, f3, f4
        real(dp) :: variation
        
        ! Evaluate function at corners
        z1 = cmplx(rect%xmin, rect%ymin, dp)
        z2 = cmplx(rect%xmax, rect%ymin, dp)
        z3 = cmplx(rect%xmax, rect%ymax, dp)
        z4 = cmplx(rect%xmin, rect%ymax, dp)
        
        f1 = func(z1)
        f2 = func(z2)
        f3 = func(z3)
        f4 = func(z4)
        
        ! Simple variation measure
        variation = abs(f2 - f1) + abs(f3 - f2) + abs(f4 - f3) + abs(f1 - f4)
        
        ! Subdivide if variation is large or function values are small (possible zero)
        should_subdivide = (variation > 0.1_dp) .or. &
                          (min(abs(f1), abs(f2), abs(f3), abs(f4)) < 0.1_dp)
        
    end subroutine check_subdivision_needed
    
    !> Count leaf rectangles in tree
    recursive subroutine rect_tree_count_leaves(tree, n_leaves)
        type(rectangle_tree_t), intent(in) :: tree
        integer, intent(out) :: n_leaves
        
        n_leaves = 0
        call count_leaves_recursive(tree%root, n_leaves)
        
    end subroutine rect_tree_count_leaves
    
    !> Recursive leaf counting helper
    recursive subroutine count_leaves_recursive(rect, n_leaves)
        type(rectangle_t), intent(in) :: rect
        integer, intent(inout) :: n_leaves
        integer :: i
        
        if (.not. rect%is_subdivided) then
            n_leaves = n_leaves + 1
        else
            do i = 1, 4
                call count_leaves_recursive(rect%sub_rects(i), n_leaves)
            end do
        end if
        
    end subroutine count_leaves_recursive
    
    !---------------------------------------------------------------------------
    ! Winding number and argument calculations
    !---------------------------------------------------------------------------
    
    !> Calculate winding number around rectangle
    subroutine zerofind_calc_winding_number(func, xmin, xmax, ymin, ymax, &
                                           settings, winding_num, ierr)
        procedure(complex_func_interface) :: func
        real(dp), intent(in) :: xmin, xmax, ymin, ymax
        type(zerofind_settings_t), intent(in) :: settings
        integer, intent(out) :: winding_num
        integer, intent(out) :: ierr
        real(dp) :: arg_change
        
        call zerofind_calc_arg_change(func, xmin, xmax, ymin, ymax, &
                                     settings, arg_change, ierr)
        if (ierr /= 0) return
        
        ! Winding number is arg_change / 2π
        winding_num = nint(arg_change / (2.0_dp * PI))
        
    end subroutine zerofind_calc_winding_number
    
    !> Calculate argument change around rectangle
    subroutine zerofind_calc_arg_change(func, xmin, xmax, ymin, ymax, &
                                        settings, arg_change, ierr)
        procedure(complex_func_interface) :: func
        real(dp), intent(in) :: xmin, xmax, ymin, ymax
        type(zerofind_settings_t), intent(in) :: settings
        real(dp), intent(out) :: arg_change
        integer, intent(out) :: ierr
        integer :: n_points, i
        real(dp), allocatable :: path_arg(:)
        complex(dp), allocatable :: z_path(:), f_path(:)
        real(dp) :: t, dt
        
        ierr = 0
        
        ! Number of points along contour
        n_points = max(100, 4 * settings%min_recursion_level)
        allocate(z_path(n_points), f_path(n_points), path_arg(n_points), stat=ierr)
        if (ierr /= 0) return
        
        ! Parametrize rectangular contour
        dt = 4.0_dp / real(n_points - 1, dp)
        
        do i = 1, n_points
            t = real(i - 1, dp) * dt
            
            if (t <= 1.0_dp) then
                ! Bottom edge: (xmin,ymin) to (xmax,ymin)
                z_path(i) = cmplx(xmin + t * (xmax - xmin), ymin, dp)
            else if (t <= 2.0_dp) then
                ! Right edge: (xmax,ymin) to (xmax,ymax)
                z_path(i) = cmplx(xmax, ymin + (t - 1.0_dp) * (ymax - ymin), dp)
            else if (t <= 3.0_dp) then
                ! Top edge: (xmax,ymax) to (xmin,ymax)
                z_path(i) = cmplx(xmax - (t - 2.0_dp) * (xmax - xmin), ymax, dp)
            else
                ! Left edge: (xmin,ymax) to (xmin,ymin)
                z_path(i) = cmplx(xmin, ymax - (t - 3.0_dp) * (ymax - ymin), dp)
            end if
            
            f_path(i) = func(z_path(i))
            
            if (abs(f_path(i)) < tiny(1.0_dp)) then
                ierr = -1  ! Function too small on contour
                deallocate(z_path, f_path, path_arg)
                return
            end if
            
            path_arg(i) = atan2(aimag(f_path(i)), real(f_path(i), dp))
        end do
        
        ! Calculate total argument change
        arg_change = 0.0_dp
        do i = 2, n_points
            arg_change = arg_change + unwrap_angle(path_arg(i) - path_arg(i-1))
        end do
        
        deallocate(z_path, f_path, path_arg)
        
    end subroutine zerofind_calc_arg_change
    
    !> Unwrap angle difference to [-π, π]
    function unwrap_angle(angle) result(unwrapped)
        real(dp), intent(in) :: angle
        real(dp) :: unwrapped
        
        unwrapped = angle
        
        do while (unwrapped > PI)
            unwrapped = unwrapped - 2.0_dp * PI
        end do
        
        do while (unwrapped < -PI)
            unwrapped = unwrapped + 2.0_dp * PI
        end do
        
    end function unwrap_angle
    
    !---------------------------------------------------------------------------
    ! Newton's method
    !---------------------------------------------------------------------------
    
    !> Solve for zero using Newton's method
    subroutine newton_solve(func, dfunc, z0, settings, z_root, n_iter, ierr)
        procedure(complex_func_interface) :: func
        procedure(complex_func_interface) :: dfunc
        complex(dp), intent(in) :: z0
        type(newton_settings_t), intent(in) :: settings
        complex(dp), intent(out) :: z_root
        integer, intent(out) :: n_iter
        integer, intent(out) :: ierr
        complex(dp) :: z, z_new, f, df
        logical :: converged
        
        ierr = 0
        z = z0
        n_iter = 0
        converged = .false.
        
        do n_iter = 1, settings%max_iter
            f = func(z)
            
            ! Check function convergence
            if (abs(f) < settings%eps_func) then
                converged = .true.
                exit
            end if
            
            
            ! Calculate derivative
            if (settings%use_numerical_deriv) then
                df = numerical_derivative(func, z, settings%delta)
            else
                df = dfunc(z)
            end if
            
            if (abs(df) < tiny(1.0_dp)) then
                ierr = -2  ! Derivative too small
                exit
            end if
            
            ! Newton step
            z_new = z - f / df
            
            ! Check convergence
            call check_newton_convergence(z, z_new, settings%eps_abs, &
                                         settings%eps_rel, converged)
            
            if (converged) exit
            
            z = z_new
        end do
        
        if (.not. converged .and. ierr == 0) then
            ierr = -1  ! Max iterations reached
        end if
        
        z_root = z
        
    end subroutine newton_solve
    
    !> Calculate numerical derivative
    function numerical_derivative(func, z, h) result(df)
        procedure(complex_func_interface) :: func
        complex(dp), intent(in) :: z
        real(dp), intent(in) :: h
        complex(dp) :: df
        complex(dp) :: f1, f2, f3, f4
        complex(dp) :: dz
        
        ! 4-point central difference formula
        dz = cmplx(h, 0.0_dp, dp)
        f1 = func(z + dz)
        f2 = func(z - dz)
        
        dz = cmplx(0.0_dp, h, dp)
        f3 = func(z + dz)
        f4 = func(z - dz)
        
        df = cmplx((real(f1,dp) - real(f2,dp)) / (2.0_dp * h), &
                   (aimag(f3) - aimag(f4)) / (2.0_dp * h), dp)
        
    end function numerical_derivative
    
    !> Check Newton convergence
    subroutine check_newton_convergence(z_old, z_new, eps_abs, eps_rel, converged)
        complex(dp), intent(in) :: z_old, z_new
        real(dp), intent(in) :: eps_abs, eps_rel
        logical, intent(out) :: converged
        real(dp) :: diff
        
        diff = abs(z_new - z_old)
        converged = (diff < eps_abs) .or. (diff < eps_rel * abs(z_old))
        
    end subroutine check_newton_convergence
    
    !---------------------------------------------------------------------------
    ! Main zero-finding algorithm
    !---------------------------------------------------------------------------
    
    !> Find all roots in rectangular region
    subroutine zerofind_find_roots(func, xmin, xmax, ymin, ymax, settings, &
                                   roots, n_roots, ierr)
        procedure(complex_func_interface) :: func
        real(dp), intent(in) :: xmin, xmax, ymin, ymax
        type(zerofind_settings_t), intent(in) :: settings
        complex(dp), intent(out) :: roots(:)
        integer, intent(inout) :: n_roots
        integer, intent(out) :: ierr
        type(rectangle_tree_t) :: tree
        type(newton_settings_t) :: newton_settings
        complex(dp), allocatable :: candidate_roots(:)
        integer :: max_roots, n_candidates
        integer :: i, j, winding_num
        logical :: is_new_root
        
        ierr = 0
        n_roots = 0
        max_roots = size(roots)
        
        ! Validate settings
        call zerofind_settings_validate(settings, ierr)
        if (ierr /= 0) return
        
        ! Create subdivision tree
        call rect_tree_create(tree, xmin, xmax, ymin, ymax, ierr)
        if (ierr /= 0) return
        
        if (settings%use_winding) then
            ! Check winding number for whole region
            call zerofind_calc_winding_number(func, xmin, xmax, ymin, ymax, &
                                             settings, winding_num, ierr)
            
            if (ierr /= 0) then
                call rect_tree_destroy(tree, ierr)
                return
            end if
            
            if (winding_num == 0) then
                ! No zeros detected by winding number - try Newton anyway
            end if
            
            ! Adaptive subdivision based on function behavior
            call rect_tree_adaptive_subdivide(tree, func, settings%max_partition_level, ierr)
            if (ierr /= 0) then
                call rect_tree_destroy(tree, ierr)
                return
            end if
        end if
        
        
        ! Allocate space for candidate roots
        allocate(candidate_roots(max_roots * 2), stat=ierr)
        if (ierr /= 0) then
            call rect_tree_destroy(tree, ierr)
            return
        end if
        
        ! Setup Newton settings
        call newton_settings_create(newton_settings, ierr)
        newton_settings%eps_abs = settings%eps_abs
        newton_settings%eps_rel = settings%eps_rel
        newton_settings%eps_func = settings%eps_residual
        newton_settings%use_numerical_deriv = .true.
        
        ! Generate starting points and run Newton iterations
        n_candidates = 0
        call generate_newton_starts(tree%root, func, newton_settings, &
                                   candidate_roots, n_candidates, &
                                   settings%n_split_x, settings%n_split_y)
        
        
        ! Filter out duplicates and points outside region
        do i = 1, n_candidates
            is_new_root = .true.
            
            ! Check if already found
            do j = 1, n_roots
                if (abs(candidate_roots(i) - roots(j)) < settings%eps_abs) then
                    is_new_root = .false.
                    exit
                end if
            end do
            
            ! Check if inside region
            if (is_new_root) then
                if (real(candidate_roots(i), dp) < xmin .or. &
                    real(candidate_roots(i), dp) > xmax .or. &
                    aimag(candidate_roots(i)) < ymin .or. &
                    aimag(candidate_roots(i)) > ymax) then
                    is_new_root = .false.
                end if
            end if
            
            ! Add to roots list
            if (is_new_root .and. n_roots < max_roots) then
                n_roots = n_roots + 1
                roots(n_roots) = candidate_roots(i)
            end if
        end do
        
        ! Clean up
        deallocate(candidate_roots)
        call newton_settings_destroy(newton_settings, ierr)
        call rect_tree_destroy(tree, ierr)
        
    end subroutine zerofind_find_roots
    
    !> Generate Newton starting points in rectangle
    recursive subroutine generate_newton_starts(rect, func, newton_settings, &
                                               candidates, n_candidates, nx, ny)
        type(rectangle_t), intent(in) :: rect
        procedure(complex_func_interface) :: func
        type(newton_settings_t), intent(in) :: newton_settings
        complex(dp), intent(inout) :: candidates(:)
        integer, intent(inout) :: n_candidates
        integer, intent(in) :: nx, ny
        integer :: i, j, n_iter, ierr
        real(dp) :: dx, dy, x, y
        complex(dp) :: z0, z_root
        
        if (rect%is_subdivided) then
            ! Process subdivisions
            do i = 1, 4
                call generate_newton_starts(rect%sub_rects(i), func, newton_settings, &
                                          candidates, n_candidates, nx, ny)
            end do
        else
            ! Generate grid of starting points
            dx = (rect%xmax - rect%xmin) / real(nx + 1, dp)
            dy = (rect%ymax - rect%ymin) / real(ny + 1, dp)
            
            
            do i = 1, nx
                x = rect%xmin + real(i, dp) * dx
                do j = 1, ny
                    y = rect%ymin + real(j, dp) * dy
                    z0 = cmplx(x, y, dp)
                    
                    ! Run Newton iteration
                    call newton_solve(func, dummy_deriv, z0, newton_settings, &
                                     z_root, n_iter, ierr)
                    
                    if (ierr == 0 .and. n_candidates < size(candidates)) then
                        n_candidates = n_candidates + 1
                        candidates(n_candidates) = z_root
                    end if
                end do
            end do
        end if
        
    end subroutine generate_newton_starts
    
    !> Check convergence between two complex numbers
    subroutine zerofind_check_convergence(z1, z2, settings, converged)
        complex(dp), intent(in) :: z1, z2
        type(zerofind_settings_t), intent(in) :: settings
        logical, intent(out) :: converged
        real(dp) :: diff
        
        diff = abs(z2 - z1)
        converged = (diff < settings%eps_abs) .or. &
                   (diff < settings%eps_rel * max(abs(z1), abs(z2)))
        
    end subroutine zerofind_check_convergence
    
    !> Dummy derivative function (not used when numerical derivatives enabled)
    function dummy_deriv(z) result(df)
        complex(dp), intent(in) :: z
        complex(dp) :: df
        df = cmplx(0.0_dp, 0.0_dp, dp)
    end function dummy_deriv
    
end module kilca_zerofind_m