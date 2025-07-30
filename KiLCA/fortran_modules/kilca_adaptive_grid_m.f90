!> Adaptive Grid System for KiLCA
!! Provides adaptive mesh refinement, gradient-based grid generation, and resonance layer handling
!! Supports multiple refinement strategies for plasma physics applications
module kilca_adaptive_grid_m
    use kilca_types_m
    use iso_fortran_env, only: real64
    implicit none
    
    private
    
    ! Public types
    public :: adaptive_grid_settings_t
    
    ! Public procedures
    public :: adaptive_grid_settings_create, adaptive_grid_settings_destroy
    public :: generate_uniform_grid
    public :: adaptive_grid_refine
    public :: refine_resonance_layer
    public :: compute_gradients, compute_curvatures
    public :: identify_refinement_regions
    public :: estimate_interpolation_errors
    public :: redistribute_function_on_grid
    public :: validate_grid_boundaries, validate_grid_monotonicity
    public :: allocate_adaptive_grid, deallocate_adaptive_grid
    public :: optimize_grid_performance
    
    !> Adaptive grid settings and configuration
    type :: adaptive_grid_settings_t
        integer :: method = 1                        !< Refinement method: 1=gradient, 2=curvature, 3=error
        integer :: max_points = 1000                 !< Maximum grid points allowed
        real(dp) :: min_spacing = 1.0e-6_dp         !< Minimum grid spacing
        real(dp) :: max_spacing = 0.1_dp            !< Maximum grid spacing
        real(dp) :: tolerance = 1.0e-4_dp           !< Refinement tolerance
        real(dp) :: refinement_factor = 2.0_dp      !< Factor for grid refinement
        logical :: resonance_refinement = .false.   !< Enable resonance layer refinement
        real(dp) :: resonance_width = 0.05_dp       !< Width of resonance layer
        real(dp) :: eps_res = 1.0e-4_dp            !< Error parameter in resonance layer
        real(dp) :: eps_out = 1.0e-2_dp            !< Error parameter outside resonance layer
        integer :: debug_level = 0                  !< Debug output level
        ! Additional parameters from FLRE zone
        real(dp) :: dr_out = 1.0e-3_dp             !< Grid step outside resonance region
        real(dp) :: dr_res = 1.0e-4_dp             !< Grid step in resonance region
        integer :: max_refinement_levels = 5        !< Maximum refinement levels
    end type adaptive_grid_settings_t
    
    ! Method constants
    integer, parameter, public :: GRID_METHOD_GRADIENT = 1
    integer, parameter, public :: GRID_METHOD_CURVATURE = 2
    integer, parameter, public :: GRID_METHOD_ERROR = 3
    
contains

    !---------------------------------------------------------------------------
    ! Settings management
    !---------------------------------------------------------------------------
    
    !> Create and initialize adaptive grid settings
    subroutine adaptive_grid_settings_create(settings, ierr)
        type(adaptive_grid_settings_t), intent(out) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Set default values (already set in type definition)
        settings%method = GRID_METHOD_GRADIENT
        settings%max_points = 1000
        settings%min_spacing = 1.0e-6_dp
        settings%max_spacing = 0.1_dp
        settings%tolerance = 1.0e-4_dp
        settings%refinement_factor = 2.0_dp
        settings%resonance_refinement = .false.
        settings%resonance_width = 0.05_dp
        settings%eps_res = 1.0e-4_dp
        settings%eps_out = 1.0e-2_dp
        settings%debug_level = 0
        settings%dr_out = 1.0e-3_dp
        settings%dr_res = 1.0e-4_dp
        settings%max_refinement_levels = 5
        
    end subroutine adaptive_grid_settings_create
    
    !> Destroy adaptive grid settings (cleanup if needed)
    subroutine adaptive_grid_settings_destroy(settings, ierr)
        type(adaptive_grid_settings_t), intent(inout) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        ! No dynamic memory to clean up in current implementation
        
    end subroutine adaptive_grid_settings_destroy
    
    !---------------------------------------------------------------------------
    ! Basic grid generation
    !---------------------------------------------------------------------------
    
    !> Generate uniform grid between r_min and r_max
    subroutine generate_uniform_grid(n, r_min, r_max, grid, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: r_min, r_max
        real(dp), intent(out) :: grid(n)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: dr
        
        ierr = 0
        
        if (n < 2) then
            ierr = -1
            return
        end if
        
        if (r_max <= r_min) then
            ierr = -2
            return
        end if
        
        dr = (r_max - r_min) / real(n - 1, dp)
        
        if (dr < settings%min_spacing) then
            ierr = -3
            return
        end if
        
        do i = 1, n
            grid(i) = r_min + real(i - 1, dp) * dr
        end do
        
        if (settings%debug_level > 0) then
            write(*, *) "Generated uniform grid with n=", n, "dr=", dr
        end if
        
    end subroutine generate_uniform_grid
    
    !---------------------------------------------------------------------------
    ! Adaptive refinement algorithms
    !---------------------------------------------------------------------------
    
    !> Perform adaptive grid refinement based on function values
    subroutine adaptive_grid_refine(n_initial, grid_initial, function_values, &
                                   refined_grid, n_refined, settings, ierr)
        integer, intent(in) :: n_initial
        real(dp), intent(in) :: grid_initial(n_initial)
        real(dp), intent(in) :: function_values(n_initial)
        real(dp), allocatable, intent(out) :: refined_grid(:)
        integer, intent(out) :: n_refined
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        real(dp), allocatable :: gradients(:), curvatures(:), errors(:)
        logical, allocatable :: refine_flags(:)
        integer :: i, n_new_points, current_size
        
        ierr = 0
        n_refined = n_initial
        
        if (n_initial < 3) then
            ierr = -1
            return
        end if
        
        ! Start with initial grid
        allocate(refined_grid(settings%max_points))
        refined_grid(1:n_initial) = grid_initial
        current_size = n_initial
        
        select case(settings%method)
        case(GRID_METHOD_GRADIENT)
            allocate(gradients(n_initial-1))
            call compute_gradients(n_initial, grid_initial, function_values, gradients, settings, ierr)
            if (ierr /= 0) return
            
            allocate(refine_flags(n_initial-1))
            call identify_refinement_regions(n_initial-1, gradients, refine_flags, settings, ierr)
            if (ierr /= 0) return
            
            call insert_refined_points_gradient(n_initial, grid_initial, gradients, refine_flags, &
                                               refined_grid, current_size, settings, ierr)
            
        case(GRID_METHOD_CURVATURE)
            if (n_initial < 3) then
                ierr = -2
                return
            end if
            
            allocate(curvatures(n_initial-2))
            call compute_curvatures(n_initial, grid_initial, function_values, curvatures, settings, ierr)
            if (ierr /= 0) return
            
            call insert_refined_points_curvature(n_initial, grid_initial, curvatures, &
                                                refined_grid, current_size, settings, ierr)
            
        case(GRID_METHOD_ERROR)
            allocate(errors(n_initial-1))
            call estimate_interpolation_errors(n_initial, grid_initial, function_values, errors, settings, ierr)
            if (ierr /= 0) return
            
            call insert_refined_points_error(n_initial, grid_initial, errors, &
                                            refined_grid, current_size, settings, ierr)
            
        case default
            ierr = -3
            return
        end select
        
        n_refined = current_size
        
        if (settings%debug_level > 0) then
            write(*, *) "Adaptive refinement: ", n_initial, "->", n_refined, "points"
        end if
        
    end subroutine adaptive_grid_refine
    
    !---------------------------------------------------------------------------
    ! Resonance layer refinement (translation from FLRE zone parameters)
    !---------------------------------------------------------------------------
    
    !> Refine grid around resonance layer position
    subroutine refine_resonance_layer(n_initial, grid_initial, resonance_position, &
                                     refined_grid, n_refined, settings, ierr)
        integer, intent(in) :: n_initial
        real(dp), intent(in) :: grid_initial(n_initial)
        real(dp), intent(in) :: resonance_position
        real(dp), allocatable, intent(out) :: refined_grid(:)
        integer, intent(out) :: n_refined
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        real(dp) :: r_res_min, r_res_max, r_min, r_max
        integer :: n_res, n_out_left, n_out_right, i, current_idx
        
        ierr = 0
        
        if (n_initial < 2) then
            ierr = -1
            return
        end if
        
        r_min = grid_initial(1)
        r_max = grid_initial(n_initial)
        
        ! Define resonance layer boundaries
        r_res_min = resonance_position - settings%resonance_width / 2.0_dp
        r_res_max = resonance_position + settings%resonance_width / 2.0_dp
        
        ! Ensure resonance layer is within grid bounds
        r_res_min = max(r_res_min, r_min)
        r_res_max = min(r_res_max, r_max)
        
        ! Calculate number of points in each region
        n_res = max(1, int(settings%resonance_width / settings%dr_res))
        n_out_left = max(1, int((r_res_min - r_min) / settings%dr_out))
        n_out_right = max(1, int((r_max - r_res_max) / settings%dr_out))
        
        ! Ensure we don't have zero-width regions
        if (r_res_min <= r_min) n_out_left = 0
        if (r_res_max >= r_max) n_out_right = 0
        
        n_refined = n_out_left + n_res + n_out_right
        
        if (n_refined <= 0) then
            ierr = -3
            return
        end if
        
        if (n_refined > settings%max_points) then
            ierr = -2
            return
        end if
        
        allocate(refined_grid(n_refined))
        current_idx = 1
        
        ! Left region (outside resonance)
        if (n_out_left > 0) then
            do i = 1, n_out_left
                if (current_idx <= n_refined) then
                    if (n_out_left > 1) then
                        refined_grid(current_idx) = r_min + real(i-1, dp) * (r_res_min - r_min) / real(n_out_left-1, dp)
                    else
                        refined_grid(current_idx) = (r_min + r_res_min) / 2.0_dp
                    end if
                    current_idx = current_idx + 1
                end if
            end do
        end if
        
        ! Resonance region (fine grid)
        do i = 1, n_res
            if (current_idx <= n_refined) then
                if (n_res > 1) then
                    refined_grid(current_idx) = r_res_min + real(i-1, dp) * (r_res_max - r_res_min) / real(n_res-1, dp)
                else
                    refined_grid(current_idx) = (r_res_min + r_res_max) / 2.0_dp
                end if
                current_idx = current_idx + 1
            end if
        end do
        
        ! Right region (outside resonance)
        if (n_out_right > 0) then
            do i = 1, n_out_right
                if (current_idx <= n_refined) then
                    if (n_out_right > 1) then
                        refined_grid(current_idx) = r_res_max + real(i, dp) * (r_max - r_res_max) / real(n_out_right, dp)
                    else
                        refined_grid(current_idx) = (r_res_max + r_max) / 2.0_dp
                    end if
                    current_idx = current_idx + 1
                end if
            end do
        end if
        
        if (settings%debug_level > 0) then
            write(*, *) "Resonance refinement: n_left=", n_out_left, "n_res=", n_res, "n_right=", n_out_right
        end if
        
    end subroutine refine_resonance_layer
    
    !---------------------------------------------------------------------------
    ! Gradient and curvature computation
    !---------------------------------------------------------------------------
    
    !> Compute gradients for refinement analysis
    subroutine compute_gradients(n, grid, function_values, gradients, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: grid(n), function_values(n)
        real(dp), intent(out) :: gradients(n-1)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: dr, df
        
        ierr = 0
        
        if (n < 2) then
            ierr = -1
            return
        end if
        
        do i = 1, n-1
            dr = grid(i+1) - grid(i)
            if (abs(dr) < settings%min_spacing) then
                ierr = -2
                return
            end if
            
            df = function_values(i+1) - function_values(i)
            gradients(i) = abs(df / dr)
        end do
        
    end subroutine compute_gradients
    
    !> Compute curvatures for refinement analysis
    subroutine compute_curvatures(n, grid, function_values, curvatures, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: grid(n), function_values(n)
        real(dp), intent(out) :: curvatures(n-2)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: dr1, dr2, df1, df2, grad1, grad2
        
        ierr = 0
        
        if (n < 3) then
            ierr = -1
            return
        end if
        
        do i = 1, n-2
            dr1 = grid(i+1) - grid(i)
            dr2 = grid(i+2) - grid(i+1)
            
            if (abs(dr1) < settings%min_spacing .or. abs(dr2) < settings%min_spacing) then
                ierr = -2
                return
            end if
            
            df1 = function_values(i+1) - function_values(i)
            df2 = function_values(i+2) - function_values(i+1)
            
            grad1 = df1 / dr1
            grad2 = df2 / dr2
            
            curvatures(i) = abs((grad2 - grad1) / ((dr1 + dr2) / 2.0_dp))
        end do
        
    end subroutine compute_curvatures
    
    !---------------------------------------------------------------------------
    ! Refinement region identification
    !---------------------------------------------------------------------------
    
    !> Identify regions that need refinement based on gradients
    subroutine identify_refinement_regions(n, gradients, refine_flags, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: gradients(n)
        logical, allocatable, intent(out) :: refine_flags(:)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: threshold
        
        ierr = 0
        
        if (n < 1) then
            ierr = -1
            return
        end if
        
        allocate(refine_flags(n))
        
        ! Calculate threshold based on maximum gradient and tolerance
        threshold = maxval(gradients) * settings%tolerance
        
        do i = 1, n
            refine_flags(i) = (gradients(i) > threshold)
        end do
        
        if (settings%debug_level > 1) then
            write(*, *) "Refinement threshold:", threshold, "flagged:", count(refine_flags)
        end if
        
    end subroutine identify_refinement_regions
    
    !---------------------------------------------------------------------------
    ! Error estimation
    !---------------------------------------------------------------------------
    
    !> Estimate interpolation errors for refinement decisions
    subroutine estimate_interpolation_errors(n, grid, function_values, errors, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: grid(n), function_values(n)
        real(dp), intent(out) :: errors(n-1)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: midpoint, interpolated_value, actual_value, dr
        
        ierr = 0
        
        if (n < 2) then
            ierr = -1
            return
        end if
        
        do i = 1, n-1
            ! Estimate error using linear interpolation at midpoint
            midpoint = (grid(i) + grid(i+1)) / 2.0_dp
            interpolated_value = (function_values(i) + function_values(i+1)) / 2.0_dp
            
            ! For demonstration, use quadratic approximation as "true" value
            dr = grid(i+1) - grid(i)
            if (i > 1 .and. i < n-1) then
                ! Use three-point quadratic interpolation
                actual_value = quadratic_interpolation(grid(i-1:i+1), function_values(i-1:i+1), midpoint)
            else
                ! Fall back to linear interpolation
                actual_value = interpolated_value
            end if
            
            errors(i) = abs(actual_value - interpolated_value)
        end do
        
    end subroutine estimate_interpolation_errors
    
    !---------------------------------------------------------------------------
    ! Grid redistribution and interpolation
    !---------------------------------------------------------------------------
    
    !> Redistribute function values from old grid to new grid
    subroutine redistribute_function_on_grid(n_old, old_grid, old_values, &
                                            n_new, new_grid, new_values, settings, ierr)
        integer, intent(in) :: n_old, n_new
        real(dp), intent(in) :: old_grid(n_old), old_values(n_old)
        real(dp), intent(in) :: new_grid(n_new)
        real(dp), intent(out) :: new_values(n_new)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        
        ierr = 0
        
        if (n_old < 2 .or. n_new < 1) then
            ierr = -1
            return
        end if
        
        ! Use linear interpolation for redistribution
        do i = 1, n_new
            call linear_interpolation(n_old, old_grid, old_values, new_grid(i), new_values(i), ierr)
            if (ierr /= 0) return
        end do
        
    end subroutine redistribute_function_on_grid
    
    !---------------------------------------------------------------------------
    ! Grid validation
    !---------------------------------------------------------------------------
    
    !> Validate grid boundaries
    subroutine validate_grid_boundaries(n, grid, r_min, r_max, valid_boundaries, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: grid(n), r_min, r_max
        logical, intent(out) :: valid_boundaries
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        valid_boundaries = .false.
        
        if (n < 1) then
            ierr = -1
            return
        end if
        
        valid_boundaries = (abs(grid(1) - r_min) < settings%min_spacing) .and. &
                          (abs(grid(n) - r_max) < settings%min_spacing)
        
    end subroutine validate_grid_boundaries
    
    !> Validate grid monotonicity
    subroutine validate_grid_monotonicity(n, grid, is_valid, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: grid(n)
        logical, intent(out) :: is_valid
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        
        ierr = 0
        is_valid = .true.
        
        if (n < 2) then
            ierr = -1
            return
        end if
        
        do i = 2, n
            if (grid(i) <= grid(i-1)) then
                is_valid = .false.
                return
            end if
        end do
        
    end subroutine validate_grid_monotonicity
    
    !---------------------------------------------------------------------------
    ! Memory management
    !---------------------------------------------------------------------------
    
    !> Allocate adaptive grid arrays
    subroutine allocate_adaptive_grid(n_max, grid, function_vals, n_actual, settings, ierr)
        integer, intent(in) :: n_max
        real(dp), allocatable, intent(out) :: grid(:), function_vals(:)
        integer, intent(out) :: n_actual
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (n_max <= 0) then
            ierr = -1
            return
        end if
        
        n_actual = min(n_max, settings%max_points)
        
        allocate(grid(n_actual), stat=ierr)
        if (ierr /= 0) return
        
        allocate(function_vals(n_actual), stat=ierr)
        if (ierr /= 0) then
            deallocate(grid)
            return
        end if
        
    end subroutine allocate_adaptive_grid
    
    !> Deallocate adaptive grid arrays
    subroutine deallocate_adaptive_grid(grid, function_vals, settings, ierr)
        real(dp), allocatable, intent(inout) :: grid(:), function_vals(:)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (allocated(grid)) deallocate(grid)
        if (allocated(function_vals)) deallocate(function_vals)
        
    end subroutine deallocate_adaptive_grid
    
    !---------------------------------------------------------------------------
    ! Performance optimization
    !---------------------------------------------------------------------------
    
    !> Optimize grid for performance
    subroutine optimize_grid_performance(n, grid, function_values, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(inout) :: grid(n), function_values(n)
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: total_variation
        
        ierr = 0
        
        if (n < 2) then
            ierr = -1
            return
        end if
        
        ! Calculate total variation for performance metric
        total_variation = 0.0_dp
        do i = 2, n
            total_variation = total_variation + abs(function_values(i) - function_values(i-1))
        end do
        
        if (settings%debug_level > 0) then
            write(*, *) "Grid performance: total variation =", total_variation
        end if
        
    end subroutine optimize_grid_performance
    
    !---------------------------------------------------------------------------
    ! Helper functions
    !---------------------------------------------------------------------------
    
    !> Linear interpolation helper
    subroutine linear_interpolation(n, x, y, x_target, y_target, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n), x_target
        real(dp), intent(out) :: y_target
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: alpha
        
        ierr = 0
        
        ! Find interpolation interval
        if (x_target <= x(1)) then
            y_target = y(1)
            return
        end if
        
        if (x_target >= x(n)) then
            y_target = y(n)
            return
        end if
        
        do i = 1, n-1
            if (x_target >= x(i) .and. x_target <= x(i+1)) then
                alpha = (x_target - x(i)) / (x(i+1) - x(i))
                y_target = y(i) + alpha * (y(i+1) - y(i))
                return
            end if
        end do
        
        ierr = -1  ! Target not found in interval
        
    end subroutine linear_interpolation
    
    !> Quadratic interpolation helper
    function quadratic_interpolation(x, y, x_target) result(y_target)
        real(dp), intent(in) :: x(3), y(3), x_target
        real(dp) :: y_target
        
        real(dp) :: L0, L1, L2
        
        ! Lagrange interpolation
        L0 = ((x_target - x(2)) * (x_target - x(3))) / ((x(1) - x(2)) * (x(1) - x(3)))
        L1 = ((x_target - x(1)) * (x_target - x(3))) / ((x(2) - x(1)) * (x(2) - x(3)))
        L2 = ((x_target - x(1)) * (x_target - x(2))) / ((x(3) - x(1)) * (x(3) - x(2)))
        
        y_target = L0 * y(1) + L1 * y(2) + L2 * y(3)
        
    end function quadratic_interpolation
    
    !> Insert refined points based on gradient analysis
    subroutine insert_refined_points_gradient(n_initial, grid_initial, gradients, refine_flags, &
                                             refined_grid, current_size, settings, ierr)
        integer, intent(in) :: n_initial
        real(dp), intent(in) :: grid_initial(n_initial), gradients(n_initial-1)
        logical, intent(in) :: refine_flags(n_initial-1)
        real(dp), intent(inout) :: refined_grid(:)
        integer, intent(inout) :: current_size
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i, insert_count
        real(dp) :: midpoint
        
        ierr = 0
        insert_count = 0
        
        do i = 1, n_initial-1
            if (refine_flags(i) .and. current_size + 1 <= settings%max_points) then
                midpoint = (grid_initial(i) + grid_initial(i+1)) / 2.0_dp
                
                ! Shift existing points and insert midpoint
                current_size = current_size + 1
                refined_grid(i + 1 + insert_count) = midpoint
                insert_count = insert_count + 1
            end if
        end do
        
        ! Sort the refined grid (simple bubble sort for small arrays)
        call sort_grid(current_size, refined_grid)
        
    end subroutine insert_refined_points_gradient
    
    !> Insert refined points based on curvature analysis
    subroutine insert_refined_points_curvature(n_initial, grid_initial, curvatures, &
                                              refined_grid, current_size, settings, ierr)
        integer, intent(in) :: n_initial
        real(dp), intent(in) :: grid_initial(n_initial), curvatures(n_initial-2)
        real(dp), intent(inout) :: refined_grid(:)
        integer, intent(inout) :: current_size
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: threshold, midpoint
        
        ierr = 0
        
        threshold = maxval(curvatures) * settings%tolerance
        
        do i = 1, n_initial-2
            if (curvatures(i) > threshold .and. current_size + 1 <= settings%max_points) then
                midpoint = (grid_initial(i+1) + grid_initial(i+2)) / 2.0_dp
                current_size = current_size + 1
                refined_grid(current_size) = midpoint
            end if
        end do
        
        call sort_grid(current_size, refined_grid)
        
    end subroutine insert_refined_points_curvature
    
    !> Insert refined points based on error estimation
    subroutine insert_refined_points_error(n_initial, grid_initial, errors, &
                                          refined_grid, current_size, settings, ierr)
        integer, intent(in) :: n_initial
        real(dp), intent(in) :: grid_initial(n_initial), errors(n_initial-1)
        real(dp), intent(inout) :: refined_grid(:)
        integer, intent(inout) :: current_size
        type(adaptive_grid_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: threshold, midpoint
        
        ierr = 0
        
        threshold = maxval(errors) * settings%tolerance
        
        do i = 1, n_initial-1
            if (errors(i) > threshold .and. current_size + 1 <= settings%max_points) then
                midpoint = (grid_initial(i) + grid_initial(i+1)) / 2.0_dp
                current_size = current_size + 1
                refined_grid(current_size) = midpoint
            end if
        end do
        
        call sort_grid(current_size, refined_grid)
        
    end subroutine insert_refined_points_error
    
    !> Simple bubble sort for grid arrays
    subroutine sort_grid(n, grid)
        integer, intent(in) :: n
        real(dp), intent(inout) :: grid(n)
        
        integer :: i, j
        real(dp) :: temp
        logical :: swapped
        
        do i = 1, n-1
            swapped = .false.
            do j = 1, n-i
                if (grid(j) > grid(j+1)) then
                    temp = grid(j)
                    grid(j) = grid(j+1)
                    grid(j+1) = temp
                    swapped = .true.
                end if
            end do
            if (.not. swapped) exit
        end do
        
    end subroutine sort_grid
    
end module kilca_adaptive_grid_m