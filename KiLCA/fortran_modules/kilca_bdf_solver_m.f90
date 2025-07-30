!> BDF (Backward Differentiation Formula) solver module for stiff ODEs
!! Implements variable-order BDF methods (orders 1-5) with adaptive step size control
module kilca_bdf_solver_m
    use kilca_types_m
    use iso_fortran_env, only: real64
    implicit none
    
    private
    
    ! Public types
    public :: bdf_solver_t
    
    ! Public procedures
    public :: bdf_solver_create, bdf_solver_destroy
    public :: bdf_solver_init, bdf_solver_reset
    public :: bdf_solver_step, bdf_solver_step_adaptive
    public :: bdf_compute_coefficients
    public :: bdf_newton_solve
    public :: bdf_error_estimate
    
    !> BDF solver type with state information
    type :: bdf_solver_t
        integer :: order = 2                        !< Current BDF order (1-5)
        integer :: max_order = 5                    !< Maximum allowed order
        integer :: n_steps = 0                      !< Number of steps taken
        integer :: n_fails = 0                      !< Number of failed steps
        real(real64) :: rtol = 1.0e-6_dp               !< Relative tolerance
        real(real64) :: atol = 1.0e-12_dp              !< Absolute tolerance
        real(real64) :: h_current = 0.0_dp             !< Current step size
        real(real64) :: h_next = 0.0_dp                !< Next step size
        real(real64) :: t_current = 0.0_dp             !< Current time
        real(real64), allocatable :: alpha(:)          !< BDF alpha coefficients
        real(real64) :: beta = 0.0_dp                  !< BDF beta coefficient
        real(real64), allocatable :: gamma(:)          !< Error estimation coefficients
        real(real64), allocatable :: y_history(:,:)    !< Solution history [neq, order+1]
        real(real64), allocatable :: h_history(:)      !< Step size history [order]
        real(real64), allocatable :: work(:)           !< Work array
        logical :: initialized = .false.           !< Initialization flag
        integer :: neq = 0                         !< Number of equations
        
        ! Newton solver parameters
        integer :: max_newton_iter = 10            !< Max Newton iterations
        real(real64) :: newton_tol = 1.0e-10_dp       !< Newton convergence tolerance
        
        ! Statistics
        integer :: n_jac_evals = 0                !< Jacobian evaluations
        integer :: n_lin_solves = 0               !< Linear system solves
        integer :: n_newton_iters = 0             !< Total Newton iterations
    end type bdf_solver_t
    
    ! BDF coefficients for orders 1-5
    real(real64), parameter :: BDF_ALPHA(6,5) = reshape([ &
        ! Order 1 (Backward Euler)
        1.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        ! Order 2
        3.0_dp/2.0_dp, -2.0_dp, 1.0_dp/2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        ! Order 3
        11.0_dp/6.0_dp, -3.0_dp, 3.0_dp/2.0_dp, -1.0_dp/3.0_dp, 0.0_dp, 0.0_dp, &
        ! Order 4
        25.0_dp/12.0_dp, -4.0_dp, 3.0_dp, -4.0_dp/3.0_dp, 1.0_dp/4.0_dp, 0.0_dp, &
        ! Order 5
        137.0_dp/60.0_dp, -5.0_dp, 5.0_dp, -10.0_dp/3.0_dp, 5.0_dp/4.0_dp, -1.0_dp/5.0_dp &
        ], [6, 5])
    
    real(real64), parameter :: BDF_BETA(5) = [ &
        1.0_dp,           & ! Order 1
        2.0_dp/3.0_dp,    & ! Order 2
        6.0_dp/11.0_dp,   & ! Order 3
        12.0_dp/25.0_dp,  & ! Order 4
        60.0_dp/137.0_dp  & ! Order 5
        ]
    
    ! Error coefficients for order selection
    real(real64), parameter :: ERROR_CONST(5) = [ &
        1.0_dp/2.0_dp,    & ! Order 1
        2.0_dp/9.0_dp,    & ! Order 2
        3.0_dp/22.0_dp,   & ! Order 3
        12.0_dp/125.0_dp, & ! Order 4
        10.0_dp/137.0_dp  & ! Order 5
        ]
    
    ! Interface for RHS function
    abstract interface
        subroutine rhs_func_interface(neq, t, y, ydot)
            import :: real64
            integer, intent(in) :: neq
            real(real64), intent(in) :: t
            real(real64), intent(in) :: y(*)
            real(real64), intent(out) :: ydot(*)
        end subroutine rhs_func_interface
        
        subroutine bdf_jac_func_interface(neq, t, y, jac, ldj)
            import :: real64
            integer, intent(in) :: neq, ldj
            real(real64), intent(in) :: t
            real(real64), intent(in) :: y(*)
            real(real64), intent(out) :: jac(ldj, *)
        end subroutine bdf_jac_func_interface
    end interface
    
contains
    
    !---------------------------------------------------------------------------
    ! BDF solver management
    !---------------------------------------------------------------------------
    
    !> Create and allocate BDF solver
    subroutine bdf_solver_create(solver, neq, order, ierr)
        type(bdf_solver_t), intent(out) :: solver
        integer, intent(in) :: neq
        integer, intent(in) :: order
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (order < 1 .or. order > 5) then
            ierr = -1
            return
        end if
        
        if (neq < 1) then
            ierr = -2
            return
        end if
        
        solver%neq = neq
        solver%order = order
        solver%max_order = min(5, max(order, 2))
        
        ! Allocate arrays
        allocate(solver%alpha(0:solver%max_order), stat=ierr)
        if (ierr /= 0) return
        
        allocate(solver%gamma(0:solver%max_order), stat=ierr)
        if (ierr /= 0) return
        
        allocate(solver%y_history(neq, 0:solver%max_order), stat=ierr)
        if (ierr /= 0) return
        
        allocate(solver%h_history(solver%max_order), stat=ierr)
        if (ierr /= 0) return
        
        allocate(solver%work(neq * 3), stat=ierr)
        if (ierr /= 0) return
        
        solver%initialized = .false.
        solver%n_steps = 0
        solver%n_fails = 0
        
    end subroutine bdf_solver_create
    
    !> Destroy BDF solver and deallocate memory
    subroutine bdf_solver_destroy(solver, ierr)
        type(bdf_solver_t), intent(inout) :: solver
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (allocated(solver%alpha)) deallocate(solver%alpha)
        if (allocated(solver%gamma)) deallocate(solver%gamma)
        if (allocated(solver%y_history)) deallocate(solver%y_history)
        if (allocated(solver%h_history)) deallocate(solver%h_history)
        if (allocated(solver%work)) deallocate(solver%work)
        
        solver%initialized = .false.
        
    end subroutine bdf_solver_destroy
    
    !> Initialize BDF solver with initial conditions
    subroutine bdf_solver_init(solver, t0, y0, h0, ierr)
        type(bdf_solver_t), intent(inout) :: solver
        real(real64), intent(in) :: t0
        real(real64), intent(in) :: y0(:)
        real(real64), intent(in) :: h0
        integer, intent(out) :: ierr
        integer :: i
        
        ierr = 0
        
        if (size(y0) /= solver%neq) then
            ierr = -1
            return
        end if
        
        ! Store initial condition
        solver%t_current = t0
        solver%h_current = h0
        solver%h_next = h0
        solver%y_history(:, 0) = y0
        
        ! Initialize history with initial values
        do i = 1, solver%max_order
            solver%y_history(:, i) = y0
            solver%h_history(i) = h0
        end do
        
        ! Start with order 1
        solver%order = 1
        
        ! Compute initial BDF coefficients
        call bdf_compute_coefficients(solver, ierr)
        
        solver%initialized = .true.
        solver%n_steps = 0
        solver%n_fails = 0
        
    end subroutine bdf_solver_init
    
    !> Reset solver statistics
    subroutine bdf_solver_reset(solver)
        type(bdf_solver_t), intent(inout) :: solver
        
        solver%n_jac_evals = 0
        solver%n_lin_solves = 0
        solver%n_newton_iters = 0
        solver%n_steps = 0
        solver%n_fails = 0
        
    end subroutine bdf_solver_reset
    
    !---------------------------------------------------------------------------
    ! BDF coefficients
    !---------------------------------------------------------------------------
    
    !> Compute BDF coefficients for current order and step sizes
    subroutine bdf_compute_coefficients(solver, ierr)
        type(bdf_solver_t), intent(inout) :: solver
        integer, intent(out) :: ierr
        integer :: i, k
        real(real64) :: xi(0:5), c(0:5)
        
        ierr = 0
        
        if (solver%order < 1 .or. solver%order > 5) then
            ierr = -1
            return
        end if
        
        ! For constant step size, use standard coefficients
        if (solver%n_steps < solver%order) then
            ! Use standard coefficients
            solver%alpha(0:solver%order) = BDF_ALPHA(1:solver%order+1, solver%order)
            solver%beta = BDF_BETA(solver%order)
        else
            ! Variable step size - compute coefficients
            ! This is a simplified version; full implementation would use
            ! divided differences or Lagrange interpolation
            
            ! Build xi array (step size ratios)
            xi(0) = 0.0_dp
            do k = 1, solver%order
                if (k <= solver%n_steps) then
                    xi(k) = xi(k-1) + solver%h_history(k) / solver%h_current
                else
                    xi(k) = xi(k-1) + 1.0_dp
                end if
            end do
            
            ! Compute coefficients (simplified for constant step)
            do i = 0, solver%order
                solver%alpha(i) = BDF_ALPHA(i+1, solver%order)
            end do
            solver%beta = BDF_BETA(solver%order)
        end if
        
        ! Set error coefficients
        solver%gamma(0:solver%order) = 0.0_dp
        solver%gamma(solver%order) = ERROR_CONST(solver%order)
        
    end subroutine bdf_compute_coefficients
    
    !---------------------------------------------------------------------------
    ! BDF step
    !---------------------------------------------------------------------------
    
    !> Take one BDF step
    subroutine bdf_solver_step(solver, f, jac, t_end, y_new, ierr)
        type(bdf_solver_t), intent(inout) :: solver
        procedure(rhs_func_interface) :: f
        procedure(bdf_jac_func_interface) :: jac
        real(real64), intent(in) :: t_end
        real(real64), intent(out) :: y_new(:)
        integer, intent(out) :: ierr
        real(real64) :: h, t_new
        integer :: i
        
        ierr = 0
        
        if (.not. solver%initialized) then
            ierr = -1
            return
        end if
        
        ! Determine step size
        h = min(solver%h_current, t_end - solver%t_current)
        t_new = solver%t_current + h
        
        ! Shift history
        do i = solver%order, 1, -1
            solver%y_history(:, i) = solver%y_history(:, i-1)
        end do
        
        ! Predict new solution using polynomial extrapolation
        call bdf_predict(solver, y_new)
        
        ! Solve nonlinear system using Newton's method
        call bdf_newton_solve_internal(solver, f, jac, t_new, h, y_new, ierr)
        
        if (ierr /= 0) then
            ! Step failed - restore history
            do i = 1, solver%order
                solver%y_history(:, i) = solver%y_history(:, i-1)
            end do
            solver%n_fails = solver%n_fails + 1
            return
        end if
        
        ! Store new solution
        solver%y_history(:, 0) = y_new
        solver%t_current = t_new
        solver%h_history(1) = h
        
        ! Shift step size history
        do i = solver%max_order, 2, -1
            solver%h_history(i) = solver%h_history(i-1)
        end do
        
        solver%n_steps = solver%n_steps + 1
        
    end subroutine bdf_solver_step
    
    !> Take adaptive BDF step with error control
    subroutine bdf_solver_step_adaptive(solver, f, jac, t_end, y_new, err_est, ierr)
        type(bdf_solver_t), intent(inout) :: solver
        procedure(rhs_func_interface) :: f
        procedure(bdf_jac_func_interface) :: jac
        real(real64), intent(in) :: t_end
        real(real64), intent(out) :: y_new(:)
        real(real64), intent(out) :: err_est
        integer, intent(out) :: ierr
        real(real64) :: h_new, safety_factor
        integer :: new_order
        
        ! Take step
        call bdf_solver_step(solver, f, jac, t_end, y_new, ierr)
        
        if (ierr /= 0) then
            ! Reduce step size on failure
            solver%h_current = solver%h_current * 0.5_dp
            return
        end if
        
        ! Estimate error
        call bdf_error_estimate(solver, y_new, err_est)
        
        ! Compute new step size
        safety_factor = 0.9_dp
        if (err_est > 0.0_dp) then
            h_new = solver%h_current * safety_factor * &
                    (1.0_dp / err_est) ** (1.0_dp / real(solver%order + 1, dp))
        else
            h_new = solver%h_current * 2.0_dp
        end if
        
        ! Limit step size change
        h_new = max(0.1_dp * solver%h_current, &
                    min(10.0_dp * solver%h_current, h_new))
        
        solver%h_next = h_new
        
        ! Consider order change (simplified)
        if (solver%n_steps > solver%order + 1) then
            if (err_est < 0.1_dp .and. solver%order < solver%max_order) then
                ! Increase order
                new_order = solver%order + 1
            else if (err_est > 10.0_dp .and. solver%order > 1) then
                ! Decrease order
                new_order = solver%order - 1
            else
                new_order = solver%order
            end if
            
            if (new_order /= solver%order) then
                solver%order = new_order
                call bdf_compute_coefficients(solver, ierr)
            end if
        end if
        
        solver%h_current = h_new
        
    end subroutine bdf_solver_step_adaptive
    
    !---------------------------------------------------------------------------
    ! Newton solver for implicit BDF equation
    !---------------------------------------------------------------------------
    
    !> Solve BDF nonlinear system using Newton's method
    subroutine bdf_newton_solve(solver, f, jac, t_new, h, y, ierr)
        type(bdf_solver_t), intent(inout) :: solver
        procedure(rhs_func_interface) :: f
        procedure(bdf_jac_func_interface) :: jac
        real(real64), intent(in) :: t_new, h
        real(real64), intent(inout) :: y(:)
        integer, intent(out) :: ierr
        real(real64), allocatable :: jacobian(:,:), rhs(:), dy(:)
        real(real64), allocatable :: ydot(:), y_pred(:)
        real(real64) :: conv_rate, conv_test, old_conv_test
        integer :: iter, i, j, info
        integer, allocatable :: ipiv(:)
        
        ierr = 0
        
        ! Allocate work arrays
        allocate(jacobian(solver%neq, solver%neq))
        allocate(rhs(solver%neq))
        allocate(dy(solver%neq))
        allocate(ydot(solver%neq))
        allocate(y_pred(solver%neq))
        allocate(ipiv(solver%neq))
        
        ! Save predicted value
        y_pred = y
        
        ! Newton iteration
        conv_rate = 1.0_dp
        old_conv_test = 1.0e10_dp
        
        do iter = 1, solver%max_newton_iter
            ! Evaluate f(t_new, y)
            call f(solver%neq, t_new, y, ydot)
            
            ! Build residual: r = sum(alpha_i * y_{n-i}) - h * beta * f(t_n, y_n)
            rhs = 0.0_dp
            do i = 0, solver%order
                rhs = rhs + solver%alpha(i) * solver%y_history(:, i)
            end do
            rhs = rhs - h * solver%beta * ydot
            
            ! Check convergence
            conv_test = maxval(abs(rhs) / (solver%atol + solver%rtol * abs(y)))
            
            if (conv_test < solver%newton_tol) then
                ! Converged
                exit
            end if
            
            ! Check for divergence
            if (iter > 1) then
                conv_rate = conv_test / old_conv_test
                if (conv_rate > 0.9_dp) then
                    ! Not converging
                    ierr = -1
                    exit
                end if
            end if
            old_conv_test = conv_test
            
            ! Evaluate Jacobian
            call jac(solver%neq, t_new, y, jacobian, solver%neq)
            solver%n_jac_evals = solver%n_jac_evals + 1
            
            ! Form Newton matrix: M = alpha_0 * I - h * beta * J
            do j = 1, solver%neq
                do i = 1, solver%neq
                    if (i == j) then
                        jacobian(i, j) = solver%alpha(0) - h * solver%beta * jacobian(i, j)
                    else
                        jacobian(i, j) = -h * solver%beta * jacobian(i, j)
                    end if
                end do
            end do
            
            ! Solve M * dy = -r
            dy = -rhs
            call dgesv(solver%neq, 1, jacobian, solver%neq, ipiv, dy, solver%neq, info)
            solver%n_lin_solves = solver%n_lin_solves + 1
            
            if (info /= 0) then
                ierr = -2
                exit
            end if
            
            ! Update solution
            y = y + dy
            
            solver%n_newton_iters = solver%n_newton_iters + 1
        end do
        
        if (iter > solver%max_newton_iter) then
            ierr = -3
        end if
        
        ! Clean up
        deallocate(jacobian, rhs, dy, ydot, y_pred, ipiv)
        
    end subroutine bdf_newton_solve
    
    !---------------------------------------------------------------------------
    ! Helper procedures
    !---------------------------------------------------------------------------
    
    !> Predict solution using polynomial extrapolation
    subroutine bdf_predict(solver, y_pred)
        type(bdf_solver_t), intent(in) :: solver
        real(real64), intent(out) :: y_pred(:)
        integer :: i
        
        ! Simple prediction: extrapolate from history
        ! For order 1: y_pred = y_{n-1}
        ! For order 2: y_pred = 2*y_{n-1} - y_{n-2}
        ! etc.
        
        y_pred = 0.0_dp
        
        select case(solver%order)
        case(1)
            y_pred = solver%y_history(:, 1)
        case(2)
            y_pred = 2.0_dp * solver%y_history(:, 1) - solver%y_history(:, 2)
        case(3)
            y_pred = 3.0_dp * solver%y_history(:, 1) - &
                     3.0_dp * solver%y_history(:, 2) + &
                     solver%y_history(:, 3)
        case default
            ! General formula
            y_pred = solver%y_history(:, 1)
        end select
        
    end subroutine bdf_predict
    
    !> Estimate local truncation error
    subroutine bdf_error_estimate(solver, y_new, err_est)
        type(bdf_solver_t), intent(in) :: solver
        real(real64), intent(in) :: y_new(:)
        real(real64), intent(out) :: err_est
        real(real64), allocatable :: err_vec(:)
        integer :: i
        
        allocate(err_vec(solver%neq))
        
        ! Estimate error using difference between predicted and corrected values
        ! This is a simplified error estimate
        err_vec = 0.0_dp
        
        ! Use Milne's error estimate for BDF methods
        do i = 0, solver%order
            err_vec = err_vec + solver%gamma(i) * solver%y_history(:, i)
        end do
        
        ! Normalize error
        err_vec = abs(err_vec) / (solver%atol + solver%rtol * abs(y_new))
        
        ! RMS norm
        err_est = sqrt(sum(err_vec**2) / real(solver%neq, dp))
        
        deallocate(err_vec)
        
    end subroutine bdf_error_estimate
    
    !> Internal wrapper to avoid interface issues
    subroutine bdf_newton_solve_internal(solver, f, jac, t_new, h, y, ierr)
        type(bdf_solver_t), intent(inout) :: solver
        procedure(rhs_func_interface) :: f
        procedure(bdf_jac_func_interface) :: jac
        real(real64), intent(in) :: t_new, h
        real(real64), intent(inout) :: y(:)
        integer, intent(out) :: ierr
        
        ! Just call the original function
        call bdf_newton_solve(solver, f, jac, t_new, h, y, ierr)
        
    end subroutine bdf_newton_solve_internal
    
end module kilca_bdf_solver_m