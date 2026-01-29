module quadpack_integration_m
    ! QUADPACK integration wrapper for KIM Fokker-Planck kernels
    ! Provides thread-safe context passing and unified interface

    use KIM_kinds_m, only: dp
    use integrands_rkf45_m, only: rkf45_integrand_context_t

    implicit none
    private

    ! Public procedures
    public :: quadpack_integrate_F1
    public :: quadpack_integrate_F2
    public :: quadpack_integrate_F3
    public :: init_quadpack_module
    public :: cleanup_quadpack_module

    ! QUADPACK algorithm selection (finite-interval theta integrals)
    integer, parameter, public :: QUADPACK_QAG = 1
    integer, parameter, public :: QUADPACK_QAGS = 2

    ! QUADPACK Gauss-Kronrod rule keys
    integer, parameter, public :: QK15 = 1  ! 7-15 points
    integer, parameter, public :: QK21 = 2  ! 10-21 points
    integer, parameter, public :: QK31 = 3  ! 15-31 points
    integer, parameter, public :: QK41 = 4  ! 20-41 points
    integer, parameter, public :: QK51 = 5  ! 25-51 points
    integer, parameter, public :: QK61 = 6  ! 30-61 points

    ! Thread-local context storage (single threadprivate scalar)
    type(rkf45_integrand_context_t), save :: current_context
    !$omp threadprivate(current_context)

    ! Module initialization flag
    logical :: module_initialized = .false.

    ! Per-thread reusable QUADPACK workspaces to avoid heavy per-call alloc/dealloc
    integer, allocatable, save :: iwork_tl(:)
    real(dp), allocatable, save :: work_tl(:)
    integer, save :: work_limit_tl = 0
    integer, save :: work_lenw_tl  = 0
    !$omp threadprivate(iwork_tl, work_tl, work_limit_tl, work_lenw_tl)

    contains

    subroutine ensure_workspace(limit)
        ! Ensure per-thread QUADPACK workspace has at least the requested capacity
        implicit none
        integer, intent(in) :: limit
        integer :: need_lenw

        need_lenw = 4 * max(limit, 1)

        if (work_limit_tl < limit) then
            if (allocated(iwork_tl)) deallocate(iwork_tl)
            allocate(iwork_tl(limit))
            work_limit_tl = limit
        end if

        if (work_lenw_tl < need_lenw) then
            if (allocated(work_tl)) deallocate(work_tl)
            allocate(work_tl(need_lenw))
            work_lenw_tl = need_lenw
        end if
    end subroutine ensure_workspace

    subroutine init_quadpack_module()
        ! Initialize QUADPACK module and thread-local storage
        implicit none
        if (module_initialized) return
        module_initialized = .true.
    end subroutine init_quadpack_module

    subroutine cleanup_quadpack_module()
        ! Cleanup module resources
        implicit none

        module_initialized = .false.

        ! Deallocate any thread-local workspaces
        !$omp parallel
        if (allocated(iwork_tl)) deallocate(iwork_tl)
        if (allocated(work_tl))  deallocate(work_tl)
        work_limit_tl = 0
        work_lenw_tl  = 0
        !$omp end parallel

    end subroutine cleanup_quadpack_module

    subroutine set_thread_context(context)
        ! Store context for current thread
        implicit none
        type(rkf45_integrand_context_t), intent(in) :: context
        current_context = context
    end subroutine set_thread_context

    function get_thread_context() result(context)
        ! Retrieve context for current thread
        implicit none
        type(rkf45_integrand_context_t) :: context
        context = current_context
    end function get_thread_context

    ! Wrapper functions for QUADPACK - these match QUADPACK's expected signature
    function quadpack_wrapper_F1(theta) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F1
        implicit none
        real(dp), intent(in) :: theta
        real(dp) :: val
        type(rkf45_integrand_context_t) :: context

        ! Retrieve context for this thread
        context = get_thread_context()

        ! Call the actual integrand with context
        val = rkf45_integrand_F1(theta, context)

    end function quadpack_wrapper_F1

    function quadpack_wrapper_F2(theta) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F2
        implicit none
        real(dp), intent(in) :: theta
        real(dp) :: val
        type(rkf45_integrand_context_t) :: context

        context = get_thread_context()
        val = rkf45_integrand_F2(theta, context)

    end function quadpack_wrapper_F2

    function quadpack_wrapper_F3(theta) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F3
        implicit none
        real(dp), intent(in) :: theta
        real(dp) :: val
        type(rkf45_integrand_context_t) :: context

        context = get_thread_context()
        val = rkf45_integrand_F3(theta, context)

    end function quadpack_wrapper_F3

    ! U-substitution wrapper functions (u = sin(theta/2))
    function quadpack_wrapper_F1_u(u) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F1
        implicit none
        real(dp), intent(in) :: u
        real(dp) :: val
        real(dp) :: theta, jac, uu
        type(rkf45_integrand_context_t) :: context

        context = get_thread_context()

        ! Transform from u to theta with Jacobian
        uu = min(max(u, -1.0d0), 1.0d0)  ! Ensure u in [-1, 1]
        theta = 2.0d0 * asin(uu)
        jac = 2.0d0 / max(sqrt(1.0d0 - uu*uu), 1.0d-300)

        val = rkf45_integrand_F1(theta, context) * jac

    end function quadpack_wrapper_F1_u

    function quadpack_wrapper_F2_u(u) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F2
        implicit none
        real(dp), intent(in) :: u
        real(dp) :: val
        real(dp) :: theta, jac, uu
        type(rkf45_integrand_context_t) :: context

        context = get_thread_context()

        uu = min(max(u, -1.0d0), 1.0d0)
        theta = 2.0d0 * asin(uu)
        jac = 2.0d0 / max(sqrt(1.0d0 - uu*uu), 1.0d-300)

        val = rkf45_integrand_F2(theta, context) * jac

    end function quadpack_wrapper_F2_u

    function quadpack_wrapper_F3_u(u) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F3
        implicit none
        real(dp), intent(in) :: u
        real(dp) :: val
        real(dp) :: theta, jac, uu
        type(rkf45_integrand_context_t) :: context

        context = get_thread_context()

        uu = min(max(u, -1.0d0), 1.0d0)
        theta = 2.0d0 * asin(uu)
        jac = 2.0d0 / max(sqrt(1.0d0 - uu*uu), 1.0d-300)

        val = rkf45_integrand_F3(theta, context) * jac

    end function quadpack_wrapper_F3_u

    subroutine quadpack_integrate_F1(result, epsabs, epsrel, context, &
                                     theta_min, theta_max, use_u_substitution)
        ! Integrate F1 using QUADPACK
        use grid_m, only: quadpack_key, quadpack_limit, quadpack_algorithm
        use config_m, only: fdebug
        implicit none

        real(dp), intent(out) :: result
        real(dp), intent(in) :: epsabs, epsrel
        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp), intent(in) :: theta_min, theta_max
        logical, intent(in) :: use_u_substitution

        ! QUADPACK variables
        real(dp) :: abserr
        integer :: neval, ier, last
        integer :: limit, lenw
        real(dp) :: a, b
        integer :: key

        ! External QUADPACK entry points
        external :: dqag, dqags, dqagi

        ! Get thread ID and set context
        call set_thread_context(context)

        ! Set integration parameters
        key = quadpack_key
        if (key < 1 .or. key > 6) key = QK61  ! Default to highest order

        limit = quadpack_limit
        if (limit < 1) limit = 500

        lenw = 4 * limit
        call ensure_workspace(limit)

        select case (trim(quadpack_algorithm))
        case('QAG')
            if (use_u_substitution) then
                a = sin(0.5d0 * theta_min)
                b = sin(0.5d0 * theta_max)
            call dqag(quadpack_wrapper_F1_u, a, b, epsabs, epsrel, key, &
                     result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            else
                a = theta_min
                b = theta_max
            call dqag(quadpack_wrapper_F1, a, b, epsabs, epsrel, key, &
                     result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            end if
        case('QAGS')
            ! Endpoint-singular integrands on finite intervals
            if (use_u_substitution) then
                a = sin(0.5d0 * theta_min)
                b = sin(0.5d0 * theta_max)
                call dqags(quadpack_wrapper_F1_u, a, b, epsabs, epsrel, &
                           result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            else
                a = theta_min
                b = theta_max
                call dqags(quadpack_wrapper_F1, a, b, epsabs, epsrel, &
                           result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            end if
        case default
            print *, 'QUADPACK Error: Unsupported algorithm for finite-interval theta integrals: ', trim(quadpack_algorithm)
            print *, 'Supported: QAG, QAGS'
            error stop 'Invalid QUADPACK algorithm'
        end select

        ! Check for errors (throttle prints to higher debug levels)
        if (ier /= 0) then
            if (ier == 6) then
                print *, "QUADPACK Error: Invalid input parameters"
                error stop "QUADPACK integration failed"
            else
                if (fdebug >= 2) then
                    select case(ier)
                        case(1)
                            print *, "QUADPACK Warning: Maximum subdivisions reached"
                        case(2)
                            print *, "QUADPACK Warning: Roundoff error detected"
                        case(3)
                            print *, "QUADPACK Warning: Extremely bad integrand behavior"
                        case(4)
                            print *, "QUADPACK Warning: Roundoff error in extrapolation"
                        case(5)
                            print *, "QUADPACK Warning: Divergent integral"
                    end select
                end if
            end if
        end if

        ! Workspaces are reused per-thread; do not deallocate here

    end subroutine quadpack_integrate_F1

    subroutine quadpack_integrate_F2(result, epsabs, epsrel, context, &
                                     theta_min, theta_max, use_u_substitution)
        ! Integrate F2 using QUADPACK
        use grid_m, only: quadpack_key, quadpack_limit, quadpack_algorithm
        use config_m, only: fdebug
        implicit none

        real(dp), intent(out) :: result
        real(dp), intent(in) :: epsabs, epsrel
        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp), intent(in) :: theta_min, theta_max
        logical, intent(in) :: use_u_substitution

        ! QUADPACK variables
        real(dp) :: abserr
        integer :: neval, ier, last
        integer :: limit, lenw
        real(dp) :: a, b
        integer :: key

        ! External QUADPACK entry points
        external :: dqag, dqags, dqagi

        ! Get thread ID and set context
        call set_thread_context(context)

        ! Set integration parameters
        key = quadpack_key
        if (key < 1 .or. key > 6) key = QK61

        limit = quadpack_limit
        if (limit < 1) limit = 500

        lenw = 4 * limit
        call ensure_workspace(limit)

        select case (trim(quadpack_algorithm))
        case('QAG')
            if (use_u_substitution) then
                a = sin(0.5d0 * theta_min)
                b = sin(0.5d0 * theta_max)
            call dqag(quadpack_wrapper_F2_u, a, b, epsabs, epsrel, key, &
                     result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            else
                a = theta_min
                b = theta_max
            call dqag(quadpack_wrapper_F2, a, b, epsabs, epsrel, key, &
                     result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            end if
        case('QAGS')
            if (use_u_substitution) then
                a = sin(0.5d0 * theta_min)
                b = sin(0.5d0 * theta_max)
                call dqags(quadpack_wrapper_F2_u, a, b, epsabs, epsrel, &
                           result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            else
                a = theta_min
                b = theta_max
                call dqags(quadpack_wrapper_F2, a, b, epsabs, epsrel, &
                           result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            end if
        case default
            print *, 'QUADPACK Error: Unsupported algorithm for finite-interval theta integrals: ', trim(quadpack_algorithm)
            print *, 'Supported: QAG, QAGS'
            error stop 'Invalid QUADPACK algorithm'
        end select

        ! Check for errors (same as F1)
        if (ier /= 0) then
            if (ier == 6) then
                print *, "QUADPACK Error: Invalid input parameters"
                error stop "QUADPACK integration failed"
            else
                if (fdebug >= 2) then
                    select case(ier)
                        case(1)
                            print *, "QUADPACK Warning: Maximum subdivisions reached"
                        case(2)
                            print *, "QUADPACK Warning: Roundoff error detected"
                        case(3)
                            print *, "QUADPACK Warning: Extremely bad integrand behavior"
                        case(4)
                            print *, "QUADPACK Warning: Roundoff error in extrapolation"
                        case(5)
                            print *, "QUADPACK Warning: Divergent integral"
                    end select
                end if
            end if
        end if

        ! Workspaces are reused per-thread; do not deallocate here

    end subroutine quadpack_integrate_F2

    subroutine quadpack_integrate_F3(result, epsabs, epsrel, context, &
                                     theta_min, theta_max, use_u_substitution)
        ! Integrate F3 using QUADPACK
        use grid_m, only: quadpack_key, quadpack_limit, quadpack_algorithm
        use config_m, only: fdebug
        implicit none

        real(dp), intent(out) :: result
        real(dp), intent(in) :: epsabs, epsrel
        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp), intent(in) :: theta_min, theta_max
        logical, intent(in) :: use_u_substitution

        ! QUADPACK variables
        real(dp) :: abserr
        integer :: neval, ier, last
        integer :: limit, lenw
        real(dp) :: a, b
        integer :: key

        ! External QUADPACK entry points
        external :: dqag, dqags, dqagi

        ! Get thread ID and set context
        call set_thread_context(context)

        ! Set integration parameters
        key = quadpack_key
        if (key < 1 .or. key > 6) key = QK61

        limit = quadpack_limit
        if (limit < 1) limit = 500

        lenw = 4 * limit
        call ensure_workspace(limit)

        select case (trim(quadpack_algorithm))
        case('QAG')
            if (use_u_substitution) then
                a = sin(0.5d0 * theta_min)
                b = sin(0.5d0 * theta_max)
            call dqag(quadpack_wrapper_F3_u, a, b, epsabs, epsrel, key, &
                     result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            else
                a = theta_min
                b = theta_max
            call dqag(quadpack_wrapper_F3, a, b, epsabs, epsrel, key, &
                     result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            end if
        case('QAGS')
            if (use_u_substitution) then
                a = sin(0.5d0 * theta_min)
                b = sin(0.5d0 * theta_max)
                call dqags(quadpack_wrapper_F3_u, a, b, epsabs, epsrel, &
                         result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            else
                a = theta_min
                b = theta_max
                call dqags(quadpack_wrapper_F3, a, b, epsabs, epsrel, &
                         result, abserr, neval, ier, limit, lenw, last, iwork_tl, work_tl)
            end if
        case default
            print *, 'QUADPACK Error: Unsupported algorithm for finite-interval theta integrals: ', trim(quadpack_algorithm)
            print *, 'Supported: QAG, QAGS'
            error stop 'Invalid QUADPACK algorithm'
        end select

        ! Check for errors (same as F1)
        if (ier /= 0) then
            if (ier == 6) then
                print *, "QUADPACK Error: Invalid input parameters"
                error stop "QUADPACK integration failed"
            else
                if (fdebug >= 2) then
                    select case(ier)
                        case(1)
                            print *, "QUADPACK Warning: Maximum subdivisions reached"
                        case(2)
                            print *, "QUADPACK Warning: Roundoff error detected"
                        case(3)
                            print *, "QUADPACK Warning: Extremely bad integrand behavior"
                        case(4)
                            print *, "QUADPACK Warning: Roundoff error in extrapolation"
                        case(5)
                            print *, "QUADPACK Warning: Divergent integral"
                    end select
                end if
            end if
        end if

        ! Workspaces are reused per-thread; do not deallocate here

    end subroutine quadpack_integrate_F3

end module quadpack_integration_m
