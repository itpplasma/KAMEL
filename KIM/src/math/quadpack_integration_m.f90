module quadpack_integration_m
    ! Adaptive Gauss-Kronrod integration wrapper for KIM Fokker-Planck kernels.
    ! Backed by the fortnum numerical core (fortnum_integrate); the
    ! rkf45_integrand_context_t flows through the integrand callback's class(*)
    ! ctx, so no module-scope thread context is needed.

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

    ! Gauss-Kronrod rule keys (config selector values)
    integer, parameter, public :: QK15 = 1  ! 7-15 points
    integer, parameter, public :: QK21 = 2  ! 10-21 points
    integer, parameter, public :: QK31 = 3  ! 15-31 points
    integer, parameter, public :: QK41 = 4  ! 20-41 points
    integer, parameter, public :: QK51 = 5  ! 25-51 points
    integer, parameter, public :: QK61 = 6  ! 30-61 points

    ! Module initialization flag
    logical :: module_initialized = .false.

    contains

    subroutine init_quadpack_module()
        ! Initialize QUADPACK module
        implicit none
        if (module_initialized) return
        module_initialized = .true.
    end subroutine init_quadpack_module

    subroutine cleanup_quadpack_module()
        ! Cleanup module resources
        implicit none
        module_initialized = .false.
    end subroutine cleanup_quadpack_module

    integer function fortnum_gk_key(config_key) result(key)
        ! Map the config rule selector (QK15..QK61) to a fortnum GK rule key.
        ! fortnum supports the GK rules with 15, 21, 31, and 61 points; the
        ! 41/51 selectors round up to the next available higher-order rule.
        implicit none
        integer, intent(in) :: config_key
        select case (config_key)
        case (QK15)
            key = 15
        case (QK21)
            key = 21
        case (QK31)
            key = 31
        case (QK41, QK51, QK61)
            key = 61
        case default
            key = 61
        end select
    end function fortnum_gk_key

    ! Integrand callbacks on the fortnum integrate_integrand_t ABI. The
    ! rkf45_integrand_context_t is forwarded unchanged through ctx; the _u
    ! variants apply the u = sin(theta/2) substitution with its Jacobian.
    function quadpack_wrapper_F1(theta, ctx) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F1
        implicit none
        real(dp), intent(in) :: theta
        class(*), intent(in), optional :: ctx
        real(dp) :: val
        val = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (rkf45_integrand_context_t)
                val = rkf45_integrand_F1(theta, ctx)
            end select
        end if
    end function quadpack_wrapper_F1

    function quadpack_wrapper_F2(theta, ctx) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F2
        implicit none
        real(dp), intent(in) :: theta
        class(*), intent(in), optional :: ctx
        real(dp) :: val
        val = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (rkf45_integrand_context_t)
                val = rkf45_integrand_F2(theta, ctx)
            end select
        end if
    end function quadpack_wrapper_F2

    function quadpack_wrapper_F3(theta, ctx) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F3
        implicit none
        real(dp), intent(in) :: theta
        class(*), intent(in), optional :: ctx
        real(dp) :: val
        val = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (rkf45_integrand_context_t)
                val = rkf45_integrand_F3(theta, ctx)
            end select
        end if
    end function quadpack_wrapper_F3

    function quadpack_wrapper_F1_u(u, ctx) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F1
        implicit none
        real(dp), intent(in) :: u
        class(*), intent(in), optional :: ctx
        real(dp) :: val
        real(dp) :: theta, jac, uu
        val = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (rkf45_integrand_context_t)
                uu = min(max(u, -1.0d0), 1.0d0)  ! Ensure u in [-1, 1]
                theta = 2.0d0 * asin(uu)
                jac = 2.0d0 / max(sqrt(1.0d0 - uu*uu), 1.0d-300)
                val = rkf45_integrand_F1(theta, ctx) * jac
            end select
        end if
    end function quadpack_wrapper_F1_u

    function quadpack_wrapper_F2_u(u, ctx) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F2
        implicit none
        real(dp), intent(in) :: u
        class(*), intent(in), optional :: ctx
        real(dp) :: val
        real(dp) :: theta, jac, uu
        val = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (rkf45_integrand_context_t)
                uu = min(max(u, -1.0d0), 1.0d0)
                theta = 2.0d0 * asin(uu)
                jac = 2.0d0 / max(sqrt(1.0d0 - uu*uu), 1.0d-300)
                val = rkf45_integrand_F2(theta, ctx) * jac
            end select
        end if
    end function quadpack_wrapper_F2_u

    function quadpack_wrapper_F3_u(u, ctx) result(val)
        use integrands_rkf45_m, only: rkf45_integrand_F3
        implicit none
        real(dp), intent(in) :: u
        class(*), intent(in), optional :: ctx
        real(dp) :: val
        real(dp) :: theta, jac, uu
        val = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (rkf45_integrand_context_t)
                uu = min(max(u, -1.0d0), 1.0d0)
                theta = 2.0d0 * asin(uu)
                jac = 2.0d0 / max(sqrt(1.0d0 - uu*uu), 1.0d-300)
                val = rkf45_integrand_F3(theta, ctx) * jac
            end select
        end if
    end function quadpack_wrapper_F3_u

    subroutine run_quadpack(f, fu, result, epsabs, epsrel, context, &
                            theta_min, theta_max, use_u_substitution)
        ! Shared driver: select QAG/QAGS via grid config, apply the optional
        ! u-substitution, and route to the matching fortnum adaptive routine.
        use fortnum_integrate, only: integrate_integrand_t, integrate_qag, &
            integrate_qags, integrate_workspace_t, integrate_epstab_t, &
            integrate_result_t
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        use grid_m, only: quadpack_key, quadpack_limit, quadpack_algorithm
        use logger_m, only: log_debug, log_error
        implicit none

        procedure(integrate_integrand_t) :: f, fu
        real(dp), intent(out) :: result
        real(dp), intent(in) :: epsabs, epsrel
        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp), intent(in) :: theta_min, theta_max
        logical, intent(in) :: use_u_substitution

        type(integrate_workspace_t) :: workspace
        type(integrate_epstab_t) :: epstab
        type(integrate_result_t) :: res
        type(fortnum_status_t) :: status
        real(dp) :: a, b
        integer :: key, limit

        key = fortnum_gk_key(quadpack_key)

        limit = quadpack_limit
        if (limit < 1) limit = 500

        if (use_u_substitution) then
            a = sin(0.5d0 * theta_min)
            b = sin(0.5d0 * theta_max)
        else
            a = theta_min
            b = theta_max
        end if

        select case (trim(quadpack_algorithm))
        case ('QAG')
            if (use_u_substitution) then
                call integrate_qag(fu, a, b, epsabs, epsrel, workspace, res, &
                                   status, key=key, limit=limit, ctx=context)
            else
                call integrate_qag(f, a, b, epsabs, epsrel, workspace, res, &
                                   status, key=key, limit=limit, ctx=context)
            end if
        case ('QAGS')
            if (use_u_substitution) then
                call integrate_qags(fu, a, b, epsabs, epsrel, workspace, &
                                    epstab, res, status, limit=limit, ctx=context)
            else
                call integrate_qags(f, a, b, epsabs, epsrel, workspace, &
                                    epstab, res, status, limit=limit, ctx=context)
            end if
        case default
            call log_error('QUADPACK: Unsupported algorithm: ' // &
                trim(quadpack_algorithm) // '. Supported: QAG, QAGS')
            result = 0.0_dp
            return
        end select

        result = res%value

        if (status%code /= FORTNUM_OK) then
            call log_debug('QUADPACK: integration did not fully converge')
        end if
    end subroutine run_quadpack

    subroutine quadpack_integrate_F1(result, epsabs, epsrel, context, &
                                     theta_min, theta_max, use_u_substitution)
        implicit none
        real(dp), intent(out) :: result
        real(dp), intent(in) :: epsabs, epsrel
        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp), intent(in) :: theta_min, theta_max
        logical, intent(in) :: use_u_substitution

        call run_quadpack(quadpack_wrapper_F1, quadpack_wrapper_F1_u, result, &
                          epsabs, epsrel, context, theta_min, theta_max, &
                          use_u_substitution)
    end subroutine quadpack_integrate_F1

    subroutine quadpack_integrate_F2(result, epsabs, epsrel, context, &
                                     theta_min, theta_max, use_u_substitution)
        implicit none
        real(dp), intent(out) :: result
        real(dp), intent(in) :: epsabs, epsrel
        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp), intent(in) :: theta_min, theta_max
        logical, intent(in) :: use_u_substitution

        call run_quadpack(quadpack_wrapper_F2, quadpack_wrapper_F2_u, result, &
                          epsabs, epsrel, context, theta_min, theta_max, &
                          use_u_substitution)
    end subroutine quadpack_integrate_F2

    subroutine quadpack_integrate_F3(result, epsabs, epsrel, context, &
                                     theta_min, theta_max, use_u_substitution)
        implicit none
        real(dp), intent(out) :: result
        real(dp), intent(in) :: epsabs, epsrel
        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp), intent(in) :: theta_min, theta_max
        logical, intent(in) :: use_u_substitution

        call run_quadpack(quadpack_wrapper_F3, quadpack_wrapper_F3_u, result, &
                          epsabs, epsrel, context, theta_min, theta_max, &
                          use_u_substitution)
    end subroutine quadpack_integrate_F3

end module quadpack_integration_m
