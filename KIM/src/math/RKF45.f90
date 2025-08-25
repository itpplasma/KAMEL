module RKF45_mod

    use KIM_kinds, only: dp

    implicit none

    real(dp), parameter :: c41 = 0.115740740740741d0
    real(dp), parameter :: c42 = 0.548927875243665d0
    real(dp), parameter :: c43 = 0.535722994391612d0
    real(dp), parameter :: c44 = 0.2d0

    real(dp), parameter :: c51 = 0.118518518518519d0
    real(dp), parameter :: c52 = 0.518986354775828d0
    real(dp), parameter :: c53 = 0.506131490342017d0
    real(dp), parameter :: c54 = 0.18d0
    real(dp), parameter :: c55 = 0.036363636363636d0

    real(dp), parameter :: one_over_cubic_root_2 = 0.840896415253715d0

    integer(dp), parameter :: MAX_INT_COUNT = 1000000

    contains

    subroutine RKF45_1D(f, y0, x0, xmax, h0, tol, sol)

        use KIM_kinds, only: dp
        implicit none

        ! Arguments
        real(dp), intent(in)  :: y0, x0, xmax, h0, tol
        real(dp), intent(out) :: sol

        ! Locals
        real(dp) :: yk, xk, hk
        real(dp) :: ytrial, xtrial, hnew
        integer  :: step_count

        interface
            function f(x)
                use KIM_kinds, only: dp
                implicit none
                real(dp), intent(in) :: x
                real(dp) :: f
            end function
        end interface

        ! Initialization
        yk = y0
        xk = x0
        hk = h0
        step_count = 0

        ! Integration loop
        do while (xk < xmax)

            if (step_count > MAX_INT_COUNT) then
                print *, "Integration error: too many steps (", MAX_INT_COUNT, ")"
                print *, "Current xk: ", xk, " Current yk: ", yk, " Current hk: ", hk
                print *, "xmax: ", xmax, " tol: ", tol, " h0: ", h0
                stop
            end if

            if (xk + hk > xmax) hk = xmax - xk

            call RKF45_step_1D(f, yk, xk, hk, tol, xtrial, ytrial, hnew)

            if (xtrial > xk) then
                xk = xtrial
                yk = ytrial
            end if

            hk = hnew
            step_count = step_count + 1

        end do

        sol = yk

    end subroutine RKF45_1D


    subroutine RKF45_step_1D(f, yk, xk, h, tol, xkp1, ykp1, hnew)
        ! Classical Runge-Kutta-Fehlberg 4(5) method (adaptive stepper)

        use KIM_kinds, only: dp
        implicit none

        real(dp), intent(in)  :: yk, xk, h, tol
        real(dp), intent(out) :: xkp1, ykp1, hnew

        real(dp) :: k1, k2, k3, k4, k5, k6
        real(dp) :: y4, y5, err, s

        interface
            function f(x)
                use KIM_kinds, only: dp
                implicit none
                real(dp), intent(in) :: x
                real(dp) :: f
            end function
        end interface

        k1 = h * f(xk)
        k2 = h * f(xk + 0.25d0*h)
        k3 = h * f(xk + 0.375d0*h)
        k4 = h * f(xk + 0.923076923076923d0*h)
        k5 = h * f(xk + h)
        k6 = h * f(xk + 0.5d0*h)

        y4 = yk + 0.115740740740741d0*k1 + 0.548927875243665d0*k3 &
                + 0.535722994391612d0*k4 - 0.2d0*k5

        y5 = yk + 0.118518518518519d0*k1 + 0.518986354775828d0*k3 &
                + 0.506131490342017d0*k4 - 0.18d0*k5 &
                + 0.036363636363636d0*k6

        err = abs(y5 - y4)

        if (err > 0.0d0) then
            s = (tol / (2.0d0*err))**0.2d0
        else
            s = 2.0d0
        end if

        hnew = max(0.1d0, min(5.0d0, s)) * h

        if (err <= tol) then
            ykp1 = y5   ! use the higher-order solution
            xkp1 = xk + h
        else
            ykp1 = yk
            xkp1 = xk
        end if

    end subroutine RKF45_step_1D

    subroutine RKF45_1D_with_context(f, y0, x0, xmax, h0, tol, sol, context)

        use KIM_kinds, only: dp
        implicit none

        ! Arguments
        real(dp), intent(in)  :: y0, x0, xmax, h0, tol
        real(dp), intent(out) :: sol
        class(*), intent(in) :: context

        ! Locals
        real(dp) :: yk, xk, hk
        real(dp) :: ytrial, xtrial, hnew
        integer  :: step_count

        interface
            function f(x, context)
                use KIM_kinds, only: dp
                implicit none
                real(dp), intent(in) :: x
                class(*), intent(in) :: context
                real(dp) :: f
            end function
        end interface

        ! Initialization
        yk = y0
        xk = x0
        hk = h0
        step_count = 0

        ! Integration loop
        do while (xk < xmax)

            if (step_count > MAX_INT_COUNT) then
                print *, "Integration error: too many steps (", MAX_INT_COUNT, ")"
                print *, "Current xk: ", xk, " Current yk: ", yk, " Current hk: ", hk
                print *, "xmax: ", xmax, " tol: ", tol, " h0: ", h0
                stop
            end if

            if (xk + hk > xmax) hk = xmax - xk

            call RKF45_step_1D_with_context(f, yk, xk, hk, tol, xtrial, ytrial, hnew, context)

            if (xtrial > xk) then
                xk = xtrial
                yk = ytrial
            end if

            hk = hnew
            step_count = step_count + 1

        end do

        sol = yk

    end subroutine RKF45_1D_with_context


    subroutine RKF45_step_1D_with_context(f, yk, xk, h, tol, xkp1, ykp1, hnew, context)
        ! Classical Runge-Kutta-Fehlberg 4(5) method (adaptive stepper)

        use KIM_kinds, only: dp
        implicit none

        real(dp), intent(in)  :: yk, xk, h, tol
        real(dp), intent(out) :: xkp1, ykp1, hnew
        class(*), intent(in) :: context

        real(dp) :: k1, k2, k3, k4, k5, k6
        real(dp) :: y4, y5, err, s

        interface
            function f(x, context)
                use KIM_kinds, only: dp
                implicit none
                real(dp), intent(in) :: x
                class(*), intent(in) :: context
                real(dp) :: f
            end function
        end interface

        k1 = h * f(xk, context)
        k2 = h * f(xk + 0.25d0*h, context)
        k3 = h * f(xk + 0.375d0*h, context)
        k4 = h * f(xk + 0.923076923076923d0*h, context)
        k5 = h * f(xk + h, context)
        k6 = h * f(xk + 0.5d0*h, context)

        y4 = yk + 0.115740740740741d0*k1 + 0.548927875243665d0*k3 &
                + 0.535331384015595d0*k4 - 0.2d0*k5

        y5 = yk + 0.118518518518519d0*k1 + 0.518986354775828d0*k3 &
                + 0.506131490342017d0*k4 - 0.18d0*k5 &
                + 0.036363636363636d0*k6

        err = abs(y5 - y4)

        if (err > 0.0d0) then
            s = (tol / (2.0d0*err))**0.2d0
        else
            s = 2.0d0
        end if

        hnew = max(0.1d0, min(5.0d0, s)) * h

        if (err <= tol) then
            ykp1 = y5   ! use the higher-order solution
            xkp1 = xk + h
        else
            ykp1 = yk
            xkp1 = xk
            !hnew = 0.5d0 * h
        end if

    end subroutine RKF45_step_1D_with_context

end module
