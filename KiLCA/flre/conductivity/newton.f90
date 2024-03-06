!------------------------------------------------------------------------------
!Copyright (C) 2012 Ivan Ivanov (navi.adler@gmail.com)
!------------------------------------------------------------------------------

! The file containes the implementation of the Newton iterations for the solution of a nonlinear equation

!------------------------------------------------------------------------------

subroutine solve_system_newton(func, jacob, zmin, zmax, n_start, z_start, max_iter_num, &
                               abs_eps_z, rel_eps_z, abs_eps_f, rel_eps_f, z, f, status)

! The Newton iterations for the solution of a nonlinear equation

implicit none

external func, jacob  ! external subroutines that specify the nonlinear equation and its derivative

integer, parameter :: prcsn = 8

real(prcsn), intent(in)                      :: zmin, zmax           ! 1-dim box for the allowed solution
integer, intent(in)                          :: n_start              ! number of starting points
real(prcsn), dimension(n_start), intent(in)  :: z_start              ! starting 1-dim points
integer, intent(in)                          :: max_iter_num         ! maximum number of iterations
real(prcsn), intent(in)                      :: abs_eps_z, rel_eps_z ! tolerances for the solution
real(prcsn), intent(in)                      :: abs_eps_f, rel_eps_f ! tolerances for the function values
real(prcsn), intent(out)                     :: z                    ! the solution point
real(prcsn), intent(out)                     :: f                    ! the function value f(z)
integer, intent(out)                         :: status               ! status of iterations

integer     :: i, debug_level = 0, iter
real(prcsn) :: z_prev, f_prev, z_next, f_next
real(prcsn) :: j_prev
real(prcsn) :: mof

status = 1

do i = 1,n_start  ! loop over starting points

    ! z_prev, f_prev - current values of z and f(z)
    ! z_next, f_next - next values of z and f(z) after one Newton iteration

    z_prev = z_start(i)

    call func(z_prev, f_prev)

    mof = abs(f_prev)  ! function norm at the starting guess point

    if (abs(f_prev) < max(abs_eps_f, rel_eps_f * mof)) then ! z-prev already satisfies the function tolerances

        z = z_prev
        f = f_prev
        status = 0
        return

    end if

    iter = -1
    do while (iter < max_iter_num)  ! loop over Newtonian iterations

        iter = iter + 1

        call jacob(z_prev, j_prev)  ! evaluate the jacobian matrix at z_prev

        z_next = z_prev - f_prev / j_prev  ! next approximation to the solution

        ! check if the point z_next is still inside the 1-dim rectangle:
        if ( .not.(z_next >= zmin .and. (z_next <= zmax)) ) exit

        call func(z_next, f_next)  ! the function value at new point

        if (debug_level > 0) then
            call print_iteration_info(iter, z_prev, f_prev, z_next, f_next)
        end if

        ! check convergence criteria:
        if ( (abs(z_next - z_prev) < max(abs_eps_z, rel_eps_z * abs(z_next))) .and. &
                      (abs(f_next) < max(abs_eps_f, rel_eps_f * mof)) ) then

             z = z_next
             f = f_next
             status = 0
             return

        end if

        z_prev = z_next  ! prepare for the next iteration
        f_prev = f_next  ! prepare for the next iteration

    end do

end do

end subroutine

!------------------------------------------------------------------------------

subroutine print_iteration_info(iter, z_prev, f_prev, z_next, f_next)

implicit none

integer, parameter :: prcsn = 8

integer, intent(in)     :: iter
real(prcsn), intent(in) :: z_prev, f_prev, z_next, f_next

write(*,*) "iter=", iter
write(*,*) "z_prev=", z_prev
write(*,*) "z_next=", z_next
write(*,*) "f_prev=", f_prev
write(*,*) "f_next=", f_next

end subroutine

!------------------------------------------------------------------------------
