!------------------------------------------------------------------------------
!Copyright (C) 2012 Ivan Ivanov (navi.adler@gmail.com)
!------------------------------------------------------------------------------

! The file containes the implementation of the bisection iterations for the solution of a nonlinear equation

!------------------------------------------------------------------------------

subroutine solve_system_bisection(func, zmin, zmax, max_iter_num, abs_eps_z, rel_eps_z, abs_eps_f, rel_eps_f, z, f, status)

! The bisection iterations for the solution of a nonlinear equation

implicit none

external func ! external subroutines that specify the nonlinear equation

integer, parameter :: prcsn = 8

real(prcsn), intent(in)                      :: zmin, zmax           ! 1-dim box for the allowed solution
integer, intent(in)                          :: max_iter_num         ! maximum number of iterations
real(prcsn), intent(in)                      :: abs_eps_z, rel_eps_z ! tolerances for the solution
real(prcsn), intent(in)                      :: abs_eps_f, rel_eps_f ! tolerances for the function values
real(prcsn), intent(out)                     :: z                    ! the solution point
real(prcsn), intent(out)                     :: f                    ! the function value f(z)
integer, intent(out)                         :: status               ! status of iterations

integer     :: debug_level = 0, iter
real(prcsn) :: z1, f1, z2, f2, zc, fc
real(prcsn) :: mof

status = 1

z1 = zmin
z2 = zmax

call func(z1, f1)
call func(z2, f2)

if (f1 * f2 > 0.0d0) then

    print *, 'error: solve_system_bisection: the function has the same sign on both edges:', z1, f1, z2, f2
    return

end if

mof = (abs(f1) + abs(f2)) / 2.0d0  ! function norm at the beginning

iter = -1
do while (iter < max_iter_num)  ! loop over bisection iterations

    iter = iter + 1

    zc = (z1 + z2) / 2.0d0

    call func(zc, fc)

    if (debug_level > 0) then
        call print_bisection_info(iter, z1, f1, z2, f2, zc, fc)
    end if

    ! check convergence criteria:
    if ( (abs(z2 - z1) < max(abs_eps_z, rel_eps_z * abs(zc))) .and. &
         (abs(fc) < max(abs_eps_f, rel_eps_f * mof)) ) then

         z = zc
         f = fc
         status = 0
         return

    end if

    if (f1 * fc <= 0.0d0) then

        z2 = zc
        f2 = fc

    else if (f2 * fc <= 0.0d0) then

        z1 = zc
        f1 = fc

    else

        print *, 'error: solve_system_bisection: internal problem:', z1, f1, z2, f2, zc, fc

    end if

end do

end subroutine

!------------------------------------------------------------------------------

subroutine print_bisection_info(iter, z1, f1, z2, f2, zc, fc)

implicit none

integer, parameter :: prcsn = 8

integer, intent(in)     :: iter
real(prcsn), intent(in) :: z1, f1, z2, f2, zc, fc

write(*,*) "iter=", iter
write(*,*) "z1 =", z1
write(*,*) "f1 =", f1
write(*,*) "z2 =", z2
write(*,*) "f2 =", f2
write(*,*) "zc =", zc
write(*,*) "fc =", fc

end subroutine

!------------------------------------------------------------------------------
