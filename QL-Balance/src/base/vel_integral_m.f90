!> Velocity-space integral for the QL-Balance drift-kinetic response.
!>
!> Fortran port of the former vel_integral.cpp. Keeps the C symbol
!> calc_velocity_integral_ and drives the same fortnum semi-infinite adaptive
!> quadrature (fortnum_integrate_qagiu) with a bind(C) integrand, so the result
!> is numerically identical to the C++ version.
module vel_integral_m
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex, &
        c_ptr, c_funptr, c_loc, c_funloc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    private

    public :: calc_velocity_integral

    type, bind(C) :: quad_func_params
        integer(c_int) :: ind
        real(c_double) :: a
        real(c_double) :: omE
        real(c_double) :: nu
        complex(c_double_complex) :: c
        complex(c_double_complex) :: d
    end type quad_func_params

    interface
        function fortnum_integrate_qagiu(f, bound, inf, epsabs, epsrel, &
                                         value, abserr, ctx) result(code) &
            bind(C, name="fortnum_integrate_qagiu")
            import :: c_funptr, c_double, c_int, c_ptr
            type(c_funptr), value :: f
            real(c_double), value :: bound, epsabs, epsrel
            integer(c_int), value :: inf
            real(c_double) :: value, abserr
            type(c_ptr), value :: ctx
            integer(c_int) :: code
        end function fortnum_integrate_qagiu
    end interface

contains

    function vi_func(x, ctx) result(fx) bind(C)
        real(c_double), value :: x
        type(c_ptr), value :: ctx
        real(c_double) :: fx
        type(quad_func_params), pointer :: P
        complex(c_double_complex) :: field
        real(c_double) :: ind_fac, re_field

        call c_f_pointer(ctx, P)
        field = P%c + P%d*x
        re_field = real(field*conjg(field), c_double)

        select case (P%ind)
        case (1)
            ind_fac = 1.0_c_double
        case (2)
            ind_fac = 1.0_c_double + x*x
        case (3)
            ind_fac = 1.0_c_double + (1.0_c_double + x*x)*(1.0_c_double + x*x)
        case default
            ind_fac = 0.0_c_double
            write (error_unit, '(a)') "unknown index!"
        end select

        fx = exp(-x*x)/((P%a*x + P%omE)**2 + P%nu*P%nu)*re_field*ind_fac
    end function vi_func

    subroutine calc_velocity_integral(ind, vT, ks, kp, omE, nu, Es, Ep, res) &
        bind(C, name="calc_velocity_integral_")
        integer(c_int), value :: ind
        real(c_double), value :: vT, ks, kp, omE, nu
        real(c_double) :: Es(*), Ep(*), res(*)
        type(quad_func_params), target :: P
        real(c_double) :: err
        integer(c_int) :: code

        P%ind = ind
        P%a = sqrt(2.0_c_double)*vT*kp
        P%omE = omE
        P%nu = nu
        P%c = (omE/ks)*cmplx(Es(1), Es(2), c_double_complex)
        P%d = sqrt(2.0_c_double)*vT*cmplx(Ep(1), Ep(2), c_double_complex)

        code = fortnum_integrate_qagiu(c_funloc(vi_func), 0.0_c_double, 2_c_int, &
                                       1.0e-6_c_double, 1.0e-6_c_double, &
                                       res(1), err, c_loc(P))
    end subroutine calc_velocity_integral

end module vel_integral_m
