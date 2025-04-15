
module gsl_mod

    interface
        function gsl_sf_bessel_In(n, x) bind(C, name="gsl_sf_bessel_In")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int), value :: n
            real(c_double), value :: x
            real(c_double) :: gsl_sf_bessel_In
        end function
    end interface

    interface
        function gsl_sf_erf(x) bind(C, name="gsl_sf_erf")
            use, intrinsic :: iso_c_binding
            implicit none
            real(c_double), value :: x
            real(c_double) :: gsl_sf_erf
        end function gsl_sf_erf
    end interface

end module 