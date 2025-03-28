
module gsl_mod

    interface
        function gsl_sf_bessel_In(n, x) bind(C, name="gsl_sf_bessel_In")
            use, intrinsic :: iso_c_binding
            integer(c_int), value :: n
            real(c_double), value :: x
            real(c_double) :: gsl_sf_bessel_In
        end function
    end interface

end module 