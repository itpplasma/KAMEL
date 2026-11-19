module local_dispersion_geometry_m

    use KIM_kinds_m, only: dp

    implicit none
    private

    public :: local_bessel_argument
    public :: local_laplacian_wavenumber_squared

contains

    pure complex(dp) function local_bessel_argument(kr_squared, ks, rho_t) result(b)
        complex(dp), intent(in) :: kr_squared
        real(dp), intent(in) :: ks, rho_t

        b = (kr_squared + ks**2) * rho_t**2
    end function local_bessel_argument

    pure complex(dp) function local_laplacian_wavenumber_squared( &
        kr_squared, ks, kparallel) result(k_squared)
        complex(dp), intent(in) :: kr_squared
        real(dp), intent(in) :: ks, kparallel

        k_squared = kr_squared + ks**2 + kparallel**2
    end function local_laplacian_wavenumber_squared

end module local_dispersion_geometry_m
