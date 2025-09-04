double complex function plasma_Z(z) result (res)

    use use_libcerf_m
    use constants_m, only: pi
    
    implicit none
    double complex, intent(in) :: z

    res = cmplx(0.0d0, 1.0d0) * sqrt(pi) * w_of_z_F(z)!exp(-z**2) * cerfc_F(z)

end function