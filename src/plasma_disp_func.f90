double complex function plasma_Z(z) result (res)

    use use_libcerf
    
    implicit none
    double complex, intent(in) :: z

    res = exp(-z**2) * cerfc_F(z)

end function