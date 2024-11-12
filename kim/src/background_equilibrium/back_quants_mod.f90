module back_quants

    implicit none

    double precision, dimension(:), allocatable :: A1e !thermodynamic force electrons
    double precision, dimension(:), allocatable :: A2e !thermodynamic force electrons
    double precision, dimension(:,:), allocatable :: A1i !thermodynamic force ions
    double precision, dimension(:,:), allocatable :: A2i ! thermodynamic force ions

    double precision, dimension(:), allocatable :: nue !collision frequency electrons
    double precision, dimension(:,:), allocatable :: nui !collision frequency ions

    double precision, dimension(:), allocatable :: vTe !thermal velocity electrons
    double precision, dimension(:,:), allocatable :: vTi !thermal velocity ions

    double precision, dimension(:), allocatable :: omce ! cyclotron frequency electrons
    double precision, dimension(:,:), allocatable :: omci !cyclotron frequency ions

    double precision, dimension(:), allocatable :: lambda_De ! Debye length electrons
    double precision, dimension(:,:), allocatable :: lambda_Di ! Deby length ions

    double precision, dimension(:), allocatable :: ks ! "senkrecht" wavenumber
    double precision, dimension(:), allocatable :: kp ! parallel wavenumber

    double precision, dimension(:), allocatable :: om_E ! ExB rotation frequency

    double complex, dimension(:), allocatable :: z0e ! parameter z for m_\phi=0 for electrons
    double complex, dimension(:,:), allocatable :: z0i ! parameter z for m_\phi=0 for electrons
end module