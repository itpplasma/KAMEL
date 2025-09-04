module back_quants_m

    use KIM_kinds_m, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: A1e !thermodynamic force electrons
    real(dp), dimension(:), allocatable :: A2e !thermodynamic force electrons
    real(dp), dimension(:,:), allocatable :: A1i !thermodynamic force ions
    real(dp), dimension(:,:), allocatable :: A2i ! thermodynamic force ions

    real(dp), dimension(:), allocatable :: nue !collision frequency electrons
    real(dp), dimension(:,:), allocatable :: nui !collision frequency ions

    real(dp), dimension(:), allocatable :: vTe !thermal velocity electrons
    real(dp), dimension(:,:), allocatable :: vTi !thermal velocity ions

    real(dp), dimension(:), allocatable :: omce ! cyclotron frequency electrons
    real(dp), dimension(:,:), allocatable :: omci !cyclotron frequency ions

    real(dp), dimension(:), allocatable :: lambda_De ! Debye length electrons
    real(dp), dimension(:,:), allocatable :: lambda_Di ! Deby length ions

    real(dp), dimension(:), allocatable :: rho_Le ! Debye length electrons
    real(dp), dimension(:,:), allocatable :: rho_Li ! Deby length ions

    real(dp), dimension(:), allocatable :: ks ! "senkrecht" wavenumber
    real(dp), dimension(:), allocatable :: kp ! parallel wavenumber

    real(dp), dimension(:), allocatable :: om_E ! ExB rotation frequency

    complex(dp), dimension(:), allocatable :: z0e ! parameter z for m_\phi=0 for electrons
    complex(dp), dimension(:,:), allocatable :: z0i ! parameter z for m_\phi=0 for electrons


end module