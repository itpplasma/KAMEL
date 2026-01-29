module constants_m

    use KIM_kinds_m, only: dp

    implicit none

    real(dp), parameter :: pi       = 3.14159265358979d0
    real(dp), parameter :: sol      = 29979245800.0
    real(dp), parameter :: e_mass   = 9.1094d-28
    real(dp), parameter :: p_mass   = 1.6726d-24
    real(dp), parameter :: e_charge = 4.8032d-10
    real(dp), parameter :: ev       = 1.6022d-12
    real(dp), parameter :: kB       = 1.380649e-16

    complex(dp), parameter :: com_unit = (0.0d0, 1.0d0)

end module
