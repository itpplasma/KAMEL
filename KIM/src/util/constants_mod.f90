module constants

    use KIM_kinds, only: dp
    
    implicit none

    real(dp) :: pi       = 3.14159265358979d0
    real(dp) :: sol      = 29979245800.0
    real(dp) :: e_mass   = 9.1094d-28
    real(dp) :: p_mass   = 1.6726d-24
    real(dp) :: e_charge = 4.8032d-10
    real(dp) :: ev       = 1.6022d-12
    real(dp) :: kB       = 1.380649e-16

    complex(dp) :: com_unit = (0.0, 1.0)

end module