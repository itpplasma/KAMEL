
module baseparam_mod

    use QLBalance_kinds, only: dp
    
    implicit none

    real(dp), parameter :: pi=3.14159265358979d0
    real(dp), parameter  :: c = 29979245800.0;
    real(dp), parameter  :: e_charge=4.8032d-10
    real(dp), parameter  :: e_mass=9.1094d-28
    real(dp), parameter  :: p_mass=1.6726d-24
    real(dp), parameter  :: ev=1.6022d-12
    real(dp) :: btor,rtor,dperp,Z_i,am,pertamp,omega,rsepar

    real(dp) :: urelax ! under relaxation factor
    real(dp) :: tol_max = 3.d-2 !3.d-4 !3.d-3 !3.d-2
    real(dp) :: factolmax = 3.d0 ! keep
    real(dp) :: factolred = 0.5d0 ! keep

end module baseparam_mod