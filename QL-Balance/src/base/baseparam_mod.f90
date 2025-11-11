module baseparam_mod
    use QLBalance_kinds, only: dp

    implicit none

    real(dp), parameter :: pi = 3.14159265358979_dp
    real(dp), parameter :: c = 29979245800.0_dp  ! cm/s
    real(dp), parameter :: e_charge = 4.8032e-10_dp  ! statC
    real(dp), parameter :: e_mass = 9.1094e-28_dp  ! g
    real(dp), parameter :: p_mass = 1.6726e-24_dp  ! g
    real(dp), parameter :: ev = 1.6022e-12_dp  ! eV in erg
    real(dp) :: btor, rtor, dperp, Z_i, am, pertamp, omega, rsepar

    real(dp) :: urelax  ! under relaxation factor
    real(dp) :: tol_max = 3e-2_dp  ! or 3e-2, 3e-3, 3e-4
    real(dp) :: factolmax = 3.0_dp  ! keep
    real(dp) :: factolred = 0.5_dp  ! keep

end module baseparam_mod
