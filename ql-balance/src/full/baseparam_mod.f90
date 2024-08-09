
module baseparam_mod
    double precision, parameter :: pi=3.14159265358979d0
    double precision,parameter  :: c = 29979245800.0;
    double precision,parameter  :: e_charge=4.8032d-10
    double precision,parameter  :: e_mass=9.1094d-28
    double precision,parameter  :: p_mass=1.6726d-24
    double precision,parameter  :: ev=1.6022d-12
    double precision :: btor,rtor,dperp,Z_i,am,pertamp,omega,rsepar

    double precision :: urelax = 0.5e0 !0.5d0  !0.9d0
    double precision :: tol_max = 3.d-2 !3.d-4 !3.d-3 !3.d-2
    double precision :: factolmax = 3.d0 ! keep
    double precision :: factolred = 0.5d0 ! keep

end module baseparam_mod