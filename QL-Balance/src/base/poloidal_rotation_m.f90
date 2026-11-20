module poloidal_rotation_m
    !! Neoclassical poloidal rotation in the radial force balance.
    !!
    !! Markl et al., Nucl. Fusion 63 (2023) 126007, eq. (26):
    !!     E0r = 1/(e_i n_i) d(n_i T_i)/dr + sqrt(g) B0^theta/c (V^phi - q V^theta),
    !!     V^theta = c k B0_phi dT_i/dr / (e_i sqrt(g) <B0^2>),
    !! with the neoclassical coefficient k (NEO-2). QL-Balance carries the poloidal
    !! velocity as V_pol = r V^theta (input profile Vth, cm/s) on the boundary grid.
    !! The equilibrium is fixed, so k and the geometric factor are constant in time
    !! and V_pol is proportional to dT_i/dr:
    !!     V_pol(r, t) = C(r) dT_i/dr(r, t),   C = V_pol,input / dT_i/dr(initial state).
    !! C encodes k from the input Vth (cylinder: k = C e B0^2 / (c B0z), see
    !! rmp_integral_ion_tensor_kamel/scripts/neoclassical_k.py). It is set once, from
    !! the (unsmoothed) ion temperature gradient of the first state the code
    !! evaluates, so V_pol reproduces the input exactly at t = 0.
    !! Where dT_i/dr is exactly zero C is set to zero (V_pol = 0 there).
    !! Restart: C is formed from the restart state's stored Vth and T_i; since the
    !! stored Vth is the evolved C dT_i/dr, this recovers C up to the change of
    !! dT_i/dr within one step. A restart from output written with a fixed Vth
    !! uses that Vth as reference instead.
    use QLBalance_kinds, only: dp
    implicit none
    private
    public :: init_poloidal_rotation_coefficient, poloidal_velocity, &
              poloidal_rotation_ready, poloidal_rotation_coefficient, &
              reset_poloidal_rotation

    real(dp), allocatable, save :: coef(:)  !! C = V_pol / (dT_i/dr) [cm/s per (erg/cm)]

contains

    subroutine init_poloidal_rotation_coefficient(vpol_ref, dTi_dr)
        !! Set C from a reference V_pol profile and dT_i/dr [erg/cm] on the same grid.
        real(dp), intent(in) :: vpol_ref(:), dTi_dr(:)
        if (size(vpol_ref) /= size(dTi_dr)) &
            error stop 'init_poloidal_rotation_coefficient: size mismatch'
        if (allocated(coef)) deallocate (coef)
        allocate (coef(size(dTi_dr)))
        where (dTi_dr /= 0.0_dp)
            coef = vpol_ref / dTi_dr
        elsewhere
            coef = 0.0_dp
        end where
    end subroutine init_poloidal_rotation_coefficient

    logical function poloidal_rotation_ready()
        poloidal_rotation_ready = allocated(coef)
    end function poloidal_rotation_ready

    function poloidal_velocity(dTi_dr) result(vpol)
        !! V_pol = C dT_i/dr on the grid of C.
        real(dp), intent(in) :: dTi_dr(:)
        real(dp) :: vpol(size(dTi_dr))
        if (.not. allocated(coef)) error stop 'poloidal_velocity: coefficient not initialised'
        if (size(dTi_dr) /= size(coef)) error stop 'poloidal_velocity: size mismatch'
        vpol = coef * dTi_dr
    end function poloidal_velocity

    function poloidal_rotation_coefficient() result(c)
        real(dp), allocatable :: c(:)
        if (.not. allocated(coef)) error stop 'poloidal_rotation_coefficient: not initialised'
        c = coef
    end function poloidal_rotation_coefficient

    subroutine reset_poloidal_rotation()
        !! Forget C (tests; a new run initialises it again).
        if (allocated(coef)) deallocate (coef)
    end subroutine reset_poloidal_rotation

end module poloidal_rotation_m
