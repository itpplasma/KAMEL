module getIfunc_config_m
    ! Compatibility flag for legacy callers of getIfunc. New callers can use
    ! getIfunc_model to select the conservation model explicitly.
    implicit none
    logical :: boole_energy_conservation = .true.
end module getIfunc_config_m

subroutine getIfunc(x1, x2, symbI)
    use getIfunc_config_m, only: boole_energy_conservation
    implicit none
    integer, parameter :: mnmax = 3
    integer :: conservation_model
    double precision, intent(in) :: x1, x2
    double complex, dimension(0:mnmax, 0:mnmax), intent(out) :: symbI

    conservation_model = merge(1, 0, boole_energy_conservation)
    call getIfunc_model(x1, x2, conservation_model, symbI)

end subroutine getIfunc

subroutine getIfunc_model(x1, x2, conservation_model, symbI)
    ! KiLCA FPGEN conservation corrections, including the coupled
    ! energy/momentum form of Leitner et al., Eq. (24).
    implicit none
    integer, parameter :: mnmax = 3
    integer, intent(in) :: conservation_model
    integer :: m, n
    double precision, intent(in) :: x1, x2
    double complex :: energy_denom, momentum_denom, coupling, determinant
    double complex, dimension(0:mnmax, 0:mnmax), intent(out) :: symbI
    double complex, dimension(0:mnmax, 0:mnmax) :: Imn

    call W2_arr(x1, x2, Imn)

    select case (conservation_model)
    case (0)
        ! Particle-number conservation only.
        symbI = Imn
    case (1)
        ! Particle-number and energy conservation.
        energy_denom = (1.d0, 0.d0) - Imn(0, 0) &
            + (2.d0, 0.d0)*Imn(2, 0) - Imn(2, 2)
        do m = 0, mnmax
            do n = 0, mnmax
                symbI(m, n) = Imn(m, n) &
                    + (Imn(m, 0) - Imn(m, 2)) &
                    * (Imn(n, 0) - Imn(n, 2))/energy_denom
            end do
        end do
    case (2)
        ! Particle-number and parallel-momentum conservation.
        momentum_denom = (1.d0, 0.d0) - Imn(1, 1)
        do m = 0, mnmax
            do n = 0, mnmax
                symbI(m, n) = Imn(m, n) &
                    + Imn(m, 1)*Imn(n, 1)/momentum_denom
            end do
        end do
    case (3)
        ! Particle-number, energy, and parallel-momentum conservation.
        energy_denom = (1.d0, 0.d0) - Imn(0, 0) &
            + (2.d0, 0.d0)*Imn(2, 0) - Imn(2, 2)
        momentum_denom = (1.d0, 0.d0) - Imn(1, 1)
        coupling = Imn(1, 2) - Imn(1, 0)
        determinant = energy_denom*momentum_denom - coupling*coupling
        do m = 0, mnmax
            do n = 0, mnmax
                symbI(m, n) = Imn(m, n) + ( &
                    coupling*Imn(m, 1)*(Imn(n, 2) - Imn(n, 0)) &
                    + energy_denom*Imn(m, 1)*Imn(n, 1) &
                    + momentum_denom*(Imn(m, 2) - Imn(m, 0)) &
                    * (Imn(n, 2) - Imn(n, 0)) &
                    + coupling*(Imn(m, 2) - Imn(m, 0))*Imn(n, 1) &
                    )/determinant
            end do
        end do
    case default
        error stop 'getIfunc conservation model must be between 0 and 3'
    end select

end subroutine getIfunc_model
