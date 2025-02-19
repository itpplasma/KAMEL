subroutine get_D_one_over_nu


    use libneo_transport, only: calc_D_one_over_nu_e, init_gauss_laguerre_integration, calc_D_one_over_nu
    use libneo_species, only: species_t, init_deuterium_plasma
    use libneo_collisions, only: fill_species_arr_coulomb_log
    use math_constants, only: c, e, ev_to_cgs
    use libneo_kinds, only: real_kind

    implicit none

    integer, parameter :: num_species = 2
    integer :: i
    type(species_t) :: species_array(num_species)
    real(kind=real_kind) :: D11
    real(kind=real_kind) :: R0 = 550.0d0 ! values for W-7X
    real(kind=real_kind) :: B0 = 2.5d4 

    call init_deuterium_plasma(3.0d3, 1.5d3, 1.0d14, species_array)

    do i = 1, num_species
        species_array(i)%rho_L = species_array(i)%mass * c * &
            sqrt(species_array(i)%temp * ev_to_cgs / species_array(i)%mass)/ e / B0
    end do
        
    call fill_species_arr_coulomb_log(2, species_array)

    call init_gauss_laguerre_integration(5.0d0/2.0d0)
        
    call calc_D_one_over_nu(1, 2, species_array, R0, D11)
    print *, "D11 = ", D11


end subroutine
