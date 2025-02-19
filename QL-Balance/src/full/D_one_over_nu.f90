subroutine get_D_one_over_nu


    use libneo_transport, only: calc_D_one_over_nu_e, init_gauss_laguerre_integration,&
        calc_D_one_over_nu, w, x, gauss_laguerre_order
    use libneo_species, only: species_t, init_deuterium_plasma
    use libneo_collisions, only: fill_species_arr_coulomb_log
    use libneo_kinds, only: real_kind

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22, npoib
    
    use plasma_parameters, only: params
    use baseparam_mod, only: c, Btor, Rtor, ev, e_charge

    implicit none

    integer, parameter :: num_species = 2
    integer :: i, ipoi
    type(species_t) :: species_array(num_species)
    double precision, dimension(gauss_laguerre_order) :: w52, x52, w72, x72, w92, x92
    real(kind=real_kind) :: diff_coeff

    if (.not. allocated(Donue11)) then
        allocate(Donue11(npoib), Donue12(npoib), Donue21(npoib), Donue22(npoib))
        allocate(Donui11(npoib), Donui12(npoib), Donui21(npoib), Donui22(npoib))
    end if

    call init_gauss_laguerre_integration(5.0d0/2.0d0)
    w52 = w
    x52 = x
    call init_gauss_laguerre_integration(7.0d0/2.0d0)
    w72 = w
    x72 = x
    call init_gauss_laguerre_integration(9.0d0/2.0d0)
    w92 = w
    x92 = x

    do ipoi=1, npoib

        ! set up species, (Te, Ti, ne)
        call init_deuterium_plasma(params(3,ipoi), params(4,ipoi), params(1,ipoi), species_array)

        do i = 1, num_species
            species_array(i)%rho_L = species_array(i)%mass * c * &
                sqrt(species_array(i)%temp * ev / species_array(i)%mass)/ e_charge / Btor
        end do
        
        call fill_species_arr_coulomb_log(2, species_array)
        
        w = w52
        x = x52
        call calc_D_one_over_nu(1, 2, species_array, Rtor, diff_coeff)
        Donue11(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, diff_coeff)
        Donui11(ipoi) = diff_coeff

        w = w72
        x = x72
        call calc_D_one_over_nu(1, 2, species_array, Rtor, diff_coeff)
        Donue21(ipoi) = diff_coeff
        Donue21(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, diff_coeff)
        Donui21(ipoi) = diff_coeff
        Donui21(ipoi) = diff_coeff

        w = w92
        x = x92
        call calc_D_one_over_nu(1, 2, species_array, Rtor, diff_coeff)
        Donue22(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, diff_coeff)
        Donui22(ipoi) = diff_coeff

    end do

end subroutine
