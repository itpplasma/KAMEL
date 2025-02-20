subroutine get_D_one_over_nu


    use libneo_transport, only: calc_D_one_over_nu_e, init_gauss_laguerre_integration,&
        calc_D_one_over_nu, w, x, gauss_laguerre_order
    use libneo_species, only: species_t, init_deuterium_plasma
    use libneo_collisions, only: fill_species_arr_coulomb_log
    use libneo_kinds, only: real_kind

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22, npoic, rb
    
    use plasma_parameters, only: params
    use baseparam_mod, only: c, Btor, Rtor, ev, e_charge

    implicit none

    integer, parameter :: num_species = 2
    integer :: i, ipoi
    type(species_t) :: species_array(num_species)
    double precision, dimension(gauss_laguerre_order) :: w52, x52, w72, x72, w92, x92
    real(kind=real_kind) :: diff_coeff

    if (.not. allocated(Donue11)) then
        allocate(Donue11(npoic), Donue12(npoic), Donue21(npoic), Donue22(npoic))
        allocate(Donui11(npoic), Donui12(npoic), Donui21(npoic), Donui22(npoic))
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

    do ipoi=1, npoic

        ! set up species, (Te, Ti, ne)
        call init_deuterium_plasma(params(3,ipoi)/ev, params(4,ipoi)/ev, params(1,ipoi), species_array)
        do i = 1, num_species
            species_array(i)%rho_L = species_array(i)%mass * c * &
                sqrt(species_array(i)%temp * ev / species_array(i)%mass)/ e_charge / abs(Btor)
        end do

        !print *, "ipoi = ", ipoi,  " Te = ", species_array(1)%temp, "Ti = ", species_array(2)%temp
        call fill_species_arr_coulomb_log(2, species_array)

        
        w = w52
        x = x52
        call calc_D_one_over_nu(1, 2, species_array, Rtor, diff_coeff)
        Donue11(ipoi) = diff_coeff
        print *, "ipoi = ", ipoi, " Donue11 = ", diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, diff_coeff)
        Donui11(ipoi) = diff_coeff

        w = w72
        x = x72
        call calc_D_one_over_nu(1, 2, species_array, Rtor, diff_coeff)
        Donue12(ipoi) = diff_coeff
        Donue21(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, diff_coeff)
        Donui12(ipoi) = diff_coeff
        Donui21(ipoi) = diff_coeff

        w = w92
        x = x92
        call calc_D_one_over_nu(1, 2, species_array, Rtor, diff_coeff)
        Donue22(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, diff_coeff)
        Donui22(ipoi) = diff_coeff

    end do

    print *, "D_one_over_nu calculated"

    open(10, file='D_one_over_nu.dat')
    do ipoi=1, npoic
        write(10,*) rb(ipoi), Donue11(ipoi), Donue12(ipoi), Donue21(ipoi), Donue22(ipoi), &
            Donui11(ipoi), Donui12(ipoi), Donui21(ipoi), Donui22(ipoi)
    end do
    close(10)

end subroutine
