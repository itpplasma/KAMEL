subroutine get_D_one_over_nu

    use libneo_transport, only: init_gauss_laguerre_integration, calc_D_one_over_nu, gauss_laguerre_order
    use libneo_species, only: species_t, init_deuterium_plasma
    use libneo_collisions, only: fill_species_arr_coulomb_log
    use QLBalance_kinds, only: dp

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22, npoic, rb
    
    use plasma_parameters, only: params
    use baseparam_mod, only: c, Btor, Rtor, ev, e_charge

    implicit none

    integer, parameter :: num_species = 2
    integer :: i, ipoi
    type(species_t) :: species_array(num_species)
    real(dp), dimension(gauss_laguerre_order) :: w52, x52, w72, x72, w92, x92
    real(dp) :: diff_coeff

    if (.not. allocated(Donue11)) then
        allocate(Donue11(npoic), Donue12(npoic), Donue21(npoic), Donue22(npoic))
        allocate(Donui11(npoic), Donui12(npoic), Donui21(npoic), Donui22(npoic))
    end if

    call init_gauss_laguerre_integration(5.0d0/2.0d0, w52, x52)
    call init_gauss_laguerre_integration(7.0d0/2.0d0, w72, x72)
    call init_gauss_laguerre_integration(9.0d0/2.0d0, w92, x92)
    
    !$omp parallel do private(ipoi, species_array, i, diff_coeff)
    do ipoi=1, npoic

        ! set up species, (Te, Ti, ne)
        call init_deuterium_plasma(params(3,ipoi)/ev, params(4,ipoi)/ev, params(1,ipoi), species_array)
        do i = 1, num_species
            species_array(i)%rho_L = species_array(i)%mass * c * &
                sqrt(species_array(i)%temp * ev / species_array(i)%mass)/ e_charge / abs(Btor)
        end do

        call fill_species_arr_coulomb_log(2, species_array)

        call calc_D_one_over_nu(1, 2, species_array, Rtor, w52, x52, diff_coeff)
        Donue11(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, w52, x52, diff_coeff)
        Donui11(ipoi) = diff_coeff

        call calc_D_one_over_nu(1, 2, species_array, Rtor, w72, x72, diff_coeff)
        Donue12(ipoi) = diff_coeff
        Donue21(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, w72, x72, diff_coeff)
        Donui12(ipoi) = diff_coeff
        Donui21(ipoi) = diff_coeff

        call calc_D_one_over_nu(1, 2, species_array, Rtor, w92, x92, diff_coeff)
        Donue22(ipoi) = diff_coeff
        call calc_D_one_over_nu(2, 2, species_array, Rtor, w92, x92, diff_coeff)
        Donui22(ipoi) = diff_coeff

        if (.false.) print *, "ipoi = ", ipoi, " npoic = ", npoic, " e11 = ", Donue11(ipoi), " e12 = ", Donue12(ipoi), " e21 = ", &
            Donue21(ipoi), " e22 = ", Donue22(ipoi), " i11 = ", Donui11(ipoi), " i12 = ", Donui12(ipoi), &
            " i21 = ", Donui21(ipoi), " i22 = ", Donui22(ipoi)

    end do
    !$omp end parallel do

    if (.false.) then
        open(10, file='D_one_over_nu.dat')
        do ipoi=1, npoic
            write(10,*) rb(ipoi), Donue11(ipoi), Donue12(ipoi), Donue21(ipoi), Donue22(ipoi), &
                Donui11(ipoi), Donui12(ipoi), Donui21(ipoi), Donui22(ipoi)
        end do
        close(10)
    end if

end subroutine

subroutine rescale_D_one_over_nu_to_new_n_and_T

    use QLBalance_kinds, only: dp

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22, npoic, rb, &
        init_Donue11, init_Donue12, init_Donue21, init_Donue22, &
        init_Donui11, init_Donui12, init_Donui21, init_Donui22
    
    use plasma_parameters, only: params, init_params
    use baseparam_mod, only: c, Rtor, ev, e_charge

    implicit none

    integer :: ipoi

    do ipoi=1, npoic

        Donue11(ipoi) = init_Donue11(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(3, ipoi) / init_params(3, ipoi)
        Donui11(ipoi) = init_Donui11(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(4, ipoi) / init_params(4, ipoi)

        Donue12(ipoi) = init_Donue12(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(3, ipoi) / init_params(3, ipoi)
        Donue21(ipoi) = init_Donue21(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(3, ipoi) / init_params(3, ipoi)
        Donui12(ipoi) = init_Donui12(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(4, ipoi) / init_params(4, ipoi)
        Donui21(ipoi) = init_Donui21(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(4, ipoi) / init_params(4, ipoi)

        Donue22(ipoi) = init_Donue22(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(3, ipoi) / init_params(3, ipoi)
        Donui22(ipoi) = init_Donui22(ipoi) * init_params(1, ipoi) / params(1, ipoi) &
            * params(4, ipoi) / init_params(4, ipoi)

        if (.false.) print *, "ipoi = ", ipoi, " npoic = ", npoic, " e11 = ", Donue11(ipoi), " e12 = ", Donue12(ipoi), " e21 = ", &
            Donue21(ipoi), " e22 = ", Donue22(ipoi), " i11 = ", Donui11(ipoi), " i12 = ", Donui12(ipoi), &
            " i21 = ", Donui21(ipoi), " i22 = ", Donui22(ipoi)

    end do

end subroutine

subroutine initialize_D_one_over_nu

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22, &
        init_Donue11, init_Donue12, init_Donue21, init_Donue22, &
        init_Donui11, init_Donui12, init_Donui21, init_Donui22

    implicit none

    call get_D_one_over_nu

    init_Donue11 = Donue11
    init_Donue12 = Donue12
    init_Donue21 = Donue21
    init_Donue22 = Donue22

    init_Donui11 = Donui11
    init_Donui12 = Donui12
    init_Donui21 = Donui21
    init_Donui22 = Donui22

end subroutine
