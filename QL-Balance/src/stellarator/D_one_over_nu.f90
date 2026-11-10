subroutine get_D_one_over_nu

    use libneo_transport, only: init_gauss_laguerre_integration, calc_D_one_over_nu, gauss_laguerre_order
    use libneo_species, only: species_t, init_deuterium_plasma
    use libneo_collisions, only: fill_species_arr_coulomb_log
    use QLBalance_kinds, only: dp

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22, npoib, rb, ipbeg, ipend, reint_coef

    use plasma_parameters, only: params
    use baseparam_mod, only: c, Btor, Rtor, ev, e_charge

    implicit none

    integer, parameter :: num_species = 2
    real(dp) :: dens, temp_e, temp_i
    integer :: i, ipoi
    type(species_t) :: species_array(num_species)
    real(dp), dimension(gauss_laguerre_order) :: w52, x52, w72, x72, w92, x92
    real(dp) :: diff_coeff

    if (.not. allocated(Donue11)) then
        allocate(Donue11(npoib), Donue12(npoib), Donue21(npoib), Donue22(npoib))
        allocate(Donui11(npoib), Donui12(npoib), Donui21(npoib), Donui22(npoib))
    end if

    call init_gauss_laguerre_integration(5.0d0/2.0d0, w52, x52)
    call init_gauss_laguerre_integration(7.0d0/2.0d0, w72, x72)
    call init_gauss_laguerre_integration(9.0d0/2.0d0, w92, x92)

    ! The fluxes that use these coefficients live on the cell boundaries.
    !$omp parallel do private(ipoi, species_array, i, diff_coeff, dens, temp_e, temp_i)
    do ipoi=1, npoib

        dens = sum(params(1, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
        temp_e = sum(params(3, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
        temp_i = sum(params(4, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))

        ! set up species, (Te, Ti, ne)
        call init_deuterium_plasma(temp_e/ev, temp_i/ev, dens, species_array)
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

        if (.false.) print *, "ipoi = ", ipoi, " npoib = ", npoib, " e11 = ", Donue11(ipoi), " e12 = ", Donue12(ipoi), " e21 = ", &
            Donue21(ipoi), " e22 = ", Donue22(ipoi), " i11 = ", Donui11(ipoi), " i12 = ", Donui12(ipoi), &
            " i21 = ", Donui21(ipoi), " i22 = ", Donui22(ipoi)

    end do
    !$omp end parallel do

    if (.false.) then
        open(10, file='D_one_over_nu.dat')
        do ipoi=1, npoib
            write(10,*) rb(ipoi), Donue11(ipoi), Donue12(ipoi), Donue21(ipoi), Donue22(ipoi), &
                Donui11(ipoi), Donui12(ipoi), Donui21(ipoi), Donui22(ipoi)
        end do
        close(10)
    end if

end subroutine

subroutine rescale_D_one_over_nu_to_new_n_and_T

    use QLBalance_kinds, only: dp

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22, npoib, ipbeg, ipend, reint_coef, &
        init_Donue11, init_Donue12, init_Donue21, init_Donue22, &
        init_Donui11, init_Donui12, init_Donui21, init_Donui22

    use plasma_parameters, only: params, init_params

    implicit none

    integer :: ipoi
    real(dp) :: scale_e, scale_i

    do ipoi=1, npoib

        ! 1/nu scaling T/n at fixed remaining dependence, on the cell boundary
        scale_e = at_boundary(init_params, 1)/at_boundary(params, 1) &
            *at_boundary(params, 3)/at_boundary(init_params, 3)
        scale_i = at_boundary(init_params, 1)/at_boundary(params, 1) &
            *at_boundary(params, 4)/at_boundary(init_params, 4)

        Donue11(ipoi) = init_Donue11(ipoi)*scale_e
        Donue12(ipoi) = init_Donue12(ipoi)*scale_e
        Donue21(ipoi) = init_Donue21(ipoi)*scale_e
        Donue22(ipoi) = init_Donue22(ipoi)*scale_e

        Donui11(ipoi) = init_Donui11(ipoi)*scale_i
        Donui12(ipoi) = init_Donui12(ipoi)*scale_i
        Donui21(ipoi) = init_Donui21(ipoi)*scale_i
        Donui22(ipoi) = init_Donui22(ipoi)*scale_i

    end do

contains

    real(dp) function at_boundary(p, ieq)
        real(dp), intent(in) :: p(:, :)
        integer, intent(in) :: ieq
        at_boundary = sum(p(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
    end function

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
