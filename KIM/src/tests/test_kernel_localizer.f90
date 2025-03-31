program test_kernel_localizer

    use kernels, only: fill_rho_kernels, K_rho_phi_of_rg
    use KIM_kinds, only: dp
    use grid, only: rg_grid
    use plasma_parameter, only: r_prof, n_prof, ni_prof, Te_prof, Ti_prof, iprof_length, Er_prof

    implicit none

    integer :: i, ikr, ikrp
    real(dp) :: ne_core
    real(dp) :: delta_r
    integer :: r_ind = 50

    call kim_init_for_test
    ne_core = n_prof(1)
    delta_r = r_prof(iprof_length) - r_prof(1)
    do i=2, iprof_length
        n_prof(i) = ne_core * (1.0d0 - r_prof(i) / delta_r)
    end do

    call calculate_backs(.false.)
    call generate_grids
    call fill_rho_kernels

    ikr = 1
    ikrp = 50
    open(unit=10, file='data.txt', status='replace')
    do i = 1, rg_grid%npts_b
        write(10,*) rg_grid%xb(i), real(K_rho_phi_of_rg(ikr,ikrp,i)), dimag(K_rho_phi_of_rg(ikr,ikrp,i))
        print *, rg_grid%xb(i), real(K_rho_phi_of_rg(ikr,ikrp,i)), dimag(K_rho_phi_of_rg(ikr,ikrp,i))
    end do
    close(10)

    call system('gnuplot -e "p ''data.txt'' using 1:2 w l title ''real'', '''' using 1:3 w l title ''imag''; pause -1"')

end program