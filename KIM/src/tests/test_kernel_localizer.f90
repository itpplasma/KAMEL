program test_kernel_localizer

    use kernels, only: fill_rho_kernels, K_rho_phi_of_rg, kernel_rho_phi_of_kr_krp_rg
    use KIM_kinds, only: dp
    use grid, only: rg_grid
    use plasma_parameter, only: r_prof, n_prof, ni_prof, Te_prof, Ti_prof, iprof_length, Er_prof
    use resonances_mod, only: r_res
    use regularization_funcs, only: theta_middle, theta_right, theta_left

    implicit none

    integer :: i, ikr, ikrp
    real(dp) :: ne_core
    real(dp) :: delta_r
    integer :: r_ind = 50
    complex(dp) :: kernel_value
    real(dp) :: kr, krp, rg
    real(dp), allocatable :: localizer_middle(:)
    real(dp), allocatable :: localizer_right(:)
    real(dp), allocatable :: localizer_left(:)

    call kim_init_for_test
    ne_core = n_prof(1)
    delta_r = r_prof(iprof_length) + 3.0d0
    
    do i=2, iprof_length
        n_prof(i) = ne_core * (1.0d0 - r_prof(i) / delta_r)
    end do

    call calculate_backs(.false.)
    call generate_grids

    allocate(localizer_right(rg_grid%npts_b), localizer_middle(rg_grid%npts_b), &
        localizer_left(rg_grid%npts_b))

    r_res = 45.0d0

    open(unit=10, file='localizer.txt', status='replace')
    do i = 1, rg_grid%npts_b
        localizer_right(i) = theta_right(rg_grid%xb(i))
        localizer_middle(i) = theta_middle(rg_grid%xb(i))
        localizer_left(i) = theta_left(rg_grid%xb(i))
        write(10,*) rg_grid%xb(i), localizer_right(i), localizer_middle(i), localizer_left(i)
    end do
    close(10)
    call system('gnuplot -p -e "&
                set xlabel ''r_g [cm]''; &
                set ylabel ''localizer''; &
                set grid;&
                p ''localizer.txt'' using 1:2 w linespoints title ''right'',&
                '''' using 1:3 w linespoints title ''middle'', &
                '''' using 1:4 w linespoints title ''left''; &
                "')

    kr = 1.0d0
    krp = 1.0d0

    call fill_rho_kernels

    ikr = 1
    ikrp = 100
    open(unit=10, file='data.txt', status='replace')
    do i = 1, rg_grid%npts_b
        write(10,*) rg_grid%xb(i), real(K_rho_phi_of_rg(ikr,ikrp,i)) * theta_middle(rg_grid%xb(i)),&
            dimag(K_rho_phi_of_rg(ikr,ikrp,i)) * theta_middle(rg_grid%xb(i))
        !print *, rg_grid%xb(i), real(K_rho_phi_of_rg(ikr,ikrp,i)), dimag(K_rho_phi_of_rg(ikr,ikrp,i))
    end do
    close(10)

    call system('gnuplot -p -e "&
                set xlabel ''r_g [cm]''; &
                set ylabel ''localizer''; &
                set grid;&
                p ''data.txt'' using 1:2 w linespoints title ''real'',&
                '''' using 1:3 w linespoints title ''imag'';"')

end program