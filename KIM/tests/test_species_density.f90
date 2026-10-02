program test_species_density
    use KIM_kinds_m, only: dp
    use kim_solver_m, only: kim_solver_t, kim_profiles_t, KIM_OK
    use fields_m, only: EBdat, postprocess_species_density
    use kernel_m, only: kernel_spl_t, FP_fill_kernels
    use grid_m, only: xl_grid
    use constants_m, only: e_charge, ev
    use species_m, only: plasma, plasma_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(kim_solver_t) :: solver
    type(kim_profiles_t) :: profiles
    type(kernel_spl_t) :: rho_phi, rho_B, j_phi, j_B
    complex(dp), allocatable :: electron_exact(:), ion_exact(:)
    real(dp) :: error, previous, x
    integer :: i, level, stat, n
    integer :: template_unit, config_unit, ios
    character(len=1024) :: line

    allocate(profiles%r(40), profiles%n(40), profiles%Te(40), &
        profiles%Ti(40), profiles%q(40), profiles%Er(40))
    do i = 1, 40
        x = real(i-1, dp)/39.0_dp
        profiles%r(i) = 3.0_dp+64.0_dp*x
        profiles%n(i) = 5.0e13_dp
        profiles%Te(i) = 1000.0_dp
        profiles%Ti(i) = 800.0_dp
        profiles%q(i) = 1.0_dp+3.0_dp*x
        profiles%Er(i) = 0.0_dp
    end do
    previous = huge(1.0_dp)
    do level = 1, 2
        ! Independent constructor fixtures use fresh global background storage.
        plasma = plasma_t(n_species=0, grid_size=0)
        open(newunit=template_unit, file='KIM_config_em_small.nml', status='old')
        open(newunit=config_unit, file='KIM_config_density.nml', status='replace')
        do
            read(template_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'l_space_dim =') > 0) line = '    l_space_dim = 16'
            if (index(line, 'rg_space_dim =') > 0) &
                write(line, '(A,I0)') '    rg_space_dim = ', 64*4**(level-1)
            if (index(line, 'gauss_int_nodes_Nx =') > 0) &
                write(line, '(A,I0)') '    gauss_int_nodes_Nx = ', 21*4**(level-1)
            write(config_unit, '(A)') trim(line)
        end do
        close(template_unit)
        close(config_unit)
        call solver%init('KIM_config_density.nml', profiles=profiles, stat=stat)
        if (stat /= KIM_OK) error stop 'Density fixture initialization failed'
        n = xl_grid%npts_b
        EBdat%r_grid = xl_grid%xb
        allocate(EBdat%Phi(n), EBdat%Br(n))
        rho_phi = kernel_spl_t(npts_l=0, npts_lp=0)
        rho_B = kernel_spl_t(npts_l=0, npts_lp=0)
        j_phi = kernel_spl_t(npts_l=0, npts_lp=0)
        j_B = kernel_spl_t(npts_l=0, npts_lp=0)
        call rho_phi%init_kernel(n, n)
        call rho_B%init_kernel(n, n)
        call j_phi%init_kernel(n, n)
        call j_B%init_kernel(n, n)
        call FP_fill_kernels(rho_phi, rho_B, j_phi, j_B)
        ! Exact static Boltzmann response to a prescribed electric potential;
        ! this reference uses physical n, signed charge, and temperature only.
        ! Compact nodal support buffers the finite guiding-centre domain;
        ! the exact reference concerns the interior Boltzmann limit.
        EBdat%Phi = cmplx(0.0_dp, 0.0_dp, dp)
        EBdat%Phi(3:n-2) = cmplx(0.02_dp, -0.01_dp, dp)* &
            sin(acos(-1.0_dp)*(EBdat%r_grid(3:n-2)-EBdat%r_grid(3))/ &
            (EBdat%r_grid(n-2)-EBdat%r_grid(3)))
        EBdat%Br = cmplx(0.0_dp, 0.0_dp, dp)
        call postprocess_species_density(EBdat, rho_phi, rho_B)
        electron_exact = e_charge*5.0e13_dp/(1000.0_dp*ev)*EBdat%Phi
        ion_exact = -e_charge*5.0e13_dp/(800.0_dp*ev)*EBdat%Phi
        call require_finite(EBdat%delta_n_e)
        call require_finite(EBdat%delta_n_i(:, 1))
        call require_finite(EBdat%rho)
        call require_finite(electron_exact)
        call require_finite(ion_exact)
        error = max(maxval(abs(EBdat%delta_n_e-electron_exact))/ &
            maxval(abs(electron_exact)), &
            maxval(abs(EBdat%delta_n_i(:, 1)-ion_exact))/ &
            maxval(abs(ion_exact)))
        print *, 'Actual kernel Boltzmann level, relative error: ', level, error
        print *, 'Electron/ion maximum error:', &
            maxval(abs(EBdat%delta_n_e-electron_exact))/maxval(abs(electron_exact)), &
            maxval(abs(EBdat%delta_n_i(:, 1)-ion_exact))/maxval(abs(ion_exact))
        print *, 'Interior electron/ion error:', &
            maxval(abs(EBdat%delta_n_e(3:n-2)-electron_exact(3:n-2)))/ &
            maxval(abs(electron_exact)), &
            maxval(abs(EBdat%delta_n_i(3:n-2, 1)-ion_exact(3:n-2)))/ &
            maxval(abs(ion_exact))
        if (level > 1) then
            if (error > 0.6_dp*previous) error stop 'Boltzmann density does not converge'
            if (error > 2.0e-3_dp) error stop 'Boltzmann density accuracy failed'
        end if
        previous = error
        if (maxval(abs(EBdat%rho-e_charge* &
            (EBdat%delta_n_i(:, 1)-EBdat%delta_n_e))) > &
            1.0e-13_dp*maxval(abs(EBdat%rho))) error stop 'Species charge sum failed'
        call solver%finalize()
    end do
    print *, 'Actual species kernels passed static Boltzmann density oracle'
contains
    subroutine require_finite(values)
        complex(dp), intent(in) :: values(:)
        if (.not. all(ieee_is_finite(real(values, dp)))) &
            error stop 'Nonfinite real Boltzmann quantity'
        if (.not. all(ieee_is_finite(aimag(values)))) &
            error stop 'Nonfinite imaginary Boltzmann quantity'
    end subroutine require_finite
end program test_species_density
