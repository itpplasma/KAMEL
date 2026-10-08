program test_frozen_fourier_gyroaverage
    ! Frozen local Taylor coefficients, not a globally constant equilibrium.
    ! Only mphi=0, T'=Er=omega=0 and an independently specified n' are tested.
    use KIM_kinds_m, only: dp
    use constants_m, only: pi, sol, e_charge, e_mass, p_mass, ev, com_unit
    use species_m, only: plasma_t, nmmax, evaluate_susceptibility
    use config_m, only: artificial_debye_case
    use flr2_fourier_kernel_m, only: flr_arg_pair_sp, core_rho_B_sp, &
        core_j_phi_sp, kern_include_ks2, kern_zero_flr_electrons
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), parameter :: density = 8.0e13_dp, gradient = -1.0e12_dp
    real(dp), parameter :: temperature = 2000.0_dp*ev, magnetic = 53000.0_dp
    real(dp), parameter :: b_values(8) = &
        [0.0_dp, 0.001_dp, 0.01_dp, 0.1_dp, 0.5_dp, 1.0_dp, 3.0_dp, 10.0_dp]
    real(dp), parameter :: x_values(6) = &
        [-100.0_dp, -3.0_dp, -0.3_dp, 0.3_dp, 3.0_dp, 100.0_dp]
    complex(dp), parameter :: drive = (0.7_dp, -0.4_dp)
    type(plasma_t) :: local
    complex(dp) :: moments(0:nmmax, 0:nmmax)
    complex(dp) :: dn, current, dn_fluid, j_fluid, expected_dn, expected_j
    real(dp) :: gamma(8), quad_bound(8), coarse, fine, charge, mass, thermal
    real(dp) :: omega, rho, debye, kp, kr, ks, bplus, bcross, x1
    real(dp) :: ward_error, max_ward, error, max_error, current_scale
    integer :: sp, model, ix, ib, orientation, ks_sign, kr_sign, rows, unit
    character(len=1024) :: filename
    logical :: write_csv

    artificial_debye_case = 0
    kern_include_ks2 = .true.
    kern_zero_flr_electrons = .false.
    local%n_species = 2
    local%grid_size = 1
    allocate(local%spec(0:1), local%ks(1))
    do sp = 0, 1
        allocate(local%spec(sp)%A1(1), local%spec(sp)%A2(1), &
            local%spec(sp)%nu(1), local%spec(sp)%vT(1), &
            local%spec(sp)%omega_c(1), local%spec(sp)%lambda_D(1), &
            local%spec(sp)%rho_L(1), local%spec(sp)%I01(1, 0:0), &
            local%spec(sp)%I21(1, 0:0))
    end do
    do ib = 1, size(b_values)
        coarse = velocity_gyroaverage(b_values(ib), 16384)
        fine = velocity_gyroaverage(b_values(ib), 32768)
        if (.not. ieee_is_finite(coarse)) error stop 'Nonfinite coarse quadrature'
        if (.not. ieee_is_finite(fine)) error stop 'Nonfinite fine quadrature'
        ! Simpson Richardson correction; the independent quadrature uses J0,
        ! never the production modified Bessel routine or its implementation.
        gamma(ib) = fine+(fine-coarse)/15.0_dp
        quad_bound(ib) = 2.0_dp*abs(fine-coarse)/15.0_dp+1.0e-11_dp
        if (abs(fine-coarse) > 1.0e-7_dp) &
            error stop 'Velocity quadrature refinement failed'
    end do
    if (abs(gamma(1)-1.0_dp) > 1.0e-11_dp) &
        error stop 'Maxwellian perpendicular normalization failed'

    write_csv = command_argument_count() == 1
    if (command_argument_count() > 1) error stop 'Expected at most one CSV path'
    if (write_csv) then
        call get_command_argument(1, filename)
        open(newunit=unit, file=trim(filename), status='new', action='write')
        write(unit, '(A)') 'species,model,x1,b,orientation,ks_sign,kr_sign,'// &
            'rho_cm,kr_cm_inv,ks_cm_inv,gamma_reference,quadrature_bound,'// &
            'dn_real,dn_imag,dn_fluid_real,dn_fluid_imag,'// &
            'j_real,j_imag,j_fluid_real,j_fluid_imag'
    end if
    rows = 0
    max_ward = 0.0_dp
    max_error = 0.0_dp
    do sp = 0, 1
        charge = -e_charge
        mass = e_mass
        if (sp == 1) then
            charge = e_charge
            mass = 2.0_dp*p_mass
        end if
        thermal = sqrt(temperature/mass)
        omega = charge*magnetic/(mass*sol)
        rho = thermal/abs(omega)
        debye = sqrt(temperature/(4.0_dp*pi*density*charge**2))
        local%spec(sp)%A1 = gradient/density
        local%spec(sp)%A2 = 0.0_dp
        local%spec(sp)%vT = thermal
        local%spec(sp)%omega_c = omega
        local%spec(sp)%lambda_D = debye
        local%spec(sp)%rho_L = rho
        do model = 0, 1
            do ix = 1, size(x_values)
                x1 = x_values(ix)
                kp = sign(0.004_dp, x1)
                local%spec(sp)%nu = abs(kp)*thermal/abs(x1)
                call evaluate_susceptibility(x1, 0.0_dp, model, moments)
                if (.not. all(ieee_is_finite(real(moments, dp)))) &
                    error stop 'Nonfinite susceptibility real part'
                if (.not. all(ieee_is_finite(aimag(moments)))) &
                    error stop 'Nonfinite susceptibility imaginary part'
                ward_error = max(abs(com_unit*x1*moments(0, 1)-1.0_dp), &
                    abs(com_unit*x1*moments(1, 0)-1.0_dp))
                if (.not. ieee_is_finite(ward_error)) &
                    error stop 'Nonfinite particle conservation moment'
                max_ward = max(max_ward, ward_error)
                if (ward_error > 1.0e-10_dp) &
                    error stop 'Independent particle conservation identity failed'
                local%spec(sp)%I01(1, 0) = moments(0, 1)
                local%spec(sp)%I21(1, 0) = moments(2, 1)
                do ib = 1, size(b_values)
                    do orientation = 1, 3
                        do ks_sign = -1, 1, 2
                            do kr_sign = -1, 1, 2
                                kr = real(kr_sign, dp)*sqrt(b_values(ib))/rho
                                ks = 0.0_dp
                                if (orientation == 2) then
                                    kr = kr/sqrt(2.0_dp)
                                    ks = real(ks_sign, dp)*sqrt(b_values(ib)/2.0_dp)/rho
                                else if (orientation == 3) then
                                    kr = 0.0_dp
                                    ks = real(ks_sign, dp)*sqrt(b_values(ib))/rho
                                end if
                                local%ks = ks
                                call flr_arg_pair_sp(local, sp, kr, kr, 1, &
                                    bplus, bcross)
                                dn = core_rho_B_sp(local, sp, bplus, bcross, 1) &
                                    /(4.0_dp*pi*charge)*drive
                                current = core_j_phi_sp(local, sp, bplus, bcross, 1) &
                                    /(4.0_dp*pi)*drive
                                if (.not. ieee_is_finite(real(dn, dp))) &
                                    error stop 'Nonfinite density real part'
                                if (.not. ieee_is_finite(aimag(dn))) &
                                    error stop 'Nonfinite density imaginary part'
                                if (.not. ieee_is_finite(real(current, dp))) &
                                    error stop 'Nonfinite current real part'
                                if (.not. ieee_is_finite(aimag(current))) &
                                    error stop 'Nonfinite current imaginary part'
                                ! Physical signed density and current coefficients;
                                ! no fitted magnitude, phase, or collision scaling.
                                dn_fluid = -gradient/(com_unit*kp*magnetic)*drive
                                j_fluid = charge*sol*ks*gradient/(kp*magnetic)*drive
                                expected_dn = dn_fluid*gamma(ib)
                                expected_j = j_fluid*gamma(ib)
                                current_scale = abs(charge*sol*gradient &
                                    /(kp*magnetic*rho)*drive)
                                error = max(abs(dn-expected_dn)/abs(dn_fluid), &
                                    abs(current-expected_j)/current_scale)
                                if (.not. ieee_is_finite(error)) &
                                    error stop 'Nonfinite Fourier response'
                                max_error = max(max_error, error)
                                if (error > quad_bound(ib)+1.0e-10_dp) then
                                    print *, 'Failed sp,model,x1,b,orientation,error:', &
                                        sp, model, x1, b_values(ib), orientation, error
                                    error stop 'Independent gyroaverage response failed'
                                end if
                                if (write_csv) write(unit, &
                                    '(I0,",",I0,2(",",ES24.16E3),3(",",I0),'// &
                                    '13(",",ES24.16E3))') &
                                    sp, model, x1, b_values(ib), orientation, &
                                    ks_sign, kr_sign, rho, kr, ks, gamma(ib), &
                                    quad_bound(ib), real(dn, dp), aimag(dn), &
                                    real(dn_fluid, dp), aimag(dn_fluid), &
                                    real(current, dp), aimag(current), &
                                    real(j_fluid, dp), aimag(j_fluid)
                                rows = rows+1
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    if (write_csv) close(unit)
    print *, 'Frozen full Fourier response cases:', rows
    print *, 'Maximum Ward and independent gyroaverage errors:', max_ward, max_error
    print *, 'Does not validate the different global FEM kernel or an EM field solve.'
contains
    real(dp) function velocity_gyroaverage(b, panels) result(value)
        real(dp), intent(in) :: b
        integer, intent(in) :: panels
        real(dp) :: t, weight, step
        integer :: it

        ! Omitted positive tail <= exp(-48), since |J0(real)| <= 1.
        step = 48.0_dp/real(panels, dp)
        value = 0.0_dp
        do it = 0, panels
            t = step*real(it, dp)
            weight = 2.0_dp
            if (mod(it, 2) == 1) weight = 4.0_dp
            if (it == 0 .or. it == panels) weight = 1.0_dp
            value = value+weight*exp(-t)*bessel_j0(sqrt(2.0_dp*b*t))**2
        end do
        value = value*step/3.0_dp
    end function velocity_gyroaverage
end program test_frozen_fourier_gyroaverage
