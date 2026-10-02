program test_manufactured_species
    use KIM_kinds_m, only: dp
    use constants_m, only: e_charge
    use fields_m, only: EBdat_t, postprocess_species_density, calculate_current_density
    use species_m, only: plasma, plasma_t
    use grid_m, only: xl_grid
    use kernel_m, only: kernel_spl_t
    use config_m, only: hdf5_output
    use IO_collection_m, only: h5id
    use KAMEL_hdf5_tools, only: h5_create, h5_close, h5_deinit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer, parameter :: n = 7
    integer, parameter :: charges(0:2) = [-1, 2, 6]
    real(dp), parameter :: abscissae(3) = [-sqrt(3.0_dp/5.0_dp), &
        0.0_dp, sqrt(3.0_dp/5.0_dp)]
    real(dp), parameter :: weights(3) = [5.0_dp/9.0_dp, 8.0_dp/9.0_dp, &
        5.0_dp/9.0_dp]
    complex(dp), parameter :: response_phi(0:2) = [ &
        cmplx(2.0e10_dp, -0.3e10_dp, dp), cmplx(-0.7e10_dp, 1.1e10_dp, dp), &
        cmplx(1.7e10_dp, 0.4e10_dp, dp)]
    complex(dp), parameter :: response_B(0:2) = [ &
        cmplx(1.4e10_dp, 0.2e10_dp, dp), cmplx(0.2e10_dp, -0.8e10_dp, dp), &
        cmplx(-1.3e10_dp, 0.5e10_dp, dp)]
    complex(dp), parameter :: current_phi = (1.2_dp, -0.8_dp)
    complex(dp), parameter :: current_B = (0.4_dp, 2.1_dp)
    real(dp) :: r(n), shape(2), h, t, qmass(n, n)
    complex(dp) :: exact(n, 0:2), charge_sum(n), current_exact(n)
    complex(dp) :: current_phi_kernel(n, n), current_B_kernel(n, n)
    complex(dp), allocatable :: current(:)
    type(EBdat_t) :: fields
    type(kernel_spl_t) :: kernel_phi, kernel_B
    integer :: sp, i, j, k, q

    r = [3.0_dp, 3.1_dp, 3.6_dp, 4.2_dp, 5.0_dp, 6.6_dp, 7.0_dp]
    plasma = plasma_t(n_species=3, grid_size=n)
    allocate(plasma%spec(0:2))
    plasma%spec(0)%name = 'e'
    plasma%spec(1)%name = 'He'
    plasma%spec(2)%name = 'C'
    do sp = 0, 2
        plasma%spec(sp)%Zspec = charges(sp)
    end do
    xl_grid%xb = r
    xl_grid%npts_b = n
    fields%r_grid = r
    fields%Phi = cmplx(0.02_dp*sin(r), -0.01_dp*cos(r), dp)
    fields%Br = cmplx(0.03_dp*cos(2.0_dp*r), 0.015_dp*sin(r), dp)

    ! Quadrature of the physical P1 basis, independent of projector assembly.
    qmass = 0.0_dp
    do i = 1, n-1
        h = r(i+1)-r(i)
        do q = 1, 3
            t = (abscissae(q)+1.0_dp)/2.0_dp
            shape(1) = 1.0_dp-t
            shape(2) = t
            do k = 1, 2
                do j = 1, 2
                    qmass(i+j-1, i+k-1) = qmass(i+j-1, i+k-1)+ &
                        shape(j)*shape(k)*h*weights(q)/2.0_dp
                end do
            end do
        end do
    end do
    call kernel_phi%init_kernel(n, n)
    call kernel_B%init_kernel(n, n)
    kernel_phi%Kllp_e = charges(0)*e_charge*response_phi(0)*qmass
    kernel_B%Kllp_e = charges(0)*e_charge*response_B(0)*qmass
    do sp = 1, 2
        kernel_phi%Kllp_i(:, :, sp) = charges(sp)*e_charge*response_phi(sp)*qmass
        kernel_B%Kllp_i(:, :, sp) = charges(sp)*e_charge*response_B(sp)*qmass
    end do
    kernel_phi%Kllp = kernel_phi%Kllp_e+sum(kernel_phi%Kllp_i, dim=3)
    kernel_B%Kllp = kernel_B%Kllp_e+sum(kernel_B%Kllp_i, dim=3)
    charge_sum = cmplx(0.0_dp, 0.0_dp, dp)
    do sp = 0, 2
        exact(:, sp) = response_phi(sp)*fields%Phi+response_B(sp)*fields%Br
        charge_sum = charge_sum+charges(sp)*e_charge*exact(:, sp)
    end do

    hdf5_output = .true.
    call h5_create('manufactured-species-density.h5', h5id)
    call postprocess_species_density(fields, kernel_phi, kernel_B)
    call h5_close(h5id)
    call h5_deinit()
    call assert_vector('signed electron Phi plus Br', fields%delta_n_e, exact(:, 0))
    call assert_vector('Z=2 ion Phi plus Br', fields%delta_n_i(:, 1), exact(:, 1))
    call assert_vector('Z=6 ion Phi plus Br', fields%delta_n_i(:, 2), exact(:, 2))
    call assert_vector('physical total charge', fields%rho, charge_sum)
    current_exact = current_phi*fields%Phi+current_B*fields%Br
    current_phi_kernel = current_phi*qmass
    current_B_kernel = current_B*qmass
    call calculate_current_density(current, fields%Phi, fields%Br, &
        current_phi_kernel, current_B_kernel)
    call assert_vector('raw physical parallel current', current, current_exact)
    print *, 'Manufactured species and current projection passed'

contains

    subroutine assert_vector(label, actual, expected)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: actual(:), expected(:)
        real(dp) :: relative_error
        if (.not. all(ieee_is_finite(real(actual, dp)))) &
            error stop 'Nonfinite real manufactured result'
        if (.not. all(ieee_is_finite(aimag(actual)))) &
            error stop 'Nonfinite imaginary manufactured result'
        if (.not. all(ieee_is_finite(real(expected, dp)))) &
            error stop 'Nonfinite real manufactured reference'
        if (.not. all(ieee_is_finite(aimag(expected)))) &
            error stop 'Nonfinite imaginary manufactured reference'
        relative_error = maxval(abs(actual-expected))/maxval(abs(expected))
        if (.not. ieee_is_finite(relative_error)) &
            error stop 'Nonfinite manufactured relative error'
        print *, trim(label), ': relative error = ', relative_error
        if (relative_error > 2.0e-12_dp) then
            print *, 'FAIL: ', trim(label)
            error stop 'Manufactured projection oracle failed'
        end if
    end subroutine assert_vector

end program test_manufactured_species
