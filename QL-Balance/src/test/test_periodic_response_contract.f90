program test_periodic_response_contract
    use QLBalance_kinds, only: dp
    use periodic_embedding_m, only: embed_complex_profile, embed_tensor_profile
    use periodic_embedding_m, only: resolved_embedding_width
    implicit none

    integer, parameter :: nl = 6, ng = 11
    real(dp), parameter :: local_r(nl) = [-1.0_dp, -0.5_dp, 0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp]
    real(dp), parameter :: global_r(ng) = [-2.0_dp, -1.0_dp, -0.5_dp, -0.25_dp, 0.0_dp, &
                                             0.5_dp, 1.0_dp, 1.25_dp, 1.5_dp, 2.0_dp, 2.5_dp]
    complex(dp) :: local_z(nl), global_z(ng)
    real(dp) :: weights(ng)
    real(dp) :: local_d(2, 2, nl), global_d(2, 2, ng)
    integer :: i, j, k

    local_z = cmplx(2.0_dp, -3.0_dp, dp)
    local_d = 0.0_dp
    do k = 1, nl
        local_d(:, :, k) = reshape([1.0_dp, 0.25_dp, 0.25_dp, 2.0_dp], [2, 2])
    end do

    call embed_complex_profile(local_r, local_z, global_r, 0.0_dp, 1.0_dp, 0.5_dp, &
                               global_z, weights)
    call embed_tensor_profile(local_r, local_d, global_r, 0.0_dp, 1.0_dp, 0.5_dp, &
                              global_d)

    do i = 1, ng
        if (global_r(i) < -0.5_dp .or. global_r(i) > 1.5_dp) then
            if (abs(global_z(i)) > 1.0e-14_dp) error stop 'field not zero outside compact support'
            if (maxval(abs(global_d(:, :, i))) > 1.0e-14_dp) &
                error stop 'tensor not zero outside compact support'
        else
            if (abs(global_z(i) - weights(i)*local_z(3)) > 1.0e-12_dp) &
                error stop 'field interpolation/weight mismatch'
            do j = 1, 2
                do k = 1, 2
                    if (abs(global_d(j, k, i) - weights(i)**2 * local_d(j, k, 3)) > 1.0e-12_dp) &
                        error stop 'tensor does not use squared field weight'
                end do
            end do
        end if
    end do

    if (abs(weights(1)) > 1.0e-14_dp .or. abs(weights(ng)) > 1.0e-14_dp) &
        error stop 'weight outside support'
    if (abs(weights(5) - 1.0_dp) > 1.0e-14_dp .or. abs(weights(7) - 1.0_dp) > 1.0e-14_dp) &
        error stop 'weight not one in trusted core'
    if (.not. (weights(4) > 0.0_dp .and. weights(4) < 1.0_dp)) &
        error stop 'left transition weight'
    if (.not. (weights(8) > 0.0_dp .and. weights(8) < 1.0_dp)) &
        error stop 'right transition weight'

    call test_solver_grid()
    print *, 'periodic response contract tests passed'

contains

    subroutine test_solver_grid()
        integer, parameter :: n = 64
        real(dp), parameter :: rm = 30.0_dp, core = 2.0_dp, requested = 1.0_dp
        real(dp), parameter :: length = 2.0_dp * (core + requested)
        real(dp) :: grid(n), points(7), effective, step, weights(7)
        complex(dp) :: field(n), embedded(7)
        real(dp) :: tensor(2, 2, n), embedded_tensor(2, 2, 7)
        integer :: i

        step = length / real(n, dp)
        do i = 1, n
            grid(i) = rm - length / 2.0_dp + real(i - 1, dp) * step
            field(i) = cmplx(grid(i), -2.0_dp * grid(i), dp)
            tensor(:, :, i) = grid(i)
        end do
        effective = resolved_embedding_width(grid, rm - core, rm + core, requested)
        points = [rm - core, rm, rm + core, rm + core + effective / 2.0_dp, &
            grid(n), grid(n) + step / 2.0_dp, rm + length / 2.0_dp]
        call embed_complex_profile(grid, field, points, rm - core, rm + core, &
            effective, embedded, weights)
        call embed_tensor_profile(grid, tensor, points, rm - core, rm + core, &
            effective, embedded_tensor)
        if (abs(effective - (requested - step)) > 1.0e-12_dp) &
            error stop 'transition must stay inside the sampled periodic domain'
        do i = 1, 3
            if (abs(embedded(i) - cmplx(points(i), -2.0_dp * points(i), dp)) > 1.0e-12_dp) &
                error stop 'periodic embedding changed a trusted-core field'
        end do
        if (abs(weights(4) - 0.5_dp) > 1.0e-12_dp) &
            error stop 'contracted transition is not smooth and centered'
        if (maxval(abs(embedded_tensor(:, :, 4) - 0.25_dp * points(4))) > 1.0e-12_dp) &
            error stop 'contracted tensor transition must use squared field weight'
        if (any(abs(embedded(5:7)) > 0.0_dp)) error stop 'field extrapolated past sampled support'
        if (any(abs(embedded_tensor(:, :, 5:7)) > 0.0_dp)) &
            error stop 'tensor extrapolated past sampled support'
    end subroutine test_solver_grid
end program test_periodic_response_contract
