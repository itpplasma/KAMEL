program test_periodic_response_contract
    use QLBalance_kinds, only: dp
    use periodic_embedding_m, only: embed_complex_profile, embed_tensor_profile
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

    print *, 'periodic response contract tests passed'
end program test_periodic_response_contract
