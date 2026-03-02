program test_ampere_matrices

    use KIM_kinds_m, only: dp
    use grid_m, only: xl_grid, calc_mass_matrix, M_mat
    use ampere_matrices_m, only: calc_potential_matrix, &
        calc_weighted_mass_matrix, calc_weighted_potential_matrix

    implicit none

    integer, parameter :: n = 10
    real(dp) :: r_min = 10.0d0   ! cm
    real(dp) :: r_max = 50.0d0   ! cm
    real(dp), allocatable :: Q_mat_loc(:,:), M_hth(:,:), Q_hz(:,:)
    real(dp), allocatable :: hth_const(:), hz_const(:)
    real(dp) :: h_const_val, tol
    integer :: i, j

    tol = 1.0d-12

    ! Set up an equidistant xl_grid
    xl_grid%npts_b = n
    allocate(xl_grid%xb(n))
    do i = 1, n
        xl_grid%xb(i) = r_min + (r_max - r_min) * dble(i-1) / dble(n-1)
    end do

    ! Compute mass matrix
    allocate(M_mat(n, n))
    call calc_mass_matrix(M_mat)

    ! Compute potential matrix
    allocate(Q_mat_loc(n, n))
    call calc_potential_matrix(Q_mat_loc)

    ! Test 1: Q_mat symmetry
    do i = 2, n-1
        do j = 2, n-1
            if (abs(Q_mat_loc(i,j) - Q_mat_loc(j,i)) > tol) then
                print *, 'FAIL: Q_mat not symmetric at (', i, ',', j, ')'
                print *, '  Q(i,j) =', Q_mat_loc(i,j), ' Q(j,i) =', Q_mat_loc(j,i)
                error stop
            end if
        end do
    end do
    print *, 'PASS: Q_mat is symmetric'

    ! Test 2: Q_mat entries are positive on diagonal (for r > 0)
    do i = 2, n-1
        if (Q_mat_loc(i,i) <= 0.0d0) then
            print *, 'FAIL: Q_mat diagonal entry non-positive at i =', i
            error stop
        end if
    end do
    print *, 'PASS: Q_mat diagonal entries are positive'

    ! Test 3: Weighted mass matrix with constant h_theta = h_const
    !         should give M_hth = h_const * M_mat
    h_const_val = 0.05d0
    allocate(hth_const(n), hz_const(n))
    hth_const = h_const_val
    hz_const  = h_const_val

    allocate(M_hth(n, n), Q_hz(n, n))
    call calc_weighted_mass_matrix(M_hth, M_mat, hth_const)

    do i = 1, n
        do j = 1, n
            if (abs(M_hth(i,j) - h_const_val * M_mat(i,j)) > tol) then
                print *, 'FAIL: M_hth /= h * M_mat at (', i, ',', j, ')'
                print *, '  M_hth =', M_hth(i,j), ' expected =', h_const_val * M_mat(i,j)
                error stop
            end if
        end do
    end do
    print *, 'PASS: M_hth = h_const * M_mat for constant h_theta'

    ! Test 4: Weighted potential matrix with constant h_z = h_const
    !         should give Q_hz = h_const * Q_mat
    call calc_weighted_potential_matrix(Q_hz, Q_mat_loc, hz_const)

    do i = 1, n
        do j = 1, n
            if (abs(Q_hz(i,j) - h_const_val * Q_mat_loc(i,j)) > tol) then
                print *, 'FAIL: Q_hz /= h * Q_mat at (', i, ',', j, ')'
                print *, '  Q_hz =', Q_hz(i,j), ' expected =', h_const_val * Q_mat_loc(i,j)
                error stop
            end if
        end do
    end do
    print *, 'PASS: Q_hz = h_const * Q_mat for constant h_z'

    print *, 'All Ampere matrix tests passed.'

end program
