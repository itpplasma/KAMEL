program test_sparse_singular

    use sparse_mod

    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: nrow = 3, ncol = 3, nz = 3
    integer, dimension(nz) :: irow
    integer, dimension(ncol+1) :: pcol
    real(kind=dp), dimension(nz) :: val
    real(kind=dp), dimension(nrow) :: b

    ! 3x3 matrix with an all-zero row 2:
    ! [ 2 1 0 ]
    ! [ 0 0 0 ]
    ! [ 0 0 3 ]
    irow = [1, 1, 3]
    pcol = [1, 2, 3, 4]
    val = [2.0_dp, 1.0_dp, 3.0_dp]
    b = [1.0_dp, 1.0_dp, 1.0_dp]

    print *, 'solving singular 3x3 system with zero row 2'
    call sparse_solve(nrow, ncol, nz, irow, pcol, val, b)

    print *, 'NO_ERROR_STOP: singular system was accepted, b = ', b

end program test_sparse_singular
