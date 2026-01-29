module numerics_utils_m

    use KIM_kinds_m, only: dp
    use constants_m, only: pi
    use gsl_mod, only: erf => gsl_sf_erf

    implicit none

    contains

    function erf_diff(z1, z2) result(d)
        ! Numerically stable difference of error functions.
        ! Uses linearized form for small deltas to avoid cancellation.
        real(dp), intent(in) :: z1, z2
        real(dp) :: d
        real(dp) :: dz, zm

        dz = z1 - z2

        if (abs(dz) < 1.0e-6_dp) then
            zm = 0.5d0 * (z1 + z2)
            d = 2.0d0 / sqrt(pi) * exp(-zm*zm) * dz
        else
            d = erf(z1) - erf(z2)
        end if

    end function erf_diff

    ! Invert a square complex(dp) matrix using LAPACK (zgetrf + zgetri)
    subroutine invert_complex_matrix(A_in, A_inv)

        use KIM_kinds_m, only: dp

        implicit none

        complex(dp), intent(in)  :: A_in(:,:)
        complex(dp), allocatable, intent(out) :: A_inv(:,:)

        complex(dp), allocatable :: A(:,:), work(:)
        integer, allocatable     :: ipiv(:)
        integer :: n, lda, info, lwork
        complex(dp) :: work_query(1)

        interface
            subroutine zgetrf(m, n, a, lda, ipiv, info)
                use KIM_kinds_m, only: dp
                integer, intent(in) :: m, n, lda
                integer, intent(out) :: ipiv(*)
                integer, intent(out) :: info
                complex(dp) :: a(lda,*)
            end subroutine zgetrf
            subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
                use KIM_kinds_m, only: dp
                integer, intent(in) :: n, lda, lwork
                integer, intent(in) :: ipiv(*)
                integer, intent(out) :: info
                complex(dp) :: a(lda,*), work(*)
            end subroutine zgetri
        end interface

        if (size(A_in,1) /= size(A_in,2)) then
            stop "invert_complex_matrix: input must be square"
        end if

        n = size(A_in,1)
        lda = max(1, n)

        allocate(A(n,n))
        A = A_in

        allocate(ipiv(n))

        call zgetrf(n, n, A, lda, ipiv, info)
        if (info /= 0) then
            stop "invert_complex_matrix: LU factorization failed"
        end if

        ! Workspace query
        lwork = -1
        call zgetri(n, A, lda, ipiv, work_query, lwork, info)
        lwork = max(1, int(real(work_query(1))))
        allocate(work(lwork))

        call zgetri(n, A, lda, ipiv, work, lwork, info)
        if (info /= 0) then
            stop "invert_complex_matrix: inversion failed"
        end if

        allocate(A_inv(n,n))
        A_inv = A

        deallocate(A, ipiv, work)

    end subroutine invert_complex_matrix

    ! Invert a square real(dp) matrix using LAPACK (dgetrf + dgetri)
    subroutine invert_real_matrix(A_in, A_inv)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), intent(in)  :: A_in(:,:)
        real(dp), allocatable, intent(out) :: A_inv(:,:)

        real(dp), allocatable :: A(:,:), work(:)
        integer, allocatable  :: ipiv(:)
        integer :: n, lda, info, lwork
        real(dp) :: work_query(1)

        interface
            subroutine dgetrf(m, n, a, lda, ipiv, info)
                use KIM_kinds_m, only: dp
                integer, intent(in) :: m, n, lda
                integer, intent(out) :: ipiv(*)
                integer, intent(out) :: info
                real(dp) :: a(lda,*)
            end subroutine dgetrf
            subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
                use KIM_kinds_m, only: dp
                integer, intent(in) :: n, lda, lwork
                integer, intent(in) :: ipiv(*)
                integer, intent(out) :: info
                real(dp) :: a(lda,*), work(*)
            end subroutine dgetri
        end interface

        if (size(A_in,1) /= size(A_in,2)) then
            stop "invert_real_matrix: input must be square"
        end if

        n = size(A_in,1)
        lda = max(1, n)

        if (allocated(A)) deallocate(A)
        allocate(A(n,n))
        A = A_in

        allocate(ipiv(n))

        call dgetrf(n, n, A, lda, ipiv, info)
        if (info /= 0) then
            stop "invert_real_matrix: LU factorization failed"
        end if

        ! Workspace query
        lwork = -1
        call dgetri(n, A, lda, ipiv, work_query, lwork, info)
        lwork = max(1, int(work_query(1)))
        allocate(work(lwork))

        call dgetri(n, A, lda, ipiv, work, lwork, info)
        if (info /= 0) then
            stop "invert_real_matrix: inversion failed"
        end if

        allocate(A_inv(n,n))
        A_inv = A

        deallocate(A, ipiv, work)

    end subroutine invert_real_matrix

end module numerics_utils_m
