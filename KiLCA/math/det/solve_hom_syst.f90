!------------------------------------------------------------------------------

subroutine solve_homegenious_linear_system (N, A, C)

implicit none;

integer, intent(in) :: N;
complex(8), dimension(N, N), intent(in) :: A;
complex(8), dimension(N), intent(out) :: C;

complex(8), dimension(N, N) :: Ac, VT;
real(8), dimension(N) :: S;

complex(8), allocatable, dimension(:) :: WORK;
complex(8), dimension(1) :: TEST;
real(8), dimension(5*N) :: RWORK;

real :: eps = 2.2d-16, err;

integer :: k, info, LWORK, rank;

Ac = A; !matrix copy

!query for optimal LWORK:
LWORK = -1;
call zgesvd ('N', 'A', N, N, Ac, N, S, Ac, 1, VT, N, TEST, LWORK, RWORK, info);

if (info /= 0) then
    print *, 'solve_hom_sys(): zgesvd: info=', info;
end if

LWORK = TEST(1);

allocate (WORK(LWORK));

call zgesvd ('N', 'A', N, N, Ac, N, S, Ac, 1, VT, N, WORK, LWORK, RWORK, info);

if (info /= 0) then
    print *, 'solve_hom_sys(): zgesvd: info=', info;
end if

deallocate (WORK);

rank = N;

!matrix rank:
do k = 1,N
    if (S(k) < N*S(1)*eps) then
        rank = k-1;
        exit;
    end if
end do

if (rank == N) then
    print *, 'warning: solve_homegenious_linear_system: determinant is not zero.'
    C = cmplx(0.0d0, 0.0d0, 8);
    return;
end if

!the right singular vector with the smallest singular value:
C = conjg(VT(N,:));

!check:
err = maxval(abs(matmul(A, C)));

if (err > N*S(1)*eps) then
    print *, 'warning: solve_homegenious_linear_system: suspicious error estimation =', err;
end if

end subroutine

!------------------------------------------------------------------------------
