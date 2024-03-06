!------------------------------------------------------------------------------

subroutine calc_det_lu (N, A, det, ipiv)

!The routine computes det(A) via LU decomposition
!Attention: the subroutine destroys input matrix!!!

implicit none;

integer, intent(in) :: N;
complex(8), dimension(N, N), intent(inout) :: A;
complex(8), intent(out) :: det;
integer, dimension(N), intent(out) :: ipiv;

integer :: k, info;

complex(8), dimension(N) :: diag;
integer, dimension(N) :: iperm;
integer :: ier, L, i1, i2;

call zgetrf (N, N, A, N, ipiv, info);
!call zgetf2 (N, N, tmp, N, ipiv, info);

if (info /= 0) then
    print *, 'warning: calc_det_lu(): zgetrf failed: info=', info;
    if (info > 0) then
        det = cmplx(0.0d0, 0.0d0, 8);
        return;
    end if
end if

do k=1,N
    diag(k) = abs(A(k,k));
end do

!sorting:
call dpsort (diag, N, iperm, 1, ier);
if (ier /= 0) then
    print *, 'determinant: error: sorting failed!'
end if

L = N / 2;

!determinant:
det = cmplx(1.0d0, 0.0d0, 8);

do k = 1,L

    i1 = iperm(k);
    i2 = iperm(N+1-k);

    det = det * (A(i1, i1) * A(i2, i2));

    if (ipiv(i1) /= i1) then
        det = -det;
    end if

    if (ipiv(i2) /= i2) then
        det = -det;
    end if

end do

if (2*L /= N) then
    i1 = iperm(L+1);

    det = det * A(i1, i1);

    if (ipiv(i1) /= i1) then
        det = -det;
    end if
end if

! simple code but overflows are possible:
! det = cmplx(1.0d0, 0.0d0, 8);
!
! do k=1,N
!     det = det*A(k,k);
!     if (ipiv(k) /= k) then
!         det = -det;
!     end if
! end do

!print *, 'A=', A

end subroutine

!------------------------------------------------------------------------------
