!<The subroutines used to build and solve linear system of stitching equations (at zones boundaries).

!------------------------------------------------------------------------------

subroutine update_system_matrix_and_rhs_vector (Nc, A, B, neq, ieq, nvar, ivar, M, J)

integer, intent(in) :: Nc;                       !number of unknown coeficients
complex(8), dimension(Nc,Nc), intent(out) :: A;  !total superposition matrix
complex(8), dimension(Nc), intent(out)    :: B;  !total rhs vector

integer, intent(in) :: neq;                      !number of eqs for the boundary
integer, intent(in) :: ieq;                      !current index of equations
integer, intent(in) :: nvar;                     !number of unknowns for the boundary
integer, intent(in) :: ivar;                     !current index of unknowns

complex(8), dimension(neq,nvar), intent(in) :: M; !superposition matrix for the boundary
complex(8), dimension(neq), intent(in)      :: J; !rhs vector for the boundary

integer :: ie1, ie2, iv1, iv2;

!first call:
if (ieq == 0) then
    A = cmplx(0.0d0, 0.0d0, 8);
    B = cmplx(0.0d0, 0.0d0, 8);
    return;
end if

!further calls:
ie1 = ieq;
ie2 = ie1 + neq - 1;

iv1 = ivar;
iv2 = iv1 + nvar - 1;

A(ie1:ie2, iv1:iv2) = M;

B(ie1:ie2) = J;

end subroutine

!------------------------------------------------------------------------------

subroutine find_superposition_coeffs (NoC, A, B, S)

integer, intent(in) :: NoC;                       !Number of Coefficients
complex(8), dimension(NoC,NoC), intent(in) :: A;  !system matrix
complex(8), dimension(NoC), intent(in) :: B;      !system rhs vector
complex(8), dimension(NoC), intent(out) :: S;     !coefficients vector

complex(8), dimension(NoC,NoC) :: Ac; !a copy of superposition matrix

integer, dimension(NoC) :: ipiv;

integer :: info;

Ac = A;
S = B;

!solve general complex equations: Ac matrix is destroyed after the call!!!
call zgesv (NoC, 1, Ac, NoC, ipiv, S, NoC, info);

if (info /= 0) then
    print *, 'warning: find_superposition_coeffs: failed to solve the system: ierr=', info

    if (info > 0) then
        print *, 'warning: find_superposition_coeffs: homogenious system solver is called...'

        Ac = A;
        call solve_homegenious_linear_system (NoC, Ac, S);
    end if
end if

end subroutine

!------------------------------------------------------------------------------

subroutine calc_system_determinant (NoC, A, det)

implicit none;

integer, intent(in) :: NoC;
complex(8), dimension(NoC,NoC), intent(in) :: A;
complex(8), intent(out) :: det

integer, dimension(NoC) :: ipiv;

complex(8), dimension(NoC,NoC) :: Ac; !a copy of A

Ac = A;

!Ac matrix is destroyed after the call!!!
call calc_det_lu (NoC, Ac, det, ipiv);

end subroutine

!------------------------------------------------------------------------------

subroutine eval_superposition_of_basis_functions (D, Nw, dim, basis, S, EB)

implicit none;

integer, intent(in) :: D, Nw, dim;
complex(8), dimension(D,Nw,dim), intent(in) :: basis;
complex(8), dimension(Nw), intent(in) :: S;
complex(8), dimension(D,dim), intent(out) :: EB;

integer :: k, i;

do k = 1,dim

    EB(:,k) = cmplx(0.0d0, 0.0d0, 8);

    do i = 1,Nw
        EB(:,k) = EB(:,k) + S(i)*basis(:,i,k);
    end do

end do

end subroutine

!------------------------------------------------------------------------------
