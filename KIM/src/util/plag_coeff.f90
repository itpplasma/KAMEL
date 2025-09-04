
subroutine binsrc(p,nmin,nmax,xi,i)
    !
    ! Finds the index  i  of the array of increasing numbers   p  with dimension  n
    ! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.
    !
    use KIM_kinds_m, only: dp
    implicit none

    integer                                :: n,nmin,nmax,i,imin,imax,k
    real(dp)                       :: xi
    real(dp), dimension(nmin:nmax) :: p

    imin=nmin
    imax=nmax
    n=nmax-nmin

    do k=1,n
      i=(imax-imin)/2+imin
      if(p(i).gt.xi) then
        imax=i
      else
        imin=i
      endif
      if(imax.eq.imin+1) exit
    enddo

    i=imax

    return
end subroutine binsrc


SUBROUTINE plag_coeff(npoi, nder, x, xp, coef)
    !
    ! npoi  - number of stencil points
    ! nder  - number of derivatives requested (0: value, 1: +1st deriv, 2: +2nd deriv)
    ! x     - evaluation point
    ! xp    - stencil abscissae (size npoi)
    ! coef  - weights: coef(0,:) for value, coef(1,:) for first deriv, coef(2,:) for second deriv (if nder>=2)
    !
    USE KIM_kinds_m, ONLY: dp
    IMPLICIT NONE

    INTEGER, INTENT(in)                                :: npoi, nder
    REAL(dp), INTENT(in)                                :: x
    REAL(dp), DIMENSION(npoi), INTENT(in)               :: xp
    REAL(dp), DIMENSION(0:nder, npoi), INTENT(out)      :: coef

    REAL(dp), DIMENSION(:), ALLOCATABLE                 :: dummy
    INTEGER                                             :: i, k, j
    REAL(dp)                                            :: fac
    REAL(dp)                                            :: s1, s2, denom

    ! -------- value weights (Lagrange basis values) --------
    DO i = 1, npoi
      coef(0, i) = 1.0_dp
      DO k = 1, npoi
          IF (k == i) CYCLE
          coef(0, i) = coef(0, i) * (x - xp(k)) / (xp(i) - xp(k))
      END DO
    END DO

    IF (nder == 0) RETURN

    ! -------- first-derivative weights (original logic, preserved) --------
    ALLOCATE(dummy(npoi))
    DO i = 1, npoi
      dummy = 1.0_dp
      dummy(i) = 0.0_dp
      DO k = 1, npoi
          IF (k == i) CYCLE
          fac = (x - xp(k)) / (xp(i) - xp(k))
          DO j = 1, npoi
            IF (j == k) THEN
                dummy(j) = dummy(j) / (xp(i) - xp(k))
            ELSE
                dummy(j) = dummy(j) * fac
            END IF
          END DO
      END DO
      coef(1, i) = SUM(dummy)
    END DO
    DEALLOCATE(dummy)

    ! -------- second-derivative weights (new; uses L_i, S1, S2 identities) --------
    IF (nder >= 2) THEN
      DO i = 1, npoi
          s1 = 0.0_dp
          s2 = 0.0_dp
          DO k = 1, npoi
            IF (k == i) CYCLE
            denom = x - xp(k)
            s1 = s1 + 1.0_dp / denom
            s2 = s2 + 1.0_dp / (denom * denom)
          END DO
          ! L_i''(x) = L_i(x) * (s1^2 - s2)
          coef(2, i) = coef(0, i) * (s1 * s1 - s2)
      END DO
    END IF

END SUBROUTINE plag_coeff
