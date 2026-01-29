!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE sample_matrix(iab,ierr)
!
  USE sample_matrix_mod
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: nder=0
  DOUBLE PRECISION, PARAMETER :: symm_break=0.01d0
  INTEGER :: i,j,iter,npoi_old,iold,inew,ibeg,iend,nshift,npoilag,n1,n2,ierr,iab
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: isplit
!
  DOUBLE PRECISION :: h,hh
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xold
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef,amat_maxmod
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: amat1,amat2
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: amat_old,wmat_old
!
  ierr=0
!
  npoilag=nlagr+1
  nshift=nlagr/2
  npoi=npoilag+1
!
  ALLOCATE(coef(0:nder,npoilag))
!
  h=(xend-xbeg)/npoilag/(1.d0+symm_break)
  hh=symm_break*h/npoilag
!
  x=xbeg
  CALL get_matrix(iab)
!
  n1=nsize
  n2=nsize
!
  ALLOCATE(xarr(npoi),amat_arr(n1,n2,npoi))
  IF(iab.EQ.3) ALLOCATE(wmat(nsize,nrhs,npoi))
!
  xarr(1)=x
  amat_arr(:,:,1)=amat
  IF(iab.EQ.3) wmat(:,:,1)=brhs
!
  x=xend
  CALL get_matrix(iab)
  xarr(npoi)=x
  amat_arr(:,:,npoi)=amat
  IF(iab.EQ.3) wmat(:,:,npoi)=brhs
!
  DO i=2,npoi-1
    x=xbeg+h*(i-1)+hh*(i-1)**2
    xarr(i)=x
  ENDDO
  if(iflag_ant.eq.1.and.x_ant.gt.xarr(1).and.x_ant.lt.xarr(npoi)) then
    DO i=2,npoi-1
      if(x_ant.lt.xarr(i)) then
        if(xarr(i)-x_ant.lt.(xarr(i)-xarr(i-1))*antshift) then
          iflag_ant=0
          xarr(i)=x_ant
        endif
      endif
    ENDDO
    DO i=npoi-1,2,-1
      if(x_ant.gt.xarr(i)) then
        if(x_ant-xarr(i).lt.(xarr(i+1)-xarr(i))*antshift) then
          iflag_ant=0
          xarr(i)=x_ant
        endif
      endif
    ENDDO
  endif
  DO i=2,npoi-1
    x=xarr(i)
    CALL get_matrix(iab)
    amat_arr(:,:,i)=amat
    IF(iab.EQ.3) wmat(:,:,i)=brhs
  ENDDO
!
  ALLOCATE(amat1(n1,n2),amat2(n1,n2),amat_maxmod(n1,n2))
!
! first check which intervals should be splitted
  ALLOCATE(isplit(npoi))
  isplit=0
  DO inew=1,npoi-1
    x=0.5d0*(xarr(inew)+xarr(inew+1))
    ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        amat_maxmod(i,j)=MAXVAL(ABS(amat_arr(i,j,ibeg:iend)))
      ENDDO
    ENDDO
    ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      amat2(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        amat_maxmod(i,j)=MAX(amat_maxmod(i,j),ABS(amat_arr(i,j,iend)))
      ENDDO
    ENDDO
    DO i=1,n1
      DO j=1,n2
        IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*amat_maxmod(i,j)) isplit(inew)=1
      ENDDO
    ENDDO
  ENDDO
  IF(MAXVAL(isplit).GT.0) THEN
    npoi_old=npoi
    ALLOCATE(xold(npoi),amat_old(n1,n2,npoi))
    xold=xarr
    amat_old=amat_arr
    IF(iab.EQ.3) THEN
      ALLOCATE(wmat_old(nsize,nrhs,npoi))
      wmat_old=wmat
    ENDIF
  ELSE
    IF(iab.EQ.3) THEN
      npoi_rhs=npoi
      ALLOCATE(xarr_rhs(npoi))
      xarr_rhs=xarr
    ENDIF
    RETURN
  ENDIF
!
  iter=0
  DO
    iter=iter+1
    IF(iter.GT.itermax) THEN
      ierr=2
      PRINT *,'sample_matrix : maximum number of iterations exceeded'
      RETURN
    ENDIF
!
! determine the dimension of new arrays and re-allocate them:
    DO iold=1,npoi_old-1
      IF(isplit(iold).EQ.1) npoi=npoi+1
    ENDDO
    IF(ALLOCATED(xarr)) THEN
      DEALLOCATE(xarr,amat_arr)
      IF(iab.EQ.3) DEALLOCATE(wmat)
    ENDIF
    ALLOCATE(xarr(npoi),amat_arr(n1,n2,npoi))
    IF(iab.EQ.3) ALLOCATE(wmat(nsize,nrhs,npoi))
!
! fill new arrays:
    inew=0
    DO iold=1,npoi_old-1
      inew=inew+1
      xarr(inew)=xold(iold)
      amat_arr(:,:,inew)=amat_old(:,:,iold)
      IF(iab.EQ.3) wmat(:,:,inew)=wmat_old(:,:,iold)
      IF(isplit(iold).EQ.1) THEN
        inew=inew+1
        x=0.5d0*(xold(iold)+xold(iold+1))
        if(iflag_ant.eq.1) then
          if(x.gt.xold(iold).and.x.lt.xold(iold+1)) then
            if(abs(x-x_ant).lt.abs(x-xold(iold))*antshift) then
              iflag_ant=0
              x=x_ant
            endif
          endif
        endif
        CALL get_matrix(iab)
        xarr(inew)=x
        amat_arr(:,:,inew)=amat
        IF(iab.EQ.3) wmat(:,:,inew)=brhs
      ENDIF
    ENDDO
    inew=inew+1
    xarr(inew)=xold(npoi_old)
    amat_arr(:,:,inew)=amat_old(:,:,npoi_old)
    IF(iab.EQ.3) wmat(:,:,inew)=wmat_old(:,:,npoi_old)
    DEALLOCATE(isplit)
!
! check which intervals should be splitted
    ALLOCATE(isplit(npoi))
    isplit=0
    DO inew=1,npoi-1
      x=0.5d0*(xarr(inew)+xarr(inew+1))
      ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          amat_maxmod(i,j)=MAXVAL(ABS(amat_arr(i,j,ibeg:iend)))
        ENDDO
      ENDDO
      ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        amat2(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          amat_maxmod(i,j)=MAX(amat_maxmod(i,j),ABS(amat_arr(i,j,iend)))
        ENDDO
      ENDDO
      DO i=1,n1
        DO j=1,n2
          IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*amat_maxmod(i,j)) isplit(inew)=1
        ENDDO
      ENDDO
    ENDDO
    IF(MAXVAL(isplit).GT.0) THEN
      npoi_old=npoi
      DEALLOCATE(xold,amat_old)
      ALLOCATE(xold(npoi),amat_old(n1,n2,npoi))
      xold=xarr
      amat_old=amat_arr
      IF(iab.EQ.3) THEN
        DEALLOCATE(wmat_old)
        ALLOCATE(wmat_old(nsize,nrhs,npoi))
        wmat_old=wmat
      ENDIF
    ELSE
      IF(iab.EQ.3) THEN
        npoi_rhs=npoi
        ALLOCATE(xarr_rhs(npoi))
        xarr_rhs=xarr
      ENDIF
      EXIT
    ENDIF
  ENDDO
!
  END SUBROUTINE sample_matrix
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE sample_eigvecs(iregime,iab,ierr)
!
! Computes egenvectors and eigenvalues for given matrix array sapled on
! the grid. Continues refinement of the grid using eigenvectors and eigenvalues
! as a measure of sampling quality.
!
! Input parameters:
!           Formal: iregime - regime switch, iregime=0 - allocate and and
!                             compute eigenvalues and eigenvector arrays,
!                             then proceed with the refinement;
!                             iregime=1 - do the refiement only.
! Output parameters:
!           Formal: ierr    - error message, ierr=0 - normal work
!
  USE sample_matrix_mod
!
  IMPLICIT NONE
!
  INTEGER,          PARAMETER :: nder=0
!
  INTEGER :: i,j,iter,npoi_old,iold,inew,ibeg,iend,nshift,npoilag,n1,n2
  INTEGER :: info,ierr,nsmall,keep,iregime,iab
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: isplit,ipoi
!
  DOUBLE PRECISION :: h,Deltax,philowlim
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xold
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef,phi_maxmod
  DOUBLE COMPLEX,   DIMENSION(:),     ALLOCATABLE :: alam_tmp
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: alam_old,phi_tmp
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: phi1,phi2
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: amat_old,phi_old,wmat_old
!
  ierr=0
  philowlim=1.d2*epsilon(1.d0)/eps
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  n1=nsize
  n2=nsize
!
  ALLOCATE(alam_tmp(nsize),phi_tmp(nsize,nsize),ipoi(nsize))
!
  IF(iregime.EQ.0) THEN
!
! Compute and align eigenvalues and eigenvectors:
!
    ALLOCATE(alam(nsize,npoi),phi(n1,n2,npoi))
!
    DO i=1,npoi
      CALL eigen_cmat(nsize,amat_arr(:,:,i),alam(:,i),phi(:,:,i),info)
      IF(info.NE.0) THEN
        ierr=1
        PRINT *,'sample_eigvecs error: eigenvectors failed'
        RETURN
      ENDIF
    ENDDO
!
    CALL order_eigvecs(nsize,npoi,xarr,alam,phi,ierr)
!
  ENDIF
!
! Mesh refinement:
!
  ALLOCATE(coef(0:nder,npoilag))
  ALLOCATE(phi1(n1,n2),phi2(n1,n2),phi_maxmod(n1,n2))
!
! first check which intervals should be splitted
  ALLOCATE(isplit(npoi))
  isplit=0
  DO inew=1,npoi-1
    x=0.5d0*(xarr(inew)+xarr(inew+1))
    ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      phi1(i,:)=MATMUL(phi(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        phi_maxmod(i,j)=MAXVAL(ABS(phi(i,j,ibeg:iend)))
      ENDDO
    ENDDO
    ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      phi2(i,:)=MATMUL(phi(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        phi_maxmod(i,j)=MAX(phi_maxmod(i,j),ABS(phi(i,j,iend)),philowlim)
      ENDDO
    ENDDO
    DO i=1,n1
      DO j=1,n2
        IF(ABS(phi1(i,j)-phi2(i,j)).GT.eps*phi_maxmod(i,j)) isplit(inew)=1
      ENDDO
    ENDDO
  ENDDO
  IF(MAXVAL(isplit).GT.0) THEN
    npoi_old=npoi
    ALLOCATE(xold(npoi),amat_old(n1,n2,npoi),      &
             phi_old(n1,n2,npoi),alam_old(nsize,npoi))
    xold=xarr
    amat_old=amat_arr
    phi_old=phi
    alam_old=alam
    IF(iab.EQ.3) THEN
      ALLOCATE(wmat_old(nsize,nrhs,npoi))
      wmat_old=wmat
    ENDIF
  ELSE
    IF(iab.EQ.3) THEN
      DEALLOCATE(xarr_rhs)
      npoi_rhs=npoi
      ALLOCATE(xarr_rhs(npoi))
      xarr_rhs=xarr
    ENDIF
    RETURN
  ENDIF
!
  iter=0
  DO
    iter=iter+1
    IF(iter.GT.itermax) THEN
      ierr=2
      PRINT *,'sample_eigvecs : maximum number of iterations exceeded'
      RETURN
    ENDIF
!
! determine the dimension of new arrays and re-allocate them:
    DO iold=1,npoi_old-1
      IF(isplit(iold).EQ.1) npoi=npoi+1
    ENDDO
    IF(ALLOCATED(xarr)) THEN
      DEALLOCATE(xarr,amat_arr,phi,alam)
      IF(iab.EQ.3) DEALLOCATE(wmat)
    ENDIF
    ALLOCATE(xarr(npoi),amat_arr(n1,n2,npoi))
    ALLOCATE(alam(nsize,npoi),phi(n1,n2,npoi))
    IF(iab.EQ.3) ALLOCATE(wmat(nsize,nrhs,npoi))
!
! fill new arrays:
    inew=0
    DO iold=1,npoi_old-1
      inew=inew+1
      xarr(inew)=xold(iold)
      alam(:,inew)=alam_old(:,iold)
      amat_arr(:,:,inew)=amat_old(:,:,iold)
      phi(:,:,inew)=phi_old(:,:,iold)
      IF(iab.EQ.3) wmat(:,:,inew)=wmat_old(:,:,iold)
      IF(isplit(iold).EQ.1) THEN
        inew=inew+1
        x=0.5d0*(xold(iold)+xold(iold+1))
        if(iflag_ant.eq.1) then
          if(x.gt.xold(iold).and.x.lt.xold(iold+1)) then
            if(abs(x-x_ant).lt.abs(x-xold(iold))*antshift) then
              iflag_ant=0
              x=x_ant
            endif
          endif
        endif
        CALL get_matrix(iab)
        xarr(inew)=x
        amat_arr(:,:,inew)=amat
        IF(iab.EQ.3) wmat(:,:,inew)=brhs
        CALL eigen_cmat(nsize,amat,alam(:,inew),phi(:,:,inew),info)
        IF(info.NE.0) THEN
          ierr=1
          PRINT *,'sample_eigvecs error: eigenvectors failed'
          RETURN
        ENDIF
        alam_tmp=0.5d0*(alam_old(:,iold)+alam_old(:,iold+1))
        CALL align_lambdas(nsize,nsize,alam_tmp,alam(:,inew),ipoi,info)
        IF(info.GT.0) THEN
          ierr=2
          PRINT *,'sample_eigvecs error: cannot align eigenvalues'
        ELSEIF(info.LT.0) THEN
          alam_tmp=alam(ipoi,inew)
          alam(:,inew)=alam_tmp
          phi_tmp=phi(:,ipoi,inew)
          phi(:,:,inew)=phi_tmp
        ENDIF
        phi_tmp=0.5d0*(phi_old(:,:,iold)+phi_old(:,:,iold+1))
        DO j=1,nsize
          IF(REAL(SUM(CONJG(phi_tmp(:,j))*phi(:,j,inew))).LT.0.d0) THEN
            phi(:,j,inew)=-phi(:,j,inew)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    inew=inew+1
    xarr(inew)=xold(npoi_old)
    alam(:,inew)=alam_old(:,npoi_old)
    amat_arr(:,:,inew)=amat_old(:,:,npoi_old)
    phi(:,:,inew)=phi_old(:,:,npoi_old)
    IF(iab.EQ.3) wmat(:,:,inew)=wmat_old(:,:,npoi_old)
    DEALLOCATE(isplit)
!
    CALL order_eigvecs(nsize,npoi,xarr,alam,phi,ierr)
!
! check which intervals should be splitted
    ALLOCATE(isplit(npoi))
    isplit=0
    DO inew=1,npoi-1
      x=0.5d0*(xarr(inew)+xarr(inew+1))
      ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        phi1(i,:)=MATMUL(phi(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          phi_maxmod(i,j)=MAXVAL(ABS(phi(i,j,ibeg:iend)))
        ENDDO
      ENDDO
      ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        phi2(i,:)=MATMUL(phi(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          phi_maxmod(i,j)=MAX(phi_maxmod(i,j),ABS(phi(i,j,iend)),philowlim)
        ENDDO
      ENDDO
      DO i=1,n1
        DO j=1,n2
          IF(ABS(phi1(i,j)-phi2(i,j)).GT.eps*phi_maxmod(i,j)) isplit(inew)=1
        ENDDO
      ENDDO
    ENDDO
    IF(MAXVAL(isplit).GT.0) THEN
      npoi_old=npoi
      DEALLOCATE(xold,amat_old,phi_old,alam_old)
      ALLOCATE(xold(npoi),amat_old(n1,n2,npoi))
      ALLOCATE(alam_old(nsize,npoi),phi_old(n1,n2,npoi))
      xold=xarr
      alam_old=alam
      amat_old=amat_arr
      phi_old=phi
      IF(iab.EQ.3) THEN
        DEALLOCATE(wmat_old)
        ALLOCATE(wmat_old(nsize,nrhs,npoi))
        wmat_old=wmat
      ENDIF
    ELSE
      IF(iab.EQ.3) THEN
        DEALLOCATE(xarr_rhs)
        npoi_rhs=npoi
        ALLOCATE(xarr_rhs(npoi))
        xarr_rhs=xarr
      ENDIF
      EXIT
    ENDIF
  ENDDO
!
  END SUBROUTINE sample_eigvecs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE align_lambdas(nsize,nstiff,alam_old,alam_new,ipoi,ierr)
!
! Finds indices j=ipoi(i) of the array alam_new so that their values are
! closest by module to elements i of array alam_old. i.e abs(alam_new(j))
! is closest to abs(alam(i)). If more than one element of alam_new are
! closer to some element of alam_old then to the rest elements, ierr=1.
! For normal work ierr < 1. If ierr=0, ipoi(i)=i, othervice ierr=-1.
!
  IMPLICIT NONE
!
  INTEGER :: nsize,nstiff,i,j,ierr
  DOUBLE PRECISION :: dist,dist_new,dist_beg
!
  INTEGER,        DIMENSION(nsize) :: ipoi
  DOUBLE COMPLEX, DIMENSION(nsize) :: alam_old,alam_new
!
  ierr=0
!
  dist_beg=MAX(MAXVAL(ABS(alam_old(:))),MAXVAL(ABS(alam_new(:))))
  DO i=1,nstiff
    dist=dist_beg
    DO j=1,nstiff
      dist_new=ABS(alam_old(i)-alam_new(j))
      IF(dist_new.LT.dist) THEN
        dist=dist_new
        ipoi(i)=j
      ENDIF
    ENDDO
  ENDDO
!
  DO i=1,nstiff
    DO j=1,nstiff
      IF(i.NE.j.AND.ipoi(i).EQ.ipoi(j)) THEN
        ierr=1
        RETURN
      ELSEIF(i.EQ.j.AND.ipoi(i).NE.j) THEN
        ierr=-1
      ENDIF
    ENDDO
  ENDDO
!
  END SUBROUTINE align_lambdas
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE sample_brhs(ierr)
!
  USE sample_matrix_mod
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: nder=0,iab=3
  INTEGER :: i,j,iter,npoi_old,iold,inew,ibeg,iend,nshift,npoilag,n1,n2,ierr
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: isplit
!
  DOUBLE PRECISION :: h
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xold
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef,amat_maxmod,w_maxmod
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: amat1,amat2
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: amat_old
!
  ierr=0
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  ALLOCATE(coef(0:nder,npoilag))
!
  x=xarr_rhs(1)
  CALL get_matrix(iab)
!
  n1=nsize
  n2=nrhs
!
  ALLOCATE(amat1(n1,n2),amat2(n1,n2),amat_maxmod(n1,n2),w_maxmod(n1,n2))
!
  wmat(:,:,1)=brhs
  w_maxmod=ABS(brhs)
!
  DO i=2,npoi_rhs
    x=xarr_rhs(i)
    CALL get_matrix(iab)
    wmat(:,:,i)=brhs
    w_maxmod=MAX(w_maxmod,ABS(brhs))
  ENDDO
!
! first check which intervals should be splitted
  ALLOCATE(isplit(npoi_rhs))
  isplit=0
  DO inew=1,npoi_rhs-1
    x=0.5d0*(xarr_rhs(inew)+xarr_rhs(inew+1))
    ibeg=MAX(1,MIN(npoi_rhs-nlagr-1,inew-nshift-1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr_rhs(ibeg:iend),coef)
    DO i=1,n1
      amat1(i,:)=MATMUL(wmat(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        amat_maxmod(i,j)=MAXVAL(ABS(wmat(i,j,ibeg:iend)))
      ENDDO
    ENDDO
    ibeg=MAX(2,MIN(npoi_rhs-nlagr,inew-nshift+1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr_rhs(ibeg:iend),coef)
    DO i=1,n1
      amat2(i,:)=MATMUL(wmat(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        amat_maxmod(i,j)=MAX(amat_maxmod(i,j),ABS(wmat(i,j,iend)),w_maxmod(i,j))
      ENDDO
    ENDDO
    DO i=1,n1
      DO j=1,n2
        IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*amat_maxmod(i,j)) isplit(inew)=1
      ENDDO
    ENDDO
  ENDDO
  IF(MAXVAL(isplit).GT.0) THEN
    npoi_old=npoi_rhs
    ALLOCATE(xold(npoi_rhs),amat_old(n1,n2,npoi_rhs))
    xold=xarr_rhs
    amat_old=wmat
  ELSE
    RETURN
  ENDIF
!
  iter=0
  DO
    iter=iter+1
    IF(iter.GT.itermax) THEN
      ierr=2
      PRINT *,'sample_brhs : maximum number of iterations exceeded'
      RETURN
    ENDIF
!
! determine the dimension of new arrays and re-allocate them:
    DO iold=1,npoi_old-1
      IF(isplit(iold).EQ.1) npoi_rhs=npoi_rhs+1
    ENDDO
    IF(ALLOCATED(xarr_rhs)) THEN
      DEALLOCATE(xarr_rhs,wmat)
    ENDIF
    ALLOCATE(xarr_rhs(npoi_rhs),wmat(n1,n2,npoi_rhs))
!
! fill new arrays:
    inew=0
    DO iold=1,npoi_old-1
      inew=inew+1
      xarr_rhs(inew)=xold(iold)
      wmat(:,:,inew)=amat_old(:,:,iold)
      IF(isplit(iold).EQ.1) THEN
        inew=inew+1
        x=0.5d0*(xold(iold)+xold(iold+1))
        CALL get_matrix(iab)
        xarr_rhs(inew)=x
        wmat(:,:,inew)=brhs
      ENDIF
    ENDDO
    inew=inew+1
    xarr_rhs(inew)=xold(npoi_old)
    wmat(:,:,inew)=amat_old(:,:,npoi_old)
    DEALLOCATE(isplit)
!
! check which intervals should be splitted
    ALLOCATE(isplit(npoi_rhs))
    isplit=0
    DO inew=1,npoi_rhs-1
      x=0.5d0*(xarr_rhs(inew)+xarr_rhs(inew+1))
      ibeg=MAX(1,MIN(npoi_rhs-nlagr-1,inew-nshift-1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr_rhs(ibeg:iend),coef)
      DO i=1,n1
        amat1(i,:)=MATMUL(wmat(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          amat_maxmod(i,j)=MAXVAL(ABS(wmat(i,j,ibeg:iend)))
        ENDDO
      ENDDO
      ibeg=MAX(2,MIN(npoi_rhs-nlagr,inew-nshift+1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr_rhs(ibeg:iend),coef)
      DO i=1,n1
        amat2(i,:)=MATMUL(wmat(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          amat_maxmod(i,j)=MAX(amat_maxmod(i,j),ABS(wmat(i,j,iend)),    &
                               w_maxmod(i,j))
        ENDDO
      ENDDO
      DO i=1,n1
        DO j=1,n2
          IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*amat_maxmod(i,j)) isplit(inew)=1
        ENDDO
      ENDDO
    ENDDO
    IF(MAXVAL(isplit).GT.0) THEN
      npoi_old=npoi_rhs
      DEALLOCATE(xold,amat_old)
      ALLOCATE(xold(npoi_rhs),amat_old(n1,n2,npoi_rhs))
      xold=xarr_rhs
      amat_old=wmat
    ELSE
      EXIT
    ENDIF
  ENDDO
!
  END SUBROUTINE sample_brhs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE convert_to_w
!
  USE sample_matrix_mod
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: nder=0
  INTEGER :: i,j,ibeg,iend,nshift,npoilag,ind
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coef
  DOUBLE COMPLEX,   DIMENSION(:,:), ALLOCATABLE :: phi_inv_loc,wnew
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  ALLOCATE(coef(0:nder,npoilag),phi_inv_loc(nsize,nsize),wnew(nsize,nrhs))
!
  DO i=1,npoi_rhs
!
    CALL binsrc(xarr,1,npoi,xarr_rhs(i),ind)
!
    ibeg=MAX(1,MIN(npoi-nlagr,ind-nshift-1))
    iend=ibeg+nlagr
!
    CALL plag_coeff(npoilag,nder,xarr_rhs(i),xarr(ibeg:iend),coef)
!
    DO j=1,nsize
      phi_inv_loc(j,:)=MATMUL(phi_inv(j,:,ibeg:iend),coef(0,:))
    ENDDO
!
    wnew=MATMUL(phi_inv_loc,wmat(:,:,i))
!
    wmat(:,:,i)=wnew
!
  ENDDO
!
  DEALLOCATE(coef,phi_inv_loc,wnew)
!
  END SUBROUTINE convert_to_w
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE order_eigvecs(nsize,npoi,xarr,alam,phi,ierr)
!
  IMPLICIT NONE
!
  INTEGER :: nsize,npoi,keep,ierr,i,j,info
!
  DOUBLE PRECISION, DIMENSION(npoi)               :: xarr
  DOUBLE COMPLEX,   DIMENSION(nsize,npoi)         :: alam
  DOUBLE COMPLEX,   DIMENSION(nsize,nsize,npoi)   :: phi
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: ipoi
  DOUBLE COMPLEX,   DIMENSION(:),     ALLOCATABLE :: alam_tmp
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: phi_tmp
!
  ALLOCATE(alam_tmp(nsize),phi_tmp(nsize,nsize),ipoi(nsize))
!
  keep=0
  DO i=2,npoi
    CALL align_lambdas(nsize,nsize,alam(:,i-1),alam(:,i),ipoi,info)
    IF(info.GT.0) THEN
      ierr=2
      PRINT *,'order_eigvecs error: cannot align eigenvalues'
    ELSEIF(info.EQ.0) THEN
      keep=i
      EXIT
    ENDIF
  ENDDO
!
  IF(keep.EQ.0) THEN
    ierr=2
    PRINT *,'order_eigvecs error: cannot align eigenvalues, jumpy'
  ENDIF
!
  DO i=keep+1,npoi
    alam_tmp=(alam(:,i-1)*(xarr(i)-xarr(i-2))-alam(:,i-2)        &
            *(xarr(i)-xarr(i-1)))/(xarr(i-1)-xarr(i-2))
!
    CALL align_lambdas(nsize,nsize,alam_tmp,alam(:,i),ipoi,info)
!
    IF(info.GT.0) THEN
      ierr=2
      PRINT *,'sample_eigvecs error: cannot align eigenvalues'
    ELSEIF(info.LT.0) THEN
      alam_tmp=alam(ipoi,i)
      alam(:,i)=alam_tmp
      phi_tmp=phi(:,ipoi,i)
      phi(:,:,i)=phi_tmp
    ENDIF
  ENDDO
!
  DO i=keep-2,1,-1
    alam_tmp=(alam(:,i+1)*(xarr(i)-xarr(i+2))-alam(:,i+2)        &
            *(xarr(i)-xarr(i+1)))/(xarr(i+1)-xarr(i+2))
    CALL align_lambdas(nsize,nsize,alam_tmp,alam(:,i),ipoi,info)
    IF(info.GT.0) THEN
      ierr=2
      PRINT *,'sample_eigvecs error: cannot align eigenvalues'
    ELSEIF(info.LT.0) THEN
      alam_tmp=alam(ipoi,i)
      alam(:,i)=alam_tmp
      phi_tmp=phi(:,ipoi,i)
      phi(:,:,i)=phi_tmp
    ENDIF
  ENDDO
!
  DO i=2,npoi
    DO j=1,nsize
      IF(REAL(SUM(CONJG(phi(:,j,i-1))*phi(:,j,i))).LT.0.d0) THEN
        phi(:,j,i)=-phi(:,j,i)
      ENDIF
    ENDDO
  ENDDO
!
  DEALLOCATE(alam_tmp,phi_tmp,ipoi)
!
  END SUBROUTINE order_eigvecs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
