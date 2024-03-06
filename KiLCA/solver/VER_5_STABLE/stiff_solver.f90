!
  SUBROUTINE stiff_solver(nsize_in,nlagr_in,nrhs_in,itermax_in,npoi_out,   &
                          eps_in,xbeg_in,xend_in,x_a,boumat,bouvec,antvec, &
                          n_inner,iswout,iblablabla,npoi_out_adapt,x,fun,ierr)
!
  USE sample_matrix_mod, ONLY : nlagr,nsize,npoi,itermax,eikonals,nstiff,  &
                                xbeg,xend,eps,xarr,amat_arr,bmat,phi,alam, &
                                phi_inv,npoi_rhs,xarr_rhs,wmat,nrhs,       &
                                x_ant,iflag_ant,antshift
!
  IMPLICIT NONE
!
  INTEGER,          PARAMETER :: nder=0,niter_max=100
  DOUBLE PRECISION, PARAMETER :: tolfix=2.d0, bignumber=1.d5
  INTEGER            :: nsize_in,npoi_out,nlagr_in,itermax_in,nstep,ierr,nrhs_in
  INTEGER            :: i,j,k,l,ind,nshift,npoilag,ibeg,iend,niter,iab,nstep_rhs
  INTEGER            :: nstencil,nreg,iregime,npoi_old,npoi_tot,ireg,iren
  INTEGER            :: istart,istop,npoi_rhs_tot,nstep_tot,nstep_rhs_tot
  INTEGER            :: n_inner,nren_tot,iren_tot,iren_rhs,n_outer,iblablabla
  INTEGER            :: npoi_out_adapt,iswout,j_ant
  DOUBLE PRECISION   :: eps_in,xbeg_in,xend_in,h_out,x_a,h,tol,tol_res
  DOUBLE PRECISION   :: alam_mild_glob,psi_mild,delpsi_mild,psi_renorm
  DOUBLE PRECISION   :: hh,bignumlog
  DOUBLE COMPLEX :: deltapsi,psibeg
  INTEGER,          DIMENSION(4) :: ind_ant
  DOUBLE PRECISION, DIMENSION(npoi_out)                           :: x
  DOUBLE COMPLEX,   DIMENSION(nsize_in)                           :: bouvec
  DOUBLE COMPLEX,   DIMENSION(nsize_in)                           :: antvec
  DOUBLE COMPLEX,   DIMENSION(nsize_in,nsize_in)                  :: boumat
  DOUBLE COMPLEX,   DIMENSION(nsize_in,nsize_in+nrhs_in,npoi_out) :: fun
!
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: nstiff_loc,ipoi
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: ibounds,nstiff_reg
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xfine,xfine_rhs,alam_mild
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef
  DOUBLE COMPLEX,   DIMENSION(:),     ALLOCATABLE :: psi_loc
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: Fbig_loc,phi_loc
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: Fbig,Fbig_rhs
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: Fbig_tot,Fbig_rhs_tot
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: Fbig_sol,fun_sol
!
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: ibeg_tot,iend_tot
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: ibeg_tmp,iend_tmp
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: ibeg_rhs,iend_rhs
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: ibst,iest,nrenorm
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: ibst_rhs,iest_rhs
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: iflag_a
  INTEGER,          DIMENSION(:,:),   ALLOCATABLE :: irenorm,irenorm_rhs
  INTEGER,          DIMENSION(:,:),   ALLOCATABLE :: ipoiren
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xarr_tot,xarr_tmp
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xarr_rhs_tot,xarr_rhs_tmp
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xfine_tot,xfine_rhs_tot
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: alam_mild_tot,xarr_sol
  DOUBLE COMPLEX,   DIMENSION(:),     ALLOCATABLE :: bouvec_eig,antvec_eig
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: alam_tot,alam_tmp
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: eikonals_tot
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: boumat_eig
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: dummy_in,dummy_out
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: dum_in,dum_out
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: startvecs,endvecs
  DOUBLE COMPLEX,   DIMENSION(:,:),   ALLOCATABLE :: roam_rhs
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: amat_tot,amat_tmp
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: phi_tot,phi_tmp
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: phi_inv_tot,Bmat_tot
  DOUBLE COMPLEX,   DIMENSION(:,:,:), ALLOCATABLE :: wmat_tot,wmat_tmp
  DOUBLE COMPLEX,   DIMENSION(:,:,:),   ALLOCATABLE :: F_in,F_out
  DOUBLE COMPLEX,   DIMENSION(:,:,:),   ALLOCATABLE :: Frhs_in,Frhs_out
  DOUBLE COMPLEX,   DIMENSION(:,:,:),   ALLOCATABLE :: c_forw,c_back
  DOUBLE COMPLEX,   DIMENSION(:,:,:),   ALLOCATABLE :: amplF,alpha,beta
  DOUBLE COMPLEX,   DIMENSION(:,:,:),   ALLOCATABLE :: add_rhs
!
  bignumlog=log(bignumber)
!
  nsize=nsize_in
  nlagr=nlagr_in
  nrhs=nrhs_in
  itermax=itermax_in
  xbeg=xbeg_in
  xend=xend_in
  eps=eps_in
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  ALLOCATE(coef(0:nder,npoilag))
!
  tol_res=eps**(-1.d0/(1.d0+nlagr))
  psi_renorm=LOG(1.d5)
!
! If iab=1 - sample_matrix samples only A,  adaptive source integration
! If iab=2 - sample_matrix samples A and B, source integrated on the same
!            grid as the homogeneous system
! If iab=3 - sample_matrix does not run (not allowed here)
!
  iab=1
!  iab=2
!
! Initial sampling - grid is adjusted for ODE matrix A :
!
  x_ant=x_a
  iflag_ant=1
  nrhs=0
!
  CALL sample_matrix(iab,ierr)
!
  IF(ierr.NE.0) then
    print *,'stiff_solver error: ierror in sample_matrix'
    call stiff_solver_deallocate
    RETURN
  ENDIF
!
  if(iblablabla.eq.1) print *,'sampling matrix done'
!
! Initial refinement : eigenvectors and eigenvalues are computed, grid
! is refined to sample eigenvectors properly :
!
  iregime=0
!
  CALL sample_eigvecs(iregime,iab,ierr)
!
  IF(ierr.NE.0) then
    print *,'stiff_solver error: error in sample_eigvecs'
    call stiff_solver_deallocate
    RETURN
  endif
!
  if(iblablabla.eq.1) print *,'sampling eigenvectors done'
!
  ALLOCATE(nstiff_loc(npoi),alam_mild(npoi))
  ALLOCATE(ipoi(nsize))
  ALLOCATE(ibounds(npoi),nstiff_reg(npoi))
!
! Splitting into subintervals (regions) with constant number of stiff modes :
!
  nstiff=nsize
!  tol=tolfix*2.d0
  tol=tolfix*4.d0
!
  CALL stiffness_estimator(nsize,npoi,nstiff,tol,xarr,alam, &
                           nstiff_loc,alam_mild)
!
  nstiff=MINVAL(nstiff_loc)
!
  nstencil=nlagr+2
!
  CALL  detstiffreg(npoi,nstiff_loc,nstencil,nreg,ibounds,nstiff_reg,ierr)
!
  npoi_tot=npoi+nreg-1
!
  ALLOCATE(ibeg_tot(nreg),iend_tot(nreg),ibeg_tmp(nreg),iend_tmp(nreg))
  ALLOCATE(ibst(nreg),iest(nreg))
  ALLOCATE(xarr_tot(npoi_tot),amat_tot(nsize,nsize,npoi_tot))
  ALLOCATE(phi_tot(nsize,nsize,npoi_tot),alam_tot(nsize,npoi_tot))
!
  IF(nrhs.NE.0) THEN
    ALLOCATE(ibeg_rhs(nreg),iend_rhs(nreg))
    ALLOCATE(ibst_rhs(nreg),iest_rhs(nreg))
    npoi_rhs_tot=npoi_tot
    ALLOCATE(xarr_rhs_tot(npoi_rhs_tot),wmat_tot(nsize,nrhs,npoi_rhs_tot))
  ENDIF
!
  DO i=1,nreg
    ibeg_tot(i)=ibounds(i)+i-1
    iend_tot(i)=ibounds(i+1)+i-1
    xarr_tot(ibeg_tot(i):iend_tot(i))=xarr(ibounds(i):ibounds(i+1))
    amat_tot(:,:,ibeg_tot(i):iend_tot(i))=amat_arr(:,:,ibounds(i):ibounds(i+1))
    CALL order_lambdas(nsize,alam(:,ibounds(i)),ipoi)
    DO j=1,nsize
      phi_tot(:,j,ibeg_tot(i):iend_tot(i))                        &
                  =phi(:,ipoi(j),ibounds(i):ibounds(i+1))
      alam_tot(j,ibeg_tot(i):iend_tot(i))                         &
                  =alam(ipoi(j),ibounds(i):ibounds(i+1))
    ENDDO
    IF(nrhs.ne.0.and.iab.eq.2) THEN
       wmat_tot(:,:,ibeg_tot(i):iend_tot(i))=wmat(:,:,ibounds(i):ibounds(i+1))
    ELSEIF(nrhs.ne.0.and.iab.eq.1) THEN
       wmat_tot(:,:,ibeg_tot(i):iend_tot(i))=(0.d0,0.d0)
    ENDIF
  ENDDO
!
  IF(nrhs.NE.0) THEN
    xarr_rhs_tot=xarr_tot
    ibeg_rhs=ibeg_tot
    iend_rhs=iend_tot
    IF(iab.EQ.1) ALLOCATE(xarr_rhs(1),wmat(1,1,1))
  ENDIF
!
! Additional refinement of the grid after splitting into subintervals, 
! formation of adaptive grid for the source problem in the case iab=1 :
!
  DO i=1,nreg
!
    DEALLOCATE(xarr,amat_arr,phi,alam)
    npoi=iend_tot(i)-ibeg_tot(i)+1
    ALLOCATE(xarr(npoi),amat_arr(nsize,nsize,npoi))
    ALLOCATE(alam(nsize,npoi),phi(nsize,nsize,npoi))
    xarr=xarr_tot(ibeg_tot(i):iend_tot(i))
    amat_arr=amat_tot(:,:,ibeg_tot(i):iend_tot(i))
    phi=phi_tot(:,:,ibeg_tot(i):iend_tot(i))
    alam=alam_tot(:,ibeg_tot(i):iend_tot(i))
!
    IF(nrhs.NE.0) THEN
      DEALLOCATE(xarr_rhs,wmat)
    ENDIF
!
    IF(nrhs.ne.0.and.iab.eq.2) THEN
      npoi_rhs=npoi
      ALLOCATE(xarr_rhs(npoi_rhs),wmat(nsize,nrhs,npoi_rhs))
      xarr_rhs=xarr
      wmat=wmat_tot(:,:,ibeg_tot(i):iend_tot(i))
    ENDIF
!
    npoi_old=npoi
    iregime=1
!
    CALL sample_eigvecs(iregime,iab,ierr)
!
    IF(nrhs.ne.0.and.iab.eq.1) THEN
      npoi_rhs=npoi
      ALLOCATE(xarr_rhs(npoi_rhs),wmat(nsize,nrhs,npoi_rhs))
      xarr_rhs=xarr
!
      CALL sample_brhs(ierr)
!
      IF(ierr.NE.0) then
        print *,'stiff_solver error: error in sample_brhs'
        call stiff_solver_deallocate
        RETURN
      endif
!
    ENDIF
!
    IF(npoi.NE.npoi_old) THEN
      ALLOCATE(xarr_tmp(npoi_tot),amat_tmp(nsize,nsize,npoi_tot))
      ALLOCATE(phi_tmp(nsize,nsize,npoi_tot),alam_tmp(nsize,npoi_tot))
      xarr_tmp=xarr_tot
      amat_tmp=amat_tot
      phi_tmp=phi_tot
      alam_tmp=alam_tot
      ibeg_tmp=ibeg_tot
      iend_tmp=iend_tot
      DEALLOCATE(xarr_tot,amat_tot,phi_tot,alam_tot)
      npoi_tot=npoi_tot+npoi-npoi_old
      ALLOCATE(xarr_tot(npoi_tot),amat_tot(nsize,nsize,npoi_tot))
      ALLOCATE(phi_tot(nsize,nsize,npoi_tot),alam_tot(nsize,npoi_tot))
      iend_tot(i)=iend_tot(i)+npoi-npoi_old
      DO j=i+1,nreg
        ibeg_tot(j)=ibeg_tot(j)+npoi-npoi_old
        iend_tot(j)=iend_tot(j)+npoi-npoi_old
      ENDDO
      DO j=1,i-1
        xarr_tot(ibeg_tot(j):iend_tot(j))=xarr_tmp(ibeg_tmp(j):iend_tmp(j))
        amat_tot(:,:,ibeg_tot(j):iend_tot(j))                                  &
                                         =amat_tmp(:,:,ibeg_tmp(j):iend_tmp(j))
        phi_tot(:,:,ibeg_tot(j):iend_tot(j))                                   &
                                         =phi_tmp(:,:,ibeg_tmp(j):iend_tmp(j))
        alam_tot(:,ibeg_tot(j):iend_tot(j))=alam_tmp(:,ibeg_tmp(j):iend_tmp(j))
      ENDDO
      xarr_tot(ibeg_tot(i):iend_tot(i))=xarr
      amat_tot(:,:,ibeg_tot(i):iend_tot(i))=amat_arr
      phi_tot(:,:,ibeg_tot(i):iend_tot(i))=phi
      alam_tot(:,ibeg_tot(i):iend_tot(i))=alam
      DO j=i+1,nreg
        xarr_tot(ibeg_tot(j):iend_tot(j))=xarr_tmp(ibeg_tmp(j):iend_tmp(j))
        amat_tot(:,:,ibeg_tot(j):iend_tot(j))                                  &
                                         =amat_tmp(:,:,ibeg_tmp(j):iend_tmp(j))
        phi_tot(:,:,ibeg_tot(j):iend_tot(j))                                   &
                                         =phi_tmp(:,:,ibeg_tmp(j):iend_tmp(j))
        alam_tot(:,ibeg_tot(j):iend_tot(j))=alam_tmp(:,ibeg_tmp(j):iend_tmp(j))
      ENDDO
      DEALLOCATE(xarr_tmp,amat_tmp,phi_tmp,alam_tmp)
    ENDIF
!
    IF(nrhs.NE.0) THEN
      ALLOCATE(xarr_tmp(npoi_rhs_tot),wmat_tmp(nsize,nrhs,npoi_rhs_tot))
      xarr_tmp=xarr_rhs_tot
      wmat_tmp=wmat_tot
      ibeg_tmp=ibeg_rhs
      iend_tmp=iend_rhs
      DEALLOCATE(xarr_rhs_tot,wmat_tot)
      npoi_rhs_tot=npoi_rhs_tot+npoi_rhs-npoi_old
      ALLOCATE(xarr_rhs_tot(npoi_rhs_tot),wmat_tot(nsize,nrhs,npoi_rhs_tot))
      iend_rhs(i)=iend_rhs(i)+npoi_rhs-npoi_old
      DO j=i+1,nreg
        ibeg_rhs(j)=ibeg_rhs(j)+npoi_rhs-npoi_old
        iend_rhs(j)=iend_rhs(j)+npoi_rhs-npoi_old
      ENDDO
      DO j=1,i-1
        xarr_rhs_tot(ibeg_rhs(j):iend_rhs(j))=xarr_tmp(ibeg_tmp(j):iend_tmp(j))
        wmat_tot(:,:,ibeg_rhs(j):iend_rhs(j))                                  &
                                         =wmat_tmp(:,:,ibeg_tmp(j):iend_tmp(j))
      ENDDO
      xarr_rhs_tot(ibeg_rhs(i):iend_rhs(i))=xarr_rhs
      wmat_tot(:,:,ibeg_rhs(i):iend_rhs(i))=wmat
      DO j=i+1,nreg
        xarr_rhs_tot(ibeg_rhs(j):iend_rhs(j))=xarr_tmp(ibeg_tmp(j):iend_tmp(j))
        wmat_tot(:,:,ibeg_rhs(j):iend_rhs(j))                                  &
                                         =wmat_tmp(:,:,ibeg_tmp(j):iend_tmp(j))
      ENDDO
      DEALLOCATE(xarr_tmp,wmat_tmp)
    ENDIF
!
  ENDDO
!
  DEALLOCATE(xarr,amat_arr,phi,alam,alam_mild)
  IF(nrhs.NE.0) THEN
    DEALLOCATE(xarr_rhs,wmat)
  ENDIF
!
! Transformation of the ODE set to the normal mode basis, integration 
! of eikonals, calcualtion of grid sizes for the residual ODE set:
!
  ALLOCATE(phi_inv_tot(nsize,nsize,npoi_tot))
  ALLOCATE(Bmat_tot(nsize,nsize,npoi_tot),alam_mild_tot(npoi_tot))
  ALLOCATE(eikonals_tot(nsize,npoi_tot),iflag_a(nreg))
!
  tol=tolfix
  nstep_tot=0
  nstep_rhs_tot=0
!
  if(iblablabla.eq.1) PRINT *,' '
!
  DO ireg=1,nreg
!
    niter=niter_max
    ibeg=ibeg_tot(ireg)
    iend=iend_tot(ireg)
    npoi=iend-ibeg+1
!
    CALL split_normal_modes(nlagr,nsize,npoi,niter,eps,xarr_tot(ibeg:iend),    &
                            amat_tot(:,:,ibeg:iend),nstiff_reg(ireg),          &
                            phi_tot(:,:,ibeg:iend),phi_inv_tot(:,:,ibeg:iend), &
                            alam_tot(:,ibeg:iend),Bmat_tot(:,:,ibeg:iend),ierr)
!
    DEALLOCATE(nstiff_loc)
    ALLOCATE(nstiff_loc(npoi))
!
    IF(ierr.EQ.0.and.iblablabla.eq.1) THEN
      PRINT *,'split_normal_modes: '
      PRINT *,'region ',ireg,': ',sngl(xarr_tot(ibeg)),' < x < ',              &
               sngl(xarr_tot(iend)),', number of points = ',npoi
      PRINT *,'system size = ',nsize,', number of stiff modes = ',             &
               nstiff_reg(ireg),', number of iterations = ',niter
      PRINT *,' '
    ELSEIF(ierr.NE.0) then
      PRINT *,'split_normal_modes: error = ',ierr
      PRINT *,' '
      call stiff_solver_deallocate
      RETURN
    ENDIF
!
    CALL stiffness_estimator(nsize,npoi,nstiff_reg(ireg),tol,                  &
                             xarr_tot(ibeg:iend),alam_tot(:,ibeg:iend),        &
                             nstiff_loc,alam_mild_tot(ibeg:iend))
!
    IF(MINVAL(nstiff_loc).LT.nstiff_reg(ireg)) THEN
      PRINT *,'stiffness_estimator : nstiff changed after refinement'
      ierr=1
      call stiff_solver_deallocate
      RETURN
    ENDIF
!
    nstiff=nstiff_reg(ireg)
    ALLOCATE(xarr(npoi),bmat(nsize,nsize,npoi),eikonals(nstiff,npoi))
    ALLOCATE(phi_inv(nsize,nsize,npoi),alam_mild(npoi))
    xarr=xarr_tot(ibeg:iend)
    bmat=bmat_tot(:,:,ibeg:iend)
    phi_inv=phi_inv_tot(:,:,ibeg:iend)
    alam_mild=alam_mild_tot(ibeg:iend)
!
    IF(nstiff.GT.0) THEN
!
      CALL integrate_eikonals(eps)
!
      eikonals_tot(1:nstiff,ibeg:iend)=eikonals(1:nstiff,:)
!
    ENDIF
!
    IF(nstiff.LT.nsize) THEN
      nstep=0
      iflag_a(ireg)=0
      DO i=2,npoi
        delpsi_mild=(xarr(i)-xarr(i-1))*MAX(alam_mild(i-1),alam_mild(i))
        k=CEILING(delpsi_mild*tol_res)
        if(k.gt.1) then
          j=2
          do while(j.lt.k)
            j=j*2
          enddo
          k=j
        endif
        nstep=nstep+k
        if(iflag_ant.eq.1) then
          if(i.eq.2.and.x_ant.gt.xarr(1).and.                             &
             x_ant-xarr(1).lt.antshift*(xarr(2)-xarr(1))/k) then
            nstep=nstep+1
            iflag_a(ireg)=-1
          elseif(i.eq.npoi.and.x_ant.lt.xarr(npoi).and.                   &
                 xarr(npoi)-x_ant.lt.antshift*(xarr(i)-xarr(i-1))/k) then
            nstep=nstep+1
            iflag_a(ireg)=1
          elseif(k.eq.1.and.x_ant.gt.xarr(i-1).and.x_ant.lt.xarr(i)) then
            nstep=nstep+1
          endif
        endif
      ENDDO
    ELSE
      if(nrhs.eq.0.and.x_a.ge.xarr(1).and.x_a.le.xarr(npoi) &
                  .and.iflag_ant.eq.1) then
        nstep=npoi
      else
        nstep=npoi-1
      endif
    ENDIF
    ibst(ireg)=nstep_tot+1
    nstep_tot=nstep_tot+nstep+1
    iest(ireg)=nstep_tot
!
    IF(nrhs.NE.0) THEN
      ibeg=ibeg_rhs(ireg)
      iend=iend_rhs(ireg)
      npoi_rhs=iend-ibeg+1
      ALLOCATE(xarr_rhs(npoi_rhs),wmat(nsize,nrhs,npoi_rhs))
      xarr_rhs=xarr_rhs_tot(ibeg:iend)
      wmat=wmat_tot(:,:,ibeg:iend)
!
      CALL convert_to_w
!
      wmat_tot(:,:,ibeg:iend)=wmat
      IF(nstiff.LT.nsize) THEN
        nstep_rhs=0
        i=2
        DO j=2,npoi_rhs
          IF(xarr_rhs(j).GT.xarr(i)) i=i+1
          delpsi_mild=(xarr(i)-xarr(i-1))*MAX(alam_mild(i-1),alam_mild(i))
          k=CEILING(delpsi_mild*tol_res)
          if(k.gt.1) then
            l=2
            do while(l.lt.k)
              l=l*2
            enddo
            k=l
          endif
          l=nint((xarr(i)-xarr(i-1))/(xarr_rhs(j)-xarr_rhs(j-1)))
          k=max(1,k/l)
          nstep_rhs=nstep_rhs+k
        ENDDO
      ELSE
        nstep_rhs=npoi_rhs-1
      ENDIF
      ibst_rhs(ireg)=nstep_rhs_tot+1
      nstep_rhs_tot=nstep_rhs_tot+nstep_rhs+1
      iest_rhs(ireg)=nstep_rhs_tot
      DEALLOCATE(xarr_rhs,wmat)
    ENDIF
!
    DEALLOCATE(xarr,Bmat,eikonals,phi_inv,alam_mild)
!
  ENDDO
!
! Formation of grids for the residual ODE set, splitting the grids
! into renormalization subintevals:
!
  ALLOCATE(xfine_tot(nstep_tot))
  ALLOCATE(irenorm(0:MAXVAL(iest-ibst),nreg),nrenorm(nreg))
  irenorm(0,:)=-1
  if(nrhs.ne.0) then
    ALLOCATE(xfine_rhs_tot(nstep_rhs_tot))
    ALLOCATE(irenorm_rhs(0:MAXVAL(iest-ibst),nreg))
    irenorm_rhs(0,:)=-1
  endif
!
  DO ireg=1,nreg
!
    ibeg=ibeg_tot(ireg)
    iend=iend_tot(ireg)
    npoi=iend-ibeg+1
    nstiff=nstiff_reg(ireg)
    ALLOCATE(xarr(npoi),alam_mild(npoi))
    xarr=xarr_tot(ibeg:iend)
    alam_mild=alam_mild_tot(ibeg:iend)
!
    IF(nstiff.LT.nsize) THEN
!
      nstep=iest(ireg)-ibst(ireg)
      ALLOCATE(xfine(0:nstep))
!
      nstep=0
      nrenorm(ireg)=0
      iren=1
      psi_mild=0.d0
      xfine(0)=xarr(1)
      if(iflag_a(ireg).eq.-1) then
        nstep=nstep+1
        xfine(nstep)=x_ant
        nrenorm(ireg)=nrenorm(ireg)+1
        irenorm(nrenorm(ireg),ireg)=nstep
        ind_ant(1)=ireg
        ind_ant(2)=nrenorm(ireg)
        ind_ant(3)=nstep
      endif
      DO i=2,npoi
        delpsi_mild=(xarr(i)-xarr(i-1))*MAX(alam_mild(i-1),alam_mild(i))
        k=CEILING(delpsi_mild*tol_res)
        if(k.gt.1) then
          j=2
          do while(j.lt.k)
            j=j*2
          enddo
          k=j
        endif
        h=(xarr(i)-xarr(i-1))/k
        delpsi_mild=delpsi_mild/k
        if(iflag_ant.eq.1.and.x_a.gt.xarr(i-1).and.x_a.lt.xarr(i)) then
          if(i.eq.npoi.and.iflag_a(ireg).eq.1) then
            do j=1,k-1
              nstep=nstep+1
              xfine(nstep)=xarr(i-1)+h*j
              psi_mild=psi_mild+delpsi_mild
              IF(psi_mild.GT.psi_renorm) THEN
                nrenorm(ireg)=nrenorm(ireg)+1
                irenorm(nrenorm(ireg),ireg)=nstep
                psi_mild=0.d0
              ENDIF
            enddo
            nstep=nstep+1
            xfine(nstep)=x_ant
            nrenorm(ireg)=nrenorm(ireg)+1
            irenorm(nrenorm(ireg),ireg)=nstep
            psi_mild=0.d0
            ind_ant(1)=ireg
            ind_ant(2)=nrenorm(ireg)
            ind_ant(3)=nstep
            nstep=nstep+1
            xfine(nstep)=xarr(i)
            psi_mild=psi_mild+delpsi_mild
            IF(psi_mild.GT.psi_renorm) THEN
              nrenorm(ireg)=nrenorm(ireg)+1
              irenorm(nrenorm(ireg),ireg)=nstep
              psi_mild=0.d0
            ENDIF
          elseif(k.eq.1) then
            nstep=nstep+1
            xfine(nstep)=x_ant
            nrenorm(ireg)=nrenorm(ireg)+1
            irenorm(nrenorm(ireg),ireg)=nstep
            psi_mild=0.d0
            ind_ant(1)=ireg
            ind_ant(2)=nrenorm(ireg)
            ind_ant(3)=nstep
            nstep=nstep+1
            xfine(nstep)=xarr(i)
            psi_mild=psi_mild+delpsi_mild
          else
            l=max(1,int((x_ant-xarr(i-1))/h))
            hh=(x_ant-xarr(i-1))/l
            do j=1,l
              nstep=nstep+1
              xfine(nstep)=xarr(i-1)+hh*j
              psi_mild=psi_mild+delpsi_mild*hh/h
              IF(psi_mild.GT.psi_renorm) THEN
                nrenorm(ireg)=nrenorm(ireg)+1
                irenorm(nrenorm(ireg),ireg)=nstep
                psi_mild=0.d0
              ENDIF
            enddo
            if(irenorm(nrenorm(ireg),ireg).ne.nstep) then
              nrenorm(ireg)=nrenorm(ireg)+1
              irenorm(nrenorm(ireg),ireg)=nstep
              psi_mild=0.d0
            endif
            ind_ant(1)=ireg
            ind_ant(2)=nrenorm(ireg)
            ind_ant(3)=nstep
            k=k-l
            hh=(xarr(i)-x_ant)/k
            do j=1,k-1
              nstep=nstep+1
              xfine(nstep)=x_ant+h*j
              psi_mild=psi_mild+delpsi_mild*hh/h
              IF(psi_mild.GT.psi_renorm) THEN
                nrenorm(ireg)=nrenorm(ireg)+1
                irenorm(nrenorm(ireg),ireg)=nstep
                psi_mild=0.d0
              ENDIF
            enddo
            nstep=nstep+1
            xfine(nstep)=xarr(i)
            psi_mild=psi_mild+delpsi_mild*hh/h
            IF(psi_mild.GT.psi_renorm) THEN
              nrenorm(ireg)=nrenorm(ireg)+1
              irenorm(nrenorm(ireg),ireg)=nstep
              psi_mild=0.d0
            ENDIF
          endif
        else
          do j=1,k-1
            nstep=nstep+1
            xfine(nstep)=xarr(i-1)+h*j
            psi_mild=psi_mild+delpsi_mild
            IF(psi_mild.GT.psi_renorm) THEN
              nrenorm(ireg)=nrenorm(ireg)+1
              irenorm(nrenorm(ireg),ireg)=nstep
              psi_mild=0.d0
            ENDIF
          enddo
          nstep=nstep+1
          xfine(nstep)=xarr(i)
          psi_mild=psi_mild+delpsi_mild
          if(xfine(nstep).eq.x_ant) then
            nrenorm(ireg)=nrenorm(ireg)+1
            irenorm(nrenorm(ireg),ireg)=nstep
            psi_mild=0.d0
            ind_ant(1)=ireg
            ind_ant(2)=nrenorm(ireg)
            ind_ant(3)=nstep
          endif
          IF(psi_mild.GT.psi_renorm) THEN
            nrenorm(ireg)=nrenorm(ireg)+1
            irenorm(nrenorm(ireg),ireg)=nstep
            psi_mild=0.d0
          ENDIF
        endif
      ENDDO
      if(irenorm(nrenorm(ireg),ireg).lt.nstep) then
        nrenorm(ireg)=nrenorm(ireg)+1
        irenorm(nrenorm(ireg),ireg)=nstep
      endif
!
      xfine_tot(ibst(ireg):iest(ireg))=xfine
!
    ELSE
!
      if(nrhs.eq.0.and.x_a.ge.xarr(1).and.x_a.le.xarr(npoi)) then
        ind_ant(1)=ireg
        ind_ant(2)=1
        nrenorm(ireg)=2
        if(iflag_ant.eq.0) then
          do i=1,npoi
            if(xarr(i).eq.x_a) then
              irenorm(1,ireg)=i-1
              ind_ant(3)=i-1
              cycle
            endif
          enddo
          irenorm(2,ireg)=npoi-1
          xfine_tot(ibst(ireg):iest(ireg))=xarr
        else
          k=0
          do i=1,npoi-1
            xfine_tot(ibst(ireg)+k)=xarr(i)
            k=k+1
            if(x_a.gt.xarr(i).and.x_a.lt.xarr(i+1)) then
              xfine_tot(ibst(ireg)+k)=x_a
              irenorm(1,ireg)=k
              ind_ant(3)=k
              k=k+1
            endif
          enddo
          irenorm(2,ireg)=k
          xfine_tot(ibst(ireg)+k)=xarr(npoi)
        endif
      else
        xfine_tot(ibst(ireg):iest(ireg))=xarr
        nrenorm(ireg)=1
        irenorm(1,ireg)=npoi-1
      endif
!
    ENDIF
!
    IF(nrhs.NE.0) THEN
!
      ibeg=ibeg_rhs(ireg)
      iend=iend_rhs(ireg)
      npoi_rhs=iend-ibeg+1
      ALLOCATE(xarr_rhs(npoi_rhs))
      nstep_rhs=iest_rhs(ireg)-ibst_rhs(ireg)
      ALLOCATE(xfine_rhs(0:nstep_rhs))
      xarr_rhs=xarr_rhs_tot(ibeg:iend)
!
      IF(nstiff.LT.nsize) THEN
        nstep_rhs=0
        xfine_rhs(0)=xarr_rhs(1)
        iren=1
        i=2
        DO j=2,npoi_rhs
          IF(xarr_rhs(j).GT.xarr(i)) i=i+1
          delpsi_mild=(xarr(i)-xarr(i-1))*MAX(alam_mild(i-1),alam_mild(i))
          k=CEILING(delpsi_mild*tol_res)
          if(k.gt.1) then
            l=2
            do while(l.lt.k)
              l=l*2
            enddo
            k=l
          endif
          l=nint((xarr(i)-xarr(i-1))/(xarr_rhs(j)-xarr_rhs(j-1)))
          k=max(1,k/l)
          h=(xarr_rhs(j)-xarr_rhs(j-1))/k
          DO l=1,k-1
            nstep_rhs=nstep_rhs+1
            xfine_rhs(nstep_rhs)=xarr_rhs(j-1)+h*l
            if(abs(xfine_rhs(nstep_rhs)-xfine(irenorm(iren,ireg))) .lt. &
               100.d0*epsilon(1.d0)*abs(xfine_rhs(nstep_rhs))) then
              irenorm_rhs(iren,ireg)=nstep_rhs
              iren=iren+1
            endif
          ENDDO
          nstep_rhs=nstep_rhs+1
          xfine_rhs(nstep_rhs)=xarr_rhs(j)
          if(abs(xfine_rhs(nstep_rhs)-xfine(irenorm(iren,ireg))) .lt. &
             100.d0*epsilon(1.d0)*abs(xfine_rhs(nstep_rhs))) then
            irenorm_rhs(iren,ireg)=nstep_rhs
            iren=iren+1
          endif
        ENDDO
        if(iren-1.ne.nrenorm(ireg)) then
          print *,'stiff_solver error: mismatch between the grids'
          ierr=1
          call stiff_solver_deallocate
          return
        endif
      ELSE
        xfine_rhs=xarr_rhs
        irenorm_rhs(1,ireg)=npoi_rhs
      ENDIF
!
      xfine_rhs_tot(ibst_rhs(ireg):iest_rhs(ireg))=xfine_rhs
!
      DEALLOCATE(xarr_rhs,xfine_rhs)
!
    ENDIF
!
    if(nstiff.LT.nsize) DEALLOCATE(xfine)
    DEALLOCATE(xarr,alam_mild)
!
    if(iblablabla.eq.1) then
      print *,'grid formation:'
      if(nrhs.ne.0) then
        print *,'region ',ireg,' : nstep = ',nstep,',  nstep_rhs = ',nstep_rhs &
               ,', nrenorm = ',nrenorm(ireg)
        print *,' '
      else
        print *,'region ',ireg,' : nstep = ',nstep,',  nstep_rhs = ',nstep     &
               ,', nrenorm = ',nrenorm(ireg)
        print *,' '
      endif
    endif
!
  ENDDO
!
! Joint renormalization pointer:
!
  nren_tot=0
  do ireg=1,nreg
    nren_tot=nren_tot+nrenorm(ireg)
  enddo
  allocate(ipoiren(2,nren_tot))
  nren_tot=0
  do ireg=1,nreg
    do iren=1,nrenorm(ireg)
      nren_tot=nren_tot+1
      ipoiren(1,nren_tot)=ireg
      ipoiren(2,nren_tot)=iren
      if(ind_ant(1).eq.ireg.and.ind_ant(2).eq.iren) ind_ant(4)=nren_tot
    enddo
  enddo
!
! Integration of the residual ODE set:
!
  ALLOCATE(Fbig_tot(nsize,nsize,nstep_tot))
  allocate(F_in(nsize,nsize,nren_tot))
  allocate(F_out(nsize,nsize,nren_tot))
  if(nrhs.ne.0) then
    ALLOCATE(Fbig_rhs_tot(nsize,nrhs,nstep_rhs_tot))
    allocate(Frhs_in(nsize,nrhs,nren_tot))
    allocate(Frhs_out(nsize,nrhs,nren_tot))
  endif
  allocate(boumat_eig(nsize,nsize),bouvec_eig(nsize),antvec_eig(nsize))
!
  iren_tot=0
  iren_rhs=0
!
  DO ireg=1,nreg
!
    ibeg=ibeg_tot(ireg)
    iend=iend_tot(ireg)
    npoi=iend-ibeg+1
    nstiff=nstiff_reg(ireg)
    ALLOCATE(xarr(npoi),bmat(nsize,nsize,npoi),eikonals(nstiff,npoi))
    ALLOCATE(phi_inv(nsize,nsize,npoi),alam_mild(npoi))
    xarr=xarr_tot(ibeg:iend)
    bmat=bmat_tot(:,:,ibeg:iend)
    phi_inv=phi_inv_tot(:,:,ibeg:iend)
    eikonals(1:nstiff,:)=eikonals_tot(1:nstiff,ibeg:iend)
    alam_mild=alam_mild_tot(ibeg:iend)
!
    nstep=iest(ireg)-ibst(ireg)
    ALLOCATE(xfine(0:nstep),Fbig(nsize,nsize,0:nstep))
    xfine=xfine_tot(ibst(ireg):iest(ireg))
!
    if(ind_ant(1).eq.ireg) then
!
! Translation of jump conditions at the antenna into jump conditions for 
! F-functions:
!
      CALL binsrc(xarr,1,npoi,x_a,ind)
!
      ibeg=MAX(1,MIN(npoi-nlagr,ind-nshift-1))
      iend=ibeg+nlagr
!
      CALL plag_coeff(npoilag,nder,x_a,xarr(ibeg:iend),coef)
!
      DO k=1,nsize
        antvec_eig(k)=sum(MATMUL(phi_inv(k,:,ibeg:iend),coef(0,:))*antvec)
      ENDDO
!
    endif
!
    istart=0
!
    DO iren=1,nrenorm(ireg)
!
      istop=irenorm(iren,ireg)
!
!
      CALL residual_system(nlagr,nsize,npoi,nstiff,nstep,eps,xarr,     &
                           Bmat,eikonals,xfine(0:nstep),istart,istop,  &
                             iflag_a(ireg),Fbig(:,1:nsize,0:nstep),ierr)
!
      IF(ierr.NE.0)  then
        print *,'stiff_solver error: error in residual_system'
        call stiff_solver_deallocate
        RETURN
      ENDIF
!
      iren_tot=iren_tot+1
      F_in(:,:,iren_tot)=Fbig(:,:,istart)
!
      do i=1,nstiff
        deltapsi=Fbig(i,i,istop)-Fbig(i,i,istart)
        if(real(deltapsi).gt.bignumlog) deltapsi=bignumlog
        psibeg=Fbig(i,i,istop)-deltapsi
        Fbig(i,i,istart:istop)=exp(Fbig(i,i,istart:istop)-psibeg)
        F_in(i,i,iren_tot)=(1.d0,0.d0)
      enddo
!
      F_out(:,:,iren_tot)=Fbig(:,:,istop)
      Fbig_tot(:,:,ibst(ireg)+istart:ibst(ireg)+istop)=Fbig(:,:,istart:istop)
!
      istart=istop
!
    ENDDO
!
    DEALLOCATE(xfine,Fbig)
!
    IF(nrhs.NE.0) THEN
!
      ibeg=ibeg_rhs(ireg)
      iend=iend_rhs(ireg)
      npoi_rhs=iend-ibeg+1
      ALLOCATE(xarr_rhs(npoi_rhs),wmat(nsize,nrhs,npoi_rhs))
      nstep_rhs=iest_rhs(ireg)-ibst_rhs(ireg)
      ALLOCATE(xfine_rhs(0:nstep_rhs),Fbig_rhs(nsize,nrhs,0:nstep_rhs))
      xarr_rhs=xarr_rhs_tot(ibeg:iend)
      wmat=wmat_tot(:,:,ibeg:iend)
!
      xfine_rhs=xfine_rhs_tot(ibst_rhs(ireg):iest_rhs(ireg))
!
      istart=0
!
      DO iren=1,nrenorm(ireg)
!
        istop=irenorm_rhs(iren,ireg)
!
        CALL residual_system_inhom(nlagr,nsize,nrhs,npoi,npoi_rhs,nstiff,      &
                                   nstep_rhs,eps,xarr,xarr_rhs,Bmat,wmat,      &
                                   eikonals,xfine_rhs,istart,istop,            &
                                   Fbig_rhs(:,:,:),ierr)
!
        IF(ierr.NE.0) then
          print *,'stiff_solver error: error in residual_system_inhom'
          call stiff_solver_deallocate
          RETURN
        endif
!
        Fbig_rhs_tot(:,:,ibst_rhs(ireg)+istart:ibst_rhs(ireg)+istop)           &
                                                 =Fbig_rhs(:,:,istart:istop)
!
        iren_rhs=iren_rhs+1
        Frhs_in(:,:,iren_rhs)=Fbig_rhs(:,:,istart)
        Frhs_out(:,:,iren_rhs)=Fbig_rhs(:,:,istop)
!
        istart=istop
!
      ENDDO
!
      DEALLOCATE(xarr_rhs,wmat,xfine_rhs,Fbig_rhs)
!
    ENDIF
!
    DEALLOCATE(xarr,bmat,eikonals,phi_inv,alam_mild)
!
  ENDDO
!
! Translation of boundary conditions to the conditions for solution amplitudes:
!
  n_outer=nsize-n_inner
!
  allocate(startvecs(nsize,n_outer),endvecs(nsize,n_inner))
  allocate(dummy_in(nsize,nsize),dummy_out(nsize,nsize))
  allocate(dum_in(nsize,n_outer),dum_out(nsize,n_outer))
!
! Starting vectors at the inner boundary:
!
  do i=1,n_inner
    boumat_eig(i,:)=matmul(boumat(i,:),phi_tot(:,:,1))
  enddo
!
  call null_space(n_inner,nsize,boumat_eig(1:n_inner,:),dummy_out,ierr)
!
  if(ierr.ne.0) then
    print *,'stiff_solver error: cannot form the starting vector'
    call stiff_solver_deallocate
    return
  endif
!
  startvecs=dummy_out(:,1:n_outer)
!
  dummy_in=F_in(:,:,1)
!
  call invert_cmat(nsize,dummy_in,dummy_out,ierr)
!
  if(ierr.ne.0) then
    print *,'stiff_solver error: cannot invert the solution set'
    call stiff_solver_deallocate
    return
  endif
!
  startvecs=matmul(dummy_out,startvecs)
!
  if(nrhs.eq.0) then
!
! Inhomogeneous solution in case of drive by the antenna:
!
    nrhs=1
    nstep_rhs_tot=nstep_tot
    npoi_rhs_tot=npoi_tot
    ALLOCATE(ibeg_rhs(nreg),iend_rhs(nreg))
    ALLOCATE(ibst_rhs(nreg),iest_rhs(nreg))
    ALLOCATE(xarr_rhs_tot(npoi_rhs_tot),xfine_rhs_tot(nstep_rhs_tot))
    ALLOCATE(Fbig_rhs_tot(nsize,nrhs,nstep_rhs_tot))
    allocate(Frhs_in(nsize,nrhs,nren_tot))
    allocate(Frhs_out(nsize,nrhs,nren_tot))
    ibeg_rhs=ibeg_tot
    iend_rhs=iend_tot
    ibst_rhs=ibst
    iest_rhs=iest
    xarr_rhs_tot=xarr_tot
    xfine_rhs_tot=xfine_tot
    Fbig_rhs_tot=(0.d0,0.d0)
    Frhs_in=(0.d0,0.d0)
    Frhs_out=(0.d0,0.d0)
!
! Antenna jump:
!
    iren_tot=ind_ant(4)
    ireg=ind_ant(1)
    Frhs_out(:,1,iren_tot)=antvec_eig
!
  endif
!
! Coupling of renormalization intervals:
!
  allocate(c_forw(nsize,nsize,nren_tot-1))
!
  do iren_tot=1,nren_tot-1
!
    dummy_in=F_in(:,:,iren_tot+1)
!
    if(ipoiren(1,iren_tot+1).ne.ipoiren(1,iren_tot)) then
      ibeg=ibeg_tot(ipoiren(1,iren_tot+1))
      iend=iend_tot(ipoiren(1,iren_tot))
      dummy_in=matmul(matmul(phi_inv_tot(:,:,iend),phi_tot(:,:,ibeg)),dummy_in)
    endif
!
    call invert_cmat(nsize,dummy_in,dummy_out,ierr)
!
    if(ierr.ne.0) then
      print *,'stiff_solver error: cannot invert the solution set'
      call stiff_solver_deallocate
      return
    endif
!
    c_forw(:,:,iren_tot)=dummy_out
!
  enddo
!
! Renormalization:
!
  allocate(amplF(nsize,n_outer,nren_tot),alpha(n_outer,n_outer,nren_tot-1))
  allocate(beta(n_outer,nrhs,nren_tot-1))
  allocate(roam_rhs(nsize,nrhs),add_rhs(nsize,nrhs,nren_tot))
!
  dummy_in=F_in(:,:,1)
!
  call invert_cmat(nsize,dummy_in,dummy_out,ierr)
!
  if(ierr.ne.0) then
    print *,'stiff_solver error: cannot invert the solution set'
    call stiff_solver_deallocate
    return
  endif
!
  add_rhs(:,:,1)=-matmul(dummy_out,Frhs_in(:,:,1))
!
! Forward run:
!
  amplF(:,:,1)=startvecs
!
  do iren_tot=1,nren_tot-1
!
    dum_in=matmul(F_out(:,:,iren_tot),amplF(:,:,iren_tot))
!
    call qr_factorize(nsize,n_outer,dum_in,dummy_out,dum_out,ierr)
!
    if(ierr.ne.0) then
      print *,'stiff_solver error: QR factorization failed'
      call stiff_solver_deallocate
      return
    endif
!
    alpha(:,:,iren_tot)=dum_out(1:n_outer,:)
    amplF(:,:,iren_tot+1)=matmul(c_forw(:,:,iren_tot),dummy_out(:,1:n_outer))
    roam_rhs=Frhs_out(:,:,iren_tot)                                          &
            +matmul(F_out(:,:,iren_tot),add_rhs(:,:,iren_tot))
    beta(:,:,iren_tot)=matmul(transpose(conjg(dummy_out(:,1:n_outer))),roam_rhs)
    roam_rhs=roam_rhs-matmul(dummy_out(:,1:n_outer),beta(:,:,iren_tot)) 
!
    dummy_in=F_in(:,:,iren_tot+1)
!
    call invert_cmat(nsize,dummy_in,dummy_out,ierr)
!
    if(ierr.ne.0) then
      print *,'stiff_solver error: cannot invert the solution set'
      call stiff_solver_deallocate
      return
    endif
!
    add_rhs(:,:,iren_tot+1)=matmul(c_forw(:,:,iren_tot),roam_rhs)           &
                           -matmul(dummy_out,Frhs_in(:,:,iren_tot+1))
!
  enddo
!
  deallocate(dum_in,dum_out)
  allocate(dum_in(n_outer,n_outer),dum_out(n_outer,nrhs))
!
! Conditions at the outer boundary:
!
  do i=n_inner+1,nsize
    boumat_eig(i,:)=matmul(boumat(i,:),phi_tot(:,:,npoi_tot))
  enddo
!
  roam_rhs=Frhs_out(:,:,nren_tot)                                           &
          +matmul(F_out(:,:,nren_tot),add_rhs(:,:,nren_tot))
!
  dum_in=matmul(boumat_eig(n_inner+1:nsize,:),                              &
         matmul(F_out(:,:,nren_tot),amplF(:,:,nren_tot)))
!
  call invert_cmat(n_outer,dum_in,dum_in,ierr)
!
  if(ierr.ne.0) then
    print *,'stiff_solver error: cannot solve the outer boundary condition'
    call stiff_solver_deallocate
    return
  endif
!
  dum_out=-matmul(dum_in,matmul(boumat_eig(n_inner+1:nsize,:),roam_rhs))
!
  add_rhs(:,:,nren_tot)=add_rhs(:,:,nren_tot) &
                       +matmul(amplF(:,:,nren_tot),dum_out)
!
! Backward run
!
  do iren_tot=nren_tot-1,1,-1
!
    dum_out=dum_out-beta(:,:,iren_tot)
!
    dum_in=alpha(:,:,iren_tot)
!
    call solve_cutriasys(n_outer,nrhs,dum_in,dum_out,ierr)
!
    if(ierr.ne.0) then
      print *,'stiff_solver error: inversion of LT matrix failed'
      call stiff_solver_deallocate
      return
    endif
!
    add_rhs(:,:,iren_tot)=add_rhs(:,:,iren_tot) &
                         +matmul(amplF(:,:,iren_tot),dum_out)
!
  enddo
!
  allocate(Fbig_sol(nsize,nrhs,nstep_tot),phi_loc(nsize,nsize))
  allocate(xarr_sol(nstep_tot+1),fun_sol(nsize,nrhs,nstep_tot+1))
!
  j=0
  ireg=0
  do iren_tot=1,nren_tot
    if(ipoiren(1,iren_tot).ne.ireg) then
      ireg=ipoiren(1,iren_tot)
      nstiff=nstiff_reg(ireg)
      istart=0
    endif
    iren=ipoiren(2,iren_tot)
    istop=irenorm(iren,ireg)-1
    if(iren_tot.eq.nren_tot) istop=istop+1
    do i=ibst(ireg)+istart,ibst(ireg)+istop
      Fbig_sol(:,:,i)=matmul(Fbig_tot(:,:,i),add_rhs(:,:,iren_tot))
      npoi=iend_tot(ireg)-ibeg_tot(ireg)+1
!
      CALL binsrc(xarr_tot(ibeg_tot(ireg):iend_tot(ireg)),     &
                  1,npoi,xfine_tot(i),ind)
!
      ibeg=MAX(1,MIN(npoi-nlagr,ind-nshift-1))+ibeg_tot(ireg)-1
      iend=ibeg+nlagr
!
      CALL plag_coeff(npoilag,nder,xfine_tot(i),xarr_tot(ibeg:iend),coef)
!
      DO k=1,nsize
        phi_loc(k,:)=MATMUL(phi_tot(k,:,ibeg:iend),coef(0,:))
      ENDDO
!
      j=j+1
      xarr_sol(j)=xfine_tot(i)
      fun_sol(:,:,j)=matmul(phi_loc,Fbig_sol(:,:,i))
!
!write(1000,*) xfine_tot(i),dble(Fbig_sol(:,:,i)) &
!                          ,dimag(Fbig_sol(:,:,i))
!
    enddo
    istart=istop+1
    if(iren_tot.eq.ind_ant(4)) then
      i=ibst(ireg)+istart
      Fbig_sol(:,:,i)=matmul(F_out(:,:,iren_tot),add_rhs(:,:,iren_tot))
!
      CALL binsrc(xarr_tot(ibeg_tot(ireg):iend_tot(ireg)),     &
                  1,npoi,xfine_tot(i),ind)
!
      ibeg=MAX(1,MIN(npoi-nlagr,ind-nshift-1))+ibeg_tot(ireg)-1
      iend=ibeg+nlagr
!
      CALL plag_coeff(npoilag,nder,xfine_tot(i),xarr_tot(ibeg:iend),coef)
!
      DO k=1,nsize
        phi_loc(k,:)=MATMUL(phi_tot(k,:,ibeg:iend),coef(0,:))
      ENDDO
!
      j=j+1
      j_ant=j
      xarr_sol(j)=xfine_tot(i)
      fun_sol(:,:,j)=matmul(phi_loc,Fbig_sol(:,:,i))
!
!write(1000,*) xfine_tot(i),dble(Fbig_sol(:,:,i)) &
!                          ,dimag(Fbig_sol(:,:,i))
    endif
  enddo
!
  npoi_out_adapt=j
!
! Processing output :
!
  if(iswout.eq.1.and.npoi_out_adapt.le.npoi_out) then
    x(1:npoi_out_adapt)=xarr_sol(1:npoi_out_adapt)
    fun(:,1:nrhs,1:npoi_out_adapt)=fun_sol(:,:,1:npoi_out_adapt)
  else
    h_out=(xend-xbeg)/(npoi_out-1)
!
    DO i=1,npoi_out
!
      x(i)=xbeg+h_out*(i-1)
      if(x(i).le.xarr_sol(j_ant)) then
!
        CALL binsrc(xarr_sol(1:j_ant),1,j_ant,x(i),ind)
!
        ibeg=MAX(1,MIN(j_ant-nlagr,ind-nshift-1))
        iend=MIN(j_ant,ibeg+nlagr)
!
      else
!
        CALL binsrc(xarr_sol(1:npoi_out_adapt),j_ant+1,npoi_out_adapt,x(i),ind)
!
        ibeg=MAX(j_ant+1,MIN(npoi_out_adapt-nlagr,ind-nshift-1))
        iend=MIN(npoi_out_adapt,ibeg+nlagr)
!
      endif
!
      CALL plag_coeff(iend-ibeg+1,nder,x(i),xarr_sol(ibeg:iend), &
                      coef(:,1:iend-ibeg+1))
!
      DO k=1,nrhs
        fun(:,k,i)=MATMUL(fun_sol(:,k,ibeg:iend),coef(0,1:iend-ibeg+1))
      ENDDO
!
    ENDDO
!
  endif
!
  DEALLOCATE(coef,phi_loc)
  IF(ALLOCATED(xfine)) DEALLOCATE(xfine,Fbig)
  call stiff_solver_deallocate
!
  END SUBROUTINE stiff_solver
