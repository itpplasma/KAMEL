!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine split_normal_modes(nlagr,nsize,npoi,niter,eps,t,amat,   &
                                nstiff,phi,phi_inv,alam,Bmat,ierr)
!
  implicit none
!
  integer,          parameter :: nder=1
!
  integer :: nlagr,nsize,npoi,nstiff,niter,ierr
  double precision :: eps,philowlim
  double precision, dimension(npoi)             :: t
  double complex,   dimension(nsize,npoi)       :: alam
  double complex,   dimension(nsize,nsize,npoi) :: amat,phi,phi_inv,Bmat
!
  integer :: nshift,npoilag,ibeg,iend
  integer :: info,i,j,k,iter,iflag
  integer,          dimension(:),     allocatable :: ipoi
  double precision, dimension(:,:),   allocatable :: coef
  double complex,   dimension(:),     allocatable :: alam_tmp
  double complex,   dimension(:,:),   allocatable :: Abig,chi_tmp,alam_old
  double complex,   dimension(:,:),   allocatable :: chi,chi_inv,c
  double complex,   dimension(:,:),   allocatable :: der_phi
  double complex,   dimension(:,:,:), allocatable :: phi_old,c_inv,chi_old
!
  ierr=0
  philowlim=100.d0*epsilon(1.d0)/eps
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  allocate(coef(0:nder,npoilag))
  allocate(Abig(nsize,nsize))
  allocate(chi(nsize,nsize),chi_inv(nsize,nsize))
  allocate(der_phi(nsize,nsize))
  allocate(c(nsize,nsize),c_inv(nsize,nsize,npoi))
  allocate(phi_old(nsize,nsize,npoi),alam_old(nsize,npoi))
  allocate(chi_old(nsize,nsize,npoi),ipoi(nsize))
  allocate(alam_tmp(nsize),chi_tmp(nsize,nsize))
!
  chi_old=phi
!
  do iter=1,niter
!
    alam_old=alam
    phi_old=phi
    iflag=0
!
    do i=1,npoi
!
! derivative of varphi :
!
      ibeg=max(1,min(npoi-nlagr,i-nshift))
      iend=ibeg+nlagr
      call plag_coeff(npoilag,nder,t(i),t(ibeg:iend),coef)
      do j=1,nsize
        der_phi(:,j)=matmul(phi_old(:,j,ibeg:iend),coef(1,:))
      enddo
!
! inverse varphi :
!
      call invert_cmat(nsize,phi_old(:,:,i),phi_inv(:,:,i),info)
!
      if(info.ne.0) then
        ierr=2
        print *,'split_normal_modes error: inversion of varphi failed'
        return
      endif
!
! matrix A :
!
      Abig=amat(:,:,i)-matmul(der_phi,phi_inv(:,:,i))
!
      if(nstiff.eq.0) then
!
! compute matrix B immediately :
!
        Bmat(:,:,i)=matmul(matmul(phi_inv(:,:,i),Abig),phi_old(:,:,i))
        cycle
!
      endif
!
! chi and lambda :
!
      call eigen_cmat(nsize,Abig,alam(:,i),chi,info)
!
      if(info.ne.0) then
        ierr=3
        print *,'split_normal_modes error: eigenvectors failed'
        return
      endif
!
      call align_lambdas(nsize,nstiff,alam_old(:,i),alam(:,i),ipoi,info)
!
      if(info.gt.0) then
        ierr=2
        print *,'split_normal_modes error: cannot align eigenvalues'
        return
      elseif(info.lt.0) then
        alam_tmp(1:nstiff)=alam(ipoi(1:nstiff),i)
        alam(1:nstiff,i)=alam_tmp(1:nstiff)
        chi_tmp(:,1:nstiff)=chi(:,ipoi(1:nstiff))
        chi(:,1:nstiff)=chi_tmp(:,1:nstiff)
      endif
!
      do j=1,nstiff
        if(real(sum(conjg(chi(:,j))*chi_old(:,j,i))).lt.0.d0) then
          chi(:,j)=-chi(:,j)
        endif
      enddo
!
! update varphi :
!
      phi(:,1:nstiff,i)=chi(:,1:nstiff)
!
      c_inv(:,:,i)=matmul(phi_inv(:,:,i),chi)
!
      do k=1,nstiff
        if(abs(alam_old(k,i)-alam(k,i))       .gt.                            &
           (abs(alam_old(k,i))+abs(alam(k,i)))*eps) iflag=1
        do j=1,nsize
          if( abs(chi(j,k))+abs(chi_old(j,k,i)).gt.philowlim .and.            &
              abs(chi(j,k)-chi_old(j,k,i))     .gt.                           &
             (abs(chi(j,k))+abs(chi_old(j,k,i)))*eps) iflag=1
        enddo
      enddo
!
      chi_old(:,:,i)=chi
!
    enddo
!
    if(iflag.eq.0) then
      niter=iter
      exit
    endif
!
  enddo
!
  if(nstiff.eq.0) then
    deallocate(coef,Abig,chi,chi_inv,der_phi,phi_old,c,c_inv,alam_old)
    deallocate(chi_old,ipoi,alam_tmp,chi_tmp)
    niter=0
    return
  endif
!
  if(iflag.eq.1) then
    ierr=4
    print *,'split_normal_modes error: maximum number of iterations exceeded'
    return
  endif
!
  do i=1,npoi
    Bmat(1:nsize,1:nstiff,i)=(0.d0,0.d0)
    do j=1,nstiff
      Bmat(j,j,i)=alam(j,i)
    enddo
!
    call invert_cmat(nsize,c_inv(:,:,i),c,info)
!
    if(info.ne.0) then
      ierr=5
      print *,'split_normal_modes error: inversion of c_inv failed'
      return
    endif
!
    do j=1,nsize
      do k=nstiff+1,nsize
        Bmat(j,k,i)=sum(c_inv(j,nstiff+1:nsize,i)*alam(nstiff+1:nsize,i) &
                        *c(nstiff+1:nsize,k))
        if(j.le.nstiff)  Bmat(j,k,i)= Bmat(j,k,i)+alam(j,i)*c(j,k)
      enddo
    enddo
!
  enddo
!
  deallocate(coef,Abig,chi,chi_inv,der_phi,phi_old,c,c_inv,alam_old)
  deallocate(chi_old,ipoi,alam_tmp,chi_tmp)
!
  end subroutine split_normal_modes
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine stiffness_estimator(nsize,npoi,nstiff_min,tol,t,alam, &
                                 nstiff,alam_mild)
!
  implicit none
!
  integer :: nsize,npoi,nstiff_min,i,j
  double precision :: tol,Deltat
!
  integer,          dimension(npoi)       :: nstiff
  double precision, dimension(npoi)       :: t,alam_mild
  double complex,   dimension(nsize,npoi) :: alam
!
  integer,          dimension(:),   allocatable :: ipoi
  double complex,   dimension(:),   allocatable :: alam_ord
!
  allocate(ipoi(nsize),alam_ord(nsize))
!
  do i=2,npoi-1
    nstiff(i)=nsize
    Deltat=min(t(i+1)-t(i),t(i)-t(i-1))
    call order_lambdas(nsize,alam(:,i),ipoi)
    alam_ord=alam(ipoi(:),i)
    do j=1,nsize-1
      if(minval(abs(alam_ord(j)-alam_ord(j+1:nsize)))*Deltat.lt.tol) then
        nstiff(i)=min(nstiff(i),j-1)
        exit
      endif
    enddo
    nstiff(i)=min(nstiff(i),nstiff_min)
    alam_mild(i)=maxval(abs(alam_ord(nstiff(i)+1:nsize)))
  enddo
!
  nstiff(1)=nstiff(2)
  call order_lambdas(nsize,alam(:,1),ipoi)
  alam_ord=alam(ipoi(:),1)
  alam_mild(1)=maxval(abs(alam_ord(nstiff(1)+1:nsize)))
  nstiff(npoi)=nstiff(npoi-1)
  call order_lambdas(nsize,alam(:,npoi),ipoi)
  alam_ord=alam(ipoi(:),npoi)
  alam_mild(npoi)=maxval(abs(alam_ord(nstiff(npoi)+1:nsize)))
!
  deallocate(ipoi,alam_ord)
!
  end subroutine stiffness_estimator
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine detstiffreg(npoi,nstiff_loc,nstencil,nreg,ibounds,nstiff_reg,ierr)
!
  implicit none
!
  integer :: npoi,nstencil,nreg,nstiff,ierr,nreg_new
  integer :: i,j,k,iflag,nstiff_min,ibeg,iend,nstiff_m,nstiff_p,nstiff_absmin
  integer, dimension(npoi) :: nstiff_loc,ibounds,nstiff_reg
  integer, dimension(:), allocatable :: nstiff_old,nstiff_new,nvtul
!
  ierr=0
!
  if(nstencil.lt.2) then
    print *,'detstifreg error: stencil size < 2'
    ierr=1
    return
  endif
!
  allocate(nstiff_old(npoi),nstiff_new(npoi),nvtul(npoi))
!
  nstiff_old=nstiff_loc
  nstiff_new=nstiff_old
!
! Enlarge all "less stiff modes" intervals by one at the ends:
  do i=2,npoi-1
    nstiff_new(i)=min(nstiff_old(i-1),nstiff_old(i),nstiff_old(i+1))
  enddo
  nstiff_old=nstiff_new
!
! Enlarge "less stiff modes" intervals to the stencil size:
  do
    iflag=0
    nstiff=nstiff_old(1)
    nvtul(1)=1
    k=1
    nstiff_min=maxval(nstiff_old)
    do i=2,npoi
      if(nstiff_old(i).ne.nstiff) then
        if(k.lt.nstencil) nstiff_min=min(nstiff_min,nstiff)
        nstiff=nstiff_old(i)
        k=1
      else
        k=k+1
      endif
      nvtul(i)=k
    enddo
    do i=1,npoi-1
      if(nstiff_old(i).eq.nstiff_min.and.nvtul(i).lt.nstencil) then
        nstiff_new(i+1)=min(nstiff_old(i),nstiff_old(i+1))
        if(nstiff_new(i+1).ne.nstiff_old(i+1)) iflag=1
      endif
    enddo
    nstiff_old=nstiff_new
    nstiff_absmin=nstiff_min
!
    nstiff=nstiff_old(npoi)
    nvtul(npoi)=1
    k=1
    nstiff_min=maxval(nstiff_old)
    do i=npoi-1,1,-1
      if(nstiff_old(i).ne.nstiff) then
        if(k.lt.nstencil) nstiff_min=min(nstiff_min,nstiff)
        nstiff=nstiff_old(i)
        k=1
      else
        k=k+1
      endif
      nvtul(i)=k
    enddo
    do i=npoi,2,-1
      if(nstiff_old(i).eq.nstiff_min.and.nvtul(i).lt.nstencil) then
        nstiff_new(i-1)=min(nstiff_old(i),nstiff_old(i-1))
        if(nstiff_new(i-1).ne.nstiff_old(i-1)) iflag=1
      endif
    enddo
    nstiff_old=nstiff_new
    nstiff_min=min(nstiff_absmin,nstiff_min)
!
    if(iflag.eq.0) then
      ibeg=0
      iend=0
      do i=npoi,1,-1
        if(nstiff_old(i).eq.nstiff_min) then
          if(iend.eq.0) then
            iend=i
            ibeg=i
          else
            ibeg=ibeg-1
          endif
        elseif(iend.ne.0) then
          if(iend-ibeg+1.lt.nstencil) then
            if(iend.eq.npoi) then
              nstiff_new(ibeg:iend)=nstiff_old(ibeg-1)
            else
              nstiff_m=nstiff_old(ibeg-1)
              nstiff_p=nstiff_old(iend+1)
              nstiff_new(ibeg:iend)=max(nstiff_m,nstiff_p)
            endif
            iflag=1
          endif
          ibeg=0
          iend=0
        endif
      enddo
      if(ibeg.eq.1.and.iend.lt.nstencil) then
        nstiff_new(ibeg:iend)=nstiff_old(iend+1)
        iflag=1
      endif
      nstiff_old=nstiff_new
    endif
    if(iflag.eq.0) exit
  enddo
!
! Determine the number of regions and boundaries:
  nreg=1
  ibounds(1)=1
  nstiff=nstiff_old(1)
  nstiff_reg(1)=nstiff
  do i=2,npoi
    if(nstiff_old(i).ne.nstiff) then
      nreg=nreg+1
      ibounds(nreg)=i
      nstiff=nstiff_old(i)
      nstiff_reg(nreg)=nstiff
    endif
  enddo
  ibounds(nreg+1)=npoi
!
  deallocate(nstiff_old,nstiff_new,nvtul)
!
  end subroutine detstiffreg
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
