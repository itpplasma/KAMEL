!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs_stiff(x,y,dy)
!
  use sample_matrix_mod, only : nlagr,npoi,nstiff,xarr,bmat,i_int
!
  implicit none
!
  integer, parameter :: nder=0
  integer :: ibeg,iend,nshift,npoilag,i,k
  double precision  :: x
  double complex    :: alam
  double precision, dimension(2*nstiff) :: y,dy
  double precision, dimension(:,:),   allocatable :: coef
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  allocate(coef(0:nder,npoilag))
!
  ibeg=max(1,min(npoi-nlagr,i_int-nshift))
  iend=ibeg+nlagr
!
  call plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
!
  k=0
!
  do i=1,nstiff
    alam=sum(bmat(i,i,ibeg:iend)*coef(0,:))
    k=k+1
    dy(k)=dble(alam)
    k=k+1
    dy(k)=dimag(alam)
  enddo
!
  deallocate(coef)
!
  return
  end subroutine rhs_stiff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine integrate_eikonals(eps)
!
  use sample_matrix_mod, only : npoi,nstiff,xarr,i_int,eikonals
!
  implicit none
!
  integer :: nvar,i
  double precision :: eps
  double precision, dimension(:), allocatable :: y
!
  external :: rhs_stiff
!
  nvar=2*nstiff
!
  eikonals(:,1)=(0.d0,0.d0)
!
  allocate(y(nvar))
  y=0.d0
!
  do i_int=2,npoi
!
    call odeint_allroutines(y,nvar,xarr(i_int-1),xarr(i_int),eps,rhs_stiff)
!
    do i=1,nstiff
      eikonals(i,i_int)=dcmplx(y(2*i-1),y(2*i))
    enddo
!
  enddo
!
  deallocate(y)
!
  do i=1,nstiff
    eikonals(i,:)=eikonals(i,:)-maxval(real(eikonals(i,:)))
  enddo
!
  return
  end subroutine integrate_eikonals
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs_mild(x,y,dy)
!
  use sample_matrix_mod, only : nlagr,npoi,nsize,nstiff,xarr,bmat
!
  implicit none
!
  integer, parameter :: nder=0
  integer :: ibeg,iend,nshift,npoilag,i,j,k,nblock
  double precision  :: x
  double precision, dimension(2*(nsize-nstiff)**2) :: y,dy
  double precision, dimension(:,:), allocatable :: coef
  double complex,   dimension(:,:), allocatable :: z,dz,bmat_block
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  nblock=nsize-nstiff
!
  allocate(coef(0:nder,npoilag),bmat_block(nblock,nblock))
  allocate(z(nblock,nblock),dz(nblock,nblock))
!
  call binsrc(xarr,1,npoi,x,i)
!
  ibeg=max(1,min(npoi-nlagr,i-nshift-1))
  iend=ibeg+nlagr
!
  call plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
!
  do i=1,nblock
    bmat_block(i,:)=matmul(bmat(nstiff+i,nstiff+1:nsize,ibeg:iend),coef(0,:))
  enddo
!
  k=0
!
  do i=1,nblock
    do j=1,nblock
      k=k+1
      z(i,j)=dcmplx(y(k),y(k+1))
      k=k+1
    enddo
  enddo
!
  dz=matmul(bmat_block,z)
!
  k=0
!
  do i=1,nblock
    do j=1,nblock
      k=k+1
      dy(k)=dble(dz(i,j))
      k=k+1
      dy(k)=dimag(dz(i,j))
    enddo
  enddo
!
  deallocate(coef,bmat_block,z,dz)
!
  return
  end subroutine rhs_mild
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhs_mild_inhom(x,y,dy)
!
  use sample_matrix_mod, only : nlagr,npoi,nsize,nstiff,xarr,bmat, &
                                nrhs,npoi_rhs,xarr_rhs,wmat
!
  implicit none
!
  integer, parameter :: nder=0
  integer :: ibeg,iend,nshift,npoilag,i,j,k,nblock
  double precision  :: x
  double precision, dimension(2*(nsize-nstiff)*nrhs) :: y,dy
  double precision, dimension(:,:), allocatable :: coef
  double complex,   dimension(:,:), allocatable :: z,dz,bmat_block,wmat_block
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  nblock=nsize-nstiff
!
  allocate(coef(0:nder,npoilag),bmat_block(nblock,nblock))
  allocate(z(nblock,nrhs),dz(nblock,nrhs),wmat_block(nblock,nrhs))
!
  call binsrc(xarr,1,npoi,x,i)
!
  ibeg=max(1,min(npoi-nlagr,i-nshift-1))
  iend=ibeg+nlagr
!
  call plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
!
  do i=1,nblock
    bmat_block(i,:)=matmul(bmat(nstiff+i,nstiff+1:nsize,ibeg:iend),coef(0,:))
  enddo
!
  call binsrc(xarr_rhs,1,npoi_rhs,x,i)
!
  ibeg=max(1,min(npoi_rhs-nlagr,i-nshift-1))
  iend=ibeg+nlagr
!
  call plag_coeff(npoilag,nder,x,xarr_rhs(ibeg:iend),coef)
!
  do i=1,nblock
    wmat_block(i,:)=matmul(wmat(nstiff+i,:,ibeg:iend),coef(0,:))
  enddo
!
  k=0
!
  do i=1,nblock
    do j=1,nrhs
      k=k+1
      z(i,j)=dcmplx(y(k),y(k+1))
      k=k+1
    enddo
  enddo
!
  dz=matmul(bmat_block,z)+wmat_block
!
  k=0
!
  do i=1,nblock
    do j=1,nrhs
      k=k+1
      dy(k)=dble(dz(i,j))
      k=k+1
      dy(k)=dimag(dz(i,j))
    enddo
  enddo
!
  deallocate(coef,bmat_block,z,dz,wmat_block)
!
  return
  end subroutine rhs_mild_inhom
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine step_mild(eps,x1,x2,z)
!
  use sample_matrix_mod, only : nsize,nstiff
!
  implicit none
!
  integer :: nblock,i,j,k,nvar
  double precision :: eps,x1,x2
  double complex,   dimension(nsize-nstiff,nsize-nstiff) :: z
  double precision, dimension(:), allocatable            :: y
!
  external :: rhs_mild
!
  nblock=nsize-nstiff
  nvar=2*nblock*nblock
!
  allocate(y(nvar))
!
  k=0
!
  do i=1,nblock
    do j=1,nblock
      k=k+1
      y(k)=dble(z(i,j))
      k=k+1
      y(k)=dimag(z(i,j))
    enddo
  enddo
!
  call odeint_allroutines(y,nvar,x1,x2,eps,rhs_mild)
!
  k=0
!
  do i=1,nblock
    do j=1,nblock
      k=k+1
      z(i,j)=dcmplx(y(k),y(k+1))
      k=k+1
    enddo
  enddo
!
  deallocate(y)
!
  end subroutine step_mild
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine step_mild_inhom(eps,x1,x2,z)
!
  use sample_matrix_mod, only : nsize,nstiff,nrhs
!
  implicit none
!
  integer :: nblock,i,j,k,nvar
  double precision :: eps,x1,x2
  double complex,   dimension(nsize-nstiff,nrhs)   :: z
  double precision, dimension(:), allocatable      :: y
!
  external :: rhs_mild_inhom
!
  nblock=nsize-nstiff
  nvar=2*nblock*nrhs
!
  allocate(y(nvar))
!
  k=0
!
  do i=1,nblock
    do j=1,nrhs
      k=k+1
      y(k)=dble(z(i,j))
      k=k+1
      y(k)=dimag(z(i,j))
    enddo
  enddo
!
  call odeint_allroutines(y,nvar,x1,x2,eps,rhs_mild_inhom)
!
  k=0
!
  do i=1,nblock
    do j=1,nrhs
      k=k+1
      z(i,j)=dcmplx(y(k),y(k+1))
      k=k+1
    enddo
  enddo
!
  deallocate(y)
!
  end subroutine step_mild_inhom
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine step_varconstint(npol,psi,psi_b,psi_e,weight,ierr)
!
  implicit none
!
  integer :: npol,i,nsize,ierr,info
!
  double complex :: psi_b,psi_e,del,term
  double complex, dimension(0:npol) :: psi,weight
!
  double complex, dimension(:),   allocatable :: delpsi,sumser
  double complex, dimension(:,:), allocatable :: tailcoef,dercoef
!
  ierr=0
!
  allocate(delpsi(0:npol),sumser(0:npol))
  allocate(tailcoef(0:npol,0:npol),dercoef(0:npol,0:npol))
!
  delpsi=psi-psi_b
  del=psi_e-psi_b
!
  tailcoef(:,0)=1.d0
!
  do i=1,npol
     tailcoef(:,i)=tailcoef(:,i-1)*delpsi/dble(i)
  enddo
!
  nsize=npol+1
!
  call invert_cmat(nsize,tailcoef(0:npol,0:npol),dercoef(0:npol,0:npol),info)
!
  if(info.ne.0) then
    ierr=1
    print *,'varconsts : error in finding derivatives'
    return
  endif
!
  term=exp(-del)
  sumser(0)=term
!
  do i=1,npol
    term=term*del/i
    sumser(i)=sumser(i-1)+term
  enddo
!
  weight=matmul((1.d0,0.d0)-sumser,dercoef)
!
  deallocate(delpsi,sumser,tailcoef,dercoef)
!
  end subroutine step_varconstint
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine residual_system(nlagr,nsize,npoi,nstiff,nstep,eps,xarr,Bmat,      &
                             eikonals,xfine,istart,istop,iflag,Fbig,ierr)
!
  implicit none
!
  integer, parameter :: nder=0
!
  integer :: nlagr,nsize,nstiff,npoi,nstep,nblock,npoilag,nshift,ind,ibeg,iend
  integer :: j,k,ierr,istart,istop,istart_margin,istop_margin,imn,imx,iflag
  double precision :: eps
  double precision,  dimension(npoi)             :: xarr
  double complex,    dimension(nstiff,npoi)      :: eikonals
  double complex,    dimension(nsize,nsize,npoi) :: Bmat
  double precision,  dimension(0:nstep)             :: xfine
  double complex,    dimension(nsize,nsize,0:nstep) :: Fbig
!
  integer :: i
  double precision :: x1,x2
  double precision, dimension(:,:),   allocatable :: coef
  double complex,   dimension(:),     allocatable :: weight,expdpsi
  double complex,   dimension(:,:),   allocatable :: z,psi_arr
  double complex,   dimension(:,:),   allocatable :: del_psi,psi_int
  double complex,   dimension(:,:,:), allocatable :: subint
!
  ierr=0
!
  Fbig=(0.d0,0.d0)
!
  npoilag=nlagr+1
  nshift=nlagr/2
  allocate(coef(0:nder,npoilag))
!
  if(nstiff.eq.nsize) then
!
    do i=istart,istop
!
      call binsrc(xarr,1,npoi,xfine(i),ind)
!
      ibeg=max(1,min(npoi-nlagr,ind-nshift-1))
      iend=ibeg+nlagr
!
      call plag_coeff(npoilag,nder,xfine(i),xarr(ibeg:iend),coef)
!
      do k=1,nstiff
        Fbig(k,k,i)=sum(eikonals(k,ibeg:iend)*coef(0,:))
      enddo
!
    enddo
!
    deallocate(coef)
    return
!
  endif
!
  nblock=nsize-nstiff
!
  allocate(z(nblock,nblock))
!
  z=(0.d0,0.d0)
!
  do i=1,nblock
    z(i,i)=(1.d0,0.d0)
  enddo
!
  Fbig(nstiff+1:nsize,nstiff+1:nsize,istart)=z
!
  x1=xfine(istart)
!
  imn=0
  imx=nstep
  if(iflag.eq.-1) then
    imn=1
  elseif(iflag.eq.1) then
    imx=nstep-1
  endif
  istart_margin=max(imn,istart-nshift)
  istop_margin=min(imx,istop+nlagr-nshift)
  if(istop_margin-istart_margin.lt.nlagr) then
    if(istart_margin.eq.imn) then
      istop_margin=nlagr
    elseif(istop_margin.eq.imx) then
      istart_margin=imx-nlagr
    else
      print *,'residual_system error : npoi < nlagr'
      ierr=2
      return
    endif
  endif
!
  do i=istart-1,istart_margin,-1
    x2=xfine(i)
    call step_mild(eps,x1,x2,z)
    x1=x2
    Fbig(nstiff+1:nsize,nstiff+1:nsize,i)=z
  enddo
!
  x1=xfine(istart)
  z=Fbig(nstiff+1:nsize,nstiff+1:nsize,istart)
!
  do i=istart+1,istop_margin
    x2=xfine(i)
    call step_mild(eps,x1,x2,z)
    x1=x2
    Fbig(nstiff+1:nsize,nstiff+1:nsize,i)=z
  enddo
!
  if(nstiff.eq.0) then
    deallocate(z)
    return
  endif
!
  allocate(psi_arr(nstiff,0:nstep),weight(npoilag),subint(nstiff,nsize,0:nstep))
  allocate(del_psi(0:nstep,0:nstep),psi_int(0:nstep,nsize),expdpsi(0:nstep))
!
  subint=(0.d0,0.d0)
!
  do i=istart_margin,istop_margin
!
    call binsrc(xarr,1,npoi,xfine(i),ind)
!
    ibeg=max(1,min(npoi-nlagr,ind-nshift-1))
    iend=ibeg+nlagr
!
    call plag_coeff(npoilag,nder,xfine(i),xarr(ibeg:iend),coef)
!
    do j=1,nstiff
      subint(j,nstiff+1:nsize,i)=matmul(matmul(coef(0,:),   &
             transpose(Bmat(j,nstiff+1:nsize,ibeg:iend))),  &
             Fbig(nstiff+1:nsize,nstiff+1:nsize,i))         &
            /sum(Bmat(j,j,ibeg:iend)*coef(0,:))
    enddo
!
    psi_arr(:,i)=matmul(eikonals(:,ibeg:iend),coef(0,:))
!
  enddo
!
  do k=1,nstiff
    Fbig(k,k,istart:istop)=psi_arr(k,istart:istop)
    if(dble(psi_arr(k,istop)).gt.dble(psi_arr(k,istart))) then
      do i=istart,istop-1
        ibeg=max(istart_margin,min(istop_margin-nlagr,i-nshift))
        iend=ibeg+nlagr
        call step_varconstint(nlagr,psi_arr(k,ibeg:iend),                    &
                              psi_arr(k,i),psi_arr(k,i+1),weight,ierr)
        if(ierr.ne.0) return
        psi_int(i,nstiff+1:nsize)=matmul(subint(k,nstiff+1:nsize,ibeg:iend), &
                                         weight)
        expdpsi(i)=exp(psi_arr(k,i)-psi_arr(k,i+1))
      enddo
      Fbig(k,nstiff+1:nsize,istop)=(0.d0,0.d0)
      do i=istop-1,istart,-1
        Fbig(k,nstiff+1:nsize,i)=Fbig(k,nstiff+1:nsize,i+1)*expdpsi(i)       &
                                -psi_int(i,nstiff+1:nsize)
      enddo
    else
      do i=istart+1,istop
        ibeg=max(istart_margin,min(istop_margin-nlagr,i-nshift))
        iend=ibeg+nlagr
        call step_varconstint(nlagr,psi_arr(k,ibeg:iend),                    &
                              psi_arr(k,i),psi_arr(k,i-1),weight,ierr)
        if(ierr.ne.0) return
        psi_int(i,nstiff+1:nsize)=matmul(subint(k,nstiff+1:nsize,ibeg:iend), &
                                         weight)
        expdpsi(i)=exp(psi_arr(k,i)-psi_arr(k,i-1))
      enddo
      Fbig(k,nstiff+1:nsize,istart)=(0.d0,0.d0)
      do i=istart+1,istop
        Fbig(k,nstiff+1:nsize,i)=Fbig(k,nstiff+1:nsize,i-1)*expdpsi(i)       &
                                -psi_int(i,nstiff+1:nsize)
      enddo
    endif
  enddo
!
  deallocate(z,subint,coef,psi_arr,weight,del_psi,psi_int,expdpsi)
!
  end subroutine residual_system
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine residual_system_inhom(nlagr,nsize,nrhs,npoi,npoi_rhs,nstiff,      &
                                   nstep,eps,xarr,xarr_rhs,Bmat,wmat,eikonals, &
                                   xfine,istart,istop,Fbig_rhs,ierr)
!
  implicit none
!
  integer, parameter :: nder=0
!
  integer :: nlagr,nsize,nstiff,npoi,nstep,nblock,npoilag,nshift,ind,ibeg,iend
  integer :: j,k,ierr,nrhs,npoi_rhs,istart,istop,istart_margin,istop_margin
  double precision :: eps
  double precision,  dimension(npoi)                :: xarr
  double precision,  dimension(npoi_rhs)            :: xarr_rhs
  double complex,    dimension(nstiff,npoi)         :: eikonals
  double complex,    dimension(nsize,nsize,npoi)    :: Bmat
  double complex,    dimension(nsize,nrhs,npoi_rhs) :: wmat
  double precision,  dimension(0:nstep)             :: xfine
  double complex,    dimension(nsize,nrhs,0:nstep)  :: Fbig_rhs
!
  integer :: i
  double precision :: x1,x2
  double precision, dimension(:,:),   allocatable :: coef
  double complex,   dimension(:),     allocatable :: weight,expdpsi,alam_loc
  double complex,   dimension(:,:),   allocatable :: z,psi_arr
  double complex,   dimension(:,:),   allocatable :: del_psi,psi_int
  double complex,   dimension(:,:,:), allocatable :: subint
!
  ierr=0
!
  Fbig_rhs=(0.d0,0.d0)
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  nblock=nsize-nstiff
!
  allocate(z(nblock,nrhs),alam_loc(nstiff))
!
  z=(0.d0,0.d0)
!
  Fbig_rhs(nstiff+1:nsize,:,0)=z
!
  x1=xfine(istart)
!
  istart_margin=max(0,istart-nshift-2)
  istop_margin=min(nstep,istop+nlagr-nshift)
  if(istop_margin-istart_margin.lt.nlagr) then
    if(istart_margin.eq.0) then
      istop_margin=nlagr
    elseif(istop_margin.eq.nstep) then
      istart_margin=nstep-nlagr
    else
      print *,'residual_system_inhom error : npoi < nlagr'
      ierr=2
      return
    endif
  endif
!
  do i=istart-1,istart_margin,-1
    x2=xfine(i)
    call step_mild_inhom(eps,x1,x2,z)
    x1=x2
    Fbig_rhs(nstiff+1:nsize,:,i)=z
  enddo
!
  x1=xfine(istart)
  z=Fbig_rhs(nstiff+1:nsize,:,istart)
!
  do i=istart+1,istop_margin
    x2=xfine(i)
    call step_mild_inhom(eps,x1,x2,z)
    x1=x2
    Fbig_rhs(nstiff+1:nsize,:,i)=z
  enddo
!
  if(nstiff.eq.0) then
    deallocate(z)
    return
  endif
!
  allocate(coef(0:nder,npoilag))
  allocate(psi_arr(nstiff,0:nstep),weight(npoilag),subint(nstiff,nrhs,0:nstep))
  allocate(del_psi(0:nstep,0:nstep),psi_int(0:nstep,nrhs),expdpsi(0:nstep))
!
  subint=(0.d0,0.d0)
!
  do i=istart_margin,istop_margin
!
    call binsrc(xarr,1,npoi,xfine(i),ind)
!
    ibeg=max(1,min(npoi-nlagr,ind-nshift-1))
    iend=ibeg+nlagr
!
    call plag_coeff(npoilag,nder,xfine(i),xarr(ibeg:iend),coef)
!
    psi_arr(:,i)=matmul(eikonals(:,ibeg:iend),coef(0,:))
!
    do j=1,nstiff
      alam_loc(j)=sum(Bmat(j,j,ibeg:iend)*coef(0,:))
    enddo
!
    call binsrc(xarr_rhs,1,npoi_rhs,xfine(i),ind)
!
    ibeg=max(1,min(npoi_rhs-nlagr,ind-nshift-1))
    iend=ibeg+nlagr
!
    call plag_coeff(npoilag,nder,xfine(i),xarr_rhs(ibeg:iend),coef)
!
    do j=1,nstiff
      subint(j,:,i)=matmul(wmat(j,:,ibeg:iend),coef(0,:))/alam_loc(j)
    enddo
!
  enddo
!
  do k=1,nstiff
    if(dble(psi_arr(k,istop)).gt.dble(psi_arr(k,istart))) then
      do i=istart,istop-1
        ibeg=max(istart_margin,min(istop_margin-nlagr,i-nshift))
        iend=ibeg+nlagr
        call step_varconstint(nlagr,psi_arr(k,ibeg:iend),                    &
                              psi_arr(k,i),psi_arr(k,i+1),weight,ierr)
        if(ierr.ne.0) return
        psi_int(i,:)=matmul(subint(k,:,ibeg:iend),weight)
        expdpsi(i)=exp(psi_arr(k,i)-psi_arr(k,i+1))
      enddo
      Fbig_rhs(k,:,istop)=(0.d0,0.d0)
      do i=istop-1,istart,-1
        Fbig_rhs(k,:,i)=Fbig_rhs(k,:,i+1)*expdpsi(i)-psi_int(i,:)
      enddo
    else
      do i=istart+1,istop
        ibeg=max(istart_margin,min(istop_margin-nlagr,i-nshift))
        iend=ibeg+nlagr
        call step_varconstint(nlagr,psi_arr(k,ibeg:iend),                    &
                              psi_arr(k,i),psi_arr(k,i-1),weight,ierr)
        if(ierr.ne.0) return
        psi_int(i,:)=matmul(subint(k,:,ibeg:iend),weight)
        expdpsi(i)=exp(psi_arr(k,i)-psi_arr(k,i-1))
      enddo
      Fbig_rhs(k,:,istart)=(0.d0,0.d0)
      do i=istart+1,istop
        Fbig_rhs(k,:,i)=Fbig_rhs(k,:,i-1)*expdpsi(i)-psi_int(i,:)
      enddo
    endif
  enddo
!
  deallocate(z,subint,coef,psi_arr,weight,del_psi,psi_int,expdpsi)
!
  end subroutine residual_system_inhom
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
