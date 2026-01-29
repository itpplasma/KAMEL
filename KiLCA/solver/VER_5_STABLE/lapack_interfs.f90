!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine invert_cmat(n,a,a_inv,info)
!
! Inverts double complex matrix a of the size n. Puts the result in a_inv.
!
  implicit none
!
  integer :: n,lwork,info
!
  double complex, dimension(n,n) :: a,a_inv
!
  integer,        dimension(:), allocatable :: ipiv
  double complex, dimension(:), allocatable :: work
!
  a_inv=a
!
  allocate(ipiv(n))
!
  call zgetrf(n,n,a_inv,n,ipiv,info)
!
  if(info.ne.0) then
    info=1
    return
  endif
!
  allocate(work(1))
  lwork=-1
!
  call zgetri(n,a_inv,n,ipiv,work,lwork,info)
!
  if(info.ne.0) then
    info=2
    return
  endif
!
  lwork=int(work(1))
  deallocate(work)
  allocate(work(lwork))
!
  call zgetri(n,a_inv,n,ipiv,work,lwork,info)
!
  if(info.ne.0) then
    info=3
    return
  endif
!
  deallocate(work,ipiv)
!
  end subroutine invert_cmat
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine eigen_cmat(n,a_in,alam,vecs,info)
!
! Find eigenvalues alam and eigenvectors vecs of the double complex matrix a_in
! of the size n. Eigenvalues are numbered in the descending order by their
! modules. If modules are equal, descending order is by comlex argument.
! Corresponding eigenvectors are numbered by the second index. The are
! normalized by module to one. Phase factor is fixed so that the square of
! eigenvector is positive (fixed up to the factor -1). If square of eigenvector
! is a numerical 0, then first in occurence component whose module is larger
! than half of largest component module, is fixed to be real and positive
!
  implicit none
!
  character(len=1), parameter :: jobvl='n',jobvr='v'
!
  integer :: n,lwork,info,i,j,isave
!
  double complex   :: theta
!
  double complex, dimension(n)   :: alam
  double complex, dimension(n,n) :: a_in,vecs
!
  integer,          dimension(:),   allocatable :: ipoi
  double precision, dimension(:),   allocatable :: rwork
  double complex,   dimension(:),   allocatable :: work
  double complex,   dimension(:,:), allocatable :: a,vl
!
  allocate(a(n,n),vl(n,n),rwork(2*n),work(1),ipoi(n))
!
  a=a_in
  lwork=-1
!
  call zgeev(jobvl,jobvr,n,a,n,alam,vl,n,vecs,n,work,lwork,rwork,info)
!
  if(info.ne.0) then
    info=1
    return
  endif
!
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
!
  call zgeev(jobvl,jobvr,n,a,n,alam,vl,n,vecs,n,work,lwork,rwork,info)
!
  if(info.ne.0) then
    info=2
    return
  endif
!
  call order_lambdas(n,alam,ipoi)
!
  deallocate(work)
  allocate(work(n))
!
  work=alam
  vl=vecs
!
  do i=1,n
    alam(i)=work(ipoi(i))
    theta=sum(vl(:,ipoi(i))**2)
    if(abs(theta).gt.1000.d0*epsilon(1.d0)) then
      theta=sqrt(theta/abs(theta))
    else
      theta=(0.d0,0.d0)
      do j=1,n
        if(abs(vl(j,ipoi(i))).gt.2.d0*abs(theta)) theta=vl(j,ipoi(i))
      enddo
      theta=theta/abs(theta)
    endif
    vecs(:,i)=vl(:,ipoi(i))/theta
  enddo
!
  deallocate(a,vl,rwork,work,ipoi)
!
  end subroutine eigen_cmat
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine qr_factorize(m,n,amat,umat,rmat,info)
!
  implicit none
!
  integer :: m,n,lda,mda,lwork,info,i,k
  double complex, dimension(m,n) :: amat
  double complex, dimension(m,m) :: umat
  double complex, dimension(m,n) :: rmat
  double complex, dimension(:),   allocatable :: tau,work
  double complex, dimension(:,:), allocatable :: atmp,v,one
!
  lda=max(m,n)
  mda=min(m,n)
!
  allocate(atmp(m,n),tau(mda),work(n),v(m,mda),one(m,m))
!
  atmp=amat
  one=(0.d0,0.d0)
  do i=1,m
    one(i,i)=(1.d0,0.d0)
  enddo
  lwork=-1
!
  call zgeqrf(m,n,atmp,lda,tau,work,lwork,info)
!
  if(info.ne.0) then
    info=1
    return
  endif
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
!
  call zgeqrf(m,n,atmp,lda,tau,work,lwork,info)
!
  if(info.ne.0) then
    info=2
    return
  endif
  rmat=(0.d0,0.d0)
  do i=1,m
    rmat(i,i:n)=atmp(i,i:n)
  enddo
  do k=1,mda
    v(1:k-1,k)=(0.d0,0.d0)
    v(k,k)=(1.d0,0.d0)
    v(k+1:m,k)=atmp(k+1:m,k)
  enddo
  umat=one
  do k=1,mda
    do i=1,m
      umat(i,:)=umat(i,:)-tau(k)*sum(umat(i,:)*v(:,k))*conjg(v(:,k))
    enddo
  enddo
!
  deallocate(atmp,tau,work,v,one)
!
  end subroutine qr_factorize
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine solve_cutriasys(n,nrhs,a,b,info)
!
  implicit none
!
  character(len=1), parameter :: uplo='u', trans='n', diag='n'
!
  integer :: n,nrhs,info
  double complex, dimension(n,n)    :: a
  double complex, dimension(n,nrhs) :: b
!
  call ztrtrs(uplo,trans,diag,n,nrhs,a,n,b,n,info)
!
  end subroutine solve_cutriasys
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine order_lambdas(n,alam,ipoi)
!
  implicit none
!
  integer :: n,i,j,isave
  double precision :: aboveone,belowone,abs1,abs2,arg1,arg2
  integer,        dimension(n)   :: ipoi
  double complex, dimension(n)   :: alam
!
  aboveone=1.d0+100.d0*epsilon(1.d0)
  belowone=1.d0-100.d0*epsilon(1.d0)
!
  do i=1,n
    ipoi(i)=i
  enddo
!
  do i=1,n-1
    do j=i+1,n
      abs1=abs(alam(ipoi(j)))
      abs2=abs(alam(ipoi(i)))
      if(abs1.gt.abs2*aboveone) then
        isave=ipoi(i)
        ipoi(i)=ipoi(j)
        ipoi(j)=isave
      elseif(abs1.gt.abs2*belowone) then
        arg1=atan2(dimag((alam(ipoi(j)))),dble((alam(ipoi(j)))))
        arg2=atan2(dimag((alam(ipoi(i)))),dble((alam(ipoi(i)))))
        if(arg1.gt.arg2) then
          isave=ipoi(i)
          ipoi(i)=ipoi(j)
          ipoi(j)=isave
        endif
      endif
    enddo
  enddo
!
  end subroutine order_lambdas
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine null_space(m,n,amat,zerospace,ierr)
!
  implicit none
!
  integer :: m,n,ierr,info
  double complex, dimension(m,n) :: amat
  double complex, dimension(n,n) :: zerospace
  double complex, dimension(:,:), allocatable :: umat,rmat
!
  if(n.le.m) then
    ierr=1
    return
  else
    ierr=0
  endif
!
  allocate(umat(n,n),rmat(n,m))
  call qr_factorize(n,m,transpose(amat),umat,rmat,info)
  if(info.ne.0) then
    ierr=2
    return
  endif
  zerospace(:,1:n-m)=conjg(umat(:,m+1:n))
  zerospace(:,n-m+1:n)=(0.d0,0.d0)
  deallocate(umat,rmat)
!
  end subroutine null_space
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
