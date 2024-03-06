!
  program driver_separate_eigenmodes
!
  use dieltens_mod, only : vecy,vecz
  use sample_matrix_mod, only : x
!
  implicit none
!
  integer :: i,j,nsize,nlagr,nrhs,itermax,npoiout,info,ndim,nant,ierr,iunit
  integer :: n_inner,iblablabla,npoiout_adaptive,iswout
  double precision :: eps,xbeg,xend,x_A
  double precision, dimension(:),     allocatable :: xx
  double complex,   dimension(:),     allocatable :: bouvec,antvec
  double complex,   dimension(:,:,:), allocatable :: fun
  double complex,   dimension(:,:),   allocatable :: amat,amat_inv,boumat
  double complex,   dimension(:,:),   allocatable :: fun_finsol
!
  iunit=1
  open(iunit,file='wave_driver.inp')
  read(iunit,*)
  read(iunit,*)
  read(iunit,*)
  read(iunit,*) npoiout
  read(iunit,*) xbeg
  read(iunit,*) xend
  read(iunit,*) x_A
  close(iunit)
  npoiout=npoiout+1
!
  nant=nint((npoiout-1)*(x_A-xbeg)/(xend-xbeg))
!
  nlagr=15
  itermax=49 !20
!  eps=1.d-10
  eps=1.d-7
!
  nsize=4
  n_inner=2
  nrhs=0
!
  allocate(xx(npoiout),fun(nsize,nsize+nrhs,npoiout))
  allocate(boumat(nsize,nsize),bouvec(nsize),antvec(nsize))
!
  x=1.d0
  call get_matrix(1)
!
  boumat=(0.d0,0.d0)
  boumat(1,1)=(1.d0,0.d0)
  boumat(2,2)=(1.d0,0.d0)
  boumat(3,1)=(1.d0,0.d0)
  boumat(4,2)=(1.d0,0.d0)
  bouvec=(0.d0,0.d0)
  antvec=(0.d0,0.d0)
  antvec(3)=vecy/vecz
  antvec(4)=(1.d0,0.d0)
!
  iblablabla=1
  iswout=1
!
  call stiff_solver(nsize,nlagr,nrhs,itermax,npoiout,eps,xbeg,xend,x_a, &
                    boumat,bouvec,antvec,n_inner,iswout,iblablabla,     &
                    npoiout_adaptive,xx,fun,ierr)
!
  if(ierr.ne.0) stop
!
  if(iswout.eq.0) npoiout_adaptive=npoiout
!
  do i=1,npoiout_adaptive
    write(100,*) xx(i),real(fun(:,1,i)),dimag(fun(:,1,i))
  enddo
!
  end program driver_separate_eigenmodes
