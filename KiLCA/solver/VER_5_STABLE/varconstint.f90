!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine varconstint
  end subroutine varconstint
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  implicit none
!
  integer, parameter :: npol=10
  integer :: i,nsize,ierr
!
  double complex :: psi_b,psi_e,imun,alam,alam1
  double complex, dimension(0:npol) :: psi,weight
!
  imun=(0.d0,1.d0)
!
  do i=0,npol
    psi(i)=(3.14d0+imun)*i/dble(npol)
  enddo
!
  psi_b=psi(3)
  psi_e=psi(7)
!
  call step_varconstint(npol,psi,psi_b,psi_e,weight,ierr)
!
  alam=imun-(1.d0,0.d0)
  alam1=-imun-(1.d0,0.d0)
!
  print *,sum(weight*sin(psi)), &
  0.5d0*((exp(alam*psi_e)-exp(alam*psi_b))*exp(psi_b)/alam &
  -(exp(alam1*psi_e)-exp(alam1*psi_b))*exp(psi_b)/alam1)/imun
!
  end
