!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_matrix(iab)
!
! Computes the right-hand side matrix A and vector B
! of the linear ODE set f' = A f + B
! B assumes many right hand sides, i.e. B is a matrix itself
!
  use sample_matrix_mod, only : nsize,nrhs,x,amat,brhs
!
  implicit none
!
  integer :: iab
  double complex :: vec1,vec2
!
  nsize=4
  nrhs=0
!
  if(.not.allocated(amat)) allocate(amat(nsize,nsize))
  if(.not.allocated(brhs).and.nrhs.gt.0) allocate(brhs(nsize,nrhs))
!
  if(iab.ne.3) then
!
    call wave_matrix(x,amat)
    vec2=(0.d0,0.1d0)
    vec1=(10.d0,0.1d0)
!    amat=(0.d0,0.d0)
!    amat(1,3)=vec1
!    amat(2,4)=vec2
!    amat(3,1)=-vec1
!    amat(4,2)=-vec2
!
  endif
!
  if(iab.eq.1) return
!
  brhs=(0.d0,0.d0)
  vec2=(0.d0,0.1d0)
  brhs(4,1)=exp(-(x-120.d0)**2)*(vec2+(4.d0*(x-120.d0)**2-2.d0)/vec2)
!
  end subroutine get_matrix
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
