  module sample_matrix_mod
    double precision, parameter :: antshift=0.3d0
    integer :: nlagr,nsize,nrhs,npoi,itermax,nstiff,i_int,npoi_rhs
    integer :: iflag_ant=0
    double precision :: x,xbeg,xend,eps,x_ant
    double precision, dimension(:),     allocatable :: xarr,xarr_rhs
    double complex,   dimension(:,:),   allocatable :: amat,brhs
    double complex,   dimension(:,:),   allocatable :: eikonals,alam
    double complex,   dimension(:,:,:), allocatable :: amat_arr,bmat,phi
    double complex,   dimension(:,:,:), allocatable :: phi_inv,wmat
  end module sample_matrix_mod
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine stiff_solver_deallocate
!
  use sample_matrix_mod
!
  if(allocated(xarr)) deallocate(xarr)
  if(allocated(xarr_rhs)) deallocate(xarr_rhs)
  if(allocated(amat)) deallocate(amat)
  if(allocated(brhs)) deallocate(brhs)
  if(allocated(eikonals)) deallocate(eikonals)
  if(allocated(alam)) deallocate(alam)
  if(allocated(amat_arr)) deallocate(amat_arr)
  if(allocated(bmat)) deallocate(bmat)
  if(allocated(phi)) deallocate(phi)
  if(allocated(phi_inv)) deallocate(phi_inv)
  if(allocated(wmat)) deallocate(wmat)
!
  end subroutine stiff_solver_deallocate
