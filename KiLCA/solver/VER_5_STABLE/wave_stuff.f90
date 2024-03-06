!
  module dieltens_mod
    logical :: prop=.true.
    integer :: iunit=71,iswitch_cold_warm=1
    double precision :: R_0,singam,cosgam,ompe2,omph2,ompd2,omce0,omch0,omcd0
    double precision :: anu_e,anu_h,anu_d,omega,vec0,vecy,vecz,vpar,B0_x,B0_z
    double precision :: v_Te,v_Th,v_Td,sq2,sqpi
    double complex   :: omover4pii,eps_1,eps_2,eps_3,vecx
  end module dieltens_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module fields_mod
    double complex :: Ex,Ey,Ez,Bx,By,Bz,curx,cury,curz
  end module fields_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module nueff_mod
    double precision :: xmid,sigma,ampnu
  end module nueff_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module keeproots_mod
    logical :: prop_keep=.true.
    double complex, dimension(:), allocatable :: vecs_keep
  end module keeproots_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module baseparam_mod
    double precision, parameter :: pi=3.14159265358979d0
    double precision,parameter  :: c=2.9979d10
    double precision,parameter  :: e_charge=4.8032d-10
    double precision,parameter  :: e_mass=9.1094d-28
    double precision,parameter  :: p_mass=1.6726d-24
    double precision,parameter  :: ev=1.6022d-12
  end module baseparam_mod
!
  module rhs_mod
    integer :: nrhs_vac,nrhs_plas
  end module rhs_mod
!
  module fulltens_mod
    integer, parameter :: nacc=20
    integer :: nx_max
    double precision :: vecx0,xmargin,delta_band
  end module fulltens_mod
!
  module mhd_sol_mod
    integer :: nrad_mhd,nsols
    double precision, dimension(:),     allocatable :: rad_arr
    double complex,   dimension(:,:,:), allocatable :: efield_arr,curr_arr
  end module mhd_sol_mod
!
  module equilpar_mod
! particle sort, 0-electrons, 1-first sort of ions, 2-second sort of ions:
    integer :: isort
! ion charge numbers, atomic weights, relative density of sort #2 :
    double precision :: Z1,Z2,am1,am2,conc2
! density and its derivative at the resonance zone:
    double precision :: den0,den0_pr
! electron temperature, its derivative, parallel velocity and its derivative:
    double precision :: temp0_e,temp_pr_e, V0_e,V0_pr_e
! first ion temperature, its derivative, parallel velocity and its derivative:
    double precision :: temp0_i1,temp_pr_i1, V0_i1,V0_pr_i1
! second ion temperature, its derivative, parallel velocity and its derivative:
    double precision :: temp0_i2,temp_pr_i2, V0_i2,V0_pr_i2
! collision frequencies of electrons and ions
    double precision :: collfr_e,collfr_i1,collfr_i2
! gardient of electrostatic potential (electric field) and its derivative:
    double precision :: Phi_pr_0,Phi_pr_pr
  end module equilpar_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine dieltens(R,eps)
!
  use dieltens_mod, only : iswitch_cold_warm
!
  implicit none
!
  double precision :: R
  double complex, dimension(3,3) :: eps
!
  call dieltens_cold(R,eps)
!
  end subroutine dieltens
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine dieltens_cold(R,eps)
!
  use baseparam_mod
  use dieltens_mod
  use nueff_mod
!
  implicit none
!
  double complex, parameter :: imun=(0.d0,1.d0)
!
  integer :: isort
  double precision :: R,omce,omch,omcd,fac,B0,gam,den,conc
  double complex :: eps_1e,eps_2e,eps_3e
  double complex :: eps_1h,eps_2h,eps_3h
  double complex :: eps_1d,eps_2d,eps_3d
  double complex, dimension(3,3) :: eps
!
  if(prop) then
    prop=.false.
    open(iunit,file='wave_driver.inp')
    read(iunit,*) vecy    ! k_y, 1/cm
    read(iunit,*) vecz    ! k_z, 1/cm
    close(iunit)
    open(iunit,file='plasparam.inp')
    read(iunit,*) R_0
    read(iunit,*) B0
    read(iunit,*) gam
    read(iunit,*) isort
    read(iunit,*) den
    read(iunit,*) conc
    read(iunit,*) anu_e
    read(iunit,*) anu_h
    read(iunit,*) anu_d
    close(iunit)
    singam=sin(gam)
    cosgam=cos(gam)
    B0_x=singam
    B0_z=cosgam
    omce0=-e_charge*B0/(e_mass*c)
    omch0=e_charge*B0/(p_mass*c)
    omcd0=omch0/2.d0
    ompe2=4.d0*pi*den*e_charge**2/e_mass
    if(isort.eq.1) then
      omega=omch0
      omph2=4.d0*pi*den*conc*e_charge**2/p_mass
      ompd2=4.d0*pi*den*(1.d0-conc)*e_charge**2/p_mass/2.d0
    else
      omega=omcd0
      omph2=4.d0*pi*den*(1.d0-conc)*e_charge**2/p_mass
      ompd2=4.d0*pi*den*conc*e_charge**2/p_mass/2.d0
    endif
    vec0=omega/c
    omover4pii=omega/(4.d0*pi*imun)
  endif
!
  fac=sqrt((B0_z*R_0/R)**2+B0_x**2)
  singam=B0_x/fac
  cosgam=B0_z*R_0/(R*fac)
!
  omce=omce0*fac
  omch=omch0*fac
  omcd=omcd0*fac
!
  eps_1e=ompe2/(omce**2-(omega+imun*anu_e)**2)
  eps_2e=eps_1e*omce/omega
  eps_1e=eps_1e*(omega+imun*anu_e)/omega
  eps_3e=-ompe2/omega/(omega+imun*anu_e)
!
  eps_1h=omph2/(omch**2-(omega+imun*anu_h)**2)
  eps_2h=eps_1h*omch/omega
  eps_1h=eps_1h*(omega+imun*anu_h)/omega
  eps_3h=-omph2/omega/(omega+imun*anu_h)
!
  eps_1d=ompd2/(omcd**2-(omega+imun*anu_d)**2)
  eps_2d=eps_1d*omcd/omega
  eps_1d=eps_1d*(omega+imun*anu_d)/omega
  eps_3d=-ompd2/omega/(omega+imun*anu_d)
!
  eps_1=1.d0+eps_1e+eps_1h+eps_1d
  eps_2=eps_2e+eps_2h+eps_2d
  eps_3=1.d0+eps_3e+eps_3h+eps_3d
!
!  eps_1=eps_1+(0.d0,1.d0)*ampnu*exp(-((R-xmid)/sigma)**2)
!
  eps(1,1)=eps_1*cosgam**2+eps_3*singam**2
  eps(1,2)=imun*eps_2*cosgam
  eps(1,3)=(eps_3-eps_1)*singam*cosgam
  eps(2,1)=-eps(1,2)
  eps(2,2)=eps_1
  eps(2,3)=imun*eps_2*singam
  eps(3,1)=eps(1,3)
  eps(3,2)=-eps(2,3)
  eps(3,3)=eps_1*singam**2+eps_3*cosgam**2
!
  return
  end subroutine dieltens_cold
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine wave_matrix(x,amat)
!
! Right hand side for the 4-th order ODE set ("full wave")
!
  use dieltens_mod
!
  implicit none
!
  double complex, parameter :: imun=(0.d0,1.d0)
!
  double precision :: x
  double complex, dimension(3,3) :: eps
  double complex, dimension(4,4) :: amat
  double complex, dimension(4) :: Ex,Ey,Ez,Bx,By,Bz
!
  Ey=(0.d0,0.d0)
  Ey(1)=(1.d0,0.d0)
  Ez=(0.d0,0.d0)
  Ez(2)=(1.d0,0.d0)
  By=(0.d0,0.d0)
  By(3)=(1.d0,0.d0)
  Bz=(0.d0,0.d0)
  Bz(4)=(1.d0,0.d0)
!
  call dieltens(x,eps)
!
  Ex=((vecz*By-vecy*Bz)/vec0-eps(1,2)*Ey-eps(1,3)*Ez)/eps(1,1)
  Bx=(vecy*Ez-vecz*Ey)/vec0
!
  amat(1,:)=imun*(vecy*Ex+vec0*Bz)
  amat(2,:)=imun*(vecz*Ex-vec0*By)
  amat(3,:)=imun*vecy*Bx-imun*vec0*(eps(3,1)*Ex+eps(3,2)*Ey+eps(3,3)*Ez)
  amat(4,:)=imun*vecz*Bx+imun*vec0*(eps(2,1)*Ex+eps(2,2)*Ey+eps(2,3)*Ez)
!
  end subroutine wave_matrix
