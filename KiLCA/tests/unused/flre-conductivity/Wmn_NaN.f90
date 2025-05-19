!<The definitions of some constants used throughout the library.

module constants

implicit none

integer, parameter :: pp = 8 !integers size for C pointers, 32 bit - 4, 64 bit - 8 and so on.

integer, parameter :: sp = 4
integer, parameter :: dp = 8
integer, parameter :: qp = 16 !quad precision
integer, parameter :: vp = dp !

integer, parameter :: spc = 4
integer, parameter :: dpc = 8
integer, parameter :: qpc = 16 !quad precision
integer, parameter :: vpc = dpc !
integer, parameter :: lgt = kind (.true.)

real(dp), parameter :: pi    = 3.141592653589793238462643383279502884197_dp
real(dp), parameter :: pio2  = 1.57079632679489661923132169163975144209858_dp
real(dp), parameter :: twopi = 6.283185307179586476925286766559005768394_dp
real(dp), parameter :: sqrt2 = 1.41421356237309504880168872420969807856967_dp
real(dp), parameter :: euler = 0.5772156649015328606065120900824024310422_dp
real(dp), parameter :: boltz = 1.60216428d-12 ! erg/ev

real(dp), parameter :: srpi = sqrt(pi)
real(dp), parameter :: sr2 = sqrt(2.0d0)
real(dp), parameter :: sqrt2p = sqrt(2.0d0*pi)
real(dp), parameter :: isqrt2pi = 1.0d0/sqrt2p
real(dp), parameter :: tppoh = (2.0d0*pi)**(1.5d0)

complex(dpc), parameter :: im = (0.0d0, 1.0d0)

real(dp), parameter :: c  = 29979245800.0       !speed of light in vacuum
real(dp), parameter :: mp = 1.67262158d-24      !proton mass
real(dp), parameter :: me = mp/1.8361526675d3 !electron mass
real(dp), parameter :: e  = 4.8032d-10        !elementary charge

integer :: PInf, NaN

data Pinf /B'01111111100000000000000000000000'/      ! +Infinity
data NaN  /B'01111111100000100000000000000000'/      ! NaN

end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module params

  integer :: nsize, lmax, jmax, jpmax, mnmax
  double complex :: bx1, bx20, gam, bx4, theta

end module params

module help_igamma_mod
  double complex :: zbeg,zend
end module help_igamma_mod

!-----------------------------------------------------------------------------!

subroutine incomplete_gamma(s, z, abs_err, rel_err, G, est_err, Nterms, status)

! Evaluates G(s, z) by continued fractions method.
! abs_err, rel_err - desired tolerances
! G - G(s, z) value
! est_err - estimated relative error
! Nterms - number of terms used
! status == 0 - Ok.

implicit none;

complex(8), intent(in)   :: s, z;
real(8),    intent(in)   :: abs_err, rel_err;
complex(8), intent(out)  :: G;
real(8),    intent(out)  :: est_err;
integer,    intent(out)  :: Nterms;
integer,    intent(out)  :: status;

integer,       parameter :: prec = 8;
integer,       parameter :: Nmin = 2, Nmax = 4194304;
complex(prec), parameter :: NaN = cmplx(1.0d100, 1.0d100, prec);
complex(prec), parameter :: c0 = cmplx(0.0d0, 0.0d0, prec), c1 = cmplx(1.0d0, 0.0d0, prec);

complex(prec), dimension(Nmax) :: a, b;

complex(prec) :: term, F, Fp;

integer :: k, j, M, N;

logical :: goon;

M = 0;
N = 4;

F  = NaN;
Fp = NaN;

goon = .true.;

do while (N <= Nmax .and. goon)

    term = c0; ! tail value

    do k = N,Nmin,-1

        if (k > M) then

            j = k / 2;

            a(k) = cmplx(j, 0, prec);

            if (j * 2 .eq. k) a(k) = a(k) - s;

            b(k) = c1;

        end if

        term = a(k) / z / (b(k) + term);

    end do

    Fp = F;

    F = c1 / z / (c1 + term); ! top level term

    goon = ( abs(F - Fp) > abs_err + rel_err * abs(F) );

    if (N > M) M = N;

    N = N * 2;

end do

est_err = abs(c1 - Fp/F);

if (goon) then
    status = 1;
else
    status = 0;
end if

Nterms = N / 2;

!G = z**s * exp(-z) * F;
G = z**s * F;

end subroutine

!-----------------------------------------------------------------------------!

!--- subroutines ---------------------------------------------------------------

subroutine dwdtau(x, y, f)

   use params, only: bx1, bx20, gam, bx4, nsize, lmax, jmax, jpmax, mnmax, theta
   implicit none
   double complex, parameter :: imun=(0.d0,1.d0)
   integer :: j, jp, l, mn, k, jjpmax
   double precision :: x
   double complex   :: tau, t
   double precision, dimension(nsize) :: y, f
   double complex :: expon, subint, ax4l, argexp, exponpow, exponax4, exponj
   double complex :: aj, ajp, ajmn

   tau=x*theta
   t=exp(-tau)
!
   ajmn=(1.d0-t)/(1.d0+gam*(1.d0+t))
   aj  =(1.d0+t)/(1.d0+gam*(1.d0+t))
   ajp =(1.d0-t)/(1.d0+gam*(1.d0-t))
   ax4l=1.d0/sqrt(1.d0+imun*bx4*tau)
   argexp=(imun*bx20-bx1**2)*tau+bx1**2*(1.d0+2.d0*gam)*ajmn
   aj  =aj/ajmn**2
   ajp =ajp/ajmn**2
!
   expon = exp(argexp)/sqrt((1.d0+gam)**2-(gam*t)**2)*ax4l
   expon = expon*theta
!
   k=0
   do mn=0,mnmax
     jjpmax=mn/2
     exponax4=expon
     do l=1,lmax
       expon=expon*ax4l
       exponj=expon
       do j=0,jjpmax
         exponpow=exponj
         do jp=0,jjpmax-j
           k=k+1
           f(k)=real(exponpow)
           k=k+1
           f(k)=dimag(exponpow)
           exponpow=exponpow*ajp
         end do
         exponj=exponj*aj
       end do
     end do
     expon=exponax4*ajmn
   end do

end subroutine dwdtau

subroutine help_igamma(x, y, f)

  use params, only: lmax
  use help_igamma_mod, only : zbeg,zend
  implicit none
  integer :: l, k
  double precision :: x
  double complex   :: z,dzdx
  double precision, dimension(2*lmax) :: y, f
  double complex :: expon, sqinv
  double complex :: aj, ajp, ajmn
!
  dzdx=zend-zbeg
  z=zbeg+dzdx*x
  sqinv=(1.d0,0.d0)/sqrt(z)
  expon=exp(zbeg-z)*sqinv*dzdx
  k=0
  do l=1,lmax
    expon=expon*sqinv
    k=k+1
    f(k)=real(expon)
    k=k+1
    f(k)=dimag(expon)
  enddo

end subroutine help_igamma

subroutine evaluate_integral(wintegral)

  use params
  use help_igamma_mod, only : zbeg,zend
  implicit none
  double precision, parameter :: eps = 1.0d-10, taumax=30.d0, zshift=0.1d0
  logical :: switch
  integer :: k, mn, l, j, jp, jjpmax, Nterms, ierr
  double precision, dimension(:), allocatable :: y
  double precision :: x_ini, x_end, est_err
  double complex :: argexp,expon,z,sqxbfac,sqix4fac,oneovopg,znumer,zdenom
  double complex :: G,s,sqonepix4t,Gscale
  double complex, dimension(0:mnmax,lmax,0:jmax,0:jpmax) :: wintegral
  double complex, dimension(:), allocatable :: delG
!
  external :: dwdtau,help_igamma
!
  k=0
  do mn=0,mnmax
    jjpmax=mn/2
    do l=1,lmax
      do j=0,jjpmax
        do jp=0,jjpmax-j
          k=k+2
        end do
      end do
    end do
  end do
!
  nsize=k
  allocate(y(nsize))
!
!
  x_ini=1.d-14
  x_end = 20.0d0/abs(bx1) + 400.0d0/abs(bx1)**2
  switch=x_end.gt.taumax
  if(switch) x_end=taumax
  y = 0.0d0
!
  call odeint_allroutines(y, nsize, x_ini, x_end, eps, dwdtau)
!
  k=0
  do mn=0,mnmax
    jjpmax=mn/2
    do l=1,lmax
      do j=0,jjpmax
        do jp=0,jjpmax-j
          k=k+1
          wintegral(mn,l,j,jp)=dcmplx(y(k),y(k+1))
          k=k+1
        end do
      end do
    end do
  end do

  deallocate(y)
!
  if(switch) then
    oneovopg=1.d0/(1.d0+gam)
    sqxbfac=bx1**2-(0.d0,1.d0)*bx20
    argexp=-sqxbfac*taumax*theta+bx1**2*(1.d0+2.d0*gam)/(1.d0+gam)
    expon=exp(argexp)
    sqonepix4t=1.d0/sqrt(1.d0+(0.d0,1.d0)*bx4*taumax*theta)
    znumer=sqxbfac*(1.d0+(0.d0,1.d0)*bx4*taumax*theta)
    sqxbfac=sqrt(sqxbfac)
    zdenom=(0.d0,1.d0)*bx4
    if(abs(zdenom).gt.abs(znumer)*eps) then
      allocate(delG(lmax))
      z=znumer/zdenom
      sqix4fac=1.d0/sqrt(zdenom)
      if(real(z).lt.0.d0.and.abs(dimag(z)).lt.zshift) then
        zbeg=z
        zend=z+dcmplx(0.d0,sign(zshift,dimag(z)))
        nsize=k
        allocate(y(nsize))
        x_ini=0.d0
        x_end=1.d0
        y=0.d0
!
        call odeint_allroutines(y, nsize, x_ini, x_end, eps, help_igamma)
!
        k=0
        do l=1,lmax
          delG(l)=dcmplx(y(k+1),y(k+2))
          k=k+2
        enddo
        deallocate(y)
        z=zend
        Gscale=exp(zbeg-zend)
      else
        delG=(0.d0,0.d0)
        Gscale=(1.d0,0.d0)
      endif
!
      do l=1,lmax
        s=1.d0-0.5d0*dfloat(l+1)
!
        call incomplete_gamma(s,z,eps,eps,G,est_err,Nterms,ierr)
!
        G=G*Gscale+delG(l)
        do mn=0,mnmax
          jjpmax=mn/2
          do j=0,jjpmax
            do jp=0,jjpmax-j
              wintegral(mn,l,j,jp)=wintegral(mn,l,j,jp)+expon*G         &
                                  *oneovopg**(mn+1-j-jp)*sqxbfac**(l-1) &
                                  *sqix4fac**(l+1)
            end do
          end do
        end do
      end do
!
      deallocate(delG)
    else
!
      do l=1,lmax
        do mn=0,mnmax
          jjpmax=mn/2
          do j=0,jjpmax
            do jp=0,jjpmax-j
              wintegral(mn,l,j,jp)=wintegral(mn,l,j,jp)+expon         &
                                  *oneovopg**(mn+1-j-jp)/znumer       &
                                  *sqonepix4t**(l-1)
            end do
          end do
        end do
      end do
!
    endif
  endif
!
end subroutine evaluate_integral

subroutine GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl)

  use params

  implicit none
  double precision, parameter :: addimpart=0.1d0,epscontour=0.03d0
  integer :: mpar, lperp, m, k, kp, n, j, jp, l, lp
  double precision :: x1,x2,x3,x4
  double complex :: onemin4ix3, onemin4ix3pow, bx1ovi, omix3, sqomix3
  double complex :: oldarg
  double complex, dimension(0:mpar,0:mpar,1:lperp) :: Imnl

  double precision, dimension(:),       allocatable :: doubfac
  double precision, dimension(:,:),     allocatable :: bincoef
  double precision, dimension(:,:,:,:), allocatable :: symbm
  double complex,   dimension(:,:,:,:), allocatable :: calsymbm
  double complex,   dimension(:,:,:,:), allocatable :: wintegral

  mnmax=2*mpar
  lmax=lperp
  jmax=mpar
  jpmax=mpar

  allocate(doubfac(0:mpar))
  allocate(bincoef(0:mpar,0:mpar))
  allocate(symbm(0:mpar,0:mpar,0:mpar,0:mpar))
  allocate(calsymbm(0:mpar,0:mpar,0:mpar,0:mpar))
  allocate(wintegral(0:mnmax,lmax,0:jmax,0:jpmax))
!
  doubfac(0)=1.d0
  do m=1,mpar
    doubfac(m)=doubfac(m-1)*dfloat(2*m-1)
  enddo
!
  bincoef=0.d0
  do m=0,mpar
    bincoef(m,0)=1.d0
    do k=1,m
      bincoef(m,k)=bincoef(m,k-1)*dfloat(m-k+1)/dfloat(k)
    enddo
  enddo
!
  symbm=0.d0
!
  do m=0,mpar
    do n=0,mpar
      do j=0,mpar
        do jp=0,mpar
          symbm(m,n,j,jp)=0.d0
          do k=0,m
            do kp=0,n
              do l=0,k
                do lp=0,kp
                  if (2*j.eq.l+lp.and.2*jp.eq.k+kp-l-lp) then
                    symbm(m,n,j,jp)=symbm(m,n,j,jp)+(-1.d0)**(k-l) &
                                   *bincoef(m,k)*bincoef(n,kp)     &
                                   *bincoef(k,l)*bincoef(kp,lp)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
!
  omix3=1.d0-(0.d0,4.d0)*x3
  sqomix3=sqrt(omix3)
  gam=(0.d0,2.d0)*x3/(omix3+sqomix3)
  onemin4ix3=1.d0/sqrt(sqomix3)
  bx1=x1/sqomix3*onemin4ix3
  bx1ovi=(0.d0,-1.d0)*bx1
  bx20=(x2+2.d0*x3/(sqomix3+1.d0))/sqomix3
  bx4=x4/sqomix3
  oldarg=(0.d0,1.d0)*bx20-bx1**2
  theta=dcmplx(1.d0,sign(addimpart                                            &
       +max(0.d0,real(oldarg)/max(abs(dimag(oldarg)),abs(oldarg)*epscontour)) &
       ,dimag(oldarg)))
!
  do m=0,mpar
    do n=0,mpar
      onemin4ix3pow=onemin4ix3**(m+n+3)
      do j=0,mpar
        do jp=0,mpar
          calsymbm(m,n,j,jp)=symbm(m,n,j,jp)*doubfac(j)*doubfac(jp) &
                            *0.5d0**(j+jp)*bx1ovi**(m+n-2*(j+jp))   &
                            *onemin4ix3pow
        enddo
      enddo
    enddo
  enddo
!
  call  evaluate_integral(wintegral)
!
  Imnl=(0.d0,0.d0)
  do m=0,mpar
    do n=0,mpar
      do l=1,lperp
        do j=0,mpar
          do jp=0,mpar
            Imnl(m,n,l)=Imnl(m,n,l)+wintegral(m+n,l,j,jp)*calsymbm(m,n,j,jp)
          enddo
        enddo
      enddo
    enddo
  enddo


end subroutine GreenMom

subroutine getIfunc_drift(consenergy,mpar,mperp,x1,x2,x3,x4,Imnkl)
!
  implicit none
!
  logical :: consenergy
  integer :: mpar,mperp,lperp,l,m,n,k
  double precision :: x1, x2, x3, x4, arggam
  double complex, dimension(0:mpar,0:mpar,0:mperp,0:mperp)  :: Imnkl
  double precision, dimension(:),     allocatable :: gamkl
  double complex,   dimension(:,:,:), allocatable :: Imnl
!
  lperp=2*mperp+1
  allocate(Imnl(0:mpar,0:mpar,1:lperp))
  if(consenergy) then
    allocate(gamkl(0:2*mperp))
    do k=0,2*mperp
      arggam=0.5d0*dfloat(k)+1
      gamkl(k)=gamma(arggam)
    enddo
  endif
!
  call GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl)
!
  do k=0,mperp
    do l=0,mperp
      do m=0,mpar
        do n=0,mpar
          if(consenergy) then
            Imnkl(m,n,k,l)=Imnl(m,n,k+l+1)+gamkl(k)*gamkl(l)/gamkl(k+l)      &
                *(Imnl(m,2,k+1)-Imnl(m,0,k+1))*(Imnl(n,2,l+1)-Imnl(n,0,l+1)) &
                /(1.d0-Imnl(0,0,1)+2.d0*Imnl(2,0,1)-Imnl(2,2,1))
          else
            Imnkl(m,n,k,l)=Imnl(m,n,k+l+1)
          endif
        enddo
      enddo
    enddo
  enddo
!
  deallocate(Imnl)
  if(consenergy) deallocate(gamkl)
end subroutine getIfunc_drift
!--------------------------------------------
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module get_matrix_tmp
    integer :: ispec
  end module get_matrix_tmp
!
!-----------------------------------------------------------------------------
!

program test

   use params
   implicit none
   logical :: consenergy
   integer :: i, nmax, mpar, mperp, lperp, m, j, jp, l, n, l1, l2
   double precision :: hx, x1, x2, x3, x4
   double complex, parameter :: imun=(0.d0,1.d0)
   double complex :: z1,z2
   double complex, dimension(0:3,0:3) :: Imn
   double complex, dimension(:,:,:), allocatable  :: Imnl
   double complex, dimension(:,:,:), allocatable  :: IImnl
   double complex, dimension(:,:,:,:), allocatable  :: IImnkl
   double complex, dimension(:,:,:,:), allocatable  :: SImnkl
   double precision, dimension(:), allocatable :: doubfac
   double precision :: null = 0.0d0
   real(8), parameter :: pi    = 3.141592653589793238462643383279502884197
   real(8), parameter :: sqrt2pi = sqrt(2.0d0*pi)
   double complex :: x1c, x2c, x3c, x4c;
   integer :: status
!
   mpar=5
   mperp=9
   lperp=2*mperp+1

   allocate(Imnl(0:mpar,0:mpar,1:lperp))
   allocate(IImnkl(0:mpar,0:mpar,0:mperp,0:mperp))
   allocate(SImnkl(0:mpar,0:mpar,0:mperp,0:mperp))
   allocate(doubfac(-2:2*mpar-1))

   allocate(IImnl(0:mpar,0:mpar,1:lperp))
!
   doubfac=0.d0
   doubfac(-1)=1.d0
   do m=1,mpar
     doubfac(2*m-1)=doubfac(2*m-3)*dfloat(2*m-1)
   enddo

!lc = -1
 x1 = -189.02963851415251
 x2 = -46753.511704959514
 x3 = -0.54539986769432813
 x4 =  0.56254880990715950

!call GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl) ! comment here to get correct values!!!
!print *, 'first:';
!print *, Imnl;

!lc = 0
 x1 = -189.02963851415251
 x2 = -48.812002536152967
 x3 = -0.54539986769432813
 x4 =  0.56254880990715950

 x1c = -189.02963851415251
 x2c = -48.812002536152967
 x3c = -0.54539986769432813
 x4c =  0.56254880990715950

 x1 = -189.02963851415251
 x2 = -46753.511704959514
 x3 = -0.54539986769432813
 x4 =  0.56254880990715950

 x1c= -189.02963851415251
 x2c= -46753.511704959514
 x3c = -0.54539986769432813
 x4c =  0.56254880990715950


x1 = 24.963007776733463
x2 = 1.1620797243229703
x3 = 7.6855794633391580e-2
x4 = 4.4250943880273737e-2

x1c = 24.963007776733463
x2c = 1.1620797243229703
x3c = 7.6855794633391580e-2
x4c = 4.4250943880273737e-2

x1 = -3.4050112348510937
x2 = -1.2197902823671228E-002
x3 = 1.9895084243796345E-004
x4 = 2.0549011000562675E-004

! x1 = 5.8921523459531224
! x2 = 5955.9014672252561
! x3 = -8.6872069670694129E-004
! x4 = -7.5712515574091113E-003

 x1 = 1.6287848200204527
 x2 = 116442.47103820778
 x3 = 6.3885834695745553E-006
 x4 = 2.0305545715023568E-005

! x1 = 1.6d0
! x2 = 116442
! x3 = 64.0d-7
! x4 = 2.0d-5
!
!  x1 = -140.47226325221121
!  x2 = -19.160769732482777
!  x3 = -0.39004192054644837
!  x4 = 3.5518191772447097E-002

 x1 = -0.15548575383155275
 x2 = -1.3981937829735058
 x3 = -2.0838042516747982E-002
 x4 = -2.6022827172207063E-002

  x1 = -189.02963851415251
  x2 = -48.812002536152967
  x3 = -0.54539986769432813
  x4 =  0.56254880990715950
!
!  x1 =  -54.373205158510601
!  x2 = -913179.34786194994
!  x3 = 2.4474542564594847E-003
!  x4 = -9.7688538680105260E-004

! x1 = 0.39733225332036942
! x2 = 0.37656301615958365
! x3 = -1.1162421484418730E-005
! x4 = -1.8188303353704633E-004

! x1 = -189.02963851415251
! x2 = -46851.135710031806
! x3 = -0.54539986769432813
! x4 = 0.56254880990715950
!
! x1 = -1.4952627531160615E-002
! x2 = -1.3923571119795708
! x3 = -2.0539733101908705E-002
! x4 = -2.5716584715356681E-002

!-0.10345896294773248       -1.3960314619430063       -2.0727622585503425E-002  -2.5909647944190106E-002
!-1.4952627531160615E-002  -1.3923571119795708       -2.0539733101908705E-002  -2.5716584715356681E-002
!7.2768354702616181E-002  -1.3887218497684983       -2.0353449425747985E-002  -2.5524561381443375E-002

x1 = -57.495244072264008
x2 = -5.1868652101341199
x3 = -0.15436544344311581
x4 = -8.8030256224353190E-002

x1 = -1259.6976813440651
x2 = -93.483351834082541
x3 = -5.2292762150680376
x4 = 4.4567822695718826

x1 = -181.56927999358916
x2 =  15.255808887885236
x3 = -0.51968185069744288
x4 = -0.12449524385361799


x1 =  -77.928847343509034
x2 = -0.11351659097169517
x3 =   3.4814783754558397E-003
x4 =   1.0988952323586173E-004

x1 = -7.6015941985375151E-003
x2 =  1.7417852256132210E-002
x3 =  1.7082185155978689E-004
x4 =  1.7813164197335554E-004

! x1 =  -1227.1450167640457
! x2 =  -71.897352402962127
! x3 =  -5.0570739471295116
! x4 =   3.8893149644463398
!
! x1 =   -0.10345896294773248
! x2 =   -1.3960314619430063
! x3 =   -2.0727622585503425E-002
! x4 =   -2.5909647944190106E-002

! x1 =  7.1801224079603185E-004
! x2 = -6.6402554522991602E-003
! x3 =  8.2209915089726310E-005
! x4 =  1.0297186722561036E-004


!x1 = -2.1990920564759989E-003
!x2 = -6.6425899639394248E-003
!x3 =  8.2310206769808418E-005
!x4 =  1.0307662062406069E-004

x1 =   1.2350118204346608E-002
x2 =  -6.6310020584000257E-003
x3 =   8.1809921768705203E-005
x4 =   1.0255368120469605E-004

x1 = -148.67656526649230
x2 = -0.71005670428972312
x3 = 3.5543917292194943E-003
x4 = -2.3154367418953105E-003

x1 = -254.9545520199325
x2 = 62933.62227644278
x3 = -0.4749892745764920
x4 = 0.3094271751497384

x1 = -254.9545520199325
x2 = 62.93362227644278
x3 = -0.4749892745764920
x4 = 0.3094271751497384

! x1 = 2.5118619905876738E-004
! x2 = 7.6126906346330096E-005
! x3 = 1.0443495244172603E-005
! x4 = 7.9445099866432262E-006
!
! x1 = -0.001477524271151
! x2 = -0.078933966058232
! x3 = -0.002946930602055
! x4 = -0.003382656999166

! x1 = 2.3780830605437656E-004
! x2 = -3.5482644166362005E-004
! x3 = 1.1048107011420668E-005
! x4 = 1.2696539201333587E-005

x1c = x1
x2c = x2
x3c = x3
x4c = x4

!call GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl)

call init_spec_data(mpar, lperp);
!call calc_imnl_quad(mpar, lperp, x1c, x2c, x3c, x4c, IImnl, status);

! do m = 0,mpar
!     do n = 0,mpar
!         do l = 1,lperp
!
!             !print *, '(',m, n, l,') =', IImnl(m,n,l), Imnl(m,n,l);
!
!         end do
!     end do
! end do

consenergy=.true.
consenergy=.false.

call calc_imnkl_quad(consenergy, mpar, mperp, x1c, x2c, x3c, x4c, IImnkl, status);

call getIfunc_drift(consenergy,mpar,mperp,x1,x2,x3,x4,SImnkl)

do m = 0,mpar
    do n = 0,mpar
        do l1 = 0,mperp
            do l2 = 0,mperp

            print *, '(',m, n, l1, l2,') =', IImnkl(m,n,l1,l2), SImnkl(m,n,l1,l2);

            end do
        end do
    end do
end do

print *, maxval(abs((IImnkl-SImnkl)/IImnkl));

call close_spec_data();

end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
