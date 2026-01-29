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

program test

   use params
   implicit none
   logical :: consenergy
   integer :: i, nmax, mpar, mperp, lperp, m, j, jp, l, n
   double precision :: hx, x1, x2, x3, x4
   double complex, parameter :: imun=(0.d0,1.d0)
   double complex :: z1,z2
   double complex, dimension(0:3,0:3) :: Imn
   double complex, dimension(:,:,:), allocatable  :: Imnl
   double complex, dimension(:,:,:,:), allocatable  :: Imnkl
   double precision, dimension(:), allocatable :: doubfac
!
   mpar=3
   mperp=2
   lperp=2*mperp+1
   allocate(Imnl(0:mpar,0:mpar,1:lperp))
   allocate(Imnkl(0:mpar,0:mpar,0:mperp,0:mperp))
   allocate(doubfac(-2:2*mpar-1))
!
   doubfac=0.d0
   doubfac(-1)=1.d0
   do m=1,mpar
     doubfac(2*m-1)=doubfac(2*m-3)*dfloat(2*m-1)
   enddo

!lc = -1
 x1 = -189.02963851415251
 x2 =  46753.511704959514
 x3 = -0.54539986769432813
 x4 =  0.56254880990715950

call GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl) ! comment here to get correct values!!!
print *, 'first:';
print *, Imnl;

!lc = 0
 x1 = -189.02963851415251
 x2 = -48.812002536152967
 x3 = -0.54539986769432813
 x4 =  0.56254880990715950

call GreenMom(mpar,lperp,x1,x2,x3,x4,Imnl)
print *, 'second:';
print *, Imnl;

     write(100,*) x1, real(Imnl(0,0,1)), dimag(Imnl(0,0,1))
     write(101,*) x1, real(Imnl(0,1,1)), dimag(Imnl(0,1,1))
     write(102,*) x1, real(Imnl(0,2,1)), dimag(Imnl(0,2,1))
     write(103,*) x1, real(Imnl(0,3,1)), dimag(Imnl(0,3,1))
     write(110,*) x1, real(Imnl(1,0,1)), dimag(Imnl(1,0,1))
     write(111,*) x1, real(Imnl(1,1,1)), dimag(Imnl(1,1,1))
     write(112,*) x1, real(Imnl(1,2,1)), dimag(Imnl(1,2,1))
     write(113,*) x1, real(Imnl(1,3,1)), dimag(Imnl(1,3,1))
     write(120,*) x1, real(Imnl(2,0,1)), dimag(Imnl(2,0,1))
     write(121,*) x1, real(Imnl(2,1,1)), dimag(Imnl(2,1,1))
     write(122,*) x1, real(Imnl(2,2,1)), dimag(Imnl(2,2,1))
     write(123,*) x1, real(Imnl(2,3,1)), dimag(Imnl(2,3,1))
     write(130,*) x1, real(Imnl(3,0,1)), dimag(Imnl(3,0,1))
     write(131,*) x1, real(Imnl(3,1,1)), dimag(Imnl(3,1,1))
     write(132,*) x1, real(Imnl(3,2,1)), dimag(Imnl(3,2,1))
     write(133,*) x1, real(Imnl(3,3,1)), dimag(Imnl(3,3,1))
!
     consenergy=.true.
!
     call getIfunc_drift(consenergy,mpar,mperp,x1,x2,x3,x4,Imnkl)
!
     z1=Imnkl(2,1,0,0)
     z2=Imnkl(3,0,0,0)
     write(1001,*) x1,abs(z1-z2)/(abs(z1)+abs(z2)),real(z1-z2),dimag(z1-z2)
!
     write(200,*) x1, real(Imnkl(0,0,0,0)), dimag(Imnkl(0,0,0,0))
     write(201,*) x1, real(Imnkl(0,1,0,0)), dimag(Imnkl(0,1,0,0))
     write(202,*) x1, real(Imnkl(0,2,0,0)), dimag(Imnkl(0,2,0,0))
     write(203,*) x1, real(Imnkl(0,3,0,0)), dimag(Imnkl(0,3,0,0))
     write(210,*) x1, real(Imnkl(1,0,0,0)), dimag(Imnkl(1,0,0,0))
     write(211,*) x1, real(Imnkl(1,1,0,0)), dimag(Imnkl(1,1,0,0))
     write(212,*) x1, real(Imnkl(1,2,0,0)), dimag(Imnkl(1,2,0,0))
     write(213,*) x1, real(Imnkl(1,3,0,0)), dimag(Imnkl(1,3,0,0))
     write(220,*) x1, real(Imnkl(2,0,0,0)), dimag(Imnkl(2,0,0,0))
     write(221,*) x1, real(Imnkl(2,1,0,0)), dimag(Imnkl(2,1,0,0))
     write(222,*) x1, real(Imnkl(2,2,0,0)), dimag(Imnkl(2,2,0,0))
     write(223,*) x1, real(Imnkl(2,3,0,0)), dimag(Imnkl(2,3,0,0))
     write(230,*) x1, real(Imnkl(3,0,0,0)), dimag(Imnkl(3,0,0,0))
     write(231,*) x1, real(Imnkl(3,1,0,0)), dimag(Imnkl(3,1,0,0))
     write(232,*) x1, real(Imnkl(3,2,0,0)), dimag(Imnkl(3,2,0,0))
     write(233,*) x1, real(Imnkl(3,3,0,0)), dimag(Imnkl(3,3,0,0))


end
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
