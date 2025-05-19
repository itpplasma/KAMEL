module params

  integer :: nsize, lmax, jmax, jpmax, mnmax
  double complex :: bx1, bx20, gam, bx4

end module params


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

   use params, only: bx1, bx20, gam, bx4, nsize, lmax, jmax, jpmax, mnmax
   implicit none
   double complex, parameter :: imun=(0.d0,1.d0)
   integer :: j, jp, l, mn, k, jjpmax
   double precision :: x, tau, t
   double precision, dimension(nsize) :: y, f
   double complex :: expon, subint, ax4l, argexp, exponpow, exponax4, exponj
   double complex :: aj, ajp, ajmn

   tau=x
   t=exp(-tau)
!
   ajmn=(1.d0-t)/(1.d0+gam*(1.d0+t))
   aj  =(1.d0+t)/(1.d0+gam*(1.d0+t))
   ajp =(1.d0-t)/(1.d0+gam*(1.d0-t))
   ax4l=1.d0/sqrt(1.d0+imun*bx4*tau)
   argexp=(imun*bx20-bx1**2)*tau+bx1**2*(1+2.d0*gam)*ajmn
   aj  =aj/ajmn**2
   ajp =ajp/ajmn**2
!
   expon = exp(argexp)/sqrt((1.d0+gam)**2-(gam*t)**2)*ax4l
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



subroutine evaluate_integral(wintegral)

  use params
  implicit none
  double precision, parameter :: eps = 1.0d-12, taumax=30.d0
  logical :: switch
  integer :: k, mn, l, j, jp, jjpmax, Nterms, ierr
  double precision, dimension(:), allocatable :: y
  double precision :: x_ini, x_end, est_err
  double complex :: argexp,expon,z,sqxbfac,sqix4fac,oneovopg,znumer,zdenom,G,s,sqonepix4t
  double complex, dimension(0:mnmax,lmax,0:jmax,0:jpmax) :: wintegral
!
  external :: dwdtau
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
    argexp=-sqxbfac*taumax+bx1**2*(1.d0+2.d0*gam)/(1.d0+gam)
    expon=exp(argexp)
    sqonepix4t=1.d0/sqrt(1.d0+(0.d0,1.d0)*bx4*taumax)
    znumer=sqxbfac*(1.d0+(0.d0,1.d0)*bx4*taumax)
    sqxbfac=sqrt(sqxbfac)
    zdenom=(0.d0,1.d0)*bx4
    if(abs(zdenom).gt.abs(znumer)*eps) then
      z=znumer/zdenom
      sqix4fac=1.d0/sqrt(zdenom)
!
      do l=1,lmax
        s=1.d0-0.5d0*dfloat(l+1)
!
        call incomplete_gamma(s,z,eps,eps,G,est_err,Nterms,ierr)
!
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
  integer :: mpar, lperp, m, k, kp, n, j, jp, l, lp
  double precision :: x1,x2,x3,x4
  double complex :: onemin4ix3, onemin4ix3pow, bx1ovi, omix3, sqomix3
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
