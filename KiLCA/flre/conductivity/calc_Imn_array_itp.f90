!<Calculates array of Imn special functions.

!------------------------------------------------------------------------------

subroutine calc_Imn_array(x1_in, x2_in, Imn)

implicit none;

integer, parameter :: dp = 8
integer, parameter :: dpc = 8
complex(8), parameter :: I = cmplx(0.0d0, 1.0d0, 8), one = cmplx(1.0d0, 0.0d0, 8);

complex(dpc) :: x1_in, x2_in;
complex(dpc), dimension(0:3,0:3) :: Imn

complex(dpc) :: t1, t2, F11m, x1, x2
real(dp) :: F_im, F_re;
integer :: l, m, n;
integer, parameter :: nmax = 2

complex(dpc), dimension(0:3,0:3,0:nmax) :: W2;

x1 = x1_in;

t1 = x1*x1;

do l = 0,nmax

    x2 = x2_in + I*l;

    t2 = - I*x2 + t1;

    call hypergeometric1f1_cont_fract_1_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);

    F11m = F_re + I*F_im;

!asymptotic form of expressions: has no fake singularity at kp=0

W2(0,0,l) =  &
(-I)*(-x2**2 + 2 + 5*x1**2 -  &
   (3*I)*x2*(1 + x1**2) + (3 + F11m)*x1**4)/ &
 ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
  (x2 + I*(2 + x1**2)))

W2(0,1,l) =  &
(-I)*x1*(x2*(-1 + F11m*x1**2) -  &
   I*(2 + 3*x1**2+ x1**4))/ &
 ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
  (x2 + I*(2 + x1**2)))

W2(0,2,l) =  &
-((x2 + I)*(2 + 3*x1**2 -  &
    I*x2*(1 - F11m*x1**2) + x1**4))/ &
  ((x2 + I*x1**2)*(x2 + I*(1 + x1**2))* &
   (x2 + I*(2 + x1**2)))

W2(0,3,l) =  &
(-I)*x1*(F11m*x2**3- (3 + 2*F11m)*x2 +  &
   I*x2**2*(3*F11m - x1**2) -  &
   I*(6 + (7 + 2*F11m)*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,0,l) =  &
(-I)*x1*(x2*(-1 + F11m*x1**2) -  &
   I*(2 + 3*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,1,l) =  &
(-I)*x2*(x2*(-1 + F11m*x1**2) -  &
   I*(2 + 3*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,2,l) =  &
(-I)*x1*(F11m*x2**3+ I*x2**2*(F11m - x1**2) -  &
   x2*(3 + 2*x1**2) -  &
   I*(2 + 3*x1**2 + x1**4))/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

W2(1,3,l) =  &
x2*((-I)*F11m*x2**3+ I*(3 + 2*F11m)*x2 -  &
   6 - (7 + 2*F11m)*x1**2 +  &
   x2**2*(3*F11m - x1**2) - x1**4)/ &
 ((x2+ I*x1**2)*(x2+ I*(1 + x1**2))* &
  (x2+ I*(2 + x1**2)))

end do

W2(2,2,0) = I*x1*W2(1,2,1) + 2.d0*W2(1,1,1) - I*x1*W2(1,2,0) + W2(0,2,0)
W2(2,2,1) = I*x1*W2(1,2,2) + 2.d0*W2(1,1,2) - I*x1*W2(1,2,1) + W2(0,2,1)
W2(2,3,0) = I*x1*W2(2,2,1) + 2.d0*W2(1,2,1) - I*x1*W2(2,2,0) + 2.0d0*W2(1,2,0)
W2(2,3,1) = I*x1*W2(1,3,2) + 3.d0*W2(1,2,2) - I*x1*W2(1,3,1) + W2(0,3,1)
W2(3,3,0) = I*x1*W2(2,3,1) + 3.d0*W2(2,2,1) - I*x1*W2(2,3,0) + 2.0d0*W2(1,3,0)

Imn(0,0) = W2(0,0,0);
Imn(0,1) = W2(0,1,0);
Imn(0,2) = W2(0,2,0);
Imn(0,3) = W2(0,3,0);

Imn(1,0) = W2(1,0,0);
Imn(1,1) = W2(1,1,0);
Imn(1,2) = W2(1,2,0);
Imn(1,3) = W2(1,3,0);

Imn(2,0) = Imn(0,2);
Imn(2,1) = Imn(1,2);
Imn(2,2) = W2(2,2,0);
Imn(2,3) = W2(2,3,0);

Imn(3,0) = Imn(0,3);
Imn(3,1) = Imn(1,3);
Imn(3,2) = Imn(2,3);
Imn(3,3) = W2(3,3,0);

!Imn = W2(:,:,0);

!do m = 0,3
!    do n = 0,m-1
!        Imn(m,n) = Imn(n,m)
!    end do
!end do

end subroutine

!------------------------------------------------------------------------------
