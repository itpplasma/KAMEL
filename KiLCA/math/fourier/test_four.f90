program test_fourier

implicit none;

real(8), parameter :: pi = 3.1415926
complex(8), parameter :: I = cmplx(0.0d0, 1.0d0)

integer, parameter :: dimx = 240, dimxf = 500
complex(8), dimension(dimx) :: J, S1, S2, E, fac
real(8), dimension(dimx) :: r

real(8) :: s, f, r0, d
integer, parameter :: M = 50;
complex(8), dimension(-M:M) :: S1t, S2t, Et, Jt
complex(8), dimension(-2*M:2*M) :: Jtt
complex(8), dimension(-M:M) :: Jttt

complex(8), dimension(dimxf) :: Jf, Jff, Jfff
real(8), dimension(dimxf) :: rf

integer :: k, k2

external :: calc_direct_fourier, calc_inverse_fourier

r0 = 35.0;
d = 3.0;
s = 2.0;
f = -log(0.01d0/abs(s-1));

do k = 1,dimx
    r(k) = (k-1)*60.0/(dimx-1);
    fac(k) = 1.0 + (s-1)*exp(-f*(r(k)-r0)**2/(d**2));
    S1(k) = 100 + (r(k)-r0+10)**2 + I*((r(k)-r0-10) - 500 + 2.3*(r(k)-r0+10)**2);
    S2(k) = S1(k)*fac(k);
    E(k) = r(k)*sin(r(k)*pi/60.0)*sin(0.1*(r(k)-r0)**2)*exp(0.1*(r(k)-r0)) + &
           I*r(k)**2*(cos(r(k)*2*pi/60.0)-1.0)*sin(0.1*(r(k)-r0)**2)*exp(0.1*(r(k)-r0));
    J(k) = (S1(k) - S2(k))*E(k);
    write (10, *) r(k), real(S1(k)), aimag(S1(k)), real(S2(k)), aimag(S2(k)), real(E(k)), aimag(E(k)), real(J(k)), aimag(J(k))
end do

call calc_direct_fourier (%val(dimx), r, S1, %val(M), %val(r0), %val(d), S1t);
call calc_direct_fourier (%val(dimx), r, S2, %val(M), %val(r0), %val(d), S2t);
call calc_direct_fourier (%val(dimx), r, E, %val(M), %val(r0), %val(d), Et);
call calc_direct_fourier (%val(dimx), r, J, %val(M), %val(r0), %val(d), Jt);

do k = -M, M
    write (30, *) k, real(S1t(k)), aimag(S1t(k)), real(S2t(k)), aimag(S2t(k)), real(Et(k)), aimag(Et(k)), real(Jt(k)), aimag(Jt(k))
end do;

do k2 = -2*M,2*M
    Jtt(k2) = cmplx(0.0d0,0.0d0)
    do k = max(-M,k2-M),min(k2+M,M)
        Jtt(k2) = Jtt(k2) + (S1t(k)-S2t(k))*Et(k2-k);
    end do
    write (40, *) k2, real(Jtt(k2)), aimag(Jtt(k2))
end do

Jttt = Jtt(-M:M)

do k = 1,dimxf
    rf(k) = (k-1)*60.0/(dimxf-1);
end do

call calc_inverse_fourier (%val(dimxf), rf, Jt, %val(M), %val(r0), %val(d), Jf);

call calc_inverse_fourier (%val(dimxf), rf, Jtt, %val(2*M), %val(r0), %val(d), Jff);

call calc_inverse_fourier (%val(dimxf), rf, Jttt, %val(M), %val(r0), %val(d), Jfff);

do k = 1,dimxf
    write (20, *) rf(k), real(Jf(k)), aimag(Jf(k)), real(Jff(k)), aimag(Jff(k)), real(Jfff(k)), aimag(Jfff(k))
end do;

stop
end program
