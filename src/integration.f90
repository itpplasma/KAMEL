module integration

    implicit none

    double precision :: eps = 1d-3
    double precision :: h1 = 0.1d0
    double precision :: hmin = 0.0d0
    integer :: num_ok, num_bad

    integer, parameter :: max_steps = 10000
    double precision, parameter :: TINY=1.0d-30


    contains

        subroutine integrate_krp(ystart, nvar, kr1, kr2, j, k1, lp, derivs)!,rkqs) 

            implicit none
        
            integer, intent(in) :: nvar
            integer, intent(in) :: j     ! index for r profile
            integer, intent(in) :: k1    ! index for kr profile
            integer, intent(in) :: lp ! indices for spline basis functions
            double precision, intent(in) :: kr1, kr2 ! limits of integration
            double complex, dimension(nvar), intent(inout) :: ystart

            integer :: i, nstp
            double precision :: h,hdid,hnext,x,xsav
            double complex, dimension(size(ystart)) :: dydx, y, yscal

            interface
                subroutine derivs(r,y, dydr, j, k1, lp)
                    implicit none
                    double precision, intent(in) :: r
                    integer, intent(in) :: j, k1, lp
                    double complex, dimension(:), intent(in) :: y
                    double complex, dimension(:), intent(out) :: dydr
                end subroutine
            end interface

            x = kr1
            h = sign(h1, kr2 - kr1) 
            num_ok  = 0 
            num_bad = 0

            y(:) = ystart(:)

            do nstp = 1, max_steps
                call derivs(x,y,dydx, j, k1, lp)
                !write(*,*) 'y = ', y, '; x = ', x

                yscal(:) = abs(y(:)) + abs(h * dydx(:)) + TINY 

                if((x + h - kr2) * (x+h-kr1) > 0.0) h = kr2-x
                

                call rkqs_c(y,dydx,x,h, yscal,hdid,hnext)!,derivs)

                if(hdid == h)then
                    num_ok=num_ok + 1
                else
                    num_bad=num_bad + 1
                endif
                if((x-kr2)*(kr2-kr1).ge.0.)then 
                    ystart(:)=y(:) 
                return 
                endif
        
                if(abs(hnext).lt.hmin) then 
                    write(*,*) 'stepsize smaller than minimum in odeint'
                end if
                h=hnext 
            enddo

            write(*,*) 'too many steps in odeint'

            return

            contains

                ! Fifth-order Runge-Kutta step with monitoring of local truncation error to 
                ! ensure accuracy and adjust stepsize. Input are the dependent variable vector 
                ! y(1:n) and its derivative dydx(1:n) at the starting value of the independent 
                ! variable x. Also input are the stepsize to be attempted htry, the required 
                ! accuracy eps, and the vector yscal(1:n) against which the error is scaled. 
                ! On output, y and x are replaced by their new values, hdid is the stepsize 
                ! that was actually accomplished, and hnext is the estimated next stepsize. 
                ! derivs is the user-supplied subroutine that computes the right-hand side derivatives.
                subroutine rkqs_c(y,dydx,x,htry,yscal,hdid,hnext) 
    
                    implicit none
                    double complex, dimension(:), intent(inout) :: y
                    double complex, dimension(:), intent(in) :: dydx, yscal
                    double precision, intent(inout) :: x
                    double precision, intent(in) :: htry
                    double precision, intent(out) :: hdid, hnext

 
                    integer :: i
                    double precision :: errmax,h,htemp,xnew
                    double complex, dimension(size(y)) :: yerr, ytemp
                    double precision, parameter :: SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

                    h=htry

                    do
                        call rkck_c(y, dydx, x, h, ytemp, yerr)
                        errmax = maxval(abs(yerr(:) / yscal(:))) / eps
                        if (errmax <= 1.0) exit
                        htemp = SAFETY * h * (errmax** PSHRNK)
                        h = sign(max(abs(htemp), 0.1d0*abs(h)),h)
                        xnew = x + h 
                        if (xnew == x) then 
                            write(*,*) 'stepsize underflow in integrate_krp, rkqs_c'
                            write(*,*) 'yerr = ', yerr
                            write(*,*) 'yscal = ', yscal
                            write(*,*) 'x =', x, '; xnew = ', xnew, ';h = ', h, '; errmax = ', errmax, '; htemp = ', htemp
                            stop
                        end if

                    end do

                    if(errmax > ERRCON)then 
                        hnext = SAFETY * h * (errmax**PGROW)
                    else
                        hnext = 5.0d0 * h
                    endif 

                    hdid=h 
                    x=x+h
                    y(:)=ytemp(:) 

                end subroutine

                ! Given values for n variables y and their derivatives dydx known at x, 
                ! use the fifth-order Cash-Karp Runge-Kutta method to advance the solution over 
                ! an interval h and return the incremented variables as yout. Also return an 
                ! estimate of the local truncation error in yout using the embedded fourth-order method. 
                ! The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
                subroutine rkck_c(y,dydx,x,h,yout,yerr) 
            
                    implicit none
                    double complex, dimension(:), intent(in) :: y, dydx
                    double precision, intent(in) :: x, h
                    double complex, dimension(:), intent(out) :: yout, yerr

                   !interface
                        !subroutine derivs(x,y, dydx)
                            !implicit none
                            !double precision, intent(in) :: x
                            !double complex, dimension(:), intent(in) :: y
                            !double complex, dimension(:), intent(out) :: dydx
                        !end subroutine
                    !end interface

                    double complex, dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
                    double precision, parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., B32=9./40.,&
                                B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, B53=-70./27.,&
                                B54=35./27.,B61=1631./55296.,B62=175./512., B63=575./13824.,&
                                B64=44275./110592.,B65=253./4096., C1=37./378.,C3=250./621.,&
                                C4=125./594.,C6=512./1771., DC1=C1-2825./27648.,DC3=C3-18575./48384., &
                                DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25

                    ytemp = y + B21 * h * dydx
                
                    call derivs(x+A2*h, ytemp, ak2, j, k1, lp)

                    ytemp = y + h * (B31 * dydx + B32 * ak2) 

                    call derivs(x+A3*h, ytemp, ak3, j, k1, lp) 

                    ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)

                    call derivs(x+A4*h, ytemp, ak4, j, k1, lp)

                    ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)

                    call derivs(x+A5*h, ytemp, ak5, j, k1, lp) 

                    ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)

                    call derivs(x+A6*h, ytemp, ak6, j, k1, lp)

                    yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
                    yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)
                
                end subroutine

        end subroutine


!        subroutine integrate_kr(ystart, nvar, krp_int, x1, x2, derivs)!,rkqs) 

            !implicit none
        
            !integer, intent(in) :: nvar
            !double precision, intent(in) :: x1, x2 ! limits of integration
            !double complex, dimension(nvar), intent(inout) :: ystart
            !integer , intent(in) :: r_ind

            !integer :: i, nstp
            !double precision :: h,hdid,hnext,x,xsav
            !double complex, dimension(size(ystart)) :: dydx, y, yscal

            !interface
                !subroutine derivs(r,y, dydr, c_r, c_kr)
                    !implicit none
                    !double precision, intent(in) :: r, c_r, c_kr
                    !double complex, dimension(:), intent(in) :: y
                    !double complex, dimension(:), intent(out) :: dydr
                !end subroutine
            !end interface

            !x = x1
            !h = sign(h1, x2 - x1) 
            !num_ok  = 0 
            !num_bad = 0

            !y(:) = ystart(:)

            !do nstp = 1, max_steps
                !call derivs(x,y,dydx, c_r, c_kr)

                !yscal(:) = abs(y(:)) + abs(h * dydx(:)) + TINY 

                !if((x + h - x2) * (x+h-x1) > 0.0) h = x2-x

                !call rkqs_c(y,dydx,x,h, yscal,hdid,hnext)!,derivs)

                !if(hdid == h)then
                    !num_ok=num_ok+1
                !else
                    !num_bad=num_bad+1
                !endif
                !if((x-x2)*(x2-x1).ge.0.)then 
                    !ystart(:)=y(:) 
                !return 
                !endif
        
                !if(abs(hnext).lt.hmin) then 
                    !write(*,*) 'stepsize smaller than minimum in odeint'
                !end if
                !h=hnext 
            !enddo

            !write(*,*) 'too many steps in odeint'

            !return

            !contains

                !! Fifth-order Runge-Kutta step with monitoring of local truncation error to 
                !! ensure accuracy and adjust stepsize. Input are the dependent variable vector 
                !! y(1:n) and its derivative dydx(1:n) at the starting value of the independent 
                !! variable x. Also input are the stepsize to be attempted htry, the required 
                !! accuracy eps, and the vector yscal(1:n) against which the error is scaled. 
                !! On output, y and x are replaced by their new values, hdid is the stepsize 
                !! that was actually accomplished, and hnext is the estimated next stepsize. 
                !! derivs is the user-supplied subroutine that computes the right-hand side derivatives.
                !subroutine rkqs_c(y,dydx,x,htry,yscal,hdid,hnext) 
    
                    !implicit none
                    !double complex, dimension(:), intent(inout) :: y
                    !double complex, dimension(:), intent(in) :: dydx, yscal
                    !double precision, intent(inout) :: x
                    !double precision, intent(in) :: htry
                    !double precision, intent(out) :: hdid, hnext

 
                    !integer :: i
                    !double precision :: errmax,h,htemp,xnew
                    !double complex, dimension(size(y)) :: yerr, ytemp
                    !double precision, parameter :: SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

                    !h=htry

                    !do
                        !call rkck_c(y,dydx,x,h,ytemp,yerr)
                        !errmax = maxval(abs(yerr(:) / yscal(:))) / eps
                        !if (errmax <= 1.0) exit
                        !htemp = SAFETY * h * (errmax** PSHRNK)
                        !h = sign(max(abs(htemp), 0.1d0*abs(h)),h)
                        !xnew = x + h 
                        !if (xnew == x) then 
                            !write(*,*) 'stepsize underflow in rkqs_c'
                            !write(*,*) 'yerr = ', yerr
                            !write(*,*) 'yscal = ', yscal
                            !write(*,*) 'x =', x, '; xnew = ', xnew, ';h = ', h, '; errmax = ', errmax, '; htemp = ', htemp
                            !stop
                        !end if

                    !end do

                    !if(errmax > ERRCON)then 
                        !hnext = SAFETY * h * (errmax**PGROW)
                    !else
                        !hnext = 5.0d0 * h
                    !endif 

                    !hdid=h 
                    !x=x+h
                    !y(:)=ytemp(:) 

                !end subroutine

                !! Given values for n variables y and their derivatives dydx known at x, 
                !! use the fifth-order Cash-Karp Runge-Kutta method to advance the solution over 
                !! an interval h and return the incremented variables as yout. Also return an 
                !! estimate of the local truncation error in yout using the embedded fourth-order method. 
                !! The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
                !subroutine rkck_c(y,dydx,x,h,yout,yerr) 
            
                    !implicit none
                    !double complex, dimension(:), intent(in) :: y, dydx
                    !double precision, intent(in) :: x, h
                    !double complex, dimension(:), intent(out) :: yout, yerr

                    !!interface
                        !!subroutine derivs(x,y, dydx)
                            !!implicit none
                            !!double precision, intent(in) :: x
                            !!double complex, dimension(:), intent(in) :: y
                            !!double complex, dimension(:), intent(out) :: dydx
                        !!end subroutine
                    !!end interface

                    !double complex, dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
                    !double precision, parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., B32=9./40.,&
                                !B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, B53=-70./27.,&
                                !B54=35./27.,B61=1631./55296.,B62=175./512., B63=575./13824.,&
                                !B64=44275./110592.,B65=253./4096., C1=37./378.,C3=250./621.,&
                                !C4=125./594.,C6=512./1771., DC1=C1-2825./27648.,DC3=C3-18575./48384., &
                                !DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25

                    !ytemp = y + B21 * h * dydx
                
                    !call derivs(x+A2*h, ytemp, ak2, c_r, c_kr)

                    !ytemp = y + h * (B31 * dydx + B32 * ak2) 

                    !call derivs(x+A3*h, ytemp, ak3, c_r, c_kr) 

                    !ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)

                    !call derivs(x+A4*h, ytemp, ak4, c_r, c_kr)

                    !ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)

                    !call derivs(x+A5*h, ytemp, ak5, c_r, c_kr) 

                    !ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)

                    !call derivs(x+A6*h, ytemp, ak6, c_r, c_kr)

                    !yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
                    !yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)
                
                !end subroutine

        !end subroutine
        

        subroutine integrate_kr(ystart, nvar, kr1, kr2, int_kr_r, j, l, derivs)!,rkqs) 

            use grid, only: k_space_dim, npoib, kr, krp
            implicit none
        
            integer, intent(in) :: nvar
            integer, intent(in) :: j     ! index for r profile
            integer, intent(in) :: l     ! index for spline function
            double complex, dimension(k_space_dim, npoib), intent(in) :: int_kr_r
            double precision, intent(in) :: kr1, kr2 ! limits of integration
            double complex, dimension(nvar), intent(inout) :: ystart

            integer :: i, nstp
            double precision :: h,hdid,hnext,x,xsav
            double complex, dimension(size(ystart)) :: dydx, y, yscal

            interface
                subroutine derivs(r,y, dydr, int_kr_r, j, l)
                    use grid, only: k_space_dim, npoib
                    implicit none
                    integer, intent(in) :: j
                    integer, intent(in) :: l
                    double precision, intent(in) :: r
                    double complex, dimension(k_space_dim, npoib), intent(in) :: int_kr_r
                    double complex, dimension(:), intent(in) :: y
                    double complex, dimension(:), intent(out) :: dydr
                end subroutine
            end interface

            x = kr1
            h = sign(h1, kr2 - kr1) 
            num_ok  = 0 
            num_bad = 0

            y(:) = ystart(:)

            do nstp = 1, max_steps
                call derivs(x,y,dydx, int_kr_r, j, l)

                yscal(:) = abs(y(:)) + abs(h * dydx(:)) + TINY 

                if((x + h - kr2) * (x+h-kr1) > 0.0) h = kr2-x

                call rkqs_c(y,dydx,x,h, yscal,hdid,hnext)!,derivs)

                if(hdid == h)then
                    num_ok=num_ok+1
                else
                    num_bad=num_bad+1
                endif
                if((x-kr2)*(kr2-kr1).ge.0.)then 
                    ystart(:)=y(:) 
                return 
                endif
        
                if(abs(hnext).lt.hmin) then 
                    write(*,*) 'stepsize smaller than minimum in odeint'
                end if
                h=hnext 
            enddo

            write(*,*) 'too many steps in odeint'

            return

            contains

                ! Fifth-order Runge-Kutta step with monitoring of local truncation error to 
                ! ensure accuracy and adjust stepsize. Input are the dependent variable vector 
                ! y(1:n) and its derivative dydx(1:n) at the starting value of the independent 
                ! variable x. Also input are the stepsize to be attempted htry, the required 
                ! accuracy eps, and the vector yscal(1:n) against which the error is scaled. 
                ! On output, y and x are replaced by their new values, hdid is the stepsize 
                ! that was actually accomplished, and hnext is the estimated next stepsize. 
                ! derivs is the user-supplied subroutine that computes the right-hand side derivatives.
                subroutine rkqs_c(y,dydx,x,htry,yscal,hdid,hnext) 
    
                    implicit none
                    double complex, dimension(:), intent(inout) :: y
                    double complex, dimension(:), intent(in) :: dydx, yscal
                    double precision, intent(inout) :: x
                    double precision, intent(in) :: htry
                    double precision, intent(out) :: hdid, hnext

 
                    integer :: i
                    double precision :: errmax,h,htemp,xnew
                    double complex, dimension(size(y)) :: yerr, ytemp
                    double precision, parameter :: SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

                    h=htry

                    do
                        call rkck_c(y,dydx,x,h,ytemp,yerr)
                        errmax = maxval(abs(yerr(:) / yscal(:))) / eps
                        if (errmax <= 1.0) exit
                        htemp = SAFETY * h * (errmax** PSHRNK)
                        h = sign(max(abs(htemp), 0.1d0*abs(h)),h)
                        xnew = x + h 
                        if (xnew == x) then 
                            write(*,*) 'stepsize underflow in integrate_kr, rkqs_c'
                            write(*,*) 'yerr = ', yerr
                            write(*,*) 'yscal = ', yscal
                            write(*,*) 'x =', x, '; xnew = ', xnew, ';h = ', h, '; errmax = ', errmax, '; htemp = ', htemp
                            stop
                        end if

                    end do

                    if(errmax > ERRCON)then 
                        hnext = SAFETY * h * (errmax**PGROW)
                    else
                        hnext = 5.0d0 * h
                    endif 

                    hdid=h 
                    x=x+h
                    y(:)=ytemp(:) 

                end subroutine

                ! Given values for n variables y and their derivatives dydx known at x, 
                ! use the fifth-order Cash-Karp Runge-Kutta method to advance the solution over 
                ! an interval h and return the incremented variables as yout. Also return an 
                ! estimate of the local truncation error in yout using the embedded fourth-order method. 
                ! The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
                subroutine rkck_c(y,dydx,x,h,yout,yerr) 
            
                    implicit none
                    double complex, dimension(:), intent(in) :: y, dydx
                    double precision, intent(in) :: x, h
                    double complex, dimension(:), intent(out) :: yout, yerr

                   !interface
                        !subroutine derivs(x,y, dydx)
                            !implicit none
                            !double precision, intent(in) :: x
                            !double complex, dimension(:), intent(in) :: y
                            !double complex, dimension(:), intent(out) :: dydx
                        !end subroutine
                    !end interface

                    double complex, dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
                    double precision, parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., B32=9./40.,&
                                B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, B53=-70./27.,&
                                B54=35./27.,B61=1631./55296.,B62=175./512., B63=575./13824.,&
                                B64=44275./110592.,B65=253./4096., C1=37./378.,C3=250./621.,&
                                C4=125./594.,C6=512./1771., DC1=C1-2825./27648.,DC3=C3-18575./48384., &
                                DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25

                    ytemp = y + B21 * h * dydx
                
                    call derivs(x+A2*h, ytemp, ak2, int_kr_r, j, l)

                    ytemp = y + h * (B31 * dydx + B32 * ak2) 

                    call derivs(x+A3*h, ytemp, ak3, int_kr_r, j, l) 

                    ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)

                    call derivs(x+A4*h, ytemp, ak4, int_kr_r, j, l)

                    ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)

                    call derivs(x+A5*h, ytemp, ak5, int_kr_r, j, l) 

                    ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)

                    call derivs(x+A6*h, ytemp, ak6, int_kr_r, j, l)

                    yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
                    yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)
                
                end subroutine

        end subroutine

        subroutine integrate_rg(ystart, nvar, r1, r2, int_rg, derivs)!,rkqs) 

            use grid, only: npoib
            implicit none
        
            integer, intent(in) :: nvar
            !integer, intent(in) :: j     ! index for r profile
            !integer, intent(in) :: l     ! index for spline function
            double complex, dimension(npoib), intent(in) :: int_rg
            double precision, intent(in) :: r1, r2 ! limits of integration
            double complex, dimension(nvar), intent(inout) :: ystart

            integer :: i, nstp
            double precision :: h,hdid,hnext,x,xsav
            double complex, dimension(size(ystart)) :: dydx, y, yscal

            interface
                subroutine derivs(r,y, dydr, int_rg)
                    use grid, only: npoib
                    implicit none
                    double precision, intent(in) :: r
                    double complex, dimension(npoib), intent(in) :: int_rg
                    double complex, dimension(:), intent(in) :: y
                    double complex, dimension(:), intent(out) :: dydr
                end subroutine
            end interface

            x = r1
            h = sign(h1, r2 - r1) 
            num_ok  = 0 
            num_bad = 0

            y(:) = ystart(:)

            do nstp = 1, max_steps
                call derivs(x,y,dydx, int_rg)

                yscal(:) = abs(y(:)) + abs(h * dydx(:)) + TINY 

                if((x + h - r2) * (x+h-r1) > 0.0) h = r2-x

                call rkqs_c(y,dydx,x,h, yscal,hdid,hnext)!,derivs)

                if(hdid == h)then
                    num_ok=num_ok+1
                else
                    num_bad=num_bad+1
                endif
                if((x-r2)*(r2-r1).ge.0.)then 
                    ystart(:)=y(:) 
                return 
                endif
        
                if(abs(hnext).lt.hmin) then 
                    write(*,*) 'stepsize smaller than minimum in odeint'
                end if
                h=hnext 
            enddo

            write(*,*) 'too many steps in odeint'

            return

            contains

                ! Fifth-order Runge-Kutta step with monitoring of local truncation error to 
                ! ensure accuracy and adjust stepsize. Input are the dependent variable vector 
                ! y(1:n) and its derivative dydx(1:n) at the starting value of the independent 
                ! variable x. Also input are the stepsize to be attempted htry, the required 
                ! accuracy eps, and the vector yscal(1:n) against which the error is scaled. 
                ! On output, y and x are replaced by their new values, hdid is the stepsize 
                ! that was actually accomplished, and hnext is the estimated next stepsize. 
                ! derivs is the user-supplied subroutine that computes the right-hand side derivatives.
                subroutine rkqs_c(y,dydx,x,htry,yscal,hdid,hnext) 
    
                    implicit none
                    double complex, dimension(:), intent(inout) :: y
                    double complex, dimension(:), intent(in) :: dydx, yscal
                    double precision, intent(inout) :: x
                    double precision, intent(in) :: htry
                    double precision, intent(out) :: hdid, hnext

 
                    integer :: i
                    double precision :: errmax,h,htemp,xnew
                    double complex, dimension(size(y)) :: yerr, ytemp
                    double precision, parameter :: SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

                    h=htry

                    do
                        call rkck_c(y,dydx,x,h,ytemp,yerr)
                        errmax = maxval(abs(yerr(:) / yscal(:))) / eps
                        if (errmax <= 1.0) exit
                        htemp = SAFETY * h * (errmax** PSHRNK)
                        h = sign(max(abs(htemp), 0.1d0*abs(h)),h)
                        xnew = x + h 
                        if (xnew == x) then 
                            write(*,*) 'stepsize underflow in integrate_rg, rkqs_c'
                            write(*,*) 'yerr = ', yerr
                            write(*,*) 'yscal = ', yscal
                            write(*,*) 'x =', x, '; xnew = ', xnew, ';h = ', h, '; errmax = ', errmax, '; htemp = ', htemp
                            stop
                        end if

                    end do

                    if(errmax > ERRCON)then 
                        hnext = SAFETY * h * (errmax**PGROW)
                    else
                        hnext = 5.0d0 * h
                    endif 

                    hdid=h 
                    x=x+h
                    y(:)=ytemp(:) 

                end subroutine

                ! Given values for n variables y and their derivatives dydx known at x, 
                ! use the fifth-order Cash-Karp Runge-Kutta method to advance the solution over 
                ! an interval h and return the incremented variables as yout. Also return an 
                ! estimate of the local truncation error in yout using the embedded fourth-order method. 
                ! The user supplies the subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
                subroutine rkck_c(y,dydx,x,h,yout,yerr) 
            
                    implicit none
                    double complex, dimension(:), intent(in) :: y, dydx
                    double precision, intent(in) :: x, h
                    double complex, dimension(:), intent(out) :: yout, yerr

                   !interface
                        !subroutine derivs(x,y, dydx)
                            !implicit none
                            !double precision, intent(in) :: x
                            !double complex, dimension(:), intent(in) :: y
                            !double complex, dimension(:), intent(out) :: dydx
                        !end subroutine
                    !end interface

                    double complex, dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
                    double precision, parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., B32=9./40.,&
                                B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, B53=-70./27.,&
                                B54=35./27.,B61=1631./55296.,B62=175./512., B63=575./13824.,&
                                B64=44275./110592.,B65=253./4096., C1=37./378.,C3=250./621.,&
                                C4=125./594.,C6=512./1771., DC1=C1-2825./27648.,DC3=C3-18575./48384., &
                                DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25

                    ytemp = y + B21 * h * dydx
                
                    call derivs(x+A2*h, ytemp, ak2, int_rg)

                    ytemp = y + h * (B31 * dydx + B32 * ak2) 

                    call derivs(x+A3*h, ytemp, ak3, int_rg) 

                    ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)

                    call derivs(x+A4*h, ytemp, ak4, int_rg)

                    ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)

                    call derivs(x+A5*h, ytemp, ak5, int_rg) 

                    ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)

                    call derivs(x+A6*h, ytemp, ak6, int_rg)

                    yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
                    yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)
                
                end subroutine

        end subroutine





end module