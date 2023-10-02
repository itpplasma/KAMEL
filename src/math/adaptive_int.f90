! adaptive Runge-Kutta integrator
! used for calculating integrals
module adaptive_int

    implicit none

    logical :: save_steps = .true.
    integer, parameter :: MAXSTP=10000!, NMAX=50, KMAXX=200
    integer :: kmax = 0
    integer :: kount
    double precision :: dxsav = 0.01d0
    double precision, dimension(:), pointer :: xp
    double precision, dimension(:,:), pointer :: yp
    double complex, dimension(:,:), pointer :: yp_c

    contains

    subroutine odeint(ystart, x1, x2, eps, h1, hmin, nok, nbad, derivs)!,rkqs) 

        implicit none

        integer, intent(out) :: nbad, nok
        double precision, intent(in) :: eps, h1, hmin, x1, x2
        double precision, dimension(:), intent(inout) :: ystart
        !external derivs!,rkqs
        double precision, parameter :: TINY=1.0d-30

        integer :: i, nstp
        double precision :: h,hdid,hnext,x,xsav
        double precision, dimension(size(ystart)) :: dydx, y, yscal

        interface
            subroutine derivs(x,y, dydx)
                implicit none
                double precision, intent(in) :: x
                double precision, dimension(:), intent(in) :: y
                double precision, dimension(:), intent(out) :: dydx
            end subroutine
        end interface

        x = x1
        h = sign(h1, x2 - x1) 
        nok  = 0 
        nbad = 0
        kount= 0

        y(:)=ystart(:)

        nullify(xp,yp)

        if (save_steps) then
            xsav = x-2.0d0 * dxsav
            allocate(xp(256))
            allocate(yp(size(ystart), size(xp)))
        end if

        do nstp=1,MAXSTP
            call derivs(x,y,dydx) 
            yscal(:) = abs(y(:)) + abs(h * dydx(:)) + TINY 

            if(save_steps .and. (abs(x-xsav) > abs(dxsav))) call save_a_step

            if((x + h - x2) * (x+h-x1) > 0.0) h = x2-x

            call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)

            if(hdid == h)then
                nok=nok+1
            else
                nbad=nbad+1
            endif
            if((x-x2)*(x2-x1).ge.0.)then 

                ystart(:)=y(:) 

                if(kmax.ne.0)then 
                    kount=kount+1
                    xp(kount)=x 
                    yp(:,kount)=y(:) 
                endif
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

            subroutine save_a_step
                kount = kount+1
                if(kount > size(xp)) then
                    xp=> reallocate_dv(xp, 2*size(xp))
                    yp=> reallocate_dm(yp, size(yp,1), size(xp))
                end if
                xp(kount) = x
                yp(:, kount) = y(:)
                xsav = x
            end subroutine

            function reallocate_dv(p,n)
                double precision, dimension(:), pointer :: p, reallocate_dv
                integer, intent(in) :: n
                integer :: nold, ierr
                allocate(reallocate_dv(n), stat=ierr)
                if (ierr /= 0) write(*,*) ' not possible to reallocate'
                if (.not. associated(p)) return

                nold = size(p)
                reallocate_dv(1:min(nold,n))=p(1:min(nold,n))
                deallocate(p)
            end function

            function reallocate_dm(p,n,m)
                double precision, dimension(:,:), pointer :: p, reallocate_dm
                integer, intent(in) :: n, m
                integer :: nold, mold, ierr
                allocate(reallocate_dm(n,m), stat=ierr)
                if (ierr /= 0) write(*,*) ' not possible to reallocate'
                if (.not. associated(p)) return

                nold = size(p,1)
                mold = size(p,2)
                reallocate_dm(1:min(nold,n), 1:min(mold,m)) = p(1:min(nold,n), 1:min(mold,m))
                deallocate(p)
            end function


    end subroutine

    ! Fifth-order Runge-Kutta step with monitoring of local truncation error to 
    ! ensure accuracy and adjust stepsize. Input are the dependent variable vector 
    ! y(1:n) and its derivative dydx(1:n) at the starting value of the independent 
    ! variable x. Also input are the stepsize to be attempted htry, the required 
    ! accuracy eps, and the vector yscal(1:n) against which the error is scaled. 
    ! On output, y and x are replaced by their new values, hdid is the stepsize 
    ! that was actually accomplished, and hnext is the estimated next stepsize. 
    ! derivs is the user-supplied subroutine that computes the right-hand side derivatives.
    subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs) 
    
        implicit none
        double precision, dimension(:), intent(inout) :: y
        double precision, dimension(:), intent(in) :: dydx, yscal
        double precision, intent(inout) :: x
        double precision, intent(in) :: htry, eps
        double precision, intent(out) :: hdid, hnext

        interface
            subroutine derivs(x,y, dydx)
                implicit none
                double precision, intent(in) :: x
                double precision, dimension(:), intent(in) :: y
                double precision, dimension(:), intent(out) :: dydx
            end subroutine
        end interface
 
        integer :: i
        double precision :: errmax,h,htemp,xnew
        double precision, dimension(size(y)) :: yerr, ytemp
        double precision, parameter :: SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

       
        h=htry

        do
            call rkck(y,dydx,x,h,ytemp,yerr,derivs)
            errmax = maxval(abs(yerr(:) / yscal(:))) / eps
            if (errmax <= 1.0) exit
            htemp = SAFETY * h * (errmax** PSHRNK)
            h = sign(max(abs(htemp), 0.1d0*abs(h)),h)
            xnew = x+ h 
            if (xnew == x) then
                write(*,*) 'stepsize underflow in rkqs'
                !stop
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
    subroutine rkck(y,dydx,x,h,yout,yerr,derivs) 
    
        implicit none
        double precision, dimension(:), intent(in) :: y, dydx
        double precision, intent(in) :: x, h
        double precision, dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs(x,y, dydx)
                implicit none
                double precision, intent(in) :: x
                double precision, dimension(:), intent(in) :: y
                double precision, dimension(:), intent(out) :: dydx
            end subroutine
        end interface

        double precision, dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
        double precision, parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., B32=9./40.,&
                    B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, B53=-70./27.,&
                    B54=35./27.,B61=1631./55296.,B62=175./512., B63=575./13824.,&
                    B64=44275./110592.,B65=253./4096., C1=37./378.,C3=250./621.,&
                    C4=125./594.,C6=512./1771., DC1=C1-2825./27648.,DC3=C3-18575./48384., &
                    DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25

        ytemp = y + B21 * h * dydx
        
        call derivs(x+A2*h, ytemp, ak2)  

        ytemp = y + h * (B31 * dydx + B32 * ak2) 

        call derivs(x+A3*h, ytemp, ak3) 

        ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)

        call derivs(x+A4*h, ytemp, ak4)

        ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)

        call derivs(x+A5*h, ytemp, ak5) 

        ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)

        call derivs(x+A6*h, ytemp, ak6)

        yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
        yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)
        
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! complex 
    subroutine odeint_c(ystart, x1, x2, eps, h1, hmin, nok, nbad, derivs)!,rkqs) 

        implicit none

        integer, intent(out) :: nbad, nok
        double precision, intent(in) :: eps, h1, hmin, x1, x2
        double complex, dimension(:), intent(inout) :: ystart
        !external derivs!,rkqs
        double precision, parameter :: TINY=1.0d-30

        integer :: i, nstp
        double precision :: h,hdid,hnext,x,xsav
        double complex, dimension(size(ystart)) :: dydx, y, yscal

        interface
            subroutine derivs(x,y, dydx)
                implicit none
                double precision, intent(in) :: x
                double complex, dimension(:), intent(in) :: y
                double complex, dimension(:), intent(out) :: dydx
            end subroutine
        end interface

        x = x1
        h = sign(h1, x2 - x1) 
        nok  = 0 
        nbad = 0
        kount= 0

        y(:) = ystart(:)

        nullify(xp,yp_c)

        if (save_steps) then
            xsav = x-2.0d0 * dxsav
            allocate(xp(256))
            allocate(yp_c(size(ystart), size(xp)))
        end if

        do nstp = 1, MAXSTP
            call derivs(x,y,dydx) 
            yscal(:) = abs(y(:)) + abs(h * dydx(:)) + TINY 

            if(save_steps .and. (abs(x-xsav) > abs(dxsav))) call save_a_step

            if((x + h - x2) * (x+h-x1) > 0.0) h = x2-x

            call rkqs_c(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)

            if(hdid == h)then
                nok=nok+1
            else
                nbad=nbad+1
            endif
            if((x-x2)*(x2-x1).ge.0.)then 

                ystart(:)=y(:) 

                if(kmax.ne.0)then 
                    kount=kount+1
                    xp(kount)=x 
                    yp_c(:,kount)=y(:) 
                endif
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

            subroutine save_a_step
                kount = kount+1
                if(kount > size(xp)) then
                    xp=> reallocate_dv(xp, 2*size(xp))
                    yp_c=> reallocate_cm(yp_c, size(yp_c,1), size(xp))
                end if
                xp(kount) = x
                yp_c(:, kount) = y(:)
                xsav = x
            end subroutine

            function reallocate_dv(p,n)
                double precision, dimension(:), pointer :: p, reallocate_dv
                integer, intent(in) :: n
                integer :: nold, ierr
                allocate(reallocate_dv(n), stat=ierr)
                if (ierr /= 0) write(*,*) ' not possible to reallocate'
                if (.not. associated(p)) return

                nold = size(p)
                reallocate_dv(1:min(nold,n))=p(1:min(nold,n))
                deallocate(p)
            end function

            function reallocate_cm(p,n,m)
                double complex, dimension(:,:), pointer :: p, reallocate_cm
                integer, intent(in) :: n, m
                integer :: nold, mold, ierr
                allocate(reallocate_cm(n,m), stat=ierr)
                if (ierr /= 0) write(*,*) ' not possible to reallocate'
                if (.not. associated(p)) return

                nold = size(p,1)
                mold = size(p,2)
                reallocate_cm(1:min(nold,n), 1:min(mold,m)) = p(1:min(nold,n), 1:min(mold,m))
                deallocate(p)
            end function


    end subroutine

    ! Fifth-order Runge-Kutta step with monitoring of local truncation error to 
    ! ensure accuracy and adjust stepsize. Input are the dependent variable vector 
    ! y(1:n) and its derivative dydx(1:n) at the starting value of the independent 
    ! variable x. Also input are the stepsize to be attempted htry, the required 
    ! accuracy eps, and the vector yscal(1:n) against which the error is scaled. 
    ! On output, y and x are replaced by their new values, hdid is the stepsize 
    ! that was actually accomplished, and hnext is the estimated next stepsize. 
    ! derivs is the user-supplied subroutine that computes the right-hand side derivatives.
    subroutine rkqs_c(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs) 
    
        implicit none
        double complex, dimension(:), intent(inout) :: y
        double complex, dimension(:), intent(in) :: dydx, yscal
        double precision, intent(inout) :: x
        double precision, intent(in) :: htry, eps
        double precision, intent(out) :: hdid, hnext

        interface
            subroutine derivs(x,y, dydx)
                implicit none
                double precision, intent(in) :: x
                double complex, dimension(:), intent(in) :: y
                double complex, dimension(:), intent(out) :: dydx
            end subroutine
        end interface
 
        integer :: i
        double precision :: errmax,h,htemp,xnew
        double complex, dimension(size(y)) :: yerr, ytemp
        double precision, parameter :: SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

        h=htry

        do
            call rkck_c(y,dydx,x,h,ytemp,yerr,derivs)
            errmax = maxval(abs(yerr(:) / yscal(:))) / eps
            if (errmax <= 1.0) exit
            htemp = SAFETY * h * (errmax** PSHRNK)
            h = sign(max(abs(htemp), 0.1d0*abs(h)),h)
            xnew = x + h 
            !write(*,*) 'x =', x, '; xnew = ', xnew, ';h = ', h, '; errmax = ', errmax, '; htemp = ', htemp
            if (xnew == x) then 
                write(*,*) 'stepsize underflow in rkqs_c'
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
    subroutine rkck_c(y,dydx,x,h,yout,yerr,derivs) 
    
        implicit none
        double complex, dimension(:), intent(in) :: y, dydx
        double precision, intent(in) :: x, h
        double complex, dimension(:), intent(out) :: yout, yerr

        interface
            subroutine derivs(x,y, dydx)
                implicit none
                double precision, intent(in) :: x
                double complex, dimension(:), intent(in) :: y
                double complex, dimension(:), intent(out) :: dydx
            end subroutine
        end interface

        double complex, dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
        double precision, parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., B32=9./40.,&
                    B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, B53=-70./27.,&
                    B54=35./27.,B61=1631./55296.,B62=175./512., B63=575./13824.,&
                    B64=44275./110592.,B65=253./4096., C1=37./378.,C3=250./621.,&
                    C4=125./594.,C6=512./1771., DC1=C1-2825./27648.,DC3=C3-18575./48384., &
                    DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25

        ytemp = y + B21 * h * dydx
        
        call derivs(x+A2*h, ytemp, ak2)  

        ytemp = y + h * (B31 * dydx + B32 * ak2) 

        call derivs(x+A3*h, ytemp, ak3) 

        ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)

        call derivs(x+A4*h, ytemp, ak4)

        ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)

        call derivs(x+A5*h, ytemp, ak5) 

        ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)

        call derivs(x+A6*h, ytemp, ak6)

        yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
        yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)
        
    end subroutine


end module