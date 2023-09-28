module adaptive_int

    implicit none



    contains

    subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs) 

        implicit none

        integer, intent(out) :: nbad, nok
        integer, intent(in) :: nvar
        integer :: KMAXX, MAXSTP, NMAX
        double precision, intent(in) :: eps, h1, hmin, x1, x2, ystart(nvar)
        double precision :: TINY
        external derivs,rkqs
        parameter :: MAXSTP=10000, NMAX=50, KMAXX=200, TINY=1.0d-30

        integer :: i,kmax,kount,nstp
        double precision :: dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX), &
                            yp(NMAX,KMAXX),yscal(NMAX) 
        !COMMON /path/ kmax,kount,dxsav,xp,yp
        h=sign(h1,x2-x1) 
        nok=0 
        nbad=0
        kount=0
        do i=1,nvar 
            y(i)=ystart(i)
        enddo

        if (kmax.gt.0) xsav=x-2.*dxsav

        do nstp=1,MAXSTP
            call derivs(x,y,dydx) 
            do i=1,nvar
                yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY 
            enddo

            if(kmax.gt.0)then 
                if(abs(x-xsav).gt.abs(dxsav)) then
                    if(kount.lt.kmax-1)then 
                        kount=kount+1
                        xp(kount)=x 
                        do i=1,nvar
                            yp(i,kount)=y(i)
                        enddo
                        xsav=x 
                    endif
                endif 
            endif
            if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
            call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs) 
            if(hdid.eq.h)then
                nok=nok+1
            else
                nbad=nbad+1
            endif
            if((x-x2)*(x2-x1).ge.0.)then 
                do i=1,nvar
                    ystart(i)=y(i) 
                enddo
                if(kmax.ne.0)then 
                    kount=kount+1
                    xp(kount)=x 
                    do i=1,nvar
                        yp(i,kount)=y(i) 
                    enddo
                endif
                return 
            endif
        
            if(abs(hnext).lt.hmin) pause ’stepsize smaller than minimum in odeint’
            h=hnext 
        enddo

        pause ’too many steps in odeint’ 

        return

    end subroutine

end module