module amn_of_r_mod

    use QLBalance_kinds, only: dp

    implicit none

    logical :: prop=.true.
    integer :: iunit=77, nder=0, nlagr=4
    integer :: ntor,mpol,nlabel,nsqpsi
    real(dp) :: ampl,flabel_min,flabel_max
    real(dp), dimension(:),     allocatable :: rsmall,qsaf,psisurf,flabel
    real(dp), dimension(:,:),   allocatable :: coef
    complex(dp),   dimension(:,:,:), allocatable :: apsimn,athetmn

contains

    subroutine binsrc_amn(p,nmin,nmax,xi,i)

        ! Finds the index  i  of the array of increasing numbers   p  with dimension  n 
        ! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.

        use QLBalance_kinds, only: dp

        implicit none

        integer :: n,nmin,nmax,i,imin,imax,k
        real(dp) :: xi
        real(dp), dimension(nmin:nmax) :: p
      
        imin=nmin
        imax=nmax
        n=nmax-nmin
      
        do k=1,n
            i=(imax-imin)/2+imin
            if(p(i).gt.xi) then
                imax=i
            else
                imin=i
            endif
            if(imax.eq.imin+1) exit
        enddo
      
        i=imax
      
    end subroutine binsrc_amn


    SUBROUTINE plag_coeff_amn(npoi,nder,x,xp,coef)
        !
        ! npoi - number of points (determines the order of Lagrange
        ! polynomial
        ! which is equal npoi-1)
        ! nder - number of derivatives computed 0 - function only, 1 - first
        ! derivative
        ! x - actual point where function and derivatives are evaluated
        ! xp(npoi) - array of points where function is known
        ! coef(0:nder,npoi) - weights for computation of function and
        ! derivatives,
        ! f=sum(fun(1:npoi)*coef(0,1:npoi) gives the function value
        ! df=sum(fun(1:npoi)*coef(1,1:npoi) gives the derivative value value
        !
        !

        use QLBalance_kinds, only: dp

        implicit none

        INTEGER, INTENT(in) :: npoi,nder
        real(dp), INTENT(in)                          :: x
        real(dp), DIMENSION(npoi), INTENT(in)         :: xp
        real(dp), DIMENSION(0:nder,npoi), INTENT(out) :: coef
        real(dp), DIMENSION(:), ALLOCATABLE           :: dummy
        
        INTEGER :: i,k,j
        double precision :: fac
        
        DO i=1,npoi
            coef(0,i)=1.d0
            DO k=1,npoi
                IF(k.EQ.i) CYCLE
                coef(0,i)=coef(0,i)*(x-xp(k))/(xp(i)-xp(k))
            ENDDO
        ENDDO
        
        IF(nder.EQ.0) RETURN
        
        ALLOCATE(dummy(npoi))
        
        DO i=1,npoi
            dummy=1.d0
            dummy(i)=0.d0
            DO k=1,npoi
                IF(k.EQ.i) CYCLE
                fac=(x-xp(k))/(xp(i)-xp(k))
                DO j=1,npoi
                    IF(j.EQ.k) THEN
                        dummy(j)=dummy(j)/(xp(i)-xp(k))
                    ELSE
                        dummy(j)=dummy(j)*fac
                    ENDIF
                ENDDO
            ENDDO
            coef(1,i)=SUM(dummy)
        ENDDO
        
        DEALLOCATE(dummy)
        
    END SUBROUTINE plag_coeff_amn

end module amn_of_r_mod

subroutine amn_of_r(m,n,r_in,amn_psi,amn_theta,ierr)
    !
    ! Computes Fourier amplitudes of vector potential components $A^\psi$ and
    ! $A_\vartheta$ as functions of effective radius. All quantitites are in 
    ! Gaussian units
    !
    ! Input  arguments:
    !           Formal: m         - poloidal wave number
    !                   n         - toroidal wave number
    !                   r         - effective radius
    ! Output arguments:
    !           Formal: amn_psi   - $(A^\psi)_{m,n}$
    !                   amn_theta - $(A_\vartheta)_{m,n}$
    !                   ierr      - error code: 0 - normal work
    !                                           1 - m out of range
    !                                           2 - n out of range
    !                                           3 - r out of range
    !
    use amn_of_r_mod
    use grid_mod, only : rb_cut_in, rb_cut_out,re_cut_out
    use QLBalance_kinds, only: dp
    
    implicit none
    
    character*1024 :: fluxdatapath
    complex(dp), parameter :: cutradius=0.97d0
    integer :: m,n,ierr,ibeg,iend,k,i,nsqpsi_max,nsqpsi_rq
    real(dp) :: rtor,sqpsimin,sqpsimax,sqpsi,hsqpsi_rq,sqpsimin_rq,w
    real(dp) :: r_in,r,q,dummy,weight
    real(dp), save :: rcut,rmax
    complex(dp) :: amn_psi,amn_theta
    
    r=r_in
    
    if(prop) then
        prop=.false.
        allocate(coef(0:nder,nlagr))
    
        open(iunit, file='field_divB0.inp')
        read(iunit,*)
        read(iunit,*)
        read(iunit,*) ampl         ! amplitude of perturbation, a.u.
        read(iunit,*)
        read(iunit,*)
        read(iunit,*)
        read(iunit,*) 
        read(iunit,*)
        read(iunit,*)
        read(iunit,*) fluxdatapath ! directory with data in flux coord.
        close(iunit)
    
        ! Fourier ampitudes of the original field:
    
        open(iunit,form='unformatted',file=trim(fluxdatapath)//'/amn.dat')
        read (iunit) ntor,mpol,nlabel,flabel_min,flabel_max
        allocate(apsimn(-mpol:mpol,ntor,nlabel))
        allocate(athetmn(-mpol:mpol,ntor,nlabel))
        read (iunit) apsimn,athetmn
        close(iunit)
    
        ! Flux surface label for Amn s is equal to r*q
        allocate(flabel(nlabel))
        do i=1,nlabel
            flabel(i)=flabel_min+(flabel_max-flabel_min)*dfloat(i-1)/dfloat(nlabel-1)
        enddo
    
        nsqpsi_max=10000000
        nsqpsi=0
        open(iunit,file=trim(fluxdatapath)//'/equil_r_q_psi.dat')
        read (iunit,*)
        read (iunit,*)
        read (iunit,*)
        do i=1,nsqpsi_max
            read (iunit,*,end=1) dummy
            nsqpsi=nsqpsi+1
        enddo
    1   continue
        close(iunit)
        allocate(rsmall(nsqpsi),qsaf(nsqpsi),psisurf(nsqpsi))
    
        open(iunit,file=trim(fluxdatapath)//'/equil_r_q_psi.dat')
        read (iunit,*)
        read (iunit,*)
        read (iunit,*)
        do i=1,nsqpsi
            read (iunit,*) rsmall(i),qsaf(i),psisurf(i)
        enddo
        close(iunit)
        rcut=cutradius*rsmall(nsqpsi)
        rmax=rsmall(nsqpsi)
    endif
    
    ierr=0
    if(m.lt.-mpol.or.m.gt.mpol) then
        ierr=1
        return
    elseif(n.lt.1.or.n.gt.ntor) then
        ierr=2
        return
    elseif(r.gt.rsmall(nsqpsi)) then
        !    ierr=3
        !    return
        r=rsmall(nsqpsi)
    endif
    
    call binsrc_amn(rsmall,1,nsqpsi,r,i)
    
    ibeg=max(1,i-nlagr/2)
    iend=ibeg+nlagr-1
    if(iend.gt.nsqpsi) then
        iend=nsqpsi
        ibeg=iend-nlagr+1
    endif
    
    call plag_coeff_amn(nlagr,nder,r,rsmall(ibeg:iend),coef)
    
    q=sum(coef(0,:)*qsaf(ibeg:iend))
    dummy=r*q
    
    if(dummy.gt.flabel(nlabel)) then
        ierr=3
        return
    endif
    
    call binsrc_amn(flabel,1,nlabel,dummy,i)
    
    ibeg=max(1,i-nlagr/2)
    iend=ibeg+nlagr-1
    if(iend.gt.nlabel) then
      iend=nlabel
      ibeg=iend-nlagr+1
    endif
    
    call plag_coeff_amn(nlagr,nder,dummy,flabel(ibeg:iend),coef)
    
    amn_psi=sum(coef(0,:)*apsimn(m,n,ibeg:iend))*ampl
    amn_theta=sum(coef(0,:)*athetmn(m,n,ibeg:iend))*ampl
    
    call localizer(1.d0,rb_cut_out,re_cut_out,r,weight)
    amn_psi=amn_psi*weight
    amn_theta=amn_theta*weight
    
    call localizer(-1.d0,rb_cut_in,rb_cut_in,r,weight)
    amn_psi=amn_psi*weight
    amn_theta=amn_theta*weight
    
end subroutine amn_of_r

subroutine localizer(sig,x1,x2,x,weight)

  use QLBalance_kinds, only: dp

  implicit none

  real(dp) :: sig,x1,x2,x,t,weight

  if(sig.gt.0.d0) then
      ! from 1 to 0:
      t=(x-x1)/(x2-x1)
  else
      ! from 0 to 1:
      t=(x2-x)/(x2-x1)
  endif

  if(t.le.0.d0) then
      weight=1.d0
  elseif(t.ge.1.d0) then
      weight=0.d0
  else
      weight=exp(-6.283185307179586d0/(1.d0-t)*exp(-1.414213562373095d0/t))
  endif

end subroutine localizer
