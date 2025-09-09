!
  subroutine gengrid
!
! Generates the grid for two types of boundary conditions at the outer
! boundary: iboutype=1 - fixed parameters, iboutype=2 - fixed fluxes
! At the inner boundary fixed fluxes = 0 are always assumed

    use grid_mod
    use plasma_parameters
    use control_mod, only: debug_mode
    use PolyLagrangeInterpolation

    implicit none

    integer :: ipoib, ipb, ipe
    double precision :: hrmax, r, rnext, recnsp
    double precision, dimension(:), allocatable :: x

    if (debug_mode) write(*,*) "Debug: coming in gengrid"
    nbaleqs = 4

    nder = 1
    npoi_der = 4
    allocate(x(npoi_der), coef(0:nder,npoi_der))

    hrmax = (rmax - rmin)/(npoimin + 1)

    npoib = 1
    r = rmin

    do while(r .lt. rmax)
      call recnsplit(r,recnsp)
      rnext = r + hrmax / recnsp
      call recnsplit(rnext,recnsp)
      r = 0.5d0 * (rnext + r + hrmax / recnsp)
      npoib = npoib + 1
    enddo
    npoic = npoib - 1

    allocate(rb(npoib), rc(npoic))
    allocate(Sb(npoib), Sc(npoic))

    r = rmin

    rb(1) = r

    do ipoib = 2,npoib
      call recnsplit(r,recnsp)
      rnext = r + hrmax / recnsp
      call recnsplit(rnext,recnsp)
      r = 0.5d0 * (rnext + r + hrmax / recnsp)
      rb(ipoib) = r
      rc(ipoib - 1) = 0.5d0 * (rb(ipoib - 1) + rb(ipoib))
    enddo

    if(npoi_der.gt.npoic) then
      stop 'gengrid : not enough grid points for derivatives'
    endif

    allocate(deriv_coef(npoi_der,npoib),ipbeg(npoib),ipend(npoib))
    allocate(reint_coef(npoi_der,npoib))

    do ipoib = 1, npoib
      ipb = ipoib - npoi_der / 2
      ipe = ipb + npoi_der - 1
      if(ipb .lt. 1) then
        ipb = 1
        ipe = ipb + npoi_der - 1
      elseif(ipe.gt.npoic) then
        ipe = npoic
        ipb = ipe - npoi_der + 1
      endif
      ipbeg(ipoib) = ipb
      ipend(ipoib) = ipe
      call plag_coeff(npoi_der, nder, rb(ipoib), rc(ipb:ipe), coef)
      deriv_coef(:,ipoib) = coef(1,:)
      reint_coef(:,ipoib) = coef(0,:)
    enddo

    deallocate(coef)

    allocate(params(nbaleqs,npoic),dot_params(nbaleqs,npoic))
    allocate(init_params(nbaleqs,npoic))
    allocate(params_b(nbaleqs,npoib),ddr_params(nbaleqs,npoib),ddr_params_nl(nbaleqs,npoib))
    allocate(params_lin(nbaleqs,npoic),params_b_lin(nbaleqs,npoib))
    allocate(fluxes_dif(nbaleqs,npoib),fluxes_con(nbaleqs,npoib),fluxes_con_nl(nbaleqs,npoib))

    if(iboutype.eq.1) then
      neqset=nbaleqs*(npoic-1)
    else
      neqset=nbaleqs*npoic
    endif
    allocate(y(neqset),dery(neqset),dery_equisource(neqset))
    !allocate(alpha(neqset*neqset))
    allocate(source_term(neqset))

    allocate(dae11(npoib),dae12(npoib),dae22(npoib))
    allocate(dai11(npoib),dai12(npoib),dai22(npoib))
    allocate(dni22(npoib),visca(npoib))
    allocate(dqle11(npoib),dqle12(npoib),dqle21(npoib),dqle22(npoib))
    allocate(dqli11(npoib),dqli12(npoib),dqli21(npoib),dqli22(npoib))
    allocate(de11(npoib),de12(npoib),de21(npoib),de22(npoib))
    allocate(di11(npoib),di12(npoib),di21(npoib),di22(npoib))
    allocate(T_EM_phi_e(npoib), T_EM_phi_i(npoib))
    allocate(T_EM_phi_e_source(npoib), T_EM_phi_i_source(npoib))
    allocate(polforce(npoib), polforce_ql(npoib), qlheat_e(npoib), qlheat_i(npoib))

    dni22=0.d0

    allocate(cneo(npoib),gpp_av(npoib))
    allocate(qsafb(npoib),qsaf(npoic))
    allocate(sqrt_g_times_B_theta_over_c(npoib),Ercov(npoib),sqg_bthet_overcavg(npoib), &
        Ercovavg(npoib))
    allocate(Ercov_lin(npoib))


    if (debug_mode) write(*,*) "Debug: going out in gengrid"
    return
end subroutine gengrid
  

subroutine calc_geometric_parameter_profiles

    use grid_mod
    use plasma_parameters
    use baseparam_mod

    implicit none

    integer :: ipoi
    double precision :: cneo_0,coullog,om_ci

    Sb=rb
    Sc=rc
    gpp_av=rtor**2

    coullog=15.d0
    om_ci=Z_i*e_charge*btor/(am*p_mass*c)
    cneo_0=1.32*4.d0*sqrt(pi)*Z_i**3*e_charge**4*coullog &
          /(3.d0*(am*p_mass)**1.5d0*om_ci**2)

    do ipoi=1,npoib
        qsafb(ipoi)=sum(qsaf(ipbeg(ipoi):ipend(ipoi))*reint_coef(:,ipoi))
        cneo(ipoi)=(rtor/rb(ipoi))**1.5d0*qsafb(ipoi)**2*cneo_0
        sqrt_g_times_B_theta_over_c(ipoi)=btor*rb(ipoi)/qsafb(ipoi)/c
    enddo

end subroutine calc_geometric_parameter_profiles


subroutine recnsplit(r,recnsp)

    use resonances_mod
    implicit none;

    integer :: k
    double precision :: r, recnsp

    if(prop) then
        prop=.false.
        call prepare_resonances
    endif

    !  recnsp = 1.d0 + gg_factor*exp(-((r-gg_r_res)/gg_width)**2);
    recnsp = 1.d0
    do k=1,numres
      recnsp = recnsp + ampl_res(k)*exp(-((r-r_res(k))/width_res(k))**2)
    enddo

end subroutine recnsplit


subroutine prepare_resonances

    use resonances_mod
    use grid_mod, only: gg_width, gg_factor,r_resonant
    use control_mod, only: ihdf5IO, debug_mode
    use h5mod
    use mpi
    use parallelTools, only: irank

    implicit none

    integer :: m,n,i,j,k,nr,jj,numres_orig
    double precision :: qres,qmin,qmax
    complex :: a
    integer, dimension(:), allocatable :: m_a,n_a
    integer, dimension(:), allocatable :: m_aa,n_aa
    double precision, dimension(:), allocatable :: r,q
    integer :: lb, ub

    iunit_res=157

    if (ihdf5IO .eq. 1) then
        CALL h5_init()
        CALL h5_open_rw(path2inp, h5_id)
        CALL h5_open_group(h5_id, '/preprocprof', group_id_1)
        CALL h5_get_bounds_1(group_id_1, 'q', lb, ub)
        allocate(r(ub),q(ub))
        CALL h5_get_double_1(group_id_1, 'q', q)
        !print *, 'q: ', q
        CALL h5_get_double_1(group_id_1, 'r_out', r)
        !print *, 'r: ', r
        CALL h5_close_group(group_id_1)
        CALL h5_close(h5_id)

        CALL h5_deinit()
        nr = ub

    else
        if (debug_mode) print *, "get number of q data points"
        nr=0
        open(iunit_res,file='profiles/q.dat')
        do
            read(iunit_res,*,end=1)
            nr=nr+1
        enddo
        1 continue
        close(iunit_res)
        allocate(r(nr),q(nr))

        if (debug_mode) print *, "reading profiles/q.dat"
        open(iunit_res,file='profiles/q.dat')
        do i=1,nr
            read(iunit_res,*) r(i),q(i)
        enddo
        close(iunit_res)
    end if
    q=abs(q)
    qmin=minval(q)
    qmax=maxval(q)

    open(iunit_res,file='flre/antenna.in')
    read(iunit_res,*)
    read(iunit_res,*)
    read(iunit_res,*)
    read(iunit_res,*)
    read(iunit_res,*)
    read(iunit_res,*) numres
    close(iunit_res)
    if (debug_mode) write(*,*) "numres = ", numres

    allocate(r_res(numres),width_res(numres),ampl_res(numres))
    allocate(r_resonant(numres))
    width_res=gg_width
    ampl_res=gg_factor
    allocate(m_a(numres),n_a(numres))
    allocate(m_aa(numres),n_aa(numres))
    numres_orig=numres

    open(iunit_res,file='flre/modes.in')
    k=1
    jj=1
    read(iunit_res,*) a
    m=nint(real(a))
    n=abs(nint(imag(a)))
    m_a(k)=m
    n_a(k)=n
    m_aa(jj)=m
    n_aa(jj)=n
    r_res(k)=abs(dfloat(m)/dfloat(n))
    outer: do i=2,numres
        read(iunit_res,*) a
        m=-abs(nint(real(a)))
        n=abs(nint(imag(a)))
        qres=abs(dfloat(m)/dfloat(n))
        !check for existence of resonant point:
        if(qres.lt.qmin.or.qres.gt.qmax) cycle
        !check for repeated resonance radii:
        jj=jj+1
        m_aa(jj)=m
        n_aa(jj)=n
        do j=1,k
          if(m*n_a(j)-n*m_a(j).eq.0) cycle outer
        enddo
        k=k+1
        m_a(k)=m
        n_a(k)=n
        r_res(k)=qres
    enddo outer
    close(iunit_res)
    numres=k

    do i=1,numres
        qres=r_res(i)
        do j=2,nr
            if((qres .ge. min(q(j-1), q(j))) .and. (qres .le. max(q(j-1), q(j)))) then
                r_res(i)=(r(j-1)*(q(j)-qres)+r(j)*(qres-q(j-1)))/(q(j)-q(j-1))
                exit
            endif
        enddo
    enddo

    if (irank .eq. 0 ) then
        if (debug_mode) print *,'Debug: gengrid: number of resonance points = ',numres
        do i = 1, numres
            ! maximum width for resonant radius is 999.999 cm with this format
            ! adjust if necessary
            write (*, "(a, i0, a, i0, a, f7.3, a)") &
                  "For mode (m,n) = (", m, ",", n, ") the resonant radius is at ", r_res(i), " cm."
        end do
    endif

    do i=1,numres_orig
        do j=1,numres
          if(m_aa(i)*n_a(j)-n_aa(i)*m_a(j).eq.0) then
              r_resonant(i)=r_res(j)
          end if
        enddo
    enddo
  

    deallocate(m_a,n_a,r,q)

  end subroutine prepare_resonances
