!
  subroutine gengrid(npoimin, write_out)

! Generates the grid for two types of boundary conditions at the outer
! boundary: iboutype=1 - fixed parameters, iboutype=2 - fixed fluxes
! At the inner boundary fixed fluxes = 0 are always assumed

    use grid
    use config, only: fdebug, output_path
    use plas_parameter, only: r_prof

    implicit none

    integer, intent(in) :: npoimin
    logical, intent(in) :: write_out
    integer :: ipoib, nder, ipb, ipe
    double precision :: hrmax, r, rnext, recnsp, rscale
    !double precision, dimension(:),   allocatable :: x
    double precision, dimension(:,:), allocatable :: coef

    if (fdebug == 1) write(*,*) "Debug: coming in gengrid"
    !nbaleqs=4

    nder=1
    npoi_der=4
    allocate(coef(0:nder,npoi_der))
    !allocate(x(npoi_der))

    rmin = minval(r_prof)
    rmax = maxval(r_prof)

    hrmax = (rmax - rmin)/(npoimin+1)

    npoib = 1
    r = rmin

    do while(r .lt. rmax)
        call recnsplit(r,recnsp)
        rnext = r + hrmax / recnsp
        call recnsplit(rnext,recnsp)
        r = 0.5d0 * (rnext + r + hrmax / recnsp)
        npoib= npoib + 1
    enddo

    npoic = npoib - 1

    allocate(rb(npoib), rc(npoic))
    allocate(Sb(npoib), Sc(npoic))

    r = rmin
    rb(1) = r

    do ipoib=2,npoib
        call recnsplit(r,recnsp)
        rnext = r + hrmax / recnsp
        call recnsplit(rnext,recnsp)
        r = 0.5d0 * (rnext + r + hrmax / recnsp)
        rb(ipoib) = r
        rc(ipoib-1) = 0.5 * (rb(ipoib-1) + rb(ipoib))
    enddo

    if(iboutype .eq. 1) then
        rscale = (rmax - rmin) / (rc(npoic) - rmin)
    else
        rscale = (rmax - rmin) / (rb(npoib) - rmin)
    endif
  !rb=rmin+rscale*(rb-rmin)
  !rc=rmin+rscale*(rc-rmin)

    if(npoi_der .gt. npoic) then
        write(*,*) '! Error : not enough grid points for derivatives'
        stop
    endif

  allocate(deriv_coef(npoi_der, npoib), ipbeg(npoib), ipend(npoib))
  allocate(reint_coef(npoi_der, npoib))

    do ipoib = 1, npoib
        ipb = ipoib - npoi_der / 2
        ipe = ipb + npoi_der - 1
        if(ipb .lt. 1) then
            ipb = 1
            ipe = ipb + npoi_der - 1
        elseif(ipe .gt. npoic) then
          ipe = npoic
          ipb = ipe - npoi_der + 1
        endif
        ipbeg(ipoib) = ipb
        ipend(ipoib) = ipe
        call plag_coeff(npoi_der, nder, rb(ipoib), rc(ipb:ipe), coef)
        deriv_coef(:, ipoib) = coef(1,:)
        reint_coef(:, ipoib) = coef(0,:)
    enddo

    deallocate(coef)

    if (write_out) call write_new_grid

    if (fdebug == 1) write(*,*) "Debug: going out in gengrid"


    contains

    subroutine write_new_grid

        implicit none
        integer :: i
        logical :: ex

        inquire(file = trim(output_path)//'grid', exist = ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'grid')
        end if

        open(unit = 77, file=trim(output_path)//'grid/rb.dat')
        do i = 1, npoib
            write(77,*) i, rb(i)
        end do
        close(77)

        open(unit = 78, file=trim(output_path)//'grid/rc.dat')
        do i = 1, npoic
            write(78,*) i, rc(i)
        end do
        close(78)

    end subroutine

end subroutine gengrid

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine recnsplit(r,recnsp)

    use resonances_mod
!  use grid_mod, only: gg_width, gg_factor, gg_r_res;
!
    implicit none;
!
    logical :: prop=.true.
    double precision :: r, recnsp;
!
    if(prop) then
        prop=.false.
        call prepare_resonances
    endif
!
!  recnsp = 1.d0 + gg_factor*exp(-((r-gg_r_res)/gg_width)**2);
    recnsp = 1.d0
  !do k=1,numres
    recnsp = recnsp +  ampl_res * exp(-((r - r_res) / width_res)**2)
  !recnsp = recnsp + ampl_res(k)*exp(-((r-r_res(k))/width_res(k))**2)
  !enddo
!
    return
end subroutine recnsplit


subroutine recnsplit_kr(r,recnsp)

    implicit none;
!
    double precision :: r, recnsp;
    double precision :: kr_res = 0.0d0
    double precision :: width_res = 2.0d0
    double precision :: ampl_res = 40.0d0
!
!  recnsp = 1.d0 + gg_factor*exp(-((r-gg_r_res)/gg_width)**2);
    recnsp = 1.d0
    recnsp = recnsp +  ampl_res * exp(-((r - kr_res) / width_res)**2)

end subroutine recnsplit_kr

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine prepare_resonances

    use resonances_mod
    use grid, only: gg_width, gg_factor,r_resonant
    use config, only: hdf5_output, fdebug
    use setup, only: m_mode, n_mode
    use plas_parameter, only: iprof_length, r_prof, q_prof
    !use h5mod
    !use mpi

    implicit none

    integer :: j
    double precision :: qres,qmin,qmax
    double precision, dimension(:), allocatable :: q
    integer :: lb, ub

    iunit_res=157

    allocate(q(iprof_length))

    q = abs(q_prof)
    qmin = minval(q)
    qmax = maxval(q)

    width_res = gg_width

    ampl_res = gg_factor

    qres = abs(dfloat(m_mode)/dfloat(n_mode))
    if(qres.lt.qmin.or.qres.gt.qmax) write(*,*) "Resonance location not found in q"

    r_res = qres

    do j=2,iprof_length
      if(qres .gt. q(j-1) .and. qres .lt. q(j)) then
        r_res = (r_prof(j-1) * (q(j) - qres) + r_prof(j) * (qres-q(j-1))) / (q(j)-q(j-1))
        exit
      endif
    enddo

    write(*,*) 'resonant radius: ',r_res
    deallocate(q)

end subroutine prepare_resonances


! generate grid for spline function space
!subroutine generate_l_space_grid
!
!    use plas_parameter
!    use back_quants
!    use setup
!
!    implicit none
!    integer :: ind, ibeg, iend
!    integer :: nlagr = 4
!    integer :: nder = 0
!    double precision, dimension(:,:), allocatable :: coef
!    double precision :: r_res
!    double precision :: q_res
!!
!    if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))
!
!    q_res = -dble(m_mode) / dble(n_mode)
!    write(*,*) q_res
!    ! find center, i.e. rational surface
!    call binsrc(abs(q_prof), 1, iprof_length, abs(q_res), ind)
!    ibeg = max(1, ind - nlagr/2)
!    iend = ibeg + nlagr - 1
!    if (iend .gt. iprof_length) then
!        iend = iprof_length
!        ibeg = iend - nlagr + 1
!    end if
!
!    call plag_coeff(nlagr, nder, q_res, q_prof(ibeg:iend), coef)
!
!    r_res = sum(coef(0,:) * r_prof(ibeg:iend))
!    write(*,*) r_res
!
!end subroutine

subroutine generate_k_space_grid(npoi_min, write_out)
    
    use grid
    use setup
    use config, only: output_path, fstatus
    use constants, only: pi

    implicit none

    integer :: i
    logical, intent(in) :: write_out
    double precision :: h, krmin, krmax, hrmax, kr_val, krnext, recnsp
    integer :: k_grid_mode = 2

    integer :: npoi_kr, npoi_min, ipoib

    if (fstatus == 1) write(*,*) 'Status: Generating k-space grid, write out=', write_out

    if (k_grid_mode==1) then
        allocate(kr(k_space_dim), krp(k_space_dim))
        do i=1, k_space_dim
            kr(i) = i
            krp(i) = i
        end do
        kr = kr - k_space_dim / 2d0
        krp = (krp+0.1d0) - k_space_dim / 2d0
    else if (k_grid_mode == 2) then
        h = pi / (k_space_dim + 1)
        allocate(kr(k_space_dim), krp(k_space_dim))
        do i=1, k_space_dim
            kr(i) = tan(- pi / 2.0d0 + i * h)
        end do
        krp = kr +0.01
    else if (k_grid_mode == 3) then ! non-equidistant grid, similar to l grid

        krmin = -k_space_dim / 3.0d0
        krmax = k_space_dim / 3.0d0

        hrmax = (krmax - krmin)/(npoi_min+1)

        npoi_kr = 1
        kr_val = krmin

        do while(kr_val .lt. krmax)
            call recnsplit_kr(kr_val,recnsp)
            krnext = kr_val + hrmax / recnsp
            call recnsplit_kr(krnext,recnsp)
            kr_val = 0.5d0 * (krnext + kr_val + hrmax / recnsp)
            npoi_kr = npoi_kr + 1
        enddo

        allocate(kr(npoi_kr), krp(npoi_kr))

        kr_val = krmin
        kr(1) = kr_val

        do ipoib=2,npoi_kr
            call recnsplit_kr(kr_val,recnsp)
            krnext = kr_val + hrmax / recnsp
            call recnsplit_kr(krnext,recnsp)
            kr_val = 0.5d0 * (krnext + kr_val + hrmax / recnsp)
            kr(ipoib) = kr_val
            !rc(ipoib-1) = 0.5 * (rb(ipoib-1) + rb(ipoib))
            krp(ipoib) = kr_val - 0.1d0
        enddo
        k_space_dim = npoi_kr
        if (fstatus == 1) write(*,*) ' Status: new k space dim = ', k_space_dim

    end if

    if (write_out) call write_k_space

    contains

        subroutine write_k_space

            implicit none

            open(unit = 78, file = trim(output_path)//'backs/kr.dat')
            open(unit = 79, file = trim(output_path)//'backs/krp.dat')
            do i = 1, k_space_dim
                write(78,*) kr(i)
                write(79,*) krp(i)
            end do
            close(78)
            close(79)

        end subroutine

end subroutine