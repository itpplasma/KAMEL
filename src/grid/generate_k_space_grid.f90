subroutine generate_k_space_grid(npoi_min, write_out, kr_cut)
    
    use grid
    use kr_grid, only: k_space_dim, kr, krp
    use setup
    use config, only: output_path, fstatus
    use constants, only: pi

    implicit none

    integer :: i
    logical, intent(in) :: write_out
    double precision, intent(in) :: kr_cut
    double precision :: h, krmin, krmax, hrmax, kr_val, krnext, recnsp
    integer :: k_grid_mode = 3

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
        krp = kr !+0.01d0
    else if (k_grid_mode == 3) then ! non-equidistant grid, similar to l grid

        krmin = - kr_cut
        krmax = kr_cut

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
        krp(1) = kr_val

        do ipoib=2,npoi_kr
            call recnsplit_kr(kr_val,recnsp)
            krnext = kr_val + hrmax / recnsp
            call recnsplit_kr(krnext,recnsp)
            kr_val = 0.5d0 * (krnext + kr_val + hrmax / recnsp)
            kr(ipoib) = kr_val
            !rc(ipoib-1) = 0.5 * (rb(ipoib-1) + rb(ipoib))
            krp(ipoib) = kr_val !- 0.1d0
        enddo
        k_space_dim = npoi_kr
        if (fstatus == 1) write(*,*) ' Status: new k space dim = ', k_space_dim

    end if

    if (write_out) call write_k_space

    contains

        subroutine write_k_space

            implicit none

            open(unit = 78, file = trim(output_path)//'grid/kr.dat')
            open(unit = 79, file = trim(output_path)//'grid/krp.dat')
            do i = 1, k_space_dim
                write(78,*) kr(i)
                write(79,*) krp(i)
            end do
            close(78)
            close(79)

        end subroutine

end subroutine


subroutine recnsplit_kr(kr_val, recnsp)

    use kr_grid, only: kr_grid_ampl_res, kr_grid_width_res, kr_res

    implicit none

    double precision, intent(in) :: kr_val
    double precision, intent(out) :: recnsp

    recnsp = 1.0d0 +  kr_grid_ampl_res * exp(-((kr_val - kr_res) / kr_grid_width_res)**2)

end subroutine recnsplit_kr