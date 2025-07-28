subroutine generate_k_space_grid(npoi_min, write_out, kr_cut)
    
    use grid
    use kr_grid, only: k_space_dim, kr, krp
    use setup
    use config, only: output_path, fstatus
    use constants, only: pi
    use KIM_kinds, only: dp

    implicit none

    integer :: i
    logical, intent(in) :: write_out
    real(dp), intent(in) :: kr_cut
    real(dp) :: h, krmin, krmax, hrmax, kr_val, krnext, recnsp
    integer :: k_grid_mode = 3
    real(dp), allocatable, dimension(:) :: delta_kr

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
        allocate(delta_kr(npoi_kr))

        kr_val = krmin
        kr(1) = krmin
        krp(1) = krmin

        do ipoib=2,npoi_kr-1
            call recnsplit_kr(kr_val,recnsp)
            krnext = kr_val + hrmax / recnsp
            call recnsplit_kr(krnext,recnsp)
            kr_val = 0.5d0 * (krnext + kr_val + hrmax / recnsp)
            kr(ipoib) = kr_val
            !rc(ipoib-1) = 0.5 * (rb(ipoib-1) + rb(ipoib))
            krp(ipoib) = kr_val !- 0.1d0
            delta_kr(ipoib-1) = kr(ipoib) - kr(ipoib-1)
        enddo
        krp(npoi_kr) = krmax
        kr(npoi_kr) = krmax

        k_space_dim = npoi_kr

    end if
    
    if (fstatus == 1) then 
        write(*,*) '- - - k space: - - -'
        write(*,*) '    new k space dim = ', k_space_dim
        write(*,*) '    min delta kr    = ', minval(delta_kr)
        write(*,*) '    max delta kr    = ', maxval(delta_kr)
        write(*,*) '- - - - - - - - - -'
    end if

    if (allocated(delta_kr)) deallocate(delta_kr)

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
    use KIM_kinds, only: dp

    implicit none

    real(dp), intent(in) :: kr_val
    real(dp), intent(out) :: recnsp

    recnsp = 1.0d0 +  kr_grid_ampl_res * exp(-((kr_val - kr_res) / kr_grid_width_res)**2)

end subroutine recnsplit_kr