! generate grid for spline function space
subroutine generate_l_space_grid

    use plas_parameter
    use back_quants
    use setup

    implicit none
    integer :: ind, ibeg, iend
    integer :: nlagr = 4
    integer :: nder = 0
    double precision, dimension(:,:), allocatable :: coef
    double precision :: r_res
    double precision :: q_res

    if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

    q_res = -dble(m_mode) / dble(n_mode)
    write(*,*) q_res
    ! find center, i.e. rational surface
    call binsrc(abs(q_prof), 1, iprof_length, abs(q_res), ind)
    ibeg = max(1, ind - nlagr/2)
    iend = ibeg + nlagr - 1
    if (iend .gt. iprof_length) then
        iend = iprof_length
        ibeg = iend - nlagr + 1
    end if

    call plag_coeff(nlagr, nder, q_res, q_prof(ibeg:iend), coef)

    r_res = sum(coef(0,:) * r_prof(ibeg:iend))
    write(*,*) r_res

end subroutine

subroutine generate_k_space_grid(write_out)
    
    use grid
    use setup
    use config, only: output_path, fstatus

    implicit none

    integer :: i
    logical, intent(in) :: write_out

    if (fstatus == 1) write(*,*) 'Status: Generating k-space grid, write out=', write_out

    allocate(kr(k_space_dim), krp(k_space_dim))

    do i=1, k_space_dim
        kr(i) = i
        krp(i) = i
    end do

    ! center around zero
    kr = kr - k_space_dim / 2d0
    krp = (krp+0.1d0) - k_space_dim / 2d0

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