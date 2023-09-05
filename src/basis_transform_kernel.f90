subroutine basis_transform_kernel(write_out)

    use kernel
    use grid
    use plas_parameter, only: iprof_length
    use config

    implicit none
    logical, intent(in) :: write_out
    integer :: i,j
    integer :: k1, k2

    double precision :: int_fac1, int_fac2

    if (fstatus == 1) write(*,*) 'Status: Transforming basis of kernels, Fourier -> Spline, write_out=',write_out

    allocate(K_rho_phi_llp(l_space_dim, l_space_dim),&
             K_rho_B_llp(l_space_dim, l_space_dim))

    if (.not. allocated(varphi_lkr)) then
        call calculate_fourier_trans_spline_funcs(.true.)
    end if

    do i=1, l_space_dim ! l
        do j=1, l_space_dim ! l'
            do k1=1, k_space_dim ! kr
                do k2 = 1, k_space_dim !kr'

                    if (k1==1 .or. k1==k_space_dim) then
                        int_fac1 = 0.5d0
                    else
                        int_fac1 = 1.0d0
                    end if
                    if (k2==1 .or. k2==k_space_dim) then
                        int_fac2 = 0.5d0
                    else
                        int_fac2 = 1.0d0
                    end if

                    K_rho_phi_llp(i,j) = K_rho_phi_llp(i,j) + int_fac1 * int_fac2 * varphi_lkr(i,k1) *&
                                         conjg(varphi_lkr(j,k2)) * K_rho_phi(k2,k1)
                    K_rho_B_llp(i,j) = K_rho_B_llp(i,j) + int_fac1 * int_fac2 * varphi_lkr(i,k1) *&
                                         conjg(varphi_lkr(j,k2)) * K_rho_B(k2,k1)
                    !if (isnan(real(K_rho_phi_llp(i,j)))) then
                    !    write(*,*) 'NAN: i=', i, ', j=', j, ', k1=', k1, ', k2=', k2
                    !end if
                end do
            end do
        end do
    end do

    K_rho_phi_llp = K_rho_phi_llp * ((kr(size(kr)) - kr(1)) / k_space_dim)**2
    K_rho_B_llp = K_rho_B_llp * ((kr(size(kr)) - kr(1)) / k_space_dim)**2

    if (write_out) call write_basis_trans_kernel

    contains

    subroutine write_basis_trans_kernel

        use config

        implicit none
        integer :: i,j
        logical :: ex

        inquire(file=trim(output_path)//'kernel', exist=ex)
        if(.not. ex) then
            call system('mkdir -p '//trim(output_path)//'kernel')
        end if

        open(unit=77, file=trim(output_path)//'kernel/K_rho_phi_llp_re.dat')
        open(unit=78, file=trim(output_path)//'kernel/K_rho_phi_llp_im.dat')

        open(unit=79, file=trim(output_path)//'kernel/K_rho_B_llp_re.dat')
        open(unit=80, file=trim(output_path)//'kernel/K_rho_B_llp_im.dat')

        do i=1, l_space_dim
            do j=1, l_space_dim
                write(77,*) real(K_rho_phi_llp(i,j))
                write(78,*) dimag(K_rho_phi_llp(i,j))
                write(79,*) real(K_rho_B_llp(i,j))
                write(80,*) dimag(K_rho_B_llp(i,j))
            end do
        end do

        close(77)
        close(78)

    end subroutine

end subroutine