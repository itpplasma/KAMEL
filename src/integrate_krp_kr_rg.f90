module integrate_krp_kr_rg

    use grid, only: npoib, kr, varphi_lkr, k_space_dim, krp
    use integrands

    implicit none


    contains

    subroutine integrate_kernel(write_out)

        use kernel, only: K_rho_phi_llp
        use integrands
        use integration

        implicit none

        logical, intent(in) :: write_out
        integer :: l, lp, j, k1
        double complex, dimension(1) :: res

        double complex, dimension(:,:), allocatable :: int_kr_rg
        double complex, dimension(:), allocatable :: int_rg

        allocate(int_kr_rg(k_space_dim, npoib))
        allocate(int_rg(npoib))

        if (.not. allocated(K_rho_phi_llp)) allocate(K_rho_phi_llp(l_space_dim, l_space_dim))


        if (.not. allocated(varphi_lkr)) then
            if (fstatus == 1) write(*,*) 'Status: Calculate Fourier transformed spline functions'
            call calculate_fourier_trans_spline_funcs(.true.)
        end if

        do l = 1, l_space_dim
            do lp = 1, l_space_dim
                do j = 1, npoib ! over r space
                    do k1 = 1, k_space_dim ! over kr space

                        call integrate_krp(res, size(res), krp(1), krp(size(krp)), j, k1,&
                                          lp, integrand_K_rho_phi_krp)
                        int_kr_rg(k1, j) = res(1)
                        write(*,*) 'res(1) = ', res(1)

                    end do

                    call integrate_kr(res, size(res), kr(1), kr(size(kr)), int_kr_rg, j, l, integrand_K_rho_phi_kr)
                    int_rg(j) = res(1)

                end do

                call integrate_rg(res, size(res), r_prof(1), r_prof(npoib), int_rg, integrand_K_rho_phi_rg)
                K_rho_phi_llp(lp,l) = res(1)

            end do
        end do

        if (write_out) call write_kernel_phi

        contains

            subroutine write_kernel_phi

                implicit none

                integer :: i,j
                logical :: ex

                inquire(file=trim(output_path)//'kernel', exist=ex)
                if (.not. ex) then
                    call system('mkdir -p '//trim(output_path)//'kernel')
                end if
                open(unit=77, file=trim(output_path)//'kernel/K_rho_phi_llp_re.dat')
                open(unit=78, file=trim(output_path)//'kernel/K_rho_phi_llp_im.dat')
                do i=1,l_space_dim
                    do j=1,l_space_dim
                        write(77,*) real(K_rho_phi_llp(i,j))
                        write(78,*) dimag(K_rho_phi_llp(i,j))
                    end do
                end do
                close(77)
                close(78)

            end subroutine



    end subroutine

end module