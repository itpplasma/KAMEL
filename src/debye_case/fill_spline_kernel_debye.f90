subroutine fill_spline_kernel_debye(write_out)

    !use integrands, only: integrand_K_rho_phi_krp, integrand_K_rho_phi_kr
    !use integration, only: integrate_krp, integrate_kr
    use config, only: fstatus
    use debye_kernel, only: func_debye_kernel
    use grid, only: xl, varphi_lkr, l_space_dim
    use kernel, only: K_rho_phi_llp, K_rho_B_llp

    implicit none

    integer :: l, lp, i_rg, i_k
    logical, intent(in) :: write_out
    integer :: count_elem_to_calc

    if (fstatus==1) write(*,*) 'Status: Basis transformation'

    call kr_space_adjustments

    if (.not. allocated(varphi_lkr)) then
        if (fstatus == 1) write(*,*) 'Status: Calculate Fourier transformed spline functions'
        call calculate_fourier_trans_spline_funcs(.true.)
    end if

    if (.not. allocated(K_rho_phi_llp)) allocate(K_rho_phi_llp(l_space_dim, l_space_dim))
    K_rho_phi_llp = 0.0d0

    if (.not. allocated(K_rho_B_llp)) allocate(K_rho_B_llp(l_space_dim, l_space_dim))
    K_rho_B_llp = 0.0d0

    do l = 1, l_space_dim ! first index of basis transformed kernel
        K_rho_phi_llp(l,l) = func_debye_kernel(0.0d0,0.0d0)!int_debye_kernel(l,lp)
    end do

    write(*,*) ''
    write(*,*) 'finished basis trafo'

    if (write_out) call write_kernel_in_spline_space_debye

    contains

    subroutine write_kernel_in_spline_space_debye

        use config, only: output_path

        implicit none

        integer :: i,j
        logical :: ex

        inquire(file=trim(output_path)//'kernel', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'kernel')
        end if
        open(unit=77, file=trim(output_path)//'kernel/K_rho_phi_llp_re_debye.dat')
        open(unit=78, file=trim(output_path)//'kernel/K_rho_phi_llp_im_debye.dat')
        do i=1,l_space_dim
            do j=1,l_space_dim
                write(77,*) real(K_rho_phi_llp(i,j))
                write(78,*) dimag(K_rho_phi_llp(i,j))
            end do
        end do
        close(77)
        close(78)
        close(79)
        close(80)

    end subroutine

end subroutine
