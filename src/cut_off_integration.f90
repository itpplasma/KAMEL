!> module cut_off_integration
!@brief: Determines the basis transformed kernels by the cut-off integration, 
! i.e. cut-off the integrations over kr and kr' prime and do partial 
! integration. All integrals are done with the trapezoidal method.
module cut_off_integration

    !use kernel, only: K_rho_phi_of_rg, K_rho_phi_llp, K_rho_B_llp
    use grid, only: k_space_dim, l_space_dim, r_space_dim, kr, krp, &
                    varphi_lkr, npoib, rb
    !use plas_parameter, only: r_prof
    use omp_lib
    use plas_parameter, only: rho_L
    use setup, only: cut_off_fac

    implicit none

    double precision :: kr_cutoff
    integer :: closest_kr_ind_upper
    integer :: closest_kr_ind_lower
    

    contains

    subroutine basis_transformation_integration(write_out)

        !use integrands, only: integrand_K_rho_phi_krp, integrand_K_rho_phi_kr
        !use integration, only: integrate_krp, integrate_kr
        use config, only: fstatus
        use kernel, only: fill_all_kernels, K_rho_phi_llp, K_rho_B_llp, K_rho_phi_of_rg
        use grid, only: xl
        use loading_bar

        implicit none

        integer :: l, lp, i_rg, i_k
        logical, intent(in) :: write_out
        integer :: count_elem_to_calc
        integer :: element_counter = 0

        double complex, dimension(1) :: res

        !double complex, dimension(:,:), allocatable :: int_kr_rg
        !double complex, dimension(:), allocatable :: int_rg
        !double complex, dimension(:,:,:), allocatable :: K_rho_phi_llp_rg

        !allocate(int_kr_rg(k_space_dim, npoib))
        !allocate(int_rg(npoib))

        write(*,*) "max(rho_Li) = ", rho_L

        ! kr space adjustments:
        ! determine cut-off in kr and corresponding indices
        kr_cutoff = 1.5d0 !cut_off_fac / rho_L
        call generate_k_space_grid(5, .true., kr_cutoff)
        closest_kr_ind_upper = findClosestIndex(kr, kr_cutoff)
        closest_kr_ind_lower = findClosestIndex(-kr, kr_cutoff)
        write(*,*) 'kr cut-off: ', kr_cutoff
        write(*,*) ' closest index lower: ', closest_kr_ind_lower, ', closest index upper: ', closest_kr_ind_upper


        if (.not. allocated(K_rho_phi_llp)) allocate(K_rho_phi_llp(l_space_dim, l_space_dim))
        if (.not. allocated(K_rho_B_llp)) allocate(K_rho_B_llp(l_space_dim, l_space_dim))
        !if (.not. allocated(K_rho_phi_llp_rg)) allocate(K_rho_phi_llp_rg(l_space_dim, l_space_dim, npoib))

        if (.not. allocated(varphi_lkr)) then
            if (fstatus == 1) write(*,*) 'Status: Calculate Fourier transformed spline functions'
            call calculate_fourier_trans_spline_funcs(.true.)
        end if

        ! fill kernel matrix
        if (.not. allocated(K_rho_phi_of_rg)) then
            !allocate(K_rho_phi(k_space_dim, k_space_dim))
            !call fill_kernel_rho_phi
            call fill_all_kernels
        end if

        if (fstatus==1) write(*,*) 'Status: Basis transformation'

        ! count the elements to be calculated as there is a cut-off for |xl - xl'|> rho_L
        count_elem_to_calc = 0
        do l = 1, l_space_dim
            do lp = 1, l_space_dim
                if (abs(xl(l) - xl(lp)) <= cut_off_fac * rho_L ) then
                    count_elem_to_calc = count_elem_to_calc + 1
                end if
            end do
        end do

        write(*,*) 'Total elements of each kernel array in l space: ', l_space_dim**2
        write(*,*) 'Elements to calculate for each kernel         : ', count_elem_to_calc
        write(*,*) 'I.e. ', dble(count_elem_to_calc) / dble(l_space_dim)**2 * 100.0d0, '% of the total matrix'

        element_counter = 0

        ! integrate over kr and krp for every spline base element
        !$OMP PARALLEL DO PRIVATE(l, lp, i_rg) &
        !$OMP SHARED(l_space_dim, xl, cut_off_fac, rho_L, element_counter, K_rho_phi_llp, rb, npoib) schedule(guided)
        do l = 1, l_space_dim ! first index of basis transformed kernel
            do lp = 1, l_space_dim ! second index of basis transformed kernel
                if (abs(xl(l) - xl(lp)) <= cut_off_fac * rho_L ) then

                    element_counter = element_counter + 1
                    ! integrate over r_g
                    do i_rg = 2, npoib
                        !K_rho_phi_llp_rg(l, lp, i_rg) = func_trapz_int_2D(l, lp, i_rg) ! + remaining terms
                        K_rho_phi_llp(l,lp) = 0.5d0 * (func_trapz_int_2D(l,lp, i_rg) + func_trapz_int_2D(l,lp, i_rg-1)) &
                            * (rb(i_rg) - rb(i_rg - 1))
                        K_rho_B_llp(l,lp) = 0.5d0 * (func_trapz_int_2D_rho_B(l,lp, i_rg) + func_trapz_int_2D_rho_B(l,lp, i_rg-1)) &
                            * (rb(i_rg) - rb(i_rg - 1))
                    end do
                    
                    call updateLoadingBar(element_counter, count_elem_to_calc)

                end if
            end do
        end do
        !$OMP END PARALLEL DO

        write(*,*) ''
        write(*,*) 'finished basis trafo'

        if (write_out) call write_kernel_in_spline_space

        contains

        subroutine write_kernel_in_spline_space

            use config, only: output_path

            implicit none

            integer :: i,j
            logical :: ex

            inquire(file=trim(output_path)//'kernel', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'kernel')
            end if
            open(unit=77, file=trim(output_path)//'kernel/K_rho_phi_llp_re.dat')
            open(unit=78, file=trim(output_path)//'kernel/K_rho_phi_llp_im.dat')
            open(unit=79, file=trim(output_path)//'kernel/K_rho_B_llp_re.dat')
            open(unit=80, file=trim(output_path)//'kernel/K_rho_B_llp_im.dat')
            do i=1,l_space_dim
                do j=1,l_space_dim
                    write(77,*) real(K_rho_phi_llp(i,j))
                    write(78,*) dimag(K_rho_phi_llp(i,j))
                    write(79,*) real(K_rho_B_llp(i,j))
                    write(80,*) dimag(K_rho_B_llp(i,j))
                end do
            end do
            close(77)
            close(78)
            close(79)
            close(80)

        end subroutine

    end subroutine



    subroutine fill_kernel_rho_phi

        use kernel_functions, only: kernel_rho_phi_of_kr_krp_rg
        use kernel, only: K_rho_phi_of_rg
        use config, only: fstatus
        use loading_bar

        implicit none
        integer :: i_kr, i_krp, i_rg
        integer :: count_loading = 0
        integer :: total_count_loading = 0


        if (fstatus == 1) write(*,*) 'Status: Fill kernel rho phi'

        if (.not. allocated(K_rho_phi_of_rg)) allocate(K_rho_phi_of_rg(k_space_dim, k_space_dim, npoib))

        total_count_loading = k_space_dim * k_space_dim * npoib
        
        !$OMP PARALLEL DO collapse(3) default(none) schedule(guided) &
        !$OMP PRIVATE(i_krp, i_kr, i_rg, total_count_loading) &
        !$OMP SHARED(K_rho_phi_of_rg, kr, krp, rb, k_space_dim, npoib, count_loading)
        do i_krp = 1, k_space_dim
            do i_kr = 1, k_space_dim
                do i_rg = 1, npoib
                    K_rho_phi_of_rg(i_krp, i_kr, i_rg) = kernel_rho_phi_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        if (fstatus == 1) write(*,*) 'Status: Finished filling kernel rho phi'

    end subroutine fill_kernel_rho_phi

    subroutine fill_kernel_rho_B

        use kernel_functions, only: kernel_rho_B_of_kr_krp_rg
        use kernel, only: K_rho_B_of_rg
        use config, only: fstatus
        use loading_bar

        implicit none
        integer :: i_kr, i_krp, i_rg
        integer :: count_loading = 0
        integer :: total_count_loading = 0


        if (fstatus == 1) write(*,*) 'Status: Fill kernel rho B'

        if (.not. allocated(K_rho_B_of_rg)) allocate(K_rho_B_of_rg(k_space_dim, k_space_dim, npoib))

        total_count_loading = k_space_dim * k_space_dim * npoib
        
        !$OMP PARALLEL DO collapse(3) default(none) schedule(guided) &
        !$OMP PRIVATE(i_krp, i_kr, i_rg, total_count_loading) &
        !$OMP SHARED(K_rho_B_of_rg, kr, krp, rb, k_space_dim, npoib, count_loading)
        do i_krp = 1, k_space_dim
            do i_kr = 1, k_space_dim
                do i_rg = 1, npoib
                    K_rho_B_of_rg(i_krp, i_kr, i_rg) = kernel_rho_B_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        if (fstatus == 1) write(*,*) 'Status: Finished filling kernel rho B'

    end subroutine fill_kernel_rho_B

    ! integrate 2D integral over k_r and k_r' with the trapezoidal rule
    double complex function func_trapz_int_2D(l, lp, i_rg)

        use setup, only: cut_off_fac
        use grid, only: xl
        use constants, only: com_unit
        use kernel, only: K_rho_phi_of_rg

        implicit none
        integer, intent(in) :: l, lp, i_rg
        integer :: i_kr, i_krp
                
        integer :: max_threads


        func_trapz_int_2D = 0.0d0

        !!$OMP PARALLEL DO default(none) schedule(guided) &
        !!$OMP PRIVATE(i_kr, i_krp) &
        !!$OMP SHARED(func_trapz_int_2D, kr, krp, k_space_dim, kr_cutoff, l, lp, i_rg, com_unit, &
        !!$OMP xl, rb, varphi_lkr, K_rho_phi_of_rg, closest_kr_ind_upper, closest_kr_ind_lower)
        do i_kr = closest_kr_ind_lower+1, closest_kr_ind_upper
            !if (kr_cutoff - abs(kr(i_kr)) .ge. 0.0d0) then
                do i_krp = closest_kr_ind_lower+1, closest_kr_ind_upper 
                    !if (kr_cutoff - abs(krp(i_krp)) .ge. 0.0d0) then
                        func_trapz_int_2D = func_trapz_int_2D + 0.25d0 * ((kr(i_kr)- kr(i_kr-1)) * (krp(i_krp) - krp(i_krp-1)) &
                            * (integrand_w_exp_facs(l,lp, i_kr, i_krp, i_rg) + integrand_w_exp_facs(l,lp,i_kr,i_krp-1, i_rg) &
                            + integrand_w_exp_facs(l,lp, i_kr-1, i_krp, i_rg) &
                            + integrand_w_exp_facs(l,lp, i_kr-1, i_krp-1, i_rg)))
                    !end if
                end do

                ! second term of transformation (integral over kr', kr at boundary)
                ! use loop over kr for integration over krp. That's why here the i_kr index is used.
                func_trapz_int_2D = func_trapz_int_2D + 0.5d0 * com_unit / (rb(i_rg) - xl(l)) &
                    * (varphi_lkr(closest_kr_ind_upper, l) &
                        * (exp(com_unit * krp(i_kr) * (xl(lp) - rb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_phi_of_rg(i_kr, closest_kr_ind_upper, i_rg) &
                        - exp(com_unit * krp(i_kr-1) * (xl(lp) - rb(i_rg))) * conjg(varphi_lkr(i_kr-1, lp))&
                        * K_rho_phi_of_rg(i_kr-1, closest_kr_ind_upper, i_rg)) &
                    ! second term in bracket:
                    - varphi_lkr(closest_kr_ind_lower, l) &
                        * (exp(com_unit * krp(i_kr) * (xl(lp) - rb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_phi_of_rg(i_kr, closest_kr_ind_lower, i_rg) &
                        - exp(com_unit * krp(i_kr-1) * (xl(lp) - rb(i_rg))) &
                        * conjg(varphi_lkr(i_kr-1, lp)) * K_rho_phi_of_rg(i_kr-1, closest_kr_ind_lower, i_rg))) &
                    * (krp(i_kr) - krp(i_kr-1))

                ! third term of transformation (integral over kr, kr' at boundary)
                func_trapz_int_2D = func_trapz_int_2D + 0.5d0 * com_unit / (xl(lp) - rb(i_rg)) &
                    ! first term in bracket, kr' upper limit
                    * (conjg(varphi_lkr(closest_kr_ind_upper, lp)) &
                        * (exp(com_unit * kr(i_kr) * (rb(i_rg) - xl(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_phi_of_rg(closest_kr_ind_upper, i_kr, i_rg) &
                        -  exp(com_unit * kr(i_kr-1) * (rb(i_rg) - xl(lp))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_phi_of_rg(closest_kr_ind_upper, i_kr-1, i_rg)) &
                    ! second term in bracket, kr' lower limit
                    - conjg(varphi_lkr(closest_kr_ind_lower, lp)) &
                        * (exp(com_unit * kr(i_kr) * (rb(i_rg) - xl(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_phi_of_rg(closest_kr_ind_lower, i_kr, i_rg) &
                        -  exp(com_unit * kr(i_kr-1) * (rb(i_rg) - xl(lp))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_phi_of_rg(closest_kr_ind_lower, i_kr-1, i_rg))) &
                    * (kr(i_kr) - kr(i_kr-1))

            !end if
            !!$OMP critical
            !if (OMP_GET_THREAD_NUM() == 0) then
            !    write(*,*) '    ', dble(i_kr) * 100.0d0/dble(k_space_dim), '% of kernel phi filled'
            !end if
            !!$OMP end critical
        end do
        !!$OMP END PARALLEL DO

        ! fourth term of transformation (term on the boundaries in kr and kr')
        func_trapz_int_2D = func_trapz_int_2D + 1.0d0 / ((xl(lp) - rb(i_rg)) * (xl(l) - rb(i_rg))) &
            * (conjg(varphi_lkr(closest_kr_ind_upper, lp)) &
            * (varphi_lkr(closest_kr_ind_upper, l) * K_rho_phi_of_rg(closest_kr_ind_upper, closest_kr_ind_upper, i_rg) &
            - varphi_lkr(closest_kr_ind_lower, l) * K_rho_phi_of_rg(closest_kr_ind_upper, closest_kr_ind_lower, i_rg)))&
            * (conjg(varphi_lkr(closest_kr_ind_lower, lp)) &
            * (varphi_lkr(closest_kr_ind_upper, l) * K_rho_phi_of_rg(closest_kr_ind_lower, closest_kr_ind_upper, i_rg) &
            - varphi_lkr(closest_kr_ind_lower, l) * K_rho_phi_of_rg(closest_kr_ind_lower, closest_kr_ind_lower, i_rg)))

    end function func_trapz_int_2D



    double complex function func_trapz_int_2D_rho_B(l, lp, i_rg)

        use setup, only: cut_off_fac
        use grid, only: xl
        use constants, only: com_unit
        use kernel, only: K_rho_B_of_rg

        implicit none
        integer, intent(in) :: l, lp, i_rg
        integer :: i_kr, i_krp
                
        func_trapz_int_2D_rho_B = 0.0d0

        !!$OMP PARALLEL DO default(none) schedule(guided) &
        !!$OMP PRIVATE(i_kr, i_krp) &
        !!$OMP SHARED(func_trapz_int_2D, kr, krp, k_space_dim, kr_cutoff, l, lp, i_rg, com_unit, &
        !!$OMP xl, rb, varphi_lkr, K_rho_phi_of_rg, closest_kr_ind_upper, closest_kr_ind_lower)
        do i_kr = closest_kr_ind_lower+1, closest_kr_ind_upper
            !if (kr_cutoff - abs(kr(i_kr)) .ge. 0.0d0) then
                do i_krp = closest_kr_ind_lower+1, closest_kr_ind_upper 
                    !if (kr_cutoff - abs(krp(i_krp)) .ge. 0.0d0) then
                        func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 0.25d0 &
                            * ((kr(i_kr)- kr(i_kr-1)) * (krp(i_krp) - krp(i_krp-1)) &
                            * (integrand_w_exp_facs_rho_B(l,lp, i_kr, i_krp, i_rg) &
                            + integrand_w_exp_facs_rho_B(l,lp,i_kr,i_krp-1, i_rg) &
                            + integrand_w_exp_facs_rho_B(l,lp, i_kr-1, i_krp, i_rg) &
                            + integrand_w_exp_facs_rho_B(l,lp, i_kr-1, i_krp-1, i_rg)))
                    !end if
                end do

                ! second term of transformation (integral over kr', kr at boundary)
                ! use loop over kr for integration over krp. That's why here the i_kr index is used.
                func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 0.5d0 * com_unit / (rb(i_rg) - xl(l)) &
                    * (varphi_lkr(closest_kr_ind_upper, l) &
                        * (exp(com_unit * krp(i_kr) * (xl(lp) - rb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_B_of_rg(i_kr, closest_kr_ind_upper, i_rg) &
                        - exp(com_unit * krp(i_kr-1) * (xl(lp) - rb(i_rg))) * conjg(varphi_lkr(i_kr-1, lp))&
                        * K_rho_B_of_rg(i_kr-1, closest_kr_ind_upper, i_rg)) &
                    ! second term in bracket:
                    - varphi_lkr(closest_kr_ind_lower, l) &
                        * (exp(com_unit * krp(i_kr) * (xl(lp) - rb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_B_of_rg(i_kr, closest_kr_ind_lower, i_rg) &
                        - exp(com_unit * krp(i_kr-1) * (xl(lp) - rb(i_rg))) &
                        * conjg(varphi_lkr(i_kr-1, lp)) * K_rho_B_of_rg(i_kr-1, closest_kr_ind_lower, i_rg))) &
                    * (krp(i_kr) - krp(i_kr-1))

                ! third term of transformation (integral over kr, kr' at boundary)
                func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 0.5d0 * com_unit / (xl(lp) - rb(i_rg)) &
                    ! first term in bracket, kr' upper limit
                    * (conjg(varphi_lkr(closest_kr_ind_upper, lp)) &
                        * (exp(com_unit * kr(i_kr) * (rb(i_rg) - xl(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_B_of_rg(closest_kr_ind_upper, i_kr, i_rg) &
                        -  exp(com_unit * kr(i_kr-1) * (rb(i_rg) - xl(lp))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_B_of_rg(closest_kr_ind_upper, i_kr-1, i_rg)) &
                    ! second term in bracket, kr' lower limit
                    - conjg(varphi_lkr(closest_kr_ind_lower, lp)) &
                        * (exp(com_unit * kr(i_kr) * (rb(i_rg) - xl(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_B_of_rg(closest_kr_ind_lower, i_kr, i_rg) &
                        -  exp(com_unit * kr(i_kr-1) * (rb(i_rg) - xl(lp))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_B_of_rg(closest_kr_ind_lower, i_kr-1, i_rg))) &
                    * (kr(i_kr) - kr(i_kr-1))

            !end if
            !!$OMP critical
            !if (OMP_GET_THREAD_NUM() == 0) then
            !    write(*,*) '    ', dble(i_kr) * 100.0d0/dble(k_space_dim), '% of kernel phi filled'
            !end if
            !!$OMP end critical
        end do
        !!$OMP END PARALLEL DO

        ! fourth term of transformation (term on the boundaries in kr and kr')
        func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 1.0d0 / ((xl(lp) - rb(i_rg)) * (xl(l) - rb(i_rg))) &
            * (conjg(varphi_lkr(closest_kr_ind_upper, lp)) &
            * (varphi_lkr(closest_kr_ind_upper, l) * K_rho_B_of_rg(closest_kr_ind_upper, closest_kr_ind_upper, i_rg) &
            - varphi_lkr(closest_kr_ind_lower, l) * K_rho_B_of_rg(closest_kr_ind_upper, closest_kr_ind_lower, i_rg)))&
            * (conjg(varphi_lkr(closest_kr_ind_lower, lp)) &
            * (varphi_lkr(closest_kr_ind_upper, l) * K_rho_B_of_rg(closest_kr_ind_lower, closest_kr_ind_upper, i_rg) &
            - varphi_lkr(closest_kr_ind_lower, l) * K_rho_B_of_rg(closest_kr_ind_lower, closest_kr_ind_lower, i_rg)))

    end function func_trapz_int_2D_rho_B



    double complex function integrand_w_exp_facs(l, lp, i_kr, i_krp, i_rg)

        use kernel_functions, only: kernel_rho_phi_of_kr_krp_rg
        use kernel, only: K_rho_phi_of_rg
        use constants, only: com_unit
        use grid, only: varphi_lkr, xl

        implicit none
        integer, intent(in) :: l, lp
        integer, intent(in) :: i_kr, i_krp
        integer, intent(in) :: i_rg

        !double complex tilde_varphi_lkr

        !double complex :: kernel_rho_phi_of_kr_krp_rg
 
        integrand_w_exp_facs = 0.0d0

        ! return \tilde{\varphi}_l(k_r) \tilde{\varphi}_{l'}^*(k_r') K(k_r, k_r'; r_g) exp(i k_r(r_g - x_l)+i k_r'(x_l' - r_g)) 

        integrand_w_exp_facs = varphi_lkr(i_kr, l) &!rb(l), rb(l-1), rb(l+1), rb(l+2), kr(i_kr)) &
            * conjg(varphi_lkr(i_krp, lp)) &!, rb(lp-1), rb(lp+1), rb(lp+2), krp(i_krp))) &
            * exp(com_unit * kr(i_kr) * (rb(i_rg)-xl(l)) + com_unit * krp(i_krp) * (xl(lp)-rb(i_rg))) &
            * K_rho_phi_of_rg(i_krp, i_kr, i_rg)

    end function integrand_w_exp_facs

    double complex function integrand_w_exp_facs_rho_B(l, lp, i_kr, i_krp, i_rg)

        use kernel_functions, only: kernel_rho_B_of_kr_krp_rg
        use kernel, only: K_rho_B_of_rg
        use constants, only: com_unit
        use grid, only: varphi_lkr, xl

        implicit none
        integer, intent(in) :: l, lp
        integer, intent(in) :: i_kr, i_krp
        integer, intent(in) :: i_rg

        !double complex tilde_varphi_lkr

        !double complex :: kernel_rho_phi_of_kr_krp_rg
 
        integrand_w_exp_facs_rho_B = 0.0d0

        ! return \tilde{\varphi}_l(k_r) \tilde{\varphi}_{l'}^*(k_r') K(k_r, k_r'; r_g) exp(i k_r(r_g - x_l)+i k_r'(x_l' - r_g)) 

        integrand_w_exp_facs_rho_B = varphi_lkr(i_kr, l) &!rb(l), rb(l-1), rb(l+1), rb(l+2), kr(i_kr)) &
            * conjg(varphi_lkr(i_krp, lp)) &!, rb(lp-1), rb(lp+1), rb(lp+2), krp(i_krp))) &
            * exp(com_unit * kr(i_kr) * (rb(i_rg)-xl(l)) + com_unit * krp(i_krp) * (xl(lp)-rb(i_rg))) &
            * K_rho_B_of_rg(i_krp, i_kr, i_rg)

    end function integrand_w_exp_facs_rho_B


    integer function findClosestIndex(array, target)

        implicit none
        real(kind=8), intent(in) :: array(:)
        real(kind=8), intent(in) :: target
        integer :: closest_index
        real(kind=8) :: min_difference
        integer :: i

        ! Initialize with a large value
        min_difference = huge(1.0)

        ! Loop through the array to find the closest element
        do i = 1, size(array)
            if (abs(array(i) - target) < min_difference) then
                min_difference = abs(array(i) - target)
                closest_index = i
            end if
        end do

        ! Return the index of the closest element
        findClosestIndex = closest_index
  end function findClosestIndex

end module