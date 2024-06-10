!> module cut_off_integration
!@brief: Determines the basis transformed kernels by the cut-off integration, 
! i.e. cut-off the integrations over kr and kr' prime and do partial 
! integration. All integrals are done with the trapezoidal method.
module cut_off_integration

    use grid, only: varphi_lkr, rg_grid, xl_grid
    use resonances_mod, only: r_res
    use omp_lib
    use plasma_parameter, only: rho_L
    use setup, only: cut_off_fac, kr_cut_off_fac, eps_reg
    use config, only: fstatus
    use kernels, only: fill_rho_kernels, K_rho_phi_llp, K_rho_B_llp, K_rho_phi_of_rg
    use loading_bar
    use constants, only: pi, com_unit
    use use_libcerf, only: cerf_F

    implicit none

    double precision :: kr_cutoff
    
    contains

    subroutine basis_transformation_of_kernels(write_out)

        implicit none

        integer :: l, lp, i_rg, i_k
        logical, intent(in) :: write_out
        integer :: count_elem_to_calc
        integer :: element_counter = 0

        !call kr_space_adjustments

        if (.not. allocated(K_rho_phi_llp)) allocate(K_rho_phi_llp(xl_grid%npts_b, xl_grid%npts_b))
        if (.not. allocated(K_rho_B_llp)) allocate(K_rho_B_llp(xl_grid%npts_b, xl_grid%npts_b))
        K_rho_phi_llp = 0.0d0
        K_rho_B_llp = 0.0d0

        if (.not. allocated(varphi_lkr)) then
            if (fstatus == 1) write(*,*) 'Status: Calculate Fourier transformed spline functions'
            call calculate_fourier_trans_spline_funcs(.true.)
        end if

        ! fill kernel matrix
        if (.not. allocated(K_rho_phi_of_rg)) then
            call fill_rho_kernels
        end if

        if (fstatus==1) write(*,*) 'Status: Basis transformation'

        ! count the elements to be calculated as there is a cut-off for |xl - xl'|> rho_L
        count_elem_to_calc = 0
        do l = 1, xl_grid%npts_b
            do lp = 1, xl_grid%npts_b
                !if (abs(xl(l) - xl(lp)) <= cut_off_fac * rho_L ) then
                if (abs(l-lp) <= 1) then
                    count_elem_to_calc = count_elem_to_calc + 1
                end if
            end do
        end do

        write(*,*) 'Total elements of each kernel array in l space: ', xl_grid%npts_b**2
        write(*,*) 'Elements to calculate for each kernel         : ', count_elem_to_calc
        write(*,*) 'I.e. ', dble(count_elem_to_calc) / dble(xl_grid%npts_b)**2 * 100.0d0, '% of the total matrix'

        element_counter = 0

        ! integrate over kr and krp for every spline base element
        !$OMP PARALLEL DO PRIVATE(l, lp, i_rg) &
        !$OMP SHARED(xl_grid, cut_off_fac, rho_L, element_counter, K_rho_phi_llp, rg_grid) schedule(guided)
        do l = 1, xl_grid%npts_b ! first index of basis transformed kernel
            do lp = 1, xl_grid%npts_b ! second index of basis transformed kernel
                !if (abs(xl(l) - xl(lp)) <= cut_off_fac * rho_L ) then
                if (abs(l-lp) <= 1) then

                    element_counter = element_counter + 1
                    ! integrate over r_g
                    do i_rg = 2, rg_grid%npts_b
                        !K_rho_phi_llp_rg(l, lp, i_rg) = func_trapz_int_2D(l, lp, i_rg) ! + remaining terms
                        K_rho_phi_llp(l,lp) = K_rho_phi_llp(l,lp) + 0.5d0 * &
                        (func_trapz_int_2D_rho_phi(l,lp, i_rg)   * exp(- eps_reg * (rg_grid%xb(i_rg)   - r_res)**2) &
                        +func_trapz_int_2D_rho_phi(l,lp, i_rg-1) * exp(- eps_reg * (rg_grid%xb(i_rg-1) - r_res)**2) )&
                            * (rg_grid%xb(i_rg) - rg_grid%xb(i_rg - 1))  
                        !K_rho_B_llp(l,lp) = K_rho_B_llp(l,lp) + 0.5d0 * &
                        !(func_trapz_int_2D_rho_B(l,lp, i_rg)    * exp(- eps_reg * (rb(i_rg)   - r_res)**2) &
                        !+func_trapz_int_2D_rho_B(l,lp, i_rg-1)  * exp(- eps_reg * (rb(i_rg-1) - r_res)**2))  &
                        !    * (rb(i_rg) - rb(i_rg - 1)) * (sqrt(eps_reg / pi))
                    end do
                    
                    call updateLoadingBar(element_counter, count_elem_to_calc)

                end if
            end do
        end do
        !$OMP END PARALLEL DO

        ! From Fourier transformation of spline functions
        K_rho_phi_llp = K_rho_phi_llp / (2.0d0 * pi) * (sqrt(eps_reg / pi))  !* (2.0d0 * pi**2)
        !K_rho_phi_llp = K_rho_phi_llp / (2.0d0 * pi)
        ! error functions are from normalization of the regularization with a Gaussian over a finite interval
        !K_rho_phi_llp = K_rho_phi_llp / (2.0d0 * pi) * (2 * sqrt(eps_reg / pi)) &
                        !/ (cerf_F(cmplx((rb(number_points_rg) - r_res) * sqrt(eps_reg), 0.0d0, kind=kind(1.0d0))) &
                        !-  cerf_F(cmplx((rb(1) - r_res) * sqrt(eps_reg), 0.0d0, kind=kind(1.0d0))))  !&
                        !* (4.0d0 * pi**2) * exp(1.0d0/2.0d0)
        !K_rho_B_llp = K_rho_B_llp / (2.0d0 * pi) * (2 * sqrt(eps_reg / pi)) &
        !                / (cerf_F(cmplx((rb(number_points_rg) - r_res) * sqrt(eps_reg), 0.0d0, kind=kind(1.0d0))) &
        !                -  cerf_F(cmplx((rb(1) - r_res) * sqrt(eps_reg), 0.0d0, kind=kind(1.0d0))))  &
        !                * (2.0d0 * pi**2)


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
            do i=1, xl_grid%npts_b
                do j=1, xl_grid%npts_b
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


    ! integrate 2D integral over k_r and k_r' with the trapezoidal rule
    double complex function func_trapz_int_2D_rho_phi(l, lp, i_rg)

        use grid, only: xl_grid, kr_grid, krp_grid, rg_grid
        use kernels, only: K_rho_phi_of_rg
        use integrands, only: integrand_w_exp_facs_rho_phi

        implicit none
        integer, intent(in) :: l, lp, i_rg
        integer :: i_kr, i_krp
        double precision :: rg, kr, krm1, krp, krpm1, xl, xlp

        xl = xl_grid%xb(l)
        xlp = xl_grid%xb(lp)
        rg = rg_grid%xb(i_rg)
                
        func_trapz_int_2D_rho_phi = 0.0d0

        do i_kr = 2, kr_grid%npts_b

            kr = kr_grid%xb(i_kr)
            krm1 = kr_grid%xb(i_kr-1)

            do i_krp = 2, krp_grid%npts_b

                krp = krp_grid%xb(i_krp)
                krpm1 = krp_grid%xb(i_krp-1)

                func_trapz_int_2D_rho_phi = func_trapz_int_2D_rho_phi + 0.25d0 * ((kr - krm1) * (krp - krpm1) &
                    * (integrand_w_exp_facs_rho_phi(l,lp, i_kr, i_krp, i_rg)& 
                    + integrand_w_exp_facs_rho_phi(l,lp,i_kr,i_krp-1, i_rg) &
                    + integrand_w_exp_facs_rho_phi(l,lp, i_kr-1, i_krp, i_rg) &
                    + integrand_w_exp_facs_rho_phi(l,lp, i_kr-1, i_krp-1, i_rg)))
            end do

                ! second term of transformation (integral over kr', kr at boundary)
                ! use loop over kr for integration over krp. That's why here the i_kr index is used.
                func_trapz_int_2D_rho_phi = func_trapz_int_2D_rho_phi + 0.5d0 * com_unit / (rg_grid%xb(i_rg) - xl_grid%xb(l)) &
                    * (varphi_lkr(kr_grid%npts_b, l) * exp(com_unit * kr_grid%max_val * (rg_grid%xb(i_rg) - xl_grid%xb(l))) &
                        ! first integral
                        * (exp(com_unit * krp_grid%xb(i_kr) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_phi_of_rg(i_kr, kr_grid%npts_b, i_rg) &
                        - exp(com_unit * krp_grid%xb(i_kr-1) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) * conjg(varphi_lkr(i_kr-1, lp))&
                        * K_rho_phi_of_rg(i_kr-1, kr_grid%npts_b, i_rg)) &
                    ! second term in bracket:
                    - varphi_lkr(1, l) * exp(com_unit * kr_grid%min_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                        * (exp(com_unit * krp_grid%xb(i_kr) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_phi_of_rg(i_kr, 1, i_rg) &
                        - exp(com_unit * krp_grid%xb(i_kr-1) * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))* conjg(varphi_lkr(i_kr-1, lp))&
                        * K_rho_phi_of_rg(i_kr-1, 1, i_rg))) &
                    * (krp_grid%xb(i_kr) - krp_grid%xb(i_kr-1))

                !! third term of transformation (integral over kr, kr' at boundary)
                func_trapz_int_2D_rho_phi = func_trapz_int_2D_rho_phi + 0.5d0 * com_unit / (xl_grid%xb(lp) - rg_grid%xb(i_rg)) &
                    ! first term in bracket, kr' upper limit
                    * (conjg(varphi_lkr(kr_grid%npts_b, lp)) &
                    * exp(com_unit * kr_grid%max_val * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
                        * (exp(com_unit * kr_grid%xb(i_kr) * (rg_grid%xb(i_rg) &
                        - xl_grid%xb(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_phi_of_rg(kr_grid%npts_b, i_kr, i_rg) &
                        -  exp(com_unit * kr_grid%xb(i_kr-1) * (rg_grid%xb(i_rg) - xl_grid%xb(l))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_phi_of_rg(kr_grid%npts_b, i_kr-1, i_rg)) &
                    ! second term in bracket, kr' lower limit
                    - conjg(varphi_lkr(1, lp)) * exp(com_unit * kr_grid%min_val * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
                        * (exp(com_unit * kr_grid%xb(i_kr) * (rg_grid%xb(i_rg) - xl_grid%xb(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_phi_of_rg(1, i_kr, i_rg) &
                        -  exp(com_unit * kr_grid%xb(i_kr-1) * (rg_grid%xb(i_rg) - xl_grid%xb(l))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_phi_of_rg(1, i_kr-1, i_rg))) &
                    * (kr_grid%xb(i_kr) - kr_grid%xb(i_kr-1))
        end do

        ! fourth term of transformation (term on the boundaries in kr and kr')
        func_trapz_int_2D_rho_phi = func_trapz_int_2D_rho_phi + 1.0d0 / ((xl_grid%xb(lp) - rg_grid%xb(i_rg)) &
            * (xl_grid%xb(l) - rg_grid%xb(i_rg))) &
            * (conjg(varphi_lkr(kr_grid%npts_b, lp)) * exp(com_unit * kr_grid%max_val * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
            * (varphi_lkr(kr_grid%npts_b, l) * exp(com_unit * kr_grid%max_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_phi_of_rg(kr_grid%npts_b, kr_grid%npts_b, i_rg) &
            - varphi_lkr(1, l)* exp(com_unit * kr_grid%min_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_phi_of_rg(kr_grid%npts_b, 1, i_rg)))&
            * (conjg(varphi_lkr(1, lp)) * exp(com_unit * kr_grid%min_val * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
            * (varphi_lkr(kr_grid%npts_b, l) * exp(com_unit * kr_grid%max_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_phi_of_rg(1, kr_grid%npts_b, i_rg) &
            - varphi_lkr(1, l) * exp(com_unit * kr_grid%min_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_phi_of_rg(1, 1, i_rg)))

        if (isnan(real(func_trapz_int_2D_rho_phi))) then
            write(*,*) 'NaN in func_trapz_int_2D_rho_phi, l = ', l, ', lp = ', lp, ', i_rg = ', i_rg
            !stop
        end if


    end function func_trapz_int_2D_rho_phi



    double complex function func_trapz_int_2D_rho_B(l, lp, i_rg)

        use grid, only: xl_grid, kr_grid, krp_grid
        use kernels, only: K_rho_B_of_rg
        use integrands, only: integrand_w_exp_facs_rho_B

        implicit none
        integer, intent(in) :: l, lp, i_rg
        integer :: i_kr, i_krp
                
        func_trapz_int_2D_rho_B = 0.0d0

        do i_kr = 1+1, kr_grid%npts_b - 1
            do i_krp = 1+1, krp_grid%npts_b - 1 
                func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 0.25d0 &
                    * ((kr_grid%xb(i_kr)- kr_grid%xb(i_kr-1)) * (krp_grid%xb(i_krp) - krp_grid%xb(i_krp-1)) &
                    * (integrand_w_exp_facs_rho_B(l,lp, i_kr, i_krp, i_rg) &
                    + integrand_w_exp_facs_rho_B(l,lp,i_kr,i_krp-1, i_rg) &
                    + integrand_w_exp_facs_rho_B(l,lp, i_kr-1, i_krp, i_rg) &
                    + integrand_w_exp_facs_rho_B(l,lp, i_kr-1, i_krp-1, i_rg)))
            end do

                ! second term of transformation (integral over kr', kr at boundary)
                ! use loop over kr for integration over krp. That's why here the i_kr index is used.
                func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 0.5d0 * com_unit / (rg_grid%xb(i_rg) - xl_grid%xb(l)) &
                    * (varphi_lkr(kr_grid%npts_b, l) * exp(com_unit * kr_grid%max_val * (rg_grid%xb(i_rg) - xl_grid%xb(l))) &
                        * (exp(com_unit * krp_grid%xb(i_kr) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_B_of_rg(i_kr, kr_grid%npts_b, i_rg) &
                        - exp(com_unit * krp_grid%xb(i_kr-1) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) * conjg(varphi_lkr(i_kr-1, lp))&
                        * K_rho_B_of_rg(i_kr-1, kr_grid%npts_b, i_rg)) &
                    ! second term in bracket:
                    - varphi_lkr(1, l) * exp(com_unit * kr_grid%min_val * (rg_grid%xb(i_rg) - xl_grid%xb(l))) &
                        * (exp(com_unit * krp_grid%xb(i_kr) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) * conjg(varphi_lkr(i_kr, lp)) &
                        * K_rho_B_of_rg(i_kr, 1, i_rg) &
                        - exp(com_unit * krp_grid%xb(i_kr-1) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) &
                        * conjg(varphi_lkr(i_kr-1, lp)) * K_rho_B_of_rg(i_kr-1, 1, i_rg))) &
                    * (krp_grid%xb(i_kr) - krp_grid%xb(i_kr-1))

                ! third term of transformation (integral over kr, kr' at boundary)
                func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 0.5d0 * com_unit / (xl_grid%xb(lp) - rg_grid%xb(i_rg)) &
                    ! first term in bracket, kr' upper limit
                    * (conjg(varphi_lkr(kr_grid%npts_b, lp)) * exp(com_unit * kr_grid%max_val * &
                    (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
                        * (exp(com_unit * kr_grid%xb(i_kr) * (rg_grid%xb(i_rg) - xl_grid%xb(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_B_of_rg(kr_grid%npts_b, i_kr, i_rg) &
                        -  exp(com_unit * kr_grid%xb(i_kr-1) * (rg_grid%xb(i_rg) - xl_grid%xb(lp))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_B_of_rg(kr_grid%npts_b, i_kr-1, i_rg)) &
                    ! second term in bracket, kr' lower limit
                    - conjg(varphi_lkr(1, lp)) * exp(com_unit * kr_grid%min_val * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
                        * (exp(com_unit * kr_grid%xb(i_kr) * (rg_grid%xb(i_rg) - xl_grid%xb(l))) * varphi_lkr(i_kr, l) &
                        * K_rho_B_of_rg(1, i_kr, i_rg) &
                        -  exp(com_unit * kr_grid%xb(i_kr-1) * (rg_grid%xb(i_rg) - xl_grid%xb(lp))) &
                        * varphi_lkr(i_kr-1, l) * K_rho_B_of_rg(1, i_kr-1, i_rg))) &
                    * (kr_grid%xb(i_kr) - kr_grid%xb(i_kr-1))

        end do

        func_trapz_int_2D_rho_B = func_trapz_int_2D_rho_B + 1.0d0 / ((xl_grid%xb(lp) - rg_grid%xb(i_rg)) &
            * (xl_grid%xb(l) - rg_grid%xb(i_rg))) &
            * (conjg(varphi_lkr(kr_grid%npts_b, lp)) * exp(com_unit * kr_grid%max_val * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
            * (varphi_lkr(kr_grid%npts_b, l) * exp(com_unit * kr_grid%max_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_B_of_rg(kr_grid%npts_b, kr_grid%npts_b, i_rg) &
            - varphi_lkr(1, l)* exp(com_unit * kr_grid%min_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_B_of_rg(kr_grid%npts_b, 1, i_rg)))&
            * (conjg(varphi_lkr(1, lp)) * exp(com_unit * kr_grid%min_val * (xl_grid%xb(lp) - rg_grid%xb(i_rg)))&
            * (varphi_lkr(kr_grid%npts_b, l) * exp(com_unit * kr_grid%max_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_B_of_rg(1, kr_grid%npts_b, i_rg) &
            - varphi_lkr(1, l) * exp(com_unit * kr_grid%min_val * (rg_grid%xb(i_rg) - xl_grid%xb(l)))&
                * K_rho_B_of_rg(1, 1, i_rg)))

    end function func_trapz_int_2D_rho_B


end module

