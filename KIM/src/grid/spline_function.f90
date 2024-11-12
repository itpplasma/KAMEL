! Subroutine to calculate the Fourier transformed spline functions
! For the spline functions, hat functions are used, where the Fourier
! transform was determined analytically
subroutine calculate_fourier_trans_spline_funcs(write_out)

    use constants, only: com_unit
    use grid, only: rg_grid, varphi_lkr, &
        grid_spacing, spline_base, xl_grid, kr_grid
    use config, only: output_path
    !use plasma_parameter, only: r_prof, iprof_length

    implicit none
    integer :: i,n
    logical, intent(in) :: write_out

    allocate(varphi_lkr(kr_grid%npts_b, xl_grid%npts_b))

    if (spline_base == 1) then
        ! hat functions
        if (grid_spacing == 1) then
        ! equidistant
            do i=1, kr_grid%npts_b
                do n =2, xl_grid%npts_b-1
                    !varphi_lkr(i,n) = FT_hat_function_e(xl_grid%xb(n), xl_grid%xb(n+1), kr_grid%xb(i))
                    varphi_lkr(i,n) = tilde_varphi_lkr_e(xl_grid%xb(n), xl_grid%xb(n+1), kr_grid%xb(i))
                end do
                varphi_lkr(i,1) =  1.0d0 / ((xl_grid%xb(2) - xl_grid%xb(1)) * kr_grid%xb(i)**2) * &
                                  (1.0d0 - exp(-com_unit * kr_grid%xb(i) * (xl_grid%xb(2) - xl_grid%xb(1))) &
                                  - com_unit * kr_grid%xb(i) * (xl_grid%xb(2) - xl_grid%xb(1)))
                varphi_lkr(i,xl_grid%npts_b) = 1.0d0 / ((xl_grid%xb(xl_grid%npts_b) &
                                  - xl_grid%xb(xl_grid%npts_b-1)) * kr_grid%xb(i)**2) &
                                  * (1.0d0 - exp(com_unit * kr_grid%xb(i) &
                                  * (xl_grid%xb(xl_grid%npts_b) - xl_grid%xb(xl_grid%npts_b-1))) &
                                  + com_unit * kr_grid%xb(i) * (xl_grid%max_val - xl_grid%xb(xl_grid%npts_b-1)))
            end do

        elseif(grid_spacing == 2) then
        ! non-equidistant, more points around rational surface
            do i=1, kr_grid%npts_b
                do n =2, xl_grid%npts_b-2
                    !varphi_lkr(i,n) = FT_hat_function_ne(rb(n), rb(n-1), rb(n+1), &
                    !                                    rb(n+2), kr_grid%xb(i))
                    varphi_lkr(i,n) = tilde_varphi_lkr(xl_grid%xb(n), xl_grid%xb(n-1), xl_grid%xb(n+1), kr_grid%xb(i))
                end do
                varphi_lkr(i,1) = 1.0d0 / ((xl_grid%xb(2) - xl_grid%xb(1)) * kr_grid%xb(i)**2) * &
                                  (1.0d0 - exp(-com_unit * kr_grid%xb(i) &
                                  * (xl_grid%xb(2) - xl_grid%xb(1))) - com_unit &
                                  * kr_grid%xb(i) * (xl_grid%xb(2) - xl_grid%xb(1)))
                varphi_lkr(i,xl_grid%npts_b) = 1.0d0 / ((xl_grid%xb(xl_grid%npts_b) - xl_grid%xb(xl_grid%npts_b-1)) &
                                  * kr_grid%xb(i)**2)  * &
                                  (1.0d0 - exp(com_unit * kr_grid%xb(i) * (xl_grid%xb(xl_grid%npts_b) &
                                  - xl_grid%xb(xl_grid%npts_b-1))) &
                                  + com_unit * kr_grid%xb(i) * (xl_grid%xb(xl_grid%npts_b) - xl_grid%xb(xl_grid%npts_b-1)))
            end do
        end if
    end if

    if (write_out) call write_FT_varphi

    contains

    ! varphi_{l,k_r} for non-equidistant grid
    double complex function FT_hat_function_ne(xl, xlm1, xlp1, xlp2, krr) result (res)

        implicit none
        double precision, intent(in) :: xl, xlm1, xlp1, xlp2, krr

        if (krr == 0.0d0) then ! analytical limit
            res = 0.5d0 * (xlp1 - xlm1)
        else 
            res = exp(- com_unit * krr * xl) * ((xlp1-xl) * (1.0d0 - exp(com_unit * krr * (xl-xlm1))) &
                  + (xl-xlm1) * (1.0d0 - exp(-com_unit * krr * (xlp1 - xl)))) &
                  / ((xlp1 - xl) * (xl - xlm1) * krr**2d0)
        end if

    end function

    ! \tilde varphi_{l,k_r} for non-equidistant grid, i.e. without the factor exp(-i kr xl)
    double complex function tilde_varphi_lkr(xl, xlm1, xlp1, krr) result (res)

        implicit none
        double precision, intent(in) :: xl, xlm1, xlp1, krr

        if (krr == 0.0d0) then ! analytical limit
            res = 0.5d0 * (xlp1 - xlm1)
        else 
            res = ((xlp1-xl) * (1.0d0 - exp(com_unit * krr * (xl-xlm1))) &
                  + (xl-xlm1) * (1.0d0 - exp(-com_unit * krr * (xlp1 - xl)))) &
                  / ((xlp1 - xl) * (xl - xlm1) * krr**2d0)
        end if

    end function


    ! varphi_{l,k_r} for equidistant grid
    double complex function FT_hat_function_e(xl, xlp1, krr) result (res)

        implicit none
        double precision, intent(in) :: xl, xlp1, krr

        if (krr == 0.0d0) then ! analytical limit
            res = 0.0d0
        else
            res = 2d0 * exp(-com_unit * krr * xl) * (1d0 - cos((xlp1 - xl) * krr)) / ((xlp1 - xl) * krr**2)
        end if

    end function

    ! tilde varphi_{l,k_r} for equidistant grid
    double complex function tilde_varphi_lkr_e(xl, xlp1, krr) result (res)

        implicit none
        double precision, intent(in) :: xl, xlp1, krr

        if (krr == 0.0d0) then ! analytical limit
            res = 0.0d0
        else
            res = 2d0 * (1d0 - cos((xlp1 - xl) * krr)) / ((xlp1 - xl) * krr**2)
        end if

    end function

    subroutine write_FT_varphi

        implicit none
        logical :: ex

        inquire(file=trim(output_path)//'basis_transform', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'basis_transform')
        end if
        open(unit = 77, file=trim(output_path)//'basis_transform/varphi_re.dat')
        open(unit = 78, file=trim(output_path)//'basis_transform/varphi_im.dat')
        do i=1, kr_grid%npts_b
            do n=1, xl_grid%npts_b
                write(77,*) real(varphi_lkr(i,n))
                write(78,*) dimag(varphi_lkr(i,n))
            end do
        end do
        close(77)
        close(78)

    end subroutine

end subroutine



    ! \tilde varphi_{l,k_r} for non-equidistant grid, i.e. without the factor exp(-i kr xl)
    !double complex function tilde_varphi_lkr(xl, xlm1, xlp1, krr) result (res)

    !    use constants, only: com_unit        

    !    implicit none
    !    double precision, intent(in) :: xl, xlm1, xlp1, krr
!
!        if (krr == 0.0d0) then ! analytical limit
!            res = 0.5d0 * (xlp1 - xlm1)
!        else 
!            res = ((xlp1-xl) * (1.0d0 - exp(com_unit * krr * (xl-xlm1))) &
!                  + (xl-xlm1) * (1.0d0 - exp(-com_unit * krr * (xlp1 - xl)))) &
!                  / ((xlp1 - xl) * (xl - xlm1) * krr**2d0)
!        end if
!
!    end function