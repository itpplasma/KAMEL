module fft_kernel

    use ISO_C_BINDING
    use integrands
    use grid, only: npoib, kr, varphi_lkr, k_space_dim, krp, r_space_dim
    use config, only: fdebug

    implicit none

    include 'fftw3.f03'

    integer, parameter ::  N = 10
    integer(8) :: plan

    double complex, dimension(:,:,:), allocatable :: kernel_rho_for_fft ! kernel for inverse fft input
    double complex, dimension(:,:,:), allocatable :: kernel_rho_llp ! kernel for inverse fft output, i.e. in the spline space (and r_g space)

    contains 

    subroutine fill_input_rho(i_l, i_lp)
    
        implicit none
        integer, intent(in) :: i_l, i_lp
        integer :: i_kr, i_krp, i_rg

        if(.not. allocated(kernel_rho_for_fft)) allocate(kernel_rho_for_fft(k_space_dim, k_space_dim, r_space_dim))

        do i_kr= 1, k_space_dim
            do i_krp = 1, k_space_dim
                do i_rg = 1, r_space_dim
                    ! here comes calculation of integral kernel elements
                    !kernel_rho_for_fft(i_kr, i_krp, i_rg) = !
                end do
            end do
        end do


    end subroutine

    subroutine fft_kernel_rho

        use finufft_mod

        implicit none

        integer :: l, lp
        integer :: ier,iflag,ntrans,type,dim
        integer :: M,N1,N2,N3,Nk
        integer :: plan,n_modes(3)
        double precision, allocatable :: xj(:),yj(:),zj(:), sk(:),tk(:),uk(:)
        double precision :: tol
        double complex, allocatable :: cj(:), fk(:)
        type(finufft_opts) opts

        if(.not. allocated(kernel_rho_for_fft)) allocate(kernel_rho_for_fft(k_space_dim, k_space_dim, r_space_dim))

        ! allocate data

        if (fdebug == 1) write(*,*) "Debug: fft of kernel"

        do l = 1, l_space_dim
            do lp = 1, l_space_dim
                !call fill_input_rho(l, lp, kernel_in)
                !call dfftw_plan_dft_2d(plan, k_space_dim, k_space_dim, kernel_rho_for_fft, kernel_rho_for_fft, FFTW_BACKWARD, FFTW_ESTIMATE)
                !call dfftw_execute_dft(plan, inp, outp)
                !kernel_rho_llp(l, lp,:) = kernel_rho_for_fft
                call finufft2d3(M,xj,yj,cj,iflag,tol,Nk,sk,tk,fk,opts,ier)
            end do
        end do

        call dfftw_destroy_plan(plan)


    end subroutine

end module