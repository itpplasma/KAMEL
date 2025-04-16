module reduced_kernel

    use KIM_kinds, only: dp

    implicit none

    type :: kernel_spl_t
        integer :: npts_l, npts_lp
        complex(dp), allocatable :: Kllp(:,:)
        contains
            procedure :: init_kernel
    end type kernel_spl_t

    contains

    subroutine init_kernel(this, npts_l, npts_lp)

        implicit none

        class(kernel_spl_t), intent(inout) :: this
        integer, intent(in) :: npts_l, npts_lp

        this%npts_l = npts_l
        this%npts_lp = npts_lp
        allocate(this%Kllp(npts_l, npts_lp))
        this%Kllp = 0.0d0

    end subroutine init_kernel

    subroutine fill_kernel_phi(kernel_phi_llp)

        use grid, only: rg_grid, kr_grid, krp_grid
        use functions, only: varphi_l
        use KIM_kinds, only: dp

        implicit none

        type(kernel_spl_t), intent(inout) :: kernel_phi_llp

        integer :: j, l, lp

        do l = 1, kernel_phi_llp%npts_l
            do lp = 1, kernel_phi_llp%npts_lp
                if (abs(l-lp).gt. 1) cycle
                do j = 1, rg_grid%npts_b
                    kernel_phi_llp%Kllp(l, lp) = l + lp !varphi_l(kr_grid%xb(j), xl_grid%xb(l-1), xl_grid%xb(l), xl_grid%xb(l+1))
                end do
            end do
        end do

    end subroutine fill_kernel_phi

end module