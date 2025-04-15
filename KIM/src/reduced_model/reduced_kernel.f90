module reduced_kernel

    implicit none

    contains

    subroutine fill_kernel_phi(kernel_phi_llp)

        use grid, only: xl_grid, rg_grid, kr_grid, krp_grid
        use functions, only: varphi_l
        use KIM_kinds, only: dp

        implicit none

        complex(dp), intent(out) :: kernel_phi_llp(:,:)
        integer :: j, l, lp

        kernel_phi_llp = 0.0d0

        do l = 1, xl_grid%npts_b
            do lp = 1, xl_grid%npts_b
                if (abs(l-lp).gt. 1) cycle
                do j = 1, rg_grid%npts_b
                    kernel_phi_llp(l, lp) = varphi_l(kr_grid%xb(j), xl_grid%xb(l-1), xl_grid%xb(l), xl_grid%xb(l+1))
                end do
            end do
        end do

    end subroutine fill_kernel_phi

end module