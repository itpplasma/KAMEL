module periodic_electron_response_m
    ! Adapter for an independently computed electron response in KIM's basis.
    ! All blocks carry physical charge/current units and KIM's Fourier weights.
    use KIM_kinds_m, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: electron_response_t, replace_periodic_electrons
    type :: electron_response_t
        complex(dp), allocatable :: rho_phi(:, :), rho_b(:, :)
        complex(dp), allocatable :: j_phi(:, :), j_b(:, :)
    end type
contains
    subroutine replace_periodic_electrons(response, kphi, kb, kjphi, kjb, &
                                          kphi_sp, kb_sp, kjphi_sp, kjb_sp)
        type(electron_response_t), intent(in) :: response
        complex(dp), intent(inout) :: kphi(:, :), kb(:, :), kjphi(:, :), kjb(:, :)
        complex(dp), intent(inout) :: kphi_sp(:, :, 0:), kb_sp(:, :, 0:)
        complex(dp), intent(inout) :: kjphi_sp(:, :, 0:), kjb_sp(:, :, 0:)
        integer :: dim
        dim = size(kphi, 1)
        if (.not. allocated(response%rho_phi) .or. .not. allocated(response%rho_b) .or. &
            .not. allocated(response%j_phi) .or. .not. allocated(response%j_b)) &
            error stop 'external electron response requires all four blocks'
        call check_block(response%rho_phi, dim)
        call check_block(response%rho_b, dim)
        call check_block(response%j_phi, dim)
        call check_block(response%j_b, dim)
        kphi_sp(:, :, 0) = response%rho_phi
        kb_sp(:, :, 0) = response%rho_b
        kjphi_sp(:, :, 0) = response%j_phi
        kjb_sp(:, :, 0) = response%j_b
        ! Re-sum species rather than subtracting nearly cancelling old blocks.
        kphi = sum(kphi_sp, dim=3)
        kb = sum(kb_sp, dim=3)
        kjphi = sum(kjphi_sp, dim=3)
        kjb = sum(kjb_sp, dim=3)
    end subroutine

    subroutine check_block(block, dim)
        complex(dp), intent(in) :: block(:, :)
        integer, intent(in) :: dim
        if (any(shape(block) /= [dim, dim])) error stop 'external electron block shape mismatch'
        if (.not. all(ieee_is_finite(real(block))) .or. &
            .not. all(ieee_is_finite(aimag(block)))) error stop 'nonfinite external electron block'
    end subroutine
end module
