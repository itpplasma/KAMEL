module fields

    use KIM_kinds, only: dp

    implicit none

    type :: EBdat_t
        real(dp), allocatable :: r_grid(:)
        complex(dp), allocatable :: Br(:)
        complex(dp), allocatable :: E_perp_psi(:)
        complex(dp), allocatable :: E_perp(:)
        complex(dp), allocatable :: Phi(:)
    end type

    type(EBdat_t) :: EBdat

    contains

    subroutine set_Br_constant(EBdat_in, const)

        implicit none
        type(EBdat_t) , intent(inout) :: EBdat_in
        complex(dp), intent(in) :: const

        EBdat_in%Br = const

    end subroutine

    subroutine calculate_E_perp_psi(plasma_in, EBdat_in)

        use species, only: plasma_t
        use KIM_kinds, only: dp
        use equilibrium, only: B0
        use plasma_parameter, only: iprof_length, r_prof

        implicit none

        type(plasma_t) , intent(in) :: plasma_in
        type(EBdat_t) , intent(inout) :: EBdat_in
        integer :: i
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef

        real(dp) :: Er_int, ks_int, kp_int
        real(dp) :: B0_int

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        do i = 1, size(EBdat_in%r_grid)
            call binsrc(plasma_in%r_grid, 1, size(plasma_in%r_grid), EBdat_in%r_grid(i), ir) 
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(plasma_in%r_grid)) then
                iend = size(plasma_in%r_grid)
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, EBdat_in%r_grid(i), plasma_in%r_grid(ibeg:iend), coef)

            Er_int = sum(coef(0,:) * plasma_in%Er(ibeg:iend))
            ks_int = sum(coef(0,:) * plasma_in%ks(ibeg:iend))
            kp_int = sum(coef(0,:) * plasma_in%kp(ibeg:iend))

            call binsrc(r_prof, 1, iprof_length, EBdat_in%r_grid(i), ir) 
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(plasma_in%r_grid)) then
                iend = size(plasma_in%r_grid)
                ibeg = iend -nlagr + 1
            end if
            call plag_coeff(nlagr, nder, EBdat_in%r_grid(i), plasma_in%r_grid(ibeg:iend), coef)

            B0_int = sum(coef(0,:) * B0(ibeg:iend))

            EBdat_in%E_perp_psi(i) = - Er_int * EBdat_in%Br(i) * ks_int / (B0_int * kp_int)
        end do

    end subroutine

    subroutine calculate_E_perp(EBdat_in)

        implicit none

        type(EBdat_t) , intent(inout) :: EBdat_in
        integer :: i

        do i=1, size(EBdat_in%r_grid)-1
            EBdat_in%E_perp(i) = -(EBdat_in%Phi(i+1) - EBdat_in%Phi(i))/(EBdat_in%r_grid(i+1) - EBdat_in%r_grid(i))
        end do

    end subroutine

end module