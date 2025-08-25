module flr2_asymptotics_m

    implicit none


    contains

    subroutine calc_flr2_asymptotic_Phi_MA(plasma_in, EBdat)

        use KIM_kinds, only: dp
        use species_mod, only: plasma_t
        use fields_mod, only: EBdat_t
        use constants, only: pi, com_unit, sol
        use grid, only: xl_grid

        implicit none

        type(plasma_t), intent(in) :: plasma_in
        type(EBdat_t), intent(inout) :: EBdat

        integer :: sp
        complex(dp), allocatable :: F0(:), F2(:), H(:)

        allocate(F0(size(EBdat%Br)))
        allocate(F2(size(EBdat%Br)))
        allocate(H(size(EBdat%Br)))

        F0 = 0.0d0
        H = 0.0d0

        do j = 1, size(EBdat%Br)
            do sp = 0, plasma_in%Nspec-1
                F0(j) = F0(j) + 1.0d0/(plasma_in%spec(sp)%lambda_D(j)**2.0d0)
                H(j) = H(j) + plasma%spec(sp)%rho_L(j)**2.0d0 / (plasma_in%spec(sp)%lambda_D(j)**2.0d0) &
                    * (1.0d0 - (plasma_in%spec(sp)%dndr(j) * plasma_in%spec(sp)%T(j)&
                        + plasma_in%spec(sp)%n(j) * plasma_in%spec(sp)%dTdr(j)) / &
                        (plasma_in%spec(sp)%Zspec * e_charge * plasma_in%spec(sp)%n(j) &
                        * plasma_in%Er(j)) )
            end do
            F2(j) = - com_unit * plasma_in%om_E(j) * plasma_in%Er(j) / plasma_in%kp(j)**2.0d0 * H(j)
            F0(j) = F0(j) * (- com_unit) * (sol * plasma_in%om_E(j) * plasma_in%Er(j) / (4.0d0 * pi *plasma_in%kp(j)**2.0d0))
            H(j) = H(j) / (8.0d0 * pi)
        end do

    end subroutine

end module