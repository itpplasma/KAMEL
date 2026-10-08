module fields_m

    use KIM_kinds_m, only: dp

    implicit none

    type :: EBdat_t
        real(dp), allocatable :: r_grid(:)
        real(dp) :: r_resonance = 0.0_dp
        real(dp) :: dx_asis = 0.0_dp
        real(dp) :: dx_transition = 0.0_dp
        complex(dp), allocatable :: Br(:)
        complex(dp), allocatable :: Bparallel(:)
        real(dp), allocatable :: D_ion(:,:,:)
        complex(dp), allocatable :: Apar(:)
        complex(dp), allocatable :: E_perp_psi(:)
        complex(dp), allocatable :: E_perp(:)
        complex(dp), allocatable :: E_perp_MA(:)
        complex(dp), allocatable :: Er(:)
        complex(dp), allocatable :: Etheta(:)
        complex(dp), allocatable :: Ez(:)
        complex(dp), allocatable :: Es(:) ! E "senkrecht", i.e. perpendicular to radial and parallel direction
        complex(dp), allocatable :: Ep(:) ! parallel to the equilibrium magnetic field
        complex(dp), allocatable :: jpar(:)   ! total parallel current density
        complex(dp), allocatable :: jpar_e(:) ! electron parallel current density
        complex(dp), allocatable :: jpar_i(:) ! ion parallel current density (sum over ion species)
        complex(dp), allocatable :: jrad(:)   ! radial current density
        complex(dp), allocatable :: rho(:), delta_n_e(:), delta_n_i(:, :)
        complex(dp), allocatable :: Phi(:)
        complex(dp), allocatable :: Phi_e(:)
        complex(dp), allocatable :: Phi_i(:)
        complex(dp), allocatable :: Phi_aligned(:)
        complex(dp), allocatable :: Phi_MA(:)
        complex(dp), allocatable :: Phi_MA_ideal(:)
        complex(dp), allocatable :: Phi_MA_asymptotic(:)
    end type

    type(EBdat_t) :: EBdat

    contains

    subroutine set_Br_field(EBdat_in, type_Br)

        use KIM_kinds_m, only: dp
        use grid_m, only: xl_grid
        use IO_collection_m, only: write_complex_profile
        use config_m, only: output_path

        implicit none

        type(EBdat_t) , intent(inout) :: EBdat_in
        integer, intent(in) :: type_Br
        complex(dp) :: Br_const
        character(len=256) :: file_path

        if (type_Br == 0) then ! Br constant
            Br_const = (1.0d0, 0.0d0)
            call set_Br_constant(EBdat_in, Br_const)
        else if (type_Br == 1) then ! Br from file
            file_path = './inp/Br_in.dat'
            call get_Br_from_txt(EBdat_in, file_path)
        end if

        call write_complex_profile(xl_grid%xb, EBdat_in%Br, xl_grid%npts_b, "/fields/br", &
            'Radial magnetic field perturbation', 'G')

    end subroutine

    subroutine set_Br_constant(EBdat_in, const)

        implicit none
        type(EBdat_t) , intent(inout) :: EBdat_in
        complex(dp), intent(in) :: const

        EBdat_in%Br = const

    end subroutine

    subroutine calculate_E_perp_psi(plasma_in, EBdat_in)

        use species_m, only: plasma_t
        use KIM_kinds_m, only: dp
        use equilibrium_m, only: B0

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

            call binsrc(plasma_in%r_grid, 1, plasma_in%grid_size, EBdat_in%r_grid(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(plasma_in%r_grid)) then
                iend = size(plasma_in%r_grid)
                ibeg = iend -nlagr + 1
            end if
            call plag_coeff(nlagr, nder, EBdat_in%r_grid(i), plasma_in%r_grid(ibeg:iend), coef)

            B0_int = sum(coef(0,:) * B0(ibeg:iend))

            EBdat_in%E_perp_psi(i) = Er_int * EBdat_in%Br(i) * ks_int / (B0_int * kp_int)
        end do

    end subroutine

    subroutine calculate_phi_aligned(plasma_in, EBdat_in)

        use species_m, only: plasma_t
        use KIM_kinds_m, only: dp
        use grid_m, only: xl_grid
        use constants_m, only: com_unit

        implicit none

        type(plasma_t) , intent(in) :: plasma_in
        type(EBdat_t) , intent(inout) :: EBdat_in
        integer :: i
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef

        real(dp) :: Er_int, kp_int
        real(dp) :: B0_int

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))
        if (.not. allocated(EBdat_in%phi_aligned)) allocate(EBdat_in%phi_aligned(xl_grid%npts_b))

        do i = 1, xl_grid%npts_b
            call binsrc(plasma_in%r_grid, 1, size(plasma_in%r_grid), xl_grid%xb(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(plasma_in%r_grid)) then
                iend = size(plasma_in%r_grid)
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, xl_grid%xb(i), plasma_in%r_grid(ibeg:iend), coef)

            Er_int = sum(coef(0,:) * plasma_in%Er(ibeg:iend))
            kp_int = sum(coef(0,:) * plasma_in%kp(ibeg:iend))
            B0_int = sum(coef(0,:) * plasma_in%B0(ibeg:iend))

            EBdat_in%phi_aligned(i) = com_unit * Er_int * EBdat_in%Br(i) / (B0_int * kp_int)
        end do

    end subroutine

    subroutine calculate_E_perp(EBdat_in)

        use constants_m, only: com_unit
        use species_m, only: plasma
        use KIM_kinds_m, only: dp

        implicit none

        type(EBdat_t) , intent(inout) :: EBdat_in
        integer :: i, ir, ibeg, iend
        integer :: nlagr = 4
        integer :: nder = 0
        real(dp), dimension(:,:), allocatable :: coef
        real(dp) :: ks

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        do i=1, size(EBdat_in%r_grid)
            call binsrc(plasma%r_grid, 1, size(plasma%r_grid), EBdat_in%r_grid(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(plasma%r_grid)) then
                iend = size(plasma%r_grid)
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, EBdat_in%r_grid(i), plasma%r_grid(ibeg:iend), coef)

            ks = sum(coef(0,:) * plasma%ks(ibeg:iend))

            EBdat_in%E_perp(i) = - com_unit * ks * EBdat_in%Phi(i)
        end do

    end subroutine

    subroutine calculate_MA_field(plasma_in, EBdat_in)

        use species_m, only: plasma_t
        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit

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
        if (.not. allocated(EBdat_in%E_perp_psi)) allocate(EBdat_in%E_perp_psi(size(EBdat_in%r_grid)))
        if (.not. allocated(EBdat_in%E_perp)) allocate(EBdat_in%E_perp(size(EBdat_in%r_grid)))
        if (.not. allocated(EBdat_in%E_perp_MA)) allocate(EBdat_in%E_perp_MA(size(EBdat_in%r_grid)))
        if (.not. allocated(EBdat_in%Phi_MA)) allocate(EBdat_in%Phi_MA(size(EBdat_in%r_grid)))

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
            B0_int = sum(coef(0,:) * plasma_in%B0(ibeg:iend))

            EBdat_in%E_perp_psi(i) = Er_int * EBdat_in%Br(i) * ks_int / (B0_int * kp_int)
            EBdat_in%E_perp(i) = - com_unit * ks_int * EBdat_in%Phi(i)
            EBdat_in%E_perp_MA(i) = EBdat_in%E_perp(i) + EBdat_in%E_perp_psi(i)
            EBdat_in%Phi_MA(i) = - EBdat_in%E_perp_MA(i) / (com_unit * ks_int)

        end do

    end subroutine

    subroutine get_Br_from_txt(EBdat_in, file_path)

        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit

        implicit none

        type(EBdat_t) , intent(inout) :: EBdat_in
        character(len=*), intent(in) :: file_path
        real(dp), allocatable :: Br_in_re(:), Br_in_im(:), r_in(:)
        integer :: i, ios, l
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef
        character(len=256) :: buffer

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        open(10, file=file_path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file: ", file_path
            stop
        end if

        l = 0
        do while (ios == 0)
            read(10, '(A)', iostat=ios) buffer
            if (ios == 0) then
                l = l + 1;
            end if
        end do

        allocate(Br_in_re(l), Br_in_im(l), r_in(l))

        rewind(10)

        do i = 1, l
            read(10, *) r_in(i), Br_in_re(i), Br_in_im(i)
        end do

        close(10)

        do i = 1, size(EBdat_in%r_grid)
            call binsrc(r_in, 1, size(r_in), EBdat_in%r_grid(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(r_in)) then
                iend = size(r_in)
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, EBdat_in%r_grid(i), r_in(ibeg:iend), coef)

            EBdat%Br(i) = sum(coef(0,:) * Br_in_re(ibeg:iend)) + com_unit * sum(coef(0,:) * Br_in_im(ibeg:iend))
        end do

    end subroutine

    subroutine calculate_E_from_phi(EBdat)

        use constants_m, only: com_unit
        use setup_m, only: m_mode, n_mode, R0
        use grid_m, only: xl_grid
        use functions_m, only: dvarphi_l_dx

        implicit none

        type(EBdat_t), intent(inout) :: EBdat
        integer :: i, j

        if (.not. allocated(EBdat%Er)) allocate(EBdat%Er(size(EBdat%r_grid)), EBdat%Etheta(size(EBdat%r_grid)), EBdat%Ez(size(EBdat%r_grid)))

        do i = 1, size(EBdat%r_grid)

            EBdat%Etheta(i) = - com_unit * m_mode * EBdat%Phi(i) / EBdat%r_grid(i)
            EBdat%Ez(i) = - com_unit * n_mode * EBdat%Phi(i) / R0
            EBdat%Er(i) = 0.0d0
            do j = 2, size(xl_grid%xb)-1
                EBdat%Er(i) = EBdat%Er(i) - dvarphi_l_dx(EBdat%r_grid(i), xl_grid%xb(j-1), xl_grid%xb(j), xl_grid%xb(j+1)) * EBdat%Phi(j)
            end do

        end do

    end subroutine

    subroutine calculate_E_in_rsp_from_cyl(EBdat)

        use equilibrium_m, only: hz, hth, equil_grid

        implicit none

        type(EBdat_t), intent(inout) :: EBdat
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef
        real(dp) :: hz_int, hth_int
        integer :: i

        if (.not. allocated(EBdat%Es)) allocate(EBdat%Es(size(EBdat%r_grid)), EBdat%Ep(size(EBdat%r_grid)))
        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        do i = 1, size(EBdat%r_grid)
            call binsrc(equil_grid, 1, size(equil_grid), EBdat%r_grid(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(equil_grid)) then
                iend = size(equil_grid)
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, EBdat%r_grid(i), equil_grid(ibeg:iend), coef)
            hz_int = sum(coef(0, :) * hz(ibeg:iend))
            hth_int = sum(coef(0, :) * hth(ibeg:iend))

            EBdat%Es(i) = hz_int * EBdat%Etheta(i) - hth_int * EBdat%Ez(i)
            EBdat%Ep(i) = hth_int * EBdat%Etheta(i) + hz_int * EBdat%Ez(i)

        end do

    end subroutine


    subroutine postprocess_electric_field(EBdat)

        use IO_collection_m, only: write_complex_profile_abs
        use grid_m, only: xl_grid
        use config_m, only: output_path, collision_model
        use species_m, only: plasma

        implicit none

        type(EBdat_t), intent(inout) :: EBdat
        character(len=50) :: suffix

        ! Determine suffix based on collision model
        select case (trim(collision_model))
        case ("Krook")
            suffix = "krook"
        case ("FokkerPlanck")
            suffix = "fp"
        case default
            suffix = "sol"
        end select

        call calculate_MA_field(plasma, EBdat)
        call write_complex_profile_abs(xl_grid%xb, EBdat%E_perp_psi, xl_grid%npts_b, "/fields/E_perp_psi", &
            'Field from perturbed magnetic flux surfaces', 'statV/cm')
        call write_complex_profile_abs(xl_grid%xb, EBdat%E_perp, xl_grid%npts_b, "/fields/E_perp", &
            'Field from potential perturbation', 'statV/cm')
        call write_complex_profile_abs(xl_grid%xb, EBdat%E_perp_MA, xl_grid%npts_b, "/fields/E_perp_MA", &
            'Total perpendicular misalignment electric field', 'statV/cm')
        call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_MA, xl_grid%npts_b, "/fields/Phi_MA", &
            'Misalignment electrostatic potential from total perpendicular misalignment electric field', 'statV')

        call calculate_E_from_phi(EBdat)
        call calculate_E_in_rsp_from_cyl(EBdat)

        call write_complex_profile_abs(EBdat%r_grid, EBdat%Er, size(EBdat%r_grid), "/fields/Er", &
            'Radial electric field perturbation in cylindrical coordinates', 'statV/cm')
        call write_complex_profile_abs(EBdat%r_grid, EBdat%Etheta, size(EBdat%r_grid), "/fields/Etheta", &
            'Poloidal electric field perturbation in cylindrical coordinates', 'statV/cm')
        call write_complex_profile_abs(EBdat%r_grid, EBdat%Ez, size(EBdat%r_grid), "/fields/Ez", &
            'Axial electric field perturbation in cylindrical coordinates', 'statV/cm')

        call write_complex_profile_abs(EBdat%r_grid, EBdat%Es, size(EBdat%r_grid), "/fields/Es", &
            'Electric field component perpendicular to radial and parallel direction in rsp coordinates', 'statV/cm')
        call write_complex_profile_abs(EBdat%r_grid, EBdat%Ep, size(EBdat%r_grid), "/fields/Ep", &
            'Electric field component parallel to equilibrium magnetic field in rsp coordinates', 'statV/cm')

        if (allocated(EBdat%Phi_aligned)) then
            call write_complex_profile_abs(xl_grid%xb, EBdat%Phi_aligned, xl_grid%npts_b, "/fields/Phi_aligned", &
                'Electrostatic potential aligned with magnetic field perturbation', 'statV')
        end if

    end subroutine


    subroutine postprocess_electric_field_no_write(EBdat)
        !! Same as postprocess_electric_field but without HDF5 writes.
        !! Used when hdf5_output is disabled (e.g. QL-Balance adapter mode).
        use species_m, only: plasma

        implicit none

        type(EBdat_t), intent(inout) :: EBdat

        call calculate_MA_field(plasma, EBdat)
        call calculate_E_from_phi(EBdat)
        call calculate_E_in_rsp_from_cyl(EBdat)

    end subroutine

    subroutine calculate_charge_density(rho, EBdat)

        use KIM_kinds_m, only: dp
        use species_m, only: plasma
        use constants_m, only: pi, com_unit

        implicit none

        type(EBdat_t), intent(in) :: EBdat
        complex(dp), allocatable, intent(out) :: rho(:)
        complex(dp) :: phi, Br
        integer :: sp

        integer :: i
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef


        if (.not. allocated(rho)) allocate(rho(size(plasma%r_grid)))
        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        rho = 0.0d0

        do i = 1, size(plasma%r_grid)

            call binsrc(EBdat%r_grid, 1, size(EBdat%r_grid), plasma%r_grid(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(EBdat%r_grid)) then
                iend = size(EBdat%r_grid)
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, plasma%r_grid(i), EBdat%r_grid(ibeg:iend), coef)
            phi = sum(coef(0, :) * EBdat%Phi(ibeg:iend))
            Br = sum(coef(0, :) * EBdat%Br(ibeg:iend))

            do sp=0, plasma%n_species-1
                rho(i) = plasma%spec(sp)%lambda_D(i)**-2.0d0
            end do

            rho(i) = rho(i) * (- 1.0d0) / (4.0d0 * pi) * (phi + com_unit * Br * plasma%Er(i) /(plasma%kp(i) * plasma%B0(i)))
        end do

    end subroutine

    subroutine calculate_current_density(jpar, Phi, Br, K_j_phi, K_j_B)

        use KIM_kinds_m, only: dp
        use weak_moments_m, only: weak_moment_projector_t
        use grid_m, only: xl_grid

        implicit none

        complex(dp), allocatable, intent(out) :: jpar(:)
        complex(dp), intent(in) :: Phi(:), Br(:)
        complex(dp), intent(in) :: K_j_phi(:,:)
        complex(dp), intent(in) :: K_j_B(:,:)
        type(weak_moment_projector_t) :: projector
        complex(dp), allocatable :: weak_load(:)

        weak_load = matmul(K_j_phi, Phi) + matmul(K_j_B, Br)
        call projector%init(xl_grid%xb)
        call projector%project(weak_load, jpar)

    end subroutine

    subroutine postprocess_species_density(EBdat, kernel_phi, kernel_B)
        use grid_m, only: xl_grid
        use kernel_m, only: kernel_spl_t
        use weak_moments_m, only: weak_moment_projector_t
        use species_m, only: plasma
        use constants_m, only: e_charge
        use IO_collection_m, only: write_complex_profile_abs
        implicit none
        type(EBdat_t), intent(inout) :: EBdat
        type(kernel_spl_t), intent(in) :: kernel_phi, kernel_B
        type(weak_moment_projector_t) :: projector
        complex(dp), allocatable :: rho_species(:), weak_load(:)
        integer :: sp, n
        character(len=256) :: dataset
        character(len=32) :: species_index

        n = size(xl_grid%xb)
        if (size(EBdat%r_grid) /= n) error stop 'Density field grid size mismatch'
        if (any(EBdat%r_grid /= xl_grid%xb)) &
            error stop 'Density requires the charge kernel field grid'
   if (plasma%spec(0)%Zspec >= 0) error stop 'Electron density requires negative charge'
        call projector%init(EBdat%r_grid)
        weak_load = matmul(kernel_phi%Kllp_e, EBdat%Phi) + &
                    matmul(kernel_B%Kllp_e, EBdat%Br)
        call projector%project(weak_load, rho_species)
        EBdat%rho = rho_species
        EBdat%delta_n_e = rho_species/(plasma%spec(0)%Zspec*e_charge)
        if (allocated(EBdat%delta_n_i)) deallocate (EBdat%delta_n_i)
        allocate (EBdat%delta_n_i(n, plasma%n_species - 1))
        call write_complex_profile_abs(EBdat%r_grid, EBdat%delta_n_e, n, &
            '/fields/delta_n_e', 'Electron number density from charge kernels', 'cm^-3')
        do sp = 1, plasma%n_species - 1
           if (plasma%spec(sp)%Zspec == 0) error stop 'Density requires charged species'
            weak_load = matmul(kernel_phi%Kllp_i(:, :, sp), EBdat%Phi) + &
                        matmul(kernel_B%Kllp_i(:, :, sp), EBdat%Br)
            call projector%project(weak_load, rho_species)
            EBdat%rho = EBdat%rho + rho_species
            EBdat%delta_n_i(:, sp) = rho_species/(plasma%spec(sp)%Zspec*e_charge)
            dataset = '/fields/delta_n_'//trim(plasma%spec(sp)%name)
            if (plasma%n_species > 2) then
                write (species_index, '(I0)') sp
                dataset = trim(dataset)//'_'//trim(species_index)
            end if
            call write_complex_profile_abs(EBdat%r_grid, EBdat%delta_n_i(:, sp), n, &
                       trim(dataset), 'Ion number density from charge kernels', 'cm^-3')
        end do
        call write_complex_profile_abs(EBdat%r_grid, EBdat%rho, n, '/fields/rho', &
                               'Total charge density from charge kernels', 'statC/cm^3')
    end subroutine postprocess_species_density

    subroutine calc_ideal_MA_phi(EBdat, kernel_phi, kernel_B)
        ! calcualte the misalignment Phi for E_perp_MA = 0, i.e. the ideal cancellation case
        ! where the flux surface corrugated phi cancels the potential surface corrugated phi
        ! this is used to check the second order derivative from the Laplace operator

        use kernel_m, only: kernel_spl_t

        implicit none

        type(EBdat_t), intent(inout) :: EBdat
        type(kernel_spl_t), intent(in) :: kernel_phi
        type(kernel_spl_t), intent(in) :: kernel_B

        integer :: l

        allocate(EBdat%Phi_MA_ideal(size(EBdat%Br)))

        do l = 1, size(EBdat%r_grid)
            EBdat%Phi_MA_ideal(l) = kernel_B%Kllp(l,l) * EBdat%Br(l) / kernel_phi%Kllp(l,l)
        end do

    end subroutine

end module
