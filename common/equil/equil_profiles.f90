!> @file equil_profiles.f90
!> @brief Equilibrium profile computation from geqdsk files
!>
!> Provides functionality to:
!> - Read axisymmetric equilibrium from geqdsk files (via libneo)
!> - Find the magnetic axis by iterative field-line tracing
!> - Trace flux surfaces to determine the separatrix
!> - Compute radial profiles: effective radius r, safety factor q,
!>   poloidal flux psi, toroidal flux phi, volume
!>
!> The profiles map between different flux coordinate representations,
!> enabling conversion from psi_pol (poloidal flux) to psi_tor (toroidal flux)
!> and effective radius.

module equil_profiles_m

    use rk4_integrator_m, only: rk4_step
    use field_line_rhs_m, only: set_field_line_mode, get_dz_dphi, field_line_rhs
    use field_eq_mod, only: btf, rtf, psif, nrad, nzet, rad, zet
    use field_sub, only: field_eq

    implicit none
    private

    double precision, parameter :: pi = 3.14159265358979d0

    !> Equilibrium profile data type
    type, public :: equil_profiles_t
        !> Magnetic axis position
        double precision :: raxis = 0.d0
        double precision :: zaxis = 0.d0

        !> Poloidal flux at axis
        double precision :: psi_axis = 0.d0

        !> Toroidal field and major radius at axis
        double precision :: btor = 0.d0
        double precision :: rbig = 0.d0

        !> Separatrix data
        double precision :: psi_sep = 0.d0
        integer :: nsurf = 0

        !> Number of radial grid points
        integer :: nsqpsi = 0

        !> Radial profiles (allocated with dimension nsqpsi)
        double precision, allocatable :: rsmall(:)      !< Effective minor radius
        double precision, allocatable :: qsaf(:)        !< Safety factor
        double precision, allocatable :: psisurf(:)     !< Poloidal flux (relative to axis)
        double precision, allocatable :: phitor(:)      !< Toroidal flux
        double precision, allocatable :: volume(:)      !< Enclosed volume
        double precision, allocatable :: dphidpsi(:)    !< d(phi_tor)/d(psi_pol)

        !> Normalized toroidal flux for interpolation
        double precision, allocatable :: phinorm(:)     !< phi_tor / phi_tor_max

        !> Computation box bounds
        double precision :: rmn = 0.d0, rmx = 0.d0
        double precision :: zmn = 0.d0, zmx = 0.d0

        !> Initialization flag
        logical :: initialized = .false.

    contains
        procedure :: init => equil_profiles_init
        procedure :: find_axis => equil_profiles_find_axis
        procedure :: compute_profiles => equil_profiles_compute
        procedure :: cleanup => equil_profiles_cleanup
        procedure :: write_output => equil_profiles_write_output
    end type equil_profiles_t

    public :: init_equilibrium_field

contains

    !> Initialize the equilibrium field from input file
    !> This must be called before any field-line tracing
    subroutine init_equilibrium_field()
        use field_sub

        double precision :: rrr, ppp, zzz
        double precision :: Brad, Bphi, Bzet
        double precision :: dBrdR, dBrdp, dBrdZ
        double precision :: dBpdR, dBpdp, dBpdZ
        double precision :: dBzdR, dBzdp, dBzdZ

        ! Initialize field by calling with dummy coordinates
        ! This triggers libneo to read the geqdsk file specified in field_divB0.inp
        rrr = 0.d0
        ppp = 0.d0
        zzz = 0.d0

        call field(rrr, ppp, zzz, Brad, Bphi, Bzet, &
                   dBrdR, dBrdp, dBrdZ, &
                   dBpdR, dBpdp, dBpdZ, &
                   dBzdR, dBzdp, dBzdZ)

    end subroutine init_equilibrium_field

    !> Initialize equilibrium profiles computation
    subroutine equil_profiles_init(self, nsqpsi_in)
        class(equil_profiles_t), intent(inout) :: self
        integer, intent(in) :: nsqpsi_in

        ! Clean up any existing allocations
        call self%cleanup()

        self%nsqpsi = nsqpsi_in

        ! Allocate profile arrays
        allocate(self%rsmall(nsqpsi_in))
        allocate(self%qsaf(nsqpsi_in))
        allocate(self%psisurf(nsqpsi_in))
        allocate(self%phitor(nsqpsi_in))
        allocate(self%volume(nsqpsi_in))
        allocate(self%dphidpsi(nsqpsi_in))
        allocate(self%phinorm(nsqpsi_in))

        ! Initialize arrays
        self%rsmall = 0.d0
        self%qsaf = 0.d0
        self%psisurf = 0.d0
        self%phitor = 0.d0
        self%volume = 0.d0
        self%dphidpsi = 0.d0
        self%phinorm = 0.d0

        ! Get computation box from equilibrium grid
        self%rmn = rad(1)
        self%rmx = rad(nrad)
        self%zmn = zet(1)
        self%zmx = zet(nzet)

        ! Store toroidal field info
        self%btor = btf
        self%rbig = rtf

        self%initialized = .true.

    end subroutine equil_profiles_init

    !> Find the magnetic axis by iterative field-line tracing
    subroutine equil_profiles_find_axis(self, nstep, nmap, niter)
        class(equil_profiles_t), intent(inout) :: self
        integer, intent(in) :: nstep   !< Integration steps per toroidal circuit
        integer, intent(in) :: nmap    !< Number of toroidal circuits per iteration
        integer, intent(in) :: niter   !< Number of Newton iterations

        integer :: ndim, iter, i, ntotstep
        double precision :: h, phi, phibeg
        double precision, allocatable :: y(:)
        double precision :: rbeg, zbeg
        double precision :: ppp
        double precision :: Brad, Bphi, Bzet
        double precision :: dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ

        external :: rhs1

        ndim = 5
        allocate(y(ndim))

        ! Initial guess: center of computation box
        rbeg = 0.5d0 * (self%rmn + self%rmx)
        zbeg = 0.5d0 * (self%zmn + self%zmx)

        ntotstep = nstep * nmap
        h = 2.d0 * pi / nstep

        ! Set mode for axis finding
        call set_field_line_mode(1)

        ! Initialize state
        phi = 0.d0
        y = 0.d0
        y(1) = rbeg
        y(2) = zbeg

        ! Iterative refinement of axis position
        do iter = 1, niter
            phibeg = phi
            y(4:5) = 0.d0

            ! Trace field lines for nmap toroidal circuits
            do i = 1, ntotstep
                call RK4D(y, ndim, phi, h, rhs1)
            end do

            ! Update axis estimate as average position
            y(1:2) = y(4:5) / (phi - phibeg)
        end do

        self%raxis = y(1)
        self%zaxis = y(2)

        ! Get psi at axis
        ppp = 0.d0
        call field_eq(self%raxis, ppp, self%zaxis, Brad, Bphi, Bzet, &
                      dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
        self%psi_axis = psif

        deallocate(y)

    end subroutine equil_profiles_find_axis

    !> Compute equilibrium profiles by tracing flux surfaces
    subroutine equil_profiles_compute(self, nstep, nsurfmax, niter)
        class(equil_profiles_t), intent(inout) :: self
        integer, intent(in) :: nstep     !< Integration steps per circuit
        integer, intent(in) :: nsurfmax  !< Number of trial surfaces for separatrix search
        integer, intent(in) :: niter     !< Newton iterations for surface closure

        integer :: ndim, isurf, iter, i, i1
        double precision :: h, phi, hr, sig, sig0, hh
        double precision, allocatable :: y(:)
        double precision, allocatable :: psisurf_scan(:), sqpsi(:), startrind(:)
        double precision :: ppp, rrr, zzz
        double precision :: Brad, Bphi, Bzet
        double precision :: dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ
        double precision :: aiota, phitor_max, psisurf_max, rsmall_max
        double precision :: sqpsimin, sqpsimax, rindmin, rindmax
        double precision :: dz_dphi

        external :: rhs1

        ndim = 5
        allocate(y(ndim))
        allocate(psisurf_scan(nsurfmax))

        h = 2.d0 * pi / nstep
        hr = (self%rmx - self%raxis) / nsurfmax

        ! Set mode for profile computation
        call set_field_line_mode(2)

        ! Scan surfaces to find separatrix
        self%nsurf = 0
        do isurf = 1, nsurfmax
            phi = 0.d0
            y(1) = self%raxis + hr * isurf
            y(2) = self%zaxis
            y(3:5) = 0.d0

            ! Trace until returning to midplane
            call RK4D(y, ndim, phi, h, rhs1)
            sig = y(2) - self%zaxis

            do while (sig * (y(2) - self%zaxis) > 0.d0)
                call RK4D(y, ndim, phi, h, rhs1)
                ! Check if we've left the computation box
                if (y(1) < self%rmn .or. y(1) > self%rmx .or. &
                    y(2) < self%zmn .or. y(2) > self%zmx) then
                    exit
                end if
            end do

            ! Check if surface closed properly
            if (y(1) < self%rmn .or. y(1) > self%rmx .or. &
                y(2) < self%zmn .or. y(2) > self%zmx) then
                self%nsurf = isurf - 1
                exit
            end if

            ! Continue to second half of circuit
            sig = y(2) - self%zaxis
            do while (sig * (y(2) - self%zaxis) > 0.d0)
                call RK4D(y, ndim, phi, h, rhs1)
                if (y(1) < self%rmn .or. y(1) > self%rmx .or. &
                    y(2) < self%zmn .or. y(2) > self%zmx) then
                    exit
                end if
            end do

            if (y(1) < self%rmn .or. y(1) > self%rmx .or. &
                y(2) < self%zmn .or. y(2) > self%zmx) then
                self%nsurf = isurf - 1
                exit
            end if

            ! Newton iteration to close surface exactly
            do iter = 1, niter
                dz_dphi = get_dz_dphi()
                if (abs(dz_dphi) > 1.d-12) then
                    hh = (self%zaxis - y(2)) / dz_dphi
                    call RK4D(y, ndim, phi, hh, rhs1)
                end if
            end do

            ! Get psi at this surface
            rrr = y(1)
            zzz = y(2)
            ppp = 0.d0
            call field_eq(rrr, ppp, zzz, Brad, Bphi, Bzet, &
                          dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

            psisurf_scan(isurf) = psif - self%psi_axis
            rsmall_max = sqrt(abs(y(3)) / pi)
            psisurf_max = psif - self%psi_axis
            phitor_max = y(4) / (2.d0 * pi)

            self%nsurf = isurf
        end do

        ! Remove last (potentially bad) surface
        if (self%nsurf > 1) self%nsurf = self%nsurf - 1

        self%psi_sep = psisurf_scan(self%nsurf) + self%psi_axis

        ! Re-interpolate to equidistant grid in sqrt(psi)
        allocate(sqpsi(0:self%nsurf))
        allocate(startrind(0:self%nsqpsi))

        sqpsi(0) = 0.d0
        do i = 1, self%nsurf
            sqpsi(i) = sqrt(abs(psisurf_scan(i)))
        end do

        rindmin = 0.d0
        rindmax = dble(self%nsurf)
        sqpsimin = 0.d0
        sqpsimax = sqpsi(self%nsurf)

        ! Compute profiles on equidistant sqrt(psi) grid
        do isurf = 1, self%nsqpsi
            phi = 0.d0
            ! Linear interpolation for starting radius
            startrind(isurf) = rindmax * dble(isurf) / dble(self%nsqpsi)

            y(1) = self%raxis + hr * startrind(isurf)
            y(2) = self%zaxis
            y(3:5) = 0.d0

            ! Trace surface
            call RK4D(y, ndim, phi, h, rhs1)
            sig = y(2) - self%zaxis

            do while (sig * (y(2) - self%zaxis) > 0.d0)
                call RK4D(y, ndim, phi, h, rhs1)
            end do

            sig = y(2) - self%zaxis
            do while (sig * (y(2) - self%zaxis) > 0.d0)
                call RK4D(y, ndim, phi, h, rhs1)
            end do

            ! Newton iteration
            do iter = 1, niter
                dz_dphi = get_dz_dphi()
                if (abs(dz_dphi) > 1.d-12) then
                    hh = (self%zaxis - y(2)) / dz_dphi
                    call RK4D(y, ndim, phi, hh, rhs1)
                end if
            end do

            ! Compute profile values
            aiota = 2.d0 * pi / phi

            rrr = y(1)
            zzz = y(2)
            ppp = 0.d0
            call field_eq(rrr, ppp, zzz, Brad, Bphi, Bzet, &
                          dBrdR, dBrdp, dBrdZ, dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)

            self%rsmall(isurf) = sqrt(abs(y(3)) / pi)
            self%qsaf(isurf) = 1.d0 / aiota
            self%psisurf(isurf) = psif - self%psi_axis
            self%phitor(isurf) = y(4) / (2.d0 * pi)
            self%volume(isurf) = abs(y(5)) * pi
        end do

        ! Compute derivatives and normalized flux
        phitor_max = self%phitor(self%nsqpsi)
        do i = 1, self%nsqpsi
            self%phinorm(i) = self%phitor(i) / phitor_max

            ! Central difference for dphidpsi
            i1 = min(self%nsqpsi - 1, max(2, i))
            self%dphidpsi(i) = (self%phitor(i1+1) - self%phitor(i1-1)) / &
                               (self%psisurf(i1+1) - self%psisurf(i1-1))
        end do

        deallocate(y, psisurf_scan, sqpsi, startrind)

    end subroutine equil_profiles_compute

    !> Write equilibrium profiles to output files
    subroutine equil_profiles_write_output(self, output_dir)
        class(equil_profiles_t), intent(in) :: self
        character(len=*), intent(in) :: output_dir

        integer :: i, iunit
        character(len=512) :: filename

        ! Write btor_rbig.dat
        filename = trim(output_dir) // '/btor_rbig.dat'
        open(newunit=iunit, file=filename, status='replace', action='write')
        write(iunit, *) self%btor, self%rbig
        close(iunit)

        ! Write axis.dat
        filename = trim(output_dir) // '/axis.dat'
        open(newunit=iunit, file=filename, status='replace', action='write')
        write(iunit, *) 'raxis = ', self%raxis
        write(iunit, *) 'zaxis = ', self%zaxis
        close(iunit)

        ! Write equil_r_q_psi.dat
        filename = trim(output_dir) // '/equil_r_q_psi.dat'
        open(newunit=iunit, file=filename, status='replace', action='write')
        write(iunit, *) '# nrad = ', self%nsqpsi
        write(iunit, *) '# rmax ', self%rsmall(self%nsqpsi), '   psimax ', &
                        self%psisurf(self%nsqpsi), '   phimax ', self%phitor(self%nsqpsi)
        write(iunit, *) '# radius r,  safety factor q, poloidal flux psi, toroidal flux phi, ', &
                        'd phi / d psi, geom. radius r, volume'

        do i = 1, self%nsqpsi
            write(iunit, *) sqrt(2.d0 * abs(self%phitor(i) / self%btor)), &
                            self%qsaf(i), self%psisurf(i), self%phitor(i), &
                            self%dphidpsi(i), self%rsmall(i), self%volume(i)
        end do
        close(iunit)

    end subroutine equil_profiles_write_output

    !> Clean up allocated memory
    subroutine equil_profiles_cleanup(self)
        class(equil_profiles_t), intent(inout) :: self

        if (allocated(self%rsmall)) deallocate(self%rsmall)
        if (allocated(self%qsaf)) deallocate(self%qsaf)
        if (allocated(self%psisurf)) deallocate(self%psisurf)
        if (allocated(self%phitor)) deallocate(self%phitor)
        if (allocated(self%volume)) deallocate(self%volume)
        if (allocated(self%dphidpsi)) deallocate(self%dphidpsi)
        if (allocated(self%phinorm)) deallocate(self%phinorm)

        self%initialized = .false.

    end subroutine equil_profiles_cleanup

end module equil_profiles_m
