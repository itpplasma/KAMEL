module host_geometry_m
    !> Flux-surface geometry supplied by a host code that solves the toroidal
    !> problem (e.g. MEPHIT). It replaces KIM's screw-pinch field magnitude and
    !> wavenumbers and adds the radial metric of the perpendicular operators.
    !>
    !> The host maps its flux label to the radial coordinate r of the KIM
    !> profiles and evaluates, on each flux surface,
    !>
    !>   B0             field magnitude used for omega_c and rho_L [G]
    !>   kp_m, kp_n     k_parallel = m*kp_m + n*kp_n [1/cm]
    !>   ks_m, ks_n     k_s        = m*ks_m + n*ks_n [1/cm]
    !>   radial_metric  M = <|grad r|^2> with the weight of the host's
    !>                  perpendicular operator [1]
    !>
    !> for the KIM mode numbers (m, n) = (m_mode, n_mode). The physical radial
    !> wavenumber is sqrt(M) k_r, so M multiplies k_r^2 in the Larmor-radius
    !> arguments, and the radial Poisson operator becomes d/dr (M d/dr).
    !> The E x B frequency follows from Er, k_s and B0 as in the cylinder.
    !>
    !> File format: one surface per line, columns
    !>   r  B0  kp_m  kp_n  ks_m  ks_n  radial_metric
    !> with r strictly increasing; lines starting with '#' are comments.
    !> The table should cover the radial range of the KIM profiles. Outside it
    !> the geometry is held constant, as KIM holds its profiles and periodic
    !> background; a periodic window reaching past the plasma edge needs this.

    use KIM_kinds_m, only: dp

    implicit none
    private

    public :: configure_host_geometry, reset_host_geometry
    public :: host_geometry_active, host_geometry_has_metric
    public :: load_host_geometry, evaluate_host_geometry, apply_host_geometry

    integer, parameter :: n_columns = 7

    logical, save :: active = .false.
    logical, save :: nonunit_metric = .false.
    real(dp), allocatable, save :: table(:, :)  ! (n_surfaces, n_columns)

contains

    !> Load geometry_file from the KIM configuration, or reset when it is empty.
    subroutine configure_host_geometry()
        use config_m, only: geometry_file, profile_location

        character(len=:), allocatable :: path

        call reset_host_geometry()
        if (len_trim(geometry_file) == 0) return
        if (geometry_file(1:1) == '/') then
            path = trim(geometry_file)
        else
            path = trim(profile_location) // trim(geometry_file)
        end if
        call load_host_geometry(path)
    end subroutine configure_host_geometry

    subroutine reset_host_geometry()
        active = .false.
        nonunit_metric = .false.
        if (allocated(table)) deallocate(table)
    end subroutine reset_host_geometry

    logical function host_geometry_active()
        host_geometry_active = active
    end function host_geometry_active

    !> True when the radial metric differs from 1 somewhere in the table.
    logical function host_geometry_has_metric()
        host_geometry_has_metric = nonunit_metric
    end function host_geometry_has_metric

    subroutine load_host_geometry(path)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use logger_m, only: log_info

        character(*), intent(in) :: path

        character(len=1024) :: line
        real(dp) :: row(n_columns)
        integer :: unit, ios, n_rows, i
        logical :: exists

        call reset_host_geometry()
        inquire(file=path, exist=exists)
        if (.not. exists) then
            print *, 'host geometry file not found: ', trim(path)
            error stop 'host geometry file not found'
        end if

        n_rows = 0
        open(newunit=unit, file=path, status='old', action='read')
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (is_data_line(line)) n_rows = n_rows + 1
        end do
        if (n_rows < 4) error stop 'host geometry file needs at least 4 surfaces'

        allocate(table(n_rows, n_columns))
        rewind(unit)
        i = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (.not. is_data_line(line)) cycle
            read(line, *, iostat=ios) row
            if (ios /= 0) then
                print *, 'host geometry: cannot read 7 columns from: ', trim(line)
                error stop 'malformed host geometry file'
            end if
            i = i + 1
            table(i, :) = row
        end do
        close(unit)

        if (.not. all(ieee_is_finite(table))) error stop 'host geometry contains non-finite values'
        if (any(table(2:, 1) <= table(:n_rows - 1, 1))) then
            error stop 'host geometry radii must increase strictly'
        end if
        if (any(table(:, 2) <= 0.0_dp)) error stop 'host geometry B0 must be positive'
        if (any(table(:, 7) <= 0.0_dp)) error stop 'host geometry radial metric must be positive'

        active = .true.
        nonunit_metric = any(table(:, 7) /= 1.0_dp)
        call log_info('Using host flux-surface geometry from ' // trim(path))
    end subroutine load_host_geometry

    logical function is_data_line(line)
        character(*), intent(in) :: line
        character(len=len(line)) :: stripped

        stripped = adjustl(line)
        is_data_line = len_trim(stripped) > 0 .and. stripped(1:1) /= '#'
    end function is_data_line

    !> Host geometry at radius r for mode numbers (m, n), held constant beyond
    !> the ends of the table.
    subroutine evaluate_host_geometry(r, m, n, B0, kp, ks, radial_metric)
        real(dp), intent(in) :: r
        integer, intent(in) :: m, n
        real(dp), intent(out) :: B0, kp, ks, radial_metric

        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr), values(n_columns), r_lo, r_hi, rc
        integer :: n_rows, ir, ibeg, iend, c

        if (.not. active) error stop 'evaluate_host_geometry called without host geometry'
        n_rows = size(table, 1)
        r_lo = table(1, 1)
        r_hi = table(n_rows, 1)
        rc = min(max(r, r_lo), r_hi)

        call binsrc(table(:, 1), 1, n_rows, rc, ir)
        ibeg = max(1, ir - nlagr / 2)
        iend = ibeg + nlagr - 1
        if (iend > n_rows) then
            iend = n_rows
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, rc, table(ibeg:iend, 1), coef)
        do c = 2, n_columns
            values(c) = sum(coef(0, :) * table(ibeg:iend, c))
        end do

        B0 = values(2)
        kp = real(m, dp) * values(3) + real(n, dp) * values(4)
        ks = real(m, dp) * values(5) + real(n, dp) * values(6)
        radial_metric = values(7)
    end subroutine evaluate_host_geometry

    !> Overwrite B0, kp, ks and om_E of plasma_in on its radial grid with the
    !> host geometry and store the radial metric. om_E keeps KIM's definition
    !> -c Er ks / B0, so the E x B frequency and the thermodynamic forces use
    !> the same Er.
    subroutine apply_host_geometry(plasma_in, B0_equil)
        use species_m, only: plasma_t, reallocate
        use setup_m, only: m_mode, n_mode
        use wavenumber_geometry_m, only: exb_rotation_frequency

        type(plasma_t), intent(inout) :: plasma_in
        real(dp), intent(inout), optional :: B0_equil(:)

        integer :: i

        if (.not. active) return
        call reallocate(plasma_in%radial_metric, plasma_in%grid_size)
        do i = 1, plasma_in%grid_size
            call evaluate_host_geometry(plasma_in%r_grid(i), m_mode, n_mode, &
                plasma_in%B0(i), plasma_in%kp(i), plasma_in%ks(i), &
                plasma_in%radial_metric(i))
            plasma_in%om_E(i) = exb_rotation_frequency(plasma_in%ks(i), &
                plasma_in%Er(i), plasma_in%B0(i))
        end do
        if (present(B0_equil)) B0_equil(1:plasma_in%grid_size) = plasma_in%B0
    end subroutine apply_host_geometry

end module host_geometry_m
