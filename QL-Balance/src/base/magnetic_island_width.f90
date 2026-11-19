module magnetic_island_width_m

    implicit none

    private
    public :: inverse_q_derivative

contains

    !> d(1/q)/dr at every point of the grid r that carries q, from Lagrange stencils on that
    !> same grid (shifted at the ends to stay inside it).
    subroutine inverse_q_derivative(r, q, diotadr)
        use PolyLagrangeInterpolation, only: nlagr, plag_coeff, get_lagr_stencil

        double precision, intent(in) :: r(:), q(:)
        double precision, intent(out) :: diotadr(:)

        double precision :: stencil_coef(0:1, nlagr)
        integer :: ipoi, ibeg, iend

        do ipoi = 1, size(r)
            call get_lagr_stencil(ipoi, size(r), ibeg, iend)
            call plag_coeff(nlagr, 1, r(ipoi), r(ibeg:iend), stencil_coef)
            diotadr(ipoi) = sum(stencil_coef(1, :)/q(ibeg:iend))
        end do
    end subroutine inverse_q_derivative

end module magnetic_island_width_m

!> @brief subroutine magnetic_island_width Calculate the width of the magnetic island. Used to cut off the part of the perpendicular electric field coming from perturbed flux surfaces. This part would otherwise diverge at the resonant surface.
!> @author Markus Markl
!> @date 20.10.2022
subroutine magnetic_island_width(coef, nder, nlagr, ibeg, iend, mode, mi_width)

    use wave_code_data
    use h5mod
    use control_mod, only: equil_path, ihdf5IO
    use baseparam_mod, only: rsepar
    use magnetic_island_width_m, only: inverse_q_derivative

    implicit none
    double precision, intent(out) :: mi_width
    integer, intent(in) :: nder, nlagr, mode, ibeg, iend
    double precision, dimension(0:nder, nlagr), intent(in) :: coef

    double precision, dimension(:), allocatable :: diotadr
    complex(8), dimension(:), allocatable :: Delta_r
    double precision :: dummy, psi_tor_a
    complex(8) :: eye
    integer :: ipoi, iunit_equil
    character(1024) :: tempch

    allocate(diotadr(dim_r))
    allocate(Delta_r(dim_r))

    eye = (0.0, 1.0d0)

    ! read psi_tor at the separatrix. this is the last entry
    open(newunit=iunit_equil, file=trim(equil_path), status='old', action='read')
    ! skip first three text lines
    do ipoi = 1,3
        read(iunit_equil, *)
    end do
    do
        read(iunit_equil, *, end=1) dummy, dummy, dummy, psi_tor_a
    end do
1   continue
    close(iunit_equil)

    ! q, Br and r are on the wave-code grid (the balance boundary grid rb)
    call inverse_q_derivative(r, q, diotadr)

    !Delta_r = sqrt(- 4.0d0 * eye * r * Br * B0 / (diotadr * mode * psi_tor_a))
    Delta_r = abs(sqrt(- eye * Br * sqrt(antenna_factor) * rsepar**4 / (diotadr * mode * psi_tor_a * r)))

    write(*,*) "psi tor a", psi_tor_a
    write(*,*) "eye ", eye
    write(*,*) "mode ", mode
    write(*,*) "sqrt(ant fac) ", sqrt(antenna_factor)

    !open(77)
    !do ipoi = 1, npoib
    !    write(77,*) r(ipoi), real(Delta_r(ipoi)), imag(Delta_r(ipoi)), B0(ipoi), real(Br(ipoi)), imag(Br(ipoi))
    !end do
    !close(77)

    mi_width = sum(Delta_r(ibeg:iend) * coef(0,:)) * 2
    write(*,*) "magnetic island width for mode ", mode, " is ", mi_width

    deallocate(diotadr)
    deallocate(Delta_r)

    if (ihdf5IO .eq. 1) then
        write(*,*) "Write island width to hdf5"
        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        tempch = "/"//trim(h5_mode_groupname)
        write(*,*) trim(tempch)
        call h5_add_double_0(h5_id, trim(tempch)//"/mi_width", mi_width)
        CALL h5_close(h5_id)
        CALL h5_deinit()
    end if

end subroutine ! magnetic_island_width
