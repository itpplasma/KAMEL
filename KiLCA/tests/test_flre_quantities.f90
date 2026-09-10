! Exercise production accessors and writers with distinguishable complex channels.
program test_flre_quantities
    use, intrinsic :: iso_c_binding, only: c_intptr_t, c_double, c_loc, c_f_pointer
    use kilca_flre_quants_m, only: flre_quants_t, flre_quants_interp_current_density, &
                                   flre_quants_save_profiles, save_current_density_ext
    use output_data, only: flag_quants
    implicit none

    integer, parameter :: n = 8
    type(flre_quants_t), target :: qp
    real(c_double), target :: radius(n)
    real(c_double), pointer :: current(:, :, :, :, :), density(:, :, :)
    integer(c_intptr_t) :: handle
    integer :: k, part, comp, kind, spec, failures
    real(c_double) :: actual(2), position
    character(len=16) :: mode
    character(len=1), parameter :: species(3) = ['i', 'e', 't']
    character(len=1), parameter :: components(3) = ['r', 's', 'p']
    character(len=200) :: filename

    failures = 0
    call get_command_argument(1, mode)
    radius = [(real(k, c_double), k=1, n)]
    qp%dimx = n
    qp%x => radius
    qp%zone_index = 7
    qp%path2linear = './'
    allocate (qp%cdlab(0:36 * n - 1), qp%current_dens(36 * n), qp%number_dens(6 * n))
    call c_f_pointer(c_loc(qp%cdlab), current, [n, 2, 3, 2, 3])
    call c_f_pointer(c_loc(qp%number_dens), density, [n, 2, 3])
    do spec = 1, 3
        do kind = 1, 2
            do comp = 1, 3
                do part = 1, 2
                    current(:, part, comp, kind, spec) = &
                        current_value(radius, part, comp, kind, spec)
                end do
            end do
        end do
        do part = 1, 2
            density(:, part, spec) = 1000 * spec + 100 * part + part * radius
        end do
    end do
    qp%current_dens = qp%cdlab
    handle = transfer(c_loc(qp), handle)

    if (trim(mode) == 'interpolation') then
        do spec = 1, 3
            do kind = 1, 2
                do comp = 1, 3
                    ! Both endpoints, all interior grid points, and half-grid points.
                    do k = 2, 2 * n
                        position = 0.5d0 * k
                        call flre_quants_interp_current_density(handle, position, &
                                                               kind - 1, spec - 1, comp - 1, actual)
                        do part = 1, 2
                            call check('current interpolation', actual(part), &
                                       current_value(position, part, comp, kind, spec))
                        end do
                    end do
                end do
            end do
        end do
    else if (trim(mode) == 'output') then
        allocate (flag_quants(8))
        flag_quants = 0
        flag_quants(1) = 2
        flag_quants(7) = 2
        call flre_quants_save_profiles(handle)
        call save_current_density_ext(qp, qp%cdlab, 'lab', 'rsp')
        do spec = 1, 3
            do kind = 1, 2
                do comp = 1, 3
                    write (filename, '(a,a,a,i0,a,a)') &
                        'zone_7_current_dens_', components(comp), '_', kind - 1, '_', species(spec)
                    call check_file(trim(filename)//'.dat', &
                                   current(:, 1, comp, kind, spec), current(:, 2, comp, kind, spec))
                    call check_file(trim(filename)//'_lab.dat', &
                                   current(:, 1, comp, kind, spec), current(:, 2, comp, kind, spec))
                end do
            end do
            call check_file('zone_7_density_'//species(spec)//'.dat', &
                            density(:, 1, spec), density(:, 2, spec))
        end do
    else
        error stop 'expected interpolation or output argument'
    end if
    if (failures /= 0) then
        print *, 'FAILED FLRE quantities:', failures
        stop 1
    end if
    print *, 'PASS FLRE quantities ', trim(mode)

contains

    elemental real(c_double) function current_value(r, part, comp, kind, spec) result(value)
        real(c_double), intent(in) :: r
        integer, intent(in) :: part, comp, kind, spec
        value = 10000 * spec + 1000 * kind + 100 * comp + 10 * part + (-1)**part * r
    end function current_value

    subroutine check(label, actual, expected)
        character(len=*), intent(in) :: label
        real(c_double), intent(in) :: actual, expected
        if (abs(actual - expected) <= 1.0d-9) return
        failures = failures + 1
        if (failures < 10) print *, 'FAIL ', label, ': got', actual, 'expected', expected
    end subroutine check

    subroutine check_file(name, re, im)
        character(len=*), intent(in) :: name
        real(c_double), intent(in) :: re(n), im(n)
        integer :: unit, ios, row
        real(c_double) :: values(4)
        character(len=300) :: line

        open (newunit=unit, file=name, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'FAIL missing output: ', name
            failures = failures + 1
            return
        end if
        do row = 1, n
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                print *, 'FAIL missing row in ', name
                failures = failures + 1
                exit
            end if
            read (line, *, iostat=ios) values(1:3)
            if (ios /= 0) then
                print *, 'FAIL expected radius, real, imaginary columns in ', name
                failures = failures + 1
                exit
            end if
            call check(name//' radius', values(1), radius(row))
            call check(name//' real', values(2), re(row))
            call check(name//' imaginary', values(3), im(row))
            read (line, *, iostat=ios) values
            if (ios == 0) then
                print *, 'FAIL extra output column in ', name
                failures = failures + 1
            end if
        end do
        if (row > n) then
            read (unit, '(a)', iostat=ios) line
            if (ios >= 0) then
                print *, 'FAIL extra output row in ', name
                failures = failures + 1
            end if
        end if
        close (unit, status='delete')
    end subroutine check_file
end program test_flre_quantities
