module zersol_test_callbacks
    use iso_c_binding
    implicit none
    integer :: calls = 0
    complex(c_double_complex), parameter :: root = (0.25d0, -0.375d0)
contains
    function residual(z, data) result(value) bind(C)
        complex(c_double_complex), value :: z
        type(c_ptr), value :: data
        complex(c_double_complex) :: value
        calls = calls + 1
        value = z - root
    end function
    function derivative(z, data) result(value) bind(C)
        complex(c_double_complex), value :: z
        type(c_ptr), value :: data
        complex(c_double_complex) :: value
        value = (1d0, 0d0)
    end function
end module

program test_zersol_bridge
    use zersol_test_callbacks
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    interface
        function create_solver(f, df, data, description, xmin, xmax, ymin, ymax) &
            result(solver) bind(C)
            import
            type(c_funptr), value :: f, df
            type(c_ptr), value :: data
            character(c_char), intent(in) :: description(*)
            real(c_double), value :: xmin, xmax, ymin, ymax
            type(c_ptr) :: solver
        end function
        subroutine set_start_array(solver, count, starts) bind(C)
            import
            type(c_ptr), value :: solver
            integer(c_int), value :: count
            complex(c_double_complex), intent(in) :: starts(*)
        end subroutine
        subroutine set_n_target(solver, count) bind(C)
            import
            type(c_ptr), value :: solver
            integer(c_int), value :: count
        end subroutine
        function find_zeros(solver, capacity, zeros, values, count) result(rc) bind(C)
            import
            type(c_ptr), value :: solver
            integer(c_int), value :: capacity
            complex(c_double_complex), intent(out) :: zeros(*), values(*)
            integer(c_int), intent(out) :: count
            integer(c_int) :: rc
        end function
        subroutine free_solver(solver) bind(C)
            import
            type(c_ptr), value :: solver
        end subroutine
    end interface
    type(c_ptr) :: solver
    complex(c_double_complex) :: starts(1), zeros(4), values(4)
    integer(c_int) :: count, rc

    solver = create_solver(c_funloc(residual), c_funloc(derivative), c_null_ptr, &
                           'complex ABI test'//c_null_char, -1d0, 1d0, -1d0, 1d0)
    if (.not. c_associated(solver)) error stop 'Failed creating root solver'
    starts = (0.2d0, -0.3d0)
    call set_start_array(solver, 1_c_int, starts)
    call set_n_target(solver, 1_c_int)
    rc = find_zeros(solver, 4_c_int, zeros, values, count)
    if (rc /= 0 .or. count /= 1) error stop 'Complex root search failed'
    if (.not. all(ieee_is_finite([real(zeros(1)), aimag(zeros(1))]))) &
        error stop 'Nonfinite root'
    if (abs(zeros(1) - root) > 1d-8 .or. abs(values(1)) > 1d-8) &
        error stop 'Complex callback or result ABI mismatch'
    if (calls < 1) error stop 'Callback was not called'
    call free_solver(solver)
end program
