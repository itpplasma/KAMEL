!> Per-zone snapshot of the Maxwell-equations system layout, formerly the C++
!> maxwell_eqs_data class. Multiple flre_zone instances can be alive at once,
!> each with its own snapshot (the maxwell_equations Fortran module itself is
!> transient/shared, reused and cleaned between zones), so this uses the same
!> opaque-handle pattern as kilca_spline_m: each instance is heap-allocated and
!> addressed by a C pointer stored as an integer(c_intptr_t) handle.
!>
!> The data-fill calls below bind to the existing C ABI entry points
!> (copy_module_data_to_maxwell_eqs_data_struct_f_, etc.) exactly as the former
!> C++ constructor did, including der_order's C-row-major flat layout, so the
!> values are byte-identical to the original.
module kilca_maxwell_eqs_data_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr, c_loc, c_f_pointer
    implicit none
    private

    public :: maxwell_eqs_data_create, maxwell_eqs_data_destroy
    public :: get_me_num_vars, get_me_num_eqs, get_me_der_order
    public :: get_me_dim_ersp_state, get_me_iersp_state
    public :: get_me_dim_ersp_sys, get_me_iersp_sys
    public :: get_me_dim_brsp_sys, get_me_ibrsp_sys
    public :: get_me_sys_ind, get_me_nwaves

    type :: maxwell_eqs_data_t
        integer(c_int) :: Nwaves
        integer(c_int) :: num_vars, num_eqs
        integer(c_int) :: der_order_flat(0:8)
        integer(c_int) :: dim_Ersp_state(0:2), iErsp_state(0:2)
        integer(c_int) :: dim_Ersp_sys(0:2), iErsp_sys(0:2)
        integer(c_int) :: dim_Brsp_sys(0:2), iBrsp_sys(0:2)
        integer(c_int), allocatable :: sys_ind(:)
    end type maxwell_eqs_data_t

    interface
        subroutine copy_module_data_to_maxwell_eqs_data_struct_f(num_vars_p, num_eqs_p, &
            dim_Ersp_sys_p, iErsp_sys_p, dim_Brsp_sys_p, iBrsp_sys_p, der_order_p) &
            bind(C, name="copy_module_data_to_maxwell_eqs_data_struct_f_")
            import :: c_int
            integer(c_int), intent(out) :: num_vars_p, num_eqs_p
            integer(c_int), intent(out) :: dim_Ersp_sys_p(3), iErsp_sys_p(3)
            integer(c_int), intent(out) :: dim_Brsp_sys_p(3), iBrsp_sys_p(3)
            integer(c_int), intent(out) :: der_order_p(9)
        end subroutine copy_module_data_to_maxwell_eqs_data_struct_f

        subroutine get_ersp_state_indices_and_dims_f(dim_Ersp_state_p, iErsp_state_p) &
            bind(C, name="get_ersp_state_indices_and_dims_f_")
            import :: c_int
            integer(c_int), intent(out) :: dim_Ersp_state_p(3), iErsp_state_p(3)
        end subroutine get_ersp_state_indices_and_dims_f

        subroutine get_sys_ind_array_f(sys_ind_p, n) bind(C, name="get_sys_ind_array_f_")
            import :: c_int
            integer(c_int), value :: n
            integer(c_int), intent(out) :: sys_ind_p(n)
        end subroutine get_sys_ind_array_f
    end interface

contains

    !> Allocates a new instance and snapshots the current (transient) Fortran
    !> maxwell_equations module state into it, applying the same -1 index
    !> adjustment (Fortran 1-based -> C 0-based) the former C++ constructor did.
    function maxwell_eqs_data_create(Nwaves) result(handle) &
        bind(C, name="maxwell_eqs_data_create_")
        integer(c_int), value :: Nwaves
        integer(c_intptr_t) :: handle
        type(maxwell_eqs_data_t), pointer :: me
        integer(c_int) :: der_order_flat(9)
        integer(c_int) :: i

        allocate (me)
        me%Nwaves = Nwaves

        call copy_module_data_to_maxwell_eqs_data_struct_f(me%num_vars, me%num_eqs, &
            me%dim_Ersp_sys, me%iErsp_sys, me%dim_Brsp_sys, me%iBrsp_sys, der_order_flat)
        me%der_order_flat = der_order_flat

        call get_ersp_state_indices_and_dims_f(me%dim_Ersp_state, me%iErsp_state)

        allocate (me%sys_ind(Nwaves))
        call get_sys_ind_array_f(me%sys_ind, Nwaves)

        do i = 0, 2
            me%iErsp_state(i) = me%iErsp_state(i) - 1
            me%iErsp_sys(i) = me%iErsp_sys(i) - 1
            me%iBrsp_sys(i) = me%iBrsp_sys(i) - 1
        end do
        do i = 1, Nwaves
            me%sys_ind(i) = me%sys_ind(i) - 1
        end do

        handle = transfer(c_loc(me), handle)
    end function maxwell_eqs_data_create

    subroutine maxwell_eqs_data_destroy(handle) bind(C, name="maxwell_eqs_data_destroy_")
        integer(c_intptr_t), value :: handle
        type(maxwell_eqs_data_t), pointer :: me

        if (handle == 0_c_intptr_t) return
        call handle_to_me(handle, me)
        if (allocated(me%sys_ind)) deallocate (me%sys_ind)
        deallocate (me)
    end subroutine maxwell_eqs_data_destroy

    subroutine handle_to_me(handle, me)
        integer(c_intptr_t), value :: handle
        type(maxwell_eqs_data_t), pointer, intent(out) :: me
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, me)
    end subroutine handle_to_me

    integer(c_int) function get_me_nwaves(handle) bind(C, name="get_me_nwaves_")
        integer(c_intptr_t), value :: handle
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_nwaves = me%Nwaves
    end function get_me_nwaves

    integer(c_int) function get_me_num_vars(handle) bind(C, name="get_me_num_vars_")
        integer(c_intptr_t), value :: handle
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_num_vars = me%num_vars
    end function get_me_num_vars

    integer(c_int) function get_me_num_eqs(handle) bind(C, name="get_me_num_eqs_")
        integer(c_intptr_t), value :: handle
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_num_eqs = me%num_eqs
    end function get_me_num_eqs

    !> Matches the original C++ me->der_order[i][j] (0-based i,j): the C ABI
    !> call above filled der_order_flat with C row-major [3][3] layout.
    integer(c_int) function get_me_der_order(handle, i, j) bind(C, name="get_me_der_order_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: i, j
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_der_order = me%der_order_flat(i*3 + j)
    end function get_me_der_order

    integer(c_int) function get_me_dim_ersp_state(handle, k) bind(C, name="get_me_dim_ersp_state_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: k
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_dim_ersp_state = me%dim_Ersp_state(k)
    end function get_me_dim_ersp_state

    integer(c_int) function get_me_iersp_state(handle, k) bind(C, name="get_me_iersp_state_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: k
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_iersp_state = me%iErsp_state(k)
    end function get_me_iersp_state

    integer(c_int) function get_me_dim_ersp_sys(handle, k) bind(C, name="get_me_dim_ersp_sys_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: k
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_dim_ersp_sys = me%dim_Ersp_sys(k)
    end function get_me_dim_ersp_sys

    integer(c_int) function get_me_iersp_sys(handle, k) bind(C, name="get_me_iersp_sys_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: k
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_iersp_sys = me%iErsp_sys(k)
    end function get_me_iersp_sys

    integer(c_int) function get_me_dim_brsp_sys(handle, k) bind(C, name="get_me_dim_brsp_sys_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: k
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_dim_brsp_sys = me%dim_Brsp_sys(k)
    end function get_me_dim_brsp_sys

    integer(c_int) function get_me_ibrsp_sys(handle, k) bind(C, name="get_me_ibrsp_sys_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: k
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_ibrsp_sys = me%iBrsp_sys(k)
    end function get_me_ibrsp_sys

    integer(c_int) function get_me_sys_ind(handle, k) bind(C, name="get_me_sys_ind_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: k
        type(maxwell_eqs_data_t), pointer :: me
        call handle_to_me(handle, me)
        get_me_sys_ind = me%sys_ind(k + 1)
    end function get_me_sys_ind

end module kilca_maxwell_eqs_data_m
