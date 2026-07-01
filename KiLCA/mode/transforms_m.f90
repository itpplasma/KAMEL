!> Field-frame and coordinate transforms, formerly transforms.{h,cpp}.
!> bind(C) (no trailing underscore, matching the oracle's plain C++ function
!> names - transforms.h never wrapped them in extern "C", but nothing relied
!> on C++ name mangling either, so giving them C linkage here is a safe,
!> behavior-preserving choice): called both from native Fortran (flre_zone,
!> also part of this port) and from still-C++ wave_code_interface.cpp (S7).
!> bp is kept as a parameter for call-site fidelity but ignored, matching
!> every other post-singleton background handle in this port (background is
!> a Fortran singleton; see kilca_background_data_m).
module kilca_transforms_m
    use, intrinsic :: iso_c_binding, only: c_intptr_t, c_double, c_double_complex, c_null_ptr
    use constants, only: dp
    use kilca_background_data_m, only: eval_hthz_native => eval_hthz
    implicit none
    private

    !> Speed of light at full double precision, mirroring the C++ oracle's
    !> constants.h that this ported code descends from. The Fortran constants
    !> module keeps a legacy single-precision c for its own physics.
    real(dp), parameter :: c = 29979245800.0_dp

    public :: galilean_transform_of_EB_fields
    public :: transform_EB_from_cyl_to_rsp
    public :: transform_EB_from_rsp_to_cyl

contains

    !> Transforms (E, B) fields in cylindrical coordinates from a frame 1 to
    !> a frame 2 moving with relative velocity V. Landau II, p. 91:
    !> E' = E + 1/c V x B,  B' = B - 1/c V x E
    !> E = E' - 1/c V x B', B = B' + 1/c V x E'
    !> V x A = e_r (-V A_th) + e_th (V A_r) + e_z (0)
    subroutine galilean_transform_of_EB_fields(EB1, EB2, V) &
        bind(C, name="galilean_transform_of_EB_fields")
        complex(c_double_complex), intent(in) :: EB1(6)
        complex(c_double_complex), intent(out) :: EB2(6)
        real(c_double), value :: V

        EB2(1) = EB1(1) - (V/c)*EB1(5) !Er
        EB2(2) = EB1(2) + (V/c)*EB1(4) !Etheta
        EB2(3) = EB1(3) !Ez

        !for consistence, magnetic fields should not transform
        EB2(4) = EB1(4) ! + (V/c)*EB1(2); !Br
        EB2(5) = EB1(5) ! - (V/c)*EB1(1); !Btheta
        EB2(6) = EB1(6) !Bz
    end subroutine galilean_transform_of_EB_fields

    subroutine transform_EB_from_cyl_to_rsp(bp, rval, EBcyl, EBrsp) &
        bind(C, name="transform_EB_from_cyl_to_rsp")
        integer(c_intptr_t), value :: bp
        real(c_double), value :: rval
        complex(c_double_complex), intent(in) :: EBcyl(6)
        complex(c_double_complex), intent(out) :: EBrsp(6)
        real(dp) :: htz(0:1)

        call eval_hthz_native(rval, 0, 0, c_null_ptr, htz)

        !transform EB from cyl to rsp system:
        EBrsp(1) = EBcyl(1)
        EBrsp(2) = htz(1)*EBcyl(2) - htz(0)*EBcyl(3)
        EBrsp(3) = htz(0)*EBcyl(2) + htz(1)*EBcyl(3)

        EBrsp(4) = EBcyl(4)
        EBrsp(5) = htz(1)*EBcyl(5) - htz(0)*EBcyl(6)
        EBrsp(6) = htz(0)*EBcyl(5) + htz(1)*EBcyl(6)
    end subroutine transform_EB_from_cyl_to_rsp

    subroutine transform_EB_from_rsp_to_cyl(bp, rval, EBrsp, EBcyl) &
        bind(C, name="transform_EB_from_rsp_to_cyl")
        integer(c_intptr_t), value :: bp
        real(c_double), value :: rval
        complex(c_double_complex), intent(in) :: EBrsp(6)
        complex(c_double_complex), intent(out) :: EBcyl(6)
        real(dp) :: htz(0:1)

        call eval_hthz_native(rval, 0, 0, c_null_ptr, htz)

        !transform EB from rsp to cyl system:
        EBcyl(1) = EBrsp(1)
        EBcyl(2) = htz(1)*EBrsp(2) + htz(0)*EBrsp(3)
        EBcyl(3) = -htz(0)*EBrsp(2) + htz(1)*EBrsp(3)

        EBcyl(4) = EBrsp(4)
        EBcyl(5) = htz(1)*EBrsp(5) + htz(0)*EBrsp(6)
        EBcyl(6) = -htz(0)*EBrsp(5) + htz(1)*EBrsp(6)
    end subroutine transform_EB_from_rsp_to_cyl

end module kilca_transforms_m
