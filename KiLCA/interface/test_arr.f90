program test

use iso_c_binding

implicit none

integer, parameter :: pp = 8    ! pp = 4 for 32 bit machines and pp = 8 for 64 bit machines
integer, parameter :: prec = 8  ! prec = 4 for single precision and prec = 8 for double precision

type pointer2array
    complex(prec), pointer :: a(:,:,:);
end type pointer2array

type(pointer2array), allocatable, dimension(:) :: array

integer :: n1, n2, n3, i1, i2, i3;

type(c_ptr) :: cptr;
!complex(prec), pointer :: arr(:,:,:);

external fill_array;

call fill_array(cptr, n1, n2, n3);

write(*,*), n1, n2, n3;

allocate(array(10));

call c_f_pointer(cptr, array(1)%a,  [n3,n2,n1]);
!call c_f_pointer(cptr, array(1)%a,  [n3,n2,n1]);
call c_f_pointer(cptr, array(10)%a, [n3,n2,n1]);

do i1 = 1,n1
    do i2 = 1,n2
        do i3 = 1,n3
            write(*,*) i1, i2, i3, real(array(1)%a(i3,i2,i1)), imag(array(1)%a(i3,i2,i1));
            write(*,*) i1, i2, i3, real(array(10)%a(i3,i2,i1)), imag(array(10)%a(i3,i2,i1));
        end do
    end do
end do

deallocate(array);

end program

!********************************************************************!
