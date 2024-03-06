program testImn;

implicit none;

integer :: i;

complex(8) :: x1, x2;

complex(8), dimension(0:3, 0:3) :: I0, I1, I2;

do while(.true.)

    write(*,*) 'x1 =';
    read(*,*) x1;

    write(*,*) 'x2 =';
    read(*,*) x2;

    write(*,*) x1, x2;

    call calc_Imn_array0(x1, x2, I0);
    call calc_Imn_array1(x1, x2, I1);
    call calc_Imn_array2(x1, x2, I2);

    i = 0;
    write(*,*) '(0,0:3):'
    write(*,*) I0(i,0:3);
    write(*,*) I1(i,0:3);
    write(*,*) I2(i,0:3);

    i = 1;
    write(*,*) '(1,0:3):'
    write(*,*) I0(i,0:3);
    write(*,*) I1(i,0:3);
    write(*,*) I2(i,0:3);

    i = 2;
    write(*,*) '(2,0:3):'
    write(*,*) I0(i,0:3);
    write(*,*) I1(i,0:3);
    write(*,*) I2(i,0:3);

    i = 3;
    write(*,*) '(3,0:3):'
    write(*,*) I0(i,0:3);
    write(*,*) I1(i,0:3);
    write(*,*) I2(i,0:3);

end do

end program
