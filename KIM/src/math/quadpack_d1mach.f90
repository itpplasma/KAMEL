! Portable replacement for Netlib D1MACH for IEEE double
double precision function d1mach(i)
    implicit none
    integer, intent(in) :: i
    double precision :: one
    one = 1.0d0
    select case(i)
    case(1)
        d1mach = tiny(one)                 ! smallest positive magnitude
    case(2)
        d1mach = huge(one)                 ! largest magnitude
    case(3)
        d1mach = 0.5d0*epsilon(one)        ! smallest relative spacing near 1
    case(4)
        d1mach = epsilon(one)              ! largest relative spacing near 1
    case(5)
        d1mach = log10(2.0d0)              ! log10 of base
    case default
        d1mach = epsilon(one)
    end select
end function d1mach

