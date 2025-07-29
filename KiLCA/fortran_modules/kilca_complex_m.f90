module kilca_complex_m
    use iso_fortran_env, only: real64, int32
    use kilca_shared_m, only: safe_allocate_cmplx_1d, safe_deallocate_cmplx_1d, &
                              check_allocation
    implicit none
    private
    
    ! Public complex constants (matching C++ O, E, I)
    complex(real64), parameter, public :: cmplx_O = cmplx(0.0_real64, 0.0_real64, real64)
    complex(real64), parameter, public :: cmplx_E = cmplx(1.0_real64, 0.0_real64, real64)
    complex(real64), parameter, public :: cmplx_I = cmplx(0.0_real64, 1.0_real64, real64)
    
    ! Public procedures
    ! Creation functions
    public :: cmplx_make
    public :: cmplx_polar
    
    ! Utility functions
    public :: cmplx_real
    public :: cmplx_imag
    public :: cmplx_abs
    public :: cmplx_abs2
    public :: cmplx_arg
    public :: cmplx_conj
    public :: cmplx_norm
    
    ! Arithmetic operations
    public :: cmplx_add
    public :: cmplx_sub
    public :: cmplx_mult
    public :: cmplx_div
    public :: cmplx_pow_int
    public :: cmplx_pow_real
    public :: cmplx_pow_complex
    
    ! Transcendental functions
    public :: cmplx_exp
    public :: cmplx_log
    public :: cmplx_log10
    public :: cmplx_sqrt
    public :: cmplx_sin
    public :: cmplx_cos
    public :: cmplx_tan
    public :: cmplx_sinh
    public :: cmplx_cosh
    public :: cmplx_tanh
    public :: cmplx_asin
    public :: cmplx_acos
    public :: cmplx_atan
    public :: cmplx_asinh
    public :: cmplx_acosh
    public :: cmplx_atanh
    public :: cmplx_exp2
    public :: cmplx_expm1
    public :: cmplx_log2
    public :: cmplx_log1p
    public :: cmplx_cbrt
    public :: cmplx_hypot
    public :: cmplx_atan2
    
    ! Array operations
    public :: cmplx_allocate_1d
    public :: cmplx_allocate_2d
    public :: cmplx_deallocate_1d
    public :: cmplx_deallocate_2d
    public :: cmplx_array_sum
    public :: cmplx_array_product
    public :: cmplx_array_dot
    public :: cmplx_array_norm2
    
    ! Matrix operations
    public :: cmplx_matrix_trace
    public :: cmplx_matrix_mult
    
    ! I/O formatting
    public :: cmplx_to_string
    public :: cmplx_print
    public :: cmplx_read_formatted
    public :: cmplx_write_formatted
    public :: cmplx_parse_string
    
contains
    
    !============= Creation Functions =============
    
    !> Create complex number from real and imaginary parts
    elemental function cmplx_make(re, im) result(z)
        real(real64), intent(in) :: re, im
        complex(real64) :: z
        
        z = cmplx(re, im, real64)
    end function cmplx_make
    
    !> Create complex number from polar form
    elemental function cmplx_polar(r, theta) result(z)
        real(real64), intent(in) :: r, theta
        complex(real64) :: z
        
        z = cmplx(r * cos(theta), r * sin(theta), real64)
    end function cmplx_polar
    
    !============= Utility Functions =============
    
    !> Get real part
    elemental function cmplx_real(z) result(re)
        complex(real64), intent(in) :: z
        real(real64) :: re
        
        re = real(z, real64)
    end function cmplx_real
    
    !> Get imaginary part
    elemental function cmplx_imag(z) result(im)
        complex(real64), intent(in) :: z
        real(real64) :: im
        
        im = aimag(z)
    end function cmplx_imag
    
    !> Get absolute value (magnitude)
    elemental function cmplx_abs(z) result(r)
        complex(real64), intent(in) :: z
        real(real64) :: r
        
        r = abs(z)
    end function cmplx_abs
    
    !> Get squared absolute value |z|^2
    elemental function cmplx_abs2(z) result(r2)
        complex(real64), intent(in) :: z
        real(real64) :: r2
        
        r2 = real(z, real64)**2 + aimag(z)**2
    end function cmplx_abs2
    
    !> Get argument (phase angle)
    elemental function cmplx_arg(z) result(theta)
        complex(real64), intent(in) :: z
        real(real64) :: theta
        
        theta = atan2(aimag(z), real(z, real64))
    end function cmplx_arg
    
    !> Get complex conjugate
    elemental function cmplx_conj(z) result(zc)
        complex(real64), intent(in) :: z
        complex(real64) :: zc
        
        zc = conjg(z)
    end function cmplx_conj
    
    !> Get norm (same as abs for complex numbers)
    elemental function cmplx_norm(z) result(r)
        complex(real64), intent(in) :: z
        real(real64) :: r
        
        r = abs(z)
    end function cmplx_norm
    
    !============= Arithmetic Operations =============
    
    !> Complex addition
    elemental function cmplx_add(z1, z2) result(z3)
        complex(real64), intent(in) :: z1, z2
        complex(real64) :: z3
        
        z3 = z1 + z2
    end function cmplx_add
    
    !> Complex subtraction
    elemental function cmplx_sub(z1, z2) result(z3)
        complex(real64), intent(in) :: z1, z2
        complex(real64) :: z3
        
        z3 = z1 - z2
    end function cmplx_sub
    
    !> Complex multiplication
    elemental function cmplx_mult(z1, z2) result(z3)
        complex(real64), intent(in) :: z1, z2
        complex(real64) :: z3
        
        z3 = z1 * z2
    end function cmplx_mult
    
    !> Complex division
    elemental function cmplx_div(z1, z2) result(z3)
        complex(real64), intent(in) :: z1, z2
        complex(real64) :: z3
        
        z3 = z1 / z2
    end function cmplx_div
    
    !> Complex power with integer exponent
    elemental function cmplx_pow_int(z, n) result(zn)
        complex(real64), intent(in) :: z
        integer(int32), intent(in) :: n
        complex(real64) :: zn
        
        zn = z**n
    end function cmplx_pow_int
    
    !> Complex power with real exponent
    elemental function cmplx_pow_real(z, x) result(zx)
        complex(real64), intent(in) :: z
        real(real64), intent(in) :: x
        complex(real64) :: zx
        
        if (z == cmplx_O) then
            if (x > 0.0_real64) then
                zx = cmplx_O
            else
                zx = cmplx(huge(1.0_real64), 0.0_real64, real64)  ! Infinity
            end if
        else
            zx = exp(x * log(z))
        end if
    end function cmplx_pow_real
    
    !> Complex power with complex exponent
    elemental function cmplx_pow_complex(z1, z2) result(z3)
        complex(real64), intent(in) :: z1, z2
        complex(real64) :: z3
        
        if (z1 == cmplx_O) then
            if (real(z2, real64) > 0.0_real64) then
                z3 = cmplx_O
            else
                z3 = cmplx(huge(1.0_real64), 0.0_real64, real64)  ! Infinity
            end if
        else
            z3 = exp(z2 * log(z1))
        end if
    end function cmplx_pow_complex
    
    !============= Transcendental Functions =============
    
    !> Complex exponential
    elemental function cmplx_exp(z) result(ez)
        complex(real64), intent(in) :: z
        complex(real64) :: ez
        
        ez = exp(z)
    end function cmplx_exp
    
    !> Complex natural logarithm
    elemental function cmplx_log(z) result(lnz)
        complex(real64), intent(in) :: z
        complex(real64) :: lnz
        
        lnz = log(z)
    end function cmplx_log
    
    !> Complex base-10 logarithm
    elemental function cmplx_log10(z) result(log10z)
        complex(real64), intent(in) :: z
        complex(real64) :: log10z
        
        log10z = log(z) / log(10.0_real64)
    end function cmplx_log10
    
    !> Complex square root
    elemental function cmplx_sqrt(z) result(sqrtz)
        complex(real64), intent(in) :: z
        complex(real64) :: sqrtz
        
        sqrtz = sqrt(z)
    end function cmplx_sqrt
    
    !> Complex sine
    elemental function cmplx_sin(z) result(sinz)
        complex(real64), intent(in) :: z
        complex(real64) :: sinz
        
        sinz = sin(z)
    end function cmplx_sin
    
    !> Complex cosine
    elemental function cmplx_cos(z) result(cosz)
        complex(real64), intent(in) :: z
        complex(real64) :: cosz
        
        cosz = cos(z)
    end function cmplx_cos
    
    !> Complex tangent
    elemental function cmplx_tan(z) result(tanz)
        complex(real64), intent(in) :: z
        complex(real64) :: tanz
        
        tanz = tan(z)
    end function cmplx_tan
    
    !> Complex hyperbolic sine
    elemental function cmplx_sinh(z) result(sinhz)
        complex(real64), intent(in) :: z
        complex(real64) :: sinhz
        
        sinhz = sinh(z)
    end function cmplx_sinh
    
    !> Complex hyperbolic cosine
    elemental function cmplx_cosh(z) result(coshz)
        complex(real64), intent(in) :: z
        complex(real64) :: coshz
        
        coshz = cosh(z)
    end function cmplx_cosh
    
    !> Complex hyperbolic tangent
    elemental function cmplx_tanh(z) result(tanhz)
        complex(real64), intent(in) :: z
        complex(real64) :: tanhz
        
        tanhz = tanh(z)
    end function cmplx_tanh
    
    !> Complex arcsine
    elemental function cmplx_asin(z) result(asinz)
        complex(real64), intent(in) :: z
        complex(real64) :: asinz
        
        asinz = asin(z)
    end function cmplx_asin
    
    !> Complex arccosine
    elemental function cmplx_acos(z) result(acosz)
        complex(real64), intent(in) :: z
        complex(real64) :: acosz
        
        acosz = acos(z)
    end function cmplx_acos
    
    !> Complex arctangent
    elemental function cmplx_atan(z) result(atanz)
        complex(real64), intent(in) :: z
        complex(real64) :: atanz
        
        atanz = atan(z)
    end function cmplx_atan
    
    !> Complex inverse hyperbolic sine (asinh)
    elemental function cmplx_asinh(z) result(asinhz)
        complex(real64), intent(in) :: z
        complex(real64) :: asinhz
        
        ! asinh(z) = log(z + sqrt(z^2 + 1))
        asinhz = log(z + sqrt(z**2 + cmplx_E))
    end function cmplx_asinh
    
    !> Complex inverse hyperbolic cosine (acosh)
    elemental function cmplx_acosh(z) result(acoshz)
        complex(real64), intent(in) :: z
        complex(real64) :: acoshz
        
        ! acosh(z) = log(z + sqrt(z^2 - 1))
        acoshz = log(z + sqrt(z**2 - cmplx_E))
    end function cmplx_acosh
    
    !> Complex inverse hyperbolic tangent (atanh)
    elemental function cmplx_atanh(z) result(atanhz)
        complex(real64), intent(in) :: z
        complex(real64) :: atanhz
        
        ! atanh(z) = 0.5 * log((1 + z) / (1 - z))
        atanhz = 0.5_real64 * log((cmplx_E + z) / (cmplx_E - z))
    end function cmplx_atanh
    
    !> Complex base-2 exponential
    elemental function cmplx_exp2(z) result(exp2z)
        complex(real64), intent(in) :: z
        complex(real64) :: exp2z
        
        ! exp2(z) = 2^z = exp(z * log(2))
        exp2z = exp(z * log(2.0_real64))
    end function cmplx_exp2
    
    !> Complex exponential minus 1 (for better accuracy near z=0)
    elemental function cmplx_expm1(z) result(expm1z)
        complex(real64), intent(in) :: z
        complex(real64) :: expm1z
        
        real(real64) :: absz
        
        absz = abs(z)
        if (absz < 1.0e-2_real64) then
            ! Use Taylor series for small |z|
            ! expm1(z) = z + z^2/2! + z^3/3! + ...
            expm1z = z * (1.0_real64 + z * (0.5_real64 + z * (1.0_real64/6.0_real64 + &
                     z * (1.0_real64/24.0_real64 + z * (1.0_real64/120.0_real64)))))
        else
            expm1z = exp(z) - cmplx_E
        end if
    end function cmplx_expm1
    
    !> Complex base-2 logarithm
    elemental function cmplx_log2(z) result(log2z)
        complex(real64), intent(in) :: z
        complex(real64) :: log2z
        
        ! log2(z) = log(z) / log(2)
        log2z = log(z) / log(2.0_real64)
    end function cmplx_log2
    
    !> Complex logarithm of 1 plus z (for better accuracy near z=0)
    elemental function cmplx_log1p(z) result(log1pz)
        complex(real64), intent(in) :: z
        complex(real64) :: log1pz
        
        real(real64) :: absz
        
        absz = abs(z)
        if (absz < 1.0e-2_real64) then
            ! Use Taylor series for small |z|
            ! log1p(z) = z - z^2/2 + z^3/3 - z^4/4 + ...
            log1pz = z * (1.0_real64 - z * (0.5_real64 - z * (1.0_real64/3.0_real64 - &
                     z * (0.25_real64 - z * (0.2_real64)))))
        else
            log1pz = log(cmplx_E + z)
        end if
    end function cmplx_log1p
    
    !> Complex cube root
    elemental function cmplx_cbrt(z) result(cbrtz)
        complex(real64), intent(in) :: z
        complex(real64) :: cbrtz
        
        ! cbrt(z) = z^(1/3) = exp(log(z)/3)
        if (z == cmplx_O) then
            cbrtz = cmplx_O
        else
            cbrtz = exp(log(z) / 3.0_real64)
        end if
    end function cmplx_cbrt
    
    !> Complex hypotenuse (absolute value of complex number formed from two reals)
    elemental function cmplx_hypot(x, y) result(r)
        real(real64), intent(in) :: x, y
        real(real64) :: r
        
        ! hypot(x, y) = sqrt(x^2 + y^2) with overflow protection
        real(real64) :: ax, ay, w
        
        ax = abs(x)
        ay = abs(y)
        
        if (ax > ay) then
            w = ay / ax
            r = ax * sqrt(1.0_real64 + w**2)
        else if (ay /= 0.0_real64) then
            w = ax / ay
            r = ay * sqrt(1.0_real64 + w**2)
        else
            r = 0.0_real64
        end if
    end function cmplx_hypot
    
    !> Complex two-argument arctangent
    elemental function cmplx_atan2(y, x) result(theta)
        complex(real64), intent(in) :: y, x
        complex(real64) :: theta
        
        ! atan2(y, x) = -i * log((x + i*y) / sqrt(x^2 + y^2))
        complex(real64) :: w
        
        if (x == cmplx_O .and. y == cmplx_O) then
            theta = cmplx_O
        else
            w = x + cmplx_I * y
            theta = -cmplx_I * log(w / sqrt(x**2 + y**2))
        end if
    end function cmplx_atan2
    
    !============= Array Operations =============
    
    !> Allocate 1D complex array
    subroutine cmplx_allocate_1d(array, n)
        complex(real64), allocatable, intent(out) :: array(:)
        integer(int32), intent(in) :: n
        
        call safe_allocate_cmplx_1d(array, n, "complex_array_1d")
    end subroutine cmplx_allocate_1d
    
    !> Allocate 2D complex array
    subroutine cmplx_allocate_2d(array, n1, n2)
        complex(real64), allocatable, intent(out) :: array(:,:)
        integer(int32), intent(in) :: n1, n2
        
        integer :: alloc_stat
        
        allocate(array(n1, n2), stat=alloc_stat)
        call check_allocation(alloc_stat, "complex_array_2d")
    end subroutine cmplx_allocate_2d
    
    !> Deallocate 1D complex array
    subroutine cmplx_deallocate_1d(array)
        complex(real64), allocatable, intent(inout) :: array(:)
        
        call safe_deallocate_cmplx_1d(array, "complex_array_1d")
    end subroutine cmplx_deallocate_1d
    
    !> Deallocate 2D complex array
    subroutine cmplx_deallocate_2d(array)
        complex(real64), allocatable, intent(inout) :: array(:,:)
        
        if (allocated(array)) deallocate(array)
    end subroutine cmplx_deallocate_2d
    
    !> Sum of complex array elements
    function cmplx_array_sum(array) result(sum_val)
        complex(real64), intent(in) :: array(:)
        complex(real64) :: sum_val
        
        sum_val = sum(array)
    end function cmplx_array_sum
    
    !> Product of complex array elements
    function cmplx_array_product(array) result(prod_val)
        complex(real64), intent(in) :: array(:)
        complex(real64) :: prod_val
        
        prod_val = product(array)
    end function cmplx_array_product
    
    !> Dot product of two complex arrays
    function cmplx_array_dot(array1, array2) result(dot_val)
        complex(real64), intent(in) :: array1(:), array2(:)
        complex(real64) :: dot_val
        
        dot_val = dot_product(array1, array2)
    end function cmplx_array_dot
    
    !> L2 norm of complex array
    function cmplx_array_norm2(array) result(norm_val)
        complex(real64), intent(in) :: array(:)
        real(real64) :: norm_val
        
        norm_val = sqrt(real(dot_product(array, array), real64))
    end function cmplx_array_norm2
    
    !============= Matrix Operations =============
    
    !> Trace of complex matrix
    function cmplx_matrix_trace(matrix) result(trace_val)
        complex(real64), intent(in) :: matrix(:,:)
        complex(real64) :: trace_val
        
        integer :: i, n
        
        n = min(size(matrix, 1), size(matrix, 2))
        trace_val = cmplx_O
        
        do i = 1, n
            trace_val = trace_val + matrix(i, i)
        end do
    end function cmplx_matrix_trace
    
    !> Matrix multiplication
    function cmplx_matrix_mult(A, B) result(C)
        complex(real64), intent(in) :: A(:,:), B(:,:)
        complex(real64), allocatable :: C(:,:)
        
        integer :: m, n, k
        
        m = size(A, 1)
        k = size(A, 2)
        n = size(B, 2)
        
        if (k /= size(B, 1)) then
            stop "Matrix dimensions incompatible for multiplication"
        end if
        
        allocate(C(m, n))
        C = matmul(A, B)
    end function cmplx_matrix_mult
    
    !============= I/O Formatting =============
    
    !> Convert complex number to string
    function cmplx_to_string(z, fmt) result(str)
        complex(real64), intent(in) :: z
        character(len=*), intent(in), optional :: fmt
        character(len=100) :: str
        
        character(len=50) :: fmt_use
        
        if (present(fmt)) then
            fmt_use = fmt
        else
            fmt_use = '(F12.6,SP,F12.6,"i")'
        end if
        
        write(str, fmt_use) real(z, real64), aimag(z)
        str = adjustl(str)
    end function cmplx_to_string
    
    !> Print complex number
    subroutine cmplx_print(z, label)
        complex(real64), intent(in) :: z
        character(len=*), intent(in), optional :: label
        
        if (present(label)) then
            write(*, '(A,": ",F12.6,SP,F12.6,"i")') trim(label), real(z, real64), aimag(z)
        else
            write(*, '(F12.6,SP,F12.6,"i")') real(z, real64), aimag(z)
        end if
    end subroutine cmplx_print
    
    !> Read complex number from formatted string (supports multiple formats)
    function cmplx_read_formatted(str) result(z)
        character(len=*), intent(in) :: str
        complex(real64) :: z
        
        character(len=len(str)) :: work_str
        integer :: i, j, ios
        real(real64) :: re, im
        
        ! Initialize
        z = cmplx_O
        work_str = adjustl(str)
        
        ! Try to parse different formats
        
        ! Format 1: (real, imag) with parentheses
        if (work_str(1:1) == '(') then
            ! Find comma position
            i = index(work_str, ',')
            if (i > 0) then
                ! Extract real part
                read(work_str(2:i-1), *, iostat=ios) re
                if (ios == 0) then
                    ! Find closing parenthesis
                    j = index(work_str, ')')
                    if (j > i) then
                        ! Extract imaginary part
                        read(work_str(i+1:j-1), *, iostat=ios) im
                        if (ios == 0) then
                            z = cmplx(re, im, real64)
                            return
                        end if
                    end if
                end if
            end if
        end if
        
        ! Format 2: real+imagi or real-imagi
        i = index(work_str, '+', back=.true.)
        if (i == 0) i = index(work_str, '-', back=.true.)
        if (i > 1) then
            j = index(work_str, 'i', back=.true.)
            if (j > i) then
                read(work_str(1:i-1), *, iostat=ios) re
                if (ios == 0) then
                    read(work_str(i:j-1), *, iostat=ios) im
                    if (ios == 0) then
                        z = cmplx(re, im, real64)
                        return
                    end if
                end if
            end if
        end if
        
        ! Format 3: two space-separated numbers
        read(work_str, *, iostat=ios) re, im
        if (ios == 0) then
            z = cmplx(re, im, real64)
            return
        end if
        
        ! Format 4: single real number (imaginary part = 0)
        read(work_str, *, iostat=ios) re
        if (ios == 0) then
            z = cmplx(re, 0.0_real64, real64)
        end if
        
    end function cmplx_read_formatted
    
    !> Write complex number with specified format
    subroutine cmplx_write_formatted(unit, z, fmt, advance)
        integer, intent(in) :: unit
        complex(real64), intent(in) :: z
        character(len=*), intent(in), optional :: fmt
        character(len=*), intent(in), optional :: advance
        
        character(len=100) :: fmt_use
        character(len=10) :: adv_use
        
        ! Set format
        if (present(fmt)) then
            fmt_use = fmt
        else
            fmt_use = '(ES20.12,SP,ES20.12,"i")'
        end if
        
        ! Set advance option
        if (present(advance)) then
            adv_use = advance
        else
            adv_use = 'YES'
        end if
        
        ! Write complex number
        write(unit, fmt_use, advance=adv_use) real(z, real64), aimag(z)
        
    end subroutine cmplx_write_formatted
    
    !> Parse complex number from string with error handling
    function cmplx_parse_string(str, status) result(z)
        character(len=*), intent(in) :: str
        integer, intent(out), optional :: status
        complex(real64) :: z
        
        character(len=len(str)) :: work_str
        character(len=1), parameter :: delimiters(4) = ['#', char(10), char(13), char(0)]
        integer :: i, ios
        
        ! Initialize
        z = cmplx_O
        if (present(status)) status = 0
        
        ! Remove comments (everything after #)
        work_str = str
        do i = 1, size(delimiters)
            ios = index(work_str, delimiters(i))
            if (ios > 0) work_str = work_str(1:ios-1)
        end do
        
        ! Try to parse
        z = cmplx_read_formatted(work_str)
        
        ! Verify parsing succeeded by checking if result is non-zero
        ! or if the string explicitly represents zero
        work_str = adjustl(work_str)
        if (abs(z) < tiny(1.0_real64) .and. &
            index(work_str, '0') == 0 .and. &
            len_trim(work_str) > 0) then
            if (present(status)) status = -1
        end if
        
    end function cmplx_parse_string
    
end module kilca_complex_m