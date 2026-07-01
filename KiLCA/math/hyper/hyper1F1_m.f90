!> Confluent hypergeometric function 1F1(a,b,z) for a = 1 and complex b, z.
!>
!> Fortran port of the former hyper1F1.cpp. Each routine keeps its single-trailing
!> -underscore C symbol so the existing Fortran callers (KiLCA conductivity,
!> QL-Balance W2_arr) link unchanged. Inputs/outputs are real(c_double) by
!> reference, matching the former extern "C" double* signatures. The quadrature
!> variant defers to fortnum's clean-room 1F1; the others reproduce the Kummer
!> series and continued-fraction evaluations verbatim.
module kilca_hyper1f1_m
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    implicit none
    private

    complex(dp), parameter :: im = (0.0_dp, 1.0_dp)

contains

    subroutine h_quad(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_quad_")
        use fortnum_special_hypergeometric_1f1, only: hyperg_1f1_a1
        use fortnum_status, only: fortnum_status_t
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, res
        type(fortnum_status_t) :: status

        b = cmplx(b_re, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        call hyperg_1f1_a1(b, z, res, status)
        f_re = real(res, dp)
        f_im = aimag(res)
    end subroutine h_quad

    subroutine h_kummer_nmax(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_kummer_nmax_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, term
        integer(int64) :: N, n_

        b = cmplx(b_re, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        N = int(ceiling(-20.0_dp/log10(abs(z/b))), int64) + 5
        if (N < 1 .or. N > 1000000_int64) call warn_nb(b, z)

        term = z/(b + real(N, dp))
        do n_ = N - 1, 0, -1
            term = 1.0_dp + z/(b + real(n_, dp))*term
        end do

        f_re = real(term, dp)
        f_im = aimag(term)
    end subroutine h_kummer_nmax

    subroutine h_kummer_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_kummer_ada_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, term, S1, S2
        real(dp) :: err, eps
        integer(int64) :: n_, Nmax, maxNmax

        b = cmplx(b_re, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        eps = epsilon(1.0_dp)
        maxNmax = 100000000_int64
        err = 0.0_dp; S2 = (0.0_dp, 0.0_dp)

        Nmax = 4
        do while (Nmax < maxNmax)
            term = z/(b + real(Nmax, dp))
            do n_ = Nmax - 1, 0, -1
                term = 1.0_dp + z/(b + real(n_, dp))*term
            end do
            S1 = term

            term = z/(b + real(Nmax + 1, dp))
            do n_ = Nmax, 0, -1
                term = 1.0_dp + z/(b + real(n_, dp))*term
            end do
            S2 = term

            err = min(abs((S2 - S1)/S2), abs(S2 - S1))
            if (err < eps) exit
            Nmax = Nmax*2
        end do

        if (Nmax >= maxNmax) call warn_conv("kummer_ada", Nmax, err, S2, b, z)

        f_re = real(S2, dp)
        f_im = aimag(S2)
    end subroutine h_kummer_ada

    subroutine h_kummer_modified_1(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_kummer_modified_1_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, term
        integer(int64) :: N, n_

        b = cmplx(b_re, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        N = int(ceiling(-20.0_dp/log10(abs(z/b))), int64) + 5
        if (N < 1 .or. N > 1000000_int64) call warn_nb(b, z)

        term = z/(b + real(N, dp))
        do n_ = N - 1, 2, -1
            term = 1.0_dp + z/(b + real(n_, dp))*term
        end do

        f_re = real(term, dp)
        f_im = aimag(term)
    end subroutine h_kummer_modified_1

    subroutine h_kummer_modified_0_nmax(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_kummer_modified_0_nmax_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, term
        integer(int64) :: N, n_

        b = cmplx(b_re, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        N = int(ceiling(-20.0_dp/log10(abs(z/b))), int64) + 5
        if (N < 1 .or. N > 1000000_int64) call warn_nb(b, z)

        term = z/(b + real(N, dp))
        do n_ = N - 1, 3, -1
            term = 1.0_dp + z/(b + real(n_, dp))*term
        end do
        term = term*(z/(b + 2.0_dp))

        f_re = real(term, dp)
        f_im = aimag(term)
    end subroutine h_kummer_modified_0_nmax

    subroutine h_kummer_modified_0_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_kummer_modified_0_ada_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, term, S1, S2
        real(dp) :: err, eps
        integer(int64) :: n_, Nmax, maxNmax

        b = cmplx(b_re, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        eps = epsilon(1.0_dp)
        maxNmax = 100000000_int64
        err = 0.0_dp; S2 = (0.0_dp, 0.0_dp)

        Nmax = 4
        do while (Nmax < maxNmax)
            term = z/(b + real(Nmax, dp))
            do n_ = Nmax - 1, 3, -1
                term = 1.0_dp + (z/(b + real(n_, dp)))*term
            end do
            S1 = term*(z/(b + 2.0_dp))

            term = z/(b + real(Nmax + 1, dp))
            do n_ = Nmax, 3, -1
                term = 1.0_dp + (z/(b + real(n_, dp)))*term
            end do
            S2 = term*(z/(b + 2.0_dp))

            err = abs((S2 - S1)/S2)
            if (err < eps) exit
            Nmax = Nmax*2
        end do

        if (Nmax >= maxNmax) call warn_conv("kummer_modified_0_ada", Nmax, err, S2, b, z)

        f_re = real(S2, dp)
        f_im = aimag(S2)
    end subroutine h_kummer_modified_0_ada

    subroutine h_cont_fract_1_modified_0_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_cont_fract_1_modified_0_ada_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, F11m

        b = cmplx(b_re, b_im, dp)
        z = cmplx(z_re, z_im, dp)

        if (abs(z/b) < 0.1_dp) then
            call h_kummer_modified_0_ada(b_re, b_im, z_re, z_im, f_re, f_im)
        else
            call h_cont_fract_1_inv_ada(b_re, b_im, z_re, z_im, f_re, f_im)
            F11m = (cmplx(f_re, f_im, dp) - 1.0_dp - z/b)*(b/z)*((b + 1.0_dp)/z) - 1.0_dp
            f_re = real(F11m, dp)
            f_im = aimag(F11m)
        end if
    end subroutine h_cont_fract_1_modified_0_ada

    subroutine h_cont_fract_1_inv_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_cont_fract_1_inv_ada_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, term, S1, S2
        real(dp) :: err, eps
        integer(int64) :: n_, Nmax, maxNmax

        b = cmplx(b_re - 1.0_dp, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        eps = epsilon(1.0_dp)
        maxNmax = 1000000_int64
        err = 0.0_dp; S2 = (0.0_dp, 0.0_dp)

        Nmax = 4
        do while (Nmax < maxNmax)
            term = real(Nmax, dp)*z/(b - z + real(Nmax, dp))
            do n_ = Nmax - 1, 1, -1
                term = real(n_, dp)*z/(b - z + real(n_, dp) + term)
            end do
            S1 = b/(b - z + term)

            term = real(Nmax + 1, dp)*z/(b - z + real(Nmax + 1, dp))
            do n_ = Nmax, 1, -1
                term = real(n_, dp)*z/(b - z + real(n_, dp) + term)
            end do
            S2 = b/(b - z + term)

            err = abs((S2 - S1)/S2)
            if (err < eps) exit
            Nmax = Nmax*2
        end do

        if (Nmax >= maxNmax) call warn_conv("cont_fract_1_inv_ada", Nmax, err, S2, b, z)

        f_re = real(S2, dp)
        f_im = aimag(S2)
    end subroutine h_cont_fract_1_inv_ada

    subroutine h_cont_fract_1_inv_nmax(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_cont_fract_1_inv_nmax_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, term
        integer(int64), parameter :: Nmax = 1000_int64
        integer(int64) :: n_

        b = cmplx(b_re - 1.0_dp, b_im, dp)
        z = cmplx(z_re, z_im, dp)

        term = real(Nmax, dp)*z/(b - z + real(Nmax, dp))
        do n_ = Nmax - 1, 1, -1
            term = real(n_, dp)*z/(b - z + real(n_, dp) + term)
        end do
        term = b/(b - z + term)

        f_re = real(term, dp)
        f_im = aimag(term)
    end subroutine h_cont_fract_1_inv_nmax

    subroutine h_cont_fract_1_dir(b_re, b_im, z_re, z_im, f_re, f_im) &
        bind(C, name="hypergeometric1f1_cont_fract_1_dir_")
        real(c_double), intent(in) :: b_re, b_im, z_re, z_im
        real(c_double), intent(out) :: f_re, f_im
        complex(dp) :: b, z, An2, An1, An, Bn2, Bn1, Bn, acoef, bcoef, Sn1, Sn
        real(dp) :: err, eps
        integer(int64) :: n_, Nmax

        b = cmplx(b_re - 1.0_dp, b_im, dp)
        z = cmplx(z_re, z_im, dp)
        eps = epsilon(1.0_dp)
        Nmax = 1000000_int64
        err = 0.0_dp

        An2 = (1.0_dp, 0.0_dp); An1 = (0.0_dp, 0.0_dp)
        Bn2 = (0.0_dp, 0.0_dp); Bn1 = (1.0_dp, 0.0_dp)
        Sn1 = An1/Bn1
        An = An1; Bn = Bn1; Sn = Sn1

        n_ = 1
        do while (n_ < Nmax)
            acoef = real(n_, dp)*z
            bcoef = b - z + real(n_, dp)
            An = bcoef*An1 + acoef*An2
            Bn = bcoef*Bn1 + acoef*Bn2
            Sn = An/Bn
            err = min(abs((Sn - Sn1)/Sn), abs(Sn - Sn1))
            if (err < eps) exit
            An2 = An1; An1 = An
            Bn2 = Bn1; Bn1 = Bn
            Sn1 = Sn
            n_ = n_ + 1
        end do

        Sn = b/(b - z + Sn)

        if (n_ == Nmax) call warn_conv("cont_fract_1_dir", n_, err, Sn, b, z)

        f_re = real(Sn, dp)
        f_im = aimag(Sn)
    end subroutine h_cont_fract_1_dir

    subroutine warn_nb(b, z)
        complex(dp), intent(in) :: b, z
        write (error_unit, '(a)') "warning: hypergeometric1F1_kummer: N out of range"
        write (error_unit, '(a,4es12.4)') "  b,z=", real(b, dp), aimag(b), real(z, dp), aimag(z)
    end subroutine warn_nb

    subroutine warn_conv(tag, Nmax, err, S, b, z)
        character(*), intent(in) :: tag
        integer(int64), intent(in) :: Nmax
        real(dp), intent(in) :: err
        complex(dp), intent(in) :: S, b, z
        write (error_unit, '(a,a,a,i0,a,es12.4)') &
            "warning: hypergeometric1f1_", tag, ": Nmax=", Nmax, " err=", err
        write (error_unit, '(a,4es12.4)') "  S,b,z(partial)=", &
            real(S, dp), aimag(S), real(z, dp), aimag(z)
    end subroutine warn_conv

end module kilca_hyper1f1_m
