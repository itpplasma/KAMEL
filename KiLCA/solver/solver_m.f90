!> ODE solver for the stiff linear system u' = A(r)*u(r) (the FLRE basis-vector
!> integration), formerly solver.cpp + rhs_func.cpp. Translated to Fortran
!> using the F2003 SUNDIALS CVODE interface (fcvode_mod/fnvector_serial_mod/
!> fsunmatrix_dense_mod/fsunlinsol_dense_mod, enabled in S0) and native LAPACK
!> complex routines (called via implicit/external interfaces, exactly as the
!> C++ oracle treated COMPLEX*16 arrays as flat REAL*8 pairs with no type
!> checking across the call -- using `external` here instead of an explicit
!> interface block avoids re-litigating that real/complex aliasing under
!> Fortran's stricter type checking).
!>
!> solver_a.cpp (a near-duplicate of solver.cpp) and eigtransform.{h,cpp} are
!> NOT in the CMake build (KiLCA/CMakeLists.txt only lists solver.cpp) and
!> have zero live callers - dropped entirely, matching the
!> ADAPTIVE_GRID_GENERATOR=USUAL/cond_profs.cpp precedent. Within the live
!> solver.cpp itself: integrate_basis_vecs_ (trailing-underscore variant) is
!> declared in solver.h but never defined anywhere - dropped.
!> superpose_basis_vecs_ is defined but has zero callers anywhere (not even
!> within solver.cpp) - dropped. The Jacobian path is hard-coded
!> `#if false // USE_JACOBIAN_IN_ODE_SOLVER == 1` in the live solver.cpp
!> (overriding the macro itself being 1 in code_settings.h/CMakeLists.txt) -
!> CVodeSetJacFn is therefore never called, so Jacobian (rhs_func.cpp) and
!> rhs_func_coeff (declared, never defined) are both dropped too.
!>
!> integrate_basis_vecs has exactly one call site in the whole live tree
!> (flre_zone.cpp), always with f = &rhs_func - but the generic SysRHSFcn
!> callback indirection (a module-level function pointer, mirroring the
!> oracle's `SysRHSFcn rhs_mat` global) is kept rather than hard-coding the
!> call to rhs_func, since the entry point's own signature still accepts an
!> arbitrary callback.
module kilca_solver_m
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use fsundials_core_mod
    use fcvode_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use fsunlinsol_dense_mod
    implicit none
    private

    public :: integrate_basis_vecs, rhs_func

    !> Mirrors solver.h's `struct solver_settings`.
    type, bind(C) :: solver_settings_t
        integer(c_int) :: Nort
        real(c_double) :: eps_rel
        real(c_double) :: eps_abs
        real(c_double) :: norm_fac
        integer(c_int) :: debug
    end type solver_settings_t

    !> Mirrors rhs_func.h's `struct rhs_func_params`.
    type, bind(C) :: rhs_func_params_t
        integer(c_int) :: Nwaves
        integer(c_int) :: Nphys
        integer(c_int) :: Nfs
        type(c_ptr) :: Dmat
        integer(c_intptr_t) :: sp
    end type rhs_func_params_t

    abstract interface
        subroutine sys_rhs_fcn(t, y, ydot, params) bind(C)
            import :: c_double, c_ptr
            real(c_double), value :: t
            real(c_double), intent(in) :: y(*)
            real(c_double), intent(out) :: ydot(*)
            type(c_ptr), value :: params
        end subroutine sys_rhs_fcn
    end interface

    !> Mirrors the oracle's file-scope `SysRHSFcn rhs_mat;` global, set once
    !> per integrate_basis_vecs call and read back inside the CVODE RHS
    !> callback (func).
    type(c_funptr) :: rhs_mat_funptr

    external :: zgeqrf, zungqr, ztrtri, ztrmm, zgemm, zgemv

    interface
        function get_sysmat_flag_back_c(handle) result(ch) bind(C, name="get_sysmat_flag_back_")
            import :: c_intptr_t, c_char
            integer(c_intptr_t), value :: handle
            character(kind=c_char) :: ch
        end function get_sysmat_flag_back_c

        subroutine calc_diff_sys_matrix_c(r, flagback, Dmat, fb_len) &
            bind(C, name="calc_diff_sys_matrix_")
            import :: c_double, c_char, c_int
            real(c_double), intent(in) :: r
            character(kind=c_char), intent(in) :: flagback(*)
            real(c_double), intent(out) :: Dmat(*)
            integer(c_int), value :: fb_len
        end subroutine calc_diff_sys_matrix_c
    end interface

contains

    integer(c_int) function signum(x) result(s)
        real(c_double), intent(in) :: x
        if (x < 0.0d0) then
            s = -1
        else if (x == 0.0d0) then
            s = 0
        else
            s = 1
        end if
    end function signum

    !> CVODE RhsFn callback: dispatches to whatever was registered as `f` in
    !> integrate_basis_vecs (mirrors the oracle's `func`).
    integer(c_int) function cvode_rhs_func(t, sunvec_y, sunvec_ydot, user_data) &
        result(ierr) bind(C)
        real(c_double), value :: t
        type(N_Vector) :: sunvec_y
        type(N_Vector) :: sunvec_ydot
        type(c_ptr), value :: user_data
        real(c_double), pointer :: yvec(:), ydotvec(:)
        procedure(sys_rhs_fcn), pointer :: rhs_mat

        yvec => FN_VGetArrayPointer(sunvec_y)
        ydotvec => FN_VGetArrayPointer(sunvec_ydot)

        call c_f_procpointer(rhs_mat_funptr, rhs_mat)
        call rhs_mat(t, yvec, ydotvec, user_data)

        ierr = 0
    end function cvode_rhs_func

    !> f: a SysRHSFcn rhs callback; Nfs: number of fundamental solutions
    !> integrated simultaneously; Nw: dimension of the problem (number of
    !> waves); dim: dimension of the r-grid; rvec: grid points where the
    !> solution is needed; Smat: on entrance Smat(0:Neq-1) holds the starting
    !> values (re,im,re,im,...), on exit holds the basis vectors at every
    !> rvec point packed one after another.
    function integrate_basis_vecs(f, Nfs, Nw, dim, rvec, Smat, ss_ptr, params) &
        result(ret) bind(C, name="integrate_basis_vecs")
        type(c_funptr), value :: f
        integer(c_int), value :: Nfs, Nw, dim
        real(c_double), intent(in) :: rvec(0:dim - 1)
        real(c_double), target, intent(inout) :: Smat(0:2*Nfs*Nw*dim - 1)
        type(c_ptr), value :: ss_ptr
        type(c_ptr), value :: params
        integer(c_int) :: ret

        type(solver_settings_t), pointer :: ss
        type(c_ptr) :: sunctx, cvode_mem
        type(SUNMatrix), pointer :: A
        type(SUNLinearSolver), pointer :: LS
        type(N_Vector), pointer :: y, yval

        integer(c_int) :: Neq, Nort, flag, flag2
        real(c_double), allocatable, target :: mem(:)
        integer(c_int) :: rdata_i, ydata_i, taudata_i
        real(c_double) :: reltol, abstol, rf

        integer(c_int) :: i, step, ort_flag, info, lwork
        real(c_double), allocatable :: work(:)
        real(c_double) :: work_query(1)
        real(c_double) :: maxv, minv, modv
        integer(c_int) :: ind, rpind, tot_steps, dirint
        integer(c_int) :: adr_i

        ret = 0

        flag = FSUNContext_Create(SUN_COMM_NULL, sunctx)
        if (flag /= 0) then
            write (error_unit, '(a)') 'Error creating SUNContext'
            ret = 1
            return
        end if

        call c_f_pointer(ss_ptr, ss)

        rhs_mat_funptr = f

        Neq = 2*Nfs*Nw
        Nort = ss%Nort

        allocate (mem(0:Nort*(1 + Neq + 2*Nfs) - 1))

        rdata_i = 0
        ydata_i = Nort
        taudata_i = Nort + Neq*Nort

        mem(rdata_i) = rvec(0)
        do i = 0, Neq - 1
            mem(ydata_i + i) = Smat(i)
        end do

        y => FN_VMake_Serial(int(Neq, c_int64_t), mem(ydata_i:ydata_i + Neq - 1), sunctx)
        if (.not. associated(y)) then
            write (error_unit, '(a)') 'error: int_basis_vecs: y vector allocation failed!..'
            ret = 1
            return
        end if

        yval => FN_VMake_Serial(int(Neq, c_int64_t), mem(ydata_i:ydata_i + Neq - 1), sunctx)
        if (.not. associated(yval)) then
            write (error_unit, '(a)') 'error: int_basis_vecs: yval vector allocation failed!..'
            ret = 1
            return
        end if

        reltol = ss%eps_rel
        abstol = ss%eps_abs

        cvode_mem = FCVodeCreate(CV_ADAMS, sunctx)
        if (.not. c_associated(cvode_mem)) then
            write (error_unit, '(a)') 'error: int_basis_vecs: cvodecreate failed!..'
            ret = 1
            return
        end if

        flag = FCVodeInit(cvode_mem, c_funloc(cvode_rhs_func), rvec(0), y)
        if (flag /= 0) then
            write (error_unit, '(a)') 'error: int_basis_vecs: CVodeInit failed!..'
            ret = 1
            return
        end if

        flag = FCVodeSetMaxOrd(cvode_mem, 12)
        if (flag /= 0) then
            write (error_unit, '(a)') 'error: int_basis_vecs: CVodeSetMaxOrd failed!..'
            ret = 1
            return
        end if

        flag = FCVodeSStolerances(cvode_mem, reltol, abstol)
        if (flag /= 0) then
            write (error_unit, '(a)') 'error: int_basis_vecs: CVodeSStolerances failed!..'
            ret = 1
            return
        end if

        A => FSUNDenseMatrix(int(Neq, c_int64_t), int(Neq, c_int64_t), sunctx)
        if (.not. associated(A)) then
            write (error_unit, '(a)') 'error: int_basis_vecs: SUNDenseMatrix failed!..'
            ret = 1
            return
        end if

        LS => FSUNLinSol_Dense(y, A, sunctx)
        if (.not. associated(LS)) then
            write (error_unit, '(a)') 'error: int_basis_vecs: SUNLinSol_Dense failed!..'
            ret = 1
            return
        end if

        flag = FCVodeSetLinearSolver(cvode_mem, LS, A)
        if (flag /= 0) then
            write (error_unit, '(a)') 'error: int_basis_vecs: CVodeSetLinearSolver failed!..'
            ret = 1
            return
        end if

        rf = rvec(dim - 1)

        flag = FCVodeSetStopTime(cvode_mem, rf)
        if (flag /= 0) then
            write (error_unit, '(a)') 'error: int_basis_vecs: cvodestoptime failed!..'
            ret = 1
            return
        end if

        flag = FCVodeSetUserData(cvode_mem, params)
        if (flag /= 0) then
            write (error_unit, '(a)') 'error: int_basis_vecs: CVodeSetUserData failed!..'
            ret = 1
            return
        end if

        ! Jacobian setup is hard-coded dead (#if false) in the oracle - omitted.

        ! Determine the optimal LWORK; ydata must not change here.
        lwork = -1
        call zgeqrf(Nw, Nfs, mem(ydata_i:), Nw, mem(taudata_i:), work_query, lwork, info)
        if (info /= 0) then
            write (error_unit, '(a,i0)') 'error: int_basis_vecs: zgeqrf_ failed!: ', info
            ret = 1
            return
        end if

        lwork = int(work_query(1))
        allocate (work(0:2*lwork - 1))

        rpind = 0
        tot_steps = 0
        dirint = signum(rvec(dim - 1) - rvec(0))
        step = 0

        do
            if (step == Nort - 1) then
                write (error_unit, '(a,i0)') &
                    'error: int_basis_vecs: maximum number of ortonormalization steps is reached: ', step
                exit
            end if

            flag = FCVode(cvode_mem, rf, y, mem(rdata_i:rdata_i), CV_ONE_STEP)
            tot_steps = tot_steps + 1

            if (.not. (flag == CV_SUCCESS .or. flag == CV_TSTOP_RETURN)) then
                write (error_unit, '(a,es12.4,a,i0)') &
                    'error: int_basis_vecs: cvode failed!: t=', mem(rdata_i), ' flag=', flag
                exit
            end if

            do i = rpind + 1, dim - 1
                if (dirint*(rvec(i) - mem(rdata_i)) <= 0.0d0) then
                    call FN_VSetArrayPointer(Smat(i*Neq:i*Neq + Neq - 1), yval)
                    flag2 = FCVodeGetDky(cvode_mem, rvec(i), 0, yval)

                    if (flag2 /= CV_SUCCESS) then
                        write (error_unit, '(a,es12.4,a,i0)') &
                            'error: int_basis_vecs_: cvodegetdky failed!: t=', mem(rdata_i), ' flag=', flag2
                        ret = 1
                        return
                    end if
                    rpind = i
                else
                    exit
                end if
            end do

            if (flag == CV_TSTOP_RETURN) then
                if (rpind /= dim - 1) then
                    write (error_unit, '(a,i0,a,i0)') &
                        'error: int_basis_vecs: a wrong index is detected: rpind=', rpind, ' dim=', dim
                end if
                exit
            end if

            ! QR to check orthogonality; if still OK, continue.
            adr_i = ydata_i + Neq
            do i = 0, Neq - 1
                mem(adr_i + i) = mem(ydata_i + i)
            end do

            call zgeqrf(Nw, Nfs, mem(adr_i:), Nw, mem(taudata_i:), work, lwork, info)
            if (info /= 0) then
                write (error_unit, '(a,i0)') 'error: int_basis_vecs: zgeqrf_ failed!: ', info
                exit
            end if

            maxv = sqrt(mem(adr_i)**2 + mem(adr_i + 1)**2)
            minv = maxv

            do i = 1, Nfs - 1
                ind = 2*i*(Nw + 1)
                modv = sqrt(mem(adr_i + ind)**2 + mem(adr_i + ind + 1)**2)
                if (modv > maxv) maxv = modv
                if (modv < minv) minv = modv
            end do

            if (maxv > ss%norm_fac*minv) then
                ort_flag = 0
            else
                ort_flag = 1
            end if

            if (ort_flag == 1) cycle

            if (ss%debug > 1) then
                write (output_unit, '(a,i5,a,f9.6,a,es8.3,a,es8.3)') &
                    'int_basis_vecs: QRstep = ', step, '  r = ', mem(rdata_i), &
                    '  min = ', minv, '  max = ', maxv
            end if

            do i = 0, Neq - 1
                mem(ydata_i + i) = mem(adr_i + i)
            end do

            rdata_i = rdata_i + 1
            mem(rdata_i) = mem(rdata_i - 1)

            ydata_i = ydata_i + Neq
            call zungqr(Nw, Nfs, Nfs, mem(ydata_i:), Nw, mem(taudata_i:), work, lwork, info)
            if (info /= 0) then
                write (error_unit, '(a,i0)') 'error: int_basis_vecs: zungqr_ failed!: ', info
                exit
            end if

            call FN_VSetArrayPointer(mem(ydata_i:ydata_i + Neq - 1), y)

            taudata_i = taudata_i + 2*Nfs
            step = step + 1

            flag = FCVodeReInit(cvode_mem, mem(rdata_i), y)
            if (flag /= 0) then
                write (error_unit, '(a,i0)') 'error: int_basis_vecs: cvodereinit failed!: flag=', flag
                exit
            end if
        end do

        if (ss%debug > 0) then
            write (output_unit, '(a,i0,a,i0,a,f0.6)') &
                'information: int_basis_vecs: Nsteps = ', tot_steps, ' Norts = ', step, ' r = ', mem(rdata_i)
        end if

        if ((rdata_i /= step) .or. (ydata_i /= Nort + Neq*step) .or. &
            (taudata_i /= Nort + Neq*Nort + 2*Nfs*step)) then
            write (error_unit, '(a)') 'error: int_basis_vecs: wrong pointers to renormalization info:'
            ret = 1
            return
        end if

        Nort = step

        call renorm_basis_vecs(Nfs, Nw, dim, rvec, Smat, (dim - 1)*Neq, Nort, &
                               mem, rdata_i, ydata_i, taudata_i)

        call FN_VDestroy(y)
        call FN_VDestroy(yval)
        flag = FSUNLinSolFree(LS)
        call FSUNMatDestroy(A)
        call FCVodeFree(cvode_mem)
        flag = FSUNContext_Free(sunctx)

        deallocate (work)
        deallocate (mem)
    end function integrate_basis_vecs

    !> All complex arrays are stored as flat double arrays of double length
    !> (re,im interleaved), exactly as the oracle. rdata_i/ydata_i/taudata_i
    !> are the 0-based offsets into mem at the LAST orthonormalization step
    !> (mirrors the oracle's rdata-1/ydata-Neq/taudata-2*Nfs pointers).
    !>
    !> udata_i0 is the 0-based offset (into the FULL udata array, which the
    !> caller now passes in its entirety rather than a section starting at
    !> the last block) of the last grid point's block - (dim-1)*Neq. The
    !> oracle walks a raw pointer `udata -= Neq` from that last block back
    !> to the first (absolute position 0); an earlier translation attempt
    !> passed the caller's Smat((dim-1)*Neq:) SECTION as the actual argument
    !> and let udata_i go NEGATIVE relative to that section's own start to
    !> reach the same absolute positions - legal C pointer arithmetic, but
    !> gfortran's -fcheck=all (enabled by this project's build) flags any
    !> negative index into an assumed-size dummy regardless of whether the
    !> underlying memory access is actually in-bounds. Passing the FULL
    !> array and tracking udata_i as an ABSOLUTE 0-based offset (starting at
    !> udata_i0, walking down to 0, never negative) reaches the exact same
    !> memory locations without ever forming an out-of-declared-range index.
    subroutine renorm_basis_vecs(Nfs, Nw, dim, rvec, udata, udata_i0, Nort, mem, &
                                  rdata_i0, ydata_i0, taudata_i0)
        integer(c_int), intent(in) :: Nfs, Nw, dim, Nort
        real(c_double), intent(in) :: rvec(0:dim - 1)
        real(c_double), intent(inout) :: udata(*)
        integer(c_int), intent(in) :: udata_i0
        real(c_double), intent(inout) :: mem(0:)
        integer(c_int), intent(in) :: rdata_i0, ydata_i0, taudata_i0

        integer(c_int) :: Neq, i, k, ind
        real(c_double) :: tmp(0:2*Nfs*Nw - 1), buf(0:2*Nfs*Nw - 1)
        real(c_double) :: alpha(2), beta(2)
        character(len=1) :: side, trans, uplo, diag
        integer(c_int) :: info, dirint, step, ropind
        integer(c_int) :: rdata_i, ydata_i, taudata_i, udata_i

        Neq = 2*Nfs*Nw

        alpha = [1.0d0, 0.0d0]
        beta = [0.0d0, 0.0d0]
        side = 'L'; trans = 'N'; uplo = 'U'; diag = 'N'

        do k = 0, Nw - 1
            do i = 0, Nfs - 1
                ind = 2*Nw*i + 2*k
                buf(ind) = 0.0d0
                buf(ind + 1) = 0.0d0
                if (i == k) buf(ind) = 1.0d0
            end do
        end do

        dirint = signum(rvec(dim - 1) - rvec(0))
        ropind = Nort

        rdata_i = rdata_i0
        ydata_i = ydata_i0
        taudata_i = taudata_i0
        udata_i = udata_i0

        do k = dim - 1, 0, -1
            do step = ropind - 1, 0, -1
                if (dirint*(rvec(k) - mem(rdata_i)) > 0.0d0) exit

                ropind = ropind - 1

                call ztrtri(uplo, diag, Nfs, mem(ydata_i:), Nw, info)
                if (info /= 0) then
                    write (error_unit, '(a,i0,a,i0,a,f0.6)') &
                        'error: renorm_basis_vecs_: ztrtri_ failed!: info=', info, ' k=', k, ' r=', rvec(k)
                end if

                call ztrmm(side, uplo, trans, diag, Nfs, Nfs, alpha, mem(ydata_i:), Nw, buf, Nw)

                rdata_i = rdata_i - 1
                ydata_i = ydata_i - Neq
                taudata_i = taudata_i - 2*Nfs
            end do

            call zgemm(trans, trans, Nw, Nfs, Nfs, alpha, udata(udata_i + 1), Nw, buf, Nw, beta, tmp, Nw)

            do i = 0, Neq - 1
                udata(udata_i + i + 1) = tmp(i)
            end do

            udata_i = udata_i - Neq
        end do
    end subroutine renorm_basis_vecs

    !> Mirrors rhs_func.cpp's rhs_func exactly: fills Dmat via the existing
    !> (Fortran-resident) calc_diff_sys_matrix_, then multiplies it onto y to
    !> get ydot. USE_SPLINES_IN_RHS_EVALUATION's spline branch is dead
    !> (hard-coded 0), matching the precedent established for sysmat_profiles.
    subroutine rhs_func(r, y, ydot, params) bind(C, name="rhs_func")
        real(c_double), value :: r
        real(c_double), intent(in) :: y(*)
        real(c_double), intent(out) :: ydot(*)
        type(c_ptr), value :: params

        type(rhs_func_params_t), pointer :: fp
        real(c_double), pointer :: Dmat(:)
        character(kind=c_char) :: flag_back_buf(1)
        real(c_double) :: alpha(2), beta(2)
        character(len=1) :: trans
        integer(c_int) :: Nw, Nfs

        call c_f_pointer(params, fp)
        call c_f_pointer(fp%Dmat, Dmat, [2*fp%Nwaves*(fp%Nwaves + 2*fp%Nfs)])

        flag_back_buf(1) = get_sysmat_flag_back_c(fp%sp)
        call calc_diff_sys_matrix_c(r, flag_back_buf, Dmat, 1_c_int)

        alpha = [1.0d0, 0.0d0]
        beta = [0.0d0, 0.0d0]
        trans = 'N'

        Nw = fp%Nwaves
        Nfs = fp%Nfs

        call zgemm(trans, trans, Nw, Nfs, Nw, alpha, Dmat, Nw, y, Nw, beta, ydot, Nw)
    end subroutine rhs_func

end module kilca_solver_m
