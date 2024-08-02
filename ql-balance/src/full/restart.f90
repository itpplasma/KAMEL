module restart_mod

    use control_mod

    implicit none

    logical :: opnd, redostep, scratch
    integer :: iunit_redo

    contains
    subroutine InquiryToRestart

        use parallelTools, only: irank
        use grid_mod, only: npoi, nbaleqs, params, y
        use time_evolution, only: timstep
        use recstep_mod, only: tol
        use baseparam_mod, only: ev
        use control_mod, only: debug_mode

        implicit none

        integer :: k, ipoi, ieq

        inquire(file='restart.dat', exist=opnd)
        if (opnd) then
            if (irank .eq. 0) then
                print *, 'restart'
            end if
            open (201, file='final.restart')
            do ipoi = 1, npoi
                read (201, *) timstep, params(:, ipoi)
                params(3:4, ipoi) = params(3:4, ipoi)*ev
                do ieq = 1, nbaleqs
                    k = nbaleqs*(ipoi - 1) + ieq
                    y(k) = params(ieq, ipoi)
                end do
            end do
            close (201)
            open (201, file='restart.dat')
            read (201, *) timstep
            close (201)
            scratch = .false.
        else
            timstep = timstep*tol
            scratch = .true.
            if (irank .eq. 0) then
                if (debug_mode) write(*,*) 'Debug: start from scratch'
                write(*,*) "timstep = ", timstep
            end if
        end if


    end subroutine

    subroutine redoTimeStep

        use parallelTools, only: irank
        use time_evolution, only: savePrevTranspCoefficients, timstep
        use recstep_mod, only: timstep_arr
        use grid_mod, only: npoic, rc, params, Ercov, params_begbeg
        use baseparam_mod, only: eV, factolmax

        implicit none

        integer :: ipoi

        if (irank .eq. 0) then
            print *, 'redo step with old DQL'
        end if

        call savePrevTranspCoefficients
        iunit_redo = 137

        if (irank .eq. 0) then
            open (iunit_redo, file='params_redostep.after')
            do ipoi = 1, npoic
                write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                    , params(3, ipoi)/ev &
                    , params(4, ipoi)/ev &
                    , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            close (iunit_redo)
        end if
        params = params_begbeg
        if (irank .eq. 0) then
            open (iunit_redo, file='params_redostep.before')
            do ipoi = 1, npoic
                write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                    , params(3, ipoi)/ev &
                    , params(4, ipoi)/ev &
                    , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            close (iunit_redo)
        end if

        timstep = timstep/factolmax
        timstep_arr = timstep

    end subroutine

end module