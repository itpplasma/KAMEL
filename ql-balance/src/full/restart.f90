module restart_mod

    implicit none

    logical :: opnd, dostep, scratch

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
end module