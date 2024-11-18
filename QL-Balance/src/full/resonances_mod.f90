
module resonances_mod

    implicit none

    integer :: numres,iunit_res
    double precision, dimension(:), allocatable :: r_res,width_res,ampl_res

    contains 

    subroutine add_mode_group_to_h5_mode_groupname

        use h5mod
        use wave_code_data, only: m_vals, n_vals

        implicit none

        if (numres .eq. 1) then
            if (m_vals(1) <= 9) then
                write (h5_mode_groupname, "(A,I1,A,I1)") "f_", m_vals(1), "_", n_vals(1)
            else
                write (h5_mode_groupname, "(A,I2,A,I1)") "f_", m_vals(1), "_", n_vals(1)
            end if
        else
            write (h5_mode_groupname, "(A)") "multi_mode"
        end if

    end subroutine

end module resonances_mod