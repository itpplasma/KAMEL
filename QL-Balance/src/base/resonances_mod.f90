
module resonances_mod

    implicit none

    integer :: numres,iunit_res
    double precision, dimension(:), allocatable :: r_res,width_res,ampl_res
    logical :: prop=.true.

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

    subroutine write_resonant_radii_to_hdf5
        use KAMEL_hdf5_tools, only: h5_init, h5_deinit, h5_open_rw, h5_close, h5_add_double_1
        use h5mod, only: h5_id, path2out, h5_mode_groupname

        call h5_init()
        call h5_open_rw(path2out, h5_id)
        call h5_add_double_1(h5_id, trim("/"//trim(h5_mode_groupname)//"/r_res"), r_res, &
                            lbound(r_res), ubound(r_res), comment="resonant radii", unit="cm")
        call h5_close(h5_id)
        call h5_deinit()

    end subroutine write_resonant_radii_to_hdf5

end module resonances_mod
