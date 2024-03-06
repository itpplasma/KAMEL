!<Saves form-factors to a file

subroutine save_form_facs_fortran (fname, modes_dim, modes_mn, m_min, m_max, n_min, n_max,&
modes_index, dimg, sqr_psi_grid, ffpsi, format_flag)

implicit none;

character(1024), intent(in) :: fname;
integer, intent(in) :: modes_dim;
integer, dimension(2, modes_dim), intent(in) :: modes_mn;
integer, intent(in) :: m_min, m_max, n_min, n_max;
integer, dimension(m_min:m_max, n_min:n_max), intent(in) :: modes_index
integer, intent(in) :: dimg;
real(8), dimension(1:dimg), intent(in) :: sqr_psi_grid;
complex(8), dimension(1:dimg,1:modes_dim), intent(in) :: ffpsi;
integer, intent(in) :: format_flag;

if (format_flag == 0) then

    open(10, form='unformatted', file=trim(fname))

    write (10) modes_dim, dimg, m_min, m_max, n_min, n_max;
    write (10) sqr_psi_grid(1), sqr_psi_grid(dimg);
    !write (10) sqr_psi_grid;
    write (10) modes_index;
    write (10) transpose(ffpsi);

else

    open(10, form='formatted', file=trim(fname))

    write (10,*) modes_dim, dimg, m_min, m_max, n_min, n_max;
    write (10,*) sqr_psi_grid(1), sqr_psi_grid(dimg);
    !write (10,*) sqr_psi_grid;
    write (10,*) modes_index;
    write (10,*) transpose(ffpsi);

end if

close(10);

end subroutine
