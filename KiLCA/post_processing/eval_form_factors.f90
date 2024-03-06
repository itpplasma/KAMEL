subroutine eval_form_factors (m, n, s, ff)

implicit none;

!parameters:
integer, intent(in) :: m, n; !poloidal, toridal number

real(8), intent(in) :: s; !radial variable (as stored in formfactors file). Can be: psi_poloidal, psi_toroidal, normalized fluxes, etc...

complex(8), intent(out) :: ff; !form factor for given mode at given radial point

!for fromfactors reading:
integer :: modes_dim;
integer, save :: m_min, m_max, n_min, n_max;
integer, save :: dimg;
integer, allocatable, dimension(:,:), save :: modes_index
real(8), allocatable, dimension(:), save :: sqr_psi_grid;
complex(8), allocatable, dimension(:,:), save :: ffpsi;

real(8), save :: del_sqr_psi;

integer, parameter :: deg = 5; !polynom degree

real(8) :: spsi;

complex(8) :: fac;

integer :: k, l, i, ind;

integer :: format_flag = 1; !0 - unformatted, 1 - formatted.
!Unformatted output is compiler dependent,
!but reading and writing unformatted files is much faster

integer, save :: first_time = 1;

if (first_time == 1) then !load data

    if (format_flag == 0) then

        open(10, form='unformatted', file='./formfactors.dat')

        read (10) modes_dim, dimg, m_min, m_max, n_min, n_max;

        allocate(modes_index(m_min:m_max,n_min:n_max));
        allocate(sqr_psi_grid(1:dimg));
        allocate(ffpsi(1:modes_dim,1:dimg));

        read (10) sqr_psi_grid(1), sqr_psi_grid(dimg);
        read (10) modes_index;
        read (10) ffpsi;

    else

        open(10, form='formatted', file='./formfactors.dat')

        read (10,*) modes_dim, dimg, m_min, m_max, n_min, n_max;

        allocate(modes_index(m_min:m_max,n_min:n_max));
        allocate(sqr_psi_grid(1:dimg));
        allocate(ffpsi(1:modes_dim,1:dimg));

        read (10,*) sqr_psi_grid(1), sqr_psi_grid(dimg);
        read (10,*) modes_index;
        read (10,*) ffpsi;

    end if

    close(10);

    !for interpolation:
    del_sqr_psi = (sqr_psi_grid(dimg)-sqr_psi_grid(1))/(dimg-1);

    do k = 2,dimg-1
        sqr_psi_grid(k) = sqr_psi_grid(1) + (k-1)*del_sqr_psi;
    end do

    print *, modes_dim, dimg, m_min, m_max, n_min, n_max;

    first_time = 0;

end if

if (m<m_min .or. m>m_max .or. n<n_min .or. n>n_max) then
    !(m,n) harmonics is outside of the form-factors array
    ff = cmplx(1.0d0,0.0d0,8);
    return;
end if;

ind = modes_index(m,n);

if (ind<0) then
    !(m,n) harmonics is absent in the form-factors array
    ff = cmplx(1.0d0,0.0d0,8);
    return;
end if

!spsi = sqrt(s)*sqr_psi_grid(dimg); !sqrt(psi) radial variable
spsi = s;

k = floor((spsi-sqr_psi_grid(1))/del_sqr_psi) + 1; !<=node index

if (k < 1 .or. k > dimg) then
    print *, 'warning: s value is outside the data range: extrapolation:', s
end if

k = min(max(1, k-(deg-1)/2), dimg-deg); !starting index for interpolation

ff = cmplx(0.0d0,0.0d0,8);

!Lagrange polynom:
do l = 0,deg
    fac = cmplx(1.0d0,0.0d0,8);
    do i = 0,deg
        if (i /= l) then
            fac = fac*(spsi - sqr_psi_grid(k+i))/(sqr_psi_grid(k+l)-sqr_psi_grid(k+i));
        end if
    end do
    ff = ff + ffpsi(ind,k+l)*fac;
end do

end subroutine
