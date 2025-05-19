program test_interface

use wave_code_data;

implicit none;

integer, parameter :: nrad = 1000;
real(8), dimension(nrad) :: r_grid;

integer :: mid = 1; ! mode id
integer :: i, j, n1, n2;

real(8) :: dr = (70.0 - 1.0) / (nrad - 1);

integer :: k;

real(8) :: kz, omega_mov_re, omega_mov_im;

complex(8), dimension(1:3,1:3,0:2) :: ctti, ctte, cttt;

complex(8), dimension(1:3,1:3,0:2) :: cct;

character (1) :: flag_back = 'f';

flre_path = 'flre/';
vac_path = 'vac/';

do k = 1,nrad
    r_grid(k) = 1.0 + dr * (k-1);
end do

call initialize_wave_code_interface(nrad, r_grid);

call get_conductivity_matrices(mid, mid); ! conductivity for the first mode

! save matrix: k_mat_* include factor r in (39): k_mat = r * sigma
i = 3;  ! -> index k  in (39)
j = 3;  ! -> index j  in (39)
n1 = 1; ! -> index n  in (39)
n2 = 1; ! -> index n' in (39)

! k_mat_* indices: (node=1:rdim, part={re,im}=1:2, j=1:3, k=1:3, n'=1:flre_order+1, n=1:flre_order+1)

do k = 1,rdim(mid)

    write (100,*) rgrd(mid)%a(k), k_mat_i(mid)%a(k,1,j,i,n2+1,n1+1), k_mat_i(mid)%a(k,2,j,i,n2+1,n1+1); ! ions

    write (200,*) rgrd(mid)%a(k), k_mat_e(mid)%a(k,1,j,i,n2+1,n1+1), k_mat_e(mid)%a(k,2,j,i,n2+1,n1+1); ! electrons

    write (300,*) rgrd(mid)%a(k), k_mat_i(mid)%a(k,1,j,i,n2+1,n1+1) + k_mat_e(mid)%a(k,1,j,i,n2+1,n1+1), &
                                  k_mat_i(mid)%a(k,2,j,i,n2+1,n1+1) + k_mat_e(mid)%a(k,2,j,i,n2+1,n1+1); ! total

end do

!c_mat_* indices: (k=1:3, j=1:3, n=0:1, node=1:rdim)

open (10, file='cct_0_re.dat');
open (20, file='cct_0_im.dat');
open (30, file='cct_1_re.dat');
open (40, file='cct_1_im.dat');

do k = 1,rdim(mid)

    write (10,*) transpose(real(c_mat_i(mid)%a(:,:,0,k) + c_mat_e(mid)%a(:,:,0,k))), rgrd(mid)%a(k);
    write (20,*) transpose(imag(c_mat_i(mid)%a(:,:,0,k) + c_mat_e(mid)%a(:,:,0,k))), rgrd(mid)%a(k);
    write (30,*) transpose(real(c_mat_i(mid)%a(:,:,1,k) + c_mat_e(mid)%a(:,:,1,k))), rgrd(mid)%a(k);
    write (40,*) transpose(imag(c_mat_i(mid)%a(:,:,1,k) + c_mat_e(mid)%a(:,:,1,k))), rgrd(mid)%a(k);

end do

close(10); close(20); close(30); close(40);

open (10, file='cct_0_re_rsp_t.dat');
open (20, file='cct_0_im_rsp_t.dat');
open (30, file='cct_1_re_rsp_t.dat');
open (40, file='cct_1_im_rsp_t.dat');

do k = 1,rdim(mid)

!     cct(:,:,0) = c_mat_i(mid)%a(:,:,0,k) + c_mat_e(mid)%a(:,:,0,k);
!     cct(:,:,1) = c_mat_i(mid)%a(:,:,1,k) + c_mat_e(mid)%a(:,:,1,k);

    do i = 1,3
        do j = 1,3
            cct(j,i,0) = c_mat_i(mid)%a(i,j,0,k) + c_mat_e(mid)%a(i,j,0,k);
            cct(j,i,1) = c_mat_i(mid)%a(i,j,1,k) + c_mat_e(mid)%a(i,j,1,k);
            cct(j,i,2) = cmplx(0.0d0,0.0d0,8);
        end do
    end do

    call transform_c_matrices_to_rsp (rgrd(mid)%a(k), cct);

    write (10,*) real(cct(:,:,0)), rgrd(mid)%a(k);
    write (20,*) imag(cct(:,:,0)), rgrd(mid)%a(k);
    write (30,*) real(cct(:,:,1)), rgrd(mid)%a(k);
    write (40,*) imag(cct(:,:,1)), rgrd(mid)%a(k);

!     write (10,*) transpose(real(cct(:,:,0))), rgrd(mid)%a(k);
!     write (20,*) transpose(imag(cct(:,:,0))), rgrd(mid)%a(k);
!     write (30,*) transpose(real(cct(:,:,1))), rgrd(mid)%a(k);
!     write (40,*) transpose(imag(cct(:,:,1))), rgrd(mid)%a(k);

end do

close(10); close(20); close(30); close(40);

! test against old expressions:
open (10, file='cct_0_re_rsp.dat');
open (20, file='cct_0_im_rsp.dat');
open (30, file='cct_1_re_rsp.dat');
open (40, file='cct_1_im_rsp.dat');

call set_wave_parameters(flre_cd_ptr(1), m_vals(1), n_vals(1));

do k = 1,rdim(mid)

    ctti = cmplx(0.0d0,0.0d0,8);
    ctte = cmplx(0.0d0,0.0d0,8);

    call calc_and_add_galilelian_correction(rgrd(mid)%a(k), 0, flag_back, ctti);

    call calc_and_add_galilelian_correction(rgrd(mid)%a(k), 1, flag_back, ctte);

    cttt = ctti + ctte;

    write (10,*) real(cttt(:,:,0)), rgrd(mid)%a(k);
    write (20,*) imag(cttt(:,:,0)), rgrd(mid)%a(k);
    write (30,*) real(cttt(:,:,1)), rgrd(mid)%a(k);
    write (40,*) imag(cttt(:,:,1)), rgrd(mid)%a(k);

end do

close(10); close(20); close(30); close(40);

! end test

close(100); close(200); close(300); close(400); close(500); close(600); close(700);

call get_mode_parameters(flre_cd_ptr(mid), m_vals(mid), n_vals(mid), kz, omega_mov_re, omega_mov_im);

print *, kz, omega_mov_re, omega_mov_im;

call deallocate_wave_code_data();

end program

!********************************************************************!
