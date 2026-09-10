
!--------------------------------------------------------------------

subroutine setup_flre_data_module (zone)
    use kilca_flre_zone_m, only: set_collisions_settings_c => set_collisions_settings_c_
    use kilca_flre_zone_m, only: set_conductivity_settings_c => set_conductivity_settings_c_
    use kilca_flre_zone_m, only: set_equations_settings_c => set_equations_settings_c_

use constants, only: pp;

use flre_sett;

implicit none;

integer(pp), intent(in) :: zone;

Nbmax = 2;

call set_conductivity_settings_c (zone, flre_order, Nmax, gal_corr, rsp);

call set_equations_settings_c (zone, hom_sys, Nwaves, Nfs, Nphys, flag_debug);

call set_collisions_settings_c (zone, collmod);

if (flag_debug > 0) then
    call print_flre_settings ();
end if

end subroutine

!--------------------------------------------------------------------

subroutine print_flre_settings ()

use flre_sett;

implicit none;

print *;

print *, "order of FLR expansion: ", flre_order;
print *, "highest cyclotron harmonic: ", Nmax;
print *, "flag if to use correction term in conductivity: ", gal_corr;
print *, "Expansion degree for Bessels: ", Nbmax;
print *, "rsp flag: ", rsp;

print *, "flag for homogenious system: ", hom_sys;
print *, "number of waves: ", Nwaves;
print *, "number of fundamental solutions: ", Nfs;
print *, "number of physical solutions: ", Nphys;

print *, "flag for debugging mode: ", flag_debug;

print *, "collisions model settings: ", collmod;

end subroutine

!--------------------------------------------------------------------

subroutine clean_flre_data_module ()

use flre_sett;

implicit none;

flre_order = 0;
Nmax = 0;
gal_corr = 0;
Nbmax = 0;
rsp = 0;

hom_sys = 0;
Nwaves = 0;
Nfs = 0;
Nphys = 0;

flag_debug = 0;

collmod = 0;

end subroutine

!--------------------------------------------------------------------

subroutine get_flre_order (flre_order_p)

use flre_sett;

implicit none;

integer, intent(out) :: flre_order_p;

flre_order_p = flre_order;

end subroutine

!--------------------------------------------------------------------

subroutine get_gal_corr (gal_corr_p)

use flre_sett;

implicit none;

integer, intent(out) :: gal_corr_p;

gal_corr_p = gal_corr;

end subroutine

!--------------------------------------------------------------------
