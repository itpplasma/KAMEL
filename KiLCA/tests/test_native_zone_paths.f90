program test_native_zone_paths
    use, intrinsic :: iso_c_binding, only: c_intptr_t
    use kilca_wave_data_m, only: wave_data_create, wave_data_destroy
    use kilca_zone_m, only: zone_t, handle_to_zone, zone_destroy_c
    use kilca_hmedium_zone_m, only: hmedium_zone_create
    use kilca_imhd_zone_m, only: imhd_zone_create
    use kilca_flre_zone_m, only: flre_zone_create_
    implicit none
    character(len=32) :: storage
    integer(c_intptr_t) :: wave, zones(3)
    class(zone_t), pointer :: zone
    integer :: i

    ! Bytes after the substring are deliberately nonzero: a C-string scan would
    ! copy them into the path instead of respecting the native character length.
    storage = 'short/'//repeat('x', 26)
    wave = wave_data_create(1, 1, 1.0d0, 0.0d0, 1.0d0, 0.0d0)
    zones(1) = hmedium_zone_create(0_c_intptr_t, 1_c_intptr_t, wave, storage(:6), 0)
    zones(2) = imhd_zone_create(0_c_intptr_t, 1_c_intptr_t, wave, storage(:6), 1)
    zones(3) = flre_zone_create_(0_c_intptr_t, 1_c_intptr_t, wave, storage(:6), 2)
    do i = 1, size(zones)
        call handle_to_zone(zones(i), zone)
        if (zone%path /= 'short/') error stop 'Zone factory ignored native path length'
        if (zone%index /= i - 1) error stop 'Zone factory changed the zone index'
        if (zone%wd%m /= 1) error stop 'Zone factory lost the shared wave description'
        call zone_destroy_c(zones(i))
    end do
    call wave_data_destroy(wave)
    print *, 'PASS: native paths and shared wave ownership in all zone factories'
end program test_native_zone_paths
