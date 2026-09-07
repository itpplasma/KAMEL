program test_zone_registry
    use, intrinsic :: iso_c_binding, only: c_intptr_t
    use kilca_wave_data_m, only: wave_data_create, wave_data_destroy
    use kilca_zone_m, only: zone_t, handle_to_zone, zone_destroy_c
    use kilca_hmedium_zone_m, only: hmedium_zone_create
    use kilca_imhd_zone_m, only: imhd_zone_create
    use kilca_flre_zone_m, only: flre_zone_create_
    implicit none

    integer(c_intptr_t) :: wave, first, last, transient
    class(zone_t), pointer :: first_zone, last_zone, zone
    integer :: i

    wave = wave_data_create(1, 1, 1.0d0, 0.0d0, 1.0d0, 0.0d0)
    first = hmedium_zone_create(0_c_intptr_t, 1_c_intptr_t, wave, 'first/', -1)
    transient = hmedium_zone_create(0_c_intptr_t, 1_c_intptr_t, wave, 'hole/', 0)
    last = hmedium_zone_create(0_c_intptr_t, 1_c_intptr_t, wave, 'last/', -2)
    call handle_to_zone(first, first_zone)
    call handle_to_zone(last, last_zone)
    call zone_destroy_c(transient)

    ! Repeated transport updates must reuse released entries, including holes
    ! between live zones. Exercise all concrete types and their finalizers.
    do i = 1, 5000
        select case (mod(i, 3))
        case (0)
            transient = hmedium_zone_create(0_c_intptr_t, 1_c_intptr_t, wave, 'new/', i)
        case (1)
            transient = imhd_zone_create(0_c_intptr_t, 1_c_intptr_t, wave, 'new/', i)
        case (2)
            transient = flre_zone_create_(0_c_intptr_t, 1_c_intptr_t, wave, 'new/', i)
        end select
        if (transient == first .or. transient == last) error stop 'Reused a live zone handle'
        call handle_to_zone(transient, zone)
        if (zone%index /= i .or. zone%path /= 'new/') error stop 'Incorrect registered zone'
        if (.not. associated(zone%wd, first_zone%wd)) error stop 'Lost shared wave ownership'
        call zone_destroy_c(transient)

        call handle_to_zone(first, zone)
        if (.not. associated(zone, first_zone)) error stop 'First live handle changed'
        if (zone%index /= -1 .or. zone%path /= 'first/') error stop 'First live zone changed'
        call handle_to_zone(last, zone)
        if (.not. associated(zone, last_zone)) error stop 'Last live handle changed'
        if (zone%index /= -2 .or. zone%path /= 'last/') error stop 'Last live zone changed'
    end do

    call zone_destroy_c(first)
    call zone_destroy_c(last)
    call wave_data_destroy(wave)
    print *, 'PASS: 5000 zone lifecycles preserve both live handles and shared wave ownership'
end program test_zone_registry
