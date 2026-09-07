program test_flr2_factory
    use kim_base_m, only: kim_t
    use kim_mod_m, only: from_kim_factory_get_kim
    use rt_flr2_m, only: flr2_t
    implicit none

    class(kim_t), allocatable :: run

    call from_kim_factory_get_kim('flr2', run)
    select type (run)
        type is (flr2_t)
            write(*,*) 'FLR2 factory test PASSED'
        class default
            write(*,*) 'FAIL: factory did not create flr2_t'
            stop 1
    end select
end program test_flr2_factory
