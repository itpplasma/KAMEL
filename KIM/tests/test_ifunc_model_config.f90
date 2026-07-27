program test_ifunc_model_config

    use config_m, only: IFUNC_MODEL_INHERIT, ifunc_model_for_species, &
                        ifunc_model_is_valid, resolve_ifunc_model, &
                        resolved_electron_ifunc_conservation_model, &
                        resolved_ion_ifunc_conservation_model

    implicit none

    integer :: model

    call assert_equal(resolve_ifunc_model(IFUNC_MODEL_INHERIT, .false.), 0, &
        'legacy energy-off fallback')
    call assert_equal(resolve_ifunc_model(IFUNC_MODEL_INHERIT, .true.), 1, &
        'legacy energy-on fallback')

    do model = 0, 3
        call assert_equal(resolve_ifunc_model(model, .false.), model, &
            'explicit model with legacy energy off')
        call assert_equal(resolve_ifunc_model(model, .true.), model, &
            'explicit model with legacy energy on')
        if (.not. ifunc_model_is_valid(model)) then
            print *, 'FAIL: rejected valid model ', model
            error stop
        end if
    end do

    if (.not. ifunc_model_is_valid(IFUNC_MODEL_INHERIT)) then
        error stop 'FAIL: rejected inherit sentinel'
    end if
    if (ifunc_model_is_valid(-2)) error stop 'FAIL: accepted model -2'
    if (ifunc_model_is_valid(4)) error stop 'FAIL: accepted model 4'

    resolved_electron_ifunc_conservation_model = 1
    resolved_ion_ifunc_conservation_model = 3
    call assert_equal(ifunc_model_for_species(0), 1, 'electron species routing')
    call assert_equal(ifunc_model_for_species(1), 3, 'first ion species routing')
    call assert_equal(ifunc_model_for_species(2), 3, 'second ion species routing')

    print *, 'PASS: per-species I-function model configuration resolves correctly'

contains

    subroutine assert_equal(actual, expected, label)

        integer, intent(in) :: actual, expected
        character(*), intent(in) :: label

        if (actual /= expected) then
            print *, 'FAIL: ', trim(label), ': expected ', expected, ', got ', actual
            error stop
        end if

    end subroutine assert_equal

end program test_ifunc_model_config
