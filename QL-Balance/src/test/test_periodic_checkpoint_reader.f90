program test_periodic_checkpoint_reader
    use QLBalance_kinds, only: dp
    use KAMEL_hdf5_tools
    use h5mod, only: path2inp, path2time, path2out
    use baseparam_mod, only: rtor
    use wave_code_data, only: m_vals, n_vals, rq, iq, rn, in, rTi, iTi, rTe, iTe, &
                              rVth, iVth, rVz, iVz, rep, idPhi0
    use periodic_amplitude_state_m, only: periodic_amplitudes
    use periodic_checkpoint_m, only: pending_periodic_restart
    implicit none
    integer(HID_T) :: file_id
    integer :: case_index
    character(128) :: group
    real(dp), parameter :: radius(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: density(3) = [1.1e13_dp, 1.2e13_dp, 1.3e13_dp]
    real(dp), parameter :: electric(3) = [0.11_dp, -0.22_dp, 0.33_dp]

    path2inp = 'checkpoint_reader_input.h5'
    path2time = 'checkpoint_reader_profiles.h5'
    path2out = 'checkpoint_reader_output.h5'
    rtor = 100.0_dp
    call h5_init()
    call h5_create(trim(path2inp), file_id)
    call h5_add(file_id, '/preprocprof/q', [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [1], [4])
    call h5_add(file_id, '/preprocprof/r_out', [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1], [4])
    call h5_close(file_id)
    call h5_create(trim(path2out), file_id)
    call h5_close(file_id)
    call h5_deinit()

    do case_index = 1, 3
        if (case_index == 1) then
            m_vals = [-6, -7]
            n_vals = [2, 2]
            group = 'multi_mode/KinProfiles/1004/'
        elseif (case_index == 2) then
            m_vals = [-16]
            n_vals = [12]
            group = 'f_-16_12/KinProfiles/1004/'
        else
            m_vals = [-6]
            n_vals = [2]
            group = 'f_*_2/fort.1000/1004/'
        end if
        call h5_init()
        call h5_create(trim(path2time), file_id)
        call h5_add(file_id, trim(group)//'n', density, [1], [3])
        call h5_add(file_id, trim(group)//'Ti', radius * 100.0_dp, [1], [3])
        call h5_add(file_id, trim(group)//'Te', radius * 200.0_dp, [1], [3])
        call h5_add(file_id, trim(group)//'Vz', radius, [1], [3])
        call h5_add(file_id, trim(group)//'rc', radius, [1], [3])
        ! The historical writer saves cell-averaged Er but boundary-point Vth.
        call h5_add(file_id, trim(group)//'Er', electric, [1], [3])
        call h5_add(file_id, trim(group)//'Vth', [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp], &
                    [1], [4])
        call h5_close(file_id)
        call h5_deinit()
        call periodic_amplitudes%initialize([(2.0_dp, 3.0_dp)])
        pending_periodic_restart%active = .true.
        pending_periodic_restart%params = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [4, 1])

        call read_background_profiles_h5_timeevol(4)

        call require(.not. periodic_amplitudes%initialized, &
                     'legacy reader retained stale amplitudes')
        call require(.not. pending_periodic_restart%active, &
                     'legacy reader retained pending restart')
        call require(.not. allocated(pending_periodic_restart%params), &
                     'legacy pending profiles remain')
        call require(all(in == density) .and. all(rn == radius), 'legacy density or grid changed')
        call require(all(idPhi0 == -electric), 'legacy Er sign or cell-length reader is wrong')
        call require(all(iVz == radius * rtor), &
                     'legacy angular to linear velocity conversion changed')
        call require(all(iVth == [10.0_dp, 20.0_dp, 30.0_dp]), 'legacy boundary Vth changed')
        call require(all(iq == [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]), 'legacy input q was not loaded')
        deallocate (rq, iq, rn, in, rTi, iTi, rTe, iTe, rVth, iVth, rVz, iVz, rep, idPhi0)
    end do
    print *, 'periodic checkpoint legacy reader passed'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) error stop message
    end subroutine require
end program test_periodic_checkpoint_reader
