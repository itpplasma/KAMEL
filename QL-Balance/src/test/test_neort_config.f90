program test_neort_config
    use iso_fortran_env, only: dp => real64
    use neort_interface, only: meta_config_neort_t, read_neort_meta_config, &
                               meta_config_neort_problem
    use wave_code_data, only: m_vals, n_vals

    implicit none

    character(len=*), parameter :: nml_file = "test_neort_config.nml"
    type(meta_config_neort_t) :: meta_config
    integer :: num_failed, iunit
    logical :: unit_open

    num_failed = 0
    allocate (m_vals(1), n_vals(1))
    m_vals = [-6]
    n_vals = [2]

    ! A namelist that only names the input files must fall back to NEO-RT's defaults.
    open (newunit=iunit, file=nml_file, status="replace", action="write")
    write (iunit, '(A)') "&NEORT"
    write (iunit, '(A)') "  boozer_file = 'equil.bc'"
    write (iunit, '(A)') "  boozer_pert_file = 'pert.bc'"
    write (iunit, '(A)') "/"
    close (iunit)

    call read_neort_meta_config(nml_file, meta_config)

    call assert_true(meta_config%config%pertfile, "pertfile set from pert file name")
    call assert_true(meta_config%config%comptorque, "comptorque forced on")
    call assert_true(.not. meta_config%config%nonlin, "nonlin defaults off")
    call assert_true(meta_config%config%m0 == 6, "m0 from mode list")
    call assert_true(meta_config%config%mph == 2, "mph from mode list")
    call assert_true(meta_config%config%epsmn == 1.0_dp, "epsmn NEO-RT default")
    call assert_true(meta_config%config%vmax_over_vth == 4.0_dp, "vmax_over_vth default")
    call assert_true(meta_config%amount_of_s == 100, "amount_of_s default")
    ! antenna_factor is already the squared field ratio, quasilinear torque ~ dB^2
    call assert_true(meta_config%antenna_factor_exponent == 1.0_dp, &
                     "antenna_factor_exponent defaults to quasilinear scaling")
    call assert_true(len_trim(meta_config_neort_problem(meta_config)) == 0, &
                     "complete config has no problem")

    inquire (unit=22, opened=unit_open)
    call assert_true(.not. unit_open, "legacy unit 22 not left open")

    ! Validation: a missing perturbation file and a negative exponent are rejected.
    meta_config%boozer_pert_file = ""
    call assert_true(len_trim(meta_config_neort_problem(meta_config)) > 0, &
                     "missing pert file reported")
    meta_config%boozer_pert_file = "pert.bc"
    meta_config%antenna_factor_exponent = -1.0_dp
    call assert_true(len_trim(meta_config_neort_problem(meta_config)) > 0, &
                     "negative exponent reported")

    open (newunit=iunit, file=nml_file, status="old")
    close (iunit, status="delete")

    if (num_failed > 0) then
        print '(A,I0)', "FAILED assertions: ", num_failed
        stop 1
    end if
    print *, "ALL TESTS PASSED"

contains

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(A,A)', "FAIL: ", label
            num_failed = num_failed + 1
        end if
    end subroutine assert_true

end program test_neort_config
