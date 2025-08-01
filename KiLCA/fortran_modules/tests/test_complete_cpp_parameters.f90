program test_complete_cpp_parameters
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: settings_t, settings_read_namelist
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    
    print *, "==========================================="
    print *, "Testing ALL C++ Settings Parameters [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_complete_cpp_parameter_reading()
    
    if (test_status == 0) then
        print *, ""
        print *, "Complete C++ parameter tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Complete C++ parameter tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test reading of ALL C++ parameters from namelist
    subroutine test_complete_cpp_parameter_reading()
        print *, "Testing complete C++ parameter reading functionality..."
        
        ! Create comprehensive test configuration file
        call create_complete_cpp_config()
        
        ! This should FAIL - not all parameters implemented yet
        print *, "Attempting to read ALL C++ parameters from namelist..."
        call settings_read_namelist("complete_cpp_config.conf", settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read complete C++ parameters, error:", ierr
            test_status = test_status + 1
            return
        else
            print *, "SUCCESS: Complete C++ parameters read correctly"
        end if
        
        ! Test EVERY SINGLE parameter from C++ implementation
        call validate_all_antenna_parameters()
        call validate_all_background_parameters()
        call validate_all_output_parameters()
        call validate_all_eigenmode_parameters()
        
    end subroutine test_complete_cpp_parameter_reading
    
    !> Validate ALL antenna parameters from C++ antenna_t
    subroutine validate_all_antenna_parameters()
        print *, "Validating ALL antenna parameters..."
        
        ! Basic parameters (already implemented)
        if (abs(settings%antenna_settings%ra - 85.5_dp) > 1.0e-12_dp) then
            print *, "FAIL: Antenna ra not read correctly:", settings%antenna_settings%ra
            test_status = test_status + 1
        else
            print *, "PASS: Antenna ra parameter"
        end if
        
        if (abs(settings%antenna_settings%wa - 2.5_dp) > 1.0e-12_dp) then
            print *, "FAIL: Antenna wa not read correctly:", settings%antenna_settings%wa
            test_status = test_status + 1
        else
            print *, "PASS: Antenna wa parameter"
        end if
        
        if (abs(settings%antenna_settings%I0 - 5.0e12_dp) > 1.0e-12_dp) then
            print *, "FAIL: Antenna I0 not read correctly:", settings%antenna_settings%I0
            test_status = test_status + 1
        else
            print *, "PASS: Antenna I0 parameter"
        end if
        
        if (abs(settings%antenna_settings%flab - cmplx(1.5e6_dp, 0.1e6_dp, dp)) > 1.0e-12_dp) then
            print *, "FAIL: Antenna flab not read correctly:", settings%antenna_settings%flab
            test_status = test_status + 1
        else
            print *, "PASS: Antenna flab parameter"
        end if
        
        if (settings%antenna_settings%dma /= 3) then
            print *, "FAIL: Antenna dma not read correctly:", settings%antenna_settings%dma
            test_status = test_status + 1
        else
            print *, "PASS: Antenna dma parameter"
        end if
        
        ! Missing parameters that should fail
        if (settings%antenna_settings%flag_eigmode /= 1) then
            print *, "FAIL: Antenna flag_eigmode not read correctly:", settings%antenna_settings%flag_eigmode
            test_status = test_status + 1
        else
            print *, "PASS: Antenna flag_eigmode parameter"
        end if
        
        ! Modes array (missing - should fail)
        if (.not. allocated(settings%antenna_settings%modes)) then
            print *, "FAIL: Antenna modes array not allocated"
            test_status = test_status + 1
        else if (size(settings%antenna_settings%modes) /= 6) then
            print *, "FAIL: Antenna modes array wrong size:", size(settings%antenna_settings%modes)
            test_status = test_status + 1
        else if (settings%antenna_settings%modes(1) /= 1 .or. settings%antenna_settings%modes(6) /= 3) then
            print *, "FAIL: Antenna modes array values incorrect"
            test_status = test_status + 1
        else
            print *, "PASS: Antenna modes array"
        end if
        
    end subroutine validate_all_antenna_parameters
    
    !> Validate ALL background parameters from C++ back_sett_t  
    subroutine validate_all_background_parameters()
        print *, "Validating ALL background parameters..."
        
        ! Basic parameters (already implemented)
        if (abs(settings%background_settings%rtor - 175.0_dp) > 1.0e-12_dp) then
            print *, "FAIL: Background rtor not read correctly:", settings%background_settings%rtor
            test_status = test_status + 1
        else
            print *, "PASS: Background rtor parameter"
        end if
        
        ! Missing string parameters that should fail
        if (.not. allocated(settings%background_settings%path2profiles)) then
            print *, "FAIL: Background path2profiles not allocated"
            test_status = test_status + 1
        else if (settings%background_settings%path2profiles /= "../test_profiles/") then
            print *, "FAIL: Background path2profiles incorrect:", settings%background_settings%path2profiles
            test_status = test_status + 1
        else
            print *, "PASS: Background path2profiles parameter"
        end if
        
        if (.not. allocated(settings%background_settings%flag_back)) then
            print *, "FAIL: Background flag_back not allocated"
            test_status = test_status + 1
        else if (settings%background_settings%flag_back /= "normal") then
            print *, "FAIL: Background flag_back incorrect:", settings%background_settings%flag_back
            test_status = test_status + 1
        else
            print *, "PASS: Background flag_back parameter"
        end if
        
        ! Missing numeric parameters that should fail
        if (settings%background_settings%N /= 7) then
            print *, "FAIL: Background N not read correctly:", settings%background_settings%N
            test_status = test_status + 1
        else
            print *, "PASS: Background N parameter"
        end if
        
        if (abs(settings%background_settings%V_scale - 1.2_dp) > 1.0e-12_dp) then
            print *, "FAIL: Background V_scale not read correctly:", settings%background_settings%V_scale
            test_status = test_status + 1
        else
            print *, "PASS: Background V_scale parameter"
        end if
        
        if (abs(settings%background_settings%zele - 0.8_dp) > 1.0e-12_dp) then
            print *, "FAIL: Background zele not read correctly:", settings%background_settings%zele
            test_status = test_status + 1
        else
            print *, "PASS: Background zele parameter"
        end if
        
        if (abs(settings%background_settings%zion - 1.1_dp) > 1.0e-12_dp) then
            print *, "FAIL: Background zion not read correctly:", settings%background_settings%zion
            test_status = test_status + 1
        else
            print *, "PASS: Background zion parameter"
        end if
        
    end subroutine validate_all_background_parameters
    
    !> Validate ALL output parameters from C++ output_sett_t
    subroutine validate_all_output_parameters()
        print *, "Validating ALL output parameters..."
        
        ! Basic parameters already implemented should work
        if (settings%output_settings%flag_background /= 1) then
            print *, "FAIL: Output flag_background not read correctly:", settings%output_settings%flag_background
            test_status = test_status + 1
        else
            print *, "PASS: Output flag_background parameter"
        end if
        
        ! Missing flag_debug parameter that should fail
        if (settings%output_settings%flag_debug /= 1) then
            print *, "FAIL: Output flag_debug not read correctly:", settings%output_settings%flag_debug
            test_status = test_status + 1
        else
            print *, "PASS: Output flag_debug parameter"
        end if
        
    end subroutine validate_all_output_parameters
    
    !> Validate ALL eigenmode parameters from C++ eigmode_sett_t
    subroutine validate_all_eigenmode_parameters()
        print *, "Validating ALL eigenmode parameters..."
        
        ! Basic parameters already implemented should work
        if (settings%eigmode_settings%search_flag /= 2) then
            print *, "FAIL: Eigenmode search_flag not read correctly:", settings%eigmode_settings%search_flag
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode search_flag parameter"
        end if
        
        ! Missing string parameter that should fail
        if (.not. allocated(settings%eigmode_settings%fname)) then
            print *, "FAIL: Eigenmode fname not allocated"
            test_status = test_status + 1
        else if (settings%eigmode_settings%fname /= "test_roots.dat") then
            print *, "FAIL: Eigenmode fname incorrect:", settings%eigmode_settings%fname
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode fname parameter"
        end if
        
        ! Missing numeric parameters that should fail
        if (settings%eigmode_settings%stop_flag /= 1) then
            print *, "FAIL: Eigenmode stop_flag not read correctly:", settings%eigmode_settings%stop_flag
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode stop_flag parameter"
        end if
        
        if (abs(settings%eigmode_settings%eps_res - 1.5e-7_dp) > 1.0e-15_dp) then
            print *, "FAIL: Eigenmode eps_res not read correctly:", settings%eigmode_settings%eps_res
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode eps_res parameter"
        end if
        
        if (abs(settings%eigmode_settings%eps_abs - 2.0e-9_dp) > 1.0e-15_dp) then
            print *, "FAIL: Eigenmode eps_abs not read correctly:", settings%eigmode_settings%eps_abs
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode eps_abs parameter"
        end if
        
        if (abs(settings%eigmode_settings%eps_rel - 1.2e-7_dp) > 1.0e-15_dp) then
            print *, "FAIL: Eigenmode eps_rel not read correctly:", settings%eigmode_settings%eps_rel
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode eps_rel parameter"
        end if
        
        if (abs(settings%eigmode_settings%delta - 5.0e-7_dp) > 1.0e-15_dp) then
            print *, "FAIL: Eigenmode delta not read correctly:", settings%eigmode_settings%delta
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode delta parameter"
        end if
        
        if (settings%eigmode_settings%test_roots /= 1) then
            print *, "FAIL: Eigenmode test_roots not read correctly:", settings%eigmode_settings%test_roots
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode test_roots parameter"
        end if
        
        if (settings%eigmode_settings%Nguess /= 2) then
            print *, "FAIL: Eigenmode Nguess not read correctly:", settings%eigmode_settings%Nguess
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode Nguess parameter"
        end if
        
        if (settings%eigmode_settings%kmin /= 2) then
            print *, "FAIL: Eigenmode kmin not read correctly:", settings%eigmode_settings%kmin
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode kmin parameter"
        end if
        
        if (settings%eigmode_settings%kmax /= 8) then
            print *, "FAIL: Eigenmode kmax not read correctly:", settings%eigmode_settings%kmax
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode kmax parameter"
        end if
        
        if (settings%eigmode_settings%n_zeros /= 5) then
            print *, "FAIL: Eigenmode n_zeros not read correctly:", settings%eigmode_settings%n_zeros
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode n_zeros parameter"
        end if
        
        if (settings%eigmode_settings%use_winding /= 1) then
            print *, "FAIL: Eigenmode use_winding not read correctly:", settings%eigmode_settings%use_winding
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode use_winding parameter"
        end if
        
        ! Complex array parameter that should fail
        if (.not. allocated(settings%eigmode_settings%fstart)) then
            print *, "FAIL: Eigenmode fstart array not allocated"
            test_status = test_status + 1
        else if (size(settings%eigmode_settings%fstart) /= 2) then
            print *, "FAIL: Eigenmode fstart array wrong size:", size(settings%eigmode_settings%fstart)
            test_status = test_status + 1
        else if (abs(settings%eigmode_settings%fstart(1) - cmplx(1.0e6_dp, 0.1e6_dp, dp)) > 1.0e-12_dp) then
            print *, "FAIL: Eigenmode fstart(1) value incorrect:", settings%eigmode_settings%fstart(1)
            test_status = test_status + 1
        else
            print *, "PASS: Eigenmode fstart array"
        end if
        
    end subroutine validate_all_eigenmode_parameters
    
    !> Create comprehensive test configuration file with ALL C++ parameters
    subroutine create_complete_cpp_config()
        integer :: unit
        integer :: iostat
        
        open(newunit=unit, file="complete_cpp_config.conf", status="replace", action="write", iostat=iostat)
        if (iostat /= 0) then
            print *, "ERROR: Could not create complete C++ config file, iostat =", iostat
            return
        end if
        
        write(unit, '(a)') "! Complete C++ Parameter Configuration - ALL parameters from C++ types"
        write(unit, '(a)') ""
        
        ! Antenna section with ALL parameters
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 85.5"
        write(unit, '(a)') "  wa = 2.5"
        write(unit, '(a)') "  I0 = 5.0e12"
        write(unit, '(a)') "  flab = (1.5e6, 0.1e6)"
        write(unit, '(a)') "  dma = 3"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1, 3, 3"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        
        ! Background section with ALL parameters
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 175.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "  path2profiles = '../test_profiles/'"
        write(unit, '(a)') "  calc_back = 1"
        write(unit, '(a)') "  flag_back = 'normal'"
        write(unit, '(a)') "  N = 7"
        write(unit, '(a)') "  V_gal_sys = 2.0e9"
        write(unit, '(a)') "  V_scale = 1.2"
        write(unit, '(a)') "  m_i = 2.5"
        write(unit, '(a)') "  zele = 0.8"
        write(unit, '(a)') "  zion = 1.1"
        write(unit, '(a)') "  flag_debug_bg = 1"
        write(unit, '(a)') "  mass = 2.0, 1.0, 4.0"
        write(unit, '(a)') "  charge = 1.0, -1.0, 2.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        
        ! Output section with ALL parameters
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "  flag_emfield = 2"
        write(unit, '(a)') "  flag_additional = 1"
        write(unit, '(a)') "  flag_dispersion = 0"
        write(unit, '(a)') "  flag_debug_out = 1"
        write(unit, '(a)') "  num_quants = 6"
        write(unit, '(a)') "  flag_quants = 1, 0, 1, 1, 0, 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        
        ! Eigenmode section with ALL parameters
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  fname = 'test_roots.dat'"
        write(unit, '(a)') "  search_flag = 2"
        write(unit, '(a)') "  rdim = 150"
        write(unit, '(a)') "  rfmin = 1.0e5"
        write(unit, '(a)') "  rfmax = 2.0e6"
        write(unit, '(a)') "  idim = 80"
        write(unit, '(a)') "  ifmin = -5.0e4"
        write(unit, '(a)') "  ifmax = 5.0e4"
        write(unit, '(a)') "  stop_flag = 1"
        write(unit, '(a)') "  eps_res = 1.5e-7"
        write(unit, '(a)') "  eps_abs = 2.0e-9"
        write(unit, '(a)') "  eps_rel = 1.2e-7"
        write(unit, '(a)') "  delta = 5.0e-7"
        write(unit, '(a)') "  test_roots = 1"
        write(unit, '(a)') "  flag_debug_eig = 0"
        write(unit, '(a)') "  Nguess = 2"
        write(unit, '(a)') "  kmin = 2"
        write(unit, '(a)') "  kmax = 8"
        write(unit, '(a)') "  n_zeros = 5"
        write(unit, '(a)') "  use_winding = 1"
        write(unit, '(a)') "  fstart = (1.0e6, 0.1e6), (1.5e6, -0.05e6)"
        write(unit, '(a)') "/"
        
        close(unit)
        print *, "Complete C++ config file created successfully"
    end subroutine create_complete_cpp_config

end program test_complete_cpp_parameters