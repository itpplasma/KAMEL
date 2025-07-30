!> KiLCA Fortran Main Program
!! Fortran implementation of the main KiLCA driver program
!! Equivalent to the C++ main_linear.cpp functionality
!!
!! This program:
!! 1. Parses command line arguments for project path
!! 2. Initializes core data structures (settings, background, modes)
!! 3. Performs mode-independent calculations
!! 4. Performs mode-dependent calculations (antenna or eigmode)
!! 5. Cleans up and exits
!!
!! Usage:
!!   kilca_main [project_path]
!!
!! If no project_path is provided, uses current directory
program kilca_main
    use iso_fortran_env, only: real64, error_unit
    use kilca_types_m
    use kilca_settings_m
    use kilca_background_m
    use kilca_mode_m
    use kilca_core_m
    use kilca_main_utils_m
    implicit none
    
    ! Main program variables
    character(len=MAX_PATH_LEN) :: project_path
    type(core_data_t), pointer :: core_data => null()
    integer :: argc, ierr
    character(len=MAX_PATH_LEN), allocatable :: argv(:)
    
    ! Get command line arguments
    call get_command_line_arguments(argc, argv)
    
    ! Parse arguments to get project path
    call kilca_main_parse_arguments(argc, argv, project_path, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A)') "Error: Failed to parse command line arguments"
        stop 1
    end if
    
    ! Validate and fix project path format
    call kilca_main_validate_project_path(project_path, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A)') "Error: Invalid project path: " // trim(project_path)
        stop 1
    end if
    
    write(*, '(A)') "KiLCA Fortran Main Program"
    write(*, '(A)') "Project path: " // trim(project_path)
    
    ! Initialize core data structure
    call kilca_main_initialize_core_data(core_data, project_path, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A)') "Error: Failed to initialize core data"
        stop 1
    end if
    
    ! Perform mode-independent calculations
    write(*, '(A)') "Performing mode-independent calculations..."
    call kilca_main_calc_mode_independent_data(core_data, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A)') "Error: Mode-independent calculations failed"
        call kilca_main_cleanup_core_data(core_data, ierr)
        stop 1
    end if
    
    ! Check eigmode flag and perform appropriate calculations
    if (associated(core_data%sd)) then
        if (core_data%sd%as%flag_eigmode == 0) then
            ! Normal operation with antenna
            write(*, '(A)') "Performing mode-dependent antenna calculations..."
            call kilca_main_calc_mode_dependent_data_antenna(core_data, ierr)
            if (ierr /= 0) then
                write(error_unit, '(A)') "Error: Mode-dependent antenna calculations failed"
                call kilca_main_cleanup_core_data(core_data, ierr)
                stop 1
            end if
        else
            ! Special operation for instability hunting
            write(*, '(A)') "Performing mode-dependent eigmode calculations..."
            call kilca_main_calc_mode_dependent_data_eigmode(core_data, ierr)
            if (ierr /= 0) then
                write(error_unit, '(A)') "Error: Mode-dependent eigmode calculations failed"
                call kilca_main_cleanup_core_data(core_data, ierr)
                stop 1
            end if
        end if
    else
        write(error_unit, '(A)') "Error: Settings not properly initialized"
        call kilca_main_cleanup_core_data(core_data, ierr)
        stop 1
    end if
    
    write(*, '(A)') "KiLCA calculations completed successfully"
    
    ! Clean up memory
    call kilca_main_cleanup_core_data(core_data, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A)') "Warning: Cleanup had errors"
    end if
    
    ! Clean up command line arguments
    if (allocated(argv)) deallocate(argv)
    
    write(*, '(A)') "Program completed successfully"

contains

    !> Get command line arguments
    subroutine get_command_line_arguments(argc, argv)
        integer, intent(out) :: argc
        character(len=MAX_PATH_LEN), allocatable, intent(out) :: argv(:)
        integer :: i
        
        argc = command_argument_count()
        
        if (argc > 0) then
            allocate(argv(argc))
            do i = 1, argc
                call get_command_argument(i, argv(i))
            end do
        else
            allocate(argv(0))  ! Empty array
        end if
    end subroutine get_command_line_arguments

end program kilca_main