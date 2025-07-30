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

    !> Parse command line arguments to extract project path
    subroutine kilca_main_parse_arguments(argc, argv, project_path, ierr)
        implicit none
        
        integer, intent(in) :: argc
        character(len=*), intent(in) :: argv(*)
        character(len=*), intent(out) :: project_path
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (argc == 0) then
            ! No arguments - use current directory
            call getcwd(project_path)
            if (len_trim(project_path) == 0) then
                project_path = "."
            end if
        else if (argc >= 1) then
            ! Use first argument as project path
            project_path = trim(argv(1))
        end if
        
    end subroutine kilca_main_parse_arguments

    !> Validate and fix project path format
    subroutine kilca_main_validate_project_path(project_path, ierr)
        implicit none
        
        character(len=*), intent(inout) :: project_path
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Check for empty path
        if (len_trim(project_path) == 0) then
            ierr = -1
            return
        end if
        
        ! Ensure path ends with slash
        if (project_path(len_trim(project_path):len_trim(project_path)) /= '/') then
            project_path = trim(project_path) // '/'
        end if
        
    end subroutine kilca_main_validate_project_path

    !> Initialize core data structure
    subroutine kilca_main_initialize_core_data(core_data, project_path, ierr)
        implicit none
        
        type(core_data_t), pointer, intent(out) :: core_data
        character(len=*), intent(in) :: project_path
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Initialize core data structure
        call core_data_create(core_data, project_path, ierr)
        if (ierr /= 0) return
        
    end subroutine kilca_main_initialize_core_data

    !> Perform mode-independent calculations
    subroutine kilca_main_calc_mode_independent_data(core_data, ierr)
        implicit none
        
        type(core_data_t), pointer, intent(inout) :: core_data
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Call core module function equivalent to calc_and_set_mode_independent_core_data
        call core_calc_mode_independent_data(core_data, ierr)
        
    end subroutine kilca_main_calc_mode_independent_data

    !> Perform mode-dependent antenna calculations
    subroutine kilca_main_calc_mode_dependent_data_antenna(core_data, ierr)
        implicit none
        
        type(core_data_t), pointer, intent(inout) :: core_data
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Call core module function equivalent to calc_and_set_mode_dependent_core_data_antenna
        call core_calc_mode_dependent_data_antenna(core_data, ierr)
        
    end subroutine kilca_main_calc_mode_dependent_data_antenna

    !> Perform mode-dependent eigmode calculations
    subroutine kilca_main_calc_mode_dependent_data_eigmode(core_data, ierr)
        implicit none
        
        type(core_data_t), pointer, intent(inout) :: core_data
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Call core module function equivalent to calc_and_set_mode_dependent_core_data_eigmode
        call core_calc_mode_dependent_data_eigmode(core_data, ierr)
        
    end subroutine kilca_main_calc_mode_dependent_data_eigmode

    !> Clean up core data and free memory
    subroutine kilca_main_cleanup_core_data(core_data, ierr)
        implicit none
        
        type(core_data_t), pointer, intent(inout) :: core_data
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Clean up core data structure
        if (associated(core_data)) then
            call core_data_destroy(core_data, ierr)
            nullify(core_data)
        end if
        
    end subroutine kilca_main_cleanup_core_data

    !> Get current working directory (simple implementation)
    subroutine getcwd(path)
        character(len=*), intent(out) :: path
        integer :: ierr
        
        ! Try to get current directory using intrinsic if available
        call get_environment_variable('PWD', path)
        if (len_trim(path) == 0) then
            path = "."
        end if
    end subroutine getcwd

end program kilca_main