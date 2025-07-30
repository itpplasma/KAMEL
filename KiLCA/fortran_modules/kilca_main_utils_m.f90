!> @file kilca_main_utils_m.f90
!> @brief Utility functions for KiLCA main program
!> @details This module provides shared utility functions used by
!>          the KiLCA main program and its tests.

module kilca_main_utils_m
    use iso_fortran_env, only: real64, error_unit
    use kilca_types_m
    use kilca_settings_m
    use kilca_background_m
    use kilca_mode_m
    use kilca_core_m
    implicit none
    private
    
    ! Public procedures
    public :: kilca_main_parse_arguments
    public :: kilca_main_validate_project_path
    public :: kilca_main_initialize_core_data
    public :: kilca_main_calc_mode_independent_data
    public :: kilca_main_calc_mode_dependent_data_antenna
    public :: kilca_main_calc_mode_dependent_data_eigmode
    public :: kilca_main_cleanup_core_data
    
contains

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
        
        ! Try to get current directory using intrinsic if available
        call get_environment_variable('PWD', path)
        if (len_trim(path) == 0) then
            path = "."
        end if
    end subroutine getcwd

end module kilca_main_utils_m