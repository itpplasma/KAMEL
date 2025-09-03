module config_display_m
    ! Module for displaying configuration with formatted tables, colors, and status indicators
    
    use KIM_kinds_m, only: dp
    
    implicit none
    
    private
    public :: display_kim_configuration, display_kim_banner
    
    ! ANSI color codes
    character(len=*), parameter :: RESET = char(27)//'[0m'
    character(len=*), parameter :: BOLD = char(27)//'[1m'
    character(len=*), parameter :: RED = char(27)//'[31m'
    character(len=*), parameter :: GREEN = char(27)//'[32m'
    character(len=*), parameter :: YELLOW = char(27)//'[33m'
    character(len=*), parameter :: BLUE = char(27)//'[34m'
    character(len=*), parameter :: MAGENTA = char(27)//'[35m'
    character(len=*), parameter :: CYAN = char(27)//'[36m'
    character(len=*), parameter :: WHITE = char(27)//'[37m'
    
    ! Box drawing characters (UTF-8)
    character(len=*), parameter :: TL = '╔'  ! Top left
    character(len=*), parameter :: TR = '╗'  ! Top right
    character(len=*), parameter :: BL = '╚'  ! Bottom left
    character(len=*), parameter :: BR = '╝'  ! Bottom right
    character(len=*), parameter :: H = '═'   ! Horizontal
    character(len=*), parameter :: V = '║'   ! Vertical
    character(len=*), parameter :: ML = '╠'  ! Middle left
    character(len=*), parameter :: MR = '╣'  ! Middle right
    character(len=*), parameter :: MC = '╬'  ! Middle cross
    character(len=*), parameter :: MT = '╦'  ! Middle top
    character(len=*), parameter :: MB = '╩'  ! Middle bottom
    
    ! Status indicators
    character(len=*), parameter :: CHECK = '✓'
    character(len=*), parameter :: CROSS = '✗'
    character(len=*), parameter :: ARROW = '▸'
    
contains
    
    subroutine display_kim_banner()
        implicit none
        integer :: i
        integer, parameter :: rep_num = 10
        character(len=100) :: gradient_line
        
        ! Clear screen effect
        write(*,*)
        write(*,*)
        
        ! Top decorative border with gradient effect
        write(*,'(A)') BLUE//BOLD//repeat('░', rep_num)//CYAN//repeat('▒', rep_num)//&
                       MAGENTA//repeat('▓', rep_num)//CYAN//repeat('▒', rep_num)//&
                       BLUE//repeat('░', rep_num)//RESET
        write(*,*)
        
        ! Large KIM ASCII art with gradient colors
        write(*,'(A)') BLUE//BOLD//'         ██╗  ██╗    ██╗    ███╗   ███╗'//RESET
        write(*,'(A)') BLUE//BOLD//'         ██║ ██╔╝    ██║    ████╗ ████║'//RESET
        write(*,'(A)') CYAN//BOLD//'         █████╔╝     ██║    ██╔████╔██║'//RESET  
        write(*,'(A)') CYAN//BOLD//'         ██╔═██╗     ██║    ██║╚██╔╝██║'//RESET
        write(*,'(A)') MAGENTA//BOLD//'         ██║  ██╗    ██║    ██║ ╚═╝ ██║'//RESET
        write(*,'(A)') MAGENTA//BOLD//'         ╚═╝  ╚═╝    ╚═╝    ╚═╝     ╚═╝'//RESET
        
        write(*,*)
        
        ! Title box with double borders
        write(*,'(A)') YELLOW//'     ╔═══════════════════════════════════════╗'//RESET
        write(*,'(A)') YELLOW//'     ║║                                     ║║'//RESET
        write(*,'(A)') YELLOW//'     ║║  '//WHITE//BOLD//'     KiLCA INTEGRAL MODEL'//YELLOW//'          ║║'//RESET
        write(*,'(A)') YELLOW//'     ║║  '//CYAN//'  Kinetic Plasma Response Solver'//YELLOW//'   ║║'//RESET
        write(*,'(A)') YELLOW//'     ║║                                     ║║'//RESET
        write(*,'(A)') YELLOW//'     ╚═══════════════════════════════════════╝'//RESET
        
        write(*,*)
        
        ! Feature highlights with icons
        !write(*,'(A)') GREEN//'              ⚛  '//WHITE//'Fast Integral Formalism'//RESET
        !write(*,'(A)') GREEN//'              ⚡ '//WHITE//'Non-local Plasma Response'//RESET
        !write(*,'(A)') GREEN//'              🔬 '//WHITE//'MHD Perturbation Analysis'//RESET
        !write(*,'(A)') GREEN//'              📊 '//WHITE//'High Performance Computing'//RESET
        
        write(*,*)
        
        ! Bottom decorative border  
        write(*,'(A)') BLUE//BOLD//repeat('░', rep_num)//CYAN//repeat('▒', rep_num)//&
                       MAGENTA//repeat('▓', rep_num)//CYAN//repeat('▒', rep_num)//&
                       BLUE//repeat('░', rep_num)//RESET
        
        write(*,*)
        write(*,*)
        
        ! Startup message
        write(*,'(A)') YELLOW//BOLD//'  ▸ Initializing KIM...'//RESET
        write(*,*)
        
    end subroutine display_kim_banner
    
    subroutine display_kim_configuration()
        use config_m
        use setup_m
        use grid_m
        use species_m, only: plasma
        
        implicit none
        integer :: width = 70
        character(len=50) :: value_str
        
        ! Display header
        call print_header(width)
        
        ! Display KIM_CONFIG section
        call print_section_header('INPUT/OUTPUT CONFIGURATION', width)
        call print_config_line('Profile Location', trim(profile_location), width)
        call print_config_line('Output Path', trim(output_path), width)
        call print_bool_line('HDF5 Input', hdf5_input, width)
        call print_bool_line('HDF5 Output', hdf5_output, width)
        call print_config_line('Debug Level', get_debug_string(fdebug), width)
        call print_config_line('Status Level', get_status_string(fstatus), width)
        
        ! Display Physics Configuration
        call print_section_header('PHYSICS CONFIGURATION', width)
        call print_config_line('Run Type', trim(type_of_run), width)
        call print_config_line('Collision Model', trim(collision_model), width)
        call print_bool_line('Collisions', .not. collisions_off, width)
        call print_bool_line('Artificial Debye Case', artificial_debye_case, width)
        call print_bool_line('Turn Off Ions', turn_off_ions, width)
        call print_bool_line('Turn Off Electrons', turn_off_electrons, width)
        write(value_str, '(I0)') number_of_ion_species
        call print_config_line('Ion Species', trim(value_str), width)
        
        ! Display KIM_SETUP section
        call print_section_header('PLASMA PARAMETERS', width)
        write(value_str, '(F12.2,A)') btor, ' G'
        call print_config_line('Toroidal B Field', trim(adjustl(value_str)), width)
        write(value_str, '(F8.2,A)') r0, ' cm'
        call print_config_line('Major Radius R₀', trim(adjustl(value_str)), width)
        write(value_str, '(F8.2,A)') r_plas, ' cm'
        call print_config_line('Minor Radius', trim(adjustl(value_str)), width)
        write(value_str, '(A,I0,A,I0,A)') '(', m_mode, ', ', n_mode, ')'
        call print_config_line('Mode Numbers (m,n)', trim(value_str), width)
        if (abs(omega) > 1.0e-10) then
            write(value_str, '(ES12.3,A)') omega, ' rad/s'
        else
            value_str = '0 (static)'
        end if
        call print_config_line('Perturbation ω', trim(adjustl(value_str)), width)
        write(value_str, '(ES10.2)') eps_reg
        call print_config_line('Regularization ε', trim(adjustl(value_str)), width)
        
        ! Display Grid Configuration
        call print_section_header('GRID CONFIGURATION', width)
        write(value_str, '(I0)') rg_space_dim
        call print_config_line('rg-space Dimension', trim(value_str), width)
        write(value_str, '(I0)') l_space_dim
        call print_config_line('l-space Dimension', trim(value_str), width)
        call print_config_line('rg Grid Spacing', trim(adjustl(grid_spacing_rg)), width)
        call print_config_line('xl Grid Spacing', trim(adjustl(grid_spacing_xl)), width)
        write(value_str, '(F6.2)') Larmor_skip_factor
        call print_config_line('Larmor Skip Factor', trim(adjustl(value_str)), width)
        call print_config_line('Theta integration method: ', trim(theta_integration), width)
        if (trim(theta_integration) == "RKF45") then
            write(value_str, '(ES10.3)') rkf45_atol
            call print_config_line('RKF45 abs tol', trim(adjustl(value_str)), width)
            write(value_str, '(ES10.3)') rkf45_rtol
            call print_config_line('RKF45 rel tol', trim(adjustl(value_str)), width)
            write(value_str, '(ES10.3)') kernel_taper_skip_threshold
            call print_config_line('Kernel taper skip thr', trim(adjustl(value_str)), width)
            write(value_str, '(A,I0,A,I0,A)') '(', gauss_int_nodes_Nx, ', ', &
                                                    gauss_int_nodes_Nxp, ')'
            call print_config_line('Gauss Nodes (x,x′)', trim(value_str), width)
        else if (trim(theta_integration) == "GaussLegendre") then
            write(value_str, '(A,I0,A,I0,A,I0,A)') '(', gauss_int_nodes_Nx, ', ', &
                                                    gauss_int_nodes_Nxp, ', ', gauss_int_nodes_Ntheta, ')'
            call print_config_line('Gauss Nodes (x,x′,θ)', trim(value_str), width)
        end if
        
        ! Display Ion Species if configured
        if (allocated(plasma%spec)) then
            call print_section_header('ION SPECIES', width)
            call print_species_table(width)
        end if
        
        ! Display footer
        call print_footer(width)
        
    end subroutine display_kim_configuration
    
    subroutine print_header(width)
        implicit none
        integer, intent(in) :: width
        character(len=50) :: title
        integer :: padding
        
        title = 'KIM CONFIGURATION'
        padding = (width - len_trim(title) - 2) / 2
        
        ! Top border
        write(*,'(A,A,A)') CYAN//BOLD, repeat('═', width), RESET
        
        ! Title line
        write(*,'(A,A,A,A,A,A,A)') CYAN//BOLD//V, &
            repeat(' ', padding), WHITE//BOLD, trim(title), CYAN, &
            repeat(' ', width - padding - len_trim(title) - 2), V//RESET
        
        ! Bottom border of header
        write(*,'(A,A,A)') CYAN//BOLD, repeat('═', width), RESET
        
    end subroutine print_header
    
    subroutine print_footer(width)
        implicit none
        integer, intent(in) :: width
        
        write(*,'(A,A,A)') CYAN//BOLD, repeat('═', width), RESET
        
    end subroutine print_footer
    
    subroutine print_section_header(title, width)
        implicit none
        character(len=*), intent(in) :: title
        integer, intent(in) :: width
        integer :: dash_count
        
        dash_count = max(1, width - len_trim(title) - 6)
        
        write(*,*)  ! Empty line before section
        write(*,'(A,A,A,A,A,A,A,A)') YELLOW//BOLD, '  ', ARROW, ' ', &
            trim(title), ' ', repeat('─', dash_count), RESET
        
    end subroutine print_section_header
    
    subroutine print_config_line(label, value, width)
        implicit none
        character(len=*), intent(in) :: label, value
        integer, intent(in) :: width
        integer :: label_width = 25
        integer :: dot_count
        
        dot_count = max(3, label_width - len_trim(label) + 3)
        
        write(*,'(A,A,A,A,A,A,A)') '    ', CYAN, trim(adjustl(label)), &
            repeat('.', dot_count), ' ', WHITE, trim(value)//RESET
        
    end subroutine print_config_line
    
    subroutine print_bool_line(label, value, width)
        implicit none
        character(len=*), intent(in) :: label
        logical, intent(in) :: value
        integer, intent(in) :: width
        integer :: label_width = 25
        integer :: dot_count
        
        dot_count = max(3, label_width - len_trim(label) + 3)
        
        if (value) then
            write(*,'(A,A,A,A,A,A,A,A)') '    ', CYAN, trim(adjustl(label)), &
                repeat('.', dot_count), ' ', GREEN, '['//CHECK//'] Enabled', RESET
        else
            write(*,'(A,A,A,A,A,A,A,A)') '    ', CYAN, trim(adjustl(label)), &
                repeat('.', dot_count), ' ', RED, '['//CROSS//'] Disabled', RESET
        end if
        
    end subroutine print_bool_line
    
    subroutine print_species_table(width)
        use species_m, only: plasma
        implicit none
        integer, intent(in) :: width
        integer :: i
        
        ! Table header
        write(*,'(A,A,A,A,A,A,A,A,A)') '    ', &
            CYAN//BOLD, 'Species', WHITE, '   Z   ', CYAN, '   A   ', WHITE, &
            '   Name'//RESET
        write(*,'(A,A)') '    ', repeat('─', 35)
        
        ! Species data
        do i = 0, plasma%n_species - 1
            write(*,'(A,I5,A,I3,A,I5,A,A15)') '      ', i+1, '    ', &
                plasma%spec(i)%Zspec, '    ', plasma%spec(i)%Aspec, '    ', &
                trim(get_species_name(plasma%spec(i)%Zspec, plasma%spec(i)%Aspec))
        end do
        
    end subroutine print_species_table
    
    ! Helper functions
    
    function get_debug_string(level) result(str)
        implicit none
        integer, intent(in) :: level
        character(len=20) :: str
        
        select case(level)
        case(0)
            str = 'Off'
        case(1)
            str = 'Basic'
        case(2)
            str = 'Detailed'
        case(3)
            str = 'Verbose'
        case default
            write(str, '(A,I0)') 'Level ', level
        end select
        
    end function get_debug_string
    
    function get_status_string(level) result(str)
        implicit none
        integer, intent(in) :: level
        character(len=20) :: str
        
        select case(level)
        case(0)
            str = 'Silent'
        case(1)
            str = 'Normal'
        case(2)
            str = 'Verbose'
        case default
            write(str, '(A,I0)') 'Level ', level
        end select
        
    end function get_status_string
    
    function get_grid_type(gtype) result(str)
        implicit none
        integer, intent(in) :: gtype
        character(len=20) :: str
        
        select case(gtype)
        case(1)
            str = ''
        case(2)
            str = 'Equidistant'
        case(3)
            str = 'Non-equidistant'
        case default
            write(str, '(A,I0)') 'Type ', gtype
        end select
        
    end function get_grid_type
    
    function get_species_name(z, a) result(name)
        implicit none
        integer, intent(in) :: z, a
        character(len=15) :: name
        
        select case(z)
        case(1)
            if (a == 1) then
                name = 'Hydrogen'
            else if (a == 2) then
                name = 'Deuterium'
            else if (a == 3) then
                name = 'Tritium'
            else
                write(name, '(A,I0)') 'H-', a
            end if
        case(2)
            if (a == 4) then
                name = 'Helium-4'
            else if (a == 3) then
                name = 'Helium-3'
            else
                write(name, '(A,I0)') 'He-', a
            end if
        case(6)
            name = 'Carbon'
        case(7)
            name = 'Nitrogen'
        case(8)
            name = 'Oxygen'
        case(10)
            name = 'Neon'
        case(18)
            name = 'Argon'
        case default
            write(name, '(A,I0,A,I0)') 'Z=', z, ' A=', a
        end select
        
    end function get_species_name
    
end module config_display_m
