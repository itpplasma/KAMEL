program demo_config_display
    implicit none
    
    ! ANSI color codes
    character(len=*), parameter :: RESET = char(27)//'[0m'
    character(len=*), parameter :: BOLD = char(27)//'[1m'
    character(len=*), parameter :: RED = char(27)//'[31m'
    character(len=*), parameter :: GREEN = char(27)//'[32m'
    character(len=*), parameter :: YELLOW = char(27)//'[33m'
    character(len=*), parameter :: BLUE = char(27)//'[34m'
    character(len=*), parameter :: CYAN = char(27)//'[36m'
    character(len=*), parameter :: WHITE = char(27)//'[37m'
    
    ! Status indicators
    character(len=*), parameter :: CHECK = '✓'
    character(len=*), parameter :: CROSS = '✗'
    character(len=*), parameter :: ARROW = '▸'
    
    integer :: width = 70
    character(len=200) :: line
    
    ! Display header
    write(line, '(A,A,A)') CYAN//BOLD, repeat('═', width), RESET
    write(*,'(A)') trim(line)
    write(*,'(A,A,A,A,A,A,A)') CYAN//BOLD//'║', &
        repeat(' ', 26), WHITE//BOLD, 'KIM CONFIGURATION', CYAN, &
        repeat(' ', 26), '║'//RESET
    write(line, '(A,A,A)') CYAN//BOLD, repeat('═', width), RESET
    write(*,'(A)') trim(line)
    
    ! Display section
    write(*,*)
    write(line, '(A,A,A,A,A,A,A,A)') YELLOW//BOLD, '  ', ARROW, ' ', &
        'INPUT/OUTPUT CONFIGURATION', ' ', repeat('─', 39), RESET
    write(*,'(A)') trim(line)
    
    ! Sample configuration lines
    write(line, '(A,A,A,A,A,A)') '    ', CYAN, 'Profile Location', &
        repeat('.', 12), WHITE, ' ./profiles/'//RESET
    write(*,'(A)') trim(line)
    
    write(line, '(A,A,A,A,A)') '    ', CYAN, 'HDF5 Input', &
        repeat('.', 18), ' '//GREEN//'['//CHECK//'] Enabled'//RESET
    write(*,'(A)') trim(line)
    
    write(line, '(A,A,A,A,A)') '    ', CYAN, 'HDF5 Output', &
        repeat('.', 17), ' '//RED//'['//CROSS//'] Disabled'//RESET
    write(*,'(A)') trim(line)
    
    ! Physics section
    write(*,*)
    write(line, '(A,A,A,A,A,A,A,A)') YELLOW//BOLD, '  ', ARROW, ' ', &
        'PLASMA PARAMETERS', ' ', repeat('─', 47), RESET
    write(*,'(A)') trim(line)
    
    write(line, '(A,A,A,A,A,A)') '    ', CYAN, 'Toroidal B Field', &
        repeat('.', 11), WHITE, ' -17977.41 G'//RESET
    write(*,'(A)') trim(line)
    
    write(line, '(A,A,A,A,A,A)') '    ', CYAN, 'Mode Numbers (m,n)', &
        repeat('.', 9), WHITE, ' (-6, 2)'//RESET
    write(*,'(A)') trim(line)
    
    ! Species table
    write(*,*)
    write(line, '(A,A,A,A,A,A,A,A)') YELLOW//BOLD, '  ', ARROW, ' ', &
        'ION SPECIES', ' ', repeat('─', 53), RESET
    write(*,'(A)') trim(line)
    
    write(line, '(A,A,A,A,A,A,A,A,A)') '    ', &
        CYAN//BOLD, 'Species', WHITE, '   Z   ', CYAN, '   A   ', WHITE, &
        '   Name'//RESET
    write(*,'(A)') trim(line)
    write(*,'(A,A)') '    ', repeat('─', 35)
    write(*,'(A)') '        1      1      2      Deuterium'
    write(*,'(A)') '        2      2      4      Helium-4'
    
    ! Footer
    write(*,*)
    write(line, '(A,A,A)') CYAN//BOLD, repeat('═', width), RESET
    write(*,'(A)') trim(line)
    
end program demo_config_display