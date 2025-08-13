program demo_display_simple
    implicit none
    
    ! ANSI color codes
    character(len=10) :: RESET = char(27)//'[0m'
    character(len=10) :: BOLD = char(27)//'[1m'
    character(len=10) :: RED = char(27)//'[31m'
    character(len=10) :: GREEN = char(27)//'[32m'
    character(len=10) :: YELLOW = char(27)//'[33m'
    character(len=10) :: CYAN = char(27)//'[36m'
    character(len=10) :: WHITE = char(27)//'[37m'
    
    ! Display header
    print *, CYAN//BOLD//'══════════════════════════════════════════════════════════'//RESET
    print *, CYAN//BOLD//'║                   '//WHITE//BOLD//'KIM CONFIGURATION'//CYAN//'                   ║'//RESET
    print *, CYAN//BOLD//'══════════════════════════════════════════════════════════'//RESET
    
    ! Display section
    print *, ''
    print *, YELLOW//BOLD//'  ▸ INPUT/OUTPUT CONFIGURATION ─────────────────'//RESET
    print *, '    '//CYAN//'Profile Location........... '//WHITE//'./profiles/'//RESET
    print *, '    '//CYAN//'HDF5 Input................. '//GREEN//'[✓] Enabled'//RESET
    print *, '    '//CYAN//'HDF5 Output................ '//RED//'[✗] Disabled'//RESET
    print *, '    '//CYAN//'Debug Level................ '//WHITE//'Basic'//RESET
    
    ! Physics section
    print *, ''
    print *, YELLOW//BOLD//'  ▸ PLASMA PARAMETERS ───────────────────────────'//RESET
    print *, '    '//CYAN//'Toroidal B Field........... '//WHITE//'-17977.41 G'//RESET
    print *, '    '//CYAN//'Major Radius R₀............ '//WHITE//'165.00 cm'//RESET
    print *, '    '//CYAN//'Mode Numbers (m,n)......... '//WHITE//'(-6, 2)'//RESET
    print *, '    '//CYAN//'Collisions................. '//GREEN//'[✓] Enabled'//RESET
    
    ! Grid section
    print *, ''
    print *, YELLOW//BOLD//'  ▸ GRID CONFIGURATION ──────────────────────────'//RESET
    print *, '    '//CYAN//'k-space Dimension.......... '//WHITE//'100'//RESET
    print *, '    '//CYAN//'l-space Dimension.......... '//WHITE//'1000'//RESET
    print *, '    '//CYAN//'Grid Spacing............... '//WHITE//'Adaptive'//RESET
    
    ! Species table
    print *, ''
    print *, YELLOW//BOLD//'  ▸ ION SPECIES ─────────────────────────────────'//RESET
    print *, '    '//CYAN//BOLD//'Species   Z     A     Name'//RESET
    print *, '    ────────────────────────────'
    print *, '        1     1     2     Deuterium'
    print *, '        2     2     4     Helium-4'
    
    ! Footer
    print *, ''
    print *, CYAN//BOLD//'══════════════════════════════════════════════════════════'//RESET
    
end program demo_display_simple