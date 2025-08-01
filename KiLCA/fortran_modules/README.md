# KiLCA Fortran Modules

This directory contains the Fortran translation of KiLCA core components.

## Directory Structure

```
fortran_modules/
├── *.f90                 # Source modules
├── Makefile              # Build system
├── .gitignore            # Git ignore patterns
├── tests/                # Test files organized by module
│   ├── complex/          # Complex number tests
│   │   ├── test_kilca_complex*.f90
│   │   └── ...
│   └── test_*.f90        # General module tests
├── build/                # Build artifacts (ignored by git)
│   ├── obj/              # Object files and modules
│   └── tests/            # Compiled test executables
└── examples/             # Usage examples
```

## Building

```bash
# Build all modules
make all

# Run all tests
make test

# Run only complex number tests
make test-complex

# Run only general tests
make test-general

# Clean build artifacts
make clean
```

## Modules

### Core Modules
- `kilca_types_m.f90` - Basic type definitions
- `kilca_constants_m.f90` - Physical and mathematical constants
- `kilca_shared_m.f90` - Shared utility functions
- `kilca_complex_m.f90` - Complex number operations
- `kilca_settings_m.f90` - Configuration management
- `kilca_core_m.f90` - Core data structures

### Complex Number Module Features
- Complete set of complex arithmetic operations
- Transcendental functions (sin, cos, exp, log, etc.)
- Advanced functions (asinh, acosh, expm1, log1p, etc.)
- Safe functions with overflow/underflow protection
- Array and matrix operations
- Comprehensive I/O formatting

## Testing

All modules have comprehensive test suites that follow the RED-GREEN-REFACTOR TDD methodology:

1. **RED Phase**: Write failing tests first
2. **GREEN Phase**: Implement functionality to pass tests
3. **REFACTOR Phase**: Improve code while maintaining tests

Tests are organized by module type:
- `tests/complex/` - Complex number function tests
- `tests/test_*.f90` - Individual module tests

## KiLCA Fortran Main Program Manual

### Overview

The KiLCA Fortran main program (`kilca_main`) is a complete standalone implementation of the KiLCA (Kinetic plAsma response ModEL - Cylindrical) solver. It provides the same functionality as the original C++ `main_linear.cpp` program but is implemented entirely in Fortran for better integration with the existing Fortran codebase.

### Program Flow

The main program follows this execution sequence:

1. **Command Line Parsing**: Reads project path from command line arguments
2. **Project Path Validation**: Ensures the path exists and is properly formatted
3. **Core Data Initialization**: Creates and initializes the main data structures
4. **Mode-Independent Calculations**: Performs background plasma calculations
5. **Mode-Dependent Calculations**: Performs either antenna or eigenmode calculations
6. **Cleanup**: Releases memory and exits

### Usage

#### Basic Usage

```bash
# Run with project directory
./kilca_main /path/to/project

# Run in current directory (uses current directory as project path)
./kilca_main

# Show help (any unrecognized argument shows basic info)
./kilca_main --help
```

#### Building the Main Program

The main program is built automatically when you build the KiLCA Fortran modules:

```bash
# Build all modules including main program
make all

# The executable will be created as:
# build/bin/kilca_main
```

#### Program Options

The program accepts these command line arguments:

- **No arguments**: Uses current working directory as project path
- **Single path argument**: Uses specified directory as project path
- **Any other arguments**: Shows program information and exits with error

### Input Files Required

The KiLCA Fortran main program requires specific input files in the project directory to function correctly. These files configure the plasma parameters, antenna settings, and calculation options.

#### Required Project Directory Structure

```
project_directory/
├── settings.conf          # Main configuration file (or settings.nml)
├── profiles/              # Background plasma profiles
│   ├── density.dat
│   ├── temperature.dat
│   └── magnetic_field.dat
├── antenna/               # Antenna configuration (if flag_eigmode = 0)
│   └── antenna_spec.dat
└── output/                # Output directory (created automatically)
```

#### Main Configuration File: `settings.conf`

This is the primary configuration file that controls all aspects of the KiLCA calculation. The file uses a namelist format or key-value pairs:

**Format Example:**
```fortran
&antenna_settings
  ra = 18.5         ! Antenna radius (cm)
  wa = 0.8          ! Antenna current layer width (cm)  
  I0 = 800.0        ! Antenna current (statamps)
  flab = (30.0e6, 0.0)  ! Antenna frequency (Hz) - complex format
  dma = 2           ! Number of mode pairs
  flag_debug = 1    ! Enable debug output
/

&background_settings
  rtor = 165.0      ! Major radius (cm)
  rp = 50.0         ! Minor radius (cm)
  B0 = 25000.0      ! Toroidal magnetic field (G)
  calc_back = 1     ! Calculate background profiles
  N = 5             ! Spline order for interpolation
  path2profiles = "./profiles/"  ! Path to profile files
/

&output_settings
  flag_background = 1     ! Output background quantities
  flag_emfield = 1        ! Output electromagnetic fields
  flag_additional = 1     ! Output additional diagnostics
  num_quants = 5          ! Number of physical quantities to calculate
/

&eigmode_settings
  flag_eigmode = 0        ! 0 = antenna mode, 1 = eigenmode search
  search_flag = 0         ! Enable eigenmode search
  rdim = 100             ! Real frequency resolution
  idim = 100             ! Imaginary frequency resolution
  rfmin = 20.0e6         ! Minimum real frequency (Hz)
  rfmax = 40.0e6         ! Maximum real frequency (Hz)
  ifmin = -1000.0        ! Minimum imaginary frequency (Hz)
  ifmax = 1000.0         ! Maximum imaginary frequency (Hz)
  n_zeros = 10           ! Maximum number of zeros to find
/
```

**Key Parameters:**

- **Antenna Settings**: Configure the RF antenna
  - `ra`: Antenna radial position (cm)
  - `wa`: Current layer width (cm)
  - `I0`: Antenna current amplitude (statamps)
  - `flab`: Antenna frequency (Hz, complex number)
  - `dma`: Number of (m,n) mode pairs to calculate

- **Background Settings**: Define the plasma background
  - `rtor`: Major radius of the torus (cm)
  - `rp`: Minor radius of the plasma (cm)
  - `B0`: Toroidal magnetic field strength (Gauss)
  - `calc_back`: Flag to calculate background profiles (0/1)
  - `N`: Order of spline interpolation

- **Output Settings**: Control what quantities are saved
  - `flag_background`: Save background plasma quantities
  - `flag_emfield`: Save electromagnetic field profiles
  - `flag_additional`: Save additional diagnostic quantities
  - `num_quants`: Number of physical quantities to compute

- **Eigenmode Settings**: Configure eigenmode analysis
  - `flag_eigmode`: 0 for antenna mode, 1 for eigenmode search
  - `search_flag`: Enable detailed eigenmode search
  - `rdim`, `idim`: Resolution in real/imaginary frequency
  - `rfmin`, `rfmax`: Real frequency search range (Hz)
  - `ifmin`, `ifmax`: Imaginary frequency search range (Hz)

#### Background Plasma Profile Files

Profile files in the `profiles/` directory define the background plasma:

**`density.dat`** - Plasma density profile:
```
# Radius (cm)    Electron density (cm^-3)    Ion density (cm^-3)
0.0              1.0e14                       1.0e14
10.0             8.0e13                       8.0e13
20.0             5.0e13                       5.0e13
30.0             2.0e13                       2.0e13
40.0             1.0e12                       1.0e12
50.0             1.0e11                       1.0e11
```

**`temperature.dat`** - Temperature profiles:
```
# Radius (cm)    Electron temp (keV)    Ion temp (keV)
0.0              5.0                     3.0
10.0             4.5                     2.8
20.0             3.5                     2.2
30.0             2.0                     1.5
40.0             1.0                     0.8
50.0             0.2                     0.2
```

**`magnetic_field.dat`** - Magnetic field configuration:
```
# Radius (cm)    Bpol (G)    Btor (G)    q-factor
0.0              0.0         25000.0     0.8
10.0             1500.0      24800.0     1.1
20.0             2800.0      24200.0     1.8
30.0             3500.0      23500.0     2.8
40.0             4000.0      22500.0     4.2
50.0             4200.0      21000.0     6.8
```

#### Antenna Configuration File: `antenna/antenna_spec.dat`

(Only required when `flag_eigmode = 0`)

```
# Antenna spectrum definition
# Mode number (m,n)    Current amplitude    Phase (radians)
1   15                 1.0                  0.0
1   16                 0.8                  0.5
2   30                 0.6                  1.0
2   32                 0.4                  1.5
```

#### Mode Configuration

For antenna calculations (`flag_eigmode = 0`), specify the toroidal mode numbers in the settings file:

```fortran
&antenna_settings
  dma = 3           ! Number of mode pairs
  ! This creates an array: modes(2*dma) = [m1, n1, m2, n2, m3, n3]
  ! You need to modify the main program or provide a separate modes file
/
```

### Output Files

The program generates several output files in the project directory:

- **`background_profiles.dat`**: Background plasma quantities
- **`electromagnetic_fields.dat`**: RF field profiles
- **`power_deposition.dat`**: Power absorption profiles
- **`eigenmode_results.dat`**: Eigenmode frequencies and growth rates (if `flag_eigmode = 1`)
- **`debug.log`**: Detailed calculation log (if debug flags enabled)

### Error Handling

The program includes comprehensive error handling:

1. **Command Line Errors**: Invalid arguments or missing project path
2. **File System Errors**: Missing project directory or input files
3. **Configuration Errors**: Invalid parameter values or missing settings
4. **Physics Errors**: Unphysical plasma parameters or calculation failures
5. **Memory Errors**: Allocation failures or cleanup issues

**Common Error Messages:**

- `"Error: calc_and_set_mode_independent_core_data: unknown project type"`
  - Check that settings.conf exists and has valid parameters
  - Verify that background profile files are present

- `"Error: Failed to initialize core data"`
  - Project directory may not exist or be readable
  - Settings file may have syntax errors

- `"Error: Mode-independent calculations failed"`
  - Background plasma parameters may be unphysical
  - Profile files may have incorrect format or missing data

### Examples

#### Example 1: Basic ICRF Heating Calculation

```bash
# Set up project directory
mkdir my_icrf_run
cd my_icrf_run

# Create basic settings file for 30 MHz ICRF heating
cat > settings.conf << EOF
&antenna_settings
  ra = 18.5
  wa = 0.8  
  I0 = 800.0
  flab = (30.0e6, 0.0)
  dma = 2
  flag_debug = 1
/
&background_settings
  rtor = 165.0
  rp = 50.0
  B0 = 25000.0
  calc_back = 1
  N = 5
/
&output_settings
  flag_background = 1
  flag_emfield = 1
  flag_additional = 1
  num_quants = 5
/
&eigmode_settings
  flag_eigmode = 0
/
EOF

# Create profile files (simplified example)
mkdir profiles
# ... create density.dat, temperature.dat, magnetic_field.dat ...

# Run calculation
/path/to/kilca_main .
```

#### Example 2: Eigenmode Search

```bash
# Modify settings for eigenmode search
cat > settings.conf << EOF
&antenna_settings
  ra = 18.5
  wa = 0.8
  I0 = 800.0
  flab = (30.0e6, 0.0)
  dma = 1
/
&background_settings
  rtor = 165.0
  rp = 50.0
  B0 = 25000.0
  calc_back = 1
  N = 5
/
&output_settings
  flag_background = 1
  flag_emfield = 1
/
&eigmode_settings
  flag_eigmode = 1        ! Enable eigenmode search
  search_flag = 1
  rdim = 200
  idim = 200
  rfmin = 25.0e6
  rfmax = 35.0e6
  ifmin = -5000.0
  ifmax = 5000.0
  n_zeros = 20
/
EOF

# Run eigenmode search
/path/to/kilca_main .
```

### Performance Notes

- **Memory Usage**: Large frequency meshes (`rdim` × `idim`) can require significant memory
- **Calculation Time**: Higher spline orders (`N`) and finer meshes increase computation time
- **Convergence**: For eigenmode search, use appropriate frequency ranges around expected resonances

### Integration with KAMEL Workflow

The KiLCA Fortran main program integrates with the broader KAMEL workflow:

1. **Preprocessing**: Use KAMEL template scripts to prepare input files
2. **Main Calculation**: Run `kilca_main` with prepared inputs
3. **Post-processing**: Use Python/MATLAB interfaces for analysis and visualization

See the main KAMEL documentation for complete workflow examples.

## Git Workflow

- Build artifacts (`*.o`, `*.mod`, executables) are ignored
- Test source files are version controlled
- Clean directory structure for collaboration