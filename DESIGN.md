# KiLCA C/C++ to Fortran Translation - Architectural Design

## Executive Summary

This document outlines the comprehensive architectural plan for translating all C/C++ code in the KiLCA (Kinetic plAsma response ModEL - Cylindrical) component to Fortran. The translation will maintain exact computational equivalence while leveraging Fortran's strengths in scientific computing and ensuring seamless integration with the existing Fortran codebase.

## Current State Analysis

### C/C++ Code Inventory

**Total Files**: 110 C/C++ files (.cpp and .h files)

**Core Architecture Components**:

1. **Main Program Structure** (4 files)
   - `progs/main_linear.cpp/.h` - Primary driver program
   - `progs/main_eig_param.cpp` - Eigenmode parameter analysis

2. **Core System** (12 files)
   - `core/core.cpp/.h` - Top-level core_data class
   - `core/settings.cpp/.h` - Configuration management
   - `core/shared.cpp/.h` - Shared utilities
   - `core/output_sett.cpp/.h` - Output settings management
   - `core/constants.h` - Physical constants
   - `core/typedefs.h` - Type definitions
   - `core/slatec.h` - SLATEC interface declarations

3. **Background Plasma Calculations** (10 files)
   - `background/background.cpp/.h` - Main background class
   - `background/calc_back.cpp/.h` - Background calculations
   - `background/eval_back.cpp/.h` - Background evaluation
   - `background/back_sett.cpp/.h` - Background settings
   - `background/calc_back_slatec.cpp` - SLATEC-based calculations

4. **Mode Analysis System** (12 files)
   - `mode/mode.cpp/.h` - Main mode_data class
   - `mode/calc_mode.cpp` - Mode calculations
   - `mode/zone.cpp/.h` - Zone management
   - `mode/transforms.cpp/.h` - Coordinate transformations
   - `mode/stitching.h` - Mode stitching algorithms
   - `mode/wave_data.h` - Wave data structures

5. **FLRE (Finite Larmor Radius Effects) System** (18 files)
   - `flre/conductivity/` - 8 files for conductivity calculations
   - `flre/dispersion/` - 2 files for dispersion relations
   - `flre/maxwell_eqs/` - 6 files for Maxwell equations
   - `flre/quants/` - 6 files for physical quantities
   - `flre/flre_zone.cpp/.h` - FLRE zone management

6. **Solver Framework** (8 files)
   - `solver/solver.cpp/.h` - Main solver class
   - `solver/rhs_func.cpp/.h` - Right-hand side functions
   - `solver/eigtransform.cpp/.h` - Eigenmode transformations
   - `solver/solver_a.cpp` - Alternative solver
   - `solver/method.h` - Solver method definitions

7. **Mathematical Libraries** (30+ files)
   - `math/adapt_grid/` - Adaptive grid algorithms
   - `math/fourier/` - Fourier transforms
   - `math/hyper/` - Hypergeometric functions
   - `math/zersol-0.0.0/` - Zero-finding solver (C++ template library)

8. **Supporting Systems** (16 files)
   - `antenna/antenna.cpp/.h` - Antenna modeling
   - `interface/` - External interfaces
   - `io/inout.cpp/.h` - Input/output operations
   - `interp/interp.cpp/.h` - Interpolation routines
   - `spline/spline.cpp/.h` - Spline interpolation
   - `hom_medium/`, `imhd/` - Physical model implementations

### Key Architectural Patterns

1. **Object-Oriented Design**: Heavy use of C++ classes with inheritance
2. **Pointer-based Data Management**: Extensive use of pointers and dynamic allocation
3. **Template Programming**: Generic mathematical algorithms
4. **Mixed Language Integration**: C/C++ calling Fortran subroutines
5. **Complex Data Structures**: Multi-dimensional arrays and nested objects

## Translation Strategy

### Fundamental Principles

1. **No Shortcuts or Simplifications**: Every line of C/C++ code must have a direct Fortran equivalent
2. **Computational Equivalence**: All calculations must produce identical results
3. **Memory Layout Preservation**: Array storage and memory access patterns must be maintained
4. **Interface Compatibility**: External interfaces must remain unchanged
5. **Performance Parity**: No performance degradation allowed

### Translation Methodology

#### Phase 1: Data Structure Migration
- Convert C++ classes to Fortran derived types
- Migrate pointer-based structures to allocatable arrays
- Preserve memory layout and access patterns
- Implement object lifecycle management (constructor/destructor equivalents)

#### Phase 2: Algorithm Translation
- Direct line-by-line translation of computational routines
- Preserve all mathematical operations exactly
- Maintain numerical precision and error handling
- Convert template functions to specific typed procedures

#### Phase 3: Interface Harmonization
- Replace C++ method calls with Fortran procedure calls
- Implement C++ virtual function behavior with Fortran procedure pointers
- Maintain external API compatibility
- Preserve file I/O formats and protocols

#### Phase 4: Integration and Validation
- Link with existing Fortran modules seamlessly
- Validate against reference calculations
- Performance benchmarking
- Memory usage optimization

### Technical Approach

#### Object-Oriented to Procedural Mapping

**C++ Class Structure**:
```cpp
class core_data {
    char *path2project;
    settings *sd;
    background *bp;
    mode_data **mda;
public:
    void calc_and_set_mode_independent_core_data();
};
```

**Fortran Equivalent**:
```fortran
type :: core_data_t
    character(len=:), allocatable :: path2project
    type(settings_t), pointer :: sd
    type(background_t), pointer :: bp
    type(mode_data_t), dimension(:), allocatable :: mda
end type

subroutine calc_and_set_mode_independent_core_data(this)
    type(core_data_t), intent(inout) :: this
```

#### Memory Management Strategy

1. **Dynamic Allocation**: Convert `new`/`delete` to `allocate`/`deallocate`
2. **Pointer Management**: Use Fortran pointers and targets where necessary
3. **Array Bounds**: Implement explicit bounds checking
4. **Memory Leaks**: Ensure proper cleanup in all code paths

#### Template Code Handling

C++ templates will be expanded to specific Fortran procedures:
- Generic mathematical functions → specific real/complex procedures
- Container templates → specific array handling routines
- Type-generic algorithms → overloaded procedures

## Module Structure Design

### Core Modules Hierarchy

```
kilca_core_m
├── kilca_types_m           ! Type definitions
├── kilca_constants_m       ! Physical constants
├── kilca_settings_m        ! Configuration management
├── kilca_background_m      ! Background plasma
├── kilca_mode_m           ! Mode analysis
├── kilca_flre_m           ! FLRE calculations
├── kilca_solver_m         ! Equation solvers
├── kilca_math_m           ! Mathematical utilities
├── kilca_io_m             ! Input/output
└── kilca_interface_m      ! External interfaces
```

### Data Type Definitions

**Core Types**:
```fortran
type :: core_data_t
    character(len=:), allocatable :: path2project
    type(settings_t) :: settings_data
    type(background_t) :: background_data
    integer :: n_modes
    type(mode_data_t), dimension(:), allocatable :: modes
end type

type :: mode_data_t
    integer :: m_number, n_number
    integer :: n_zones
    type(zone_t), dimension(:), allocatable :: zones
    real(dp), dimension(:), allocatable :: r_grid
    complex(dp), dimension(:,:), allocatable :: eb_fields
end type
```

### Procedure Interface Design

**Consistent Naming Convention**:
- Module procedures: `module_action_object` (e.g., `core_create_data`)
- Type-bound procedures: `action_object` (e.g., `calculate_background`)
- Mathematical functions: descriptive names (e.g., `hypergeometric_1f1`)

**Error Handling Strategy**:
- Optional integer `status` parameter for all procedures
- Standardized error codes
- Comprehensive error messages
- Graceful failure modes

## Critical Translation Challenges

### 1. Complex Template Libraries

**Challenge**: The `zersol-0.0.0` library uses advanced C++ templates
**Solution**: 
- Expand templates to specific Fortran procedures
- Maintain algorithmic complexity and performance
- Preserve numerical stability properties

### 2. Dynamic Memory Management

**Challenge**: Extensive use of pointers and dynamic allocation
**Solution**:
- Systematic conversion to Fortran allocatable arrays
- Reference counting for shared data structures
- Explicit memory management procedures

### 3. Mixed-Language Interfaces

**Challenge**: Current C++ code calls Fortran subroutines
**Solution**:
- Direct procedure calls within Fortran
- Maintain external C interface for backward compatibility
- Preserve calling conventions and data passing

### 4. Performance-Critical Sections

**Challenge**: Mathematical kernels requiring optimal performance
**Solution**:
- Profile-guided optimization during translation
- Maintain algorithmic complexity
- Leverage Fortran's array operations
- Preserve vectorization opportunities

## Quality Assurance Framework

### Validation Strategy

1. **Unit Testing**: Every translated procedure against C++ reference
2. **Integration Testing**: Complete physics calculations comparison
3. **Regression Testing**: Against known reference solutions
4. **Performance Testing**: Computational efficiency validation
5. **Memory Testing**: No memory leaks or corruption

### Verification Metrics

- **Numerical Accuracy**: Machine precision equivalence
- **Performance Ratio**: ±5% of original C++ performance
- **Memory Usage**: Comparable or improved memory footprint
- **Code Coverage**: 100% translation of source lines
- **Interface Compatibility**: All external calls preserved

### Testing Framework

**Reference Comparison**:
```fortran
program validate_translation
    ! Load reference C++ results
    ! Execute Fortran implementation
    ! Compare with specified tolerance
    ! Report discrepancies
end program
```

## Implementation Dependencies

### Required Fortran Features

- **Fortran 2008**: Object-oriented programming support
- **Fortran 2018**: Advanced parameterized derived types
- **ISO_C_BINDING**: For external interface compatibility
- **IEEE_ARITHMETIC**: Numerical precision control

### External Library Integration

- **LAPACK**: Linear algebra operations (existing integration)
- **MPI**: Parallel computing support (existing integration)
- **HDF5**: Data I/O compatibility (existing integration)
- **GSL**: Mathematical functions (requires interface)

### Build System Modifications

- Update CMake configuration for pure Fortran build
- Maintain compatibility with existing Fortran modules
- Preserve external library dependencies
- Support both static and dynamic linking

## Risk Mitigation

### High-Risk Elements

1. **Complex Mathematical Algorithms**: Extensive validation required
2. **Memory Management**: Potential for subtle bugs
3. **Floating-Point Precision**: Numerical equivalence challenges
4. **Performance Regression**: Optimization difficulties

### Mitigation Strategies

1. **Incremental Translation**: Module-by-module approach
2. **Continuous Validation**: Test after each module translation
3. **Reference Preservation**: Keep C++ version for comparison
4. **Expert Review**: Mathematical physicist code review
5. **Automated Testing**: Comprehensive CI/CD pipeline

## Success Criteria

### Technical Metrics

- [ ] 100% of C++ code translated to Fortran
- [ ] Numerical equivalence to machine precision
- [ ] Performance within ±5% of original
- [ ] Memory usage comparable or improved
- [ ] All existing interfaces preserved
- [ ] Zero memory leaks or corruption
- [ ] Complete test coverage

### Qualitative Goals

- [ ] Code maintainability improved
- [ ] Documentation enhanced
- [ ] Debugging capabilities enhanced
- [ ] Development workflow simplified
- [ ] Scientific reproducibility maintained

## Timeline Implications

This translation represents a substantial undertaking requiring:
- **Estimated Effort**: 6-12 months for complete translation
- **Resource Requirements**: Experienced Fortran/C++ developer
- **Validation Period**: 2-3 months for comprehensive testing
- **Documentation**: Technical documentation and user guides

The "no shortcuts, no simplifications" requirement ensures this will be a comprehensive, methodical translation that preserves all functionality while modernizing the codebase architecture.