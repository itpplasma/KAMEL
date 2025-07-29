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

## Git Workflow

- Build artifacts (`*.o`, `*.mod`, executables) are ignored
- Test source files are version controlled
- Clean directory structure for collaboration