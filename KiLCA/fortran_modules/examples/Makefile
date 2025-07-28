# Makefile for KiLCA Settings Examples
# Usage: make [target]

# Compiler settings
FC = mpif90
FFLAGS = -fopenmp -O2 -g -fbacktrace -fcheck=all -Wall

# Directories
BUILD_DIR = ../build
MODULE_DIR = $(BUILD_DIR)/modules
LIB_DIR = $(BUILD_DIR)

# Library
LIBRARY = $(LIB_DIR)/libkilca_fortran_modules.a

# Include directories
INCLUDES = -I$(MODULE_DIR)

# Link libraries
LDFLAGS = -L$(LIB_DIR) -lkilca_fortran_modules

# Source files
BASIC_SRC = example_basic_usage.f90
ADVANCED_SRC = example_advanced_usage.f90

# Executables
BASIC_EXE = example_basic_usage
ADVANCED_EXE = example_advanced_usage

.PHONY: all basic advanced clean help

# Default target
all: basic advanced

# Basic usage example
basic: $(BASIC_EXE)

$(BASIC_EXE): $(BASIC_SRC) $(LIBRARY)
	@echo "Building basic usage example..."
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)
	@echo "Basic example built successfully"

# Advanced usage example  
advanced: $(ADVANCED_EXE)

$(ADVANCED_EXE): $(ADVANCED_SRC) $(LIBRARY)
	@echo "Building advanced usage example..."
	$(FC) $(FFLAGS) $(INCLUDES) -o $@ $< $(LDFLAGS)
	@echo "Advanced example built successfully"

# Ensure library is built
$(LIBRARY):
	@echo "Building KiLCA Fortran modules library..."
	@cd .. && mkdir -p build && cd build && cmake .. && make kilca_fortran_modules
	@echo "Library built successfully"

# Run examples
run-basic: $(BASIC_EXE)
	@echo "Running basic usage example..."
	@echo "================================"
	./$(BASIC_EXE)

run-advanced: $(ADVANCED_EXE)
	@echo "Running advanced usage example..."
	@echo "=================================="
	./$(ADVANCED_EXE)

run-all: run-basic run-advanced

# Clean build artifacts
clean:
	@echo "Cleaning examples..."
	rm -f $(BASIC_EXE) $(ADVANCED_EXE)
	rm -f *.mod *.o
	rm -f advanced_settings_output.txt
	rm -f advanced_example_errors.log
	rm -f example_eigenmode_output.dat
	@echo "Examples cleaned"

# Clean everything including library
clean-all: clean
	@echo "Cleaning library..."
	rm -rf $(BUILD_DIR)
	@echo "Everything cleaned"

# Help target
help:
	@echo "KiLCA Settings Examples Makefile"
	@echo "================================"
	@echo ""
	@echo "Available targets:"
	@echo "  all         - Build both examples (default)"
	@echo "  basic       - Build basic usage example"
	@echo "  advanced    - Build advanced usage example"
	@echo "  run-basic   - Build and run basic example"
	@echo "  run-advanced- Build and run advanced example"
	@echo "  run-all     - Build and run both examples"
	@echo "  clean       - Clean example executables"
	@echo "  clean-all   - Clean examples and library"
	@echo "  help        - Show this help"
	@echo ""
	@echo "Examples:"
	@echo "  make basic          # Build basic example"
	@echo "  make run-advanced   # Build and run advanced example"
	@echo "  make clean-all      # Clean everything"