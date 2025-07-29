# Makefile for KiLCA Settings Module Examples
# This makefile builds the example programs that demonstrate usage of the kilca_settings_m module

# Compiler settings
FC = mpif90
FFLAGS = -fopenmp -O2 -g -fbacktrace -fcheck=all -Wall
DEBUG_FLAGS = -g -O0 -fcheck=all
RELEASE_FLAGS = -O3

# Directories  
PARENT_DIR = ..
BUILD_DIR = build
MODULE_DIR = $(PARENT_DIR)/build/modules
LIB_DIR = $(PARENT_DIR)/build

# Library
LIBNAME = kilca_fortran_modules
LIBRARY = $(LIB_DIR)/lib$(LIBNAME).a

# Source files
BASIC_SRC = example_basic_usage.f90
ADVANCED_SRC = example_advanced_usage.f90

# Executables
BASIC_EXE = $(BUILD_DIR)/example_basic_usage
ADVANCED_EXE = $(BUILD_DIR)/example_advanced_usage

# Default target
.PHONY: all clean examples help basic advanced debug

all: examples

examples: $(BASIC_EXE) $(ADVANCED_EXE)

basic: $(BASIC_EXE)

advanced: $(ADVANCED_EXE)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Check if parent library exists
$(LIBRARY):
	@echo "Building parent library first..."
	cd $(PARENT_DIR) && make

# Build basic example
$(BASIC_EXE): $(BASIC_SRC) $(LIBRARY) | $(BUILD_DIR)
	@echo "Building basic usage example..."
	$(FC) $(FFLAGS) -I$(MODULE_DIR) -o $@ $< -L$(LIB_DIR) -l$(LIBNAME)

# Build advanced example  
$(ADVANCED_EXE): $(ADVANCED_SRC) $(LIBRARY) | $(BUILD_DIR)
	@echo "Building advanced usage example..."
	$(FC) $(FFLAGS) -I$(MODULE_DIR) -o $@ $< -L$(LIB_DIR) -l$(LIBNAME)

# Debug builds
debug: FFLAGS = $(DEBUG_FLAGS) -fopenmp -Wall
debug: examples

# Release builds
release: FFLAGS = $(RELEASE_FLAGS) -fopenmp
release: examples

# Run examples
.PHONY: run-basic run-advanced run-all

run-basic: $(BASIC_EXE)
	@echo "Running basic usage example..."
	@echo "================================"
	cd $(BUILD_DIR) && ./example_basic_usage

run-advanced: $(ADVANCED_EXE)
	@echo "Running advanced usage example..."
	@echo "================================="
	cd $(BUILD_DIR) && ./example_advanced_usage

run-all: run-basic run-advanced

# Test that examples compile and link correctly
.PHONY: test

test: examples
	@echo "Testing example compilation and linking..."
	@echo "✓ Basic example compiled successfully"
	@echo "✓ Advanced example compiled successfully" 
	@echo "Examples are ready to run with 'make run-all'"

# Install examples (copy to a system location)
.PHONY: install

install: examples
	@echo "Installing examples to /usr/local/share/kilca/examples..."
	@mkdir -p /usr/local/share/kilca/examples
	@cp $(BASIC_EXE) $(ADVANCED_EXE) /usr/local/share/kilca/examples/
	@cp $(BASIC_SRC) $(ADVANCED_SRC) /usr/local/share/kilca/examples/
	@cp Makefile /usr/local/share/kilca/examples/
	@echo "Examples installed successfully"

# Clean build artifacts
clean:
	@echo "Cleaning example build artifacts..."
	rm -rf $(BUILD_DIR)
	rm -f *.mod *.o
	rm -f error_demo.log

# Help target
help:
	@echo "KiLCA Settings Module Examples - Available targets:"
	@echo ""
	@echo "Building:"
	@echo "  all          - Build all examples (default)"
	@echo "  examples     - Build all examples"
	@echo "  basic        - Build basic usage example only"
	@echo "  advanced     - Build advanced usage example only"
	@echo "  debug        - Build with debug flags"
	@echo "  release      - Build with optimized flags"
	@echo ""
	@echo "Running:"
	@echo "  run-basic    - Run basic usage example"
	@echo "  run-advanced - Run advanced usage example"  
	@echo "  run-all      - Run both examples"
	@echo ""
	@echo "Maintenance:"
	@echo "  test         - Test that examples compile"
	@echo "  clean        - Remove build artifacts"
	@echo "  install      - Install examples to system location"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - Parent kilca_fortran_modules library must be built first"
	@echo "  - Use 'cd .. && make' to build the library if needed"
	@echo ""
	@echo "Example usage:"
	@echo "  make              # Build all examples"
	@echo "  make run-all      # Build and run both examples"
	@echo "  make debug        # Build with debugging enabled"
	@echo "  make clean        # Clean up build files"

# Show variables (for debugging makefile)
.PHONY: show-vars

show-vars:
	@echo "Makefile variables:"
	@echo "  FC = $(FC)"
	@echo "  FFLAGS = $(FFLAGS)"
	@echo "  PARENT_DIR = $(PARENT_DIR)"
	@echo "  BUILD_DIR = $(BUILD_DIR)"
	@echo "  MODULE_DIR = $(MODULE_DIR)"
	@echo "  LIB_DIR = $(LIB_DIR)"
	@echo "  LIBRARY = $(LIBRARY)"
	@echo "  BASIC_EXE = $(BASIC_EXE)"
	@echo "  ADVANCED_EXE = $(ADVANCED_EXE)"