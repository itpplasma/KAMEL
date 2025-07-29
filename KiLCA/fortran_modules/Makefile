# Fortran modules Makefile
FC = gfortran
FFLAGS = -Wall -std=f2008 -g -O2 -Wno-compare-reals
LDFLAGS = -llapack -lblas
# Add bessel.f90 object file for Bessel function support
BESSEL_OBJ = ../math/bessel/bessel.o
BUILDDIR = build
OBJDIR = $(BUILDDIR)/obj
TESTDIR = $(BUILDDIR)/tests
SOURCES = kilca_types_m.f90 kilca_constants_m.f90 kilca_shared_m.f90 bessel_stubs.f90 kilca_complex_m.f90 kilca_settings_m.f90 kilca_core_m.f90 kilca_stability_m.f90 kilca_spline_m.f90 kilca_background_m.f90 kilca_mode_m.f90 kilca_eigensolver_m.f90 kilca_orthogonalization_m.f90 kilca_mode_solver_m.f90 kilca_interp_integ_m.f90 kilca_zone_m.f90
OBJECTS = $(SOURCES:%.f90=$(OBJDIR)/%.o)
MODULES = $(SOURCES:%.f90=$(OBJDIR)/%.mod)

# Test source files
COMPLEX_TESTS_SRC = $(wildcard tests/complex/test_*.f90)
GENERAL_TESTS_SRC = $(wildcard tests/test_*.f90)
COMPLEX_TESTS = $(COMPLEX_TESTS_SRC:tests/complex/%.f90=$(TESTDIR)/%)
GENERAL_TESTS = $(GENERAL_TESTS_SRC:tests/%.f90=$(TESTDIR)/%)

.PHONY: all clean test test-complex test-general directories

all: directories $(OBJECTS)

directories:
	@mkdir -p $(OBJDIR) $(TESTDIR)

# Module compilation rules
$(OBJDIR)/kilca_types_m.o: kilca_types_m.f90 | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_constants_m.o: kilca_constants_m.f90 $(OBJDIR)/kilca_types_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_shared_m.o: kilca_shared_m.f90 $(OBJDIR)/kilca_types_m.o $(OBJDIR)/kilca_constants_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/bessel_stubs.o: bessel_stubs.f90 $(OBJDIR)/kilca_types_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_complex_m.o: kilca_complex_m.f90 $(OBJDIR)/kilca_types_m.o $(OBJDIR)/kilca_shared_m.o $(OBJDIR)/bessel_stubs.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_settings_m.o: kilca_settings_m.f90 $(OBJDIR)/kilca_types_m.o $(OBJDIR)/kilca_constants_m.o $(OBJDIR)/kilca_shared_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_core_m.o: kilca_core_m.f90 $(OBJDIR)/kilca_types_m.o $(OBJDIR)/kilca_constants_m.o $(OBJDIR)/kilca_shared_m.o $(OBJDIR)/kilca_settings_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_stability_m.o: kilca_stability_m.f90 $(OBJDIR)/kilca_types_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_spline_m.o: kilca_spline_m.f90 $(OBJDIR)/kilca_types_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_background_m.o: kilca_background_m.f90 $(OBJDIR)/kilca_types_m.o $(OBJDIR)/kilca_settings_m.o $(OBJDIR)/kilca_spline_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_mode_m.o: kilca_mode_m.f90 $(OBJDIR)/kilca_types_m.o $(OBJDIR)/kilca_settings_m.o $(OBJDIR)/kilca_background_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_eigensolver_m.o: kilca_eigensolver_m.f90 $(OBJDIR)/kilca_types_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_orthogonalization_m.o: kilca_orthogonalization_m.f90 $(OBJDIR)/kilca_types_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_mode_solver_m.o: kilca_mode_solver_m.f90 $(OBJDIR)/kilca_types_m.o $(OBJDIR)/kilca_mode_m.o $(OBJDIR)/kilca_background_m.o $(OBJDIR)/kilca_eigensolver_m.o $(OBJDIR)/kilca_orthogonalization_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_interp_integ_m.o: kilca_interp_integ_m.f90 $(OBJDIR)/kilca_types_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/kilca_zone_m.o: kilca_zone_m.f90 $(OBJDIR)/kilca_constants_m.o $(OBJDIR)/kilca_settings_m.o $(OBJDIR)/kilca_background_m.o $(OBJDIR)/kilca_mode_m.o | directories
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

# Compile bessel.f90 if needed
$(BESSEL_OBJ): ../math/bessel/bessel.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Complex test compilation rules
$(TESTDIR)/test_%: tests/complex/test_%.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(OBJDIR) $(OBJECTS) $< -o $@ $(LDFLAGS)

# General test compilation rules  
$(TESTDIR)/test_%: tests/test_%.f90 $(OBJECTS) | directories
	$(FC) $(FFLAGS) -I$(OBJDIR) $(OBJECTS) $< -o $@ $(LDFLAGS)

# Test targets
test-complex: $(COMPLEX_TESTS)
	@echo "Running complex number tests..."
	@for test in $(COMPLEX_TESTS); do \
		echo "Running $$test..."; \
		$$test || exit 1; \
		echo ""; \
	done
	@echo "All complex tests passed!"

test-general: $(GENERAL_TESTS)
	@echo "Running general tests..."
	@for test in $(GENERAL_TESTS); do \
		echo "Running $$test..."; \
		$$test || exit 1; \
		echo ""; \
	done
	@echo "All general tests passed!"

test: test-general test-complex
	@echo "All tests completed successfully!"

# Clean targets
clean:
	rm -rf $(BUILDDIR)
	rm -f $(BESSEL_OBJ)

clean-tests:
	rm -rf $(TESTDIR)

clean-obj:
	rm -rf $(OBJDIR)

# Help target
help:
	@echo "Available targets:"
	@echo "  all           - Build all modules"
	@echo "  test          - Run all tests"
	@echo "  test-complex  - Run complex number tests only"
	@echo "  test-general  - Run general tests only"
	@echo "  clean         - Remove all build artifacts"
	@echo "  clean-tests   - Remove test executables"
	@echo "  clean-obj     - Remove object files"
	@echo "  help          - Show this help message"