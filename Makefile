CONFIG ?= Release
INSTALL_KIM_SYMLINK ?= OFF

# Honor LIBNEO_REF/LIBNEO_PATH only when passed on the make command line; an
# ambient value from the shell is ignored so it cannot change the libneo fetch.
unexport LIBNEO_REF LIBNEO_PATH

_LIBNEO_DEFS :=
ifeq ($(origin LIBNEO_REF),command line)
  _LIBNEO_DEFS += -DLIBNEO_REF=$(LIBNEO_REF)
endif
ifeq ($(origin LIBNEO_PATH),command line)
  _LIBNEO_DEFS += -DLIBNEO_PATH=$(LIBNEO_PATH)
endif

.PHONY: all ninja test clean KIM KiLCA QL-Balance PreProc install install-kim ctest golden pytest \
	verify-mathematica-inventory verify-mathematica-inventory-negative \
	verify-mathematica-susceptibilities verify-mathematica-susceptibilities-negative \
	regenerate-fp-susceptibility-oracle verify-mathematica-charge-kernels \
	verify-mathematica-charge-kernels-negative regenerate-rho-phi-kernel-oracle \
	verify-mathematica-field-operators verify-mathematica-field-operators-negative \
	regenerate-field-operator-oracle

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DINSTALL_KIM_SYMLINK=$(INSTALL_KIM_SYMLINK) $(_LIBNEO_DEFS)

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

KIM: build/build.ninja
	cmake --build build --config $(CONFIG) --target KIM.x

KiLCA: build/build.ninja
	cmake --build build --config $(CONFIG) --target KiLCA

QL-Balance: build/build.ninja
	cmake --build build --config $(CONFIG) --target ql-balance.x

PreProc:
	$(MAKE) -C PreProc/fourier

test: ninja
	ctest --test-dir build --stop-on-failure --output-on-failure --no-label-summary

WOLFRAM_KERNEL ?= $(firstword $(shell command -v WolframKernel 2>/dev/null) \
	$(shell command -v math 2>/dev/null) $(wildcard \
	/Applications/Wolfram.app/Contents/MacOS/WolframKernel \
	/Applications/Mathematica.app/Contents/MacOS/WolframKernel))

verify-mathematica-inventory:
	@test -n "$(WOLFRAM_KERNEL)" || { echo "WolframKernel not found; set WOLFRAM_KERNEL"; exit 1; }
	$(WOLFRAM_KERNEL) -noinit -noprompt -script verification/mathematica/01_conventions_and_inventory.wl

verify-mathematica-inventory-negative:
	@for fixture in missing-thesis missing-code duplicate-convention unclassified; do \
		if KAMEL_VERIFY_NEGATIVE_FIXTURE=$$fixture $(MAKE) --no-print-directory verify-mathematica-inventory >/dev/null 2>&1; then \
			echo "negative fixture unexpectedly passed: $$fixture"; exit 1; \
		else \
			echo "negative fixture rejected: $$fixture"; \
		fi; \
	done

verify-mathematica-susceptibilities:
	@test -n "$(WOLFRAM_KERNEL)" || { echo "WolframKernel not found; set WOLFRAM_KERNEL"; exit 1; }
	$(WOLFRAM_KERNEL) -noinit -noprompt -script verification/mathematica/02_kinetic_and_susceptibilities.wl

verify-mathematica-susceptibilities-negative:
	@for fixture in wrong-branch thermal-speed flow-sign recurrence-index; do \
		if KAMEL_VERIFY_NEGATIVE_FIXTURE=$$fixture $(MAKE) --no-print-directory verify-mathematica-susceptibilities >/dev/null 2>&1; then \
			echo "negative fixture unexpectedly passed: $$fixture"; exit 1; \
		else \
			echo "negative fixture rejected: $$fixture"; \
		fi; \
	done

regenerate-fp-susceptibility-oracle:
	@test -n "$(WOLFRAM_KERNEL)" || { echo "WolframKernel not found; set WOLFRAM_KERNEL"; exit 1; }
	KAMEL_REGENERATE_ORACLE=1 $(WOLFRAM_KERNEL) -noinit -noprompt -script verification/mathematica/02_kinetic_and_susceptibilities.wl

verify-mathematica-charge-kernels:
	@test -n "$(WOLFRAM_KERNEL)" || { echo "WolframKernel not found; set WOLFRAM_KERNEL"; exit 1; }
	@test -f verification/mathematica/04_charge_potential_kernels.wl || { echo "charge-kernel verification script not found"; exit 1; }
	$(WOLFRAM_KERNEL) -noinit -noprompt -script verification/mathematica/04_charge_potential_kernels.wl

verify-mathematica-charge-kernels-negative:
	@for fixture in four-pi fourier-measure phase bessel-index bplus-sign hidden-ks2; do \
		if KAMEL_VERIFY_NEGATIVE_FIXTURE=$$fixture $(MAKE) --no-print-directory verify-mathematica-charge-kernels >/dev/null 2>&1; then \
			echo "negative fixture unexpectedly passed: $$fixture"; exit 1; \
		else \
			echo "negative fixture rejected: $$fixture"; \
		fi; \
	done

regenerate-rho-phi-kernel-oracle:
	@test -n "$(WOLFRAM_KERNEL)" || { echo "WolframKernel not found; set WOLFRAM_KERNEL"; exit 1; }
	KAMEL_REGENERATE_ORACLE=1 $(WOLFRAM_KERNEL) -noinit -noprompt -script verification/mathematica/04_charge_potential_kernels.wl

verify-mathematica-field-operators:
	@test -n "$(WOLFRAM_KERNEL)" || { echo "WolframKernel not found; set WOLFRAM_KERNEL"; exit 1; }
	@test -f verification/mathematica/05_magnetic_current_basis_poisson.wl || { echo "field-operator verification script not found"; exit 1; }
	$(WOLFRAM_KERNEL) -noinit -noprompt -script verification/mathematica/05_magnetic_current_basis_poisson.wl

verify-mathematica-field-operators-negative:
	@for fixture in source-sign boundary-term transform-factor transposed-kernel unjustified-gauge hidden-ks2; do \
		if KAMEL_VERIFY_NEGATIVE_FIXTURE=$$fixture $(MAKE) --no-print-directory verify-mathematica-field-operators >/dev/null 2>&1; then \
			echo "negative fixture unexpectedly passed: $$fixture"; exit 1; \
		else \
			echo "negative fixture rejected: $$fixture"; \
		fi; \
	done

regenerate-field-operator-oracle:
	@test -n "$(WOLFRAM_KERNEL)" || { echo "WolframKernel not found; set WOLFRAM_KERNEL"; exit 1; }
	KAMEL_REGENERATE_ORACLE=1 $(WOLFRAM_KERNEL) -noinit -noprompt -script verification/mathematica/05_magnetic_current_basis_poisson.wl

# Golden-record regression now lives in test/golden/ and runs in the dedicated
# GitHub Actions job, NOT in `make test`/ctest. This target is for manual local
# runs only; it does the A/B double build and needs the golden-baseline tag.
golden:
	test/golden/run_golden.sh

clean:
	rm -rf build
	rm -rf test/golden/build_ref test/golden/build_cur \
	       test/golden/*/build_ref test/golden/*/build_cur \
	       test/golden/*/out test/golden/*/logs test/golden/*/claims \
	       test/golden/*/caselist_builds.txt
	-git worktree prune

install: ninja
	cmake --install build

install-kim: KIM
	@echo "Creating symbolic link for KIM..."
	@if [ ! -f build/install/bin/KIM.x ]; then \
		echo "Error: KIM.x not found. Please build KIM first with 'make KIM'"; \
		exit 1; \
	fi
	@echo ""
	@echo "To complete the installation, run the following command:"
	@echo ""
	@echo "  sudo ln -sf $(PWD)/build/install/bin/KIM.x /usr/local/bin/kim"
	@echo ""
	@echo "After running the above command, you'll be able to run 'kim' from anywhere."
	@echo ""
	@echo "Alternatively, add this alias to your shell configuration:"
	@echo "  alias kim='$(PWD)/build/install/bin/KIM.x'"
