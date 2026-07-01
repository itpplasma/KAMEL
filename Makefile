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

# Detect Ninja; fall back to Unix Makefiles if not found
NINJA := $(shell command -v ninja 2>/dev/null)
ifdef NINJA
    CMAKE_GENERATOR := Ninja
    BUILD_SENTINEL := build/build.ninja
else
    $(warning *** Ninja not found — falling back to Unix Makefiles. Install ninja for faster builds.)
    CMAKE_GENERATOR := "Unix Makefiles"
    BUILD_SENTINEL := build/Makefile
endif

.PHONY: all build test clean KIM KiLCA QL-Balance PreProc install install-kim ctest golden pytest

all: build

$(BUILD_SENTINEL):
	cmake -S . -B build -G $(CMAKE_GENERATOR) -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DINSTALL_KIM_SYMLINK=$(INSTALL_KIM_SYMLINK) $(_LIBNEO_DEFS)

build: $(BUILD_SENTINEL)
	cmake --build build --config $(CONFIG)

KIM: $(BUILD_SENTINEL)
	cmake --build build --config $(CONFIG) --target KIM.x

KiLCA: $(BUILD_SENTINEL)
	cmake --build build --config $(CONFIG) --target KiLCA

QL-Balance: $(BUILD_SENTINEL)
	cmake --build build --config $(CONFIG) --target ql-balance.x

PreProc:
	$(MAKE) -C PreProc/fourier

test: build
	ctest --test-dir build --stop-on-failure --output-on-failure --no-label-summary

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

install: build
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
