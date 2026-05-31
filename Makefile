CONFIG ?= Release
INSTALL_KIM_SYMLINK ?= OFF

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
	cmake -S . -B build -G $(CMAKE_GENERATOR) -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DINSTALL_KIM_SYMLINK=$(INSTALL_KIM_SYMLINK)

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

test: build golden
	ctest --test-dir build --stop-on-failure --output-on-failure --no-label-summary

golden: build
	$(MAKE) -C test/ql-balance/golden_record

clean:
	rm -rf build
	$(MAKE) -C test/ql-balance/golden_record clean

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
