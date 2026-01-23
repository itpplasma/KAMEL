CONFIG ?= Release
INSTALL_KIM_SYMLINK ?= OFF

.PHONY: all ninja test clean KIM KiLCA QL-Balance PreProc install install-kim ctest golden pytest

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DINSTALL_KIM_SYMLINK=$(INSTALL_KIM_SYMLINK)

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

test: ninja golden
	pytest test/
	ctest --test-dir build --stop-on-failure --output-on-failure --no-label-summary

golden: ninja
	$(MAKE) -C test/ql-balance/golden_record

clean:
	rm -rf build
	$(MAKE) -C test/ql-balance/golden_record clean

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
