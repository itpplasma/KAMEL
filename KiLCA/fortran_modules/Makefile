CONFIG ?= Release

.PHONY: all ninja test test-complex test-general clean help

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

test: ninja
	ctest --test-dir build --stop-on-failure --output-on-failure

test-complex: ninja
	ctest --test-dir build -R "_complex$$" --stop-on-failure --output-on-failure

test-general: ninja
	ctest --test-dir build -E "_complex$$" --stop-on-failure --output-on-failure

clean:
	rm -rf build

help:
	@echo "Available targets:"
	@echo "  all           - Build all modules (default)"
	@echo "  ninja         - Build using Ninja"
	@echo "  test          - Run all tests"
	@echo "  test-complex  - Run complex number tests only"
	@echo "  test-general  - Run general tests only"
	@echo "  clean         - Remove all build artifacts"
	@echo "  help          - Show this help message"
	@echo ""
	@echo "Build configurations:"
	@echo "  CONFIG=Debug   - Build with debug symbols"
	@echo "  CONFIG=Release - Build with optimizations (default)"
	@echo ""
	@echo "Example usage:"
	@echo "  make                   # Build in Release mode"
	@echo "  make CONFIG=Debug      # Build in Debug mode"
	@echo "  make test              # Build and run all tests"
	@echo "  make clean             # Clean build directory"