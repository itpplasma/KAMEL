CONFIG ?= Release

.PHONY: all ninja test clean

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

test: build
	ctest --test-dir build/tests --output-on-failure

clean:
	rm -rf build
