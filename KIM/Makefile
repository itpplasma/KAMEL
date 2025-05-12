BUILD_DIR := build

.PHONY: all ninja test install clean

all: ninja

build/build.ninja:
	cmake -G Ninja -S . -B build

ninja: build/build.ninja
	cmake --build build

test: ninja
	ctest --test-dir build

clean:
	rm -rf build
