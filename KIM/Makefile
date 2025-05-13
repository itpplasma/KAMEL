BUILD_DIR := build

.PHONY: all ninja test install clean

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -D CMAKE_EXPORT_COMPILE_COMMANDS=ON

ninja: build/build.ninja
	cmake --build build

test: ninja
	ctest --test-dir build

clean:
	rm -rf build
