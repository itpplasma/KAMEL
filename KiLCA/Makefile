CONFIG ?= Release

.PHONY: all ninja test install clean
all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -D CMAKE_EXPORT_COMPILE_COMMANDS=ON

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

clean:
	rm -rf build
