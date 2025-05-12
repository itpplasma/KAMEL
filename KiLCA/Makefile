CONFIG ?= Release

.PHONY: all ninja test install clean
all: ninja

build/build.ninja:
	cmake -G Ninja -S . -B build -DCMAKE_BUILD_TYPE=$(CONFIG)

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

clean:
	rm -rf build
