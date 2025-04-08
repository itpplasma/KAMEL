CONFIG ?= Release

.PHONY: all ninja test install clean
all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) 

ninja: build/build.ninja
	cmake --build build

test: build
	ctest --test-dir build/tests --output-on-failure

clean:
	rm -rf build
