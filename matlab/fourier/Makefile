.PHONY: all ninja test install clean
all: ninja

build/build.ninja:
	cmake -GNinja -S . -B build

ninja: build/build.ninja
	cmake --build build

clean:
	rm -rf build