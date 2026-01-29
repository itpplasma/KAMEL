.PHONY: all ninja test install clean
all: ninja

build/build.ninja:
	cmake -GNinja -S . -B build

ninja: build/build.ninja
	cmake --build build

clean:
	rm -rf build
	rm src/libneo*.f90
	rm src/*libneo.f90
