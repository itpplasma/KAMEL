CONFIG ?= Release

.PHONY: all ninja test clean KIM KiLCA QL-Balance

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

KIM: build/build.ninja
	cmake --build build --config $(CONFIG) --target KIM.x

KiLCA: build/build.ninja
	cmake --build build --config $(CONFIG) --target KiLCA

QL-Balance: build/build.ninja
	cmake --build build --config $(CONFIG) --target QL-Balance

test: ninja
	ctest --test-dir build --stop-on-failure --output-on-failure

clean:
	rm -rf build
