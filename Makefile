CONFIG ?= Release
PYTHON ?= python3

# Extra CMake arguments for the first configure, e.g. a dependency override:
# make all CMAKE_ARGS=-DITP_LIBNEO_REF=<ref>
CMAKE_ARGS ?=

.PHONY: all ninja test clean KIM KiLCA QL-Balance PreProc install ctest golden pytest

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON $(CMAKE_ARGS)

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

KIM: build/build.ninja
	cmake --build build --config $(CONFIG) --target KIM.x

KiLCA: build/build.ninja
	cmake --build build --config $(CONFIG) --target KiLCA

QL-Balance: build/build.ninja
	cmake --build build --config $(CONFIG) --target ql-balance.x

PreProc: build/build.ninja
	cmake --build build --config $(CONFIG) --target fouriermodes.x

test: ninja
	ctest --test-dir build --stop-on-failure --output-on-failure --no-label-summary

pytest:
	$(PYTHON) -m pytest test/golden/bin test/test_periodic_workflow_config.py test/test_kimpy_collision_frequency.py -q

# Golden-record regression now lives in test/golden/ and runs in the dedicated
# GitHub Actions job, NOT in `make test`/ctest. This target is for manual local
# runs only; it does the A/B double build and needs the golden-baseline tag.
golden:
	test/golden/run_golden.sh

clean:
	rm -rf build
	rm -rf test/golden/build_ref test/golden/build_cur \
	       test/golden/*/build_ref test/golden/*/build_cur \
	       test/golden/*/out test/golden/*/logs test/golden/*/claims \
	       test/golden/*/caselist_builds.txt
	-git worktree prune

install: ninja
	cmake --install build
