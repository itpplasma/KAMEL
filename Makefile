CONFIG ?= Release
PYTHON ?= python3

# Honor LIBNEO_REF/LIBNEO_PATH only when passed on the make command line; an
# ambient value from the shell is ignored so it cannot change the libneo fetch.
unexport LIBNEO_REF LIBNEO_PATH

_LIBNEO_DEFS :=
ifeq ($(origin LIBNEO_REF),command line)
  _LIBNEO_DEFS += -DLIBNEO_REF=$(LIBNEO_REF)
endif
ifeq ($(origin LIBNEO_PATH),command line)
  _LIBNEO_DEFS += -DLIBNEO_PATH=$(LIBNEO_PATH)
endif

.PHONY: all ninja test clean KIM KiLCA QL-Balance PreProc install ctest golden pytest

all: ninja

build/build.ninja:
	cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_EXPORT_COMPILE_COMMANDS=ON $(_LIBNEO_DEFS)

ninja: build/build.ninja
	cmake --build build --config $(CONFIG)

KIM: build/build.ninja
	cmake --build build --config $(CONFIG) --target KIM.x

KiLCA: build/build.ninja
	cmake --build build --config $(CONFIG) --target KiLCA

QL-Balance: build/build.ninja
	cmake --build build --config $(CONFIG) --target ql-balance.x

PreProc:
	$(MAKE) -C PreProc/fourier

test: ninja
	ctest --test-dir build --stop-on-failure --output-on-failure --no-label-summary

pytest:
	$(PYTHON) -m pytest test/golden/bin test/test_periodic_workflow_config.py -q

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
