# Layered Dependency Discovery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the KAMEL build system work on macOS without Homebrew by using standard CMake discovery first and Homebrew as a fallback.

**Architecture:** Restructure dependency detection in CMakeLists.txt as layered cascades: user-provided variables first, then standard PATH/CMAKE_PREFIX_PATH discovery, then Homebrew as fallback. Remove duplicate OpenMP detection from QL-Balance/CMakeLists.txt.

**Tech Stack:** CMake 3.24+, Fortran/C/C++ compilers, OpenMP, HDF5

---

### Task 1: Refactor Compiler Detection Cascade

**Files:**
- Modify: `CMakeLists.txt:4-22`

- [ ] **Step 1: Replace the compiler detection block**

Replace lines 4-30 of `CMakeLists.txt` with the layered cascade. The key change: search `$PATH` first (no path restriction), then fall back to Homebrew LLVM paths. Append to `CMAKE_PREFIX_PATH` instead of overwriting.

```cmake
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    if(NOT DEFINED CMAKE_C_COMPILER)
        # 1. Search PATH first (works with Nix, manual installs, etc.)
        find_program(LLVM_CLANG clang)
        if(NOT LLVM_CLANG)
            # 2. Homebrew LLVM as fallback
            find_program(LLVM_CLANG clang
                PATHS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin
                NO_DEFAULT_PATH)
        endif()
        if(LLVM_CLANG)
            set(CMAKE_C_COMPILER "${LLVM_CLANG}")
        endif()
    endif()
    if(NOT DEFINED CMAKE_CXX_COMPILER)
        find_program(LLVM_CLANGXX clang++)
        if(NOT LLVM_CLANGXX)
            find_program(LLVM_CLANGXX clang++
                PATHS /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin
                NO_DEFAULT_PATH)
        endif()
        if(LLVM_CLANGXX)
            set(CMAKE_CXX_COMPILER "${LLVM_CLANGXX}")
        endif()
    endif()
    # Append Homebrew prefix paths only if brew is available (don't overwrite)
    find_program(_BREW_EXE brew HINTS /opt/homebrew/bin /usr/local/bin)
    if(_BREW_EXE)
        list(APPEND CMAKE_PREFIX_PATH /opt/homebrew /usr/local)
    endif()
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    if(NOT DEFINED CMAKE_C_COMPILER)
        set(CMAKE_C_COMPILER gcc)
    endif()
    if(NOT DEFINED CMAKE_CXX_COMPILER)
        set(CMAKE_CXX_COMPILER g++)
    endif()
endif()
```

- [ ] **Step 2: Verify the change parses correctly**

Run: `cd /Users/markl_m/code/KAMEL && cmake -S . -B build-test --trace-expand 2>&1 | head -50`

Expected: CMake configure starts without syntax errors in the compiler detection block. The exact configure may fail due to missing deps — that's fine at this stage, we just want to confirm the new logic is syntactically valid and runs.

- [ ] **Step 3: Commit**

```bash
git add CMakeLists.txt
git commit -m "build: refactor compiler detection to search PATH before Homebrew"
```

---

### Task 2: Refactor HDF5 Discovery Cascade

**Files:**
- Modify: `CMakeLists.txt:32-48`

- [ ] **Step 1: Replace the HDF5 discovery block**

Replace lines 32-48 of `CMakeLists.txt` with a cascade that tries standard CMake discovery before Homebrew. Note: `FetchHDF5.cmake` (included via `Dependencies.cmake`) already calls `find_package(HDF5)` and builds from source if not found — so this block only needs to help CMake *find* an existing HDF5 installation by setting `HDF5_DIR`.

```cmake
if(APPLE AND NOT DEFINED HDF5_DIR)
    # 1. Check if HDF5 is findable via standard paths / CMAKE_PREFIX_PATH
    find_package(HDF5 COMPONENTS C Fortran HL QUIET CONFIG)
    if(HDF5_FOUND)
        message(STATUS "HDF5 found via standard paths: ${HDF5_VERSION}")
    else()
        # 2. Homebrew fallback
        find_program(BREW_EXECUTABLE brew HINTS /opt/homebrew/bin /usr/local/bin)
        if(BREW_EXECUTABLE)
            execute_process(
                COMMAND "${BREW_EXECUTABLE}" --prefix hdf5
                OUTPUT_VARIABLE _hdf5_prefix
                OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE _hdf5_res
            )
            if(_hdf5_res EQUAL 0 AND EXISTS "${_hdf5_prefix}/lib/cmake/hdf5")
                set(HDF5_DIR "${_hdf5_prefix}/lib/cmake/hdf5" CACHE PATH "Path to HDF5Config.cmake")
                message(STATUS "HDF5_DIR set from Homebrew: ${HDF5_DIR}")
            endif()
        endif()
        # 3. If still not found, FetchHDF5.cmake will build from source
    endif()
endif()
```

- [ ] **Step 2: Verify the change parses correctly**

Run: `cd /Users/markl_m/code/KAMEL && cmake -S . -B build-test 2>&1 | grep -i "hdf5"`

Expected: HDF5 discovery messages appear — either "found via standard paths", "set from Homebrew", or the fetch-from-source message from `FetchHDF5.cmake`.

- [ ] **Step 3: Commit**

```bash
git add CMakeLists.txt
git commit -m "build: try standard HDF5 discovery before Homebrew fallback"
```

---

### Task 3: Refactor OpenMP/libomp Detection

**Files:**
- Modify: `CMakeLists.txt:57-81`

- [ ] **Step 1: Replace the OpenMP detection block**

Replace lines 57-81 of `CMakeLists.txt` (the `if(APPLE)` block for OpenMP) with a cascade that trusts `find_package(OpenMP)` and only falls back to Homebrew libomp if the OpenMP libraries weren't fully resolved.

```cmake
if(APPLE)
    find_package(OpenMP REQUIRED C Fortran)

    # Check if find_package fully resolved OpenMP (compiler knows where it is)
    if(NOT OpenMP_C_LIBRARIES OR NOT OpenMP_Fortran_LIBRARIES)
        # Homebrew libomp fallback — only if find_package didn't populate libraries
        find_program(BREW_EXECUTABLE brew HINTS /opt/homebrew/bin /usr/local/bin)
        if(BREW_EXECUTABLE)
            execute_process(
                COMMAND "${BREW_EXECUTABLE}" --prefix libomp
                OUTPUT_VARIABLE LIBOMP_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            if(EXISTS "${LIBOMP_PREFIX}")
                message(STATUS "Using Homebrew libomp at ${LIBOMP_PREFIX}")
                include_directories(SYSTEM "${LIBOMP_PREFIX}/include")
                link_directories("${LIBOMP_PREFIX}/lib")
                set(LIBOMP_LINK_FLAGS "-L${LIBOMP_PREFIX}/lib -lomp")
            endif()
        endif()
    else()
        message(STATUS "OpenMP fully resolved by find_package (no Homebrew fallback needed)")
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -O2 -g -fbacktrace -fcheck=all -Wall -Wextra -Wconversion -Wno-external-argument-mismatch")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    if(LIBOMP_LINK_FLAGS)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBOMP_LINK_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LIBOMP_LINK_FLAGS}")
    endif()
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp -O2 -mtune=generic -g -fbacktrace -fcheck=all -Wall -Wextra -Wconversion -Wno-external-argument-mismatch")
endif()
```

- [ ] **Step 2: Verify the change parses correctly**

Run: `cd /Users/markl_m/code/KAMEL && cmake -S . -B build-test 2>&1 | grep -iE "openmp|libomp"`

Expected: Either "OpenMP fully resolved by find_package" or "Using Homebrew libomp at ..." messages appear.

- [ ] **Step 3: Commit**

```bash
git add CMakeLists.txt
git commit -m "build: trust find_package(OpenMP) before Homebrew libomp fallback"
```

---

### Task 4: Remove Duplicate OpenMP Detection from QL-Balance

**Files:**
- Modify: `QL-Balance/CMakeLists.txt:6-38`

- [ ] **Step 1: Simplify the Darwin block**

The parent `CMakeLists.txt` already configures OpenMP flags, compiler flags, and libomp paths. These are inherited via `add_subdirectory`. Replace lines 6-38 of `QL-Balance/CMakeLists.txt` with a simplified block that only sets QL-Balance-specific flags (the parent doesn't set `-cpp`, `-fPIC`, `-DDOUBLE_APPEND_FORTRAN`, etc.).

```cmake
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
    message("Building on macOS")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DDOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -g -Og -Wall -Wno-unused-function -fbacktrace")
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
    message("Building on Linux; using default Fortran compiler: ${CMAKE_Fortran_COMPILER}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -D DOUBLE_APPEND_FORTRAN")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
else()
    message("Building on another operating system might not be implemented.")
endif()
```

This removes:
- `find_package(OpenMP REQUIRED)` — already done in parent
- `find_program(BREW_EXECUTABLE brew ...)` — no longer needed
- `brew --prefix libomp` detection — handled by parent
- `-L/opt/homebrew/lib` hardcoded path — removed entirely
- `${OpenMP_C_FLAGS}` / `${OpenMP_CXX_FLAGS}` — already set in parent's `CMAKE_C_FLAGS`/`CMAKE_CXX_FLAGS`
- `CMAKE_EXE_LINKER_FLAGS` with libomp — already set in parent

- [ ] **Step 2: Verify the change parses correctly**

Run: `cd /Users/markl_m/code/KAMEL && cmake -S . -B build-test 2>&1 | grep -i "ql-balance\|Building on"`

Expected: "Building on macOS" message from QL-Balance without any Homebrew-related warnings.

- [ ] **Step 3: Commit**

```bash
git add QL-Balance/CMakeLists.txt
git commit -m "build(QL-Balance): remove duplicate Homebrew libomp detection, inherit from parent"
```

---

### Task 5: Full Build Verification

**Files:** None (verification only)

- [ ] **Step 1: Clean build directory and run full configure**

```bash
cd /Users/markl_m/code/KAMEL
rm -rf build-test
cmake -S . -B build-test 2>&1 | tee /tmp/kamel-configure.log
```

Expected: Configure completes without errors. Check that compiler, HDF5, and OpenMP detection messages appear and use the correct cascade logic.

- [ ] **Step 2: Build the project**

```bash
cd /Users/markl_m/code/KAMEL
cmake --build build-test -j$(sysctl -n hw.ncpu) 2>&1 | tail -20
```

Expected: Build completes successfully (or at least gets past the CMake configure stage without dependency-related errors).

- [ ] **Step 3: Run tests**

```bash
cd /Users/markl_m/code/KAMEL/build-test
ctest --output-on-failure 2>&1 | tail -20
```

Expected: Tests pass (or at least no regressions from this change).

- [ ] **Step 4: Clean up test build directory**

```bash
rm -rf /Users/markl_m/code/KAMEL/build-test
```
