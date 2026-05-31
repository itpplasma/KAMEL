# Layered Dependency Discovery for macOS

**Date:** 2026-03-31
**Status:** Approved

## Problem

The KAMEL build system on macOS/Darwin relies on Homebrew for compiler detection, OpenMP/libomp, and HDF5 discovery. This prevents building in environments without Homebrew, such as Nix or manual installations.

## Design Principle

Standard CMake discovery first (respects `$PATH`, `$CMAKE_PREFIX_PATH`, user-provided variables), Homebrew as fallback only when needed.

## Changes

### 1. Compiler Detection (`CMakeLists.txt:4-22`)

**Cascade:**
1. User-provided `-DCMAKE_C_COMPILER` / `-DCMAKE_CXX_COMPILER` (already works, skip detection)
2. `find_program(clang/clang++)` with no path restrictions — picks up whatever is on `$PATH`
3. Homebrew LLVM paths (`/opt/homebrew/opt/llvm/bin`, `/usr/local/opt/llvm/bin`) as fallback

`CMAKE_PREFIX_PATH`: only append Homebrew paths if `brew` is found. Do not overwrite user/environment-provided paths.

### 2. HDF5 Discovery (`CMakeLists.txt:32-48`)

**Cascade:**
1. User-provided `-DHDF5_DIR` (already works, skip)
2. Check if HDF5 is findable via standard `CMAKE_PREFIX_PATH` / system paths
3. Fall back to `brew --prefix hdf5` if not found
4. If still not found, `FetchHDF5.cmake` builds from source (existing behavior)

### 3. OpenMP/libomp (`CMakeLists.txt:57-81`)

**Cascade:**
1. `find_package(OpenMP REQUIRED)` — if the compiler knows where OpenMP is (Nix, configured LLVM), this is sufficient
2. If `find_package` fully populated flags and libraries, use those directly
3. Only if link paths are incomplete, fall back to `brew --prefix libomp`

### 4. QL-Balance OpenMP (`QL-Balance/CMakeLists.txt:6-28`)

Remove duplicate Homebrew libomp detection block. The parent `CMakeLists.txt` already sets `CMAKE_*_FLAGS` with OpenMP configuration, and these are inherited via `add_subdirectory`. Remove hardcoded `-L/opt/homebrew/lib`.

### 5. FindSuiteSparse (`cmake/FindSuiteSparse.cmake`)

No changes. Homebrew paths in the search list are harmless fallbacks alongside other standard paths. User hints already allow overriding.

## Files Modified

| File | Scope |
|---|---|
| `CMakeLists.txt` | Compiler detection, HDF5 discovery, OpenMP/libomp |
| `QL-Balance/CMakeLists.txt` | Remove duplicate Homebrew libomp logic |

## Not Changed

| File | Reason |
|---|---|
| `cmake/FindSuiteSparse.cmake` | Search paths are already broad and non-exclusive |
| `cmake/Fetch*.cmake` | These handle source-build fallbacks, unrelated to Homebrew |
