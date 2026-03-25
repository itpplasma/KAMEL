# IO Refactor Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace ad-hoc write/print statements in KIM and QL-Balance with a shared leveled logger module, and add integer data-verbosity control for HDF5 output.

**Architecture:** A new `common/logger/logger_m.f90` static library provides leveled logging (SILENT..TRACE) and format helpers. Each code reads `log_level` and `data_verbosity` from its own namelist and calls `set_log_level()`. Old boolean/integer verbosity flags are removed (clean break).

**Tech Stack:** Fortran 2008, CMake/Ninja, CTest, f90nml (Python)

**Design document:** `docs/plans/2026-03-25-io-refactor-design.md`

---

## Task 1: Create Logger Module

**Files:**
- Create: `common/logger/logger_m.f90`
- Create: `common/logger/CMakeLists.txt`

**Step 1: Write `common/logger/logger_m.f90`**

```fortran
module logger_m
    implicit none
    private

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: FMT_LEN = 256

    ! Log levels
    integer, public, parameter :: LVL_SILENT  = -1
    integer, public, parameter :: LVL_RESULT  =  0
    integer, public, parameter :: LVL_ERROR   =  1
    integer, public, parameter :: LVL_WARNING =  2
    integer, public, parameter :: LVL_INFO    =  3
    integer, public, parameter :: LVL_DEBUG   =  4
    integer, public, parameter :: LVL_TRACE   =  5

    integer :: current_level = LVL_INFO

    public :: set_log_level, get_log_level
    public :: log_error, log_warning, log_info, log_debug, log_trace, log_result

    interface fmt_val
        module procedure fmt_val_real, fmt_val_int, fmt_val_logical, fmt_val_str
    end interface fmt_val
    public :: fmt_val

contains

    subroutine set_log_level(level)
        integer, intent(in) :: level
        current_level = max(LVL_SILENT, min(LVL_TRACE, level))
    end subroutine set_log_level

    function get_log_level() result(level)
        integer :: level
        level = current_level
    end function get_log_level

    subroutine log_error(msg)
        character(*), intent(in) :: msg
        write(0, '(A)') '[ERROR] ' // trim(msg)
        error stop
    end subroutine log_error

    subroutine log_warning(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_WARNING) write(0, '(A)') '[WARN ] ' // trim(msg)
    end subroutine log_warning

    subroutine log_info(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_INFO) write(6, '(A)') '[INFO ] ' // trim(msg)
    end subroutine log_info

    subroutine log_debug(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_DEBUG) write(6, '(A)') '[DEBUG] ' // trim(msg)
    end subroutine log_debug

    subroutine log_trace(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_TRACE) write(6, '(A)') '[TRACE] ' // trim(msg)
    end subroutine log_trace

    subroutine log_result(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_RESULT) write(6, '(A)') trim(msg)
    end subroutine log_result

    ! --- Format helpers ---

    function fmt_val_real(label, value, unit) result(s)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        character(len=20) :: vbuf
        write(vbuf, '(ES15.8)') value
        if (present(unit)) then
            s = trim(label) // ' = ' // trim(adjustl(vbuf)) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // trim(adjustl(vbuf))
        end if
    end function fmt_val_real

    function fmt_val_int(label, value, unit) result(s)
        character(*), intent(in) :: label
        integer, intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        character(len=20) :: vbuf
        write(vbuf, '(I0)') value
        if (present(unit)) then
            s = trim(label) // ' = ' // trim(adjustl(vbuf)) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // trim(adjustl(vbuf))
        end if
    end function fmt_val_int

    function fmt_val_logical(label, value, unit) result(s)
        character(*), intent(in) :: label
        logical, intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        if (present(unit)) then
            s = trim(label) // ' = ' // merge('T', 'F', value) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // merge('T', 'F', value)
        end if
    end function fmt_val_logical

    function fmt_val_str(label, value, unit) result(s)
        character(*), intent(in) :: label
        character(*), intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        if (present(unit)) then
            s = trim(label) // ' = ' // trim(value) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // trim(value)
        end if
    end function fmt_val_str

end module logger_m
```

**Step 2: Write `common/logger/CMakeLists.txt`**

```cmake
set(LOGGER_MODULE_DIR ${CMAKE_BINARY_DIR}/modules/logger)
file(MAKE_DIRECTORY ${LOGGER_MODULE_DIR})

add_library(kamel_logger STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/logger_m.f90
)

set_target_properties(kamel_logger PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${LOGGER_MODULE_DIR}
    LINKER_LANGUAGE Fortran
)

target_include_directories(kamel_logger
    PUBLIC ${LOGGER_MODULE_DIR}
)
```

**Step 3: Build and verify**

Run: `make clean && make all`
Expected: Build succeeds (logger is built but not yet used)

**Step 4: Commit**

```
git add common/logger/
git commit -m "feat: add shared logger module with leveled output and format helpers"
```

---

## Task 2: Wire Logger into CMake for Both Codes

**Files:**
- Modify: `CMakeLists.txt:103-105` (add subdirectory)
- Modify: `KIM/src/CMakeLists.txt:147-156` (link library)
- Modify: `QL-Balance/CMakeLists.txt:145-167` (link library)

**Step 1: Add `common/logger` subdirectory to top-level CMakeLists.txt**

In `CMakeLists.txt`, after line 103 (`add_subdirectory(common/math)`), add:

```cmake
add_subdirectory(common/logger)
```

So lines 103-106 become:

```cmake
add_subdirectory(common/math)
add_subdirectory(common/logger)
add_subdirectory(common/hdf5_tools)
add_subdirectory(common/equil)
```

**Step 2: Link `kamel_logger` into KIM**

In `KIM/src/CMakeLists.txt`, add `kamel_logger` to the `target_link_libraries(KIM_lib PUBLIC ...)` block (after `kamel_hdf5_tools`).

**Step 3: Link `kamel_logger` into QL-Balance**

In `QL-Balance/CMakeLists.txt`, add `kamel_logger` to the `target_link_libraries(ql-balance_lib PUBLIC ...)` block (after `kamel_hdf5_tools`).

**Step 4: Build and verify**

Run: `make all`
Expected: Build succeeds, logger library is linked into both codes.

**Step 5: Commit**

```
git add CMakeLists.txt KIM/src/CMakeLists.txt QL-Balance/CMakeLists.txt
git commit -m "build: wire kamel_logger into KIM and QL-Balance"
```

---

## Task 3: Add Config Variables and Namelist Wiring for KIM

**Files:**
- Modify: `KIM/src/setup/config_mod.f90:34` (replace fdebug/fstatus/fdiagnostics)
- Modify: `KIM/src/general/read_config.f90:27-29` (update KIM_IO namelist)
- Modify: `KIM/nmls/KIM_config.nml:34-36` (update defaults)
- Modify: `KIM/src/util/IO_collection.f90:123-131` (update HDF5 output)

**Step 1: Update `config_mod.f90`**

Replace the declaration at line 34:
```fortran
integer :: fdebug, fstatus, fdiagnostics
```
with:
```fortran
integer :: log_level = 3         ! maps to LVL_INFO
integer :: data_verbosity = 1    ! standard output
```

**Step 2: Update `read_config.f90` namelist definition**

In the `KIM_IO` namelist (lines 27-29), replace `fdebug, fstatus, fdiagnostics` with `log_level, data_verbosity`.

After the namelist read (after line 69 `close(unit = 77)`), add:
```fortran
use logger_m, only: set_log_level
call set_log_level(log_level)
```

**Step 3: Update `KIM_config.nml`**

Replace lines 34-36:
```
    fdebug = 1
    fstatus = 1
    fdiagnostics = 0
```
with:
```
    log_level = 3
    data_verbosity = 1
```

**Step 4: Update `IO_collection.f90` HDF5 output**

Replace the `fdebug`/`fstatus`/`fdiagnostics` writes (lines 123, 125, 131) with:
```fortran
call h5_add(h5grpid, 'log_level', log_level, 'Log verbosity level', 'i')
call h5_add(h5grpid, 'data_verbosity', data_verbosity, 'Data output verbosity', 'i')
```

**Step 5: Build and verify**

Run: `make all`
Expected: Build succeeds. KIM reads new variables from namelist.

**Step 6: Commit**

```
git commit -am "feat(KIM): replace fdebug/fstatus/fdiagnostics with log_level and data_verbosity"
```

---

## Task 4: Add Config Variables and Namelist Wiring for QL-Balance

**Files:**
- Modify: `QL-Balance/src/base/control_mod.f90:18,22` (replace debug_mode/diagnostics_output)
- Modify: `QL-Balance/src/base/read_config.f90:3-8,23-32,37-38,72,74` (update namelist + printing)
- Modify: `QL-Balance/namelists/balance_conf.nml:41,46` (update defaults)

**Step 1: Update `control_mod.f90`**

Replace line 18 (`logical :: debug_mode`) with:
```fortran
integer :: log_level = 3         ! maps to LVL_INFO
```

Replace line 22 (`logical :: diagnostics_output`) with:
```fortran
integer :: data_verbosity = 1    ! standard output
```

**Step 2: Update `read_config.f90`**

In the `use control_mod` import (lines 3-8): replace `debug_mode` and `diagnostics_output` references with `log_level, data_verbosity`. Remove `debug_mode` from the import list entirely.

In the namelist definition (lines 23-32): replace `diagnostics_output` with `data_verbosity` and `debug_mode` with `log_level`. Remove `suppression_mode` only if it was tied to diagnostics (check first — it may be independent).

After the namelist read (after line 38), add:
```fortran
use logger_m, only: set_log_level
call set_log_level(log_level)
```

Update the print statements (lines 72, 74) to use the new variable names:
```fortran
write (*, "(A,I0)") "    log_level = ", log_level
write (*, "(A,I0)") "    data_verbosity = ", data_verbosity
```

**Step 3: Update `balance_conf.nml`**

Replace:
```
 diagnostics_output = .false. ,
```
with:
```
 data_verbosity = 1 ,
```

Replace:
```
 debug_mode = .false. ,
```
with:
```
 log_level = 3 ,
```

**Step 4: Build — expect errors**

Run: `make all`
Expected: **Compilation errors** in files that still reference `debug_mode` and `diagnostics_output`. This is expected — Task 6 will fix QL-Balance source files. For now just verify the config/namelist changes compile in `control_mod.f90` and `read_config.f90` themselves.

**Step 5: Commit (WIP)**

```
git commit -am "feat(QL-Balance): replace debug_mode/diagnostics_output with log_level and data_verbosity

WIP: source files still reference old variables, will be migrated in subsequent commits"
```

---

## Task 5: Migrate KIM Source Files to Logger

**Files to modify (~12 files with active control variable usage):**
- `KIM/src/general/config_display.f90:119-121,338-374`
- `KIM/src/kernels/kernel.f90:80,129,372,473-475,847,947-949`
- `KIM/src/electrostatic_poisson/solve_poisson.f90:11,34-35,51,168,188`
- `KIM/src/background_equilibrium/calculate_equil.f90:74,180,185`
- `KIM/src/background_equilibrium/species_mod.f90:1171,1184,1275,1373,1377`
- `KIM/src/electromagnetic/electromagnetic_solver.f90:47,222`
- `KIM/src/grid/grid_mod.f90:197,274`
- `KIM/src/grid/gengrid.f90:6,37`
- `KIM/src/math/quadpack_integration_m.f90:204,274,299,368,393,462`

Also migrate any remaining bare `write(*,*)` and `print *` statements in these files to appropriate logger calls.

**Migration rules:**
- `if (fstatus >= 1) write(*,*) 'Status: ...'` → `call log_info('...')`
- `if (fstatus >= 2) write(*,*) ...` → `call log_debug('...')`
- `if (fdebug == 1) ...` → `call log_debug('...')`
- `if (fdebug >= 2) ...` → `call log_debug('...')` or `call log_trace('...')`
- `if (fdebug == 3) ...` → `call log_trace('...')`
- `if (fdiagnostics == 3) ...` → `if (data_verbosity >= 3) ...` (keep as data guard, not logger)
- Bare `write(*,*) 'Warning: ...'` → `call log_warning('...')`
- Bare `write(*,*) 'Error: ...'` before `error stop` → `call log_error('...')`
- Bare `write(*,*) ...` status messages → `call log_info('...')`
- Formatted output like `write(*,"(A,ES15.8,A)") "label", val, "unit"` → `call log_info(fmt_val('label', val, 'unit'))`

**For each file:**
1. Replace `use config_m, only: ... fdebug ...` with `use logger_m, only: log_info, log_debug, ...` (and `use config_m, only: ... data_verbosity ...` where `fdiagnostics` was used for data guards)
2. Convert write/print statements per the rules above
3. Remove bare `if (fdebug ...)` / `if (fstatus ...)` guards — the logger checks internally

**Step 1: Migrate config_display.f90**

Update lines 119-121 to use `log_level` and `data_verbosity` display strings. Update or remove `get_debug_string`/`get_status_string` helpers (lines 338-374) — replace with a single helper or inline.

**Step 2: Migrate remaining KIM source files**

Work through each file listed above, applying the migration rules.

**Step 3: Build and verify**

Run: `make all`
Expected: Build succeeds. No references to `fdebug`, `fstatus`, or `fdiagnostics` remain.

**Step 4: Verify no old references remain**

Run: `grep -rn 'fdebug\|fstatus\|fdiagnostics' KIM/src/`
Expected: No matches (except possibly comments).

**Step 5: Run tests**

Run: `make test`
Expected: All KIM tests pass.

**Step 6: Commit**

```
git commit -am "refactor(KIM): migrate all IO to shared logger module"
```

---

## Task 6: Migrate QL-Balance Source Files to Logger

**Files to modify (~17 files referencing debug_mode, ~9 referencing diagnostics_output):**

Key files:
- `QL-Balance/src/base/wave_code_data_64bit.f90` (heaviest: ~20 debug_mode refs)
- `QL-Balance/src/base/ramp_coil.f90` (~25 debug_mode refs)
- `QL-Balance/src/base/paramscan.f90` (~12 debug_mode refs)
- `QL-Balance/src/base/SingleStep.f90` (~8 debug_mode refs)
- `QL-Balance/src/base/writeData.f90` (~8 debug_mode + diagnostics_output refs)
- `QL-Balance/src/base/get_dql.f90`
- `QL-Balance/src/base/diff_coeffs.f90`
- `QL-Balance/src/base/time_evolution.f90`
- `QL-Balance/src/base/calc_current_densities.f90`
- `QL-Balance/src/base/plasma_parameters.f90`
- `QL-Balance/src/base/h5mod.f90`
- `QL-Balance/src/base/balance_eqs_source_terms.f90`
- `QL-Balance/src/base/kim_wave_code_adapter.f90`
- `QL-Balance/src/base/integrate_parallel_current.f90`
- `QL-Balance/src/base/transp_coeffs_mod.f90`
- `QL-Balance/src/stellarator/balance_eqs_source_terms_stell.f90`
- `QL-Balance/src/stellarator/time_evol_stell.f90`

**Migration rules:**
- `if (debug_mode) write(*,*) ...` → `call log_debug('...')`
- `if (debug_mode) print *, ...` → `call log_debug('...')`
- `if (diagnostics_output) call write_..._to_hdf5(...)` → `if (data_verbosity >= 2) call write_..._to_hdf5(...)`
- Bare `write(*,*) 'WARNING: ...'` → `call log_warning('...')`
- `error stop "msg"` → `call log_error('msg')` (which calls error stop internally)
- Bare status `write(*,*)` → `call log_info('...')`
- Formatted values → `call log_info(fmt_val(...))`

**For each file:**
1. Replace `use control_mod, only: ... debug_mode ...` with `use logger_m, only: log_debug, ...` (and `use control_mod, only: ... data_verbosity ...` where `diagnostics_output` was used)
2. Convert write/print statements per the rules above
3. Remove bare `if (debug_mode)` guards — the logger checks internally

**Step 1: Migrate files in batches**

Start with the most-referenced files (wave_code_data_64bit, ramp_coil, paramscan) then work through the rest.

**Step 2: Handle read_config.f90 print block (lines 41-98)**

The large config-display block of ~50 write statements should be converted to `log_info(fmt_val(...))` calls. For example:
```fortran
! Before:
write (*, "(A,ES15.8,A)") "    B_tor = ", btor, " G"
! After:
call log_info(fmt_val('    B_tor', btor, 'G'))
```

**Step 3: Build and verify**

Run: `make all`
Expected: Build succeeds. No references to `debug_mode` or `diagnostics_output` remain.

**Step 4: Verify no old references remain**

Run: `grep -rn 'debug_mode\|diagnostics_output' QL-Balance/src/`
Expected: No matches (except possibly comments).

**Step 5: Run tests**

Run: `make test`
Expected: All tests pass.

**Step 6: Commit**

```
git commit -am "refactor(QL-Balance): migrate all IO to shared logger module"
```

---

## Task 7: Update Python Interface and Golden Record Test

**Files:**
- Modify: `python/balance_interface/balance_interface.py` — update any references to old namelist keys
- Modify: `test/ql-balance/golden_record/setup_runfolder.py` — update namelist key if `diagnostics_output` is set there

**Step 1: Update `setup_runfolder.py`**

In `setup_runfolder.py`, line 96 currently sets:
```python
bi.conf.conf["balancenml"]["diagnostics_output"] = True
```
Replace with:
```python
bi.conf.conf["balancenml"]["data_verbosity"] = 2
```

Also check if `debug_mode` is set anywhere in the test setup and replace with `log_level`.

**Step 2: Check `balance_interface.py`**

Search for any hardcoded references to `debug_mode` or `diagnostics_output` in namelist manipulation. The Python `debug` flag (constructor parameter) is a Python-side flag and should NOT be changed — it's independent.

**Step 3: Update `suppression_mode` reference**

In `setup_runfolder.py` line 97:
```python
bi.conf.conf["balancenml"]["suppression_mode"] = False
```
`suppression_mode` is independent of logging — keep it unchanged.

**Step 4: Run golden record test locally (if possible)**

Run: `make test`
Expected: All tests pass including golden record.

**Step 5: Commit**

```
git commit -am "fix(test): update golden record and Python interface for new IO variables"
```

---

## Task 8: Final Verification and Cleanup

**Step 1: Full clean build**

Run: `make clean && make all`
Expected: Build succeeds with no warnings related to unused variables.

**Step 2: Run all tests**

Run: `make test`
Expected: All tests pass.

**Step 3: Verify no old variable references remain**

Run these greps — all should return no matches (except comments):
```bash
grep -rn 'fdebug\|fstatus\|fdiagnostics' KIM/src/ --include='*.f90'
grep -rn 'debug_mode\|diagnostics_output' QL-Balance/src/ --include='*.f90'
grep -rn 'debug_mode\|diagnostics_output' python/ --include='*.py'
```

**Step 4: Verify logger is the only console output path**

Run: `grep -rn 'write(\*' KIM/src/ QL-Balance/src/base/ --include='*.f90' | grep -v '!' | wc -l`
Expected: Zero or near-zero (only file I/O to specific units should remain, no `write(*,*)` to stdout).

**Step 5: Commit any final cleanup**

```
git commit -am "chore: final cleanup of IO refactor"
```
