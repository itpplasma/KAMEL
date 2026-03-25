# IO Refactor Design: Shared Logger and Data Verbosity

> **GitHub Issue:** #116 — Refactor IO in QL-Balance and KIM
>
> **Goal:** Replace ad-hoc `write(*,*)`/`print *` statements across both codes
> with a shared leveled logger module, and add an integer data-verbosity control
> for HDF5/file output.

---

## 1. Logger Module

New shared module at `common/logger/logger_m.f90`, built as a static library
(`kamel_logger`) and linked by both KIM and QL-Balance.

### Log Levels

| Constant       | Value | Purpose                                                       |
|----------------|-------|---------------------------------------------------------------|
| `LVL_SILENT`   |    -1 | No output                                                     |
| `LVL_RESULT`   |     0 | Final results only (banner, key computed quantities)          |
| `LVL_ERROR`    |     1 | Errors before `error stop`                                    |
| `LVL_WARNING`  |     2 | Warnings (sign mismatch, missing optional files, etc.)        |
| `LVL_INFO`     |     3 | Normal progress (**default**) — config summary, solver entry  |
| `LVL_DEBUG`    |     4 | Detailed diagnostics — per-resonance values, intermediates    |
| `LVL_TRACE`    |     5 | Everything — per-gridpoint dumps, inner-loop output           |

### API

```fortran
! Level control (called once after namelist read)
call set_log_level(level)
integer function get_log_level()

! Logging subroutines (string-only, level check is internal)
call log_error(msg)      ! writes to stderr, then calls error stop
call log_warning(msg)    ! writes to stderr
call log_info(msg)       ! writes to stdout
call log_debug(msg)      ! writes to stdout
call log_trace(msg)      ! writes to stdout
call log_result(msg)     ! writes to stdout (shown at LVL_RESULT and above)

! Format helpers — return formatted strings, overloaded for real(dp),
! integer, logical, and character value types.
character(len=256) function fmt_val(label, value, unit)
```

Output format follows NEO-RT conventions:

```
[ERROR] message        → stderr
[WARN ] message        → stderr
[INFO ] message        → stdout
[DEBUG] message        → stdout
[TRACE] message        → stdout
```

`log_result` prints the message without a prefix, for banners and key results.

## 2. Data Verbosity

A separate integer variable `data_verbosity` controlling how much HDF5/file
data gets written. Lives in each code's own config module, not in the shared
logger.

### Levels

| Value | Meaning              | What gets written                                              |
|-------|----------------------|----------------------------------------------------------------|
|     0 | Minimal              | Final-state output only (the HDF5 downstream tools need)       |
|     1 | Standard (**default**) | + per-timestep profiles, transport coefficients                |
|     2 | Detailed             | + per-mode fields (Es, Br), resonance diagnostics, currents   |
|     3 | Full                 | + Jacobian data, intermediate solver states, gridpoint dumps   |

### Configuration

Each code reads `log_level` and `data_verbosity` from its own namelist:

- **QL-Balance:** `balance_conf.nml` → `BALANCENML` namelist
- **KIM:** `KIM_config.nml` → `KIM_CONFIG` namelist

After reading, each code calls `set_log_level(log_level)`.

The existing format-toggle variables (`ihdf5IO` in QL-Balance, `hdf5_output`
in KIM) are kept — they control *format*, not *amount*.

## 3. Removed Variables (Clean Break)

| Code        | Removed                                        | Replaced by                    |
|-------------|------------------------------------------------|--------------------------------|
| QL-Balance  | `debug_mode` (logical)                         | `log_level`                    |
| QL-Balance  | `diagnostics_output` (logical)                 | `data_verbosity`               |
| KIM         | `fdebug` (integer)                             | `log_level`                    |
| KIM         | `fstatus` (integer)                            | `log_level`                    |
| KIM         | `fdiagnostics` (integer)                       | `data_verbosity`               |

No backward-compatibility shims — old namelists must be updated.

## 4. Migration Strategy

### Classification of existing output (~880 statements across ~47 files)

1. **Config display** (`read_config.f90`, `config_display.f90`) →
   `log_info(fmt_val(...))`
2. **Progress/status** (`"Status: solve poisson"`) → `log_info(msg)`
3. **Warnings** (`"WARNING: ..."`) → `log_warning(msg)`
4. **Debug prints** (currently guarded by `if (debug_mode)` / `if (fdebug)`) →
   `log_debug(msg)` or `log_trace(msg)`
5. **Error messages** before `error stop` → `log_error(msg)`
6. **HDF5 data writes** → add `if (data_verbosity >= N)` guards

### What NOT to change

- File I/O to unit numbers (namelist writes, data files)
- ASCII art banner → `log_result` (shown at verbosity 0)
- Third-party code (libcerf, ddeabm, slatec)

### Order of work

1. Create `common/logger/logger_m.f90` + CMake
2. Wire into both codes' CMakeLists
3. Migrate KIM (~27 files, ~565 statements)
4. Migrate QL-Balance (~20 files, ~316 statements)
5. Update namelists and Python interface (`balance_conf.py`,
   `balance_interface.py`)
6. Update golden record test config

## 5. File Changes

### New files

- `common/logger/logger_m.f90`
- `common/logger/CMakeLists.txt`

### Config/build modifications

- `CMakeLists.txt` — add `common/logger` subdirectory
- `KIM/src/CMakeLists.txt` — link `kamel_logger`
- `QL-Balance/CMakeLists.txt` — link `kamel_logger`
- `KIM/src/setup/config_mod.f90` — replace `fdebug`/`fstatus`/`fdiagnostics`
- `KIM/src/general/read_config.f90` — read new variables, call `set_log_level`
- `QL-Balance/src/base/control_mod.f90` — replace `debug_mode`/`diagnostics_output`
- `QL-Balance/src/base/read_config.f90` — read new variables, call `set_log_level`
- `QL-Balance/namelists/balance_conf.nml`
- `python/balance_interface/balance_interface.py`

### Source migration (~47 files)

- KIM: ~27 files
- QL-Balance: ~20 files
