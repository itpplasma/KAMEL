# Golden Record Test Framework for QL-Balance

This document describes the golden record regression test framework for QL-Balance, modeled after the NEO-RT implementation at `/home/andreas/git/master/neo-rt/test/golden_record/`.

## Overview

The golden record test compares QL-Balance output from the **current branch** against a **golden reference** generated from the `main` branch. This catches regressions where code changes inadvertently alter numerical results.

### How It Works

1. **Golden Generation** (`ensure_golden.py`):
   - Clones the `main` branch to `test/golden_record_balance/main_ref/`
   - Builds `ql-balance.x` from main
   - Sets up a runfolder using `setup_runfolder.py`
   - Runs the main-branch executable
   - Copies the output HDF5 to `golden.h5`

2. **Test Execution** (`test_golden_record.py`):
   - Sets up a fresh runfolder
   - Runs the current-branch `ql-balance.x`
   - Compares specified quantities against `golden.h5` using tolerances

3. **Comparison** (`compare_hdf5.py`):
   - Compares HDF5 datasets by path
   - Supports configurable `rtol` (relative) and `atol` (absolute) tolerances
   - Reports max absolute/relative differences on failure

## Directory Structure

```
test/golden_record_balance/
├── ensure_golden.py         # Clones main, builds, generates golden.h5
├── compare_hdf5.py          # HDF5 comparison utilities
├── test_golden_record.py    # Main pytest test file
├── setup_runfolder.py       # USER TODO: Implement runfolder setup
├── quantities_to_compare.py # USER TODO: Define quantities to compare
├── conftest.py              # Pytest configuration
├── .gitignore               # Ignores generated artifacts
│
├── main_ref/                # [generated] Cloned main branch
├── runfolder/               # [generated] Working runfolder
└── golden.h5                # [generated] Golden reference data
```

## Makefile Targets

```bash
make test      # Runs: ctest + golden + pytest (full test suite)
make ctest     # Runs CTest unit tests only
make golden    # Ensures golden.h5 exists (clones/builds main if needed)
make pytest    # Runs all pytest tests
make clean     # Removes build/ and golden record artifacts
```

## What Needs to Be Implemented

### 1. `setup_runfolder.py`

This file needs two functions implemented:

#### `setup_runfolder(target_dir: Path) -> None`

Creates a complete runfolder with all necessary input files for `ql-balance.x`.

**Reference runfolder:** `/home/andreas/git/master/runfolders/runfolder_33353_2900/`

**Required files typically include:**
- `balance_conf.nml` - QL-Balance configuration namelist
- `background.in` - Background plasma parameters
- `eigmode.in` - Eigenmode configuration
- `output.in` - Output configuration
- `profiles/` directory with:
  - `n.dat` - Electron density
  - `Te.dat` - Electron temperature
  - `Ti.dat` - Ion temperature
  - `q.dat` - Safety factor
  - `Er.dat` - Radial electric field
  - `Vz.dat` - Toroidal rotation (optional)
  - `Vth.dat` - Thermal velocity (optional)
  - `Da.dat` - Anomalous diffusion (optional)
- Equilibrium data file (e.g., `*_equil_r_q_psi.dat`)
- Input HDF5 file (e.g., `*_inp.hdf5`)

**Example implementation:**
```python
import shutil
from pathlib import Path

REFERENCE_RUNFOLDER = Path("/home/andreas/git/master/runfolders/runfolder_33353_2900")

def setup_runfolder(target_dir: Path) -> None:
    # Copy all necessary files from reference
    for item in REFERENCE_RUNFOLDER.iterdir():
        if item.name in ["ql-balance.x", "out", ".claude", ".mypy_cache", ".gdb_history"]:
            continue  # Skip these
        if item.is_file():
            shutil.copy2(item, target_dir / item.name)
        elif item.is_dir():
            shutil.copytree(item, target_dir / item.name)

    # Create output directory
    (target_dir / "out").mkdir(exist_ok=True)
```

#### `get_output_hdf5_path(runfolder: Path) -> Path`

Returns the path to the output HDF5 file produced by `ql-balance.x`.

**Example implementation:**
```python
def get_output_hdf5_path(runfolder: Path) -> Path:
    # Adjust this path based on actual output location
    output_file = runfolder / "out" / "balance_output.h5"
    if not output_file.exists():
        # Try to find any .h5 file in out/
        h5_files = list((runfolder / "out").glob("*.h5"))
        if h5_files:
            return h5_files[0]
        raise FileNotFoundError(f"No HDF5 output found in {runfolder / 'out'}")
    return output_file
```

### 2. `quantities_to_compare.py`

Define the HDF5 dataset paths and tolerances for quantities to compare.

**To discover available paths in your HDF5 file:**
```bash
h5dump -H your_output.h5
```
or
```python
import h5py
with h5py.File("your_output.h5", "r") as f:
    f.visititems(lambda name, obj: print(name))
```

**Example configuration:**
```python
from compare_hdf5 import QuantitySpec

QUANTITIES_TO_COMPARE = [
    # Basic profiles
    QuantitySpec("/profiles/Te"),
    QuantitySpec("/profiles/Ti"),
    QuantitySpec("/profiles/n"),

    # Transport coefficients (may need relaxed tolerances)
    QuantitySpec("/output/D11", rtol=1e-6, atol=1e-12),
    QuantitySpec("/output/D12", rtol=1e-6, atol=1e-12),

    # Torque
    QuantitySpec("/output/torque", rtol=1e-6),

    # Add more quantities as needed...
]
```

**QuantitySpec parameters:**
- `path`: HDF5 path to the dataset (e.g., `/group/subgroup/dataset`)
- `rtol`: Relative tolerance (default: `1e-8`)
- `atol`: Absolute tolerance (default: `1e-15`)

## Running the Tests

```bash
# Build first
make QL-Balance

# Run just the golden record test
make golden && pytest test/golden_record_balance/ -v

# Run full test suite
make test
```

## Comparison with NEO-RT Implementation

| Aspect | NEO-RT | QL-Balance |
|--------|--------|------------|
| Output format | Multiple text files | Single HDF5 file |
| Test cases | 9 parametrized s-values | Single runfolder |
| Input files | Minimal template + symlinks | Full runfolder copy |
| Build config | `CONFIG=Fast` (required for FPE) | `CONFIG=Release` |
| Comparison | Text file parsing | HDF5 dataset paths |

## Key Files Reference

- **NEO-RT golden record:** `/home/andreas/git/master/neo-rt/test/golden_record/`
- **Reference runfolder:** `/home/andreas/git/master/runfolders/runfolder_33353_2900/`
- **QL-Balance source:** `/home/andreas/git/master/KAMEL/QL-Balance/`

## Next Steps

1. Implement `setup_runfolder()` in `test/golden_record_balance/setup_runfolder.py`
2. Implement `get_output_hdf5_path()` in the same file
3. Identify output HDF5 structure and populate `quantities_to_compare.py`
4. Test with `make golden && pytest test/golden_record_balance/ -v`
5. Adjust tolerances as needed for numerical stability

## Notes

- The `main_ref/` directory is cached locally to avoid unnecessary rebuilds
- Golden is regenerated only when main branch has new commits
- Tests are marked with `@pytest.mark.golden_record` for selective running
- MPI is used to run `ql-balance.x` (`mpirun -np 1`)
