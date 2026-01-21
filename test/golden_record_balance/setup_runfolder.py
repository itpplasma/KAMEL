#!/usr/bin/env python3
"""
Setup script for the QL-Balance golden record test runfolder.

==============================================================================
USER TODO: Implement this file!

This script must create a complete runfolder with all necessary input files
for running ql-balance.x. The runfolder structure should match the example at:
    ../runfolders/runfolder_33353_2900/

Required files typically include:
    - balance_conf.nml (QL-Balance configuration namelist)
    - background.in (background plasma parameters)
    - eigmode.in (eigenmode configuration)
    - output.in (output configuration)
    - profiles/ directory with:
        - n.dat (electron density)
        - Te.dat (electron temperature)
        - Ti.dat (ion temperature)
        - q.dat (safety factor)
        - Er.dat (radial electric field)
        - Vz.dat (toroidal rotation) - optional
        - Vth.dat (thermal velocity) - optional
        - Da.dat (anomalous diffusion) - optional
    - Equilibrium data file (e.g., *_equil_r_q_psi.dat)
    - Input HDF5 file (e.g., *_inp.hdf5)

You can either:
1. Copy files from a reference runfolder
2. Generate files programmatically
3. Use a combination of both

==============================================================================
"""

from pathlib import Path


def setup_runfolder(target_dir: Path) -> None:
    """
    Set up a runfolder in the target directory with all necessary input files.

    This function is called by ensure_golden.py and test_golden_record.py to
    create the runfolder before running ql-balance.x.

    Args:
        target_dir: The directory where the runfolder should be created.
                   This directory already exists when this function is called.

    Raises:
        RuntimeError: If setup fails for any reason.

    Example implementation:
        ```python
        import shutil

        # Copy from reference runfolder
        reference = Path("/path/to/reference/runfolder")
        for item in reference.iterdir():
            if item.is_file():
                shutil.copy2(item, target_dir / item.name)
            elif item.is_dir() and item.name in ["profiles", "vacuum", "flre"]:
                shutil.copytree(item, target_dir / item.name)
        ```
    """
    # TODO: Implement this function!
    raise NotImplementedError(
        "setup_runfolder() is not implemented. "
        "Please edit test/golden_record_balance/setup_runfolder.py"
    )


def get_output_hdf5_path(runfolder: Path) -> Path:
    """
    Return the path to the output HDF5 file produced by ql-balance.x.

    This function is called after ql-balance.x has finished running to locate
    the output file for comparison.

    Args:
        runfolder: The runfolder directory where ql-balance.x was run.

    Returns:
        Path to the output HDF5 file.

    Raises:
        FileNotFoundError: If the output file doesn't exist.

    Example implementation:
        ```python
        # Output file is typically in out/ subdirectory
        output_file = runfolder / "out" / "balance_output.h5"
        if not output_file.exists():
            raise FileNotFoundError(f"Output file not found: {output_file}")
        return output_file
        ```
    """
    # TODO: Implement this function!
    raise NotImplementedError(
        "get_output_hdf5_path() is not implemented. "
        "Please edit test/golden_record_balance/setup_runfolder.py"
    )
