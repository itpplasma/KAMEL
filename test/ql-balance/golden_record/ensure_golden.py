#!/usr/bin/env python3
"""
Ensure golden record data exists by building from main branch if needed.

This module handles:
1. Cloning/pulling the main branch to a reference directory
2. Building the main branch
3. Running ql-balance.x with the main branch executable
4. Storing the output as golden.h5

The reference directory is cached locally to avoid unnecessary rebuilds.
"""

import shutil
import subprocess
from pathlib import Path
from setup_runfolder import setup_runfolder as setup_runfolder_external, get_output_hdf5_path

SCRIPT_DIR = Path(__file__).resolve().parent
MAIN_REF_DIR = SCRIPT_DIR / "main_ref"
GOLDEN_H5 = SCRIPT_DIR / "golden.h5"
RUNFOLDER_DIR = SCRIPT_DIR / "runfolder"
REPO_URL = "https://github.com/itpplasma/KAMEL"
CONFIG = "Release"


def run_cmd(cmd: list[str], cwd: Path | None = None, check: bool = True) -> str:
    """Run a command and return stdout."""
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=check)
    return result.stdout.strip()


def create_conftest() -> None:
    """Create conftest.py to exclude main_ref from pytest collection."""
    conftest = MAIN_REF_DIR / "conftest.py"
    conftest.write_text('collect_ignore_glob = ["**/*"]\n')


def clone_main_ref() -> None:
    """Clone the main branch to the reference directory."""
    print(f"Cloning main branch to {MAIN_REF_DIR}...")
    subprocess.run(
        [
            "git",
            "clone",
            "--branch",
            "main",
            "--single-branch",
            REPO_URL,
            str(MAIN_REF_DIR),
        ],
        check=True,
    )
    create_conftest()


def pull_main_ref() -> bool:
    """Pull latest changes. Returns True if there were updates."""
    print("Fetching latest changes from main...")
    run_cmd(["git", "fetch", "origin", "main"], cwd=MAIN_REF_DIR)

    local_head = run_cmd(["git", "rev-parse", "HEAD"], cwd=MAIN_REF_DIR)
    remote_head = run_cmd(["git", "rev-parse", "origin/main"], cwd=MAIN_REF_DIR)

    if local_head != remote_head:
        print(f"Updating from {local_head[:8]} to {remote_head[:8]}...")
        run_cmd(["git", "reset", "--hard", "origin/main"], cwd=MAIN_REF_DIR)
        return True

    print("Main branch is up to date.")
    return False


def build_main_ref() -> Path:
    """Build the main branch and return path to ql-balance.x executable."""
    print(f"Building main branch with CONFIG={CONFIG}...")
    subprocess.run(
        ["make", f"CONFIG={CONFIG}", "QL-Balance"],
        cwd=MAIN_REF_DIR,
        check=True,
    )
    return MAIN_REF_DIR / "build" / "install" / "bin" / "ql-balance.x"


def setup_runfolder() -> Path:
    """Set up the runfolder using the user-provided setup script."""
    # Clean up any existing runfolder
    if RUNFOLDER_DIR.exists():
        shutil.rmtree(RUNFOLDER_DIR)
    RUNFOLDER_DIR.mkdir(parents=True)

    setup_runfolder_external(RUNFOLDER_DIR)
    return RUNFOLDER_DIR


def run_ql_balance(executable: Path, runfolder: Path) -> Path:
    """Run ql-balance.x and return path to output HDF5 file."""
    print(f"Running {executable}...")
    result = subprocess.run(
        [str(executable)],
        cwd=runfolder,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(f"ql-balance.x failed with return code {result.returncode}")
        print(f"stdout:\n{result.stdout}")
        print(f"stderr:\n{result.stderr}")
        raise RuntimeError("ql-balance.x execution failed")

    return get_output_hdf5_path(runfolder)


def copy_to_golden(output_h5: Path) -> None:
    """Copy the output HDF5 to golden.h5."""
    print(f"Copying {output_h5} to {GOLDEN_H5}...")
    shutil.copy2(output_h5, GOLDEN_H5)


def ensure_golden() -> Path:
    """
    Ensure golden.h5 exists and is up to date with main branch.

    The golden record is rebuilt when the remote main branch has new commits.

    Returns the path to golden.h5.
    """
    needs_clone = not MAIN_REF_DIR.exists()
    needs_build = False
    needs_regenerate = not GOLDEN_H5.exists()

    if needs_clone:
        clone_main_ref()
        needs_build = True
        needs_regenerate = True
    else:
        updated = pull_main_ref()
        if updated:
            needs_build = True
            needs_regenerate = True

    executable = MAIN_REF_DIR / "build" / "install" / "bin" / "ql-balance.x"

    if needs_build or not executable.exists():
        executable = build_main_ref()
        needs_regenerate = True

    if needs_regenerate:
        runfolder = setup_runfolder()
        output_h5 = run_ql_balance(executable, runfolder)
        copy_to_golden(output_h5)

    if not GOLDEN_H5.exists():
        raise RuntimeError(f"Failed to create {GOLDEN_H5}")

    print(f"Golden record ready: {GOLDEN_H5}")
    return GOLDEN_H5


def main() -> None:
    """CLI entry point."""
    ensure_golden()


if __name__ == "__main__":
    main()
