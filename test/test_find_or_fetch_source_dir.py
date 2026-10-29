"""find_or_fetch must honour CMake's FETCHCONTENT_SOURCE_DIR_<DEP> override."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

UTIL = Path(__file__).resolve().parents[1] / "cmake" / "Util.cmake"


@pytest.mark.skipif(shutil.which("cmake") is None, reason="cmake not found")
def test_pinned_sha_uses_fetchcontent_source_dir(tmp_path):
    libneo = tmp_path / "libneo"
    libneo.mkdir()
    (libneo / "CMakeLists.txt").write_text("project(libneo NONE)\n")
    project = tmp_path / "project"
    project.mkdir()
    (project / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.24)\n"
        "project(consumer NONE)\n"
        f'include("{UTIL.as_posix()}")\n'
        "find_or_fetch(libneo)\n"
    )
    # Forbid network transports so a fetch attempt fails instead of succeeding.
    env = dict(os.environ, GIT_ALLOW_PROTOCOL="file")
    result = subprocess.run(
        [
            "cmake",
            "-S",
            str(project),
            "-B",
            str(tmp_path / "build"),
            "-DLIBNEO_REF=d03789a1fef15e188752614328dd50ea9854919f",
            f"-DFETCHCONTENT_SOURCE_DIR_LIBNEO={libneo}",
        ],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"Using libneo from {libneo}" in result.stdout
