"""Run every tutorial in tutorials/ exactly as written.

Each tutorials/<name>/README.md is turned into one bash script: `sh` code
blocks run as they are, `python` code blocks run through `python -` in the
shell's current directory. The script starts in the repository root with
KAMEL set to it and HOME set to a temporary directory, so the tutorials write
their runs below that temporary home.

The tutorials need the built executables and the `kim` and KAMELpy Python
packages. Without them the tests are skipped, unless
KAMEL_REQUIRE_TUTORIALS=1 is set, as in CI.
"""

from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL_DIR = REPO_ROOT / "tutorials"
BIN_DIR = REPO_ROOT / "build" / "install" / "bin"
TIMEOUT_SECONDS = 900

CODE_BLOCK = re.compile(r"^```(\w*)[^\n]*\n(.*?)^```\s*$", re.MULTILINE | re.DOTALL)
PYTHON_HEREDOC = "KAMEL_TUTORIAL_PYTHON"


def tutorial_readmes() -> list[Path]:
    return sorted(TUTORIAL_DIR.glob("*/README.md"))


def build_script(readme: Path) -> str:
    lines = ["set -euo pipefail", 'cd "$KAMEL"']
    for language, body in CODE_BLOCK.findall(readme.read_text(encoding="utf-8")):
        if language == "sh":
            lines.append(body.rstrip("\n"))
        elif language == "python":
            assert PYTHON_HEREDOC not in body
            lines.append(f"python - <<'{PYTHON_HEREDOC}'")
            lines.append(body.rstrip("\n"))
            lines.append(PYTHON_HEREDOC)
    return "\n".join(lines) + "\n"


def missing_prerequisites() -> list[str]:
    missing = []
    for pattern in ("KIM.x", "ql-balance.x", "KiLCA_Normal_*"):
        if not list(BIN_DIR.glob(pattern)):
            missing.append(f"build/install/bin/{pattern}")
    for module in ("kim", "h5py", "matplotlib"):
        if importlib.util.find_spec(module) is None:
            missing.append(f"Python module {module}")
    return missing


def tutorial_environment(home: Path, shim_bin: Path) -> dict[str, str]:
    """Run `python` and the `kim` command from the interpreter running pytest."""
    environment = dict(os.environ)
    python_bin = str(Path(sys.executable).parent)
    environment.update(
        HOME=str(home),
        KAMEL=str(REPO_ROOT),
        MPLBACKEND="Agg",
        MPLCONFIGDIR=str(home / ".matplotlib"),
        PATH=os.pathsep.join([str(shim_bin), python_bin, environment.get("PATH", "")]),
        PYTHONPATH=os.pathsep.join(
            filter(None, [str(REPO_ROOT / "python"), environment.get("PYTHONPATH")])
        ),
    )
    return environment


def test_every_tutorial_has_runnable_blocks():
    readmes = tutorial_readmes()
    assert readmes, "no tutorials found"
    for readme in readmes:
        assert "```sh" in readme.read_text(encoding="utf-8"), f"{readme} has no sh blocks"


@pytest.mark.parametrize("readme", tutorial_readmes(), ids=lambda path: path.parent.name)
def test_tutorial_runs(readme: Path, tmp_path: Path):
    missing = missing_prerequisites()
    if missing:
        message = "tutorial prerequisites missing: " + ", ".join(missing)
        if os.environ.get("KAMEL_REQUIRE_TUTORIALS") == "1":
            pytest.fail(message)
        pytest.skip(message)

    script = tmp_path / "tutorial.sh"
    script.write_text(build_script(readme), encoding="utf-8")
    home = tmp_path / "home"
    home.mkdir()
    shim_bin = tmp_path / "bin"
    shim_bin.mkdir()
    python_shim = shim_bin / "python"
    python_shim.write_text(f'#!/bin/sh\nexec "{sys.executable}" "$@"\n', encoding="utf-8")
    python_shim.chmod(0o755)
    result = subprocess.run(
        ["bash", str(script)],
        cwd=REPO_ROOT,
        env=tutorial_environment(home, shim_bin),
        capture_output=True,
        text=True,
        timeout=TIMEOUT_SECONDS,
    )
    assert result.returncode == 0, (
        f"{readme.relative_to(REPO_ROOT)} failed (exit {result.returncode})\n"
        f"--- stdout ---\n{result.stdout[-4000:]}\n--- stderr ---\n{result.stderr[-4000:]}"
    )
