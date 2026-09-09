#!/usr/bin/env python3
"""Controllable stand-in for KIM.x used by orchestration tests."""

from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path


def namelist_value(text: str, name: str) -> str:
    match = re.search(rf"^\s*{name}\s*=\s*['\"]([^'\"]+)['\"]", text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"missing {name} in namelist")
    return match.group(1)


def integer_value(text: str, name: str) -> int:
    match = re.search(rf"^\s*{name}\s*=\s*(-?\d+)", text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"missing {name} in namelist")
    return int(match.group(1))


def main() -> int:
    mode = os.environ.get("FAKE_KIM_MODE", "success")
    invocation = {"argv": sys.argv, "cwd": str(Path.cwd())}
    Path("invocation.json").write_text(json.dumps(invocation), encoding="utf-8")
    print("fake KIM stdout", flush=True)
    if mode == "stderr":
        print("fake KIM diagnostic", file=sys.stderr, flush=True)
    if mode == "nonzero":
        print("fake KIM failure", file=sys.stderr, flush=True)
        return 7
    if mode == "timeout":
        print("fake KIM waiting", file=sys.stderr, flush=True)
        time.sleep(10)
        return 0

    import h5py
    import numpy as np

    namelist = Path(sys.argv[1]).read_text(encoding="utf-8")
    output_root = Path(namelist_value(namelist, "output_path"))
    output_file = namelist_value(namelist, "h5_out_file")
    mode_name = f"m{integer_value(namelist, 'm_mode')}_n{integer_value(namelist, 'n_mode')}"
    mode_dir = output_root / mode_name
    mode_dir.mkdir(parents=True, exist_ok=True)
    if mode == "missing_output":
        return 0

    with h5py.File(mode_dir / output_file, "w") as handle:
        potential_name = "fields/Phi_m" if output_file == "out_ES.h5" else "fields/Phi"
        handle.create_dataset(potential_name, data=np.ones(3))
        if mode != "partial_output":
            handle.create_dataset("fields/jpar", data=np.ones(3))
            handle.create_dataset("backs/e/r", data=np.arange(3.0))
            handle.create_dataset("setup/periodic_scale/dx_asis", data=1.0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
