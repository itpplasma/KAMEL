"""Run the production FLRE/vacuum pipeline with complex quantity output enabled."""

import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def rows(path, columns, shared_boundaries=False):
    values = [
        [float(token) for token in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("%")
    ]
    assert len(values) > 1, f"Empty or truncated output: {path}"
    assert all(len(row) == columns for row in values), f"Wrong column count: {path}"
    assert all(math.isfinite(x) for row in values for x in row), f"Nonfinite output: {path}"
    assert all(
        a[0] <= b[0] if shared_boundaries else a[0] < b[0] for a, b in zip(values, values[1:])
    ), f"Invalid grid: {path}"
    return values


def main(executable, case):
    with tempfile.TemporaryDirectory(prefix="kilca_full_output_", dir=Path.cwd()) as tmp:
        run = Path(tmp) / "flre_m6n2"
        shutil.copytree(case, run)
        settings = run / "output.in"
        settings.write_text(
            "\n".join(
                "2" + line[1:] if line.startswith("1 ") else line
                for line in settings.read_text().splitlines()
            )
            + "\n"
        )
        result = subprocess.run(
            [executable], cwd=run, capture_output=True, text=True, timeout=120, check=False
        )
        assert result.returncode == 0, result.stdout + result.stderr
        output = run / "linear-data" / "m_6_n_2_flab_[1,0]"
        expected = {
            f"zone_0_current_dens_{axis}_{kind}_{species}{frame}.dat"
            for axes, frame in [("rsp", ""), ("rsptz", "_lab")]
            for axis in axes
            for kind in (0, 1)
            for species in "iet"
        }
        assert {p.name for p in output.glob("*current_dens*.dat")} == expected
        expected.update(f"zone_0_density_{species}.dat" for species in "iet")
        for name in sorted(expected):
            data = rows(output / name, 3)
            assert any(abs(row[2]) > 0 for row in data), f"Lost imaginary component: {name}"
        for component in "rtz":
            for species in "iet":
                for quantity in ("dens", "int"):
                    rows(output / f"zone_0_torque_{quantity}_{component}_{species}.dat", 2)
        fields = rows(output / "EB.dat", 13, shared_boundaries=True)
        assert fields[0][0] == 3.0 and fields[-1][0] == 80.0, "Incomplete zone coverage"
        print("PASS: production FLRE/vacuum solve, complex output, and Lorentz torque")


if __name__ == "__main__":
    main(str(Path(sys.argv[1]).resolve()), Path(sys.argv[2]).resolve())
