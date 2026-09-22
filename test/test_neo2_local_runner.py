import sys
from pathlib import Path

import h5py
import numpy as np
from neo2_for_Er.local_runner import LocalNeo2Plan, run_staged_surfaces, stage_surfaces


def test_stages_and_runs_explicit_surfaces_locally(tmp_path: Path):
    np.savetxt(tmp_path / "surfaces.dat", [[0.2, 12, 150, 0, 100, 2e13, -0.5]])
    template = tmp_path / "template"
    template.mkdir()
    (template / "neo2.in").write_text("s=<boozer_s> k=<conl_over_mfp> R=<rbeg> Z=<zbeg>\n")
    jobs = stage_surfaces(tmp_path, template)
    executable = tmp_path / "neo2"
    executable.write_text(
        f"#!{sys.executable}\nimport h5py\n"
        "from pathlib import Path\n"
        "s=float(Path('neo2.in').read_text().split()[0].split('=')[1].replace('d','e'))\n"
        "with h5py.File('neo2_config.h5','w') as f:f.create_dataset('settings/boozer_s',data=s)\n"
        "with h5py.File('fulltransp.h5','w') as f:f.create_dataset('k_cof',data=s/2)\n"
    )
    executable.chmod(0o755)
    output = run_staged_surfaces(tmp_path, LocalNeo2Plan(executable, max_workers=1))
    assert jobs[0].is_dir()
    np.testing.assert_allclose(output[:, 1], [0.1])
