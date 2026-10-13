import sys
import textwrap
from pathlib import Path

import f90nml
import h5py
import numpy as np
import pytest
from neo2_for_Er.local_runner import (
    LocalNeo2Plan,
    Neo2LocalError,
    collect_staged_surfaces,
    run_staged_surfaces,
    stage_surfaces,
)


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


@pytest.fixture
def batch_case(tmp_path):
    rows = [
        [0.3, 30, 150, 1, 100, 2e13, -0.5],
        [0.1, 10, 160, 2, 200, 3e13, -0.7],
        [0.3, 31, 170, 3, 300, 4e13, -0.9],
    ]
    np.savetxt(tmp_path / "surfaces.dat", rows)
    template = tmp_path / "template"
    template.mkdir()
    (template / "neo2.in").write_text(
        "&settings\nboozer_s=<boozer_s>,rbeg=<rbeg>,zbeg=<zbeg>\n/\n"
        "&collision\nconl_over_mfp=<conl_over_mfp>\n/\n"
    )
    return tmp_path, template


def write_batch_output(root, job, first, last):
    rows = np.loadtxt(root / "k_surfaces.dat", ndmin=2)[first - 1 : last]
    with h5py.File(job / "k_coefficients.h5", "w") as output:
        output["surface_first"], output["surface_last"] = first, last
        output["k_coefficient_only"] = 1
        for name, column in zip(("boozer_s", "rbeg", "zbeg", "conl_over_mfp"), rows.T):
            output[name] = column
        output["k_cof"] = rows[:, 0] / 2
        output["qcurr_g"] = np.ones(len(rows))
        output["qcurr_t"] = 2.5 - rows[:, 0] / 2
        output["integral_converged"] = np.zeros(len(rows), dtype=int)


@pytest.mark.parametrize("surfaces_per_run, chunks", [(0, 1), (2, 2)])
def test_native_chunks_execute_and_preserve_row_mapping(batch_case, surfaces_per_run, chunks):
    root, template = batch_case
    jobs = stage_surfaces(root, template, surfaces_per_run)
    executable = root / "fake_neo2"
    executable.write_text(f"#!{sys.executable}\n" + textwrap.dedent("""\
        import f90nml, h5py, numpy as np
        nml = f90nml.read('neo2.in')['settings']
        first, last = nml['surface_first'], nml['surface_last']
        rows = np.loadtxt(nml['k_surface_file'], ndmin=2)[first-1:last]
        with h5py.File(nml['k_output_file'], 'w') as output:
            output['surface_first'], output['surface_last'] = first, last
            output['k_coefficient_only'] = 1
            for key, values in zip(('boozer_s','rbeg','zbeg','conl_over_mfp'), rows.T):
                output[key] = values
            output['k_cof'] = rows[:,0] / 2
            output['qcurr_g'] = np.ones(len(rows))
            output['qcurr_t'] = 2.5 - rows[:,0] / 2
            output['integral_converged'] = np.zeros(len(rows), dtype=int)
        """))
    executable.chmod(0o755)
    result = run_staged_surfaces(root, LocalNeo2Plan(executable, max_workers=1))
    np.testing.assert_allclose(result, [[10, 0.05], [30, 0.15], [31, 0.15]])
    assert len(jobs) == chunks
    assert (root / "k_surfaces.dat").is_file()
    with pytest.raises(Neo2LocalError, match="overwrite"):
        stage_surfaces(root, template, surfaces_per_run)


@pytest.mark.parametrize(
    "fault", ["missing", "nan", "identity", "shape", "range", "overlap", "coverage"]
)
def test_batch_collection_rejects_incomplete_or_wrong_results(batch_case, fault):
    root, template = batch_case
    jobs = stage_surfaces(root, template, 2)
    write_batch_output(root, jobs[0], 1, 2)
    write_batch_output(root, jobs[1], 3, 3)
    if fault == "missing":
        (jobs[1] / "k_coefficients.h5").unlink()
    elif fault == "overlap":
        write_batch_output(root, jobs[1], 2, 3)
    elif fault == "coverage":
        (root / "jobs_list.txt").write_text(jobs[0].name + "\n")
    else:
        with h5py.File(jobs[1] / "k_coefficients.h5", "r+") as output:
            if fault == "nan":
                output["k_cof"][0] = np.nan
            elif fault == "identity":
                output["rbeg"][0] += 1
            elif fault == "shape":
                del output["k_cof"]
                output["k_cof"] = []
            elif fault == "range":
                output["surface_last"][...] = 4
    with pytest.raises(Neo2LocalError):
        collect_staged_surfaces(root)


def test_profile_collector_reads_native_chunks(batch_case):
    from profile_processor import Profile_Processor

    root, template = batch_case
    jobs = stage_surfaces(root, template, 0)
    write_batch_output(root, jobs[0], 1, 3)
    processor = Profile_Processor.__new__(Profile_Processor)
    processor.collect_k_profile(str(root), str(root))
    np.testing.assert_allclose(np.loadtxt(root / "k.dat"), [[10, 0.05], [30, 0.15], [31, 0.15]])
    np.testing.assert_allclose(processor.k_prof, [0.05, 0.15, 0.15])


@pytest.mark.parametrize("value", [-1, 1.5])
def test_invalid_chunk_size_does_not_stage(batch_case, value):
    root, template = batch_case
    with pytest.raises(Neo2LocalError, match="nonnegative integer"):
        stage_surfaces(root, template, value)
    assert not (root / "k_surfaces.dat").exists()
