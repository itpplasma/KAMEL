import os
import subprocess
import sys
import textwrap


def test_neo2_import_and_staging_use_python38_dataclass_and_zip_apis():
    script = textwrap.dedent(
        """
        import builtins
        import dataclasses
        import tempfile
        from pathlib import Path

        import h5py
        import numpy as np

        original_dataclass = dataclasses.dataclass
        def python38_dataclass(*args, **kwargs):
            if "slots" in kwargs:
                raise TypeError("Python 3.8 dataclass does not accept slots")
            return original_dataclass(*args, **kwargs)
        dataclasses.dataclass = python38_dataclass

        original_zip = builtins.zip
        def python38_zip(*iterables, **kwargs):
            if kwargs:
                raise TypeError("Python 3.8 zip does not accept keyword arguments")
            return original_zip(*iterables)
        builtins.zip = python38_zip

        from neo2_for_Er import Neo2CondorPlan, stage_neo2_condor_jobs
        from neo2_for_Er.local_runner import LocalNeo2Plan, run_staged_surfaces, stage_surfaces

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            np.savetxt(root / "surfaces.dat", [[0.25, 42, 150, 0, 1000, 1e13, -0.5]])
            template = root / "template"
            template.mkdir()
            (template / "neo2.in").write_text(
                "s=<boozer_s> k=<conl_over_mfp> R=<rbeg> Z=<zbeg>\\n"
            )
            jobs = stage_surfaces(root, template)
            executable = root / "neo2.x"
            executable.write_text(
                "#!" + __import__("sys").executable + "\\n"
                "import h5py\\n"
                "with h5py.File('neo2_config.h5','w') as f:f.create_dataset('settings/boozer_s',data=0.25)\\n"
                "with h5py.File('fulltransp.h5','w') as f:f.create_dataset('k_cof',data=0.5)\\n"
            )
            executable.chmod(0o755)
            plan = Neo2CondorPlan(
                executable=executable,
                python_executable=Path(__import__("sys").executable),
                shared_filesystem_prefixes=(root,),
                parallel_runtime={"omp": False, "mpi": False},
            )
            stage_neo2_condor_jobs(root, plan=plan)
            profile = run_staged_surfaces(root, LocalNeo2Plan(executable, max_workers=1))
            assert profile.shape == (1, 2)
        """
    )
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, ["python", environment.get("PYTHONPATH")])
    )

    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, env=environment
    )

    assert result.returncode == 0, result.stderr
