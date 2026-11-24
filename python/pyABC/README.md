# pyABC — analytical bifurcation criterion (incomplete)

Python helpers intended to prepare the inputs for the analytical bifurcation
criterion. The package is unfinished: `pyABC.calc_ABC` and
`pyABC.set_Da_estimation` are empty stubs, so no criterion is calculated.

| Module | Content |
| --- | --- |
| `pyABC.py` | `pyABC` class: creates `<shot>_<time>_a2bc_inp.hdf5` in the run directory and stores shot, time slice, mode numbers, coil file and profiles |
| `da_estimator.py` | `da_estimator`: estimate of the anomalous diffusion coefficient from profiles and flux data, with text and HDF5 export |
| `field_divB0.py` | `field_divB0`: writes `field_divB0.inp` (see [`PreProc/fourier`](../../PreProc/fourier/README.md)) |
| `fouriermodes.py` | `fouriermodes()`: runs `fouriermodes.x` and stores its output in a flux-data directory |
| `GPEC_interface.py` | `GPEC_interface`: writes GPEC namelists and runs GPEC |
| `AUG_coil.py` | `AUG_coil`: reads ASDEX Upgrade coil-current files |
| `smooth2level.py` | `smooth2level()`: smoothing helper |

Known path problems: `field_divB0.py` imports `InpOut` via
`../kilca_interface/` (the package is `python/KiLCA_interface/`), and
`GPEC_interface.py` expects templates under `../matlab/balance/`, which is not
part of KAMEL, and a GPEC installation at a fixed site path.
