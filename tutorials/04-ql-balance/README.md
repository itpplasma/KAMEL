# Tutorial 4: Run QL-Balance with KiLCA

QL-Balance evolves the radial plasma profiles under quasilinear transport
driven by a resonant magnetic perturbation. It needs the perturbation fields
from a response code. This tutorial uses KiLCA, prepared through the KAMELpy
Python package, and runs a short time evolution for the mode m = 6, n = 2.
Preparation takes a few seconds and the QL-Balance run about 15 seconds.

The case uses synthetic parabolic profiles and is meant to show the workflow.

Before you start, build KAMEL and set up Python as described in
[the tutorial overview](../README.md). Working through
[Tutorial 3](../03-kilca/README.md) first helps, since QL-Balance runs KiLCA
internally.

## 1. Prepare the run directory

A QL-Balance run directory needs:

- `balance_conf.nml`, the QL-Balance settings;
- `profiles/` with the plasma profiles;
- `flre/` and `vacuum/`, KiLCA projects for the plasma and vacuum response;
- the executable, or a link to it.

Assembling this by hand is error-prone, so the script
[`setup_run.py`](setup_run.py) does it with KAMELpy. It creates parabolic
profiles around the resonant surface q = 3, runs KiLCA for both projects,
writes `balance_conf.nml` from the template in
[`QL-Balance/namelists/`](../../QL-Balance/namelists/balance_conf.nml), and
links `ql-balance.x` as `ql-balance`:

```sh
python "$KAMEL/tutorials/04-ql-balance/setup_run.py" "$HOME/kamel-tutorials/ql-balance"
cd "$HOME/kamel-tutorials/ql-balance"
ls
```

Read `setup_run.py` to see which values are set. The settings that differ
from the template are:

| Setting | Value | Meaning |
| --- | --- | --- |
| `type_of_run` | `'TimeEvolution'` | Evolve the profiles in time |
| `wave_code` | `'KiLCA'` (template default) | Response fields from KiLCA |
| `ramp_up_mode` | `3` | Switch the full coil current on at once |
| `suppression_mode` | `.false.` | Save the full profile history |
| `save_prof_time_step` | `1` | Save the profiles at every time step |

The paths to the KiLCA projects (`flre_path`, `vac_path`) and to the input
and output files are written as absolute paths.

## 2. Run QL-Balance

QL-Balance reads `balance_conf.nml` from the current directory:

```sh
cd "$HOME/kamel-tutorials/ql-balance"
./ql-balance > ql-balance.log
tail -n 1 ql-balance.log
ls out
```

The run stops when it reaches the maximum time and writes
`out/0_0_tutorial.hdf5`.

## 3. Look at the result

The output file has one group per mode, here `f_6_2`. It contains:

- `r_res`: the radius of the resonant surface in cm;
- `KinProfiles/<step>/`: snapshots of the plasma profiles (`n`, `Te`, `Ti`,
  `Vz`, `Er`) on the radial grid `rc`;
- `LinearProfiles/<step>/`: the perturbation fields and quasilinear transport
  coefficients, such as `dqle22`, after each KiLCA evaluation;
- time traces such as `br_abs_time`, the perturbation amplitude at the
  resonance.

[Data.md](../../Data.md) describes the complete file.

This script plots a transport coefficient and the change of the density
between the first and the last saved step:

```python
import h5py
import matplotlib.pyplot as plt

with h5py.File("out/0_0_tutorial.hdf5", "r") as f:
    mode = f["f_6_2"]
    r_res = mode["r_res"][0]

    snapshots = mode["KinProfiles"]
    steps = sorted(snapshots, key=int)
    rc = snapshots[steps[0]]["rc"][()]
    density_change = snapshots[steps[-1]]["n"][()] - snapshots[steps[0]]["n"][()]

    linear = mode["LinearProfiles"]
    last = linear[sorted(linear, key=int)[-1]]
    r = last["r"][()]
    dqle22 = last["dqle22"][()]

print(f"resonant surface at r = {r_res:.2f} cm, {len(steps)} saved steps")

fig, (ax_d, ax_n) = plt.subplots(2, 1, sharex=True)
ax_d.semilogy(r, dqle22)
ax_d.set_ylabel("dqle22")
ax_n.plot(rc, density_change)
ax_n.set_ylabel("n(last) - n(first) [1/cm³]")
ax_n.set_xlabel("r [cm]")
for ax in (ax_d, ax_n):
    ax.axvline(r_res, linestyle="--", color="gray")
ax_n.set_xlim(r_res - 10, r_res + 10)
fig.savefig("ql_balance.png", dpi=120)
print("saved ql_balance.png")
```

## Without Python

`ql-balance.x` itself needs no Python, but the run directory does not come
from a single input file. Besides `balance_conf.nml` and the two KiLCA
projects, QL-Balance reads an input HDF5 file (`path2inp`, here
`INPUT_0_0_tutorial.hdf5`) that KAMELpy writes during preparation.
[Data.md](../../Data.md) lists what each of these files contains. To prepare
runs without Python, start from the directory created in step 1 and replace
its contents piece by piece. Once a directory exists, rerunning
`./ql-balance` there needs no Python at all.

## Using KIM instead of KiLCA

Set `wave_code = 'KIM'` in `balance_conf.nml` to compute the response fields
with KIM, which QL-Balance calls in-process. `kim_config_path` then points to
a KIM namelist (see [Tutorial 1](../01-kim-namelist/README.md)), and the
`kim_*` options select the run type, transport models and modes. The
[QL-Balance README](../../QL-Balance/README.md) lists these options.
