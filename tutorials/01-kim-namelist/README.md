# Tutorial 1: Run KIM from a namelist

This tutorial runs the KIM executable `KIM.x` directly, without the Python
interface. It computes the electrostatic plasma response to a radial magnetic
perturbation with poloidal mode number m = -6 and toroidal mode number n = 2
in a small deuterium plasma. The run takes about one second.

The case is a demonstration: the profiles have seven points and the grid has
32 points. It shows the workflow, not a converged physics result.

Before you start, build KAMEL and set up Python as described in
[the tutorial overview](../README.md).

## 1. Create a run directory

KIM reads a namelist file and a directory of profile files. Copy both into a
new working directory:

```sh
mkdir -p "$HOME/kamel-tutorials/kim-namelist"
cd "$HOME/kamel-tutorials/kim-namelist"
cp "$KAMEL/tutorials/01-kim-namelist/KIM_config.nml" .
cp -r "$KAMEL/tutorials/data/kim-profiles" profiles
ls profiles
```

The profile directory contains four two-column text files, each with the
effective radius in cm in the first column:

| File | Quantity | Unit |
| --- | --- | --- |
| `n.dat` | Electron density | 1/cm³ |
| `Te.dat` | Electron temperature | eV |
| `Ti.dat` | Ion temperature | eV |
| `q.dat` | Safety factor | 1 |

All KIM inputs use CGS units. A density of 10¹³ 1/cm³ corresponds to
10¹⁹ 1/m³.

## 2. Read the namelist

Open `KIM_config.nml`. The groups that matter most for a first run are:

- `&kim_config`: `type_of_run = 'electrostatic'` selects the global
  electrostatic solver and `collision_model = 'FokkerPlanck'` the collision
  model.
- `&kim_io`: `profile_location` points to the profile directory and must end
  with `/`. `hdf5_output = .true.` writes the results to an HDF5 file.
- `&kim_setup`: the toroidal field `btor` in G, the major radius `r0` in cm,
  and the mode numbers `m_mode` and `n_mode`.
- `&kim_grid`: the radial domain from `r_min` to `r_plas` in cm and the number
  of grid points.

The [namelist reference](../../KIM/nmls/README.md) describes every option.
KIM stops with an error if the file contains an option it does not know.

## 3. Run KIM

Pass the namelist path as the only argument:

```sh
cd "$HOME/kamel-tutorials/kim-namelist"
"$KAMEL/build/install/bin/KIM.x" KIM_config.nml > kim.log
tail -n 3 kim.log
```

KIM prints `WARNING: Er.dat not found` and computes the equilibrium radial
electric field from radial force balance instead. It writes the result to
`profiles/Er.dat`. To use your own radial electric field, put an `Er.dat`
file in the profile directory before the run.

The results are written to `out/m<m>_n<n>/`:

```sh
cd "$HOME/kamel-tutorials/kim-namelist"
ls out/m-6_n2
```

## 4. Look at the result

`KIM_out.h5` contains the solved fields in the group `/fields`, the radial
grids in `/grid`, the background plasma in `/backs` and a copy of the input
settings. Every dataset carries `unit` and `comment` attributes. Complex
fields are stored with `real` and `imag` components. The fields are given on
the grid `/grid/xl_xb`.

This script plots the electrostatic potential and the parallel current
density. The radial magnetic perturbation `fields/br` is the prescribed
constant drive in this run (`type_br_field = 12`):

```python
import h5py
import matplotlib.pyplot as plt
import numpy as np


def read_complex(dataset):
    data = dataset[()]
    return data["real"] + 1j * data["imag"]


with h5py.File("out/m-6_n2/KIM_out.h5", "r") as f:
    r = f["grid/xl_xb"][()]
    phi = read_complex(f["fields/Phi_m"])
    jpar = read_complex(f["fields/jpar"])
    print("Phi unit:", f["fields/Phi_m"].attrs["unit"])

fig, (ax_phi, ax_j) = plt.subplots(2, 1, sharex=True)
ax_phi.plot(r, np.abs(phi))
ax_phi.set_ylabel("|Phi| [statV]")
ax_j.plot(r, np.abs(jpar))
ax_j.set_ylabel("|j_par| [statA/cm²]")
ax_j.set_xlabel("r [cm]")
fig.savefig("kim_fields.png", dpi=120)
print("saved kim_fields.png")
```

To list everything in the file, run `h5ls -r out/m-6_n2/KIM_out.h5` or
iterate over the file with `h5py`.

## Next steps

- Change `m_mode` and `n_mode` and rerun. The resonant surface lies where
  q = -m/n; with these profiles it must stay between `r_min` and `r_plas`.
- Increase `l_space_dim` and `rg_space_dim` and compare the fields.
- [Tutorial 2](../02-kim-python/README.md) runs the same case through the
  Python interface, which validates inputs and records every run.
