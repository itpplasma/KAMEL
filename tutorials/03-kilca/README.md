# Tutorial 3: Run KiLCA

KiLCA computes the linear response of a cylindrical plasma to an external
antenna current. The radial domain is split into zones, and each zone uses its
own plasma model. This tutorial runs a single mode (m = 6, n = 2) with a
finite-Larmor-radius (FLR) plasma zone, two vacuum zones and an antenna. The
run takes less than a second.

Before you start, build KAMEL and set up Python as described in
[the tutorial overview](../README.md).

## 1. Create a run directory

KiLCA reads all its inputs from one directory, called the project directory.
**The path of that directory must contain `flre` or `vacuum`**: KiLCA uses the
path to decide whether it computes a plasma or a vacuum project and stops with
`unknown project type in path` otherwise. A plasma run therefore needs a
directory such as `.../flre`.

```sh
mkdir -p "$HOME/kamel-tutorials/kilca"
cd "$HOME/kamel-tutorials/kilca"
rm -rf flre
cp -r "$KAMEL/tutorials/03-kilca/flre" flre
ls flre flre/profiles
```

## 2. Read the inputs

Each `.in` file holds one value per line, followed by a `#` comment that
explains it. The files are:

| File | Content in this example |
| --- | --- |
| `background.in` | Major radius 165 cm, plasma radius 67 cm, toroidal field -1.76·10⁴ G, profile directory `./profiles/`, ion mass 2 (deuterium) |
| `antenna.in` | Antenna at r = 70 cm with a current of 5.1·10¹² statA, frequency, number of modes |
| `modes.in` | The mode list, here `(6, 2)` |
| `zone_1_flre.in` | FLR plasma from r = 3 cm (cylinder centre) to 67 cm |
| `zone_2_vacuum.in` | Vacuum from 67 cm to the antenna at 70 cm |
| `zone_3_vacuum.in` | Vacuum from the antenna to an ideal wall at 80 cm |
| `output.in` | Which quantities are computed and saved |
| `eigmode.in` | Settings for eigenmode searches (not used here) |

The `profiles/` directory contains two-column text files with the radius in
cm and the profile value in CGS units: density `n.dat`, temperatures `Te.dat`
and `Ti.dat` in eV, safety factor `q.dat`, radial electric field `Er.dat`, and
the toroidal and poloidal rotation velocities `Vz.dat` and `Vth.dat`.
[Data.md](../../Data.md) describes all formats.

## 3. Run KiLCA

The executable name contains the version and build options. KiLCA uses the
current directory as the project directory, or the directory given as its
first argument:

```sh
cd "$HOME/kamel-tutorials/kilca/flre"
KILCA=$(ls "$KAMEL"/build/install/bin/KiLCA_Normal_*)
"$KILCA" > kilca.log
ls background-data linear-data
```

`background-data/` contains the equilibrium quantities computed from the
profiles. `linear-data/` holds one directory per mode and frequency, here
`m_6_n_2_flab_[1,0]` for m = 6, n = 2 and the laboratory frequency (1, 0)
from `antenna.in`.

## 4. Look at the result

`EB.dat` in the mode directory contains the perturbation fields. Its 13
columns are the radius in cm, followed by the real and imaginary parts of
`Er`, `Etheta`, `Ez` (statV/cm) and `Br`, `Btheta`, `Bz` (G).

```python
import matplotlib.pyplot as plt
import numpy as np

eb = np.loadtxt("linear-data/m_6_n_2_flab_[1,0]/EB.dat")
r = eb[:, 0]
br = eb[:, 7] + 1j * eb[:, 8]

fig, ax = plt.subplots()
ax.plot(r, np.abs(br))
ax.axvline(67.0, linestyle="--", color="gray", label="plasma edge")
ax.axvline(70.0, linestyle=":", color="gray", label="antenna")
ax.set_xlabel("r [cm]")
ax.set_ylabel("|Br| [G]")
ax.set_yscale("log")
ax.legend()
fig.savefig("kilca_br.png", dpi=120)
print("saved kilca_br.png")
```

## Using Python to prepare KiLCA runs

The package [`python/KiLCA_interface`](../../python/KiLCA_interface/)
reads and writes all KiLCA input files, so parameters can be changed from
Python instead of editing the files. Templates for every file are in
[`python/KiLCA_interface/blueprints/`](../../python/KiLCA_interface/blueprints/).
[Tutorial 4](../04-ql-balance/README.md) uses this package indirectly: the
QL-Balance interface prepares and runs KiLCA before the transport run.

## Next steps

- Change the antenna current in `antenna.in` and rerun. The response is
  linear, so the fields scale with the current.
- Add a mode to `modes.in` and raise the number of modes in `antenna.in`.
