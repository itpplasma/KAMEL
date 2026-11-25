# KiLCA — Kinetic Linear Cylindrical Approximation

KiLCA computes the linear response of a cylindrical plasma to external
magnetic perturbations, using a finite-Larmor-radius (FLR) expansion of the
kinetic conductivity. The radial domain is split into zones; each zone is
solved with its own plasma model (FLR, ideal MHD, homogeneous medium, or
vacuum), and the zone solutions are stitched at the boundaries. The model is
described in Heyn et al., Nucl. Fusion 46 (2006), Heyn et al., Nucl. Fusion
48 (2008) and Ivanov et al., Phys. Plasmas 18 (2011); see
[Publications](../README.md#publications) for the full references.

## Build

From the KAMEL repository root:

```sh
make KiLCA        # or `make all` for every component
```

The build installs into `build/install/bin/`:

| Executable | Purpose |
| --- | --- |
| `KiLCA_Normal_<signature>` | Linear response run (antenna drive or eigenmode search) |
| `KiLCA_EigParam_<signature>` | Eigenmode parameter scan; stages one run directory per scan value |
| `post_processor_<arch>` | Post-processing of a finished run (reads `post_proc.in`) |

`<signature>` encodes the version and compile-time options, e.g.
`V_2.4.2_..._Release_64bit`. Run `ls build/install/bin/KiLCA_*` to see the
exact name.

## Run

The executable takes the run (project) directory as its first argument, or
uses the current working directory if none is given:

```sh
build/install/bin/KiLCA_Normal_<signature> /path/to/flre
```

The path of the project directory must contain `flre` (plasma) or `vacuum`
(vacuum): KiLCA selects the project type from the path and stops with
`unknown project type in path` otherwise.

A run directory contains:

| File | Content |
| --- | --- |
| `background.in` | Background settings, including the profile directory |
| `antenna.in` | Antenna settings; selects antenna drive or eigenmode search |
| `modes.in` | List of (m, n) modes |
| `output.in` | Output selection |
| `eigmode.in` | Eigenmode search settings (eigenmode runs) |
| `*zone_<k>*.in` | One file per radial zone, `k = 1, 2, ...`; the file content selects the zone model |

Templates for all files are in
[`python/KiLCA_interface/blueprints/`](../python/KiLCA_interface/blueprints/).
The Python package [`python/KiLCA_interface`](../python/KiLCA_interface/)
reads and writes these files and is the recommended way to prepare runs.
Input profiles, file formats, and outputs (including the per-zone
`zone_<k>_EB.dat` field files) are described in [Data.md](../Data.md).

## Test

KiLCA unit and regression tests are registered with CTest and run with
`make test` from the repository root. KiLCA is also covered by the
golden-record suite in [`test/golden/`](../test/golden/README.md).

The solver was ported from C++ to Fortran in 2026; the design and
verification records are in
[`docs/plans/2026-02-16-kilca-cpp-to-fortran-design.md`](../docs/plans/2026-02-16-kilca-cpp-to-fortran-design.md)
and [`docs/verification/pr161-fortran-migration.md`](../docs/verification/pr161-fortran-migration.md).
