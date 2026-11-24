# QL-Balance

QL-Balance solves one-dimensional radial transport equations for the plasma
profiles with quasilinear transport coefficients driven by resonant magnetic
perturbations. The perturbation fields come from a linear plasma-response
code: KiLCA (default) or KIM. Relevant publications are Heyn et al., Nucl.
Fusion 54 (2014) and Markl et al., Nucl. Fusion 63 (2023).

## Build

From the KAMEL repository root:

```sh
make QL-Balance   # builds KiLCA and KIM as dependencies
```

The executable is `build/install/bin/ql-balance.x`.

## Run

QL-Balance reads `balance_conf.nml` from the current working directory:

```sh
cd /path/to/run-directory
/path/to/KAMEL/build/install/bin/ql-balance.x
```

The template [`namelists/balance_conf.nml`](namelists/balance_conf.nml)
contains the `&BALANCENML` group (run control, grid, paths, response-code
selection) and the `&NEORT` group used by NTV runs. Stellarator runs
additionally read `stell_conf.nml` (template in
[`namelists/stell_conf.nml`](namelists/stell_conf.nml)).

### Run types (`type_of_run`)

| Value | Purpose |
| --- | --- |
| `SingleStep` | One evaluation of the quasilinear transport coefficients |
| `TimeEvolution` | Time evolution of the profiles under quasilinear transport |
| `TimeEvolutionNTV` | Time evolution including neoclassical toroidal viscosity from NEO-RT (`&NEORT`) |
| `ParameterScan` | Scan over factor sets for the n, Te, Ti and V_z profiles (`fac_n`, `fac_Te`, `fac_Ti`, `fac_vz`) |
| `TimeEvolutionStellarator` | Stellarator time evolution (`stell_conf.nml`) |

### Response code (`wave_code`)

- `'KiLCA'` (default): the run directory must also contain the KiLCA input
  files (see [KiLCA/README.md](../KiLCA/README.md)); `flre_path` and
  `vac_path` point to the KiLCA plasma and vacuum projects.
- `'KIM'`: KIM is called in-process with the namelist given by
  `kim_config_path`. `kim_run_type` selects the KIM run type
  (`electromagnetic` by default, or `electrostatic_periodic`), and the
  `kim_*` options set the transport models and mode list. See the
  [KIM namelist reference](../KIM/nmls/README.md).

Input profiles, output files and units are described in
[Data.md](../Data.md).

## Test

Unit and regression tests live in [`src/test/`](src/test/) and are run
with `make test` from the repository root. QL-Balance is also covered by
the golden-record suite in [`test/golden/`](../test/golden/README.md).
