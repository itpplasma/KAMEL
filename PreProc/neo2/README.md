# NEO-2 run templates

Templates and helper scripts for running [NEO-2](https://github.com/itpplasma/neo-2)
on a set of flux surfaces. NEO-2 computes the neoclassical transport coefficients
that are needed to calculate the poloidal rotation velocity and hence the radial
electric field. NEO-2 itself is not part of KAMEL.

The Python workflow in [`python/neo2_for_Er/`](../../python/neo2_for_Er/) uses
this directory as its template location.

## Files

| File | Purpose |
| --- | --- |
| `TEMPLATE_DIR/neo2.in` | NEO-2 input with placeholders `<conl_over_mfp>`, `<boozer_s>`, `<rbeg>`, `<zbeg>` |
| `TEMPLATE_DIR/condor.submit` | HTCondor submit file for one `neo_2.x` run |
| `create_surf_realspace.py` | Creates one run directory per flux surface |
| `run_condor.sh` | Submits every run directory listed in `jobs_list.txt` to HTCondor |
| `remote_run.conf` | Hosts for `neo2_for_Er.run_remote_neo2`: `host cpus memory_MB` per line, `#` comments a line out |

## Usage

1. Provide `surfaces.dat` in the working directory with seven columns:
   `s`, `r_eff`, `rbeg`, `zbeg`, `T`, `N`, `kappa`.
2. Run `create_surf_realspace.py`. It copies `TEMPLATE_DIR` to one directory per
   surface (named after `s`), fills the placeholders in `neo2.in` from
   `surfaces.dat`, and writes the directory names to `jobs_list.txt`.
   Existing directories with the same names are deleted first.
3. Run `run_condor.sh` (or `run_condor.sh -d` for a dry run) to submit the jobs.

`condor.submit` and `remote_run.conf` contain site-specific settings (CPU and
memory requests, excluded and listed machines); adapt them to your cluster.
