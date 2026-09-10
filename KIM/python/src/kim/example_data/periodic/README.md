# Periodic parabolic reference case

This directory contains a complete, compact forced-periodicity case for the electrostatic KIM
solver. The six profile files use effective radius (`r_eff`) in centimetres and KIM's CGS profile
units: density in `1/cm^3`, temperatures in eV, radial electric field in statV/cm, and toroidal
velocity in cm/s. The request uses deuterium with `(m, n) = (7, 2)` and therefore has its signed
`q = -m/n = -3.5` resonant surface near `r_eff = 49.762 cm`.

The `request.json` values and 65 point profiles are the reference inputs used by the opt-in real
integration test. The grid and periodic-window settings are deliberately compact demonstration
settings so that the case is practical for examples and automated checks; they are not generally
converged production defaults. New physical cases should perform their own radial, periodic-window,
and quadrature convergence checks.

Run this case from a directory containing a copy of both `request.json` and `profiles/`, for
example:

```bash
kim validate request.json
kim run request.json --executable /path/to/KIM.x --runs-dir runs
```
