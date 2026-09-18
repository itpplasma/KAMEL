# MARS-F reader fixture

This synthetic fixture exercises the read-only A08 MARS-F profile boundary.
It uses the legacy four-file layout:

- `PROFDEN.IN`: density in `1/m^3`
- `PROFTE.IN`: electron temperature in `eV`
- `PROFTI.IN`: ion temperature in `eV`
- `PROFROT.IN`: toroidal velocity in `m/s`

The coordinate is `sqrt_psiN` with unit `1`. These values are test material,
not an approved physical case. The reader preserves the source units and
coordinate grids; it does not interpolate or convert them.
