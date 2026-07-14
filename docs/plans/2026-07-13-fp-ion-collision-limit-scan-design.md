# FP Ion-Collision Limit Scan Design

The scan varies only the Fokker--Planck ion collision frequency.  A positive
configuration multiplier `ion_fp_collision_scale`, defaulting to one, is
applied to every ion species after KIM computes its physical collision
frequency and before KIM constructs `z0`, cell-centred backgrounds, `x1`,
`x2`, or the I-functions.  Electron collision frequencies are never scaled.

The benchmark evaluates both `(m,n)=(6,2)` and `(7,2)` at multipliers
`1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001`.  The scale-one FP result is the collisional
reference, and the existing causal analytical collisionless-ion result is the
zero-collision comparison.  Each scale has its own immutable configuration,
output directory, log, and metadata.  The analysis reports finiteness,
relative complex L2 distance to both endpoints, peak amplitude, resonance
amplitude and phase, and writes JSON, CSV, and a PDF plot.

This is a numerical limiting scan, not an assumption that the FP I-function
implementation is well conditioned at zero collision frequency.  KIM rejects
non-positive multipliers because the FP expressions explicitly divide by
`nu`.  If a low-scale run becomes non-finite or fails, the scan records that
as the observed numerical limit rather than substituting the analytical
collisionless kernel.
