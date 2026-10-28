# Species density output

The global Fokker--Planck electrostatic and electromagnetic paths expose
electron `delta_n_e` and per-ion `delta_n_i` number-density amplitudes on the
field grid (`cm^-3`). The in-process result also supplies the signed charge
numbers in electron-then-ion order. HDF5 names are `/fields/delta_n_e` and
`/fields/delta_n_<species name>` (with a species index appended for multiple
ions); `/fields/rho` is their charge-weighted sum.

The code applies the assembled charge kernels to the solved potential and
radial magnetic field, then solves the raw unweighted `dr` mass system.
Field boundary conditions never constrain this physical moment projection.
The same projection is used for parallel current. Kernel charge responses
already include the local Debye term; no extra Boltzmann term is added.

Other run types leave these new API buffers unallocated. Their existing
approximate charge diagnostics are not validated species number densities.
These outputs do not supply a radial particle flux or a toroidal closure.

The uniform finite-interval paths require the closed-grid repair. Adaptive
field and background endpoint agreement has not been validated here. The
Boltzmann oracle uses compact potential support with a buffer from the finite
guiding-centre boundaries, and refines both background and kernel quadrature.
