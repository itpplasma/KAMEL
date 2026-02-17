# Design: Self-Consistent Br Calculation (Issue #52)

## Overview

Implement self-consistent radial magnetic field (Br) calculation in KIM by coupling
the existing electrostatic Poisson equation with the parallel Ampere equation from
Section 5 of the KIM writeup. The two equations are assembled into a single 2N x 2N
block system and solved directly.

## Architecture

### New Run Type

A new `electromagnetic_t` type extending `kim_t`, following the existing factory pattern.

- Selected via `type_of_run = "electromagnetic"` in `&kim_config`
- Requires `collision_model = "FokkerPlanck"` (all 4 kernel matrices needed)
- Registered in `from_kim_factory_get_kim` in `KIM_mod.f90`

### New Namelist Parameters (`&kim_setup`)

- `Br_boundary_re` (real, default 1.0) -- real part of Br at right boundary
- `Br_boundary_im` (real, default 0.0) -- imaginary part of Br at right boundary

### File Organization

- `electromagnetic/electromagnetic_solver.f90` -- `electromagnetic_t` type and block solve
- `electromagnetic/ampere_matrices.f90` -- weighted mass/potential matrix assembly
- Modified: `KIM_mod.f90` (factory), `setup_m.f90` + `read_config.f90` (namelist)

## Block System

The coupled system is:

```
| A_Phi    C_PhiB | | Phi | = | b_Phi |
| C_BPhi   A_B    | | Br  |   | b_B   |
```

### Block Definitions (each N x N)

- `A_Phi = Delta + 4*pi * K_rho_phi` (existing Poisson LHS)
- `C_PhiB = 4*pi * K_rho_B` (existing kernel)
- `C_BPhi = i*4*pi/c * M * K_j_phi` (mass matrix x existing kernel)
- `A_B = kz * M^h_theta - m * Q^h_z + i*4*pi/c * M * K_j_B` (new + existing)

### Weighted Matrices (Midpoint Approximation)

- `M^h_theta_{l'l} ~ h_theta(r_bar) * M_{l'l}`
- `Q^h_z_{l'l} ~ h_z(r_bar) * Q_{l'l}`

where `r_bar` is the midpoint of overlapping basis function support. The unweighted
`Q_{l'l}` uses the analytical log-expressions from the writeup.

### Solve

The 2N x 2N system is converted to sparse format via existing `dense_to_sparse` and
solved with existing UMFPACK/SuperLU solver (`sparse_solveComplex_b1`). The solution
vector is split: Phi = sol(1:N), Br = sol(N+1:2N).

## Boundary Conditions

Imposed by row/column elimination (same technique as existing `impose_bc_on_matrix_and_rhs`):

- **Phi** (rows 1, N): existing BC types (Dirichlet or zero-misalignment `bc_type=3`)
- **Br** (rows N+1, 2N): Br(left) = 0, Br(right) = Br_boundary from namelist

For zero-misalignment BCs, `Phi_aligned` uses the known boundary Br values.

## Postprocessing & Output

After the solve:

- `EBdat%Phi` <- first N elements
- `EBdat%Br` <- last N elements (self-consistent)
- All existing postprocessing applies unchanged:
  - `calculate_current_density` (uses K_j_phi*Phi + K_j_B*Br)
  - `postprocess_electric_field` (uses both Phi and Br)
  - `calculate_MA_field` (misalignment from both fields)
- Diagnostic output of Ampere blocks at `fdebug >= 2`

## Testing & Validation

1. **Unit test**: Ampere matrix assembly with constant h_theta, h_z -- weighted matrices
   must equal h * M and h * Q exactly.

2. **Regression**: With K_j_phi = K_j_B = 0, Br passes through unchanged from BC,
   Phi matches existing electrostatic result for prescribed Br.

3. **Physics**: Self-consistent Br shows screening at resonant surface
   (|Br_sc(r_res)| < |Br_vac(r_res)|). Quantitative benchmark against KiLCA.
