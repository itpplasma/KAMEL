# Jacobian Probing in rhs_balance_m.f90

## High-Level Overview

The Jacobian probing in `rhs_balance_m.f90` is a technique for efficiently computing the sparse Jacobian matrix **A** that appears in the linearized 1D transport equations:

```
dy/dt = A·y + q
```

where:
- **y** is the state vector containing plasma parameters (density n, toroidal velocity Vphi, electron temperature Te, ion temperature Ti) at all radial grid points (size 4N)
- **A** is the sparse Jacobian matrix representing how changes in each parameter affect the time derivatives
- **q** is the source vector (external heating, fueling, etc.)

The implicit time integrator (CVODE/SUNDIALS) needs this Jacobian to efficiently solve the stiff transport equations.

## Why Probing?

Instead of computing the full dense Jacobian analytically (which would be enormous and mostly zeros), the code uses **probing** to:

1. **Discover the sparsity pattern** - identify which entries are non-zero
2. **Compute only non-zero entries** - save memory and computation
3. **Use finite differences** - approximate derivatives numerically by perturbing each state variable

The Jacobian is sparse because transport is local: a change in temperature at one radial point only affects neighboring points through diffusion/convection, not distant points.

## The Two-Pass Algorithm

### Pass 1: Count Non-Zeros (isw_rhs = 0)

```fortran
call initialize_rhs(y, dy)
  └─> calls rhs_balance with isw_rhs=0
      └─> probes each column, counts non-zeros
      └─> stores count in nz
  └─> allocates sparse arrays: irow(nz), icol(nz), amat(nz)
```

**Purpose**: Determine how much memory to allocate for sparse storage before filling the matrix.

### Pass 2: Fill Matrix (isw_rhs = 1)

```fortran
call rhs_balance(x, y, dy)
  └─> loops over all columns (iprobe = 1 to neqset)
  └─> for each column:
      • perturb y_lin(iprobe) by small amount
      • compute resulting dy (time derivatives)
      • store non-zero dy(i) values as Jacobian entries
      • J[i,iprobe] = dy(i) / y_lin(iprobe)
  └─> stores in COO format: irow(k), icol(k), amat(k)
  └─> also computes source vector q via rhs_balance_source
```

## The Probing Loop in Detail

For each column `iprobe` (each state variable):

### 1. Determine Affected Range

```fortran
ibeg = iprobe / nbaleqs - nshift
iend = iprobe / nbaleqs + nshift
```

Since transport is local, perturbing variable at point `i` only affects points `i ± nshift` (typically nshift=4). This dramatically reduces computation.

### 2. Set Perturbation

```fortran
if (abs(y(iprobe)) > 1.0e-30_dp) then
    y_lin(iprobe) = 1.0e-6_dp * y(iprobe)  ! Relative perturbation
else
    y_lin(iprobe) = 1.0e-6_dp  ! Fallback for zero values
end if
```

**Key insight**: Uses relative perturbation (~10⁻⁶ × current value) instead of absolute perturbation. This prevents:
- Extreme Jacobian entries near resonances
- Numerical issues when probing near-zero values

### 3. Propagate Perturbation

The perturbation is propagated through the physics:

```fortran
! Copy to params_lin grid
params_lin(ieq, ipoi) = y_lin(i)

! Interpolate to cell boundaries
ddr_params_lin = derivatives via deriv_coef
params_b_lin = boundary values via reint_coef

! Compute linearized electric field
Ercov_lin = sqrt_g_times_B_theta_over_c * params_b_lin(2, :) + ...

! Compute linearized fluxes at boundaries
call compute_fluxes_at_boundary(...)
  └─> computes particle fluxes Gamma_e, Gamma_i
  └─> computes heat fluxes Q_e, Q_i
  └─> including both classical and quasi-linear contributions
```

### 4. Compute Response (dy)

```fortran
do ipoi = ibeg, iend
    call compute_dot_params_at_point(...)
    dot_params(:, ipoi) = dot_params_loc
end do

! Copy to dy vector
dy(i) = dot_params(ieq, ipoi)
```

This computes how the time derivatives change due to the perturbation.

### 5. Store Non-Zero Entries

```fortran
do i = ibegtot, iendtot
    if (dy(i) /= 0.0_dp) then
        k = k + 1
        if (isw_rhs == 1) then
            irow(k) = i          ! Row index
            icol(k) = iprobe     ! Column index (perturbed variable)
            amat(k) = dy(i) / y_lin(iprobe)  ! J[i,iprobe]
        end if
    end if
end do
```

**Sparse format**: Stores only non-zero entries in COO (Coordinate) format:
- `irow(k)` - row index
- `icol(k)` - column index
- `amat(k)` - matrix value

### 6. Clean Up

```fortran
y_lin(iprobe) = 0.0_dp
params_lin(:, :) = 0.0_dp
...
dy = 0.0_dp
```

Reset all linearized quantities to zero before probing next column.

## Physical Interpretation

The Jacobian entries represent **linearized coupling** between state variables:

- **J[i,j]** tells how quickly variable `i` changes when variable `j` is perturbed
- **Diagonal dominance**: Variables are most sensitive to themselves at the same radial point
- **Off-diagonal structure**: Reflects diffusive/convective coupling to neighbors
- **Block structure**: Each radial point has 4 variables (n, Vphi, Te, Ti) that couple strongly

Example couplings:
- `∂(dn/dt)/∂n` - how density evolution depends on density (diagonal, typically negative for stability)
- `∂(dn/dt)/∂Te` - how temperature affects density evolution (off-diagonal, through thermal forces)
- `∂(dn/dt)/∂n_{i+1}` - how neighboring density affects local density (super-diagonal, through diffusion)

## Integration with Time Stepper

The sparse Jacobian is used by CVODE (from SUNDIALS) to solve the implicit system:

```
y^{n+1} - Δt·A(y^n)·y^{n+1} = y^n + Δt·q(y^n)
```

CVODE uses the Jacobian to:
1. Form the iteration matrix `I - Δt·A`
2. Solve linear systems during Newton iterations
3. Adapt time step based on stiffness (eigenvalues of A)

The sparse format is crucial because:
- Full matrix would be (4N)×(4N) ≈ millions of entries
- Sparse matrix has only ~20N entries (0.5% fill for N=100)
- Enables efficient storage and linear solves

## Key Design Decisions

1. **Local perturbations**: Only perturb and compute in affected region (±nshift points)
   - Reduces work from O(N²) to O(N)
   
2. **Relative perturbations**: Scale perturbation by current value
   - Prevents numerical issues at resonances
   - Better conditioning of Jacobian

3. **Two-pass algorithm**: Count then fill
   - Enables exact memory allocation
   - Avoids dynamic resizing

4. **COO format**: Simple coordinate format
   - Easy to construct incrementally
   - Can be converted to CSR/CSC for solvers

5. **Separate source computation**: `rhs_balance_source` computes q independently
   - Cleaner separation of concerns
   - Allows checking residual: `dy_check = A·y + q`

## Performance Characteristics

For typical problem with N=100 radial points (400 equations):

- **Full Jacobian**: 160,000 entries, mostly zeros
- **Sparse Jacobian**: ~8,000 non-zero entries (~5% fill)
- **Memory savings**: 20× reduction
- **Computational savings**: Loop over 400 columns, but only ~20 entries per column on average
- **Time complexity**: O(N × nshift × nbaleqs) ≈ O(N) for fixed nshift

## Summary

The Jacobian probing in `rhs_balance_m.f90` is an elegant finite-difference scheme that:

1. **Exploits locality** of transport equations to reduce computational cost
2. **Constructs sparse representation** to save memory and enable efficient solves
3. **Uses physically-motivated perturbations** to maintain numerical stability
4. **Integrates seamlessly** with SUNDIALS implicit ODE solver
5. **Enables robust simulation** of stiff plasma transport on long time scales

The technique is essential for making implicit time integration feasible for these multi-scale, stiff transport problems where explicit methods would require prohibitively small time steps.
