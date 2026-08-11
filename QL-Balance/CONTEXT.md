# QL-Balance Domain Context

## Purpose

QL-Balance evolves global one-dimensional plasma profiles using neoclassical, anomalous, and quasilinear transport. Wave solvers provide mode-resolved perturbations and currents through adapters.

## Canonical terms

- **Global profiles:** the QL-Balance radial density, temperature, safety-factor, electric-field, and rotation profiles advanced in time.
- **Wave-code adapter:** the boundary that maps QL profiles and modes into a wave solver and maps its response back to QL data.
- **Local-to-global embedding:** the explicit placement of a periodic KIM response on the global QL radial grid.
- **Compact transition:** a smooth weight equal to one in the trusted local core and zero outside its bounded support.
- **Misalignment field:** the physical electric-field contribution associated with the mismatch between equipotential and perturbed magnetic surfaces. Its cancelling components must retain their relative phase through embedding.
- **Target-current normalization:** the accepted closure in which the per-mode magnetic amplitude is chosen so the integrated shielding current matches `I_par_toroidal`.
- **Shielding amplitude:** the persistent complex per-mode field scale \(s_{mn}\). Linear response quantities scale with \(s_{mn}\), while diffusion coefficients scale with \(|s_{mn}|^2\).
- **Species transport policy:** electrons use the existing drift-kinetic Heyn/Markl coefficients; ions use KIM's gyrokinetic integral coefficients.

## Ownership boundaries

QL-Balance owns:

- global profiles and accepted/trial time-integration state;
- the KIM wave-code adapter and local-to-global transition policy;
- the legacy drift-kinetic electron coefficient calculation;
- mode accumulation on the global radial grid;
- target-current shielding normalization, guards, rollback, and restart state;
- global output and transport provenance.

QL-Balance does not own:

- KIM susceptibility moments or the gyrokinetic integral ion contraction;
- KIM's local Fourier differentiation and field conventions;
- KiLCA's electromagnetic vacuum response calculation.

## Required invariants

- Derive local physical electric fields before applying a transition; differentiating a tapered potential is forbidden because it creates a spurious window-derivative field.
- Apply the same transition to cancelling misalignment components.
- A local tensor embedded with field weight \(w\) scales as \(w^2D_{ij}\).
- Independent mode diffusion contributions accumulate incoherently at the established assembly point.
- The shielding amplitude is committed only with an accepted time step and is restored on rejection/restart.
- The compatibility `antenna_factor`, if retained, equals \(|s_{mn}|^2\) and is never applied twice.
- Production dispatch uses drift-kinetic electrons and integral ions explicitly; benchmark evaluation does not silently change the evolution model.

## Current design source

Read the repository-root `roadmap.md`, especially the local embedding, target-current evolution, species policy, and benchmark milestones.
