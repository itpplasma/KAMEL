# KIM Domain Context

## Purpose

KIM is the KiLCA Integral Model. It computes kinetic plasma response using integral kernels and supports both global and local forced-periodicity calculations.

## Canonical terms

- **Global KIM:** the full radial electromagnetic/FEM calculation. Its results have global radial-domain semantics.
- **Periodic KIM:** the forced-periodicity calculation in a local Fourier window around one rational surface.
- **Resonant surface:** the radius selected by the signed mode and safety-factor convention at which the local window is centred.
- **As-is core:** the central part of a periodic window whose equilibrium profiles are not modified by the periodic transition construction.
- **Transition zone:** the window region used to make the local problem periodic; it is not automatically valid physical output.
- **Local response:** the periodic solution together with its grid, validity metadata, physical fields, species currents, and local transport tensor.
- **Integral moments:** the susceptibility moments \(\mathcal I_{\ell\sigma}^{p,q}\) used by the response and transport contractions.
- **Magnetic drives:** prescribed complex \(B^r\) and \(B_\parallel\) perturbations. In electrostatic periodic KIM they are inputs, not an Ampere self-consistent solution.
- **Integral ion tensor:** the four gyrokinetic coefficients \(D_{11},D_{12},D_{21},D_{22}\) obtained from the quadratic form over \((\Phi,B^r,B_\parallel)\).

## Ownership boundaries

KIM owns:

- construction and solution of the periodic local response;
- Fourier, phase, charge, cyclotron-harmonic, and susceptibility conventions;
- conversion from periodic potential coefficients to local physical electric fields;
- species-resolved parallel-current response;
- the gyrokinetic integral ion transport contraction.

KIM does not own:

- extrapolation of a local response over the complete QL-Balance radius;
- global profile time evolution;
- the target-current normalization policy;
- the legacy drift-kinetic electron coefficients.

## Required invariants

- A library call must restore process-global plasma/grid state before returning, including on failure.
- Reordering independent mode solves must not change their results.
- Complex phase and Fourier normalization must be preserved at every public boundary.
- Periodic results must distinguish the trusted as-is core from the transition zone.
- Scaling both prescribed magnetic drives by \(s\) scales linear fields/current by \(s\) and the transport tensor by \(|s|^2\).
- The \(B_\parallel=0\), zero-FLR, zero-harmonic limit must reproduce the established drift-kinetic tensor.

## Current design source

Read the repository-root `roadmap.md` for the staged periodic KIM/QL-Balance implementation. The companion mathematical derivation and Mathematica verification live in the KIM documentation repository until copied into KAMEL as a checked algebraic oracle.
