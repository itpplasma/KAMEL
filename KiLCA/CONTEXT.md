# KiLCA Domain Context

## Purpose

KiLCA is the cylindrical linear plasma-response solver. In the periodic KIM/QL-Balance workflow it remains the initial source of the electromagnetic vacuum magnetic-drive shape.

## Canonical terms

- **Vacuum response:** the KiLCA calculation used to obtain magnetic perturbations without the periodic KIM plasma response.
- **Radial magnetic drive:** the complex \(B^r\) component used to normalize the constant-\(\psi\) periodic KIM forcing.
- **Parallel magnetic drive:** the complex perturbation component parallel to the equilibrium magnetic field, mapped from KiLCA's field representation to \(B_\parallel\).
- **Native mode convention:** KiLCA's signed \((m,n)\) and Fourier phase convention. Interfaces must convert conventions once and retain provenance.

## Ownership boundaries

KiLCA owns:

- the electromagnetic vacuum calculation and its native field components;
- the physical phase and units of the vacuum magnetic drive.

The KIM/QL adapter owns:

- the documented mapping from KiLCA field components, including `Bp`, into \(B^r\) and \(B_\parallel\);
- normalization of the extracted shapes for a periodic KIM request;
- recording source mode, phase, units, and normalization in coupled output.

KiLCA is not used to replace the periodic KIM plasma response in this workflow.

## Required invariants

- Signed mode numbers and the native phase convention are never erased by absolute-value preprocessing.
- \(B_\parallel\) is not silently set to zero in production.
- Setting \(B_\parallel=0\) is an explicit benchmark configuration only.
- One common shielding amplitude scales the initial \(B^r\) and \(B_\parallel\) shapes unless a later reviewed closure establishes separate amplitudes.

## Current design source

Read the repository-root `roadmap.md` and the existing signed-resonance and magnetic-compression GitHub issues before changing this interface.
