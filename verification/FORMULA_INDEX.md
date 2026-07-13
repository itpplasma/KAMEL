# Formula Index and Convention Ledger

This ledger is the source-of-truth inventory for KAMEL mathematical
verification. The thesis source is Markus J. Markl, *Kinetic investigation of
the bifurcation of resonant magnetic perturbations in magnetically confined
plasmas and development of an integral response model*,
DOI `10.3217/efp2p-0x485` (repository file `87404.pdf`, MD5
`bc3a1f7afc861f54f15a3c6f94815bde`). “Indexed” means located and classified;
it does **not** mean independently verified.

## Locked conventions

| ID | Convention | Canonical definition | Code evidence |
| --- | --- | --- | --- |
| C-unit | Units | Gaussian CGS internally: cm, s, g, statC, gauss; profile temperature is eV and is multiplied by `ev` before energy formulas. | `KIM/src/general/constants.f90`, `species_m::calculate_plasma_backs` |
| C-charge | Species and charge | Species 0 is the electron; 1…S−1 are ions. `Zspec` is signed (electron −1), and `e_sigma=Zspec*e_charge`. | `species_m::init_*_species` |
| C-vt | Thermal speed | `vT_sigma=sqrt(T_sigma*ev/m_sigma)`; this is the one-dimensional Maxwellian scale, with no factor of two. | thesis (2.52), (5.8); `species_m::calculate_plasma_backs` |
| C-gyro | Gyrofrequency | `omega_c_sigma=Zspec*e_charge*abs(B0)/(m_sigma*c)` is signed; `rho_L=abs(vT/omega_c)` is non-negative. | thesis (4.14)–(4.16); `species_m::calculate_plasma_backs` |
| C-phase | Space/time phase | Physical perturbations are real parts of amplitudes proportional to `exp(+i(m theta+n phi)-i omega t)`, with `z=R0*phi` and `k_z=n/R0`. | thesis (4.1), (4.3), (4.26); `calculate_equil` |
| C-reality | Complex reality | A real full field has `F[-m,-n]=conjugate(F[m,n])`; KIM may solve one complex harmonic independently. | consequence of C-phase |
| C-radial-ft | Radial Fourier pair | `Phi(r)=integral dk_r exp(+i k_r r) Phi(k_r)/(2 pi)`; the forward transform uses `exp(-i k_r r)`. A period of length `L` uses `k_r=2 pi ell/L`. | thesis (13.3); `fourier_grid_m` |
| C-coords | Coordinates | Right-handed cylindrical `(r,theta,z)`; `z=R0 phi`. `h=B0/|B0|=(0,h_theta,h_z)`; radial is outward. | thesis Ch. 3–4; `calculate_equil` |
| C-field | Field components | `B_r` is contravariant radial perturbation; `E_perp=E·(h×grad r)`. `V_z` is cylindrical toroidal velocity. | thesis (5.10)–(5.11), (5.22); KIM fields/profile input |
| C-q | Safety factor | `q=r B_z/(R0 B_theta)`; the configured resonance is `q_res=abs(m/n)`. Signed `m,n` still enter `k_s,k_parallel`. | thesis (3.5)–(3.6); `prepare_resonances`, `calculate_equil` |
| C-wavevectors | Wavevectors | `k_s=(m h_z-n h_theta/R0)/r`; `k_parallel=m h_theta/r+n h_z/R0`. | `calculate_equil`; thesis (5.19)–(5.20) |
| C-flow | Parallel flow | File `Vz` is projected once: `V_parallel=Vz/h_z`. Ion FP `x2=-(omega_E+k_parallel V_parallel+m_phi omega_c-omega)/nu`; electrons are stationary under the shared ion-flow policy. | `profile_input_m`, `species_m::calculate_thermodynamic_forces_and_susc` |
| C-bessel | Bessel function | `J_n(x)=(1/(2 pi)) integral[-pi,pi] exp(i(n tau-x sin tau)) d tau`. | thesis (13.7); intrinsic/library Bessel calls |
| C-pdf | Plasma dispersion | `Z(z)=i sqrt(pi) w(z)` with analytic continuation fixed by the `exp(-i omega t)` convention. Susceptibility `I^{kl}` follows thesis (5.12)–(5.17) and Appendix A. | `getIfunc`, libcerf |
| C-collision | Collision scale | `nu_sigma` is the total species collision frequency used in `x1=k_parallel vT/nu` and `x2`; it is positive. | thesis (5.16)–(5.17); `species_m` |
| C-velocity-measure | Velocity measure | Moments use `d^3p` (or the exactly transformed velocity measure) and species charge once in the outer species sum. | thesis (4.8), (11.4)–(11.8) |
| C-poisson | Electrostatic equation | Gaussian CGS Poisson equation is `nabla^2 deltaPhi=-4 pi delta rho`; kernel assembly must not absorb a second `4 pi`. | thesis (12.25); global/periodic Poisson solvers |
| C-boundary | Boundary and gauge | Current global `bc_type=3` fixes the aligned-potential boundary. The accepted periodic unknown/gauge is not yet derived; it remains explicitly pending #197/#188. | `poisson.f90`, `periodic_solve.f90` |
| C-parity | Mode/harmonic parity | Cyclotron harmonics `m_phi` are signed integers; no implicit replacement by `abs(m_phi)` is allowed. | thesis Ch. 13–14; FP susceptibility loops |

### Assumption classes

`C1` context; `C2` kinetic premises; `C3` cylindrical equilibrium;
`C4` linear response/action-angle conventions; `C5` differential FLRE;
`C6` quasilinear transport; `C7` realistic-geometry correction; `C8`
profile/equilibrium preparation; `C9` bifurcation analysis; `C10` LHD
application; `C11` constitutive kernels; `C12` unmagnetized toy problem;
`C13` magnetized Krook integral model; `C14` magnetized Fokker–Planck
model; `C15` variational/basis/Poisson formulation; `CA` susceptibility
identities. Consumer keys are `EQ=common/equil`, `PRE=PreProc`,
`KIM=KIM`, `KiLCA=KiLCA`, and `QLB=QL-Balance`.

## Thesis displayed-equation inventory

The extraction contains 494 unique numbered displayed equations: every integer
label from (1.1), chapter-contiguous ranges through (15.33), and (A.1)–(A.29).
Repeated references/reprints are counted once. Printed page and physical PDF
page are both recorded so the inventory is reproducible.

| Equation | Thesis page | Chapter | Assumptions | Consumers | Verification |
| --- | ---: | ---: | --- | --- | --- |
| (1.1) | 4 (PDF 15) | 1 | C1 | — | premise |
| (2.1) | 12 (PDF 23) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.2) | 12 (PDF 23) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.3) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.4) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.5) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.6) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.7) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.8) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.9) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.10) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.11) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.12) | 13 (PDF 24) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.13) | 14 (PDF 25) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.14) | 14 (PDF 25) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.15) | 15 (PDF 26) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.16) | 15 (PDF 26) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.17) | 15 (PDF 26) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.18) | 15 (PDF 26) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.19) | 15 (PDF 26) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.20) | 15 (PDF 26) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.21) | 15 (PDF 26) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.22) | 16 (PDF 27) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.23) | 16 (PDF 27) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.24) | 16 (PDF 27) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.25) | 16 (PDF 27) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.26) | 17 (PDF 28) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.27) | 17 (PDF 28) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.28) | 17 (PDF 28) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.29) | 17 (PDF 28) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.30) | 18 (PDF 29) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.31) | 18 (PDF 29) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.32) | 18 (PDF 29) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.33) | 18 (PDF 29) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.34) | 18 (PDF 29) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.35) | 19 (PDF 30) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.36) | 19 (PDF 30) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.37) | 19 (PDF 30) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.38) | 19 (PDF 30) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.39) | 20 (PDF 31) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.40) | 20 (PDF 31) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.41) | 20 (PDF 31) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.42) | 20 (PDF 31) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.43) | 20 (PDF 31) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.44) | 20 (PDF 31) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.45) | 20 (PDF 31) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.46) | 21 (PDF 32) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.47) | 21 (PDF 32) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.48) | 21 (PDF 32) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.49) | 21 (PDF 32) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.50) | 21 (PDF 32) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.51) | 21 (PDF 32) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.52) | 22 (PDF 33) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.53) | 22 (PDF 33) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (2.54) | 22 (PDF 33) | 2 | C2 | EQ/KIM/KiLCA/QLB | indexed; convention checks |
| (3.1) | 23 (PDF 34) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.2) | 23 (PDF 34) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.3) | 24 (PDF 35) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.4) | 24 (PDF 35) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.5) | 25 (PDF 36) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.6) | 25 (PDF 36) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.7) | 25 (PDF 36) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.8) | 25 (PDF 36) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.9) | 25 (PDF 36) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.10) | 26 (PDF 37) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.11) | 26 (PDF 37) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.12) | 26 (PDF 37) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.13) | 26 (PDF 37) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.14) | 26 (PDF 37) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.15) | 26 (PDF 37) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (3.16) | 26 (PDF 37) | 3 | C3 | EQ/KIM/KiLCA | indexed; oracle #198 |
| (4.1) | 27 (PDF 38) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.2) | 27 (PDF 38) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.3) | 27 (PDF 38) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.4) | 27 (PDF 38) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.5) | 27 (PDF 38) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.6) | 28 (PDF 39) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.7) | 28 (PDF 39) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.8) | 28 (PDF 39) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.9) | 28 (PDF 39) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.10) | 28 (PDF 39) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.11) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.12) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.13) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.14) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.15) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.16) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.17) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.18) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.19) | 29 (PDF 40) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.20) | 30 (PDF 41) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.21) | 30 (PDF 41) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.22) | 30 (PDF 41) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.23) | 30 (PDF 41) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.24) | 30 (PDF 41) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.25) | 30 (PDF 41) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.26) | 31 (PDF 42) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.27) | 31 (PDF 42) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.28) | 31 (PDF 42) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.29) | 31 (PDF 42) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.30) | 31 (PDF 42) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.31) | 31 (PDF 42) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.32) | 31 (PDF 42) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.33) | 32 (PDF 43) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.34) | 32 (PDF 43) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.35) | 32 (PDF 43) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.36) | 32 (PDF 43) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.37) | 32 (PDF 43) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.38) | 32 (PDF 43) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.39) | 33 (PDF 44) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.40) | 33 (PDF 44) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.41) | 33 (PDF 44) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.42) | 33 (PDF 44) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.43) | 33 (PDF 44) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.44) | 33 (PDF 44) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (4.45) | 33 (PDF 44) | 4 | C4 | KIM/KiLCA | indexed; oracles #194/#196 |
| (5.1) | 35 (PDF 46) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.2) | 35 (PDF 46) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.3) | 35 (PDF 46) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.4) | 35 (PDF 46) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.5) | 36 (PDF 47) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.6) | 36 (PDF 47) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.7) | 36 (PDF 47) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.8) | 37 (PDF 48) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.9) | 37 (PDF 48) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.10) | 37 (PDF 48) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.11) | 37 (PDF 48) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.12) | 37 (PDF 48) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.13) | 37 (PDF 48) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.14) | 37 (PDF 48) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.15) | 38 (PDF 49) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.16) | 38 (PDF 49) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.17) | 38 (PDF 49) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.18) | 38 (PDF 49) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.19) | 38 (PDF 49) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.20) | 38 (PDF 49) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.21) | 39 (PDF 50) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.22) | 39 (PDF 50) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.23) | 39 (PDF 50) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.24) | 39 (PDF 50) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.25) | 40 (PDF 51) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.26) | 40 (PDF 51) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.27) | 41 (PDF 52) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.28) | 41 (PDF 52) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.29) | 41 (PDF 52) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.30) | 41 (PDF 52) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.31) | 41 (PDF 52) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.32) | 41 (PDF 52) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.33) | 41 (PDF 52) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.34) | 42 (PDF 53) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.35) | 42 (PDF 53) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.15) | 42 (PDF 53) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.36) | 42 (PDF 53) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.37) | 42 (PDF 53) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.38) | 42 (PDF 53) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.39) | 43 (PDF 54) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.40) | 44 (PDF 55) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.41) | 44 (PDF 55) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.42) | 44 (PDF 55) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (5.43) | 45 (PDF 56) | 5 | C5 | KIM/KiLCA/QLB | indexed; oracle #194 |
| (6.1) | 46 (PDF 57) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.2) | 47 (PDF 58) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.3) | 47 (PDF 58) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.4) | 47 (PDF 58) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.5) | 47 (PDF 58) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.6) | 47 (PDF 58) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.7) | 47 (PDF 58) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.8) | 47 (PDF 58) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.9) | 48 (PDF 59) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.10) | 48 (PDF 59) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.11) | 48 (PDF 59) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.12) | 48 (PDF 59) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.13) | 48 (PDF 59) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.14) | 48 (PDF 59) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.15) | 49 (PDF 60) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.16) | 49 (PDF 60) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.17) | 50 (PDF 61) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.18) | 50 (PDF 61) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (6.19) | 50 (PDF 61) | 6 | C6 | QLB | indexed; broader oracle #199 |
| (7.1) | 51 (PDF 62) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.2) | 51 (PDF 62) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.3) | 51 (PDF 62) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.4) | 52 (PDF 63) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.5) | 52 (PDF 63) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.6) | 52 (PDF 63) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.7) | 52 (PDF 63) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.8) | 53 (PDF 64) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (7.9) | 54 (PDF 65) | 7 | C7 | PRE/EQ/KiLCA/QLB | indexed; oracle #198 |
| (8.1) | 56 (PDF 67) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.2) | 56 (PDF 67) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.3) | 56 (PDF 67) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.4) | 56 (PDF 67) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.5) | 56 (PDF 67) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.6) | 56 (PDF 67) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.7) | 56 (PDF 67) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.8) | 57 (PDF 68) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.9) | 57 (PDF 68) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.10) | 58 (PDF 69) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.11) | 59 (PDF 70) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (8.12) | 59 (PDF 70) | 8 | C8 | PRE/EQ/QLB | indexed; oracle #198 |
| (9.1) | 64 (PDF 75) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.2) | 65 (PDF 76) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.3) | 65 (PDF 76) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.4) | 65 (PDF 76) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.5) | 65 (PDF 76) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.6) | 65 (PDF 76) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.7) | 66 (PDF 77) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.8) | 66 (PDF 77) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.9) | 66 (PDF 77) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.10) | 66 (PDF 77) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.11) | 69 (PDF 80) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.12) | 70 (PDF 81) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.13) | 72 (PDF 83) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.14) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.15) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.16) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.17) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.18) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.19) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.20) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.21) | 73 (PDF 84) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.22) | 74 (PDF 85) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.23) | 74 (PDF 85) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.24) | 75 (PDF 86) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.25) | 75 (PDF 86) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.26) | 76 (PDF 87) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.27) | 76 (PDF 87) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.28) | 77 (PDF 88) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (9.29) | 77 (PDF 88) | 9 | C9 | QLB | indexed; broader oracle #199 |
| (10.1) | 83 (PDF 94) | 10 | C10 | — | indexed; published criterion |
| (10.2) | 83 (PDF 94) | 10 | C10 | — | indexed; published criterion |
| (11.1) | 86 (PDF 97) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.2) | 86 (PDF 97) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.3) | 86 (PDF 97) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.4) | 86 (PDF 97) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.5) | 86 (PDF 97) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.6) | 86 (PDF 97) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.7) | 87 (PDF 98) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.8) | 87 (PDF 98) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (11.9) | 87 (PDF 98) | 11 | C11 | KIM/KiLCA | indexed; oracles #196/#197 |
| (12.1) | 88 (PDF 99) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.2) | 88 (PDF 99) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.3) | 88 (PDF 99) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.4) | 88 (PDF 99) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.5) | 88 (PDF 99) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.6) | 88 (PDF 99) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.7) | 88 (PDF 99) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.8) | 89 (PDF 100) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.9) | 89 (PDF 100) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.10) | 89 (PDF 100) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.11) | 89 (PDF 100) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.12) | 89 (PDF 100) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.13) | 90 (PDF 101) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.14) | 90 (PDF 101) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.15) | 90 (PDF 101) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.16) | 90 (PDF 101) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.17) | 90 (PDF 101) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.18) | 90 (PDF 101) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.19) | 91 (PDF 102) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.20) | 91 (PDF 102) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.21) | 91 (PDF 102) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.22) | 91 (PDF 102) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.23) | 92 (PDF 103) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.24) | 92 (PDF 103) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.25) | 92 (PDF 103) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.26) | 92 (PDF 103) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.27) | 92 (PDF 103) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (12.28) | 92 (PDF 103) | 12 | C12 | KIM tests | indexed; oracle #197 |
| (13.1) | 93 (PDF 104) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.2) | 93 (PDF 104) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.3) | 93 (PDF 104) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.4) | 94 (PDF 105) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.5) | 94 (PDF 105) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.6) | 94 (PDF 105) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.7) | 94 (PDF 105) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.8) | 94 (PDF 105) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.9) | 94 (PDF 105) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.10) | 95 (PDF 106) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.11) | 95 (PDF 106) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.12) | 95 (PDF 106) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.13) | 95 (PDF 106) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.14) | 95 (PDF 106) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.15) | 95 (PDF 106) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.16) | 95 (PDF 106) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.17) | 96 (PDF 107) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.18) | 96 (PDF 107) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.19) | 96 (PDF 107) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.20) | 96 (PDF 107) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.21) | 96 (PDF 107) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.22) | 96 (PDF 107) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.23) | 96 (PDF 107) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.24) | 97 (PDF 108) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.25) | 97 (PDF 108) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.26) | 97 (PDF 108) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.27) | 97 (PDF 108) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.28) | 97 (PDF 108) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.29) | 97 (PDF 108) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.30) | 98 (PDF 109) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.31) | 98 (PDF 109) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.32) | 98 (PDF 109) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.33) | 98 (PDF 109) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.34) | 99 (PDF 110) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.35) | 99 (PDF 110) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.36) | 99 (PDF 110) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.37) | 99 (PDF 110) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.38) | 99 (PDF 110) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.39) | 99 (PDF 110) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.40) | 100 (PDF 111) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.41) | 100 (PDF 111) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.42) | 100 (PDF 111) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.43) | 101 (PDF 112) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.44) | 101 (PDF 112) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.45) | 101 (PDF 112) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.46) | 101 (PDF 112) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.47) | 101 (PDF 112) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.48) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.49) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.50) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.51) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.52) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.53) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.54) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.55) | 102 (PDF 113) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.56) | 103 (PDF 114) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.57) | 103 (PDF 114) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.58) | 103 (PDF 114) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.59) | 103 (PDF 114) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.60) | 103 (PDF 114) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.61) | 103 (PDF 114) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.62) | 104 (PDF 115) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.63) | 104 (PDF 115) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.64) | 104 (PDF 115) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.65) | 104 (PDF 115) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.66) | 104 (PDF 115) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.67) | 104 (PDF 115) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.68) | 105 (PDF 116) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.69) | 105 (PDF 116) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.70) | 105 (PDF 116) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.71) | 105 (PDF 116) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.72) | 105 (PDF 116) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.73) | 105 (PDF 116) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.74) | 106 (PDF 117) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.75) | 106 (PDF 117) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.76) | 106 (PDF 117) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.77) | 106 (PDF 117) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.78) | 106 (PDF 117) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.79) | 106 (PDF 117) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.80) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.81) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.82) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.83) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.84) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.85) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.86) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.87) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.88) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.89) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.90) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.91) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.92) | 107 (PDF 118) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.93) | 108 (PDF 119) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.94) | 108 (PDF 119) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.95) | 108 (PDF 119) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.96) | 108 (PDF 119) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.97) | 108 (PDF 119) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.98) | 108 (PDF 119) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.99) | 108 (PDF 119) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.100) | 109 (PDF 120) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.101) | 109 (PDF 120) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.102) | 109 (PDF 120) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.103) | 109 (PDF 120) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.104) | 109 (PDF 120) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.105) | 109 (PDF 120) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.106) | 109 (PDF 120) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.107) | 110 (PDF 121) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.108) | 110 (PDF 121) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.109) | 110 (PDF 121) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.110) | 110 (PDF 121) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.111) | 110 (PDF 121) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.112) | 110 (PDF 121) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.113) | 110 (PDF 121) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.114) | 111 (PDF 122) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.115) | 111 (PDF 122) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.116) | 111 (PDF 122) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.117) | 111 (PDF 122) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.118) | 111 (PDF 122) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.119) | 111 (PDF 122) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.120) | 111 (PDF 122) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.121) | 112 (PDF 123) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.122) | 112 (PDF 123) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.123) | 112 (PDF 123) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.124) | 112 (PDF 123) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.125) | 112 (PDF 123) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.126) | 112 (PDF 123) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.127) | 113 (PDF 124) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.128) | 113 (PDF 124) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.129) | 113 (PDF 124) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.130) | 113 (PDF 124) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.131) | 114 (PDF 125) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.132) | 114 (PDF 125) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.133) | 114 (PDF 125) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.134) | 114 (PDF 125) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.135) | 115 (PDF 126) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.136) | 115 (PDF 126) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.137) | 115 (PDF 126) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.138) | 115 (PDF 126) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.139) | 115 (PDF 126) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.140) | 115 (PDF 126) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.141) | 115 (PDF 126) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.142) | 116 (PDF 127) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.143) | 116 (PDF 127) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.144) | 116 (PDF 127) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.145) | 116 (PDF 127) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.146) | 116 (PDF 127) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.147) | 116 (PDF 127) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.148) | 117 (PDF 128) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.149) | 117 (PDF 128) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.150) | 117 (PDF 128) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.151) | 118 (PDF 129) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.152) | 118 (PDF 129) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.153) | 118 (PDF 129) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.154) | 118 (PDF 129) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.155) | 118 (PDF 129) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.156) | 119 (PDF 130) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.157) | 119 (PDF 130) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.158) | 119 (PDF 130) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (13.159) | 119 (PDF 130) | 13 | C13 | KIM/KiLCA | indexed; oracles #196/#197 |
| (14.1) | 120 (PDF 131) | 14 | C14 | KIM | indexed; oracles #194/#196 |
| (14.2) | 120 (PDF 131) | 14 | C14 | KIM | indexed; oracles #194/#196 |
| (14.3) | 120 (PDF 131) | 14 | C14 | KIM | indexed; oracles #194/#196 |
| (14.4) | 121 (PDF 132) | 14 | C14 | KIM | indexed; oracles #194/#196 |
| (14.5) | 121 (PDF 132) | 14 | C14 | KIM | indexed; oracles #194/#196 |
| (14.6) | 121 (PDF 132) | 14 | C14 | KIM | indexed; oracles #194/#196 |
| (15.1) | 122 (PDF 133) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.2) | 122 (PDF 133) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.3) | 122 (PDF 133) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.4) | 123 (PDF 134) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.5) | 123 (PDF 134) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.6) | 123 (PDF 134) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.7) | 124 (PDF 135) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.8) | 124 (PDF 135) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.9) | 124 (PDF 135) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.10) | 124 (PDF 135) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.11) | 124 (PDF 135) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.12) | 125 (PDF 136) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.13) | 125 (PDF 136) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.14) | 125 (PDF 136) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.15) | 126 (PDF 137) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.16) | 126 (PDF 137) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.17) | 126 (PDF 137) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.18) | 127 (PDF 138) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.19) | 127 (PDF 138) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.20) | 127 (PDF 138) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.21) | 127 (PDF 138) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.22) | 127 (PDF 138) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.23) | 128 (PDF 139) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.24) | 128 (PDF 139) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.25) | 128 (PDF 139) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.26) | 128 (PDF 139) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.27) | 128 (PDF 139) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.28) | 128 (PDF 139) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.29) | 129 (PDF 140) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.30) | 129 (PDF 140) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.31) | 129 (PDF 140) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.32) | 130 (PDF 141) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (15.33) | 130 (PDF 141) | 15 | C15 | KIM periodic | indexed; oracle #197 |
| (A.1) | 136 (PDF 147) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.2) | 136 (PDF 147) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.3) | 136 (PDF 147) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.4) | 136 (PDF 147) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.5) | 136 (PDF 147) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.6) | 137 (PDF 148) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.7) | 137 (PDF 148) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.8) | 137 (PDF 148) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.9) | 137 (PDF 148) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.10) | 137 (PDF 148) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.11) | 137 (PDF 148) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.12) | 138 (PDF 149) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.13) | 138 (PDF 149) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.14) | 138 (PDF 149) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.16) | 138 (PDF 149) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.17) | 138 (PDF 149) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.18) | 138 (PDF 149) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.19) | 139 (PDF 150) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.20) | 139 (PDF 150) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.21) | 139 (PDF 150) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.22) | 139 (PDF 150) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.23) | 139 (PDF 150) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.24) | 139 (PDF 150) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.25) | 140 (PDF 151) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.26) | 140 (PDF 151) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.27) | 140 (PDF 151) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.28) | 141 (PDF 152) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |
| (A.29) | 141 (PDF 152) | A | CA | KIM/KiLCA/QLB | indexed; oracle #194 |

## Physics-bearing code formula sites

This table is the reviewable semantic inventory. A “site” is a stable routine
or symbol, not every continuation line. Generic interpolation, indexing, and
I/O arithmetic are excluded. The Mathematica gate verifies that every listed
path exists and that each code family has a source-to-check chain.

<!-- CODE-SITES-BEGIN -->
| Site | Code symbol | Class | Thesis mapping | Verification |
| --- | --- | --- | --- | --- |
| KIM-01 | `KIM/src/background_equilibrium/species_mod.f90::calculate_plasma_backs` | vT, gyrofrequency, Larmor/Debye lengths, collision scale | (4.14)–(4.16), (5.8) | conventions script; #198 pending |
| KIM-02 | `KIM/src/background_equilibrium/species_mod.f90::calculate_thermodynamic_forces_and_susc` | A1/A2, x1/x2, FP susceptibilities | (5.8), (5.12)–(5.20) | flow test; independent #194 oracle |
| KIM-03 | `KIM/src/background_equilibrium/calculate_equil.f90::calculate_equil` | h_theta/h_z, k_s, k_parallel, omega_E | (3.5), (5.18)–(5.20) | conventions script; #198 pending |
| KIM-04 | `KIM/src/asymptotics/flr2_fourier_kernel.f90` | b_plus, b_cross, phases, charge/current kernels | Ch. 14 | #196 pending |
| KIM-05 | `KIM/src/kernels/FP_kernel_plasma_prefacs.f90` | finite-radius FP kernel prefactors | Ch. 14 | #196 pending |
| KIM-06 | `KIM/src/kernels/kernel.f90` | real-space charge/current kernels and 1/(4 pi) normalization | Ch. 13–14 | #196/#197 pending |
| KIM-07 | `KIM/src/electrostatic_poisson/poisson.f90` | global Gaussian-CGS weak Poisson problem | (12.25), Ch. 15 | #197 pending |
| KIM-08 | `KIM/src/electrostatic_poisson/periodic_solve.f90` | periodic Fourier Poisson matrix/RHS | (12.25), Ch. 15 | #197/#188 pending |
| KIM-09 | `KIM/src/asymptotics/FLR2_asymptotics.f90` | aligned potential and FLR2 limits | Ch. 5 | #197 pending |
| KIM-10 | `KIM/src/grid/prepare_resonances.f90` | q_res=abs(m/n) | (3.5)–(3.6) | conventions script |
| KIM-11 | `KIM/src/background_equilibrium/profile_input_m.f90` | radial force balance and Vz-to-Vparallel projection | (8.8)–(8.9) | profile/flow tests |
| KIM-12 | `KIM/src/background_equilibrium/periodic_background.f90` | periodic primitive state and derivative reconstruction | Ch. 3, 8 | #198 pending |
| KILCA-01 | `KiLCA/flre/conductivity/calc_I_array.f90::calc_Imn_array` | susceptibility-function definition | (5.12)–(5.17), App. A | not exercised by KIM #194 oracle |
| KILCA-02 | `KiLCA/solver/VER_5_STABLE/wave_stuff.f90` | plasma/cyclotron frequencies and dielectric response | Ch. 4–5 | not exercised by KIM #194 oracle; #197 pending |
| KILCA-03 | `KiLCA/flre/conductivity/calc_I_array_drift_serg.f90` | drift/FP conductivity integrals | Ch. 5, App. A | not exercised by KIM #194 oracle |
| KILCA-04 | `PreProc/fourier/src/rhs_flt.f90` | exp(-i m iota phi) field-line Fourier phase | (4.3) | conventions script |
| QLB-01 | `QL-Balance/src/base/getIfunc.f90` | susceptibility functions | App. A | independent #194 oracle |
| QLB-02 | `QL-Balance/src/base/get_dql.f90` | quasilinear transport coefficients | (6.11)–(6.16) | broader #199 |
| QLB-03 | `QL-Balance/src/base/rhs_balance_m.f90::compute_particle_fluxes` | particle fluxes | (6.2), (6.11)–(6.16) | unit tests; broader #199 |
| QLB-04 | `QL-Balance/src/base/rhs_balance_m.f90::compute_total_heat_fluxes` | heat fluxes | (6.9), (6.11)–(6.16) | unit tests; broader #199 |
| QLB-05 | `QL-Balance/src/base/calc_current_densities.f90` | current moments and thermal speeds | (4.8), (5.8) | broader #199 |
| QLB-06 | `QL-Balance/src/base/W2_arr.f90::W2_arr` | production differential FP moments consumed by KIM | (5.14)–(5.15), (A.1)–(A.3) | independent #194 oracle |
| EQ-01 | `common/equil/mag_wrapper.f90::magfie` | cylindrical metric and field unit vector | Ch. 3 | #198 pending |
| EQ-02 | `common/equil/equil_profiles.f90` | flux-surface equilibrium profiles | Ch. 3, 8 | #198 pending |
| PRE-01 | `PreProc/fourier/src/fouriermodes.f90` | straight-field-line angle and equilibrium Fourier modes | (4.3), Ch. 7–8 | #198 pending |
<!-- CODE-SITES-END -->

## Coverage policy

### Fokker–Planck susceptibility oracle (#194)

`verification/mathematica/02_kinetic_and_susceptibilities.wl` independently
evaluates the generating integral (5.14)–(5.15), applies the energy-conserving
rank-one correction (5.13), and traces `x1/x2` through (4.1), (4.26),
(5.16)–(5.20), (2.52)–(2.54), and Appendix (A.1)–(A.29). The generated
`verification/oracles/fp_susceptibilities.dat` covers every `I00`, `I20`,
`I02`, `I01`, `I21`, `I22`, `I11`, and `I13` consumed by KIM, including
finite flow, `m_phi=-1,0,+1`, near-resonance, strong-collision, and
collisionless-scaled points. `test_fp_susceptibility_oracle` reads this file
without requiring Mathematica and compares it to the production implementation.

The inventory script requires all 494 thesis labels, all 19 convention rows,
and all 25 semantic code sites. It rejects a missing equation, missing code
site/path, incompatible duplicate convention, or unclassified convention.
Subsequent oracle issues replace “pending” statuses with deterministic
script/test links; copying a formula here never upgrades its status.
