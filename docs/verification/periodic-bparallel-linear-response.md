# Periodic KIM parallel-magnetic linear response

This note records the field mapping and the executable reduction used for
issue #288. It covers only the linear prescribed-drive columns; the quadratic
transport tensor remains the responsibility of the existing #299 stack.

## KiLCA field mapping

KiLCA returns the six complex field amplitudes in cylindrical coordinates and
then calls `transform_EB_from_cyl_to_rsp`. The interface copies `EBrsp[5]`
directly to `wave_code_data%Bp`. The underlying RSP transform is

\[
B_s=h_z B_\theta-h_\theta B_z,\qquad
B_p=h_\theta B_\theta+h_z B_z,
\]

where \((h_\theta,h_z)=\mathbf B_0/|\mathbf B_0|\). Consequently KiLCA's
`Bp` is the physical \(B_\parallel=\mathbf B_1\mathbin{\cdot}\mathbf B_0/|\mathbf
B_0|\), with no sign or normalization conversion. `Br` and `Bp` are both CGS
magnetic-field amplitudes in Gauss.

KiLCA differentiates a mode as \(i m/r\) and \(i n/R_0\), so the retained
spatial convention is \(\exp[i(m\theta+n z/R_0)]\). Its time convention is
\(\operatorname{Re}[F\exp(-i\omega t)]\). The adapter asks KiLCA for the exact
signed `(m,n)` entry, copies both complex amplitudes without conjugation, and
linearly samples both profiles at the same signed \(q=-m/n\) resonance. It records
that radius, `(m,n)`, and both sampled Gauss amplitudes before the solve. The
later QL-Balance current normalization multiplies the complete response by one
common complex scalar, preserving the relative phase of `Br` and `Bparallel`.

## Source-only gyrogeometry derivative

For the supported Fokker--Planck, \(m_\phi=0\), full-periodic model define

\[
k_o=\sqrt{k_s^2+k_r^2},\quad
k_g=\sqrt{k_s^2+k_r'^2},\quad
b_+=\frac{\rho^2}{2}(k_o^2+k_g^2),\quad
b_\times=\rho^2 k_o k_g.
\]

The revised paper differentiates only source gyrogeometry. Thus \(k_o\), the
radial phase, detuning, equilibrium forces, and susceptibility moments are held
fixed, while

\[
\partial_{k_s^{(g)}}b_+=\rho^2 k_s,\qquad
\partial_{k_s^{(g)}}b_\times
 =\rho^2 k_o\frac{k_s}{k_g}.
\]

For

\[
g_0=e^{-b_+}I_0(b_\times),\qquad
g_1=e^{-b_+}b_\times I_1(b_\times),
\]

the implementation uses

\[
\partial g_0=-(\partial b_+)g_0+(\partial b_\times)e^{-b_+}I_1,
\qquad
\partial g_1=-(\partial b_+)g_1
 +(\partial b_\times)b_\times g_0,
\]

where the second identity follows from
\(d[xI_1(x)]/dx=xI_0(x)\). Applying this derivative to the paper's common
Gaussian--Bessel bracket gives the two per-species cores

\[
K_a^{\rho B_\parallel}=-i\frac{v_{Ta}^2}
 {\lambda_{Da}^2\nu_a c}\,\partial_{k_s^{(g)}}
 \mathcal G_a[I_{00},I_{02}],
\]

\[
K_a^{j_\parallel B_\parallel}=-i\frac{v_{Ta}^3}
 {\lambda_{Da}^2\nu_a c}\,\partial_{k_s^{(g)}}
 \mathcal G_a[I_{10},I_{12}].
\]

Assembly supplies the unchanged radial Fourier phase and \(1/(8\pi^2)\)
normalization, then sums species. The charge column enters the Poisson
right-hand side as
\(-4\pi B_\parallel K^{\rho B_\parallel}_{:,0}\); the current column enters
the reconstructed current as
\(B_\parallel K^{j_\parallel B_\parallel}_{:,0}\).

## Executable oracle and support boundary

`test_bparallel_linear_response` evaluates the Gaussian--Bessel bracket by a
direct 20,000-point gyrophase quadrature and differentiates it with a centered
source-only finite difference. This is independent of the production Bessel
identity. It also constructs an off-diagonal assembled matrix element directly,
including the radial Fourier phase, quadrature weight, and $1/(8\pi^2)$
normalization. Deliberately reversing that phase, changing the normalization,
or swapping observation and source arguments makes the test fail. Sentinel
values in unrelated moment routes detect index transposition, and a
full-geometry derivative mutation is required to differ. The solve test
independently checks a manufactured Bparallel-only matrix, exact zero-drive
regression, and common-complex-scalar linearity.

The public contracts are tested separately. An optional `KIM_PERIODIC`
namelist with an unknown field or a missing closing delimiter is rejected
instead of silently restoring or accepting a partial drive;
per-solve drive overrides are restored before the next solve; and the
QL-Balance adapter test compares a Bparallel-only run with an exact zero-drive
run at public fields and electron/ion currents. Disconnecting the Bparallel
argument at that adapter seam makes the behavioral test fail.

Before a periodic run redirects KIM onto its local window, the solver handle
keeps a deep snapshot of the authoritative global plasma and radial grid. Each
later mode restores that snapshot, invalidates the periodic-background cache,
and recomputes mode-dependent geometry. A two-resonance ordering test omits any
intervening profile reset and compares both ordered results with fresh
single-mode solves; retaining the first local window or its geometry fails the
test.

Profile updates likewise restore the global state before injection and clear
all grid-sized derived arrays. A 128-point global-grid fixture followed by the
96-point periodic window compares an updated solve with a fresh solve; keeping
window-sized cell-center or susceptibility storage fails with bounds checking.

Adaptive and non-equidistant global meshes are rebuilt around the current
signed resonance after either a mode change or profile update. Independent
mode-order and shifted-profile fixtures compare a reused handle with a genuinely
fresh handle whose configuration starts at the target resonance; retaining the
first resonance-centred mesh makes both comparisons fail.
The same lifecycle reset invalidates the electromagnetic FEM mass matrix and
cell-centred kinetic prefactors before the regenerated field mesh is used. An
independent element-by-element linear hat-function assembly checks the cached
mass matrix against the current adaptive nodes after a mode change; the
bounds-checked two-mode test also exercises prefactor reallocation when the
cell count changes.

Sequential configuration reads reset every optional periodic-window value
before parsing, so a file without `KIM_PERIODIC` receives the documented
defaults rather than values from a preceding handle. A separate lifecycle
fixture transitions from an adaptive in-memory profile centred at one known
crossing to an adaptive file-backed profile centred at another, then repeats
the file-backed initialization. It checks the independently known signed
crossings and the grid generator's cached resonance, and fails if profile-source
mode, primitive allocations, or adaptive-grid state leak.

KIM and the KiLCA adapter use the same signed $q=-m/n$ resonance locator for
both increasing and decreasing $q$ profiles. Wrong-helicity modes with no
signed crossing and profiles with multiple crossings are rejected. The
adapter locates the crossing on the
authoritative QL-Balance profile, samples both vacuum drives there, and passes
that physical radius to KIM as a solve-scoped override. KIM validates the
radius against its local profile and reports the exact supplied value; the
override is restored after the solve. This avoids moving the sampled drive to
a nearby crossing of KIM's independently discretized profile. A curved
reversed-shear fixture uses nonconstant complex Br and Bparallel profiles and
checks the analytic crossing radius, both interpolated drives, and KIM's
reported radius. Omitting the override or restoring the former
increasing-only locator makes the test fail.

A nonzero drive is accepted only for Fokker--Planck species kernels, an active
FP response (`artificial_debye_case=0` or `2`), `mphi_max=0`, and the full
periodic gyrogeometry. Other combinations fail validation. An exactly zero
Bparallel drive bypasses the new assembly and arithmetic, preserving the
pre-existing Phi/Br path bit for bit.
