(* ::Package:: *)

(*
  Reproducible m_phi = 0 reduction of the collisionless Fourier kernels.

  Source: PhD_thesis.pdf, equations (13.53)--(13.59), (13.63),
  (13.67), (13.151), and (13.153).

  KIM convention:
    kAbs  = Sqrt[kParallel^2 + epsilon^2]
    kPole = kParallel + I epsilon
    z0    = -(omegaE - omega)/(Sqrt[2] vT kAbs), kept real

  The common continuous-Fourier normalization 1/(8 Pi^2) is included
  in the four G expressions below; the radial Fourier phase is not.
*)

ClearAll["Global`*"];

plasmaZ[z_] := I Sqrt[Pi] Exp[-z^2] Erfc[-I z];

kAbs = Sqrt[kParallel^2 + epsilon^2];
kPole = kParallel + I epsilon;
z0 = -(omegaE - omega)/(Sqrt[2] vT kAbs);

bPlus = rhoL^2 (2 ks^2 + kr^2 + krp^2)/2;
bCross = rhoL^2 Sqrt[ks^2 + kr^2] Sqrt[ks^2 + krp^2];
s0 = Exp[-bPlus] BesselI[0, bCross];
s1 = Exp[-bPlus] BesselI[-1, bCross];

(* Equations (13.53)--(13.55), m_phi = 0.  The kAbs replacement in a1
   is the same even-denominator regularization used by KIM's existing
   collisionless global charge kernel. *)
a0 = s0 (-omegaE/omegaC + ks vT^2/omegaC^2 (A1 + (1 + bPlus) A2)) +
  ks vT^2/omegaC^2 A2 bCross s1;
a1 = -kAbs/omegaC s0;
a2 = ks A2/(2 omegaC^2) s0;

(* Equations (13.58) and (13.59), with their kernel prefactors rewritten
   so that 1/(8 Pi^2) is visibly common to every Fourier-space kernel. *)
rhoPhiMoment =
  plasmaZ[z0] (a0 + Sqrt[2] vT z0 a1 + 2 vT^2 z0^2 a2) +
  Sqrt[2] vT (a1 + Sqrt[2] vT z0 a2);
jPhiMoment =
  (z0 plasmaZ[z0] + 1) a0 + Sqrt[2] vT z0 a1 +
  2 vT^2 z0^2 a2 + 2 vT^2 a2;

rhoPhiCore = omegaC/(Sqrt[2] lambdaD^2 vT kAbs) rhoPhiMoment;
jPhiCore = omegaC/(lambdaD^2 kPole) jPhiMoment;

(* Equations (13.151) and (13.153), reduced to one common bracket. *)
magneticBracket =
  A2 s0/2 + (z0 plasmaZ[z0] + 1) ((A1 + A2 (1 + bPlus + z0^2)) s0 +
    A2 bCross s1);
rhoBCore = I vT^2/(c lambdaD^2 omegaC kPole) magneticBracket;
jBCore = I Sqrt[2] vT^3 z0/(c lambdaD^2 omegaC kPole) magneticBracket;

rhoPhiG = rhoPhiCore/(8 Pi^2);
rhoBG = rhoBCore/(8 Pi^2);
jPhiG = jPhiCore/(8 Pi^2);
jBG = jBCore/(8 Pi^2);

Print["magnetic recurrence check = ",
  FullSimplify[jBCore - Sqrt[2] vT z0 rhoBCore]];
Print["homogeneous static Debye check = ",
  FullSimplify[rhoPhiG /. {A1 -> 0, A2 -> 0, omegaE -> 0, omega -> 0,
      kr -> 0, krp -> 0, ks -> 0},
    Assumptions -> {vT > 0, omegaC > 0, lambdaD > 0, epsilon > 0,
      kParallel ∈ Reals}]];

oracleRules = {
  kr -> 0.73, krp -> -0.41, rhoL -> 0.62, ks -> 0.37,
  kParallel -> -0.23, epsilon -> 0.05,
  omegaE -> 0.41, omega -> 0.17,
  vT -> 1.3, omegaC -> 2.1, lambdaD -> 0.8,
  A1 -> 0.12, A2 -> -0.07, c -> 2.99792458*^10
};

Print["rho_phi = ", N[rhoPhiG /. oracleRules, 18]];
Print["rho_B   = ", N[rhoBG /. oracleRules, 18]];
Print["j_phi   = ", N[jPhiG /. oracleRules, 18]];
Print["j_B     = ", N[jBG /. oracleRules, 18]];
