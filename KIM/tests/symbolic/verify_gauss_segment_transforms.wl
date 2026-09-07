(* Exact audit of the standard-Gauss global-FEM F1/F2 cell transforms. *)

ClearAll["Global`*"];

assertTrue[value_, label_] := If[TrueQ[value],
  Print["PASS: " <> label],
  Print["FAIL: " <> label <> " -> " <> ToString[value]]; Exit[1]];

assumptions = rho > 0 && 0 < theta < Pi && lo < hi &&
  Element[{rho, theta, lo, hi, x, xp, ks, rg}, Reals];

c = Cos[theta];
s = Sin[theta];
mu = (x + xp)/2;
delta = x - xp;
kappa = ks^2 rho^2;
zlo = (mu - lo)/(rho Sqrt[1 + c]);
zhi = (mu - hi)/(rho Sqrt[1 + c]);
erfDifference = Erf[zlo] - Erf[zhi];
common = Exp[-kappa - delta^2/(4 rho^2 (1 - c))];
radialGaussian[r_] := Exp[-(r - mu)^2/(rho^2 (1 + c))];

q[r_, xv_, xpv_] := Exp[-ks^2 rho^2 -
  (r - (xv + xpv)/2)^2/(rho^2 (1 + c)) -
  (xv - xpv)^2/(4 rho^2 (1 - c))]/(rho^2 s);

j1 = Sqrt[Pi]/2 Sqrt[1 + c] erfDifference;
j23 = (1 + c)^(3/2) (
  Sqrt[Pi] (delta^2/(4 rho^2 (1 + c)) + 1/2) erfDifference
  + Exp[-zlo^2] (lo - mu)/(rho Sqrt[1 + c])
  - Exp[-zhi^2] (hi - mu)/(rho Sqrt[1 + c]));
j4 = 1/2 (1 + c)^(3/2) (
  Sqrt[Pi] (1/2 - delta^2/(4 rho^2 (1 + c))) erfDifference
  + Exp[-zlo^2] (lo - mu)/(rho Sqrt[1 + c])
  - Exp[-zhi^2] (hi - mu)/(rho Sqrt[1 + c]));

f1Direct = Integrate[q[rg, x, xp], {rg, lo, hi},
  Assumptions -> assumptions, GenerateConditions -> False];
f1Corrected = common j1/(rho s);
assertTrue[FullSimplify[f1Direct == f1Corrected, assumptions],
  "F1 corrected form equals the direct Gaussian cell integral"];

assertTrue[
  FullSimplify[D[rho j1, hi] == radialGaussian[hi], assumptions] &&
  FullSimplify[(rho j1 /. hi -> lo) == 0, assumptions],
  "Jrg1 is the zeroth Gaussian cell moment"];
assertTrue[
  FullSimplify[D[rho^3 j23, hi] ==
    ((hi - x)^2 + (hi - xp)^2) radialGaussian[hi], assumptions] &&
  FullSimplify[(rho^3 j23 /. hi -> lo) == 0, assumptions],
  "Jrg23 is the summed quadratic Gaussian cell moment"];
assertTrue[
  FullSimplify[D[rho^3 j4, hi] ==
    (hi - x) (hi - xp) radialGaussian[hi], assumptions] &&
  FullSimplify[(rho^3 j4 /. hi -> lo) == 0, assumptions],
  "Jrg4 is the mixed quadratic Gaussian cell moment"];

f2Density = rho^2 ks^2 q[rg, x, xp] -
  rho^2/2 (D[q[rg, x, xp], {x, 2}] +
    D[q[rg, x, xp], {xp, 2}]);

preRefactorCoefficient = 4 Cos[2 theta] (kappa + 1) -
  kappa Cos[4 theta] - (3 kappa + 4);
preRefactorDensity = -common radialGaussian[rg]/(8 rho^4 s^5) (
  rho^2 preRefactorCoefficient +
  (2 Cos[2 theta] + 6) ((rg - x)^2 + (rg - xp)^2) -
  16 c (rg - x) (rg - xp));
assertTrue[FullSimplify[
  TrigExpand[f2Density - preRefactorDensity] == 0, assumptions],
  "F2 pointwise density follows from the b+ differential operator"];

preRefactorIntegrated = -common/(8 rho^4 s^5) (
  rho^3 j1 preRefactorCoefficient +
  (2 Cos[2 theta] + 6) rho^3 j23 -
  16 c rho^3 j4);
f2Corrected = common/(2 rho s^5) (
  2 s^2 (1 + kappa s^2) j1 -
  (c^2 + 1) j23 + 4 c j4);
assertTrue[FullSimplify[
  TrigExpand[preRefactorIntegrated - f2Corrected] == 0, assumptions],
  "F2 corrected form equals the Gaussian-moment reduction"];

f1Stale = common Sqrt[Pi] (1 + c) erfDifference/(2 rho s);
f2Stale = common/(2 rho s^5) (
  (kappa + 2) s^2 j1 - (c^2 + 1) j23 + 4 c j4);
f1StaleRelativeError = FullSimplify[
  (f1Stale - f1Direct)/f1Direct, assumptions];
f2StaleResidual = FullSimplify[TrigExpand[
  (f2Stale - f2Corrected)/(common j1/(2 rho s^5))], assumptions];
assertTrue[f1StaleRelativeError === -1 + Sqrt[2] Cos[theta/2],
  "stale F1 mutation has the expected nonzero residual"];
assertTrue[f2StaleResidual === kappa Cos[2 theta] s^2,
  "stale F2 mutation has the expected nonzero residual"];

Print["GAUSS_SEGMENT_TRANSFORMS_SYMBOLIC_OK"];
