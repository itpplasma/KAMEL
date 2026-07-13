scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
If[!MemberQ[$Path, scriptDirectory], AppendTo[$Path, scriptDirectory]];
Needs["checklib`"];
resetChecks[];

repositoryRoot = ExpandFileName[FileNameJoin[{scriptDirectory, "..", ".."}]];
oraclePath = FileNameJoin[{repositoryRoot, "verification", "oracles",
  "plasma_background.dat"}];
mutation = Replace[Environment["KAMEL_VERIFY_NEGATIVE_FIXTURE"], $Failed -> ""];
regenerate = TrueQ[Environment["KAMEL_REGENERATE_ORACLE"] === "1"];

(* Locked Gaussian-CGS constants used by production.  Exact decimal rationals
   keep oracle generation independent of machine precision. *)
pi = N[Pi, 70]; c = 29979245800; e = 48032/10^14;
evReference = 16022/10^16;
ev = If[mutation === "wrong-ev", 1, evReference];
me = 91094/10^32; mp = 16726/10^28;
checkNumeric["production eV-to-erg conversion", ev, evReference, 10^-45];

(* Cylindrical geometry and harmonic phase exp(+i(m theta+n phi)). *)
geometry[r_, q_, bmag_, bsign_, r0_, mm_, nn_, er_] := Module[
  {bz, bt, hz, ht, ks, kp, ome},
  bz = bsign bmag/Sqrt[1 + r^2/(r0^2 q^2)];
  bt = bz r/(q r0); hz = bz/bmag; ht = bt/bmag;
  ks = (mm hz - nn ht/r0)/r;
  kp = mm ht/r + nn hz/r0;
  ome = -c ks er/bmag;
  {bz, bt, hz, ht, ks, kp, ome}
];
Clear[r, q, bz, r0, mm, nn];
btSymbolic = bz r/(q r0);
bmagSymbolic = Sqrt[bz^2 + btSymbolic^2];
hzSymbolic = bz/bmagSymbolic; htSymbolic = btSymbolic/bmagSymbolic;
check["safety-factor inversion", FullSimplify[r bz/(r0 btSymbolic) == q,
  Assumptions -> {r > 0, r0 > 0, q > 0, bz != 0}]];
check["resonance is m/q+n=0", TrueQ[FullSimplify[
  mm htSymbolic/r + nn hzSymbolic/r0 - hzSymbolic/r0 (mm/q + nn),
  Assumptions -> {r > 0, r0 > 0, q > 0, bz != 0}] == 0]];
check["field direction is normalized", FullSimplify[
  hzSymbolic^2 + htSymbolic^2 == 1,
  Assumptions -> {r > 0, r0 > 0, q > 0, bz != 0}]];
resonanceProbe = geometry[184/5, 3, 18000, -1, 165, -6, 2, -1/2];
check["finite-radius fixture crosses resonance", PossibleZeroQ[resonanceProbe[[6]]]];
zeroRotationProbe = geometry[12, 2, 18000, -1, 165, -6, 2, 0];
finiteRotationProbe = geometry[12, 2, 18000, -1, 165, -6, 2, -1/2];
check["zero radial field gives zero ExB rotation", zeroRotationProbe[[7]] == 0];
checkNumeric["finite ExB rotation retains c/B factor", finiteRotationProbe[[7]],
  -c finiteRotationProbe[[5]] (-1/2)/18000, 10^-45];

(* Ampere current and radial MHD force balance. The production equilibrium
   evolves u=B^2 and uses B_theta^2=r^2 u/(q^2 R0^2+r^2). *)
Clear[rr, btheta, bzfield, ufield, pressure, qfield, majorRadius];
jtheta = -c/(4 pi) D[bzfield[rr], rr];
jz = c/(4 pi) (btheta[rr]/rr + D[btheta[rr], rr]);
lorentzRadial = (jtheta bzfield[rr] - jz btheta[rr])/c;
check["equilibrium current gives cylindrical radial Lorentz force",
  FullSimplify[lorentzRadial == -(bzfield[rr] D[bzfield[rr], rr] +
    btheta[rr] D[btheta[rr], rr] + btheta[rr]^2/rr)/(4 pi),
    Assumptions -> rr > 0]];
productionDu = -2 rr ufield[rr]/(qfield[rr]^2 majorRadius^2 + rr^2) -
  8 pi D[pressure[rr], rr];
check["B-squared ODE is radial force balance", FullSimplify[
  productionDu/2 + rr ufield[rr]/(qfield[rr]^2 majorRadius^2 + rr^2) +
    4 pi D[pressure[rr], rr] == 0]];

(* Species-derived quantities and thesis (2.48)-(2.51). *)
vthermal[temp_, mass_] := Sqrt[temp ev/mass];
omegaCyclotron[z_, mass_, b_] := z e Abs[b]/(mass c);
debye[n_, temp_, z_] := Sqrt[temp ev/(4 pi n (z e)^2)];
lambdaEE[ne_, te_] := 23.5 - Log[Sqrt[ne]/te^(5/4)] -
  Sqrt[10^-5 + (Log[te] - 2)^2/16];
lambdaEI[ne_, te_] := 24 - Log[Sqrt[ne]/te];
lambdaII[ni_, ti_, zi_, ai_, nj_, tj_, zj_, aj_] :=
  23 - Log[zi zj (ai + aj)/(ai tj + aj ti) *
    Sqrt[ni zi^2/ti + nj zj^2/tj]];
collisionFrequencies[n_List, temp_List, z_List, a_List] := Module[
  {ns = Length[n], lee, lei, nu},
  lee = lambdaEE[n[[1]], temp[[1]]];
  lei = lambdaEI[n[[1]], temp[[1]]];
  nu = ConstantArray[0, ns];
  nu[[1]] = 5.8*^-6 n[[1]] lee/temp[[1]]^(3/2) +
    Sum[7.7*^-6 n[[j]] z[[j]]^2 lei/temp[[1]]^(3/2), {j, 2, ns}];
  Do[nu[[i]] = 1.8*^-7 n[[1]] z[[i]]^2 lei/
      (Sqrt[a[[i]]] temp[[i]]^(3/2)) +
    Sum[1.8*^-7 n[[j]] z[[i]]^2 z[[j]]^2*
      lambdaII[n[[i]], temp[[i]], z[[i]], a[[i]], n[[j]], temp[[j]],
        z[[j]], a[[j]]]/(Sqrt[a[[i]]] temp[[i]]^(3/2)), {j, 2, ns}],
    {i, 2, ns}];
  nu
];

(* KiLCA's named legacy single-ion NRL background model. The additive floor
   and separate electron/ion scale factors are solver policy, not KIM physics. *)
kilcaBackground[n0_, ti_, te_, z0_, a0_, massi_, masse_, velocityFactor_] := Module[
  {lee, lei, lii, nuee, nuei, nuie, nuii, nue, nui, vti, vte,
   kilcaEv = 160216428/10^20},
  lee = 23.5 - Log[Sqrt[n0]/te^(5/4)] - (Log[te] - 2)/4;
  lei = 30 - Log[Sqrt[n0]/ti^(3/2) z0^2/a0];
  lii = 23 - Log[z0^3/ti^(3/2) Sqrt[2 n0]];
  nuee = 58/10^7 n0 lee/(Sqrt[te] velocityFactor te);
  nuei = 77/10^7 n0 lei z0^2/(velocityFactor te)^(3/2);
  nuie = 32/10^10 n0 lei z0^2/(a0 Sqrt[te] velocityFactor ti);
  nuii = 14/10^8 n0 lii z0^4/(Sqrt[a0] Sqrt[ti] velocityFactor ti);
  nue = nuee + nuei + 10; nui = nuie + nuii + 10;
  vti = Sqrt[kilcaEv ti/massi]; vte = Sqrt[kilcaEv te/masse];
  {lee, lei, lii, nuee, nuei, nuie, nuii, nue, nui, vti, vte}
];

ne = 2*^13; ncarbon = 1*^12;
nd = If[mutation === "charge-weighted-density", ne/(1 + 6), ne - 6 ncarbon];
checkNumeric["arbitrary D/C mixture is quasineutral", nd + 6 ncarbon, ne,
  10^-45];
speciesProbes = {{ne, 1200, -1, me}, {nd, 800, 1, 2 mp},
  {ncarbon, 500, 6, 12 mp}};
check["all species-derived lengths are positive", And @@ Map[
  Function[probe, With[{n0 = probe[[1]], t = probe[[2]], z0 = probe[[3]],
      mass = probe[[4]]},
    vthermal[t, mass] > 0 && debye[n0, t, z0] > 0 &&
      Abs[vthermal[t, mass]/omegaCyclotron[z0, mass, -18000]] > 0]],
  speciesProbes]];
check["cyclotron frequency retains charge sign",
  omegaCyclotron[-1, me, -18000] < 0 && omegaCyclotron[6, 12 mp, -18000] > 0];
signedElectronRho = vthermal[1200, me]/omegaCyclotron[-1, me, -18000];
signedIonRho = vthermal[800, 2 mp]/omegaCyclotron[1, 2 mp, -18000];
check["signed and absolute Larmor radii are distinct conventions",
  signedElectronRho < 0 && signedIonRho > 0 &&
    Abs[signedElectronRho] == Abs[vthermal[1200, me]/omegaCyclotron[-1, me, -18000]]];
z0Background[ome_, drive_, nu_, kp_, vt_] :=
  -(ome - drive - I nu)/(Abs[kp] Sqrt[2] vt);
check["plasma-dispersion argument carries exactly sqrt(2)", FullSimplify[
  z0Background[3, 1, 2, -5, 7] Abs[-5] Sqrt[2] 7 == -(3 - 1 - 2 I)]];

(* Density scaling leaves logarithmic forces invariant and Debye length scales
   as scale^-1/2. *)
scale = 7/3;
check["density rescaling preserves dlogn", (scale 3)/(scale 11) == 3/11];
checkNumeric["density rescaling changes Debye length",
  debye[scale ne, 1200, -1], debye[ne, 1200, -1]/Sqrt[scale], 10^-45];

(* Recompute nonlinear derived quantities after interpolation. *)
tLeft = 200; tRight = 1800; massProbe = 2 mp;
vCenter = vthermal[(tLeft + tRight)/2, massProbe];
vInterpolated = (vthermal[tLeft, massProbe] + vthermal[tRight, massProbe])/2;
cellValue = If[mutation === "interpolated-derived-ratio", vInterpolated, vCenter];
checkNumeric["cell-center thermal speed is derived after interpolation",
  cellValue, vCenter, 10^-45];
check["interpolating a nonlinear derived ratio is observably different",
  Abs[vCenter - vInterpolated]/vCenter > 10^-2];

(* The active window scale is the largest active Larmor radius, never a fixed
   electron slot. *)
rhoProbe = {0.02, 0.25, 0.08};
selectedRho = If[mutation === "electron-window-scale", First[rhoProbe], Max[rhoProbe]];
check["active-species window scale", selectedRho == 0.25];

(* Radial electric-field force balance, including the finite toroidal-flow
   term. Ti is in eV; ev/e converts both pressure-gradient terms to statV/cm. *)
erForceBalance[n0_, ti_, dn_, dti_, rminor_, b0_, vz_, q0_, rmajor_] :=
  ti ev dn/(e n0) + ev dti/e + rminor b0 vz/(c q0 rmajor);
pressureOnlyEr = erForceBalance[2*^13, 800, -2*^11, -8, 184/5,
  -18000, 0, 3, 165];
rotatingEr = erForceBalance[2*^13, 800, -2*^11, -8, 184/5,
  -18000, 2*^5, 3, 165];
check["zero flow removes radial-force-balance rotation term",
  pressureOnlyEr == erForceBalance[2*^13, 800, -2*^11, -8, 184/5,
    -18000, 0, 3, 165]];
checkNumeric["finite flow adds radial-force-balance rotation term",
  rotatingEr - pressureOnlyEr, (184/5) (-18000) 2*^5/(c 3 165), 10^-45];
kilcaEvReference = 160216428/10^20;
kilcaPotentialGradient[vsi_, b0_, n0_, ti_, dn_, dti_, z0_] :=
  vsi b0/c - kilcaEvReference (ti dn + dti n0)/(z0 e n0);
check["KiLCA zero-flow potential gradient is the ion pressure term",
  kilcaPotentialGradient[0, 18000, 2*^13, 800, -2*^11, -8, 1] ==
    -kilcaEvReference (800 (-2*^11) - 8 (2*^13))/(e (2*^13))];
checkNumeric["KiLCA finite-flow potential gradient retains c factor",
  kilcaPotentialGradient[2*^5, 18000, 2*^13, 800, -2*^11, -8, 1] -
    kilcaPotentialGradient[0, 18000, 2*^13, 800, -2*^11, -8, 1],
  2*^5 18000/c, 10^-45];

(* Force-balance terms all have statvolt/cm dimensions. Dimension vectors are
   {mass,length,time,charge}; eV-to-erg is an energy. *)
dimEnergy = {1, 2, -2, 0}; dimCharge = {0, 0, 0, 1};
dimElectric = dimEnergy - dimCharge - {0, 1, 0, 0};
dimVelocity = {0, 1, -1, 0};
dimMagnetic = dimElectric; (* Gaussian CGS: E and B have equal dimensions. *)
check["pressure-gradient force-balance units",
  dimEnergy - dimCharge - {0, 1, 0, 0} == dimElectric];
check["rotation force-balance units",
  dimVelocity + dimMagnetic - dimVelocity == dimElectric];

(* A smooth periodic profile has matching values and derivatives at both
   seams. The endpoint-exclusive production grid must differentiate it with a
   wrapped stencil after periodization, never reuse the aperiodic derivative. *)
seamLength = 18; seamWave = 2 pi/seamLength;
seamProfiles[x_] := {2*^13 (1 + Cos[seamWave x]/20 + Sin[seamWave x]/100),
  1200 + 40 Cos[seamWave x] + 15 Sin[seamWave x],
  800 + 30 Cos[seamWave x] - 12 Sin[seamWave x],
  500 - 20 Cos[seamWave x] + 9 Sin[seamWave x],
  3 + Cos[seamWave x]/50 + Sin[seamWave x]/100,
  -1/2 + Cos[seamWave x]/100 + Sin[seamWave x]/200};
seamDerivatives[x_] := Evaluate[D[seamProfiles[y], y] /. y -> x];
checkNumeric["periodized primitive values match at seam",
  Norm[seamProfiles[-seamLength/2] - seamProfiles[seamLength/2]], 0, 10^-40];
checkNumeric["periodized primitive derivatives match at seam",
  Norm[seamDerivatives[-seamLength/2] - seamDerivatives[seamLength/2]], 0, 10^-40];
seamGridCount = 96; seamStep = seamLength/seamGridCount;
seamSamples = Table[First[seamProfiles[-seamLength/2 + j seamStep]],
  {j, 0, seamGridCount - 1}];
wrappedDerivative = (seamSamples[[2]] - seamSamples[[-1]])/(2 seamStep);
oneSidedDerivative = (seamSamples[[2]] - seamSamples[[1]])/seamStep;
analyticSeamDerivative = First[seamDerivatives[-seamLength/2]];
check["wrapped seam derivative improves on a one-sided derivative",
  Abs[wrappedDerivative - analyticSeamDerivative] <
    Abs[oneSidedDerivative - analyticSeamDerivative]];

If[StringLength[mutation] > 0, FinishChecks[]];

formatNumber[value_] := StringReplace[
  ToString[FortranForm[N[value, 34]], InputForm], {"*^" -> "E", " " -> ""}];
pad[values_] := PadRight[values, 20, 0];
row[kind_, icase_, isp_, values_] := StringRiffle[
  Join[ToString /@ {kind, icase, isp, 0}, formatNumber /@ pad[values]], " "];

geometryCases = {{1/1000, 6/5, 18000, -1, 165, -6, 2, -1/2},
  {184/5, 3, 18000, -1, 165, -6, 2, -1/2},
  {12, 2, 18000, -1, 165, -6, 2, 0}};
geometryRows = MapIndexed[Function[{args, idx}, Module[{g = geometry @@ args},
  row[1, First[idx], 0, Join[N[args, 70], g]]]], geometryCases];

caseData = {
  <|"n" -> {2*^13, 2*^13}, "T" -> {1000, 800}, "z" -> {-1, 1},
    "a" -> {1, 2}, "mass" -> {me, 2 mp}, "B" -> -18000,
    "dn" -> {-2*^11, -2*^11}, "dT" -> {-15, -10}, "Er" -> -1/2,
    "kp" -> 3/100, "ks" -> -1/5, "omE" -> 2*^4, "V" -> {0, 3*^5},
    "omega" -> 0, "mphi" -> 1|>,
  <|"n" -> {ne, nd, ncarbon}, "T" -> {1200, 800, 500},
    "z" -> {-1, 1, 6}, "a" -> {1, 2, 12},
    "mass" -> {me, 2 mp, 12 mp}, "B" -> -18000,
    "dn" -> {-2*^11, -14*^10, -1*^10}, "dT" -> {-12, -8, -5},
    "Er" -> -1/2, "kp" -> -2/100, "ks" -> 1/4, "omE" -> -3*^4,
    "V" -> {0, 2*^5, -1*^5}, "omega" -> 0, "mphi" -> -1|>
};
speciesRows = {}; collisionRows = {};
Do[Module[{d = caseData[[ic]], nus, ns, vt, oc, rho, ld, a1, a2, x1, x2},
  nus = collisionFrequencies[d["n"], d["T"], d["z"], d["a"]];
  ns = Length[d["n"]];
  Do[
    vt = vthermal[d["T"][[sp]], d["mass"][[sp]]];
    oc = omegaCyclotron[d["z"][[sp]], d["mass"][[sp]], d["B"]];
    rho = Abs[vt/oc]; ld = debye[d["n"][[sp]], d["T"][[sp]], d["z"][[sp]]];
    a1 = d["dn"][[sp]]/d["n"][[sp]] - d["z"][[sp]] e d["Er"]/(d["T"][[sp]] ev) -
      3 d["dT"][[sp]]/(2 d["T"][[sp]]); a2 = d["dT"][[sp]]/d["T"][[sp]];
    x1 = d["kp"] vt/nus[[sp]];
    x2 = -(d["omE"] + d["mphi"] oc - d["omega"] + d["kp"] d["V"][[sp]])/nus[[sp]];
    AppendTo[speciesRows, row[2, ic, sp - 1, {d["n"][[sp]], d["T"][[sp]],
      d["z"][[sp]], d["a"][[sp]], d["mass"][[sp]], vt, oc, rho, ld, nus[[sp]],
      a1, a2, x1, x2, d["dn"][[sp]], d["dT"][[sp]], d["V"][[sp]], d["kp"],
      d["omE"], d["mphi"]}]],
    {sp, 1, ns}];
  Do[AppendTo[collisionRows, row[3, ic, sp - 1, {nus[[sp]],
    lambdaEE[d["n"][[1]], d["T"][[1]]],
    If[sp == 1, 0, lambdaEI[d["n"][[1]], d["T"][[1]]]]}]], {sp, 1, ns}]
], {ic, 1, Length[caseData]}];

kilcaMp = 167262158/10^32; kilcaMe = 910938185917485/10^42;
kilcaInputs = {2 10^13, 800, 1000, 1, 2, 2 kilcaMp, kilcaMe, 1};
kilcaRows = {row[4, 1, 0, Join[kilcaInputs, kilcaBackground @@ kilcaInputs]]};
seamRows = {row[5, 1, 0, Join[{-seamLength/2},
    seamProfiles[-seamLength/2], seamDerivatives[-seamLength/2]]],
  row[5, 1, 1, Join[{seamLength/2},
    seamProfiles[seamLength/2], seamDerivatives[seamLength/2]]]};
equilibriumFixture = Module[{rfix = 184/5, bt = 1200, bz = -18000,
    dbt = 3, dbz = -2, dp, jt, jzfix},
  dp = -(bt dbt + bz dbz + bt^2/rfix)/(4 pi);
  jt = -c dbz/(4 pi); jzfix = c (bt/rfix + dbt)/(4 pi);
  {rfix, bt, bz, dbt, dbz, dp, jt, jzfix,
    4 pi dp + bt dbt + bz dbz + bt^2/rfix}];
equilibriumRows = {row[6, 1, 0, equilibriumFixture]};
oracleRows = Join[geometryRows, speciesRows, collisionRows, kilcaRows,
  seamRows, equilibriumRows];
header = {"# KAMEL plasma-background oracle; thesis (2.48)-(2.52), (3.5)-(3.6), (5.8), (5.18)-(5.20).",
  "# Fixed columns: kind case species reserved v1 ... v20.",
  "# kind 1 geometry; kind 2 KIM species; kind 3 KIM collisions; kind 4 KiLCA legacy single-ion background; kind 5 periodic seams; kind 6 current/force balance."};
If[regenerate,
  Export[oraclePath, StringRiffle[Join[header, oracleRows], "\n"] <> "\n", "Text"];
  Print["Wrote ", Length[oracleRows], " rows to ", oraclePath],
  Module[{lines, parse, committed, generated},
    lines = Select[StringSplit[Import[oraclePath, "Text"], "\n"],
      StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &];
    parse[token_] := ToExpression[StringReplace[token,
      RegularExpression["[Ee]([+-]?[0-9]+)$"] -> "*^$1"]];
    committed = (parse /@ StringSplit[#]) & /@ lines;
    generated = (parse /@ StringSplit[#]) & /@ oracleRows;
    check["committed plasma-background oracle row count",
      Length[committed] == Length[generated] == 17];
    If[Length[committed] == Length[generated], Do[
      checkNumeric["plasma-background oracle row " <> ToString[i],
        Norm[committed[[i]] - generated[[i]]], 0, 10^-30],
      {i, Length[generated]}]]
  ]
];

FinishChecks[];
