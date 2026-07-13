scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
If[!MemberQ[$Path, scriptDirectory], AppendTo[$Path, scriptDirectory]];
Needs["checklib`"];
resetChecks[];

repositoryRoot = ExpandFileName[FileNameJoin[{scriptDirectory, "..", ".."}]];
oraclePath = FileNameJoin[{repositoryRoot, "verification", "oracles",
  "field_operator_matrices.dat"}];
workingPrecision = 70;
$MaxExtraPrecision = 5000;
mutation = Replace[Environment["KAMEL_VERIFY_NEGATIVE_FIXTURE"], $Failed -> ""];
regenerate = TrueQ[Environment["KAMEL_REGENERATE_ORACLE"] === "1"];

(* The full cylindrical harmonic operator precedes the local one-dimensional
   reduction used by the chapter-15 integral-model implementation. *)
Clear[r, theta, z, mm, kz, radialField];
cylindricalLaplacian[expr_] := D[expr, {r, 2}] + D[expr, r]/r +
  D[expr, {theta, 2}]/r^2 + D[expr, {z, 2}];
harmonicField = radialField[r] Exp[I (mm theta + kz z)];
harmonicLaplacian = FullSimplify[
  cylindricalLaplacian[harmonicField]/Exp[I (mm theta + kz z)]];
check["cylindrical harmonic Laplacian", TrueQ[FullSimplify[
  harmonicLaplacian - (radialField''[r] + radialField'[r]/r -
    (mm^2/r^2 + kz^2) radialField[r])] == 0]];

fullWKBLabel[kr_, ks_, kpar_] := -(kr^2 + ks^2 + kpar^2);
radialSymbol[kr_, ks_] := If[mutation === "hidden-ks2", -(kr^2 + ks^2), -kr^2];
check["full cylindrical WKB symbol retains helical terms",
  fullWKBLabel[3/5, 2/7, 1/4] == -((3/5)^2 + (2/7)^2 + (1/4)^2)];
check["chapter-15 reduced KIM operator is radial only",
  radialSymbol[3/5, 2/7] == -(3/5)^2];

(* Radial Fourier pair: f(r)=Sum[f_m exp(+i k_m r)] and
   f_m=(1/L) Integral[f exp(-i k_m r)]. *)
kOfM[m_, length_] := 2 Pi m/length;
transformFactor[n_] := If[mutation === "transform-factor", 2/n, 1/n];
periodicCoefficient[values_, nodes_, m_, length_] :=
  transformFactor[Length[values]] Total[
    values Exp[-I kOfM[m, length] nodes]];
checkNumeric["inverse radial Fourier phase",
  D[Exp[I k r], {r, 2}], -k^2 Exp[I k r], 10^-60];

driveL = 4;
driveNodes = Range[0, 3];
driveValues = 2 + Exp[I 2 Pi driveNodes/driveL];
driveCoefficients = Association@Table[
  m -> periodicCoefficient[driveValues, driveNodes, m, driveL], {m, -2, 2}];
checkNumeric["constant magnetic drive coefficient",
  driveCoefficients[0], 2, 10^-60];
checkNumeric["nonconstant magnetic drive coefficient",
  driveCoefficients[1], 1, 10^-60];
check["unexcited magnetic drive modes vanish",
  And @@ (Abs[driveCoefficients[#]] < 10^-60 & /@ {-2, -1, 2})];

(* Thesis (15.26)-(15.30): variable-grid hat overlap and weak Laplacian. *)
hatNodes = {0, 1, 3};
hatMassInterior = {1/6, 1, 1/3};
hatLaplaceInterior = {1, -3/2, 1/2};
check["variable-grid hat overlap row",
  hatMassInterior == {1/6, (1 + 2)/3, 2/6}];
check["variable-grid weak Laplace row",
  hatLaplaceInterior == {1/1, -(1/1 + 1/2), 1/2}];

(* Green's identity requires its boundary contribution when the test function
   does not vanish. u=x and v=1-x expose a missing term exactly. *)
uProbe[x_] := x;
vProbe[x_] := 1 - x;
strongProbe = Integrate[vProbe[x] uProbe''[x], {x, 0, 1}];
weakProbe = -Integrate[vProbe'[x] uProbe'[x], {x, 0, 1}];
boundaryProbe = If[mutation === "boundary-term", 0,
  vProbe[1] uProbe'[1] - vProbe[0] uProbe'[0]];
check["weak form includes boundary term", strongProbe == weakProbe + boundaryProbe];

(* Poisson convention Delta Phi=-4 pi rho and non-symmetric kernel action. *)
sourceSign = If[mutation === "source-sign", 1, -1];
magneticSource[kb_, br_] := sourceSign 4 Pi kb.br;
checkNumeric["magnetic Poisson source sign",
  Norm[magneticSource[IdentityMatrix[2], {2, -3}] - (-4 Pi {2, -3})],
  0, 10^-60];

(* With the locked exp(+i k_parallel s) phase, the electrostatic parallel
   field is -i k_parallel Phi.  The ideal edge lifting cancels the equilibrium
   perpendicular-field projection along the perturbed magnetic field. *)
alignedPotential[eperp_, br_, b0_, kpar_] := I eperp br/(b0 kpar);
parallelIdealResidual[phi_, eperp_, br_, b0_, kpar_] :=
  -I kpar phi - eperp br/b0;
checkNumeric["aligned-potential phase and sign",
  parallelIdealResidual[alignedPotential[2 - I, 3 + I/2, 5, 7/4],
    2 - I, 3 + I/2, 5, 7/4], 0, 10^-60];

kernelProbe = {{1 + I, 2 - I}, {3 + 2 I, 4 - I}};
phiProbe = {2 - I, -1 + 3 I};
kernelAction = If[mutation === "transposed-kernel",
  Transpose[kernelProbe].phiProbe, kernelProbe.phiProbe];
check["kernel row is outgoing/test index", kernelAction == kernelProbe.phiProbe];

(* Screening fixes mode zero. Only the vacuum operator has a constant null
   space; deleting the screened zero mode is an unjustified gauge. *)
screenedZero[kappa_] := If[mutation === "unjustified-gauge", 0, -kappa^2];
check["screened mode zero is physical and nonsingular", screenedZero[3/4] == -9/16];
check["vacuum mode zero requires a declared gauge", radialSymbol[0, 0] == 0];

(* A same-space shift is algebraically a no-op. Phi_ref must remain a
   physical-space lifting field during reconstruction. *)
aProbe = {{3 + I, 1/5}, {-2/7 I, 2 - I}};
gProbe = {1 + 2 I, -3 + I};
pProbe = {2/3 - I/4, -1/5 + 3 I/7};
checkNumeric["same-basis deviation shift is a no-op",
  Norm[LinearSolve[aProbe, gProbe - aProbe.pProbe] + pProbe -
    LinearSolve[aProbe, gProbe]], 0, 10^-55];
finiteFourier[rp_] := Sum[(1/(1 + m^2)) Exp[I kOfM[m, 4] rp], {m, -3, 3}];
checkNumeric["finite Fourier reconstruction is periodic",
  finiteFourier[1/3], finiteFourier[1/3 + 4], 10^-55];
oneSidedAmplitude = 2/3 - 4 I/5;
realFieldFromPair[phase_] := (oneSidedAmplitude Exp[I phase] +
  Conjugate[oneSidedAmplitude] Exp[-I phase])/2;
checkNumeric["one-sided complex amplitude reconstructs a real field",
  realFieldFromPair[11/13], Re[oneSidedAmplitude Exp[I 11/13]], 10^-55];

If[StringLength[mutation] > 0, FinishChecks[]];

(* Independent FP mixed-moment check for current kernels. #197 depends on the
   arbitrary-precision #194 susceptibility oracle, so consume one complete
   verified moment set instead of duplicating that expensive integration. The
   magnetic drive supplies one parallel velocity and parallel current supplies
   the other: j-B therefore requires I11/I13. The printed I02/I22 in thesis
   (14.6) is not the mixed response moment. *)
susceptibilityPath = FileNameJoin[{repositoryRoot, "verification", "oracles",
  "fp_susceptibilities.dat"}];
parseNumber[token_] := ToExpression[StringReplace[token,
  RegularExpression["[Ee]([+-]?[0-9]+)$"] -> "*^$1"]];
susceptibilityLines = Select[StringSplit[Import[susceptibilityPath, "Text"], "\n"],
  StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &];
susceptibilityRows = (parseNumber /@ StringSplit[#]) & /@ susceptibilityLines;
currentRecords = Select[susceptibilityRows, #[[1]] == 1 && #[[8]] == -1 &];
moment[k_, l_] := FirstCase[currentRecords,
  row_ /; row[[11]] == k && row[[12]] == l :> row[[15]] + I row[[16]]];
currentX1 = currentRecords[[1, 13]];
currentX2 = currentRecords[[1, 14]];
i01 = moment[0, 1]; i10 = i01;
i21 = moment[2, 1]; i12 = i21;
i11 = moment[1, 1]; i13 = moment[1, 3];
i02 = moment[0, 2]; i22 = moment[2, 2];
fpGeneratingSymmetric[al_, be_, x1_, x2_, tau_] :=
  Exp[(I x2 - x1^2) tau +
    (al + I x1) (be + I x1) (Exp[-tau] - 1) + (al + be)^2/2];
check["FP generator is symmetric in incoming/outgoing moment variables",
  fpGeneratingSymmetric[al, be, x1, x2, tau] ==
    fpGeneratingSymmetric[be, al, x1, x2, tau]];
checkNumeric["j-Phi transpose-equivalent moment I10=I01", i10, i01, 10^-36];
checkNumeric["j-Phi energy moment I12=I21", i12, i21, 10^-36];
check["j-B mixed moment differs from thesis I02", Abs[i11 - i02] > 10^-6];

a1Current = 2/7; a2Current = -1/5;
f0Current = 4/5; f2Current = 3/10;
jBCorrect = i11 (a1Current f0Current + a2Current f2Current) +
  a2Current i13 f0Current/2;
jBPrinted = i02 (a1Current f0Current + a2Current f2Current) +
  a2Current i22 f0Current/2;
check["j-B production mixed moments reject printed I02/I22",
  Abs[jBCorrect - jBPrinted] > 10^-6];

(* Deterministic matrix/field oracle. All row kinds have 20 value columns:
   kind 1 screened Fourier mode; 2 hat row; 3 magnetic-drive coefficient;
   4 current-moment discrepancy. *)
formatNumber[value_] := StringReplace[
  ToString[FortranForm[N[value, 34]], InputForm], {"*^" -> "E", " " -> ""}];
padValues[values_] := PadRight[values, 20, 0];
rowString[kind_, case_, row_, col_, values_] := StringRiffle[
  Join[ToString /@ {kind, case, row, col}, formatNumber /@ padValues[values]], " "];

fourierL = 5;
fourierKappa = 3/4;
fourierRows = Table[Module[{kk, aa, pref, psi, total, source, bpsi, er},
  kk = kOfM[m, fourierL];
  aa = -kk^2 - fourierKappa^2;
  pref = (1/5 + I/7)/(1 + m^2);
  psi = (m + 3)/20 + I (2 - m)/30;
  total = pref + psi;
  source = aa total;
  bpsi = source - aa pref;
  er = -I kk total;
  rowString[1, 1, m, m, {fourierL, fourierKappa, kk, Re[aa], Im[aa],
    Re[pref], Im[pref], Re[psi], Im[psi], Re[source], Im[source],
    Re[bpsi], Im[bpsi], Re[total], Im[total], Re[er], Im[er],
    2/7, 1/4, fullWKBLabel[kk, 2/7, 1/4]}]
], {m, -2, 2}];

hatRows = Table[rowString[2, 1, 2, j, {hatNodes[[2]], hatNodes[[j]],
  hatMassInterior[[j]], hatLaplaceInterior[[j]], boundaryProbe}], {j, 1, 3}];
driveRows = Table[rowString[3, 1, m, 0, {driveL,
  Re[driveCoefficients[m]], Im[driveCoefficients[m]]}], {m, -2, 2}];
currentRows = {rowString[4, 1, 1, 1, {currentX1, currentX2,
  Re[i11], Im[i11], Re[i02], Im[i02], Re[i13], Im[i13], Re[i22], Im[i22],
  Re[jBCorrect], Im[jBCorrect], Re[jBPrinted], Im[jBPrinted]}]};
oracleRows = Join[fourierRows, hatRows, driveRows, currentRows];

oracleHeader = {
  "# KAMEL field/operator oracle; Markl thesis (14.4)-(14.6), (15.10)-(15.31).",
  "# Fixed columns: kind case row col v1 ... v20. Kind meanings are documented in the Mathematica source.",
  "# Kind 1 includes screened mode zero, physical-space reference, deviation RHS, total field, radial E, and full WKB label.",
  "# Kind 2 records the variable-grid hat mass/Laplace row; kind 3 drive DFT; kind 4 proves j-B uses mixed I11/I13."
};

If[regenerate,
  Export[oraclePath, StringRiffle[Join[oracleHeader, oracleRows], "\n"] <> "\n", "Text"];
  Print["Wrote ", Length[oracleRows], " rows to ", oraclePath],
  Module[{lines, parseNumber, committed, generated},
    lines = Select[StringSplit[Import[oraclePath, "Text"], "\n"],
      StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &];
    parseNumber[token_] := ToExpression[StringReplace[token,
      RegularExpression["[Ee]([+-]?[0-9]+)$"] -> "*^$1"]];
    committed = (parseNumber /@ StringSplit[#]) & /@ lines;
    generated = (parseNumber /@ StringSplit[#]) & /@ oracleRows;
    check["14 committed field/operator oracle rows",
      Length[committed] == Length[generated] == 14];
    If[Length[committed] == Length[generated],
      Do[checkNumeric["field/operator oracle row " <> ToString[i],
        Norm[committed[[i]] - generated[[i]]], 0, 10^-30],
        {i, 1, Length[generated]}]]
  ]
];

FinishChecks[];
