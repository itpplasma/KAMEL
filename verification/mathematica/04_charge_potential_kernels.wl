scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
If[!MemberQ[$Path, scriptDirectory], AppendTo[$Path, scriptDirectory]];
Needs["checklib`"];
resetChecks[];

repositoryRoot = ExpandFileName[FileNameJoin[{scriptDirectory, "..", ".."}]];
oraclePath = FileNameJoin[{repositoryRoot, "verification", "oracles", "rho_phi_kernels.dat"}];
workingPrecision = 70;
$MaxExtraPrecision = 5000;
mutation = Replace[Environment["KAMEL_VERIFY_NEGATIVE_FIXTURE"], $Failed -> ""];
regenerate = TrueQ[Environment["KAMEL_REGENERATE_ORACLE"] === "1"];

(* Fourier convention fixed by thesis (13.3), (15.12). *)
inverseMeasure = If[mutation === "fourier-measure", 1/Sqrt[2 Pi], 1/(2 Pi)];
poissonMeasure = If[mutation === "four-pi", 1, 1/(4 Pi)];
kernelNormalization = poissonMeasure inverseMeasure/(1/(2 Pi));
check["radial inverse-Fourier measure", inverseMeasure == 1/(2 Pi)];
check["Gaussian-CGS kernel normalization", poissonMeasure == 1/(4 Pi)];

kPerp[ks_, kr_] := Sqrt[ks^2 + kr^2];
alpha[ks_, kr_] := ArcTan[ks, kr];
bPlus[ks_, kr_, krp_, rho_] := If[
  mutation === "hidden-ks2",
  rho^2 (kr^2 + krp^2)/2,
  rho^2 (2 ks^2 + kr^2 + krp^2)/2
];
bCross[ks_, kr_, krp_, rho_] := If[
  mutation === "hidden-ks2",
  rho^2 Abs[kr krp],
  rho^2 kPerp[ks, kr] kPerp[ks, krp]
];
fourierGyroPhase[ks_, kr_, krp_, rg_, m_Integer] := Module[{sign = 1},
  If[mutation === "phase", sign = -1];
  Exp[I sign ((kr - krp) rg + m (alpha[ks, krp] - alpha[ks, kr]))]
];

checkNumeric[
  "finite-radius b+ retains ks^2",
  bPlus[2/5, -3/7, 5/6, 4/3],
  (4/3)^2 (2 (2/5)^2 + (-3/7)^2 + (5/6)^2)/2,
  10^-60
];
checkNumeric[
  "finite-radius bx",
  bCross[2/5, -3/7, 5/6, 4/3],
  (4/3)^2 Sqrt[(2/5)^2 + (-3/7)^2] Sqrt[(2/5)^2 + (5/6)^2],
  10^-60
];
checkNumeric[
  "Fourier and gyro phase",
  fourierGyroPhase[2/5, -3/7, 5/6, 7/10, 2],
  Exp[I (((-3/7) - 5/6) 7/10 +
    2 (ArcTan[2/5, 5/6] - ArcTan[2/5, -3/7]))],
  10^-60
];

(* Upstream perpendicular-velocity/gyro integral. The Weber integral is the
   Maxwellian average of the two independently gyro-averaged plane waves. *)
weberGenerating[s_, a_, c_, m_] :=
  Exp[-(a^2 + c^2)/(4 s)] BesselI[m, a c/(2 s)]/(2 s);
flrDensity[bp_, bx_, m_Integer] := Exp[-bp] BesselI[
  If[mutation === "bessel-index", m + 1, m], bx
];
flrEnergy[bp_, bx_, m_Integer] := Exp[-bp] (
  (1 + If[mutation === "bplus-sign", bp, -bp] - m) BesselI[m, bx] +
  bx BesselI[m - 1, bx]
);

symbolicEnergy = FullSimplify[
  -D[weberGenerating[s, a, c, m], s]/2 /. s -> 1/2,
  Assumptions -> {a > 0, c > 0, Element[m, Integers]}
];
derivedEnergy = Exp[-(a^2 + c^2)/2] (
  (1 - (a^2 + c^2)/2 - m) BesselI[m, a c] + a c BesselI[m - 1, a c]
);
check[
  "upstream Gaussian derivative fixes 1-b+-m",
  TrueQ[FullSimplify[symbolicEnergy == derivedEnergy,
    Assumptions -> {a > 0, c > 0, Element[m, Integers]}]]
];

directMoment[power_Integer, m_Integer, ain_, cin_] := NIntegrate[
  u^(1 + 2 power)/2^power Exp[-u^2/2] BesselJ[m, ain u] BesselJ[m, cin u],
  {u, 0, Infinity}, WorkingPrecision -> 55, AccuracyGoal -> 38,
  PrecisionGoal -> 34, MaxRecursion -> 30
];
directGyroBessel[m_Integer, x_] := NIntegrate[
  Exp[I (m tau - x Sin[tau])]/(2 Pi), {tau, -Pi, Pi},
  WorkingPrecision -> 55, AccuracyGoal -> 38, PrecisionGoal -> 34
];

probeA = SetPrecision[7/10, workingPrecision];
probeC = SetPrecision[11/10, workingPrecision];
probeBp = (probeA^2 + probeC^2)/2;
probeBx = probeA probeC;
checkNumeric[
  "density Bessel product from velocity integral",
  directMoment[0, 1, probeA, probeC], flrDensity[probeBp, probeBx, 1], 10^-32
];
checkNumeric[
  "energy Bessel product proves minus b+ sign",
  directMoment[1, 1, probeA, probeC], flrEnergy[probeBp, probeBx, 1], 10^-32
];
checkNumeric[
  "unreduced gyro integral for incoming plane wave",
  directGyroBessel[1, probeA], BesselJ[1, probeA], 10^-34
];
checkNumeric[
  "unreduced gyro integral for outgoing plane wave",
  directGyroBessel[1, probeC], BesselJ[1, probeC], 10^-34
];

If[StringLength[mutation] > 0, FinishChecks[]];

(* Krook and FP assumptions remain deliberately separate. *)
plasmaZ[z_] := I Sqrt[Pi] Exp[-z^2] Erfc[-I z];
krookParallelIntegral[z_?NumericQ] := NIntegrate[
  Exp[-u^2]/(u - z)/Sqrt[Pi], {u, -Infinity, Infinity},
  WorkingPrecision -> 55, AccuracyGoal -> 38, PrecisionGoal -> 34
];
zProbe = SetPrecision[3/5 + 4 I/5, workingPrecision];
checkNumeric[
  "Krook parallel integral gives plasma dispersion function",
  krookParallelIntegral[zProbe], plasmaZ[zProbe], 10^-34
];
harmonicProbe = SetPrecision[7/10, workingPrecision];
checkNumeric[
  "Krook cyclotron sum identity",
  Sum[BesselI[m, harmonicProbe] Exp[I m 2/5], {m, -80, 80}],
  Exp[harmonicProbe Cos[2/5]], 10^-40
];
checkNumeric[
  "homogeneous gyro sum reduces to radial Gaussian",
  -bPlus[2/5, -3/7, 5/6, 4/3] +
    bCross[2/5, -3/7, 5/6, 4/3] Cos[
      alpha[2/5, 5/6] - alpha[2/5, -3/7]],
  -(4/3)^2 ((-3/7) - 5/6)^2/2,
  10^-55
];
check[
  "Horton small-FLR Gamma0 series",
  Normal[Series[Exp[-b] BesselI[0, b], {b, 0, 3}]] ==
    1 - b + 3 b^2/4 - 5 b^3/12
];

(* Thesis (5.13)-(5.15), independently integrated here so the full kernel
   fixture does not inherit production getIfunc/W2_arr recurrence algebra. *)
generatingIntegrand[x1_, x2_, al_, be_, y_] := Module[
  {scale = 1 + x1^2, tau},
  tau = y/scale;
  Exp[(I x2 - x1^2) tau +
      (al + I x1) (be + I x1) (Exp[-tau] - 1) +
      (al + be)^2/2]/scale
];
Clear[fpMomentMatrix];
fpMomentMatrix[x1in_?NumericQ, x2in_?NumericQ] :=
  fpMomentMatrix[x1in, x2in] = Module[
    {x1 = SetPrecision[x1in, workingPrecision],
     x2 = SetPrecision[x2in, workingPrecision], expressions},
    expressions = Flatten[Table[
      D[generatingIntegrand[x1, x2, al, be, y], {al, k}, {be, l}] /.
        {al -> 0, be -> 0},
      {k, 0, 2}, {l, 0, 2}
    ]];
    Partition[
      Quiet[NIntegrate[
        Evaluate[expressions], {y, 0, Infinity},
        WorkingPrecision -> workingPrecision, AccuracyGoal -> 46,
        PrecisionGoal -> 40, MaxRecursion -> 40,
        Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
      ], {NIntegrate::precw, NIntegrate::slwcon}],
      3
    ]
  ];
fpDifferential[k_, l_, x1_, x2_] := fpMomentMatrix[x1, x2][[k + 1, l + 1]];
fpConserving[k_, l_, x1_, x2_] := Module[{denominator},
  denominator = 1 - fpDifferential[0, 0, x1, x2] +
    2 fpDifferential[2, 0, x1, x2] - fpDifferential[2, 2, x1, x2];
  fpDifferential[k, l, x1, x2] +
    (fpDifferential[k, 0, x1, x2] - fpDifferential[k, 2, x1, x2]) *
    (fpDifferential[l, 0, x1, x2] - fpDifferential[l, 2, x1, x2]) /
    denominator
];

caseRecord[id_, charge_, lambda_, vt_, omegaC_, nu_, ks_, kr_, krp_, rg_,
    mphi_, a1_, a2_, kpar_, omegaE_, vpar_, omega_, debye_] := <|
  "id" -> id, "charge" -> charge, "lambda" -> lambda, "vT" -> vt,
  "omegaC" -> omegaC, "nu" -> nu, "ks" -> ks, "kr" -> kr,
  "krp" -> krp, "rg" -> rg, "mphi" -> mphi, "A1" -> a1,
  "A2" -> a2, "kpar" -> kpar, "omegaE" -> omegaE,
  "Vpar" -> vpar, "omega" -> omega, "debye" -> debye
|>;

cases = {
  caseRecord[1,  1, 2, 11/10, 3, 4/5, 2/5, 1/5, 1/5, 7/10,  0, 3/10, -1/5, 1/2, 1/10, 1/4, 7/10, 1],
  caseRecord[2, -1, 3/2, 9/10, -4, 1, 1/3, -2/5, 7/10, 4/5, -1, -1/4, 2/5, 3/5, -1/5, 1/3, 1/2, 0],
  caseRecord[3,  1, 5/2, 6/5, 5/2, 3/4, 1/2, -3/4, 2/3, 3/5, 1, 1/5, 1/3, 2/5, 1/6, -1/5, 3/4, 0],
  caseRecord[4,  1, 7/4, 1, 2, 1, 1/4, 3/5, -4/5, 9/10, 2, -2/5, 1/2, 1/3, 1/10, 1/5, 2/3, 0],
  caseRecord[5, -1, 9/5, 4/5, -3, 6/5, 2/7, -1/2, 5/6, 2/5, -2, 1/6, -3/10, 1/2, -1/8, 1/4, 3/5, 0],
  caseRecord[6,  1, 2, 1, 1, 1, 1/10, Sqrt[600 - 1/100], Sqrt[600 - 1/100], 1/3, 0, 1/4, 1/5, 1/2, 0, 0, 1/3, 0],
  caseRecord[7,  1, 2, 1, 10^8, 1, 1/3, 2/5, -3/5, 1/2, 0, 1/3, -1/4, 1/2, 0, 0, 1/4, 1],
  caseRecord[8, -1, 3, 1, -5, 1, 2/5, 1/3, -1/4, 3/4, 0, 0, 0, 1/2, 0, 0, 1/3, 1]
};
Do[
  AppendTo[cases,
    caseRecord[9 + m + 2, 1, 2, 1, 5, 1/2, 3/20, 4/5, -1/2,
      11/20, m, 1/5, -1/6, 1/4, 1/10, 1/5, 3/5, 0]
  ],
  {m, -2, 2}
];

x1FromCase[c_] := c["kpar"] c["vT"]/c["nu"];
x2FromCase[c_] := -(c["omegaE"] + c["kpar"] c["Vpar"] +
  c["mphi"] c["omegaC"] - c["omega"])/c["nu"];
kernelParts[c_] := Module[
  {rho, bp, bx, ph, x1, x2, i00, i20, i01, i21, f0, f2, dynamic,
   rhoBDynamic, debye, total, speedOfLight = 29979245800},
  rho = Abs[c["vT"]/c["omegaC"]];
  bp = bPlus[c["ks"], c["kr"], c["krp"], rho];
  bx = bCross[c["ks"], c["kr"], c["krp"], rho];
  ph = fourierGyroPhase[c["ks"], c["kr"], c["krp"], c["rg"], c["mphi"]];
  x1 = x1FromCase[c]; x2 = x2FromCase[c];
  i00 = fpConserving[0, 0, x1, x2];
  i20 = fpConserving[2, 0, x1, x2];
  i01 = fpConserving[0, 1, x1, x2];
  i21 = fpConserving[2, 1, x1, x2];
  f0 = flrDensity[bp, bx, c["mphi"]];
  f2 = flrEnergy[bp, bx, c["mphi"]];
  dynamic = ph kernelNormalization I c["vT"]^2 c["ks"] /
    (c["lambda"]^2 c["omegaC"] c["nu"]) *
    (i00 (c["A1"] f0 + c["A2"] f2) + c["A2"] i20 f0/2);
  rhoBDynamic = -ph kernelNormalization c["vT"]^3 /
    (c["lambda"]^2 c["omegaC"] c["nu"] speedOfLight) *
    (i01 (c["A1"] f0 + c["A2"] f2) + c["A2"] i21 f0/2);
  debye = If[c["debye"] == 1 && c["mphi"] == 0,
    -ph kernelNormalization/c["lambda"]^2, 0];
  total = debye + dynamic;
  <|"rho" -> rho, "bp" -> bp, "bx" -> bx, "phase" -> ph,
    "x1" -> x1, "x2" -> x2, "I00" -> i00, "I20" -> i20,
    "I01" -> i01, "I21" -> i21,
    "F0" -> f0, "F2" -> f2, "dynamic" -> dynamic,
    "debyeTerm" -> debye, "total" -> total, "rhoBDynamic" -> rhoBDynamic|>
];

(* Geometry and limit proofs do not assume a collision model. *)
geom[kr_, krp_] := fourierGyroPhase[2/5, kr, krp, 7/10, 0] *
  flrDensity[bPlus[2/5, kr, krp, 3/4], bCross[2/5, kr, krp, 3/4], 0];
checkNumeric["exchange conjugacy of geometric kernel", geom[-3/7, 5/6],
  Conjugate[geom[5/6, -3/7]], 10^-55];
checkNumeric["diagonal b+ reduction",
  bPlus[2/5, 3/7, 3/7, 4/3], (4/3)^2 ((2/5)^2 + (3/7)^2), 10^-55];
checkNumeric["diagonal bx reduction",
  bCross[2/5, 3/7, 3/7, 4/3], (4/3)^2 ((2/5)^2 + (3/7)^2), 10^-55];
checkNumeric["real-field Fourier reconstruction",
  fourierGyroPhase[-2/5, 3/7, -5/6, 7/10, -2],
  Conjugate[fourierGyroPhase[2/5, -3/7, 5/6, 7/10, 2]], 10^-14];
check["zero-FLR density selects mphi=0",
  Limit[flrDensity[bp rho^2, bx rho^2, 0], rho -> 0] == 1 &&
  Limit[flrDensity[bp rho^2, bx rho^2, 1], rho -> 0] == 0];
check["Debye-only case has no dynamic source", cases[[8]]["A1"] == 0 && cases[[8]]["A2"] == 0];
check["electron and ion charge signs agree with signed gyrofrequency",
  And @@ ((Sign[#1["charge"]] == Sign[#1["omegaC"]]) & /@ cases)];

largeB = SetPrecision[600, workingPrecision];
largeI0Asymptotic = (1 + 1/(8 largeB) + 9/(128 largeB^2) +
  225/(3072 largeB^3))/Sqrt[2 Pi largeB];
checkNumeric["large-bx stabilized I0 expansion",
  Exp[-largeB] BesselI[0, largeB], largeI0Asymptotic, 10^-12];

(* Conservative geometric truncation bound: at equal gyro angles the full
   positive harmonic sum is one. This excludes no susceptibility damping. *)
truncationB = SetPrecision[1/10, workingPrecision];
retainedWeight = Exp[-truncationB] Sum[BesselI[m, truncationB], {m, -2, 2}];
truncationTailBound = 1 - retainedWeight;
check["mphi=-2..2 geometric tail below 1e-4", 0 < truncationTailBound < 10^-4];

formatNumber[value_] := StringReplace[
  ToString[FortranForm[N[value, 34]], InputForm], {"*^" -> "E", " " -> ""}
];

oracleRows = Table[
  Module[{p = kernelParts[c]},
    Join[
      {ToString[c["id"]], ToString[c["charge"]]},
      formatNumber /@ {c["lambda"], c["vT"], c["omegaC"], c["nu"],
        c["ks"], c["kr"], c["krp"], c["rg"]},
      {ToString[c["mphi"]]},
      formatNumber /@ {c["A1"], c["A2"], c["kpar"], c["omegaE"],
        c["Vpar"], c["omega"]},
      {ToString[c["debye"]]},
      formatNumber /@ {p["x1"], p["x2"], p["rho"], p["bp"], p["bx"],
        Re[p["phase"]], Im[p["phase"]], Re[p["I00"]], Im[p["I00"]],
        Re[p["I20"]], Im[p["I20"]], Re[p["I01"]], Im[p["I01"]],
        Re[p["I21"]], Im[p["I21"]], p["F0"], p["F2"],
        Re[p["dynamic"]], Im[p["dynamic"]], Re[p["debyeTerm"]],
        Im[p["debyeTerm"]], Re[p["total"]], Im[p["total"]],
        Re[p["rhoBDynamic"]], Im[p["rhoBDynamic"]]}
    ]
  ],
  {c, cases}
];

standardHarmonics = kernelParts /@ cases[[9 ;; 13]];
standardTotal = Total[#1["dynamic"] & /@ standardHarmonics];
m0RelativeError = Abs[standardTotal - standardHarmonics[[3]]["dynamic"]]/Abs[standardTotal];
m2ShellRelative = (Abs[standardHarmonics[[1]]["dynamic"]] +
  Abs[standardHarmonics[[5]]["dynamic"]])/Abs[standardTotal];
check["representative mphi=0 truncation error below 0.5 percent", m0RelativeError < 5/1000];
check["representative |mphi|=2 shell below 5e-5", m2ShellRelative < 5 10^-5];

oracleHeader = {
  "# KAMEL rho-Phi oracle; Markl thesis (13.1)-(13.73), (14.1)-(14.4); independent Weber/FP integrals.",
  "# Full finite-radius b+/bx include ks^2. The upstream Gaussian moment proves 1-b+-mphi.",
  "# Case 6 has bx=600; case 7 is zero-FLR; case 8 is Debye-only; cases 9-13 record mphi=-2..2.",
  "# Values generated at 70-digit working precision and committed with 34 digits.",
  "# id charge lambdaD vT omegaC nu ks kr krp rg mphi A1 A2 kpar omegaE Vpar omega debye x1 x2 rhoL bplus bcross phaseRe phaseIm I00Re I00Im I20Re I20Im I01Re I01Im I21Re I21Im F0 F2 rhoPhiDynamicRe rhoPhiDynamicIm debyeRe debyeIm rhoPhiTotalRe rhoPhiTotalIm rhoBDynamicRe rhoBDynamicIm"
};

If[regenerate,
  Export[oraclePath,
    StringRiffle[Join[oracleHeader, StringRiffle[#, " "] & /@ oracleRows], "\n"] <> "\n",
    "Text"];
  Print["Wrote ", Length[oracleRows], " rows to ", oraclePath],
  Module[{lines, parseNumber, committedRows},
    lines = Select[StringSplit[Import[oraclePath, "Text"], "\n"],
      StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &];
    parseNumber[token_] := ToExpression[StringReplace[token,
      RegularExpression["[Ee]([+-]?[0-9]+)$"] -> "*^$1"]];
    committedRows = (parseNumber /@ StringSplit[#]) & /@ lines;
    check["13 committed kernel oracle rows", Length[committedRows] == Length[oracleRows] == 13];
    If[Length[committedRows] == Length[oracleRows],
      Do[
        checkNumeric["kernel oracle row " <> ToString[row],
          committedRows[[row, -4]] + I committedRows[[row, -3]],
          parseNumber[oracleRows[[row, -4]]] + I parseNumber[oracleRows[[row, -3]]],
          10^-30];
        checkNumeric["rho-B oracle row " <> ToString[row],
          committedRows[[row, -2]] + I committedRows[[row, -1]],
          parseNumber[oracleRows[[row, -2]]] + I parseNumber[oracleRows[[row, -1]]],
          10^-30],
        {row, 1, Length[oracleRows]}
      ]
    ]
  ]
];

FinishChecks[];
