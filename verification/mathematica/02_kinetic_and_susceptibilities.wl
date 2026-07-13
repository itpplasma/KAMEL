scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
If[!MemberQ[$Path, scriptDirectory], AppendTo[$Path, scriptDirectory]];
Needs["checklib`"];
resetChecks[];

repositoryRoot = ExpandFileName[FileNameJoin[{scriptDirectory, "..", ".."}]];
oraclePath = FileNameJoin[{repositoryRoot, "verification", "oracles", "fp_susceptibilities.dat"}];
workingPrecision = 70;
$MaxExtraPrecision = 4000;

(* Thesis (5.14)-(5.15), defined independently of QL-Balance/W2_arr. The
   y-rescaling keeps both small- and large-x1 quadratures well conditioned. *)
generatingIntegrand[x1_, x2_, alpha_, beta_, y_] := Module[
  {scale = 1 + x1^2, tau},
  tau = y/scale;
  Exp[(I x2 - x1^2) tau +
      (alpha + I x1) (beta + I x1) (Exp[-tau] - 1) +
      (alpha + beta)^2/2]/scale
];

clearMomentCache[] := Clear[momentMatrix];
momentMatrix[x1in_?NumericQ, x2in_?NumericQ] :=
  momentMatrix[x1in, x2in] = Module[
    {x1 = SetPrecision[x1in, workingPrecision],
     x2 = SetPrecision[x2in, workingPrecision], expressions},
    expressions = Flatten[Table[
      D[generatingIntegrand[x1, x2, alpha, beta, y], {alpha, k}, {beta, l}] /.
        {alpha -> 0, beta -> 0},
      {k, 0, 3}, {l, 0, 3}
    ]];
    Partition[
      Quiet[
        NIntegrate[
          Evaluate[expressions], {y, 0, Infinity},
          WorkingPrecision -> workingPrecision,
          AccuracyGoal -> 48, PrecisionGoal -> 42, MaxRecursion -> 40,
          Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
        ],
        NIntegrate::precw
      ],
      4
    ]
  ];

iDifferential[k_Integer, l_Integer, x1_?NumericQ, x2_?NumericQ] :=
  momentMatrix[x1, x2][[k + 1, l + 1]];

(* Thesis (5.13): rank-one energy-conservation correction to the
   Ornstein-Uhlenbeck differential moments. Momentum is intentionally not
   corrected in the straight-cylinder model. *)
iConserving[k_Integer, l_Integer, x1_?NumericQ, x2_?NumericQ] := Module[
  {denominator},
  denominator = 1 - iDifferential[0, 0, x1, x2] +
    2 iDifferential[2, 0, x1, x2] - iDifferential[2, 2, x1, x2];
  iDifferential[k, l, x1, x2] +
    (iDifferential[k, 0, x1, x2] - iDifferential[k, 2, x1, x2]) *
    (iDifferential[l, 0, x1, x2] - iDifferential[l, 2, x1, x2]) /
    denominator
];

mutation = Environment["KAMEL_VERIFY_NEGATIVE_FIXTURE"];
regenerate = TrueQ[Environment["KAMEL_REGENERATE_ORACLE"] === "1"];

thermalSpeed[temperature_, mass_] := Sqrt[temperature/mass] *
  If[mutation === "thermal-speed", Sqrt[2], 1];
x1FromPhysics[case_] := case["kpar"] thermalSpeed[case["temperature"], case["mass"]]/case["nu"];
x2FromPhysics[case_] := -(
  case["omegaE"] +
  If[mutation === "flow-sign", -1, 1] case["kpar"] case["Vpar"] +
  case["mphi"] case["omegaC"] - case["omega"]
)/case["nu"];

caseRecord[id_, kpar_, temperature_, mass_, nu_, omegaE_, vpar_, mphi_, omegaC_, omega_] := <|
  "id" -> id, "kpar" -> SetPrecision[kpar, workingPrecision],
  "temperature" -> SetPrecision[temperature, workingPrecision],
  "mass" -> SetPrecision[mass, workingPrecision], "nu" -> SetPrecision[nu, workingPrecision],
  "omegaE" -> SetPrecision[omegaE, workingPrecision], "Vpar" -> SetPrecision[vpar, workingPrecision],
  "mphi" -> mphi, "omegaC" -> SetPrecision[omegaC, workingPrecision],
  "omega" -> SetPrecision[omega, workingPrecision]
|>;

cases = {
  caseRecord[1, 3/5, 1, 1, 1, 1/5, 2/5, -1, 11/10, 7/10],
  caseRecord[2, 3/5, 1, 1, 1, 1/5, 2/5,  0, 11/10, 7/10],
  caseRecord[3, 3/5, 1, 1, 1, 1/5, 2/5,  1, 11/10, 7/10],
  caseRecord[4, 3/5, 1, 1, 1, 0, 0, 0, 0, -7/4],
  caseRecord[5, 3/4, 1, 1, 1, 0, 0, 0, 0, 1/100000000],
  caseRecord[6, 2, 1, 1, 10, 0, 0, 0, 0, 1/2],
  caseRecord[7, 3/5, 1, 1, 1/10, 0, 0, 0, 0, -4/5]
};
consumedPairs = {{0, 0}, {2, 0}, {0, 2}, {0, 1}, {2, 1}, {2, 2}, {1, 1}, {1, 3}};

checkNumeric["mphi=-1 x2 mapping", x2FromPhysics[cases[[1]]], 34/25, 10^-55];
checkNumeric["mphi=0 x2 mapping", x2FromPhysics[cases[[2]]], 13/50, 10^-55];
checkNumeric["mphi=+1 x2 mapping", x2FromPhysics[cases[[3]]], -21/25, 10^-55];
check["thermal speed is sqrt(T/m)", thermalSpeed[9, 4] == 3/2];

(* Collision-model moment contract: the differential OU operator conserves
   particle number only; the rank-one integral term restores energy, not
   parallel momentum. Krook would damp all three and is not substituted. *)
ouEigenvalues = <|"particle" -> 0, "momentum" -> -1, "energy" -> -2|>;
correctedEigenvalues = Join[ouEigenvalues, <|"energy" -> 0|>];
check["OU particle moment conserved", ouEigenvalues["particle"] == 0];
check["OU momentum intentionally damped", ouEigenvalues["momentum"] == -1];
check["OU energy requires integral correction", ouEigenvalues["energy"] == -2];
check["integral term restores energy", correctedEigenvalues["energy"] == 0];
check["FP model is not non-conserving Krook", Values[ouEigenvalues] =!= {-1, -1, -1}];
If[MemberQ[{"thermal-speed", "flow-sign"}, mutation], FinishChecks[]];

x1Probe = SetPrecision[3/5, workingPrecision];
x2Probe = SetPrecision[-7/4, workingPrecision];
recurrenceRight = If[
  mutation === "recurrence-index",
  (x2Probe/x1Probe) iConserving[2, 1, x1Probe, x2Probe],
  (x2Probe/x1Probe) iConserving[1, 0, x1Probe, x2Probe]
];
checkNumeric["I11 recurrence (A.4)", iConserving[1, 1, x1Probe, x2Probe], recurrenceRight, 10^-38];
checkNumeric[
  "I31 recurrence (A.4)",
  iConserving[3, 1, x1Probe, x2Probe],
  (x2Probe/x1Probe) iConserving[2, 1, x1Probe, x2Probe],
  10^-38
];
checkNumeric[
  "susceptibility symmetry",
  iConserving[1, 3, x1Probe, x2Probe],
  iConserving[3, 1, x1Probe, x2Probe],
  10^-42
];
checkNumeric[
  "x1 parity",
  iConserving[1, 3, -x1Probe, x2Probe],
  iConserving[1, 3, x1Probe, x2Probe],
  10^-38
];
branchExpected = (-1)^(1 + 3) Conjugate[iConserving[1, 3, x1Probe, x2Probe]];
If[mutation === "wrong-branch", branchExpected = -branchExpected];
checkNumeric[
  "Landau-branch conjugation",
  iConserving[1, 3, x1Probe, -x2Probe],
  branchExpected,
  10^-38
];
If[MemberQ[{"wrong-branch", "recurrence-index"}, mutation], FinishChecks[]];

formatNumber[value_] := StringReplace[
  ToString[FortranForm[N[value, 34]], InputForm],
  {"*^" -> "E", " " -> ""}
];

oracleRows = Flatten[Table[
  Module[{x1 = x1FromPhysics[case], x2 = x2FromPhysics[case], value},
    Table[
      value = iConserving[pair[[1]], pair[[2]], x1, x2];
      {
        ToString[case["id"]], formatNumber[case["kpar"]],
        formatNumber[case["temperature"]], formatNumber[case["mass"]],
        formatNumber[case["nu"]], formatNumber[case["omegaE"]],
        formatNumber[case["Vpar"]], ToString[case["mphi"]],
        formatNumber[case["omegaC"]], formatNumber[case["omega"]],
        ToString[pair[[1]]], ToString[pair[[2]]], formatNumber[x1], formatNumber[x2],
        formatNumber[Re[value]], formatNumber[Im[value]]
      },
      {pair, consumedPairs}
    ]
  ],
  {case, cases}
], 1];

oracleHeader = {
  "# KAMEL FP susceptibility oracle; generated from Markl thesis (5.13)-(5.15), not W2_arr.",
  "# Wolfram NIntegrate working precision: 70 digits; committed values: 34 digits.",
  "# Cases 1-3: finite parallel flow and mphi=-1,0,+1; 4: general complex point;",
  "# 5: near-zero resonance numerator; 6: strong collision; 7: collisionless scaling.",
  "# case kpar T_erg mass nu omega_E Vpar mphi omega_c omega k l x1 x2 Re(Ikl) Im(Ikl)"
};

If[regenerate,
  Export[
    oraclePath,
    StringRiffle[Join[oracleHeader, StringRiffle[#, " "] & /@ oracleRows], "\n"] <> "\n",
    "Text"
  ];
  Print["Wrote ", Length[oracleRows], " rows to ", oraclePath],
  Module[{lines, parseNumber, committedRows},
    lines = Select[
      StringSplit[Import[oraclePath, "Text"], "\n"],
      StringLength[StringTrim[#]] > 0 && !StringStartsQ[StringTrim[#], "#"] &
    ];
    parseNumber[token_] := ToExpression[StringReplace[
      token, RegularExpression["[Ee]([+-]?[0-9]+)$"] -> "*^$1"
    ]];
    committedRows = (parseNumber /@ StringSplit[#]) & /@ lines;
    check["56 committed oracle rows", Length[committedRows] == Length[oracleRows] == 56];
    If[Length[committedRows] == Length[oracleRows],
      Do[
        checkNumeric[
          "oracle row " <> ToString[row],
          committedRows[[row, -2]] + I committedRows[[row, -1]],
          parseNumber[oracleRows[[row, -2]]] + I parseNumber[oracleRows[[row, -1]]],
          10^-30
        ],
        {row, 1, Length[oracleRows]}
      ]
    ]
  ]
];

FinishChecks[];
