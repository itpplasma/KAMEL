scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
If[!MemberQ[$Path, scriptDirectory], AppendTo[$Path, scriptDirectory]];
Needs["checklib`"];
resetChecks[];

repositoryRoot = ExpandFileName[FileNameJoin[{scriptDirectory, "..", ".."}]];
indexPath = FileNameJoin[{repositoryRoot, "verification", "FORMULA_INDEX.md"}];
indexText = Import[indexPath, "Text"];

chapterMaxima = <|
  "1" -> 1, "2" -> 54, "3" -> 16, "4" -> 45, "5" -> 43,
  "6" -> 19, "7" -> 9, "8" -> 12, "9" -> 29, "10" -> 2,
  "11" -> 9, "12" -> 28, "13" -> 159, "14" -> 6, "15" -> 33,
  "A" -> 29
|>;
expectedEquations = Flatten[KeyValueMap[
  Function[{chapter, maximum},
    Table[chapter <> "." <> ToString[number], {number, 1, maximum}]
  ],
  chapterMaxima
]];
indexedEquations = StringCases[
  indexText,
  RegularExpression["(?m)^\\| \\(((?:[0-9]+|A)\\.[0-9]+)\\) \\|"] -> "$1"
];

codeSiteRows = StringCases[
  indexText,
  RegularExpression["(?m)^\\| ((?:KIM|KILCA|QLB|EQ|PRE)-[0-9]+) \\|"] -> "$1"
];
conventionRows = StringCases[
  indexText,
  RegularExpression["(?m)^\\| (C-[a-z0-9-]+) \\|"] -> "$1"
];

negativeFixture = Environment["KAMEL_VERIFY_NEGATIVE_FIXTURE"];
If[negativeFixture === "missing-thesis", indexedEquations = Rest[indexedEquations]];
If[negativeFixture === "missing-code", codeSiteRows = Rest[codeSiteRows]];
If[negativeFixture === "unclassified", conventionRows = Rest[conventionRows]];

check["formula index exists", FileExistsQ[indexPath]];
check["494 unique thesis equations indexed", Length[DeleteDuplicates[indexedEquations]] == 494];
check["complete thesis equation labels", Sort[indexedEquations] === Sort[expectedEquations]];
check["no duplicate thesis equation rows", DuplicateFreeQ[indexedEquations]];
check["24 semantic code formula sites classified", Length[codeSiteRows] == 24];
check["19 conventions classified", Length[conventionRows] == 19];

codePaths = {
  "KIM/src/background_equilibrium/species_mod.f90",
  "KIM/src/background_equilibrium/calculate_equil.f90",
  "KIM/src/asymptotics/flr2_fourier_kernel.f90",
  "KIM/src/kernels/FP_kernel_plasma_prefacs.f90",
  "KIM/src/kernels/kernel.f90",
  "KIM/src/electrostatic_poisson/poisson.f90",
  "KIM/src/electrostatic_poisson/periodic_solve.f90",
  "KIM/src/asymptotics/FLR2_asymptotics.f90",
  "KIM/src/grid/prepare_resonances.f90",
  "KIM/src/background_equilibrium/profile_input_m.f90",
  "KIM/src/background_equilibrium/periodic_background.f90",
  "KiLCA/flre/conductivity/calc_I_array.f90",
  "KiLCA/solver/VER_5_STABLE/wave_stuff.f90",
  "KiLCA/flre/conductivity/calc_I_array_drift_serg.f90",
  "PreProc/fourier/src/rhs_flt.f90",
  "QL-Balance/src/base/getIfunc.f90",
  "QL-Balance/src/base/get_dql.f90",
  "QL-Balance/src/base/rhs_balance_m.f90",
  "QL-Balance/src/base/calc_current_densities.f90",
  "common/equil/mag_wrapper.f90",
  "common/equil/equil_profiles.f90",
  "PreProc/fourier/src/fouriermodes.f90"
};
check["all indexed code paths exist", And @@ (FileExistsQ[FileNameJoin[{repositoryRoot, #}]] & /@ codePaths)];

registerConvention["units", "Gaussian CGS; profile temperature in eV", "FORMULA_INDEX C-unit"];
registerConvention["charge", "e_sigma=Zspec e; species 0 electron", "FORMULA_INDEX C-charge"];
registerConvention["vT", "sqrt(T ev/m)", "FORMULA_INDEX C-vT"];
registerConvention["gyro", "Zspec e abs(B)/(m c)", "FORMULA_INDEX C-gyro"];
registerConvention["phase", "exp(+i(m theta+n phi)-i omega t)", "FORMULA_INDEX C-phase"];
registerConvention["radial Fourier", "inverse exp(+i k r)/(2 pi)", "FORMULA_INDEX C-radial-FT"];
registerConvention["q resonance", "q=r Bz/(R0 Btheta); qres=abs(m/n)", "FORMULA_INDEX C-q"];
registerConvention["Poisson", "Laplacian Phi=-4 pi rho", "FORMULA_INDEX C-poisson"];
If[negativeFixture === "duplicate-convention",
  registerConvention["phase", "exp(-i(m theta+n phi)-i omega t)", "negative fixture"]
];
check["eight executable convention definitions registered", Length[registeredConventions[]] == 8];

massDimension = {1, 0, 0};
lengthDimension = {0, 1, 0};
timeDimension = {0, 0, 1};
chargeDimensionCGS = {1/2, 3/2, -1};
magneticFieldDimensionCGS = {1/2, -1/2, -1};
speedDimension = lengthDimension - timeDimension;
gyroDimension = chargeDimensionCGS + magneticFieldDimensionCGS - massDimension - speedDimension;
check["CGS gyrofrequency has inverse-time units", gyroDimension === -timeDimension];
energyDimension = massDimension + 2 lengthDimension - 2 timeDimension;
thermalSpeedDimension = (energyDimension - massDimension)/2;
check["thermal speed has velocity units", thermalSpeedDimension === speedDimension];

phase[mm_, nn_] := Exp[I (mm theta + nn phi)];
check[
  "Fourier conjugate reality",
  TrueQ[Simplify[
    phase[-mm, -nn] == Conjugate[phase[mm, nn]],
    Element[{mm, nn, theta, phi}, Reals]
  ]]
];
check[
  "radial Fourier signs are inverse pairs",
  Simplify[Exp[I kr r] Exp[-I kr r], Element[{kr, r}, Reals]] === 1
];

FinishChecks[];
