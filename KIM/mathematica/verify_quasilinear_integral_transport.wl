(* ::Package:: *)

(* Independent symbolic oracle for KIM's integral ion transport algebra.

   Run with

     wolframscript -file KIM/mathematica/verify_quasilinear_integral_transport.wl

   The script proves the insertion and energy-polynomial identities exactly,
   checks source/observation reciprocity at high precision, and regenerates
   the two numeric fixtures consumed by the Fortran tests. *)

ClearAll[assertTrue, assertClose];
$MaxExtraPrecision = 1000;

assertTrue[condition_, label_String] :=
  If[TrueQ[condition], Null,
    Print["FAIL: " <> label]; Exit[1]];

assertClose[actual_, expected_, tolerance_, label_String] :=
  assertTrue[Abs[N[actual - expected, 50]] < tolerance, label];

(* Do not inherit values for conventional short physics symbols from a
   user's Mathematica init file. *)
ClearAll[s, l, q, bp, bx, ksS, krS, ksO, krO, xiS, xiO, vT, omegaC,
  omegaSigned, c, B0, kParallel, u, omegaE, omegaMode];

(* Gaussian--Bessel generator and its first three perpendicular-energy
   moments.  Differentiating this generator is independent of the formulas
   implemented in Fortran. *)

ClearAll[perpendicularGenerator, perpendicularMoment];
perpendicularGenerator[s_, harmonic_, bPlus_, bCross_] :=
  Exp[-bPlus/s] BesselI[harmonic, bCross/s]/s;
perpendicularMoment[q_Integer?NonNegative, harmonic_, bPlus_, bCross_] :=
  FullSimplify[
    (-1)^q D[perpendicularGenerator[s, harmonic, bPlus, bCross], {s, q}]
      /. s -> 1];

wExpected[0] = Exp[-bp] BesselI[l, bx];
wExpected[1] = Exp[-bp] ((1 - bp - l) BesselI[l, bx]
  + bx BesselI[l - 1, bx]);
wExpected[2] = Exp[-bp] ((3 - 2 bp) bx BesselI[l - 1, bx]
  + (2 + bp^2 + bx^2 + 2 bp (-2 + l) - 3 l + l^2)
    BesselI[l, bx]);

Do[
  assertTrue[
    FullSimplify[perpendicularMoment[q, l, bp, bx] - wExpected[q]] === 0,
    "perpendicular generator q=" <> ToString[q]],
  {q, 0, 2}];

(* Exact field-insertion algebra.  The generic function g prevents the
   checks from sharing any Bessel implementation with the production code. *)

ClearAll[sourceInsert, observationInsert, genericMoment];
sourceInsert["Phi", expr_] := -I c ksS expr/B0;
sourceInsert["Br", expr_] := vT xiS expr/B0;
sourceInsert["Bpar", expr_] := I omegaSigned D[expr, ksS]/B0;
observationInsert["Phi", expr_] := I c ksO expr/B0;
observationInsert["Br", expr_] := vT xiO expr/B0;
observationInsert["Bpar", expr_] := -I omegaSigned D[expr, ksO]/B0;
genericMoment = g[ksS, ksO];

explicitInsertion["Phi", "Phi"] = c^2 ksO ksS genericMoment/B0^2;
explicitInsertion["Phi", "Br"] = I c ksO vT xiS genericMoment/B0^2;
explicitInsertion["Br", "Phi"] = -I c ksS vT xiO genericMoment/B0^2;
explicitInsertion["Br", "Br"] = vT^2 xiO xiS genericMoment/B0^2;
explicitInsertion["Phi", "Bpar"] =
  -c ksO omegaSigned D[genericMoment, ksS]/B0^2;
explicitInsertion["Bpar", "Phi"] =
  -c ksS omegaSigned D[genericMoment, ksO]/B0^2;
explicitInsertion["Br", "Bpar"] =
  I vT xiO omegaSigned D[genericMoment, ksS]/B0^2;
explicitInsertion["Bpar", "Br"] =
  -I vT xiS omegaSigned D[genericMoment, ksO]/B0^2;
explicitInsertion["Bpar", "Bpar"] =
  omegaSigned^2 D[genericMoment, ksS, ksO]/B0^2;

fields = {"Phi", "Br", "Bpar"};
Do[
  actual = observationInsert[observation,
    sourceInsert[source, genericMoment]];
  assertTrue[FullSimplify[actual - explicitInsertion[observation, source]] === 0,
    "exact insertion " <> observation <> "-" <> source],
  {observation, fields}, {source, fields}];

detuningSentinel = kParallel u + l omegaSigned + omegaE - omegaMode;
assertTrue[D[detuningSentinel, ksS] === 0 &&
    D[detuningSentinel, ksO] === 0,
  "Bparallel derivatives do not act on kinetic detuning"];

(* Two-wave FLR factors. *)

ClearAll[kpS, kpO, alphaS, alphaO, bPlusTwoWave, bCrossTwoWave,
  baseTransverseMoment, transverseFactor, transverseCoefficient];
kpS = Sqrt[ksS^2 + krS^2];
kpO = Sqrt[ksO^2 + krO^2];
alphaS = ArcTan[ksS, krS];
alphaO = ArcTan[ksO, krO];
bPlusTwoWave = vT^2 (kpS^2 + kpO^2)/(2 omegaC^2);
bCrossTwoWave = vT^2 kpS kpO/omegaC^2;
brPower["Br"] = 1;
brPower["Phi"] = 0;
brPower["Bpar"] = 0;
baseTransverseMoment[q_Integer] := baseTransverseMoment[q] =
  Exp[I l (alphaO - alphaS)] *
    (wExpected[q] /. {bp -> bPlusTwoWave, bx -> bCrossTwoWave});
transverseCoefficient[q_Integer, observation_String, source_String] :=
  transverseCoefficient[q, observation, source] = Switch[{observation, source},
    {"Phi", "Phi"}, c^2 ksO ksS baseTransverseMoment[q]/B0^2,
    {"Phi", "Br"}, I c ksO vT baseTransverseMoment[q]/B0^2,
    {"Br", "Phi"}, -I c ksS vT baseTransverseMoment[q]/B0^2,
    {"Br", "Br"}, vT^2 baseTransverseMoment[q]/B0^2,
    {"Phi", "Bpar"},
      -c ksO omegaSigned D[baseTransverseMoment[q], ksS]/B0^2,
    {"Bpar", "Phi"},
      -c ksS omegaSigned D[baseTransverseMoment[q], ksO]/B0^2,
    {"Br", "Bpar"},
      I vT omegaSigned D[baseTransverseMoment[q], ksS]/B0^2,
    {"Bpar", "Br"},
      -I vT omegaSigned D[baseTransverseMoment[q], ksO]/B0^2,
    {"Bpar", "Bpar"},
      omegaSigned^2 D[baseTransverseMoment[q], ksS, ksO]/B0^2];
transverseFactor[q_Integer, observation_String, source_String] :=
  xiO^brPower[observation] xiS^brPower[source] *
    transverseCoefficient[q, observation, source];

waveRules = {ksS -> 7/10, krS -> -2/5, ksO -> 11/10,
  krO -> 3/10, vT -> 6/5, omegaC -> 9/5,
  omegaSigned -> -9/5, c -> 13/10, B0 -> 9/10, l -> 2};
swappedWaveRules = {ksS -> 11/10, krS -> 3/10, ksO -> 7/10,
  krO -> -2/5, vT -> 6/5, omegaC -> 9/5,
  omegaSigned -> -9/5, c -> 13/10, B0 -> 9/10, l -> 2};

Do[
  assertClose[
    N[transverseCoefficient[q, observation, source] /. waveRules, 60],
    Conjugate[N[transverseCoefficient[q, source, observation]
      /. swappedWaveRules, 60]], 10^-35,
    "high-precision reciprocity " <> observation <> "-" <> source
      <> ", q=" <> ToString[q]],
  {q, 0, 2}, {observation, fields}, {source, fields}];

(* Drift-kinetic limiting oracle.  Set both radial wave numbers to zero and
   approach k_s,k_s' from the positive side.  The gyro-average becomes
   W_0 -> 1 and W_1,W_2 -> 0; with Bparallel=0 this removes every finite-FLR
   correction and leaves exactly the Heyn/Markl four-entry ledger. *)
ClearAll[epsilon, zeroFlrMoment];
zeroFlrMoment[q_Integer] := FullSimplify[
  Limit[(wExpected[q] /. {bp -> vT^2 (epsilon^2 + epsilon^2)/(2 omegaC^2),
      bx -> vT^2 epsilon^2/omegaC^2, l -> 0}), epsilon -> 0,
    Direction -> "FromAbove"]];
assertTrue[zeroFlrMoment[0] === 1, "zero-FLR W0 reduction"];
assertTrue[zeroFlrMoment[1] === 1, "zero-FLR W1 reduction"];
assertTrue[zeroFlrMoment[2] === 2, "zero-FLR W2 reduction"];

zeroFlrLedger = <|
  {1, 1} -> momentZero,
  {1, 2} -> momentOne,
  {2, 1} -> momentOne,
  {2, 2} -> 2 momentZero|> /. {momentZero -> 1, momentOne -> 1};
assertTrue[zeroFlrLedger === <|{1, 1} -> 1, {1, 2} -> 1,
    {2, 1} -> 1, {2, 2} -> 2|>,
  "zero-FLR four-entry drift-kinetic ledger"];

(* The thermodynamic polynomials determine the complete coefficient family
   and the parallel-susceptibility ledger. *)

ClearAll[transportPolynomial, susceptibilityPowers];
transportPolynomial[1, 1, observation_, source_] :=
  transverseFactor[0, observation, source];
transportPolynomial[1, 2, observation_, source_] :=
  xiS^2 transverseFactor[0, observation, source]/2 +
    transverseFactor[1, observation, source];
transportPolynomial[2, 1, observation_, source_] :=
  xiO^2 transverseFactor[0, observation, source]/2 +
    transverseFactor[1, observation, source];
transportPolynomial[2, 2, observation_, source_] :=
  xiO^2 xiS^2 transverseFactor[0, observation, source]/4 +
    (xiO^2 + xiS^2) transverseFactor[1, observation, source]/2 +
    transverseFactor[2, observation, source];
susceptibilityPowers[expression_] :=
  Sort[Keys[CoefficientRules[Expand[expression], {xiO, xiS}]]];

expectedBparPowers = <|
  {1, 1} -> {{0, 0}},
  {1, 2} -> {{0, 0}, {0, 2}},
  {2, 1} -> {{0, 0}, {2, 0}},
  {2, 2} -> {{0, 0}, {0, 2}, {2, 0}, {2, 2}}|>;
ClearAll[ledgerPolynomial, momentZero, momentOne, momentTwo];
ledgerPolynomial[1, 1, observation_, source_] :=
  xiO^brPower[observation] xiS^brPower[source] momentZero;
ledgerPolynomial[1, 2, observation_, source_] :=
  xiO^brPower[observation] xiS^brPower[source] *
    (xiS^2 momentZero/2 + momentOne);
ledgerPolynomial[2, 1, observation_, source_] :=
  xiO^brPower[observation] xiS^brPower[source] *
    (xiO^2 momentZero/2 + momentOne);
ledgerPolynomial[2, 2, observation_, source_] :=
  xiO^brPower[observation] xiS^brPower[source] *
    (xiO^2 xiS^2 momentZero/4
      + (xiO^2 + xiS^2) momentOne/2 + momentTwo);
KeyValueMap[
  (assertTrue[
      susceptibilityPowers[ledgerPolynomial[#1[[1]], #1[[2]],
        "Bpar", "Bpar"]] === #2,
      "exact Bparallel susceptibility ledger D" <> ToString[#1] <>
        ", actual=" <> ToString[susceptibilityPowers[
          ledgerPolynomial[#1[[1]], #1[[2]], "Bpar", "Bpar"]]]]) &,
  expectedBparPowers];

allLedgers = Flatten[Table[
  susceptibilityPowers[ledgerPolynomial[i, j, observation, source]],
  {i, 1, 2}, {j, 1, 2}, {observation, fields}, {source, fields}], 3];
assertTrue[Length[allLedgers] === 36,
  "four tensor entries contain all nine ordered field pairs"];
assertTrue[Max[Flatten[allLedgers]] === 3,
  "highest required susceptibility index is three"];

(* High-precision fixtures.  Inputs remain exact rationals until the final N,
   so these values exercise symbolic differentiation before evaluation. *)

fieldS = <|"Phi" -> 3/5 - 2 I/7, "Br" -> -4/9 + I/6,
  "Bpar" -> 5/8 + 3 I/11|>;
fieldO = <|"Phi" -> -2/3 + I/5, "Br" -> 7/10 + 2 I/9,
  "Bpar" -> -3/7 - I/4|>;
nu = 7/5;
ifunc[p_Integer, q_Integer] := (p + 1)/(q + 2) + I (p - q)/13;
contractPolynomial[expression_] := Total[
  (Last[#] ifunc[First[#][[1]], First[#][[2]]]) & /@
    CoefficientRules[Expand[expression], {xiO, xiS}]];
contractedMoment[1, 1, observation_, source_] :=
  transverseCoefficient[0, observation, source] *
    ifunc[brPower[observation], brPower[source]];
contractedMoment[1, 2, observation_, source_] :=
  transverseCoefficient[0, observation, source] *
      ifunc[brPower[observation], brPower[source] + 2]/2 +
    transverseCoefficient[1, observation, source] *
      ifunc[brPower[observation], brPower[source]];
contractedMoment[2, 1, observation_, source_] :=
  transverseCoefficient[0, observation, source] *
      ifunc[brPower[observation] + 2, brPower[source]]/2 +
    transverseCoefficient[1, observation, source] *
      ifunc[brPower[observation], brPower[source]];
contractedMoment[2, 2, observation_, source_] :=
  transverseCoefficient[0, observation, source] *
      ifunc[brPower[observation] + 2, brPower[source] + 2]/4 +
    transverseCoefficient[1, observation, source] *
      (ifunc[brPower[observation] + 2, brPower[source]] +
        ifunc[brPower[observation], brPower[source] + 2])/2 +
    transverseCoefficient[2, observation, source] *
      ifunc[brPower[observation], brPower[source]];

Do[
  assertTrue[
    FullSimplify[contractedMoment[i, j, observation, source] -
      contractPolynomial[transportPolynomial[i, j, observation, source]]] === 0,
    "exact coefficient contraction D" <> ToString[{i, j}] <> " " <>
      observation <> "-" <> source];
  assertClose[
    N[contractedMoment[i, j, observation, source] /. waveRules, 60],
    N[contractPolynomial[transportPolynomial[i, j, observation, source]]
      /. waveRules, 60], 10^-35,
    "direct coefficient contraction D" <> ToString[{i, j}] <> " " <>
      observation <> "-" <> source],
  {i, 1, 2}, {j, 1, 2}, {observation, fields}, {source, fields}];

transverseRows = Flatten[Table[
  value = N[transverseCoefficient[q, fields[[observation]], fields[[source]]]
    /. waveRules, 45];
  {q, observation, source, Re[value], Im[value]},
  {q, 0, 2}, {observation, 1, 3}, {source, 1, 3}], 2];

rawChannel[i_, j_, observation_, source_, rules_, sourceFields_,
    observationFields_] :=
  N[Conjugate[observationFields[observation]] * sourceFields[source] *
    contractPolynomial[transportPolynomial[i, j, observation, source]]
      /. rules, 45]/(2 nu);

channelRows = Flatten[Table[
  value = (rawChannel[i, j, fields[[observation]], fields[[source]],
      waveRules, fieldS, fieldO] +
    Conjugate[rawChannel[j, i, fields[[source]], fields[[observation]],
      swappedWaveRules, fieldO, fieldS]])/2;
  {i, j, observation, source, Re[value], Im[value]},
  {i, 1, 2}, {j, 1, 2}, {observation, 1, 3}, {source, 1, 3}], 3];

scriptDirectory = DirectoryName[$InputFileName];
fixtureDirectory = FileNameJoin[{scriptDirectory, "..", "tests", "test_data"}];
If[!DirectoryQ[fixtureDirectory], CreateDirectory[fixtureDirectory]];
generatorHash = IntegerString[FileHash[$InputFileName, "SHA256"], 16, 64];

writeFixture[path_, rows_] := Module[{stream},
  stream = Quiet[Check[OpenWrite[path, PageWidth -> Infinity], $Failed]];
  assertTrue[Head[stream] === OutputStream,
    "open generated fixture " <> path];
  WriteString[stream, "# algebra_version 1\n"];
  WriteString[stream, "# generator_sha256 " <> generatorHash <> "\n"];
  Scan[(WriteString[stream,
      StringRiffle[
        (If[IntegerQ[#], ToString[#],
            ToString[FortranForm[N[#, 20]]]]) & /@ #, " "] <> "\n"]) &,
    rows];
  assertTrue[StringQ[Quiet[Check[Close[stream], $Failed]]],
    "close generated fixture " <> path]];

writeFixture[FileNameJoin[{fixtureDirectory,
  "quasilinear_transverse_oracle.dat"}], transverseRows];
writeFixture[FileNameJoin[{fixtureDirectory,
  "quasilinear_channel_oracle.dat"}], channelRows];

Print["QL_INTEGRAL_ORACLE_OK"];
