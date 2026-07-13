BeginPackage["checklib`"];

check::usage = "check[name, condition] records one Boolean verification.";
checkNumeric::usage = "checkNumeric[name, actual, expected, tolerance] checks a scaled numerical error.";
registerConvention::usage = "registerConvention[name, definition, source] rejects incompatible duplicates.";
registeredConventions::usage = "registeredConventions[] returns the convention registry.";
resetChecks::usage = "resetChecks[] clears checks and registered conventions.";
FinishChecks::usage = "FinishChecks[] prints the summary and exits nonzero after any failure.";

Begin["`Private`"];

$checkCount = 0;
$failureCount = 0;
$conventionDefinitions = <||>;
$conventionSources = <||>;

check[name_String, condition_] := Module[{passed = TrueQ[condition]},
  $checkCount++;
  If[passed,
    Print["[PASS] ", name],
    $failureCount++;
    Print["[FAIL] ", name]
  ];
  passed
];

checkNumeric[name_String, actual_, expected_, tolerance_] := Module[
  {scale = Max[1, Abs[N[expected]]], error},
  error = Abs[N[actual] - N[expected]]/scale;
  check[name <> " (scaled error " <> ToString[ScientificForm[error, 3]] <> ")", error <= tolerance]
];

registerConvention[name_String, definition_String, source_String] := Module[{},
  If[StringLength[StringTrim[name]] == 0 ||
     StringLength[StringTrim[definition]] == 0 ||
     StringLength[StringTrim[source]] == 0,
    check["classified convention " <> name, False];
    Return[False]
  ];
  If[KeyExistsQ[$conventionDefinitions, name],
    Return[check[
      "compatible duplicate convention " <> name,
      SameQ[$conventionDefinitions[name], definition]
    ]]
  ];
  AssociateTo[$conventionDefinitions, name -> definition];
  AssociateTo[$conventionSources, name -> source];
  True
];

registeredConventions[] := Association[$conventionDefinitions];

resetChecks[] := (
  $checkCount = 0;
  $failureCount = 0;
  $conventionDefinitions = <||>;
  $conventionSources = <||>;
);

FinishChecks[] := Module[{exitCode = If[$failureCount == 0, 0, 1]},
  Print["Checks: ", $checkCount, "; failures: ", $failureCount];
  Exit[exitCode]
];

End[];
EndPackage[];
