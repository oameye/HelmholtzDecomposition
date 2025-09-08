BeginTestSection["Verbose Setting Tests"];
Needs["HelmholtzDecomposition`"];

(* Test default verbose setting *)
VerificationTest[
  GetVerbose[],
  False,
  TestID -> "Default verbose setting is False"
  ]

(* Test setting verbose to True *)
VerificationTest[
  (SetVerbose[True]; GetVerbose[]),
  True,
  TestID -> "SetVerbose[True] works correctly"
  ]

(* Test setting verbose to False *)
VerificationTest[
  (SetVerbose[False]; GetVerbose[]),
  False,
  TestID -> "SetVerbose[False] works correctly"
  ]

(* Test that verbose setting affects output - test that both modes return same result *)
VerificationTest[
  Module[{result1, result2, f, coords},
    f = {x^2, y};
    coords = {x, y};
    SetVerbose[False];
    result1 = HelmholtzDecomposition[f, coords];
    SetVerbose[True];
    result2 = HelmholtzDecomposition[f, coords];
    SetVerbose[False]; (* Reset to default *)
    result1 == result2
  ],
  True,
  TestID -> "Verbose setting doesn't change calculation results"
  ]

(* Test that SetVerbose only accepts Boolean values *)
VerificationTest[
  SetVerbose["invalid"],
  $Failed,
  {SetVerbose::bool},
  TestID -> "SetVerbose rejects non-Boolean input"
  ]

EndTestSection[]
