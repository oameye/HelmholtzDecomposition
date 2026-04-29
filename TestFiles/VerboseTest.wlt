BeginTestSection["Verbose option"];
Needs["HelmholtzDecomposition`"];

(* Verbose -> True/False is accepted; computation result is independent of it *)
VerificationTest[
  Module[{a, b},
    a = HelmholtzDecomposition[{x^2, y}, {x, y}, Verbose -> False];
    b = Quiet @ HelmholtzDecomposition[{x^2, y}, {x, y}, Verbose -> True];
    a["Association"] === b["Association"]
  ],
  True,
  TestID -> "Verbose option does not change result"
];

(* Non-boolean Verbose is rejected with a message *)
VerificationTest[
  HelmholtzDecomposition[{x^2, y}, {x, y}, Verbose -> "yes"],
  $Failed,
  {HelmholtzDecomposition::badverbose},
  TestID -> "Verbose -> non-boolean -> $Failed with message"
];

EndTestSection[]
