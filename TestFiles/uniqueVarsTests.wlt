BeginTestSection["uniqueVars helper"];
Needs["HelmholtzDecomposition`"];

VerificationTest[
  HelmholtzDecomposition`Private`uniqueVars[Subscript[x, 2] Subscript[x, 3]],
  {Subscript[x, 2], Subscript[x, 3]},
  TestID -> "uniqueVars: product"
];

VerificationTest[
  HelmholtzDecomposition`Private`uniqueVars[Exp[Subscript[x, 2] Subscript[x, 3]]],
  {Subscript[x, 2], Subscript[x, 3]},
  TestID -> "uniqueVars: exponential"
];

VerificationTest[
  HelmholtzDecomposition`Private`uniqueVars[
    Sin[Subscript[x, 3]] Exp[Subscript[x, 3] Subscript[x, 2]]],
  {Subscript[x, 2], Subscript[x, 3]},
  TestID -> "uniqueVars: trig+exp composite"
];

EndTestSection[]
