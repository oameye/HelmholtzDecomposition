BeginTestSection["getUniqueVariables Function Tests"];
Needs["HelmholtzDecomposition`"];

(* Test the getUniqueVariables helper function *)
VerificationTest[
  HelmholtzDecomposition`Private`getUniqueVariables[Subscript[x, 2] Subscript[x, 3]],
  {Subscript[x, 2], Subscript[x, 3]},
  TestID -> "getUniqueVariables product test"
  ]

VerificationTest[
  HelmholtzDecomposition`Private`getUniqueVariables[Exp[Subscript[x, 2] Subscript[x, 3]]],
  {Subscript[x, 2], Subscript[x, 3]},
  TestID -> "getUniqueVariables exponential test"
  ]

VerificationTest[
  HelmholtzDecomposition`Private`getUniqueVariables[Sin[Subscript[x, 3]]Exp[Subscript[x, 3] Subscript[x, 2]]],
  {Subscript[x, 2], Subscript[x, 3]},
  TestID -> "getUniqueVariables complex expression test"
  ]

EndTestSection[]
