BeginTestSection["Failure modes"];
Needs["HelmholtzDecomposition`"];

(* Length mismatch *)
VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y, z}],
  $Failed,
  {HelmholtzDecomposition::nomatch},
  TestID -> "Length mismatch -> $Failed with message"
];

(* Auto-coord with too few/many symbols *)
VerificationTest[
  HelmholtzDecomposition[{x^2 + y + z, x y}],
  $Failed,
  {HelmholtzDecomposition::autocoord},
  TestID -> "Auto-coord fails when symbols != length"
];

(* Undecomposable atom -> graceful $Failed *)
VerificationTest[
  HelmholtzDecomposition[{Tan[x y], 0}, {x, y}],
  $Failed,
  {HelmholtzDecomposition::nodec},
  TestID -> "Undecomposable atom -> $Failed via Throw, not Abort"
];

(* Unknown property *)
VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["NotAProperty"],
  $Failed,
  {HelmholtzDecomposition::nokey},
  TestID -> "Unknown property -> $Failed with message"
];

(* Non-boolean Caching *)
VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}, Caching -> "yes"],
  $Failed,
  {HelmholtzDecomposition::badcache},
  TestID -> "Caching -> non-boolean -> $Failed with message"
];

VerificationTest[
  HelmholtzGaugeShift[HelmholtzDecomposition[{x^2 + y, x y}, {x, y}], x^2],
  $Failed,
  {HelmholtzGaugeShift::badharm},
  TestID -> "Non-harmonic scalar gauge -> $Failed with message"
];

VerificationTest[
  HelmholtzGaugeShift[HelmholtzDecomposition[{x^2 + y, x y}, {x, y}], {1, 1, 1}],
  $Failed,
  {HelmholtzGaugeShift::badlen},
  TestID -> "Gauge vector length mismatch -> $Failed with message"
];

VerificationTest[
  HelmholtzGaugeShift[HelmholtzDecomposition[{x^2 + y, x y}, {x, y}], {x, 0}],
  $Failed,
  {HelmholtzGaugeShift::badharm},
  TestID -> "Non-harmonic vector gauge -> $Failed with message"
];

(* MakeBoxes on a malformed bare Association does not crash *)
VerificationTest[
  Head @ ToBoxes[HelmholtzDecomposition[<|"foo" -> 1|>], StandardForm],
  RowBox,
  TestID -> "MakeBoxes on malformed wrapper falls through to default"
];

(* MakeBoxes on a real result evaluates without error *)
VerificationTest[
  Head @ ToBoxes[HelmholtzDecomposition[{x^2 + y, x y}, {x, y}], StandardForm],
  InterpretationBox,
  TestID -> "MakeBoxes on a real result produces an InterpretationBox"
];

(* Field shape: a 2x2 matrix is not a vector field *)
VerificationTest[
  Head @ HelmholtzDecomposition[{{1, 2}, {3, 4}}, {x, y}],
  HelmholtzDecomposition,
  TestID -> "Matrix in place of f is rejected (left unevaluated)"
];

EndTestSection[]
