BeginTestSection["API"];
Needs["HelmholtzDecomposition`"];

(* HelmholtzDecomposition still returns {g, r} (backward compat) *)
VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}],
  {{x^2 + y^2/2, x y}, {-(1/2) (-2 + y) y, 0}},
  TestID -> "HelmholtzDecomposition list form preserved"
];

VerificationTest[
  HelmholtzGradient[{x^2 + y, x y}, {x, y}],
  {x^2 + y^2/2, x y},
  TestID -> "HelmholtzGradient component"
];

VerificationTest[
  HelmholtzRotational[{x^2 + y, x y}, {x, y}],
  {-(1/2) (-2 + y) y, 0},
  TestID -> "HelmholtzRotational component"
];

(* HelmholtzPotential is a primitive of HelmholtzGradient *)
VerificationTest[
  Module[{P, g},
    P = HelmholtzPotential[{x^2 + y, x y}, {x, y}];
    g = HelmholtzGradient[{x^2 + y, x y}, {x, y}];
    Simplify[Grad[P, {x, y}] - g]
  ],
  {0, 0},
  TestID -> "HelmholtzPotential is a primitive of HelmholtzGradient"
];

VerificationTest[
  Sort@Keys@HelmholtzAssociation[{x^2 + y, x y}, {x, y}],
  Sort@{"Gradient", "Rotational", "Potential",
        "VectorPotentialMatrix", "PotentialMatrix",
        "Residual", "Coordinates", "Field"},
  TestID -> "HelmholtzAssociation keys"
];

VerificationTest[
  HelmholtzAssociation[{x^2 + y, x y}, {x, y}]["Residual"],
  {0, 0},
  TestID -> "Residual is zero"
];

VerificationTest[
  Head@PrintHelmholtz[{x^2 + y, x y}, {x, y}],
  Association,
  TestID -> "PrintHelmholtz returns Association"
];

(* Auto-detect coordinates from the field *)
VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}],
  {{x^2 + y^2/2, x y}, {-(1/2) (-2 + y) y, 0}},
  TestID -> "Auto-detect coordinates"
];

VerificationTest[
  HelmholtzPotential[{x^2 + y, x y}],
  HelmholtzPotential[{x^2 + y, x y}, {x, y}],
  TestID -> "Auto-detect agrees with explicit coords"
];

(* PrintHelmholtz on a precomputed Association reuses it *)
VerificationTest[
  Module[{a},
    a = HelmholtzAssociation[{x^2 + y, x y}, {x, y}];
    PrintHelmholtz[a] === a
  ],
  True,
  TestID -> "PrintHelmholtz[assoc] reuses precomputed data"
];

(* Assumptions option propagates to the simplifier *)
VerificationTest[
  Module[{P},
    P = HelmholtzPotential[{a x^2, b y}, {x, y},
          Assumptions -> a > 0 && b > 0];
    Simplify[Grad[P, {x, y}] - HelmholtzGradient[{a x^2, b y}, {x, y}],
            a > 0 && b > 0]
  ],
  {0, 0},
  TestID -> "Assumptions option does not break correctness"
];

(* Simplify -> Identity skips simplification (faster, larger output) *)
VerificationTest[
  AssociationQ@HelmholtzAssociation[{x^2 + y, x y}, {x, y}, Simplify -> Identity],
  True,
  TestID -> "Simplify -> Identity returns valid Association"
];

(* Memoization: second call should be (near-)instantaneous *)
VerificationTest[
  Module[{f, xv, t1, t2},
    f = {x^3 y + y^2, x y^2 - x^2};
    xv = {x, y};
    t1 = First @ AbsoluteTiming[HelmholtzAssociation[f, xv]];
    t2 = First @ AbsoluteTiming[HelmholtzAssociation[f, xv]];
    t2 < 0.01 || t2 < t1 / 5
  ],
  True,
  TestID -> "Memoization: repeat call is much faster"
];

EndTestSection[]
