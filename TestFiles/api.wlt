BeginTestSection["API v2"];
Needs["HelmholtzDecomposition`"];

(* ---- shape of the result ---- *)

VerificationTest[
  Head @ HelmholtzDecomposition[{x^2 + y, x y}, {x, y}],
  HelmholtzDecomposition,
  TestID -> "Result has head HelmholtzDecomposition"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Properties"],
  {"Properties", "Association", "Field", "Coordinates",
   "Gradient", "Rotational",
   "Potential", "AntisymmetricPotential", "PotentialMatrix",
   "VectorPotentialA",
   "GaugeShift", "GaugePotential", "CanonicalAssociation",
   "Residual", "Verify"},
  TestID -> "Properties list"
];

(* ---- property access ---- *)

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Gradient"],
  {x^2 + y^2/2, x y},
  TestID -> "[\"Gradient\"]"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Rotational"],
  {-(1/2) (-2 + y) y, 0},
  TestID -> "[\"Rotational\"]"
];

VerificationTest[
  Module[{hd, P, g},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    P = hd["Potential"]; g = hd["Gradient"];
    Simplify[Grad[P, {x, y}] - g]
  ],
  {0, 0},
  TestID -> "[\"Potential\"] is a primitive of [\"Gradient\"]"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Field"],
  {x^2 + y, x y},
  TestID -> "[\"Field\"] returns input"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Coordinates"],
  {x, y},
  TestID -> "[\"Coordinates\"] returns input"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["GaugeShift"],
  {0, 0},
  TestID -> "[\"GaugeShift\"] defaults to zero"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["GaugePotential"],
  0,
  TestID -> "[\"GaugePotential\"] defaults to zero"
];

VerificationTest[
  AssociationQ @ HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["CanonicalAssociation"],
  True,
  TestID -> "[\"CanonicalAssociation\"] is available on canonical objects"
];

(* ---- Verify ---- *)

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Verify"],
  <|"GradientCurlFree" -> True,
    "RotationalDivergenceFree" -> True,
    "ResidualZero" -> True|>,
  TestID -> "[\"Verify\"] reports all checks"
];

(* ---- 3D vector potential A ---- *)

VerificationTest[
  Module[{hd, A, r},
    hd = HelmholtzDecomposition[{-y - z, x, 0}, {x, y, z}];
    A  = hd["VectorPotentialA"];
    r  = hd["Rotational"];
    Simplify[Curl[A, {x, y, z}] - r]
  ],
  {0, 0, 0},
  TestID -> "[\"VectorPotentialA\"]: Curl[A] == r in 3D"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["VectorPotentialA"],
  $Failed,
  {HelmholtzDecomposition::baddim},
  TestID -> "[\"VectorPotentialA\"] in 2D -> $Failed with message"
];

(* ---- legacy List @@ form ---- *)

VerificationTest[
  List @@ HelmholtzDecomposition[{x^2 + y, x y}, {x, y}],
  {{x^2 + y^2/2, x y}, {-(1/2) (-2 + y) y, 0}},
  TestID -> "List @@ hd returns {g, r}"
];

VerificationTest[
  AssociationQ @ Normal @ HelmholtzDecomposition[{x^2 + y, x y}, {x, y}],
  True,
  TestID -> "Normal[hd] returns the Association"
];

(* ---- coordinate auto-detection ---- *)

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}]["Gradient"],
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Gradient"],
  TestID -> "Auto-coord agrees with explicit coords"
];

(* ---- options ---- *)

VerificationTest[
  Module[{P},
    P = HelmholtzDecomposition[{a x^2, b y}, {x, y},
          Assumptions -> a > 0 && b > 0]["Potential"];
    Simplify[Grad[P, {x, y}] - HelmholtzDecomposition[{a x^2, b y}, {x, y}]["Gradient"],
            a > 0 && b > 0]
  ],
  {0, 0},
  TestID -> "Assumptions option round-trips"
];

VerificationTest[
  Head @ HelmholtzDecomposition[{x^2 + y, x y}, {x, y},
            SimplificationFunction -> Identity],
  HelmholtzDecomposition,
  TestID -> "SimplificationFunction -> Identity returns valid object"
];

(* ---- gauge shifting ---- *)

VerificationTest[
  Module[{hd, shifted},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    shifted = HelmholtzGaugeShift[hd, x + y];
    Simplify[shifted["Gradient"] - (hd["Gradient"] + {1, 1})]
  ],
  {0, 0},
  TestID -> "HelmholtzGaugeShift[hd, phi]: gradient shifts by Grad[phi]"
];

VerificationTest[
  Module[{hd, shifted},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    shifted = HelmholtzGaugeShift[hd, x + y];
    Simplify[shifted["Rotational"] - (hd["Rotational"] - {1, 1})]
  ],
  {0, 0},
  TestID -> "HelmholtzGaugeShift[hd, phi]: rotational shifts oppositely"
];

VerificationTest[
  Module[{hd, shifted},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    shifted = HelmholtzGaugeShift[hd, x + y];
    Simplify[Grad[shifted["Potential"], {x, y}] - shifted["Gradient"]]
  ],
  {0, 0},
  TestID -> "HelmholtzGaugeShift[hd, phi]: shifted potential matches shifted gradient"
];

VerificationTest[
  Module[{hd, shifted},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    shifted = HelmholtzGaugeShift[hd, {1, 1}];
    Simplify[shifted["GaugePotential"] - (x + y)]
  ],
  0,
  TestID -> "HelmholtzGaugeShift[hd, h]: reconstructs scalar potential from harmonic vector"
];

VerificationTest[
  Module[{hd, shifted},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    shifted = HelmholtzGaugeShift[hd, x + y];
    shifted["Verify"]
  ],
  <|"GradientCurlFree" -> True,
    "RotationalDivergenceFree" -> True,
    "ResidualZero" -> True|>,
  TestID -> "HelmholtzGaugeShift preserves verification checks"
];

VerificationTest[
  Module[{hd, shifted},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    shifted = HelmholtzGaugeShift[hd, x + y];
    Simplify[shifted["CanonicalAssociation"]["Gradient"] - hd["Gradient"]]
  ],
  {0, 0},
  TestID -> "[\"CanonicalAssociation\"] preserves the original decomposition"
];

(* ---- caching ---- *)

VerificationTest[
  Module[{f = {x^3 y + y^2, x y^2 - x^2}, xv = {x, y}, t1, t2},
    ClearHelmholtzCache[];
    t1 = First @ AbsoluteTiming[HelmholtzDecomposition[f, xv]];
    t2 = First @ AbsoluteTiming[HelmholtzDecomposition[f, xv]];
    t2 < t1 / 5
  ],
  True,
  TestID -> "Caching: repeat call is much faster"
];

VerificationTest[
  Module[{f = {x^2, y^2}, xv = {x, y}, n},
    ClearHelmholtzCache[];
    HelmholtzDecomposition[f, xv];
    n = ClearHelmholtzCache[];
    n >= 1
  ],
  True,
  TestID -> "ClearHelmholtzCache returns count"
];

VerificationTest[
  Module[{f = {x^2, y^2}, xv = {x, y}},
    ClearHelmholtzCache[];
    HelmholtzDecomposition[f, xv, Caching -> False];
    HelmholtzDecomposition[f, xv, Caching -> False];
    ClearHelmholtzCache[]   (* should be 0 — nothing was stored *)
  ],
  0,
  TestID -> "Caching -> False does not write to the cache"
];

(* Cache key includes the simplifier: switching produces a fresh entry, not a stale hit *)
VerificationTest[
  Module[{f = {x^2 + y, x y}, xv = {x, y}, hdFull, hdId},
    ClearHelmholtzCache[];
    hdId = HelmholtzDecomposition[f, xv, SimplificationFunction -> Identity];
    hdFull = HelmholtzDecomposition[f, xv, SimplificationFunction -> FullSimplify];
    Length[HelmholtzDecomposition`Private`$cache] === 2
      && hdId =!= hdFull
  ],
  True,
  TestID -> "Cache key includes SimplificationFunction"
];

(* TimeConstraint is honored *)
VerificationTest[
  Head @ HelmholtzDecomposition[{x^2 + y, x y}, {x, y}, TimeConstraint -> 30],
  HelmholtzDecomposition,
  TestID -> "TimeConstraint -> finite N is accepted"
];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}, TimeConstraint -> -1],
  $Failed,
  {HelmholtzDecomposition::badtc},
  TestID -> "TimeConstraint -> negative -> $Failed"
];

(* ---- introspection ---- *)

VerificationTest[
  Length @ HelmholtzDecomposition[{x^2 + y, x y, z}, {x, y, z}]["Coordinates"],
  3,
  TestID -> "Field dimension via Length @ hd[\"Coordinates\"]"
];

VerificationTest[
  Keys @ HelmholtzDecomposition[{x^2 + y, x y}, {x, y}],
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]["Properties"],
  TestID -> "Keys[hd] equals hd[\"Properties\"]"
];

(* Regression: rational-in-parameter coefficients on linear monomial.
   Pre-fix the package threw HelmholtzDecomposition::nodec because
   PolynomialQ[x/(1+a^2), {}] returned False (auto-detect mode rejects
   1/(1+a^2) as non-polynomial in {a}, even though we'd want a treated
   as a parameter). *)
VerificationTest[
  Module[{hd, P},
    hd = HelmholtzDecomposition[{x/(1 + a^2), 0}, {x, y}];
    P = hd["Potential"];
    Simplify[Grad[P, {x, y}] - hd["Gradient"]]
  ],
  {0, 0},
  TestID -> "Linear term with rational-in-parameter coefficient"
];

VerificationTest[
  Head @ HelmholtzDecomposition[
    {(x^2 + y)/(1 + a^2), x y/(1 + b^2)}, {x, y}],
  HelmholtzDecomposition,
  TestID -> "Mixed polynomial with parameter-rational coefficients"
];

EndTestSection[]
