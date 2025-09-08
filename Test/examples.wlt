BeginTestSection["Examples"];
Needs["HelmholtzDecomposition`"];

VerificationTest[
  HelmholtzDecomposition[{x^2 + y, x y}, {x, y}],
  {{x^2 + y^2/2, x y},
  {-(1/2) (-2 + y) y, 0}},
  TestID -> "2d polynomial field"
  ]

VerificationTest[
  HelmholtzDecomposition[{Cos["w" x] Exp["a" y], 0}, {x, y}],
{{(("w")^2 
\!\(\*SuperscriptBox[\(\[ExponentialE]\), \("a"\ y\)]\) Cos["w" x])/(-(
    "a")^2 + ("w")^2), ("a" "w" 
\!\(\*SuperscriptBox[\(\[ExponentialE]\), \("a"\ y\)]\) Sin["w" x])/(-(
    "a")^2 + ("w")^2)}, {(("a")^2 
\!\(\*SuperscriptBox[\(\[ExponentialE]\), \("a"\ y\)]\) Cos["w" x])/((
   "a")^2 - ("w")^2), ("a" "w" 
\!\(\*SuperscriptBox[\(\[ExponentialE]\), \("a"\ y\)]\) Sin["w" x])/((
   "a")^2 - ("w")^2)}},
  TestID -> "2d Cos Exp field"
  ]

VerificationTest[
  HelmholtzDecomposition[{-y - z, x + "a" y, "b" + (x - "c") z}, {x, y, z}], 
  {{z^2/2, "a" y, "b" + (-"c" + x) z}, {-y - 1/2 z (2 + z), x, 0}},
  TestID -> "Rössler attractor"
  ]

VerificationTest[
  HelmholtzDecomposition[{"a" (y - x), x ("b" - z) - y, x y - "c" z}, {x, y, z}], 
  {{-"a" x, -y, -"c" z}, {"a" y, x ("b" - z), x y}},
  TestID -> "Lorenz Attractor"
  ]

dim = 3;
xvec = Array[Subscript[x, #] &, dim];
rvec = Array[Subscript["b", #] &, dim];
\[Alpha]ij = Array[Subscript["a", #1, #2] &, {dim, dim}];
f = Table[rvec[[i]] Subscript[x, i] (1 - Sum[\[Alpha]ij[[i, j]] Subscript[x, j], {j, dim}]), {i, dim}];

LotkaVolterraResult= {{1/2 (-2 Subscript["b", 1] Subscript[x, 
      1] (-1 + Subscript[x, 1] Subscript["a", 1, 1] + 
        Subscript[x, 2] Subscript["a", 1, 2] + 
        Subscript[x, 3] Subscript["a", 1, 3]) - Subscript["b", 2] 
\!\(\*SubsuperscriptBox[\(x\), \(2\), \(2\)]\) Subscript["a", 2, 1] - 
     Subscript["b", 3] 
\!\(\*SubsuperscriptBox[\(x\), \(3\), \(2\)]\) Subscript["a", 3, 1]), 
  1/2 (-Subscript["b", 1] 
\!\(\*SubsuperscriptBox[\(x\), \(1\), \(2\)]\) Subscript["a", 1, 2] - 
     2 Subscript["b", 2] Subscript[x, 
      2] (-1 + Subscript[x, 1] Subscript["a", 2, 1] + 
        Subscript[x, 2] Subscript["a", 2, 2] + 
        Subscript[x, 3] Subscript["a", 2, 3]) - Subscript["b", 3] 
\!\(\*SubsuperscriptBox[\(x\), \(3\), \(2\)]\) Subscript["a", 3, 2]), 
  1/2 (-Subscript["b", 1] 
\!\(\*SubsuperscriptBox[\(x\), \(1\), \(2\)]\) Subscript["a", 1, 3] - 
     Subscript["b", 2] 
\!\(\*SubsuperscriptBox[\(x\), \(2\), \(2\)]\) Subscript["a", 2, 3] - 
     2 Subscript["b", 3] Subscript[x, 
      3] (-1 + Subscript[x, 1] Subscript["a", 3, 1] + 
        Subscript[x, 2] Subscript["a", 3, 2] + 
        Subscript[x, 3] Subscript["a", 3, 3]))}, {1/
   2 (Subscript["b", 2] 
\!\(\*SubsuperscriptBox[\(x\), \(2\), \(2\)]\) Subscript["a", 2, 1] + 
     Subscript["b", 3] 
\!\(\*SubsuperscriptBox[\(x\), \(3\), \(2\)]\) Subscript["a", 3, 1]), 
  1/2 (Subscript["b", 1] 
\!\(\*SubsuperscriptBox[\(x\), \(1\), \(2\)]\) Subscript["a", 1, 2] + 
     Subscript["b", 3] 
\!\(\*SubsuperscriptBox[\(x\), \(3\), \(2\)]\) Subscript["a", 3, 2]), 
  1/2 (Subscript["b", 1] 
\!\(\*SubsuperscriptBox[\(x\), \(1\), \(2\)]\) Subscript["a", 1, 3] + 
     Subscript["b", 2] 
\!\(\*SubsuperscriptBox[\(x\), \(2\), \(2\)]\) Subscript["a", 2, 3])}}
VerificationTest[
  HelmholtzDecomposition[f, xvec], 
LotkaVolterraResult,
  TestID -> "Competitive Lotka-Volterra equations with n species" ]

EndTestSection[]
