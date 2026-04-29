BeginTestSection["Examples"];
Needs["HelmholtzDecomposition`"];

(* These tests assert the algorithmic answer matches the paper, robust to
   FullSimplify's canonical-form drift across Mathematica versions:
   Simplify[actual - expected] === 0 instead of literal === . *)

zeroQ[expr_] := PossibleZeroQ[expr] || (Simplify[expr] === 0);
zerosQ[list_List] := AllTrue[Flatten[list], zeroQ];
zerosQ[other_] := zeroQ[other];

(* 2D polynomial *)
VerificationTest[
  Module[{hd, expectedG, expectedR},
    hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}];
    expectedG = {x^2 + y^2/2, x y};
    expectedR = {-(1/2) (-2 + y) y, 0};
    zerosQ[hd["Gradient"] - expectedG] && zerosQ[hd["Rotational"] - expectedR]
      && AllTrue[hd["Verify"], TrueQ]
  ],
  True,
  TestID -> "2d polynomial field"
];

(* 2D Cos*Exp *)
VerificationTest[
  Module[{hd, w, a, expectedG, expectedR},
    hd = HelmholtzDecomposition[{Cos[w x] Exp[a y], 0}, {x, y}];
    expectedG = {(w^2 Exp[a y] Cos[w x])/(-a^2 + w^2),
                 (a w Exp[a y] Sin[w x])/(-a^2 + w^2)};
    expectedR = {(a^2 Exp[a y] Cos[w x])/(a^2 - w^2),
                 (a w Exp[a y] Sin[w x])/(a^2 - w^2)};
    zerosQ[hd["Gradient"] - expectedG] && zerosQ[hd["Rotational"] - expectedR]
      && AllTrue[hd["Verify"], TrueQ]
  ],
  True,
  TestID -> "2d Cos*Exp field"
];

(* Rossler-style 3D field *)
VerificationTest[
  Module[{hd, a, b, c, expectedG, expectedR},
    hd = HelmholtzDecomposition[{-y - z, x + a y, b + (x - c) z}, {x, y, z}];
    expectedG = {z^2/2, a y, b + (-c + x) z};
    expectedR = {-y - 1/2 z (2 + z), x, 0};
    zerosQ[hd["Gradient"] - expectedG] && zerosQ[hd["Rotational"] - expectedR]
      && AllTrue[hd["Verify"], TrueQ]
  ],
  True,
  TestID -> "Rossler attractor"
];

(* Lorenz attractor *)
VerificationTest[
  Module[{hd, a, b, c, expectedG, expectedR},
    hd = HelmholtzDecomposition[{a (y - x), x (b - z) - y, x y - c z}, {x, y, z}];
    expectedG = {-a x, -y, -c z};
    expectedR = {a y, x (b - z), x y};
    zerosQ[hd["Gradient"] - expectedG] && zerosQ[hd["Rotational"] - expectedR]
      && AllTrue[hd["Verify"], TrueQ]
  ],
  True,
  TestID -> "Lorenz attractor"
];

(* Competitive Lotka-Volterra equations with n species *)
VerificationTest[
  Module[{dim, xvec, rvec, alphaij, f, hd, expG, expR},
    dim = 3;
    xvec = Array[Subscript[x, #] &, dim];
    rvec = Array[Subscript[bb, #] &, dim];
    alphaij = Array[Subscript[aa, #1, #2] &, {dim, dim}];
    f = Table[
      rvec[[i]] xvec[[i]] (1 - Sum[alphaij[[i, j]] xvec[[j]], {j, dim}]),
      {i, dim}
    ];
    hd = HelmholtzDecomposition[f, xvec];
    (* Just assert verification + that the components add to f *)
    zerosQ[hd["Gradient"] + hd["Rotational"] - f]
      && AllTrue[hd["Verify"], TrueQ]
  ],
  True,
  TestID -> "Competitive Lotka-Volterra (n=3) decomposes consistently"
];

(* 1D edge case: trivially gradient *)
VerificationTest[
  Module[{hd},
    hd = HelmholtzDecomposition[{x^3 + 2 x}, {x}];
    zerosQ[hd["Gradient"] - {x^3 + 2 x}] && zerosQ[hd["Rotational"]]
  ],
  True,
  TestID -> "1D field is purely gradient"
];

EndTestSection[]
