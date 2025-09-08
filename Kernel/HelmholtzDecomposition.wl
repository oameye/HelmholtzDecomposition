(* ::Package:: *)
BeginPackage["HelmholtzDecomposition`"]

(*
This package provides functions to perform a Helmholtz decomposition of a vector field in n dimensions, following the approach in Richters and Glotz, "Analytical Helmholtz Decomposition of n-Dimensional Vector Fields", https://doi.org/10.1016/j.jmaa.2023.127138. The code has been taken from the corresponding Mathematica notebook (published in https://zenodo.org/records/7680297) and adapted to be used as a package.
Any intellectual property rights of the original authors are fully respected. Please cite the original paper when using this code.
*)

(* Exported symbols added here with SymbolName::usage *)

PrintHelmholtz::usage = "PrintHelmholtz[f, xvec] performs and displays the Helmholtz decomposition of vector field f with coordinate variables xvec. PrintHelmholtz[f, xvec, options] allows additional options.";

HelmholtzGradient::usage = "HelmholtzGradient[f, xvec] returns the gradient component g of the Helmholtz decomposition f = g + r.";

HelmholtzRotational::usage = "HelmholtzRotational[f, xvec] returns the rotational component r of the Helmholtz decomposition f = g + r.";

HelmholtzDecomposition::usage = "HelmholtzDecomposition[f, xvec] returns {g, r} where f = g + r is the Helmholtz decomposition.";

SetVerbose::usage = "SetVerbose[True/False] enables or disables verbose output during decomposition calculations.";

GetVerbose::usage = "GetVerbose[] returns the current verbose setting.";

Begin["`Private`"]

(* Global verbose setting *)
$PrintInfo = False;

(* Functions to control verbose output *)
SetVerbose[value_?BooleanQ] := ($PrintInfo = value);
SetVerbose[value_] := (Message[SetVerbose::bool, value]; $Failed);
SetVerbose::bool = "The argument `1` is not True or False.";

GetVerbose[] := $PrintInfo;

(* Integrals and Incomplete Laplacians *)

(* Riemann-Liouville Integral, integrate func p-times with respect to xvec[[m]] *)
IntegrateXmP[func_, 0, m_, xvec_] := func;
IntegrateXmP[func_, p_Integer?Positive, m_, xvec_] := 
  Integrate[IntegrateXmP[func, p-1, m, xvec], {xvec[[m]], 0, xvec[[m]]}];

(* Incomplete Laplacian: Subscript[\[CapitalDelta], \m] *)
LaplacianIncomplete[func_, m_, xvec_] := Laplacian[func, xvec] - D[func, {xvec[[m]], 2}];

(* Incomplete Laplacian to the power of p: Subscript[\[CapitalDelta], \m]^p *)
LaplacianIncompleteP[func_, 0, m_, xvec_] := func;
LaplacianIncompleteP[func_, p_Integer?Positive, m_, xvec_] := 
  LaplacianIncomplete[LaplacianIncompleteP[func, p-1, m, xvec], m, xvec];

(* Get all variables *)

(* helper function getAllVariables from https://mathematica.stackexchange.com/questions/21257/extracting-variables-from-an-expression *)
headlist = {Or, And, Equal, Unequal, Less, LessEqual, Greater, GreaterEqual, Inequality};

getAllVariables[f_?NumericQ] := List[]
getAllVariables[{}] := List[]
getAllVariables[t_] /; MemberQ[headlist, t] := List[]
getAllVariables[ll_List] := Flatten[Union[Map[getAllVariables[#] &, ll]]]
getAllVariables[Derivative[n_Integer][f_][arg__]] := getAllVariables[{arg}]
getAllVariables[f_Symbol[arg__]] := Module[{fvars},
  If[MemberQ[Attributes[f], NumericFunction] || MemberQ[headlist, f],
    fvars = getAllVariables[{arg}],
    (* else *)
    fvars = f[arg]
  ];
  fvars
]
getAllVariables[other_] := other
getUniqueVariables[all_] := Sort[DeleteDuplicates[Flatten[List[getAllVariables[all]]]]]

(* Calculate Potential Matrix *)
(* === CalculateFij returns Fij with all the inputs manually given === *)
CalculateFij[f_, k_, m_, \[Lambda]_, u_, W_, xvec_, dim_] := Module[{Fij},
  (* Check whether conditions of Theorem 6.1 are fulfilled *)
  On[Assert];
  Assert[Length[f] == Length[xvec] == dim, "Length[f], Length[xvec] and dim must be identical"];
  Assert[1 <= k <= dim, "k must be between 1 and dim"];
  Assert[1 <= m <= dim, "m must be between 1 and dim"];
  Assert[u != 0 || u == 1 - "a"^2/"w"^2 (* added to avoid assertion errors *), "u must be non-zero"];
  If[k != m, 
    For[i = 1, i <= dim, i++, 
      Assert[D[f[[k]], xvec[[i]]] === 0 || i == k || i == m]
    ]
  ];
  Assert[FullSimplify[D[W, {xvec[[m]], 2\[Lambda]}] - f[[k]]] == 0,
    StringJoin["Differentiating W=,", W, " with respect to \!\(\*SubscriptBox[\(x\), \(m\)]\) with m=", m, " 2\[Lambda]-times with \[Lambda]=,", \[Lambda], " must yield \!\(\*SubscriptBox[\(f\), \(k\)]\)=", f[[k]]]
  ];
  Assert[FullSimplify[(-1)^\[Lambda] LaplacianIncompleteP[W, \[Lambda], m, xvec] - (1 - u) f[[k]] == 0], "Terminal condition not fulfilled"];
  (* Calculate Fij *)
  Fij = Simplify[Table[
    KroneckerDelta[k, i]*Sum[
      (-1)^p/u*D[LaplacianIncompleteP[W, p, m, xvec], xvec[[j]], {xvec[[m]], 2\[Lambda] - 2p - 2}], 
      {p, 0, \[Lambda] - 1}
    ], 
    {i, dim}, {j, dim}
  ]];
  Return[Fij]
]; (* end of function CalculateFij *)

(* === AutocalculateSingleFij returns Fij given a simple expression (no sums!) M and the vector coordinate k === *)
(* === Don't call this function directly, use AutocalculateFij[f] === *)
AutocalculateSingleFij[S_, k_, xvec_, dim_] := Module[{i, j, Singlef, xiWithoutxk, Q1, Q2, Q3, Q4, m, \[Lambda], u, W, v1, v2, Fij},
  If[S === 0, Return[Table[0, {i, dim}, {j, dim}]]];
  Singlef = Table[S KroneckerDelta[i, k], {i, dim}];
  xiWithoutxk = Intersection[xvec, DeleteCases[getUniqueVariables[S], xvec[[k]]]];
  (* find suitable monomial to differentiate *)
  Q1 = PolynomialQ[S, xiWithoutxk];
  Q2 = PolynomialQ[S, xvec] && Length[xiWithoutxk] == 1 && Exponent[S, First[xiWithoutxk]] > Exponent[S, xvec[[k]]];
  Q3 = PolynomialQ[S, xvec[[k]]] && Length[xiWithoutxk] == 1;
  If[Q1 || Q3,
    m = If[Q1 && !Q2, k, (* i *) First[Flatten[Position[xvec, First[xiWithoutxk]]]]];
    u = 1;
    \[Lambda] = Ceiling[(Total[Exponent[S, xvec]] - Exponent[S, xvec[[m]]] + 1)/2]; (* nochmal checken ob M ein polynom ist *)
    W = IntegrateXmP[S, 2\[Lambda], m, xvec];
    If[$PrintInfo, Print[S, " in f", k, " is partly monomial, use m=", m, ", \[Lambda]=", \[Lambda]]];
    Fij = CalculateFij[Singlef, k, m, \[Lambda], u, W, xvec, dim];
    Return[Fij]
  ];
  (* try to catch functions exp, sin, cos *)
  v1 = FullSimplify[D[S, {xvec[[k]], 2}]/S];
  v2 = FullSimplify[LaplacianIncomplete[S, k, xvec]/S];
  xiInv1 = Intersection[xvec, getUniqueVariables[v1]];
  xiInv2 = Intersection[xvec, getUniqueVariables[v2]];
  Q4 = !v1 === 0 && Length[xiInv1] == 0 && Length[xiInv2] == 0;
  If[Q4,
    \[Lambda] = 1;
    m = k;
    W = S/v1;
    u = 1 + v2/v1;
    If[$PrintInfo, Print[S, " in f", k, " can be solved with m=", m, ", \[Lambda]=", \[Lambda], ", W=", W, ", u=", u]];
    Fij = CalculateFij[Singlef, k, m, \[Lambda], u, W, xvec, dim];
    Return[Fij]
  ];
  Print["No decomposition found for expression ", S, " in f", k, "."];
  Abort[];
]; (* end of function AutocalculateSingleFij *)

(* === AutocalculateFij[f] returns Fij for a given vector field f, if possible === *)
AutocalculateFij[f_, xvec_, dim_] := Module[{M, Fij, Poly, pol, k, i, j},
  Fij = Table[0, {i, dim}, {j, dim}]; (* initialize with 0 *)
  Do[
    Poly = Expand[f[[k]]];
    Poly = If[Head[Poly] === Plus, List @@ Poly, List[Poly]]; (* Make sure Poly is a proper list, even if f[[k]] has only one summand *)
    If[$PrintInfo, Print["k=", k, ": ", Poly]];
    Do[
      M = Poly[[pol]];
      Fij = Fij + AutocalculateSingleFij[M, k, xvec, dim];
      , {pol, Length[Poly]}
    ];
    , {k, dim}
  ];
  Return[Fij]
]; (* end of function AutocalculateFij *)

ROTij[v_, xvec_, dim_] := FullSimplify[Table[D[v[[i]], xvec[[j]]] - D[v[[j]], xvec[[i]]], {i, dim}, {j, dim}]];
GF[Fij_, dim_] := FullSimplify[Sum[Fij[[k, k]], {k, dim}]];
RijF[Fij_] := FullSimplify[Fij - Transpose[Fij]];
gF[Fij_, xvec_, dim_] := FullSimplify[Grad[GF[Fij, dim], xvec]];
rF[Fij_, xvec_, dim_] := FullSimplify[Table[Sum[D[RijF[Fij][[j, i]], xvec[[i]]], {i, dim}], {j, dim}]];

(* PrintHelmholtz *)

PrintHelmholtz[f_, xvec_, Fij_, SF_] := Module[{PrintCurl, G, Rij, g, r, dim},
  dim = Length[f];
  PrintCurl = False;
  Print["f = ", MatrixForm[f]];
  G = GF[Fij, dim];
  Rij = RijF[Fij];
  g = gF[Fij, xvec, dim];
  r = rF[Fij, xvec, dim];
  Print["Div f = ", Div[f, xvec]];
  If[dim < 4 && PrintCurl, Print["Curl f = ", Curl[f, xvec]]];
  Print[OverBar["ROT"], " f = ", If[MemberQ[SF, "ROTf"], ROTij[f, xvec, dim], MatrixForm[ROTij[f, xvec, dim], TableSpacing -> {1, 3}]]];
  Print["Fij = ", If[MemberQ[SF, "Fij"], Fij, MatrixForm[Fij, TableSpacing -> {1, 4}]]];
  Print["G = ", G];
  Print["Rij = ", If[MemberQ[SF, "Rij"], Rij, MatrixForm[Rij, TableSpacing -> {1, 4}]]];
  Print["g = ", If[MemberQ[SF, "g"], g, MatrixForm[g]]];
  Print["r = ", If[MemberQ[SF, "r"], r, MatrixForm[r]]];
  If[dim < 4 && PrintCurl, Print["Curl g = ", Curl[g, xvec]]];
  Print[OverBar["ROT"], " g = ", MatrixForm[ROTij[g, xvec, dim]]];
  Print["Div r = ", FullSimplify[Div[r, xvec]]];
  Assert[FullSimplify[f - g - r] === Array[0 &, dim], "f not equal g + r"];
  Print["g + r = ", MatrixForm[Simplify[g + r]]];
  Print["f = g + r: ", FullSimplify[Expand[f]] == FullSimplify[Expand[g + r]]];
];

PrintHelmholtz[f_, xvec_, Fij_] := Module[{dim},
  dim = Length[f];
  If[SubsetQ[{"ROTf", "Fij", "Rij", "g", "r"}, Fij],
    PrintHelmholtz[f, xvec, AutocalculateFij[f, xvec, dim], Fij],
    PrintHelmholtz[f, xvec, Fij, {}]
  ];
];

PrintHelmholtz[f_, xvec_] := Module[{dim},
  dim = Length[f];
  On[Assert];
  Assert[Length[f] == Length[xvec] == dim, "Length of vector field f must match length of coordinate vector xvec"];
  PrintHelmholtz[f, xvec, AutocalculateFij[f, xvec, dim], {}];
];

(* Backward compatibility: if only vector field is provided, use default coordinate names *)
PrintHelmholtz[f_] := Module[{xvec, dim},
  dim = Length[f];
  xvec = Table[Subscript[x, i], {i, dim}];
  PrintHelmholtz[f, xvec];
];

(* Additional utility functions for direct access to components *)
HelmholtzGradient[f_, xvec_] := Module[{dim, Fij},
  dim = Length[f];
  Fij = AutocalculateFij[f, xvec, dim];
  gF[Fij, xvec, dim]
];

HelmholtzRotational[f_, xvec_] := Module[{dim, Fij},
  dim = Length[f];
  Fij = AutocalculateFij[f, xvec, dim];
  rF[Fij, xvec, dim]
];

HelmholtzDecomposition[f_, xvec_] := Module[{dim, Fij, g, r},
  dim = Length[f];
  Fij = AutocalculateFij[f, xvec, dim];
  g = gF[Fij, xvec, dim];
  r = rF[Fij, xvec, dim];
  {g, r}
];

End[] (* End `Private` *)

EndPackage[]
