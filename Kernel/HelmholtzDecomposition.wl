(* ::Package:: *)
BeginPackage["HelmholtzDecomposition`"]

(*
This package provides functions to perform a Helmholtz decomposition of a vector field in n dimensions, following the approach in Richters and Glotz, "Analytical Helmholtz Decomposition of n-Dimensional Vector Fields", https://doi.org/10.1016/j.jmaa.2023.127138. The code has been taken from the corresponding Mathematica notebook (published in https://zenodo.org/records/7680297) and adapted to be used as a package.
Any intellectual property rights of the original authors are fully respected. Please cite the original paper when using this code.
*)

HelmholtzDecomposition::usage = "HelmholtzDecomposition[f, xvec] returns {g, r} where f = g + r, g is curl-free and r is divergence-free. HelmholtzDecomposition[f] auto-detects the coordinates from the symbols appearing in f.";

HelmholtzGradient::usage = "HelmholtzGradient[f, xvec] returns the curl-free (gradient) component g of the Helmholtz decomposition f = g + r.";

HelmholtzRotational::usage = "HelmholtzRotational[f, xvec] returns the divergence-free (rotational) component r of the Helmholtz decomposition f = g + r.";

HelmholtzPotential::usage = "HelmholtzPotential[f, xvec] returns the scalar potential P with Grad[P, xvec] == HelmholtzGradient[f, xvec].";

HelmholtzVectorPotential::usage = "HelmholtzVectorPotential[f, xvec] returns the antisymmetric matrix Rij such that the rotational component r satisfies r[[j]] == Sum[D[Rij[[j,i]], xvec[[i]]], i].";

HelmholtzAssociation::usage = "HelmholtzAssociation[f, xvec] returns an Association with keys \"Gradient\", \"Rotational\", \"Potential\", \"VectorPotentialMatrix\", \"PotentialMatrix\", \"Residual\", \"Coordinates\", \"Field\". Computes the underlying potential matrix once and is cached on the (f, xvec) pair.";

PrintHelmholtz::usage = "PrintHelmholtz[f, xvec] (or PrintHelmholtz[f]) prints a full diagnostic of the Helmholtz decomposition and returns the same Association as HelmholtzAssociation. PrintHelmholtz[assoc] reuses an already-computed Association without redoing any work.";

SetVerbose::usage = "SetVerbose[True|False] enables or disables verbose output of the AutocalculateFij internals.";

GetVerbose::usage = "GetVerbose[] returns the current verbose setting.";

(* Options used by all surface functions *)
Options[HelmholtzDecomposition]    = {Assumptions :> $Assumptions, Simplify -> FullSimplify};
Options[HelmholtzGradient]         = Options[HelmholtzDecomposition];
Options[HelmholtzRotational]       = Options[HelmholtzDecomposition];
Options[HelmholtzPotential]        = Options[HelmholtzDecomposition];
Options[HelmholtzVectorPotential]  = Options[HelmholtzDecomposition];
Options[HelmholtzAssociation]      = Options[HelmholtzDecomposition];
Options[PrintHelmholtz]            = Options[HelmholtzDecomposition];

Begin["`Private`"]

(* ============================================================ *)
(* Verbose toggle                                                *)
(* ============================================================ *)

$PrintInfo = False;

SetVerbose[value_?BooleanQ] := ($PrintInfo = value);
SetVerbose[value_] := (Message[SetVerbose::bool, value]; $Failed);
SetVerbose::bool = "The argument `1` is not True or False.";

GetVerbose[] := $PrintInfo;

(* ============================================================ *)
(* Configurable simplifier                                       *)
(* Internal default is FullSimplify with $Assumptions; the public *)
(* surface lets the user override both via options.              *)
(* ============================================================ *)

$simplify  = FullSimplify;
$assumptions = True;

withSimp[expr_] := $simplify[expr, $assumptions];

(* ============================================================ *)
(* Iterated integrals and incomplete Laplacians                  *)
(* ============================================================ *)

IntegrateXmP[func_, 0, m_, xvec_] := func;
IntegrateXmP[func_, p_Integer?Positive, m_, xvec_] :=
  Integrate[IntegrateXmP[func, p - 1, m, xvec], {xvec[[m]], 0, xvec[[m]]}];

LaplacianIncomplete[func_, m_, xvec_] := Laplacian[func, xvec] - D[func, {xvec[[m]], 2}];

LaplacianIncompleteP[func_, 0, m_, xvec_] := func;
LaplacianIncompleteP[func_, p_Integer?Positive, m_, xvec_] :=
  LaplacianIncomplete[LaplacianIncompleteP[func, p - 1, m, xvec], m, xvec];

(* ============================================================ *)
(* Variable extraction                                           *)
(* ============================================================ *)

headlist = {Or, And, Equal, Unequal, Less, LessEqual, Greater, GreaterEqual, Inequality};

getAllVariables[f_?NumericQ] := {};
getAllVariables[{}] := {};
getAllVariables[t_] /; MemberQ[headlist, t] := {};
getAllVariables[ll_List] := Flatten[Union[Map[getAllVariables[#] &, ll]]];
getAllVariables[Derivative[_Integer][_][arg__]] := getAllVariables[{arg}];
getAllVariables[f_Symbol[arg__]] :=
  If[MemberQ[Attributes[f], NumericFunction] || MemberQ[headlist, f],
    getAllVariables[{arg}],
    f[arg]
  ];
getAllVariables[other_] := other;
getUniqueVariables[all_] := Sort[DeleteDuplicates[Flatten[List[getAllVariables[all]]]]];

(* Auto-detect coordinate vector. Used when xvec is omitted. *)
autoDetectCoordinates[f_List] :=
  Module[{vars = getUniqueVariables[f]},
    If[Length[vars] =!= Length[f],
      Message[HelmholtzDecomposition::autocoord, Length[f], vars];
      $Failed,
      vars
    ]
  ];

HelmholtzDecomposition::autocoord =
  "Cannot auto-detect coordinates: vector field has length `1` but the symbols `2` were found in the expression. Pass an explicit coordinate list as the second argument.";

(* ============================================================ *)
(* Potential matrix Fij (Theorem 6.1 of Richters & Glotzl)       *)
(* ============================================================ *)

CalculateFij[f_, k_, m_, lambda_, u_, W_, xvec_, dim_] :=
  Module[{Fij, i, p, j},
    On[Assert];
    Assert[Length[f] == Length[xvec] == dim, "Length[f], Length[xvec] and dim must be identical"];
    Assert[1 <= k <= dim, "k must be between 1 and dim"];
    Assert[1 <= m <= dim, "m must be between 1 and dim"];
    If[k != m,
      Do[Assert[D[f[[k]], xvec[[i]]] === 0 || i == k || i == m], {i, dim}]
    ];
    Assert[
      FullSimplify[D[W, {xvec[[m]], 2 lambda}] - f[[k]]] == 0,
      "Differentiating W 2*lambda times wrt x_m must yield f_k"
    ];
    Assert[
      FullSimplify[(-1)^lambda LaplacianIncompleteP[W, lambda, m, xvec] - (1 - u) f[[k]] == 0],
      "Terminal condition not fulfilled"
    ];
    Fij = Simplify[Table[
      KroneckerDelta[k, i] * Sum[
        (-1)^p / u * D[LaplacianIncompleteP[W, p, m, xvec], xvec[[j]], {xvec[[m]], 2 lambda - 2 p - 2}],
        {p, 0, lambda - 1}
      ],
      {i, dim}, {j, dim}
    ]];
    Fij
  ];

AutocalculateSingleFij[S_, k_, xvec_, dim_] :=
  Module[{i, j, Singlef, xiWithoutxk, Q1, Q2, Q3, Q4, m, lambda, u, W,
          v1, v2, xiInv1, xiInv2},
    If[S === 0, Return[Table[0, {i, dim}, {j, dim}]]];
    Singlef = Table[S KroneckerDelta[i, k], {i, dim}];
    xiWithoutxk = Intersection[xvec, DeleteCases[getUniqueVariables[S], xvec[[k]]]];
    Q1 = PolynomialQ[S, xiWithoutxk];
    Q2 = PolynomialQ[S, xvec] && Length[xiWithoutxk] == 1 &&
         Exponent[S, First[xiWithoutxk]] > Exponent[S, xvec[[k]]];
    Q3 = PolynomialQ[S, xvec[[k]]] && Length[xiWithoutxk] == 1;
    If[Q1 || Q3,
      m = If[Q1 && !Q2, k, First[Flatten[Position[xvec, First[xiWithoutxk]]]]];
      u = 1;
      lambda = Ceiling[(Total[Exponent[S, xvec]] - Exponent[S, xvec[[m]]] + 1) / 2];
      W = IntegrateXmP[S, 2 lambda, m, xvec];
      If[$PrintInfo, Print[S, " in f", k, " is partly monomial, use m=", m, ", lambda=", lambda]];
      Return[CalculateFij[Singlef, k, m, lambda, u, W, xvec, dim]]
    ];
    v1 = FullSimplify[D[S, {xvec[[k]], 2}] / S];
    v2 = FullSimplify[LaplacianIncomplete[S, k, xvec] / S];
    xiInv1 = Intersection[xvec, getUniqueVariables[v1]];
    xiInv2 = Intersection[xvec, getUniqueVariables[v2]];
    Q4 = !(v1 === 0) && Length[xiInv1] == 0 && Length[xiInv2] == 0;
    If[Q4,
      lambda = 1;
      m = k;
      W = S / v1;
      u = 1 + v2 / v1;
      If[$PrintInfo, Print[S, " in f", k, " can be solved with m=", m, ", lambda=", lambda, ", W=", W, ", u=", u]];
      Return[CalculateFij[Singlef, k, m, lambda, u, W, xvec, dim]]
    ];
    Message[HelmholtzDecomposition::nodec, S, k];
    Throw[$Failed, HelmholtzDecomposition]
  ];

HelmholtzDecomposition::nodec =
  "No analytic decomposition found for term `1` in component f_`2`. The Richters-Glotzl construction handles polynomial, exponential, sine and cosine atoms; mixed transcendentals may need to be split by hand.";

AutocalculateFij[f_, xvec_, dim_] :=
  Module[{M, Fij, Poly, pol, k},
    Fij = Table[0, {dim}, {dim}];
    Do[
      Poly = Expand[f[[k]]];
      Poly = If[Head[Poly] === Plus, List @@ Poly, {Poly}];
      If[$PrintInfo, Print["k=", k, ": ", Poly]];
      Do[
        M = Poly[[pol]];
        Fij = Fij + AutocalculateSingleFij[M, k, xvec, dim],
        {pol, Length[Poly]}
      ],
      {k, dim}
    ];
    Fij
  ];

(* ============================================================ *)
(* Pieces extracted from Fij                                     *)
(* ============================================================ *)

ROTij[v_, xvec_, dim_] := withSimp[
  Table[D[v[[i]], xvec[[j]]] - D[v[[j]], xvec[[i]]], {i, dim}, {j, dim}]
];

scalarPotentialFromFij[Fij_, dim_] := withSimp[Sum[Fij[[k, k]], {k, dim}]];

antisymmetricFromFij[Fij_] := withSimp[Fij - Transpose[Fij]];

(* ============================================================ *)
(* Single-shot data builder, memoized on (f, xvec, opts)         *)
(* ============================================================ *)

ClearAll[iHelmholtzData];

iHelmholtzData[f_List, xvec_List, simp_, asm_] :=
  iHelmholtzData[f, xvec, simp, asm] =
    Module[{dim, Fij, P, Rij, g, r, residual, $simplify = simp, $assumptions = asm},
      On[Assert];
      Assert[Length[f] == Length[xvec],
        "Length of vector field f must match length of coordinate vector xvec"];
      dim = Length[f];
      Fij = Catch[AutocalculateFij[f, xvec, dim], HelmholtzDecomposition];
      If[Fij === $Failed, Return[$Failed]];
      P = scalarPotentialFromFij[Fij, dim];
      Rij = antisymmetricFromFij[Fij];
      g = withSimp[Grad[P, xvec]];
      r = withSimp[Table[Sum[D[Rij[[j, i]], xvec[[i]]], {i, dim}], {j, dim}]];
      residual = withSimp[f - g - r];
      <|
        "Gradient"              -> g,
        "Rotational"            -> r,
        "Potential"             -> P,
        "VectorPotentialMatrix" -> Rij,
        "PotentialMatrix"       -> Fij,
        "Residual"              -> residual,
        "Coordinates"           -> xvec,
        "Field"                 -> f
      |>
    ];

(* Resolve options into (simp, asm) *)
resolveOpts[opts_List, headSym_] := {
  OptionValue[headSym, opts, Simplify],
  OptionValue[headSym, opts, Assumptions]
};

(* ============================================================ *)
(* Public API                                                    *)
(* ============================================================ *)

HelmholtzAssociation[f_List, xvec_List, opts:OptionsPattern[]] :=
  Module[{simp, asm},
    {simp, asm} = resolveOpts[{opts}, HelmholtzAssociation];
    iHelmholtzData[f, xvec, simp, asm]
  ];

HelmholtzDecomposition[f_List, xvec_List, opts:OptionsPattern[]] :=
  Module[{d = HelmholtzAssociation[f, xvec, opts]},
    If[AssociationQ[d], {d["Gradient"], d["Rotational"]}, $Failed]
  ];

HelmholtzGradient[f_List, xvec_List, opts:OptionsPattern[]] :=
  Replace[HelmholtzAssociation[f, xvec, opts], a_Association :> a["Gradient"]];

HelmholtzRotational[f_List, xvec_List, opts:OptionsPattern[]] :=
  Replace[HelmholtzAssociation[f, xvec, opts], a_Association :> a["Rotational"]];

HelmholtzPotential[f_List, xvec_List, opts:OptionsPattern[]] :=
  Replace[HelmholtzAssociation[f, xvec, opts], a_Association :> a["Potential"]];

HelmholtzVectorPotential[f_List, xvec_List, opts:OptionsPattern[]] :=
  Replace[HelmholtzAssociation[f, xvec, opts], a_Association :> a["VectorPotentialMatrix"]];

(* Auto-detect coordinates when xvec is omitted *)
autoCall[head_, f_List, opts___] :=
  Module[{xv = autoDetectCoordinates[f]},
    If[xv === $Failed, $Failed, head[f, xv, opts]]
  ];

HelmholtzDecomposition[f_List, opts:OptionsPattern[]]      := autoCall[HelmholtzDecomposition,     f, opts];
HelmholtzGradient[f_List, opts:OptionsPattern[]]           := autoCall[HelmholtzGradient,          f, opts];
HelmholtzRotational[f_List, opts:OptionsPattern[]]         := autoCall[HelmholtzRotational,        f, opts];
HelmholtzPotential[f_List, opts:OptionsPattern[]]          := autoCall[HelmholtzPotential,         f, opts];
HelmholtzVectorPotential[f_List, opts:OptionsPattern[]]    := autoCall[HelmholtzVectorPotential,   f, opts];
HelmholtzAssociation[f_List, opts:OptionsPattern[]]        := autoCall[HelmholtzAssociation,       f, opts];

(* ============================================================ *)
(* Diagnostic printer — returns the Association                  *)
(* ============================================================ *)

printAssociation[data_Association] :=
  Module[{xvec = data["Coordinates"], f = data["Field"], dim},
    dim = Length[f];
    Print["f   = ", MatrixForm[f]];
    Print["Div f = ", Div[f, xvec]];
    Print[OverBar["ROT"], " f = ", MatrixForm[ROTij[f, xvec, dim], TableSpacing -> {1, 3}]];
    Print["Fij = ", MatrixForm[data["PotentialMatrix"], TableSpacing -> {1, 4}]];
    Print["G (scalar potential) = ", data["Potential"]];
    Print["Rij = ", MatrixForm[data["VectorPotentialMatrix"], TableSpacing -> {1, 4}]];
    Print["g   = ", MatrixForm[data["Gradient"]]];
    Print["r   = ", MatrixForm[data["Rotational"]]];
    Print[OverBar["ROT"], " g = ", MatrixForm[ROTij[data["Gradient"], xvec, dim]]];
    Print["Div r = ", FullSimplify[Div[data["Rotational"], xvec]]];
    Print["f - g - r = ", MatrixForm[data["Residual"]]];
    Print["f = g + r: ", data["Residual"] === ConstantArray[0, dim]];
    data
  ];

PrintHelmholtz[f_List, xvec_List, opts:OptionsPattern[]] :=
  Replace[HelmholtzAssociation[f, xvec, opts], a_Association :> printAssociation[a]];

PrintHelmholtz[f_List, opts:OptionsPattern[]] :=
  autoCall[PrintHelmholtz, f, opts];

PrintHelmholtz[data_Association] /;
  SubsetQ[Keys[data], {"Gradient", "Rotational", "Potential",
    "VectorPotentialMatrix", "PotentialMatrix", "Residual",
    "Coordinates", "Field"}] := printAssociation[data];

End[] (* End `Private` *)

EndPackage[]
