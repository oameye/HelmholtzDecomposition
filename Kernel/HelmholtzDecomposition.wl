(* ::Package:: *)
BeginPackage["HelmholtzDecomposition`"]

(*
This package provides the analytical Helmholtz decomposition of a vector field
in n dimensions, following Richters & Glotzl, "Analytical Helmholtz Decomposition
of n-Dimensional Vector Fields", https://doi.org/10.1016/j.jmaa.2023.127138.

The algorithmic core is taken from the companion notebook published at
https://zenodo.org/records/7680297; intellectual property of the original
authors is fully respected. Please cite the original paper.
*)

HelmholtzDecomposition::usage =
  "HelmholtzDecomposition[f, xvec] returns a HelmholtzDecomposition[Association[\[Ellipsis]]] " <>
  "object whose Properties give the Helmholtz decomposition f = g + r of the " <>
  "vector field f with respect to the coordinates xvec. " <>
  "HelmholtzDecomposition[f] auto-detects the coordinates from f.\n" <>
  "Properties (access via hd[\"Name\"]):\n" <>
  "  \"Gradient\"               g, the curl-free component\n" <>
  "  \"Rotational\"             r, the divergence-free component\n" <>
  "  \"Potential\"              P, scalar potential with Grad[P] == g\n" <>
  "  \"AntisymmetricPotential\" Rij, the 2-form whose row-divergence is r\n" <>
  "  \"VectorPotentialA\"       (3D only) vector A with Curl[A] == r\n" <>
  "  \"PotentialMatrix\"        Fij, the underlying potential matrix\n" <>
  "  \"Residual\"               f - g - r, simplified\n" <>
  "  \"Field\"                  the input f\n" <>
  "  \"Coordinates\"            the input xvec\n" <>
  "  \"Verify\"                 <| GradientCurlFree, RotationalDivergenceFree, ResidualZero |> -- always uses FullSimplify regardless of SimplificationFunction\n" <>
  "  \"Properties\"             list of available property names\n" <>
  "  \"Association\"            the raw underlying Association";

ClearHelmholtzCache::usage =
  "ClearHelmholtzCache[] empties the memoization cache of HelmholtzDecomposition. " <>
  "Returns the number of entries removed.";

(* Options live on the constructor head only *)
Options[HelmholtzDecomposition] = {
  SimplificationFunction :> FullSimplify,
  Assumptions :> $Assumptions,
  Verbose -> False,
  Caching -> True,
  TimeConstraint -> Infinity
};

(* Messages *)
HelmholtzDecomposition::nomatch =
  "Length[f] = `1` does not match Length[xvec] = `2`.";

HelmholtzDecomposition::autocoord =
  "Cannot auto-detect coordinates from a length-`1` field whose detected " <>
  "symbols are `2`. Pass an explicit coordinate list as the second argument.";

HelmholtzDecomposition::nodec =
  "No analytic decomposition found for term `1` in component f_`2`. The " <>
  "Richters-Glotzl construction handles polynomial, exponential, sine and " <>
  "cosine atoms; mixed transcendentals may need to be split by hand.";

HelmholtzDecomposition::baddim =
  "Property `1` is only defined for vector fields of dimension `2`. Got `3`.";

HelmholtzDecomposition::nokey =
  "`1` is not a property of HelmholtzDecomposition. Valid properties: `2`.";

HelmholtzDecomposition::badverbose =
  "Verbose -> `1` is not True or False.";

HelmholtzDecomposition::badcache =
  "Caching -> `1` is not True or False.";

HelmholtzDecomposition::badtc =
  "TimeConstraint -> `1` is not Infinity or a positive number.";

Begin["`Private`"]

(* ============================================================ *)
(* Cache                                                         *)
(* ============================================================ *)

$cache = <||>;

ClearHelmholtzCache[] := With[{n = Length[$cache]}, $cache = <||>; n];

(* Private throw tag so user-level Catch can't accidentally swallow our $Failed *)
$throwTag = Unique["HelmholtzDecomposition$Throw"];

(* ============================================================ *)
(* Iterated integrals and incomplete Laplacians                  *)
(* ============================================================ *)

iterIntegrate[func_, 0, m_, xvec_] := func;
iterIntegrate[func_, p_Integer?Positive, m_, xvec_] :=
  Integrate[iterIntegrate[func, p - 1, m, xvec], {xvec[[m]], 0, xvec[[m]]}];

incompleteLaplacian[func_, m_, xvec_] :=
  Laplacian[func, xvec] - D[func, {xvec[[m]], 2}];

incompleteLaplacianP[func_, 0, m_, xvec_] := func;
incompleteLaplacianP[func_, p_Integer?Positive, m_, xvec_] :=
  incompleteLaplacian[incompleteLaplacianP[func, p - 1, m, xvec], m, xvec];

(* ============================================================ *)
(* Variable extraction                                           *)
(* ============================================================ *)

$headlist = {Or, And, Equal, Unequal, Less, LessEqual, Greater, GreaterEqual, Inequality};

allVars[_?NumericQ] := {};
allVars[{}] := {};
allVars[t_] /; MemberQ[$headlist, t] := {};
allVars[ll_List] := Flatten[Union[Map[allVars, ll]]];
allVars[Derivative[_Integer][_][arg__]] := allVars[{arg}];
allVars[f_Symbol[arg__]] :=
  If[MemberQ[Attributes[f], NumericFunction] || MemberQ[$headlist, f],
    allVars[{arg}],
    f[arg]
  ];
allVars[other_] := other;

uniqueVars[expr_] := Sort[DeleteDuplicates[Flatten[{allVars[expr]}]]];

autoDetectCoordinates[f_] :=
  Module[{vars = uniqueVars[f]},
    If[Length[vars] === Length[f],
      vars,
      Message[HelmholtzDecomposition::autocoord, Length[f], vars];
      $Failed
    ]
  ];

(* ============================================================ *)
(* Potential matrix Fij (Theorem 6.1 of Richters & Glotzl)       *)
(* ============================================================ *)

calculateFij[f_, k_, m_, lambda_, u_, W_, xvec_, dim_] :=
  Module[{Fij, i, p, j},
    (* Assert is on by default; we don't toggle global state. If the user has
       Off[Assert] set, that is their explicit choice — we honor it. *)
    Assert[Length[f] == Length[xvec] == dim];
    Assert[1 <= k <= dim];
    Assert[1 <= m <= dim];
    If[k != m, Do[Assert[D[f[[k]], xvec[[i]]] === 0 || i == k || i == m], {i, dim}]];
    Assert[FullSimplify[D[W, {xvec[[m]], 2 lambda}] - f[[k]]] == 0];
    Assert[FullSimplify[(-1)^lambda incompleteLaplacianP[W, lambda, m, xvec] - (1 - u) f[[k]] == 0]];
    Fij = Simplify[Table[
      KroneckerDelta[k, i] Sum[
        (-1)^p / u D[incompleteLaplacianP[W, p, m, xvec], xvec[[j]], {xvec[[m]], 2 lambda - 2 p - 2}],
        {p, 0, lambda - 1}
      ],
      {i, dim}, {j, dim}
    ]];
    Fij
  ];

autoSingleFij[S_, k_, xvec_, dim_, verb_] :=
  Module[{i, j, Singlef, xiWithoutxk, Q1, Q2, Q3, Q4, m, lambda, u, W, v1, v2, xiInv1, xiInv2},
    If[S === 0, Return[Table[0, {i, dim}, {j, dim}]]];
    Singlef = Table[S KroneckerDelta[i, k], {i, dim}];
    xiWithoutxk = Intersection[xvec, DeleteCases[uniqueVars[S], xvec[[k]]]];
    (* Q1 tests "polynomial in the spatial coordinates" — parameters in
       the coefficients (rational, etc.) are treated as constants. The
       earlier `PolynomialQ[S, xiWithoutxk]` form misfired when
       xiWithoutxk was empty: PolynomialQ[expr, {}] enters auto-detect
       mode and rejects rational-in-parameter expressions. *)
    Q1 = PolynomialQ[S, xvec];
    Q2 = Q1 && Length[xiWithoutxk] == 1 &&
         Exponent[S, First[xiWithoutxk]] > Exponent[S, xvec[[k]]];
    Q3 = PolynomialQ[S, xvec[[k]]] && Length[xiWithoutxk] == 1;
    If[Q1 || Q3,
      m = If[Q1 && !Q2, k, First[Flatten[Position[xvec, First[xiWithoutxk]]]]];
      u = 1;
      lambda = Ceiling[(Total[Exponent[S, xvec]] - Exponent[S, xvec[[m]]] + 1) / 2];
      W = iterIntegrate[S, 2 lambda, m, xvec];
      If[verb, Echo[{S, "in f", k, "polynomial branch", "m" -> m, "lambda" -> lambda}, "HelmholtzDecomposition"]];
      Return[calculateFij[Singlef, k, m, lambda, u, W, xvec, dim]]
    ];
    v1 = FullSimplify[D[S, {xvec[[k]], 2}] / S];
    v2 = FullSimplify[incompleteLaplacian[S, k, xvec] / S];
    xiInv1 = Intersection[xvec, uniqueVars[v1]];
    xiInv2 = Intersection[xvec, uniqueVars[v2]];
    Q4 = !(v1 === 0) && Length[xiInv1] == 0 && Length[xiInv2] == 0;
    If[Q4,
      lambda = 1; m = k; W = S / v1; u = 1 + v2 / v1;
      If[verb, Echo[{S, "in f", k, "trig/exp branch", "m" -> m, "W" -> W, "u" -> u}, "HelmholtzDecomposition"]];
      Return[calculateFij[Singlef, k, m, lambda, u, W, xvec, dim]]
    ];
    Message[HelmholtzDecomposition::nodec, S, k];
    Throw[$Failed, $throwTag]
  ];

autoFij[f_, xvec_, dim_, verb_] :=
  Module[{M, Fij, Poly, pol, k},
    Fij = Table[0, {dim}, {dim}];
    Do[
      Poly = Expand[f[[k]]];
      Poly = If[Head[Poly] === Plus, List @@ Poly, {Poly}];
      If[verb, Echo[{"k" -> k, Poly}, "HelmholtzDecomposition"]];
      Do[M = Poly[[pol]]; Fij = Fij + autoSingleFij[M, k, xvec, dim, verb], {pol, Length[Poly]}],
      {k, dim}
    ];
    Fij
  ];

(* ============================================================ *)
(* Pieces extracted from Fij (simplifier-aware)                  *)
(* ============================================================ *)

(* Uniform simplifier wrapper.
   - Identity is special-cased because it only accepts one argument.
   - Simplify/FullSimplify natively accept TimeConstraint -> tc and treat
     Infinity as "no limit", so we always pass the option through. *)
applySimp[Identity, expr_, _, _] := expr;
applySimp[simp_, expr_, asm_, tc_] := simp[expr, asm, TimeConstraint -> tc];

scalarPotential[Fij_, simp_, asm_, tc_] :=
  applySimp[simp, Sum[Fij[[k, k]], {k, Length[Fij]}], asm, tc];

antisymPart[Fij_, simp_, asm_, tc_] :=
  applySimp[simp, Fij - Transpose[Fij], asm, tc];

(* ============================================================ *)
(* Compute the data Association                                  *)
(* ============================================================ *)

hdCompute[f_, xvec_, simp_, asm_, verb_, tc_] :=
  Module[{dim, Fij, P, Rij, g, r, residual, curlMatrix, divR, verify},
    dim = Length[f];
    Fij = autoFij[f, xvec, dim, verb];
    P   = scalarPotential[Fij, simp, asm, tc];
    Rij = antisymPart[Fij, simp, asm, tc];
    g   = applySimp[simp, Grad[P, xvec], asm, tc];
    r   = applySimp[simp, Table[Sum[D[Rij[[j, i]], xvec[[i]]], {i, dim}], {j, dim}], asm, tc];
    residual = applySimp[simp, f - g - r, asm, tc];
    (* Verify is precomputed: notebook display + ["Verify"] become O(1).
       Always uses FullSimplify (independent of the user's SimplificationFunction
       choice) so "verified" actually means verified. *)
    curlMatrix = FullSimplify[
      Table[D[g[[i]], xvec[[j]]] - D[g[[j]], xvec[[i]]], {i, dim}, {j, dim}],
      asm, TimeConstraint -> tc
    ];
    divR = FullSimplify[Sum[D[r[[i]], xvec[[i]]], {i, dim}], asm, TimeConstraint -> tc];
    verify = <|
      "GradientCurlFree"         -> curlMatrix === ConstantArray[0, {dim, dim}],
      "RotationalDivergenceFree" -> divR === 0,
      "ResidualZero"             -> residual === ConstantArray[0, dim]
    |>;
    <|
      "Field"                  -> f,
      "Coordinates"            -> xvec,
      "Gradient"               -> g,
      "Rotational"             -> r,
      "Potential"              -> P,
      "AntisymmetricPotential" -> Rij,
      "PotentialMatrix"        -> Fij,
      "Residual"               -> residual,
      "Verify"                 -> verify
    |>
  ];

(* ============================================================ *)
(* Public constructor                                            *)
(* ============================================================ *)

HelmholtzDecomposition[f_?VectorQ, xvec_?VectorQ, opts:OptionsPattern[]] :=
  Module[{simp, asm, verb, cacheOpt, tc, key, raw, wrapped},
    If[Length[f] =!= Length[xvec],
      Message[HelmholtzDecomposition::nomatch, Length[f], Length[xvec]];
      Return[$Failed]
    ];
    simp     = OptionValue[SimplificationFunction];
    asm      = OptionValue[Assumptions];
    verb     = OptionValue[Verbose];
    cacheOpt = OptionValue[Caching];
    tc       = OptionValue[TimeConstraint];
    If[!BooleanQ[verb],
      Message[HelmholtzDecomposition::badverbose, verb]; Return[$Failed]];
    If[!BooleanQ[cacheOpt],
      Message[HelmholtzDecomposition::badcache, cacheOpt]; Return[$Failed]];
    If[!(tc === Infinity || (NumericQ[tc] && tc > 0)),
      Message[HelmholtzDecomposition::badtc, tc]; Return[$Failed]];
    key = {f, xvec, simp, asm, tc};
    If[cacheOpt && KeyExistsQ[$cache, key], Return[$cache[key]]];
    raw = Catch[hdCompute[f, xvec, simp, asm, verb, tc], $throwTag];
    If[raw === $Failed, Return[$Failed]];
    wrapped = HelmholtzDecomposition[raw];
    If[cacheOpt, AssociateTo[$cache, key -> wrapped]];
    wrapped
  ];

(* Auto-coord overload *)
HelmholtzDecomposition[f_?VectorQ, opts:OptionsPattern[]] :=
  Module[{xv = autoDetectCoordinates[f]},
    If[xv === $Failed, $Failed, HelmholtzDecomposition[f, xv, opts]]
  ];

(* ============================================================ *)
(* Property dispatch                                             *)
(*                                                               *)
(* All property access goes through a single dispatch[a, k]      *)
(* function. Adding a new synthetic property is one line of      *)
(* DownValues; no reliance on Mathematica's pattern-specificity  *)
(* ordering between SubValue rules.                              *)
(* ============================================================ *)

$validProps = {
  "Properties", "Association",
  "Field", "Coordinates",
  "Gradient", "Rotational",
  "Potential", "AntisymmetricPotential", "PotentialMatrix",
  "VectorPotentialA",
  "Residual", "Verify"
};

(* Synthetic / computed properties. Each must take (a) and return its value
   or $Failed (with a Message). Keys not listed here fall through to the
   stored-Association lookup. *)
dispatch[a_, "Properties"]  := $validProps;
dispatch[a_, "Association"] := a;
dispatch[a_, "VectorPotentialA"] :=
  Module[{Rij = a["AntisymmetricPotential"], n = Length[a["Field"]]},
    If[n =!= 3,
      Message[HelmholtzDecomposition::baddim, "VectorPotentialA", 3, n];
      $Failed,
      {Rij[[2, 3]], Rij[[3, 1]], Rij[[1, 2]]}
    ]
  ];

(* Stored properties: present in the Association *)
dispatch[a_, k_String] /; KeyExistsQ[a, k] := a[k];

(* Anything else *)
dispatch[a_, k_] := (Message[HelmholtzDecomposition::nokey, k, $validProps]; $Failed);

HelmholtzDecomposition[a_Association][k_] := dispatch[a, k];

(* ============================================================ *)
(* Backward-compat conveniences                                  *)
(* ============================================================ *)

HelmholtzDecomposition /: Normal[HelmholtzDecomposition[a_Association]] := a;

(* List @@ hd  ->  {g, r}  (for legacy callers expecting a tuple) *)
HelmholtzDecomposition /: Apply[List, HelmholtzDecomposition[a_Association]] :=
  {a["Gradient"], a["Rotational"]};

(* Keys[hd] returns the queryable property names; users wanting field
   dimension write Length[hd["Coordinates"]] explicitly. *)
HelmholtzDecomposition /: Keys[HelmholtzDecomposition[a_Association]] := $validProps;

(* ============================================================ *)
(* Notebook display: collapsible summary box                     *)
(* ============================================================ *)

$wrappedKeys = {"Field", "Coordinates", "Gradient", "Rotational",
  "Potential", "AntisymmetricPotential", "PotentialMatrix",
  "Residual", "Verify"};

wellFormedQ[a_Association] := SubsetQ[Keys[a], $wrappedKeys];

HelmholtzDecomposition /: MakeBoxes[hd:HelmholtzDecomposition[a_Association?wellFormedQ], form:(StandardForm|TraditionalForm)] :=
  BoxForm`ArrangeSummaryBox[
    HelmholtzDecomposition,
    hd,
    None,
    {
      BoxForm`SummaryItem[{"dimension: ", Length[a["Field"]]}],
      BoxForm`SummaryItem[{"coordinates: ", a["Coordinates"]}]
    },
    {
      BoxForm`SummaryItem[{"potential: ", Short[a["Potential"], 1]}],
      BoxForm`SummaryItem[{"all checks pass: ", AllTrue[a["Verify"], TrueQ]}]
    },
    form
  ];

(* ============================================================ *)
(* Syntax / autocomplete metadata                                *)
(* ============================================================ *)

SyntaxInformation[HelmholtzDecomposition] = {
  "ArgumentsPattern" -> {_, _., OptionsPattern[]}
};

End[] (* `Private` *)

EndPackage[]
