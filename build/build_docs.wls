#!/usr/bin/env wolframscript
(* Generates the paclet documentation notebooks from plain-text specs.
   Run with:
     wolframscript -code '$CharacterEncoding="UTF-8"; Get["/abs/path/build/build_docs.wls"]'
*)

$CharacterEncoding = "UTF-8";

$pacletRoot = FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]]}];
$docRoot = FileNameJoin[{$pacletRoot, "Documentation", "English"}];

(* ===== shared cell helpers ===== *)

styleDefs = FrontEnd`FileName[{"Wolfram"}, "Reference.nb", CharacterEncoding -> "UTF-8"];

usageRow[lhs_String, rhs_String] :=
  Cell[BoxData[GridBox[{{
      Cell[lhs, "TableText"],
      Cell[rhs, "TableText"]
    }}, ColumnsEqual -> False]], "Usage"];

textCell[s_String] := Cell[s, "Text"];
notesCell[s_String] := Cell[s, "Notes"];
inputCell[code_String] := Cell[BoxData[code], "Input"];

primaryExamplesHeader[] := Cell["Basic Examples", "PrimaryExamplesSection"];
moreExamplesHeader[]    := Cell["Scope", "ExtendedExamplesSection"];
optionsHeader[]         := Cell["Options", "NotesSection"];
notesHeader[]           := Cell["Details", "NotesSection"];
sourceHeader[]          := Cell["Source", "NotesSection"];

(* ===== HelmholtzDecomposition reference page ===== *)

mainPage = Notebook[{
  Cell["HelmholtzDecomposition", "ObjectName"],

  usageRow[
    "HelmholtzDecomposition[f, xvec]",
    "returns a HelmholtzDecomposition object whose properties give the analytical Helmholtz decomposition f = g + r of the vector field f with respect to the coordinates xvec."],
  usageRow[
    "HelmholtzDecomposition[f]",
    "auto-detects the coordinates from the symbols appearing in f."],
  usageRow[
    "hd[\"prop\"]",
    "returns the named property of a HelmholtzDecomposition object hd."],

  notesHeader[],

  notesCell["The decomposition is constructed via the algorithm of Richters and Glotzl (J. Math. Anal. Appl. 525, 127138, 2023). The algorithmic core handles polynomial, exponential, sine and cosine atomic terms; mixed transcendentals may need to be split by hand."],
  notesCell["The result is a wrapped Association with head HelmholtzDecomposition. It carries enough state to recover all derived quantities without recomputation."],
  notesCell["Available properties (queried via hd[\"name\"]):"],
  notesCell["  \"Gradient\"               g, the curl-free component"],
  notesCell["  \"Rotational\"             r, the divergence-free component"],
  notesCell["  \"Potential\"              scalar P with Grad[P] == g"],
  notesCell["  \"AntisymmetricPotential\" the antisymmetric matrix Rij whose row-divergence is r"],
  notesCell["  \"VectorPotentialA\"       (3D only) vector A with Curl[A] == r"],
  notesCell["  \"PotentialMatrix\"        the underlying potential matrix Fij"],
  notesCell["  \"Residual\"               f - g - r, simplified"],
  notesCell["  \"Field\", \"Coordinates\"   the inputs"],
  notesCell["  \"Verify\"                 <| GradientCurlFree, RotationalDivergenceFree, ResidualZero |>"],
  notesCell["  \"Properties\"             list of the keys above"],
  notesCell["  \"Association\"            the raw underlying Association"],
  notesCell["List @@ hd returns {g, r}; Normal[hd] returns the underlying Association."],

  optionsHeader[],

  notesCell["SimplificationFunction -> FullSimplify"],
  notesCell["    Function applied to the components. Use Simplify for speed, Identity to skip simplification entirely."],
  notesCell["Assumptions :> $Assumptions"],
  notesCell["    Forwarded to the simplifier."],
  notesCell["Verbose -> False"],
  notesCell["    When True, traces which atoms hit which branch of the algorithm via Echo."],
  notesCell["Caching -> True"],
  notesCell["    When True (default), memoizes results on (f, xvec, simp, asm). ClearHelmholtzCache[] empties the cache."],

  primaryExamplesHeader[],

  textCell["Decompose a 2D polynomial field:"],
  inputCell["hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]"],
  inputCell["hd[\"Gradient\"]"],
  inputCell["hd[\"Rotational\"]"],
  inputCell["hd[\"Potential\"]"],
  inputCell["hd[\"Verify\"]"],

  textCell["Auto-detect coordinates and recover the legacy {g, r} list:"],
  inputCell["List @@ HelmholtzDecomposition[{x^2 + y, x y}]"],

  textCell["Pass assumptions to the simplifier:"],
  inputCell["HelmholtzDecomposition[{a x^2, b y}, {x, y},\n  Assumptions -> a > 0 && b > 0][\"Potential\"]"],

  textCell["3D vector potential A:"],
  inputCell["HelmholtzDecomposition[{-y - z, x, 0}, {x, y, z}][\"VectorPotentialA\"]"],

  sourceHeader[],
  notesCell["Erhard Glotzl, Oliver Richters, \"Helmholtz decomposition and potential functions for n-dimensional analytic vector fields\", J. Math. Anal. Appl. 525, 127138 (2023). https://doi.org/10.1016/j.jmaa.2023.127138"]
}, StyleDefinitions -> styleDefs];

mainPagePath = FileNameJoin[{$docRoot, "ReferencePages", "Symbols", "HelmholtzDecomposition.nb"}];
Print["Writing ", mainPagePath];
Put[mainPage, mainPagePath];

(* ===== ClearHelmholtzCache reference page ===== *)

clearCachePage = Notebook[{
  Cell["ClearHelmholtzCache", "ObjectName"],

  usageRow[
    "ClearHelmholtzCache[]",
    "empties the memoization cache used by HelmholtzDecomposition. Returns the number of entries that were removed."],

  notesHeader[],

  notesCell["HelmholtzDecomposition memoizes results on the tuple (f, xvec, SimplificationFunction, Assumptions). For long-running notebooks computing many distinct fields, this can grow without bound. ClearHelmholtzCache[] is the manual escape hatch."],
  notesCell["Pass Caching -> False to a particular HelmholtzDecomposition call to skip the cache for that call only."],

  primaryExamplesHeader[],

  inputCell["HelmholtzDecomposition[{x^2, y^2}, {x, y}];"],
  inputCell["ClearHelmholtzCache[]"]
}, StyleDefinitions -> styleDefs];

clearCachePath = FileNameJoin[{$docRoot, "ReferencePages", "Symbols", "ClearHelmholtzCache.nb"}];
Print["Writing ", clearCachePath];
Put[clearCachePage, clearCachePath];

(* ===== Top-level guide page ===== *)

guidePage = Notebook[{
  Cell["HelmholtzDecomposition", "GuideTitle"],
  Cell["Analytical Helmholtz decomposition of n-dimensional vector fields, following Richters and Glotzl (2023).", "GuideAbstract"],

  Cell["Decomposing a vector field", "GuideFunctionsSection"],
  Cell["HelmholtzDecomposition", "InlineFunctionSans"],
  Cell["    construct the decomposition object f = g + r and access derived quantities via properties", "GuideText"],

  Cell["Cache management", "GuideFunctionsSection"],
  Cell["ClearHelmholtzCache", "InlineFunctionSans"],
  Cell["    drop all memoized decomposition results", "GuideText"],

  Cell["Tutorials", "GuideTutorialsSection"],
  Cell[TextData[ButtonBox["Helmholtz decomposition tutorial",
    BaseStyle -> "Link",
    ButtonData -> "paclet:HelmholtzDecomposition/tutorial/HelmholtzDecompositionTutorial"]], "GuideTutorial"]
}, StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "Guide.nb", CharacterEncoding -> "UTF-8"]];

guidePath = FileNameJoin[{$docRoot, "Guides", "HelmholtzDecomposition.nb"}];
Print["Writing ", guidePath];
Put[guidePage, guidePath];

(* ===== Tutorial ===== *)

tutorialPage = Notebook[{
  Cell["HelmholtzDecomposition tutorial", "Title"],

  Cell["Introduction", "Section"],
  textCell["Every smooth vector field f on a star-shaped domain admits a decomposition"],
  textCell["    f = g + r,    Curl[g] == 0,    Div[r] == 0,"],
  textCell["with g a gradient and r divergence-free. The Richters-Glotzl construction makes this analytic for fields built from polynomial, exponential, sine and cosine atoms in any number of dimensions."],

  Cell["Quick start", "Section"],
  inputCell["Needs[\"HelmholtzDecomposition`\"]"],
  inputCell["hd = HelmholtzDecomposition[{x^2 + y, x y}, {x, y}]"],
  textCell["The result is a wrapped object. Query it via properties:"],
  inputCell["hd[\"Gradient\"]"],
  inputCell["hd[\"Rotational\"]"],
  inputCell["hd[\"Potential\"]"],

  Cell["Verifying a decomposition", "Section"],
  textCell["The Verify property reports the three sanity checks:"],
  inputCell["hd[\"Verify\"]"],

  Cell["Working with parameters", "Section"],
  textCell["Pass Assumptions to control simplification:"],
  inputCell["HelmholtzDecomposition[{a x^2, b y}, {x, y},\n  Assumptions -> a > 0 && b > 0][\"Potential\"]"],

  Cell["3D vector potential", "Section"],
  textCell["In three dimensions you can recover the conventional vector potential A with Curl[A] == r:"],
  inputCell["hd3 = HelmholtzDecomposition[{-y - z, x, 0}, {x, y, z}];"],
  inputCell["A = hd3[\"VectorPotentialA\"]"],
  inputCell["Curl[A, {x, y, z}] - hd3[\"Rotational\"]"],

  Cell["Cache management", "Section"],
  textCell["Repeated calls on the same field reuse cached results:"],
  inputCell["First @ AbsoluteTiming @ HelmholtzDecomposition[{x^3 y, x y^2}, {x, y}]"],
  inputCell["First @ AbsoluteTiming @ HelmholtzDecomposition[{x^3 y, x y^2}, {x, y}]"],
  textCell["Drop the cache when memory matters:"],
  inputCell["ClearHelmholtzCache[]"]
}, StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "Tutorial.nb", CharacterEncoding -> "UTF-8"]];

tutorialPath = FileNameJoin[{$docRoot, "Tutorials", "HelmholtzDecompositionTutorial.nb"}];
Print["Writing ", tutorialPath];
Put[tutorialPage, tutorialPath];

Print["\nDocumentation generated under ", $docRoot];
