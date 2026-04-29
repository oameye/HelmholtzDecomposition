BeginTestSection["Documentation"];

(* Parse-check the generated paclet doc notebooks. We can't render them
   without a frontend, but we can confirm they're syntactically valid
   Notebook[...] expressions with the expected top-level cell structure. *)

$docRoot = FileNameJoin[{DirectoryName[$TestFileName, 2],
                         "Documentation", "English"}];

$docFiles = {
  FileNameJoin[{$docRoot, "ReferencePages", "Symbols", "HelmholtzDecomposition.nb"}],
  FileNameJoin[{$docRoot, "ReferencePages", "Symbols", "ClearHelmholtzCache.nb"}],
  FileNameJoin[{$docRoot, "Guides", "HelmholtzDecomposition.nb"}],
  FileNameJoin[{$docRoot, "Tutorials", "HelmholtzDecompositionTutorial.nb"}]
};

VerificationTest[
  AllTrue[$docFiles, FileExistsQ],
  True,
  TestID -> "All four doc notebooks present"
];

VerificationTest[
  AllTrue[$docFiles, MatchQ[Get[#], _Notebook] &],
  True,
  TestID -> "Each doc notebook parses as a Notebook expression"
];

VerificationTest[
  AllTrue[$docFiles,
    Function[file,
      With[{nb = Get[file]},
        MatchQ[nb, Notebook[{__Cell}, ___]]
      ]
    ]
  ],
  True,
  TestID -> "Each doc notebook contains at least one Cell"
];

VerificationTest[
  Module[{nb = Get[FileNameJoin[{$docRoot, "ReferencePages", "Symbols",
    "HelmholtzDecomposition.nb"}]]},
    MemberQ[Cases[nb, Cell[_, style_, ___] :> style, Infinity], "ObjectName"]
  ],
  True,
  TestID -> "Reference page has an ObjectName cell"
];

EndTestSection[]
