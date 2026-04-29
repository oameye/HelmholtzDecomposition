#!/usr/bin/env wolframscript
(* CLI-runnable test driver. Usage:
     wolframscript -code '$CharacterEncoding="UTF-8"; Get["/abs/path/TestFiles/run_tests.wls"]'
   Exits with non-zero status if any test fails.
*)

$CharacterEncoding = "UTF-8";
Needs["PacletTools`"];

$pacletRoot = ParentDirectory[DirectoryName[$InputFileName]];
PacletDirectoryLoad[$pacletRoot];

$testDir   = DirectoryName[$InputFileName];
$testFiles = File /@ FileNames["*.wlt", $testDir];

rep = TestReport[$testFiles];

Print["Tests run:   ", rep["TestsSucceededCount"] + rep["TestsFailedCount"]];
Print["Succeeded:   ", rep["TestsSucceededCount"]];
Print["Failed:      ", rep["TestsFailedCount"]];

If[rep["TestsFailedCount"] > 0,
  Print["\n--- Failures ---"];
  Scan[
    Function[tro,
      If[tro["Outcome"] =!= "Success",
        Print["TestID:   ", tro["TestID"]];
        Print["Outcome:  ", tro["Outcome"]];
        Print["Input:    ", tro["Input"]];
        Print["Expected: ", tro["ExpectedOutput"]];
        Print["Actual:   ", tro["ActualOutput"]];
        Print["Messages: ", tro["ActualMessages"]];
        Print[""]
      ]
    ],
    Values @ rep["TestResults"]
  ];
  Exit[1]
];

Exit[0]
