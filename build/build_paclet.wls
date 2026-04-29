#!/usr/bin/env wolframscript
(* Build the paclet archive. Stages in /tmp to avoid recursive self-include
   (CreatePacletArchive walks the whole source tree, including build/).
   Run with:
     wolframscript -code '$CharacterEncoding="UTF-8"; Get["/abs/path/build/build_paclet.wls"]'
*)

$CharacterEncoding = "UTF-8";
Needs["PacletTools`"];

$pacletRoot = ParentDirectory[DirectoryName[$InputFileName]];
$stageDir   = FileNameJoin[{$TemporaryDirectory, "HelmholtzDecomposition_build"}];

If[DirectoryQ[$stageDir], DeleteDirectory[$stageDir, DeleteContents -> True]];
CreateDirectory[$stageDir];

archive = CreatePacletArchive[$pacletRoot, $stageDir];
Print["Built ", archive, " (", FileByteCount[archive], " bytes)"];

target = FileNameJoin[{$pacletRoot, "build", FileNameTake[archive]}];
CopyFile[archive, target, OverwriteTarget -> True];
Print["Copied to ", target];
