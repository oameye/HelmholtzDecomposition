Needs["PacletTools`"];
PacletDirectoryLoad[FileNameJoin[{NotebookDirectory[], ".."}]];

testDir = FileNameJoin[{NotebookDirectory[]}];  (*adjust as needed*)
testFiles = 
 File /@ FileNames["*.wlt" | "*.mt" | "*.nb", testDir, Infinity]

rep = TestReport[testFiles]