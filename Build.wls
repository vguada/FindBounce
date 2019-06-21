#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsubsection:: *)
(*Git functions*)


(* Get human readable description of current state of repository.
For more info see: https://git-scm.com/docs/git-describe *)
getGitRevision[dir_]:=Module[
	{description},
	SetDirectory[dir];
	description=First[
		ReadList["!git describe --long --tags --dirty --always",String],
		Return[$Failed]
	];
	ResetDirectory[];
	description
];


(* Get repository commit count. It is useful to mark "build number". *)
getGitCommitCount[dir_]:=Module[
	{value},
	SetDirectory[dir];
	value=First[
		ReadList["!git rev-list --count HEAD",Number],
		Return[$Failed]
	];
	ResetDirectory[];
	value
];


(* ::Subsubsection:: *)
(*Build procedure*)


$name="FindBounce";
(* Get the distribution directory by using this notebook\[CloseCurlyQuote]s directory. *)
$rootDir=If[$Notebooks,NotebookDirectory[],Directory[]];
$source=FileNameJoin[{$rootDir,$name}];
(* Construct the target directory from this notebooks base directory. *)
$target=FileNameJoin[{$rootDir,"build",$name}];


$helpMessage=("  Usage: Build.wls [options]

  Options:
    -h, --help                   Output usage information
    -r, --release                Build public release");


Module[
	{directories},
	If[
		MemberQ[Rest@$ScriptCommandLine,"/?"|"-h"|"--help"],
		Print[$helpMessage];Quit[1]
	];
	
	If[
		Not@DirectoryQ[$target],
		CreateDirectory[$target,CreateIntermediateDirectories->True]
	];
	(* Copy package files from the source to the target directory. *)
	Map[
		CopyFile[
			FileNameJoin[{$source,#}],
			FileNameJoin[{$target,#}],
			OverwriteTarget->True
		]&,
		{"FindBounce.m","PacletInfo.m"}
	];
	
	(* Delete old and copy the new source directories. 
	Documentation should be built with a separate procedure. *)
	directories={"FrontEnd","Kernel"};
	Map[
		If[
			FileExistsQ[FileNameJoin[{$target,#}]],
			DeleteDirectory[FileNameJoin[{$target,#}],DeleteContents->True]
		]&,
		directories
	];
	Map[
		CopyDirectory[FileNameJoin[{$source,#}],FileNameJoin[{$target,#}]]&,
		directories
	];
];


Module[
	{original,modified,noCommits,revision,publicQ,version},

	original=List@@Import[FileNameJoin[{$source,"PacletInfo.m"}]];
	noCommits=getGitCommitCount[$rootDir];
	revision=getGitRevision[$rootDir];
	publicQ=MemberQ[Rest@$ScriptCommandLine,"-r"|"--release"];
	version=Version/.original;
	(* Internal releases have number of commits appended. Public releases conform to 
	"semantic versioning" convention with 3 version numbers. *)
	If[Not@publicQ,version=version<>"."<>ToString@noCommits];
		
	(* Description of git repository (git describe ...) returns a string that cannot 
	be used as BuildNumber because only integers are accepted. 
	Otherwise PackPaclet returns $Failed and procedure fails. Instead commit count is used. *)
	modified=Paclet@@Normal[Merge[{
		Association@@original,
		<|Version->version,Internal->Not@publicQ,BuildNumber->noCommits,"Hash"->revision|>
		},
		Last
	]];
	
	If[
		StringContainsQ[revision,"dirty"],
		Print["  Warning! Repository contains uncommited changes."]
	];
	Export[
		FileNameJoin[{$target,"PacletInfo.m"}],
		modified,
		"Comments"->{"Paclet Info File","Created on "<>DateString[]}
	];
];


Module[
	{pacletPath,newPaclet},
	
	Needs["PacletManager`"];
	(* Change into the build directory and create the paclet. *)	
	SetDirectory[FileNameJoin[{$rootDir,"build"}]];
	pacletPath=PackPaclet[$name];
	ResetDirectory[];
	
	(* Uninstall older versions first. *)
	If[PacletFind[$name]=!={},PacletUninstall[$name]];
	newPaclet=PacletInstall[pacletPath];
	RebuildPacletData[];
	
	If[
		Length@PacletFind["HeatTrans"]===1,
		Print["  ",FileNameTake[pacletPath]," built and installed succesfully."];Quit[0],
		Print["  Paclet build failed!"];Quit[1]
	]
];
