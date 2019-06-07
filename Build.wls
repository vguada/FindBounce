#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsubsection:: *)
(*Git functions*)


getGitRevision[dir_]:=Module[
	{description},
	SetDirectory[dir];
	(* For meaning of this see: https://git-scm.com/docs/git-describe *)
	description=First[
		ReadList["!git describe --long --tags --dirty --always",String],
		Return[$Failed]
	];
	ResetDirectory[];
	description
];


getGitCommitCount[dir_]:=Module[
	{value},
	SetDirectory[dir];
	(* For meaning of this see: https://git-scm.com/docs/git-describe *)
	value=First[
		ReadList["!git rev-list --count HEAD",Number],
		Return[$Failed]
	];
	ResetDirectory[];
	value
]


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
	(* Delete old and copy the new source directories. Documentation notebook stays unmodified. *)
	directories={"Documentation","FrontEnd","Kernel"};
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
	{original, modified,noCommits,revision,internalQ},

	original=Import[FileNameJoin[{$source,"PacletInfo.m"}]];
	noCommits=getGitCommitCount[$rootDir];
	revision=getGitRevision[$rootDir];
	internalQ=Not@MemberQ[Rest@$ScriptCommandLine,"-r"|"--release"];
		
	(* Description of git repository (git describe ...) returns a string that cannot 
	be used as BuildNumber because only integers are accepted. 
	Otherwise PackPaclet returns $Failed and procedure fails. Instead now commit count is used. *)
	modified=Fold[
			Insert[#1,#2,4]&,
			original,
			{Internal->internalQ,BuildNumber->noCommits,"Hash"->revision}
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
	
	(*Uninstall a possibly old version. *)
	If[PacletFind[$name]=!={},PacletUninstall[$name]];
	newPaclet=PacletInstall[pacletPath];
	RebuildPacletData[];
	
	If[
		Length@PacletFind[$name]===1,
		Print["  Paclet built and installed succesfully: ",FileNameTake[pacletPath]];Quit[0],
		Print["  Paclet build failed."];Quit[1]
	]
];
