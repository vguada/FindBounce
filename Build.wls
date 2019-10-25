#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Build procedure*)


(* ::Subsubsection:: *)
(*Initialization*)


$name="FindBounce";
(* Get the distribution directory by using this notebook directory. *)
$root=If[$Notebooks,NotebookDirectory[],Directory[]];
(* Directory of main source code. *)
$source=FileNameJoin[{$root,$name}];
(* Construct the target directory from this notebooks base directory. *)
$target=FileNameJoin[{$root,"build",$name}];


If[$VersionNumber<10.0, Print["Mathematica 10.0 or later required."];Quit[1]];


$helpMessage=("  Usage: wolframscript -script Build.wls [options]

  Options:
    -h, --help                   Output usage information
    --release                    Build public release
    --docs                       Build documentation");

(* This should support text based interface with "wolframscript" mechanism and 
"wolfram -script" mechanism. See tutorial/WolframLanguageScripts for more info. *)
$args=Join[$ScriptCommandLine,$CommandLine];
If[
	MemberQ[$args,"/?"|"-h"|"--help"],
	Print[$helpMessage];Quit[0]
];


(* Verify that we have git *)
With[
	{rg=RunProcess[{"git","--version"}]},
	If[
		rg===$Failed||rg["ExitCode"]!=0,
		Print["git software is not available."];Quit[1]
	]
];


If[
	Not@DirectoryQ[$target],
	CreateDirectory[$target,CreateIntermediateDirectories->True]
];


(* ::Subsubsection:: *)
(*Process documentation notebooks*)


(* Notebook name should be relative to root directory. *)
processDocumentation[nbName_]:=Module[
	{original,copy,copyDir,nb,tStart,time},
	tStart=AbsoluteTime[];
	original=FileNameJoin[{$source,nbName}];
	If[
		Not@FileExistsQ[original],
		Return[$Failed];Print[original," file does not exist!"]
	];
	copy=FileNameJoin[{$target,nbName}];
	copyDir=DirectoryName[copy];
	If[
		Not@DirectoryQ[copyDir],
		CreateDirectory[copyDir,CreateIntermediateDirectories->True]
	];
	CopyFile[original,copy,OverwriteTarget->True];

	UsingFrontEnd[
		NotebookEvaluate[copy,InsertResults->True];
		nb=NotebookOpen[copy];
		SetOptions[
			nb,
			WindowSize->{Scaled[0.5], Scaled[1]},
			WindowMargins->{{0,Automatic},{Automatic,0}},
			Saveable->False
		];
		NotebookSave[nb];
		NotebookClose[nb];
	];
	time=Round[AbsoluteTime[]-tStart,0.1];
	Print["  Processed ",nbName," in ",time," seconds."]
];


If[
	MemberQ[$args,"--docs"|"--release"],

	(* NotebookEvaluate doesn't work properly (unrendered graphic output)
	if $EvaluationEnvironment is not "Session". *)
	If[
		$EvaluationEnvironment=!="Session",
		Print["For correct processing of documentation run script with: \"wolframscript -script Build.wls [opts]\"."];
		Quit[1]
	];
	Print["  Building documentation. This may take a few minutes..."];
	(* Documentation should be evaluated with development version of the package. 
	Actuall loading of the package should happen inside each documentation notebook. *)
	Needs["PacletManager`"];
	PacletDirectoryAdd[$root];

	SetOptions[First[$Output],FormatType->StandardForm];

	(* Add more function calls here if documentation contains more notebooks. *)
	processDocumentation[
		FileNameJoin[{"Documentation","English","Tutorials","Examples_FindBounce.nb"}]
	];

	(* Remove path to development version to enble correct paclet uninstallation 
	before installing a newly produced .paclet file. *)
	PacletDirectoryRemove[$root];
	Print["\n  Documentation build finished."];
];


(* ::Subsubsection:: *)
(*Copy files*)


Module[
	{},
	(* Copy package files from the source to the target directory. *)
	Map[
		CopyFile[
			FileNameJoin[{$source,#}],
			FileNameJoin[{$target,#}],
			OverwriteTarget->True
		]&,
		{"FindBounce.m","PacletInfo.m"}
	];

	(* Delete old and copy the new source directories. *)
	Map[
		(
		If[
			FileExistsQ[FileNameJoin[{$target,#}]],
			DeleteDirectory[FileNameJoin[{$target,#}],DeleteContents->True]
		];
		CopyDirectory[FileNameJoin[{$source,#}],FileNameJoin[{$target,#}]]
		)&,
		{"FrontEnd","Kernel"}
	];

	(* Documentation should be built with a separate procedure before. *)
	If[
		Not@DirectoryQ[FileNameJoin[{$target,"Documentation"}]],
		Print["  Warning! Documentation directory is missing in build folder."]
	];
];


Module[
	{original,modified,noCommits,revision,publicQ,version},

	original=List@@Import[FileNameJoin[{$source,"PacletInfo.m"}]];
	noCommits=getGitCommitCount[$root];
	revision=getGitRevision[$root];
	publicQ=MemberQ[$args,"--release"];

	(* Description of git repository (git describe ...) returns a string that cannot 
	be used as BuildNumber because only integers are accepted. 
	Otherwise PackPaclet returns $Failed and procedure fails. Instead commit count is used. 
	Version cannot be appended with build number, because otherwise we break procedure
	of loading development version with PacletDirectoryAdd. *)
	modified=Paclet@@Normal[Merge[{
		Association@@original,
		<|Internal->Not@publicQ,BuildNumber->noCommits,"Hash"->revision|>
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


(* ::Subsubsection:: *)
(*Create .paclet file*)


Module[
	{pacletPath,newPaclet,buildNo},
	Needs["PacletManager`"];
	(* Change into the build directory and create the paclet. *)	
	SetDirectory@FileNameJoin[{$root,"build"}];
	pacletPath=PackPaclet[$name];
	ResetDirectory[];

	(* Uninstall older versions first. *)
	If[PacletFind[$name]=!={},PacletUninstall[$name]];
	newPaclet=PacletInstall[pacletPath];
	RebuildPacletData[];

	If[
		Length@PacletFind[$name]===1
		,
		buildNo=getGitCommitCount[$root];
		Print["  ",FileNameTake[pacletPath]," (",buildNo,") built and installed succesfully."];
		Quit[0]
		,
		Print["  Paclet build failed!"];
		Quit[1]
	]
];
