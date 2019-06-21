#!/usr/bin/env wolframscript
(* ::Package:: *)

Print["Building documentation. This may take a few minutes..."];


$name="FindBounce";
(* Get the distribution directory by using this notebook\[CloseCurlyQuote]s directory. *)
$rootDir=If[$Notebooks,NotebookDirectory[],Directory[]];
$source=FileNameJoin[{$rootDir,$name,"Documentation","English","Tutorials"}];
(* Construct the target directory from this notebooks base directory. *)
$target=FileNameJoin[{$rootDir,"build",$name,"Documentation","English","Tutorials"}];

(* Documentation should be evaluated with development version of the package. *)
Needs["PacletManager`"];
PacletDirectoryAdd[$rootDir];

If[Not@DirectoryQ[$target],CreateDirectory[$target,CreateIntermediateDirectories->True];]
SetOptions[First[$Output],FormatType->StandardForm];


processDocumentation[nbName_]:=Module[
	{original,copy,nb,tStart,time},
	tStart=AbsoluteTime[];
	original=FileNameJoin[{$source,nbName}];
	If[
		Not@FileExistsQ[original],
		Return[$Failed];Print[original," file does not exist!"]
	];
	copy=FileNameJoin[{$target,nbName}];
	If[FileExistsQ@copy,DeleteFile@copy];
	CopyFile[original,copy];
	
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


processDocumentation["Examples_FindBounce.nb"];


Print["\nDocumentation build finished."];
Quit[];
