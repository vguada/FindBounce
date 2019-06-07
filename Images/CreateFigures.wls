#!/usr/bin/env wolframscript
(* ::Package:: *)

(* $rootDir should be repository root directory. *)
$rootDir=If[$Notebooks,ParentDirectory@NotebookDirectory[],Directory[]];


(* ::Subsubsection:: *)
(*Create package icon*)


(* Icon should have 46x46 size. *)
icon=Plot[
	{LogisticSigmoid[x],LogisticSigmoid[-x-2]+0.2},
	{x,-10,10},
	Axes->False,
	AspectRatio->1,
	PlotStyle->Thickness[0.05],
	ImageSize->46,
	Frame->True,
	FrameTicks->None,
	FrameStyle->Gray
]


Export[
	FileNameJoin[{$rootDir,"FindBounce","FrontEnd","Icon.png"}],
	icon
]
