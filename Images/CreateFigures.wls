#!/usr/bin/env wolframscript
(* ::Package:: *)

(* $rootDir should be repository root directory. *)
$rootDir=If[$Notebooks,ParentDirectory@NotebookDirectory[],Directory[]];


(* ::Subsection:: *)
(*Package icon*)


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


(* ::Subsection:: *)
(*README examples*)


(* We want to load the development version. *)
PacletDirectoryAdd[$rootDir];
Needs["FindBounce`"]


(* ::Subsubsection:: *)
(*One field usage example*)


potential[x_]:=0.5 x^2-0.5 x^3+0.1 x^4;

extrema=Block[{x},x/.NSolve[D[potential[x],x]==0,x]]


pts=Transpose[{extrema,potential/@extrema}];
Plot[
	potential[x],
	{x,-1,4},
	AxesLabel->{"\[Phi]","V(\[Phi])"},
	LabelStyle->12,
	Epilog->{Red,PointSize[Large],Point[pts]}
]
(* We use this naming for .png file, because there might be more than on usage example. *)
Export[
	FileNameJoin[{$rootDir,"Images","UsageExample_1-1.png"}],
	%,
	ImageSize->{600,200}
];


bf=FindBounce[potential[x],{x},{extrema[[1]],extrema[[3]]}]


bf["Action"]


bf["Dimension"]


bf["Properties"]


BouncePlot[bf]

Export[
	FileNameJoin[{$rootDir,"Images","UsageExample_1-2.png"}],
	%,
	ImageSize->{600,200}
];
