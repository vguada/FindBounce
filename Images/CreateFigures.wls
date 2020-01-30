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
	LabelStyle->14,
	Epilog->{Red,PointSize[Large],Point[pts]},
	ImageSize->400,
	ImageMargins->{100*{1,1},{0,0}}
]
(* We use this naming for .png file, because there might be more than on usage example. 
Print shows the path of exported file if we used this script from command line interface. *)
Print@Export[
	FileNameJoin[{$rootDir,"Images","UsageExample_1-1.png"}],
	%
]


bf=FindBounce[potential[x],{x},{extrema[[1]],extrema[[3]]}]


bf["Action"]


bf["Dimension"]


bf["Properties"]


BouncePlot[bf,
	LabelStyle->14,
	ImageSize->400,
	ImageMargins->{100*{1,1},{0,0}}
]

Print@Export[
	FileNameJoin[{$rootDir,"Images","UsageExample_1-2.png"}],
	%
]


(* ::Subsubsection:: *)
(*Two fields example*)


Clear[V]
V[h_,s_]:= -100 h^2 + 0.1 h^4 - 60 s^2 + 0.3 s^4 + 3 h^2 s^2;
minima = {{0.,10.},{22.4,0.}};


bf = FindBounce[V[h,s],{h,s}, minima, "MidFieldPoint"-> {6.,6.}]


(* Some less important options can be deleted in README.md description. *)
{ri,rf} = MinMax@bf["Radii"];
Show[
	ContourPlot[
		V[h,s],{h,-1,25},{s,-1,11},
		Contours->40,
		ContourShading->None,
		ContourStyle->Gray,
		LabelStyle->14,
		Epilog->{PointSize[Large],Red,Point[minima]},
		FrameLabel->{"h","s"},
		ImageSize->400,
		ImageMargins->{100*{1,1},{0,0}}
	],
	ParametricPlot[Through@bf["Bounce"][r],{r,ri,rf},PlotStyle->Thick]
]

Print@Export[
	FileNameJoin[{$rootDir,"Images","UsageExample_2-1.png"}],
	%
]


BouncePlot[bf,
	PlotLabel-> Row[{"Action = ",bf["Action"]}],
	PlotStyle-> {Purple, Orange},
	LabelStyle->14,
	FrameLabel->{"\[Rho]","h,s"},
	ImageSize->400,
	ImageMargins->{100*{1,1},{0,0}}
]

Print@Export[
	FileNameJoin[{$rootDir,"Images","UsageExample_2-2.png"}],
	%
]
