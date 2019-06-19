(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*Description*)


(* "FindBounce" package must be loaded before running tests, otherwise procedure is aborted. *)
If[
	Not@MemberQ[$Packages,"FindBounce`"],
	Print["Error: Package is not loaded, try again"];Abort[];
];


(* ::Subsection::Closed:: *)
(*Begin Test Section*)


BeginTestSection["Tests"];


(* ::Subsection::Closed:: *)
(*FindBounce*)


(* ::Subsubsection::Closed:: *)
(*1 field*)


(* Potential of this type (x^4-x^2+x/4) doesn't have exact solution.
Expected output (action value) is also computed with the same function.
Actual and expected output are compared with chosen precision. Therefore we avoid 
some problems of comparing machine precison numbers. Comparing by accuracy 
is not reccomended, because magnitude of action can vary significantly. *)

VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518}]["Action"],
	483.857,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F default (Case B)"
];


VerificationTest[
	FindBounce[x^4-x^2+x/2,x,{-0.809017,0.5}]["Action"],
	16.856,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F default (Case A)"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"Segments"->10]["Action"],
	486.801,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 10 segments"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"Segments"->100]["Action"],
	483.062,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 100 segments"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.12x^4,x,{0.,2.160892},"Dimension"->3]["Action"],
	1067.706,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 3 dimensions (Case B)"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.1x^4,x,{0.,2.882782},"Dimension"->3]["Action"],
	97.156,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 3 dimensions (Case A)"
];


(* ::Subsubsection::Closed:: *)
(*2 fields*)


VerificationTest[
	FindBounce[
		0.1x^4+0.3y^4+2.x^2*y^2-80.x^2-100.y^2,
		{x,y},
		{{0.,12.91},{20.,0.}}
	]["Action"],
	4959.39,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 2F default (Case B)"
];


VerificationTest[
	FindBounce[
		0.1x^4+0.3y^4+2.x^2*y^2-80.x^2-100.y^2+800.y,
		{x,y},
		{{0.,10.},{19.917,-0.577}}
	]["Action"],
	101.195,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 2F default (Case A)"
];


VerificationTest[
	FindBounce[
		0.1x^4+0.3y^4+2.x^2*y^2-80.x^2-100.y^2,
		{x,y},
		{{0.,12.91},{20.,0.}},"Dimension"->3
	]["Action"],
	2336.512,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 2F 3 dimensions (Case B)"
];


VerificationTest[
	FindBounce[
		0.1x^4+0.3y^4+2.x^2*y^2-80.x^2-100.y^2+800.y,
		{x,y},
		{{0.,10.},{19.917,-0.577}},"Dimension"->3]["Action"],
	127.042,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 2F 3 dimensions (Case A)"
];


(* ::Subsubsection::Closed:: *)
(*3 fields*)


VerificationTest[
	FindBounce[
		100.x-40.x^2+0.1x^4+150.y-100.y^2+2.x^2*y^2+0.2y^4+100.z-120.z^2+2.1x^2*z^2+1.2y^2*z^2+0.3z^4,
		{x,y,z},
		{{-0.103556,-16.166160,-0.258163},{13.4220429,-0.2880399,-0.193469}}
	]["Action"],
	275.556,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 3F default"
]


(* ::Subsubsection::Closed:: *)
(*Fail checks*)


(* Slightly moved minima. *)
VerificationTest[
	FindBounce[x^4-x^2+x/4,{x},{-0.762844,0.633518}+0.1],
	$Failed,
	{SingleFieldBounce::extrema},
	TestID->"FindBounce - Error: wrong position of minima"
];


(* ::Subsection::Closed:: *)
(*BouncePlot*)


(* We just check if plot (Graphics expression) is produced without messages. 
Only BouncePlot is inside VerificationTest for relevant timing. *)
With[{
	bf=FindBounce[x^4-x^2+x/4,{x},{-0.762844,0.633518},"Segments"->10]
	},
	VerificationTest[
		BouncePlot[bf,PerformanceGoal->"Speed"],
		_Graphics,
		SameTest->MatchQ,
		TestID->"BouncePlot - 1F default"
	]
];


(* ::Subsection::Closed:: *)
(*End Test Section*)


EndTestSection[];
