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


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518}]["Action"],
	482.955,
	SameTest->(Abs[(#1-#2)/#2]<10^(-2)&),
	TestID->"FindBounce - default"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"Segments"->10]["Action"],
	482.955,
	SameTest->(Abs[(#1-#2)/#2]<10^(-1)&),
	TestID->"FindBounce - 10 segments"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"Segments"->100]["Action"],
	482.955,
	SameTest->(Abs[(#1-#2)/#2]<10^(-3)&),
	TestID->"FindBounce - 100 segments"
];


(* ::Subsubsection::Closed:: *)
(*2 field*)


VerificationTest[
	FindBounce[
		0.1x^4+0.3y^4+2.x^2*y^2-80.x^2-100.y^2,
		{x,y},
		{{0.,12.91},{20.,0.}}
	]["Action"],
	4945.,
	SameTest->(Abs[(#1-#2)/#2]<5.&),
	TestID->"FindBounce - 2F default"
];


(* ::Subsubsection::Closed:: *)
(*Fail checks*)


(* Slightly moved minima. *)
VerificationTest[
	FindBounce[x^4-x^2+x/4,{x},{-0.762844,0.633518}+0.1],
	$Failed,
	{SingleFieldBounce::extrema},
	TestID->"FindBounce - wrong minima"
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
		TestID->"BouncePlot - 1 field default"
	]
];


(* ::Subsection::Closed:: *)
(*End Test Section*)


EndTestSection[];
