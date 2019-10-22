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
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"FieldPoints"->11]["Action"],
	486.801,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 10 segments"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"FieldPoints"->101]["Action"],
	483.062,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 100 segments"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.12x^4,x,{0.,2.160892},"Dimension"->3]["Action"],
	1051.704,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 3 dimensions (Case B)"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.1x^4,x,{0.,2.882782},"Dimension"->3]["Action"],
	95.116,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 3 dimensions (Case A)"
];


(* ::Subsubsection::Closed:: *)
(*1 field - bi-quartic potential*)


(* This type of biquartic potential has exact solution. *)
Module[
	{Vt=5,\[Phi]m=10,Vm=-20,\[Phi]p=-5,\[Epsilon]=2,Vp,\[CapitalDelta]Vp,\[CapitalDelta]Vm,\[Alpha],\[CapitalDelta]},
	Vp = -18+(\[Epsilon]-1)*3;
	\[CapitalDelta]Vm=Vt-Vm;
	\[CapitalDelta]Vp= Vt - Vp;
	\[Alpha] = -\[Phi]p/\[Phi]m;
	\[CapitalDelta] = \[CapitalDelta]Vp/\[CapitalDelta]Vm;
	biQuarticPotential[\[Phi]_]:=Piecewise[{
		{(Vp+\[CapitalDelta]Vp/\[Phi]p^4*(\[Phi]-\[Phi]p)^4)-Vp,\[Phi]<0},
		{(Vm+\[CapitalDelta]Vm/\[Phi]m^4*(\[Phi]-\[Phi]m)^4)-Vp,\[Phi]>=0}
	}];
	biQuarticAction[]:=((2*\[Pi]^2*\[Phi]m^4)/(3\[CapitalDelta]Vm))*(4\[Alpha]^3 + 6\[Alpha]^2*\[CapitalDelta]+4\[Alpha]*\[CapitalDelta]^2+\[CapitalDelta]^3+\[Alpha]^4*(3+\[CapitalDelta](\[CapitalDelta]-3)))/(1-\[CapitalDelta])^3;
];


(* Actual output (action value) is compared to exact analytical value of expected output.
Precision is chosen to satisfy the test. *)
VerificationTest[
	FindBounce[biQuarticPotential[x],{x},{-5,10}]["Action"],
	biQuarticAction[],
	{SingleFieldBounceImprovement::dVFailed},
	SameTest->(Abs[(#1-#2)/#2]<10^(-1)&),
	TestID->"FindBounce - 1F biquartic (default)"
];


(* TODO: Why does the message SingleFieldBounceImprovement::dVFailed not appear in this case?*)
VerificationTest[
	FindBounce[
		biQuarticPotential[x],{x},{-5,10},
		"FieldPoints"->71
	]["Action"],
	biQuarticAction[],
	SameTest->(Abs[(#1-#2)/#2]<10^(-2)&),
	TestID->"FindBounce - 1F biquartic (70 segments)"
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
	275.521,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 3F default"
];


(* ::Subsubsection::Closed:: *)
(*Fail checks*)


(* Slightly moved minima. *)
VerificationTest[
	FindBounce[x^4-x^2+x/4,{x},{-0.762844,0.633518}+0.1],
	$Failed,
	{SingleFieldBounce::extrema},
	TestID->"FindBounce - Error: wrong position of minima"
];


(* Both minima of potential have the same value. This is degenerated case. *)
VerificationTest[
	FindBounce[x^4-x^2,{x},{-0.707107,0.707107}]["Action"],
	\[Infinity],
	{FindBounce::degeneracy},
	TestID->"FindBounce - Error: degenerated minima"
];


(* ::Subsubsection::Closed:: *)
(*Fail checks - wrong option values*)


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"InitialRadiusAccuracyGoal"->0],
	$Failed,
	{FindBounce::posint},
	TestID->"FindBounce - wrong InitialRadiusAccuracyGoal"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"PathTolerance"->-0.1],
	$Failed,
	{FindBounce::posreal},
	TestID->"FindBounce - wrong PathTolerance"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"ActionTolerance"->"string"],
	$Failed,
	{FindBounce::posreal},
	TestID->"FindBounce - wrong ActionTolerance"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"Dimension"->2],
	$Failed,
	{FindBounce::dim},
	TestID->"FindBounce - wrong Dimension"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"InitialRadius"->0],
	$Failed,
	{FindBounce::posreal},
	TestID->"FindBounce - wrong InitialRadius"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"MaxRadiusIterations"->0],
	$Failed,
	{FindBounce::posint},
	TestID->"FindBounce - wrong MaxRadiusIterations"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"MaxPathIterations"->0],
	$Failed,
	{FindBounce::posint},
	TestID->"FindBounce - wrong MaxPathIterations"
];


(* ::Subsection::Closed:: *)
(*BouncePlot*)


(* We just check if plot (Graphics expression) is produced without messages. 
Only BouncePlot is inside VerificationTest for relevant timing. We also check
preventing wrong options to reach Plot. *)
With[{
	bf=FindBounce[x^4-x^2+x/4,{x},{-0.762844,0.633518},"FieldPoints"->11]
	},
	VerificationTest[
		BouncePlot[bf,PerformanceGoal->"Speed","BadOption"->1],
		_Graphics,
		SameTest->MatchQ,
		TestID->"BouncePlot - 1F default"
	]
];


(* Check that function returns unevaluated if the argument is not BounceFunction or a list of them. *)
VerificationTest[
	BouncePlot["badValue"],
	_BouncePlot,
	SameTest->MatchQ,
	TestID->"BouncePlot - returns unevaluated"
];


(* Check that BounceFunction of 2 fields is plotted with 2 colors by default. *)
With[{
	bf2=FindBounce[0.1x^4+0.3y^4+2.x^2*y^2-80.x^2-100.y^2,{x,y},{{0.,12.91},{20.,0.}},"FieldPoints"->11]
	},
	VerificationTest[
		Length@Cases[
			BouncePlot[bf2,PerformanceGoal->"Speed"],
			_RGBColor,
			Infinity
		],
		2,
		TestID->"BouncePlot - 2F default"
	]
];


(* Check if two BounceFunctions are plotted with specified colors. *)
With[{
	bf=FindBounce[x^4-x^2+x/4,{x},{-0.762844,0.633518},"FieldPoints"->11]
	},
	VerificationTest[
		Cases[
			BouncePlot[{bf,bf},PlotStyle->{Red,Directive[Dashed,Blue]},PerformanceGoal->"Speed"],
			_RGBColor,
			Infinity
		],
		{Red,Blue},
		TestID->"BouncePlot - 2 functions"
	]
];


(* ::Subsection::Closed:: *)
(*End Test Section*)


EndTestSection[];
