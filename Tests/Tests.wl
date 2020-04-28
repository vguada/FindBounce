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
(*Basic syntax*)


(* This potential has minima at 0. and 1. Parameter values (c1) 0.47 and  0.2 represent 
thin-walled and thick-walled bouce respectively. This potential is taken from paper
Chigusa et. al. 2019, https://arxiv.org/abs/1906.10829 *)
singleField[x_,c1_]:=1/4*x^4-(c1+1)/3*x^3+c1/2*x^2 ;

(* Two field potential from the same paper, which has minima at {0,0} and {1,1}.
Parameter values (c2) 2 and  80 represent thin-walled and thick-walled bouce respectively. *)
twoField[x_,y_,c2_]:=(x^2+5y^2)(5(x-1)^2+(y-1)^2)+c2(1/4y^4-1/3y^3) ;


(* In the following tests we check for combinations of special case with single field. 
Only dimensions of the output are checked, but not is values, because we are testing only
for correct processing of input arguments. *)
VerificationTest[
	FindBounce[singleField[x,0.2],x,{0.,1.}]["Coefficients"]//Dimensions,
	{3,_,1},
	SameTest->MatchQ,
	TestID->"FindBounce - 1F basic syntax 1"
];


VerificationTest[
	FindBounce[singleField[x,0.2],x,{{0.},{1.}}]["Coefficients"]//Dimensions,
	{3,_,1},
	SameTest->MatchQ,
	TestID->"FindBounce - 1F basic syntax 2"
];


VerificationTest[
	FindBounce[singleField[x,0.2],{x},{0.,1.}]["Coefficients"]//Dimensions,
	{3,_,1},
	SameTest->MatchQ,
	TestID->"FindBounce - 1F basic syntax 3"
];


VerificationTest[
	FindBounce[singleField[x,0.2],{x},{{0.},{1.}}]["Coefficients"]//Dimensions,
	{3,_,1},
	SameTest->MatchQ,
	TestID->"FindBounce - 1F basic syntax 4"
];


VerificationTest[
	FindBounce[twoField[x,y,80.],{x,y},{{0.,0.},{1.,1.}}]["Coefficients"]//Dimensions,
	{3,_,2},
	SameTest->MatchQ,
	TestID->"FindBounce - 2F basic syntax 1"
];


VerificationTest[
	FindBounce[twoField[x,y,80.],{x,y},{{0.},{1.}}],
	$Failed,
	{FindBounce::mins},
	TestID->"FindBounce - 2F wrong minima dimensions"
];


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
	FindBounce[x^4-x^2+x/2,{x},{{-0.809017},{0.5}}]["Action"],
	16.856,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F default second form (multi-field like)"
];


VerificationTest[
	FindBounce[x^4-x^2+x/2,x,{1/2,1/4 (-1-Sqrt[5])},"FieldPoints"->3,"MidFieldPoint"->1/4 (-1+Sqrt[5])]["Action"],
	2.911,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 2 segments"
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


(* We are using Range instead of Subdivide because the latter is introduced in Mma 10.1 . *)
VerificationTest[
  FindBounce[1/4 x^4 - 4/10 x^3 + 1/10 x^2, x, {0, 1}, "FieldPoints" -> 11]["Action"],
  FindBounce[1/4 x^4 - 4/10 x^3 + 1/10 x^2, x, {0, 1}, "FieldPoints" -> Range[0.,1.,0.1]]["Action"],
  SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
  TestID -> "FindBounce - integer vs. list of \"FieldPoints\""
];


VerificationTest[
  FindBounce[1/4 x^4 - 4/10 x^3 + 1/10 x^2, x, {0, 1}, "FieldPoints" -> 11,"Gradient"->None]["Action"],
  FindBounce[1/4 x^4 - 4/10 x^3 + 1/10 x^2, x, {0, 1}, "FieldPoints" -> Range[0.,1.,0.1],"Gradient"->None]["Action"],
  SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
  TestID -> "FindBounce - integer vs. list of \"FieldPoints\", \"Gradient\"->None"
];


VerificationTest[
	FindBounce[0.5x^2+0.5x^3+0.1x^4,x,{0.,-2.883}]["CoefficientsExtension"][[1,1,1]],
	-0.05486,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F 3 CofficientsExtension (Case A,D=4)"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.1x^4,x,{0.,2.883},"Dimension"->3]["CoefficientsExtension"][[1,1,1]],
	0.0025433,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F CofficientsExtension (Case A,D=3)"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.12x^4,x,{0.,2.160892}]["CoefficientsExtension"][[1,1,1]],
	1610.0432764,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F CofficientsExtension (Case B,D=4)"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.12x^4,x,{0.,2.160892},"Dimension"->3]["CoefficientsExtension"][[1,1,1]],
	131.40002,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F CofficientsExtension (Case B,D=3)"
];


VerificationTest[
	FindBounce[0.5x^2-0.5x^3+0.1x^4,x,{0.,2.882782}]["CoefficientsExtension"][[1,1,1]],
	-FindBounce[0.5x^2+0.5x^3+0.1x^4,x,{0.,-2.882782}]["CoefficientsExtension"][[1,1,1]],
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - 1F CofficientsExtension Symmetry"
];


(* ::Subsubsection::Closed:: *)
(*1 field - thin vs. thick wall regime*)


With[{a=0.5},
	VerificationTest[
		FindBounce[0.5*x^2+0.5*x^3+0.125*a*x^4,x,{-5.236,0.}]["Action"],
		269.693,
		SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
		TestID->"FindBounce - 1F thick wall"
	]
];


With[{a=0.99},
	VerificationTest[
		FindBounce[0.5*x^2+0.5*x^3+0.125*a*x^4,x,{-2.04,0.}]["Action"],
		3.6305*10^6,
		SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
		TestID->"FindBounce - 1F thin wall"
	]
];


(* ::Subsubsection::Closed:: *)
(*1 field - bi-quartic potential*)


(* This type of biquartic potential has exact solution. *)
biQuarticPotential[{\[Phi]p_,\[Phi]m_},{Vt_,Vm_},\[Epsilon]_]:=Module[
	{Vp,\[CapitalDelta]Vp,\[CapitalDelta]Vm,\[Alpha],\[CapitalDelta]},
	Vp=-18+(\[Epsilon]-1)*3;
	\[CapitalDelta]Vm=Vt-Vm;
	\[CapitalDelta]Vp=Vt-Vp;
	\[Alpha]=-\[Phi]p/\[Phi]m;
	\[CapitalDelta]=\[CapitalDelta]Vp/\[CapitalDelta]Vm;
	Function[
		\[Phi],
		Evaluate@Piecewise[
			{{(Vp+\[CapitalDelta]Vp/\[Phi]p^4*(\[Phi]-\[Phi]p)^4)-Vp,\[Phi]<0}},
			(Vm+\[CapitalDelta]Vm/\[Phi]m^4*(\[Phi]-\[Phi]m)^4)-Vp]
	]
];


biQuarticAction[{\[Phi]p_,\[Phi]m_},{Vt_,Vm_},\[Epsilon]_]:=Module[
	{Vp,\[CapitalDelta]Vp,\[CapitalDelta]Vm,\[Alpha],\[CapitalDelta]},
	Vp=-18+(\[Epsilon]-1)*3;
	\[CapitalDelta]Vm=Vt-Vm;
	\[CapitalDelta]Vp=Vt-Vp;
	\[Alpha]=-\[Phi]p/\[Phi]m;
	\[CapitalDelta]=\[CapitalDelta]Vp/\[CapitalDelta]Vm;
	((2*\[Pi]^2*\[Phi]m^4)/(3\[CapitalDelta]Vm))*(4\[Alpha]^3 + 6\[Alpha]^2*\[CapitalDelta]+4\[Alpha]*\[CapitalDelta]^2+\[CapitalDelta]^3+\[Alpha]^4*(3+\[CapitalDelta](\[CapitalDelta]-3)))/(1-\[CapitalDelta])^3
];


VerificationTest[
	FindBounce[biQuarticPotential[{-5,10},{5,-20},2][x],x,{-5,10}],
	$Failed,
	{FindBounce::gradval},
	TestID->"FindBounce - 1F biquartic (nonexistent gradient)"
];


(* Actual output (action value) is compared to exact analytical value of expected output.
Precision is chosen to satisfy the test. *)
VerificationTest[
	FindBounce[biQuarticPotential[{-5,10},{5,-20},2][x],x,{-5,10},"Gradient"->None]["Action"],
	biQuarticAction[{-5,10},{5,-20},2],
	SameTest->(Abs[(#1-#2)/#2]<10^(-1)&),
	TestID->"FindBounce - 1F biquartic (default)"
];


(* TODO: Why does the message FindBounce::gradfail not appear in this case?*)
VerificationTest[
	FindBounce[
		biQuarticPotential[{-5,10},{5,-20},2][x],x,{-5,10},
		"FieldPoints"->71
	]["Action"],
	biQuarticAction[{-5,10},{5,-20},2],
	SameTest->(Abs[(#1-#2)/#2]<10^(-2)&),
	TestID->"FindBounce - 1F biquartic (70 segments)"
];


(* ::Subsubsection::Closed:: *)
(*1 field - potential defined by points*)


(* The simplest potential defined by points. *)
VerificationTest[
	FindBounce[{{1,0},{2,2},{3,1}}]["Action"],
	3576.28,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - simplest polygonal potential"
];


(* The simplest potential defined by points and D=3. *)
VerificationTest[
	FindBounce[{{1,0},{2,2},{3,1}},"Dimension"->3]["Action"],
	196.262,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - simplest polygonal potential, D=3"
];


(* Polygonal potential with 4 points and D=4. *)
VerificationTest[
	FindBounce[{{0,0},{2,2},{2.5,3},{3,1}}]["Action"],
	32937.13,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - N = 4 polygonal potential"
];


(* Polygonal potential with 4 points and D=3. *)
VerificationTest[
	FindBounce[{{0,0},{2,2},{2.5,3},{3,1}},"Dimension"->3]["Action"],
	1039.70,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - N = 4 polygonal potential, D=3"
];


(* Polygonal potential with intermediate minima. *)
VerificationTest[
	FindBounce[{{1,0},{2,2},{3,1},{4,2},{5,1}}]["Action"],
	82767.,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - intermediate minima polygonal potential"
];


(* Not enough points to define the potential. *)
VerificationTest[
	FindBounce[{{1,0},{2,2}}],
	$Failed,
	{FindBounce::points},
	TestID->"FindBounce - polygonal, not enough points"
];


(* Wrong dimensions of the array defining the potential. *)
VerificationTest[
	FindBounce[{{1,0},{2,2},{3,1},{4,2},{5,1,100}}],
	$Failed,
	{FindBounce::points},
	TestID->"FindBounce - polygonal, wrong array dimensions"
];


(* Array of points should not include Complex or non-numerical values. *)
VerificationTest[
	FindBounce[{{1,0},{2,2},{3,1},{4,2},{5,1+I}}],
	$Failed,
	{FindBounce::points},
	TestID->"FindBounce - polygonal, non-numerical values"
];


(* ::Subsubsection::Closed:: *)
(*1 field - potential unbounded from below*)


VerificationTest[
	FindBounce[0.5 x^2+0.05 x^3-0.125 x^4,x,{0,-5},"BottomlessPotential"->True]["Action"],
	51.64,
	SameTest->(Abs[(#1-#2)/#2]<10^(-4)&),
	TestID->"FindBounce - bottomles potential"
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


(* Field symbols have some value, which they should't have. *)
VerificationTest[
	Block[{x=1},FindBounce[x^4-x^2+x/4,{x},{-0.762,0.633}]],
	$Failed,
	{FindBounce::syms},
	TestID->"FindBounce - error: field symbol value"
];


(* Slightly moved minima. *)
VerificationTest[
	FindBounce[x^4-x^2+x/4,{x},{-0.762844,0.633518}+0.1],
	$Failed,
	{FindBounce::extrema},
	TestID->"FindBounce - Error: wrong position of minima"
];


(* Both minima of potential have the same value. This is degenerated case. *)
VerificationTest[
	FindBounce[x^4-x^2,{x},{-0.707107,0.707107}]["Action"],
	\[Infinity],
	{FindBounce::degeneracy},
	TestID->"FindBounce - Error: degenerated minima"
]


(* The minima are the same. *)
VerificationTest[
	FindBounce[x^4-x^2,{x},{0.707107,0.707107}],
	$Failed,
	{FindBounce::mins},
	TestID->"FindBounce - Error: minima are the same"
];


(* The minima are complex - we want that function explicitly fails. *)
VerificationTest[
	FindBounce[x^4-x^2+x/4,{x},{-0.762,0.633+0.01 I}],
	$Failed,
	{FindBounce::mins},
	TestID->"FindBounce - Error: minima are complex"
];


(* ::Subsubsection::Closed:: *)
(*Fail checks - wrong option values*)


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"PathTolerance"->-1],
	$Failed,
	{FindBounce::posreal},
	TestID->"FindBounce - wrong PathTolerance"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"ActionTolerance"->-1],
	$Failed,
	{FindBounce::posreal},
	TestID->"FindBounce - wrong InitialRadiusAccuracyGoal"
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
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"FieldPoints"->2],
	$Failed,
	{FindBounce::fieldpts},
	TestID->"FindBounce - wrong number of FieldPoints"
];


VerificationTest[
	FindBounce[
		singleField[x,0.2],x,{0.,1.},
		"FieldPoints"->{0.,0.2,0.4,0.6,0.8,1.}
	]["Coefficients"]//Dimensions,
	{3,_,1},
	SameTest->MatchQ,
	TestID->"FindBounce - FieldPoints as vector"
];


VerificationTest[
	FindBounce[
		singleField[x,0.2],x,{0.,1.},
		"FieldPoints"->{{0.},{0.2},{0.4},{0.6},{0.8},{1.}}
	]["Coefficients"]//Dimensions,
	{3,_,1},
	SameTest->MatchQ,
	TestID->"FindBounce - FieldPoints as matrix"
];


VerificationTest[
	FindBounce[
		singleField[x,0.2],x,{0.,1.},
		"FieldPoints"->{{0.},{0.2},{0.4},{0.6},{0.8},{1.+I}}
	],
	$Failed,
	{FindBounce::fieldpts},
	TestID->"FindBounce - Error: FieldPoints with complex number"
]


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"MaxRadiusIterations"->0],
	$Failed,
	{FindBounce::posint},
	TestID->"FindBounce - wrong MaxRadiusIterations"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762844,0.633518},"MaxPathIterations"->0.1],
	$Failed,
	{FindBounce::nonnegint},
	TestID->"FindBounce - wrong MaxPathIterations"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762,0.633},"FieldPoints"->{0.762,0.633}],
	$Failed,
	{FindBounce::fieldpts},
	TestID->"FindBounce - wrong array length FieldPoints"
];


VerificationTest[
	FindBounce[x^4-x^2+x/4,x,{-0.762,0.633},"FieldPoints"->RandomReal[{0,1},{5,2}]],
	$Failed,
	{FindBounce::fieldpts},
	TestID->"FindBounce - wrong array shape FieldPoints"
];


VerificationTest[
	FindBounce[singleField[x,0.2],x,{0.,1.},"Gradient"->"badValue"],
	$Failed,
	{FindBounce::gradmtd},
	TestID->"FindBounce - wrong \"Gradient\" option value"
];


VerificationTest[
	FindBounce[twoField[x,y,80.],{x,y},{{0.,0.},{1.,1.}},"Hessian"->"badValue"],
	$Failed,
	{FindBounce::hessmtd},
	TestID->"FindBounce - wrong \"Hessian\" option value"
];


(* ::Subsection::Closed:: *)
(*BounceFunction*)


(* We can use the same BounceFunction for mutiple tests of properties/keys. *)
With[{
	bf=FindBounce[x^4-x^2+x/4,x,{-0.762,0.633}]
	},
	VerificationTest[
		bf["Properties"],
		{"Action","BottomlessPotential","Bounce","Coefficients","CoefficientsExtension","Dimension","Path","PathIterations","Radii"},
		TestID->"BounceFunction - all available properties"
	];
	
	VerificationTest[
		bf["BadValue"],
		_Missing,
		{BounceFunction::noprop},
		SameTest->MatchQ,
		TestID->"BounceFunction - non-existent property"
	];

	VerificationTest[
		bf["BottomlessPotential"],
		Missing["NotAvailable"],
		TestID->"BounceFunction - valid Missing value"
	]
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
		BouncePlot[bf,"BadOption"->1],
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


(* Test for arbitrary plot range if flat parts of curve are ploted fully.
Negative values are not plotted, becasue function is not supported there. *)
With[{
	bf=FindBounce[{{1,0},{2,2},{3,1}}]
	},
	VerificationTest[
		(Join@@Cases[
			BouncePlot[bf,PlotRange->{{-5,10},Automatic}],
			Line[x_]:>x,
			Infinity
		])[[{1,-1},1]],
		{0.,10.},
		SameTest->(Norm[#1-#2]<10^-4&),
		TestID->"BouncePlot - arbitrary plot range"
	]
];


(* Test for case of degenerated vacua. *)
With[{
	bf=Quiet[
			FindBounce[-1/2*x^2+1/8*x^4 ,{x},{-1.414,1.414}],
			FindBounce::degeneracy
		]
	},
	VerificationTest[
		BouncePlot[bf],
		_Graphics,
		SameTest->MatchQ,
		TestID->"BouncePlot - degenerated vacua"
	]
];


(* ::Subsection::Closed:: *)
(*End Test Section*)


EndTestSection[];
