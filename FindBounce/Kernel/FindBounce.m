(* ::Package:: *)

(* ::Section::Closed:: *)
(*Header*)


(* :Title: FindBounce *)
(* :Context: FindBounce` *)
(* :Author: Victor Guada, Miha Nemevsek and Matevz Pintar *)
(* :Summary: Computes decay of the false vacuum in models with multiple scalar. *)
(* :Keywords: tunneling, first order phase transitions, bubble nucleation *)


(*  Copyright (C) 2019

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License aM published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*)


(* ::Section::Closed:: *)
(*Begin Package*)


BeginPackage["FindBounce`"];


(* ::Subsection::Closed:: *)
(*Available public functions*)


(* The main function of the package. 
Messages from all private functions it contains are attached to this symbol, 
so that the user is not confused by the names of internal functions appearing. *)
FindBounce;
(* Symbol which acts as a container for results of FindBounce. *)
BounceFunction;
(* A handy function for plotting the contents of BounceFunction. *)
BouncePlot;


(* Clear definitions from package symbols in public and private context. *)
ClearAll["`*","`*`*"];


(* ::Section::Closed:: *)
(*Code*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Utilities*)


(* ::Subsubsection::Closed:: *)
(*Autocomplete arguments*)


(* 
https://resources.wolframcloud.com/FunctionRepository/resources/AddCodeCompletion
https://mathematica.stackexchange.com/questions/56984
*)
addCodeCompletion[function_String][args___]:=Module[
	{processed},
	processed={args}/.{
		None -> 0,
		"AbsoluteFileName" -> 2,
		"RelativeFileName" -> 3,
		"Color" -> 4,
		"PackageName" -> 7,
		"DirectoryName" -> 8,
		"InterpreterType" -> 9
	};
	(FE`Evaluate[FEPrivate`AddSpecialArgCompletion[#1]]&)[function -> processed]
];


(* ::Subsubsection::Closed:: *)
(*Version compatibility*)


(* Manual implementation of MinMax function, because it has been only added in Mathematica 10.1
It matches the original in all examples from documentation. *)
minMax[list_List]:={Min[list],Max[list]};

minMax[list_List,d_]:={Min[list]-d,Max[list]+d};

minMax[list_List,Scaled[s_]]:=With[
	{d=s*(Max[list]-Min[list])},
	{Min[list]-d,Max[list]+d}
];

minMax[list_List,{dMin_,dMax_}]:={Min[list]-dMin,Max[list]+dMax};


subdivide[x1_,x2_,n_Integer]:=Table[x1+(x2-x1)*k/n,{k,0,n}];


(* Use user-defined function only in earlier versions. *)
If[
	$VersionNumber<10.1,
	MinMax=minMax;
	Subdivide=subdivide
];


(* ::Subsection::Closed:: *)
(*Initial calculations*)


(* ::Subsubsection::Closed:: *)
(*Utilites for working with field points*)


(* This efficient implementation keeps a PackedArray packed. *)
segmentsLength[fieldPoints_]:=Sqrt@Total[Differences[fieldPoints]^2,{2}];


(* Project field points from multi field to single field, called "longitudinal bounce" and return a vector. *)
longitudinalProjection[fieldPoints_]:=Join[{0.},Accumulate@segmentsLength[fieldPoints]];


(* This is used to project coefficents from longitudinal bounce back to true field space. *)
unitVectors[fieldPoints_]:=Differences[fieldPoints]/segmentsLength[fieldPoints];


(* ::Subsubsection::Closed:: *)
(*Segmentation*)


(* Homogeneous and bi-homogeneous segmentation between 2 or 3 points for any dimension. *)
homogeneousSegmentation[{p1_,p2_},n_Integer/;n>=2]:=Subdivide[p1,p2,n];

homogeneousSegmentation[{p1_,p2_,p3_},n_Integer/;n>=2]:=Module[
	{lengths,n1,n2},
	lengths=Norm/@Differences[{p1,p2,p3}];
	n1=Clip[
		Round[First[lengths]/Total[lengths]*n],
		{1,n-1}
	];
	n2=n-n1;
	Join[Subdivide[p1,p2,n1],Rest@Subdivide[p2,p3,n2]]
];


FindBounce::fieldpts = "\"FieldPoints\" should be an integer (n>2) or array of non-complex numbers longer than 2.";
FindBounce::midpt="\"MidFieldPoint\" should be a vector of length equal to the number of fields or symbol None.";
FindBounce::minpos="The first and the last point in given list of field points have to match the minima.";

(* Returns either a matrix representing field segmentation or $Failed. *)
fieldSegmentation[{min1_,min2_},opts:OptionsPattern[]]:=Module[
	{noFields,fieldPoints,midPoint},
	(* We assume that dimension of minima is already processed before. *)
	noFields=Length[min1];

	midPoint=N@OptionValue[FindBounce,{opts},"MidFieldPoint"];
	If[midPoint=!=None&&Length[midPoint]==0,midPoint={midPoint}];
	(* "MidFieldPoint" value has to be either None or numeric vector. *)
	If[
		midPoint=!=None,
		If[
			Not@And[
				ArrayQ[N@midPoint,1,(MatchQ[#,_Real]&)],
				Length[midPoint]===noFields
			],
			Message[FindBounce::midpt];Return[$Failed,Module]
		]
	];
	fieldPoints=OptionValue[FindBounce,{opts},"FieldPoints"];

	Which[
		IntegerQ[fieldPoints],
		If[fieldPoints<3,Message[FindBounce::fieldpts];Return[$Failed,Module]];
		fieldPoints=If[
			midPoint===None,
			homogeneousSegmentation[N@{min1,min2},fieldPoints-1],
			homogeneousSegmentation[N@{min1,midPoint,min2},fieldPoints-1]
		]
		,
		ListQ[fieldPoints],
		fieldPoints=N@fieldPoints;
		If[
			Depth[fieldPoints]==2&&noFields==1,
			fieldPoints=Partition[fieldPoints,1]
		];
		If[
			Not@And[
				ArrayQ[fieldPoints,2,(MatchQ[#,_Real]&)],
				MatchQ[Dimensions@fieldPoints,{x_/;x>=3,noFields}]
			],
			Message[FindBounce::fieldpts];Return[$Failed,Module]
		];
		(* The first and last point in given list of field points have to be the same as
		minimum 1 and 2 respectively. *)
		If[
			fieldPoints[[{1,-1}]]!={min1,min2},
			Message[FindBounce::minpos];Return[$Failed,Module]
		]
		,
		True,
		Message[FindBounce::fieldpts];Return[$Failed,Module]
	];

	Developer`ToPackedArray@fieldPoints
];


(* ::Subsubsection::Closed:: *)
(*Potential values*)


(* Elegant implementation of helper function for replacement rules, which can be used 
throughout the whole package. *)
replaceValues[fields_,fieldPoints_]:=Transpose@MapThread[
	Thread@*Rule,
	{fields,Transpose[fieldPoints]}
];


FindBounce::potvals="Potential values should be real numbers for all field points.";
FindBounce::extrema="Wrong position of the extrema.
Check the minima or use \"MidFieldPoint\" option to include the maximum/saddle point of the potential.";

(* We assume that we already get consistent dimensions of fields and fieldPoints. *)
getPotentialValues[V_,fields_,fieldPoints_]:=Module[
	{values},
	values=V/.replaceValues[fields,fieldPoints];
	If[
		Not@VectorQ[values,(MatchQ[#,_Real]&)],
		Message[FindBounce::potvals];Return[$Failed,Module]
	];
	If[
		Or[values[[1]]>=values[[2]],values[[-1]]>=values[[-2]]],
		Message[FindBounce::extrema];Return[$Failed,Module]
	];
	Developer`ToPackedArray@values
];


(* ::Subsubsection::Closed:: *)
(*Potenital gradient*)


symbolicGradient[V_,fields_,fieldPoints_]:=Module[
	{},
	Grad[V,fields]/.replaceValues[fields,fieldPoints]
];


numericGradient//Options={"Scale"->10^-4};

numericGradient[V_,fields_,fieldPoints_,opts:OptionsPattern[]]:=Module[
	{pathLength,gradients,noFields,noPts,eps,n,Vt},
	pathLength=Total@segmentsLength[fieldPoints];
	eps = OptionValue["Scale"]*pathLength;
	noFields =Length[fields];
	noPts = Length[fieldPoints];
	n = IdentityMatrix[noFields];

	Vt = V/.replaceValues[fields,fieldPoints];
	gradients = ConstantArray[0.,{noPts,noFields}];
	Do[
		gradients[[s,i]]=((V/.Thread[fields->fieldPoints[[s]]+eps*n[[i]]])-Vt[[s]])/eps,
		{s,noPts},{i,noFields}
	];
	gradients
];


(* This processes Gradient option and should work similarly to Gradient option in FindMinimum
and related functions, see
http://reference.wolfram.com/language/tutorial/UnconstrainedOptimizationSpecifyingDerivatives.html.
User can choose between symbolic gradient calculation, finite differences or custom function.
Hessian calculation is part of another independent function, because the FindBounce options
are independent. *)

(* TODO: There is a problem with potentials given as functions with ?NumericQ argument
check. Function D or Grad can actually calculate numerical derivatives
(probably some FDM in the background) but this can be quite slow. It would be better
if "FiniteDifference" option value would be automatically chosen in such cases. *)

FindBounce::gradmtd="Option value for Gradient should be Automatic, None, \"Symbolic\", \"FiniteDifference\" of custom function.";
FindBounce::gradval=(
	"The gradient of the potential is not well defined at some field point."<> 
	"Redefine the potential, choose option \"Gradient\"->None or \"Gradient\"->\"FiniteDifference\".");

getPotentialGradient[V_,fields_,fieldPoints_,opts:OptionsPattern[]]:=Module[
	{optValue,method,gradient},
	optValue=OptionValue[FindBounce,{opts},"Gradient"]/.Automatic->"Symbolic";
	method=Which[
		optValue==="Symbolic","Symbolic",
		optValue==="FiniteDifference","Numeric",
		MatchQ[Head@optValue,_Symbol]&&Not@StringQ[optValue],"Custom",
		True,Message[FindBounce::gradmtd];Return[$Failed,Module]
	];
	gradient=Switch[method,
		"Symbolic",symbolicGradient[V,fields,fieldPoints],
		"Numeric",numericGradient[V,fields,fieldPoints],
		"Custom",optValue/.replaceValues[fields,fieldPoints]
	];

	If[
		Not@ArrayQ[gradient,2,NumericQ],
		Message[FindBounce::gradval];Return[$Failed,Module]
	];
	gradient
];


(* ::Subsubsection::Closed:: *)
(*Potenital hessian*)


symbolicHessian[V_,fields_,fieldPoints_]:=Module[
	{},
	D[V,{fields,2}]/.replaceValues[fields,fieldPoints]
];


numericHessian//Options={"Scale"->10^-4};

numericHessian[V_,fields_,fieldPoints_,opts:OptionsPattern[]]:=Module[
	{pathLength,hessian,hessianM,noFields,noPts,eps,n,Vt},
	pathLength=Total@segmentsLength[fieldPoints];
	eps = OptionValue["Scale"]*pathLength;
	noFields =Length[fields];
	noPts = Length[fieldPoints];
	n = IdentityMatrix[noFields];

	Vt = V/.replaceValues[fields,fieldPoints];
	hessianM = ConstantArray[0.,{noPts,noFields,noFields}];

	(* diagonal elements *)
	Do[
		hessianM[[s,i,i]] = 
			((V/.Thread[fields->fieldPoints[[s]]+eps(2n[[i]])])- 
			2 Vt[[s]]+ 
			(V/.Thread[fields->fieldPoints[[s]]+eps(-2n[[i]])]))/(8 eps^2),
		{s,noPts},{i,noFields}
	];

	(* upper-triangle elements *)
	Do[
		hessianM[[s,i,j]] = (
			(V/.Thread[fields->fieldPoints[[s]]+eps(n[[i]]+n[[j]])])- 
			(V/.Thread[fields->fieldPoints[[s]]+eps(-n[[i]]+n[[j]])])- 
			(V/.Thread[fields->fieldPoints[[s]]+eps(n[[i]]-n[[j]])])+  
			(V/.Thread[fields->fieldPoints[[s]]+eps(-n[[i]]-n[[j]])]) )/(4 eps^2),
		{s,noPts},{i,noFields},{j,i+1,noFields}
	];
	(* full matrix assembly *)
	hessian = Table[(hessianM[[s]]+Transpose[hessianM[[s]]]),{s,noPts}];
	hessian
];


(* See comments for function calculate gradients of the potential. *)

FindBounce::hessmtd="Option value for Hessian should be Automatic, None, \"Symbolic\", \"FiniteDifference\" of custom function.";
FindBounce::hessval=(
	"The hessian of the potential is not well defined at some field point."<> 
	"Redefine the potential or choose option \"Hessian\"->\"FiniteDifference\".");

getPotentialHessian[V_,fields_,fieldPoints_,opts:OptionsPattern[]]:=Module[
	{optValue,method,hessians},
	(* Inclusion of None is a quick fix for default options of FindBounce
	 and should be eventually removed. *)
	optValue=OptionValue[FindBounce,{opts},"Hessian"]/.(Automatic|None)->"Symbolic";
	method=Which[
		optValue==="Symbolic","Symbolic",
		optValue==="FiniteDifference","Numeric",
		MatchQ[Head@optValue,_Symbol]&&Not@StringQ[optValue],"Custom",
		True,Message[FindBounce::hessmtd];Return[$Failed,Module]
	];

	hessians=Switch[method,
		"Symbolic",symbolicHessian[V,fields,fieldPoints],
		"Numeric",numericHessian[V,fields,fieldPoints],
		"Custom",optValue/.replaceValues[fields,fieldPoints]
	];

	If[
		Not@ArrayQ[hessians,3,NumericQ],
		Message[FindBounce::hessval];Return[$Failed,Module]
	];
	hessians
];


(* ::Subsubsection::Closed:: *)
(*Estimate initial radius*)


(* Estimates the initial radius with the N=2 closed form solution. *)
initialRadiusClosedForm[{{x1_,y1_},{x2_,y2_},{x3_,y3_}},dim_]:=Module[
	{a1,a2,c,f0},
	If[y1==y3,Return[Infinity,Module]];
	{a1,a2}={(y2-y1)/(x2-x1),(y3-y2)/(x3-x2)}/8;
	c = dim/(dim-2)*(a2-a1)/a1*(1 -(a2/(a2-a1))^(1-2/dim));
	f0 =(x3+c*x2)/(1+c);
	If[
		Re[f0]>=0&&Im[f0]==0,
		Sqrt[dim/4(x2-f0)/a1],	
		1/2(x3-x1)/(Sqrt[a1*(x2-x1)]-Sqrt[-a2*(x3-x2)])
	]
];


(* Estimate initial radius with closed form solution (for 3 points) of projection to single field. *)
initialRadiusEstimate[fieldPoints_,potentialValues_,dim_]:=Module[
	{projection,maxPos,x3,y3,pts,list},
	projection=longitudinalProjection[fieldPoints];
	maxPos=Position[potentialValues,Max@potentialValues[[2;;-2]]][[1,1]];
	x3=projection[[{1,maxPos,-1}]];
	y3=potentialValues[[{1,maxPos,-1}]];
	(*We want to have the local minimum always on the right hand side.*) 
	pts=If[
		y3[[3]]>y3[[1]],
		Transpose[{x3,y3}],
		Transpose[{x3,Reverse@y3}]
	];

	initialRadiusClosedForm[pts,dim]
];


(* ::Subsection::Closed:: *)
(*SingleFieldBounce*)


(* ::Subsubsection::Closed:: *)
(*FindSegment*)


FindSegment[a_,\[Phi]L_,d_,Ns_]:=
Module[{pos = 1,R},
	R = BounceParameterRvb[0,\[Phi]L,a,d,Ns,False,pos][[1]];
	
	While[Im[R[[-1]]] == 0 && pos< Ns,
		pos++;
		R = BounceParameterRvb[0,\[Phi]L,a,d,Ns,False,pos][[1]]; 
	];
	
	pos
];


(* ::Subsubsection::Closed:: *)
(*RInitial*)


(*This function return the second or second last Radii, i.e R[[pos+1]] or R[[-2]] if D = 4*)
RInitial[4,initialR_?NumericQ,a_,\[Phi]L_,pos_,backward_]:=
	If[backward,
		Sqrt[initialR^2 +  Sqrt[(\[Phi]L[[-2]]- \[Phi]L[[-1]])initialR^2 /a[[-2]]]],
		(*else---------*)
		Sqrt[initialR^2 + (a[[pos]]initialR^2 - Sqrt[a[[pos]]^2 initialR^4+4(a[[pos+1]] - 
		a[[pos]])(\[Phi]L[[pos+1]] - \[Phi]L[[pos]])initialR^2   ]  )/(2(a[[pos+1]] -a[[pos]])) ] 				        
	 ];
RInitial[3,initialR_?NumericQ,a_,\[Phi]L_,pos_,backward_]:= initialR;


(* ::Subsubsection::Closed:: *)
(*BounceParameterRvb*)


(*See eqs. 18-20*)
Rs[4,c1_?NumericQ,a_,b_] := Sqrt[ 1/2 (Sqrt[ c1^2 - 4*a*b ] + c1)/a  ];
Rs[3,c1_?NumericQ,a_,b_] := Module[{\[Xi]},  
	\[Xi]= ( Sqrt[ 36*a*b^2 - c1^3 ] - 6*b*a^(1/2) )^(1/3) /a^(1/2);
	If[Re@\[Xi]==0,
		\[Xi] = (Sqrt[ 36*a*b^2 - c1^3 ] + 6*b*a^(1/2) )^(1/3) /a^(1/2)
	]; 

	1/2 (c1/a/\[Xi] + \[Xi])	 
];
 
(*See eqs. 15-16*)
BounceParameterRvb[initialR_?NumericQ,\[Phi]L_,a_,d_,Ns_,backward_,pos_]:=
Module[{R,b,v,\[Alpha],x,y,z,Rvb},
	(*-------Backward--------------------*)
	If[backward,
		\[Alpha] = Join[a,{0.}];
		R = RInitial[d,initialR,\[Alpha],\[Phi]L,pos,backward];
		b = 0.; 
		v = \[Phi]L[[-1]];
		Rvb = Reap[
			Sow[b,x];Sow[v,y];Sow[R,z];
			Do[ b +=-(4/d)(\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^(d);Sow[b,x];
				v +=(4/(d-2))(\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^2 ;Sow[v,y];
				R  = Rs[d,(\[Phi]L[[-i-1]]-v) ,\[Alpha][[-i-1]],b];Sow[R,z];
				,{i,1,Ns}
			];
		][[2]];
		
		Return[Reverse[Rvb,{1,2}]],     
	(*--------Else-Forward---------------*)
		\[Alpha] = Join[{0.}, a ]; 
		R = RInitial[d,initialR,\[Alpha],\[Phi]L,pos,backward]; 
		b = 0.;
		v = \[Phi]L[[pos]] -(4/d)R^2\[Alpha][[pos]];
		Rvb = Reap[
			If[pos>1,Do[Sow[0,x];Sow[v,y];Sow[b,z];,{i,1,pos-1}]];
			Sow[R,x];
			Do[ v+= -( 4/(d-2)) (\[Alpha][[i+1]]-\[Alpha][[i]]) R^2 ;Sow[v,y];
				b+= (4/d) (\[Alpha][[i+1]]-\[Alpha][[i]]) R^(d);Sow[b,z];
				R = Rs[d,(\[Phi]L[[i+1]]-v) ,\[Alpha][[i+1]],b];Sow[R,x];
			,{i,pos,Ns}];
			v+= -( 4/(d-2))(-\[Alpha][[Ns+1]]) R^2 ;Sow[v,y];
			b+= (4/d) (-\[Alpha][[Ns+1]]) R^(d);Sow[b,z]; 
		][[2]];
		
		Chop@Return[Rvb]        
	]; 
];


(* ::Subsubsection::Closed:: *)
(*FindInitialRadius*)


FindBounce::errbc = "Large error in solving the boundary conditions. This potential may not have a bounce solution.";

(*Find the solution of eq. 25 or 26.*)
FindInitialRadius[d_,VL_,\[Phi]L_,a_,Ns_,maxIteR_,actionTolerance_,ansatzInitialR_,aRinitial_,pos_,switchMessage_]:= 
Module[{radii,initialR,ite,complexR,realR,switch,k,initialR0,complexityR,rw,findInitialR,
	\[Lambda],lambda,Rvb,kinetic,potential,errorAction,action,actionPrev = 0,actionToleranceB,minimum\[Lambda] =.01},
(*The initial Radius can be also found simply with FindRoot[R[x],{x,x0}]. However, this approach fails 
if it hits a singularity in R[x] while looking for x0. Thereby, in order to have a robust mechanism, 
FindInitialRadius checks the output of FindRoot after each iterations and uses the bisection method 
whenever FindRoot fails.*)
	
	(*Defines \[Lambda] of eq. 26 and other internal functions related to the validity of the solution*)
	lambda[initialR_?NumericQ] := lambda[initialR] = (
		Rvb = BounceParameterRvb[initialR,\[Phi]L,a,d,Ns,True,pos];
		kinetic = \[ScriptCapitalT][Sequence@@Rvb,a,d,Ns,pos];
		potential = \[ScriptCapitalV][Sequence@@Rvb,a,d,VL,\[Phi]L,Ns];
		(*The action and actionTolerance*)
		action = kinetic + potential;
		errorAction = Abs[(action-actionPrev)/action];
		actionToleranceB = If[Not[action == actionPrev], errorAction > actionTolerance,True];
		(*Discriminates between undershooting and overshooting*)
		complexityR = Abs@Im@Chop@Rvb[[1,1]];
		
		\[Lambda] = Sqrt[(2-d)*kinetic/(d*potential)]
	);
	
	(*Looks for the initial radius*)
	findInitialR[initialR_]:= 
		Quiet[Abs[rw/.FindRoot[Abs[Re@lambda[rw]-1],{rw,initialR}, 
					MaxIterations->1,
					PrecisionGoal->0,
					AccuracyGoal->0]
				],
		{FindRoot::lstol,FindRoot::cvmit}];   
	
	If[Ns==2&&Not[d==3],
		(*Exact radius from the close form solution N=3 field points and D=4 dimensions, see appendix B*)
		initialR = ansatzInitialR;
		lambda[initialR]
		,	
		(*Picks up the best estimate*)
		If[NumericQ[aRinitial],
			initialR = aRinitial
			,
			radii = BounceParameterRvb[0.,\[Phi]L,a,d,Ns,False,1][[1]];
			If[ Abs[Im[radii[[-1]]]] < 0,
				initialR = Abs[radii[[-2]]],
				initialR = Abs[ansatzInitialR] 
			];   
		];
		initialR0 = initialR;
	
		(*Finds the interval of the solution Sol \[Element] [realR, complexR] or reduces the interval*)
		ite = 0; 
		switch = True; 
		k = 1;
		While[ ite <= maxIteR && switch,   
			complexR = Infinity;
			realR = 0;
			switch = False;
			lambda[initialR]; (*lambda[initialR] also re-evaluates the values of Rvb,kinetic,potential and complexityR.*)
			If[ complexityR > 0,
				(*Overshooting*)
				While[ complexityR > 0 &&ite <= maxIteR&&(actionToleranceB||Abs[\[Lambda]-1]>minimum\[Lambda])&&Chop[Re@\[Lambda]]!=0,
					If[initialR < complexR,
						complexR = initialR;
						initialR = findInitialR[initialR];
						,
						(*Perturbs 10% down the ansatz since initialR>complexR*) 
						initialR = 0.9*complexR
					];
					actionPrev = action;lambda[initialR];
					ite++
				]; 
				realR = initialR;
				,(*Undershooting*)
				While[ complexityR == 0&&ite <= maxIteR&&(actionToleranceB||Abs[\[Lambda]-1]>minimum\[Lambda])&&Chop[Re@\[Lambda]]!=0,
					If[initialR > realR,
						realR = initialR;  
						initialR = findInitialR[initialR],
						(*Perturbs 10% up the ansatz since initialR<realR*)  
						initialR = 1.1*realR
					];
					actionPrev = action;lambda[initialR];
					ite++
				]; 
				complexR = initialR;    
			];  
			(*One the interval is found, reduces the interval and use bisection method*)	
			While[ite <= maxIteR&&(actionToleranceB||Abs[\[Lambda]-1]>minimum\[Lambda])&&Chop[Re@\[Lambda]]!=0, 
				
				If[ complexityR > 0,
					If[ (*Overshooting*)
						initialR < complexR,
						complexR = initialR;
						initialR = findInitialR[initialR];
						,
						initialR = Abs[complexR+realR]/2;
					];
					,
					If[ (*Undershooting*)
						initialR > realR, 
						realR = initialR; 
						initialR = findInitialR[initialR];
						,
						initialR = Abs[complexR+realR]/2;
					];
				];
				actionPrev = action;lambda[initialR];
				ite++   
			];
	   
			If[ Chop[Re@\[Lambda]] == 0,
				k++;
				switch=True;
				initialR = (1+k)*Abs[initialR0]  
			];   
		];    

		If[ ite > maxIteR&&switchMessage&&errorAction>actionTolerance*1.1,
			Message[FindBounce::cvmit,maxIteR] 
		];

		If[ Abs[\[Lambda]-1]>minimum\[Lambda]&&switchMessage, 
			Message[FindBounce::errbc];
			Return[$Failed,Module]
		];	
	];	
	Clear[lambda];
			
	Re@{initialR,Sequence@@Rvb,potential,kinetic}
]; 


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalT]*)


(*See eqs. 10*)
\[ScriptCapitalT][R_,v_,b_,a_,d_,Ns_,pos_]:= 
Module[{T},
	If[pos>1,
		T = 2\[Pi]^(d/2)/Gamma[d/2](
			32 a[[pos-1]]^2/(d^2(d+2)) (R[[pos]]^(2+d)) -
			8 a[[pos-1]]*b[[pos-1]]/d (R[[pos]]^2)-
			(2/(d-2))*b[[pos-1]]^2*(R[[pos]]^(2-d))), 
		T = 0
	];	
	T += 2\[Pi]^(d/2)/Gamma[d/2]Sum[
		 32*a[[i]]^2/(d^2(d+2)) (R[[i+1]]^(2+d) - R[[i]]^(2+d))  -
		 8*a[[i]]*b[[i]]/d  (R[[i+1]]^2 - R[[i]]^2) -
		(2/(d-2))*b[[i]]^2 (R[[i+1]]^(2-d)-R[[i]]^(2-d))
	,{i,pos,Ns}]
];


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalV]*)


 (*See eq. 11*)
\[ScriptCapitalV][R_,v_,b_,a_,d_,VL_,\[Phi]L_,Ns_]:= 
	2\[Pi]^(d/2)/Gamma[d/2](Sum[   
	32*a[[i]]^2/(d(d+2))  (R[[i+1]]^(2+d) - R[[i]]^(2+d)) + 
	8*a[[i]]*b[[i]]/(d-2) (R[[i+1]]^2 -R[[i]]^2) +
	( VL[[i]]-VL[[-1]]+ 8*a[[i]]( v[[i]] - \[Phi]L[[i]]))(R[[i+1]]^d-R[[i]]^d)/d ,{i,1,Ns}] +
		1/d R[[1]]^d (VL[[1]] - VL[[-1]]) );


(* ::Subsubsection::Closed:: *)
(*SingleFieldBounce*)


FindBounce::extrema = "Wrong position of the extrema, check the minima or use \"MidFieldPoint\" to include the maximum/saddle point of the potential.";
FindBounce::pathdef = "The path is deformed irregularly on the potential. Verify that the vacuum is a minimum of the potential (not a saddle point) or changes the number of segements.";
FindBounce::nosol = "Solution not found, increase the number of segments or accuracy.";

SingleFieldBounce[V_,potentialPoints_,Ns_,noFields_,\[Phi]L_,dim_,maxIteR_,actionTolerance_,
	ansatzInitialR_,aRinitial_,rule_,iter_,switchMessage_]:= 
Module[{a,VL,pos,initialR,R,v,b,T1,V1},
	
	If[
		potentialPoints===None,
		VL = Table[V/.rule[[s]],{s,Ns+1}]
		,
		If[potentialPoints[[1]]<potentialPoints[[-1]],
			VL = potentialPoints,
			VL = Reverse[potentialPoints]
		];
	];
	
	If[VL[[1]]>=VL[[2]]||VL[[-1]]>=VL[[-2]],
		If[iter === 0,
			Message[FindBounce::extrema];
			Return[$Failed,Module]
			,
			Message[FindBounce::pathdef];
			Return[$Failed,Module]
		]
	];
	
	a  = Table[ ((VL[[s+1]]-VL[[s]])/(\[Phi]L[[s+1]]-\[Phi]L[[s]]))/8 ,{s,Ns} ]; 
	pos = FindSegment[a,\[Phi]L,dim,Ns];
	{initialR,R,v,b,V1,T1} = FindInitialRadius[dim,VL,\[Phi]L,a,Ns,maxIteR,actionTolerance,ansatzInitialR,
		aRinitial,pos,switchMessage]/.x_/;FailureQ[x]:>Return[$Failed,Module];
		
	(*Checks if we got a consistent answer.*)
	If[V1+T1<0&&switchMessage,
		Message[FindBounce::nosol];
		Return[$Failed,Module]
	];
	
	If[pos>1, R[[pos-1]]=0 ];
	a = Join[a,{0}];	

	{V1+T1,VL,v,a,b,pos,R,initialR}
];


(* ::Subsection::Closed:: *)
(*SingleFieldBounceExtension*)


(* ::Subsubsection::Closed:: *)
(*Find\[ScriptCapitalI]*)


(*See eqs. 47-48*)
Find\[ScriptCapitalI][v_,\[Phi]L_,a_,b_,Ns_,pos_,R_,ddVL_,d_]:=
Module[{\[ScriptCapitalI],d\[ScriptCapitalI],v0,\[Phi]L0,a0,b0,ddVL0},
	v0 =Join[{0},v,{0}]; 
	\[Phi]L0 =Join[{0},\[Phi]L,{0}]; 
	a0 =Join[{0},a,{0}];
	b0 =Join[{0},b,{0}]; 
	ddVL0 =Join[{0},ddVL,{0}];
	If[d==4,
		\[ScriptCapitalI] = Table[Join[ConstantArray[0,pos-1],
			Table[   ddVL0[[s+m]] (   ( v0[[s+m]]-\[Phi]L0[[s+m]] )/8 *R[[s]]^2 +a0[[s+m]]/24*R[[s]]^4   + b0[[s+m]]/2 Log[R[[s]]] )  ,{s,pos,Ns+1}]  ],{m,0,1}];
		d\[ScriptCapitalI] = Table[Join[ConstantArray[0,pos-1], 
			Table[   ddVL0[[s+m]] (   ( v0[[s+m]]-\[Phi]L0[[s+m]] )/4 *R[[s]] +a0[[s+m]]/6*R[[s]]^3 + b0[[s+m]]/2 /R[[s]]   )   ,{s,pos,Ns+1}]  ],{m,0,1}];          
	];
	If[d==3,
		\[ScriptCapitalI] = Table[Join[ConstantArray[0,pos-1],
			Table[   ddVL0[[s+m]] (   ( v0[[s+m]]-\[Phi]L0[[s+m]] )/6*R[[s]]^2 +
			a0[[s+m]]/15*R[[s]]^4 + b0[[s+m]]*R[[s]])   ,{s,pos,Ns+1}]  ],{m,0,1}];
		d\[ScriptCapitalI] = Table[Join[ConstantArray[0,pos-1], 
			Table[   ddVL0[[s+m]] (   ( v0[[s+m]]-\[Phi]L0[[s+m]] )/3*R[[s]] +
			a0[[s+m]] 4/15 R[[s]]^3 + b0[[s+m]] )   ,{s,pos,Ns+1}]  ],{m,0,1}];
	];
	
	{\[ScriptCapitalI],d\[ScriptCapitalI]}   
];


(* ::Subsubsection::Closed:: *)
(*BounceParameterr\[Beta]\[Nu]*)


(*See eqs. 50-52*)
BounceParameterr\[Beta]\[Nu][rw_?NumericQ,a_,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_] := 
Module[{r\[Beta]\[Nu]M,r,\[Beta],\[Nu],x,y,z,\[Beta]prev,\[Alpha]0,a0,c0,b0},
	a0 = Join[{0},a]; 
	\[Alpha]0 = Join[{0},\[Alpha]];
	b0 = Join[{0},b];
	c0 = Table[2 (b0[[i]] - 4/d a0[[i]] R[[i]]^d),{i,1,Ns+1}];
	r\[Beta]\[Nu]M=Reap[  
		r = 0;
		Sow[r,x]; 
		\[Beta] = \[Beta]prev = 0;
		Sow[\[Beta],y];
		\[Nu] = (rw c0[[pos]] )R[[pos]]^(2-d)-4/d \[Alpha]0[[pos]] R[[pos]]^(2)-\[ScriptCapitalI][[1,pos]]; 
		Sow[\[Nu],z];
		r = rw; 
		Sow[r,x]; 
		\[Beta] += (4/d (\[Alpha]0[[pos+1]]-\[Alpha]0[[pos]])+4*r(a0[[pos+1]]-a0[[pos]]) )   R[[pos]]^d + 
			(d\[ScriptCapitalI][[2,pos]] -  d\[ScriptCapitalI][[1,pos]]  )R[[pos]]^(d-1) /2; 
		Sow[\[Beta],y];
		\[Nu] += -2/(d-2)(\[Beta]-\[Beta]prev)R[[pos]]^(2-d)-4/d (\[Alpha]0[[pos+1]]-\[Alpha]0[[pos]])R[[pos]]^2 - 
			(\[ScriptCapitalI][[2,pos]]-\[ScriptCapitalI][[1,pos]]); 
		Sow[\[Nu],z];
		Do[ \[Beta]prev=\[Beta];   
			r =(2/(d-2)\[Beta]+(\[Nu]+\[ScriptCapitalI][[1,i]]+4/d \[Alpha]0[[i]] R[[i]]^2)R[[i]]^(d-2))/c0[[i]]; 
			Sow[r,x]; 
			\[Beta]+=4/d (\[Alpha]0[[i+1]]-\[Alpha]0[[i]])R[[i]]^d+4r(a0[[i+1]]-a0[[i]])R[[i]]^d + 
				(  d\[ScriptCapitalI][[2,i]] -  d\[ScriptCapitalI][[1,i]]  )R[[i]]^(d-1) /2; 
			Sow[\[Beta],y];
			\[Nu]+=-2/(d-2)(\[Beta]-\[Beta]prev)/R[[i]]^(d-2)-4/d (\[Alpha]0[[i+1]]-\[Alpha]0[[i]])R[[i]]^2 -
				 (\[ScriptCapitalI][[2,i]]-\[ScriptCapitalI][[1,i]]); 
			Sow[\[Nu],z];
			,{i,pos+1,Ns}
		];
		r = ( 2 \[Beta] / R[[Ns+1]]^(d-1)  - 8/d \[Alpha]0[[Ns+1]] R[[Ns+1]] -
			d\[ScriptCapitalI][[1,Ns+1]])/(8 a0[[Ns+1]] R[[Ns+1]]);
		Sow[r,x];   
		\[Beta] = 0;
		Sow[\[Beta],y];
		\[Nu] = 0; 
		Sow[\[Nu],z];               
		][[2]];
		
	r\[Beta]\[Nu]M 
];


(* ::Subsubsection::Closed:: *)
(*FindInitialRadiusExtension*)


(*Solves Boundary Conditions \[Xi][N-1] = 0*)
FindInitialRadiusExtension[rw_?NumericQ,a_ ,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_]:=
Module[{r,\[Beta],\[Nu]},
	{r,\[Beta],\[Nu]} = BounceParameterr\[Beta]\[Nu][rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos] ;
		
	\[Nu][[-2]] + 2/(d-2) \[Beta][[-2]]*R[[-1]]^(2-d)+4/d \[Alpha][[Ns]]*R[[-1]]^2+\[ScriptCapitalI][[1,-1]]   
];


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalT]\[Xi]*)


(*Kinetic term from Int[\[Rho]^(D-1)(1/2 d\[Xi]^2)]*)
\[ScriptCapitalT]\[Xi]4[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=
	1/(24 \[Rho]^2) (a \[Rho]^4 (-48 \[Beta]+\[Rho]^4 (2 ddVL v+16 \[Alpha]+a ddVL \[Rho]^2-2 ddVL \[Phi]L))+
	b (-48 \[Beta]+2 \[Rho]^4 (-3 ddVL v-24 \[Alpha]+2 a ddVL \[Rho]^2+3 ddVL \[Phi]L))
	-24 b^2 ddVL \[Rho]^2 Log[\[Rho]])+ (r (4 b-4 a \[Rho]^4)^2)/(8 \[Rho]^2);

\[ScriptCapitalT]\[Xi]3[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=
	-((4 b \[Beta])/\[Rho])-2 b^2 ddVL \[Rho]+8/15 a b ddVL \[Rho]^4+32/315 a^2 ddVL \[Rho]^7-1/3 \[Rho]^2 (8 a \[Beta]+b (8 \[Alpha]+ddVL (v-\[Phi]L)))+
	8/45 a \[Rho]^5 (8 \[Alpha]+ddVL (v-\[Phi]L))+(2 r (3 b-4 a \[Rho]^3)^2)/(9 \[Rho]);

\[ScriptCapitalT]\[Xi][d_?NumericQ,a_,R_,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,ddVL_,VL_,\[Phi]L_,Ns_,pos_,r_] :=
Module[{\[ScriptCapitalT],VN,\[ScriptCapitalT]\[Xi]D,p},
	VN = VL[[-1]];
	If[d==4,\[ScriptCapitalT]\[Xi]D = \[ScriptCapitalT]\[Xi]4];
	If[d==3,\[ScriptCapitalT]\[Xi]D = \[ScriptCapitalT]\[Xi]3];	
	If[pos>1,
		p = pos-1;
		\[ScriptCapitalT]= \[ScriptCapitalT]\[Xi]D[a[[p]],R[[pos]],b[[p]],v[[p]],\[Alpha][[p]],\[Beta][[p]],\[Nu][[p]],VL[[p]],VN,ddVL[[p]],\[Phi]L[[p]],r[[p]]],
		\[ScriptCapitalT]= \[ScriptCapitalT]\[Xi]D[0,R[[pos]],0,0,0,0,0,VL[[pos]],VN,ddVL[[pos]],0,r[[pos]]]
	];
	
	\[ScriptCapitalT] += 2\[Pi]^(d/2)/Gamma[d/2] Sum[\[ScriptCapitalT]\[Xi]D[a[[i]],R[[i+1]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i+1]]]-
			\[ScriptCapitalT]\[Xi]D[a[[i]],R[[i]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i]]]   
	,{i,pos,Ns}]
];


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalV]\[Xi]*)


(*Potential term from Int[\[Rho]^(D-1)(Vperturbation[\[CurlyPhi]]-8a)]*)
\[ScriptCapitalV]\[Xi]4[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=
	1/48 \[Rho]^2 (5 a^2 ddVL \[Rho]^6+ 16 a (12 \[Beta]+6 \[Nu] \[Rho]^2+\[Rho]^4 (8 \[Alpha]+ddVL (v-\[Phi]L)))+
	24 b (8 \[Alpha]+ddVL (v-\[Phi]L))+6 \[Rho]^2 (16 \[Alpha]+ddVL (v-\[Phi]L)) (v-\[Phi]L))+
	1/2 b ddVL (b+2 a \[Rho]^4) Log[\[Rho]]+8 a b r \[Rho]^2+1/4 r \[Rho]^4 (4 (VL-VN)+32 a^2 \[Rho]^2+32 a (v-\[Phi]L));

\[ScriptCapitalV]\[Xi]3[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=
	1/630 \[Rho] (1260 b^2 ddVL+210 b \[Rho] (24 \[Alpha]+ddVL (3 v+8 a \[Rho]^2-3 \[Phi]L))+
	\[Rho] (128 a^2 ddVL \[Rho]^5+336 a (15 \[Beta]+5 \[Nu] \[Rho]+\[Rho]^3 (8 \[Alpha]+ddVL (v-\[Phi]L)))+105 \[Rho] (16 \[Alpha]+
	ddVL (v-\[Phi]L)) (v-\[Phi]L)))+16 a b r \[Rho]^2+1/3 r \[Rho]^3 (3 (VL-VN)+32 a^2 \[Rho]^2+24 a (v-\[Phi]L));

\[ScriptCapitalV]\[Xi][d_?NumericQ,a_,R_,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,ddVL_,VL_,\[Phi]L_,Ns_,pos_,r_] := 
Module[{\[ScriptCapitalV],\[ScriptCapitalV]\[Xi]D,p},
	If[d==4, \[ScriptCapitalV]\[Xi]D = \[ScriptCapitalV]\[Xi]4];
	If[d==3, \[ScriptCapitalV]\[Xi]D = \[ScriptCapitalV]\[Xi]3];
	If[pos>1,
		p = pos-1;
		\[ScriptCapitalV]= 2\[Pi]^(d/2)/Gamma[d/2]\[ScriptCapitalV]\[Xi]D[a[[p]],R[[pos]],b[[p]],v[[p]],\[Alpha][[p]],\[Beta][[p]],\[Nu][[p]],VL[[p]],VL[[-1]],ddVL[[p]],\[Phi]L[[p]],r[[pos]]];,
		\[ScriptCapitalV]= 2\[Pi]^(d/2)/Gamma[d/2]\[ScriptCapitalV]\[Xi]D[0,R[[pos]],0,0,0,0,0,VL[[1]],VL[[-1]],ddVL[[1]],0,r[[pos]]];
	];
	
	\[ScriptCapitalV] += 2\[Pi]^(d/2)/Gamma[d/2]Sum[\[ScriptCapitalV]\[Xi]D[a[[i]],R[[i+1]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VL[[-1]],ddVL[[i]],\[Phi]L[[i]],r[[i+1]]]-
		\[ScriptCapitalV]\[Xi]D[a[[i]],R[[i]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VL[[-1]],ddVL[[i]],\[Phi]L[[i]],r[[i]]]   
	,{i,pos,Ns}]
];


(* ::Subsubsection::Closed:: *)
(*SingleFieldBounceExtension*)


FindBounce::gradfail = "The gradient of the potential is not well defined at some field value. \"Gradient\"->None was taken.";

SingleFieldBounceExtension[VL_,Ns_,v_,a_,b_,R_,\[Phi]L_,pos_,dim_,eL_,gradient_]:=
Module[{dVL,\[Alpha],\[ScriptCapitalI],d\[ScriptCapitalI],r1,rInitial,r,\[Beta],\[Nu],eL0,ddVL,extensionPB=True,T\[Xi]=0.,V\[Xi]=0.},

	dVL=MapThread[Dot[#1,#2]&,{gradient,Join[eL,{Last@eL}]}];

	If[And@@(NumericQ[#]&/@dVL),
		\[Alpha] = Join[a[[1;;Ns]] - dVL[[2;;Ns+1]]/8 ,{0}];
		ddVL = Table[ (dVL[[s+1]]-8(a[[s]]+\[Alpha][[s]]))/(\[Phi]L[[s+1]]-\[Phi]L[[s]]),{s,Ns}];
		{\[ScriptCapitalI],d\[ScriptCapitalI]} = Find\[ScriptCapitalI][v,\[Phi]L,a,b,Ns,pos,R,ddVL,dim];
		r1 = rInitial/.FindRoot[FindInitialRadiusExtension[rInitial,a,b,dim,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos],{rInitial,-1}]//Quiet;
		{r,\[Beta],\[Nu]} = If[
			pos>1,
			Join[ConstantArray[0,{pos-2,3}],
			BounceParameterr\[Beta]\[Nu][r1,a,b,dim,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos]//Transpose]//Transpose,
			BounceParameterr\[Beta]\[Nu][r1,a,b,dim,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos][[All,2;;-1]]
		];
		T\[Xi] = \[ScriptCapitalT]\[Xi][dim,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r];
		V\[Xi] = \[ScriptCapitalV]\[Xi][dim,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r];
		,
		Message[FindBounce::gradfail];
		extensionPB=False
	];
	
	{Re[V\[Xi]+T\[Xi]],\[Alpha],\[Beta],\[Nu],ddVL,extensionPB}
];


(* ::Subsection::Closed:: *)
(*ParameterInFieldSpace*)


ParameterInFieldSpace[vs_,as_,bs_,\[Alpha]s_,\[Beta]s_,\[Nu]s_,\[Phi]_,eL_,l_,\[Phi]L_,Ns_,noFields_,pos_,dim_,bottomless_,actionOld_,actionNew_,actionTolerance_,switchPathOld_,extensionPB_]:=
Module[{v,a,b,\[Alpha],\[Beta],\[Nu],path=\[Phi],switchPath = switchPathOld},
	v = Table[\[Phi][[s+1]]+eL[[s]]*(vs[[s]]-(l[[s]]+\[Phi]L[[s]])),{s,Ns}];
	a = Table[eL[[s]]*as[[s]],{s,Ns}];
	b = Table[eL[[s]]*bs[[s]],{s,Ns}];
	
	If[extensionPB,
		\[Alpha] = Table[eL[[s]]*\[Alpha]s[[s]],{s,Ns}];
		\[Beta] = Table[eL[[s]]*\[Beta]s[[s]],{s,Ns}];
		\[Nu] = Table[eL[[s]]*\[Nu]s[[s]],{s,Ns}];
	];
	
	If[bottomless,
		path[[1]] = v[[1]]+b[[1]];
		a[[1]] = Abs@a[[1]]
	];
	
	If[Abs[(actionOld-actionNew)/actionNew] < actionTolerance,
		switchPath = True
	];
	
	{v,a,b,\[Alpha],\[Beta],\[Nu],path,actionNew,switchPath}
];


(* ::Subsection::Closed:: *)
(*MultiFieldBounce*)


MultiFieldBounce[fieldPoints_,gradient_,hessian_,noFields_,pos_,d_,R0_,\[ScriptV]0_,\[ScriptA]0_,\[ScriptB]0_,lengthPath_,pathTolerance_]:=Module[
	{Ns,\[Nu],\[Beta],rI,a,R,\[Zeta]ts,\[Nu]\[Beta],x,y,d\[CurlyPhi],rF,DV,D2V,
	\[Xi]Mc,M,c,\[Nu]0,\[Beta]0,\[Nu]\[Xi]p,\[Nu]\[Xi]m,\[Beta]\[Xi]p,\[Beta]\[Xi]m,fLowT,fD,fD1,frI,n1,p,
	\[Phi]0 = fieldPoints, v0 = \[ScriptV]0, a0 = \[ScriptA]0, b0 = \[ScriptB]0,switchPath = False,n,m},

	Ns=Length[fieldPoints]-1;	
	If[pos>1,
		p = pos-1;
		\[Phi]0[[p]] = v0[[p]]
		,
		p = pos
	];

	R = Transpose@Table[R0,{i,noFields}];

	DV = gradient;
	D2V = hessian;
	d\[CurlyPhi] = Threshold@Join[
		Table[ConstantArray[0.,{2,noFields}],{s,1,p-1}],
		Table[8./d a0[[s+m]]*R[[s+1]]- 2 b0[[s+m]]/R[[s+1]]^(d-1),{s,p,Ns-1},{m,0,1}]
	];

	(*c*)
	If[pos>1, 
		\[Nu]0[p] = ConstantArray[0,noFields]; 
		\[Beta]0[p] = ConstantArray[0,noFields];
		, 
		\[Nu]0[p] = ((16 a0[[1]]-DV[[1]]-DV[[2]]) R[[1]]^2)/(4 (-2+d));  
		\[Beta]0[p] = ((-16 a0[[1]]+DV[[1]]+DV[[2]]) R[[1]]^d)/(4 d);  
	];
	
	Do[ \[Nu]0[s] = \[Nu]0[s-1]+ 1/(4 (-2+d)) R[[s]] (4 d\[CurlyPhi][[-1+s,1]]-4 d\[CurlyPhi][[-1+s,2]]+
			(-16 a0[[-1+s]]+16 a0[[s]]+DV[[-1+s]]-DV[[1+s]]) R[[s]]);
		\[Beta]0[s] = \[Beta]0[s-1]+ 1/(4 d) R[[s]]^(-1+d) (-2 d d\[CurlyPhi][[-1+s,1]]+2 d d\[CurlyPhi][[-1+s,2]]+
			(16 a0[[-1+s]]-16 a0[[s]]-DV[[-1+s]]+DV[[1+s]]) R[[s]]);
	,{s,p+1,Ns}];
	
	c  = -Flatten[ {If[pos>1,{},ConstantArray[0,noFields]],
		Table[ ((-16 a0[[s]]+DV[[s]]+DV[[1+s]]) R[[1+s]]^2)/(4 d)+
			(2 (R[[1+s]]^(2-d)) )/(-2+d) \[Beta]0[s]+\[Nu]0[s],{s,p,Ns}]
		,ConstantArray[0,noFields]}    
	];
	
	(*M*)
	\[Nu]\[Xi]p[s_]:= \[Nu]\[Xi]p[s]= -( (R[[s]]^2) /(4 (-2+d)))D2V[[1+s]];(*\[Zeta]ts[1+s]*)
	\[Beta]\[Xi]p[s_]:= \[Beta]\[Xi]p[s]= (D2V[[1+s]] (R[[s]]^d) )/(4 d);(*\[Zeta]ts[1+s]*)
	\[Nu]\[Xi]m[s_]:= \[Nu]\[Xi]m[s]= (D2V[[s-1]] (R[[s]]^2) )/(4 (-2+d));(*\[Zeta]ts[-1+s]*)
	\[Beta]\[Xi]m[s_]:= \[Beta]\[Xi]m[s]= -((D2V[[s-1]] (R[[s]]^d) )/(4 d));(*\[Zeta]ts[-1+s]*)
	fLowT[s_,j_]:= (2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]m[j]+\[Beta]\[Xi]p[j-2])+\[Nu]\[Xi]m[j]+\[Nu]\[Xi]p[j-2];(*[Eq[s],\[Xi][j-1]]*) (* 2 \[LessEqual] j \[LessEqual] s*)

	\[Nu]\[Xi]p[p-1]= If[pos>1,IdentityMatrix[noFields],-((D2V[[p]]*(R[[p]]^2) )/(4 (-2+d))) ];
	\[Beta]\[Xi]p[p-1]= If[pos>1,ConstantArray[0,{noFields,noFields}],(D2V[[p]]*R[[p]]^d)/(4 d)];
	fD[s_]:= (D2V[[s]]*(R[[1+s]]^2) )/(4 d) +(2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]p[s-1])+\[Nu]\[Xi]p[s-1];(*[Eq[s],\[Xi][s]]*)

	\[Nu]\[Xi]p[p] = If[pos>1,ConstantArray[0,{noFields,noFields}],-((D2V[[1+p]]*R[[p]]^2)/(4 (-2+d)))];
	\[Beta]\[Xi]p[p] = If[pos>1,ConstantArray[0,{noFields,noFields}],(D2V[[1+p]] (R[[p]]^d) )/(4 d)];
	fD1[s_]:= (-IdentityMatrix[noFields]+(D2V[[1+s]]*R[[1+s]]^2)/(4 d))+(2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]p[s])+\[Nu]\[Xi]p[s] ;(*[Eq[s],\[Xi][s+1]]*)
	frI[s_]:= ((2 (R[[1+s]]^(2-d)) )/(-2+d) (4 a0[[p]]*R[[p]]^d)-(8  a0[[p]]*R[[p]]^2)/(-2+d)  )IdentityMatrix[noFields]; 
	
	n1 = If[pos>1,0,1];
	
	(*M*)
	n = (Ns+1-p)+1+n1;
	m = SparseArray[{{n,n}->0}];
	Do[
		m[[s+n1,s+n1]]=fD[s+p-1]
		,{s,1,n-1-n1}
	];
	Do[
		m[[s+n1,s+1+n1]]=fD1[s+p-1]
		,{s,1,n-1-n1}
	];
	m[[n,n]] = IdentityMatrix[noFields];
	If[pos==1,
		m[[1,2]] = IdentityMatrix[noFields];
		Do[
			m[[s+n1,1]]=frI[s]
			,{s,1,n-1-n1}
		];   
	];
	Do[
		m[[s+n1,j-1+n1]]=fLowT[s+p-1,j+p-1]
		,{s,2,n-n1-1}
		,{j,2,s}
	];
	M = ArrayFlatten[m];
	
	(*Solving M.\[Xi] = c*)
	\[Xi]Mc = LinearSolve[M,c];  
	
	(*\[Zeta]ts,a,r,R*)
	\[Zeta]ts = Join[ConstantArray[0,{p-1,noFields}],Partition[\[Xi]Mc[[ 1+n1 noFields;;-1]],noFields]];
	a  = Join[ConstantArray[0,{p-1,noFields}],Table[ 1/8( (  DV[[s]]+DV[[s+1]]+
		D2V[[s]].\[Zeta]ts[[s]]+D2V[[s+1]].\[Zeta]ts[[s+1]]  )/2)-a0[[s]],{s,p,Ns}] ];
	rI  = If[pos>1,
			ConstantArray[0,noFields], 
			\[Xi]Mc[[1;;noFields]] 
			];
	rF  = ConstantArray[0,noFields];
	
	(*\[Nu]\[Beta]*)
	If[pos>1,
		\[Nu] = \[Zeta]ts[[p]]; 
		\[Beta] = ConstantArray[0.,noFields];
		,
		\[Nu] = -4/(d-2) (a[[1]]+2 a0[[1]] rI)R[[1]]^2;
		\[Beta] = 4/d (a[[1]]+d a0[[1]] rI)R[[1]]^d;    
	];
	     
	\[Nu]\[Beta]=Reap[
		If[ pos>1, 
			Do[
				Sow[ConstantArray[0,noFields],x];
				Sow[ConstantArray[0,noFields],y];
			,{s,1,p-1}];
		];
		Sow[\[Nu],x];
		Sow[\[Beta],y];
		Do[ \[Nu]+= -4/(d-2)(a[[s+1]]-a[[s]])R[[s+1]]^2-1/(d-2) (  d\[CurlyPhi][[s,2]]-d\[CurlyPhi][[s,1]] )R[[s+1]] ;
			\[Beta]+=  4/d(a[[s+1]]-a[[s]])R[[s+1]]^d+1/2 (  d\[CurlyPhi][[s,2]]-d\[CurlyPhi][[s,1]] )R[[s+1]]^(d-1);
			Sow[\[Nu],x];
			Sow[\[Beta],y]; 
		,{s,p,Ns-1} ];  
	][[2]];
	 
	{\[Phi]0[[1]],\[Zeta]ts[[1]]} = {fieldPoints[[1]],ConstantArray[0,noFields]};

	If[Max[Abs@\[Zeta]ts]/lengthPath < pathTolerance,
		switchPath = True
	];
	Clear[\[Nu]\[Xi]p,\[Beta]\[Xi]p,\[Nu]\[Xi]m,\[Beta]\[Xi]m,\[Nu]0,\[Beta]0];
	
	{\[Phi]0+\[Zeta]ts,v0+\[Nu]\[Beta][[1]],a0+a,b0+\[Nu]\[Beta][[2]],R0,pos,switchPath}   	
];


(* ::Subsection::Closed:: *)
(*BottomlessPotential*)


(* ::Subsubsection::Closed:: *)
(*BottomlessPotential*)


BottomlessPotential[initialR_?NumericQ,a_,\[Phi]L_,\[Phi]m_] :=
BottomlessPotential[initialR,a,\[Phi]L,\[Phi]m] = 
Module[{\[Phi]0,b1,v1,bB,vB,p = 2},
	vB = \[Phi]m + \[Phi]L[[p]];
	\[Phi]0 = \[Phi]m - (1 + Sqrt[1 - 2 a[[p-1]] initialR^2 (vB-\[Phi]L[[p]])^2])/(a[[p-1]]*initialR^2 (vB-\[Phi]L[[p]]));
	bB = (\[Phi]0 - \[Phi]m);
	v1 = vB - 2 a[[p]] initialR^2+4 bB/(2+bB^2 a[[p-1]] initialR^2)^2;
	b1 = initialR^4 (2*bB^3 a[[p-1]]+a[[p]] (2+bB^2*a[[p-1]]*initialR^2)^2)/(2+bB^2*a[[p-1]]*initialR^2)^2;		
	
	{v1,b1,bB,vB}
];


(* ::Subsubsection::Closed:: *)
(*BottomlessParameterRvb*)


Rs[4,c1_?NumericQ,a_,b_] := Sqrt[ 1/2 (Sqrt[ c1^2 - 4 a b ] + c1)/a  ];

BottomlessParameterRvb[initialR_?NumericQ,\[Phi]L_,a_,d_,Ns_,backward_,\[Phi]m_]:=
BottomlessParameterRvb[initialR,\[Phi]L,a,d,Ns,backward,\[Phi]m] =
Module[{R,b,v,\[Alpha],x,y,z,Rvb,p=2,b4,v4},
    (*-------Backward--------------------*)
	If[backward,
		\[Alpha] = Join[a,{0.}]; 
		R = RInitial[d,initialR,\[Alpha],\[Phi]L,None,backward]; 
		b = 0.; 
		v=\[Phi]L[[-1]];
		Rvb = Reap[
			Sow[b,x];Sow[v,y];Sow[R,z];
			Do[ b +=-(4/d)(\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^(d);Sow[b,x];
				v +=( 4/(d-2))(\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^2 ;Sow[v,y]; 
				R  = Rs[d,(\[Phi]L[[-i-1]]-v) ,\[Alpha][[-i-1]],b];Sow[R,z];
				,{i,1,Ns-1}
			];
			{v,b,b4,v4} = BottomlessPotential[R,a,\[Phi]L,\[Phi]m];
			Sow[b4,x]; Sow[v4,y]; Sow[0,z];	 
		][[2]];
			
		Return[Reverse[Rvb,{1,2}]]
	,     
	\[Alpha] = Join[{0.},a];
	Rvb = Reap[
		{v,b,b4,v4} = BottomlessPotential[initialR,a,\[Phi]L,\[Phi]m];
		Sow[0,x]; Sow[v4,y]; Sow[b4,z];
		Sow[initialR,x]; Sow[v,y]; Sow[b,z];
		R = Rs[d, (\[Phi]L[[p+1]]-v) ,\[Alpha][[p+1]],b];Sow[R,x];
		
		Do[ v+= -( 4/(d-2)) (\[Alpha][[i+1]]-\[Alpha][[i]]) R^2 ; Sow[v,y];
			b+= (4/d) (\[Alpha][[i+1]]-\[Alpha][[i]]) R^(d); Sow[b,z];
			R = Rs[d, (\[Phi]L[[i+1]]-v) ,\[Alpha][[i+1]],b]; Sow[R,x];
		,{i,p+1,Ns}];
		v += -( 4/(d-2))(-\[Alpha][[Ns+1]]) R^2 ;Sow[v,y];
		b += (4/d) (-\[Alpha][[Ns+1]]) R^(d);Sow[b,z]; 
		][[2]];
		
		Chop@Return[Rvb]        
	];	         
];


(* ::Subsubsection::Closed:: *)
(*FindInitialRadiusB*)


(*Find the solution of eq. 25 or 26.*)
FindInitialRadiusB[d_,VL_,\[Phi]L_,a_,Ns_,maxIteR_,actionTolerance_,ansatzInitialR_,aRinitial_,\[Phi]m_]:= 
Module[{R,initialR,ite,findInitialR,\[Lambda],complexR,realR,switch,k,initialR0},
	(*Defines \[Lambda] of the bottomless potential*)
	\[Lambda][initialR_?NumericQ] := \[Lambda][initialR] = Chop[ Sqrt[\[CapitalLambda]B[initialR,d,VL,\[Phi]L,a,Ns,True,\[Phi]m]] ];
	(*Looks for the initial radius*)
	findInitialR[initialR_?NumericQ]:= findInitialR[initialR] = 
		Module[{rw}, 
			rw =rw/.Quiet[FindRoot[Abs[\[Lambda][rw]-1],{rw,initialR}, 
				MaxIterations->1,
				PrecisionGoal->0,
				AccuracyGoal->0],
				{FindRoot::lstol,FindRoot::cvmit}]; 
				  
			Abs[rw]
	]; 
	
	(*Picks up the best estimate*)
	If[NumericQ[aRinitial],
		initialR = aRinitial
		,
		R = BounceParameterRvb[0,\[Phi]L,a,d,Ns,False,1][[1]]//Chop;
		If[ Abs[ Im[R[[-1]]]]<10^(-12),
			initialR = Abs[R[[-2]]],
			initialR = Abs[ansatzInitialR] 
		];   
	];
	initialR0 = initialR;
	
	(*Finds the interval of the solution Sol \[Element] [realR, complexR] or reduces the interval*)
	ite = 0; switch = True; k = 1;
	While[ ite <= maxIteR && switch &&k<10,    
		complexR = Infinity; 
		realR=10^(-1); 
		switch = False;
		If[ Abs[Im[\[Lambda][initialR]]]>10^(-12),
			(*-Overshooting*)
			While[(Abs[Im[\[Lambda][initialR]]]>10^(-12)||(Abs@Im[\[Lambda][initialR]]<=  10^(-12)&&\[Lambda][initialR]>1))&&ite <= maxIteR&&Abs[\[Lambda][initialR]-1]>actionTolerance,
				If[initialR < complexR,
					complexR = initialR;
					initialR = findInitialR[initialR], 
					initialR = 0.9*complexR]; 
			ite++]; 
			,
			(*-Undershooting*)
			(*----else----------*)
			While[ Abs@\[Lambda][initialR]< 1 &&Abs@Im[\[Lambda][initialR]]<= 10^(-12)&&ite <= maxIteR&&Abs[\[Lambda][initialR]-1]>actionTolerance&&Chop@Re[\[Lambda][initialR]]!=0,
				If[initialR > realR,
					realR = initialR;
					initialR= findInitialR[initialR], 
					initialR = 1.1*realR;
					k++
				];
			ite++]; 
			complexR = initialR;    
		];   
		
		(*One the interval is found, reduces the interval and use bisection method*)	
		While[ite <= maxIteR &&Abs[\[Lambda][initialR]-1]>actionTolerance&&Chop[Re[\[Lambda][initialR]]]!=0, 
		
			If[Abs@Im[\[Lambda][initialR]]>10^(-12)||(Abs@Im[\[Lambda][initialR]]<=  10^(-12)&&\[Lambda][initialR]>1),
				If[ initialR <  complexR,
					complexR = initialR;
					initialR = findInitialR[initialR],
					(*-Overshooting*)
					initialR = Abs[complexR+realR]/2];,
				(*-----else-------------*)
				If[ initialR > realR, 
					realR = initialR; 
					initialR = findInitialR[initialR],
					(*-Undershooting*)
					initialR = Abs[complexR +realR]/2  
				]
			];
			ite++;     	
		];
		    
		If[ Re[Chop[\[Lambda][initialR]]] ==0, 
			k++;
			switch=True;
			initialR = (1+k)*Abs[initialR0] 
		];     
	];     
	
	If[ ite > maxIteR, 
		Message[FindInitialRadii::cvmit,maxIteR] 
	];
	
	 If[ Re[\[Lambda][initialR]-1] >0.5, 
		Message[FindBounce::errbc];
		Return[$Failed,Module]
	];  
	Clear[findInitialR,\[Lambda]];
	
	Re@initialR  
]; 


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalT]B,\[ScriptCapitalV]B,\[CapitalLambda]B*)


(*Kinetic term from Int[\[Rho]^(D-1)(1/2 d\[CurlyPhi]4^2)]*)
\[ScriptCapitalT]Bs[a4_,\[Rho]_?NumericQ,b4_] := 
	-((4 (4+3 b4^2 a4 \[Rho]^2 (2+b4^2 a4 \[Rho]^2)))/(3 a4 (2+b4^2 a4 \[Rho]^2)^3));

\[ScriptCapitalT]B[R_,Ns_,a_,b_] := 
Module[{\[ScriptCapitalT],p=2,d=4},
	\[ScriptCapitalT] = 2\[Pi]^(d/2)/Gamma[d/2] (
		\[ScriptCapitalT]Bs[a[[p-1]],R[[p]],b[[p-1]]] - 
		\[ScriptCapitalT]Bs[a[[p-1]],0.,b[[p-1]]]);
	
	\[ScriptCapitalT] += 2\[Pi]^(d/2)/Gamma[d/2] Sum[32 a[[i]]^2/(d^2(d+2)) (R[[i+1]]^(2+d)- 
		 R[[i]]^(2+d)) -8 a[[i]]*b[[i]]/d  (R[[i+1]]^2 -R[[i]]^2) -
		(2/(d-2))b[[i]]^2 (R[[i+1]]^(2-d)-R[[i]]^(2-d)),{i,p,Ns}] 
];

(*Potential term from Int[\[Rho]^3(VL[[2]] + \[Lambda] \[Phi]m^4- \[Lambda](\[CurlyPhi]\[Rho]-\[Phi]L[[2]]-\[Phi]m)^4)] &&\[Lambda]= a[[p-1]]*)
\[ScriptCapitalV]Bs[a4_,\[Rho]_?NumericQ,b4_,\[Phi]T_,v4_,\[Phi]m_,VT_] := 
	(VT \[Rho]^4)/4+(4 (2+3 b4^2 a4 \[Rho]^2))/(3 a4 (2+b4^2 a4 \[Rho]^2)^3)+1/4 a4 \[Rho]^4 \[Phi]m^4;

\[ScriptCapitalV]B[R_,v_,VL_,\[Phi]L_,Ns_,a_,b_,\[Phi]m_] := 
Module[{\[ScriptCapitalV],p=2,d=4},
	\[ScriptCapitalV] = 2\[Pi]^(d/2)/Gamma[d/2](
		\[ScriptCapitalV]Bs[a[[p-1]],R[[p]],b[[p-1]],\[Phi]L[[p]],v[[p-1]],\[Phi]m,VL[[p]]]-
		\[ScriptCapitalV]Bs[a[[p-1]],0.,b[[p-1]],\[Phi]L[[p]],v[[p-1]],\[Phi]m,VL[[p]]]);

	\[ScriptCapitalV] += 2\[Pi]^(d/2)/Gamma[d/2](Sum[
		32 a[[i]]^2/(d(d+2))  (R[[i+1]]^(2+d) - R[[i]]^(2+d))  + 
		8 a[[i]]*b[[i]]/(d-2) (R[[i+1]]^2 -R[[i]]^2)  +
		( VL[[i]]-VL[[Ns+1]]+ 8 a[[i]]( v[[i]] - \[Phi]L[[i]]))(R[[i+1]]^d-R[[i]]^d)/d  ,{i,p,Ns}]   )
];

\[CapitalLambda]B[initialR_?NumericQ,d_,VL_,\[Phi]L_,a4_,Ns_,backward_,\[Phi]m_]:=
Module[{R4,v4,b4},
	{R4,v4,b4} = BottomlessParameterRvb[initialR,\[Phi]L,a4,d,Ns,backward,\[Phi]m]; 

	(2-d)*\[ScriptCapitalT]B[R4,Ns,a4,b4]/(d*\[ScriptCapitalV]B[R4,v4,VL,\[Phi]L,Ns,a4,b4,\[Phi]m]) 
];


(* ::Subsubsection::Closed:: *)
(*BottomlessPotentialBounce*)


FindBounce::blminpos = "Wrong position of the minima in bottomless potential.";
FindBounce::blpathdef = "The path is deformed irregularly on the potential, try changing number of segments.";
FindBounce::blinitr = "Trivial solution founded, increase the number of segments or accuracy.";
FindBounce::blnrm = "The potential should be a polynomial of order 4.";

BottomlessPotentialBounce[V_,potentialPoints_,Ns_,noFields_,\[Phi]L_,dim_,maxIteR_,actionTolerance_,
ansatzInitialR_,aRinitial_,rule_,iter_,fields_]:= 
Module[{a,VL,initialR,R,v,b,T1,V1,\[Phi]m,cList,\[Lambda],v0},

	cList = CoefficientList[Expand@Normal@Series[V,{fields[[1]],Infinity,4}],fields[[1]]];
	If[Length[cList]===5,
		\[Lambda] = Abs@cList[[5]];
		v0 = cList[[4]]/(4*\[Lambda]);
		,
		Message[FindBounce::blnrm];
		Return[$Failed,Module]
	];
	
	VL = If[potentialPoints===None,Table[V/.rule[[s]],{s,Ns+1}],potentialPoints];	
	
	If[VL[[-1]]>=VL[[-2]],
		If[iter === 0,
			Message[FindBounce::blminpos];
			Return[$Failed,Module],
			Message[FindBounce::blpathdef];
			Return[$Failed,Module]
		]
	];
	
	a  = Table[ ( (VL[[s+1]]-VL[[s]])/(\[Phi]L[[s+1]]-\[Phi]L[[s]]) )/8,{s,Ns}]; 
	a[[1]] = \[Lambda];
	\[Phi]m = v0 + \[Phi]L[[-1]];
	initialR = FindInitialRadiusB[dim,VL,\[Phi]L,a,Ns,maxIteR,actionTolerance,ansatzInitialR,aRinitial,\[Phi]m]//Re;
	
	If[initialR<10^(-5.),
		Message[FindBounce::blinitr];
		Return[$Failed,Module]
	];
	
	{R,v,b} = BottomlessParameterRvb[initialR,\[Phi]L,a,dim,Ns,True,\[Phi]m]//Re;
			
	a = Join[a,{0.}];	
	T1 = \[ScriptCapitalT]B[R,Ns,a,b]; 
	V1 = \[ScriptCapitalV]B[R,v,VL,\[Phi]L,Ns,a,b,\[Phi]m];
	VL[[1]] = VL[[2]]+\[Lambda]*\[Phi]m^4;

	{V1+T1,VL,v,a,b,2,R,initialR}
];


(* ::Subsection::Closed:: *)
(*FindBounce*)


(* ::Subsubsection::Closed:: *)
(*BounceFunction*)


summaryBoxGraphics[bf_BounceFunction]:= BouncePlot[
	bf,
	(* Small plot has to render fast. *)
	PerformanceGoal->"Speed",
	FrameLabel->None, 
	FrameTicks->None,
	(* To avoid gray background before mouse-over. *)
	Background->White,
	GridLines->None,
	Axes->None,
	(* Set standard image size *)
	ImageSize -> Dynamic[{Automatic,3.5*CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]
];


BounceFunction::usage="BounceFunction object represents results from FindBounce function.";
BounceFunction::noprop="Property `1` of BounceFunction is not available. It should be one of `2`.";

(* In Mma 10. AssociationQ already works, even though it is undocumented. *)
BounceFunction[asc_?AssociationQ]["Properties"]:=Sort@Keys[asc];

BounceFunction[asc_?AssociationQ][key_]:=With[{
	value=Lookup[asc,key],
	supported=Sort@Keys[asc]
	},
	(* Check for explicit Missing["KeyAbsent",_] because other reasons like
	Missing["NotAvailable"] could be valid results. *)
	If[
		MatchQ[value,Missing["KeyAbsent",_]],
		Message[BounceFunction::noprop,Style[key,ShowStringCharacters->True],supported]
	];
	value
];

(* Nice styling of output, see https://mathematica.stackexchange.com/questions/77658 *)
BounceFunction/:MakeBoxes[obj:BounceFunction[asc_?AssociationQ],form:(StandardForm|TraditionalForm)]:=Module[
	{above,below,icon},
	
	 (* column *)
	above = {
		BoxForm`SummaryItem[{"Action: ", obj["Action"]}],
		BoxForm`SummaryItem[{"Dimension: ", obj["Dimension"]}]
	};
	below = {
		BoxForm`SummaryItem[{"Domain: ", obj["Radii"][[{1,-1}]]}],		
		BoxForm`SummaryItem[{"Path iterations: ", obj["PathIterations"]}],
		BoxForm`SummaryItem[{"Field points: ", Length@obj["Path"]}]
	};
	icon = summaryBoxGraphics[obj];
	
	BoxForm`ArrangeSummaryBox[
		BounceFunction, (* head *)
		obj,      (* interpretation *)
		Deploy@icon,    (* icon, use None if not needed *)
		(* above and below must be in a format suitable for Grid or Column *)
		above,    (* always shown content *)
		below,    (* expandable content *)
		form,
		(* "Interpretable"->Automatic works only in version 11.2+ *)
		"Interpretable" -> False
	]
];


(* ::Subsubsection::Closed:: *)
(*piecewiseBounce*)


piecewiseBounce[{v_,a_,b_,R_,\[Nu]_,\[Alpha]_,\[Beta]_,ddVL_},{\[Phi]_,dim_,pos_,noSegs_,noFields_,bottomless_,extensionPB_}]:=
Module[{\[CurlyPhi]0,\[CurlyPhi],min1=\[Phi][[1]],min2=\[Phi][[-1]],\[ScriptCapitalI],\[ScriptCapitalI]0,\[Xi],\[Xi]0},
	(*Polygonal bounce, as in eq. 4*)
	\[CurlyPhi][s_,i_,\[Rho]_]:= v[[s,i]]+4/dim*a[[s,i]]*\[Rho]^2+2/(dim-2)*b[[s,i]]/\[Rho]^(dim-2);
		
	(*If the extension the polygonal is considered or not*)
	If[extensionPB,				
		(*see equations 41 for the definition of \[Xi][\[Rho]]*)
		\[Xi][s_,i_,\[Rho]_]:= \[Nu][[s,i]]+4/dim*\[Alpha][[s,i]]*\[FormalR]^2+2/(dim-2)*\[Beta][[s,i]]/\[FormalR]^(dim-2);
		(*See eq 47,48 for \[ScriptCapitalI][\[Rho]] in D=3,4 repectively.*)
		If[
			dim==3,
			\[ScriptCapitalI][s_,i_,\[Rho]_]:= ddVL[[s]]( (v[[s,i]]-\[Phi][[s,i]])/6 \[Rho]^2 + a[[s,i]]/15 \[Rho]^4 + b[[s,i]]\[Rho]),
			(*dim\[Equal]4*)
			\[ScriptCapitalI][s_,i_,\[Rho]_]:= ddVL[[s]]( (v[[s,i]]-\[Phi][[s,i]])/8 \[Rho]^2 + a[[s,i]]/24 \[Rho]^4 + b[[s,i]]/2*Log[\[Rho]])
		];
		(*Note that consistently b is zero when \[Rho]=0, so the bounce is smooth. However, if b is computed numerical one should get rid of b by hand as:*)
		If[
			pos>1,
			(*Case A*)
			\[Xi]0[\[Rho]_]:=\[Nu][[pos-1]]+4/dim*\[Alpha][[pos-1]]*\[Rho]^2;
			If[
				dim==3,
				\[ScriptCapitalI]0[\[Rho]_]:= ddVL[[pos-1]]( (v[[pos-1]]-\[Phi][[pos-1]])/6 \[Rho]^2 + a[[pos-1]]/15 \[Rho]^4),
				(*dim\[Equal]4*)
				\[ScriptCapitalI]0[\[Rho]_]:= ddVL[[pos-1]]( (v[[pos-1]]-\[Phi][[pos-1]])/8 \[Rho]^2 + a[[pos-1]]/24 \[Rho]^4)
			];
			,
			(*Case B*)
			\[Xi]0[\[Rho]_]:=0;
			\[ScriptCapitalI]0[\[Rho]_]:=0;
		];
		,
		\[Xi]0[\[Rho]_]:=0;
		\[ScriptCapitalI]0[\[Rho]_]:=0;
		\[Xi][s_,i_,\[Rho]_]:=0.;
		\[ScriptCapitalI][s_,i_,\[Rho]_]:=0.
	];
	(*If the bottomless polygonal potential is considered or not*)
	If[bottomless,
		\[CurlyPhi]0[\[Rho]_] := v[[1]] + b[[1]]/(1 + 1/2 Norm[a[[1]]]*Norm[b[[1]]]^2\[Rho]^2)
		,
		If[
			pos>1,
			(*Case A*)
			\[CurlyPhi]0[\[Rho]_] := v[[pos-1]]+ 4/dim*a[[pos-1]]*\[Rho]^2+\[Xi]0[\[Rho]]+\[ScriptCapitalI]0[\[Rho]],
			(*Case B*)
			\[CurlyPhi]0[\[Rho]_] := min1 +\[Xi]0[\[Rho]]+\[ScriptCapitalI]0[\[Rho]]
		];
	];
	
	(* System symbol FormalR (\[FormalR]) is used as Function argument, because it has attribute
	Protected and cannot be assigned any value. Syntax with explicit Function is clearer than
	using pure functions. *)
	
	Table[
		Function[
			Evaluate@\[FormalR],
			Evaluate@Piecewise[
				Join[
					{{\[CurlyPhi]0[\[FormalR]][[i]],\[FormalR]<R[[pos]]}},
					Table[{
						\[CurlyPhi][s,i,\[FormalR]]+\[Xi][s,i,\[FormalR]]+\[ScriptCapitalI][s,i,\[FormalR]],R[[s]]<=\[FormalR]<R[[s+1]]},
						{s,pos,noSegs}
					]
				],
				min2[[i]] (* default value of Piecewise *)
			]
		],
		{i,noFields}
	]
];


(* ::Subsubsection::Closed:: *)
(*createGenericResults*)


(* Utility function to create generic data structure for results, before calculation is
even started. Values of certain fileds are updated in the end. *)
createGenericResults[fieldPoints_,dim_]:=Association[
	"Action"->Infinity,
	"Bounce"->Function[Evaluate@First@fieldPoints],
	"BottomlessPotential"->Missing["NotAvailable"],
	"Coefficients"->Missing["NotAvailable"],
	"Dimension"->dim,
	"PathIterations"->0,
	"Path"->fieldPoints,
	"Radii"->{0,Infinity}
];


(* ::Subsubsection::Closed:: *)
(*FindBounce*)


FindBounce::usage = 
	"FindBounce[V[x],x,{min1, min2}] computes false vacuum decay for potential V[x] for field x.\n"<>
	"FindBounce[V[x,y,...],{x,y,...},{min1, min2}] works with multiple scalar fields.\n"<>
	"FindBounce[{{x1,y1},{x2,y1},...}] works with single field potential given as a list of points.";
FindBounce::dim = "Only supported \"Dimension\"s are 3 and 4, default value was taken.";
FindBounce::posreal = "Value of option \"`1`\" should be a positive real number.";
FindBounce::posint = "Value of option \"`1`\" should be a positive integer.";
FindBounce::nonnegint = "Value of option \"`1`\" should be a non-negative integer.";
FindBounce::degeneracy = "Not vacuum decay, the vacua are degenerated.";
FindBounce::points = "Single field potential defined by points should be a n by 2 matrix of non-complex numbers, with n>=3.";
FindBounce::fieldpts = "\"FieldPoints\" should be an integer (n>2) or array of non-complex numbers longer than 2.";
FindBounce::syms = "Field symbols should not have any value.";
FindBounce::mins = "Dimensions of minima should be consistent to the number of fields, different and not complex.";

(* Option names are written in canonical (alphabetical) order. *)
Options[FindBounce] = {
	"ActionTolerance" -> 10.^(-6),
	"BottomlessPotential" -> False,
	"Dimension" -> 4,
	"FieldPoints" -> 31,
	"Gradient" -> Automatic,
	"Hessian" -> Automatic,
	"MaxPathIterations" -> 3,
	"MaxRadiusIterations" -> 100,
	"MidFieldPoint" -> None,
	"PathTolerance" -> 0.01
};

(* Autocomplete option names *)
With[
	{keys=ToString/@Sort@Keys@Options@FindBounce},
	addCodeCompletion["FindBounce"][Sequence@@Join[{0},ConstantArray[keys,Length@keys]]]
];

FindBounce//SyntaxInformation={
	"ArgumentsPattern"->{_,_.,_.,OptionsPattern[]},
	"LocalVariables"->{"Solve",{2,2}}
};

(* This definition transforms single field case to multi-field case with one field. *)
FindBounce[V_,field_Symbol,{minimum1_,minimum2_},opts:OptionsPattern[]]:=FindBounce[
	V,{field},{minimum1,minimum2},opts
];

(* Definition for a potential defined by a list of points. *)
FindBounce[points_List,opts:OptionsPattern[]]:=(
	If[
		Not@And[
			ArrayQ[N@points,2,(MatchQ[#,_Real]&)],
			MatchQ[Dimensions[points],{x_/;x>=3,2}]
		],
		Message[FindBounce::points];Return[$Failed]
	];
	FindBounce[points,{True},points[[{1,-1},1]],opts]
);	

FindBounce[V_,fields_List,{minimum1_,minimum2_},opts:OptionsPattern[]]:=
Module[{Ns,a,\[Phi]L,ansatzInitialR,b,v,\[Alpha],\[Beta],\[Nu],dim,noFields,VL,midPoint,fieldPoints,
	maxItePath,maxIteR,R,extensionPB,rule,pos,l,eL,dV,d2V,RM,
	actionP,action\[Xi]=0.,action,vM,aM,bM,posM,ddVL,bottomless,p,pathTolerance,actionTolerance,
	min1,min2,iter=0,potentialPoints=None,switchPath=False,initialR=None,results,potentialValues,
	gradient,hessian},
	
	(* First we check correctness of arguments and options, then we start to process them. *)
	(* Checks if field variables do not have any values.*)
	If[Not@ArrayQ[fields,1,(Head[#]===Symbol&)],Message[FindBounce::syms];Return[$Failed,Module]];
	
	(* Minima are transformed to a matrix of Dimensions {2,noFields}. *)
	noFields = Length[fields];
	{min1,min2}=N[Flatten/@{{minimum1},{minimum2}}];
	If[
		Not@And[
			ArrayQ[{min1,min2},2,(MatchQ[#,_Real]&)],
			Length[min1]==noFields,
			min1 != min2
		],
		Message[FindBounce::mins];Return[$Failed,Module]
	];
	
	(* Checking of acceptable option values. If they are wrong function immediately returns $Failed.*)
	bottomless = TrueQ@OptionValue["BottomlessPotential"];
	pathTolerance = OptionValue["PathTolerance"]/.Except[_?NonNegative]:>(Message[FindBounce::posreal,"PathTolerance"];Return[$Failed,Module]);
	actionTolerance = OptionValue["ActionTolerance"]/.Except[_?NonNegative]:>(Message[FindBounce::posreal,"ActionTolerance"];Return[$Failed,Module]);
	dim = OptionValue["Dimension"]/.Except[3|4]:>(Message[FindBounce::dim];Return[$Failed,Module]);
	maxIteR = OptionValue["MaxRadiusIterations"]/.Except[_Integer?Positive]:>(Message[FindBounce::posint,"MaxRadiusIterations"];Return[$Failed,Module]);
	maxItePath = OptionValue["MaxPathIterations"]/.Except[_Integer?NonNegative]:>(Message[FindBounce::nonnegint,"MaxPathIterations"];Return[$Failed,Module]);
	midPoint = N@OptionValue["MidFieldPoint"];

	fieldPoints=fieldSegmentation[{min1,min2},opts]/.($Failed:>Return[$Failed]);
	(* TODO: Decision if extension method is used should be taken in one place only.
	Review the following check if it catches all cases. *)
	extensionPB=And[noFields==1, OptionValue["Gradient"]=!=None, Not@bottomless];

	(* Special case  if the potential is given as a list of points.*)
	If[fields[[1]]===True,
		{fieldPoints,potentialPoints} = Transpose@V;
		fieldPoints=Partition[fieldPoints,1];
		potentialValues=potentialPoints;
		extensionPB=False,
		(* Otherwise we ask for value of potential at field points. *)
		potentialValues=getPotentialValues[V,fields,fieldPoints]/.($Failed:>Return[$Failed]);
	];

	(* Segmentation should always go from lower to higher minimum. *)
	If[
		First[potentialValues]>Last[potentialValues],
		fieldPoints=Reverse[fieldPoints];
		potentialValues=Reverse[potentialValues]
	];

	ansatzInitialR=initialRadiusEstimate[fieldPoints,potentialValues,dim];
	Ns=Length[fieldPoints]-1;

	results=createGenericResults[fieldPoints,dim];
	(*If the vacua are degenerated, return generic/initial results. *)
	If[
		Head[ansatzInitialR]===DirectedInfinity&&Not[bottomless],
		Message[FindBounce::degeneracy];Return[BounceFunction@results]
	];

	(* For single field there is no path deformation calculation. *)
	If[noFields==1,maxItePath = 0];
	(*Bounce and path deformation.*)
	While[iter <= maxItePath||switchPath,
		
		(*Rule.*)
		rule = Table[fields[[i]]->fieldPoints[[s,i]],{s,Ns+1},{i,noFields}];
		(* TODO: Eventually move this intermediate calculation inside each major function below. *)
		\[Phi]L=longitudinalProjection[fieldPoints];
		eL=unitVectors[fieldPoints];
		l=segmentsLength[fieldPoints];
		
		(*Single Field Bounce.*)
		{actionP,VL,v,a,b,pos,R,initialR} = If[
			Not[bottomless],
			SingleFieldBounce[V,potentialPoints,Ns,noFields,\[Phi]L,dim,maxIteR,
				actionTolerance,ansatzInitialR,initialR,rule,iter,
				switchPath||iter==maxItePath
			]/.(x_/;Not@FreeQ[x,$Failed]:>Return[$Failed,Module])
			,
			BottomlessPotentialBounce[V,potentialPoints,Ns,noFields,\[Phi]L,
				dim,maxIteR,actionTolerance,ansatzInitialR,initialR,rule,iter,fields
			]/.(x_/;Not@FreeQ[x,$Failed]:>Return[$Failed,Module])
		];
		
		If[
			extensionPB,
			gradient=getPotentialGradient[V,fields,fieldPoints,opts]/.($Failed:>Return[$Failed]);
			{action\[Xi],\[Alpha],\[Beta],\[Nu],ddVL,extensionPB} = SingleFieldBounceExtension[VL,Ns,v,a,b,R,\[Phi]L,pos,dim,eL,gradient];
		];
				
		(*Transforms \[Phi]L,v,a,b (logitudinal) into \[Phi] (field space) and its bounce parameters.*)
		{v,a,b,\[Alpha],\[Beta],\[Nu],fieldPoints,action,switchPath} = ParameterInFieldSpace[v,a,b,\[Alpha],\[Beta],\[Nu],fieldPoints,eL,l,\[Phi]L,Ns,noFields,pos,dim,
			bottomless,action,actionP+action\[Xi],actionTolerance,switchPath,extensionPB];

		(*Breaks the interations of path deformation.*)
		If[switchPath||iter == maxItePath||bottomless, 
			p = If[pos>1,pos-1,pos]; 
			Break[]
		];

		(*Multi Field Bounce.*)
		gradient=getPotentialGradient[V,fields,fieldPoints,opts]/.($Failed:>Return[$Failed]);
		hessian=getPotentialHessian[V,fields,fieldPoints,opts]/.($Failed:>Return[$Failed]);
		{fieldPoints,vM,aM,bM,RM,posM,switchPath} = MultiFieldBounce[fieldPoints,gradient,hessian,noFields,pos,
			dim,R,v,a,b,\[Phi]L[[-1]],pathTolerance];

		iter++
	];

	results["Action"]=action;
	results["Bounce"]=piecewiseBounce[{v,a,b,R,\[Nu],\[Alpha],\[Beta],ddVL},{fieldPoints,dim,pos,Ns,noFields,bottomless,extensionPB}];
	results["BottomlessPotential"]=If[bottomless,VL[[1]],Missing["NotAvailable"]];
	results["Coefficients"]={v,a,b}[[All,p;;-1]];
	results["CoefficientsExtension"]=If[extensionPB,{\[Nu],\[Alpha],\[Beta],ddVL}[[All,p;;-1]],Missing["NotAvailable"]];
	results["PathIterations"]=iter;
	results["Path"]=fieldPoints;
	results["Radii"]=R[[p;;-1]];

	BounceFunction@results
];


(* ::Subsection::Closed:: *)
(*BouncePlot*)


BouncePlot::usage="BouncePlot[bf] plots the content of BounceFunction bf.\n"<>
"BouncePlot[{bf1,bf2}] plots several functions bfi.";

BouncePlot//Options=Options@Plot;

BouncePlot//SyntaxInformation={"ArgumentsPattern"->{_,OptionsPattern[]}};

BouncePlot[bf_BounceFunction,opts:OptionsPattern[]]:= BouncePlot[{bf},opts];
	
BouncePlot[{bf__BounceFunction},opts:OptionsPattern[]]:= Module[
	{bounce,radii,defaultPlotRange,plotRange,rangeMin,rangeMax,bounceExtension},
	(* In case of degenerated vacua, "Action" is Infinity and empty plot is returned. 
	This is suitable form for FindBounce summary box plot.*)
	If[
		Not@AllTrue[Flatten@Through[{bf}["Action"]],NumberQ],
		Return[
			Graphics[
				{Text[Style["n/a",Large,Gray]]},
				AspectRatio->1/GoldenRatio,
				Frame->True,
				Evaluate@FilterRules[{opts},Options@Graphics]
			],
			Module
		]
	];
	
	(* Piecewise bounces of consecutive BounceFunction(s) are effectively flattened. *)
	bounce = Flatten@Through[{bf}["Bounce"]];
	(* This helps to draw discrete radii. *)
	radii = Through[{bf}["Radii"]];
	(* Clip estimated plot range to non-negative values. *)
	defaultPlotRange = MinMax[radii,Scaled[0.25]];
	(* With this flat parts of curve are ploted if explicit value for option PlotRange is given.*)
	(* Clip is used for backward comparibility, because Ramp is introduced in Mma 11. *)
	plotRange=Quiet[
		Clip[#,{0.,Infinity}]&@If[
			MatchQ[OptionValue[PlotRange],{{_?NumberQ,_?NumberQ},_}],
			First@OptionValue[PlotRange],
			defaultPlotRange
		],
		OptionValue::nodef
	];
	(* This fiddling is neccesary for compatibility with Mma 10. becase evaluation sequence 
	for Plot has changed. *)
	{rangeMin,rangeMax}=plotRange;

	Plot[
		Evaluate@Through[bounce[r]],
		{r,rangeMin,rangeMax},
		(* Exclusions->None is very important because it dramatically speeds up plotting. 
		We know that function is continious and there is no need to search for discontinuities. *)
		Exclusions->None,
		Evaluate@FilterRules[{opts},Options@Plot],
		(* Default options come here. They can be overridden by other options given explicitly. *)
		Frame->True,
		FrameLabel->{"\[Rho]","\[CurlyPhi](\[Rho])"},
		LabelStyle->Directive[Black, FontSize->17,FontFamily->"Times New Roman"],
		GridLines->Automatic
	]
];


(* ::Section::Closed:: *)
(*End Package*)


End[]; (*"`Private`"*)


(* ReadProtected attribute on public symbols prevents rendering of huge box with all 
definitions (DownValues) when they are called in Information or with shortcut ?FindBounce. *)
SetAttributes[Evaluate@Names["`*"],{ReadProtected}];


EndPackage[];
