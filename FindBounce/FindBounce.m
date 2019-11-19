(* ::Package:: *)

(* ::Chapter::Closed:: *)
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


(* ::Chapter::Closed:: *)
(*Begin Package*)


BeginPackage["FindBounce`"];


(* ::Section::Closed:: *)
(*Available public functions*)


(* The main function of the package. *)
FindBounce;
(* Symbol which acts as a container for results of FindBounce. *)
BounceFunction;
(* A handy function for plotting the contents of BounceFunction. *)
BouncePlot;


(* Clear definitions from package symbols in public and private context. *)
ClearAll["`*","`*`*"];


(* ::Section::Closed:: *)
(*Begin Private*)


Begin["`Private`"];


(* ::Chapter::Closed:: *)
(*Code*)


(* ::Section::Closed:: *)
(*Utilities*)


(* ::Subsection::Closed:: *)
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


If[
	$VersionNumber<10.1,
	MinMax=minMax
];


(* ::Section::Closed:: *)
(*InitialValue*)


(* ::Subsection::Closed:: *)
(*Segmentation*)


Options[Segmentation] = {"FieldPoints" -> 31,"Method" -> "H"};

Segmentation//SyntaxInformation={
	"ArgumentsPattern"->{{_,_,_},OptionsPattern[]}
};

Segmentation[\[Phi]3_,OptionsPattern[]]:=
Module[{fieldPoints,\[Phi],\[Delta]\[Phi],\[Delta]\[Phi]1,\[Delta]\[Phi]2,\[CapitalPhi],np,n1,n2}, 
	fieldPoints = OptionValue["FieldPoints"];
	If[ fieldPoints <= 3, 
		
		Return[\[Phi]3], 
		(*Homogeneous_Segmentation*)
		If[OptionValue["Method"] === "H" , 
			\[Delta]\[Phi]= Abs[ (\[Phi]3[[3]] - \[Phi]3[[1]]) ]/(fieldPoints-1);
			\[CapitalPhi] = Table[ i , { i, \[Phi]3[[1]], \[Phi]3[[3]], \[Delta]\[Phi]} ];
		
		Return[\[CapitalPhi]] 
		];
		
		(*Bi-Homogeneous_Segmentation*)
		If[OptionValue["Method"] === "2H", 
			\[Delta]\[Phi] = Abs[ (\[Phi]3[[3]] - \[Phi]3[[1]]) ]/(fieldPoints-1);
			n1 = Quotient[Abs[ (\[Phi]3[[2]] - \[Phi]3[[1]]) ], \[Delta]\[Phi] ];
			n2 = IntegerPart[ (fieldPoints-1) - n1 ];
			\[Delta]\[Phi]1 = Abs[ (\[Phi]3[[2]] - \[Phi]3[[1]]) ]/n1;
			\[Delta]\[Phi]2 = Abs[ (\[Phi]3[[3]] - \[Phi]3[[2]]) ]/n2;
			\[CapitalPhi] = Join[   Table[ i , { i, \[Phi]3[[1]]  , \[Phi]3[[2]] - \[Delta]\[Phi]1 ,\[Delta]\[Phi]1  }  ],  
				Table[ i , { i, \[Phi]3[[2]], \[Phi]3[[3]] ,\[Delta]\[Phi]2  } ]    ];
		
		Return[\[CapitalPhi]] 
		];       
	]    
];


(* ::Subsection::Closed:: *)
(*NewAnsatz*)


NewAnsatz[\[Phi]_,Ns_]:= 
Module[{l,\[Phi]L,eL},
	l  = Table[Norm[{\[Phi][[s+1]]-\[Phi][[s]]}],{s,1,Ns}];
	\[Phi]L = Table[Sum[l[[s1]],{s1,1,s-1}],{s,1,Ns+1}];
	eL = (\[Phi][[2;;-1]]-\[Phi][[1;;-2]])/l;
	
	{Length[\[Phi]L]-1,\[Phi]L,eL,l} 
];


(* ::Subsection::Closed:: *)
(*DerivativePotential*)


DerivativePotential::gradient = "\"Gradient\" is not a vector, default value was taken.";
DerivativePotential::hessian = "\"Hessian\" is not a matrix, default value was taken."; 

DerivativePotential[V_,fields_,noFields_,gradient_,hessian_]:=
Module[{dV,d2V},
	If[
		ArrayQ[gradient,1],
		dV = gradient,
		If[Not[gradient === None ||gradient === Automatic],
		Message[DerivativePotential::gradient]
	];
		dV = D[V,{fields}]
	];

	If[
		ArrayQ[hessian,2]&&noFields>1,
		d2V = hessian
	,
		If[hessian =!= None,
		Message[DerivativePotential::hessian]
	];
		d2V = D[V,{fields},{fields}]
	];

	{dV,d2V}
];


(* ::Subsection::Closed:: *)
(*InitialValue*)


InitialValue::wrongInput = "Wrong \"`1`\".";
InitialValue::dimArray = "The array dimention of min1, min2 and fields are inconsistent.";
InitialValue::mpts = "\"MidFieldPoint\" should be a vector of lenght equal to the number of fields.";

InitialValue[V_,fields_,noFields_,min1_,Point_,min2_,potentialPoints_,
gradient_,hessian_,dim_,bottomless_,fieldpoints_]:=
Module[{VL,Ns,\[Phi],\[Phi]L,eL,l,rule,\[Phi]3,VL3,L3,\[Phi]L3,eL3,a3,initialR,
fieldPoints,dV = None,d2V=None,improvePB=False,c,\[CurlyPhi]0,point = Point,methodSeg,path },

If[point === None,
		point = (min1+min2)/2;
		methodSeg = "H"
		,
		If[Not[Length[N@point]===noFields||(noFields===1&&Length[N@point]===0)],
			Message[InitialValue::mpts];
			Return[$Failed,Module]
		];
		methodSeg = "2H"
	];

	If[ Head[fieldpoints] === Integer,
		Ns = fieldpoints-1;
		path = None,
		Ns = Length[fieldpoints]-1;
		path = fieldpoints
	];

	If[path === None, 
		\[Phi] = N@{min1,point,min2},
		\[Phi] = N@path
	];
	fieldPoints = Length[\[Phi]];

	rule = If[
		noFields == 1,
		Table[fields[[1]]->\[Phi][[s]],{s,fieldPoints}], 
		Table[fields[[i]]->\[Phi][[s,i]],{s,fieldPoints},{i,noFields}]
	];
	
	VL = If[
		potentialPoints===None, 
		Table[V/.rule[[s]],{s,fieldPoints}], 
		{potentialPoints[[1]],Max[potentialPoints],potentialPoints[[-1]]}
	];

	(*Checks if the values of the potential are well definited*)
	If[!NumericQ[VL[[1]]] || !NumericQ[VL[[2]]] || !NumericQ[VL[[-1]]],
	Message[InitialValue::wrongInput,"Potential"];
	Return[$Failed,Module]
	];
	 
(*Checks the dimension of the field values*)
	If[Length[\[Phi][[1]]] =!= Length[\[Phi][[2]]] || Length[\[Phi][[2]]] =!= Length[\[Phi][[-1]]],
	Message[InitialValue::dimArray];
	Return[$Failed,Module]   
	];

	(*The number 3 stands for Number of Field Values: Segments+1 = 3*) 
	VL3 = MinMax[VL[[{1,-1}]]];
	VL3 = Insert[VL3,Max[VL[[2;;-2]]],2];
	(*Sets V[1]<V[3]<V[2], i.e the local minimum on the right.*) 
	If[ 
		VL3[[1]] != VL[[1]] ,
		\[Phi] = Reverse[\[Phi]]
	];

	(*Finds the Logitudinal field values \[Phi]L3, fields \[Phi]3, distance l3, and direction eL3 for 3 Segments *)
	\[Phi]3 = \[Phi][[{1,Position[ VL,VL3[[2]] ][[1,1]],-1}]];
	L3 = Table[Norm[\[Phi]3[[s+1]]-\[Phi]3[[s]]],{s,2}];
	\[Phi]L3 = Table[Sum[L3[[s1]],{s1,s-1}],{s,3}];
	eL3 = Table[(\[Phi]3[[s+1]]-\[Phi]3[[s]])/L3[[s]],{s,2}];
	a3 = Table[(VL3[[s+1]]-VL3[[s]])/(\[Phi]L3[[s+1]]-\[Phi]L3[[s]]) 1/8,{s,2}]  ;

	(*Estimates the initial radii initialR with the N=2 close form solution*)
	initialR = 
If[VL3[[1]]!= VL3[[3]] ,
	c = dim/(dim-2) (a3[[2]]-a3[[1]])/a3[[1]] (1 -(a3[[2]]/(a3[[2]]-a3[[1]]))^(1-2/dim));
	\[CurlyPhi]0 =( ( \[Phi]L3[[3]] + c \[Phi]L3[[2]] )/(1+c)  );
	If[Re@\[CurlyPhi]0>=0&&Im@\[CurlyPhi]0==0,
		Sqrt[dim/4(\[Phi]L3[[2]]-\[CurlyPhi]0)/a3[[1]]]
				,
				1/2(\[Phi]L3[[3]]-\[Phi]L3[[1]])/(Sqrt[a3[[1]] (\[Phi]L3[[2]]-\[Phi]L3[[1]])]-Sqrt[-a3[[2]] (\[Phi]L3[[3]]-\[Phi]L3[[2]])])
	]
	,
	Infinity
];

	(*Finds the Logitudinal field values \[Phi]L, fields \[Phi], distance l, and direction eL*)
If[path === None,
		\[Phi]L = Segmentation[\[Phi]L3,"FieldPoints" ->Ns +1 ,"Method"->methodSeg];
		\[Phi] = Table[ 
			If[ \[Phi]L[[s]]< \[Phi]L3[[2]],
				eL3[[1]](\[Phi]L[[s]]-L3[[1]])+ \[Phi]3[[2]],
				eL3[[2]](\[Phi]L[[s]]-(L3[[1]] + L3[[2]]))+\[Phi]3[[3]]
			]
		,{s,1,Length[\[Phi]L]} ];
		l = \[Phi]L[[2;;-1]]-\[Phi]L[[1;;-2]];
	,
		l  = Table[Norm[{\[Phi][[s+1]]-\[Phi][[s]]}],{s,fieldPoints-1}];
		\[Phi]L = Table[Sum[l[[s1]],{s1,s-1}],{s,fieldPoints}];
	];

	If[potentialPoints===None&&Not[bottomless],
		If[Not[gradient === None &&noFields==1],
	{dV,d2V} = DerivativePotential[V,fields,noFields,gradient,hessian];
	If[Ns>3&&Not[gradient === None],
		improvePB = True
	];
		];
	];
	eL = (\[Phi][[2;;-1]]-\[Phi][[1;;-2]])/l;
	
	{initialR,Length[\[Phi]L]-1,\[Phi],\[Phi]L,eL,l,dV,d2V,improvePB,path}
];


(* ::Section::Closed:: *)
(*SingleFieldBounce*)


(* ::Subsection::Closed:: *)
(*FindSegment*)


FindSegment[a_,\[Phi]L_,d_,Ns_]:=
Module[{pos = 1,R,estimateRmin,estimateRmax,ps},
	R = BounceParameterRvb[0,\[Phi]L,a,d,Ns,False,pos][[1]];
	
	While[Im[R[[-1]]] == 0 && pos< Ns,
		pos++;
		R = BounceParameterRvb[0,\[Phi]L,a,d,Ns,False,pos][[1]]; 
	];
	
	pos
];


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
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
Module[{R,b,v,\[Alpha],v1,b1,x,y,z,Rvb,position=pos},
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


(* ::Subsection::Closed:: *)
(*FindInitialRadius*)


FindInitialRadius::noSolution = "Large error solving the boundaries conditions. This potential may not have a bounce solution.";

(*Find the solution of eq. 25 or 26.*)
FindInitialRadius[d_,VL_,\[Phi]L_,a_,Ns_,maxIteR_,accuracyRadius_,ansatzInitialR_,aRinitial_,pos_,switchMessage_]:= 
Module[{R,v,b,radii,initialR,ite,complexR,realR,switch,k,initialR0,complexityR,rw,findInitialR,\[Lambda],lambda,Rvb,kinetic,potential},
(*The initial Radius can be found simply with FindRoot[R[x],{x,x0}]. However, this is not always the case since FindRoot
use Newton's Method and it fails if it hits a singularity. Thereby, in order to have a robust mechanism one check the output
of FindRoot after each iterations and use the bisection method in case FindRoot fails.*)
	
	(*Definitios of internal functions*)
	(*Saves values of BounceParameter,kinetic and Potentail energy and simplifies the notation*)
	
	(*Defines \[Lambda] of eq. 26*)
	lambda[initialR_?NumericQ] := lambda[initialR] = (
		Rvb = BounceParameterRvb[initialR,\[Phi]L,a,d,Ns,True,pos];
		kinetic = \[ScriptCapitalT][Sequence@@Rvb,a,d,Ns,pos];
		potential = \[ScriptCapitalV][Sequence@@Rvb,a,d,VL,\[Phi]L,Ns];
		(*complexityR: Discriminates between undershooting and overshooting*)
		complexityR = Abs@Im@Chop@Rvb[[1,1]];
		
		\[Lambda] = Sqrt[(2-d)*kinetic/(d*potential)]
	);
	
	(*Looks for the initial radius*)
	findInitialR[initialR_]:= 
		Quiet[Abs[rw/.FindRoot[Abs[Re@lambda[rw]-1],{rw,initialR}, 
					MaxIterations->1,
					PrecisionGoal->0,
					AccuracyGoal->accuracyRadius]
				],
		{FindRoot::lstol,FindRoot::cvmit}];   
	
If[Ns==2&&Not[d==3],(*initialR = exact radius from the close form solution N=3*)
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
			While[ complexityR > 0 &&ite <= maxIteR&&Abs[\[Lambda]-1]> 0.5*10^(-accuracyRadius)&&Chop[Re@\[Lambda]]!=0,
				If[initialR < complexR,
					complexR = initialR;
					initialR = findInitialR[initialR];
					,
					(*Perturbs 10% down the ansatz since initialR>complexR*) 
					initialR = 0.9*complexR
				];
				lambda[initialR];
				ite++
			]; 
			realR = initialR;
			,(*Undershooting*)
			While[ complexityR == 0&&ite <= maxIteR&&Abs[\[Lambda]-1] > 0.5*10^(-accuracyRadius)&&Chop[Re@\[Lambda]]!=0,
				If[initialR > realR,
					realR = initialR;  
					initialR = findInitialR[initialR],
					(*Perturbs 10% up the ansatz since initialR<realR*)  
					initialR = 1.1*realR
				];
				lambda[initialR];
				ite++
			]; 
			complexR = initialR;    
		];  
		(*One the interval is found, reduces the interval and use bisection method*)	
		While[ite <= maxIteR&&Abs[\[Lambda]-1]> 0.5*10^(-accuracyRadius)&&Chop[Re@\[Lambda]]!=0, 
			
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
			lambda[initialR];
			ite++   
		];
	   
		If[ Chop[Re@\[Lambda]] == 0,
			k++;
			switch=True;
			initialR = (1+k)*Abs[initialR0]  
		];   
	];    
	
	If[ ite > maxIteR&&switchMessage&& Abs[\[Lambda]-1]> 0.5*10^(-accuracyRadius*0.9),
		Message[FindInitialRadius::cvmit,maxIteR] 
	];

	If[ Re[\[Lambda]-1] > 0.5*10^(-accuracyRadius*0.1)&&switchMessage, 
		Message[FindInitialRadius::noSolution];
		Return[$Failed,Module]
	];	
];	
	Clear[lambda];
			
	Re@{initialR,Sequence@@Rvb,potential,kinetic}
]; 


(* ::Subsection::Closed:: *)
(*\[ScriptCapitalT]*)


(*See eqs. 10*)
\[ScriptCapitalT][R_,v_,b_,a_,d_,Ns_,pos_]:= 
Module[{p,T},
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


(* ::Subsection::Closed:: *)
(*\[ScriptCapitalV]*)


 (*See eq. 11*)
\[ScriptCapitalV][R_,v_,b_,a_,d_,VL_,\[Phi]L_,Ns_]:= 
	2\[Pi]^(d/2)/Gamma[d/2](Sum[   
	32*a[[i]]^2/(d(d+2))  (R[[i+1]]^(2+d) - R[[i]]^(2+d)) + 
	8*a[[i]]*b[[i]]/(d-2) (R[[i+1]]^2 -R[[i]]^2) +
	( VL[[i]]-VL[[-1]]+ 8*a[[i]]( v[[i]] - \[Phi]L[[i]]))(R[[i+1]]^d-R[[i]]^d)/d ,{i,1,Ns}] +
		1/d R[[1]]^d (VL[[1]] - VL[[-1]]) );


(* ::Subsection::Closed:: *)
(*SingleFieldBounce*)


SingleFieldBounce::extrema = "Wrong position of the minima.";
MultiFieldBounce::pathDeformation = "The path is deformed irregularly on the potential. Verifies that the vacuum is a minimum of the potential (not a saddle point) or changes the number of segements.";
SingleFieldBounce::noSolution = "Solution not found, increase the number of segments or accuracy.";

SingleFieldBounce[V_,potentialPoints_,Ns_,noFields_,\[Phi]L_,dim_,maxIteR_,accuracyRadius_,
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
			Message[SingleFieldBounce::extrema];
			Return[$Failed,Module]
			,
			Message[MultiFieldBounce::pathDeformation];
			Return[$Failed,Module]
		]
	];
	
	a  = Table[ ((VL[[s+1]]-VL[[s]])/(\[Phi]L[[s+1]]-\[Phi]L[[s]]))/8 ,{s,Ns} ]; 
	pos = FindSegment[a,\[Phi]L,dim,Ns];
	{initialR,R,v,b,V1,T1} = FindInitialRadius[dim,VL,\[Phi]L,a,Ns,maxIteR,accuracyRadius,ansatzInitialR,
		aRinitial,pos,switchMessage]/.x_/;FailureQ[x]:>Return[$Failed,Module];
		
	(*Checks if we got a consistent answer.*)
	If[V1+T1<0&&switchMessage,
		Message[SingleFieldBounce::noSolution];
		Return[$Failed,Module]
	];
	
	If[pos>1, R[[pos-1]]=0 ];
	a = Join[a,{0}];	

	{V1+T1,VL,v,a,b,pos,R,initialR}
];


(* ::Section::Closed:: *)
(*SingleFieldBounceImprovement*)


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*BounceParameterr\[Beta]\[Nu]*)


(*See eqs. 50-52*)
BounceParameterr\[Beta]\[Nu][rw_?NumericQ,a_,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_] := 
Module[{r\[Beta]\[Nu]M,r,\[Beta],\[Nu],x,y,z,\[Beta]prev,\[Alpha]0,a0,c0,b0,c},
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


(* ::Subsection::Closed:: *)
(*FindInitialRadiiImprovement*)


(*Solves Boundary Conditions \[Xi][N-1] = 0*)
FindInitialRadiiImprovement[rw_?NumericQ,a_ ,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_]:=
Module[{r,\[Beta],\[Nu]},
	{r,\[Beta],\[Nu]} = BounceParameterr\[Beta]\[Nu][rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos] ;
		
	\[Nu][[-2]] + 2/(d-2) \[Beta][[-2]]*R[[-1]]^(2-d)+4/d \[Alpha][[Ns]]*R[[-1]]^2+\[ScriptCapitalI][[1,-1]]   
];


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*SingleFieldBounceImprovement*)


SingleFieldBounceImprovement::dVFailed = "The first derivative of the potential is not well defined at some field value. \"Gradient\"->None was taken.";

SingleFieldBounceImprovement[VL_,dV_,noFields_,rule_,Ns_,v_,a_,b_,R_,\[Phi]L_,pos_,dim_,eL_,improvePB_]:=
Module[{dVL,\[Alpha],\[ScriptCapitalI],d\[ScriptCapitalI],r1,rInitial,r,\[Beta],\[Nu],eL0,T\[Xi]=0,V\[Xi]=0,ddVL = Missing["NotAvailable"]},
	eL0 = Join[eL,{eL[[-1]]}];
	If[improvePB,	
		dVL= If[
			noFields==1,
			Table[(dV[[1]]/.rule[[s]])*eL0[[s]],{s,Ns+1}],
			Table[(dV/.rule[[s]]).eL0[[s]],{s,Ns+1}]
		];
		If[And@@(NumericQ[#]&/@dVL),
			\[Alpha]  = Join[a[[1;;Ns]] - dVL[[2;;Ns+1]]/8 ,{0}];
			ddVL = Table[ (dVL[[s+1]]-8(a[[s]]+\[Alpha][[s]]))/(\[Phi]L[[s+1]]-\[Phi]L[[s]]),{s,Ns}];
			{\[ScriptCapitalI],d\[ScriptCapitalI]} = Find\[ScriptCapitalI][v,\[Phi]L,a,b,Ns,pos,R,ddVL,dim];
			r1 = rInitial/.FindRoot[FindInitialRadiiImprovement[rInitial,a,b,dim,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos],{rInitial,-1}]//Quiet;
			{r,\[Beta],\[Nu]} = If[
				pos>1,
				Join[ConstantArray[0,{pos-2,3}],
				BounceParameterr\[Beta]\[Nu][r1,a,b,dim,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos]//Transpose]//Transpose,
				BounceParameterr\[Beta]\[Nu][r1,a,b,dim,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos][[All,2;;-1]]
			];
			T\[Xi] = \[ScriptCapitalT]\[Xi][dim,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r];
			V\[Xi] = \[ScriptCapitalV]\[Xi][dim,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r];
			,
			Message[SingleFieldBounceImprovement::dVFailed]
		]
	];
		
	{V\[Xi]+T\[Xi],ddVL}
];


(* ::Section::Closed:: *)
(*ParameterInFieldSpace*)


ParameterInFieldSpace[vs_,as_,bs_,R_,\[Phi]_,eL_,l_,\[Phi]L_,Ns_,noFields_,pos_,dim_,bottomless_,actionOld_,actionNew_,dAction_,switchPathOld_]:=
Module[{v,a,b,path=\[Phi],switchPath = switchPathOld},
	v = Table[\[Phi][[s+1]]+eL[[s]]*(vs[[s]]-(l[[s]]+\[Phi]L[[s]])),{s,Ns}];
	a = Table[eL[[s]]*as[[s]],{s,Ns}];
	b = Table[eL[[s]]*bs[[s]],{s,Ns}]; 
	
	If[bottomless,
		path[[1]] = v[[1]]+b[[1]]
	];
	
	If[Abs[(actionOld-actionNew)/actionNew] < dAction,
		switchPath = True
	];
	
	{v,a,b,path,actionNew,switchPath}
];


(* ::Section::Closed:: *)
(*MultiFieldBounce*)


MultiFieldBounce[fields_,dV_,d2V_,Ns_,noFields_,pos_,d_,R0_,\[CapitalPhi]0_,\[ScriptV]0_,\[ScriptA]0_,\[ScriptB]0_,lengthPath_,dPath_]:=
Module[{\[Nu],\[Beta],rI,a,\[Zeta]t,R,\[Zeta]ts,\[Phi]M,\[Nu]\[Beta],x,y,d\[CurlyPhi],rF,rN,DV,D2V,
	\[Xi]Mc,M,c,\[Nu]0,\[Beta]0,\[Nu]\[Xi]p,\[Nu]\[Xi]m,\[Beta]\[Xi]p,\[Beta]\[Xi]m,fLowT,fD,fD1,frI,n1,rules,p,rulesQuartic, 
	\[Phi]0 = \[CapitalPhi]0, v0 = \[ScriptV]0, a0 = \[ScriptA]0, b0 = \[ScriptB]0,switchPath = False,n,m},
		
	If[pos>1,
		p = pos-1;
		\[Phi]0[[p]] = v0[[p]]
		,
		p = pos
	];

	R = Transpose@Table[R0,{i,noFields}];
	rules = Table[fields[[i]]->\[Phi]0[[s,i]],{s,p,Ns+1},{i,noFields}];

	(*Derivative of Poential and PB*)
	DV = Chop@Join[ConstantArray[0,{p-1,noFields}],Table[dV/.rules[[s-p+1]],{s,p,Ns+1}]];
	D2V = Chop@Join[Table[ConstantArray[0,{noFields,noFields}],{s,1,p-1}], 
		Table[d2V/.rules[[s-p+1]],{s,p,Ns+1}]];	
	d\[CurlyPhi] = Chop@Join[Table[ConstantArray[0,{2,noFields}],{s,1,p-1}],
		Table[8/d a0[[s+m]]*R[[s+1]]- 2 b0[[s+m]]/R[[s+1]]^(d-1),{s,p,Ns-1},{m,0,1}]];
		
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
	 
	{\[Phi]0[[1]],\[Zeta]ts[[1]]} = {\[CapitalPhi]0[[1]],ConstantArray[0,noFields]};

	If[Max[Abs@\[Zeta]ts]/lengthPath < dPath,
		switchPath = True
	];
	Clear[\[Nu]\[Xi]p,\[Beta]\[Xi]p,\[Nu]\[Xi]m,\[Beta]\[Xi]m,\[Nu]0,\[Beta]0];
	
	{\[Phi]0+\[Zeta]ts,v0+\[Nu]\[Beta][[1]],a0+a,b0+\[Nu]\[Beta][[2]],R0,pos,switchPath}   	
];


(* ::Section::Closed:: *)
(*FindBounce*)


(* ::Subsection::Closed:: *)
(*BounceFunction: Summary Box*)


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

BounceFunction[asc_?AssociationQ][property_]:=With[{
	value=Lookup[asc,property],
	supported=Sort@Keys[asc]
	},
	If[
		MatchQ[value,_Missing],
		Message[BounceFunction::noprop,property,supported]
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
		BoxForm`SummaryItem[{"IterationsPath: ", obj["PathIterations"]}],
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


(* ::Subsection::Closed:: *)
(*Process results: piecewiseBounce*)


piecewiseBounce[{v_,a_,b_,R_},{min1_,min2_},{dim_,pos_,noSegs_,noFields_,bottomless_}]:=Module[
	{\[CurlyPhi]0},

	If[bottomless,
		\[CurlyPhi]0[\[Rho]_] := v[[1]] + b[[1]]/(1 + 1/2 Norm[a[[1]]]*Norm[b[[1]]]^2\[Rho]^2),
		If[
			pos>1,
			(*Case A*)
			\[CurlyPhi]0[\[Rho]_] := v[[pos-1]]+ 4/dim*a[[pos-1]]*\[Rho]^2,
			(*Case B*)
			\[CurlyPhi]0[\[Rho]_] := min1 
		];
	];
	(* System symbol FormalR (\[FormalR]) is used as Function argument, because it has attribute
	Protected and cannot be assigned any value. Syntax with explicit Function is clearer than
	using pure functions. *)
	If[noFields==1,
		Function[
			Evaluate@\[FormalR],
			Evaluate@Piecewise[
				Join[
					{{\[CurlyPhi]0[\[FormalR]],\[FormalR]<R[[pos]]}},
					Table[{
						v[[s]]+4/dim*a[[s]]*\[FormalR]^2+2/(dim-2)*b[[s]]/\[FormalR]^(dim-2),R[[s]]<=\[FormalR]<R[[s+1]]},
						{s,pos,noSegs}
					]
				],
				min2 (* default value of Piecewise *)
			]
		],(* else (multifield case) *)
		Table[
		Function[
			Evaluate@\[FormalR],
			Evaluate@Piecewise[
				Join[
					{{\[CurlyPhi]0[\[FormalR]][[i]],\[FormalR]<R[[pos]]}},
					Table[{
						v[[s,i]]+4/dim*a[[s,i]]*\[FormalR]^2+2/(dim-2)*b[[s,i]]/\[FormalR]^(dim-2),R[[s]]<=\[FormalR]<R[[s+1]]},
						{s,pos,noSegs}
					]
				],
				min2[[i]] (* default value of Piecewise *)
			]
		],
		{i,noFields}
		]	
	]
];


(* ::Subsection::Closed:: *)
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
FindBounce::points = "Single field potential defined by points should be a n by 2 matrix of reals or integers, with n>=3.";
FindBounce::fieldpts = "\"FieldPoints\" should be an integer (n>2) or array of numbers longer than 2.";
FindBounce::syms = "Field symbols should not have any value.";
FindBounce::mins = "Dimensions of minima should be consistent to number of fields and not complex.";

Options[FindBounce] = {
	"BottomlessPotential" -> False,
	"PathTolerance" -> 0.01,
	"ActionTolerance" -> 0.01,
	"Dimension" -> 4,
	"FieldPoints" -> 31,
	Gradient -> Automatic,
	Hessian -> None,
	"InitialRadiusAccuracyGoal" -> 10,
	"MaxPathIterations" -> 3,
	"MaxRadiusIterations" -> 100,
	"MidFieldPoint" -> None
};

FindBounce//SyntaxInformation={
	"ArgumentsPattern"->{_,_.,_.,OptionsPattern[]},
	"LocalVariables"->{"Solve",{2,2}}
};

(* This definition should take care of single field case. *)
FindBounce[V_,fields_/;Length[fields]==0,{min1_,min2_},opts:OptionsPattern[]]:=
	FindBounce[V,{fields},{min1,min2},opts];	

(* Definition for a potential defined by a list of points. *)
FindBounce[points_List,opts:OptionsPattern[]]:=(
	If[
		Not@And[
			ArrayQ[N@points,2,(MatchQ[#,_Real]&)],
			MatchQ[Dimensions[points],{x_/;x>=3,2}]
		],
		Message[FindBounce::points];Return[$Failed]
	];
	FindBounce[points,{True},points[[1,{1,-1}]],opts]
);	
	
FindBounce[V_,fields_List,{min1_,min2_},opts:OptionsPattern[]]:=
Module[{Ns(*Number of segments*),a,path,\[Phi]L,ansatzInitialR,b,v,\[Phi],dim,accuracyRadius,
	noFields,VL,d\[Phi]L,point,fieldpoints,maxItePath,maxIteR,R,improvePB,
	rule,improvementPB,pos,l,eL,dV,d2V,\[Phi]l,RM,actionP,action\[Xi],action,vM,aM,bM,posM,
	ddVL,dPath,bottomless,p,dAction,iter=0,potentialPoints=None,switchPath=False,initialR=None},
	
	(*Checks if field variables do not have any values. *)
	If[Not@ArrayQ[fields,1,(Head[#]===Symbol&)],Message[FindBounce::syms];Return[$Failed,Module]];
	
	(* Checking of acceptable option values. If they are wrong function immediately returns $Failed. *)
	bottomless = TrueQ@OptionValue["BottomlessPotential"];
	accuracyRadius = OptionValue["InitialRadiusAccuracyGoal"]/.Except[_Integer?Positive]:>(Message[FindBounce::posint,"InitialRadiusAccuracyGoal"];Return[$Failed,Module]);
	dPath = OptionValue["PathTolerance"]/.Except[_?NonNegative]:>(Message[FindBounce::posreal,"PathTolerance"];Return[$Failed,Module]);
	dAction = OptionValue["ActionTolerance"]/.Except[_?NonNegative]:>(Message[FindBounce::posreal,"ActionTolerance"];Return[$Failed,Module]);
	dim = OptionValue["Dimension"]/.Except[3|4]:>(Message[FindBounce::dim];Return[$Failed,Module]);
	maxIteR = OptionValue["MaxRadiusIterations"]/.Except[_Integer?Positive]:>(Message[FindBounce::posint,"MaxRadiusIterations"];Return[$Failed,Module]);
	maxItePath = OptionValue["MaxPathIterations"]/.Except[_Integer?NonNegative]:>(Message[FindBounce::nonnegint,"MaxPathIterations"];Return[$Failed,Module]);
	point = OptionValue["MidFieldPoint"];
	fieldpoints = OptionValue["FieldPoints"];
	
	noFields = Length[fields];
	If[noFields == 1, maxItePath = 0];
	
	(*Checks if the minima are real.*)
	If[
		Not@And[
			ArrayQ[N@{min1,min2},1|2,(MatchQ[#,_Real]&)],
			(noFields==1&&Length[N@min1]==0)|(noFields>1&&Length[N@min1]==noFields)
		],
		Message[FindBounce::mins];Return[$Failed,Module]
	];
	
	If[
		And[
			Not[IntegerQ[fieldpoints]&&fieldpoints>2],
			Not[ArrayQ[N@fieldpoints,noFields,(MatchQ[#,_Real]&)]&&Length[N@fieldpoints]>2]
		],	
		Message[FindBounce::fieldpts];Return[$Failed,Module]
	];
	
	If[fields[[1]]===True,
		{fieldpoints,potentialPoints} = Transpose@V
	];

	(*InitialValue*)
	{ansatzInitialR,Ns,\[Phi],\[Phi]L,eL,l,dV,d2V,improvePB,path} = 
		InitialValue[V,fields,noFields,min1,point,min2,potentialPoints,
		OptionValue[Gradient],OptionValue[Hessian],dim,
		bottomless,fieldpoints]/.x_/;FailureQ[x]:>Return[$Failed,Module];
		
	(*If the vacua are degenerated*)
	If[
		Head[ansatzInitialR]===DirectedInfinity&&Not[bottomless],
		Message[FindBounce::degeneracy];
		Return[BounceFunction@Association[
				"Action"->Infinity,
				"Bounce"->Function[\[Phi][[1]]],
				"BottomlessPotential"->Missing["NotAvailable"],
				"Coefficients"->Null,
				"Dimension"->dim,
				"PathIterations"->0,
				"Path"->\[Phi],
				"Radii"->{0,Infinity}
			]
		]
	];

	(*Bounce and path deformation*)
	While[iter <= maxItePath||switchPath,
	
		(*Rule*)
		rule = If[
			Length[fields] == 1,
			Table[fields[[1]]->\[Phi][[s]],{s,Ns+1}],
			Table[fields[[i]]->\[Phi][[s,i]],{s,Ns+1},{i,noFields}]
		];
		
		(*Single Field Bounce*)
		{actionP,VL,v,a,b,pos,R,initialR} = 
			If[Not[bottomless],
				SingleFieldBounce[V,potentialPoints,Ns,noFields,\[Phi]L,dim,maxIteR,
					accuracyRadius,ansatzInitialR,initialR,rule,iter,
					switchPath||iter==maxItePath
					]/.(x_/;Not@FreeQ[x,$Failed]:>Return[$Failed,Module])
				,
				BottomlessPotentialBounce[V,potentialPoints,Ns,noFields,\[Phi]L,
					dim,maxIteR,accuracyRadius,ansatzInitialR,initialR,rule,
					iter,fields
					]/.(x_/;Not@FreeQ[x,$Failed]:>Return[$Failed,Module])
			];

		{action\[Xi],ddVL} = SingleFieldBounceImprovement[VL,dV,noFields,rule,Ns,v,a,b,R,\[Phi]L,pos,dim,eL,
			improvePB&&path===None];
				
		(*Transforms \[Phi]L,v,a,b (logitudinal) into \[Phi] (field space) and its bounce parameters.*)
		{v,a,b,\[Phi],action,switchPath} = ParameterInFieldSpace[v,a,b,R,\[Phi],eL,l,\[Phi]L,Ns,noFields,pos,dim,
			bottomless,action,actionP+action\[Xi],dAction,switchPath];
	
		(*Breaks the interations of path deformation*)
		If[switchPath||iter == maxItePath||bottomless, 
			If[noFields> 1&& maxItePath >0, {pos,v,a,b,R} = {posM,vM,aM,bM,RM}];
			p = If[pos>1,pos-1,pos]; 
			Break[]
		];

		(*Multi Field Bounce*)
		{\[Phi],vM,aM,bM,RM,posM,switchPath} = MultiFieldBounce[fields,dV,d2V,Ns,noFields,pos,
			dim,R,\[Phi],v,a,b,\[Phi]L[[-1]],dPath];
		{Ns,\[Phi]L,eL,l} = NewAnsatz[\[Phi],Ns];
		iter++
	];
	
	BounceFunction@Association[
		"Action"->action,
		"Bounce"->piecewiseBounce[{v,a,b,R},{\[Phi][[1]],\[Phi][[-1]]},{dim,pos,Ns,noFields,bottomless}],
		"BottomlessPotential"->If[bottomless,VL[[1]],Missing["NotAvailable"]],
		"Coefficients"->{v,a,b}[[All,p;;-1]],
		"Dimension"->dim,
		"PathIterations"->iter,
		"Path"->\[Phi],
		"Radii"->R[[p;;-1]]
	]

];


(* ::Section::Closed:: *)
(*Visualization: BouncePlot*)


BouncePlot::usage="BouncePlot[bf] plots the content of BounceFunction bf.\n"<>
"BouncePlot[{bf1,bf2}] plots several functions bfi.";

BouncePlot//Options=Options@Plot;

BouncePlot//SyntaxInformation={"ArgumentsPattern"->{_,OptionsPattern[]}};

BouncePlot[bf_BounceFunction,opts:OptionsPattern[]]:= BouncePlot[{bf},opts];
	
BouncePlot[{bf__BounceFunction},opts:OptionsPattern[]]:= Module[
	{bounce,radii,defaultPlotRange,plotRange,rangeMin,rangeMax},
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
(*BottomlessPotential*)


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*BottomlessParameterRvb*)


Rs[4,c1_?NumericQ,a_,b_] := Sqrt[ 1/2 (Sqrt[ c1^2 - 4 a b ] + c1)/a  ];

BottomlessParameterRvb[initialR_?NumericQ,\[Phi]L_,a_,d_,Ns_,backward_,\[Phi]m_]:=
BottomlessParameterRvb[initialR,\[Phi]L,a,d,Ns,backward,\[Phi]m] =
Module[{R,b,v,\[Alpha],v1,b1,x,y,z,Rvb,p=2,\[Phi]0,b4,v4},
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


(* ::Subsection::Closed:: *)
(*FindInitialRadiusB*)


(*Find the solution of eq. 25 or 26.*)
FindInitialRadiusB[d_,VL_,\[Phi]L_,a_,Ns_,maxIteR_,accuracyRadius_,ansatzInitialR_,aRinitial_,\[Phi]m_]:= 
Module[{R,initialR,timeRinitial,ite,findInitialR,\[Lambda],complexR,realR,switch,k,initialR0},

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
	(*Defines the method to use*)
	\[Lambda][initialR_] := \[Lambda][initialR] = Chop[ Sqrt[\[CapitalLambda]B[initialR,d,VL,\[Phi]L,a,Ns,True,\[Phi]m]] ];
	findInitialR[initialR_]:= findInitialR[initialR] = 
		Module[{rw}, 
			rw =rw/.Quiet[FindRoot[Abs[\[Lambda][rw]-1],{rw,initialR}, 
				MaxIterations->1,PrecisionGoal->0,AccuracyGoal->accuracyRadius ],
				{FindRoot::lstol,FindRoot::cvmit}]; 
				  
			Abs[rw]
	]; 
	(*Finds the interval of the solution Sol \[Element] [realR, complexR] or reduces the interval*)
	ite = 0; switch = True; k = 1;
	While[ ite <= maxIteR && switch &&k<10,    
		complexR = Infinity; realR=10^(-1); switch = False;
		If[ Abs[Im[\[Lambda][initialR]]]>10^(-12),
			(*-Overshooting*)
			While[(Abs[Im[\[Lambda][initialR]]]>10^(-12)||(Abs@Im[\[Lambda][initialR]]<=  10^(-12)&&\[Lambda][initialR]>1))&&ite <= maxIteR&&Abs[\[Lambda][initialR]-1]>0.5*10^(-accuracyRadius),
				If[initialR < complexR,
					complexR = initialR;
					initialR = findInitialR[initialR], 
					initialR = 0.9*complexR]; 
			ite++]; 
			,
			(*-Undershooting*)
			(*----else----------*)
			While[ Abs@\[Lambda][initialR]< 1 &&Abs@Im[\[Lambda][initialR]]<= 10^(-12)&&ite <= maxIteR&&Abs[\[Lambda][initialR]-1]>0.5*10^(-accuracyRadius)&&Chop@Re[\[Lambda][initialR]]!=0,
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
		While[ite <= maxIteR &&Abs[\[Lambda][initialR]-1]>0.5*10^(-accuracyRadius)&&Chop[Re[\[Lambda][initialR]]]!=0, 
		
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
	
	 If[ Re[\[Lambda][initialR]-1] >0.5*10^(-1), 
		Message[FindInitialRadius::noSolution];
		Return[$Failed,Module]
	];  
	Clear[findInitialR,\[Lambda]];
	
	Re@initialR  
]; 


(* ::Subsection::Closed:: *)
(*\[ScriptCapitalT]B,\[ScriptCapitalV]B,\[CapitalLambda]B*)


(*Kinetic term from Int[\[Rho]^(D-1)(1/2 d\[CurlyPhi]4^2)]*)
\[ScriptCapitalT]Bs[a4_,\[Rho]_?NumericQ,b4_] := 
	-((4 (4+3 b4^2 a4 \[Rho]^2 (2+b4^2 a4 \[Rho]^2)))/(3 a4 (2+b4^2 a4 \[Rho]^2)^3));

\[ScriptCapitalT]B[R_,Ns_,a_,b_] := 
Module[{\[ScriptCapitalT],VT,p=2,d=4},
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


(* ::Subsection::Closed:: *)
(*BottomlessPotentialBounce*)


BottomlessPotentialBounce::extrema = "Wrong position of the minima.";
BottomlessPotentialBounce::pathDeformation = "The path is deformed irregularly on the potential, try changing number of segments.";
BottomlessPotentialBounce::Rinitial0 = "Trivial solution founded, increase the number of segments or accuracy.";
BottomlessPotentialBounce::nrm = "The potential should be a polynomial of order 4.";

BottomlessPotentialBounce[V_,potentialPoints_,Ns_,noFields_,\[Phi]L_,dim_,maxIteR_,accuracyRadius_,
ansatzInitialR_,aRinitial_,rule_,iter_,fields_]:= 
Module[{a,VL,pos,initialR,R,v,b,T1,V1,\[CurlyPhi],\[Phi]m,cList,\[Lambda],v0},

	cList = CoefficientList[Expand@Normal@Series[V,{fields[[1]],Infinity,4}],fields[[1]]];
	If[Length[cList]===5,
		\[Lambda] = Abs@cList[[5]];
		v0 = cList[[4]]/(4*\[Lambda]);
		,
		Message[BottomlessPotentialBounce::nrm];
		Return[$Failed,Module]
	];
	
	VL = If[potentialPoints===None,Table[V/.rule[[s]],{s,Ns+1}],potentialPoints];	
	
	If[VL[[-1]]>=VL[[-2]],
		If[iter === 0,
			Message[BottomlessPotentialBounce::extrema];
			Return[$Failed,Module],
			Message[BottomlessPotentialBounce::pathDeformation];
			Return[$Failed,Module]
		]
	];
	
	a  = Table[ ( (VL[[s+1]]-VL[[s]])/(\[Phi]L[[s+1]]-\[Phi]L[[s]]) )/8,{s,Ns}]; 
	a[[1]] = \[Lambda];
	\[Phi]m = v0 + \[Phi]L[[-1]];
	initialR = FindInitialRadiusB[dim,VL,\[Phi]L,a,Ns,maxIteR,accuracyRadius,ansatzInitialR,aRinitial,\[Phi]m]//Re;
	
	If[initialR<10^(-5.),
		Message[BottomlessPotentialBounce::Rinitial0];
		Return[$Failed,Module]
	];
	
	{R,v,b} = BottomlessParameterRvb[initialR,\[Phi]L,a,dim,Ns,True,\[Phi]m]//Re;
		
	a = Join[a,{0.}];	
	T1 = \[ScriptCapitalT]B[R,Ns,a,b]; 
	V1 = \[ScriptCapitalV]B[R,v,VL,\[Phi]L,Ns,a,b,\[Phi]m];
	VL[[1]] = VL[[2]]+\[Lambda]*\[Phi]m^4;

	{V1+T1,VL,v,a,b,2,R,initialR}
];


(* ::Chapter::Closed:: *)
(*End Package*)


End[]; (*"`Private`"*)

EndPackage[];
