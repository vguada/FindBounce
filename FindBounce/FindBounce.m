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


FindBounce;
BounceFunction;
BouncePlot;
SingleFieldBounce;
SingleFieldBounceImprovement;
MultiFieldBounce;


(* ::Section::Closed:: *)
(*Begin Private*)


Begin["`Private`"];


(* ::Chapter::Closed:: *)
(*Code*)


(* ::Section::Closed:: *)
(*InitialValue*)


(* ::Subsection::Closed:: *)
(*Segmentation*)


Options[Segmentation] = {"Segments" -> 30,"Method" -> "2H"};

Segmentation//SyntaxInformation={
	"ArgumentsPattern"->{{_,_,_},OptionsPattern[]}
};

Segmentation[\[Phi]N3_,OptionsPattern[]]:=
Module[{Nfv,\[Phi]3,\[Phi],\[Delta]\[Phi],\[Delta]\[Phi]1,\[Delta]\[Phi]2,\[CapitalPhi],np,n1,n2,infP,pos1,pos2,
	\[Psi],\[Beta],p1,p2,p3,p4,p,case1,case2,case3,case4}, 
	\[Phi]3 = N[\[Phi]N3]; 
	Nfv = OptionValue["Segments"]+1; (*Nfv = field values = "number of segement" + 1*)
	If[ Nfv <= 3, 
		
		Return[\[Phi]3], 
		(*Homogeneous_Segmentation*)
		If[OptionValue["Method"] === "H" , 
			\[Delta]\[Phi]= Abs[ (\[Phi]3[[3]] - \[Phi]3[[1]]) ]/(Nfv-1);
			\[CapitalPhi] = Table[ i , { i, \[Phi]3[[1]], \[Phi]3[[3]], \[Delta]\[Phi]} ];
		
		Return[\[CapitalPhi]] 
		];
		
		(*Bi-Homogeneous_Segmentation*)
		If[OptionValue["Method"] === "2H", 
			\[Delta]\[Phi] = Abs[ (\[Phi]3[[3]] - \[Phi]3[[1]]) ]/(Nfv-1);
			n1 = Quotient[Abs[ (\[Phi]3[[2]] - \[Phi]3[[1]]) ], \[Delta]\[Phi] ];
			n2 = IntegerPart[ (Nfv-1) - n1 ];
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


(* ::Input::Initialization:: *)
DerivativePotential::gradient = "\"Gradient\" is not a vector, default value was taken.";
DerivativePotential::hessian = "\"Hessian\" is not a matrix, default value was taken."; 

DerivativePotential[V_,fields_,noFields_,gradient_,hessian_]:= Module[{dV,d2V},
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


(* ::Input::Initialization:: *)
InitialValue::wrongInput = "Wrong \"`1`\".";
InitialValue::dimArray = "The array dimention of min1, min2 and fields are inconsistent.";

InitialValue[V_,fields_,noFields_,min1_,point_,min2_,segments_,Vinitial_,methodSeg_,Path_,gradient_,hessian_,dim_,precision_]:=
Module[{VL,Ns,\[Phi],\[Phi]L,eL,l,rule,\[Phi]3,VL3,L3,\[Phi]L3,eL3,a3,Rinitial,
Length\[Phi]0,dV = None,d2V=None,improvePB=False,c,\[CurlyPhi]0,R0 },
SetPrecision[

If[Path === None, 
	\[Phi] = N[{min1,point,min2}],
	\[Phi] = N[Path]
];
Length\[Phi]0 = Length[\[Phi]];

rule = If[
		noFields == 1,
		Table[fields[[1]]->\[Phi][[s]],{s,Length\[Phi]0}], 
		Table[fields[[i]]->\[Phi][[s,i]],{s,Length\[Phi]0},{i,noFields}]
	];
	
	VL = If[
		Vinitial===None, 
		Table[V/.rule[[s]],{s,Length\[Phi]0}], 
		{Vinitial[[1]],Max[Vinitial],Vinitial[[-1]]}
	];

(*Checks for errors*)
	If[!NumericQ[VL[[1]]] || !NumericQ[VL[[2]]] || !NumericQ[VL[[-1]]],
		Message[InitialValue::wrongInput,"InitialPotential"];Return[$Failed,Module]
	];
	 
	If[Length[\[Phi][[1]]] =!= Length[\[Phi][[2]]] || Length[\[Phi][[2]]] =!= Length[\[Phi][[-1]]],
		Message[InitialValue::dimArray];Return[$Failed,Module]   
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
a3 = Table[(VL3[[s+1]]-VL3[[s]])/(\[Phi]L3[[s+1]]-\[Phi]L3[[s]]) 1/8. ,{s,2}]  ;

(*Estimates the initial radii Rinitial with the N=2 close form solution*)
Rinitial =
 If[ VL3[[1]]!= VL3[[3]] ,
	R0 =1/2(\[Phi]L3[[3]]-\[Phi]L3[[1]])/(Sqrt[a3[[1]] (\[Phi]L3[[2]]-\[Phi]L3[[1]])]-Sqrt[-a3[[2]] (\[Phi]L3[[3]]-\[Phi]L3[[2]])]);
	If[Im@R0>0,
		c = dim/(dim-2) (a3[[2]]-a3[[1]])/a3[[1]] (1 -(a3[[2]]/(a3[[2]]-a3[[1]]))^(1-2/dim));
		\[CurlyPhi]0 =( ( \[Phi]L3[[3]] + c \[Phi]L3[[2]] )/(1+c)  );
		Sqrt[dim/4 (\[Phi]L3[[2]]-\[CurlyPhi]0)/a3[[1]]]
		,
		R0
	]
,
		Infinity
];

	(*Finds the Logitudinal field values \[Phi]L, fields \[Phi], distance l, and direction eL*)
	If[Path === None,
	\[Phi]L = Segmentation[\[Phi]L3, "Segments" -> segments ,"Method" -> methodSeg];
		\[Phi] = Chop@Table[ 
		If[ \[Phi]L[[s]]< \[Phi]L3[[2]],
			eL3[[1]](\[Phi]L[[s]]-L3[[1]])+ \[Phi]3[[2]],
			eL3[[2]](\[Phi]L[[s]]-(L3[[1]] + L3[[2]]))+\[Phi]3[[3]]
		]
	,{s,1,Length[\[Phi]L]} ];
	l = \[Phi]L[[2;;-1]]-\[Phi]L[[1;;-2]];
	,
	l  = Table[Norm[{\[Phi][[s+1]]-\[Phi][[s]]}],{s,Length\[Phi]0-1}];
	\[Phi]L = Table[Sum[l[[s1]],{s1,s-1}],{s,Length\[Phi]0}];
];

If[Vinitial===None,
	If[Not[gradient === False &&noFields==1],
		{dV,d2V} = DerivativePotential[V,fields,noFields,gradient,hessian];
	];
	If[gradient===Automatic&&Ns>3,
		improvePB = True
	];
];
eL = (\[Phi][[2;;-1]]-\[Phi][[1;;-2]])/l;
,precision];

{Rinitial,Length[\[Phi]L]-1,\[Phi],\[Phi]L,eL,l,dV,d2V,improvePB}
];


(* ::Section::Closed:: *)
(*SingleFieldBounce*)


(* ::Subsection::Closed:: *)
(*FindSegment*)


FindSegment[a_,\[Phi]L_,d_,Ns_,precision_]:=Module[{pos = 1,R,estimateRmin,estimateRmax,ps},
	R = BounceParameterRvb[0,\[Phi]L,a,d,Ns,False,pos,precision][[1]];
	
	While[Im[R[[-1]]] == 0 && pos< Ns,
		pos++;
		R = BounceParameterRvb[0,\[Phi]L,a,d,Ns,False,pos,precision][[1]]; 
	];
	
	pos
];


(* ::Subsection::Closed:: *)
(*RInitial*)


(*This function return the second or second last Radii, i.e R[[pos+1]] or R[[-2]] if D = 4*)
RInitial[4,Rinitial_?NumericQ,a_,\[Phi]L_,pos_,backward_]:=
	If[backward,
		Sqrt[Rinitial^2 +  Sqrt[(\[Phi]L[[-2]]- \[Phi]L[[-1]])Rinitial^2 /a[[-2]]]],
		(*else---------*)
		Sqrt[Rinitial^2 + (a[[pos]]Rinitial^2 - Sqrt[a[[pos]]^2 Rinitial^4+4(a[[pos+1]] - 
		a[[pos]])(\[Phi]L[[pos+1]] - \[Phi]L[[pos]])Rinitial^2   ]  )/(2(a[[pos+1]] -a[[pos]])) ] 				        
	 ];
RInitial[3,Rinitial_?NumericQ,a_,\[Phi]L_,pos_,backward_]:= Rinitial;


(* ::Subsection::Closed:: *)
(*BounceParameterRvb*)


(*=========== Rs(D) ====================*) 
(*See eqs. 18-20*)
Rs[4,c1_?NumericQ,a_,b_] := Sqrt[ 1/2 (Sqrt[ c1^2 - 4 a b ] + c1)/a  ];
Rs[3,c1_?NumericQ,a_,b_] := Module[{\[Xi]},  
	\[Xi]= ( Sqrt[ 36 a b^2 -  c1^3 ] - 6 b a^(1/2)  )^(1/3) /a^(1/2); 
	
	If[Re@\[Xi]!=0, 
		1/2 (c1/a/\[Xi] + \[Xi]),
		( Sqrt[ 36 a b^2 -  c1^3 ] + 6 b a^(1/2)  )^(1/3) /a^(1/2)
	] 
];
(*=========== BounceParameterRvb ======================*) 
(*See eqs. 15-16*)
BounceParameterRvb[Rinitial_?NumericQ,\[Phi]L_,a_,d_,Ns_,backward_,pos_,precision_]:= BounceParameterRvb[Rinitial,\[Phi]L,a,d,Ns,backward,pos,precision] =
Module[{R,b,v,\[Alpha],v1,b1,x,y,z,Rvb,position=pos},
SetPrecision[
	(*-------Backward--------------------*)
	If[backward,
		\[Alpha] = Join[a,{0.}]; 
		R = RInitial[d,Rinitial,\[Alpha],\[Phi]L,pos,backward]; 
		b=0.;
		v=\[Phi]L[[-1]];
		Rvb = Reap[
			Sow[b,x];Sow[v,y];Sow[R,z];
			Do[ b +=-(4/d)(\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^(d);Sow[b,x];
				v +=(4/(d-2))(\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^2 ;Sow[v,y];
				R  = Rs[d,(\[Phi]L[[-i-1]]-v) ,\[Alpha][[-i-1]],b];Sow[R,z];
				,{i,1,Ns}
			];
		][[2]];
		
		Chop@Return[Reverse[Rvb,{1,2}]],     
	(*--------Else-Forward---------------*)
		\[Alpha] = Join[{0}, a ]; 
		R = RInitial[d,Rinitial,\[Alpha],\[Phi]L,pos,backward]; 
		b=0.;
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
, precision]; 
];


(* ::Subsection::Closed:: *)
(*FindInitialRadii*)


FindInitialRadii::noSolution = "Large error solving the boundaries conditions. This potential may not have a bounce solution.";

(*Find the solution of eq. 25 or 26.*)
FindInitialRadii[d_,VL_,\[Phi]L_,a_,Ns_,maxIteR_,accuracyRadii_,ansatzRinitial_,aRinitial_,pos_,switchMessage_,precision_]:= 
Module[{R,Rinitial,timeRinitial,ite,RW,\[Lambda],Rcomplex,Rreal,switch,k,Rinitial0,RComplexity},
SetPrecision[
	If[Ns==2&&Not[d==3&&pos==1],
		Return[ansatzRinitial]
	];
	
	(*Picks up the best estimate*)
	If[NumericQ[aRinitial],
		Rinitial = aRinitial
		,
		R = BounceParameterRvb[0.,\[Phi]L,a,d,Ns,False,1,precision][[1]];
		If[ Abs[Im[R[[-1]]]]<0,
			Rinitial = Abs[R[[-2]]],
			Rinitial = Abs[ansatzRinitial] 
		];   
	];
	Rinitial0 = Rinitial;
	
	(*Defines \[Lambda] of eq. 26*)
	\[Lambda][Rinitial_] := \[Lambda][Rinitial] = Sqrt[\[CapitalLambda][Rinitial,d,VL,\[Phi]L,a,Ns,True,pos,precision]];
	(*Discriminates between undershooting and overshooting*)
	RComplexity[Rinitial_] := RComplexity[Rinitial] = Abs@Im@BounceParameterRvb[Rinitial,\[Phi]L,a,d,Ns,True,pos,precision][[1,1]];

	RW[Rinitial_]:= RW[Rinitial] = Module[{rw}, 
		rw =rw/.Quiet[FindRoot[Abs[Re@\[Lambda][rw]-1],{rw,Rinitial}, 
			MaxIterations->1,PrecisionGoal->0,AccuracyGoal->accuracyRadii],
			{FindRoot::lstol,FindRoot::cvmit}];   

		Abs[rw]
	];      
	
	(*Finds the interval of the solution Sol \[Element] [Rreal, Rcomplex] or reduces the interval*)
	ite = 0; 
	switch = True; 
	k = 1;
	While[ ite <= maxIteR && switch,   
		Rcomplex = Infinity;
		Rreal = 0;
		switch = False;
		If[ RComplexity[Rinitial]>0,
			(*Overshooting*)
			While[ RComplexity[Rinitial]>0&&ite <= maxIteR&&Abs[\[Lambda][Rinitial]-1]>.5*10^(-accuracyRadii)&&Chop[Re[\[Lambda][Rinitial]]]!=0.,
				If[Rinitial < Rcomplex,
					Rcomplex = Rinitial;
					Rinitial = RW[Rinitial];
					,
					(*Perturbs 10% down the ansatz since Rinitial>Rcomplex*) 
					Rinitial = .9*Rcomplex
				]; 
				ite++
			]; 
			Rreal = Rinitial;,
			(*Undershooting*)
			(*----else----------*)
			While[ RComplexity[Rinitial] == 0&&ite <= maxIteR&&Abs[\[Lambda][Rinitial]-1]>.5*10^(-accuracyRadii)&&Chop@Re[\[Lambda][Rinitial]]!=0.,
				If[Rinitial > Rreal,
					Rreal = Rinitial;  
					Rinitial = RW[Rinitial],
					(*Perturbs 10% up the ansatz since Rinitial<Rreal*)  
					Rinitial = 1.1*Rreal
				];
				ite++
			]; 
			Rcomplex = Rinitial;    
		];  
		(*One the interval is found, reduces the interval and use bisection method*)	
		While[ite <= maxIteR&&Abs[\[Lambda][Rinitial]-1]>.5*10^(-accuracyRadii)&&Chop@Re@\[Lambda][Rinitial]!=0., 
			
			If[ RComplexity[Rinitial]>0,
				If[ (*Overshooting*)
					Rinitial < Rcomplex,
					Rcomplex = Rinitial;
					Rinitial = RW[Rinitial];
					,
					Rinitial = Abs[Rcomplex+Rreal]/2.;
				];
				,
				If[ (*Undershooting*)
					Rinitial > Rreal, 
					Rreal = Rinitial; 
					Rinitial = RW[Rinitial];
					,
					Rinitial = Abs[Rcomplex+Rreal]/2.;
				];
			];

			ite++   
		];
	   
		If[ Re@Chop@\[Lambda][Rinitial] == 0.,
			k++;
			switch=True;
			Rinitial = (1+k)*Abs[Rinitial0]  
		];   
	];    
	
	If[ ite > maxIteR&&switchMessage,
		Message[FindInitialRadii::cvmit,maxIteR] 
	];
	
	If[ Abs[\[Lambda][Rinitial]-1] >.5*10^(-1)&&switchMessage, 
		Message[FindInitialRadii::noSolution];
		Return[$Failed,Module]
	];  
, precision];	
			
	Re@Rinitial
]; 


(* ::Subsection::Closed:: *)
(*\[ScriptCapitalT],\[ScriptCapitalV],\[CapitalLambda]*)


(*See eqs. 10*)
\[ScriptCapitalT][v_,a_,b_,R_,d_,Ns_,pos_,precision_]:= \[ScriptCapitalT][v,a,b,R,d,Ns,pos,precision] = 
Module[{p,T},
SetPrecision[
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
,precision]
];
 
 (*See eq. 11*)
\[ScriptCapitalV][v_,a_,b_,R_,d_,VL_,\[Phi]L_,Ns_,precision_]:= \[ScriptCapitalV][v,a,b,R,d,VL,\[Phi]L,Ns,precision] = 
SetPrecision[2\[Pi]^(d/2)/Gamma[d/2](Sum[   
	32 a[[i]]^2/(d(d+2))  (R[[i+1]]^(2+d) - R[[i]]^(2+d)) + 
	8 a[[i]]b[[i]]/(d-2) (R[[i+1]]^2 -R[[i]]^2) +
	( VL[[i]]-VL[[-1]]+ 8 a[[i]]( v[[i]] - \[Phi]L[[i]]))(R[[i+1]]^d-R[[i]]^d)/d ,{i,1,Ns}] +
		1/d R[[1]]^d (VL[[1]] - VL[[-1]]) )
,precision];
 
 (*See eq. 26*)
\[CapitalLambda][Rinitial_?NumericQ,d_,VL_,\[Phi]L_,a_,Ns_,backward_,pos_,precision_]:= \[CapitalLambda][Rinitial,d,VL,\[Phi]L,a,Ns,backward,pos,precision] = 
Module[{R,v,b}, 
SetPrecision[
	{R,v,b} = BounceParameterRvb[Rinitial,\[Phi]L,a,d,Ns,backward,pos,precision]; 
	
	(2-d)*\[ScriptCapitalT][v,a,b,R,d,Ns,pos,precision]/(d*\[ScriptCapitalV][v,a,b,R,d,VL,\[Phi]L,Ns,precision])
,precision]
];


(* ::Subsection::Closed:: *)
(*SingleFieldBounce*)


SingleFieldBounce::extrema = "Wrong position of the minima.";
MultiFieldBounce::pathDeformation = "The path is deformed irregularly on the potential. Verifies that the vacuum is a minimum of the potential (not a saddle point) or changes the number of segements.";
SingleFieldBounce::Rinitial0 = "Trivial solution founded, increase the number of segments or accuracy.";

SingleFieldBounce[V_,Vinitial_,Ns_,noFields_,\[Phi]L_,dim_,maxIteR_,accuracyRadii_,ansatzRinitial_,aRinitial_,rule_,iter_,switchMessage_,precision_]:= 
Module[{a,VL,pos,Rinitial,R,v,b,T1,V1},
SetPrecision[	
	If[
		Vinitial===None,
		VL = Table[V/.rule[[s]],{s,Ns+1}]
		,
		If[Vinitial[[1]]<Vinitial[[-1]],
			VL = Vinitial,
			VL = Reverse[Vinitial]
		];
	];
	
	
	If[VL[[1]]>=VL[[2]]||VL[[-1]]>=VL[[-2]],
		If[iter === 0,
			Message[SingleFieldBounce::extrema];Return[$Failed,Module],
			Message[MultiFieldBounce::pathDeformation];Return[$Failed,Module]
		]
	];
	
	a  = Table[ ((VL[[s+1]]-VL[[s]])/(\[Phi]L[[s+1]]-\[Phi]L[[s]]))/8 ,{s,Ns} ]; 
	pos = FindSegment[a,\[Phi]L,dim,Ns,precision];
	Rinitial = FindInitialRadii[dim,VL,\[Phi]L,a,Ns,maxIteR,accuracyRadii,ansatzRinitial,aRinitial,pos,switchMessage,precision]/.x_/;FailureQ[x]:>Return[$Failed,Module];
	{R,v,b} = Re@BounceParameterRvb[Rinitial,\[Phi]L,a,dim,Ns,True,pos,precision];	
	
	If[R[[-2]]<10^(-5),
		Message[SingleFieldBounce::Rinitial0];
		Return[$Failed,Module]
	];
		
	a = Join[a,{0}];	
	T1 = \[ScriptCapitalT][v,a,b,R,dim,Ns,pos,precision]; 
	V1 = \[ScriptCapitalV][v,a,b,R,dim,VL,\[Phi]L,Ns,precision]; 
,precision];

	{V1+T1,VL,v,a,b,pos,R,Rinitial}
];


(* ::Section::Closed:: *)
(*SingleFieldBounceImprovement*)


(* ::Subsection::Closed:: *)
(*Find\[ScriptCapitalI]*)


(*See eqs. 47-48*)
Find\[ScriptCapitalI][v_,\[Phi]L_,a_,b_,Ns_,pos_,R_,ddVL_,d_]:= Find\[ScriptCapitalI][v,\[Phi]L,a,b,Ns,pos,R,ddVL,d]=
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
			Table[   ddVL0[[s+m]] (   ( v0[[s+m]]-\[Phi]L0[[s+m]] )/6 R[[s]]^2 +
			a0[[s+m]]/15*R[[s]]^4 + b0[[s+m]]*R[[s]])   ,{s,pos,Ns+1}]  ],{m,0,1}];
		d\[ScriptCapitalI] = Table[Join[ConstantArray[0,pos-1], 
			Table[   ddVL0[[s+m]] (   ( v0[[s+m]]-\[Phi]L0[[s+m]] )/3 R[[s]] +
			a0[[s+m]] 4/15 R[[s]]^3 + b0[[s+m]] )   ,{s,pos,Ns+1}]  ],{m,0,1}];
	];
	
	{\[ScriptCapitalI],d\[ScriptCapitalI]}   
];


(* ::Subsection::Closed:: *)
(*BounceParameterr\[Beta]\[Nu]*)


(*See eqs. 50-52*)
BounceParameterr\[Beta]\[Nu][rw_?NumericQ,a_,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_] :=BounceParameterr\[Beta]\[Nu][rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos] = 
Module[{r\[Beta]\[Nu]M,r,\[Beta],\[Nu],x,y,z,\[Beta]prev,\[Alpha]0,a0,c0,b0,c},
	a0 = Join[{0},a]; 
	\[Alpha]0 = Join[{0},\[Alpha]];
	b0 = Join[{0},b];
	c0 = Table[2 (b0[[i]] - 4/d a0[[i]] R[[i]]^d),{i,1,Ns+1}];
	r\[Beta]\[Nu]M=Reap[  
		r = 0.;Sow[r,x]; 
		\[Beta] = \[Beta]prev = 0.;Sow[\[Beta],y];
		\[Nu] = (rw c0[[pos]] )R[[pos]]^(2-d)-4/d \[Alpha]0[[pos]] R[[pos]]^(2)-\[ScriptCapitalI][[1,pos]]; Sow[\[Nu],z];
		r = rw; Sow[r,x]; 
		\[Beta] += (4/d (\[Alpha]0[[pos+1]]-\[Alpha]0[[pos]])+4r(a0[[pos+1]]-a0[[pos]]) )   R[[pos]]^d + 
			(d\[ScriptCapitalI][[2,pos]] -  d\[ScriptCapitalI][[1,pos]]  )R[[pos]]^(d-1) /2; Sow[\[Beta],y];
		\[Nu] += -2/(d-2)(\[Beta]-\[Beta]prev)R[[pos]]^(2-d)-4/d (\[Alpha]0[[pos+1]]-\[Alpha]0[[pos]])R[[pos]]^2 - 
			(\[ScriptCapitalI][[2,pos]]-\[ScriptCapitalI][[1,pos]]); Sow[\[Nu],z];
		Do[ \[Beta]prev=\[Beta];   
			r =(2/(d-2)\[Beta]+(\[Nu]+\[ScriptCapitalI][[1,i]]+4/d \[Alpha]0[[i]] R[[i]]^2)R[[i]]^(d-2))/c0[[i]]; Sow[r,x]; 
			\[Beta]+=4/d (\[Alpha]0[[i+1]]-\[Alpha]0[[i]])R[[i]]^d+4r(a0[[i+1]]-a0[[i]])R[[i]]^d + 
				(  d\[ScriptCapitalI][[2,i]] -  d\[ScriptCapitalI][[1,i]]  )R[[i]]^(d-1) /2; Sow[\[Beta],y];
			\[Nu]+=-2/(d-2)(\[Beta]-\[Beta]prev)/R[[i]]^(d-2)-4/d (\[Alpha]0[[i+1]]-\[Alpha]0[[i]])R[[i]]^2 -
				 (\[ScriptCapitalI][[2,i]]-\[ScriptCapitalI][[1,i]]); Sow[\[Nu],z];
			,{i,pos+1,Ns}
		];
		r = ( 2 \[Beta] / R[[Ns+1]]^(d-1)  - 8/d \[Alpha]0[[Ns+1]] R[[Ns+1]] -
			d\[ScriptCapitalI][[1,Ns+1]])/(8 a0[[Ns+1]] R[[Ns+1]]);Sow[r,x];   
		\[Beta] = 0; Sow[\[Beta],y];
		\[Nu] = 0; Sow[\[Nu],z];               
		][[2]];
		
	r\[Beta]\[Nu]M 
];


(* ::Subsection::Closed:: *)
(*FindInitialRadiiImprovement*)


(*Solves Boundary Conditions \[Xi][N-1] = 0*)
FindInitialRadiiImprovement[rw_?NumericQ,a_ ,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_]:= FindInitialRadiiImprovement[rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos]=
Module[{bc,r,\[Beta],\[Nu]},
		{r,\[Beta],\[Nu]} = BounceParameterr\[Beta]\[Nu][rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos] ;
		
		bc = \[Nu][[-2]] + 2/(d-2) \[Beta][[-2]]*R[[-1]]^(2-d)+4/d \[Alpha][[Ns]]*R[[-1]]^2+\[ScriptCapitalI][[1,-1]]   
];


(* ::Subsection::Closed:: *)
(* \[ScriptCapitalT]\[Xi], \[ScriptCapitalV]\[Xi]*)


(* ::Subsubsection::Closed:: *)
(*\[ScriptCapitalT]\[Xi]*)


(*Kinetic term from Int[\[Rho]^(D-1)(1/2 d\[Xi]^2)]*)
\[ScriptCapitalT]\[Xi]4[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] := (*\[ScriptCapitalT]\[Xi]4[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r] = *)
	1/(24 \[Rho]^2) (a \[Rho]^4 (-48 \[Beta]+\[Rho]^4 (2 ddVL v+16 \[Alpha]+a ddVL \[Rho]^2-2 ddVL \[Phi]L))+
	b (-48 \[Beta]+2 \[Rho]^4 (-3 ddVL v-24 \[Alpha]+2 a ddVL \[Rho]^2+3 ddVL \[Phi]L))
	-24 b^2 ddVL \[Rho]^2 Log[\[Rho]])+ (r (4 b-4 a \[Rho]^4)^2)/(8 \[Rho]^2);

\[ScriptCapitalT]\[Xi]3[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] := (*\[ScriptCapitalT]\[Xi]3[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r] = *)
	-((4 b \[Beta])/\[Rho])-2 b^2 ddVL \[Rho]+8/15 a b ddVL \[Rho]^4+32/315 a^2 ddVL \[Rho]^7-1/3 \[Rho]^2 (8 a \[Beta]+b (8 \[Alpha]+ddVL (v-\[Phi]L)))+
	8/45 a \[Rho]^5 (8 \[Alpha]+ddVL (v-\[Phi]L))+(2 r (3 b-4 a \[Rho]^3)^2)/(9 \[Rho]);

\[ScriptCapitalT]\[Xi][d_?NumericQ,a_,R_,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,ddVL_,VL_,\[Phi]L_,Ns_,pos_,r_] := (*\[ScriptCapitalT]\[Xi][d,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r] =*)
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
\[ScriptCapitalV]\[Xi]4[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=(*\[ScriptCapitalV]\[Xi]4[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r] =*) 
	1/48 \[Rho]^2 (5 a^2 ddVL \[Rho]^6+ 16 a (12 \[Beta]+6 \[Nu] \[Rho]^2+\[Rho]^4 (8 \[Alpha]+ddVL (v-\[Phi]L)))+
	24 b (8 \[Alpha]+ddVL (v-\[Phi]L))+6 \[Rho]^2 (16 \[Alpha]+ddVL (v-\[Phi]L)) (v-\[Phi]L))+
	1/2 b ddVL (b+2 a \[Rho]^4) Log[\[Rho]]+8 a b r \[Rho]^2+1/4 r \[Rho]^4 (4 (VL-VN)+32 a^2 \[Rho]^2+32 a (v-\[Phi]L));

\[ScriptCapitalV]\[Xi]3[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=\[ScriptCapitalV]\[Xi]3[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r]= 1/630 \[Rho] (1260 b^2 ddVL+
	210 b \[Rho] (24 \[Alpha]+ddVL (3 v+8 a \[Rho]^2-3 \[Phi]L))+\[Rho] (128 a^2 ddVL \[Rho]^5+336 a (15 \[Beta]+5 \[Nu] \[Rho]+\[Rho]^3 (8 \[Alpha]+ddVL (v-\[Phi]L)))+105 \[Rho] (16 \[Alpha]+
	ddVL (v-\[Phi]L)) (v-\[Phi]L)))+16 a b r \[Rho]^2+1/3 r \[Rho]^3 (3 (VL-VN)+32 a^2 \[Rho]^2+24 a (v-\[Phi]L));

\[ScriptCapitalV]\[Xi][d_?NumericQ,a_,R_,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,ddVL_,VL_,\[Phi]L_,Ns_,pos_,r_] := (*\[ScriptCapitalV]\[Xi][d,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r]=*)
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


SingleFieldBounceImprovement::dVFailed = "The first derivative of the Potential is not well define at some field value. \"Gradient\"->None was taken.";

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


ParameterInFieldSpace[vs_,as_,bs_,R_,\[Phi]_,eL_,l_,\[Phi]L_,VL_,Ns_,noFields_,pos_,dim_]:=
Module[{v,a,b},
	v = Table[\[Phi][[s+1]]+eL[[s]]*(vs[[s]]-(l[[s]]+\[Phi]L[[s]])),{s,Ns}];
	a = Table[eL[[s]]*as[[s]],{s,Ns}];
	b = Table[eL[[s]]*bs[[s]],{s,Ns}]; 
	
	{v,a,b}
];


(* ::Section::Closed:: *)
(*MultiFieldBounce*)


MultiFieldBounce[fields_,dV_,d2V_,Ns_,noFields_,pos_,d_,R0_,\[CapitalPhi]0_,\[ScriptV]0_,\[ScriptA]0_,\[ScriptB]0_,precision_]:=
Module[{\[Nu],\[Beta],rI,a,\[Zeta]t,R,\[Zeta]ts,\[Phi]M,\[Nu]\[Beta],x,y,d\[CurlyPhi],rF,rN,DV,D2V,
	\[Xi]Mc,M,c,\[Nu]0,\[Beta]0,\[Nu]\[Xi]p,\[Nu]\[Xi]m,\[Beta]\[Xi]p,\[Beta]\[Xi]m,fLowT,fD,fD1,frI,n1,rules,p,rulesQuartic, 
	\[Phi]0 = \[CapitalPhi]0, v0 = \[ScriptV]0, a0 = \[ScriptA]0, b0 = \[ScriptB]0},
		
SetPrecision[
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
		\[Beta]\[Xi]p[s_]:= \[Beta]\[Xi]p[s]= (D2V[[1+s]] (R[[s]]^d) )/(4 d);  (*\[Zeta]ts[1+s]*)
		\[Nu]\[Xi]m[s_]:= \[Nu]\[Xi]m[s]= (D2V[[s-1]] (R[[s]]^2) )/(4 (-2+d));  (*\[Zeta]ts[-1+s]*)
		\[Beta]\[Xi]m[s_]:= \[Beta]\[Xi]m[s]=-((D2V[[s-1]] (R[[s]]^d) )/(4 d));  (*\[Zeta]ts[-1+s]*)
	fLowT[s_,j_]:= (2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]m[j]+\[Beta]\[Xi]p[j-2])+\[Nu]\[Xi]m[j]+\[Nu]\[Xi]p[j-2];(*[Eq[s],\[Xi][j-1]]*) (* 2 \[LessEqual] j \[LessEqual] s*)

		\[Nu]\[Xi]p[p-1]= If[pos>1,IdentityMatrix[noFields],-((D2V[[p]]*(R[[p]]^2) )/(4 (-2+d))) ];
		\[Beta]\[Xi]p[p-1]= If[pos>1,ConstantArray[0,{noFields,noFields}],(D2V[[p]]*R[[p]]^d)/(4 d)];
	fD[s_] := (D2V[[s]]*(R[[1+s]]^2) )/(4 d) +(2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]p[s-1])+ \[Nu]\[Xi]p[s-1];(*[Eq[s],\[Xi][s]]*)

		\[Nu]\[Xi]p[p] = If[pos>1,ConstantArray[0,{noFields,noFields}],-((D2V[[1+p]] R[[p]]^2)/(4 (-2+d)))];
		\[Beta]\[Xi]p[p] = If[pos>1,ConstantArray[0,{noFields,noFields}],(D2V[[1+p]] (R[[p]]^d) )/(4 d)];
	fD1[s_]:= (-IdentityMatrix[noFields]+(D2V[[1+s]] R[[1+s]]^2)/(4 d))+(2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]p[s])+\[Nu]\[Xi]p[s] ;(*[Eq[s],\[Xi][s+1]]*)
	frI[s_]:= ((2 (R[[1+s]]^(2-d)) )/(-2+d) (4 a0[[p]] R[[p]]^d)-(8  a0[[p]] R[[p]]^2)/(-2+d)  )IdentityMatrix[noFields]; 
	
	n1 = If[pos>1,0,1];
	M = Module[{n,m},
		n = (Ns+1-p)+1+n1;
		m = SparseArray[{{n,n}->0}];
		Do[m[[s+n1,s+n1]]=fD[s+p-1],{s,1,n-1-n1}];
		Do[m[[s+n1,s+1+n1]]=fD1[s+p-1],{s,1,n-1-n1}];
		m[[n,n]]=IdentityMatrix[noFields];
		If[pos==1,
			m[[1,2]]=IdentityMatrix[noFields];
			Do[m[[s+n1,1]]=frI[s],{s,1,n-1-n1}];   
		];
		Do[m[[s+n1,j-1+n1]]=fLowT[s+p-1,j+p-1]; ,{s,2,n-n1-1},{j,2,s}];
		
		ArrayFlatten[m]
	];
	
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
		If[ pos>1, Do[Sow[ConstantArray[0,noFields],x];Sow[ConstantArray[0,noFields],y];,{s,1,p-1}];];
		Sow[\[Nu],x];Sow[\[Beta],y];
		Do[ \[Nu]+= -4/(d-2)(a[[s+1]]-a[[s]])R[[s+1]]^2-1/(d-2) (  d\[CurlyPhi][[s,2]]-d\[CurlyPhi][[s,1]] )R[[s+1]] ;
			\[Beta]+=   4/d(a[[s+1]]-a[[s]])R[[s+1]]^d+1/2 (  d\[CurlyPhi][[s,2]]-d\[CurlyPhi][[s,1]] )R[[s+1]]^(d-1);
			Sow[\[Nu],x];Sow[\[Beta],y]; ,{s,p,Ns-1} ];  
	][[2]];
	
	R[[p]] += R[[p]]*rI; 
	R[[Ns+1]] += R[[Ns+1]]*rF; 
	{\[Phi]0[[1]],\[Zeta]ts[[1]]} = {\[CapitalPhi]0[[1]],ConstantArray[0,noFields]};
,precision];

	{\[Phi]0+\[Zeta]ts,v0+\[Nu]\[Beta][[1]],a0+a,b0+\[Nu]\[Beta][[2]],R0,pos}   	
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
	(* Set standard image size *)
	ImageSize -> Dynamic[{Automatic,3.5*CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]
];


BounceFunction::usage="BounceFunction object represents results from FindBounce function.";

BounceFunction[asc_?AssociationQ]["Properties"]:=Sort@Keys[asc];

(* A message about missing property could be issued if neccesary (three argument Lookup).  *)
BounceFunction[asc_?AssociationQ][property_]:=Lookup[asc,property];
	
(* Nice styling of output, see https://mathematica.stackexchange.com/questions/77658 *)
BounceFunction/:MakeBoxes[obj:BounceFunction[asc_?AssociationQ],form:(StandardForm|TraditionalForm)]:=Module[
	{above,below,icon},
	
	 (* column *)
	above = {
		BoxForm`SummaryItem[{"Action: ", obj["Action"]}],
		BoxForm`SummaryItem[{"Dimension: ", obj["Dimension"]}]
	};
	below = {
		BoxForm`SummaryItem[{"Domain: ", obj["Domain"]}],		
		BoxForm`SummaryItem[{"InitialSegment: ", obj["InitialSegment"]}],
		BoxForm`SummaryItem[{"Iterations: ", obj["IterationsPath"]}],
		BoxForm`SummaryItem[{"Segments: ", obj["Segments"]}]
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


(* Returns pure Function - could be something else in the future. *)
piecewiseBounce[{v_,a_,b_,R_},{min1_,min2_},{dim_,pos_,Ns_,noFields_}]:=
Module[{\[CurlyPhi]0,MultiFieldPiecewise},

	If[
		pos>1,
		(*Case A*)
		\[CurlyPhi]0[\[Rho]_] := v[[pos-1]]+ 4/dim*a[[pos-1]]*\[Rho]^2
		,
		(*Case B*)
		\[CurlyPhi]0[\[Rho]_] := min1 
	];
	
	If[noFields==1,
		Return[Evaluate@Piecewise[
			Join[
				{{\[CurlyPhi]0[#],#<R[[pos]]}},
				Table[{
					v[[s]]+ 4/dim*a[[s]]*#^2+2/(dim-2)*b[[s]]/#^(dim-2),
					R[[s]]<=#<R[[s+1]]},
					{s,pos,Ns}
				]
			],
			min2 (* default value of Piecewise *)
			]&
		]
	];
	
	Table[Piecewise[
			Join[
				{{\[CurlyPhi]0[#][[i]],#<R[[pos]]}},
				Table[{
					v[[s,i]]+ 4/dim*a[[s,i]]*#^2+2/(dim-2)*b[[s,i]]/#^(dim-2),
					R[[s]]<=#<R[[s+1]]},
					{s,pos,Ns}
				]
			],
			min2[[i]] (* default value of Piecewise *)
			]
	,{i,noFields}]&		
	
];


(* ::Subsection::Closed:: *)
(*FindBounce*)


FindBounce::usage = "FindBounce[potential,fields,{min1, min2}] computes false vacuum decay in potential with multiple scalar fields.";
FindBounce::syms = "Field symbols should not have any value.";
FindBounce::dim = "Only supported \"Dimension\"s are 3 and 4, default value was taken.";
FindBounce::iter = "Maximum number of iterations should be a positive integer, default value `1` was taken.";
FindBounce::optionValue = "Wrong \"`1`\", default value `2` was taken.";
FindBounce::degeneracy = "Not vacuum decay, the vacua are degenerated.";

Options[FindBounce] = {
	"AccuracyInitialRadii" -> 10,
	"AccuracyPath" -> 10,
	"Dimension" -> 4,
	Gradient -> Automatic,
	Hessian -> None,
	"InitialLocalMaximum" -> Automatic,
	"InitialRadii" -> None,
	"InitialPath" -> None,
	"InitialPotential" -> None,
	"MethodSegmentation" -> "2H",
	"MaxIterationsPath" -> 3,
	"MaxIterationsRadii" -> 100,
	"Segments" -> 30,
	"SetPrecision" -> 20
};

FindBounce//SyntaxInformation={
	"ArgumentsPattern"->{_,_,_,OptionsPattern[]},
	"LocalVariables"->{"Solve",{2,2}}
};

(* This definition should take care of single field case. *)
FindBounce[V_,fields_/;Length[fields]==0,{min1_,min2_},opts:OptionsPattern[]]:=
	FindBounce[V,{fields},{min1,min2},opts];
	
FindBounce[V_,fields_List,{min1_,min2_},opts:OptionsPattern[]]:=
Module[{a,Path,\[Phi]L,ansatzRinitial,b,v,\[Phi],Ns,dim,Rinitial,accuracyRadii,accuracyPath,
	noFields,VL,d\[Phi]L,point,itePath,maxIteR,R,methodSeg,improvePB,Vinitial,
	rule,improvementPB,pos,l,eL,dV,d2V,\[Phi]l,RM,Action,Action\[Xi],vM,aM,bM,posM,
	ddVL,precision,iter=0},
	(*OptionValues*)
	Rinitial = N[OptionValue["InitialRadii"]]/.{Except[_Real?Positive|None]:>(Message[FindBounce::optionValue,"InitialRadii",None];None)};
	Path = OptionValue["InitialPath"];
	Vinitial = OptionValue["InitialPotential"];
	accuracyRadii = OptionValue["AccuracyInitialRadii"]/.{Except[_Integer?Positive]:>(Message[FindBounce::optionValue,"AccuracyInitialRadii",10];10)};
	accuracyPath = OptionValue["AccuracyPath"]/.{Except[_Integer?Positive]:>(Message[FindBounce::optionValue,"AccuracyPath",10];10)};
	dim = OptionValue["Dimension"]/.{Except[3|4]:>(Message[FindBounce::dim];4)};
	Ns = OptionValue["Segments"]/.x_/;(x<2||x \[NotElement] Integers):>(Message[FindBounce::optionValue,"Segments",30];30);
	maxIteR = OptionValue["MaxIterationsRadii"]/.{Except[_Integer?Positive]:>(Message[FindBounce::optionValue,"MaxIterationsRadii",100];100)};
	methodSeg = OptionValue["MethodSegmentation"]/.{Except["H"|"2H"]:>(Message[FindBounce::optionValue,"MethodSegmentation","\"H\""];"H")};
	noFields = Length[fields];
	point = OptionValue["InitialLocalMaximum"]/.Automatic->(min1+min2)/2;
	itePath = OptionValue["MaxIterationsPath"]/.Except[_Integer?NonNegative]:>(Message[FindBounce::iter,3];3);
	precision = OptionValue["SetPrecision"]/.{Except[_Integer?Positive|MachinePrecision]:>(Message[FindBounce::optionValue,"SetPrecision",20];20)};
	If[noFields == 1, itePath = 0];
		
	(*Checks if field variables do not have any values.*)
	If[
		Not@ArrayQ[fields,1,(Head[#]===Symbol&)],
		Message[FindBounce::syms];Return[$Failed,Module]
	];
	(*Checks if InitialPath is valid*)
	If[
		Not@ArrayQ[Path,noFields]&&Path=!=None,
		Message[FindBounce::optionValue,"InitialPath",{min1,point,min2}]
	];
	(*Checks if InitialPotential is valid*)
	If[
		Not@ArrayQ[Vinitial,1]&&Vinitial=!=None,
		Message[FindBounce::optionValue,"InitialPotential",""]
	];

	(*InitialValue*)
	{ansatzRinitial,Ns,\[Phi],\[Phi]L,eL,l,dV,d2V,improvePB} = InitialValue[V,fields,noFields,min1,point,min2,Ns,Vinitial,
		methodSeg,Path,OptionValue[Gradient],OptionValue[Hessian],dim,precision]/.x_/;FailureQ[x]:>Return[$Failed,Module];
	
	(*If the vacua are degenerated*)
	If[
		Head[ansatzRinitial]===DirectedInfinity,
		Message[FindBounce::degeneracy];
		Return[BounceFunction@Association[
				"Action"->Infinity,
				"Bounce"->Function[\[Phi][[1]]],
				"Coefficients"->Null,
				"Dimension"->dim,
				"Domain"->{0,Infinity},
				"InitialSegment"->1,
				"IterationsPath"->Missing["NotAvailable"],
				"Segments"->Ns,
				"Path"->\[Phi],
				"PathLongitudinal"->\[Phi]L,
				"Potential"->VL,
				"PotentialSecondDerivative"->Missing["NotAvailable"],
				"Radii"->{0,Infinity}
			]
		]
	];

	(*Bounce and Path deformation*)
	While[iter <= itePath,
	
		(*Rule*)
		rule = If[
			Length[fields] == 1,
			Table[fields[[1]]->\[Phi][[s]],{s,Ns+1}],
			Table[fields[[i]]->\[Phi][[s,i]],{s,Ns+1},{i,noFields}]
		];
		
		(*Single Field Bounce*)
		{Action,VL,v,a,b,pos,R,Rinitial} = SingleFieldBounce[V,Vinitial,Ns,noFields,\[Phi]L,dim,maxIteR,
				accuracyRadii,ansatzRinitial,Rinitial,rule,iter,iter == itePath,precision]/.x_/;FailureQ[x]:>Return[$Failed,Module];
		{Action\[Xi],ddVL} = SingleFieldBounceImprovement[VL,dV,noFields,rule,Ns,v,a,b,R,\[Phi]L,pos,dim,eL,improvePB&&iter==itePath&&Path===None];
	
		(*Transforms \[Phi]L,v,a,b (logitudinal) into \[Phi] (field space) and its bounce parameters.*)
		{v,a,b} = ParameterInFieldSpace[v,a,b,R,\[Phi],eL,l,\[Phi]L,VL,Ns,noFields,pos,dim];
	
		(*Breaks the interations of path deformation*)
		If[iter == itePath, 
			If[noFields>1&&itePath>0, {pos,v,a,b,R} = {posM,vM,aM,bM,RM}]; 
			Break[]
		];

		(*Multi Field Bounce*)
		{\[Phi],vM,aM,bM,RM,posM} = MultiFieldBounce[fields,dV,d2V,Ns,noFields,pos,dim,R,\[Phi],v,a,b,precision];
		{Ns,\[Phi]L,eL,l} = NewAnsatz[\[Phi],Ns];
		ansatzRinitial = Rinitial;
		iter++
	];

	BounceFunction@Association[
		"Action"-> SetPrecision[Action+Action\[Xi],MachinePrecision],
		"Bounce"->piecewiseBounce[{v,a,b,R},{\[Phi][[1]],\[Phi][[-1]]},{dim,pos,Ns,noFields}],
		"Coefficients"->{v,a,b},
		"Dimension"->dim,
		"Domain"->SetPrecision[{If[pos>1,0,R[[pos]]],R[[-1]]},MachinePrecision],
		"InitialSegment"->pos,
		"IterationsPath"->iter,
		"Segments"->Ns,
		"Path"->\[Phi],
		"PathLongitudinal"->\[Phi]L,
		"Potential"->VL,
		"PotentialSecondDerivative"->ddVL,
		"Radii"->R
	]
];


(* ::Section::Closed:: *)
(*Visualization: BouncePlot*)


BouncePlot::usage="BouncePlot[bf] plots the content of BounceFunction bf.";
BouncePlot::bfs="Argument should be a list of function(s).";

BouncePlot//Options=Options@Plot;

BouncePlot//SyntaxInformation={"ArgumentsPattern"->{_,OptionsPattern[]}};

BouncePlot[bf__BounceFunction,opts:OptionsPattern[]]:=
	BouncePlot[{bf},opts];
	
BouncePlot[bf_List,opts:OptionsPattern[]]:= Module[
	{bounce,R,markers,plotRange},
	
	(* Check if bf is a BounceFunction*)
	If[
		Not@ArrayQ[bf,1,(Head[#]===BounceFunction&)],
		Message[BouncePlot::bfs];Return[$Failed,Module]
	];
	
	bounce = #["Bounce"]&/@bf;
	(* This helps to draw discrete radii R. *)
	R = #["Radii"]&/@bf;
	
	plotRange = Clip[
			MinMax[R,Scaled[0.25]],
			{0,10^100}
	];	
	
	Plot[
		Evaluate[#[\[Rho]]&/@bounce],
		{\[Rho],Sequence@@plotRange},
		opts,
		(* Default options come here. Keep them as short as possible, for flexibiltiy.
		Explicitly given options (above) have precedence. *)
		Frame->True,
		FrameLabel->{"\[Rho]","\[CurlyPhi](\[Rho])"},
		LabelStyle->Directive[Black,FontSize->17, FontFamily->"Times New Roman",FontSlant->Plain],
		GridLines->Automatic
	]
];


(* ::Chapter::Closed:: *)
(*End Package*)


End[]; (*"`Private`"*)

EndPackage[];
