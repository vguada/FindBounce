(* ::Package:: *)

(* ::Chapter::Closed:: *)
(*Header Comments*)


(*False Vacuum Decay with Polygonal Bounce*)

(* :Title: FindBounce *)
(* :Context: Tunneling, First order first transitions, bubble nucleation*)
(* :Author: Victor Guada *)
(* :Summary: Compues decay of the false vacuum in models with multiple scalar*)
(* :Copyright: Victor Guada, Miha Nemev\[SHacek]ek and Alessio Maiezza, 2019 *)

(* This program is free software...*)


(* ::Chapter:: *)
(*BeginPackage*)


BeginPackage["FindBounce`"];


(* ::Section::Closed:: *)
(*Available Functions*)


FindBounce::usage = "FindBounce[ V[{phi1, phi2,\[Ellipsis]}],{phi1, phi2,\[Ellipsis]}, {min1, min2} ]
	computes false vacuum decay in potential with multiple scalar fields.";
BounceAction::usage = "";
BounceDimension::usage = "";
PlotBounce::usage = "";
PlotBounce2D::usage = "";
PlotBounce2DProjection::usage = "";
PlotBounce3D::usage = "";
ListPointBounce::usage = "";
Segmentation;
FindSegment;
NewAnsatz;
Rvb;
FindRw;
Find\[ScriptCapitalI];
r\[Beta]\[Nu];
Findrw;
\[Phi]vabRs;
PathDeformation;
\[ScriptCapitalT];
\[ScriptCapitalV];
\[ScriptCapitalT]\[Xi];
\[ScriptCapitalV]\[Xi];


(* ::Section::Closed:: *)
(*Options*)


Options[FindBounce] = {"AnsatzRadii" -> None,
				"AnsatzPath" -> None,
				"AnsatzV1" -> None,
				"AccuracyBounce" -> 6,
				"AccuracyPathDeformation" -> 10,
				Gradient-> None,
				Hessian-> None,
				"Dimension" -> 4,
				"ForwardBackward" -> "backward",
				"NumberSegments" -> 29,
				"ImprovementPolygonalBounce" ->True,
				"InitialFieldValue"  -> Null,
				"MaxIterationsR" -> 100,
				"MaxIterationsPathDeformation" -> 3,
				"MethodBounce" -> "DerrickFindRoot",
				"MethodSegmentation" -> "HS"};
Options[Segmentation] = {"NumberFieldValues" -> 30,
				"Method" -> "HS"};


(* ::Section::Closed:: *)
(*Messages*)


(*"Dimension"::usage = "Dimension: the dimension of the spacetime; defaul imput is 4.";
"NumberSegments"::usage = "NumberSegements: numer of segments, defult value is 29";
"AnsatzRadii"::usage = "AnsatzRadii: gives a estimation by hand.";
"MaxIterationsR"::usage = "MaxIterationsR: maximum number of iterations to get R";
"MaxIterationsPathDeformation"::usage = "MaxIterationsPathDeformation: number of iterations in path deformation (defaul value 2).";
"ForwardBackward"::usage = "ForwardBackward: set Polygonal Bounce either forward or backward.";
"InitialFieldValue"::usage = "InitialFieldValue: initial field Value to start the trayectory in N=3.";
"AccuracyBounce"::usage = "AccuracyBounce: accuracy in FindRw";
(*========= Segmentation ==========*)	
"NumberFieldValues"::usage = "NumberFieldValues: number of field values(defaul imput is 200).";
"Method"::usage = "Method: specify what method should be used {HS,biHS,HSPlus}, where H:Homogeneous, S:Segementation, bi: Double HS splitted in the Saddle Point and Plus: additional field values close to the extrema";*)
(*========== Comments ==========*)
FindBounce::Error = "Wrong input, PB has aborted.";
FindBounce::ErrorExtrema = "Wrong position of the minima, PB has aborted.";
FindBounce::ErrorPathDeformation = "The path is deformed irregularly on the potential, try changing number of segments.";
FindBounce::ErrorIte = "Wrong number of interation, PB has aborted.";
FindBounce::AnsatzPath = "Wrong AnsatzPath, PB has aborted.";
Segmentation::Error = "Warning: wrong number of field values.";
AnsatzN3::Degeneracy = "There is not tunneling decay since the vacua are degenerated.";
AnsatzN3::Error = "Wrong input, PB has aborted.";
AnsatzN3::ErrorPotential = "Wrong input Potential, PB has aborted.";
FindRw::ErrorMethod = "Wrong method, PB has aborted.";
FindRw::ErrorTolerance = "Failed to converge to the requested accuracy";


(* ::Section::Closed:: *)
(*Begin Private*)


Begin["`Private`"];


(* ::Chapter:: *)
(*Code*)


(* ::Section::Closed:: *)
(*Segmentation*)


Segmentation[\[Phi]N3_,OptionsPattern[]]:=
Block[{Nfv,\[Phi]3,\[Phi],\[Delta]\[Phi],\[Delta]\[Phi]1,\[Delta]\[Phi]2,\[CapitalPhi],np,n1,n2,infP,pos1,pos2,
	\[Psi],\[Beta],p1,p2,p3,p4,p,case1,case2,case3,case4}, 
\[Phi]3 = N[\[Phi]N3]; Nfv = OptionValue["NumberFieldValues"]; 
	If[ Nfv < 2 , Message[Segmentation::Error] ]; 
	If[ Nfv <= 3, Return[\[Phi]3], 
(*========= Homogeneous_Segmentation ===========================*)
	If[OptionValue["Method"] === "HS" , 
\[Delta]\[Phi]= Abs[ (\[Phi]3[[3]] - \[Phi]3[[1]]) ]/(Nfv-1);
\[CapitalPhi] = Table[ i , { i, \[Phi]3[[1]], \[Phi]3[[3]], \[Delta]\[Phi]} ];
	Return[ \[CapitalPhi] ] ];
(*========= Bi-Homogeneous_Segmentation ==========================*)
	If[OptionValue["Method"] === "biHS", 
\[Delta]\[Phi] = Abs[ (\[Phi]3[[3]] - \[Phi]3[[1]]) ]/(Nfv-1);
n1 = Quotient[Abs[ (\[Phi]3[[2]] - \[Phi]3[[1]]) ], \[Delta]\[Phi] ];
n2 = IntegerPart[ (Nfv-1) - n1 ];
\[Delta]\[Phi]1 = Abs[ (\[Phi]3[[2]] - \[Phi]3[[1]]) ]/n1;
\[Delta]\[Phi]2 = Abs[ (\[Phi]3[[3]] - \[Phi]3[[2]]) ]/n2;
\[CapitalPhi] = Join[   Table[ i , { i, \[Phi]3[[1]]  , \[Phi]3[[2]] - \[Delta]\[Phi]1 ,\[Delta]\[Phi]1  }  ],  
			Table[ i , { i, \[Phi]3[[2]], \[Phi]3[[3]] ,\[Delta]\[Phi]2  } ]    ];
	Return[ \[CapitalPhi] ] ];
(*======= Homogeneous_Segmentation_plus_Additional_Segmentation_in_Between ============*)
	If[OptionValue["Method"] === "HSPlus",
np = Quotient[Nfv,20];
Nfv = Nfv-1. - 4.(np+1);
\[Delta]\[Phi] = Abs[ (\[Phi]3[[3]] - \[Phi]3[[1]]) ]/(Nfv);
n1 = Quotient[Abs[ (\[Phi]3[[2]] - \[Phi]3[[1]]) ], \[Delta]\[Phi] ];
n2 = IntegerPart[ (Nfv) - n1];
\[Delta]\[Phi]1 = Abs[ (\[Phi]3[[1]] - \[Phi]3[[2]]) ]/n1;
\[Delta]\[Phi]2 = Abs[ (\[Phi]3[[2]] - \[Phi]3[[3]]) ]/n2;
\[CapitalPhi] = Join[ Table[ i , { i, \[Phi]3[[1]]         , \[Phi]3[[2]] - \[Delta]\[Phi]1 ,\[Delta]\[Phi]1  }  ],
          Table[ i , { i, \[Phi]3[[2]]         , \[Phi]3[[3]]       ,\[Delta]\[Phi]2  }  ],
          Table[ i , { i, \[Phi]3[[1]] + \[Delta]\[Phi]1/2 , \[Phi]3[[1]] + \[Delta]\[Phi]1*np + \[Delta]\[Phi]1/2 , \[Delta]\[Phi]1 }  ],
          Table[ i , { i, \[Phi]3[[2]] + \[Delta]\[Phi]2/2 , \[Phi]3[[2]] + \[Delta]\[Phi]2*np + \[Delta]\[Phi]2/2, \[Delta]\[Phi]2 }  ],
          Table[ i , { i, \[Phi]3[[2]] - \[Delta]\[Phi]1/2 - \[Delta]\[Phi]1*np  , \[Phi]3[[2]] - \[Delta]\[Phi]1/2, \[Delta]\[Phi]1 }  ],
          Table[ i , { i, \[Phi]3[[3]] - \[Delta]\[Phi]2/2 - \[Delta]\[Phi]2*np  , \[Phi]3[[3]] - \[Delta]\[Phi]2/2, \[Delta]\[Phi]2 }  ] ]//Sort//DeleteDuplicates;
	Return[\[CapitalPhi]]     ]        
	]    ];


(* ::Section::Closed:: *)
(*FindSegment*)


FindSegment[a_,\[Phi]L_,d_,Ns_]:=Block[{pos,R,estimateRmin,estimateRmax,ps},
pos = 1;
R = Rvb[0,\[Phi]L,a,d,Ns,"forward",pos][[1]]//Chop;
	While[Im[R[[-1]]]== 0.&&pos< Ns,
		pos++;
		R = Rvb[0,\[Phi]L,a,d,Ns,"forward",pos][[1]]//Chop; 
		];
Return[pos] ];


(* ::Section::Closed:: *)
(*AnsatzN3*)


AnsatzN3[V_,\[CurlyPhi]_,\[Phi]N3_,N\[CurlyPhi]_,Nfv_,aV1_,methodSeg_]:=Block[{VL3,a3,\[Phi]3,\[Phi]L3,L3,d\[Phi]L3,Rw,R30,\[CapitalDelta]2,\[Phi],\[Phi]L,Vext,\[Psi],\[Epsilon],S1,R,l,eL,rule3},
(*============\[Equal] VL3_\[Phi]3 ==============*)
\[Phi]3[2]= \[Phi]N3[[2]];
rule3 = If[Length[\[CurlyPhi]] == 1,Table[\[CurlyPhi][[1]]->\[Phi]N3[[\[Alpha]]],{\[Alpha],1,3}],If[Length[\[CurlyPhi]] == 0,Table[\[CurlyPhi] ->\[Phi]N3[[\[Alpha]]],{\[Alpha],1,3}], Table[\[CurlyPhi][[i]]->\[Phi]N3[[\[Alpha],i]],{\[Alpha],1,3},{i,N\[CurlyPhi]}]]];
Vext = If[aV1===None, Table[V/.rule3[[\[Alpha]]],{\[Alpha],1,3}], {aV1[[1]],Max[aV1],aV1[[-1]]}];	
{VL3[1], VL3[3], VL3[2]} = Vext //Sort;
{\[Phi]3[1], \[Phi]3[3]} = If[ VL3[1] == Vext[[1]] ,{\[Phi]N3[[1]], \[Phi]N3[[3]]},{\[Phi]N3[[3]], \[Phi]N3[[1]]} ] ;
(*============ Errors Test ============*)
	If[!NumericQ[VL3[1]] || !NumericQ[VL3[2]] || !NumericQ[VL3[3]],Message[AnsatzN3::ErrorPotential];Abort[]];
	If[Length[\[Phi]3[1]] =!= Length[\[Phi]3[2]] || Length[\[Phi]3[2]] =!= Length[\[Phi]3[3]] , Message[AnsatzN3::Error];Abort[]];
	If[VL3[1] == VL3[3], Message[AnsatzN3::Degeneracy];Abort[];]; 
	If[N\[CurlyPhi] < 1, Message[AnsatzN3::Error];Abort[]; ];
(*========== L3_\[Phi]L3_d\[Phi]L3_a3 =============*)
Do[L3[\[Alpha]] = Norm[\[Phi]3[\[Alpha]+1]-\[Phi]3[\[Alpha]]],{\[Alpha],1,2}];
Do[\[Phi]L3[\[Alpha]] =Sum[L3[i],{i,1,\[Alpha]-1}],{\[Alpha],1,3}];
Do[d\[Phi]L3[\[Alpha]]=(\[Phi]3[\[Alpha]+1] - \[Phi]3[\[Alpha]])/L3[\[Alpha]],{\[Alpha],1,2}];	
Do[a3[\[Alpha]]  = (VL3[\[Alpha]+1]-VL3[\[Alpha]])/(\[Phi]L3[\[Alpha]+1]-\[Phi]L3[\[Alpha]]) 1/8. ,{\[Alpha],1,2}]  ;   
(*========== Segmentation_\[Phi][0] ==============*)
\[Phi]L = Segmentation[Table[\[Phi]L3[\[Alpha]],{\[Alpha],1,3}], "NumberFieldValues" -> Nfv ,"Method" -> methodSeg];
\[Phi]  = Table[ If[ \[Phi]L[[\[Alpha]]] < \[Phi]L3[2], d\[Phi]L3[1](\[Phi]L[[\[Alpha]]]-L3[1])+ \[Phi]3[2],
		d\[Phi]L3[2](\[Phi]L[[\[Alpha]]]-(L3[1] + L3[2]))+\[Phi]3[3]], {\[Alpha],1,Length[\[Phi]L]} ]//Chop;
l  = \[Phi]L[[2;;-1]]-\[Phi]L[[1;;-2]];
eL = (\[Phi][[2;;-1]]-\[Phi][[1;;-2]])/l;
(*========== Estimated_Value_(N=3) =============*)
Rw = 1/2.(\[Phi]L3[3]-\[Phi]L3[1])/(Sqrt[a3[1] (\[Phi]L3[2]-\[Phi]L3[1])]-Sqrt[-a3[2] (\[Phi]L3[3]-\[Phi]L3[2])]);
(*==============================================*)
Return[{Rw,Length[\[Phi]L]-1,\[Phi],\[Phi]L,eL,l}]     ];


(* ::Section::Closed:: *)
(*NewAnsatz*)


NewAnsatz[\[Phi]s_,Ns_,N\[CurlyPhi]_]:=Block[{\[Phi],l,\[Phi]L,eL},
\[Phi]  =\[Phi]s;
l  = Table[Norm[{\[Phi][[s+1]]-\[Phi][[s]]}],{s,1,Ns}];
\[Phi]L = Table[Sum[l[[s1]],{s1,1,s-1}],{s,1,Ns+1}];
eL = ( \[Phi][[2;;-1]]-\[Phi][[1;;-2]])/l;
Return[{Length[\[Phi]L]-1,\[Phi],\[Phi]L,eL,l}]  ];


(* ::Section::Closed:: *)
(*>Polygonal Bounce*)


(* ::Subsection::Closed:: *)
(*RI*)


RI[4,Rw_?NumericQ,a_,\[Phi]L_,pos_,forBack_]:=
	If[forBack=="forward",Sqrt[Rw^2 + (a[[pos]]Rw^2-  Sqrt[a[[pos]]^2 Rw^4+4(a[[pos+1]] - a[[pos]])(\[Phi]L[[pos+1]] - \[Phi]L[[pos]])Rw^2   ]  )/(2(a[[pos+1]] -a[[pos]])) ], 
						  Sqrt[Rw^2 +  Sqrt[(\[Phi]L[[-2]]- \[Phi]L[[-1]])Rw^2 /a[[-2]]  ]    ]          ];
RI[3,Rw_?NumericQ,a_,\[Phi]L_,pos_,forBack_]:= Rw;


(* ::Subsection::Closed:: *)
(*Rvb*)


(*=========== Rn(D) ====================*) 
Rn[4,c1_?NumericQ,a_,b_] := Sqrt[ 1/2 (Sqrt[ c1^2 - 4 a b ] + c1)/a  ];
Rn[3,c1_?NumericQ,a_,b_] := Block[{\[Xi]},  \[Xi]= ( Sqrt[ 36 a b^2 -  c1^3 ] - 6 b a^(1/2)  )^(1/3) /a^(1/2); Return[  1/2 (c1/a/\[Xi] + \[Xi]) ] ];
(*=========== Rvb ======================*) 
Rvb[Rw_?NumericQ,\[Phi]L_,a_,d_,Ns_,forBack_,pos_]:=Rvb[Rw,\[Phi]L,a,d,Ns,forBack,pos] =
Block[{R,b,v,\[Alpha],v1,b1,x,y,z,rvb},
	If[forBack == "backward",
		\[Alpha]=Join[a,{0.}]; R = RI[d,Rw,\[Alpha],\[Phi]L,pos,forBack]; b=0.; v=\[Phi]L[[-1]];
		rvb = Reap[
			Sow[b,x];Sow[v,y];Sow[R,z];
				Do[ b +=-(4./d) (\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^(d);Sow[b,x];
					v +=( 4./(d-2))(\[Alpha][[-i]]-\[Alpha][[-i-1]]) R^2 ;Sow[v,y]; 
					R  = Rn[d,(\[Phi]L[[-i-1]]-v) ,\[Alpha][[-i-1]],b];Sow[R,z];
					,{i,1,Ns}] ;][[2]];
		Return[Reverse[rvb,{1,2}]]      ];
	If[forBack == "forward",
		\[Alpha]=Join[{0.},a ]; R = RI[d,Rw,\[Alpha],\[Phi]L,pos,forBack];b=0.;
		v=\[Phi]L[[pos]]; v+= -(4/d)R^2\[Alpha][[pos]];
		rvb =Reap[
			If[pos>1,Do[Sow[0,x];Sow[v,y];Sow[b,z];,{i,1,pos-1}]];
			Sow[R,x];
			Do[ v+= -( 4./(d-2)) (\[Alpha][[i+1]]-\[Alpha][[i]]) R^2 ;Sow[v,y];
				b+= (4./d) (\[Alpha][[i+1]]-\[Alpha][[i]]) R^(d);Sow[b,z];
				R = Rn[d,(\[Phi]L[[i+1]]-v) ,\[Alpha][[i+1]],b];Sow[R,x];
			,{i,pos,Ns}];
			v+= -( 4./(d-2))(-\[Alpha][[Ns+1]]) R^2 ;Sow[v,y];
			b+= (4./d) (-\[Alpha][[Ns+1]]) R^(d);Sow[b,z]; ][[2]];
		Return[rvb]        ];       ];


(* ::Subsection::Closed:: *)
(*FindRw*)


FindRw[d_,VL_,\[Phi]L_,a_,Ns_, method_,maxIteR_,accuracyB_,ansatzRw_,aRw_]:= 
	Quiet[Block[{R,Rw,timeRw,ite,RW,\[Lambda],Rcomplex,Rreal,\[Lambda]prev,switch,k,Rw0},
		If[method != "DerrickFindRoot"&&method != "DerrickRescaling", Message[FindRw::ErrorMethod]; Abort[];];
(*Picks up the best estimate*)
		If[NumericQ[aRw],Rw =aRw,
			R = Rvb[0.,\[Phi]L,a,d,Ns,"forward",1][[1]]//Chop;
			If[ Abs[ Im[R[[-1]]]]<10^(-12),
				Rw = Abs[R[[-2]]],
				Rw =  Abs[ansatzRw] ];   ];
		Rw0 =Rw;
(*Defines the method to use*)
		\[Lambda][Rw_] := \[Lambda][Rw] = Chop[ Sqrt[\[CapitalLambda][Rw,d,VL,\[Phi]L,a,Ns,"backward",None]] ];
		If[method == "DerrickRescaling",
			RW[Rw_]:= RW[Rw] = Block[{rw} ,rw =Abs[\[Lambda][Rw]] Abs[Rw]; Return[Abs[rw]]  ];          ];
		If[method == "DerrickFindRoot", (*Can be improved by the explicit computation of \[Lambda]*)
			RW[Rw_]:= RW[Rw] = Block[{rw}, Quiet[rw =rw/.FindRoot[Abs[\[Lambda][rw]-1],{rw, Rw}, MaxIterations->1,PrecisionGoal->0,AccuracyGoal->accuracyB ];    ]; 
									Return[Abs[rw]]];      ];
(*Finds interval of the solution Sol \[Element] [Rreal, Rcomplex] or reduces the interval*)
		ite = 0; switch = True; k = 1;
		While[ ite <= maxIteR &&switch &&k<10,    
			Rcomplex = Infinity; Rreal=0; switch = False;
			If[Abs[Im[\[Lambda][Rw]]]>10^(-12),(*Overshooting*)
				While[Abs[Im[\[Lambda][Rw]]]>10^(-12)&&ite <= maxIteR&&Abs[\[Lambda][Rw]-1]>.5 10^(-accuracyB)&&Chop[Re[\[Lambda][Rw]]]!=0.,
					If[Rw < Rcomplex, Rcomplex = Rw;   Rw=RW[Rw] , Rw =.8 Rcomplex];  ite++       ]; 
				Rreal = Rw;,(*Undershooting*)
				(*--------------*)
				While[ Abs@Im[\[Lambda][Rw]]<= 10^(-12)&&ite <= maxIteR&&Abs[\[Lambda][Rw]-1]>.5 10^(-accuracyB)&&Chop@Re[\[Lambda][Rw]]!=0.,
					If[Rw > Rreal,   Rreal = Rw;  Rw=RW[Rw] , Rw = 1.2 Rreal];  ite++      ]; 
				Rcomplex = Rw;    ] ;   
(*Reduces the interval and use bisection method*)	
			While[ite <= maxIteR &&Abs[\[Lambda][Rw]-1]>.5 10^(-accuracyB)&&Chop[Re[\[Lambda][Rw]]]!=0., 
				\[Lambda]prev = \[Lambda][Rw];
				If[Abs@Im[\[Lambda][Rw]]>10^(-12),
					If[ Rw <  Rcomplex,Rcomplex = Rw;Rw =RW[Rw]  ,Rw =Abs[Rcomplex +Rreal]/2.] ;  (*Overshooting*),
					(*---------------*)
					If[ Rw >  Rreal       ,Rreal        = Rw;Rw =RW[Rw]  ,Rw =Abs[Rcomplex +Rreal]/2.](*Undershooting*)];
				ite++; 
				If[Abs[\[Lambda][Rw] -\[Lambda]prev]<.5 10^(-accuracyB),Message[FindRw::ErrorTolerance];Break[];  ];      	];    
			If[ Re[Chop[\[Lambda][Rw]]] ==0., k++;switch=True ;  Rw = 0.5 k  Abs[Rw0]  ];     ];     
			If[ ite > maxIteR||k >10, Message[FindRw::ErrorExceedIteration] ];  
	Return[Rw]	  ];];


(* ::Section::Closed:: *)
(*>Polygonal Bounce Improvement*)


(* ::Subsection::Closed:: *)
(*Find\[ScriptCapitalI]*)


Find\[ScriptCapitalI][v_,\[Phi]L_,a_,b_,Ns_,pos_,R_,ddVL_,d_]:=Find\[ScriptCapitalI][v,\[Phi]L,a,b,Ns,pos,R,ddVL,d]=
	Block[{\[ScriptCapitalI],d\[ScriptCapitalI],v0,\[Phi]L0,a0,b0,ddVL0},
		v0 =Join[{0},v,{0}]; \[Phi]L0 =Join[{0},\[Phi]L,{0}]; a0 =Join[{0},a,{0}];
		b0 =Join[{0},b,{0}];ddVL0 =Join[{0},ddVL,{0}];
If[d==4,
	\[ScriptCapitalI] = Table[Join[ConstantArray[0,pos-1],
	Table[   ddVL0[[i +m]] (   ( v0[[i+m]]-\[Phi]L0[[i+m]] )/8 R[[i]]^2 +a0[[i+m]]/24 R[[i]]^4   + b0[[i+m]]/2 Log[R[[i]]] )  ,{i,pos,Ns+1}]  ],{m,0,1}];
	d\[ScriptCapitalI] = Table[Join[ConstantArray[0,pos-1], 
	Table[   ddVL0[[i +m]] (   ( v0[[i+m]]-\[Phi]L0[[i+m]] )/4 R[[i]] +a0[[i+m]]/6 R[[i]]^3 + b0[[i+m]]/2 /R[[i]]   )   ,{i,pos,Ns+1}]  ],{m,0,1}];          ];
If[d==3,
	\[ScriptCapitalI]=Table[Join[ConstantArray[0,pos-1],
		Table[   ddVL0[[i +m]] (   ( v0[[i+m]]-\[Phi]L0[[i+m]] )/6 R[[i]]^2 +a0[[i+m]]/15 R[[i]]^4   + b0[[i+m]]R[[i]])   ,{i,pos,Ns+1}]  ],{m,0,1}];
	d\[ScriptCapitalI]= Table[Join[ConstantArray[0,pos-1], 
		Table[   ddVL0[[i +m]] (   ( v0[[i+m]]-\[Phi]L0[[i+m]] )/3 R[[i]] +a0[[i+m]] 4/15 R[[i]]^3 + b0[[i+m]] )   ,{i,pos,Ns+1}]  ],{m,0,1}];            ];
Return[{\[ScriptCapitalI],d\[ScriptCapitalI]}];   ];


(* ::Subsection::Closed:: *)
(*r\[Beta]\[Nu]*)


r\[Beta]\[Nu][rw_?NumericQ,a_ ,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_] :=r\[Beta]\[Nu][rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos] = Block[{r\[Beta]\[Nu]M,r,\[Beta],\[Nu],x,y,z,\[Beta]prev,\[Alpha]0,a0,c0,b0,c},
	a0 =Join[{0},a];\[Alpha]0 =Join[{0},\[Alpha]];b0 =Join[{0},b];
	c0 = Table[2 (b0[[i]] - 4/d a0[[i]] R[[i]]^d),{i,1,Ns+1}];
r\[Beta]\[Nu]M=Reap[  
	r=0.;Sow[r,x]; \[Beta] = \[Beta]prev = 0.;Sow[\[Beta],y];
	\[Nu]= (rw c0[[pos]] )R[[pos]]^(2-d)-4/d \[Alpha]0[[pos]] R[[pos]]^(2)-\[ScriptCapitalI][[1,pos]]  ;Sow[\[Nu],z];
	r =rw;Sow[r,x]; 
	\[Beta]+=(4/d (\[Alpha]0[[pos+1]]-\[Alpha]0[[pos]])+4r(a0[[pos+1]]-a0[[pos]]) )   R[[pos]]^d + (  d\[ScriptCapitalI][[2,pos]] -  d\[ScriptCapitalI][[1,pos]]  )R[[pos]]^(d-1) /2;Sow[\[Beta],y];
	\[Nu]+=-2/(d-2)(\[Beta]-\[Beta]prev)R[[pos]]^(2-d)-4/d (\[Alpha]0[[pos+1]]-\[Alpha]0[[pos]])R[[pos]]^2   - (\[ScriptCapitalI][[2,pos]]-\[ScriptCapitalI][[1,pos]]);Sow[\[Nu],z];
	Do[     \[Beta]prev=\[Beta];   
		r =(2/(d-2)\[Beta]+(\[Nu]+\[ScriptCapitalI][[1,i]]+4/d \[Alpha]0[[i]] R[[i]]^2)R[[i]]^(d-2))/c0[[i]];Sow[r,x]; 
		\[Beta]+=4/d (\[Alpha]0[[i+1]]-\[Alpha]0[[i]])R[[i]]^d+4r(a0[[i+1]]-a0[[i]])R[[i]]^d + (  d\[ScriptCapitalI][[2,i]] -  d\[ScriptCapitalI][[1,i]]  )R[[i]]^(d-1) /2;Sow[\[Beta],y];
		\[Nu]+=-2/(d-2)(\[Beta]-\[Beta]prev)/R[[i]]^(d-2)-4/d (\[Alpha]0[[i+1]]-\[Alpha]0[[i]])R[[i]]^2   - (\[ScriptCapitalI][[2,i]]-\[ScriptCapitalI][[1,i]]);Sow[\[Nu],z];
		,{i,pos+1,Ns}];
	r =( 2 \[Beta] / R[[Ns+1]]^(d-1)  - 8/d \[Alpha]0[[Ns+1]] R[[Ns+1]] - d\[ScriptCapitalI][[1,Ns+1]])/(8 a0[[Ns+1]] R[[Ns+1]]);Sow[r,x];   
	\[Beta] = 0;Sow[\[Beta],y];
	\[Nu] = 0; Sow[\[Nu],z];               ][[2]];
Return[r\[Beta]\[Nu]M] ];


(* ::Subsection::Closed:: *)
(*Findrw*)


Findrw[rw_?NumericQ,a_ ,b_,d_,Ns_,\[Alpha]_,R_,\[ScriptCapitalI]_,d\[ScriptCapitalI]_,pos_]:= Findrw[rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos]=
	Block[{bc,r,\[Beta],\[Nu]},
	{r,\[Beta],\[Nu]} = r\[Beta]\[Nu][rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos] ;
	bc = \[Nu][[-2]] + 2/(d-2) \[Beta][[-2]] R[[-1]]^(2-d)+4/d \[Alpha][[Ns]]R[[-1]]^2+\[ScriptCapitalI][[1,-1]];
Return[bc]];


(* ::Section::Closed:: *)
(*>Multi-Field Polygonal Bounce*)


(* ::Subsection::Closed:: *)
(*PathDeformation*)


PathDeformation[\[CurlyPhi]_,dV_,d2V_,Ns_,N\[CurlyPhi]_,pos_,d_,R_,\[Phi]0_,v0_,a0_,b0_,accuracyPath_]:=
Block[{\[Nu],\[Beta],vs,bs,Rs,rI,as,\[Zeta]t,r0,\[Zeta]ts,\[Phi]s,\[Zeta]eqs,\[Nu]\[Beta],x,y,d\[CurlyPhi],rF,rN,\[Zeta]qN,d\[Zeta]qN,p,DV,D2V,
	\[Xi]Mc,M,c,\[Nu]0,\[Beta]0,\[Nu]\[Xi]p,\[Nu]\[Xi]m,\[Beta]\[Xi]p,\[Beta]\[Xi]m,fLowT,fD,fD1,frI,n1,rules},(*Fix rF*)
	p = If[pos >1,pos-1,pos];
	rules = Table[\[CurlyPhi][[i]]->\[Phi]0[[s,i]],{s,p,Ns+1},{i,N\[CurlyPhi]}];
	DV = Chop[Join[ConstantArray[0.,{p-1,N\[CurlyPhi]}],Table[dV/.rules[[s-p+1]],{s,p,Ns+1}]]];
	D2V = Chop[Join[Table[ConstantArray[0.,{N\[CurlyPhi],N\[CurlyPhi]}],{s,1,p-1}], Table[d2V/.rules[[s-p+1]],{s,p,Ns+1}]]];
	d\[CurlyPhi] = Join[Table[ConstantArray[0.,{2,N\[CurlyPhi]}],{s,1,p-1}],Table[   8/d a0[[s+m]] R[[s+1]]- 2 b0[[s+m]]/R[[s+1]]^(d-1)  ,{s,p,Ns-1},{m,0,1}]];
	If[pos>1, \[Nu]0[p] = ConstantArray[0.,N\[CurlyPhi]]; \[Beta]0[p]=ConstantArray[0.,N\[CurlyPhi]];, 
		\[Nu]0[p] = ((16 a0[[1]]-DV[[1]]-DV[[2]]) R[[1]]^2)/(4 (-2+d));  
		\[Beta]0[p] = ((-16 a0[[1]]+DV[[1]]+DV[[2]]) R[[1]]^d)/(4 d);  ];
	Do[ \[Nu]0[s] = \[Nu]0[s-1]+ 1/(4 (-2+d)) R[[s]] (4 d\[CurlyPhi][[-1+s,1]]-4 d\[CurlyPhi][[-1+s,2]]+(-16 a0[[-1+s]]+16 a0[[s]]+DV[[-1+s]]-DV[[1+s]]) R[[s]]);
		\[Beta]0[s] = \[Beta]0[s-1]+ 1/(4 d) R[[s]]^(-1+d) (-2 d d\[CurlyPhi][[-1+s,1]]+2 d d\[CurlyPhi][[-1+s,2]]+(16 a0[[-1+s]]-16 a0[[s]]-DV[[-1+s]]+DV[[1+s]]) R[[s]]);,{s,p+1,Ns}];
	\[Nu]\[Xi]p[s_]:=\[Nu]\[Xi]p[s]= -( (R[[s]]^2) /(4 (-2+d)))D2V[[1+s]];(*\[Zeta]ts[1+s]*)
	\[Beta]\[Xi]p[s_]:=\[Beta]\[Xi]p[s]=(D2V[[1+s]] (R[[s]]^d) )/(4 d);  (*\[Zeta]ts[1+s]*)
	\[Nu]\[Xi]m[s_]:=\[Nu]\[Xi]m[s]=(D2V[[s-1]] (R[[s]]^2) )/(4 (-2+d));  (*\[Zeta]ts[-1+s]*)
	\[Beta]\[Xi]m[s_]:=\[Beta]\[Xi]m[s]=-((D2V[[s-1]] (R[[s]]^d) )/(4 d));  (*\[Zeta]ts[-1+s]*)
	fLowT[s_,j_]:= (2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]m[j]+\[Beta]\[Xi]p[j-2])+\[Nu]\[Xi]m[j]+\[Nu]\[Xi]p[j-2];(*[Subscript[Eq, s],Subscript[\[Xi], j-1]]*) (* 2 \[LessEqual] j \[LessEqual] s*)
	\[Nu]\[Xi]p[p-1]= If[pos>1,IdentityMatrix[N\[CurlyPhi]],-((D2V[[p]] (R[[p]]^2) )/(4 (-2+d))) ];
	\[Beta]\[Xi]p[p-1]= If[pos>1,ConstantArray[0.,{N\[CurlyPhi],N\[CurlyPhi]}],(D2V[[p]] R[[p]]^d)/(4 d)];
	fD[s_] := (D2V[[s]] (R[[1+s]]^2) )/(4 d) +(2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]p[s-1])+ \[Nu]\[Xi]p[s-1];(*[Subscript[Eq, s],Subscript[\[Xi], s]]*)
	\[Nu]\[Xi]p[p] = If[pos>1,ConstantArray[0.,{N\[CurlyPhi],N\[CurlyPhi]}],-((D2V[[1+p]] R[[p]]^2)/(4 (-2+d)))];
	\[Beta]\[Xi]p[p] = If[pos>1,ConstantArray[0.,{N\[CurlyPhi],N\[CurlyPhi]}],(D2V[[1+p]] (R[[p]]^d) )/(4 d)];
	fD1[s_]:= (-IdentityMatrix[N\[CurlyPhi]]+(D2V[[1+s]] R[[1+s]]^2)/(4 d))+(2 (R[[1+s]]^(2-d)) )/(-2+d) (\[Beta]\[Xi]p[s])+\[Nu]\[Xi]p[s] ;(*[Subscript[Eq, s],Subscript[\[Xi], s+1]]*)
	frI[s_]:= ((2 (R[[1+s]]^(2-d)) )/(-2+d) (4 a0[[p]] R[[p]]^d)-(8  a0[[p]] R[[p]]^2)/(-2+d)  )IdentityMatrix[N\[CurlyPhi]]; 
	c  = -Flatten[ {If[pos>1,{},ConstantArray[0.,N\[CurlyPhi]]],Table[ ((-16 a0[[s]]+DV[[s]]+DV[[1+s]]) R[[1+s]]^2)/(4 d)+(2 (R[[1+s]]^(2-d)) )/(-2+d) \[Beta]0[s]  +\[Nu]0[s]  ,{s,p,Ns}],ConstantArray[0.,N\[CurlyPhi]]}    ];
	n1 = If[pos>1,0,1];
	M = Block[{n,m},
		n=(Ns+1-p)+1+n1;
		m = SparseArray[{{n,n}->0}];
		Do[m[[s+n1,s+n1]]=fD[s+p-1],{s,1,n-1-n1}];
		Do[m[[s+n1,s+1+n1]]=fD1[s+p-1],{s,1,n-1-n1}];
		m[[n,n]]=IdentityMatrix[N\[CurlyPhi]];
		If[pos==1,
			m[[1,2]]=IdentityMatrix[N\[CurlyPhi]];
			Do[m[[s+n1,1]]=frI[s],{s,1,n-1-n1}];   ];
		Do[m[[s+n1,j-1+n1]]=fLowT[s+p-1,j+p-1]; ,{s,2,n-n1-1},{j,2,s}];
		ArrayFlatten[m]   ];
	\[Xi]Mc = LinearSolve[M,c];  
	\[Zeta]ts = Join[ConstantArray[0.,{p-1,N\[CurlyPhi]}],Partition[\[Xi]Mc[[ 1+n1 N\[CurlyPhi];;-1]],N\[CurlyPhi]]];
	as  = Join[ConstantArray[0.,{p-1,N\[CurlyPhi]}],Table[   1/8(      (  DV[[s]]+DV[[s+1]]+D2V[[s]].\[Zeta]ts[[s]]+D2V[[s+1]].\[Zeta]ts[[s+1]]     )/2)-a0[[s]],{s,p,Ns}  ]];
	rI  = If[pos>1,ConstantArray[0,N\[CurlyPhi]],\[Xi]Mc[[1;;N\[CurlyPhi]]] ];
	rF  = ConstantArray[0,N\[CurlyPhi]];
	If[pos>1,
		\[Nu] = \[Zeta]ts[[p]];\[Beta] = ConstantArray[0.,N\[CurlyPhi]];,
		\[Nu] =-4/(d-2) (as[[1]]+2 a0[[1]] rI)R[[1]]^2;
		\[Beta] =   4/d (as[[1]]+d a0[[1]] rI)R[[1]]^d;    ];     
	\[Nu]\[Beta]=Reap[
		If[ pos>1, Do[Sow[ConstantArray[0,N\[CurlyPhi]],x];Sow[ConstantArray[0,N\[CurlyPhi]],y];,{s,1,p-1}];];
		Sow[\[Nu],x];Sow[\[Beta],y];
		Do[ \[Nu]+= -4/(d-2)(as[[s+1]]-as[[s]])R[[s+1]]^2-1/(d-2) (  d\[CurlyPhi][[s,2]]-d\[CurlyPhi][[s,1]] )R[[s+1]] ;
			\[Beta]+=   4/d(as[[s+1]]-as[[s]])R[[s+1]]^d+1/2 (  d\[CurlyPhi][[s,2]]-d\[CurlyPhi][[s,1]] )R[[s+1]]^(d-1);
			Sow[\[Nu],x];Sow[\[Beta],y]; ,{s,p,Ns-1} ];  ][[2]];
Return[ Chop[{\[Zeta]ts,\[Nu]\[Beta][[1]],as,\[Nu]\[Beta][[2]],R[[p]] rI,R[[Ns+1]] rF}] ];   ];


(* ::Subsection::Closed:: *)
(*\[Phi]vabRs*)


\[Phi]vabRs[\[CurlyPhi]_,dV_,d2V_,\[Phi]_,eL_,l_,\[Phi]L_,v_,a_,b_,Ns_,R_,N\[CurlyPhi]_,pos_,d_,accuracyPath_]:=
Block[{\[Phi]vabR,\[Phi]s,vs,as,bs,Rs,\[Zeta]s,x1,x2,x3,x4,x5,ps},
	ps = If[pos >1,pos-1,pos];
	{\[Phi]s,vs,as,bs,Rs} = {\[Phi],Table[\[Phi][[s+1]]+eL[[s]](v[[s]]-(l[[s]]+\[Phi]L[[s]])),{s,1,Ns}],Table[eL[[s]]a[[s]],{s,1,Ns}],Table[eL[[s]]b[[s]],{s,1,Ns}],Transpose[Table[R,{i,1,N\[CurlyPhi]}]]};
	{ \[Phi]s,vs,as,bs,Rs[[ps]],Rs[[Ns+1]] } += PathDeformation[\[CurlyPhi],dV,d2V,Ns,N\[CurlyPhi],pos,d,Rs,\[Phi]s ,vs,as,bs,accuracyPath];
Return[{\[Phi]s,vs,as,bs,Rs}] ];


(* ::Section::Closed:: *)
(*\[ScriptCapitalT],\[ScriptCapitalV],\[ScriptCapitalS],\[CapitalLambda]*)


\[ScriptCapitalT][v_,a_,b_,R_,d_,Ns_]:= \[ScriptCapitalT][v,a,b,R,d,Ns] = 
	2\[Pi]^(d/2)/Gamma[d/2]Sum[32 a[[i]]^2/(d^2(d+2)) (R[[i+1]]^(2+d) - R[[i]]^(2+d)) -8 a[[i]]b[[i]]/d  (R[[i+1]]^2 -R[[i]]^2) -
	If[ Re[R[[i]]] <10^(-1.) , 0.,  (2/(d-2))b[[i]]^2 (R[[i+1]]^(2-d)-R[[i]]^(2-d)) ],{i,1,Ns}]  ;
\[ScriptCapitalV][v_,a_,b_,R_,d_,VL_,\[Phi]L_,Ns_]:= \[ScriptCapitalV][v,a,b,R,d,VL,\[Phi]L,Ns] = 2\[Pi]^(d/2)/Gamma[d/2](
	Sum[ 32 a[[i]]^2/(d(d+2))  (R[[i+1]]^(2+d) - R[[i]]^(2+d))  + 
		8 a[[i]]b[[i]]/(d-2) (R[[i+1]]^2 -R[[i]]^2)  +
		( VL[[i]]-VL[[-1]]+ 8 a[[i]]( v[[i]] - \[Phi]L[[i]]))(R[[i+1]]^d-R[[i]]^d)/d ,{i,1,Ns}]+
	1/d R[[1]]^d (VL[[1]] - VL[[-1]]) ) ;
\[CapitalLambda][Rw_?NumericQ,d_,VL_,\[Phi]L_,a_,Ns_,forBack_,pos_]:=\[CapitalLambda][Rw,d,VL,\[Phi]L,a,Ns,forBack,pos] = 
	Block[{R,v,b}, {R,v,b} = Rvb[Rw,\[Phi]L,a,d,Ns,forBack,pos]; 
	Return[(2-d)\[ScriptCapitalT][v,a,b,R,d,Ns] /(d \[ScriptCapitalV][v,a,b,R,d,VL,\[Phi]L,Ns])] ];
\[ScriptCapitalS][v_,a_,b_,R_,d_,Ns_,VL_,\[Phi]L_]:= \[ScriptCapitalS][v,a,b,R,d,Ns,VL,\[Phi]L] = \[ScriptCapitalT][v,a,b,R,d,Ns] + \[ScriptCapitalV][v,a,b,R,d,VL,\[Phi]L,Ns];


(* ::Section::Closed:: *)
(*>\[ScriptCapitalT]\[Xi], \[ScriptCapitalV]\[Xi], \[ScriptCapitalS]\[Xi]*)


(* ::Subsection::Closed:: *)
(*\[ScriptCapitalT]\[Xi]*)


\[ScriptCapitalT]\[Xi]4[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] := \[ScriptCapitalT]\[Xi]4[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r] = 1/(24 \[Rho]^2) (a \[Rho]^4 (-48 \[Beta]+\[Rho]^4 (2 ddVL v+16 \[Alpha]+a ddVL \[Rho]^2-2 ddVL \[Phi]L))+b (-48 \[Beta]+2 \[Rho]^4 (-3 ddVL v-24 \[Alpha]+2 a ddVL \[Rho]^2+3 ddVL \[Phi]L))-24 b^2 ddVL \[Rho]^2 Log[\[Rho]])+ (r (4 b-4 a \[Rho]^4)^2)/(8 \[Rho]^2);
\[ScriptCapitalT]\[Xi]3[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] := \[ScriptCapitalT]\[Xi]3[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r] = -((4 b \[Beta])/\[Rho])-2 b^2 ddVL \[Rho]+8/15 a b ddVL \[Rho]^4+32/315 a^2 ddVL \[Rho]^7-1/3 \[Rho]^2 (8 a \[Beta]+b (8 \[Alpha]+ddVL (v-\[Phi]L)))+8/45 a \[Rho]^5 (8 \[Alpha]+ddVL (v-\[Phi]L))+(2 r (3 b-4 a \[Rho]^3)^2)/(9 \[Rho]);
\[ScriptCapitalT]\[Xi][d_?NumericQ,a_,R_,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,ddVL_,VL_,\[Phi]L_,Ns_,pos_,r_] := \[ScriptCapitalT]\[Xi][d,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r] =
	Block[{\[ScriptCapitalT],VN,v0,\[Phi]L0,a0,b0,VL0,\[Alpha]0,\[Nu]0,\[Beta]0,ddVL0,\[ScriptCapitalT]\[Xi]D},VN = VL[[-1]];
		a0 =Join[{0},a,{0}];VL0 =Join[{VL[[1]]},VL,{0}];b0 =Join[{0},b,{0}];v0 =Join[{0},v,{0}];\[Phi]L0 =Join[{0},\[Phi]L,{0}];
		\[Alpha]0 =Join[{0},\[Alpha],{0}];\[Nu]0 =Join[{0},\[Nu],{0}];\[Beta]0 =Join[{0},\[Beta],{0}];ddVL0 =Join[{ddVL[[1]]},ddVL,{0}];
		If[d==4,\[ScriptCapitalT]\[Xi]D = \[ScriptCapitalT]\[Xi]4];
		If[d==3,\[ScriptCapitalT]\[Xi]D = \[ScriptCapitalT]\[Xi]3];
		\[ScriptCapitalT]= \[ScriptCapitalT]\[Xi]D[a0[[pos]],R[[pos]],b0[[pos]],v0[[pos]],\[Alpha]0[[pos]],\[Beta]0[[pos]],\[Nu]0[[pos]],VL0[[pos]],VN,ddVL0[[pos]],\[Phi]L0[[pos]],r[[pos]]];
		Do[ \[ScriptCapitalT] += \[ScriptCapitalT]\[Xi]D[a[[i]],R[[i+1]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i+1]]]-
			\[ScriptCapitalT]\[Xi]D[a[[i]],R[[i]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i]]]   ,{i,pos,Ns}];
	Return[2\[Pi]^(d/2)/Gamma[d/2]\[ScriptCapitalT]] ];


(* ::Subsection::Closed:: *)
(*\[ScriptCapitalV]\[Xi]*)


\[ScriptCapitalV]\[Xi]4[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=\[ScriptCapitalV]\[Xi]4[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r] =1/48 \[Rho]^2 (5 a^2 ddVL \[Rho]^6+16 a (12 \[Beta]+6 \[Nu] \[Rho]^2+\[Rho]^4 (8 \[Alpha]+ddVL (v-\[Phi]L)))+24 b (8 \[Alpha]+ddVL (v-\[Phi]L))+6 \[Rho]^2 (16 \[Alpha]+ddVL (v-\[Phi]L)) (v-\[Phi]L))+1/2 b ddVL (b+2 a \[Rho]^4) Log[\[Rho]]+8 a b r \[Rho]^2+1/4 r \[Rho]^4 (4 (VL-VN)+32 a^2 \[Rho]^2+32 a (v-\[Phi]L));
\[ScriptCapitalV]\[Xi]3[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=\[ScriptCapitalV]\[Xi]3[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r]=1/630 \[Rho] (1260 b^2 ddVL+210 b \[Rho] (24 \[Alpha]+ddVL (3 v+8 a \[Rho]^2-3 \[Phi]L))+\[Rho] (128 a^2 ddVL \[Rho]^5+336 a (15 \[Beta]+5 \[Nu] \[Rho]+\[Rho]^3 (8 \[Alpha]+ddVL (v-\[Phi]L)))+105 \[Rho] (16 \[Alpha]+ddVL (v-\[Phi]L)) (v-\[Phi]L)))+16 a b r \[Rho]^2+1/3 r \[Rho]^3 (3 (VL-VN)+32 a^2 \[Rho]^2+24 a (v-\[Phi]L));
\[ScriptCapitalV]\[Xi][d_?NumericQ,a_,R_,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,ddVL_,VL_,\[Phi]L_,Ns_,pos_,r_] := \[ScriptCapitalV]\[Xi][d,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r]=
	Block[{\[ScriptCapitalV],VN,v0,\[Phi]L0,a0,b0,VL0,\[Alpha]0,\[Nu]0,\[Beta]0,ddVL0,\[ScriptCapitalV]\[Xi]D},
		VN = VL[[-1]];a0 =Join[{0},a,{0}];VL0 =Join[{VL[[1]]},VL,{0}];
		b0 =Join[{0},b,{0}]; v0 =Join[{0},v,{0}]; \[Phi]L0 =Join[{0},\[Phi]L,{0}];
		\[Alpha]0 =Join[{0},\[Alpha],{0}]; \[Nu]0 =Join[{0},\[Nu],{0}]; \[Beta]0 =Join[{0},\[Beta],{0}]; ddVL0 =Join[{ddVL[[1]]},ddVL,{0}];
		If[d==4,\[ScriptCapitalV]\[Xi]D = \[ScriptCapitalV]\[Xi]4];
		If[d==3,\[ScriptCapitalV]\[Xi]D = \[ScriptCapitalV]\[Xi]3];
		\[ScriptCapitalV]= \[ScriptCapitalV]\[Xi]D[a0[[pos]],R[[pos]],b0[[pos]],v0[[pos]],\[Alpha]0[[pos]],\[Beta]0[[pos]],\[Nu]0[[pos]],VL0[[pos]],VN,ddVL0[[pos]],\[Phi]L0[[pos]],r[[pos]]];
		Do[\[ScriptCapitalV]+=\[ScriptCapitalV]\[Xi]D[a[[i]],R[[i+1]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i+1]]]-
			\[ScriptCapitalV]\[Xi]D[a[[i]],R[[i]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i]]]   ,{i,pos,Ns}];
	Return[2\[Pi]^(d/2)/Gamma[d/2]\[ScriptCapitalV]] ];


(* ::Subsection::Closed:: *)
(*\[ScriptCapitalS]\[Xi]*)


\[ScriptCapitalS]\[Xi]4[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=\[ScriptCapitalS]\[Xi]4[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r] = (\[Rho]^4 (7 a^2 ddVL \[Rho]^6+4 a (24 \[Beta]+24 \[Nu] \[Rho]^2+5 \[Rho]^4 (8 \[Alpha]+ddVL (v-\[Phi]L)))+6 \[Rho]^2 (16 \[Alpha]+ddVL (v-\[Phi]L)) (v-\[Phi]L))+4 b (-24 \[Beta]+\[Rho]^4 (3 ddVL v+24 \[Alpha]+2 a ddVL \[Rho]^2-3 ddVL \[Phi]L))-24 b ddVL \[Rho]^2 (b-2 a \[Rho]^4) Log[\[Rho]])/(48 \[Rho]^2)+(r (64 b^2 \[Rho]^2+128 a b \[Rho]^6+2 \[Rho]^8 (16 (VL-VN)+160 a^2 \[Rho]^2+128 a (v-\[Phi]L))))/(32 \[Rho]^4);
\[ScriptCapitalS]\[Xi]3[a_,\[Rho]_?NumericQ,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,VL_,VN_,ddVL_,\[Phi]L_,r_] :=\[ScriptCapitalS]\[Xi]3[a,\[Rho],b,v,\[Alpha],\[Beta],\[Nu],VL,VN,ddVL,\[Phi]L,r]= (\[Rho]^3 (192 a^2 ddVL \[Rho]^5+112 a (30 \[Beta]+15 \[Nu] \[Rho]+4 \[Rho]^3 (8 \[Alpha]+ddVL (v-\[Phi]L)))+105 \[Rho] (16 \[Alpha]+ddVL (v-\[Phi]L)) (v-\[Phi]L))+84 b (-30 \[Beta]+\[Rho]^3 (5 ddVL v+40 \[Alpha]+24 a ddVL \[Rho]^2-5 ddVL \[Phi]L))  )/(630 \[Rho])+(r (18 b^2 \[Rho]^2+96 a b \[Rho]^5+\[Rho]^6 (9 (VL-VN)+128 a^2 \[Rho]^2+72 a (v-\[Phi]L))))/(9 \[Rho]^3);
\[ScriptCapitalS]\[Xi][d_?NumericQ,a_,R_,b_,v_,\[Alpha]_,\[Beta]_,\[Nu]_,ddVL_,VL_,\[Phi]L_,Ns_,pos_,r_] :=\[ScriptCapitalS]\[Xi][d,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r] = 
Block[{\[ScriptCapitalS],VN,v0,\[Phi]L0,a0,b0,VL0,\[Alpha]0,\[Nu]0,\[Beta]0,ddVL0,\[ScriptCapitalS]\[Xi]D},
	VN = VL[[-1]]; a0 =Join[{0},a,{0}];VL0 =Join[{VL[[1]]},VL,{0}];b0 =Join[{0},b,{0}];
	v0 =Join[{0},v,{0}];\[Phi]L0 =Join[{0},\[Phi]L,{0}];\[Alpha]0 =Join[{0},\[Alpha],{0}];\[Nu]0 =Join[{0},\[Nu],{0}];
	\[Beta]0 =Join[{0},\[Beta],{0}];ddVL0 =Join[{0},ddVL,{0}];
	If[d==4,\[ScriptCapitalS]\[Xi]D = \[ScriptCapitalS]\[Xi]4];
	If[d==3,\[ScriptCapitalS]\[Xi]D = \[ScriptCapitalS]\[Xi]3];
	\[ScriptCapitalS]=\[ScriptCapitalS]\[Xi]D[a0[[pos]],R[[pos]],b0[[pos]],v0[[pos]],\[Alpha]0[[pos]],\[Beta]0[[pos]],\[Nu]0[[pos]],VL0[[pos]],VN,ddVL0[[pos]],\[Phi]L0[[pos]],r[[pos]]];
	Do[\[ScriptCapitalS]+=\[ScriptCapitalS]\[Xi]D[a[[i]],R[[i+1]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i+1]]]-
		\[ScriptCapitalS]\[Xi]D[a[[i]],R[[i]],b[[i]],v[[i]],\[Alpha][[i]],\[Beta][[i]],\[Nu][[i]],VL[[i]],VN,ddVL[[i]],\[Phi]L[[i]],r[[i]]]   ,{i,pos,Ns}];
Return[2\[Pi]^(d/2)/Gamma[d/2]\[ScriptCapitalS]] ];


(* ::Section::Closed:: *)
(*FindBounce*)


FindBounce[V_,\[CurlyPhi]_,{min1_,min2_},OptionsPattern[]]:=
Block[{a,Rw,aPath,\[Phi]L,ansatzRw,b,v,\[Phi],Ns,d,Nfv,aRw,accuracyB,accuracyPath,
	N\[CurlyPhi],VL,d\[Phi]L,point,itePath,k=0,maxIteR,ps,R,forBack,methodRw,methodSeg,improvePB,ImprovePB,
	dVL,ddVL,\[Alpha],\[Beta],\[Nu],r,pos,rw,r1,\[ScriptCapitalI],d\[ScriptCapitalI],l,eL,dV,d2V,\[Phi]l,\[Phi]s=Null,vs=Null,as=Null,bs=Null,Rs=Null,V\[Xi],T\[Xi],V1,T1,aV1,\[Phi]N3,Action,rule,PB,MultiPB,improvementPB},
(*=OPTION========== OptionValues ===========*)
	aRw = OptionValue["AnsatzRadii"];
	aPath = OptionValue["AnsatzPath"];
	aV1 = OptionValue["AnsatzV1"];
	accuracyB = OptionValue["AccuracyBounce"];
	accuracyPath = OptionValue["AccuracyPathDeformation"];	
	d = IntegerPart[OptionValue["Dimension"]]; If[d<3||d>4, Message[FindBounce::Error];Abort[];];  
	forBack = OptionValue["ForwardBackward"];
	Nfv = OptionValue["NumberSegments"]+1; (*field values = "number of segement" + 1*)
	maxIteR = OptionValue["MaxIterationsR"];
	itePath =OptionValue["MaxIterationsPathDeformation"];
	methodRw = OptionValue["MethodBounce"];
	methodSeg =OptionValue["MethodSegmentation"];
	improvePB = OptionValue["ImprovementPolygonalBounce"];
	N\[CurlyPhi] = Length[min1];If[N\[CurlyPhi]==0,N\[CurlyPhi]=1;]; (*Number of Fields*)
	point= If[OptionValue["InitialFieldValue"] === Null,
			(min1+min2)/2,OptionValue["InitialFieldValue"] ];
	If[itePath<0,Message[FindBounce::ErrorIte];Abort[]; ];
	If[N\[CurlyPhi]<=1,itePath=0;];
(*=========== Derivative of the Potential ===*)
	If[OptionValue[Gradient] === None,
		dV = If[Length[\[CurlyPhi]] == 1, D[V,\[CurlyPhi]], D[V,{\[CurlyPhi]}]],
		dV = OptionValue[Gradient]  ];
	If[OptionValue[Hessian] === None && N\[CurlyPhi]>1,
		d2V = D[V,{\[CurlyPhi]},{\[CurlyPhi]}],
		d2V = OptionValue[Hessian]  ];
(*=BEGIN ========= Estimations ==============*)
	\[Phi]N3 = N[{min1,point,min2}];
	If[aPath===None||aRw ===None,
		{ansatzRw,Ns,\[Phi],\[Phi]L,eL,l}= AnsatzN3[V,\[CurlyPhi],\[Phi]N3,N\[CurlyPhi],Nfv,aV1,methodSeg];   ];
	If[aPath =!= None, 
		If[Length[aPath[[1]]]==N\[CurlyPhi]||N\[CurlyPhi]==1&&Length[aPath[[1]]]==0,
			{Ns,\[Phi],\[Phi]L,eL,l}=NewAnsatz[aPath,Length[aPath]-1,N\[CurlyPhi]],
			Message[FindBounce::AnsatzPath];Abort[];]  ];
	While[k <= itePath,
(*========= Single_Field_Polygonal_bounce ==*)
		rule = If[Length[\[CurlyPhi]] == 1,Table[\[CurlyPhi][[1]]->\[Phi][[s]],{s,Ns+1}],If[Length[\[CurlyPhi]] == 0, Table[\[CurlyPhi]->\[Phi][[s]],{s,Ns+1}],Table[\[CurlyPhi][[i]]-> \[Phi][[s,i]],{s,Ns+1},{i,N\[CurlyPhi]}]]];
		VL = If[aV1===None,Table[V/.rule[[s]],{s,Ns+1}],aV1];
		If[VL[[1]]>VL[[2]]||VL[[-1]]>VL[[-2]],
			If[k === 0,Message[FindBounce::ErrorExtrema];Abort[];, Message[FindBounce::ErrorPathDeformation];Abort[];  ]];
		a  = Table[ (VL[[s+1]]-VL[[s]])/(\[Phi]L[[s+1]]-\[Phi]L[[s]]) 1/8. ,{s,Ns} ]; 
		pos = FindSegment[a,\[Phi]L,d,Ns];
		Rw = FindRw[d,VL,\[Phi]L,a,Ns,methodRw,maxIteR,accuracyB,ansatzRw,aRw]//Re;
		If[Rw<10^(-5.),Message[FindBounce::Error];Abort[]];
		{R,v,b} = Rvb[Rw,\[Phi]L,a,d,Ns,"backward",None]//Re;
		T1= \[ScriptCapitalT][v,a,b,R,d,Ns]; V1= \[ScriptCapitalV][v,a,b,R,d,VL,\[Phi]L,Ns]; a = Join[a,{0.}];
(*========= Single_Field_Improvement========*) 
		If[improvePB&&k ==itePath,
			dVL= If[N\[CurlyPhi]==1,Table[(dV/.rule[[s]]),{s,Ns+1}],Table[(dV/.rule[[s]]).eL[[If[s==Ns+1,Ns,s]]],{s,Ns+1}]];
			\[Alpha]  = Join[a[[1;;-2]] - dVL[[2;;-1]]/8 ,{0}];
			ddVL = Join[Table[  (dVL[[s+1]]-8(a[[s]] + \[Alpha][[s]]))/(\[Phi]L[[s+1]]-\[Phi]L[[s]])  ,{s,Ns}  ],{0}];
			{\[ScriptCapitalI],d\[ScriptCapitalI]} = Find\[ScriptCapitalI][v,\[Phi]L,a,b,Ns,pos,R,ddVL,d];
			r1 = rw/.FindRoot[Findrw[rw,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos],{rw,-1}]//Quiet;
			{r,\[Beta],\[Nu]} = If[pos>1,Join[ConstantArray[0,{pos-2,3}],r\[Beta]\[Nu][r1,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos]//Transpose]//Transpose,r\[Beta]\[Nu][r1,a,b,d,Ns,\[Alpha],R,\[ScriptCapitalI],d\[ScriptCapitalI],pos][[All,2;;-1]]];
			T\[Xi] = \[ScriptCapitalT]\[Xi][d,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r] ;
			V\[Xi] = \[ScriptCapitalV]\[Xi][d,a,R,b,v,\[Alpha],\[Beta],\[Nu],ddVL,VL,\[Phi]L,Ns,pos,r];,
			(*-else-----------*)
			{r,\[Beta],\[Nu],\[Alpha],ddVL}=ConstantArray[0,{5,Ns+1}];
			{\[ScriptCapitalI],d\[ScriptCapitalI]}=ConstantArray[0,{2,2,Ns+1}];
			{T\[Xi],V\[Xi]}={0,0}; r1 = 0;];
(*========= Perturbation Multi-Field ======*)
		If[k ==itePath || N\[CurlyPhi]==1, Break[];];
		{\[Phi]s,vs,as,bs,Rs} = Chop@\[Phi]vabRs[\[CurlyPhi],dV,d2V,\[Phi],eL,l,\[Phi]L,v,a,b,Ns,R,N\[CurlyPhi],pos,d,accuracyPath];
		{Ns,\[Phi],\[Phi]L,eL,l} = NewAnsatz[\[Phi]s,Ns,N\[CurlyPhi]];
		ansatzRw = Rw;
		k++];
	Action = (T1+V1)+(T\[Xi]+V\[Xi]);
	PB = {d,\[Phi]L,VL,\[Phi],pos,v,a,b,R};
	ImprovePB = {ddVL,\[Nu],\[Alpha],\[Beta]};
	MultiPB = {\[Phi]s,vs,as,bs,Rs};
Return[ {Action,PB,ImprovePB,MultiPB}]   ];


(* ::Chapter:: *)
(*Sub-Packages*)


(* ::Section::Closed:: *)
(*BounceAction*)


BounceAction[FindBounce_List]:=FindBounce[[1]];


(* ::Section::Closed:: *)
(*BounceDimension*)


BounceDimension[FindBounce_List]:=FindBounce[[2,1]];


(* ::Section::Closed:: *)
(*PlotBounce*)


PlotBounce[FindBounce_List]:=Block[{Action,PB,MultiPB,ImprovePB,\[Phi]L,VL,\[Phi],pos,v,a,b,R,ddVL,\[Nu],\[Alpha],\[Beta],r,\[Phi]s,vs,as,bs,Rs,list\[Phi]LVL,listR\[Phi]L,d,Ns},
{Action,PB,ImprovePB,MultiPB}=FindBounce;
{d,\[Phi]L,VL,\[Phi],pos,v,a,b,R}=PB;
{ddVL,\[Nu],\[Alpha],\[Beta]}=ImprovePB;
{\[Phi]s,vs,as,bs,Rs}=MultiPB;
Ns = Length[R]-1;
(*======================*)
list\[Phi]LVL = Show[ListPlot[Transpose@{\[Phi]L,VL},PlotStyle->Blue,Joined->True,Frame->True,Mesh->Full, LabelStyle->Directive[Black,FontSize->20, FontFamily->"Times New Roman",FontSlant->Plain],FrameLabel->{"\[CurlyPhi]","V(\[CurlyPhi])"},GridLines->Automatic,ImageSize->Medium]];
(*====================*)
listR\[Phi]L = Show[ 
If[pos==1,Plot[  Sign[\[Phi][[-1]]-\[Phi][[1]]](  v[[pos]]+ 4/d a[[pos]]R[[pos]]^2 +2/(d-2)b[[pos]]/R[[pos]]^(d-2)-Abs[\[Phi][[-1]]-\[Phi][[1]]])+ \[Phi][[-1]],{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False,Epilog->{Text[Style[Row[{"Action = ",Action}],10,15,"Times New Roman"],Scaled[{0.75,0.9}],Background->LightGreen]} ],Plot[ Sign[\[Phi][[-1]]-\[Phi][[1]]]( (v[[pos-1]]+ 4/d a[[pos-1]]\[Rho]^2 )-Abs[\[Phi][[-1]]-\[Phi][[1]]])+ \[Phi][[-1]],{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False]],
(*====================*)
Table[Plot[{  Sign[\[Phi][[-1]]-\[Phi][[1]]]( ( v[[s]] + 4/d a[[s]]\[Rho]^2 +If[\[Rho]>0,2/(d-2)b[[s]]/\[Rho]^(d-2),0] ) -Abs[\[Phi][[1]]])+ Abs[\[Phi][[-1]]] },{\[Rho],R[[s]],R[[s+1]]},PlotStyle->RGBColor[1,0.5,1/8]],{s,pos,Ns}],
(*====================*)
Plot[{  Sign[\[Phi][[-1]]-\[Phi][[1]]](( v[[Ns]] + 4/d a[[Ns]]R[[-1]]^2 +2/(d-2)b[[Ns]]/R[[-1]]^(d-2) ) -Abs[\[Phi][[1]]])+ Abs[\[Phi][[-1]]]   },{\[Rho],R[[-1]],R[[-1]]2},PlotStyle->RGBColor[1,0.5,1/8]],
(*====================*)
ListPlot[{R,\[Phi]}//Transpose,PlotStyle->{RGBColor[1,0.5,1/8],PointSize[.017]}],
Frame->True,
FrameLabel->{"\[Rho]","\[CurlyPhi](\[Rho])"},LabelStyle->Directive[Black,FontSize->17, FontFamily->"Times New Roman",FontSlant->Plain],GridLines->Automatic ,GridLinesStyle->GrayLevel[.85],ImageSize->Medium,PlotRange->All ] ;
(*====================*)
{list\[Phi]LVL,listR\[Phi]L}   ];


(* ::Section::Closed:: *)
(*PlotBounce2D*)


PlotBounce2D[FindBounce_List]:=Block[{Action,PB,improvementPB,MultiPB,\[Phi]L,VL,\[Phi],pos,v,a,b,R,ddVL,\[Nu],\[Alpha],\[Beta],r,\[Phi]s,vs,as,bs,Rs,list\[Phi]LVL,listR\[Phi]L,d,Ns,imagesize,fontsize,extremaV,ps,PBPath,plot},
(*======================*)
{Action,PB,improvementPB,MultiPB}=FindBounce;
{d,\[Phi]L,VL,\[Phi],pos,v,a,b,R}=PB;
{ddVL,\[Nu],\[Alpha],\[Beta]}=improvementPB;
{\[Phi]s,vs,as,bs,Rs}=MultiPB;
Ns = Length[R]-1;
(*======================*)
imagesize = 400;fontsize =20;
extremaV =ListPlot[ \[Phi][[{1,-1}]] ,PlotStyle->{PointSize[.015],Black} ,LabelStyle->Directive[Black,FontSize->20, FontFamily->"Times New Roman",FontSlant->Plain],Frame->True,FrameLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)"},GridLinesStyle->GrayLevel[.85],GridLines->Automatic ];
ps = If[pos>1,pos-1,pos];
PBPath =  Show[ Table[  ParametricPlot[
vs[[s]]+4/d as[[s]]\[Rho]^2+If[pos>1&&s==ps,0.,2/(d-2)bs[[s]]/\[Rho]^(d-2)],{\[Rho],Min[Rs[[s]]],Min[Rs[[s+1]]]},PlotStyle->Blue, LabelStyle->Directive[Black,FontSize->16, FontFamily->"Times New Roman",FontSlant->Plain],AxesLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)"},AxesOrigin->{-1,-1}],{s,ps,Ns}],ListPlot[    If[ pos>1,Join[{vs[[pos-1]]},\[Phi]s[[pos;;-1]]],\[Phi]s]  ,PlotStyle->{Lighter[Blue,.5],Dashed},Joined->True,Mesh->All],PlotRange->All,ImageSize->400] ;
list\[Phi]LVL = ListPlot[Transpose@{\[Phi]L,VL/100},PlotStyle->Blue,Joined->True,Frame->True,Mesh->Full, LabelStyle->Directive[Black,FontSize->20, FontFamily->"Times New Roman",FontSlant->Plain],FrameLabel->{"\[CurlyPhi]","V(\[CurlyPhi]) \!\(\*SuperscriptBox[\(.10\), \(2\)]\)"},GridLines->Automatic,ImageSize->Medium];
listR\[Phi]L = Show[ 
If[pos==1,Plot[  vs[[pos]]+ 4/d as[[pos]]R[[pos]]^2 +2/(d-2)bs[[pos]]/R[[pos]]^(d-2),{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False,Epilog->{Text[Style[Row[{"Action = ",Action}],10,15,"Times New Roman"],Scaled[{0.75,0.9}],Background->LightGreen]} ],Plot[(vs[[pos-1]]+ 4/d as[[pos-1]]\[Rho]^2 ),{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False,Epilog->{Text[Style[Row[{"Action = ",Action}],10,15,"Times New Roman"],Scaled[{0.75,0.9}],Background->LightGreen]} ]],
(*====================*)
Table[Plot[( vs[[s]] + 4/d as[[s]]\[Rho]^2 +If[\[Rho]>0,2/(d-2)bs[[s]]/\[Rho]^(d-2),0] ),{\[Rho],R[[s]],R[[s+1]]},PlotStyle->RGBColor[1,0.5,1/8]],{s,pos,Ns}],
(*====================*)
Plot[( vs[[Ns]] + 4/d as[[Ns]]Rs[[-1]]^2 +2/(d-2)bs[[Ns]]/Rs[[-1]]^(d-2) ),{\[Rho],R[[-1]],R[[-1]]2},PlotStyle->RGBColor[1,0.5,1/8]],
(*====================*)
Frame->True,
FrameLabel->{"\[Rho]","\[CurlyPhi](\[Rho])"},LabelStyle->Directive[Black,FontSize->17, FontFamily->"Times New Roman",FontSlant->Plain],GridLines->Automatic ,GridLinesStyle->GrayLevel[.85],ImageSize->400,PlotRange->All ] ;
plot =Show[extremaV,PBPath,ImageSize->400];
{plot,list\[Phi]LVL ,listR\[Phi]L } ];


(* ::Section::Closed:: *)
(*PlotBounce2DProjection*)


PlotBounce2DProjection[FindBounce_List,{field1_Integer,field2_Integer}]:=Block[{Action,PB,improvementPB,MultiPB,\[Phi]L,VL,\[Phi],pos,v,a,b,R,ddVL,\[Nu],\[Alpha],\[Beta],r,\[Phi]s,vs,as,bs,Rs,list\[Phi]LVL,listR\[Phi]L,d,Ns,imagesize,fontsize,ps,z,zPBPath,extremaV,plot},
(*======================*)
{Action,PB,improvementPB,MultiPB}=FindBounce;
{d,\[Phi]L,VL,\[Phi],pos,v,a,b,R}=PB;
{ddVL,\[Nu],\[Alpha],\[Beta]}=improvementPB;
{\[Phi]s,vs,as,bs,Rs}=MultiPB;
Ns = Length[R]-1;
ps = If[pos>1,pos-1,pos];
z = {field1,field2};
extremaV =ListPlot[ {\[Phi][[1,z]],\[Phi][[-1,z]]} ,PlotStyle->{PointSize[.015],Black} ,LabelStyle->Directive[Black,FontSize->20, FontFamily->"Times New Roman",FontSlant->Plain],Frame->True,FrameLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)"},GridLinesStyle->GrayLevel[.85],GridLines->Automatic ];
zPBPath =Show[ Table[  ParametricPlot[
vs[[s,z]]+4/d as[[s,z]]\[Rho]^2+If[pos>1&&s==ps,0.,2/(d-2)bs[[s,z]]/\[Rho]^(d-2)],{\[Rho],Min[Rs[[s]]],Min[Rs[[s+1]]]},PlotStyle->Cyan, LabelStyle->Directive[Black,FontSize->16, FontFamily->"Times New Roman",FontSlant->Plain],AxesLabel->{Subscript["\[CurlyPhi]",field1],Subscript["\[CurlyPhi]",field2]},ImageSize->350,AxesOrigin->{-1,-1}],{s,ps,Ns}],PlotRange->All];
plot =Show[extremaV,zPBPath,ImageSize->400]];


(* ::Section::Closed:: *)
(*PlotBounce3D*)


PlotBounce3D[FindBounce_List]:=Block[{Action,PB,improvementPB,MultiPB,\[Phi]L,VL,\[Phi],pos,v,a,b,R,ddVL,\[Nu],\[Alpha],\[Beta],r,\[Phi]s,vs,as,bs,Rs,list\[Phi]LVL,listR\[Phi]L,d,Ns,imagesize,fontsize,ps,PBPath,extremaV,plot},
(*======================*)
{Action,PB,improvementPB,MultiPB}=FindBounce;
{d,\[Phi]L,VL,\[Phi],pos,v,a,b,R}=PB;
{ddVL,\[Nu],\[Alpha],\[Beta]}=improvementPB;
{\[Phi]s,vs,as,bs,Rs}=MultiPB;
Ns = Length[R]-1;
(*======================*)
imagesize = 400;fontsize =20;
extremaV=ListPointPlot3D[ \[Phi][[{1,-1}]] ,PlotStyle->{Directive[ Orange,PointSize[.02]]},LabelStyle->Directive[Black,FontSize->20, FontFamily->"Times New Roman",FontSlant->Plain],AxesLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)","\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(3\)]\)"},ImageSize->400];
ps = If[pos>1,pos-1,pos];
PBPath =   Show[ Table[  ParametricPlot3D[
vs[[s]]+4/d as[[s]]\[Rho]^2+If[pos>1&&s==ps,0.,2/(d-2)bs[[s]]/\[Rho]^(d-2)],{\[Rho],Min[Rs[[s]]],Min[Rs[[s+1]]]},PlotStyle->Blue, LabelStyle->Directive[Black,FontSize->16, FontFamily->"Times New Roman",FontSlant->Plain],AxesLabel->{"\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)","\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)","\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(3\)]\)"}],{s,ps,Ns}],ListPointPlot3D[    If[ pos>1,Join[{vs[[pos-1]]},\[Phi]s[[pos;;-1]]],\[Phi]s]  ,PlotStyle->{Lighter[Blue,.5],Dashed}],Graphics3D[{Darker[Cyan],Thickness[.005],Line[If[ pos>1,Join[{vs[[pos-1]]},\[Phi]s[[pos;;-1]]],\[Phi]s]]}],PlotRange->All] ;
list\[Phi]LVL = ListPlot[Transpose@{\[Phi]L,VL/100},PlotStyle->Blue,Joined->True,Frame->True,Mesh->Full, LabelStyle->Directive[Black,FontSize->20, FontFamily->"Times New Roman",FontSlant->Plain],FrameLabel->{"\[CurlyPhi]","V(\[CurlyPhi]) \!\(\*SuperscriptBox[\(.10\), \(2\)]\)"},GridLines->Automatic,ImageSize->Medium];
listR\[Phi]L = Show[ 
If[pos==1,Plot[  vs[[pos]]+ 4/d as[[pos]]R[[pos]]^2 +2/(d-2)bs[[pos]]/R[[pos]]^(d-2),{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False,Epilog->{Text[Style[Row[{"Action = ",Action}],10,15,"Times New Roman"],Scaled[{0.75,0.9}],Background->LightGreen]} ],Plot[(vs[[pos-1]]+ 4/d as[[pos-1]]\[Rho]^2 ),{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False,Epilog->{Text[Style[Row[{"Action = ",Action}],10,15,"Times New Roman"],Scaled[{0.75,0.85}],Background->LightGreen]} ]],
(*====================*)
Table[Plot[( vs[[s]] + 4/d as[[s]]\[Rho]^2 +If[\[Rho]>0,2/(d-2)bs[[s]]/\[Rho]^(d-2),0] ),{\[Rho],R[[s]],R[[s+1]]},PlotStyle->RGBColor[1,0.5,1/8]],{s,pos,Ns}],
(*====================*)
Plot[( vs[[Ns]] + 4/d as[[Ns]]Rs[[-1]]^2 +2/(d-2)bs[[Ns]]/Rs[[-1]]^(d-2) ),{\[Rho],R[[-1]],R[[-1]]2},PlotStyle->RGBColor[1,0.5,1/8]],
(*====================*)
Frame->True,
FrameLabel->{"\[Rho]","\[CurlyPhi](\[Rho])"},LabelStyle->Directive[Black,FontSize->17, FontFamily->"Times New Roman",FontSlant->Plain],GridLines->Automatic ,GridLinesStyle->GrayLevel[.85],ImageSize->Medium,PlotRange->All ] ;
plot =Show[extremaV,PBPath];
{plot,list\[Phi]LVL ,listR\[Phi]L }];


(* ::Section::Closed:: *)
(*ListPointBounce*)


ListPointBounce[pts_,d_]:=ListPointBounce[pts,d]=
Block[{Action,PB,dim,improvementPB,MultiPB,\[Phi]L,VL,\[Phi],pos,v,a,b,R,ddVL,\[Nu],\[Alpha],\[Beta],r,\[Phi]s,vs,as,bs,Rs,list\[Phi]LVL,listR\[Phi]L,Ns,imagesize,fontsize,contourV,extremaV,ps,PBPath,V,\[CurlyPhi]V,\[CurlyPhi],plot},
\[CurlyPhi]V = Transpose[pts];
CheckAbort[{Action,PB,improvementPB,MultiPB}=FindBounce[V,\[CurlyPhi],{\[CurlyPhi]V[[1,1]],\[CurlyPhi]V[[1,-1]]},"Dimension"->d,"AnsatzPath"->\[CurlyPhi]V[[1]],"AnsatzV1"->\[CurlyPhi]V[[2]]];
{dim,\[Phi]L,VL,\[Phi],pos,v,a,b,R}=PB;
{ddVL,\[Nu],\[Alpha],\[Beta]}=improvementPB;
{\[Phi]s,vs,as,bs,Rs}=MultiPB;
Ns = Length[R]-1;
plot =Show[ 
If[pos==1,Plot[  Sign[\[Phi][[-1]]-\[Phi][[1]]](  v[[pos]]+ 4/d a[[pos]]R[[pos]]^2 +2/(d-2)b[[pos]]/R[[pos]]^(d-2)-Abs[\[Phi][[-1]]-\[Phi][[1]]])+ \[Phi][[-1]],{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False],Plot[ Sign[\[Phi][[-1]]-\[Phi][[1]]]( (v[[pos-1]]+ 4/d a[[pos-1]]\[Rho]^2 )-Abs[\[Phi][[-1]]-\[Phi][[1]]])+ \[Phi][[-1]],{\[Rho],0.,R[[pos]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False]],
(*====================*)
Table[Plot[{  Sign[\[Phi][[-1]]-\[Phi][[1]]]( ( v[[s]] + 4/d a[[s]]\[Rho]^2 +If[\[Rho]>10^-5,2/(d-2)b[[s]]/\[Rho]^(d-2),0] )-Abs[\[Phi][[-1]]-\[Phi][[1]]])+  \[Phi][[-1]]  },{\[Rho],R[[s]],R[[s+1]]},PlotStyle->RGBColor[1,0.5,1/8],Axes->False],{s,pos,Ns}]
(*====================*)
,Plot[{ Sign[\[Phi][[-1]]-\[Phi][[1]]]( ( v[[Ns]] + 4/d a[[Ns]]R[[-1]]^2 +2/(d-2)b[[Ns]]/R[[-1]]^(d-2) )-Abs[\[Phi][[-1]]-\[Phi][[1]]])+ \[Phi][[-1]]  },{\[Rho],R[[-1]],R[[-1]]2},PlotStyle->RGBColor[1,0.5,1/8],Axes->False]
(*====================*)
,ListPlot[Join[If[pos==1,{},{{0, Sign[\[Phi][[-1]]-\[Phi][[1]]]( v[[pos-1]]-Abs[\[Phi][[-1]]-\[Phi][[1]]])+ \[Phi][[-1]]}}],Transpose[{R[[pos;;-1]], \[Phi][[pos;;-1]]}]],PlotStyle->{RGBColor[1,0.5,1/8],PointSize[.017]},Axes->False]
,Frame->True,
FrameLabel->{"\[Rho]","\[CurlyPhi](\[Rho])"},LabelStyle->Directive[Black,FontSize->17, FontFamily->"Times New Roman",FontSlant->Plain],GridLines->Automatic ,GridLinesStyle->GrayLevel[.85],ImageSize->Medium,PlotLabel->Style[Row[{"Action = ",Action}],Background->LightGreen,FontSize->17],PlotRange->{{0,10},{Min[\[Phi]],Max[\[Phi]]}}  ]; ,
plot =ListPlot[{{0,0}},PlotStyle->Red,Epilog->{Inset[Text[Style["$Aborted",Background->LightRed,FontSize->15, FontFamily->"Times New Roman",FontSlant->Plain]],Scaled[{.8,.8}]]}];  ];
Return[plot]];


(* ::Chapter::Closed:: *)
(*End Package*)


End[]; (*"`Private`"*)

EndPackage[];
