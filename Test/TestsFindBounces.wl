(* ::Package:: *)

(* ::Section:: *)
(*Begin Test Section*)


(* ::Subsection:: *)
(*Description*)


(* Test of Functions*)

(* "FindBounce" package must be loaded before running these tests, otherwise testing is aborted. *)


If[
	Not@MemberQ[$Packages,"FindBounces`"],
	Print["Error: Package is not loaded, try again"];Abort[];
];
BeginTestSection["Tests"];

Begin["Private`"];


(* ::Subsection:: *)
(*Global Variables*)


\[Mu]\[Lambda]\[Alpha]\[Epsilon]={{80.,100.},{.1,.3},2.,{0.,0.}} ;
V2[\[Phi]_]:=Sum[-\[Mu]\[Lambda]\[Alpha]\[Epsilon][[1,i]] (\[Phi][[i]])^2 +\[Mu]\[Lambda]\[Alpha]\[Epsilon][[2,i]](\[Phi][[i]])^4+\[Mu]\[Lambda]\[Alpha]\[Epsilon][[4,i]]\[Phi][[i]] ,{i,1,2}]+ \[Mu]\[Lambda]\[Alpha]\[Epsilon][[3]](\[Phi][[1]])^2 (\[Phi][[2]])^2; 
dV[\[CurlyPhi]_]:={-160. \[CurlyPhi][[1]]+0.4 \[CurlyPhi][[1]]^3+4 \[CurlyPhi][[1]] \[CurlyPhi][[2]]^2,-200. \[CurlyPhi][[2]]+4.\[CurlyPhi][[1]]^2 \[CurlyPhi][[2]]+1.2 \[CurlyPhi][[2]]^3};
d2V[\[CurlyPhi]_]:={{-160. +1.2 \[CurlyPhi][[1]]^2+4. \[CurlyPhi][[2]]^2,8. \[CurlyPhi][[1]] \[CurlyPhi][[2]]},{8. \[CurlyPhi][[1]] \[CurlyPhi][[2]],-200.+4. \[CurlyPhi][[1]]^2+3.6 \[CurlyPhi][[2]]^2}};
Get@FileNameJoin[{NotebookDirectory[],"variablesTest.wl"}];
If[
	Head[variableTest1]=!=List||Head[variableTest1]=!=List,
	Print["Error: Needs variablesTest.wl in the NotebookDirectory[], try again"];Abort[];
];
{aRw,\[Phi]N3,a,R,v,b,d,\[Phi]L,VL,dVL,ddVL,\[Alpha],\[ScriptCapitalI],d\[ScriptCapitalI],r,\[Beta],\[Nu],
r1,Ns,pos,forBack,timeRw,\[Phi]t,\[Phi],N\[CurlyPhi],l,eL,\[Phi]s,vs,as,bs,Rs,Ts,
Vs,T\[Xi],V\[Xi],T1,V1,Rw,methodRw,maxIteR,accuracyB,ansatzRw,accuracyPath,ite\[Zeta],methodSeg,Nfv,aV1} = variableTest1;
{\[Phi]s2,vs2,as2,bs2,Rs2,\[Phi]2,\[Phi]L2,eL2,l2,v2,a2,b2,R2,VL2,ddVL2,\[ScriptCapitalI]2,d\[ScriptCapitalI]2,\[Alpha]2,pos2,r2,\[Beta]2,\[Nu]2,r12,T\[Xi]2,V\[Xi]2} = variableTest2; 


(* ::Subsection:: *)
(*VerificationTest*)


VerificationTest[Round[Segmentation[{1,2,6},numberFieldValues->4],0.01],{1.,2.67,4.33,6.},
	TestID->"Segmentation"
];
VerificationTest[findSegment[a,\[Phi]L,d,Ns],pos,
	TestID->"findSegment"
];
(*VerificationTest[ ,
	TestID->"ansatzN3"
];*)
VerificationTest[Round[newAnsatz[\[Phi]s[[-1]],Ns,N\[CurlyPhi]],0.01],Round[{Ns,\[Phi],\[Phi]L,eL,l},0.01],
	TestID->"newAnsatz"
];
VerificationTest[Round[Re[Rvb[Rw,\[Phi]L,a,d,Ns,forBack,pos]],0.01],Round[{R,v,b},0.01],
	TestID->"Rvb"
];
VerificationTest[Round[findRw[d,VL,\[Phi]L,a,Ns, methodRw,maxIteR,accuracyB,ansatzRw,aRw][[-1]],0.0001],Round[Rw,0.0001],
	TestID->"findRw"
];
VerificationTest[Round[find\[ScriptCapitalI][v2,\[Phi]L2,a2,b2,Ns,pos,R2,ddVL2,d],0.01],Round[{\[ScriptCapitalI]2,d\[ScriptCapitalI]2},0.01],
	TestID->"find\[ScriptCapitalI]"
];
VerificationTest[Round[r\[Beta]\[Nu][r12,a2,b2,d,Ns,\[Alpha]2,R2,\[ScriptCapitalI]2,d\[ScriptCapitalI]2,pos2],0.01][[All,2;;-1]],Round[{r2,\[Beta]2,\[Nu]2},0.01],
	TestID->"r\[Beta]\[Nu]"
];
VerificationTest[Round[rw/.FindRoot[findrw[rw,a2,b2,d,Ns,\[Alpha]2,R2,\[ScriptCapitalI]2,d\[ScriptCapitalI]2,pos2],{rw,-1}],0.01],Round[r12,0.01],
	TestID->"findrw"
];
VerificationTest[Round[\[Phi]vabRs[dV,d2V,\[Phi],eL,l,\[Phi]L,v,a,b,Ns,R,N\[CurlyPhi],pos,d,ite\[Zeta],accuracyPath],0.01],Round[{\[Phi]s2,vs2,as2,bs2,Rs2},0.01],
	TestID->"\[Phi]vabRs"
];
VerificationTest[Round[PathDeformation[dV,d2V,Ns,N\[CurlyPhi],pos,d,Rs[[1]],\[Phi]s[[1]],vs[[1]],as[[1]],bs[[1]],accuracyPath],0.01],Round[{\[Phi]s[[2]]-\[Phi]s[[1]],vs[[2]]-vs[[1]],as[[2]]-as[[1]],bs[[2]]-bs[[1]],Rs[[2,pos]]-Rs[[1,pos]],Rs[[2,-1]]-Rs[[1,-1]]},0.01],
	TestID->"PathDeformation"
];
VerificationTest[Round[\[ScriptCapitalT][v,a,b,R,d,Ns],0.01],Round[T1,0.01],
	TestID->"\[ScriptCapitalT]"
];
VerificationTest[Round[\[ScriptCapitalV][v,a,b,R,d,VL,\[Phi]L,Ns],0.01],Round[V1,0.01],
	TestID->"\[ScriptCapitalV]"
];
VerificationTest[Round[\[ScriptCapitalT]\[Xi][d,a2,R2,b2,v2,\[Alpha]2,\[Beta]2,\[Nu]2,ddVL2,VL2,\[Phi]L2,Ns,pos2,r2],0.01],Round[T\[Xi]2,0.01],
	TestID->"\[ScriptCapitalT]\[Xi]"
];
VerificationTest[Round[\[ScriptCapitalV]\[Xi][d,a2,R2,b2,v2,\[Alpha]2,\[Beta]2,\[Nu]2,ddVL2,VL2,\[Phi]L2,Ns,pos2,r2],0.01],Round[V\[Xi]2,0.01],
	TestID->"\[ScriptCapitalV]\[Xi]"
];


(* ::Subsection:: *)
(*EndTestSection*)


End[];

EndTestSection[];
