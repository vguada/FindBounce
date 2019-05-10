(* ::Package:: *)

(* ::Subsection:: *)
(*Description*)


(* "FindBounce" package must be loaded before running tests, otherwise procedure is aborted. *)
If[
	Not@MemberQ[$Packages,"FindBounce`"],
	Print["Error: Package is not loaded, try again"];Abort[];
];


(* ::Subsection:: *)
(*Begin Test Section*)


BeginTestSection["Tests"];


(* ::Subsection:: *)
(*End Test Section*)


EndTestSection[];
