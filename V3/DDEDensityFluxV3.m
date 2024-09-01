
(* File: DDEDensityFluxV3.m *)
(* Third version of DDEDensityFlux.m *)

(* Last updated: Wednesday, August 5, 2009, at 17:05, by WH in Boulder *)

(* Previously updated: Saturday, October 4, 2003, at 22:00, by WH in Boulder *)

(* Blended version of DHolly0610.m and DDEDensityFluxV2.m *)
(* All changes signed by WH 10/04/2003 and WH 10/05/2003 *)

(* DDEDensityFluxV2.m : *)
(* Includes the latest version of the discrete homotopy operator *)
(* Code inserted in this file, see DISCRETE HOMOTOPY OPERATOR *)
(* Code works under Mathematica 5.0 and earlier versions *)
(* Also minor changes in the print statements *)
(* All changes signed by WH 10/04/2003 *)

(* DHolly0610.m: *)
(* Last updated: Tuesday, June 10, 2003 at 16:30 by WH at home *)
(* Sheltered some local variables a bit better. See WH 06/10/03 *)
(* Added more cases to the menu, including: *)
(* Blaszak-Marciniak three field lattices *)
(* Blaszak-Marciniak four field lattices *)
(* Blaszak-Marciniak higher-order three field lattice Case I and II *)
(* Belov-Chaltikian lattice *)

(* History DDEDensityFluxV2.m:                                              *)
(* Time-stamp: <Tue Jul  1 06:27:34 MDT 2003>                               *)
(* File: DDEDensityFlux.m                                                   *)
(* Last updated: July 1, 2003 at 05:25 by Mike C.                           *)
(* July 1, 2003: Integration of the Homotopy Operator code.  Search for     *)
(* "discreteHomotopyOperator", "Integration", and "Mike C." for changes     *)

(****************************************************************************)
(*                                                                          *)
(*          *** M A T H E M A T I C A   P R O G R A M ***                   *)
(*                                                                          *)
(*                         DDEDensityFluxV3.m                               *)
(*                                                                          *)
(*        SYMBOLIC COMPUTATION of CONSERVED DENSITIES and FLUXES            *)
(*   for SYSTEMS of EVOLUTIONARY TYPE DIFFERENTIAL-DIFFERENCE EQUATIONS     *)
(*                                                                          *)
(* program name: DDEDensityFluxV3.m                                         *)
(*                                                                          *)
(* purpose: computation of conserved densities and fluxes with possible     *)
(*          compatibility conditions                                        *)
(*                                                                          *)
(* features: use multiscale approach, and possible shifts in density,       *)
(*           linear system for the coefficients can be computed via a       *)
(*           a shifting algorithm or by applying the discrete Euler         *)
(*           operator (variational derivative)                              *)
(*                                                                          *)
(* input: system of evolution type differential difference equations of     *)
(*        any order, any degree, polynomial type, only constant parameters, *)
(*        variable coefficients are NOT allowed as parameters               *)
(*                                                                          *)
(* output: density and flux of desired rank (if it exists),                 *)
(*         and compatibility conditions for the parameters (if applicable)  *)
(*                                                                          *)
(* tested: on PC's (desktops and laptops) running Pentium II and III        *)
(*                                                                          *)
(* language: Mathematica 5.0 and also v. 4.0 and 4.1                        *)
(*                                                                          *)
(* authors: Unal Goktas, Willy Hereman, and Holly Eklund                    *)
(*          Department of Mathematical and Computer Sciences                *)
(*          Colorado School of Mines                                        *)
(*          Golden, CO 80401-1887, USA                                      *)
(*                                                                          *)
(* email:whereman@mines.edu                                                 *)
(*                                                                          *)
(* Version 3: October 5, 2003                                               *)
(*                                                                          *)
(* Copyright 1998-2003                                                      *)
(*                                                                          *)
(****************************************************************************)

Clear["Global`*"];

(* FLAGS FOR FORCING OPTIONS *)
(* Force options must be given in the data file *)

(* If forcesinglescale=True, only scale1 with w(d/dt)=1 is used to  *) 
(* construct the form of rho. The density is not split into pieces. *)

(* If forcemultiplescale=True, rho is constructed based on scale1,  *)
(* but split into pieces according to scale0 with w(d/dt) = 0.      *)
(* NOTE: If both are false, the code uses scale0 when appropriate.  *)

(* forcesinglescale   = False; *)
(* forcemultiplescale = False; *)

(* If forceshiftsinrho=True, the code will generate the form of rho *)
(* based on {u_n, u_{n+1}, ..., u_{n+maximumshift}}.                *)
(* To compute the maximumshift there are two options (see below).   *)
(* If forceshiftsinrho=False, no shifts on the dependent variable   *)
(* u_n will be used to generate the form of rho.                    *)

(* forceshiftsinrho = False; *)

(* If forcemaximumshiftrhsdde=True, then the maximumshift is based  *)
(* on the maximum of the shifts occuring in the rhs of the DDEs.    *)

(* If forcemaximumshiftpowers=True, then maximumshift = p-1 where   *)
(* p = (rhorank/weightu[1]), that is the highest power in rho.      *)
(* Rho will be generated from the list {u_n, u_{n+1},...,u_{n+p-1}}.*)
(* If u_n^p is the highest-power term in rho, then e.g. the term    *)
(* u_n u_{n+1} u_{n+2} ... u_{n+p-1} will occur in the form of rho. *)
(* For example, if p=4 then maximumshift = p-1 =3. The terms in rho *)
(* are based on {u_n, u_{n+1}, u_{n+2}}. So, rho has terms like     *) 
(* u_n^4, u_n^2 u_{n+3}^2, ..., u_n u_{n+1} u_{n+2} u_{n+3}.        *)
(* NOTE: forcemaximumshiftrhsdde and forcemaximumshiftpowers must   *)
(* opposite Boolean values. *)

(* forcemaximumshiftrhsdde = False; *)
(* forcemaximumshiftpowers = False; *)

(* If forcediscreteeuler=True, then the linear system for the c[i] is *)
(* computed via the Discrete Variational Derivative (Euler) Operator. *) 
(* If forceshifting=True, then the original shifting algorithm is     *)
(* used to compute the linear system for the c[i] and the fluxes J_n. *)

(* forcediscreteeuler = True; *)
(* forceshifting      = False; *)

(* If forcestripparameters=True, then power of the nonzero parameters, *)
(* are removed in the factored form of the linear system for the c[i]. *) 
(* Example: aa, bb^3, aa*bb, ... are removed during simplification,    *)
(* but not factors like (aa^2-bb), (aa+bb)^4, etc. *)

(* forcestripparameters = True; *)

(* If forceextrasimplifications=True, then while generating the linear *)
(* system, equations of type c[i] == 0 are automatically applied to    *)
(* the rest of the system for the c[i]. *)

(* forceextrasimplifications = True; *)

(* DEBUGGING FLAGS *)

printflagcommonfactor     = False;
debugsubroutine4          = False;
debugsubroutine5          = False;
debugmaxshift             = False;
debugconstructoriginalrho = False;
debugconstructformrho     = False; 
debugcheckpoints          = False;
debugsubscriptform        = False;
debugformrho              = False;
debugformrhoL             = False;
debugweightsK             = False;
debugweightsL             = False;
debugweightsM             = False;
debugweightsX             = False;
debugweightsY             = False;
debugweightsZ             = False;
debugconstructformrhoK    = False;
debugconstructformrhoOM   = False;
debugconstructformrhoNM   = False;
debugroutine5N            = False;
debugdiscrete             = False;
debugdeulerd              = False;
debugdiscreteeuler        = False;
debugconstructformjn      = False;
debugevaluate             = False;
debugconstructeqlist      = False;
debugstripsystem          = False;
debugstripparameters      = False;
debugmysimplify1          = False;
debugmysimplify2          = False;
debugstripadjust          = False;
debugsystemsolver         = False;
debuganalyzer             = False;

(* ######################## B1 ################################# *)
(*****************************************************************************)
(* commentinter[]: prints a welcome message                                  *)
(*****************************************************************************)
(* WH 10/04/2003 Version is now 2.0, dated October 4, 2003 *)
commentinter[] := Block[{},
Print["*********************************************************"];
Print["  WELCOME TO THE MATHEMATICA PROGRAM DDEDensityFlux.m    "];
Print["     by UNAL GOKTAS, WILLY HEREMAN, AND HOLLY EKLUND     "];
Print[" FOR THE COMPUTATION OF CONSERVED DENSITIES AND FLUXES.  "];
Print["      Version 2.0 released on October 4, 2003            "];
Print["          Last updated on August 5, 2009                 "];
Print["             Copyright 1998-2003                         "];
Print["*********************************************************"]
]; (* end block commentinter *)
(* ######################## E1 ################################# *)

(* ######################## B2 ################################# *)
(*****************************************************************************)
(* cls: clears the screen                                                    *)
(*****************************************************************************)
cls := Print["\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"];
(* ######################## E2 ################################# *)

(* ######################## B3 ################################# *)
(*****************************************************************************)
(* printpage[]: a subroutine for menu                                       *)
(*****************************************************************************)
printpage[n_,page_] := Module[{i,choice,lenpage},
(* cls taken out to save space and paper *)
(* cls; *)
lenpage = Length[page[n]];
Print[" "];
Print["  *** MENU INTERFACE ***  (page: ",n,")"];
Print["-------------------------------------------"];
For[i=1,i <= lenpage,i++,
    Print[Part[Part[page[n],i],1]]];
Print[" nn) Next Page"];
Print[" tt) Your System"];
Print[" qq) Exit the Program"];
Print["------------------------------------------"];
choice = Input["ENTER YOUR CHOICE: "];
Return[choice]
]; (* end module printpage *)
(* ######################## E3 ################################# *)

(* ######################## B4 ################################# *)
(*****************************************************************************)
(* menu: creates the menu                                                    *)
(*****************************************************************************)
menu := Module[
{counterpage = 1,menulist,numpages,page,
choice1,control,choice2,choice3,lenmenulist,i},
menulist = {
{"  1) Kac-van Moerbeke (or Volterra) Lattice (d_kdv.m)"},
{"  2) Modified KdV Lattice (with parameter) (d_mkdv.m)"},
{"  3) Modified (quadratic) Volterra Lattice (d_modvol.m)"},
{"  4) Ablowitz-Ladik Discretization of NLS Equation (d_ablnl1.m)"},
{"  5) Toda Lattice (d_toda.m)"},
{"  6) Standard Discretization of NLS Equation (d_stdnls.m)"},
{"  7) Herbst/Taha Discretization of KdV Equation (d_diskdv.m)"},
{"  8) Herbst/Taha Discretization of mKdV Equation (d_dimkdv.m)"}, 
{"  9) Herbst/Taha Discretization of combined KdV-mKdV Equation (d_herbs1.m)"},
{"  10) Self-Dual Network Lattice (d_dual.m)"}, 
{"  11) Belov-Chaltikian Lattice (d_bc1.m)"}, 
{"  12) Bogoyavlenskii Lattice (d_bogoya.m)"}, 
{"  13) Blaszak-Marciniak Three Field Lattice I (d_blmc1.m)"}, 
{"  14) Blaszak-Marciniak Three Field Lattice II (d_blmc5.m)"}, 
{"  15) Blaszak-Marciniak Four Field Lattice I (d_blmc2.m)"}, 
{"  16) Blaszak-Marciniak Four Field Lattice II (d_blmc6.m)"}, 
{"  17) Blaszak-Marciniak Four Field Lattice III (d_blmc7.m)"}, 
{"  18) Blaszak-Marciniak Higher-order Three Field Lattice I (d_blmc3.m)"}, 
{"  19) Blaszak-Marciniak Higher-order Three Field Lattice II (d_blmc4.m)"}, 
{"  20) Non-isospectral Blaszak-Marciniak Three Field Lattice (d_blmc1n.m)"}, 
{"  21) Parameterized Toda Lattice (d_ptoda.m)"},
{"  22) Generalized Toda Lattice-1 (d_gtoda1.m)"},
{"  23) Generalized Toda Lattice-2 (d_gtoda2.m)"}, 
{"  24) Henon System (d_henon.m)"} 
}; (* closes menulist *)
lenmenulist = Length[menulist];
numpages = Ceiling[lenmenulist/10];
For[i = 1,i <= numpages,i++,
   page[i] = If[lenmenulist >= (i*10),
               menulist[[Table[k,{k,(i-1)*10+1,i*10}]]],
               menulist[[Table[k,{k,(i-1)*10+1,lenmenulist}]]]]];

choice1 = printpage[counterpage,page];

control := (
Switch[choice1,
tt,Print["Make sure that you have prepared the data file for the system"];
   Print["you want to test (similar to the data files we supplied)."];
   choice2 = Input["If your file is ready, press 1, else 2: "];
   If[choice2 === 1,
    choice3 = InputString["Enter the name of your data file: "];
    Get[choice3],
    Print["Aborting the computations!"];
    Print["Prepare your data file first, then start over."];
    Abort[]
     ],
nn,If[counterpage < numpages,counterpage++;
     choice1 = printpage[counterpage,page]; control,
     counterpage = 1; choice1 = printpage[1,page]; control],
qq,Print["Aborting the computations!"];Abort[],
1,<<d_kdv.m,
2,<<d_mkdv.m,
3,<<d_modvol.m,
4,<<d_ablnl1.m,
5,<<d_toda.m,
6,<<d_stdnls.m,
7,<<d_diskdv.m,
8,<<d_dimkdv.m,
9,<<d_herbs1.m,
10,<<d_dual.m,
11,<<d_bc1.m, 
12,<<d_bogoya.m,
13,<<d_blmc1.m, 
14,<<d_blmc5.m, 
15,<<d_blmc2.m, 
16,<<d_blmc6.m, 
17,<<d_blmc7.m, 
18,<<d_blmc3.m, 
19,<<d_blmc4.m, 
20,<<d_blmc1n.m,
21,<<d_ptoda.m,
22,<<d_gtoda1.m,
23,<<d_gtoda2.m,
24,<<d_henon.m,
_,Print["Aborting the computations!"];Abort[]
]; (* closes Switch *)
); (* control *)

control;
If[Not[StringQ[myfile]],
   myfile = InputString["Enter the name of the output file: "]]; 
]; (* end module menu *)
(* ######################## E4 ################################# *)

(* ######################## B5 ################################# *)
(****************************************************************************)
(* upShift[f_]: takes variable u[i][n,t] and shifts the subscript up        *)
(*                  one place.  u[i][n,t]->u[i][n+1,t]                      *)
(****************************************************************************)
upShift[f_] := Module[{nf},
      nf = f /. {n -> n + 1};
      Return[nf]
]; (* end module upShift *)
(* ######################## E5 ################################# *)

(* ######################## B6 ################################# *)
(*****************************************************************************)
(* downShift[f_]:takes variable u[i][n,t] and shifts the subscript down      *)
(*                  one place.  u[i][n,t]->u[i][n-1,t]                       *)
(*****************************************************************************)
downShift[f_] := Module[{nf},
        nf = f /. {n -> n - 1};
        Return[nf]
]; (* end module downShift *)
(* ######################## E6 ################################# *)

(* ######################## B7 ################################# *)
(*****************************************************************************)
(* subscriptform[]: takes an expression and prints it with subscripted form  *)
(*                  and removes the functional dependence                    *)
(*****************************************************************************)
subscriptform[expr_] := Module[{tempexpres},
tempexpres = expr;
If[debugsubscriptform, 
   Print["At Pt. PDE1, ENTERING is tempexpres= "];
   Print[tempexpres]
   ];
tempexpres = tempexpres /. {u[i_][n_][x__] :>
             SequenceForm[u,Subscript[i],Subscript[","],Subscript[n]]};
If[debugsubscriptform, 
   Print["At Pt. PDE2, LEAVING is tempexpres= "];
   Print[tempexpres]
   ];
Return[tempexpres]
]; (* end module subscriptform *)
(* ######################## E7 ################################# *)

(* ######################## B8 ################################# *)
(*****************************************************************************)
(* splitrho[]: if free coefficients c[i] remain in the density, density      *)
(*          is further split in independent pieces corresponding to each     *)
(*          coefficient c[i]                                                 *)
(*****************************************************************************)
splitrho[expr_] := Module[{listofc,lenlistofc,ci,rhoi,rest,finalformrest}, 
  rest=Expand[expr];
  If[FreeQ[rest,c[__]], (* start if 12, then part *)
If[debugcheckpoints,
  Print["Passing checkpoint 1000, then part of if 12"]
  ];
  Print["The density has no free coefficients"];
(* WH 10/05/2003 *)
(*  Print["The density is: "]; *)
(*  Print[rest]; *)
   rest = Factor[rest];
   rest = Numerator[rest];
   finalformrho = subscriptform[rest];
       Print["Saving the density in an output file."];
       Save["DDEdf.out", finalformrho], (* else for if 12 *)
If[debugcheckpoints,
  Print["Passing checkpoint 1000, else part of if 12"]
  ];
    listofc=Union[Cases[rest,c[__],12]];
    lenlistofc=Length[listofc];
    Print[" "];
    Print["There is/are ", lenlistofc, " free coefficient(s) in the density."];
    Print["These free coefficients are: ", listofc, "."];
    Print["Splitting the density in independent pieces."];
     Do[ (* begin do 13 *)
       ci=Part[listofc,i];
       rhoi=Expand[Coefficient[rest,ci]]; 
       rest=Expand[rest-ci*rhoi];
       Print[" "];
       Print["Part of the density rho with coefficient ", ci, ":"];
       Print[" "];
       Print[rhoi],{i,1,lenlistofc}]; (* end do 13 *)
     Print[" "];
     Print["Part of the density rho that had no free coefficient: "];
     Print[" "];
     Print[rest];
     Print[" "] 
] (* end if 12 *) 
]; (* end module splitrho *)
(* ######################## E8 ################################# *)

(* ######################## B9 ################################# *)
(*****************************************************************************)
(* splitflux[]: if free coefficients c[i] remain in the flux, the flux       *)
(*          is further split in independent pieces corresponding to each     *)
(*          coefficient c[i]                                                 *)
(*****************************************************************************)
splitflux[expr_] := Module[{listofc,lenlistofc,ci,fluxi,rest}, 
  rest=Expand[expr];
  If[FreeQ[rest,c[__]], (* start if 14 *)
    Print["The flux has no free coefficients."];
(* WH 10/05/2003 *)
(*    Print["The flux is: "]; *)
(*    Print[rest]; *)
    finalformflux = subscriptform[rest];
       Print["Saving the flux in an output file."];
       Save["DDEdf.out", finalformflux], (* else of if 14 *)
    listofc=Union[Cases[rest,c[__],12]];
    lenlistofc=Length[listofc];
    Print[" "];
    Print["There is/are ", lenlistofc, " free coefficient(s) in the flux."];
    Print["These free coefficients are: ", listofc, "."];
    Print["Splitting the flux in independent pieces."];
     Do[
       ci=Part[listofc,i];
       fluxi=Expand[Coefficient[rest,ci]]; 
       rest=Expand[rest-ci*fluxi];
       Print[" "];
       Print["Part of the flux J_n with coefficient ", ci, ":"];
       Print[" "];
       Print[fluxi],{i,1,lenlistofc}]; (* end do *)
     Print[" "];
     Print["Part of the flux J_n that had no free coefficient: "];
     Print[" "];
     Print[rest];
     Print[" "] 
] (* end if 14 *) 
]; (* end module splitflux *)
(* ######################## E9 ################################# *)

(* ######################## B10 ################################# *)
(*************************************************************************** *)
(* stripper of nonzero parameters, used in the simplification of the linear  *)
(* system for the c[i].                                                      *)
(* Nonzero parameters is the union of parameters and weightedparameters.     *)
(*************************************************************************** *)
stripparameters[equatlist_List,param_List]:=Module[{tempequatlist,result,rules,
listcoefstripped={},coefk,parttempequatlistk,partresultk},
If[param =!= {}, (* start if 1000 *)
 tempequatlist=equatlist;
 If[debugstripparameters, 
  Print["At Pt. STR0, inside stripparameters."];
  Print["Starting removal of nonzero parameters in the system for the c[i]."];
  Print["The unsimplified system has ",Length[tempequatlist]," equations."]
  ];
 If[debugstripparameters, 
   Print["At Pt. STR1, entering stripparameters, tempequatlist ="];
   Print[tempequatlist]
   ];
 If[debugstripparameters, 
   Print["At Pt. STR2, entering stripparameters, param ="];
   Print[param]
   ];
 (* Sets up rules to be applied *)
 rules = #^_ : 1 :> Sequence[] & /@ param;
 If[debugstripparameters, 
   Print["At Pt. STR3, rules ="];
   Print[rules]
   ];
 (* Maps factor to every term, so as to have a constant base. *)
 result = MapAll[Factor, tempequatlist];
 If[debugstripparameters, 
   Print["At Pt. STR4, factoring tempequatlist, result ="];
   Print[result]
   ];
 (* If it is x*b it converts it to {x, b} *)
 result =  If[Head[#] === Times, List @@ #, {#}] & /@ result;
 If[debugstripparameters, 
   Print["At Pt. STR5, making list of coefficients, result ="];
   Print[result]
   ];
 (* Removes parameters to powers: b^2, a, ... *)
 (* apply to depth level 2 only as to not remove e.g. (a^2 - b) *)
 result = Replace[result, rules, 2];
 If[debugstripparameters, 
   Print["At Pt. STR6, after applying the rules, result ="];
   Print[result]
   ];
 (* Puts it back into standard form, {a, b} to a*b *)
 result = Times @@ # & /@ result; 
 If[debugstripparameters, 
   Print["At Pt. STR7, back to products, result ="];
   Print[result]
   ];
 (* test by making a list of coefficients that have been stripped off *)
 Do[ (* start do 11 *)
    parttempequatlistk = Part[tempequatlist,k];
    If[debugstripparameters, 
       Print["At Pt. STRx0, in do loop, for k=",k," parttempequatlistk ="];
       Print[parttempequatlistk]
       ]; 
    partresultk = Part[result,k];
    If[debugstripparameters, 
       Print["At Pt. STRx1, in do loop, for k=",k," partresultk ="];
       Print[partresultk]
       ]; 
    coefk = Factor[parttempequatlistk/partresultk]; 
    If[debugstripparameters, 
       Print["At Pt. STRx2, in do loop, for k=",k," coefk ="];
       Print[coefk]
       ]; 
 (* reporting on accidentally stripping of a-b, and the like *)
    If[Head[coefk]==Plus, 
       Print["Accidentally stripped coefficient ", coefk]
       ];
 (* reporting on accidentally stripping of (a-b)^5, and the like *)
    If[Head[coefk]==Power,  
       If[Head[Part[coefk,1]]==Plus,  
          Print["Accidentally stripped coefficient ", coefk]
         ]
    ];
 (* reporting on accidentally stripping of c[1], c[1]^2, and the like *)
    If[FreeQ[coefk,c[i_]]===False,  
       Print["Accidentally stripped coefficient ", coefk]
       ];
    listcoefstripped = Append[listcoefstripped, coefk];
    If[debugstripparameters, 
       Print["At Pt. STR8, in do loop, for k=",k,", listcoefstripped ="];
       Print[listcoefstripped]
       ], {k,1,Length[tempequatlist]}
    ]; (* end do 11 *)
 listcoefstripped = Union[listcoefstripped];
 If[debugstripparameters, 
    Print["At Pt. STR9, list all coefficients stripped, listcoefstripped ="];
    Print[listcoefstripped]
    ]; 
 result = Union[result];
 (* result = Sort[result, Length[Part[#1,1]] < Length[Part[#2,1]]&]; *)
 If[debugstripparameters, 
  Print["At Pt. STR10, inside stripparameters."];
  Print["Finished removal of nonzero parameters in the system for the c[i]."];
  Print["The simplified system has ",Length[result]," equations."]
  ], (* else *)
 result = equatlist]; (* end if 1000 *)
 Return[result];
]; (* end module stripparameters *)
(* ######################## E10 ################################# *)

(* ######################## B11 ################################# *)
(*************************************************************************** *)
(* stripadjust[]: cancels the numerical factors of each term of the argument *)
(* and puts these into a list after cancellation                             *)
(* If called with flag2=3, it also computes the obstruction for the shifting *)
(* routine (does more than stripping since it modifies the labels n also.    *)
(*****************************************************************************)
stripadjust[givenexpr_,flag1_,flag2_] := Module[
{lenexpr1,list={},expr1,expr2,expr3,expr3s,coefficient,
expr4,lenexpr3,iexpr,k,s},
If[debugstripadjust, 
  Print["At Pt. STR0, givenexpr="];
  Print[givenexpr]
  ];
If[debugstripadjust, 
   Print["At Pt. STR1, flag1="];
   Print[flag1]
  ];
If[debugstripadjust, 
   Print["At Pt. STR2, flag2="];
   Print[flag2]
  ];
expr1 = Expand[givenexpr];
If[debugstripadjust, 
   Print["At Pt. STR3, expanded given expression, expr1="];
   Print[expr1]
   ];
If[Head[expr1] === Plus, lenexpr1 = Length[expr1], lenexpr1 = 1];
If[debugstripadjust, 
   Print["At Pt. STR4, Head[expr1]="];
   Print[Head[expr1]]
   ];
If[debugstripadjust, 
   Print["At Pt. STR5, length of expr1, lenexpr1="];
   Print[lenexpr1]
   ];
If[debugstripadjust, 
   Print["At Pt. STR6, starting the for loop 100."]
   ];
For[k = 1,k <= lenexpr1, k++, (* start for 100, loop iterator is k *)
   If[lenexpr1 == 1, (* start if 35, then part *)
      expr2 = expr1;
      If[debugstripadjust, 
         Print["At Pt. STR6a, then since lenexpr1=1, for k=",k,", expr2="];
         Print[expr2]
         ], (* else if 35 *)
      expr2 = Part[expr1,k];
      If[debugstripadjust, 
         Print["At Pt. STR6b, then since lenexpr1 not 1,for k=",k,",expr2="];
         Print[expr2]
         ]
      ]; (* end if 35 *)
   coefficient = 1;
   If[debugstripadjust, 
      Print["At Pt. STR7, for k=",k,",coefficient="];
      Print[coefficient]
      ];
   expr3 = FactorList[expr2];
   If[debugstripadjust, 
      Print["At Pt. STR8, for k=",k,", using factorlist on expr2, expr3="];
      Print[expr3]
      ];
   lenexpr3 = Length[expr3];
   If[debugstripadjust, 
      Print["At Pt. STR9, for k=",k,", the length of expr3 is lenexpr3="];
      Print[lenexpr3]
      ];
If[debugstripadjust, 
   Print["At Pt. STR10, starting the for loop 200."]
   ];
   For[s = 1, s <= lenexpr3, s++, (* start for loop 200, loop iterator is s *)
      expr3s = Part[expr3,s];
   If[debugstripadjust, 
      Print["At Pt. STR11, for k=",k," and s=",s," part s of expr3, expr3s="]; 
      Print[expr3s]
      ];
      If[FreeQ[expr3s,u], (* start if 50, then part *)
         If[debugstripadjust, 
            Print["At Pt. STR12, for k=",k," and s=",s," then part if 50"]; 
            Print["FreeQ[expr3s,u] should be true, check:",FreeQ[expr3s,u]]
            ];
           If[flag1, (* start if 60, then *)
              If[debugstripadjust, 
               Print["At Pt. STR13, for k=",k," and s=",s," then part if 60"]; 
               Print["do because flag1 should be true, check:",flag1]
                ];
              If[FreeQ[weightpars,Part[expr3s,1]], (* start if 70, then part *)
              If[debugstripadjust, 
               Print["At Pt. STR14, for k=",k," and s=",s," then part if 70"]; 
         Print["do because FreeQ[weightpars,Part[expr3s,1]] should be true,"<>
                     " check:",FreeQ[weightpars,Part[expr3s,1]]]
               ];
(* added or flag2==4 *)
                  If[((flag2 == 3) || (flag2 == 4)), 
                     (* start if 80, then part  *)
              If[debugstripadjust, 
              Print["At Pt. STR15, for k=",k," and s=",s," then part if 80"]; 
              Print["do because flag1=True, and flag2 is 3 or 4, check:",flag2]
              ];
                    coefficient = coefficient*Part[expr3s,1]^Part[expr3s,2];
              If[debugstripadjust, 
              Print["At Pt. STR16, for k=",k," and s=",s," coefficient="]; 
              Print[coefficient]
              ]
                    ]; (* end if 80 *)
                 expr2 = expr2/Part[expr3s,1]^Part[expr3s,2];
              If[debugstripadjust, 
                 Print["At Pt. STR17, for k=",k," and s=",s," "<>
                       " still in then part if 70, expr2=", expr2]
                 ]
                 ], (* end if 70, else if 60 *)
(* added also or flag2 == 4 *)
                If[((flag2 == 3) || (flag2 == 4)), 
                   (* start if 90, then part *)
                  If[debugstripadjust, 
                     Print["do because flag3 should be 3 or 4, check:",flag2]
                     ];
                  coefficient = coefficient*Part[expr3s,1]^Part[expr3s,2];
               If[debugstripadjust, 
              Print["At Pt. STR17a, for k=",k," and s=",s," then part if 80"]; 
              Print["coefficient=", coefficient]
                 ]
                  ]; (* end if 90 *)
                  expr2 = expr2/Part[expr3s,1]^Part[expr3s,2];
              If[debugstripadjust, 
              Print["At Pt. STR17b, for k=",k," and s=",s," then part if 80"];
              Print["expr2=", expr2]
                ]
               ], (* end if 60, else if 50 *)
        expr4 = Part[expr3s,1];
        If[debugstripadjust, 
           Print["At Pt. STR18, for k=",k," and s=",s," "<>
                 " in else part if 50, expr4=", expr4]
                 ];
        s = lenexpr3 + 1;
        If[debugstripadjust, 
           Print["At Pt. STR19, for k=",k," and s=",s," "<>
                 " in else part if 50, s=lenexpr3+1, now s=", s]
                 ]
        ] (* end if 50 *)
      ]; (* end for loop 200 loop iterator is s *)
   If[debugstripadjust, 
      Print["At Pt. STR20, after for loop 200, k=",k," "<>
            " Start switch part."]
     ];
(* format: switch[expr,form1,value1,form2,value2,...] *)
(* goal: compare expr to form1 if they match it returns value1, etc *)
   Switch[flag2, (* start switch 00, with flag2 *)
         1,
         If[debugstripadjust, 
            Print["At Pt. STR21, in switch function with k=",k," "<>
                  " flag2 value should be 1, check:", flag2]
           ];
         If[debugstripadjust, 
            Print["At Pt. STR21b, in switch function with k=",k," "<> 
                  " before appending list, list="];
            Print[list]
           ];
         list = Append[list,expr2];
         If[debugstripadjust, 
            Print["At Pt. STR22, in switch function with k=",k," "<>
                  " after appending list, list="];
            Print[list]
           ], (* next, option 2 *)
         2,
         If[debugstripadjust, 
            Print["At Pt. STR23, in switch function with k=",k," " <>
                  " flag2 value should be 2, check:",flag2]
         ];
         iexpr = Part[Head[expr4],1];
            If[debugstripadjust, 
               Print["At Pt. STR24, in switch function with k=",k," iexpr="];
               Print[iexpr]
               ];
(* Shifting happens always under flag2=2 *)
            expr2 = expr2 /. n->(2*n-iexpr);
            If[debugstripadjust, 
               Print["At Pt. STR25, in switch function with k=",k," expr2="];
               Print[expr2]
               ];
         list = Union[list,{expr2}];
         If[debugstripadjust, 
            Print["At Pt. STR26, in switch function with k=",k," list="];
            Print[list]
            ], (* next, option 3 *)
         3,
         If[debugstripadjust, 
            Print["At Pt. STR27, in switch function with k=",k," "<>
                  " flag2 value should be 3, check:", flag2]
         ];
(* Shifting happens always under flag2=3 *)
         iexpr = Part[Head[expr4],1];
         If[debugstripadjust, 
            Print["At Pt. STR28, in switch function with k=",k," iexpr="];
            Print[iexpr]
            ];
         expr2 = expr2 /. n->(2*n-iexpr);
           If[debugstripadjust, 
             Print["At Pt. STR29, in switch function with k=",k," expr2="];
             Print[expr2]
             ];
         list = Union[list,{{coefficient,expr2}}];
           If[debugstripadjust, 
             Print["At Pt. STR30, in switch function with k=",k," list="];
             Print[list]
            ], (* next, option 4 *)
         4, 
         If[debugstripadjust, 
            Print["At Pt. STR31, in switch function with k=",k," "<>
            " flag2 value should be 4, check:", flag2]
            ];
(* No shifting under flag2=4 *)
         expr2=expr2;
         If[debugstripadjust, 
            Print["At Pt. STR32, in switch function with k=",k," expr2="];
            Print[expr2]
            ];
         list = Union[list,{{coefficient,expr2}}];
         If[debugstripadjust, 
            Print["At Pt. STR33, in switch function with k=",k," list="];
            Print[list]
            ]
         ]; (* end switch 00 *)
If[debugstripadjust, 
   Print["At Pt. STR34, after switch is finished, k=",k," and iexpr="];
   Print[iexpr]
   ];
If[debugstripadjust, 
   Print["At Pt. STR35, after switch is finished, k=",k," and expr2="];
   Print[expr2]
   ];
If[debugstripadjust, 
   Print["At Pt. STR36, after switch is finished, k=",k," and list="];
   Print[list]
   ]
]; (* end for 1300, loop iterator is k *)
Return[list]
]; (* end module stripadjust *)
(* ######################## E11 ################################# *)

(* ######################## B12 ################################# *)
(*****************************************************************************)
(* setuniformrank[]: a subroutine for the scaling function                   *)
(*****************************************************************************)
setuniformrank[list_List] := Module[{i,lengthlist,syslist = {}},
lengthlist = Length[list];
For[i = 1,i < lengthlist,i++,
   syslist = Append[syslist,Part[list,i] == Part[list,lengthlist]];
   ];
Return[syslist]
]; 
(* ######################## E12 ################################# *)

(* ######################## B13 ################################# *)
(*****************************************************************************)
(* scaling[]: determines the scaling properties, if the equations are        *)
(*            incompatible it will print the appropriate messages            *)
(*****************************************************************************)
scaling[eqlist_List] := Module[
{pointsym,expr,list,syslist,msyslist0={},msyslist={},i,j,k,scalesol,
lenmsyslist,tempscalesol = {},tempmsyslist,trouble,troubleleft,
troubleright,posleft,posright,troublelist},
If[weightpars =!= {}, (* if 1 *)
  pointsym0 = Union[Table[weightpars[[i]] -> E^zeroweight[weightpars[[i]]],
                                             {i,1,Length[weightpars]}
                         ],{u[i_Integer][_][t] :> E^zeroweightu[i], t -> E^(0)}
                   ]; 
  If[debugweightsZ, 
     Print["At Pt. Z00, Scale0 w(d/dt)=0, pointsym0= "];
     Print[pointsym0]
     ], (* else if 1 *)
  pointsym0 = {u[i_Integer][_][t] :> E^zeroweightu[i], t -> E^(0)};
  If[debugweightsZ, 
     Print["At Pt. Z00, Scale0 w(d/dt)=0, pointsym0= "];
     Print[pointsym0]
     ];
  ]; (* end if 1 *)
For[i = 1, i <= noeqs,i++,
   expr[i] = Part[eqlist,i];
If[debugweightsZ, 
   Print["At Pt. Z01, Scale0 w(d/dt)=0, for i=",i,", ENTERING stripadjust"<>
         " expr[i]= "];
   Print[expr[i]]
   ];
(* First call of stripadjust in scaling function for Scale0 with w(d/dt)=0 *)
If[debugstripadjust, 
   Print["At Pt. Z02,"];
   Print["FIRST CALL OF STRIPADJUST --> IN SCALING FUNCTION, SCALE W(D/DT)=0"];
   Print["Usage: stripadjust[expr[i],True,1]"];
   Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
  ];
   list[i] = stripadjust[expr[i],True,1];
If[debugweightsZ, 
   Print["At Pt. Z03, Scale0 w(d/dt)=0, for i=",i,", LEAVING stripadjust"<>
         " list[i]= "];
   Print[list[i]]
   ];
   list[i] = list[i] /. pointsym0;
If[debugweightsZ, 
   Print["At Pt. Z04, Scale0 w(d/dt)=0, for i=",i,", list[i]= "];
   Print[list[i]]
   ];
   list[i] = PowerExpand[Log[list[i]]];
If[debugweightsZ, 
   Print["At Pt. Z05, Scale0 w(d/dt)=0, for i=",i,", list[i]= "];
   Print[list[i]]
   ];
(* Here we assign w(d/dt)=0 for Scale0 *)
   list[i] = Prepend[list[i],zeroweightu[i]+0];
If[debugweightsZ, 
   Print["At Pt. Z06, Scale0 w(d/dt)=0, for i=",i,", list[i]= "];
   Print[list[i]]
   ];
(* Added: remove the dubliplaces in list *)
   list[i] = Union[list[i]];
If[debugweightsZ, 
   Print["At Pt. Z07, Scale0 w(d/dt)=0, for i=",i,", list[i]= "];
   Print[list[i]];
   ];
   syslist[i] = setuniformrank[list[i]];
If[debugweightsZ, 
   Print["At Pt. Z08, Scale0 w(d/dt)=0, for i=",i,", syslist[i]= "];
   Print[syslist[i]]
   ];
   msyslist0 = Union[msyslist0,syslist[i]];
   ]; (* end for *)
Print[" "];
Print["LINEAR SYSTEM FOR THE WEIGHTS corresponding to Scale0 with w(d/dt)=0:"];
Print[msyslist0]; (* was blank *)
scalesol0 = Flatten[Solve[msyslist0]];
Print["SOLUTION OF THE SCALING EQUATIONS for Scale0 with w(d/dt)=0:"];
Print[scalesol0];
(* Print["START NEXT CASE!"]; *) (* was blank *)

(* Assign to Scale1 with w(d/dt)=1 *)
If[weightpars =!= {},
  pointsym = Union[Table[weightpars[[i]] -> E^weight[weightpars[[i]]],
                         {i,1,Length[weightpars]}
                        ], {u[i_Integer][_][t] :> E^weightu[i], t -> E^(-1)}
                  ];
  If[debugweightsZ, 
     Print["At Pt. Z300, Scale1 w(d/dt)=1, pointsym= "];
     Print[pointsym]
     ], (* else if 1 *)
  pointsym = {u[i_Integer][_][t] :> E^weightu[i], t -> E^(0)};
If[debugweightsZ, 
   Print["At Pt. Z301, Scale1 w(d/dt)=1, pointsym= "];
   Print[pointsym]
   ];
  ]; (* end if 1 *)
For[i = 1, i <= noeqs,i++, (* start for 78 *)
   expr[i] = Part[eqlist,i];
If[debugweightsZ, 
   Print["At Pt. Z302, Scale1 w(d/dt)=1, for i=",i,", ENTERING stripadjust"<>
         " expr[i]= "];
   Print[expr[i]]
   ];
(* Second call of stripadjust ---> in scaling function Scale1 with w(d/dt)=1 *)
If[debugstripadjust, 
   Print["At Pt. Z303."];
   Print["SECOND CALL OF STRIPADJUST --> IN SCALING FUNCTION WITH W(D/DT)=1."];
   Print["Usage: stripadjust[expr[i],True,1]"];
   Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
  ];
   list[i] = stripadjust[expr[i],True,1];
If[debugweightsZ, 
   Print["At Pt. Z304, Scale1 w(d/dt)=1, for i=",i,", LEAVING stripadjust"<>
         " list[i]= "];
   Print[list[i]]
   ];
   list[i] = list[i] /. pointsym;
If[debugweightsZ, 
   Print["At Pt. Z305, Scale1 w(d/dt)=1, for i=",i,", list[i]= "];
   Print[list[i]]
   ];
   list[i] = PowerExpand[Log[list[i]]];
If[debugweightsZ, 
   Print["At Pt. Z306, Scale1 w(d/dt)=1, for i=",i,", list[i]= "];
   Print[list[i]]
   ];
(* Here we assign w(d/dt)=1 for Scale1 *)
   list[i] = Prepend[list[i],weightu[i]+1];
If[debugweightsZ, 
   Print["At Pt. Z307, Scale1 w(d/dt)=1, for i=",i,", list[i]= "];
   Print[list[i]]
   ];
(* Added: remove the dubliplaces in list *)
   list[i] = Union[list[i]];
If[debugweightsZ, 
   Print["At Pt. Z308, Scale1 w(d/dt)=1, for i=",i,", list[i]= "];
   Print[list[i]];
   ];
   syslist[i] = setuniformrank[list[i]];
If[debugweightsZ, 
   Print["At Pt. Z309, Scale1 w(d/dt)=1, for i=",i,", syslist[i]= "];
   Print[syslist[i]]
   ];
(* original msyslist, not yet renamed into msyslist1 *)
   msyslist = Union[msyslist,syslist[i]];
   ]; (* end for 78 *)
Print["LINEAR SYSTEM FOR THE WEIGHTS corresponding to Scale1 with w(d/dt)=1:"];
Print[msyslist]; (* was blank *)
scalesol = Flatten[Solve[msyslist]];
Print["SOLUTION OF THE SCALING EQUATIONS for Scale1 with w(d/dt)=1: "];
Print[scalesol]; (* was blank *)
(* We continue building up the form of rho for Scale1 with w(d/dt)=1 *)

If[MemberQ[msyslist,False], (* start if 345 *)
  Print["Fatal Error! Weights are incorrect."];
  Print["Check your choice(s) for the weights in the data file."];
  Print["Aborting the computations!"];
  Closelog[], (* was ; *)
  (* Abort[], *)
  If[scalesol === {} && msyslist =!= {True}, (* start if 67 *)
    Print["In the given system there is at least one equation with terms"];
    Print["of unequal rank. Scaling properties can not be determined for"];
    Print["this system. The program will try to find the conflict, and,"];
    Print["if successful, provide suggestions to help resolve the conflict."];
    lenmsyslist = Length[msyslist];
    i = 1;
    While[tempscalesol === {} && i <= lenmsyslist,
         tempmsyslist = Complement[msyslist,{Part[msyslist,i]}];
         tempscalesol = Flatten[Solve[tempmsyslist]];
         i++;
         ];
    i--;
    If[tempscalesol =!= {},
      trouble = Part[msyslist,i];
      troubleleft = Part[trouble,1];
      troubleright = Part[trouble,2];
      For[j = 1, j <= noeqs, j++,
         list[j] = Drop[list[j],1];
         posleft = Flatten[Position[list[j],troubleleft,{1}]];
         posright = Flatten[Position[list[j],troubleright,{1}]];
         If[posleft =!= {} && posright =!= {},
           troublelist = Union[
             Table[Part[Expand[u[j][n]'[t]],posright[[k]]],
                     {k,1,Length[posright]}],
             Table[Part[Expand[u[j][n]'[t]],posleft[[k]]],
                     {k,1,Length[posleft]}]
                ];
           Print["The terms"];
           Print[troublelist];
           Print["in equation ",j," are incompatible. Try to introduce an"];
           Print["auxiliary parameter with weight as coefficient of some of"];
           Print["these terms. Aborting the computations!"];
           ];
         ],
    Print["Try to introduce auxiliary parameters with weight as coefficients"];
    Print["into the system. Aborting the computations!"];
    ];
    CloseLog[];
    Abort[],
    Return[scalesol]
    ]; (* end if 67 *)
  ] (* end if 345 *)
]; (* end module scaling *)
(* ######################## E13 ################################# *)

(* ######################## B14 ################################# *)
(*****************************************************************************)
(* subroutine5[]: a subroutine for constructformrho function                 *)
(*****************************************************************************)
(* subroutine5: nodims is first argument, varscalelist is second argument *)
(* maximumshift in third argument *)
(* orginally we used nodims = noeqs + Length[weightpars]; *)
subroutine5[nodims_,givenlist_,maximumshift_] := Module[
{m,boundary,scale,var,len,tempwei,temppar,list,pair,tempnodims,
ctm = 0,k,s,tempmaximumshift},
  tempmaximumshift=maximumshift;
  If[debugsubroutine5, 
     Print["Pt. CAM0300, before up-shifts, tempmaximumshift= "];
     Print[tempmaximumshift]
     ]; 
  If[debugsubroutine5, 
     Print["Pt. CAM0301, before up-shifts, noeqs= "];
     Print[noeqs]
     ]; 
  If[debugsubroutine5, 
     Print["Pt. CAM04, before up-shifts, nodims= "];
     Print[nodims]
     ]; 
If[forceshiftsinrho, 
  Print["USING SHIFTS in the construction of RHO!"];
  tempnodims = nodims + noeqs*tempmaximumshift, (* else *)
  tempnodims = nodims, (* else undecided *)
  tempnodims = nodims
  ]; 
  If[debugsubroutine5, 
     Print["Pt. CAM05, after up-shifts, tempnodims= "];
     Print[tempnodims]
     ]; 
list[0] = {};
For[ (* start for 66, loop for m *)
   m = 1, m <= tempnodims, m++,
   pair = Part[givenlist,m];
   If[debugroutine5N, 
      Print["At Pt. N1, for m=",m," pair= "];
      Print[pair]
      ];
   var = Part[pair,1];
   If[debugroutine5N, 
      Print["At Pt. N2 and CAM06, for m=",m," var= "];
      Print[var]
      ];
   scale = Part[pair,2];
   If[debugroutine5N, 
      Print["At Pt. N3 and CAM07, for m=",m," scale= "];
      Print[scale]
      ];
   If[scale =!= 0, (* start if 98, for scale *)
(* replaced m == 1 by m === 1 *)
     If[m === 1, (* begin if 47,  m=1 *)
       boundary = Floor[rhorank/scale];
       If[debugroutine5N, 
         Print["At Pt. N4, for m=1 only, boundary= "];
         Print[boundary]
         ];
(* WH 06/10/2003, replaced i by iii *)
       list[m] = Table[{iii*scale,var^iii},{iii,0,boundary}];
       If[debugroutine5N, 
         Print["At Pt. N5 and CAM08, for m=1 only, list[1]= "];
         Print[list[m]]
         ];
       ctm = m,
       If[debugroutine5N, 
         Print["At Pt. N6, for m=1 only, ctm= "];
         Print[ctm]
         ];
       len = Length[list[m-1]];
       If[debugroutine5N, 
         Print["At Pt. N7, for m=1 only, len= "];
         Print[len]
         ];
       list[m] = {};
       If[debugroutine5N, 
         Print["At Pt. N8, for m=1 only, list[1]= "];
         Print[list[m]]
         ];
       For[k = 1,k <= len,k++, (* start for 77, loop for k *)
          tempwei = Part[Part[list[m-1],k],1];
          If[debugroutine5N, 
            Print["At Pt. N9 and CAM0800, for k=",k," and m=1, tempwei= "];
            Print[tempwei]
            ];
          temppar = Part[Part[list[m-1],k],2];
          If[debugroutine5N, 
            Print["At Pt. N10 and CAM0801, for k=",k," and m=1, temppar= "];
            Print[temppar]
            ];
          boundary = Floor[(rhorank-tempwei)/scale];
          If[debugroutine5N, 
            Print["At Pt. N11, for k=",k," and m=1, boundary= "];
            Print[boundary]
            ];
          For[s = 0,s <= boundary,s++, (* start for 32  loop for s *)
             list[m] = Union[list[m],{{tempwei+s*scale,temppar*var^s}}];
             If[debugroutine5N, 
             Print["At Pt. N12 and CAM09, for k=",k," and s=",s,", list[m]= "];
               Print[list[m]]
               ];
             ctm = m;
             If[debugroutine5N, 
               Print["At Pt. N13, for k=",k," and s=",s,", updated ctm= "];
               Print[ctm]
               ];
             ]; (* end for 32, loop for s *)
          ]; (* end for 77, loop for k *)
       ]; (* end if 47, m=1 *)
     ]; (* end if 98, for scale *)
   ]; (* end for 66, loop for m *)
Return[list[ctm]]
]; (* end module subroutine5 *)
(* ######################## E15 ################################# *)

(* ######################## B16 ################################# *)
(*****************************************************************************)
(* originalconstructformrho[]: returns the form of rho                       *)
(*****************************************************************************)
originalconstructformrho[nodims_,varscalelist_,maximumshift_] := Module[
{list,len,tempwei,temppar,inttest,difflist,lendifflist,temp,
formrholist = {},i,j},
Print["Starting the construction of the form of the density"];
Print["with the single-scale (Scale1) routine based on w(d/dt)=1."];
   If[debugconstructformrhoOM, 
      Print["At Pt. OM1, first parameter of routine5 is nodims= "];
      Print[nodims]
      ];
   If[debugconstructformrhoOM, 
      Print["At Pt. OM2, 2nd parameter of routine5 is varscalelist= "];
      Print[varscalelist]
      ];
list = subroutine5[nodims,varscalelist,maximumshift];
   If[debugconstructformrhoOM, 
      Print["At Pt. OM3, LEAVING routine5, list= "];
      Print[list]
      ];
len = Length[list];
For[i = 1,i <= len,i++, (* for 1 *)
   tempwei = Part[Part[list,i],1];
   If[debugconstructformrhoOM, 
      Print["At Pt. OM4, for i=",i," tempwei= "];
      Print[tempwei]
      ];
   temppar = Part[Part[list,i],2];
   If[debugconstructformrhoOM, 
    Print["At Pt. OM5, for i=",i," VARIABLE before adjusted "<>
          " differentiation, temppar= "];
      Print[temppar]
      ];
   inttest = rhorank-tempwei;
   If[debugconstructformrhoOM, 
      Print["At Pt. OM5, for i=",i," inttest= "];
      Print[inttest];
      Print["If integer and has u then differentiate inttest times!"]
      ];
   If[IntegerQ[inttest] && Not[FreeQ[temppar,u]], (* if 3 *)
     difflist = D[temppar,{t,inttest}];
      If[debugconstructformrhoOM, 
         Print["At Pt. OM6, for i=",i," ENTERING stripadjust, difflist= "];
         Print[difflist]
         ];
(* Third call of stripadjust in old constructformrho (single scale) *)
If[debugstripadjust, 
   Print["At Pt. OM7."];
   Print["THIRD CALL OF STRIPADJUST --> IN CONSTRUCTFORMRHO (single scale)."];
   Print["Usage: stripadjust[difflist,True,2]"];
   Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
  ];
     difflist = stripadjust[difflist,True,2];
      If[debugconstructformrhoOM, 
         Print["At Pt. OM8, for i=",i," LEAVING stripadjust, difflist= "];
         Print[difflist]
         ];
     lendifflist = Length[difflist];
       If[debugconstructformrhoOM, 
         Print["At Pt. OM9, for i=",i," lendifflist= "];
         Print[lendifflist]
         ];
     For[j = 1,j <= lendifflist,j++, (* for 2 *)
        temp = Part[difflist,j];
        If[debugconstructformrhoOM, 
           Print["At Pt. OM10, for i=",i," and j=",j," temp= "];
           Print[temp]
          ];
        If[Intersection[formrholist,{temp}] === {}, (* if 5 *)
          formrholist = Union[formrholist,{temp}]; 
          ]; (* end if 5 *)
        ]; (* end for 2 *)
     ]; (* end if 3 *)
   ]; (* end for 1 *)
lenformrho = Length[formrholist];
formrholist = Union[formrholist]; 
formrholist = Complement[formrholist, {0}];  
lenformrho = Length[formrholist];
For[i = 1,i <= lenformrho,i++,
   formrho = formrho + c[i]*Part[formrholist,i];
   ];
Print[" "];
Print["For RANK = ",rhorank,", this is the form of the density rho: "];
Print[" "];
Print[subscriptform[formrho]];
Print[" "];
If[formrho === 0,
  Print["The only density is the trivial density rho=0."];
  Print[" "];
  Print["Aborting!"]; 
  CloseLog[]; Abort[] ] 
]; (* end module original constructformrho *)               
(* ######################## E16 ################################# *)

(* ######################## B17 ################################# *)
(*****************************************************************************)
(* newconstructformrho[]: returns the form of rho                            *)
(*****************************************************************************)
newconstructformrho[nodims_,varscalelist_,maximumshift_] := Module[
{list,len,tempwei,temppar,inttest,difflist,lendifflist,tempterm,term,
specificweightlist={},weightplustermlist={},i,j,mm,specificformrholist,
densitychoice,templist,weightnow,specificformrho,formrhochosen,ii},
(* Construction of pieces for formrho happens in routine5 *)
(* which builds all pieces compatible with the original scaling. *)
Print["Starting the construction of the form of the density"];
Print["with the multiple scale routine with both "];
Print["Scale1 with w(d/dt)=1, and Scale0 with w(d/dt)=0."];
   If[debugconstructformrhoNM, 
      Print["At Pt. NM1, first parameter of routine5 is nodims= "];
      Print[nodims]
      ];
   If[debugconstructformrhoNM, 
      Print["At Pt. NM2, 2nd parameter of routine5 is varscalelist= "];
      Print[varscalelist]
      ];
list = subroutine5[nodims,varscalelist,maximumshift];
   If[debugconstructformrhoNM, 
      Print["At Pt. NM3, LEAVING routine5 is list= "];
      Print[list]
      ];
len = Length[list];
For[i = 1,i <= len,i++, (* start for loop 1 *)
   tempwei = Part[Part[list,i],1];
   If[debugconstructformrhoNM, 
      Print["At Pt. NM4, for i=",i," tempwei= "];
      Print[tempwei]
      ];
   temppar = Part[Part[list,i],2];
   If[debugconstructformrhoNM, 
    Print["At Pt. NM5, for i=",i," VARIABLE before "<>
          " adjusted differentiation, temppar= "];
      Print[temppar]
      ];
   inttest = rhorank-tempwei;
   If[debugconstructformrhoNM, 
      Print["At Pt. NM5, for i=",i," inttest= "];
      Print[inttest]
      ];
   If[ (* start if 1 *)
     IntegerQ[inttest] && Not[FreeQ[temppar,u]], 
     difflist = D[temppar,{t,inttest}];
      If[debugconstructformrhoNM, 
         Print["At Pt. NM6, for i=",i," ENTERING stripadjust, difflist= "];
         Print[difflist]
         ];
(* Fourth call of stripadjust in constructformrho (new) *)
If[debugstripadjust, 
 Print["At Pt. NM7."];
 Print["FOURTH CALL OF STRIPADJUST --> IN CONSTRUCTFORMRHO (multiple scale)."];
 Print["Usage: stripadjust[difflist,True,2]"];
 Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
 ];
     difflist = stripadjust[difflist,True,2];
      If[debugconstructformrhoNM, 
         Print["At Pt. NM8, for i=",i," LEAVING stripadjust, difflist= "];
         Print[difflist]
         ];
     lendifflist = Length[difflist];
       If[debugconstructformrhoNM, 
         Print["At Pt. NM9, for i=",i," lendifflist= "];
         Print[lendifflist]
         ];
     For[j = 1,j <= lendifflist,j++, (* start for loop 2 *)
        term = Part[difflist,j];
        If[debugconstructformrhoNM, 
           Print["At Pt. NM10, for i=",i," and j=",j," term= "];
           Print[term]
          ];

(* Here we can apply a second scale to separate the terms of different *)
(* weights in different lists. *)
(* Apply separate routine that computes the weights based on given scalings *)
If[debugformrhoL, 
   Print["At Pt. L0, start separate forms of density based on two scales!"]
   ];
   If[debugweightsL, 
      Print["At Pt. L1, term= "];
      Print[term]
      ];
   tempterm = term /. pointsym0;
   If[debugweightsL,
      Print["At Pt. L2, after applying pointsym0= "];
      Print[tempterm]
     ];
   tempterm = PowerExpand[Log[tempterm]]; 
   If[debugweightsL, 
      Print["At Pt. L3, after powerexpand pointsym0= "];
      Print[tempterm]
     ];
(* Apply the alternate scaling rules, here scalesol0 *)
   tempterm = tempterm /. scalesol0;
   If[debugweightsL, 
      Print["At Pt. L4, after applying scalesol0, tempterm= "];
      Print[tempterm]
      ];
(* tempterm below is weight under scalesol0 *)
(* We now select the terms that have weight zero under scalesol0 *)
(* Building a list of all the occuring weights under scalesol0 *)
specificweightlist = Append[specificweightlist,tempterm];
specificweightlist = Union[specificweightlist];
If[debugweightsL, 
  Print["At Pt. L4.1, specificweightlist= "];
  Print[specificweightlist]
  ];
weightplustermlist = Append[weightplustermlist,{tempterm,term}];
weightplustermlist = Union[weightplustermlist];
If[debugweightsL, 
  Print["At Pt. L4.2, weightplustermlist= "];
  Print[weightplustermlist]
  ]
        ]; (* end for loop 2 *)
     ]; (* end if 1 *)
   ]; (* end for loop 1 *)

(* Build the various forms for rho for the various specific weights *)
(* that occur in the list specificweightlist *)
If[debugweightsL, 
  Print["At Pt. L8, START NEW DO LOOP for splitting the rhos!"]
  ];
Do[ (* start do loop 9, loops over kk *)
  templist={};
  specificformrholist[kk]={};
  specificformrho[kk]=0;
  weightnow = Part[specificweightlist,kk];
  If[debugweightsM, 
    Print["At Pt. M10.1, working with weightnow="];
    Print[weightnow]
    ];
  For[mm=1, mm <= Length[weightplustermlist], mm++, (* start for loop 11 *)
      If[Part[weightplustermlist[[mm]],1]===weightnow, (* if test 12, then *)
         specificformrholist[kk]=
            Append[specificformrholist[kk],Part[weightplustermlist[[mm]],2]];
         If[debugweightsM, 
            Print["At Pt. M10.2, for kk=",kk," specificformrholist[kk]= "];
            Print[specificformrholist[kk]]
            ]; 
         templist=
            Append[templist,{weightnow,Part[weightplustermlist[[mm]],2]}];
         If[debugweightsM, 
            Print["At Pt. M10.3, for kk=",kk," templist= "];
            Print[templist]
            ] 
        ] (* end if test 12 *)
     ]; (* end for loop 11 *)
(* now update the weightplustermlist: remove pieces already used *)
   weightplustermlist = Complement[weightplustermlist, templist];
   If[debugweightsM, 
     Print["At Pt. M10.4, for kk=",kk," updated weightplustermlist= "];
     Print[weightplustermlist]
     ];
(* now build the form for specificformrho[kk] *)
For[ii=1,ii <= Length[specificformrholist[kk]],ii++, (* start for loop 13 *)
    specificformrho[kk] = specificformrho[kk]+
                         c[ii]*Part[specificformrholist[kk],ii]
   ]; (* end for loop 13 *)
 If[debugweightsM, 
    Print["At Pt. M10.5, for kk=",kk," specificformrho[kk]= "];
    Print[specificformrho[kk]]
    ], {kk,1,Length[specificweightlist]}
]; (* end do loop 9, loops over kk *)

(* Added weight to the density. *)
(* Start of do loop to display the separate densities rho that were built *)
Print["There are ",Length[specificweightlist]," candidates for the density."];
Print["LIST OF ALL THE CANDIDATE DENSITIES:"];
Do[ (* start do loop 14, loops over jj *) 
  Print["For RANK=",rhorank,", CHOICE ",jj,", scale-invariant density rho: "];
  Print[subscriptform[specificformrho[jj]]];
  Print["based on Scale0 with w(d/dt)=0. This density has weight: ",
        specificweightlist[[jj]]], 
  {jj,1,Length[specificweightlist]}
  ]; (* end do loop 14, loops over jj *) 

(* If there is more than one rho, we ask the user which rho to continue with.*)
(* Otherwise we would have to loop the rest of the code over all candidates. *)
(* May not be easy if the variables are not all local! *)
If[ (* start if 0, if there is only one piece then *)
  Length[specificweightlist]===1, 
  If[debugweightsK, 
     Print["then part of if 0"]
     ];
  formrho = specificformrho[1];
  If[debugweightsK, 
     Print["At Pt. K8.1, program continues with DENSITY, formrho= "];
     Print[formrho]
     ];
  lenformrho = Length[specificformrholist[1]],  (* else for if 0 *)
  (* if there are more pieces *)
   If[debugweightsK, 
     Print["else part of if 0"]
     ];
 If[ (* start if 1 *)
   densitychoice === 0, 
   If[debugweightsK, 
     Print["then part of if 1"]
     ];
   Print["Select a number from 1 to ", Length[specificweightlist]];
   densitychoice = Input["ENTER YOUR CHOICE: "];
   formrho=specificformrho[densitychoice];
   Print[" "];
   Print[subscriptform[formrho]];
   formrhochosen = subscriptform[formrho];
   Print["Continuing with DENSITY chosen by the user: "];
   Print[formrhochosen];
   If[debugweightsK, 
      Print["At Pt. K8.2then1, program continues with DENSITY, formrho= "];
      Print[formrho]
      ];
   lenformrho = Length[specificformrholist[densitychoice]], 
   (* else part of if 1 *)
      If[debugweightsK, 
        Print["else part of if 1"]
        ];
    If[ (* start if 2 *)
      (densitychoice >= 1 && densitychoice <= Length[specificweightlist]),
      (* then part of if 2 *)
      If[debugweightsK, 
        Print["then part of if 2"]
        ];
      formrho = specificformrho[densitychoice];
      If[debugweightsK, 
        Print["At Pt. K8.3then2, program continues with DENSITY: "];
        Print[formrho]
        ];
      Print[" "];
      Print["Continuing with DENSITY (number selected in data file): "];
      Print[formrho];
      lenformrho = Length[specificformrholist[densitychoice]], 
      (* else part of if 2 *)
      If[debugweightsK, 
         Print["else part of if 2"]
         ];
      Print["Select a number from 1 to ", Length[specificweightlist]];
      densitychoice = Input["ENTER YOUR CHOICE: "];
      formrho = specificformrho[densitychoice];
      formrhochosen = subscriptform[formrho];
      Print[" "];
      Print["Continuing with DENSITY chosen by the user: "];
      Print[formrhochosen];
        If[debugweightsK, 
           Print["At Pt. K8.4else2, program continues with DENSITY: "];
           Print[formrho]
           ];
      lenformrho = Length[specificformrholist[densitychoice]], 
      (* undecided part of if 2 *)
      If[debugweightsK, 
        Print["undecided part of if 2"]
        ];
      Print["Select a number from 1 to ", Length[specificweightlist]];
      densitychoice = Input["ENTER YOUR CHOICE: "];
      formrho = specificformrho[densitychoice];
      formrhochosen = subscriptform[formrho];
      Print[" "];
      Print["Continuing the DENSITY chosen by the user: "];
      Print[formrhochosen];
        If[debugweightsK, 
           Print["At Pt. K8.4und2, program continues with DENSITY: "];
           Print[formrho]
           ];
      lenformrho = Length[specificformrholist[densitychoice]]
      ], (* end if 2 *)
   (* undediced part of if 1 *)
      If[debugweightsK, 
        Print["undecided part of if 1"]
        ];
   Print["Select a number from 1 to ",Length[specificweightlist]];
   densitychoice = Input["ENTER YOUR CHOICE: "];
   formrho = specificformrho[densitychoice]; 
   formrhochosen = subscriptform[formrho];
   Print[" "];
   Print["Continuing with the DENSITY chosen by the user: "];
   Print[formrhochosen];
     If[debugweightsK, 
        Print["At Pt. K8.4und1, program continues with DENSITY: "];
        Print[formrho]
        ];
   lenformrho = Length[specificformrholist[densitychoice]]
   ], (* end if 1 *) 
   If[debugweightsK, 
     Print["undecided part of if 0"]
     ]
  ]; (* end of if 0 *)
If[debugweightsK, 
   Print["At Pt. K9, program continues with DENSITY: "];
   Print[formrho]
  ];
If[formrho === 0,
  Print["The only density is the trivial density rho=0."];
  Print[" "];
  Print["Aborting!"]; 
  CloseLog[]; Abort[] 
  ] 
]; (* end module new constructformrho *)               
(* ######################## E17 ################################# *)

(* ######################## B18 ################################# *)
(*****************************************************************************)
(* subroutine4[]: a subroutine for evaluate function                         *)
(*****************************************************************************)
subroutine4[givenexpr_] := Module[
{expr1,expr2,freerule={},part1,lenexpr2,i},
expr1 = Expand[givenexpr];
If[debugsubroutine4, 
  Print["At Pt. SR1, expanded version of entering expr, expr1=",expr1]
  ];
expr2 = Numerator[Factor[expr1]];
If[debugsubroutine4, 
  Print["At Pt. SR2, numerator of expr1, expr2= ",expr2]
  ];
expr2 = FactorList[expr2];
If[debugsubroutine4, 
  Print["At Pt. SR3, factorlist on expr2, expr1= ",expr2]
  ];
lenexpr2 = Length[expr2];
If[debugsubroutine4, 
  Print["At Pt. SR4, length of expr2, lenexpr2= ",lenexpr2]
  ];
For[i = 1,i<= lenexpr2,i++, (* start for 999, loop over i *)
   part1 = Part[Part[expr2,i],1];
If[debugsubroutine4, 
  Print["At Pt. SR5, in for loop, part1= ",part1]
  ];
   If[MemberQ[unknownlist,part1],
     freerule = Union[freerule,{part1 -> 1}];
     If[debugsubroutine4, 
        Print["At Pt. SR6, in for loop, freerule= ",freerule]
        ]
     ];
   ]; (* end for loop 999 *)
If[debugsubroutine4, 
  Print["At Pt. SR7, after for loop, freerule, freerule= ",freerule]
  ];
If[freerule =!= {} && printflagcommonfactor,
  Print[" "];
  Print["Caution! There is a common (free) factor in the density."];
  Print[" "];
  Print["Setting ",freerule,"."];
  ];
If[debugsubroutine4, 
  Print["At Pt. SR8, returning freerule= ",freerule]
  ];
Return[freerule] 
]; (* end module subroutine4 *)
(* ######################## E18 ################################# *)

(* ######################## B19 ################################# *)

(* WH 10/04/2003 *)
(*****************************************************************************)
(* evaluate[]: computes the density and flux with verification               *)
(* Last modified: July 1, 2003 by Mike C. at 06:13.                          *)
(* Added call to Discrete Homotopy Operator, removed old way of              *)
(* calculating J_n                                                           *)
(*****************************************************************************)
evaluate[expr1_,expr3_,expr2_,solrules_] := Module[
{newformrho,formjn = 0,test,factortest,flux,evformjn,jnplus1,
diffevjnjnplus,freerule,rho,dtrho},
newformrho = expr1 /. solrules;

(* Mike C., July 1, 2003 *)
(* Don't calculate Jn this way, use the Discrete Homotopy Operator below *)

(* WH 10/04/2003 *)
(* BEGIN TAKEN OUT BY Mike C.
(* BEGIN computation of final form of jn *)
If[debugevaluate, 
   Print["At Pt. HOL100, solrules= "];
   Print[solrules]
   ];
evformjn = Expand[expr3 /. solrules];
If[debugevaluate, 
   Print["At Pt. HOL200, evaluated jn, ENTERING subroutine4, evformjn= "];
   Print[evformjn]
   ];
(* First call of subroutine4 *)
(* added the else part below *)
If[evformjn =!= 0,
  freerule = subroutine4[evformjn];
  flux = evformjn /. freerule;
  flux = Expand[flux], (* else *)
  flux = evformjn
  ];
If[debugevaluate, 
   Print["At Pt. HOL300, evaluated jn, LEAVING subroutine4, flux= "];
   Print[flux]
   ];
  flux = flux/. backwardtranslation;
jnplus1 = Part[Map[upShift,{flux}],1];
If[debugevaluate, 
   Print["At Pt. HOL400, evaluated jn+1, jnplus1= "];
   Print[jnplus1]
   ];
diffevjnjnplus = flux-jnplus1;
If[debugevaluate, 
   Print["At Pt. HOL500, evaluated jn-jn+1, diffevjnjnplus= "];
   Print[diffevjnjnplus]
   ];
(* END computation of final form of jn *)
END TAKEN OUT by Mike C. *)

(* Mike C., July 1, 2003 *)
(* don't need any thing between BEGIN and END above *)

(* We calculate J_n after rho is computed below *).
(* WH 10/04/2003 *)

If[debugevaluate, 
   Print["At Pt. HOL501, ENTERING subroutine4, newformrho= "]; 
   Print[newformrho] 
   ]; 
(* Second call of subroutine4 *)
(* added the else part below *)
(* change the logic so that the test always happens *)
(* move matching bracked way higher up *)
If[newformrho =!= 0, (* start if 876 *)
  freerule = subroutine4[newformrho];
  rho = newformrho /. freerule;
  rho = Expand[rho];
  If[debugevaluate, 
     Print["At Pt. HOL600, LEAVING subroutine4, rho= "];
     Print[rho]
     ], (* else of if 876 *)
  rho = newformrho
  ]; (* end if 876 *)

(* Mike C., July 1, 2003 *)
(* Discrete Homotopy Operator integration checked by Mike C., July 1, 2003 *)
(* The variable "noeqs" is defined in each data file, and it specifies *)
(* the number of dependent variables of the form u[i][n+j][t] *)
(* WH 10/04/2003 HERE HERE HERE *)

(* WH 10/04/2003 Application of the homotopy operator happens below! *)

dtrho = D[rho, t] /. forwardtranslation;

(* WH 10/04/2003 *)
Print[" "];
Print["Starting the computation of J_n with the discrete homotopy operator!"]; 

flux = (-1)*discreteHomotopyOperator[dtrho, u, noeqs];

Print["Finished the computation of J_n with the discrete homotopy operator!"]; 

(* End of discreteHomotopyOperator integration in this function *)

(* Start of the final and independent test for density and flux *)
  flux = flux /. backwardtranslation;
  test = (D[rho,t] - flux + (flux /. n -> n+1)) /. solrules;
  factortest = Expand[test];
If[debugevaluate, 
   Print["At Pt. HOL700, factortest= "];
   Print[factortest]
   ];
  If[factortest =!= 0, factortest = Factor[factortest]];
  If[factortest =!= 0, (* start if 24 *)
    Print[" "];
    Print["Automatic verification of the result FAILED! 
     There is NO DENSITY of this form!"], (* else *)
If[debugevaluate,
(* WH 10/04/2003 *)
    Print["rho_n and flux J_n PASSED the test!"];
    Print[" "]
    ]
   ]; (* end if 24 *)
  Print[" "];
  Print["Result of explicit verification: D_t (rho) - J_n + J_{n+1} = ",
         factortest];
  If[factortest =!= 0, (* start if 25 *)
    Print[" "];
    Print["Automatic verification of the result FAILED!"<> 
          " There is NO DENSITY of this form!"], (* else *)
    Print["Density rho and flux PASSED the test!"];
    Print[" "]
    ];  (* end if 25 *)
  Print["***************************************************************"];
  Print[" "];
  Print["This is the density, rho: "];
  Print[" "];
  Print[subscriptform[rho]];
  Print[" "];
(* 02/04/2003 capturing the rho in list *)
  If[rho =!= 0, listallrhos = Append[listallrhos,rho]];
  If[debugevaluate,
    Print["At Pt. EVAL1, listallrhos = "];
    Print[listallrhos]
    ];
  splitrho[subscriptform[rho]];
  Print["***************************************************************"];
  Print[" "];
  Print["The corresponding flux J_n is: "];
  Print[" "];
  Print[subscriptform[flux]];
  Print[" "];
  splitflux[subscriptform[flux]]; 
(* 02/04/2003 capturing the fluxes in list *)
  If[flux =!= 0, listallfluxes = Append[listallfluxes,flux]];
  If[debugevaluate,
    Print["At Pt. EVAL2, listallfluxes = "];
    Print[listallfluxes]
    ];
(* closing square bracket was here (* end if 876 *) *)
];  (* end module evaluate *)
(* ######################## E19 ################################# *)

(* ######################## B20 ################################# *)
(* Flags were added to constructeqlist so that it does not call the *)
(* stripadjust function with the same flags. *)
(* These flag setting play role for stripadjust function and switch *)
(*****************************************************************************)
(* constructeqlist[]: forms the system of equations for the c[i]             *)
(*****************************************************************************)
(* Added flags: flag1 and flag2. *)
(* Set Flag1=False and Flag2=3 for the shifting routine. *)
(* Set Flag1=False and Flag2=4 for Discrete Euler routine. *)
(* local vars remainder, mytempterm, etc *)
constructeqlist[expr_,flag1_,flag2_] := Module[
{expr1,expr2,templist1,templist2,lentemplist2,s,eqlist={},selectsame,
coefficient,subrule,eqlistshifting},
(* WH 06/10/2003 Line move just before the call on constructeqlist *)
(* 
Print["Starting the derivation of the equations of the system for c[i]."]; 
*)
expr1=Expand[expr];
If[debugconstructeqlist, 
   Print["At Pt. DD0, ENTERING (already in) constructeqlist, flag1= "];
   Print[flag1]
   ];
(* HERE HERE HERE *)
If[debugconstructeqlist, 
   Print["At Pt. DD1, ENTERING (already in) constructeqlist, flag2= "];
   Print[flag2]
   ];
If[debugconstructeqlist, 
   Print["At Pt. DD2, ENTERING (already in) constructeqlist, expr1= "];
   Print[expr1]
   ];
coefficient = expr1 //. u[i_][n_][t] :> 0;
(* simplify the key expressions as we go along *)
If[(forceextrasimplifications && coefficient =!= 0), (* start if 1200 *)
   If[debugconstructeqlist, 
     Print["At Pt. DD2NEW1, before mysimplify1, coefficient= "];
     Print[coefficient]
     ];
   coefficient = mysimplify1[{coefficient}];
   If[debugconstructeqlist,
     Print["At Pt. DD2NEW1bis, after mysimplify1, coefficient= "];
     Print[coefficient]
     ];
   coefficient = Part[coefficient,1];
   If[debugconstructeqlist,
     Print["At Pt. DD2NEW1tris, after taking part 1, coefficient= "];
     Print[coefficient]
     ];
   If[forcestripparameters, 
      (* Print["parameters = ", parameters]; *)
      If[debugconstructeqlist, 
        Print["At Pt. DD2NEW2, before stripparameters, coefficient= "];
        Print[coefficient]
        ];
      coefficient = stripparameters[{coefficient},parameters];
      If[debugconstructeqlist, 
        Print["At Pt. DD2NEW2bis, after stripparameters, coefficient= "];
        Print[coefficient]
        ];
      coefficient = Part[coefficient,1];
      If[debugconstructeqlist, 
        Print["At Pt. DD2NEW2tris, after taking part 1, coefficient= "];
        Print[coefficient]
        ]
   ];
   If[debugconstructeqlist, 
      Print["At Pt. DD2NEW3, leaving stripparameters, coefficient= "];
      Print[coefficient]
     ];
If[(Length[coefficient] === 1 && Not[FreeQ[coefficient,c]]), (* if 1300 *)
   subrule = {coefficient -> 0};   
   If[debugconstructeqlist, 
      Print["At Pt. DD2NEW4, based on stripped coefficient, subrule= "];
      Print[subrule]
     ];
   expr1 = expr1 /. subrule;
   If[debugconstructeqlist, 
      Print["At Pt. DD2NEW5, after applying subrule, expr1= "];
      Print[expr1]
     ]
 ] (* end if 1300 *)
];(* end if 1200 *)
If[debugconstructeqlist, 
   Print["At Pt. DD3, after setting u[i][n][t], u[i][n+p][t], u[i][n-p][t]"<> 
         " all zero, `constant' coefficient, coefficient= "];
   Print[coefficient]
   ];
eqlist = Union[eqlist,{coefficient}];
If[debugconstructeqlist, 
   Print["At Pt. DD4, after union of coefficient with eqlist, eqlist= "];
   Print[eqlist]
   ];
expr2 = expr1-coefficient;
If[debugconstructeqlist, 
   Print["At Pt. DD5, after removing expr1, ENTERING stripadjust,"<>
         " expr2= "];
   Print[expr2]
   ];
If[expr2 =!= 0, (* start if 300 *)
(* Fifth call of stripadjust in constructeqlist *)
If[debugstripadjust, 
   Print["At Pt. DD6."];
   Print["FIFTH CALL OF STRIPADJUST ---> IN CONSTRUCTEQLIST."];
(* flag1 and flag2 are passed on from the function constructeqlist *)
   Print["Usage: stripadjust[expr2,flag1,flag2]."];
   Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
  ];
(* for shifting routine: we use flag1=False and flag2=3, so we used *)
(* templist1 = stripadjust[expr2,False,3]; *)
(* for discrete euler routine: we use flag1=False and flag2=4 *)
   templist1 = stripadjust[expr2,flag1,flag2]; 
If[debugconstructeqlist, 
   Print["At Pt. DD7, LEAVING stripadjust, templist1= "];
   Print[templist1]
   ];
(* Use the shifting routine to produce the linear system for c[i] *)
If[forceshifting, (* start if 1000, then part *)
  eqlistshifting = eqlist;
  templist2 = Table[Part[Part[templist1,j],2],{j,1,Length[templist1]}];
If[debugconstructeqlist, 
   Print["WORKING WITH ORIGINAL SHIFTING ROUTINE TO BUILD THE LINEAR SYSTEM!"];
   Print["At Pt. DD8, working with type, templist2= "];
   Print[templist2]
   ];
  templist2 = Union[templist2];
If[debugconstructeqlist, 
   Print["At Pt. DD9, after union, list of all the type terms, templist2= "];
   Print[templist2]
   ];
  templist2 = Reverse[templist2];
  lentemplist2 = Length[templist2];
If[debugconstructeqlist, 
   Print["At Pt. DD10, length of templist2, lentemplist2= "];
   Print[lentemplist2]
   ];
  For[s=1, s<=lentemplist2, s++, (* start for 20 *)
     selectsame = Select[templist1,SameQ[#[[2]],Part[templist2,s]]&];
If[debugconstructeqlist, 
   Print["At Pt. DD11, selectsame= "];
   Print[selectsame]
   ];
   coefficient = Sum[Part[Part[selectsame,j],1],{j,1,Length[selectsame]}];
If[debugconstructeqlist, 
   Print["At Pt. DD12, coefficient of specific type term, coefficient= "];
   Print[coefficient]
   ];
If[(forceextrasimplifications && coefficient =!= 0), (* start if 1400 *)
   (* simplify the expression as we go along *)
   If[debugconstructeqlist, 
     Print["At Pt. DD12NEW1, before mysimplify1, coefficient= "];
     Print[coefficient]
     ];
    coefficient = mysimplify1[{coefficient}];
    If[debugconstructeqlist, 
      Print["At Pt. DD2NEW1tris, after mysimplify1, coefficient= "];
      Print[coefficient]
      ];
    coefficient = Part[coefficient,1];
    If[debugconstructeqlist, 
      Print["At Pt. DD2NEW1quatro, after taking part 1, coefficient= "];
      Print[coefficient]
      ];
  If[forcestripparameters, (* start if 456 *)
      (* Print["parameters = ", parameters]; *)
      If[debugconstructeqlist, 
        Print["At Pt. DD12NEW2, before stripparameters, coefficient= "];
        Print[coefficient]
        ];
      coefficient = stripparameters[{coefficient},parameters];
      If[debugconstructeqlist, 
        Print["At Pt. DD12NEW2bis, after stripparameters, coefficient= "];
        Print[coefficient]
        ];
      coefficient = Part[coefficient,1];
      If[debugconstructeqlist, 
        Print["At Pt. DD12NEW2tris, after taking part 1, coefficient= "];
        Print[coefficient]
        ]
    ]; (* end if 456 *)
   If[debugconstructeqlist, 
      Print["At Pt. DD12NEW3, leaving stripparameters, coefficient= "];
      Print[coefficient]
     ];
If[(Length[coefficient] === 1 && Not[FreeQ[coefficient,c]]), (* if 1500 *)
   subrule = {coefficient -> 0};   
   If[debugconstructeqlist, 
      Print["At Pt. DD12NEW4, subrule= "];
      Print[subrule]
     ];
   templist1 = templist1 /. subrule;
   If[debugconstructeqlist, 
      Print["At Pt. DD12NEW5, after applying subrule, templist1= "];
      Print[templist1]
     ]
] (* end if 1500 *)
]; (* end if 1400 *)
  eqlistshifting = Union[eqlistshifting,{coefficient}];
If[debugconstructeqlist, 
   Print["At Pt. DD13, union with previous eqlist, eqlistshifting= "];
   Print[eqlistshifting]
   ]
     ]; (* end for 20 *)
eqlist = eqlistshifting;
If[debugconstructeqlist, 
   Print["At Pt. DD14, eqlistshifting is named eqlist, eqlist= "];
   Print[eqlist]
   ]
]; (* end if 1000 *)
(* Use the Discrete Euler routine to compute the linear system for c[i]. *)
If[forcediscreteeuler, (* start if 2000 then part *)
  eqlisteuler = eqlist; 
  templist2 = Table[Part[Part[templist1,j],2],{j,1,Length[templist1]}];
If[debugconstructeqlist, 
   Print["WORKING WITH DISCRETE EULER ROUTINE TO BUILD THE LINEAR SYSTEM!"];
   Print["At Pt. DD8bis, working with type, templist2= "];
   Print[templist2]
   ];
  templist2 = Union[templist2];
If[debugconstructeqlist, 
   Print["At Pt. DD9bis, after union,list of all the type terms, templist2= "];
   Print[templist2]
   ];
  lentemplist2 = Length[templist2];
If[debugconstructeqlist, 
   Print["At Pt. DD10bis, length of templist2, lentemplist2= "];
   Print[lentemplist2]
   ];
  For[s=1, s<=lentemplist2, s++, (* start for 20 *)
     selectsame = Select[templist1,SameQ[#[[2]],Part[templist2,s]]&];
If[debugconstructeqlist, 
   Print["At Pt. DD11bis, selectsame= "];
   Print[selectsame]
   ];
   coefficient = Sum[Part[Part[selectsame,j],1],{j,1,Length[selectsame]}];
If[debugconstructeqlist, 
   Print["At Pt. DD12bis, coefficient of specific type term, coefficient= "];
   Print[coefficient]
   ];
   eqlisteuler = Union[eqlisteuler,{coefficient}];
If[debugconstructeqlist, 
   Print["At Pt. DD13bis, union with previous eqlist, eqlisteuler= "];
   Print[eqlisteuler]
   ]
     ]; (* end if 20bis *)
eqlist = eqlisteuler;
If[debugconstructeqlist, 
   Print["At Pt. DD14bis, eqlisteuler is named eqlist, eqlist= "];
   Print[eqlist]
   ]
  ] (* end if 2000 *)
 ]; (* end if 300 *)
Return[eqlist] 
]; (* end module constructeqlist *)
(* ######################## E20 ################################# *)

(* ######################## B21 ################################# *)
(*****************************************************************************)
(* mysimplify1[]: simplifies the system of equations for the c[i]            *)
(*****************************************************************************)
mysimplify1[list_List] := Module[
{newlist,lennewlist,k,fterm,simplelist={}},
newlist = Complement[list,{0}];
If[debugmysimplify1, 
  Print["At Pt. FSIM1, in mysimplify1, after expanding, newlist = "];
  Print[newlist]
  ];
lennewlist = Length[newlist];
For[k = 1 ,k <= lennewlist,k++, (* start for 39 *)
   fterm = Factor[newlist[[k]]];
   If[debugmysimplify1, 
      Print["At Pt. FSIM2, in mysimplify1, fterm = "];
      Print[fterm]
      ];
   If[Part[fterm,0] === Times && NumberQ[Part[fterm,1]],
     fterm = fterm/Part[fterm,1] ];
   If[debugmysimplify1, 
      Print["At Pt. FSIM3, in mysimplify1, after division by factor, fterm= "];
      Print[fterm]
      ];
   If[Intersection[Expand[simplelist],Expand[{fterm*(-1)}]] === {},
     simplelist = Union[simplelist,{fterm}]]
   ]; (* end for 39 *)
(* simplelist = Sort[simplelist, Length[Part[#1,1]] < Length[Part[#2,1]]&]; *)
If[debugmysimplify1, 
  Print["At Pt. FSIM4, in mysimplify1, after sorting, simplelist = "];
  Print[simplelist]
  ];
Return[simplelist]
]; (* end module mysimplify1 *)
(* ######################## E21 ################################# *)

(* ######################## B22 ################################# *)
(*****************************************************************************)
(* mysimplify2[]: applies equations of type c[i] == 0 to entire system       *)
(*****************************************************************************)
mysimplify2[list_List] := Module[
{newlist,lennewlist,sublist1={},sublist2,partk,sublist1rules},
newlist = MapAll[Expand,list];
If[debugmysimplify2, 
  Print["At Pt. SIM1, in mysimplify2, after expanding, newlist = "];
  Print[newlist]
  ];
lennewlist = Length[newlist];
(* pick out single term equations involving c[i] put in sublist1 *)
Do[
   partk = Part[newlist,k];
   If[debugmysimplify2, 
     Print["At Pt. SIM2, in mysimplify2, partk = "];
     Print[partk]
     ];
   If[(Length[partk] === 1 && Not[FreeQ[partk,c]]), 
      sublist1 = Append[sublist1,partk];
      If[debugmysimplify2, 
        Print["At Pt. SIM3, in mysimplify2, sublist1 = "];
        Print[sublist1]
        ]
      ], {k,1,lennewlist}
  ];
sublist1 = Union[sublist1];
sublist2 = Complement[newlist,sublist1];
If[debugmysimplify2, 
  Print["At Pt. SIM4, in mysimplify2, sublist2 = "];
  Print[sublist2]
  ];
(* make rules from sublist1 *)
If[sublist1 =!= {}, 
   sublist1rules = Table[Part[sublist1,k] -> 0, {k,1,Length[sublist1]}];   
   If[debugmysimplify2, 
     Print["At Pt. SIM5, in mysimplify2, sublist1rules = "];
     Print[sublist1rules]
     ];
   sublist2 = sublist2 /. sublist1rules;
   sublist2 = Union[sublist2];
   sublist2 = Complement[sublist2,{0}];
   If[debugmysimplify2, 
     Print["At Pt. SIM6, in mysimplify2, after sublist1rules, sublist2 = "];
     Print[sublist2]
     ]
  ];
newlist = Flatten[Join[sublist1,sublist2]];
If[debugmysimplify2, 
   Print["At Pt. SIM7, in mysimplify2, after joining both lists, newlist = "];
   Print[newlist]
   ];
(* newlist = Sort[newlist, Length[Part[#1,1]] < Length[Part[#2,1]]&]; *)
If[debugmysimplify2, 
   Print["At Pt. SIM9, in mysimplify2, after sorting, newlist = "];
   Print[newlist]
   ];
Return[newlist]
]; (* end module mysimplify2 *)
(* ######################## E22 ################################# *)

(* ######################## B23 ################################# *)
(*****************************************************************************)
(* analyzer[]: determines the coefficients c[i] that must be zero, used in   *)
(*             the search for compatibility conditions                       *)
(*****************************************************************************)
analyzer[list_,len_] := Module[
{i,partlist,analyzelist = {}},
Print["Starting the analysis of the system."];
For[i = 1,i <= len,i++, (* start for 40 *)
   partlist = Part[list,i];
If[debuganalyzer, 
   Print["At Pt. AN1, partlist = "];
   Print[partlist]
   ];
   If[Length[partlist] === 1 && MemberQ[unknownlist,partlist], 
     (* start if 111 *)
     analyzelist = Union[analyzelist,{partlist}];
   If[debuganalyzer, 
     Print["At Pt. AN2, analyzelist = "];
     Print[analyzelist]
     ],
     If[Length[partlist] === 2 && MemberQ[unknownlist,partlist[[2]]] &&
      (MemberQ[parameters,partlist[[1]]] || MemberQ[weightpars,partlist[[1]]]),
       (* end if 110 *)
       analyzelist = Union[analyzelist,{partlist[[2]]}];
       If[debuganalyzer, 
         Print["At Pt. AN3, analyzelist = "];
         Print[analyzelist]
         ] 
     ] (* end if 110 *)
   ]; (* end if 111 *)
   ]; (* end for 40 *)
Return[analyzelist]
]; (* end module analyzer *)
(* ######################## E23 ################################# *)

(* ######################## B24 ################################# *)
(*****************************************************************************)
(* picklhs[]: makes a system from the list by setting LHS == 0               *)
(*****************************************************************************)
picklhs[list_,k_] := Part[list,k] == 0;
(* ######################## E24 ################################# *)

(* ######################## B25 ################################# *)
(*****************************************************************************)
(* coefmat[]: creates the coefficient matrix of the system for the c[i]      *)
(*****************************************************************************)
coefmat[system_List,unknowns_List] := Module[
{i,mat = {},lensys},
If[system === {} || unknowns === {}, Print["Fatal Error!"],
  lensys = Length[system];
   For[i = 1,i <= lensys,i++,
       mat = Append[mat,Table[Coefficient[Expand[Part[Part[system,i],1]],
       Part[unknowns,k]],{k,1,Length[unknowns]}]]
       ]
   ];
Return[mat] 
]; (* end module coefmat *)
(* ######################## E25 ################################# *)

(* ######################## B26 ################################# *)
(*****************************************************************************)
(* main[] :                                                                  *)
(*****************************************************************************)
main[] := Module[
{rhot,myeqlist,seqlist,lenseqlist,maineqlist,analyzelist,lengthpar,
complexanalysis,syscondlist = {},i,j,k,syscond,inputlist,
inputpart,inputval,inputrule,comcond,comcondfac,sol,lengthsol,rules,
parameter,newmaineqlist,solc,solrules,comcondfactab,myfunclist,myrhot,
myvars,myresult,myresultnew,lenmyresultnew,tempeqlist,myrhotlist,jnlist,
obstruction,term,varsinterm,depvars,newdepvars,seplist,subsinterm,
minsubsinterm,upterms,downterms,obsterm,jntermlist,formjn,formjnplus,
difference,cleanseqlist,lencleanseqlist,coefmatrix},
rhot = D[formrho,t];
rhot = Expand[rhot];
Print["Starting the construction of the system for the c[i]."];
(* Application of the Discrete Euler operator if forcediscreteeuler=True *)
If[forcediscreteeuler, (* start if forcediscreteeuler *)
 Print["APPLICATION OF DISCRETE EULER ROUTINE TO COMPUTE THE LINEAR SYSTEM!"];
  If[debugdiscreteeuler, 
     Print["At Pt. LLL-6bis, before the forward translation, rhot=: "];
     Print[rhot]
     ];
 myfunclist = Table[u[i][n][t],{i,1,noeqs}]; 
 (* function discrete variational derivative requires format u[i][n,t] *)
 myfunclist = myfunclist /. forwardtranslation;
 If[debugdiscreteeuler,
    Print["At Pt. LLL-5bis, myfunclist going into DiscreteEulerD is: "];
    Print[myfunclist]
   ];
 myrhot = rhot /. forwardtranslation;
 If[debugdiscreteeuler, 
    Print["At Pt. LLL-4bis, after translation, going into DiscreteEulerD,"<>
          " myrhot= "];
    Print[myrhot]
    ];
 myvars = {n,t};
 myresult = DiscreteEulerD[myrhot,myfunclist,myvars];
 If[debugdiscreteeuler, 
    Print["At Pt. LLL-3bis, myresult coming out of DiscreteEulerD is: "];
    Print[myresult]
   ];
 myresultnew = myresult /. backwardtranslation;
 If[debugdiscreteeuler, 
    Print["At Pt. LLL-2bis, translated back to old notation, myresultnew: "];
    Print[myresultnew]
   ];
 (* Pass myresultnew to the constructeqlist and also stripadjust. *)
 (* Format: double list of expressions corresponding to each component. *)
 (* So, we have to loop of the various obstructions and combine *)
 (* the resulting equations in one system. *)
 If[debugdiscreteeuler, 
    Print["At Pt. LLL-1bis, start of Do loop "];
    Print["(run over each component of variational derivative)"]
    ];
 myresultnew = Flatten[myresultnew];
 If[debugdiscreteeuler, 
    Print["At Pt. LLL0bis, myresultnew= "];
    Print[myresultnew]
   ];
 lenmyresultnew = Length[myresultnew];
 eqlisteuler = {};
 Do[ (* start do 111 *)
    If[debugdiscreteeuler, 
     Print["At Pt. LLL1bis, ENTERING constructeqlist, Part[myresultnew,kk]= "];
     Print[Part[myresultnew,kk]]
      ];
(* For Discrete Euler routine set flag1=False and flag2=4 *)
(* First call of constructeqlist in discrete Euler routine *)
If[debugdiscreteeuler, 
   Print["At Pt. LLL2bis"];
   Print["First call of constructeqlist ---> in Discrete Euler routine"];
   Print["Usage: constructeqlist[Part[myresultnew,kk],False,4]"];
   Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
   ];
(* WH 06/10/2003 constructeqlist is used here *)
(*    
Print["Starting the derivation of the equations of the system for c[i]."]; 
*)
    tempeqlist[kk] = constructeqlist[Part[myresultnew,kk],False,4]; 
    If[debugdiscreteeuler, 
       Print["At Pt. LLL3bis, for kk=",kk," LEAVING constructeq, "<>
             " tempeqlist[kk]= "];
       Print[tempeqlist[kk]]
      ];
 eqlisteuler = Flatten[Join[eqlisteuler,tempeqlist[kk]]], 
 {kk,1,lenmyresultnew}
 ]; (* end do 111 *)
 Print["Finished the construction of system for c[i]."];
 If[debugdiscreteeuler, 
   Print["At Pt. LLL4bis, after the Do loop, eqlisteuler= "];
   Print[eqlisteuler]
   ];
If[debugdiscreteeuler,
  Print["At Pt. LLL5bis, eqlisteuler="]; 
  Print[eqlisteuler]
  ];
(* First simplification of system *)
Print["Starting the first simplification of the system for the c[i]."];
Print["The unsimplified system has ",Length[eqlisteuler]," equations."];
If[debugdiscreteeuler,
  Print["At Pt. LLL5tris, before mysimplify1, eqlisteuler= "];
  Print[eqlisteuler]
  ];
seqlisteuler = mysimplify1[eqlisteuler];
If[debugdiscreteeuler,
  Print["At Pt. LLL5quatro, after mysimplify1, seqlisteuler= "];
  Print[seqlisteuler]
  ];
Print["Finished the first simplification of the system for the c[i]."];
Print["The simplified system has ",Length[seqlisteuler]," equations."]
(* POSTPONED TILL AFTER STRIPPARAMETERS *)
]; (* end if forcediscreteeuler *)
(* WH 10/04/2003 *)
(* original shifting routine to compute the flux J_n and the obstruction *)
If[debugdiscrete, 
   Print["At Pt. HOLO.-2, rhot after replacemant of t-derivatives, rhot: "];
   Print[rhot]
  ];

(* WH 10/03/2003 *)
(* START CONSTRUCTION OF J_n based on shifting routine *)
If[debugdiscrete, 
  Print["APPLICATION OF THE SHIFTING ROUTINE TO COMPUTE FLUX JN."]
  ];
 If[debugconstructformjn, 
   Print["At Pt. HOLO.-1"];
   Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
   ];
rhot = Expand[rhot];
If[debugconstructformjn, 
    Print["At Pt. HOL0.0, original rhot= "];
    Print[rhot];
  ];
myrhot = rhot /. forwardtranslation;
If[debugconstructformjn, 
   Print["At Pt. HOL0.2, translated form of rhot is: "];
   Print[myrhot]
  ];
myvars = {n,t};
If[debugconstructformjn,
   Print["At Pt. HOL0.3, myvars= "];
   Print[myvars]
  ];
myfunclist = Table[u[i][n][t],{i,1,noeqs}]; 
myfunclist = myfunclist /. forwardtranslation;
If[debugconstructformjn,
   Print["At Pt. HOL0.4, myfunclist after translation is: "];
   Print[myfunclist]
  ];
myrhotlist = Table[myrhot[[k]], {k,1,Length[myrhot]} ]; 
If[debugconstructformjn, 
  Print["At Pt. HOL0.5, list of term in myrhot, myrhotlist= "]; 
    Print[myrhotlist]
  ];
jnlist = {}; (* list *)
obstruction = 0; (* format: term *)
Do[ (* start do loop 21 *)
term = myrhot[[kk]];
If[debugconstructformjn, 
   Print["At Pt. HOL1.0, term= "];
   Print[term]
  ];
varsinterm = Union[Cases[{term},u[i_][Plus[___,n],t],Infinity]]; 
varsinterm = Flatten[Append[varsinterm, Union[Cases[{term}, 
                    u[i_][Plus[0, n], t], Infinity]]]];
If[debugconstructformjn, 
    Print["At Pt. HOL1.1, list of variables, varsinterm= "];
    Print[varsinterm]
  ];
depvars=Table[Head[Part[varsinterm, i]], {i, Part[Dimensions[varsinterm],1]}];
If[debugconstructformjn, 
    Print["At Pt. HOL1.1.1, list of dependent variables, depvars= "];
    Print[depvars]
  ];
newdepvars = Union[depvars];
If[debugconstructformjn, 
   Print["At Pt. HOL1.1.2, list of condensed dependent "<>
         "variables, newdepvars= "];
   Print[newdepvars]
  ];
seplist = Table[Select[varsinterm, Head[#] == Part[newdepvars, i] &], 
                {i,Part[Dimensions[newdepvars], 1]}];
If[debugconstructformjn, 
    Print["At Pt. HOL1.1.3, separated list of varsinterm, seplist= "];
    Print[seplist]
  ];
subsinterm = 
  Table[Union[Cases[{Part[seplist, i]}, Plus[___, n], Infinity]], 
        {i,Part[Dimensions[seplist], 1]}];
(* Below gives list of shifted subscripts with n added in list *) 
subsinterm = 
  Table[Append[Part[subsinterm, i], 
        Union[Cases[{varsinterm}, Plus[0, n], Infinity]]], 
        {i,Part[Dimensions[seplist], 1]}];
subsinterm =
  Table[Flatten[Part[subsinterm, i]], 
        {i,Part[Dimensions[Part[subsinterm], 1], 1]}];
(* Below gives list of shifted subscripts with n added in list dropping *)
(* the last term (n) if the dimension of variables in term is less. *)  
subsinterm = Take[Part[subsinterm, 1], Part[Dimensions[Part[seplist, 1]], 1]];
If[debugconstructformjn, 
   Print["At Pt. HOL1.2, lit of subscripts, subsinterm: "];
   Print[subsinterm]
  ];
minsubsinterm = Min[subsinterm /. n->0];
If[debugconstructformjn, 
  Print["At Pt. HOL1.3, minimum of list of subscripts, minssubsinterm= "]; 
  Print[minsubsinterm]
  ];
upterms = If[minsubsinterm < 0, 
            NestList[upShift, term, minsubsinterm*-1], (* else *)
            upterms = {} 
            ];
If[debugconstructformjn, 
   Print["At Pt. HOL1.4, list of upShifted terms, upterms= "];
   Print[upterms]
   ];
downterms = If[minsubsinterm > 0, 
            NestList[downShift, term, minsubsinterm], (* else *)
            downterms = {} 
            ];
If[debugconstructformjn, 
   Print["At Pt. HOL1.5, list of downShifted terms, downterms= "];
   Print[downterms]
   ];
If[minsubsinterm==0,
   obsterm=term];
If[minsubsinterm < 0,
   obsterm=Last[upterms]];
If[minsubsinterm>0,
  obsterm=Last[downterms]];
If[debugconstructformjn, 
   Print["At Pt. HOL1.6, list of terms to go in obstruction, obsterm= "];
   Print[obsterm]
   ];
obstruction = obstruction + obsterm;
If[debugconstructformjn, 
   Print["At Pt. HOL1.7, updated obstruction= "];
   Print[obstruction]
   ];
If[minsubsinterm < 0,
   jntermlist = Drop[upterms,-1]];
If[minsubsinterm > 0,
   jntermlist = -1*Drop[downterms,1]];
If[minsubsinterm == 0,
   jntermlist = {}];
If[debugconstructformjn, 
   Print["At Pt. HOL1.8, list of terms to go into jn, jntermlist= "];
   Print[jntermlist]
  ];
(* updating the jnlist *)
jnlist = Flatten[Join[jnlist, jntermlist]];
If[debugconstructformjn, 
   Print["At Pt. HOL1.9, updated jnlist, jnlist= "];
   Print[jnlist]
   ];
formjn = jnlist /. {List -> Plus};
If[debugconstructformjn, 
   Print["At Pt. HOL1.10, updated jnlist, formjn= "];
   Print[formjn]
   ], {kk,1,Length[myrhot]}
]; (* end of do loop 21 *)
formjnplus=Part[Map[upShift,{formjn}],1];
If[debugconstructformjn, 
   Print["At Pt. HOL1.11, jn+1, formjnplus= "];
   Print[formjnplus]
   ];
difference = formjn-formjnplus;
If[debugconstructformjn, 
   Print["At Pt. HOL1.12, jn-jn+1, difference= "];
   Print[difference]
   ];
(* WH 10/04/2003 *)
(* END CONSTRUCTION OF JN based on shifting routine *)

(* Start of the construction of the linear system for c[i] based on the *)
(* shifting routine *)
(* For the shifting routine we set flag1=False and flag2=3 *)
(* Second call of constructeqlist in shifting routine *)
If[forceshifting, (* start if 450 *)
   Print["APPLICATION OF SHIFTING ROUTINE TO COMPUTE THE LINEAR SYSTEM!"];
 If[debugconstructformjn, 
    Print["At Pt. R-1."];
    Print["Second call of constructeqlist ---> in shifting routine"];
    Print["Usage: constructeqlist[rhot,False,3]"]; 
    Print["Test flag settings: flag1= ",flag1," and flag2= ",flag2]
    ];
 If[debugdiscrete,
   Print["At Pt. R0, ENTERING constructeqlist, rhot="]; 
   Print[rhot]
   ];
(* WH 06/10/2003 constructeqlist is used here also *)
(* 
 Print["Starting the derivation of the equations of the system for c[i]."]; 
*)
 myeqlist = constructeqlist[rhot,False,3];
 Print["Finished the construction of system for c[i]."];
 If[debugdiscrete,
   Print["At Pt. R1, LEAVING constructeqlist, myeqlist="]; 
   Print[myeqlist]
   ]
]; (* end if 450 *)
If[forceshifting, (* start if 33 *)
   If[debugdiscrete,
      Print["At Pt. R2"];
      Print["USED SHIFTING ROUTINE to compute myeqlist, myeqlist="]; 
      Print[myeqlist]
     ];
 Print["Starting the first simplification of the system for the c[i]."];
 Print["The unsimplified system has ",Length[myeqlist]," equations."];
 If[debugdiscrete, 
   Print["At Pt. R2bis, before mysimplify1, myeqlist= "];
   Print[myeqlist]
   ];
 seqlist = mysimplify1[myeqlist];
 If[debugdiscrete, 
   Print["At Pt. R2tris, after mysimplify1, seqlist= "];
   Print[seqlist]
   ];
 Print["Finished the first simplification of the system for the c[i]."];
 Print["The simplified system has ",Length[seqlist]," equations."]
 ]; (* end if 33 *)
If[forcediscreteeuler, (* start if 34 *)
   seqlist = seqlisteuler;
   If[debugdiscreteeuler,
      Print["At Pt. R2bis."];
      Print["USED DISCRETE EULER ROUTINE to compute seqlisteuler, "];
      Print["seqlisteuler is renamed to seqlist,  seqlist="]; 
      Print[seqlist]
      ]
  ]; (* end if 34 *)
(* Start first simplification of system with routine stripparameters *)
If[forcestripparameters, 
   Print["The following parameters are assumed to be nonzero:"];
   Print[parameters]
   ];
If[debugstripsystem, 
   Print["At Pt. S0, the list of the parameters ="];
   Print[parameters]
   ];
If[debugstripsystem && forcestripparameters,
   Print["At Pt. S1, entering stripparameters, parameters ="];
   Print[parameters]
   ];
If[debugstripsystem && forcestripparameters,
   Print["At Pt. S2, entering stripparameters, seqlist ="];
   Print[seqlist]
   ];

(* start from cleanseqlist *)
cleanseqlist = seqlist;

(* Now adjust to system format lhs == 0 *)
lencleanseqlist = Length[cleanseqlist];
(* added: test if lencleanseqlist is nonzero *)
 If[lencleanseqlist =!= 0, 
    maineqlist = Flatten[Table[picklhs[cleanseqlist,i],
    {i,1,lencleanseqlist}]]
   ];
 If[lencleanseqlist == 0, 
   maineqlist = {0 == 0}
   ];
 maineqlist = MapAll[Factor,maineqlist];
 If[debugstripsystem,
   Print["At Pt. S6, factored maineqlist="]; 
   Print[maineqlist]
   ];
(* added sorting according to length of left hand sides *)
maineqlist = MapAll[Expand,maineqlist];
maineqlist = Sort[maineqlist, Length[Part[#1,1]] < Length[Part[#2,1]]&];
maineqlist = MapAll[Factor,maineqlist]; 
Print["This is system for the c[i] (to be solved): "];
Print[" "];
Print[maineqlist];
Print[" "];
Print["List of unknown coefficients c[i]: "];
Print[unknownlist];
Print["Starting the solution of this system with ",Length[maineqlist]," "<>
      "equations."];
 If[debugsystemsolver, 
    Print["At Pt. SOL0, debugsystemsolver flag should be true, flag ="];
    Print[debugsystemsolver]
    ];
 If[debugsystemsolver, 
    Print["At Pt. SOL0bis, list of parameters ="];
    Print[parameters]
    ];
If[parameters =!= {}, (* start if 2222, then start for parameter test *)
 If[debugsystemsolver, 
    Print["At Pt. SOL1, entering analyzer, seqlist = "];
    Print[seqlist]
    ];
 If[debugsystemsolver, 
    Print["At Pt. SOL2, entering analyzer, lenseqlist = "];
    Print[lenseqlist]
    ];
  analyzelist = analyzer[seqlist,lenseqlist];
 If[debugsystemsolver, 
    Print["At Pt. SOL3, leaving analyzer, analyzelist = "];
    Print[analyzelist]
    ];
 If[debugsystemsolver, 
    Print["At Pt. SOL4, parameters = "];
    Print[parameters]
    ];
  lengthpar = Length[parameters];
 If[debugsystemsolver, 
    Print["At Pt. SOL5, lengthpar = "];
    Print[lengthpar]
    ];
  Print["The system for the coefficients c[i] has ",
         lengthpar," parameter(s)."];
  complexanalysis=lengthpar+Length[maineqlist]+Length[unknownlist];
  If[complexanalysis > 20, (* start if 44 *)
    system = MapAll[Factor,maineqlist];
    coefmatrix = coefmat[system,unknownlist];
    Print["Rank and form of rho, the system for c[i], its"];
    Print["coefficient matrix, and the lists of unknowns and"];
    Print["parameters, will all be saved in the file worklog.m for"];
    Print["further analysis. For analysis use the Mathematica functions"];
    Print["Reduce, Solve, etc. Load the file worklog.m with the command"];
    Print["<<worklog.m ."];
    Save["worklog.m",myfile,rhorank,formrho,system,coefmatrix,
    unknownlist,parameters];
    Clear[coefmatrix,system];
   ]; (* end if 44 *)
 For[i = 1,i <= lengthpar, i++, (* start for 344 *)
     syscondlist = Union[syscondlist,{parameters[[i]] != 0}];
     If[debugsystemsolver, 
       Print["At Pt. SOL6, in for loop, with i=",i," syscondlist = "];
       Print[syscondlist]
       ]
     ]; (* end for 344 *)
  syscond = Apply[And,syscondlist];
  If[debugsystemsolver, 
     Print["At Pt. SOL7, apply after And, syscond = "];
     Print[syscond]
     ];
  Print[" "];
  Print["Starting the search for compatibility conditions."];
  Print["The program will try all possible choices for non-vanishing"];
  Print["c[i]."];
  Print[" "];
  inputlist = Complement[unknownlist,analyzelist];
  If[debugsystemsolver, 
     Print["At Pt. SOL8, complement unknownlist and analyzelist, inputlist="];
     Print[inputlist]
     ];
  If[inputlist === {},
     If[debugsystemsolver, 
        Print["At Pt. SOL9, since inputlist = {}"]
        ];
    Print["All of the coefficients c[i] in the density have to vanish."];
    Print["DENSITY AND FLUX ARE ZERO!"];
    Print["The search for compatibility conditions ends."];
    ];
  While[inputlist =!= {}, (* start while 77 *)
     If[debugsystemsolver, 
        Print["At Pt. SOL10, in while loop since inputlist =!= {}"]
        ];
   inputpart = Part[inputlist,1];
   If[debugsystemsolver, 
      Print["At Pt. SOL11, part 1 of inputlist, inputpart="];
      Print[inputpart]
      ];
   inputval = inputpart == 1;
   If[debugsystemsolver, 
      Print["At Pt. SOL12, inputval="];
      Print[inputval]
      ];
   inputlist = Complement[inputlist,{inputpart}];
   If[debugsystemsolver, 
      Print["At Pt. SOL13, complement inputlist and inputpart, inputlist="];
      Print[inputlist]
      ];
   inputrule = ToRules[inputval];
   If[debugsystemsolver, 
      Print["At Pt. SOL14, after to rules, inputrule="];
      Print[inputrule]
      ];
   Print[" "];
   Print["-------------------------------------------------------"];
   Print[" "];
   Print["* Setting ",inputrule[[1]]," :"];
   Print[" "];
   Print["Computation of the compatibility conditions."];
   comcond = Eliminate[maineqlist /. inputrule, unknownlist];
   If[debugsystemsolver, 
      Print["At Pt. SOL15, after applying eliminate, comcond="];
      Print[comcond]
      ];
   Print[" "];
   comcondfac = MapAll[Factor,comcond];
   If[debugsystemsolver, 
      Print["At Pt. SOL16, after factoring, comcondfac="];
      Print[comcondfac]
      ];
   (* Remove possible duplicates *)
   If[Head[comcondfac] === Or,
   If[debugsystemsolver, 
      Print["At Pt. SOL17, since Head[comcondfac] is Or,"]
      ];
     comcondfactab = Table[Part[comcondfac,k], {k,1,Length[comcondfac]}];
     If[debugsystemsolver, 
        Print["At Pt. SOL18, table of concondfacs, comcondfactab="];
        Print[comcondfactab]
        ];
     comcondfac = Union[comcondfactab];
     If[debugsystemsolver, 
        Print["At Pt. SOL19, after union, comcondfactab="];
        Print[comcondfactab]
        ]
     ];
   Print["This is the compatibility condition:"];
   Print[" "];
   Print[comcondfac];
   sol = Reduce[comcond && syscond,parameters];
   If[debugsystemsolver, 
      Print["At Pt. SOL20, after using reduce on comcond + syscon, sol="];
      Print[sol]
      ];
   If[comcondfac === True,
     Print[" "];
     Print["The compatibility condition is satisfied without constraints"];
     Print["on the parameters."]
     ];
   If[sol =!= False, (* start if 688 *)
     If[Not[FreeQ[sol,Or]],
       lengthsol = Length[sol],
       lengthsol = 1;
       ];
     If[debugsystemsolver, 
       Print["At Pt. SOL21, number of solutions, lengthsol="];
       Print[lengthsol]
       ];
    For[j = 1,j <= lengthsol,j++, (* start for 66 *)
        If[lengthsol == 1,
          rules = ToRules[sol],
          rules = ToRules[Part[sol,j]];
          ];
        If[debugsystemsolver, 
           Print["At Pt. SOL22, in for loop, with j=",j," rules="];
           Print[rules]
           ];
        If[comcondfac =!= True,
          For[k = 1,k <= lengthpar,k++,
             parameter = parameters[[k]] /. rules;
             If[debugsystemsolver, 
           Print["At Pt. SOL23, in for loops, k=",k," and j=",j," parameter="];
                Print[parameter]
                ];
             Print[" "];
             Print["For ",parameters[[k]]," = ",parameter];
             ]
          ];
        newmaineqlist = maineqlist /. rules;
        If[debugsystemsolver, 
           Print["At Pt. SOL24, applying rules to maineqlist, newmaineqlist="];
           Print[newmaineqlist]
           ];
        solc = Flatten[Solve[newmaineqlist /. inputrule,unknownlist]];
        If[debugsystemsolver, 
           Print["At Pt. SOL25, solc before simplification, solc="];
           Print[solc]
           ];
        solc = Simplify[Union[solc,inputrule]];
        Print[" "];
        Print["The solution of the system is:"];
        Print[" "];
        Print[solc];
        solrules = Union[rules,solc];
        If[debugsystemsolver, 
           Print["At Pt. SOL26, just before evaluate, formrho="];
           Print[formrho]
           ];
        If[debugsystemsolver, 
           Print["At Pt. SOL27, just before evaluate, formjn="];
           Print[formjn]
           ];
(* WH 10/04/2003 *)
(* evaluation of formjn happens here *)
        evaluate[formrho,formjn,rhot,solrules];
        If[debugsystemsolver, 
           Print["At Pt. SOL28, just after evaluate, formrho="];
           Print[formrho]
          ];
        If[debugsystemsolver, 
           Print["At Pt. SOl29, just after evaluate, formjn="]; 
           Print[formjn]
          ];
        inputlist = 
          Intersection[inputlist,
             Union[solc[[Flatten[Position[solc,_->0]]]] /. x:(c[a_]->0)->c[a],
                   Complement[inputlist,First[#]& /@ solc]
                  ]
                       ];
        If[debugsystemsolver, 
           Print["At Pt. SOl30, inputlist="]; 
           Print[inputlist]
          ];
        Clear[rules,solc,newmaineqlist];
        ], (* closes for 66, else of if 688 *)
     Print[" "];
     Print["The system becomes inconsistent, or the compatibility conditions"];
     Print["require that one or more of the parameters are zero."];
     Print["Not acceptable!"];
     ] (* end if 688 *)
  ], (* closes while 77, else part 222 *)
  solc = Flatten[Solve[maineqlist,unknownlist]];
  Print["Solution of the system: "];
  Print[" "];
  Print[solc];
  If[debugsystemsolver, 
     Print["At Pt. SOl31, before evaluate, solc="]; 
     Print[solc]
     ];
  If[debugsystemsolver, 
     Print["At Pt. SOl32, before evaluate, formrho="]; 
     Print[formrho]
     ];
  If[debugsystemsolver, 
     Print["At Pt. SOl33, before evaluate, formjn="]; 
     Print[formjn]
     ];
  If[debugsystemsolver, 
     Print["At Pt. SOl34, before evaluate, rhot="]; 
     Print[rhot]
     ];
(* WH 10/04/2003 *)
(* evaluation of formjn happens here again *)
  evaluate[formrho,formjn,rhot,solc];
  Clear[solc];
  ] (* end if 2222 *)
]; (* end module main *)
(* ######################## E26 ################################# *)

(* ######################## B27 ################################# *)
(*****************************************************************************)
(* set-up for collecting the data in a log file, transcript of computations  *)
(*****************************************************************************)
logfile="";
OpenLog[filename_String]:=(logfile = OpenWrite[filename];
  If[logfile === $Failed, Return[]];
  AppendTo[$Echo,logfile]; AppendTo[$Output,logfile];);
CloseLog[]:=($Echo = Complement[$Echo,{logfile}];
  $Output = Complement[$Output,{logfile}];
  Close[logfile];);
(* ######################## E27 ################################# *)

(* ######################## B28 ################################# *)
(*****************************************************************************)
(* forward and backward translation rules                                    *)
(*****************************************************************************)
(* Needed for application of Discrete Euler Operator *)
(* Translation from format u[i][n+p][t] into u[i][n+p,t] and back *)
forwardtranslation = {u[a_][b__][t] -> u[a][b,t]};
backwardtranslation = {u[a_][b__,t] -> u[a][b][t]};
(* ######################## B28 ################################# *)

(* ######################## B29 ################################# *)
(*****************************************************************************)
(* Filename: DEulerD.m *)
(* Last updated: September 24, 2002 at 11:00 *)
(* Discrete Euler Operator (Variational Derivative) for DDEs (lattices). *)
(* Function name: DiscreteEulerD[...] *)
(* Written by Holly Eklund and Willy Hereman *)
(* Used: Computation of conserved densities for DDEs *)
(* Input: testexpression, and list dependent variables, and *)
(* list independent variables *)
(* Output: single list of results for the various components in list *)
(* of dependent variables *)
(*****************************************************************************)
(* DiscreteEulerD[  ]: applies the discrete Euler operator                   *)
(* (variational derivatives).                                                *)
(* Produces a single list with results                                       *)
(*****************************************************************************)
Off[Plus::write];
Off[Set::write];
DiscreteEulerD[f_, (y_)[x_, r___], w:{x_, r___}]:= Module[
  {tempf,svars,psvars,nsvars,ssubs,pssubs,ppshift,d1,nssubs,lnshift,
   d2,result1,result30,result31, result40, result40a, result41, result41a,
   eulerresult}, 
  Print["Applying DiscreteEulerD operator (discrete variational derivative)"];
  Print["w.r.t. variable ",y[x,r]," ."];
    tempf = Expand[f];
    svars = Union[Cases[{tempf},y[Plus[___,n],t],Infinity]]; 
    psvars = Union[Cases[{tempf},y[Plus[___?Positive,n],t],Infinity]]; 
    nsvars = Union[Cases[{tempf},y[Plus[___?Negative,n],t],Infinity]]; 
    ssubs = Union[Cases[{svars},Plus[___,n],Infinity]]; 
    pssubs = Map[(#+n)&, Cases[ssubs /. {n->0},s_?Positive]];
    ppshift =(pssubs-n)/.{}->{0};
    d1=Part[Dimensions[ppshift],1];
    nssubs = Map[(#+n)&, Cases[ssubs /. {n->0},s_?Negative]];
    lnshift = -1*(nssubs-n)/.{}->{0};
    d2=Part[Dimensions[lnshift],1];
    result1 = D[tempf, y[x, r]]; 
(* Applied Flatten subexpression *)
    result30 = (D[tempf, y[#1,t]] & ) /@ pssubs;
If[debugdeulerd, 
   Print["At Pt. X60, list of PARTIAL DERIVATIVES w.r.t. " <>
         "POSITIVE shifted terms, result30= "];
   Print[result30]
  ];
    result31 = (D[tempf, y[#1,t]] & ) /@ nssubs;
If[debugdeulerd, 
   Print["At Pt. X61, list of PARTIAL DERIVATIVES w.r.t. "<>
         "NEGATIVE shifted terms, result31= "];
   Print[result31]
   ];
    result40 = Table[Switch[First[ppshift],
               0, 0,
               Last[ppshift]==1, dummyfunc[result30],_, 
               Nest[dummyfunc, Part[result30, i], Part[ppshift,i]]] 
                    /. dummyfunc -> downShift, {i, d1}];
If[debugdeulerd, 
   Print["At Pt. X70, applying appropriate number of downShifts "<>
         "to appropriate part of result 30, result40= "];
   Print[result40]
  ];
    (* Added Flatten to deal with cases where result31 had one *)
    (* or less component *)
    result40a = Apply[Plus,Flatten[result40]];
If[debugdeulerd, 
   Print["At Pt. X70a, flattening and adding result40, result40a= "];
   Print[result40a]
   ];
   result41 = Table[Switch[First[lnshift],
              0, 0,
              Last[lnshift]==1, dummyfunc[result31],_, 
              Nest[dummyfunc, Part[result31, j], Part[lnshift,j]]] 
              /. dummyfunc -> upShift, {j, d2}];
If[debugdeulerd, 
   Print["At Pt. X71, applying appropriate number of upShifts "<>
         "to appropriate part of result 31, result41= "];
   Print[result41]
  ];
    (* Added Flatten to deal with cases where result31 had one *)
    (* or less component *)
    result41a = Apply[Plus,Flatten[result41]];
If[debugdeulerd, 
   Print["At Pt. X71a, flattening and adding result41, result41a= "];
   Print[result41a]
   ]; 
    eulerresult = Expand[result1 + result40a + result41a];
If[debugdeulerd, 
   Print["At Pt. Y1, after expanding all terms, eulerresult= "]; 
   Print[eulerresult]
  ];
Print["Finished with application of DiscreteEulerD."];
If[debugdiscreteeuler, 
   Print["This is the DiscreteEulerdResult for variable ",y[x,r]," :"];
   Print[eulerresult]
   ]; 
Return[eulerresult]
]; (* end module DiscreteEulerD *)
(* ######################## E29 ################################# *)

(* ######################## B30 ################################# *)
(*****************************************************************************)
(* Allows to use DiscreteEulerD below on a list w of dependent and a list v  *)
(* of independent variables                                                  *)
(*****************************************************************************)
DiscreteEulerD[f_, v:{(y_)[x_, r___], ___}, w:{x_, r___}]:= 
   (DiscreteEulerD[f, #1, w] & ) /@ v /; 
   If[Apply[And, (MatchQ[#1, _[Apply[Sequence, w]]] & ) /@ v], 
   True, Message[DiscreteEulerD::argx, w]];
DiscreteEulerD[f_, (y_)[x_], x_]:=DiscreteEulerD[f, y[x], {x}];
DiscreteEulerD[f_, v:{(y_)[x_], ___}, x_]:=DiscreteEulerD[f, v, {x}];
(* ######################## E30 ################################# *)

(* ######################## B31 ################################# *)
(* Routine to determine the maximum shifts, based on routine DEulerD.m *)
(* BEGIN NEW ROUTINE MAX SHIFTS *)
(* Routine computemaxshiftrhsdde *)
(* Computes the extremal shift of each variable *)
(* Example: if u[i][n+p,t] or u[n-q,t] appears than the extreme shift *)
(* of u[i] is the larger(max) of the two positive numbers p and q *)
computemaxshiftrhsdde[rhs_List] := Module[
{temprhseqs,svars,psvars,nsvars,ssubs,pssubs,ppshift,maxposshift,nssubs,
lnshift,maxnegshift,extremeshift,listextremeshiftsrhsdde},
temprhseqs = rhseqs /. forwardtranslation;
If[debugmaxshift, 
   Print["At pt. NCAM1, temprhseqs= "];
   Print[temprhseqs]
   ];
Do[ (* start do 100, loops over i *)
svars[i] = Union[Cases[{temprhseqs},u[i][Plus[___,n],t],Infinity]]; 
If[debugmaxshift, 
   Print["At pt. NCAM2, list all shifted variables, i=",i,",svars[i]= "];
   Print[svars[i]]
   ]; 
psvars[i] = Union[Cases[{temprhseqs},u[i][Plus[___?Positive,n],t],Infinity]]; 
If[debugmaxshift, 
   Print["At pt. NCAM3, list positive shifted variables, i=",i,",psvars[i]= "];
   Print[psvars[i]]
   ]; 
nsvars[i] = Union[Cases[{temprhseqs},u[i][Plus[___?Negative,n],t],Infinity]]; 
If[debugmaxshift, 
   Print["At pt. NCAM4, list negative shifted variables, i=",i,",nsvars[i]= "];
   Print[nsvars[i]]
   ]; 
ssubs[i] = Union[Cases[{svars[i]},Plus[___,n],Infinity]]; 
If[debugmaxshift, 
   Print["At pt. NCAM5, list of all shifted subscripts, i=",i,",ssubs[i]= "]; 
    Print[ssubs[i]]
    ];
pssubs[i] = Map[(#+n)&, Cases[ssubs[i] /. {n->0},s_?Positive]];
If[debugmaxshift, 
  Print["At pt. NCAM6, list positive shifted subscripts, i=",i,",pssubs[i]= "];
   Print[pssubs[i]]
   ];
ppshift[i] = (pssubs[i]-n) /. {}->{0};
If[debugmaxshift, 
   Print["At pt. NCAM7, list of places shifted in POSITIVE direction"<>
         ",i=",i,", ppshift[i]= "];
   Print[ppshift[i]]
   ];
maxposshift[i] = Max[ppshift[i]];
If[debugmaxshift, 
  Print["At pt. NCAM8, maximum positive shift , for i=",i,",maxposshift[i]= "];
   Print[maxposshift[i]]
   ];
nssubs[i] = Map[(#+n)&, Cases[ssubs[i] /. {n->0},s_?Negative]];
If[debugmaxshift, 
   Print["At pt. NCAM9, list negative shifted subscripts,i=",i,",nssubs[i]= "];
   Print[nssubs[i]]
   ];
lnshift[i] = -1*(nssubs[i]-n) /. {}->{0};
If[debugmaxshift, 
   Print["At pt. NCAM10, list of places shifted in POSITIVE direction,"<>
         "i=",i,",lnshift[i]= "];
   Print[lnshift[i]]
   ];
maxnegshift[i] = Max[lnshift[i]];
If[debugmaxshift, 
 Print["At pt. NCAM11, maximum negative shift , for i=",i,",maxnegshift[i]= "];
   Print[maxnegshift[i]]
   ];
(* combine negative and positive shifts and look at the largest of the two *)
extremeshift[i] = Max[Flatten[{maxposshift[i],maxnegshift[i]}]];
If[debugmaxshift, 
 Print["At pt. NCAM12, extremeshift, for variable u[",i,"],extremeshift[i]= "];
   Print[extremeshift[i]]
   ], {i,1,noeqs} 
];  (* end do 100 *)
listextremeshiftsrhsdde = Flatten[Table[extremeshift[i], {i,1,noeqs}]];
Return[listextremeshiftsrhsdde]
]; (* end module computemaxshiftrhsdde based on Euler operator *)
(* ######################## E31 ################################# *)

(* ######################## B32 ################################# *)
(* Routine computemaxshiftpowers *)
(* Hence, p is rhorank/(weightu[1]) and maximumshift = p - 1 *)
(* LATER: should be more general in the case of multi-variables u_n and v_n *) 
computemaxshiftpowers[rhorank_,weightu1_]:=Module[{maximumshiftpowers},
If[debugmaxshift, 
   Print["At pt. MAX0, selected rank of rho, rhorank = "];
   Print[rhorank]
   ];
If[debugmaxshift, 
   Print["At pt. MAX1, weight of first variable, weightu1 = "];
   Print[weightu1]
   ];
maximumshiftpowers = (rhorank/weightu1) - 1;
If[debugmaxshift, 
   Print["At pt. MAX2, maximumshiftpowers = "];
   Print[maximumshiftpowers]
  ];
Return[maximumshiftpowers]
]; (* end module computemaxshiftpowers *)
(* ######################## E32 ################################# *)

(* ----------------------------------------------------------------------- *)
(* ---> End of the auxiliary functions and procedures *)
(* ---> Start of the executable part of the program *)

(* ######################## B33 ################################# *)

(* Change by Mike C., July 1, 2003 *)
(* location of discreteHomotopyOperator *)
(* The Discrete Homotopy Operator code will be in a different directory than *)
(* the current directory.  It is up one level in its own directory.          *)

(* WH 10/04/2003 *)
(* Old loading homotopy functions (Mike C. July 1, 2003) *)
(* WH 10/04/2003 *)
(* Homotopy code is now entered in this code *)
(*
$Path = Append[$Path, "../DDE1DFlx"];
Get["loadDiscreteHomotopyFunctions.m"];
*)

(* ***************** START DISCRETE HOMOTOPY OPERATOR CODES ******** *)

(* File DDEHOPER.M: all pieces for homotopy operator in one file! *)
(* WH 10/04/2003 *)
(* Based on code of 10/01/2003 *)

(* File DDEHOPER.M: all pieces for discrete homotopy operator in one file! *)
(* WH 10/04/2003 *)
(* Based on code of 10/01/2003 *)

(* ************************* BEGIN of all ***************************** *)

(* ************* BEGIN loadDiscreteHomotopyFunctions.m ************* *)

(* Time-stamp: <Fri May 30 14:16:14 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* loadDiscreteHomotopyFunctions                                        *)
(* Purpose: To load all the functions needed for the Discrete Homotopy  *)
(*          case.                                                       *)
(* Authors: Frances Martin, Adam T. Ringler                             *)
(* Input: none                                                          *)
(* Output: none                                                         *)
(* Created: May 28, 2003 @ 1:04 PM @ CSM                                *)
(* Last Modified: May 30, 2003 @ 2:14 PM @ CSM                          *)
(************************************************************************)

(* load all necessary functions *)
(* 
<<shiftExpression.m; 
<<shiftExtract.m; 
<<maxDownShift.m;
<<maxUpShift.m;
<<discreteEulerOperator.m;
<<discreteInteriorProduct.m;
<<discreteLambdaReplace.m;
<<discreteHomotopyOperator.m;
*)

(* ************ END loadDiscreteHomotopyFunctions.m ****************** *)

(* *************** BEGIN discreteEulerOperator.m ********************** *)

(* Time-stamp: <Tue Jun  3 19:14:20 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* discreteEulerOperator[expression_, orderOfOperator_, dependVarHead_, *)
(*                       dependVarNum_];                                *)
(* Purpose: To calculate nth order Euler Operators for DDE's            *)
(* Authors: Ingo Kabirschke, Ryan Sayers                                *)
(* Input: An expression with no negative shifts: expression,            *)
(*        the order of the Euler Operator: orderOfOperator,             *)
(*        a dependent variable: dependVar                               *)
(* Output: The result of applying the nth order Euler Operator to a DDE *)
(* Created: May 28, 2003 at CSM                                         *)
(* Last Modified: June 3, 2003 at 7:12 PM at CSM                        *)
(************************************************************************)

(* Debug Flags for the DiscreteEulerOperator *)
debugDiscreteEulerOpInput = False;
debugEulerMaxShift = False;
debugEulerSum = False;
debugDiscreteEulerResult = False;

discreteEulerOperator[expression_, orderOfOperator_, 
                      dependVarHead_, dependVarNum_] :=
  Module[{maxShift, discreteEulerSum, discreteEulerResult},

    (* Debugging input *)
    If[debugDiscreteEulerOpInput,
       Print["Computing the ", orderOfOperator,
             " order discreteEulerOperator on "];
       Print[expression];
       Print["with respect to ", dependVarHead];
    ];

    (* Checks for negative shifts - if found, returns error message *)
    (* and aborts calculation                                       *)
    If[maxDownShift[expression, dependVarHead] < 0,
       Print["The discrete Euler operator needs an expression without", 
             " negative shifts."];
       Return[False];
    ];

    (* Finds the maximum positive shift in the expression *)
    maxShift = maxUpShift[expression, dependVarHead];

    (* Debugging maximum positive shift *)
    If[debugEulerMaxShift, 
      Print["The max shift is ", maxShift, "."]
    ];
   
    (* Calculate the Euler Operator before taking the partial     *)
    (* derivative - taken from Hereman's Notes (Newest), Page 31  *)
    discreteEulerSum = Sum[Binomial[k, orderOfOperator] * 
                           shiftExpression[expression, u, -k],
                          {k, orderOfOperator, maxShift}
                       ];

    (* Debugging result before partial derivative *)
    If[debugEulerSum, 
      Print["The Euler sum is "]; 
      Print[discreteEulerSum];
    ];

    (* Take the partial derivative with respect to the dependent   *)
    (* variable to complete the calculation of the Euler Operator. *)
    (* Also taken from Hereman's notes (Newest), page 31           *)
    discreteEulerResult = Expand[ D[discreteEulerSum, 
                                   dependVarHead[dependVarNum][n, t]]
                          ];

    (* Debugging result *)
     If[debugDiscreteEulerResult, 
       Print["The result of applying the ", orderOfOperator,
             " order discrete Euler operator with respect to ",
             dependVarHead[dependVarNum], " to"];
       Print[expression];
       Print["is"];
       Print[discreteEulerResult];
       ]; 
 
    (* Return nth order Euler Operator *)
    Return[discreteEulerResult];
];

(* Print["discreteEulerResult loaded successfully."]; *)

(* *************** END discreteEulerOperator.m ********************** *)

(* ***** ********* BEGIN discreteInteriorProduct.m ************* *)

(* Time-stamp: <Fri May 30 15:08:44 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* discreteInteriorProduct[expression_, dependVarHead_, dependVarNum_]  *)
(* Purpose: To calculate the discrete interior product of a given       *)
(*          expression with respect to the dependent variable           *)
(* Authors: Frances Martin, Adam Ringler, Ryan Sayers                   *)
(* Input: An expression: expression                                     *)
(*        A dependent variable head: dependVarHead                      *)
(*          Example: for u[1], u[2], u[3], the head is u.               *)
(*        The index of the dependent variable: dependVarNum             *)
(*          Example: for u[1], the dependVarNum is 1.                   *)
(* Output: The interior product of the expression                       *)
(* Created: May 28, 2003, @ 1:43 PM @ CSM                               *)
(* Last Modified: May 30, 2003 @ 3:08 PM @ CSM                          *)
(************************************************************************)

(* Debug Flags for the DiscreteInteriorProduct *)
debugDiscreteIntProdInput = False;
debugDiscreteIntProdTerms = False;
debugDiscreteIntProdResult = False;

discreteInteriorProduct[expression_, dependVarHead_, dependVarNum_] :=
  Module[{i, l, intProdResult, maxShift},

    (* Debugging the input *)
    If[debugDiscreteIntProdInput,
       Print["Applying the Discrete Interior Product to "];
       Print[expression];
    ];

    (* Finds the maximum positive shift in the expression. This will *)
    (* be used to limit the sum in the interior product calculation. *)
    maxShift = maxUpShift[expression, dependVarHead];

    (* Calculating the interior product using page 32 of Dr. Hereman notes *)
    (* expressing (D - I)^i as a binomial expansion                        *)
    intProdResult =
      Expand[
        Sum[
          Sum[
            (-1)^(i-l) * Binomial[i, l] * 
            shiftExpression[
              dependVarHead[dependVarNum][n, t] * 
              discreteEulerOperator[expression, i + 1, dependVarHead, 
                                    dependVarNum],
              dependVarHead, l],
            {l, 0, i}],
          {i, 0, maxShift - 1}]
      ];         
                         
    (* Debugging the terms from calculating the interior product *)
    (* Outputs the terms in tabular form                         *)
    If[debugDiscreteIntProdTerms,
      Print["The terms of the interior product are: "];
      Print[
        Table[
          Sum[
            (-1)^(i - l) * Binomial[i, l] * 
            shiftExpression[dependVarHead[dependVarNum][n, t] * 
                            discreteEulerOperator[expression, 
                                                  i + 1, dependVarHead, 
                                                  dependVarNum],
                            dependVarHead, l],
             {l, 0, i}],
          {i, 0, Max[{0, maxShift - 1}]}
        ]
      ]
    ];

    (* Debugging the result of applying the interior product *)
    If[debugDiscreteIntProdResult,
       Print["The Interior Product of "];
       Print[expression];
       Print[" with respect to ", dependVarHead, "[", dependVarNum, 
             "] is: "];
       Print[intProdResult];
    ];

    (* Return the interior product *)
    Return[intProdResult];
];

(* Print["discreteInteriorProduct loaded successfully."]; *)

(* *************** END discreteInteriorProduct.m ************* *)

(* ********************* BEGIN shiftExpression.m ******************** *)

(* Time-stamp: <Fri May 30 14:39:16 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* shiftExpression[expression_, dependVarHead_, shiftAmt_]              *)
(* Purpose: To shift a discrete expression by some finite amount        *)
(* Authors: Ingo Kabirschke, Ryan Sayers                                *)
(* Input: An expression: expression,                                    *) 
(*        a dependent variable head: dependVarHead,                     *)
(*          (For example, for u[1], u[2], and u[3] the head would be u) *)
(*        the distance to shift the variables: shiftAmt                 *) 
(* Output: The expression shifted by shiftAmt                           *)
(* Created: May 28, 2003                                                *)
(* Last Modified: May 30, 2003, at 2:35 PM at CSM                       *)
(************************************************************************)

(* DEBUG FLAGS *)
debugShiftInput = False;
debugShiftResult = False;

shiftExpression[expression_, dependVarHead_, shiftAmt_] := 
  Module[{shiftReplaceRule, shiftResult},
    
    (* Debugging input *)
    If[debugShiftInput, 
      Print["Shifting "];
      Print[expression];
      Print[" by ", shiftAmt];
    ];
      
    (* Define shifting rule: replace n with n + shift *)
    shiftReplaceRule = {dependVarHead[i_][n_, t_] -> 
                        dependVarHead[i][n + shiftAmt, t]};

    (* Apply shifting rule *)
    shiftResult = expression /. shiftReplaceRule;

    (* Debugging result *)
    If[debugShiftResult,
      Print["The result of shiftExpression is "];
      Print[expression];
    ];

    (* Return shifted expression *)
    Return[shiftResult];
  ];

(* Print["shiftExpression loaded successfully."]; *)

(* ********************* END shiftExpression.m ******************** *)

(* ****************** BEGIN discreteLambdaReplace.m **************** *)

(* Time-stamp: <Fri May 30 15:25:30 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* discreteLambdaReplace[expression_, dependVarHead_]                   *)
(* Purpose: To replace the dependent variable with lambda * dependVar   *)
(* Authors: Frances Martin, Kara Namanny and Ingo Kabirschke            *)
(* Input: An expression: expression,                                    *) 
(*        a dependent variable head: dependVarHead,                     *)
(*          (For example, for u[1], u[2], and u[3], the head would be u)*)
(* Output: The expression with lambda replacement                       *)
(* Created: May 29, 2003 10:35 AM at CSM                                *)
(* Last Modified: May 30, 2003 3:25 PM at MLRC                          *)
(************************************************************************)

(* Debug Flags for DiscreteLambdaReplace *)
debugDiscreteLambdaReplaceInput = False;
debugDiscreteLambdaReplaceResult = False;

discreteLambdaReplace[expression_, dependVarHead_] :=
  Module[{i, discreteLambdaReplaceRule, discreteLambdaReplaceResult},

    (* Debugging input *)
    If[debugDiscreteLambdaReplaceInput, 
      Print["Replacing the dependent variables, ", 
            dependVarHead, "[i], in the expression, "];
      Print[expression]; 
      Print[", with lambda * ", dependVarHead, "[i]"];
    ];

    (* Define the replacement rule - replaces each occurance of the    *)
    (* dependent variable head with lambda * dependVarHead             *)
    discreteLambdaReplaceRule =
      dependVarHead[i_][n_, t] -> lambda * dependVarHead[i][n, t];

    (* Apply replacement rule *) 
    discreteLambdaReplaceResult = Expand[expression /. 
                                         discreteLambdaReplaceRule];

    (* Debugging result *)
    If[debugDiscreteLambdaReplaceResult, 
      Print["The result of replacement is: "];
      Print[discreteLambdaReplaceResult];
    ];

    (* Return new expression with lambda *)
    Return[discreteLambdaReplaceResult];
  ];

(* Print["Discrete Lambda Replace loaded successfully."]; *)

(* ****************** END discreteLambdaReplace.m **************** *)

(* ********************* BEGIN discreteHomotopyOperator.m ************ *)

(* Time-stamp: <Wed Jun 11 15:46:02 MDT 2003>        *)
(************************************************************************)
(* discreteHomotopyOperator[expression_, dependVarHead_, numDependVar_] *)
(* Purpose: To calculate the Discrete Homotopy Operator of              *)
(*          an expression                                               *)
(* Authors: Ingo Kabirschke, Frances Martin, Kara Namanny               *)
(* Input: An expression: expression,                                    *) 
(*        a dependent variable head: dependVarHead,                     *)
(*          (For example, for u[1], u[2], and u[3] the head would be u) *)
(*        the # of elements in dependent variable vector: numDependVar  *) 
(* Output: The Discrete Homotopy Operator of the expression             *)
(* Created: May 29, 2003 @ 11:42 AM @ CSM                               *)
(* Last Modified: June 3, 2003, at 8:14 PM at CSM                       *)
(************************************************************************)

(* Debug Flags for the DiscreteHomotopyOperator *)
debugDiscreteHomotopyOpInput = False;
debugShiftedExpression = False;
debugDiscreteHomotopyOpResult = False;

discreteHomotopyOperator[expression_, dependVarHead_, numDependVar_] := 
  Module[{i, maximumDownShift, shiftedExpression, fullDiscreteIntProd, 
          lambdaExpression, integrand, shiftedHomotopyOpResult, 
          homotopyOpResult},

    (* Debugging input *)
    If[debugDiscreteHomotopyOpInput,
       Print["Applying the Discrete Homotopy Operator to "];
       Print[expression];
       Print["The dependent variable vector has head ", dependVarHead,
             " and ", numDependVar, " components."];
    ];

    (* Find the maximum downshift *)
    maximumDownShift = maxDownShift[expression, dependVarHead];
    
    (* Remove any negative shifts from expression, if needed. *)
    (* Ie. Shift up the amount of the maximum negative shift. *)
    If[maximumDownShift < 0,
       shiftedExpression = shiftExpression[expression, dependVarHead, 
                                    Abs[maximumDownShift]],
       (* Else (if no negative shifts) *)
       shiftedExpression = expression;
    ];

    (* Debugging the shifted expression *)
    If[debugShiftedExpression,
       Print["The shifted expression is: "];
       Print[shiftedExpression];
    ];

    (* Checking for integrability (Zeroth Euler Operator = 0)        *)
    (* If the expression is not integrable, outputs an error message *)
    (* and aborts further caculations. (Returns false)               *)
    For[i = 1, i <= numDependVar, i++,
      If[discreteEulerOperator[shiftedExpression, 0, dependVarHead, i] =!= 0,
         (* Then *)
         Print["The expression input (after shifting), "];
         Print[shiftedExpression];
         Print[" is not integrable."];
         Print["Aborting calculations."];
         Return[False];
      ];
    ];

    (* The steps below are taken from Hereman's notes (newest) page 26   *)
    (*********************************************************************)

    (* Find the interior product (little j) of the shifted expression *)
    fullDiscreteIntProd = Sum[     
                            discreteInteriorProduct[shiftedExpression, 
                                                    dependVarHead, i],
                            {i, 1, numDependVar}
                          ];

    (* Replace dependent variable with lambda * dependent variable *)
    lambdaExpression = 
      discreteLambdaReplace[fullDiscreteIntProd, dependVarHead];

    (* Divide new expression by lambda *)
    integrand = Expand[lambdaExpression / lambda];


    (* Integrate with respect to lambda from 0 to 1 to find flux (big J) *)
    shiftedHomotopyOpResult = Expand[Integrate[integrand, {lambda, 0, 1}]];

    (* If original expression was shifted to remove negative shifts, *)
    (* shift back down the same amount to get final result.          *)
    If[maximumDownShift < 0,
       homotopyOpResult = shiftExpression[shiftedHomotopyOpResult, 
                                          dependVarHead, maximumDownShift],
       (* Else, do not shift the final result *)
       homotopyOpResult = shiftedHomotopyOpResult;
    ]; 

    (* Debugging output *)
    If[debugDiscreteHomotopyOpResult,
       Print["The result of applying the Discrete Homotopy Operator is "];
       Print[homotopyOpResult];
    ];

    (* Return the result of applying the Homotopy Operator to expression *)
    Return[homotopyOpResult];
];

(* Print["Discrete Homotopy Operator loaded successfully."]; *)

(* ********************* END discreteHomotopyOperator.m ************ *)

(* ********************* BEGIN shiftExtract.m ******************** *)

(* Time-stamp: <Fri May 30 14:20:49 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* shiftExtract[expression_, dependVarHead_]                            *)
(* Purpose: To extract a list of the shifts in expression for           *)
(*          dependVarHead                                               *)
(* Authors: Ingo Kabirschke and Ryan Sayers                             *)
(* Input: An expression: expression,                                    *)
(*        the dependent variable head: dependVarHead                    *)
(*          Example: for u[1], u[2], u[3] this would be u               *)
(* Output: The list of shifts of dependVarHead in expression            *)
(* Created: May 28, 2003 at CSM                                         *)
(* Last Modified: May 30, 2003 at 2:16 PM at CSM                        *)
(************************************************************************)

(* DEBUG FLAGS for shiftExtract *)
debugShiftExtractInput = False;
debugShiftExtractOutput = False;

shiftExtract[expression_, dependVarHead_] := Module[
  {shiftExtractRule, shiftExtractCases, shiftList},

  (* Debugging input *)  
  If[debugShiftExtractInput,
     Print["The input to shiftExtract is "];
     Print[expression];
  ];

  (* Defining extraction rules that change an occurance of the dependent*)
  (* variable to a number, representing the dependent variable's shift  *)
  (* Example: u[1][n + 3, t] would be replaced by 3                     *)
  shiftExtractRule = 
    {dependVarHead[i_][n + shift_, t_] -> shift,
     dependVarHead[i_][n, t_] -> 0};

  (* Makes a list of every occurance of the dependent variables *)
  (* in the expression                                          *)
  shiftExtractCases = 
    Cases[expression, dependVarHead[i_][n_, t_], {0, Infinity}];

  (* Applying extraction rule - this results in a list of numbers that *)
  (* represent the shifts of each occurance of the dependent variable. *)
  shiftList = shiftExtractCases /. shiftExtractRule;

  (* Debugging result *)
  If[debugShiftExtractOutput,
     Print["The output of shiftExtract is "];
     Print[shiftList];
  ];

  (* Return list of shifts *)
  Return[shiftList];
];

(* Print["shiftExtract loaded successfully."]; *)

(* ********************* END shiftExtract.m ******************** *)

(* ********************* BEGIN maxDownShift.m ******************** *)

(* Time-stamp: <Fri May 30 14:32:11 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* maxDownShift[expression_, dependVarHead_]                            *)
(* Purpose: To find the maximum negative shift in a given expression    *)
(* Authors: Frances Martin, Kara Namanny, Adam Ringler                  *)
(* Input: a discrete expression: expression                             *)
(*        a dependent variable head: dependVarHead                      *)
(*          Example: for the variables u[1], u[2], u[3] the head is u   *)
(* Output: maximum negative shift in expression                         *)
(* Created: May 28, 2003 @ 10:56 AM @ CSM                               *)
(* Last Modified: May 30, 2003 @ 2:30 PM @ CSM                          *)
(************************************************************************)

(* Debug Flags for MaxDownShift *)
debugMaxDownShiftInput = False;
debugMaxDownShiftResult = False;

maxDownShift[expression_, dependVarHead_] := 
  Module[{maxNegativeShift},

  (* Debugging input *)
  If[debugMaxDownShiftInput,
    Print["Finding the maximum negative shift in "];
    Print[expression];
  ];

  (* Find the maximum negative shift in a given expression, using  *)
  (* Mathematica's Min function                                    *)
  maxNegativeShift = Min[shiftExtract[expression, dependVarHead]];

  (* Debugging result *)
  If[debugMaxDownShiftResult,
    Print["The maximum negative shift in "];
    Print[expression];
    Print[" is: ", maxNegativeShift];
  ];

  (* Return maximum negative shift in expression *)
  Return[maxNegativeShift];
];

(* Print["maxDownShift loaded successfully."]; *)

(* ********************* END maxDownShift.m ******************** *)

(* ********************* BEGIN maxUPShift.m ******************** *)

(* Time-stamp: <Fri May 30 14:29:16 Mountain Daylight Time 2003>        *)
(************************************************************************)
(* maxUpShift[expression, dependVarHead]                                *)
(* Purpose: To find the maximum positive shift in a given expression    *)
(* Authors: Frances Martin, Kara Namanny, Adam Ringler                  *)
(* Input: a discrete expression: expression                             *)
(*        a dependent variable head: dependVarHead                      *)
(*          Example: for u[1], u[2], u[3] the head is u                 *)
(* Output: The maximum positive shift in expression                     *)
(* Created: May 28, 2003 @ 10:56 AM @ CSM                               *)
(* Last Modified: May 30, 2003 @ 2:29 PM @ CSM                          *)
(************************************************************************)

(* Debug Flags for MaxUpShift *)
debugMaxUpShiftInput = False;
debugMaxUpShiftResult = False;

maxUpShift[expression_, dependVarHead_] := Module[{maxPositiveShift},

  (* Debugging input *)
  If[debugMaxUpShiftInput,
    Print["Finding the maximum positive shift in "];
    Print[expression];
  ];

  (* Finding the maximum positive shift in the given expression, using *)
  (* Mathematica's Max function                                        *)
  maxPositiveShift = Max[shiftExtract[expression, dependVarHead]];

  (* Debugging result *)
  If[debugMaxUpShiftResult,
    Print["The maximum positive shift in "];
    Print[expression];
    Print[" is: ", maxPositiveShift];
  ];

  (* Returning the maximum positive shift in the given expression *)
  Return[maxPositiveShift];
];

(* Print["maxUpShift loaded successfully."]; *)

(* ********************* END maxUPShift.m ******************** *)

Print["Discrete Homotopy Operator (October 3, 2009) loaded successfully."];

(* ************************* END of all ***************************** *)

(* ************* END DISCRETE HOMOTOPY OPERATOR CODES ************ *)
(* WH 10/04/2003 *)
(* End of discreteHomotopyOperator integration in this function *)

(* ***************** START MAIN CODE ****************** *)

Print["DDEDensityFluxV3.m (code of August 5, 2009) loaded successfully."];

menu;
OpenLog[myfile];
commentinter[];

(* HERE INITIALLY *)
(* 02/04/2003 WH/HE setting up empty lists for listallrhos, listallfluxes *)
listallrhos = {};
listallfluxes = {};

If[formrho === {}, formrho = 0];
If[ListQ[formrho] && formrho =!= {}, formrho=Part[formrho,1]]; 
Print["Working with the data file for the ",name,"."];
(* WH 10/05/2003 *)
   Print[" "];
For[i = 1,i <= noeqs,i++, (* start for 99 *)
   Print["Equation ",i," of the system with ",noeqs," equation(s):"];
(* WH 10/05/2003 *)
   Print[" "];
   Print["  ", SequenceForm[u,Subscript[i],Subscript[","],Subscript[n]]',
         " = ", subscriptform[u[i][n]'[t]]];
   Print[" "];
   ]; (* end for 99 *)

If[formrho === 0, (* point rhogiven *)
  eqlist = Flatten[Table[u[i][n]'[t],{i,1,noeqs}]];
  scalerules = scaling[eqlist];
  counternegweight = 0;
  listfreeweights = {};
  varscalelist = {};
  allchoices = {};
  Print[" "]; 
  Print["Program determines the weights of the variables (and parameters)"];
  Print["corresponding to Scale1 with w(d/dt)=1."];
  For[ (* begin for loop 1, loops over i *)
     i = 1,i <= noeqs,i++,
     weightu[i] = weightu[i] /. scalerules;
     If[i == 1, 
       Print["For the given system:"]
     ];
     Print["* weight of u[",i,"] and all its shifts is ",weightu[i],"."];
     varscalelist = Union[varscalelist,{{u[i][n][t],weightu[i]}}];
If[debugweightsY, 
   Print["At Pt. Y1, varscalelist = ", varscalelist]
   ];
     If[Not[NumberQ[weightu[i]]],
       If[Length[weightu[i]] == 1,
         listfreeweights = Union[listfreeweights,{weightu[i]}]],
       If[weightu[i] < 0,counternegweight++];
       ];
     ]; (* end for loop 1, loops over i *)
(* weightu1 is needed later, not a local variable !!! *)
(* could be generalized later for multi-variables *)
(* for now we base the power on the variable u[1] *)
weightu1 = weightu[1];
(* NOTE: At this point varscalelist is computed *)

If[debugweightsY, 
   Print["At Pt. Y2, listfreeweights = ", listfreeweights]
  ];
  lenweightpars = Length[weightpars];
  For[i = 1,i <= lenweightpars,i++, (* start for 111 *)
     weight[Part[weightpars,i]] =
             weight[Part[weightpars,i]] /. scalerules;
     Print["* weight of ",Part[weightpars,i]," is ",
           weight[Part[weightpars,i]],"."];
     varscalelist = Union[varscalelist,
                    {{Part[weightpars,i],weight[Part[weightpars,i]]}}];

If[debugweightsY, 
   Print["At Pt. Y3, varscalelist = ", varscalelist]
   ];
     If[Not[NumberQ[weight[Part[weightpars,i]]]],
       If[Length[weight[Part[weightpars,i]]] == 1,
         listfreeweights =
            Union[listfreeweights,{weight[Part[weightpars,i]]}]],
       If[weight[Part[weightpars,i]] < 0,counternegweight++];
       ];
If[debugweightsY, 
   Print["At Pt. Y4, listfreeweights = ", listfreeweights]
  ]
  ]; (* end for 111 *)

  counterfreeweight = Length[listfreeweights];
  If[counternegweight =!= 0, (* begin if point 1 *)
    Print[" "];
    Print["One or more of the weights are negative."];
    Print["Negative weights are not allowed."];
    Print["Aborting the computations!"];
    CloseLog[]; Abort[],
    If[counterfreeweight =!= 0, (* start if 111 *)
(* Changed counterfreeweight > 1 into counterfreeweight > 0, *)
(* so that user always has to give info for the weights. *)
(* Avoids the automatic choices made by the code which did not work properly *)
(* Added the info in the data file *)
(* Note:  If{counterfreeweight > 1, replaced square bracket by { *)
        If[counterfreeweight > 0, (* start if 222 *)
        Print[" "];
        Print["Two or more of weights have freedom."];
        Print["Enter your values for the free weights or type Abort[];"];
        Print["Your values for the free weights can be entered by typing"];
        Print["`weightu[label] = value' or `weight[variable] = value' and"];
        Print["putting `;' between your entries."];
(* weights are given here by user *)
(* Added the code to update the weights after the user-info is entered *)
    Input[": "];
    Print["Program continues with the following weights of variables"];
    Print["(and parameters) corresponding to Scale1 with w(d/dt)=1, "];
    Print["and info by user."];
    For[ (* begin for loop 1 *)
       i = 1,i <= noeqs,i++,
       weightu[i] = weightu[i] /. scalerules;
       If[i == 1, 
         Print["UPDATE: for the given system:"]
         ];
       Print["* UPDATED weight of u[",i,"] is ",weightu[i],"."];
       ]; (* end for loop 1 *)
  For[ (* begin for loop 2 *)
     i = 1,i <= lenweightpars,i++,
     weight[Part[weightpars,i]] =
             weight[Part[weightpars,i]] /. scalerules;
     Print["* UPDATED weight of ",Part[weightpars,i]," is ",
           weight[Part[weightpars,i]],"."]
     ], (* end for loop 2, else corresponding to counterfreeweight not > 1 *) 
       Print[" "];
       Print["Program will try to determine CHOICES for ",
             listfreeweights[[1]],":"];
       Print[" "];
       tempvarscalelist = Sort[Table[Part[varscalelist,k][[2]],
                               {k,1,Length[varscalelist]}]];
        If[debugweightsX, 
          Print["At Pt. X1, tempvarscalelist= ",  tempvarscalelist]
          ];
        tempvarscalelist =
               Complement[tempvarscalelist,Select[tempvarscalelist,NumberQ]];
        If[debugweightsX, 
           Print["At Pt. X2, tempvarscalelist= ",  tempvarscalelist]
          ];
        lentempvarscalelist = Length[tempvarscalelist];
        Do[ (* begin do loop 112 *)
          Print["* CHOICE ",k];
          Print[" "];
          Print["Solving the equation: ",Part[tempvarscalelist,k]," = 1."];
          Print[" "];
          attempt = Flatten[Solve[Part[tempvarscalelist,k] == 1]];
          If[attempt =!= {}, (* start if 444 *)
            solfreeweight = Part[listfreeweights,1] /. attempt;
            While[solfreeweight <= 0,
                 solfreeweight = solfreeweight + 1;
                 Print["Since the weight was zero or negative,"];
                 Print["it was incremented with 1."];
                 Print[" "]
                 ];
            Print[Part[listfreeweights,1]," = ",solfreeweight];
            weightrule = Part[listfreeweights,1] -> solfreeweight;
            tempvarscalelistval = tempvarscalelist /. weightrule;
            If[MemberQ[NonNegative[tempvarscalelistval],False],
              Print["Choice is rejected!"],
              varscalelistval = varscalelist /. weightrule;
              Print[" "];
              Print["List of all the variables (and parameters) with their"];
              Print["weights:"];
              Print[" "];
              Print[varscalelistval];
              Print[" "];
              testvarscalelist = Table[Part[varscalelistval,k][[2]],
                        {k,1,Length[varscalelist]}];
              If[Union[Positive[testvarscalelist]] === {False},
                Print["Since all the weights of u[i] are zero"];
                Print["this choice is rejected!"];
                Print[" "],
                allchoices = Union[allchoices,{solfreeweight}];
                ];
              ];
            ], (* end if 444 *)
            {k,1,lentempvarscalelist} ];  (* end do 112 *)
        Print[" "];
        Print["Simple POSITIVE choices are considered."];
        Print[" "];
        intchoices = Select[allchoices,IntegerQ];
        fracchoices = Complement[allchoices,intchoices];
        If[intchoices =!= {}, (* start if 555 *)
          (* Picking the minimum integer choice *)
          intchoices = Min[intchoices];
          choicerule = Part[listfreeweights,1] -> intchoices, (* else *)
          If[fracchoices =!= {},
           (* Picking the minimum fractional choice *)
           fracchoices = Min[fracchoices];
           choicerule = Part[listfreeweights,1] -> fracchoices, (* else *)
           Print["Enter your choice for the value of ",Part[listfreeweights,1],
                 "by typing its value. NO semi-colon at the end!"];
           choicerule = Part[listfreeweights,1]-> Input[": "];
           ];
          ]; (* end if 555 *)
        varscalelistchoice = varscalelist /. choicerule;
        If[Length[allchoices] > 1, (* start if 666 *)
          Print["List of all positive choices considered: "];
          Print[" "];
          Print[Union[Table[Part[listfreeweights,1] -> Part[allchoices,k],
                             {k,1,Length[allchoices]}]]];
          Print[" "];
          Print["Some of the weights of u[i], but not all, could be zero."];
          Print["In the data file you could enter your choices in the format"];
          Print["`weightu[label] = value;' or `weight[variable] = value;'."];
          Print[" "];
          Print["Program continues with the choice: ", choicerule];
          Print[" "];
          Print["corresponding to the weights:"];
          Print[" "];
          Print[varscalelistchoice];
          Print[" "];
          ]; (* end if 666 *)
        varscalelist = varscalelistchoice;
      ] (* end if 222 *)
  ] (* end if 111 *)
]; (* end if point 1 *)

(* default value for case where both forcemaximumshiftrhsdde == False and *)
(* forcemaximumshiftpowers == False *)
If[forcemaximumshiftrhsdde == False && forcemaximumshiftpowers == False,
   maximumshift = 0
];

(* if forcemaximumshiftrhsdde then the maximumshift is based on the *)
(* computed maximum of the listextremeshiftsrhsdde *)
rhseqs = Flatten[Table[Expand[D[u[i][n][t],t]], {i,1,noeqs}]];
If[forcemaximumshiftrhsdde,
   listextremeshiftsrhsdde=computemaxshiftrhsdde[rhseqs];
   If[debugmaxshift, 
     Print["At Pt. MAX00, leaving computemaxshiftrhsdde,"<>
           " listextremeshiftsrhsdde="];
     Print[listextremeshiftsrhsdde]
     ];
  maximumshift = Max[listextremeshiftsrhsdde];
  If[debugmaxshift, 
    Print["At Pt. MAX01, maximumshift="];
    Print[maximumshift]
    ];
  Print["The forced shift of length ",maximumshift," in the density "];
  Print["is based on the maximum shifts occuring in the DDE system."]
 ];

  Print[" "];
  Print["The rank of rho should be an integer multiple of the lowest weight "];
  Print["of the DEPENDENT variable(s). Fractional weights are allowed."];
  Print[" "];
  If[
    ListQ[formrho] && formrho =!= {}, formrho=Part[formrho,1]
    ];
  If[
    formrho === 0,
    If[Not[NumberQ[rhorank]], rhorank = Input["Enter the rank of rho: "]]
    ];
  If[formrho === 0,
    Print["Computation of the density (and flux) of RANK = ",rhorank,"."];
    ];
  nodims = noeqs + Length[weightpars];

(* MAXIMUM SHIFT BASED ON POWERS *)
If[forcemaximumshiftpowers, 
  maximumshift=computemaxshiftpowers[rhorank,weightu1];
  Print["The forced shift of length ",maximumshift," in the density"];
  Print["is based on the highest degree occurring in the density."]
  ];

(* Adjusting the varscalelist! Build list of all the variables *)
(* adding the up-shifts only, upto a maximumshift level *)
(* called extremeshift[i] for variable u[i]. *)
(* One possibility:  maximumshift = maximum of all extremeshifts. *)
(* Second possibility: maximumshift = p - 1, where p is highest degree. *)

If[debugmaxshift, 
   Print["At pt. MAX31, maximumshift, for all variables, maximumshift= "];
   Print[maximumshift]
   ];
  shiftedvarscalelist = varscalelist; 
  If[debugmaxshift, 
     Print["Pt. 00, before up-shifts, shiftedvarscalelist= "];
     Print[shiftedvarscalelist]
     ]; 
  newshiftedvarscalelist = shiftedvarscalelist;
  If[forceshiftsinrho, (* start if 20 *)
    Do[ 
      newshiftedvarscalelist = 
      Append[newshiftedvarscalelist, 
           Table[{Part[shiftedvarscalelist,i][[1]] /. {n -> n+pp},
                  Part[shiftedvarscalelist,i][[2]]}, 
                   {i,1,Length[shiftedvarscalelist]}]],
      {pp,1,maximumshift}
      ];
  newshiftedvarscalelist = Partition[Flatten[newshiftedvarscalelist],2],
  (* else *)
  newshiftedvarscalelist = varscalelist
  ]; (* end if 20 *)
  If[debugmaxshift, 
     Print["Pt. 01, after up-shifts, newshiftedvarscalelist= "];
     Print[newshiftedvarscalelist]
     ]; 
(* now go back to original varscalelist *)
varscalelist = newshiftedvarscalelist;
  If[debugmaxshift, 
     Print["Pt. 02, after up-shifts, varscalelist= "];
     Print[varscalelist]
     ]; 
(* Since the weighted parameters were already present we need *)
(* to remove duplicates, added line with union below *)
  varscalelist = Union[varscalelist];
  varscalelist = Reverse[Sort[varscalelist,OrderedQ[{#1[[2]],#2[[2]]}]&]];
  If[debugmaxshift, 
     Print["Pt. 03, enters the construction of rho, varscalelist= "];
     Print[varscalelist]
     ]; 
(* We either use the extra scale OR we compute the form of rho *)
(* based on the old strategy with one scale *)
(* We only apply the extra scales for which the weights of the dependent *)
(* variables are nonzero *)
testuweights = Table[zeroweightu[i] /. scalesol0, {i,1,noeqs}];
If[debugformrhoL, 
   Print["At Pt. LL0, testuweights= "];
   Print[testuweights]   
  ];
  If[debugcheckpoints,
     Print["PASSING CHECKPOINT 0"]
     ];
(* Changed the logic, to be able to force the software *)
(* to continue with either the original or new constructformrho routine *)
If[Union[testuweights] =!= {0}, (* if 5 then *)
   (* then apply newconstructformrho routine, unless ordered otherwise *)
   Print["AN EXTRA SCALE (Scale0) BASED ON w(d/dt)=0 CAN BE USED!"];
   If[debugcheckpoints,
      Print["PASSING CHECKPOINT 1 for THEN PART of if 5 "]
     ];
   If[forcesinglescale, (* then if 6 *)
     If[debugcheckpoints,
        Print["PASSING CHECKPOINT 2 for THEN PART of if 6"]
       ];
      Print["You are forcing the use of the single scale routine anyway!"];
      originalconstructformrho[nodims,varscalelist,maximumshift], 
      (* else if 6 *)
      If[debugcheckpoints,
         Print["PASSING CHECKPOINT 3 for ELSE PART of if 6"]
        ];
(* added maximumshift *)
     newconstructformrho[nodims,varscalelist,maximumshift], 
     (* undecided if 6 *)
     If[debugcheckpoints,
        Print["PASSING CHECKPOINT 4 for UNDECIDED PART of if 6"]
        ];
     newconstructformrho[nodims,varscalelist,maximumshift]
     ], (* end if 6 *)
   (* else if 5 apply original constructformrho, unless ordered otherwise *)
   Print["AN EXTRA SCALE (Scale0) BASED ON w(d/dt)=0 CAN NOT BE USED!"];
   If[debugcheckpoints,
      Print["PASSING CHECKPOINT 5 for ELSE PART of if 5"]
      ];
      If[forcemultiplescale, (* start if 7 *)
      If[debugcheckpoints,
         Print["PASSING CHECKPOINT 6 for THEN PART of if 7"]
         ];
        Print["You are forcing the use of the multiple scale routine anyway!"];
        newconstructformrho[nodims,varscalelist,maximumshift], (* else if 7 *)
       If[debugcheckpoints,
          Print["PASSING CHECKPOINT 7 for ELSE PART of if 7 "]
          ];
       originalconstructformrho[nodims,varscalelist,maximumshift], 
       (* undecided if 7 *)
       If[debugcheckpoints,
          Print["PASSING CHECKPOINT 8 for UNDECIDED PART of if 7"]
          ];
       originalconstructformrho[nodims,varscalelist,maximumshift]
       ], (* end if 7, undecided if 5 *)
       If[debugcheckpoints,
          Print["PASSING CHECKPOINT 9 for UNDECIDED PART of if 5"]
          ]; 
       If[debugcheckpoints,       
          Print[" Could not DECIDE! Fatal Error!"]
          ]
  ], (* end if 5, else point for rhogiven *)
  If[debugcheckpoints,
     Print["PASSING CHECKPOINT 9 for ELSE PART of if rhogiven"]
     ];
  formrho = Expand[formrho];
  If[Head[formrho] === Times,
    lenformrho = 1,lenformrho = Length[formrho];
    ];
  Print[" "];
  Print["The program will only test the given form of the density rho."];
  Print["No determination of scaling properties."];
  Print["You have given this form for the density rho:"];
  Print[" "];
  Print[subscriptform[formrho]];
  Print[" "];
  ]; (* end point rhogiven *)
unknownlist = Table[c[i],{i,1,lenformrho}];

(* Brought back the two lines below from code condens.m *)
parameters = Union[parameters,weightpars];
weightpars = {};

Clear[i,j,k,lenformrho,counternegweight,counterfreeweight,scalerules,eqlist,
     lenweightpars,varscalelist,tempvarscalelist,attempt,listfreeweights,
     allchoices,varscalelistval,temparscalelistval,intchoices,fracchoices,
     choicerule,weightrule,nodims,highestordereq,weightu,name,weight,
     lentempvarscalelist,solfreeweight,varscalelistchoice,testvarscalelist];
(* ######################## E33 ################################# *)

(* ######################## B34 ################################# *)
main[];
Print[" "];
Print["*************************** SUMMARY ***************************"];
Print[" "];
Print["Total CPU time used in the current session is ",Round[TimeUsed[]],
      " seconds."];
Print[" "];
Print["To see all the densities type: listallrhos," ];
Print["or subscriptform[listallrhos]."];
Print[" "];
Print["To see the last density type: rho, or subscriptform[rho]."];
Print[" "];
Print["To split the density in independent pieces type: splitrho[rho],"];
Print["or splitrho[subscriptform[rho]]."];
Print[" "];
Print["To see all the fluxes type: listallfluxes,"];
Print["or subscriptform[listallfluxes]."];
Print[" "];
Print["To see the last flux type: flux, or subscriptform[flux]."];
Print[" "];
Print["To split the flux in independent pieces type: splitflux[flux],"];
Print["or splitflux[subscriptform[flux]]."];
Print[" "];
Print["*************************** SUMMARY ***************************"];
Print[" "];
CloseLog[];
(* ######################## E34 ################################# *)

Print["You used the code DDEDensityFluxV3.m of August 5, 2009."];

(* ******************************* END ************************************* *)

