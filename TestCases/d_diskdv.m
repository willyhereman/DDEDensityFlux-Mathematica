(* d_diskdv.m *)

(* For discretization of the combined KdV-mKdV equation *)
(* derived by Herbst and Hereman *)

(* KdV subcase !!!! *)

(* Take the discretized KdV by itself: set bb = 0 *)
bb = 0;   

(* leave aa free, you could set it to 1/6 *)

(* aa = alpha, bb = beta, cc = gamma, dd = delta *)
(* set hh = 1; *)

u[1][n_]'[t]:= -Expand[ (cc + aa*u[1][n][t]+bb*u[1][n][t]^2)*( 
dd*((1/2)*u[1][n+2][t]-u[1][n+1][t]+u[1][n-1][t]-(1/2)*u[1][n-2][t])+
(aa/2)*( 
u[1][n+1][t]^2-u[1][n-1][t]^2+u[1][n][t]*(u[1][n+1][t]-u[1][n-1][t])+
u[1][n+1][t]*u[1][n+2][t]-u[1][n-1][t]*u[1][n-2][t] )+ 
(bb/2)*( u[1][n+1][t]^2*(u[1][n+2][t]+u[1][n][t])-
u[1][n-1][t]^2*(u[1][n-2][t]+u[1][n][t]) ) ) ];

noeqs = 1;
name = "Herbst/Taha Discretization of KdV Equation";
parameters = {aa};
weightpars = {cc,dd};

formrho = 0;

(*
weightu[1] = 1/2; 
weight[aa] = 0;
weight[cc] = 1/2;
weight[dd] = 1/2;
*)

(* FORCING OPTIONS *)

(* If forcesinglescale=True, only scale1 with w(d/dt)=1 is used to  *) 
(* construct the form of rho. The density is not split into pieces. *)

(* If forcemultiplescale=True, rho is constructed based on scale1,  *)
(* but split into pieces according to scale0 with w(d/dt) = 0.      *)
(* NOTE: If both are false, the code uses scale0 when appropriate.  *)

forcesinglescale   = False;
forcemultiplescale = False;

(* If forceshiftsinrho=True, the code will generate the form of rho *)
(* based on {u_n, u_{n+1}, ..., u_{n+maximumshift}}.                *)
(* To compute the maximumshift there are two options (see below).   *)
(* If forceshiftsinrho=False, no shifts on the dependent variable   *)
(* u_n will be used to generate the form of rho.                    *)

forceshiftsinrho = False; 

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

forcemaximumshiftrhsdde = False;
forcemaximumshiftpowers = False;

(* If forcediscreteeuler=True, then the linear system for the c[i] is *)
(* computed via the Discrete Variational Derivative (Euler) Operator. *) 
(* If forceshifting=True, then the original shifting algorithm is     *)
(* used to compute the linear system for the c[i] and the fluxes J_n. *)

forcediscreteeuler = True;
forceshifting      = False;

(* If forcestripparameters=True, then power of the nonzero parameters, *)
(* are removed in the factored form of the linear system for the c[i]. *) 
(* Example: aa, bb^3, aa*bb, ... are removed during simplification,    *)
(* but not factors like (aa^2-bb), (aa+bb)^4, etc. *)

forcestripparameters = True;

(* If forceextrasimplifications=True, then while generating the linear *)
(* system, equations of type c[i] == 0 are automatically applied to    *)
(* the rest of the system for the c[i]. *)

forceextrasimplifications = True;

(* d_diskdv.m *)
