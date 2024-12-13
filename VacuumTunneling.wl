(* ::Package:: *)

(* ::Section:: *)
(*BeginPackge*)


(* ::Subsection::Closed:: *)
(*BeginPackge*)


BeginPackage["VacuumTunneling`"]


(* ::Subsection:: *)
(*Public Functions*)


(* ::Subsubsection:: *)
(*Tunneling*)


(* Declare your package's public symbols here. *)

Tunneling;


Options[Tunneling] = {
  Dimension -> 4,
  TimesToFind->50,
  TimesToDefom->20,
  RelativeAccuracy->1/100,
  NumbericalPotential->False,
  NumbericalRenormalization->False,
  PointsNumber->100,
  BarrierBetweenVacuums->Null};


(* ::Subsection::Closed:: *)
(*BeginPrivate*)


Begin["`Private`"]
Off[RuleDelayed::rhs];
Off[InterpolatingFunction::dmval];
Off[NIntegrate::inumr];


(* ::Section:: *)
(*Code*)


(* ::Subsubsection:: *)
(*Tunneling*)


(* Define your public and private symbols here. *)

Tunneling[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_, OptionsPattern[]]:=Module[{b},
Which[Length[fields]==Length[Vacuum1]==Length[Vacuum2]==0,
        Is1DTunneling=True;
        b=T1[potential,renormalization,fields,Vacuum1,Vacuum2,
        OptionValue[Dimension],OptionValue[BarrierBetweenVacuums],OptionValue[TimesToFind],OptionValue[RelativeAccuracy],OptionValue[NumbericalPotential],OptionValue[NumbericalRenormalization]];
        If[b[[4]]==False,Print["Fail to find within given times. The result may be incorrect. Please increase option 'TimesToFind' or reduce 'RelativeAccuracy'"]];
        Return[{b[[1]],b[[3]]}],
      Length[fields]==Length[Vacuum1]==Length[Vacuum2]>1,
        Is1DTunneling=False;
        Which[OptionValue[NumbericalPotential]==False && OptionValue[NumbericalRenormalization]==False,
                b=T2[potential, renormalization, fields, Vacuum1, Vacuum2,
                OptionValue[Dimension],OptionValue[BarrierBetweenVacuums],OptionValue[TimesToFind],OptionValue[RelativeAccuracy],OptionValue[NumbericalPotential],OptionValue[NumbericalRenormalization],OptionValue[PointsNumber],OptionValue[TimesToDefom]],
             OptionValue[NumbericalPotential]==True || OptionValue[NumbericalRenormalization]==True,
                b=nT2[potential, renormalization, fields, Vacuum1, Vacuum2,
                OptionValue[Dimension],OptionValue[BarrierBetweenVacuums],OptionValue[TimesToFind],OptionValue[RelativeAccuracy],OptionValue[NumbericalPotential],OptionValue[NumbericalRenormalization],OptionValue[PointsNumber],OptionValue[TimesToDefom]]
        ];
        Return[{b[[2]],b[[3]]}],
      True,Print["Number of fields has to be same as dimension of vacuums."]]]


(* ::Subsubsection::Closed:: *)
(*T1*)


(* Define your public and private symbols here. *)

T1[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_,
 Dimension_, BarrierBetweenVacuums_,
 TimesToFind_, RelativeAccuracy_, NumbericalPotential_, NumbericalRenormalization_]:=
  Which[
    NumbericalPotential===False && NumbericalRenormalization===False, 
      AnalyticalPotentialAnalyticalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, BarrierBetweenVacuums, Dimension, TimesToFind, RelativeAccuracy],
    NumbericalPotential===True && NumbericalRenormalization===False, 
      NumbericalPotentialAnalyticalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, BarrierBetweenVacuums, Dimension, TimesToFind, RelativeAccuracy],
    NumbericalPotential===False && NumbericalRenormalization===True, 
      AnalyticalPotentialNumbericalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, BarrierBetweenVacuums, Dimension, TimesToFind, RelativeAccuracy],
    NumbericalPotential===True && NumbericalRenormalization===True, 
      NumbericalPotentialNumbericalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, BarrierBetweenVacuums, Dimension, TimesToFind, RelativeAccuracy],
  True,Print[ "Options 'NumbericalPotential' and 'NumbericalRenormalization' should be 'True' or 'False'." ]
];


(* ::Subsubsection:: *)
(*T2*)


T2[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_, 
 Dimension_, BarrierBetweenVacuums_,
 TimesToFind_, RelativeAccuracy_, NumbericalPotential_, NumbericalRenormalization_, nofpoints_,TimesDefom_]:=
Module[{v,z,dov,TrueVacuum,FalseVacuum,
line,vt,zt,dv,dvdf,dz,dzdf,s,
tmaxfn,amaxfn,k,n,tnpts,tpts,tds,l,tls,tps,tp2,tpp,
lastconvergence,alwaysconvergence,q,stay,staynumber,sp,rs,fps,ts,fr,npts,segments,vs,zs},

v=Function[fields,potential];
z=Function[fields,renormalization];
If[v@@Vacuum1>v@@Vacuum2,FalseVacuum=Vacuum1;TrueVacuum=Vacuum2;,FalseVacuum=Vacuum2;TrueVacuum=Vacuum1;];
dov=Norm[Vacuum1-Vacuum2];
line=Function[t,FalseVacuum+(TrueVacuum-FalseVacuum)/dov t];

vt[t_]=v@@line[t];
zt[t_]=z@@line[t];

(*dvdf=ndgrad[potential,fields];
dzdf=ndgrad[renormalization,fields];*)

dvdf=Function[fields,Evaluate@Grad[potential,fields]];
dzdf=Function[fields,Evaluate@Grad[renormalization,fields]];

s=T1[vt[t],zt[t],t,dov,0,Dimension, BarrierBetweenVacuums,TimesToFind, RelativeAccuracy, NumbericalPotential, NumbericalRenormalization];

tmaxfn=100000;
amaxfn=\[Infinity];
k=0;
n=nofpoints;
tnpts=Table[FalseVacuum+i/(n) (TrueVacuum-FalseVacuum),{i,0,n}];
tpts=tnpts;
segments=tpts-Delete[Prepend[tpts,tpts[[1]]],n+2];
tds=ArrayReduce[Norm,segments,2];
l=0;
tls=Table[l=l+tds[[i]],{i,n+1}];

tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
q=0;
staynumber=0;
alwaysconvergence=True;
While[q<TimesDefom,q++;
lastconvergence=s[[4]];
alwaysconvergence=alwaysconvergence&&lastconvergence;

{npts,stay}=defom[tpts,tp2,tpp,z,dvdf,dzdf];

If[stay==0,staynumber=0,staynumber++];
If[staynumber>=2,
sp=Interpolation[Transpose[{tls,tpts}]];
rs=Table[i*s[[2]]/100,{i,0,100}];
ts=s[[1]][rs];
fps=sp[ts];
fr=Interpolation[Transpose[{rs,fps}]];
Break[]];

segments=npts-Delete[Prepend[npts,npts[[1]]],n+2];
tds=ArrayReduce[Norm,segments,2];

l=0;
tls=Table[l=l+tds[[i]],{i,n+1}];

vs=MapThread[v,Transpose[npts]];
zs=MapThread[z,Transpose[npts]];
vt=Interpolation[Transpose[{tls,vs}]];
zt=Interpolation[Transpose[{tls,zs}]];

s=T1[vt[t],zt[t],t,Last[tls],0,Dimension, BarrierBetweenVacuums,TimesToFind, RelativeAccuracy, NumbericalPotential, NumbericalRenormalization];

tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
tpts=npts;
];
If[q>=20,
    Print["Path deformation can not complete in given times. The result may be incorrect."]
];
If[lastconvergence==False,
    Print["Reach max steps in 1d tunneling. The result may be incorrect. "],
    If[alwaysconvergence==False,
        Print["Reach accuracy in final step, but did not in some intermediate path."]
    ]
];
WriteString["stdout","\[EmptyCircle]"];
Return[{s[[1]],fr,s[[3]]}]]


(* ::Subsubsection:: *)
(*nT2*)


nT2[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_,
 Dimension_, BarrierBetweenVacuums_,
 TimesToFind_, RelativeAccuracy_, NumbericalPotential_, NumbericalRenormalization_, nofpoints_,TimesDefom_]:=
Module[{v,z,line,dov,TrueVacuum,FalseVacuum,
xt,yt,vt,zt,s,tmaxfn,amaxfn,k,n,tnpts,tpts,segments,tds,l,tls,tps,tp2,tpp,
lastconvergence,alwaysconvergence,q,stay,staynumber,sp,rs,fps,ts,fr,npts,vs,zs},

v=Function[fields,potential];
z=Function[fields,renormalization];

If[v@@Vacuum1>v@@Vacuum2,FalseVacuum=Vacuum1;TrueVacuum=Vacuum2;,FalseVacuum=Vacuum2;TrueVacuum=Vacuum1;];
dov=Norm[Vacuum1-Vacuum2];
line=Function[t,FalseVacuum+(TrueVacuum-FalseVacuum)/dov t];

vt[t_]=v@@line[t];
zt[t_]=z@@line[t];


s=T1[vt[t],zt[t],t,dov,0,Dimension, BarrierBetweenVacuums,TimesToFind, RelativeAccuracy, NumbericalPotential, NumbericalRenormalization];

tmaxfn=100000;
amaxfn=\[Infinity];
k=0;
n=nofpoints;
tnpts=Table[FalseVacuum+i/n (TrueVacuum-FalseVacuum),{i,0,n}];
tpts=tnpts;
segments=tpts-Delete[Prepend[tpts,tpts[[1]]],n+2];
tds=ArrayReduce[Norm,segments,2];
l=0;
tls=Table[l=l+tds[[i]],{i,n+1}];
tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
q=0;
staynumber=0;
alwaysconvergence=True;
While[q<TimesDefom,q++;
lastconvergence=s[[4]];
alwaysconvergence=alwaysconvergence&&lastconvergence;
{npts,stay}=ndefom[tpts,tp2,tpp,v,z,dov];
If[stay==0,staynumber=0,staynumber++];
If[staynumber>=2,
sp=Interpolation[Transpose[{tls,tpts}]];
rs=Table[i*s[[2]]/100,{i,0,100}];
ts=s[[1]][rs];
fps=sp[ts];
fr=Interpolation[Transpose[{rs,fps}]];
Break[]];

segments=npts-Delete[Prepend[npts,npts[[1]]],n+2];
tds=ArrayReduce[Norm,segments,2];
l=0;
tls=Table[l=l+tds[[i]],{i,n+1}];
vs=Quiet[MapThread[v,Transpose[npts]]];
zs=Quiet[MapThread[z,Transpose[npts]]];
vt=Interpolation[Transpose[{tls,vs}]];
zt=Interpolation[Transpose[{tls,zs}]];

s=T1[vt[t],zt[t],t,Last[tls],0,Dimension, BarrierBetweenVacuums,TimesToFind, RelativeAccuracy, NumbericalPotential, NumbericalRenormalization];

tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
tpts=npts;
];
If[q>=20,
    Print["Path deformation can not complete in given times. The result may be incorrect."]
];
If[lastconvergence==False,
    Print["Reach max steps in 1d tunneling. The result may be incorrect. "],
    If[alwaysconvergence==False,
        Print["Reach accuracy in final step, but did not in some intermediate path."]
    ]
];
WriteString["stdout","\[EmptyCircle]"];
Return[{s[[1]],fr,s[[3]]}]]


(* ::Subsubsection::Closed:: *)
(*dbdr*)


dbdr[bounce_,ts_,maxr_]:=Module[{tmax,tmin,n,sr,db},
(*ts_ is a 1 dimension array of n numbers*)
tmax=Max[bounce[0],bounce[maxr]];
tmin=Min[bounce[0],bounce[maxr]];
n=Length[ts];
db=Transpose[Table[If[tmin<=ts[[p]]<=tmax,sr=r/.FindRoot[bounce[r]==ts[[p]],{r,0,maxr},Method->"Brent"];{Module[{r},(bounce'[r])/.{r->sr}],Module[{r},(bounce''[r])/.{r->sr}]},{0,0}],{p,n}]];

Return[db]
(*return a 2*n table,[[1]] for \[Phi]'(r), [[2]] for \[Phi]''(r) *)
]


(* ::Subsubsection:: *)
(*force*)


force[points_,dtdr2_,ddtdr_,z_,dvdf_,dzdf_]:=Module[{n,dphis,dts,dphidt,edphidt,adphidt,ddphis,ddphidt,eddphidt,addphidt,gv,zs,gz,zgv,dbsdlz2,fg,fng,gvmax,fn},
n=Length[points];
dphis=Table[points[[i+1]]-points[[i]],{i,n-1}];
dts=Table[Norm[dphis[[i]]],{i,n-1}];
dphidt=dphis/dts;
edphidt=Prepend[dphidt,dphidt[[1]]];
AppendTo[dphidt,dphidt[[n-1]]];
adphidt=(dphidt+edphidt)/2;
(*caculate the average of left and right derivative*)
ddphis=Table[adphidt[[i+1]]-adphidt[[i]],{i,n-1}];
ddphidt=ddphis/dts;
eddphidt=Prepend[ddphidt,ddphidt[[1]]];
AppendTo[ddphidt,ddphidt[[n-1]]];
addphidt=(ddphidt+eddphidt)/2;
(*caculate the average of left and right second derivative*)

gv=MapThread[dvdf,Transpose[points]];
zs=MapThread[z,Transpose[points]];

gz=MapThread[dzdf,Transpose[points]];
gvmax=Max[Abs[gv]];

fg=-ddtdr*gz/zs/2-zs*gv;
(*fng=fg-Table[fg[[i]].adphidt[[i]]*adphidt[[i]],{i,n}];*)(*not fastest*)
(*fng=fg-Diagonal[fg.Transpose[adphidt]]*adphidt;*)(*fastest in this single step, but not in tatal*)
fng=fg-ArrayReduce[Total,fg*adphidt,2]*adphidt;(*fastest in total*)
fn=addphidt*dtdr2+fng;
Return[{fn,gvmax}]]


(* ::Subsubsection::Closed:: *)
(*nforce*)


nforce[points_,dtdr2_,ddtdr_,v_,z_,dov_]:=Module[
{n,nf,dphis,dts,dphidt,edphidt,adphidt,ddphis,ddphidt,eddphidt,addphidt,deltaphi,delta,ds,pointsclone,neighbourpointsr,neighbourpointsl,
vsclone,zsclone,vsr,vsl,vs,dvr,dvl,gvr,gvl,gv,zsr,zsl,zs,dzr,dzl,gzr,gzl,gz,fg,fng,gvmax,fn},
n=Length[points];
nf=Length[points[[1]]];
dphis=Table[points[[i+1]]-points[[i]],{i,n-1}];
dts=Table[Norm[dphis[[i]]],{i,n-1}];
dphidt=dphis/dts;
edphidt=Prepend[dphidt,dphidt[[1]]];
AppendTo[dphidt,dphidt[[n-1]]];
adphidt=(dphidt+edphidt)/2;
(*caculate the average of left and right derivative*)
ddphis=Table[adphidt[[i+1]]-adphidt[[i]],{i,n-1}];
ddphidt=ddphis/dts;
eddphidt=Prepend[ddphidt,ddphidt[[1]]];
AppendTo[ddphidt,ddphidt[[n-1]]];
addphidt=(ddphidt+eddphidt)/2;
(*caculate the average of left and right second derivative*)
deltaphi=dov/n;
delta=Table[ReplacePart[Table[0,nf],i->deltaphi],{i,nf}];
ds=Transpose[Table[delta,n]];
pointsclone=Table[points,nf];
neighbourpointsr=pointsclone+ds;
neighbourpointsl=pointsclone-ds;

vsr=Table[MapThread[v,Transpose[neighbourpointsr[[i]]]],{i,nf}];
vsl=Table[MapThread[v,Transpose[neighbourpointsl[[i]]]],{i,nf}];
vs=MapThread[v,Transpose[points]];
vsclone=Table[vs,nf];
dvr=vsr-vsclone;
dvl=vsclone-vsl;
gvr=dvr/deltaphi;
gvl=dvl/deltaphi;
gv=Transpose[(gvr+gvl)/2];

zsr=Table[MapThread[z,Transpose[neighbourpointsr[[i]]]],{i,nf}];
zsl=Table[MapThread[z,Transpose[neighbourpointsl[[i]]]],{i,nf}];
zs=MapThread[z,Transpose[points]];
zsclone=Table[zs,nf];
dzr=zsr-zsclone;
dzl=zsclone-zsl;
gzr=dzr/deltaphi;
gzl=dzl/deltaphi;
gz=Transpose[(gzr+gzl)/2];


gvmax=Max[Abs[gv]];

fg=-ddtdr*gz/zs/2-zs*gv;

fng=fg-ArrayReduce[Total,fg*adphidt,2]*adphidt;
fn=addphidt*dtdr2+fng;
Return[{fn,gvmax}]]


(* ::Subsubsection::Closed:: *)
(*step*)


step[points_,dtdr2_,ddtdr_,instepsize_,z_,dvdf_,dzdf_]:=Module[{minstep,stepsize,f1,fmax,fratio,points2,f2,dfmax,nextpoints},
minstep=10^-6;
stepsize=instepsize;

f1=force[points,dtdr2,ddtdr,z,dvdf,dzdf];
fmax=Max[Abs[f1[[1]]]];
fratio=fmax/f1[[2]];
While[True,
points2=points+f1[[1]]*stepsize*0.5;
f2=force[points2,dtdr2,ddtdr,z,dvdf,dzdf][[1]];
If[stepsize<=minstep,stepsize=minstep;Break[]];
dfmax=Max[Abs[f2-f1[[1]]]];
If[dfmax<fmax*0.1,Break[]];
stepsize=stepsize*0.5;];
nextpoints=points2+f2*stepsize*0.5;
Return[{nextpoints,stepsize,fratio}]]


(* ::Subsubsection::Closed:: *)
(*nstep*)


nstep[points_,dtdr2_,ddtdr_,instepsize_,v_,z_,dov_]:=Module[{minstep,stepsize,f1,fmax,fratio,points2,f2,dfmax,nextpoints},
minstep=10^-6;
stepsize=instepsize;
f1=nforce[points,dtdr2,ddtdr,v,z,dov];
fmax=Max[Abs[f1[[1]]]];
fratio=fmax/f1[[2]];
While[True,(*Print["."];*)
points2=points+f1[[1]]*stepsize*0.5;
f2=nforce[points2,dtdr2,ddtdr,v,z,dov][[1]];
If[stepsize<=minstep,stepsize=minstep;Break[]];
dfmax=Max[Abs[f2-f1[[1]]]];
If[dfmax<fmax*0.1,Break[]];
stepsize=stepsize*0.5;];
nextpoints=points2+f2*stepsize*0.5;
Return[{nextpoints,stepsize,fratio}]]


(* ::Subsubsection::Closed:: *)
(*defom*)


defom[points_,dtdr2_,ddtdr_,z_,dvdf_,dzdf_]:=Module[{n,dov,gv,gvs,gvm,stepsize,pts,j,k,nextis,fratio},

n=Length[points];
dov=Norm[Last[points]-First[points]];

gv=MapThread[dvdf,Transpose[points]];

gvs=Table[Norm[gv[[i]]],{i,n}];
gvm=Max[gvs];
stepsize=0.1*dov/gvm/n;
pts=points;
j=0;
k=0;
fratio=1;

While[j<=500&&fratio>=0.02,j++;
nextis=step[pts,dtdr2,ddtdr,stepsize,z,dvdf,dzdf];

pts=nextis[[1]];
stepsize=1.5*nextis[[2]];
fratio=nextis[[3]];
];
If[j<=1,k++];

WriteString["stdout","\:2219"];

Return[{pts,k}]]


(* ::Subsubsection::Closed:: *)
(*ndefom*)


ndefom[points_,dtdr2_,ddtdr_,v_,z_,dov_]:=Module[{n,nf,delta,ds,pointsclone,neighbourpointsr,neighbourpointsl,vsclone,deltaphi,
vsr,vsl,vs,dvr,dvl,gvr,gvl,gv,gvs,gvm,stepsize,pts,j,k,nextis,fratio},
n=Length[points];
nf=Length[points[[1]]];
deltaphi=dov/n;
delta=Table[ReplacePart[Table[0,nf],i->deltaphi],{i,nf}];

ds=Transpose[Table[delta,n]];
pointsclone=Table[points,nf];
neighbourpointsr=pointsclone+ds;
neighbourpointsl=pointsclone-ds;

vsr=Table[MapThread[v,Transpose[neighbourpointsr[[i]]]],{i,nf}];
vsl=Table[MapThread[v,Transpose[neighbourpointsl[[i]]]],{i,nf}];
vs=MapThread[v,Transpose[points]];
vsclone=Table[vs,nf];
dvr=vsr-vsclone;
dvl=vsclone-vsl;
gvr=dvr/deltaphi;
gvl=dvl/deltaphi;
gv=Transpose[(gvr+gvl)/2];

gvs=Table[Norm[gv[[i]]],{i,n}];
gvm=Max[gvs];
stepsize=0.1*dov/gvm/n;
pts=points;
j=0;
k=0;
fratio=1;
While[j<=500&&fratio>=0.02,j++;
nextis=nstep[pts,dtdr2,ddtdr,stepsize,v,z,dov];
pts=nextis[[1]];
stepsize=1.5*nextis[[2]];
fratio=nextis[[3]];
];
If[j<=1,k++];
WriteString["stdout","\:2219"];
Return[{pts,k}]]


(* ::Subsubsection::Closed:: *)
(*AnalyticalPotentialAnalyticalRnormalization*)


AnalyticalPotentialAnalyticalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_, barrierbetweenvacuums_, dimension_, times_, accuracy_]:=
Module[{OriginalFunctionOfPotential,OriginalFunctionOfRenormalization,
  domainl,domainr,barrier,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,\[Phi]max,bigv,\[Alpha],NormalizedPotential,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  OriginalFunctionOfPotential=Function[thefield,potential/.{fieldname->thefield}];
  OriginalFunctionOfRenormalization=Function[thefield,renormalization/.{fieldname->thefield}];
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<=thefield<=domainr},thefield]],barrierbetweenvacuums];
  \[Phi]max=If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],10(barrier-FalseVacuum)+FalseVacuum,TrueVacuum];
  PotentialAtFalseVacuum=OriginalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=OriginalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  bigv=0.04*Abs[FunctionOfRelativePotential@\[Phi]max];
  (*Normalize the potential. Choosing a small factor can increase calculation speed but reduce accuracy. 0.04 makes both speed and accuracy acceptable*)
  \[Alpha]=Sqrt[bigv];
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,barrier,dimension,times, accuracy];
  FieldAndDerivative=GetFieldPoints[NormalizedPotential, OriginalFunctionOfRenormalization, \[Phi]max, FalseVacuum, barrier, dimension, times, accuracy,{InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  FieldPointsOfR=Table[{FieldAndDerivative[[i,1]]/\[Alpha],FieldAndDerivative[[i,2]]},{i,InitialFieldAndBigR[[2]]}];
  TheAction=CalculateTheAction[NormalizedPotential, OriginalFunctionOfRenormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection::Closed:: *)
(*NumbericalPotentialAnalyticalRnormalization*)


NumbericalPotentialAnalyticalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_, barrierbetweenvacuums_, dimension_, times_, accuracy_]:=
Module[{OriginalExpressionOfPotential,OriginalFunctionOfPotential,AnalyticalFunctionOfPotential,OriginalFunctionOfRenormalization,
  domainl,domainr,barrier,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,\[Phi]max,bigv,\[Alpha],NormalizedPotential,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  OriginalExpressionOfPotential[fieldname_?NumericQ]=potential;
  OriginalFunctionOfPotential=Function[thefield,OriginalExpressionOfPotential[thefield]];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<thefield<domainr},thefield]],barrierbetweenvacuums];
  \[Phi]max=If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],10(barrier-FalseVacuum)+FalseVacuum,TrueVacuum];
  AnalyticalFunctionOfPotential=InterpolationInArea[OriginalFunctionOfPotential,\[Phi]max,FalseVacuum,barrier];
  OriginalFunctionOfRenormalization=Function[thefield,renormalization/.{fieldname->thefield}];
  PotentialAtFalseVacuum=AnalyticalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=AnalyticalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  bigv=0.04*Abs[FunctionOfRelativePotential@\[Phi]max];
  \[Alpha]=Sqrt[bigv];
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,barrier,dimension,times, accuracy];
  FieldAndDerivative=GetFieldPoints[NormalizedPotential, OriginalFunctionOfRenormalization, \[Phi]max, FalseVacuum, barrier, dimension, times, accuracy,{InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  FieldPointsOfR=Table[{FieldAndDerivative[[i,1]]/\[Alpha],FieldAndDerivative[[i,2]]},{i,InitialFieldAndBigR[[2]]}];
  TheAction=CalculateTheAction[NormalizedPotential, OriginalFunctionOfRenormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection::Closed:: *)
(*AnalyticalPotentialNumbericalRnormalization*)


AnalyticalPotentialNumbericalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_, barrierbetweenvacuums_, dimension_, times_, accuracy_]:=
Module[{OriginalFunctionOfPotential,OriginalFunctionOfRenormalization,OriginalExpressionOfRnormalization,AnalyticalFunctionOfRnormalization,
  domainl,domainr,barrier,\[Phi]max,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,bigv,\[Alpha],NormalizedPotential,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  OriginalFunctionOfPotential=Function[thefield,potential/.{fieldname->thefield}];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<=thefield<=domainr},thefield]],barrierbetweenvacuums];
  \[Phi]max=If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],10(barrier-FalseVacuum)+FalseVacuum,TrueVacuum];
  FunctionOfRelativePotential=OriginalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  PotentialAtFalseVacuum=OriginalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=OriginalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  OriginalExpressionOfRnormalization[fieldname_?NumericQ]=renormalization;
  OriginalFunctionOfRenormalization=Function[thefield,OriginalExpressionOfRnormalization[thefield]];
  AnalyticalFunctionOfRnormalization=InterpolationInArea[OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,barrier];
  bigv=0.04*Abs[FunctionOfRelativePotential@\[Phi]max];
  \[Alpha]=Sqrt[bigv];
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,AnalyticalFunctionOfRnormalization,\[Phi]max,FalseVacuum,barrier,dimension,times,accuracy];
  FieldAndDerivative=GetFieldPoints[NormalizedPotential, AnalyticalFunctionOfRnormalization,\[Phi]max,FalseVacuum,barrier,dimension,times,accuracy,{InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  FieldPointsOfR=Table[{FieldAndDerivative[[i,1]]/\[Alpha],FieldAndDerivative[[i,2]]},{i,InitialFieldAndBigR[[2]]}];
  TheAction=CalculateTheAction[NormalizedPotential, AnalyticalFunctionOfRnormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));
  bounce=Interpolation[FieldPointsOfR];
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection::Closed:: *)
(*NumbericalPotentialNumbericalRnormalization*)


NumbericalPotentialNumbericalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_,barrierbetweenvacuums_, dimension_, times_, accuracy_]:=
Module[{OriginalExpressionOfpotential,OriginalFunctionOfPotential,AnalyticalFunctionOfPotential,OriginalExpressionOfRnormalization,OriginalFunctionOfRenormalization,AnalyticalFunctionOfRnormalization,
  domainl,domainr,barrier,\[Phi]max,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,bigv,\[Alpha],NormalizedPotential,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  OriginalExpressionOfpotential[fieldname_?NumericQ]=potential;
  OriginalFunctionOfPotential=Function[thefield,OriginalExpressionOfpotential[thefield]];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<=thefield<=domainr},thefield]],barrierbetweenvacuums];
  \[Phi]max=If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],10(barrier-FalseVacuum)+FalseVacuum,TrueVacuum];
  AnalyticalFunctionOfPotential=InterpolationInArea[OriginalFunctionOfPotential,\[Phi]max,FalseVacuum,barrier];
  
  OriginalExpressionOfRnormalization[fieldname_?NumericQ]=renormalization;
  OriginalFunctionOfRenormalization=Function[thefield,OriginalExpressionOfRnormalization[thefield]];
  AnalyticalFunctionOfRnormalization=InterpolationInArea[OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,barrier];
  
  PotentialAtFalseVacuum=AnalyticalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=AnalyticalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  bigv=0.04*Abs[FunctionOfRelativePotential@\[Phi]max];
  \[Alpha]=Sqrt[bigv];
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,AnalyticalFunctionOfRnormalization,\[Phi]max,FalseVacuum,barrier,dimension,times, accuracy];
  FieldAndDerivative=GetFieldPoints[NormalizedPotential, AnalyticalFunctionOfRnormalization, \[Phi]max, FalseVacuum, barrier, dimension, times, accuracy,{InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  FieldPointsOfR=Table[{FieldAndDerivative[[i,1]]/\[Alpha],FieldAndDerivative[[i,2]]},{i,InitialFieldAndBigR[[2]]}];
  TheAction=CalculateTheAction[NormalizedPotential, AnalyticalFunctionOfRnormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));
  bounce=Interpolation[FieldPointsOfR];
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection::Closed:: *)
(*InterpolationInArea*)


InterpolationInArea[NumbericalFunction_,TrueVacuum_,FalseVacuum_,barrier_]:=
  Module[{n,dx,Functionpts,\[Phi]max,ifun1},
  n=1000;
  \[Phi]max=If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],10(barrier-FalseVacuum)+FalseVacuum,TrueVacuum];
  dx=(\[Phi]max-FalseVacuum)/n;
  Functionpts=ParallelTable[Quiet[{(i-2)*dx+FalseVacuum,NumbericalFunction[(i-2)*dx+FalseVacuum]}],{i,(n+3)},DistributedContexts-> {"VacuumTunneling`Private`"}];
  ifun1=Interpolation[Cases[Functionpts,{_,_?NumericQ}]];
  
  Return[ifun1]
  ];


(* ::Subsubsection::Closed:: *)
(*FindInitialFieldAndBigR*)


FindInitialFieldAndBigR[AnalyticalPotential_, AnalyticalRenormalization_, TrueVacuum_, FalseVacuum_, barrier_, dimension_, times_, accuracy_]:=
  Module[{\[Phi]max,sign,absolutefieldaccuracy,absolutederivativeaccuracy,basicdr,dr,weib,weim,\[Phi]0,\[Phi],d\[Phi],r,c\[Phi],cd\[Phi],dd\[Phi],j,convergence},
  
  \[Phi]max=TrueVacuum;
  sign=Sign[\[Phi]max-FalseVacuum];
  
  While[Block[{y},AnalyticalPotential'[y]/.{y->\[Phi]max}]*sign>0 && Block[{y},AnalyticalPotential[y]/.{y->\[Phi]max}]<Block[{y},AnalyticalPotential[y]/.{y->FalseVacuum}],
  \[Phi]max=\[Phi]max+sign*(barrier-\[Phi]max)*10^(-6)];
  (*to prevent wrong true vacuum, which makes the field rolling far away from false vacuum.*)
  
  absolutefieldaccuracy=Abs[barrier-FalseVacuum]*accuracy;
  absolutederivativeaccuracy=Abs[barrier-FalseVacuum]*accuracy^2;
  basicdr=Abs[\[Phi]max-FalseVacuum]*accuracy*0.05;
  d\[Phi]=1;
  weib=1;weim=1;
  \[Phi]=(weib barrier+weim \[Phi]max)/(weib+weim);

  While[(Abs[\[Phi]-FalseVacuum]>=absolutefieldaccuracy||Abs[d\[Phi]]>=absolutederivativeaccuracy)&&Log[2,weib+weim]<=times,
  \[Phi]0=(weib barrier+weim \[Phi]max)/(weib+weim);
  r=0;
  \[Phi]=\[Phi]0;
  d\[Phi]=0;
  j=0;
  c\[Phi]=\[Phi];
  cd\[Phi]=d\[Phi];
  dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*Block[{y},AnalyticalPotential'[y]/.{y->c\[Phi]}];
  
  While[(d\[Phi]*sign<=0||(\[Phi]-barrier)*sign>=0)&&(\[Phi]-FalseVacuum)*sign>=0,
  (*dr=If[Abs[cd\[Phi]]<=absolutederivativeaccuracy && Abs[dd\[Phi]]<=absolutederivativeaccuracy,10basicdr,basicdr];*)
  dr=basicdr;
  d\[Phi]=cd\[Phi]+dd\[Phi]*dr;
  \[Phi]=c\[Phi]+cd\[Phi] dr (*+dd\[Phi]*dr*dr/2*);
  r=r+dr;
  c\[Phi]=\[Phi];
  cd\[Phi]=d\[Phi];
  dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*Block[{y},AnalyticalPotential'[y]/.{y->c\[Phi]}]+1/2 Block[{y},AnalyticalRenormalization'[y]/.{y->c\[Phi]}]/AnalyticalRenormalization[c\[Phi]] cd\[Phi]^2-(dimension-1)/r cd\[Phi];
  j++;
  ];
  
  weib=If[d\[Phi] sign<0,2weib+1,2weib-1];
  weim=If[d\[Phi] sign<0,2weim-1,2weim+1];
  If[Is1DTunneling==True,WriteString["stdout","\:2219"]];
  ];
  
  If[Log[2,weib+weim]>times,convergence=False,convergence=True];
  Return[{\[Phi]0,j,convergence}]
  ];


(* ::Subsubsection::Closed:: *)
(*GetFieldPoints*)


GetFieldPoints[AnalyticalPotential_, AnalyticalRenormalization_, TrueVacuum_, FalseVacuum_, barrier_, dimension_, times_, accuracy_,{TheInitialField_,NumberOfTheBigR_}]:=
  Module[{dr0,absolutefieldaccuracy,absolutederivativeaccuracy,r,dr,\[Phi],d\[Phi],dd\[Phi],c\[Phi],cd\[Phi],fieldinfo,datafieldinfo,k},
  fieldinfo=Array[datafieldinfo,{NumberOfTheBigR,5}];
  dr0=Abs[FalseVacuum-TrueVacuum]*accuracy*0.05;
  
  absolutefieldaccuracy=Abs[barrier-FalseVacuum]*accuracy;
  absolutederivativeaccuracy=Abs[barrier-FalseVacuum]*accuracy^2;
  r=0;
  \[Phi]=TheInitialField;
  d\[Phi]=0;
  k=0;
  c\[Phi]=\[Phi];
  cd\[Phi]=d\[Phi];
  dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*Block[{y},AnalyticalPotential'[y]/.{y->c\[Phi]}];
  
  While[k<NumberOfTheBigR,
    (*dr=If[Abs[cd\[Phi]]<=absolutederivativeaccuracy && Abs[dd\[Phi]]<=absolutederivativeaccuracy,10dr0,dr0];*)
    dr=dr0;
    datafieldinfo[k+1,1]=r;
    datafieldinfo[k+1,2]=c\[Phi];
    datafieldinfo[k+1,3]=cd\[Phi];
    datafieldinfo[k+1,4]=dd\[Phi]; 
    datafieldinfo[k+1,5]=dr;
    (* Record current r, field values and field's derivations. *)
    d\[Phi]=cd\[Phi]+dd\[Phi]*dr;
    \[Phi]=c\[Phi]+cd\[Phi] dr(*+dd\[Phi]*dr*dr/2*);
    r=r+dr;
    (* Evolve the field and field's first derivation by dr. *)
    c\[Phi]=\[Phi];
    cd\[Phi]=d\[Phi];
    dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*Module[{y},AnalyticalPotential'[y]/.{y->c\[Phi]}]+1/2 Module[{y},AnalyticalRenormalization'[y]/.{y->c\[Phi]}]/AnalyticalRenormalization[c\[Phi]] cd\[Phi]^2-(dimension-1)/r cd\[Phi];
    (* Then the evolved field and field's first derivation become the current ones. *)
    k++;
    
  ];
  If[Is1DTunneling==True,WriteString["stdout","\[SmallCircle]"]];
  Return[fieldinfo]
  ];


(* ::Subsubsection::Closed:: *)
(*CalculateTheAction*)


CalculateTheAction[APotential_, ARenormalization_,NumberOfTheBigR_,dimension_,FieldAndDerivative_]:=
Sum[Extract[FieldAndDerivative,{q,5}]*((Extract[FieldAndDerivative,{q,1}])^(dimension-1))*(1/(2ARenormalization[Extract[FieldAndDerivative,{q,2}]])*(Extract[FieldAndDerivative,{q,3}])^2+APotential[Extract[FieldAndDerivative,{q,2}]]),{q,NumberOfTheBigR}]*(2\[Pi]^(dimension/2))/Gamma[dimension/2];


(* ::Section::Closed:: *)
(*End*)


End[]
EndPackage[]
