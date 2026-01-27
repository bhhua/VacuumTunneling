(* ::Package:: *)

(* ::Section:: *)
(*BeginPackge*)


(* ::Subsection:: *)
(*BeginPackge*)


BeginPackage["VacuumTunneling`"]


(* ::Subsection:: *)
(*Public Functions*)


(* ::Subsubsection:: *)
(*Tunneling*)


(* Declare your package's public symbols here. *)

Tunneling;
GravityTunneling;


Options[Tunneling] = {
  Dimension -> 4,
  TimesToFind->50,
  TimesToDeform->20,
  RelativeAccuracy->1/100,
  NumericalPotential->False,
  NumericalRenormalization->False,
  PointsNumber->100,
  BarrierBetweenVacuums->Null,
  StepScale->1,
  InitialPath->"line"};


Options[GravityTunneling] = {
  (*Dimension -> 4,*)
  
  TimesToFind->100,
  Tol->10^-3,
  (*Derivative target accuracy,when this accuracy is too high and StepSize is too large,
  it will be endless because it cannot reach the judgment accuracy.*)
  StepSize->10^-5,
  Barrier->Null
  (*NumbericalPotential->False,
  NumbericalRenormalization->False,*)};


(* ::Subsection:: *)
(*BeginPrivate*)


Begin["`Private`"]
Off[RuleDelayed::rhs];
Off[InterpolatingFunction::dmval];
Off[NIntegrate::inumr];


(* ::Section:: *)
(*Parameters*)


(* ::Subsection:: *)
(*The Primordial Black Holes information*)


Print["Tips: Better not to use 'GravityTunnelingG','GravityTunnelingMp','GravityTunneling\[Mu]minus','GravityTunnelingrh' as Symbol name."];
(* Declare your package's public symbols here. *)
GravityTunnelingG=1/(8\[Pi] GravityTunnelingMp^2);
GravityTunnelingMp=2.434*10^18(*GeV*);
GravityTunneling\[Mu]minus;
GravityTunnelingrh=2 GravityTunnelingG * GravityTunneling\[Mu]minus;


(* ::Section:: *)
(*Code*)


(* ::Subsubsection:: *)
(*Tunneling*)


(* Define your public and private symbols here. *)

Tunneling[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_, OptionsPattern[]]:=Module[{b},
Which[(*Length[fields]==Length[Vacuum1]==Length[Vacuum2]==0*)ArrayQ[fields]==ArrayQ[Vacuum1]==ArrayQ[Vacuum2]==False,
        Is1DTunneling=True;
        b=T1[potential,renormalization,fields,Vacuum1,Vacuum2,
        OptionValue[Dimension],OptionValue[BarrierBetweenVacuums],
        OptionValue[TimesToFind],OptionValue[RelativeAccuracy],OptionValue[StepScale],
        OptionValue[NumericalPotential],OptionValue[NumericalRenormalization]];
        If[b[[4]]==False,Print["Fail to find within given times. The result may be incorrect. Please increase option 'TimesToFind' or reduce 'RelativeAccuracy'"]];
        Return[{b[[1]],b[[3]]}],
      (ArrayQ[fields]==ArrayQ[Vacuum1]==ArrayQ[Vacuum2]==True&&ArrayQ[renormalization]==False)&&Length[fields]==Length[Vacuum1]==Length[Vacuum2]>1,
        Is1DTunneling=False;
        Which[OptionValue[NumericalPotential]==False && OptionValue[NumericalRenormalization]==False,
                b=T2[potential, renormalization, fields, Vacuum1, Vacuum2,
                OptionValue[Dimension],OptionValue[BarrierBetweenVacuums],
                OptionValue[TimesToFind],OptionValue[RelativeAccuracy],OptionValue[StepScale],
                OptionValue[NumericalPotential],OptionValue[NumericalRenormalization],
                OptionValue[PointsNumber],OptionValue[TimesToDeform]],
             OptionValue[NumericalPotential]==True || OptionValue[NumericalRenormalization]==True,
                b=nT2[potential, renormalization, fields, Vacuum1, Vacuum2,
                OptionValue[Dimension],OptionValue[BarrierBetweenVacuums],
                OptionValue[TimesToFind],OptionValue[RelativeAccuracy],OptionValue[StepScale],
                OptionValue[NumericalPotential],OptionValue[NumericalRenormalization],
                OptionValue[PointsNumber],OptionValue[TimesToDeform]]
        ];
        Return[{b[[2]],b[[3]]}],
      ArrayQ[fields]==ArrayQ[Vacuum1]==ArrayQ[Vacuum2]==ArrayQ[renormalization]==True&&Length[fields]==Length[Vacuum1]==Length[Vacuum2]==Length[renormalization]>1,
        Is1DTunneling=False;
        If[OptionValue[NumericalRenormalization]==False,IsNumericalRenormailzation={False,False},IsNumericalRenormailzation=OptionValue[NumericalRenormalization]];
        b=Tunneling2[potential,{renormalization[[1]],renormalization[[2]]},fields,Vacuum1,Vacuum2,
                     OptionValue[NumericalRenormalization], IsNumericalRenormailzation,
                     OptionValue[Dimension],OptionValue[TimesToDeform],OptionValue[PointsNumber],
                     OptionValue[TimesToFind],OptionValue[RelativeAccuracy],OptionValue[StepScale],OptionValue[InitialPath]];
      Return[{b[[1]],b[[2]],b[[3]]}],
      True,Print["Number of fields has to be same as dimension of vacuums."]]]


(* ::Subsubsection::Closed:: *)
(*GravityTunneling*)


GravityTunneling[Vin_, fields_, \[Mu]minus_, \[Phi]tv_, \[Phi]fv_, OptionsPattern[]]:=Module[
{Vi,left,right,\[Phi]bar,InitialPoint,Bounce,dM,sp,P2,P3,P4,P5,gp,B,resultfunction,gT,gme,Gravity\[Phi],Gravity\[Mu],Gravityf,Gravity\[CapitalDelta],BE},
Vi=Function[fields,Vin];
left=Min[\[Phi]tv, \[Phi]fv];
right=Max[\[Phi]tv, \[Phi]fv];
\[Phi]bar=If[OptionValue[Barrier]===Null,thefield/.Last[FindMaximum[{Vi[thefield],left<=thefield<=right},thefield]],OptionValue[Barrier]];
GravityTunneling\[Mu]minus=\[Mu]minus;
InitialPoint=FindInitialPoint[Vi,\[Phi]tv,\[Phi]fv,\[Phi]bar,OptionValue[StepSize],OptionValue[Tol],OptionValue[TimesToFind]];
Bounce=RecordBounce[Vi,\[Phi]tv,\[Phi]fv,\[Phi]bar,OptionValue[StepSize],OptionValue[Tol],InitialPoint[[1]],InitialPoint[[2]]];
dM=Last[Bounce[[4]]][[2]];
sp=ListPlot[Bounce[[1]],(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["\[Phi]",10]}];
P2=ListPlot[Bounce[[2]],(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["\[Mu]",10]}];
P3=ListPlot[Bounce[[3]],(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["\[Delta]",10]}];
(*The mass changing of the black hole is evaluate by the SdS metric by \[Mu] parameter, dM[r]=\[Mu][r]-(8\[Pi] G Vi[\[Phi][r]])/(6 G)(r^3).*)
P4=ListPlot[Bounce[[4]],(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["dM",10]}];
P5=ListPlot[Bounce[[5]],(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["f",10]},PlotRange->All];
gp=GraphicsGrid[{{P2, P5}, {P4, P3}}];
B=ThinWallAction[dM];
Gravity\[Phi]=Interpolation[Bounce[[1]]];
Gravity\[Mu]=Interpolation[Bounce[[2]]];
Gravityf[r_]:=1-(2 GravityTunnelingG Gravity\[Mu][r])/r;
Gravity\[CapitalDelta]=Interpolation[Bounce[[3]]];
(*The Metric tensor*)
gT[r_,\[Theta]_]:=({
 {Gravityf[r]Exp[2 Gravity\[CapitalDelta][r]], 0, 0, 0},
 {0, 1/Gravityf[r], 0, 0},
 {0, 0, r^2, 0},
 {0, 0, 0, r^2 Sin[\[Theta]]^2}
});
(*The Integraal Measure of four dimensional spacetimes*)
gme[r_,\[Theta]_]:=Sqrt[Det[gT[r,\[Theta]]]];
Bounce[[1]][[Length[Bounce[[1]]],1]];
BE=8\[Pi] GravityTunnelingG GravityTunneling\[Mu]minus (NIntegrate[gme[r,\[Theta]](Gravityf[r] Gravity\[Phi]'[r]^2+3Vi[Gravity\[Phi][r]]-3Vi[\[Phi]fv]),{r,GravityTunnelingrh,Bounce[[1]][[Length[Bounce[[1]]],1]]},{\[Theta],0,\[Pi]},{\[Phi],0,2\[Pi]}]);
resultfunction[x_]=Which[x=="Action",BE,
                          x=="ThinAction",B,
                          x=="ScalarBounce",sp,
                          x=="GravityBounce",gp];
Return[resultfunction];];


(* ::Subsubsection::Closed:: *)
(*T1*)


(* Define your public and private symbols here. *)

T1[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_,
 Dimension_, BarrierBetweenVacuums_,
 TimesToFind_, RelativeAccuracy_, StepScale_,
 NumbericalPotential_, NumbericalRenormalization_]:=
  Which[
    NumbericalPotential===False && NumbericalRenormalization===False, 
      AnalyticalPotentialAnalyticalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, 
                                                  BarrierBetweenVacuums, Dimension, 
                                                  TimesToFind, RelativeAccuracy, StepScale],
    NumbericalPotential===True && NumbericalRenormalization===False, 
      NumbericalPotentialAnalyticalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, 
                                                  BarrierBetweenVacuums, Dimension, 
                                                  TimesToFind, RelativeAccuracy, StepScale],
    NumbericalPotential===False && NumbericalRenormalization===True, 
      AnalyticalPotentialNumbericalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, 
                                                  BarrierBetweenVacuums, Dimension, 
                                                  TimesToFind, RelativeAccuracy, StepScale],
    NumbericalPotential===True && NumbericalRenormalization===True, 
      NumbericalPotentialNumbericalRnormalization[potential, renormalization, fields, Vacuum1, Vacuum2, 
                                                  BarrierBetweenVacuums, Dimension, 
                                                  TimesToFind, RelativeAccuracy, StepScale],
  True,Print[ "Options 'NumericalPotential' and 'NumericalRenormalization' should be 'True' or 'False'." ]
];


(* ::Subsubsection::Closed:: *)
(*T2*)


T2[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_, 
 Dimension_, BarrierBetweenVacuums_,
 TimesToFind_, RelativeAccuracy_, StepScale_, 
 NumbericalPotential_, NumbericalRenormalization_, nofpoints_,TimesDefom_]:=
Module[{v,z,dov,TrueVacuum,FalseVacuum,
line,vt,barrier1d,zt,dv,dvdf,dz,dzdf,s,
tmaxfn,amaxfn,k,n,tnpts,tpts,tds,l,tls,tps,tp2,tpp,
lastconvergence,alwaysconvergence,q,stay,staynumber,sp,rs,fps,ts,fr,npts,segments,vs,zs},

v=Function[fields,potential];
z=Function[fields,renormalization];
If[v@@Vacuum1>v@@Vacuum2,FalseVacuum=Vacuum1;TrueVacuum=Vacuum2;,FalseVacuum=Vacuum2;TrueVacuum=Vacuum1;];
dov=Norm[Vacuum1-Vacuum2];
line=Function[t,FalseVacuum+(TrueVacuum-FalseVacuum)/dov t];
spl=Function[t,FalseVacuum+(TrueVacuum-FalseVacuum)/dov t];

vt[t_]=v@@line[t];
barrier1d = t/.Last[FindMaximum[{vt[t],0<t<dov},t]];
If[dov>10 Norm[line[barrier1d]-FalseVacuum],TrueVacuum=line[10barrier1d];dov=10barrier1d];

zt[t_]=z@@line[t];

dvdf=Function[fields,Evaluate@Grad[potential,fields]];
dzdf=Function[fields,Evaluate@Grad[renormalization,fields]];

s=T1[vt[t],zt[t],t,dov,0,
    Dimension, BarrierBetweenVacuums,
    TimesToFind, RelativeAccuracy, StepScale, 
    NumbericalPotential, NumbericalRenormalization];
Print[s[[3]]/550000];
Print[Plot[s[[1]][x],{x,0,s[[2]]}]];
Print[s[[1]]];
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
df1={};
While[q<TimesDefom,q++;
lastconvergence=s[[4]];
alwaysconvergence=alwaysconvergence&&lastconvergence;

{npts,stay}=defom[tpts,tp2,tpp,z,dvdf,dzdf];
spl=Interpolation[Transpose[{tls,tpts}]];
AppendTo[df1,npts-tpts];

If[stay==0,staynumber=0,staynumber++];
If[staynumber>=2||q>=TimesDefom,
sp=Interpolation[Transpose[{tls,tpts}]];
rs=Table[i*s[[2]]/100,{i,0,100}];
ts=s[[1]][rs];
fps=sp[ts];
fr=Interpolation[Transpose[{rs,fps}]];
Break[]];

segments=npts-Delete[Prepend[npts,npts[[1]]],n+2];
tds=ArrayReduce[Norm,segments,2];

l=0.;
tls=Table[l=l+tds[[i]],{i,n+1}];

vs=MapThread[v,Transpose[npts]];
zs=MapThread[z,Transpose[npts]];
vt=Interpolation[Transpose[{tls,vs}]];
zt=Interpolation[Transpose[{tls,zs}]];

s=T1[vt[t],zt[t],t,Last[tls],0.,
     Dimension, BarrierBetweenVacuums,
     TimesToFind, RelativeAccuracy, StepScale,
     NumbericalPotential, NumbericalRenormalization];
Print[s[[3]]/550000];
Print[Plot[s[[1]][x],{x,0,s[[2]]}]];
Print[s[[1]]];
tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
tpts=npts;
];

(*extra test*)
segments=npts-Delete[Prepend[npts,npts[[1]]],n+2];
tds=ArrayReduce[Norm,segments,2];

l=0.;
tls=Table[l=l+tds[[i]],{i,n+1}];

vs=MapThread[v,Transpose[npts]];
zs=MapThread[z,Transpose[npts]];
vt=Interpolation[Transpose[{tls,vs}]];
zt=Interpolation[Transpose[{tls,zs}]];

s=T1[vt[t],zt[t],t,Last[tls],0.,
     Dimension, BarrierBetweenVacuums,
     TimesToFind, RelativeAccuracy, StepScale,
     NumbericalPotential, NumbericalRenormalization];
Print[s[[3]]/550000];
Print[Plot[s[[1]][x],{x,0,s[[2]]}]];
Print[s[[1]]];
tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
tpts=npts;
npts=tpts+df1[[2]];

sp=Interpolation[Transpose[{tls,tpts}]];
rs=Table[i*s[[2]]/100,{i,0,100}];
ts=s[[1]][rs];
fps=sp[ts];
fr=Interpolation[Transpose[{rs,fps}]];

segments=npts-Delete[Prepend[npts,npts[[1]]],n+2];
tds=ArrayReduce[Norm,segments,2];

l=0.;
tls=Table[l=l+tds[[i]],{i,n+1}];

vs=MapThread[v,Transpose[npts]];
zs=MapThread[z,Transpose[npts]];
vt=Interpolation[Transpose[{tls,vs}]];
zt=Interpolation[Transpose[{tls,zs}]];

s=T1[vt[t],zt[t],t,Last[tls],0.,
     Dimension, BarrierBetweenVacuums,
     TimesToFind, RelativeAccuracy, StepScale,
     NumbericalPotential, NumbericalRenormalization];
Print[s[[3]]/550000];
Print[Plot[s[[1]][x],{x,0,s[[2]]}]];
Print[s[[1]]];
tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
tpts=npts;

(*end extra test*)
If[q>=TimesDefom,
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
(*nT2*)


nT2[potential_, renormalization_, fields_, Vacuum1_, Vacuum2_,
 Dimension_, BarrierBetweenVacuums_,
 TimesToFind_, RelativeAccuracy_, StepScale_,
 NumbericalPotential_, NumbericalRenormalization_, nofpoints_,TimesDefom_]:=
Module[{v,z,line,dov,TrueVacuum,FalseVacuum,lxx,
xt,yt,vt,barrier1d,zt,s,tmaxfn,amaxfn,k,n,tnpts,tpts,segments,tds,l,tls,tps,tp2,tpp,
lastconvergence,alwaysconvergence,q,stay,staynumber,sp,rs,fps,ts,fr,npts,vs,zs},

v=Function[fields,potential];
z=Function[fields,renormalization];

If[v@@Vacuum1>v@@Vacuum2,FalseVacuum=Vacuum1;TrueVacuum=Vacuum2;,FalseVacuum=Vacuum2;TrueVacuum=Vacuum1;];
dov=Norm[Vacuum1-Vacuum2];
line=Function[t,FalseVacuum+(TrueVacuum-FalseVacuum)/dov t];

lxx[x_]:=(TrueVacuum[[2]]-FalseVacuum[[2]])/(TrueVacuum[[1]]-FalseVacuum[[1]])*(x-FalseVacuum[[1]])+FalseVacuum[[2]];

vt[t_]=v@@line[t];

barrier1d = t/.Last[FindMaximum[{vt[t],0<t<dov},t]];

If[dov>10 Norm[line[barrier1d]-FalseVacuum],TrueVacuum=line[10barrier1d];dov=10barrier1d];

zt[t_]=z@@line[t];


s=T1[vt[t],zt[t],t,dov,0.,
     Dimension, BarrierBetweenVacuums,
     TimesToFind, RelativeAccuracy, StepScale,
     NumbericalPotential, NumbericalRenormalization];

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

s=T1[vt[t],zt[t],t,Last[tls],0.,
     Dimension, BarrierBetweenVacuums,
     TimesToFind, RelativeAccuracy, StepScale,
     NumbericalPotential, NumbericalRenormalization];

tps=dbdr[s[[1]],tls,s[[2]]];
tp2=tps[[1]]*tps[[1]];
tpp=tps[[2]];
tpts=npts;
];
If[q>=TimesDefom,
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
(*Tunneling2*)


Tunneling2[potential_,{renormalization1_,renormalization2_},fieldnames_,vacuum1_,vacuum2_,NumbericalPotential_, NumbericalRenormalization_,dimension_,timestodeform_,nofpoints_,timestofind_,accuracy_,StepScale_,initialpath_]:=Module[
{originv,originz1,originz2,falsevacuum,truevacuum,analyticalv,vfvto0,bigv,vscaled,dvscaled1,dvscaled2,analyticalz1,analyticalz2,dz1d1,dz1d2,dz2d1,dz2d2,dov,line1,line2,linep1,linep2,vlinep,truev,pointslist,xlist,nowdeform,staynumber,scaledbounce,dr,fieldinfo,rilist,dxdrlist,dxdr2list,nextpointslist,stay,dphilist,dxlist,integratedx,scalealpha,truebounce,mulitfieldbounce,trueaction,trueactionz1,result},
  originv=Function[Evaluate[fieldnames],potential];
  originz1=Function[Evaluate[fieldnames],renormalization1];
  originz2=Function[Evaluate[fieldnames],renormalization2];
  
  Which[
  NumericQ[originv@@vacuum1]===False||NumericQ[originv@@vacuum2]===False,Print["Check the potential values at VEVs."];Return[0],
    originv@@vacuum1>originv@@vacuum2,falsevacuum=vacuum1;truevacuum=vacuum2,
    originv@@vacuum1<originv@@vacuum2,falsevacuum=vacuum2;truevacuum=vacuum1
  ];
  
  If[NumbericalPotential===False,
    analyticalv=originv,
    Return[0]
  ];
  
  vfvto0=(analyticalv[##]-analyticalv@@falsevacuum)&;
  bigv=Abs[vfvto0@@truevacuum]/0.04;
  scalealpha=Sqrt[bigv];
  vscaled=(vfvto0[##]/bigv)&;
  dvscaled1=Function[{field1,field2},Evaluate[D[vscaled[field1,field2],field1]]];
  dvscaled2=Function[{field1,field2},Evaluate[D[vscaled[field1,field2],field2]]];
  
  
  
  If[NumbericalRenormalization[[1]]===False,
    analyticalz1=originz1,
    Return[0]
  ];
  If[NumbericalRenormalization[[2]]===False,
    analyticalz2=originz2,
    Return[0]
  ];
  
  dz1d1=Function[{field1,field2},Evaluate[D[analyticalz1[field1,field2],field1]]];
  dz1d2=Function[{field1,field2},Evaluate[D[analyticalz1[field1,field2],field2]]];
  dz2d1=Function[{field1,field2},Evaluate[D[analyticalz2[field1,field2],field1]]];
  dz2d2=Function[{field1,field2},Evaluate[D[analyticalz2[field1,field2],field2]]];
  
  dov=Norm[truevacuum-falsevacuum];
  
  If[initialpath=="gradient",
    pointslist=Table[{(truevacuum[[1]]-falsevacuum[[1]])i/nofpoints+falsevacuum[[1]],y/.FindRoot[dvscaled2[(truevacuum[[1]]-falsevacuum[[1]])i/nofpoints+falsevacuum[[1]],y]==0,{y,(falsevacuum[[2]]+truevacuum[[2]])/2}]},{i,0,nofpoints-1}];
    dphilist=Delete[pointslist,1]-Delete[pointslist,-1];
    dxlist=Sqrt[Total[Transpose[dphilist*dphilist]]];
    PrependTo[dxlist,0];
    integratedx=0;
    xlist=Table[integratedx=integratedx+dxlist[[i]],{i,nofpoints}];
    line1=Interpolation[Transpose[{xlist,Transpose[pointslist][[1]]}]];
    line2=Interpolation[Transpose[{xlist,Transpose[pointslist][[2]]}]];
    linep1=line1'[#]&;
    linep2=line2'[#]&;
    ,
    line1[parameterx_]:=(truevacuum[[1]]-falsevacuum[[1]])parameterx/dov+falsevacuum[[1]];
    line2[parameterx_]:=(truevacuum[[2]]-falsevacuum[[2]])parameterx/dov+falsevacuum[[2]];
    linep1[parameterx_]:=(truevacuum[[1]]-falsevacuum[[1]])/dov;
    linep2[parameterx_]:=(truevacuum[[2]]-falsevacuum[[2]])/dov;
    pointslist=Table[{line1[i/(nofpoints-1)dov],line2[i/(nofpoints-1)dov]},{i,0,nofpoints-1}];
    xlist=Table[i/(nofpoints-1)dov,{i,0,nofpoints-1}];
  ];
  
  
  nowdeform=0;
  staynumber=0;
  
  (*originpoints=pointslist;
  Clear[dfdf];*)
  
  While[nowdeform<timestodeform,
  nowdeform++;
  nowd=nowdeform;
  vlinep=(dvscaled1[line1[#],line2[#]]*linep1[#]+dvscaled1[line1[#],line2[#]]*linep1[#])&;
  truev=parax/.FindRoot[vlinep[parax]==0,{parax,xlist[[-1]]}];
  (*Print[truev-xlist[[-1]]];*)
  
  scaledbounce=PathBounce2[vscaled,{dvscaled1,dvscaled2},{analyticalz1,analyticalz2},{{dz1d1,dz1d2},{dz2d1,dz2d2}},{line1,line2},{linep1,linep2},truev,0,bigv,dimension,timestofind,accuracy,StepScale];
  If[scaledbounce[[3]]===False,Print["Reach max steps in 1d tunneling. Return the last result."];Break[]];
  
  dr=scaledbounce[[1]];
  fieldinfo=scaledbounce[[2]];
  
  rilist=Table[
    Which[
      xlist[[i]]>=fieldinfo[[1,2]],0,
      xlist[[i]]<=fieldinfo[[-1,2]],0,
      True,FirstPosition[fieldinfo,{_,parameterx_,_}/;parameterx<=xlist[[i]]][[1]]
    ],
  {i,nofpoints}];
  
  dxdrlist=Table[
    If[
      rilist[[i]]==0,0,
      (fieldinfo[[rilist[[i]],3]]+fieldinfo[[rilist[[i]]-1,3]])/2
    ],
  {i,nofpoints}];
  dxdr2list=dxdrlist*dxdrlist;
  
  {nextpointslist,stay}=defom2[pointslist,nofpoints,dxdr2list,{dvscaled1,dvscaled2},{analyticalz1,analyticalz2},{{dz1d1,dz1d2},{dz2d1,dz2d2}}];
  
  
  If[stay==0,staynumber=0,staynumber++];
  If[nowdeform>=timestodeform,Print["Reach max steps in path deformation. Return the last result."];Break[]];
  If[staynumber>=2,WriteString["stdout","\[SmallCircle]"];Break[]];
  (*truebounce=Interpolation[Transpose[{Transpose[fieldinfo][[1]]/scalealpha,Transpose[fieldinfo][[2]]}]];
  trueaction=CalculateTheAction2[dimension,dr,fieldinfo,vscaled,{analyticalz1,analyticalz2},{line1,line2},{linep1,linep2}]/(scalealpha^(dimension-2));*)
  
  pointslist=nextpointslist;
  dphilist=Delete[pointslist,1]-Delete[pointslist,-1];
  dxlist=Sqrt[Total[Transpose[dphilist*dphilist]]];
  PrependTo[dxlist,0];
  integratedx=0;
  xlist=Table[integratedx=integratedx+dxlist[[i]],{i,nofpoints}];
  
  line1=Interpolation[Transpose[{xlist,Transpose[pointslist][[1]]}]];
  
  line2=Interpolation[Transpose[{xlist,Transpose[pointslist][[2]]}]];
  
  linep1=line1'[#]&;
  linep2=line2'[#]&;
  
  ];
  
  (*extra test*)
  
  
  
  (*If[scaledbounce[[3]]===False,Return[0]];
  
  lastpoints=pointslist;
  pathlist={};
  dfdf=(lastpoints-originpoints)/20;
  actionlist={};
  deformtimes={};
  pointslist=originpoints;
  
  For[qqqq=0,qqqq<=40,qqqq++,
  nextpointslist=pointslist+dfdf;
  
  (*Print[ListPlot[pointslist]];*)
  pointslist=nextpointslist;
  
  dphilist=Delete[pointslist,1]-Delete[pointslist,-1];
  dxlist=Sqrt[Total[Transpose[dphilist*dphilist]]];
  PrependTo[dxlist,0];
  integratedx=0;
  xlist=Table[integratedx=integratedx+dxlist[[i]],{i,nofpoints}];
  
  line1=Interpolation[Transpose[{xlist,Transpose[pointslist][[1]]}]];
  
  line2=Interpolation[Transpose[{xlist,Transpose[pointslist][[2]]}]];
  
  linep1=line1'[#]&;
  linep2=line2'[#]&;
  vlinep=(dvscaled1[line1[#],line2[#]]*linep1[#]+dvscaled1[line1[#],line2[#]]*linep1[#])&;
  truev=parax/.FindRoot[vlinep[parax]==0,{parax,xlist[[-1]]}];
  (*Print[truev-xlist[[-1]]];*)
  scaledbounce=PathBounce2[vscaled,{dvscaled1,dvscaled2},{analyticalz1,analyticalz2},{{dz1d1,dz1d2},{dz2d1,dz2d2}},{line1,line2},{linep1,linep2},truev,0,bigv,dimension,timestofind,accuracy,StepScale];
  dr=scaledbounce[[1]];
  fieldinfo=scaledbounce[[2]];
  rilist=Table[
    Which[
      xlist[[i]]>=fieldinfo[[1,2]],0,
      xlist[[i]]<=fieldinfo[[-1,2]],0,
      True,FirstPosition[fieldinfo,{_,parameterx_,_}/;parameterx<=xlist[[i]]][[1]]
    ],
  {i,nofpoints}];
  
  dxdrlist=Table[
    If[
      rilist[[i]]==0,0,
      (fieldinfo[[rilist[[i]],3]]+fieldinfo[[rilist[[i]]-1,3]])/2
    ],
  {i,nofpoints}];
  dxdr2list=dxdrlist*dxdrlist;
  
  truebounce=Interpolation[Transpose[{Transpose[fieldinfo][[1]]/scalealpha,Transpose[fieldinfo][[2]]}]];
  trueaction=CalculateTheAction2[dimension,dr,fieldinfo,vscaled,{analyticalz1,analyticalz2},{line1,line2},{linep1,linep2}]/(scalealpha^(dimension-2));
  
  
  
  (*Print[ListPlot[(2 \[Pi]^(dimension/2))/Gamma[dimension/2]*dr*Transpose[Table[fieldinfo[[i,1]]^(dimension-1)({(1/(2analyticalz1[line1[fieldinfo[[i,2]]],line2[fieldinfo[[i,2]]]])linep1[fieldinfo[[i,2]]]^2+1/(2analyticalz2[line1[fieldinfo[[i,2]]],line2[fieldinfo[[i,2]]]])linep2[fieldinfo[[i,2]]]^2)fieldinfo[[i,3]]^2,vscaled[line1[fieldinfo[[i,2]]],line2[fieldinfo[[i,2]]]]}),{i,Length[fieldinfo]}]],PlotRange\[Rule]All]];
  Print[ListPlot[Table[vscaled[line1[fieldinfo[[i,2]]],line2[fieldinfo[[i,2]]]],{i,Length[fieldinfo]}],PlotRange\[Rule]All]];
  Print[ListPlot[Transpose[Table[{line1[fieldinfo[[i,2]]],line2[fieldinfo[[i,2]]]},{i,Length[fieldinfo]}]],PlotRange\[Rule]All]];
  Print[trueaction];*)
  If[scaledbounce[[3]]===True,
    AppendTo[pathlist,pointslist];
    AppendTo[actionlist,trueaction];,
    AppendTo[pathlist,{{,}}];
    AppendTo[actionlist,Null];
    Print["skip",qqqq];
  ];
  AppendTo[deformtimes,qqqq];
  (*Print[Plot[truebounce[x],{x,0,fieldinfo[[-1,1]]/scalealpha}]];*)
  (*Print[truebounce];*)
  
  ];
  minVal=Min[deformtimes];
  maxVal=Max[deformtimes];
  colorFunction=ColorData[{"Rainbow",{minVal,maxVal}}];
  Print[plot1=ListLinePlot[actionlist,ColorFunction->Function[{x,y},colorFunction[x]],ColorFunctionScaling->False,PlotLegends->BarLegend[{colorFunction,{minVal,maxVal}},LegendLabel->"deform times"],PlotLabel->"actions",ImageSize->400]];
  colors=colorFunction/@deformtimes;
  plot2=ListLinePlot[pathlist,PlotStyle->colors,PlotLegends->BarLegend[{colors,{minVal,maxVal}},LegendLabel->"deform times"],PlotLabel->"paths",ImageSize->400];
  plot3=ListLinePlot[lastpoints,PlotStyle\[Rule]Black,PlotLegends->Placed[LineLegend[{Black},{"path solved"}],{0.8,0.2}]];
  Print[Show[plot2,plot3]];*)
  
  
  
  (*end extra test*)
  
  scalealpha=Sqrt[bigv];
  (*truebounce=Interpolation[Transpose[{Transpose[fieldinfo][[1]]/scalealpha,Transpose[fieldinfo][[2]]}]];
  mulitfieldbounce={line1[truebounce[#]],line2[truebounce[#]]}&;*)
  mulitfieldbounce=Interpolation[Transpose[{Transpose[fieldinfo][[1]]/scalealpha,Transpose[{line1/@Transpose[fieldinfo][[2]],line2/@Transpose[fieldinfo][[2]]}]}]];
  (*Print[mulitfieldbounce];*)
  trueaction=CalculateTheAction2[dimension,dr,fieldinfo,vscaled,{analyticalz1,analyticalz2},{line1,line2},{linep1,linep2}]/(scalealpha^(dimension-2));
  trueactionz1=CalculateTheActionz1[dimension,dr,fieldinfo,vscaled,{analyticalz1,analyticalz2},{line1,line2},{linep1,linep2}]/(scalealpha^(dimension-2));

  (*result={mulitfieldbounce,trueaction,scaledbounce[[3]]};*)
  result={mulitfieldbounce,trueaction,trueactionz1,scaledbounce[[3]]};
  WriteString["stdout","\[EmptyCircle]"];
  Return[result]
];


PathBounce2[vscaled_,{dvscaled1_,dvscaled2_},{z1_,z2_},{{dz1d1_,dz1d2_},{dz2d1_,dz2d2_}},{line1_,line2_},{linep1_,linep2_},truevacuum_,falsevacuum_,bigv_,dimension_,times_,accuracy_,StepScale_]:=Module[
{mmm,sign,fieldaccuracy,derivativeaccuracy,basicdr,dr,termr,termz1,termz2,termv,thexmax,weib,weim,rho,crho,parameterx0,parameterx,cparameterx,dxdr,cdxdr,ddxdr,cddxdr,fieldinfo,bounce2z,theaction,scalealpha,convergence},
  
  (*Veff(truevacuum) must < Veff(falsevacuum).*)
  If[vscaled[line1[truevacuum],line2[truevacuum]]>=vscaled[line1[falsevacuum],line2[falsevacuum]],
    Print["Wrong potential or VEV."];
    Return[{0,False}]
  ];
  
  sign=Sign[truevacuum-falsevacuum];
  fieldaccuracy=Abs[truevacuum-falsevacuum]*accuracy/2;
  derivativeaccuracy = accuracy/(5 Sqrt[dimension]);
  basicdr= Sqrt[Abs[dimension (truevacuum-falsevacuum)^2/0.04]]/2000;
  
  
  dr=basicdr*StepScale;
  (*(D-1)/rdx/dr*)
  termr[rhoo_,dtdr_]:=Piecewise[{{0,rhoo==0},{N[(dimension-1)/rhoo*dtdr],rhoo>0}}];
  (*-\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\(
\*FractionBox[\(1\), 
SubscriptBox[\(Z\), \(i\)]]
\*FractionBox[
SubscriptBox[\(dZ\), \(i\)], \(dx\)]
\*SuperscriptBox[\((
\*FractionBox[
SubscriptBox[\(d\[Phi]\), \(i\)], \(dx\)])\), \(2\)]
\*SuperscriptBox[\((
\*FractionBox[\(dx\), \(dr\)])\), \(2\)]\)\)*)
  termz1[parametert_,dtdr_]=-(1/z1[line1[parametert],line2[parametert]] (dz1d1[line1[parametert],line2[parametert]]*linep1[parametert]+dz1d2[line1[parametert],line2[parametert]]*linep2[parametert])(linep1[parametert])^2+1/z2[line1[parametert],line2[parametert]] (dz2d1[line1[parametert],line2[parametert]]*linep1[parametert]+dz2d2[line1[parametert],line2[parametert]]*linep2[parametert])(linep2[parametert])^2)*dtdr^2;
  (*1/2\!\(
\*SubscriptBox[\(\[Sum]\), \(i, j\)]\(
\*FractionBox[
SubscriptBox[\(Z\), \(i\)], 
SuperscriptBox[
SubscriptBox[\(Z\), \(j\)], \(2\)]]
\*FractionBox[
SubscriptBox[\(dZ\), \(j\)], 
SubscriptBox[\(d\[Phi]\), \(i\)]]
\*FractionBox[
SubscriptBox[\(d\[Phi]\), \(i\)], \(dx\)]
\*SuperscriptBox[\((
\*FractionBox[
SubscriptBox[\(d\[Phi]\), \(j\)], \(dx\)])\), \(2\)]
\*SuperscriptBox[\((
\*FractionBox[\(dx\), \(dr\)])\), \(2\)]\)\)*)
  termz2[parametert_,dtdr_]=1/2 (1/z1[line1[parametert],line2[parametert]] dz1d1[line1[parametert],line2[parametert]]*linep1[parametert]^3+z1[line1[parametert],line2[parametert]]/z2[line1[parametert],line2[parametert]]^2 dz2d1[line1[parametert],line2[parametert]]*linep1[parametert]*linep2[parametert]^2+z2[line1[parametert],line2[parametert]]/z1[line1[parametert],line2[parametert]]^2 dz1d2[line1[parametert],line2[parametert]]*linep2[parametert]*linep1[parametert]^2+1/z2[line1[parametert],line2[parametert]] dz2d2[line1[parametert],line2[parametert]]*linep2[parametert]^3)*dtdr^2;
  (*-\!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\(
\*SubscriptBox[\(Z\), \(i\)]
\*FractionBox[
SubscriptBox[\(d\[Phi]\), \(i\)], \(dx\)]
\*FractionBox[\(\[PartialD]Veff\), \(\[PartialD]
\*SubscriptBox[\(\[Phi]\), \(i\)]\)]\)\)*)
  termv[parametert_]=-(z1[line1[parametert],line2[parametert]]*linep1[parametert]*dvscaled1[line1[parametert],line2[parametert]]+z2[line1[parametert],line2[parametert]]*linep2[parametert]*dvscaled2[line1[parametert],line2[parametert]]);
  
  
  thexmax=parametert/.FindRoot[termv[parametert]==0,{parametert,truevacuum},AccuracyGoal->20,PrecisionGoal->20];
  parameterx=(thexmax+falsevacuum)/2;
  dxdr=2derivativeaccuracy;
  weib=1;
  weim=1;
  
  (*Print[truevacuum];
  Print[thexmax];*)
  While[(Abs[parameterx-falsevacuum]>=fieldaccuracy||Abs[dxdr]>=derivativeaccuracy)&&Log[2,weib+weim]<=times,
  
  parameterx0=(weib*(thexmax+falsevacuum)/2+weim*thexmax)/(weib+weim);
  parameterx=parameterx0;
  rho=0;
  dxdr=0;
  ddxdr=-termv[parameterx];
  (*Print[ddxdr];*)
  (*fieldinfo1={{{rho,z1[line1[parameterx],line2[parameterx]]},{rho,z2[line1[parameterx],line2[parameterx]]}}};*)
  mmm=0;
  
  (*Print[{dxdr*sign<=0,parameterx-falsevacuum>=0,parameterx<thexmax,0.95*parameterx0<parameterx<thexmax,parameterx==thexmax,parameterx>thexmax}];*)
  (*Print[parameterx0];*)
    While[(dxdr*sign<=0&&parameterx-falsevacuum>=0&&parameterx<thexmax)||0.95*parameterx0<parameterx<thexmax,
      crho=rho;
      cparameterx=parameterx;
      cdxdr=dxdr;
      cddxdr=ddxdr;
      (*If[mmm<=5,Print["x",parameterx,"x'",dxdr,"x''",ddxdr,"tr",-termr[crho,cdxdr],"tz1",-termz1[cparameterx,cdxdr],"tz2",-termz2[cparameterx,cdxdr],"tv",-termv[cparameterx]]];*)
      rho=crho+dr;
      parameterx=cparameterx+cdxdr*dr;
      dxdr=cdxdr+cddxdr*dr;
      ddxdr=-termr[crho,cdxdr]-termz1[cparameterx,cdxdr]-termz2[cparameterx,cdxdr]-termv[cparameterx];
      (*AppendTo[fieldinfo1,{{rho,z1[line1[parameterx],line2[parameterx]]},{rho,z2[line1[parameterx],line2[parameterx]]}}];*)
      mmm++;
    ];
  (*Print[{dxdr*sign<=0,parameterx-falsevacuum>=0,parameterx<thexmax,0.95*parameterx0<parameterx<thexmax,parameterx==thexmax,parameterx>thexmax}];
  Print[ListPlot[Transpose[fieldinfo1],PlotRange\[Rule]All]];*)
  If[dxdr*sign>0,
    weib=2weib-1;weim=2weim+1;,
    weib=2weib+1;weim=2weim-1;
    ];
    WriteString["stdout","|",Log[2,weib+weim]];
  ];
  (*EmitSound[Play[Sin[440*2*Pi*t],{t,0,1}]];*)
  
  parameterx=parameterx0+0.;
  rho=0.;
  dxdr=0.;
  ddxdr=-termv[parameterx]+0.;
  fieldinfo=Table[{0.,(*parameterx0*)0.,0.},mmm];
  (*Print[parameterx0];*)
  (*Print[{dxdr*sign<=0,parameterx-falsevacuum>=0,parameterx<thexmax,0.95*parameterx0<parameterx<thexmax,parameterx==thexmax,parameterx>thexmax}];*)
  mmm=0;
  
  While[(dxdr*sign<=0&&parameterx-falsevacuum>=0&&parameterx<thexmax)||0.95*parameterx0<parameterx<thexmax,
    crho=rho;
    cparameterx=parameterx;
    cdxdr=dxdr;
    cddxdr=ddxdr;
      (*If[mmm<=5,Print["x",parameterx,"x'",dxdr,"x''",ddxdr,"tr",-termr[crho,cdxdr],"tz1",-termz1[cparameterx,cdxdr],"tz2",-termz2[cparameterx,cdxdr],"tv",-termv[cparameterx]]];*)
    rho=crho+dr;
    parameterx=cparameterx+cdxdr*dr;
    dxdr=cdxdr+cddxdr*dr;
    ddxdr=-termr[crho,cdxdr]-termz1[cparameterx,cdxdr]-termz2[cparameterx,cdxdr]-termv[cparameterx];
    
    mmm++;
    fieldinfo[[mmm]]={rho,parameterx,dxdr};
  ];
  
  (*Print[parameterx," ",dxdr," ",ddxdr," ",mmm];*)
  scalealpha=Sqrt[bigv];
  (*Print[ListPlot[Transpose[Transpose[fieldinfo][[1;;2]]]]];*)
  If[Log[2,weib+weim]>times,convergence=False,convergence=True];
  (*Print[convergence];*)
  Return[{dr,fieldinfo,convergence}]
];


CalculateTheAction2[dimension_,dr_,fieldinfos_,vscaled_,{z1_,z2_},{line1_,line2_},{linep1_,linep2_}]:=
(2 \[Pi]^(dimension/2))/Gamma[dimension/2]*dr*Sum[fieldinfos[[i,1]]^(dimension-1) ((1/(2z1[line1[fieldinfos[[i,2]]],line2[fieldinfos[[i,2]]]]) linep1[fieldinfos[[i,2]]]^2+1/(2z2[line1[fieldinfos[[i,2]]],line2[fieldinfos[[i,2]]]]) linep2[fieldinfos[[i,2]]]^2)fieldinfos[[i,3]]^2+vscaled[line1[fieldinfos[[i,2]]],line2[fieldinfos[[i,2]]]]),{i,Length[fieldinfos]}];


CalculateTheActionz1[dimension_,dr_,fieldinfos_,vscaled_,{z1_,z2_},{line1_,line2_},{linep1_,linep2_}]:=
(2 \[Pi]^(dimension/2))/Gamma[dimension/2]*dr*Sum[fieldinfos[[i,1]]^(dimension-1) ((1/(2z1[line1[fieldinfos[[i,2]]],line2[fieldinfos[[i,2]]]]) linep1[fieldinfos[[i,2]]]^2(*+1/(2z2[line1[fieldinfos[[i,2]]],line2[fieldinfos[[i,2]]]]) linep2[fieldinfos[[i,2]]]^2*))fieldinfos[[i,3]]^2+vscaled[line1[fieldinfos[[i,2]]],line2[fieldinfos[[i,2]]]]),{i,Length[fieldinfos]}];


force2[points_,nofpoints_,dtdr2_,{dvd1_,dvd2_},{z1_,z2_},{{dz1d1_,dz1d2_},{dz2d1_,dz2d2_}}]:=Module[
{dphis,dts,dphidt,adphidt,ddphis,ddphidt,eddphidt,addphidt,gv1,gv2,zs,gz,dzdt,gvmax,az,fn},

dphis=Delete[points,1]-Delete[points,-1];

dts=Sqrt[Total[Transpose[dphis*dphis]]];
dphidt=dphis/dts;

adphidt=(Append[dphidt,dphidt[[-1]]]+Prepend[dphidt,dphidt[[1]]])/2;

ddphis=Delete[adphidt,1]-Delete[adphidt,-1];
ddphidt=ddphis/dts;

addphidt=(Append[ddphidt,ddphidt[[-1]]]+Prepend[ddphidt,ddphidt[[1]]])/2;

gv1=MapThread[dvd1,Transpose[points]];
gv2=MapThread[dvd2,Transpose[points]];

(*gvmax=Max[Sqrt[gv1*gv1+gv2*gv2]];*)
gvmax=Max[Abs[{gv1,gv2}]];

fn=Table[
  {(addphidt[[i,1]]-((dz1d1@@points[[i]])*adphidt[[i,1]]+(dz1d2@@points[[i]])*adphidt[[i,2]])/z1@@points[[i]]*adphidt[[i,1]]+(((dz1d1@@points[[i]])*adphidt[[i,1]]+(dz1d2@@points[[i]])*adphidt[[i,2]])/z1@@points[[i]]*adphidt[[i,1]]^2+((dz2d1@@points[[i]])*adphidt[[i,1]]+(dz2d2@@points[[i]])*adphidt[[i,2]])/z2@@points[[i]]*adphidt[[i,2]]^2)*adphidt[[i,1]]+1/2 (((z1@@points[[i]])(dz1d1@@points[[i]]))/(z1@@points[[i]])^2 adphidt[[i,1]]^2+((z1@@points[[i]])(dz2d1@@points[[i]]))/(z2@@points[[i]])^2 adphidt[[i,2]]^2)-1/2 (((z1@@points[[i]])*(dz1d1@@points[[i]])*adphidt[[i,1]])/(z1@@points[[i]])^2*adphidt[[i,1]]^2+((z2@@points[[i]])*(dz1d2@@points[[i]])*adphidt[[i,2]])/(z1@@points[[i]])^2*adphidt[[i,1]]^2+((z1@@points[[i]])*(dz2d1@@points[[i]])*adphidt[[i,1]])/(z2@@points[[i]])^2*adphidt[[i,2]]^2+((z2@@points[[i]])*(dz2d2@@points[[i]])*adphidt[[i,2]])/(z2@@points[[i]])^2*adphidt[[i,2]]^2)adphidt[[i,1]])*dtdr2[[i]]-(z1@@points[[i]])*(dvd1@@points[[i]])+((z1@@points[[i]])*adphidt[[i,1]]*(dvd1@@points[[i]])+(z2@@points[[i]])*adphidt[[i,2]]*(dvd2@@points[[i]]))*adphidt[[i,1]],
  (addphidt[[i,2]]-((dz2d1@@points[[i]])*adphidt[[i,1]]+(dz2d2@@points[[i]])*adphidt[[i,2]])/z2@@points[[i]]*adphidt[[i,2]]+(((dz1d1@@points[[i]])*adphidt[[i,1]]+(dz1d2@@points[[i]])*adphidt[[i,2]])/z1@@points[[i]]*adphidt[[i,1]]^2+((dz2d1@@points[[i]])*adphidt[[i,1]]+(dz2d2@@points[[i]])*adphidt[[i,2]])/z2@@points[[i]]*adphidt[[i,2]]^2)*adphidt[[i,2]]+1/2 (((z2@@points[[i]])(dz1d2@@points[[i]]))/(z1@@points[[i]])^2 adphidt[[i,1]]^2+((z2@@points[[i]])(dz2d2@@points[[i]]))/(z2@@points[[i]])^2 adphidt[[i,2]]^2)-1/2 (((z1@@points[[i]])*(dz1d1@@points[[i]])*adphidt[[i,1]])/(z1@@points[[i]])^2*adphidt[[i,1]]^2+((z2@@points[[i]])*(dz1d2@@points[[i]])*adphidt[[i,2]])/(z1@@points[[i]])^2*adphidt[[i,1]]^2+((z1@@points[[i]])*(dz2d1@@points[[i]])*adphidt[[i,1]])/(z2@@points[[i]])^2*adphidt[[i,2]]^2+((z2@@points[[i]])*(dz2d2@@points[[i]])*adphidt[[i,2]])/(z2@@points[[i]])^2*adphidt[[i,2]]^2)adphidt[[i,2]])*dtdr2[[i]]-(z2@@points[[i]])*(dvd2@@points[[i]])+((z1@@points[[i]])*adphidt[[i,1]]*(dvd1@@points[[i]])+(z2@@points[[i]])*adphidt[[i,2]]*(dvd2@@points[[i]]))*adphidt[[i,2]]},
{i,nofpoints}];

Return[{fn,gvmax}]]


step2[points_,nofpoints_,dtdr2_,{dvscaled1_,dvscaled2_},{analyticalz1_,analyticalz2_},{{dz1d1_,dz1d2_},{dz2d1_,dz2d2_}},instepsize_]:=Module[{minstep,stepsize,f1,fmax,fratio,points2,f2,dfmax,nextpoints,lappoints,epsilon},
minstep=10^-6;
stepsize=instepsize;
f1=force2[points,nofpoints,dtdr2,{dvscaled1,dvscaled2},{analyticalz1,analyticalz2},{{dz1d1,dz1d2},{dz2d1,dz2d2}}];

fmax=Max[Abs[f1[[1]]]];
fratio=fmax/f1[[2]];
While[True,
  points2=points+f1[[1]]*stepsize*0.5;
  
  f2=force2[points2,nofpoints,dtdr2,{dvscaled1,dvscaled2},{analyticalz1,analyticalz2},{{dz1d1,dz1d2},{dz2d1,dz2d2}}];
  
  If[stepsize<=minstep,stepsize=minstep;Break[]];
    dfmax=Max[Abs[f2[[1]]-f1[[1]]]];
    
  If[dfmax<fmax*0.1,Break[]];
    stepsize=stepsize*0.5;
];

nextpoints=points2+f2[[1]]*stepsize*0.5;

(*Laplacian smoothing*)
epsilon=Min[5/nofpoints,0.05]*Exp[-(20*nowd+nowstep)/50];
lappoints=nextpoints+epsilon*(Prepend[Append[Delete[nextpoints,{{1},{2}}],nextpoints[[-1]]],nextpoints[[1]]]+Prepend[Append[Delete[nextpoints,{{-1},{-2}}],nextpoints[[-1]]],nextpoints[[1]]]-2nextpoints);

(*Print[ListPlot[{f1[[1]]*stepsize*0.5+f2[[1]]*stepsize*0.5,nextpoints-points}]];
Print[ListPlot[nextpoints-points]];*)
Return[{lappoints,stepsize,fratio}]]


defom2[points_,nofpoints_,dtdr2_,{dvscaled1_,dvscaled2_},{analyticalz1_,analyticalz2_},{{dz1d1_,dz1d2_},{dz2d1_,dz2d2_}}]:=Module[{dov,gv1,gv2,gvmax,stepsize,pts,j,k,nextis,fratio},

dov=Norm[Last[points]-First[points]];

gv1=MapThread[dvscaled1,Transpose[points]];
gv2=MapThread[dvscaled2,Transpose[points]];

gvmax=Max[Sqrt[gv1*gv1+gv2*gv2]];
(*gvmax=Max[Abs[{gv1,gv2}]];*)

stepsize=0.1*dov/gvmax/nofpoints;
pts=points;
j=0;
k=0;
fratio=1;
nowstep=0;
While[j<=500&&fratio>=0.02,j++;
nextis=step2[pts,nofpoints,dtdr2,{dvscaled1,dvscaled2},{analyticalz1,analyticalz2},{{dz1d1,dz1d2},{dz2d1,dz2d2}},stepsize];
nowstep++;
pts=nextis[[1]];
(*Print[pts];*)
stepsize=1.5*nextis[[2]];
fratio=nextis[[3]];

];
Print[ListPlot[pts]];
If[j<=1,k++];

WriteString["stdout","\:2219"];
Print[j];
Return[{pts,k}]]


(* ::Subsubsection:: *)
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
(*dphis=Table[points[[i+1]]-points[[i]],{i,n-1}];*)
dphis=Delete[points,1]-Delete[points,-1];
(*dts=Table[Norm[dphis[[i]]],{i,n-1}];*)
dts=Sqrt[Total[Transpose[dphis*dphis]]];
dphidt=dphis/dts;
(*edphidt=Prepend[dphidt,dphidt[[1]]];
AppendTo[dphidt,dphidt[[n-1]]];
adphidt=(dphidt+edphidt)/2;*)
adphidt=(Append[dphidt,dphidt[[-1]]]+Prepend[dphidt,dphidt[[1]]])/2;
(*caculate the average of left and right derivative*)
(*ddphis=Table[adphidt[[i+1]]-adphidt[[i]],{i,n-1}];*)
ddphis=Delete[adphidt,1]-Delete[adphidt,-1];
ddphidt=ddphis/dts;
(*eddphidt=Prepend[ddphidt,ddphidt[[1]]];
AppendTo[ddphidt,ddphidt[[n-1]]];
addphidt=(ddphidt+eddphidt)/2;*)
addphidt=(Append[ddphidt,ddphidt[[-1]]]+Prepend[ddphidt,ddphidt[[1]]])/2;
(*caculate the average of left and right second derivative*)

gv=MapThread[dvdf,Transpose[points]];
zs=MapThread[z,Transpose[points]];

gz=MapThread[dzdf,Transpose[points]];
gvmax=Max[Abs[gv]];

fg=dtdr2*gz/zs/2-zs*gv;
(*fng=fg-Table[fg[[i]].adphidt[[i]]*adphidt[[i]],{i,n}];*)(*not fastest*)
(*fng=fg-Diagonal[fg.Transpose[adphidt]]*adphidt;*)(*fastest in this single step, but not in tatal*)
fng=fg-ArrayReduce[Total,fg*adphidt,2]*adphidt;(*fastest in total*)
fn=addphidt*dtdr2+fng;

Return[{fn,gvmax}]]


(* ::Subsubsection:: *)
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

fg=dtdr2*gz/zs/2-zs*gv;

fng=fg-ArrayReduce[Total,fg*adphidt,2]*adphidt;
fn=addphidt*dtdr2+fng;

Return[{fn,gvmax}]]


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
(*nstep*)


nstep[points_,dtdr2_,ddtdr_,instepsize_,v_,z_,dov_]:=Module[{minstep,stepsize,f1,fmax,fratio,points2,f2,dfmax,nextpoints},
minstep=10^-6;
stepsize=instepsize;
f1=nforce[points,dtdr2,ddtdr,v,z,dov];
fmax=Max[Abs[f1[[1]]]];
fratio=fmax/f1[[2]];
While[True,
points2=points+f1[[1]]*stepsize*0.5;
f2=nforce[points2,dtdr2,ddtdr,v,z,dov][[1]];
If[stepsize<=minstep,stepsize=minstep;Break[]];
dfmax=Max[Abs[f2-f1[[1]]]];
If[dfmax<fmax*0.1,Break[]];
stepsize=stepsize*0.5;];
nextpoints=points2+f2*stepsize*0.5;
Return[{nextpoints,stepsize,fratio}]]


(* ::Subsubsection:: *)
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
(*Print[pts];*)
stepsize=1.5*nextis[[2]];
fratio=nextis[[3]];
];

If[j<=1,k++];

WriteString["stdout","\:2219"];

Return[{pts,k}]]


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
(*AnalyticalPotentialAnalyticalRnormalization*)


AnalyticalPotentialAnalyticalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_, 
                                            barrierbetweenvacuums_, dimension_, 
                                            times_, accuracy_, StepScale_]:=
Module[{OriginalFunctionOfPotential,OriginalFunctionOfRenormalization,
  domainl,domainr,barrier,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,\[Phi]max,bigv,\[Alpha],NormalizedPotential,supercooling,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  
  OriginalFunctionOfPotential=Function[thefield,potential/.{fieldname->thefield}];
  OriginalFunctionOfRenormalization=Function[thefield,renormalization/.{fieldname->thefield}];
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<=thefield<=domainr},thefield]],barrierbetweenvacuums];
  If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],\[Phi]max=10(barrier-FalseVacuum)+FalseVacuum;supercooling=True,\[Phi]max=TrueVacuum;supercooling=False];
  PotentialAtFalseVacuum=OriginalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=OriginalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  bigv=Abs[FunctionOfRelativePotential@\[Phi]max]/0.04;
  (*Rescale the potential. *)
  \[Alpha]=Sqrt[bigv];
  
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  
  (*Print[t0101[[1]]];*)
  
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,
                                              barrier,dimension,times, accuracy,StepScale,supercooling];
  
  
  FieldAndDerivative=(*InitialFieldAndBigR[[1]];*)GetFieldPoints[NormalizedPotential, OriginalFunctionOfRenormalization, \[Phi]max, FalseVacuum, 
                                              barrier, dimension, times, accuracy,StepScale,
                                              {InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  
  
  FieldPointsOfR=Transpose[{Transpose[FieldAndDerivative][[1]]/\[Alpha],Transpose[FieldAndDerivative][[2]]}];
  
  TheAction=ParallelSum[Extract[FieldAndDerivative,{q,5}]*((Extract[FieldAndDerivative,{q,1}])^(dimension-1))*(1/(2 OriginalFunctionOfRenormalization[Extract[FieldAndDerivative,{q,2}]])*(Extract[FieldAndDerivative,{q,3}])^2+NormalizedPotential[Extract[FieldAndDerivative,{q,2}]]),{q,InitialFieldAndBigR[[2]]},DistributedContexts-> {"VacuumTunneling`Private`"}]*(2\[Pi]^(dimension/2))/Gamma[dimension/2]/(\[Alpha]^(dimension-2));
  (*TheAction=CalculateTheAction[NormalizedPotential, OriginalFunctionOfRenormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));*)
  
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection:: *)
(*NumbericalPotentialAnalyticalRnormalization*)


NumbericalPotentialAnalyticalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_, 
                                            barrierbetweenvacuums_, dimension_, 
                                            times_, accuracy_, StepScale_]:=
Module[{OriginalExpressionOfPotential,OriginalFunctionOfPotential,AnalyticalFunctionOfPotential,OriginalFunctionOfRenormalization,
  domainl,domainr,barrier,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,\[Phi]max,bigv,\[Alpha],NormalizedPotential,supercooling,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  OriginalExpressionOfPotential[fieldname_?NumericQ]=potential;
  OriginalFunctionOfPotential=Function[thefield,OriginalExpressionOfPotential[thefield]];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<thefield<domainr},thefield]],barrierbetweenvacuums];
  If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],\[Phi]max=10(barrier-FalseVacuum)+FalseVacuum;supercooling=True,\[Phi]max=TrueVacuum;supercooling=False];
  AnalyticalFunctionOfPotential=InterpolationInArea[OriginalFunctionOfPotential,\[Phi]max,FalseVacuum,barrier];
  OriginalFunctionOfRenormalization=Function[thefield,renormalization/.{fieldname->thefield}];
  PotentialAtFalseVacuum=AnalyticalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=AnalyticalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  bigv=Abs[FunctionOfRelativePotential@\[Phi]max]/0.04;
  \[Alpha]=Sqrt[bigv];
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,
                                              barrier,dimension,times, accuracy,StepScale,supercooling];
  FieldAndDerivative=GetFieldPoints[NormalizedPotential, OriginalFunctionOfRenormalization, \[Phi]max, FalseVacuum, 
                                    barrier, dimension, times, accuracy, StepScale,
                                    {InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  FieldPointsOfR=Transpose[{Transpose[FieldAndDerivative][[1]]/\[Alpha],Transpose[FieldAndDerivative][[2]]}];
  TheAction=CalculateTheAction[NormalizedPotential, OriginalFunctionOfRenormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection:: *)
(*AnalyticalPotentialNumbericalRnormalization*)


AnalyticalPotentialNumbericalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_, 
                                            barrierbetweenvacuums_, dimension_, 
                                            times_, accuracy_, StepScale_]:=
Module[{OriginalFunctionOfPotential,OriginalFunctionOfRenormalization,OriginalExpressionOfRnormalization,AnalyticalFunctionOfRnormalization,
  domainl,domainr,barrier,\[Phi]max,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,bigv,\[Alpha],NormalizedPotential,supercooling,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  OriginalFunctionOfPotential=Function[thefield,potential/.{fieldname->thefield}];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<=thefield<=domainr},thefield]],barrierbetweenvacuums];
  If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],\[Phi]max=10(barrier-FalseVacuum)+FalseVacuum;supercooling=True,\[Phi]max=TrueVacuum;supercooling=False];
  FunctionOfRelativePotential=OriginalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  PotentialAtFalseVacuum=OriginalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=OriginalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  OriginalExpressionOfRnormalization[fieldname_?NumericQ]=renormalization;
  OriginalFunctionOfRenormalization=Function[thefield,OriginalExpressionOfRnormalization[thefield]];
  AnalyticalFunctionOfRnormalization=InterpolationInArea[OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,barrier];
  bigv=Abs[FunctionOfRelativePotential@\[Phi]max]/0.04;
  \[Alpha]=Sqrt[bigv];
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,AnalyticalFunctionOfRnormalization,\[Phi]max,FalseVacuum,
                                              barrier,dimension,times,accuracy,StepScale,supercooling];
  FieldAndDerivative=GetFieldPoints[NormalizedPotential, AnalyticalFunctionOfRnormalization,\[Phi]max,FalseVacuum,
                                    barrier,dimension,times,accuracy,StepScale,
                                    {InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  FieldPointsOfR=Transpose[{Transpose[FieldAndDerivative][[1]]/\[Alpha],Transpose[FieldAndDerivative][[2]]}];
  TheAction=CalculateTheAction[NormalizedPotential, AnalyticalFunctionOfRnormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));
  (*bounce=Interpolation[FieldPointsOfR];*)
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection:: *)
(*NumbericalPotentialNumbericalRnormalization*)


NumbericalPotentialNumbericalRnormalization[potential_, renormalization_, fieldname_, TrueVacuum_, FalseVacuum_,
                                            barrierbetweenvacuums_, dimension_, times_, accuracy_, StepScale_]:=
Module[{OriginalExpressionOfpotential,OriginalFunctionOfPotential,AnalyticalFunctionOfPotential,OriginalExpressionOfRnormalization,OriginalFunctionOfRenormalization,AnalyticalFunctionOfRnormalization,
  domainl,domainr,barrier,\[Phi]max,thefield,PotentialAtFalseVacuum,FunctionOfRelativePotential,bigv,\[Alpha],NormalizedPotential,supercooling,
  InitialFieldAndBigR,FieldAndDerivative,FieldPointsOfR,TheAction,result},
  domainl=Min[TrueVacuum, FalseVacuum];
  domainr=Max[TrueVacuum, FalseVacuum];
  OriginalExpressionOfpotential[fieldname_?NumericQ]=potential;
  OriginalFunctionOfPotential=Function[thefield,OriginalExpressionOfpotential[thefield]];
  barrier=If[barrierbetweenvacuums===Null,thefield/.Last[FindMaximum[{OriginalFunctionOfPotential[thefield],domainl<=thefield<=domainr},thefield]],barrierbetweenvacuums];
  If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],\[Phi]max=10(barrier-FalseVacuum)+FalseVacuum;supercooling=True,\[Phi]max=TrueVacuum;supercooling=False];
  AnalyticalFunctionOfPotential=InterpolationInArea[OriginalFunctionOfPotential,\[Phi]max,FalseVacuum,barrier];
  
  OriginalExpressionOfRnormalization[fieldname_?NumericQ]=renormalization;
  OriginalFunctionOfRenormalization=Function[thefield,OriginalExpressionOfRnormalization[thefield]];
  AnalyticalFunctionOfRnormalization=InterpolationInArea[OriginalFunctionOfRenormalization,\[Phi]max,FalseVacuum,barrier];
  
  PotentialAtFalseVacuum=AnalyticalFunctionOfPotential[FalseVacuum];
  FunctionOfRelativePotential=AnalyticalFunctionOfPotential[#]-PotentialAtFalseVacuum &;
  bigv=Abs[FunctionOfRelativePotential@\[Phi]max]/0.04;
  \[Alpha]=Sqrt[bigv];
  NormalizedPotential=FunctionOfRelativePotential[#]/bigv &;
  
  InitialFieldAndBigR=FindInitialFieldAndBigR[NormalizedPotential,AnalyticalFunctionOfRnormalization,\[Phi]max,FalseVacuum,
                                              barrier,dimension,times, accuracy, StepScale, supercooling];
  FieldAndDerivative=GetFieldPoints[NormalizedPotential, AnalyticalFunctionOfRnormalization, \[Phi]max, FalseVacuum, barrier, 
                                    dimension, times, accuracy, StepScale,
                                    {InitialFieldAndBigR[[1]],InitialFieldAndBigR[[2]]}];
  FieldPointsOfR=Transpose[{Transpose[FieldAndDerivative][[1]]/\[Alpha],Transpose[FieldAndDerivative][[2]]}];
  
  TheAction=CalculateTheAction[NormalizedPotential, AnalyticalFunctionOfRnormalization,InitialFieldAndBigR[[2]],dimension,FieldAndDerivative]/(\[Alpha]^(dimension-2));
  (*bounce=Interpolation[FieldPointsOfR];*)
  result={Interpolation[FieldPointsOfR],Last[FieldAndDerivative][[1]]/\[Alpha],TheAction,InitialFieldAndBigR[[3]]};
  If[Is1DTunneling==True,WriteString["stdout","\[EmptyCircle]"]];
  Return[result]
  ];


(* ::Subsubsection:: *)
(*InterpolationInArea*)


InterpolationInArea[NumbericalFunction_,TrueVacuum_,FalseVacuum_,barrier_]:=
  Module[{n,dx,Functionpts,\[Phi]max,ifun1},
  n=1000;
  \[Phi]max=If[Abs[TrueVacuum-barrier]>10Abs[barrier-FalseVacuum],10(barrier-FalseVacuum)+FalseVacuum,TrueVacuum];
  dx=(\[Phi]max-FalseVacuum)/n;
  Functionpts=ParallelTable[Off[NIntegrate::inumr];Quiet[{(i-1)*dx+FalseVacuum,NumbericalFunction[(i-1)*dx+FalseVacuum]}],{i,(n+1)},DistributedContexts-> {"VacuumTunneling`Private`"}];
  
  ifun1=Interpolation[Cases[Functionpts,{_,_?NumericQ}]];
  
  Return[ifun1]
  ];


(* ::Subsubsection:: *)
(*FindInitialFieldAndBigR*)


FindInitialFieldAndBigR[AnalyticalPotential_, AnalyticalRenormalization_, TrueVacuum_, FalseVacuum_, barrier_,
                        dimension_, times_, accuracy_,StepScale_,supercooling_]:=
  Module[{(*$MinPrecision = 50, $MaxPrecision = 50,*)afield,vp,zp,fvp,fzp,\[Phi]max,sign,precision,absolutefieldaccuracy,absolutederivativeaccuracy,basicdr,dr,weib,weim,\[Phi]0,\[Phi],d\[Phi],r,c\[Phi],cd\[Phi],dd\[Phi],j,convergence},
  
  vp[afield_]=AnalyticalPotential'[afield];
  zp[afield_]=AnalyticalRenormalization'[afield];
  fvp=Function[afield,vp[afield]];
  fzp=Function[afield,zp[afield]];
  
  \[Phi]max=TrueVacuum;
  sign=Sign[\[Phi]max-FalseVacuum];
  precision=Ceiling[Log[10,(2^times)/Abs[TrueVacuum-FalseVacuum]]];
  
  If[supercooling==False,
     \[Phi]max=thefield/.FindRoot[fvp[thefield]==0,{thefield,0.99TrueVacuum},PrecisionGoal->20,AccuracyGoal->20]];
  (*to prevent wrong true vacuum, which makes the field rolling far away from false vacuum.*)
  
  
  absolutefieldaccuracy=Abs[barrier-FalseVacuum]*accuracy;
  absolutederivativeaccuracy = accuracy/(5 Sqrt[dimension]);
  basicdr= Sqrt[Abs[dimension (\[Phi]max-FalseVacuum)^2/0.04]]/2000;
  
  d\[Phi]=1;
  weib=1;weim=1;
  \[Phi]=(weib barrier+weim \[Phi]max)/(weib+weim);
  dr=basicdr * StepScale;

  While[(Abs[\[Phi]-FalseVacuum]>=absolutefieldaccuracy||Abs[d\[Phi]]>=absolutederivativeaccuracy)&&Log[2,weib+weim]<=times,
  
  \[Phi]0=(weib barrier+weim \[Phi]max)/(weib+weim);
  r=0.;
  \[Phi]=\[Phi]0;
  d\[Phi]=0.;
  j=0;
  c\[Phi]=\[Phi];
  cd\[Phi]=d\[Phi];
  dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*fvp[c\[Phi]];
  
  While[(d\[Phi]*sign<=0||(\[Phi]-barrier)*sign>=0)&&(\[Phi]-FalseVacuum)*sign>=0,
  
  d\[Phi]=cd\[Phi]+dd\[Phi]*dr;
  \[Phi]=c\[Phi]+cd\[Phi] dr;
  r=r+dr;
  c\[Phi]=\[Phi];
  cd\[Phi]=d\[Phi];
  dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*fvp[c\[Phi]]+1/2 fzp[c\[Phi]]/AnalyticalRenormalization[c\[Phi]] cd\[Phi]^2-(dimension-1)/r cd\[Phi];
  j++;
  
  ];
  
  weib=If[d\[Phi] sign<0,2weib+1,2weib-1];
  weim=If[d\[Phi] sign<0,2weim-1,2weim+1];
  If[Is1DTunneling==True,WriteString["stdout","\:2219"]];
  ];
  
  If[Log[2,weib+weim]>times,convergence=False,convergence=True];
  
  Return[{\[Phi]0,j,convergence}]
  ];


(* ::Subsubsection:: *)
(*Finding the Initial Shooting Point/Tunneling Point*)


FindInitialPoint[Vi_,\[Phi]tv_,\[Phi]fv_,\[Phi]bar_,dr0_,tol_,times_]:=Module[
{sign,dr,d\[Phi],weib,weit,\[Phi],\[Phi]0,\[Mu],\[Delta],j,c\[Phi],c\[Mu],c\[Delta],cd\[Phi],d\[Mu],d\[Delta],dd\[Phi]},
sign=Sign[\[Phi]tv-\[Phi]fv];
dr=dr0;
d\[Phi]=1;
weib=1;weit=1;(*initial value of weight functions*)
\[Phi]=(weib \[Phi]bar+weit \[Phi]tv)/(weib+weit);
(*The first Loop to find the tunneling point, Cycle the intial value of \[Phi]0 until the result reach the accury goal or the maxiam step.*)
While[(Abs[\[Phi]-\[Phi]fv]>=tol||Abs[d\[Phi]]>=tol)&&Log[2,weib+weit]<=times,
(*The initial Value from the Bounce equation*)
\[Phi]0=(weib \[Phi]bar+weit \[Phi]tv)/(weib+weit);
\[Phi]=\[Phi]0;
d\[Phi]=(GravityTunnelingrh Vi'[\[Phi]0] )/(1-8\[Pi] GravityTunnelingG (GravityTunnelingrh^2) Vi[\[Phi]0] );
\[Mu]=0;
\[Delta]=0;
j=0;
c\[Phi]=\[Phi];
c\[Mu]=\[Mu];
c\[Delta]=\[Delta];
cd\[Phi]=d\[Phi];
(*The first evolution which getting from the EoM*)
d\[Mu]=4\[Pi] GravityTunnelingrh^2 Vi[c\[Phi]];
d\[Delta]=4\[Pi] GravityTunnelingG GravityTunnelingrh cd\[Phi]^2;
dd\[Phi]=0;
\[Mu]=c\[Mu]+d\[Mu] dr;
\[Delta]=c\[Delta]+d\[Delta] dr;
d\[Phi]=cd\[Phi]+dd\[Phi] dr;
\[Phi]=c\[Phi]+cd\[Phi] dr;
(*The second Loop to find the whether the \[Phi] evolution from \[Phi]0 will end at the false vacuum \[Phi]fv, Cycle until the derivative of \[Phi] changes sign,but at least it will pass the maximum value,and at the farthest to the false vacuum.*)
While[(sign*d\[Phi]<=0|| sign*(\[Phi]-\[Phi]bar)>=0) && sign*(\[Phi]-\[Phi]fv)>=0,
j++;
c\[Mu]=\[Mu];
c\[Delta]=\[Delta];
c\[Phi]=\[Phi];
cd\[Phi]=d\[Phi];
d\[Mu]=4\[Pi] (GravityTunnelingrh+j dr)^2 (1/2 (1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+j dr))cd\[Phi]^2+Vi[c\[Phi]]);
d\[Delta]=4\[Pi] GravityTunnelingG(GravityTunnelingrh+j dr)cd\[Phi]^2;
dd\[Phi]=(Vi'[c\[Phi]]-((2GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu])-2GravityTunnelingG(GravityTunnelingrh+j dr)d\[Mu])/(GravityTunnelingrh+j dr)^2 cd\[Phi]+2/(GravityTunnelingrh+j dr) (1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+j dr))cd\[Phi]+d\[Delta](1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+j dr))cd\[Phi]))/(1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+j dr));
\[Mu]=c\[Mu]+d\[Mu] dr;
\[Delta]=c\[Delta]+d\[Delta] dr;
d\[Phi]=cd\[Phi]+dd\[Phi]*dr;
\[Phi]=c\[Phi]+cd\[Phi] dr;];
WriteString["stdout","\:2219"];
(*If the result is overshooting, then take \[Phi]0 close to the false vacuum, undershooting then make \[Phi]0 close to ture vacuum*)
weib=If[d\[Phi]<0,2weib+1,2weib-1];weit=If[d\[Phi]<0,2weit-1,2weit+1]];
If[Log[2,weib+weit]>100,Print["The loop cannot converge within the maximum number of steps. 
Please increase the number of steps, 
decrease the step size and lower the convergence condition."],Print["The tunneling point has been found"]];
Return[{\[Phi]0,j}];
]


(* ::Subsubsection:: *)
(*GetFieldPoints*)


GetFieldPoints[AnalyticalPotential_, AnalyticalRenormalization_, TrueVacuum_, FalseVacuum_, barrier_, 
               dimension_, times_, accuracy_, StepScale_, 
               {TheInitialField_,NumberOfTheBigR_}]:=
  Module[{afield,vp,zp,fvp,fzp,dr0,absolutefieldaccuracy,absolutederivativeaccuracy,r,dr,\[Phi],d\[Phi],dd\[Phi],c\[Phi],cd\[Phi],fieldinfo,datafieldinfo,k},
  
  vp[afield_]=AnalyticalPotential'[afield];
  zp[afield_]=AnalyticalRenormalization'[afield];
  fvp=Function[afield,vp[afield]];
  fzp=Function[afield,zp[afield]];
  fieldinfo=Array[datafieldinfo,{NumberOfTheBigR,5}];
  dr0=Sqrt[Abs[dimension (TrueVacuum-FalseVacuum)^2/0.04]]/2000;
  
  absolutefieldaccuracy=Abs[barrier-FalseVacuum]*accuracy;
  
  absolutederivativeaccuracy = accuracy/(5 Sqrt[dimension]);
  
  r=0.;
  \[Phi]=TheInitialField;
  dr=dr0 * StepScale;
  d\[Phi]=0.;
  k=0;
  c\[Phi]=\[Phi];
  cd\[Phi]=d\[Phi];
  dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*fvp[c\[Phi]];
  
  While[k<NumberOfTheBigR,
    datafieldinfo[k+1,1]=r;
    datafieldinfo[k+1,2]=c\[Phi];
    datafieldinfo[k+1,3]=cd\[Phi];
    datafieldinfo[k+1,4]=dd\[Phi]; 
    datafieldinfo[k+1,5]=dr;
    (* Record current r, field values and field's derivations. *)
    d\[Phi]=cd\[Phi]+dd\[Phi]*dr;
    \[Phi]=c\[Phi]+cd\[Phi] dr;
    r=r+dr;
    (* Evolve the field and field's first derivation by dr. *)
    c\[Phi]=\[Phi];
    cd\[Phi]=d\[Phi];
    dd\[Phi]=AnalyticalRenormalization[c\[Phi]]*fvp[c\[Phi]]+1/2 fzp[c\[Phi]]/AnalyticalRenormalization[c\[Phi]] cd\[Phi]^2-(dimension-1)/r cd\[Phi];
    (* Then the evolved field and field's first derivation become the current ones. *)
    k++;
    
  ];
  If[Is1DTunneling==True,WriteString["stdout","\[SmallCircle]"]];
  
  Return[fieldinfo]
  ];


(* ::Subsubsection:: *)
(*Finding the Bounce solution and record the configuation*)


RecordBounce[Vi_,\[Phi]tv_,\[Phi]fv_,\[Phi]bar_,dr0_,tol_,\[Phi]0_,j_]:=Module[
{sign,\[Phi]info,\[Mu]info,\[Delta]info,data\[Phi]info,data\[Mu]info,data\[Delta]info,dr,d\[Phi],\[Phi],\[Mu],\[Delta],c\[Phi],c\[Mu],c\[Delta],cd\[Phi],d\[Mu],d\[Delta],dd\[Phi],k,\[Mu]r,\[Delta]r,\[Phi]r,dMr,fr},
sign=Sign[\[Phi]tv-\[Phi]fv];
\[Phi]info=Array[data\[Phi]info,{j,4}];
\[Mu]info=Array[data\[Mu]info,{j,2}];
\[Delta]info=Array[data\[Delta]info,{j,2}];
dr=dr0;
\[Mu]=0;
\[Delta]=0;
\[Phi]=\[Phi]0;
d\[Phi]=(GravityTunnelingrh Vi'[\[Phi]0] )/(1-8\[Pi] GravityTunnelingG (GravityTunnelingrh^2) Vi[\[Phi]0] );
k=0;
c\[Mu]=\[Mu];
c\[Delta]=\[Delta];
c\[Phi]=\[Phi];
cd\[Phi]=d\[Phi];
data\[Mu]info[1,1]=GravityTunnelingrh;
data\[Mu]info[1,2]=\[Mu];
data\[Delta]info[1,1]=GravityTunnelingrh;
data\[Delta]info[1,2]=\[Delta];
data\[Phi]info[1,1]=GravityTunnelingrh;
data\[Phi]info[1,2]=\[Phi];
d\[Mu]=4\[Pi] GravityTunnelingrh^2 Vi[c\[Phi]];
d\[Delta]=4\[Pi] GravityTunnelingG GravityTunnelingrh cd\[Phi]^2;
dd\[Phi]=0;
\[Mu]=c\[Mu]+d\[Mu] dr;
\[Delta]=c\[Delta]+d\[Delta] dr;
d\[Phi]=cd\[Phi]+dd\[Phi] dr;
\[Phi]=c\[Phi]+cd\[Phi] dr;
While[(sign*d\[Phi]<=0|| sign*(\[Phi]-\[Phi]bar)>=0) && sign*(\[Phi]-\[Phi]fv)>=0,
k++;
c\[Mu]=\[Mu];
c\[Delta]=\[Delta];
c\[Phi]=\[Phi];
cd\[Phi]=d\[Phi];
data\[Mu]info[k,1]=GravityTunnelingrh+k dr;
data\[Mu]info[k,2]=\[Mu];
data\[Delta]info[k,1]=GravityTunnelingrh+k dr;
data\[Delta]info[k,2]=\[Delta];
data\[Phi]info[k,1]=GravityTunnelingrh+k dr;
data\[Phi]info[k,2]=\[Phi];
data\[Phi]info[k,3]=d\[Phi];
d\[Mu]=4\[Pi] (GravityTunnelingrh+k dr)^2 (1/2 (1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+k dr))cd\[Phi]^2+Vi[c\[Phi]]);
d\[Delta]=4\[Pi] GravityTunnelingG(GravityTunnelingrh+k dr)cd\[Phi]^2;
dd\[Phi]=(Vi'[c\[Phi]]-((2GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu])-2GravityTunnelingG(GravityTunnelingrh+k dr)d\[Mu])/(GravityTunnelingrh+k dr)^2 cd\[Phi]+2/(GravityTunnelingrh+k dr) (1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+k dr))cd\[Phi]+d\[Delta](1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+k dr))cd\[Phi]))/(1-(2 GravityTunnelingG (GravityTunneling\[Mu]minus+c\[Mu]))/(GravityTunnelingrh+k dr));
data\[Phi]info[k,4]=dd\[Phi];
\[Mu]=c\[Mu]+d\[Mu] dr;
\[Delta]=c\[Delta]+d\[Delta] dr;
d\[Phi]=cd\[Phi]+dd\[Phi]*dr;
\[Phi]=c\[Phi]+cd\[Phi] dr;];
\[Mu]r=Table[{data\[Mu]info[i,1],data\[Mu]info[i,2]},{i,j-2}];
\[Delta]r=Table[{data\[Delta]info[i,1],data\[Delta]info[i,2]},{i,j-2}];
\[Phi]r=Table[{data\[Phi]info[i,1],data\[Phi]info[i,2]},{i,j-2}];
dMr=Table[{data\[Mu]info[i,1],data\[Mu]info[i,2]-(4\[Pi])/3 Vi[data\[Phi]info[i,2]]*data\[Mu]info[i,1]^3},{i,j-2}];
fr=Table[{data\[Mu]info[i,1],1-(2GravityTunnelingG (GravityTunneling\[Mu]minus+data\[Mu]info[i,2]))/data\[Mu]info[i,1]},{i,j-2}];
(*Print[ListPlot[\[Phi]r,(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["\[Phi]",10]}]];
P2=ListPlot[\[Mu]r,(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["\[Mu]",10]}];
P3=ListPlot[\[Delta]r,(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["\[Delta]",10]}];
(*The mass changing of the black hole is evaluate by the SdS metric by \[Mu] parameter, dM[r]=\[Mu][r]-(8\[Pi] G Vi[\[Phi][r]])/(6 G)(r^3).*)
P4=ListPlot[dMr,(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["dM",10]}];
P5=ListPlot[fr,(*PlotLegends\[Rule]Placed[{""},{0.2,0.7}],*)Frame->True,FrameLabel->{Style["r",10],Style["f",10]},PlotRange->All];
Print[GraphicsGrid[{{P2, P3}, {P4, P5}}]];*)
Return[{\[Phi]r,\[Mu]r,\[Delta]r,dMr,fr}];
]


(* ::Subsubsection:: *)
(*CalculateTheAction*)


CalculateTheAction[APotential_, ARenormalization_,NumberOfTheBigR_,dimension_,FieldAndDerivative_]:=
Sum[Extract[FieldAndDerivative,{q,5}]*((Extract[FieldAndDerivative,{q,1}])^(dimension-1))*(1/(2ARenormalization[Extract[FieldAndDerivative,{q,2}]])*(Extract[FieldAndDerivative,{q,3}])^2+APotential[Extract[FieldAndDerivative,{q,2}]]),{q,NumberOfTheBigR}]*(2\[Pi]^(dimension/2))/Gamma[dimension/2];


(* ::Section:: *)
(*Computing the Action*)


(* ::Subsection:: *)
(*Action in the Thin Wall limit*)


(*Tunneling action under thin-wall approximation, in this case the action can be carry out by equations B=(Subscript[A, fv]-Subscript[A, tv])/(4G),
where A is the suface area of the black hole in fv/tv and Subscript[A, fv/tv]=4Subscript[\[Pi]r^2, tv/fv]=4\[Pi](2Subscript[GM, tv/fv])^2=16\[Pi]G^2Subscript[M^2, tv/fv]. 
So, we have B=4\[Pi]G(Subscript[M, fv]+Subscript[M, tv])(Subscript[M, fv]-Subscript[M, tv])~8Subscript[\[Pi]GM, fv]dM*)
ThinWallAction[dM_]:=8\[Pi] GravityTunnelingG (GravityTunneling\[Mu]minus+dM/2)dM;


(* ::Subsection:: *)
(*Strictly calculated action quantity*)


(* ::Subsubsection:: *)
(*the GR tensor algebra *)


(*The Metric tensor*)
gT[r_,\[Theta]_]:=({
 {f[r]Exp[2 \[CapitalDelta][r]], 0, 0, 0},
 {0, 1/f[r], 0, 0},
 {0, 0, r^2, 0},
 {0, 0, 0, r^2 Sin[\[Theta]]^2}
});
cgT[r_,\[Theta]_]:=({
 {1/(f[r]Exp[2 \[CapitalDelta][r]]), 0, 0, 0},
 {0, f[r], 0, 0},
 {0, 0, 1/r^2, 0},
 {0, 0, 0, 1/(r^2 Sin[\[Theta]]^2) }
});
(*The Integraal Measure of four dimensional spacetimes*)
gme[r_,\[Theta]_]:=Det[gT[r,\[Theta]]];
hme[r_,\[Theta]_]:=Det[gT[r,\[Theta]]]/(f[r]Exp[2 \[CapitalDelta][r]]);
(*The Affine Connection*)
\[CapitalGamma]T[r_,\[Theta]_]:={{{0,f'[r]/(2f[r])+\[CapitalDelta]'[r],0,0},{f'[r]/(2f[r])+\[CapitalDelta]'[r],0,0,0},{0,0,0,0},{0,0,0,0}},
{{(-f[r]Exp[2\[CapitalDelta][r]])/2 (f'[r]+2\[CapitalDelta]'[r]f[r]),0,0,0},{0,-(f'[r]/(2f[r])),0,0},{0,0,-r f[r],0},{0,0,0,-r Sin[\[Theta]]^2 f[r]}},
{{0,0,0,0},{0,0,1/r,0},{0,1/r,0,0},{0,0,0,-Sin[\[Theta]]Cos[\[Theta]]}},
{{0,0,0,0},{0,0,0,1/r},{0,0,0,(1/Tan[\[Theta]])},{0,1/r,(1/Tan[\[Theta]]),0}}};
(*The Ricci Curvature Tensor*)
d\[CapitalGamma][r_,\[Theta]_]=({
 {-D[\[CapitalGamma]T[r,\[Theta]][[2,1,1]],r], -D[\[CapitalGamma]T[r,\[Theta]][[2,2,1]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,2,1]],\[Theta]], -D[\[CapitalGamma]T[r,\[Theta]][[2,3,1]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,3,1]],\[Theta]], -D[\[CapitalGamma]T[r,\[Theta]][[2,4,1]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,4,1]],\[Theta]]},
 {Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,1,i]],r],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,2,1]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,2,1]],\[Theta]], Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,2,i]],r],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,2,2]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,2,2]],\[Theta]], Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,3,i]],r],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,2,3]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,2,3]],\[Theta]], Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,4,i]],r],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,2,4]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,2,4]],\[Theta]]},
 {Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,1,i]],\[Theta]],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,3,1]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,3,1]],\[Theta]], Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,2,i]],\[Theta]],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,3,2]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,3,2]],\[Theta]], Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,3,i]],\[Theta]],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,3,3]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,3,3]],\[Theta]], Sum[D[\[CapitalGamma]T[r,\[Theta]][[i,4,i]],\[Theta]],{i,1,4}]-D[\[CapitalGamma]T[r,\[Theta]][[2,3,4]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,3,4]],\[Theta]]},
 {-D[\[CapitalGamma]T[r,\[Theta]][[2,4,1]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,4,1]],\[Theta]], -D[\[CapitalGamma]T[r,\[Theta]][[2,4,2]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,4,2]],\[Theta]], -D[\[CapitalGamma]T[r,\[Theta]][[2,4,3]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,4,3]],\[Theta]], -D[\[CapitalGamma]T[r,\[Theta]][[2,4,4]],r]-D[\[CapitalGamma]T[r,\[Theta]][[3,4,4]],\[Theta]]}
});
\[CapitalGamma]\[CapitalGamma][r_,\[Theta]_]:=({
 {Sum[\[CapitalGamma]T[r,\[Theta]][[i,1,j]]*\[CapitalGamma]T[r,\[Theta]][[j,1,i]]-\[CapitalGamma]T[r,\[Theta]][[i,1,1]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,1,j]]*\[CapitalGamma]T[r,\[Theta]][[j,2,i]]-\[CapitalGamma]T[r,\[Theta]][[i,1,2]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,1,j]]*\[CapitalGamma]T[r,\[Theta]][[j,3,i]]-\[CapitalGamma]T[r,\[Theta]][[i,1,3]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,1,j]]*\[CapitalGamma]T[r,\[Theta]][[j,4,i]]-\[CapitalGamma]T[r,\[Theta]][[i,1,4]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}]},
 {Sum[\[CapitalGamma]T[r,\[Theta]][[i,2,j]]*\[CapitalGamma]T[r,\[Theta]][[j,1,i]]-\[CapitalGamma]T[r,\[Theta]][[i,2,1]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,2,j]]*\[CapitalGamma]T[r,\[Theta]][[j,2,i]]-\[CapitalGamma]T[r,\[Theta]][[i,2,2]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,2,j]]*\[CapitalGamma]T[r,\[Theta]][[j,3,i]]-\[CapitalGamma]T[r,\[Theta]][[i,2,3]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,2,j]]*\[CapitalGamma]T[r,\[Theta]][[j,4,i]]-\[CapitalGamma]T[r,\[Theta]][[i,2,4]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}]},
 {Sum[\[CapitalGamma]T[r,\[Theta]][[i,3,j]]*\[CapitalGamma]T[r,\[Theta]][[j,1,i]]-\[CapitalGamma]T[r,\[Theta]][[i,3,1]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,3,j]]*\[CapitalGamma]T[r,\[Theta]][[j,2,i]]-\[CapitalGamma]T[r,\[Theta]][[i,3,2]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,3,j]]*\[CapitalGamma]T[r,\[Theta]][[j,3,i]]-\[CapitalGamma]T[r,\[Theta]][[i,3,3]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,3,j]]*\[CapitalGamma]T[r,\[Theta]][[j,4,i]]-\[CapitalGamma]T[r,\[Theta]][[i,3,4]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}]},
 {Sum[\[CapitalGamma]T[r,\[Theta]][[i,4,j]]*\[CapitalGamma]T[r,\[Theta]][[j,1,i]]-\[CapitalGamma]T[r,\[Theta]][[i,4,1]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,4,j]]*\[CapitalGamma]T[r,\[Theta]][[j,2,i]]-\[CapitalGamma]T[r,\[Theta]][[i,4,2]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,4,j]]*\[CapitalGamma]T[r,\[Theta]][[j,3,i]]-\[CapitalGamma]T[r,\[Theta]][[i,4,3]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}], Sum[\[CapitalGamma]T[r,\[Theta]][[i,4,j]]*\[CapitalGamma]T[r,\[Theta]][[j,4,i]]-\[CapitalGamma]T[r,\[Theta]][[i,4,4]]*\[CapitalGamma]T[r,\[Theta]][[j,i,j]],{i,1,4},{j,1,4}]}
});
RT[r_,\[Theta]_]:=d\[CapitalGamma][r,\[Theta]]+\[CapitalGamma]\[CapitalGamma][r,\[Theta]];
(*The Ricci Scalar*)
RS[r_,\[Theta]_]:=Sum[cgT[r,\[Theta]][[i,j]]*RT[r,\[Theta]][[i,j]],{i,1,4},{j,1,4}];


(* ::Section:: *)
(*End*)


End[]
EndPackage[]
