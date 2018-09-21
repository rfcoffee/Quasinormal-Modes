(* Calculate quasinormal modes of AdS-Schwarzschild black holes. For details, see: http://arxiv.org/abs/arXiv:1709.01641 *)

(* We are using the de Sitter unit *)
seriesorder=50;
rhlist=Range[5,10,1];
numrh=Length[rhlist];
sollist=Range[1,numrh];
Thlist=Range[1,numrh];
For[i=1,i<=numrh,i++,

  ClearAll[r,z,m,c,x,n,k,w,glist];
  h[r_]:=1;
  F[r_]:=1+r^2-2*m/r;
  f[r_]:=F[r]/h[r];
  V[r_]:=f'[r]/r+c*h[r]/r^2;
  rh=rhlist[[i]];
  m=rh/2*(1+rh^2);
  c=0;

  xh=1/rh;
  kh=1/2*f'[rh];
  Th=kh/(2*Pi);
  chopnum=10^(-17);

  fnorm[x_]=Chop[Normal[Series[f[1/x],{x,xh,seriesorder}]],chopnum];
  flist=CoefficientList[fnorm[y+xh],y];

  fprimenorm[x_]=fnorm'[x];
  fprimelist=CoefficientList[fprimenorm[y+xh],y];

  (*glist=-Drop[flist,1] *)
  x4[x_]=Chop[Normal[Series[x^4,{x,xh,seriesorder}]],chopnum];
  x4list=CoefficientList[x4[y+xh],y];
  x3[x_]=Chop[Normal[Series[x^3,{x,xh,seriesorder}]],chopnum];
  x3list=CoefficientList[x3[y+xh],y];
  x2[x_]=Chop[Normal[Series[x^2,{x,xh,seriesorder}]],chopnum];
  x2list=CoefficientList[x2[y+xh],y];
  x1[x_]=Chop[Normal[Series[x^1,{x,xh,seriesorder}]],chopnum];
  x1list=CoefficientList[x1[y+xh],y];
  hnorm[x_]=Chop[Normal[Series[h[1/x],{x,xh,seriesorder}]],chopnum];
  hlist=CoefficientList[hnorm[y+xh],y];

  (* redo the s,t,u lists for AdS BH*)
  s1[y_]=Simplify[-(y+xh)^4*1/y*f[1/(y+xh)]];
  snorm[y_]=Chop[Normal[Series[s1[y],{y,0,seriesorder}]]];
  t1[y_]=Simplify[-2*(y+xh)^3*f[1/(y+xh)]+(y+xh)^2*f'[1/(y+xh)]-2*I*w*(y+xh)^2];
  tnorm[y_]=Chop[Normal[Series[t1[y],{y,0,seriesorder}]]];
  u1[y_]=Simplify[y*(y+xh)*f'[1/(y+xh)]];
  unorm[y_]=Chop[Normal[Series[u1[y]+c*y*(y+xh)^2*h[1/(y+xh)],{y,0,seriesorder}]]];
  slist1=CoefficientList[snorm[y],y];
  slist=If[Length[slist1]>seriesorder,slist1[[1;;seriesorder]],Join[slist1,Table[0,{i,seriesorder-Length[slist1]}]]];
  tlist1=CoefficientList[tnorm[y],y];
  tlist=If[Length[tlist1]>seriesorder,tlist1[[1;;seriesorder]],Join[tlist1,Table[0,{i,seriesorder-Length[tlist1]}]]];
  ulist1=CoefficientList[unorm[y],y];
  ulist=If[Length[ulist1]>seriesorder,ulist1[[1;;seriesorder]],Join[ulist1,Table[0,{i,seriesorder-Length[ulist1]}]]];

  Slist=(-xh)^Range[0,seriesorder-1]*slist;
  Tlist=(-xh)^Range[0,seriesorder-1]*tlist;
  Ulist=(-xh)^Range[0,seriesorder-1]*ulist;
  P[n_]:=n*(n-1)*slist[[1]]+n*tlist[[1]];
  Alist=Range[1,seriesorder];

  For[n=1,n<=(seriesorder-1),n++,
    Alist[[n+1]]=Factor[-1/P[n]*Sum[(k*(k-1)*Slist[[n-k+1]]+k*Tlist[[n-k+1]]+Ulist[[n-k+1]])*Alist[[k+1]],{k,0,n-1}]];];

  psi0[w_]=Numerator[Factor[Sum[Alist[[n+1]],{n,0,seriesorder-1}]]];
  sollist[[i]]=w/.(NSolve[{psi0[w]==0,Re[w]>0,Abs[Im[w]]>chopnum},w]);
  Thlist[[i]]=Th;];

(* Export to table: Col1: BH radius, Col2: Hawking temperature: Col3-: Quasi-normal modes at increasing orders*)
Export["out_SAdS.dat",Join[Partition[rhlist,1],Join[Partition[Thlist,1],sollist,2],2],"Table"];
Exit[];
