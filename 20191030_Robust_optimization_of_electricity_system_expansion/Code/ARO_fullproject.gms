$ontext

The code for the presentation "Robust optimization of electricity system expansion" 
BTU C-S, 30.10.2019. 
Iegor Riepin

A 6-node power system and RO math I am using here is based on a great book:
Conejo, A. J., Baringo, M. L., Kazempour, S. J., & Siddiqui, A. S. (2016). 
Investment in Electricity Generation and Transmission: Decision Making under Uncertainty.
DOI: 10.1007/978-3-319-29501-5

The formulation of a column-and-constraint generation (or Berders-primal) algorithm is based on 
multiple sources. To my knowledge, there are no publicly available GAMS formulations of this monster,
thus I hope my implementation might be of use.

Feedback, bug reportings and suggestions are highly welcome: iegor.riepin@b-tu.de

$offtext

SETS
  n / n1 * n6 /
  g / g1 * g5 /
  d / d1 * d4 /
  l / l1 * l9 /
  pros(l) / l4 * l9 /
  ex(l)   / l1 * l3 /
  mapG(g,n) / g1.n1, g2.n2, g3.n3, g4.n5, g5.n6 /
  mapD(d,n) / d1.n1, d2.n4, d3.n5, d4.n6 /
  ref (n) / n1 /
  mapSL (l,n) / l1.n1, l2.n1, l3.n4, l4.n2, l5.n2, l6.n3, l7.n3, l8.n4, l9.n5 /
  mapRL (l,n) / l1.n2, l2.n3, l3.n5, l4.n3, l5.n4, l6.n4, l7.n6, l8.n6, l9.n6 /

  v iterations /v1 * v20/
  vv(v) iterations subset

  i dsens /i1*i5/
  j pesens /j1*j5/
  ;

TABLE LDATA (L ,*)
         B   FLmax IC
l1       500 150   0
l2       500 150   0
l3       500 150   0
l4       500 150   700000
l5       500 150   1400000
l6       500 200   1800000
l7       500 200   1600000
l8       500 150   800000
l9       500 150   700000;

TABLE DDATA (d ,*)
   PDmin PDmax LScost
d1 180   220   140
d2 135   165   152
d3 90    110   155
d4 180   220   165;

TABLE GDATA (g ,*)
   PEmin PEmax Gcost
g1 0     300   18
g2 0     250   25
g3 0     400   16
g4 0     300   32
g5 0     150   35;

SCALARS
ILmax     /3000000/
SIGMA     /8760/
M         /5000000/
GAMMA_D   /0/
GAMMA_G   /0/
;


VARIABLES
**************MASTER*************
Z_M
ETA

PL_M(l,v)
THETA_M(n,v)

**************SUB**************
Z
PL(l)
THETA(n)

Lam(n)
phi_l(l)
phi_lplus(l)
xi_ref(n)
;

POSITIVE VARIABLES
**************MASTER*************
PG_M(g,v)
PLS_m(d,v)

**************SUB**************
PG (g)
PLS (d)
Pdem(d)
PE(g)

phi_lmin(l)
phi_lmax(l)
phi_Emin(g)
phi_Emax(g)
phi_Dmin(d)
phi_Dmax(d)
phi_Nmin(n)
phi_Nmax(n)
;

BINARY VARIABLES
**************MASTER*************
x_m(l)


**************SUB**************
u1(l),u2(l),u3(g),u4(g),u5(d),u6(d),u7(n),u8(n);

parameters
Pdem_fixed(d,v)
PE_fixed(g,v)
x(l);

**************MASTER*************
*Pdem_fixed(d,vv) = 0;

**************SUB**************
*X(l)    = 0;
*X('l5') = 1;
*X('l7') = 1;

parameters
Tol         decomposition tolerance / 1e-6/
Z_Lower     lower bound / -inf /
Z_Upper     upper bound / inf /
;

EQUATIONS
**************MASTER*************
EQ_obj_master

EQ_invbudj
EQ_marketclear_m(n,v)

EQ_flowdef(l,v)
EQ_flow_cona(l,v)
EQ_flow_conb(l,v)
EQ_flowMILP_cona(l,v)
EQ_flowMILP_conb(l,v)
EQ_flowdefMILPa(l,v)
EQ_flowdefMILPb(l,v)

EQ_genCAP(g,v)
EQ_lsCAP(d,v)
EQ_vaCAPa(n,v)
EQ_vaCAPb(n,v)
EQ_vaREF_m(n,v)

EQ_ETAdef(v)

**************SUB**************
EQ_obj_sub

EQ_marketclear

EQ_flowdef_ex
EQ_flowdef_pros
EQ_vaREF

EQ_ARO_demlo
EQ_ARO_demup
EQ_ARO_Ubudg

EQ_ARO_gavlo
EQ_ARO_gavup
EQ_ARO_UbudgG

FOC_PG
FOC_PLS
FOC_PL_EX
FOC_PL_PROS
FOC_THETA
FOC_THETA_REF

kkt1a, kkt1b, kkt1c, kkt1d
kkt2a, kkt2b, kkt2c, kkt2d
kkt3a, kkt3b, kkt3c, kkt3d
kkt4a, kkt4b, kkt4c, kkt4d
kkt5a, kkt5b, kkt5c, kkt5d
kkt6a, kkt6b, kkt6c, kkt6d
kkt7a, kkt7b, kkt7c, kkt7d
kkt8a, kkt8b, kkt8c, kkt8d
;

**************MASTER*************

EQ_obj_master.. Z_M  =E= SUM(l$pros(l), LDATA(l,"IC")*x_m(l)) + ETA;


EQ_invbudj.. SUM(l$pros(l), LDATA(L,"IC")*x_m(l)) =L= ILmax;

EQ_marketclear_m(n,vv)..  SUM(g$mapG(g,n), PG_M(g,vv))
                        -SUM(l$mapSL(l,n), PL_M(l,vv))
                        +SUM(l$mapRL(l,n), PL_M(l,vv)) =E=
                         SUM(d$mapD(d,n),  Pdem_fixed(d,vv) - PLS_m(d,vv));

EQ_flowdef(l,vv)$EX(l).. PL_M(l,vv) =E= LDATA (l,"B")
                                *(SUM(n$mapSL(l,n), THETA_M(n,vv)) - SUM (n$mapRL(l,n), THETA_M(n,vv)));

EQ_flow_cona(l,vv)$EX(l).. -LDATA (l, "FLmax ") =L= PL_M(l,vv);
EQ_flow_conb(l,vv)$EX(l)..  PL_M(l,vv) =L= LDATA(l, "FLmax");

EQ_flowMILP_cona(l,vv)$PROS(l).. -x_m(l)*LDATA(l, "FLmax") =L= PL_M(l,vv);

EQ_flowMILP_conb(l,vv)$PROS(l).. PL_M(l,vv) =L= x_m(l) * LDATA(l, "FLmax");

EQ_flowdefMILPa(l,vv)$PROS(l).. -(1-x_m(l))*M =L= PL_M(l,vv) - LDATA(l, "B")*(SUM(n$mapSL(l,n), THETA_M(n,vv))
                                                                        - SUM(n$mapRL(l,n), THETA_M(n,vv)));

EQ_flowdefMILPb(l,vv)$PROS(l).. PL_M(l,vv) - LDATA(l, "B") * (SUM(n$mapSL (l,n), THETA_M(n,vv)) - SUM(N$mapRL(l,n), THETA_M(n,vv))) =l=
                                   (1 - x_m(l))*M;

EQ_genCAP(g,vv)..           PG_M(g,vv)    =L= PE_fixed(g,vv);
EQ_lsCAP(d,vv)..            PLS_M(d,vv)   =L= Pdem_fixed(d,vv);
EQ_vaCAPa(n,vv)..           -3.14         =L= THETA_M(n,vv);
EQ_vaCAPb(n,vv)..           THETA_M(n,vv) =L= 3.14;
EQ_vaREF_m(n,vv)$REF(n)..   THETA_M(n,vv) =L= 0;

EQ_ETAdef(vv)..  ETA =G= SIGMA*(SUM(g,GDATA (g,"Gcost")*PG_M(g,vv))
                                 +SUM(d,DDATA (d,"LSCOST")*PLS_M(d,vv)));

**************SUB**************

EQ_obj_sub.. Z =E=      SIGMA*(SUM(g,GDATA (g,"Gcost")*PG(G))
                                 +SUM(d,DDATA (d,"LSCOST")* PLS(d)));

* EQUATIONS 3C ARE DEFINITION OF BINARY VARIABLES58 2 Transmission Expansion Planning

EQ_marketclear(n)..  SUM(g$mapG(g,n), PG(G))
                    -SUM(l$mapSL(l,n), PL(l))
                    +SUM(l$mapRL(l,n), PL(l)) =E=
                     SUM(d$mapD(d,n),  Pdem(d) - PLS (d));

EQ_flowdef_ex(l)$EX(l).. PL(l) =E= LDATA (l,"B")
                                *(SUM(n$mapSL(l,n), THETA (n)) - SUM (n$mapRL(l,n),THETA (n)));

EQ_flowdef_pros(l)$PROS(l).. PL(l) =E= x(l)*LDATA (l,"B")
                                *(SUM(n$mapSL(l,n), THETA (n)) - SUM (n$mapRL(l,n),THETA (n)));

EQ_vaREF(n)$REF(n)..     THETA(n) =L= 0;

EQ_ARO_demlo(d)..   Pdem(d) =G= DDATA(d, "PDmin");
EQ_ARO_demup(d)..   Pdem(d) =L= DDATA(d, "PDmax");
EQ_ARO_Ubudg..      SUM(d, Pdem(d) - DDATA(d,"PDmin")) / SUM(d, DDATA(d,"PDmax") - DDATA(d,"PDmin")) =L= GAMMA_D;

EQ_ARO_gavlo(g)..   PE(g) =G= GDATA(g, "PEmin");
EQ_ARO_gavup(g)..   PE(g) =L= GDATA(g, "PEmax");
EQ_ARO_UbudgG..     SUM(g, GDATA(g,"PEmax") - PE(g)) / SUM(g, GDATA(g,"PEmax")) =L= GAMMA_G;

FOC_PG(g)..         SIGMA*GDATA(g,"Gcost") -  sum(n$mapG(g,n),lam(n)) + phi_Emax(g) - phi_Emin(g)  =e= 0;
FOC_PLS(d)..        SIGMA*DDATA(d,"LSCOST") - sum(n$mapD(d,n),lam(n)) + phi_Dmax(d) - phi_Dmin(d) =e= 0;

FOC_PL_EX(l)$EX(l).. sum(n$mapSL(l,n),lam(n))
                   - sum(n$mapRL(l,n),lam(n))
                   - phi_l(l)
                   + phi_lmax(l)
                   - phi_lmin(l)
                   =e= 0;

FOC_PL_PROS(L)$PROS(l).. sum(n$mapSL(l,n),lam(n))
                       - sum(n$mapRL(l,n),lam(n))
                       - phi_lplus(l)
                       + phi_lmax(l)
                       - phi_lmin(l)
                       =e= 0;

FOC_THETA(n)$(not REF(n)).. sum(l$(EX(l) and mapSL (l,n)), LDATA(l,"B")*phi_l(l))
                           + sum(l$(pros(l) and mapSL (l,n)), x(l)*LDATA(l,"B")*phi_lplus(l))
                           - sum(l$(EX(l) and mapRL (l,n)), LDATA(l,"B")*phi_l(l))
                           - sum(l$(pros(l) and mapRL (l,n)), x(l)*LDATA(l,"B")*phi_lplus(l))
                           + phi_Nmax(n) - phi_Nmin(n) =e= 0;

FOC_THETA_REF(n)$REF(n).. sum(l$(EX(l) and mapSL (l,n)), LDATA(l,"B")*phi_l(l))
                           + sum(l$(pros(l) and mapSL (l,n)), x(l)*LDATA(l,"B")*phi_lplus(l))
                           - sum(l$(EX(l) and mapRL (l,n)), LDATA(l,"B")*phi_l(l))
                           - sum(l$(pros(l) and mapRL (l,n)), x(l)*LDATA(l,"B")*phi_lplus(l))
                          + phi_Nmax(n) - phi_Nmin(n) - xi_ref(n) =e= 0;

kkt1a(l)..  phi_lmax(l)               =g= 0;
kkt1b(l)..  LDATA(l, "FLmax") - PL(l) =g= 0;
kkt1c(l)..  phi_lmax(l)               =l= M*(u1(l));
kkt1d(l)..  LDATA(l, "FLmax") - PL(l) =l= M*(1-u1(l));

kkt2a(l)..  phi_lmin(l)               =g= 0;
kkt2b(l)..  PL(l) + LDATA(l, "FLmax") =g= 0;
kkt2c(l)..  phi_lmin(l)               =l= M*(u2(l));
kkt2d(l)..  PL(l) + LDATA(l, "FLmax") =l= M*(1-u2(l));

kkt3a(g)..  phi_Emax(g)               =g= 0;
kkt3b(g)..  PE(g) - PG(g)             =g= 0;
kkt3c(g)..  phi_Emax(g)               =l= M*(u3(g));
kkt3d(g)..  PE(G) - PG(g)             =l= M*(1-u3(g));

kkt4a(g)..  phi_Emin(g)               =g= 0;
kkt4b(g)..  PG(g)                     =g= 0;
kkt4c(g)..  phi_Emin(g)               =l= M*(u4(g));
kkt4d(g)..  PG(g)                     =l= M*(1-u4(g));

kkt5a(d)..  phi_Dmax(d)               =g= 0;
kkt5b(d)..  Pdem(d) - PLS(d)          =g= 0;
kkt5c(d)..  phi_Dmax(d)               =l= M*(u5(d));
kkt5d(d)..  Pdem(d) - PLS(d)          =l= M*(1-u5(d));

kkt6a(d)..  phi_Dmin(d)               =g= 0;
kkt6b(d)..  PLS(d)                    =g= 0;
kkt6c(d)..  phi_Dmin(d)               =l= M*(u6(d));
kkt6d(d)..  PLS(d)                    =l= M*(1-u6(d));

kkt7a(n)..  phi_Nmax(n)               =g= 0;
kkt7b(n)..  3.14 - THETA(n)           =g= 0;
kkt7c(n)..  phi_Nmax(n)               =l= M*(u7(n));
kkt7d(n)..  3.14 - THETA(n)           =l= M*(1-u7(n));

kkt8a(n)..  phi_Nmin(n)               =g= 0;
kkt8b(n)..  THETA(n) + 3.14           =g= 0;
kkt8c(n)..  phi_Nmin(n)               =l= M*(u8(n));
kkt8d(n)..  THETA(n) + 3.14           =l= M*(1-u8(n));


MODEL ARO_MASTER
/
EQ_obj_master
EQ_invbudj
EQ_marketclear_m
EQ_flowdef
EQ_flow_cona
EQ_flow_conb
EQ_flowMILP_cona
EQ_flowMILP_conb
EQ_flowdefMILPa
EQ_flowdefMILPb
EQ_genCAP
EQ_lsCAP
EQ_vaCAPa
EQ_vaCAPb
EQ_vaREF_m
EQ_ETAdef
/;

MODEL ARO_MASTER_V0
/
EQ_obj_master
EQ_invbudj
/;

MODEL ARO_LVL23
/
EQ_obj_sub
EQ_marketclear

EQ_flowdef_ex
EQ_flowdef_pros
EQ_vaREF

EQ_ARO_demlo
EQ_ARO_demup
EQ_ARO_Ubudg

EQ_ARO_gavlo
EQ_ARO_gavup
EQ_ARO_UbudgG

FOC_PG
FOC_PLS
FOC_PL_EX
FOC_PL_PROS
FOC_THETA
FOC_THETA_REF

kkt1a, kkt1b, kkt1c, kkt1d
kkt2a, kkt2b, kkt2c, kkt2d
kkt3a, kkt3b, kkt3c, kkt3d
kkt4a, kkt4b, kkt4c, kkt4d
kkt5a, kkt5b, kkt5c, kkt5d
kkt6a, kkt6b, kkt6c, kkt6d
kkt7a, kkt7b, kkt7c, kkt7d
kkt8a, kkt8b, kkt8c, kkt8d
/;

* to allow CPLEX correctly detect rays in an infeasible problem
* only simplex method can be used and no preprocessing neither scaling options
* optimality and feasibility tolerances are very small to avoid primal degeneration
file COPT / cplex.opt /
put COPT putclose 'ScaInd -1' / 'LPMethod 1' / 'PreInd 0' / 'EpOpt 1e-9' / 'EpRHS 1e-9' / ;
ARO_MASTER.OptFile = 1;
ARO_LVL23.OptFile = 1;

*SOLVE ARO_MASTER USING MIP MINIMIZING Z_M;
*execute_unload "MASTER_results.gdx";

*SOLVE ARO_LVL23 USING MIP MAXIMIZING Z;
*execute_unload "SUB_results.gdx";

parameter
x_ithist(l,v)
report(j,i,*,*)
;

GAMMA_D = 0;

loop(j,

GAMMA_D = 0;
vv(v)   = no;
Pdem_fixed(d,vv) = DDATA(d, "PDmin");
PE_fixed(g,vv)   = GDATA(g, "PEmax");
Z_Lower = -inf;
Z_Upper = inf;

loop(i,

vv(v)   = no;
Pdem_fixed(d,vv) = DDATA(d, "PDmin");
PE_fixed(g,vv)   = GDATA(g, "PEmax");
Z_Lower = -inf;
Z_Upper = inf;

loop (v $((abs(1-Z_Lower/Z_Upper) > Tol)),

  if (ord(v)=1,
    SOLVE ARO_MASTER_V0 USING MIP MINIMIZING Z_M;
  else
    SOLVE ARO_MASTER USING MIP MINIMIZING Z_M;
    );

  Z_Lower = SUM(l$pros(l), LDATA(l,"IC")*x_m.l(l)) + ETA.l;
  X(l) = x_m.l(l);

*  execute_unload "MASTER_results_check.gdx";


  SOLVE ARO_LVL23 USING MIP MAXIMIZING Z;

  Z_Upper = min(Z_Upper, (SUM(l$pros(l), LDATA(l,"IC")*x_m.l(l)) + Z.l) );
  vv(v) = yes ;
  Pdem_fixed(d,vv)$(vv.last) = Pdem.l(d);
  PE_fixed(g,vv)$(vv.last) = PE.l(g);
  x_ithist(l,vv)$(vv.last) = x_m.l(l);

*  execute_unload "SUB_results_check.gdx";
  );

  report(j,i,'X(l)',l)   =X(l);
  report(j,i,'G_D','')   =GAMMA_D;
  report(j,i,'G_G','')   =GAMMA_G;

  GAMMA_D=GAMMA_D+0.2;
);


put_utility 'gdxout' / 'J-steps' j.tl:0;
execute_unload GAMMA_D, GAMMA_G, report;

  GAMMA_G=GAMMA_G+0.2;
);

execute_unload "ARO_SENS_results.gdx" report;

