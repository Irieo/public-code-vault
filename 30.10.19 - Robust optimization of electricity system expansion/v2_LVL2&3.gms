* 6 node system
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
        mapSL (l,n) / l1.n1, l2.n1, l3.n4, l4.n2, l5.n2, l6.n3, l7.n3, l8.n4, l9.n5/
        mapRL (l,n) / l1.n2, l2.n3, l3.n5, l4.n3, l5.n4, l6.n4, l7.n6, l8.n6, l9.n6 /;

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
d1 200   300   40
d2 150   250   52
d3 100   200   55
d4 200   300   65;

TABLE GDATA (g ,*)
   PEmax Gcost
g1 300   18
g2 250   25
g3 400   16
g4 300   32
g5 150   35;

SCALARS
ILmax   /3000000/
SIGMA   /8760/
M       /5000000/
GAMMA   /0.5/;

VARIABLES
Z
PL(l)
THETA(n)

Lam(n)
phi_l(l)
phi_lplus(l)
xi_ref(n)
;

POSITIVE VARIABLES
PG (g)
PLS (d)
Pdem(d)

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
*x(l)
u1(l),u2(l),u3(g),u4(g),u5(d),u6(d),u7(n),u8(n);

parameter
x(l);

X(l)    = 0;
X('l5') = 1;
X('l7') = 1;

EQUATIONS
EQ_obj_sub

*EQ_invbudj
EQ_marketclear

EQ_flowdef_ex
EQ_flowdef_pros

*EQ_flowcon_a
*EQ_flowcon_b
*EQ_genCAP
*EQ_lsCAP
*EQ_vaCAPa
*EQ_vaCAPb

EQ_vaREF

EQ_ARO_demlo
EQ_ARO_demup
EQ_ARO_Ubudg

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

*                       SUM(l$pros(l), LDATA(l,"IC")*x(l))

EQ_obj_sub.. Z =E=      SIGMA*(SUM(g,GDATA (g,"Gcost")*PG(G))
                                 +SUM(d,DDATA (d,"LSCOST")* PLS(d)));

*EQ_invbudj.. SUM(l$pros(l), LDATA(L,"IC")*x(l)) =L= ILmax;

* EQUATIONS 3C ARE DEFINITION OF BINARY VARIABLES58 2 Transmission Expansion Planning

EQ_marketclear(n)..  SUM(g$mapG(g,n), PG(G))
                    -SUM(l$mapSL(l,n), PL(l))
                    +SUM(l$mapRL(l,n), PL(l)) =E=
                     SUM(d$mapD(d,n),  Pdem(d) - PLS (d));

EQ_flowdef_ex(l)$EX(l).. PL(l) =E= LDATA (l,"B")
                                *(SUM(n$mapSL(l,n), THETA (n)) - SUM (n$mapRL(l,n),THETA (n)));

EQ_flowdef_pros(l)$PROS(l).. PL(l) =E= x(l)*LDATA (l,"B")
                                *(SUM(n$mapSL(l,n), THETA (n)) - SUM (n$mapRL(l,n),THETA (n)));

*EQ_flowcon_a(l).. -LDATA (l, "FLmax ") =L= PL (l);
*EQ_flowcon_b(l)..  PL(l) =L= LDATA(l, "FLmax");

*EQ_genCAP(g)..           PG(g)  =L= GDATA(g, "PEmax");
*EQ_lsCAP(d)..            PLS(d) =L= Pdem(d);
*EQ_vaCAPa(n)..           -3.14  =L= THETA (n);
*EQ_vaCAPb(n)..           THETA(n) =L= 3.14;

EQ_vaREF(n)$REF(n)..     THETA(n) =L= 0;

EQ_ARO_demlo(d)..   Pdem(d) =G= DDATA(d, "PDmin");
EQ_ARO_demup(d)..   Pdem(d) =L= DDATA(d, "PDmax");
EQ_ARO_Ubudg..      SUM(d, Pdem(d) - DDATA(d,"PDmin")) / SUM(d, DDATA(d,"PDmax") - DDATA(d,"PDmin")) =E= GAMMA;

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
kkt3b(g)..  GDATA(g, "PEmax") - PG(g) =g= 0;
kkt3c(g)..  phi_Emax(g)               =l= M*(u3(g));
kkt3d(g)..  GDATA(g, "PEmax") - PG(g) =l= M*(1-u3(g));

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

MODEL ARO_LVL23 /ALL/;

SOLVE ARO_LVL23 USING MIP MINIMIZING Z;

execute_unload "all_results.gdx";
